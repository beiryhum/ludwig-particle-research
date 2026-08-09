/*****************************************************************************
 *
 *  bbl_link_diag.c
 *
 *  Read-only exact moving-BBL link and implicit-force diagnostics.
 *
 *  The module deliberately separates:
 *
 *    A. Physical LINK_FLUID pass-2 exchange used for local surface loading.
 *    B. Ludwig's original implicit fhydro/Thydro assembled from pass-1 f0/t0
 *       and the zeta velocity correction.
 *
 *  No LB distribution, particle force, torque, velocity, angular velocity,
 *  position, drag matrix, or restart state is changed here.
 *
 *****************************************************************************/

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bbl_link_diag.h"

#define BBL_LINK_DIAG_DEFAULT_FREQ 1
#define BBL_LINK_DIAG_DEFAULT_TOLERANCE 1.0e-8
#define BBL_LINK_DIAG_INITIAL_CAPACITY 1024

/* Six-vector layout: force xyz followed by torque xyz. */
enum bbl_link_diag_six_enum {
  BBL_DIAG_FX = 0,
  BBL_DIAG_FY,
  BBL_DIAG_FZ,
  BBL_DIAG_TX,
  BBL_DIAG_TY,
  BBL_DIAG_TZ,
  BBL_DIAG_N6
};

enum bbl_link_diag_action_enum {
  BBL_LINK_DIAG_ACTION_OFF = 0,
  BBL_LINK_DIAG_ACTION_WARN,
  BBL_LINK_DIAG_ACTION_FATAL
};

struct bbl_link_diag_s {
  pe_t * pe;
  cs_t * cs;
  MPI_Comm comm;
  int rank;
  int size;

  int target_index;
  int frequency;
  int output_start;
  int output_end;
  int current_step;
  int active_step;
  int action;
  double tolerance;

  char link_filename[FILENAME_MAX];
  char closure_filename[FILENAME_MAX];

  FILE * link_fp;
  FILE * closure_fp;

  bbl_link_record_t * records;
  int nrecord;
  int capacity;

  /* Independent source accumulators, local to each MPI task. */
  double build_source[BBL_DIAG_N6];
  double pass1_fluid[BBL_DIAG_N6];
  double pass1_colloid[BBL_DIAG_N6];
  double pass1_boundary[BBL_DIAG_N6];
  double pass1_active[BBL_DIAG_N6];
  double pass1_total[BBL_DIAG_N6];
  int pass1_total_owner_count;
  double pass2_reference[BBL_DIAG_N6];
};

/* build.c source hooks use this process-local registration. Ludwig creates one
 * bbl object per process, so a single diagnostic registration is sufficient. */
static bbl_link_diag_t * bbl_link_diag_registered_ = NULL;

static int bbl_link_diag_record_compare(const void * a, const void * b);
static int bbl_link_diag_same_key(const bbl_link_record_t * a,
                                  const bbl_link_record_t * b);
static double bbl_link_diag_norm3(const double a[3]);
static int bbl_link_diag_finite_record(const bbl_link_record_t * r);
static int bbl_link_diag_parse_action(pe_t * pe, rt_t * rt);
static int bbl_link_diag_switch_default(rt_t * rt, const char * key,
                                        int default_value);
static int bbl_link_diag_site_is_local_interior(cs_t * cs, int index);
static void bbl_link_diag_six_zero(double a[BBL_DIAG_N6]);
static void bbl_link_diag_six_add(double dst[BBL_DIAG_N6],
                                  const double force[3],
                                  const double torque[3]);
static void bbl_link_diag_zeta_correction(const colloid_t * pc,
                                          double force[3],
                                          double torque[3]);
static double bbl_link_diag_force_scaled_error(const double a[3],
                                               const double b[3],
                                               double * abs_error);
static double bbl_link_diag_torque_scaled_error(const double a[3],
                                                const double b[3],
                                                double radius,
                                                double force_scale,
                                                double * abs_error);
static void bbl_link_diag_write_link_header(FILE * fp);
static void bbl_link_diag_write_closure_header(FILE * fp);
static void bbl_link_diag_write_record(FILE * fp,
                                       const bbl_link_record_t * r);
static void bbl_link_diag_write_vec6(FILE * fp,
                                     const double a[BBL_DIAG_N6]);

/*****************************************************************************
 *
 *  bbl_link_diag_create
 *
 *****************************************************************************/

int bbl_link_diag_create(pe_t * pe, cs_t * cs, rt_t * rt,
                         colloids_info_t * cinfo, int nstart,
                         bbl_link_diag_t ** pdiag) {

  bbl_link_diag_t * diag = NULL;
  colloid_t * pc = NULL;
  int enabled = 0;
  int local_count = 0;
  int global_count = 0;
  int local_bbl = 0;
  int global_bbl = 0;
  int owner_filter = 1;
  int include_halo = 0;
  int include_wall = 0;
  int n = 0;

  assert(pe);
  assert(cs);
  assert(rt);
  assert(cinfo);
  assert(pdiag);

  *pdiag = NULL;

  if (rt_key_present(rt, "colloid_bbl_link_output")) {
    enabled = rt_switch(rt, "colloid_bbl_link_output");
  }
  if (enabled == 0) return 0;

  if (bbl_link_diag_registered_ != NULL) {
    pe_fatal(pe, "Only one exact BBL link diagnostic may be registered\n");
  }

  diag = (bbl_link_diag_t *) calloc(1, sizeof(bbl_link_diag_t));
  if (diag == NULL) pe_fatal(pe, "calloc(bbl_link_diag_t) failed\n");

  diag->pe = pe;
  diag->cs = cs;
  cs_cart_comm(cs, &diag->comm);
  MPI_Comm_rank(diag->comm, &diag->rank);
  MPI_Comm_size(diag->comm, &diag->size);

  diag->target_index = 1;
  diag->frequency = BBL_LINK_DIAG_DEFAULT_FREQ;
  diag->output_start = nstart + 1;
  diag->output_end = nstart + 1;
  diag->current_step = nstart;
  diag->action = bbl_link_diag_parse_action(pe, rt);
  diag->tolerance = BBL_LINK_DIAG_DEFAULT_TOLERANCE;

  snprintf(diag->link_filename, FILENAME_MAX, "%s", "bbl_link_exact.csv");
  snprintf(diag->closure_filename, FILENAME_MAX, "%s",
           "bbl_link_closure.csv");

  rt_int_parameter(rt, "colloid_bbl_link_index", &diag->target_index);
  rt_int_parameter(rt, "colloid_bbl_link_output_freq", &diag->frequency);
  rt_int_parameter(rt, "colloid_bbl_link_output_start", &diag->output_start);
  rt_int_parameter(rt, "colloid_bbl_link_output_end", &diag->output_end);
  rt_double_parameter(rt, "colloid_bbl_link_closure_tolerance",
                      &diag->tolerance);

  n = rt_string_parameter(rt, "colloid_bbl_link_output_filename",
                          diag->link_filename, FILENAME_MAX);
  if (n == 0) {
    snprintf(diag->link_filename, FILENAME_MAX, "%s", "bbl_link_exact.csv");
  }

  n = rt_string_parameter(rt, "colloid_bbl_link_closure_filename",
                          diag->closure_filename, FILENAME_MAX);
  if (n == 0) {
    snprintf(diag->closure_filename, FILENAME_MAX, "%s",
             "bbl_link_closure.csv");
  }

  owner_filter = bbl_link_diag_switch_default(
      rt, "colloid_bbl_link_owner_filter", 1);
  include_halo = bbl_link_diag_switch_default(
      rt, "colloid_bbl_link_include_halo", 0);
  include_wall = bbl_link_diag_switch_default(
      rt, "colloid_bbl_link_include_wall", 0);

  if (diag->target_index <= 0) {
    pe_fatal(pe, "colloid_bbl_link_index must be positive\n");
  }
  if (diag->frequency <= 0) {
    pe_fatal(pe, "colloid_bbl_link_output_freq must be positive\n");
  }
  if (diag->output_start <= nstart) {
    pe_fatal(pe, "colloid_bbl_link_output_start must be greater than N_start\n");
  }
  if (diag->output_end < diag->output_start) {
    pe_fatal(pe, "colloid_bbl_link_output_end must be >= output_start\n");
  }
  if (!(diag->tolerance > 0.0) || !isfinite(diag->tolerance)) {
    pe_fatal(pe, "colloid_bbl_link_closure_tolerance must be finite and > 0\n");
  }
  if (owner_filter == 0) {
    pe_fatal(pe, "Exact BBL diagnostic requires owner_filter yes\n");
  }
  if (include_halo != 0) {
    pe_fatal(pe, "Exact BBL diagnostic requires include_halo no\n");
  }
  if (include_wall != 0) {
    pe_fatal(pe, "Particle physical-link CSV structurally excludes wall links\n");
  }

  colloids_info_local_head(cinfo, &pc);
  for (; pc; pc = pc->nextlocal) {
    if (pc->s.index == diag->target_index) {
      local_count += 1;
      if (pc->s.bc == COLLOID_BC_BBL) local_bbl += 1;
    }
  }

  MPI_Allreduce(&local_count, &global_count, 1, MPI_INT, MPI_SUM, diag->comm);
  MPI_Allreduce(&local_bbl, &global_bbl, 1, MPI_INT, MPI_SUM, diag->comm);

  if (global_count != 1) {
    pe_fatal(pe, "Exact BBL diagnostic target index %d has owner count %d\n",
             diag->target_index, global_count);
  }
  if (global_bbl != 1) {
    pe_fatal(pe, "Exact BBL diagnostic target index %d is not BBL\n",
             diag->target_index);
  }

  if (diag->rank == 0) {
    diag->link_fp = fopen(diag->link_filename, "w");
    if (diag->link_fp == NULL) {
      pe_fatal(pe, "Cannot open exact BBL link file %s\n",
               diag->link_filename);
    }
    diag->closure_fp = fopen(diag->closure_filename, "w");
    if (diag->closure_fp == NULL) {
      pe_fatal(pe, "Cannot open exact BBL closure file %s\n",
               diag->closure_filename);
    }
    bbl_link_diag_write_link_header(diag->link_fp);
    bbl_link_diag_write_closure_header(diag->closure_fp);
    fflush(diag->link_fp);
    fflush(diag->closure_fp);
  }

  bbl_link_diag_registered_ = diag;

  pe_info(pe, "\nExact moving-BBL two-layer diagnostic (v2)\n");
  pe_info(pe, "------------------------------------------\n");
  pe_info(pe, "Target colloid index:       %d\n", diag->target_index);
  pe_info(pe, "Output step range:          %d to %d\n",
          diag->output_start, diag->output_end);
  pe_info(pe, "Output frequency:           %d\n", diag->frequency);
  pe_info(pe, "Link output:                %s\n", diag->link_filename);
  pe_info(pe, "Closure output:             %s\n", diag->closure_filename);
  pe_info(pe, "Closure tolerance:          %.7e\n", diag->tolerance);
  pe_info(pe, "Physical force convention:  fluid acting on colloid\n");
  pe_info(pe, "Torque reference point:     current colloid centre\n");
  pe_info(pe, "Physical link semantics:    current pass2 BBL step\n");
  pe_info(pe, "Implicit force semantics:   pass1 f0/t0 minus zeta*[v,w]\n");
  pe_info(pe, "Owner rule:                 outside lattice site is local interior\n");
  pe_info(pe, "Wall links in physical CSV: structurally excluded\n");

  *pdiag = diag;

  return 0;
}

/*****************************************************************************
 *
 *  bbl_link_diag_free
 *
 *****************************************************************************/

int bbl_link_diag_free(bbl_link_diag_t ** pdiag) {

  bbl_link_diag_t * diag = NULL;

  if (pdiag == NULL || *pdiag == NULL) return 0;

  diag = *pdiag;

  if (diag->rank == 0) {
    if (diag->link_fp) {
      fflush(diag->link_fp);
      fclose(diag->link_fp);
    }
    if (diag->closure_fp) {
      fflush(diag->closure_fp);
      fclose(diag->closure_fp);
    }
  }

  if (bbl_link_diag_registered_ == diag) bbl_link_diag_registered_ = NULL;

  free(diag->records);
  free(diag);
  *pdiag = NULL;

  return 0;
}

/*****************************************************************************
 *
 *  bbl_link_diag_step_set
 *
 *  This must be called before the colloid map/build update for the step.
 *
 *****************************************************************************/

int bbl_link_diag_step_set(bbl_link_diag_t * diag, int step) {

  if (diag == NULL) return 0;

  diag->current_step = step;
  diag->active_step = 0;
  diag->nrecord = 0;
  diag->pass1_total_owner_count = 0;

  bbl_link_diag_six_zero(diag->build_source);
  bbl_link_diag_six_zero(diag->pass1_fluid);
  bbl_link_diag_six_zero(diag->pass1_colloid);
  bbl_link_diag_six_zero(diag->pass1_boundary);
  bbl_link_diag_six_zero(diag->pass1_active);
  bbl_link_diag_six_zero(diag->pass1_total);
  bbl_link_diag_six_zero(diag->pass2_reference);

  if (step >= diag->output_start && step <= diag->output_end &&
      ((step - diag->output_start) % diag->frequency) == 0) {
    diag->active_step = 1;
  }

  return 0;
}

/*****************************************************************************
 *
 *  bbl_link_diag_active_for
 *
 *****************************************************************************/

int bbl_link_diag_active_for(const bbl_link_diag_t * diag,
                             int particle_index) {

  if (diag == NULL) return 0;

  return (diag->active_step && particle_index == diag->target_index);
}

/*****************************************************************************
 *
 *  bbl_link_diag_build_source_add
 *
 *  Source hook used by build.c. Site ownership is tested here, so rebuild
 *  operations on halos cannot be counted twice.
 *
 *****************************************************************************/

int bbl_link_diag_build_source_add(cs_t * cs, int particle_index,
                                   int site_index, const double force[3],
                                   const double torque[3]) {

  bbl_link_diag_t * diag = bbl_link_diag_registered_;

  assert(cs);
  assert(force);
  assert(torque);

  if (diag == NULL || diag->active_step == 0) return 0;
  if (particle_index != diag->target_index) return 0;
  if (cs != diag->cs) return 0;
  if (!bbl_link_diag_site_is_local_interior(cs, site_index)) return 0;

  bbl_link_diag_six_add(diag->build_source, force, torque);

  return 0;
}

/*****************************************************************************
 *
 *  bbl_link_diag_pass1_component_add
 *
 *****************************************************************************/

int bbl_link_diag_pass1_component_add(
    bbl_link_diag_t * diag, int particle_index,
    bbl_link_diag_pass1_component_t component,
    const double force[3], const double torque[3]) {

  double * target = NULL;

  assert(force);
  assert(torque);

  if (diag == NULL || diag->active_step == 0) return 0;
  if (particle_index != diag->target_index) return 0;

  switch (component) {
    case BBL_LINK_DIAG_PASS1_FLUID:
      target = diag->pass1_fluid;
      break;
    case BBL_LINK_DIAG_PASS1_COLLOID:
      target = diag->pass1_colloid;
      break;
    case BBL_LINK_DIAG_PASS1_BOUNDARY:
      target = diag->pass1_boundary;
      break;
    case BBL_LINK_DIAG_PASS1_ACTIVE:
      target = diag->pass1_active;
      break;
    default:
      pe_fatal(diag->pe, "Unknown exact BBL pass1 diagnostic component\n");
  }

  bbl_link_diag_six_add(target, force, torque);

  return 0;
}

/*****************************************************************************
 *
 *  bbl_link_diag_pass1_total_capture
 *
 *  Called after COLLOID_SUM_DYNAMICS. The local colloid list contains the
 *  unique physical particle owner, whose f0/t0 now contain the actual global
 *  pass-1 totals used by the implicit update.
 *
 *****************************************************************************/

int bbl_link_diag_pass1_total_capture(bbl_link_diag_t * diag,
                                      colloids_info_t * cinfo) {

  colloid_t * pc = NULL;

  assert(cinfo);

  if (diag == NULL || diag->active_step == 0) return 0;

  colloids_info_local_head(cinfo, &pc);
  for (; pc; pc = pc->nextlocal) {
    if (pc->s.index == diag->target_index) {
      diag->pass1_total[BBL_DIAG_FX] = pc->f0[0];
      diag->pass1_total[BBL_DIAG_FY] = pc->f0[1];
      diag->pass1_total[BBL_DIAG_FZ] = pc->f0[2];
      diag->pass1_total[BBL_DIAG_TX] = pc->t0[0];
      diag->pass1_total[BBL_DIAG_TY] = pc->t0[1];
      diag->pass1_total[BBL_DIAG_TZ] = pc->t0[2];
      diag->pass1_total_owner_count += 1;
    }
  }

  return 0;
}

/*****************************************************************************
 *
 *  bbl_link_diag_pass2_reference_add
 *
 *****************************************************************************/

int bbl_link_diag_pass2_reference_add(bbl_link_diag_t * diag,
                                      int particle_index,
                                      const double force[3],
                                      const double torque[3]) {

  assert(force);
  assert(torque);

  if (diag == NULL || diag->active_step == 0) return 0;
  if (particle_index != diag->target_index) return 0;

  bbl_link_diag_six_add(diag->pass2_reference, force, torque);

  return 0;
}

/*****************************************************************************
 *
 *  bbl_link_diag_record
 *
 *****************************************************************************/

int bbl_link_diag_record(bbl_link_diag_t * diag,
                         const bbl_link_record_t * record) {

  bbl_link_record_t * tmp = NULL;
  int new_capacity = 0;

  assert(record);

  if (diag == NULL || diag->active_step == 0) return 0;
  if (record->particle_index != diag->target_index) return 0;

  if (diag->nrecord == diag->capacity) {
    new_capacity = (diag->capacity == 0)
                 ? BBL_LINK_DIAG_INITIAL_CAPACITY
                 : 2*diag->capacity;
    tmp = (bbl_link_record_t *) realloc(
        diag->records, (size_t) new_capacity*sizeof(bbl_link_record_t));
    if (tmp == NULL) {
      pe_fatal(diag->pe, "realloc exact BBL link records failed\n");
    }
    diag->records = tmp;
    diag->capacity = new_capacity;
  }

  diag->records[diag->nrecord] = *record;
  diag->records[diag->nrecord].step = diag->current_step;
  diag->nrecord += 1;

  return 0;
}

/*****************************************************************************
 *
 *  bbl_link_diag_end
 *
 *****************************************************************************/

int bbl_link_diag_end(bbl_link_diag_t * diag, colloids_info_t * cinfo) {

  enum {
    S_BUILD = 0,
    S_P1_FLUID = S_BUILD + BBL_DIAG_N6,
    S_P1_COLLOID = S_P1_FLUID + BBL_DIAG_N6,
    S_P1_BOUNDARY = S_P1_COLLOID + BBL_DIAG_N6,
    S_P1_ACTIVE = S_P1_BOUNDARY + BBL_DIAG_N6,
    S_P1_TOTAL = S_P1_ACTIVE + BBL_DIAG_N6,
    S_P2_REF = S_P1_TOTAL + BBL_DIAG_N6,
    S_ZETA = S_P2_REF + BBL_DIAG_N6,
    S_HYDRO = S_ZETA + BBL_DIAG_N6,
    S_RADIUS = S_HYDRO + BBL_DIAG_N6,
    S_NDOUBLE = S_RADIUS + 1
  };

  bbl_link_record_t * all_records = NULL;
  colloid_t * pc = NULL;
  int * record_counts = NULL;
  int * byte_counts = NULL;
  int * byte_displacements = NULL;
  int send_bytes = 0;
  int total_records = 0;
  int total_bytes = 0;
  int owner_count_local = 0;
  int owner_count = 0;
  int pass1_owner_count_local = 0;
  int pass1_owner_count = 0;
  int duplicate_count = 0;
  int unique_count = 0;
  int nonfinite_count = 0;
  int fail = 0;
  double summary_local[S_NDOUBLE];
  double summary_global[S_NDOUBLE];
  double zeta_force[3] = {0};
  double zeta_torque[3] = {0};
  double link_physical[BBL_DIAG_N6] = {0};
  double pass1_components[BBL_DIAG_N6] = {0};
  double hydro_reconstructed[BBL_DIAG_N6] = {0};
  double physical_force_abs = 0.0;
  double physical_force_scaled = 0.0;
  double physical_torque_abs = 0.0;
  double physical_torque_scaled = 0.0;
  double pass1_force_abs = 0.0;
  double pass1_force_scaled = 0.0;
  double pass1_torque_abs = 0.0;
  double pass1_torque_scaled = 0.0;
  double hydro_force_abs = 0.0;
  double hydro_force_scaled = 0.0;
  double hydro_torque_abs = 0.0;
  double hydro_torque_scaled = 0.0;
  double force_scale_physical = 0.0;
  double force_scale_hydro = 0.0;
  double radius = 0.0;

  assert(cinfo);

  if (diag == NULL || diag->active_step == 0) return 0;

  memset(summary_local, 0, sizeof(summary_local));
  memset(summary_global, 0, sizeof(summary_global));

  for (int i = 0; i < BBL_DIAG_N6; i++) {
    summary_local[S_BUILD + i] = diag->build_source[i];
    summary_local[S_P1_FLUID + i] = diag->pass1_fluid[i];
    summary_local[S_P1_COLLOID + i] = diag->pass1_colloid[i];
    summary_local[S_P1_BOUNDARY + i] = diag->pass1_boundary[i];
    summary_local[S_P1_ACTIVE + i] = diag->pass1_active[i];
    summary_local[S_P1_TOTAL + i] = diag->pass1_total[i];
    summary_local[S_P2_REF + i] = diag->pass2_reference[i];
  }

  pass1_owner_count_local = diag->pass1_total_owner_count;

  /* Unique physical owner supplies Ludwig's original diagnostic and zeta
   * correction using the same current v,w as record_force_torque(). */
  colloids_info_local_head(cinfo, &pc);
  for (; pc; pc = pc->nextlocal) {
    if (pc->s.index == diag->target_index) {
      owner_count_local += 1;
      bbl_link_diag_zeta_correction(pc, zeta_force, zeta_torque);

      summary_local[S_ZETA + BBL_DIAG_FX] = zeta_force[0];
      summary_local[S_ZETA + BBL_DIAG_FY] = zeta_force[1];
      summary_local[S_ZETA + BBL_DIAG_FZ] = zeta_force[2];
      summary_local[S_ZETA + BBL_DIAG_TX] = zeta_torque[0];
      summary_local[S_ZETA + BBL_DIAG_TY] = zeta_torque[1];
      summary_local[S_ZETA + BBL_DIAG_TZ] = zeta_torque[2];

      summary_local[S_HYDRO + BBL_DIAG_FX] = pc->diagnostic.fhydro[0];
      summary_local[S_HYDRO + BBL_DIAG_FY] = pc->diagnostic.fhydro[1];
      summary_local[S_HYDRO + BBL_DIAG_FZ] = pc->diagnostic.fhydro[2];
      summary_local[S_HYDRO + BBL_DIAG_TX] = pc->diagnostic.Thydro[0];
      summary_local[S_HYDRO + BBL_DIAG_TY] = pc->diagnostic.Thydro[1];
      summary_local[S_HYDRO + BBL_DIAG_TZ] = pc->diagnostic.Thydro[2];
      summary_local[S_RADIUS] = pc->s.a0;
    }
  }

  MPI_Reduce(&owner_count_local, &owner_count, 1, MPI_INT, MPI_SUM, 0,
             diag->comm);
  MPI_Reduce(&pass1_owner_count_local, &pass1_owner_count, 1, MPI_INT,
             MPI_SUM, 0, diag->comm);
  MPI_Reduce(summary_local, summary_global, S_NDOUBLE, MPI_DOUBLE, MPI_SUM, 0,
             diag->comm);

  send_bytes = diag->nrecord*(int) sizeof(bbl_link_record_t);

  if (diag->rank == 0) {
    record_counts = (int *) calloc((size_t) diag->size, sizeof(int));
    byte_counts = (int *) calloc((size_t) diag->size, sizeof(int));
    byte_displacements = (int *) calloc((size_t) diag->size, sizeof(int));
    if (record_counts == NULL || byte_counts == NULL ||
        byte_displacements == NULL) {
      pe_fatal(diag->pe, "calloc exact BBL gather metadata failed\n");
    }
  }

  MPI_Gather(&diag->nrecord, 1, MPI_INT, record_counts, 1, MPI_INT, 0,
             diag->comm);

  if (diag->rank == 0) {
    for (int irank = 0; irank < diag->size; irank++) {
      byte_counts[irank] = record_counts[irank]*(int) sizeof(bbl_link_record_t);
      byte_displacements[irank] = total_bytes;
      total_bytes += byte_counts[irank];
      total_records += record_counts[irank];
    }
    if (total_records > 0) {
      all_records = (bbl_link_record_t *) malloc(
          (size_t) total_records*sizeof(bbl_link_record_t));
      if (all_records == NULL) {
        pe_fatal(diag->pe, "malloc exact BBL gathered records failed\n");
      }
    }
  }

  MPI_Gatherv(diag->records, send_bytes, MPI_BYTE,
              all_records, byte_counts, byte_displacements, MPI_BYTE,
              0, diag->comm);

  if (diag->rank == 0) {
    const double * build = summary_global + S_BUILD;
    const double * p1_fluid = summary_global + S_P1_FLUID;
    const double * p1_colloid = summary_global + S_P1_COLLOID;
    const double * p1_boundary = summary_global + S_P1_BOUNDARY;
    const double * p1_active = summary_global + S_P1_ACTIVE;
    const double * p1_total = summary_global + S_P1_TOTAL;
    const double * p2_ref = summary_global + S_P2_REF;
    const double * zeta = summary_global + S_ZETA;
    const double * hydro_original = summary_global + S_HYDRO;

    radius = summary_global[S_RADIUS];

    if (total_records > 1) {
      qsort(all_records, (size_t) total_records,
            sizeof(bbl_link_record_t), bbl_link_diag_record_compare);
    }

    for (int irecord = 0; irecord < total_records; irecord++) {
      int duplicate = 0;
      bbl_link_record_t * r = all_records + irecord;

      if (irecord > 0) {
        duplicate = bbl_link_diag_same_key(r, all_records + irecord - 1);
      }
      if (duplicate) {
        duplicate_count += 1;
        continue;
      }

      unique_count += 1;
      if (!bbl_link_diag_finite_record(r)) nonfinite_count += 1;

      link_physical[BBL_DIAG_FX] += r->dF[0];
      link_physical[BBL_DIAG_FY] += r->dF[1];
      link_physical[BBL_DIAG_FZ] += r->dF[2];
      link_physical[BBL_DIAG_TX] += r->dT[0];
      link_physical[BBL_DIAG_TY] += r->dT[1];
      link_physical[BBL_DIAG_TZ] += r->dT[2];

      bbl_link_diag_write_record(diag->link_fp, r);
    }

    for (int i = 0; i < BBL_DIAG_N6; i++) {
      pass1_components[i] = build[i] + p1_fluid[i] + p1_colloid[i]
                          + p1_boundary[i] + p1_active[i];
      hydro_reconstructed[i] = pass1_components[i] - zeta[i];
    }

    physical_force_scaled = bbl_link_diag_force_scaled_error(
        link_physical, p2_ref, &physical_force_abs);
    force_scale_physical = fmax(bbl_link_diag_norm3(link_physical),
                                bbl_link_diag_norm3(p2_ref));
    physical_torque_scaled = bbl_link_diag_torque_scaled_error(
        link_physical + 3, p2_ref + 3, radius, force_scale_physical,
        &physical_torque_abs);

    pass1_force_scaled = bbl_link_diag_force_scaled_error(
        pass1_components, p1_total, &pass1_force_abs);
    force_scale_hydro = fmax(bbl_link_diag_norm3(pass1_components),
                             bbl_link_diag_norm3(p1_total));
    pass1_torque_scaled = bbl_link_diag_torque_scaled_error(
        pass1_components + 3, p1_total + 3, radius, force_scale_hydro,
        &pass1_torque_abs);

    hydro_force_scaled = bbl_link_diag_force_scaled_error(
        hydro_reconstructed, hydro_original, &hydro_force_abs);
    force_scale_hydro = fmax(bbl_link_diag_norm3(hydro_reconstructed),
                             bbl_link_diag_norm3(hydro_original));
    hydro_torque_scaled = bbl_link_diag_torque_scaled_error(
        hydro_reconstructed + 3, hydro_original + 3, radius,
        force_scale_hydro, &hydro_torque_abs);

    fprintf(diag->closure_fp, "%d,%d,", diag->current_step,
            diag->target_index);

    bbl_link_diag_write_vec6(diag->closure_fp, link_physical);
    bbl_link_diag_write_vec6(diag->closure_fp, p2_ref);
    bbl_link_diag_write_vec6(diag->closure_fp, build);
    bbl_link_diag_write_vec6(diag->closure_fp, p1_fluid);
    bbl_link_diag_write_vec6(diag->closure_fp, p1_colloid);
    bbl_link_diag_write_vec6(diag->closure_fp, p1_boundary);
    bbl_link_diag_write_vec6(diag->closure_fp, p1_active);
    bbl_link_diag_write_vec6(diag->closure_fp, pass1_components);
    bbl_link_diag_write_vec6(diag->closure_fp, p1_total);
    bbl_link_diag_write_vec6(diag->closure_fp, zeta);
    bbl_link_diag_write_vec6(diag->closure_fp, hydro_reconstructed);
    bbl_link_diag_write_vec6(diag->closure_fp, hydro_original);

    fprintf(diag->closure_fp,
            "%.17e,%.17e,%.17e,%.17e,"
            "%.17e,%.17e,%.17e,%.17e,"
            "%.17e,%.17e,%.17e,%.17e,"
            "%d,%d,%d,%d,%d,%d,"
            "fluid_on_colloid,current_colloid_centre,current_pass2_BBL_step,"
            "pass1_f0_t0_minus_zeta_current_updated_velocity\n",
            physical_force_abs, physical_force_scaled,
            physical_torque_abs, physical_torque_scaled,
            pass1_force_abs, pass1_force_scaled,
            pass1_torque_abs, pass1_torque_scaled,
            hydro_force_abs, hydro_force_scaled,
            hydro_torque_abs, hydro_torque_scaled,
            owner_count, pass1_owner_count, total_records, unique_count,
            duplicate_count, nonfinite_count);

    fflush(diag->link_fp);
    fflush(diag->closure_fp);

    fail = (owner_count != 1 || pass1_owner_count != 1 ||
            unique_count <= 0 || duplicate_count != 0 ||
            nonfinite_count != 0 ||
            !isfinite(physical_force_scaled) ||
            !isfinite(physical_torque_scaled) ||
            !isfinite(pass1_force_scaled) ||
            !isfinite(pass1_torque_scaled) ||
            !isfinite(hydro_force_scaled) ||
            !isfinite(hydro_torque_scaled) ||
            physical_force_scaled > diag->tolerance ||
            physical_torque_scaled > diag->tolerance ||
            pass1_force_scaled > diag->tolerance ||
            pass1_torque_scaled > diag->tolerance ||
            hydro_force_scaled > diag->tolerance ||
            hydro_torque_scaled > diag->tolerance);

    if (fail && diag->action == BBL_LINK_DIAG_ACTION_WARN) {
      pe_warn(diag->pe,
              "Exact BBL v2 closure failed at step %d: physical F %.7e T %.7e; pass1 F %.7e T %.7e; hydro F %.7e T %.7e; duplicates %d nonfinite %d\n",
              diag->current_step,
              physical_force_scaled, physical_torque_scaled,
              pass1_force_scaled, pass1_torque_scaled,
              hydro_force_scaled, hydro_torque_scaled,
              duplicate_count, nonfinite_count);
    }
  }

  MPI_Bcast(&fail, 1, MPI_INT, 0, diag->comm);

  free(all_records);
  free(record_counts);
  free(byte_counts);
  free(byte_displacements);

  if (fail && diag->action == BBL_LINK_DIAG_ACTION_FATAL) {
    pe_fatal(diag->pe,
             "Exact BBL two-layer closure failed at step %d\n",
             diag->current_step);
  }

  return 0;
}

/*****************************************************************************
 *
 *  Internal helpers
 *
 *****************************************************************************/

static int bbl_link_diag_record_compare(const void * a, const void * b) {

  const bbl_link_record_t * ra = (const bbl_link_record_t *) a;
  const bbl_link_record_t * rb = (const bbl_link_record_t *) b;

#define BBL_COMPARE_FIELD(field) \
  do { if (ra->field < rb->field) return -1; \
       if (ra->field > rb->field) return +1; } while (0)

  BBL_COMPARE_FIELD(particle_index);
  BBL_COMPARE_FIELD(fluid_site[0]);
  BBL_COMPARE_FIELD(fluid_site[1]);
  BBL_COMPARE_FIELD(fluid_site[2]);
  BBL_COMPARE_FIELD(p);

#undef BBL_COMPARE_FIELD

  return 0;
}

static int bbl_link_diag_same_key(const bbl_link_record_t * a,
                                  const bbl_link_record_t * b) {

  assert(a);
  assert(b);

  return (a->particle_index == b->particle_index &&
          a->fluid_site[0] == b->fluid_site[0] &&
          a->fluid_site[1] == b->fluid_site[1] &&
          a->fluid_site[2] == b->fluid_site[2] &&
          a->p == b->p);
}

static double bbl_link_diag_norm3(const double a[3]) {

  return sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
}

static int bbl_link_diag_finite_record(const bbl_link_record_t * r) {

  assert(r);

  for (int ia = 0; ia < 3; ia++) {
    if (!isfinite(r->rb[ia])) return 0;
    if (!isfinite(r->v[ia])) return 0;
    if (!isfinite(r->w[ia])) return 0;
    if (!isfinite(r->ub[ia])) return 0;
    if (!isfinite(r->dF[ia])) return 0;
    if (!isfinite(r->dT[ia])) return 0;
  }

  if (!isfinite(r->fdist_in)) return 0;
  if (!isfinite(r->dm)) return 0;
  if (!isfinite(r->df)) return 0;
  if (!isfinite(r->deltam)) return 0;
  if (!isfinite(r->dms)) return 0;
  if (!isfinite(r->sump)) return 0;

  return 1;
}

static int bbl_link_diag_parse_action(pe_t * pe, rt_t * rt) {

  char action[32] = "warn";
  int n = 0;

  assert(pe);
  assert(rt);

  n = rt_string_parameter(rt, "colloid_bbl_link_closure_action",
                          action, sizeof(action));
  if (n == 0 || strcmp(action, "warn") == 0) {
    return BBL_LINK_DIAG_ACTION_WARN;
  }
  if (strcmp(action, "fatal") == 0) {
    return BBL_LINK_DIAG_ACTION_FATAL;
  }
  if (strcmp(action, "off") == 0) {
    return BBL_LINK_DIAG_ACTION_OFF;
  }

  pe_fatal(pe,
           "colloid_bbl_link_closure_action must be off, warn, or fatal\n");

  return BBL_LINK_DIAG_ACTION_WARN;
}

static int bbl_link_diag_switch_default(rt_t * rt, const char * key,
                                        int default_value) {

  assert(rt);
  assert(key);

  if (rt_key_present(rt, key)) return rt_switch(rt, key);

  return default_value;
}

static int bbl_link_diag_site_is_local_interior(cs_t * cs, int index) {

  int ijk[3] = {0};
  int nlocal[3] = {0};

  assert(cs);

  cs_index_to_ijk(cs, index, ijk);
  cs_nlocal(cs, nlocal);

  return (ijk[0] >= 1 && ijk[0] <= nlocal[0] &&
          ijk[1] >= 1 && ijk[1] <= nlocal[1] &&
          ijk[2] >= 1 && ijk[2] <= nlocal[2]);
}

static void bbl_link_diag_six_zero(double a[BBL_DIAG_N6]) {

  assert(a);

  for (int i = 0; i < BBL_DIAG_N6; i++) a[i] = 0.0;
}

static void bbl_link_diag_six_add(double dst[BBL_DIAG_N6],
                                  const double force[3],
                                  const double torque[3]) {

  assert(dst);
  assert(force);
  assert(torque);

  dst[BBL_DIAG_FX] += force[0];
  dst[BBL_DIAG_FY] += force[1];
  dst[BBL_DIAG_FZ] += force[2];
  dst[BBL_DIAG_TX] += torque[0];
  dst[BBL_DIAG_TY] += torque[1];
  dst[BBL_DIAG_TZ] += torque[2];
}

static void bbl_link_diag_zeta_correction(const colloid_t * pc,
                                          double force[3],
                                          double torque[3]) {

  assert(pc);
  assert(force);
  assert(torque);

  force[0] = pc->zeta[0]*pc->s.v[0]
           + pc->zeta[1]*pc->s.v[1]
           + pc->zeta[2]*pc->s.v[2]
           + pc->zeta[3]*pc->s.w[0]
           + pc->zeta[4]*pc->s.w[1]
           + pc->zeta[5]*pc->s.w[2];
  force[1] = pc->zeta[1]*pc->s.v[0]
           + pc->zeta[6]*pc->s.v[1]
           + pc->zeta[7]*pc->s.v[2]
           + pc->zeta[8]*pc->s.w[0]
           + pc->zeta[9]*pc->s.w[1]
           + pc->zeta[10]*pc->s.w[2];
  force[2] = pc->zeta[2]*pc->s.v[0]
           + pc->zeta[7]*pc->s.v[1]
           + pc->zeta[11]*pc->s.v[2]
           + pc->zeta[12]*pc->s.w[0]
           + pc->zeta[13]*pc->s.w[1]
           + pc->zeta[14]*pc->s.w[2];

  torque[0] = pc->zeta[3]*pc->s.v[0]
            + pc->zeta[8]*pc->s.v[1]
            + pc->zeta[12]*pc->s.v[2]
            + pc->zeta[15]*pc->s.w[0]
            + pc->zeta[16]*pc->s.w[1]
            + pc->zeta[17]*pc->s.w[2];
  torque[1] = pc->zeta[4]*pc->s.v[0]
            + pc->zeta[9]*pc->s.v[1]
            + pc->zeta[13]*pc->s.v[2]
            + pc->zeta[16]*pc->s.w[0]
            + pc->zeta[18]*pc->s.w[1]
            + pc->zeta[19]*pc->s.w[2];
  torque[2] = pc->zeta[5]*pc->s.v[0]
            + pc->zeta[10]*pc->s.v[1]
            + pc->zeta[14]*pc->s.v[2]
            + pc->zeta[17]*pc->s.w[0]
            + pc->zeta[19]*pc->s.w[1]
            + pc->zeta[20]*pc->s.w[2];
}

static double bbl_link_diag_force_scaled_error(const double a[3],
                                               const double b[3],
                                               double * abs_error) {

  double d[3] = {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
  double scale = fmax(bbl_link_diag_norm3(a), bbl_link_diag_norm3(b));

  assert(abs_error);

  *abs_error = bbl_link_diag_norm3(d);

  return *abs_error/fmax(scale, DBL_MIN);
}

static double bbl_link_diag_torque_scaled_error(const double a[3],
                                                const double b[3],
                                                double radius,
                                                double force_scale,
                                                double * abs_error) {

  double d[3] = {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
  double scale = fmax(bbl_link_diag_norm3(a), bbl_link_diag_norm3(b));

  assert(abs_error);

  scale = fmax(scale, fabs(radius)*force_scale);
  *abs_error = bbl_link_diag_norm3(d);

  return *abs_error/fmax(scale, DBL_MIN);
}

static void bbl_link_diag_write_link_header(FILE * fp) {

  assert(fp);

  fprintf(fp,
          "step,particle_index,boundary_type,link_status,"
          "owner_cart_rank,particle_copy_is_local,"
          "fluid_gx,fluid_gy,fluid_gz,solid_gx,solid_gy,solid_gz,"
          "p,cx,cy,cz,"
          "rbx,rby,rbz,vx,vy,vz,wx,wy,wz,ubx,uby,ubz,"
          "fdist_in,dm,df,deltam,dms,sump,"
          "dFx,dFy,dFz,dTx,dTy,dTz,force_convention,"
          "torque_reference,time_semantics\n");
}

static void bbl_link_diag_write_closure_header(FILE * fp) {

  assert(fp);

#define V6(prefix) \
  prefix "Fx," prefix "Fy," prefix "Fz," \
  prefix "Tx," prefix "Ty," prefix "Tz,"

  fprintf(fp,
          "step,particle_index,"
          V6("link_physical_")
          V6("pass2_reference_")
          V6("build_source_")
          V6("f0_link_fluid_")
          V6("f0_link_colloid_")
          V6("f0_link_boundary_")
          V6("active_correction_")
          V6("f0_components_sum_")
          V6("f0_pass1_original_")
          V6("zeta_correction_")
          V6("hydro_reconstructed_")
          V6("hydro_original_")
          "physical_force_abs_error,physical_force_scaled_error,"
          "physical_torque_abs_error,physical_torque_scaled_error,"
          "pass1_component_force_abs_error,pass1_component_force_scaled_error,"
          "pass1_component_torque_abs_error,pass1_component_torque_scaled_error,"
          "hydro_reconstruction_force_abs_error,"
          "hydro_reconstruction_force_scaled_error,"
          "hydro_reconstruction_torque_abs_error,"
          "hydro_reconstruction_torque_scaled_error,"
          "particle_owner_count,pass1_owner_count,raw_record_count,"
          "physical_link_count,duplicate_count,nonfinite_count,"
          "physical_link_force_convention,physical_link_torque_reference,"
          "physical_link_time_semantics,hydro_original_semantics\n");

#undef V6
}

static void bbl_link_diag_write_record(FILE * fp,
                                       const bbl_link_record_t * r) {

  assert(fp);
  assert(r);

  fprintf(fp,
          "%d,%d,particle_bbl,LINK_FLUID,%d,%d,"
          "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,"
          "%.17e,%.17e,%.17e,"
          "%.17e,%.17e,%.17e,"
          "%.17e,%.17e,%.17e,"
          "%.17e,%.17e,%.17e,"
          "%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,"
          "%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,"
          "fluid_on_colloid,current_colloid_centre,current_pass2_BBL_step\n",
          r->step, r->particle_index, r->owner_rank,
          r->particle_copy_is_local,
          r->fluid_site[0], r->fluid_site[1], r->fluid_site[2],
          r->solid_site[0], r->solid_site[1], r->solid_site[2],
          r->p, r->c[0], r->c[1], r->c[2],
          r->rb[0], r->rb[1], r->rb[2],
          r->v[0], r->v[1], r->v[2],
          r->w[0], r->w[1], r->w[2],
          r->ub[0], r->ub[1], r->ub[2],
          r->fdist_in, r->dm, r->df, r->deltam, r->dms, r->sump,
          r->dF[0], r->dF[1], r->dF[2],
          r->dT[0], r->dT[1], r->dT[2]);
}

static void bbl_link_diag_write_vec6(FILE * fp,
                                     const double a[BBL_DIAG_N6]) {

  assert(fp);
  assert(a);

  fprintf(fp,
          "%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,",
          a[0], a[1], a[2], a[3], a[4], a[5]);
}
