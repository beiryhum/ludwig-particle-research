/*****************************************************************************
 *
 *  bbl_surface_stress_diag.c
 *
 *  Read-only local distribution / momentum-flux diagnostic on physical
 *  colloid BBL links.
 *
 *  IMPORTANT PHYSICAL SEMANTICS
 *  ----------------------------
 *  This module DOES NOT claim that a fluid-node stress tensor can be assigned
 *  a unique BBL-link surface area such that pressure + viscous traction closes
 *  exactly to the momentum-exchange force. It therefore does not manufacture
 *  dF_pressure or dF_viscous.
 *
 *  Instead it records, on the exact physical LINK_FLUID geometry:
 *    - pre-collision distributions and moments;
 *    - post-collision / pre-particle-BBL distributions and moments;
 *    - Ludwig collision velocity reconstructed with the half-force rule;
 *    - Ludwig lb_2nd_moment() quantity;
 *    - a rheology-style second-moment residual and the existing Ludwig
 *      rheology conversion factor, clearly labelled as a moment estimator.
 *
 *  This is intentionally diagnostic-only. No LB distribution, collision,
 *  BBL formula, particle state, force, torque, zeta matrix, or restart state
 *  is changed.
 *
 *****************************************************************************/

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bbl_surface_stress_diag.h"
#include "physics.h"

#define BBL_SURF_DEFAULT_FREQ 1
#define BBL_SURF_DEFAULT_TOL 1.0e-10
#define BBL_SURF_INITIAL_CAPACITY 1024
#define BBL_SURF_NSYM 6

enum bbl_surface_action_enum {
  BBL_SURF_ACTION_OFF = 0,
  BBL_SURF_ACTION_WARN,
  BBL_SURF_ACTION_FATAL
};

typedef struct bbl_surface_record_s bbl_surface_record_t;

struct bbl_surface_record_s {
  int step;
  int particle_index;
  int owner_rank;
  int particle_copy_is_local;
  int fluid_site[3];
  int solid_site[3];
  int p;
  int c[3];
  int post_seen;

  double rb[3];

  double particle_r_pre[3];
  double particle_v_pre[3];
  double particle_w_pre[3];
  double particle_r_post[3];
  double particle_v_post[3];
  double particle_w_post[3];

  double rho_pre;
  double momentum_pre[3];
  double force_total[3];
  double u_collision_pre[3];
  double p_lb_abs;
  double p_lb_gauge;

  double pi_pre[BBL_SURF_NSYM];
  double sdev_pre[BBL_SURF_NSYM];
  double sneq_collision_pre[BBL_SURF_NSYM];
  double sneq_rheology_pre[BBL_SURF_NSYM];
  double viscous_rheology_pre[BBL_SURF_NSYM];
  double cs2;
  double rho0_reference;
  double eta_shear;
  double tau_shear;
  double rheology_viscous_factor;
  double f_pre[NVELMAX];

  double rho_post;
  double momentum_post[3];
  double u_collision_from_post[3];
  double pi_post[BBL_SURF_NSYM];
  double sdev_post[BBL_SURF_NSYM];
  double f_post[NVELMAX];

  double pre_moment_identity_error;
  double post_moment_identity_error;
  double velocity_time_consistency_error;
};

struct bbl_surface_stress_diag_s {
  pe_t * pe;
  cs_t * cs;
  rt_t * rt;
  colloids_info_t * cinfo;
  lb_t * lb;
  hydro_t * hydro;
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

  char filename[FILENAME_MAX];
  char summary_filename[FILENAME_MAX];
  FILE * fp;
  FILE * summary_fp;

  bbl_surface_record_t * records;
  int nrecord;
  int capacity;
  int local_post_link_count;
  int local_unmatched_post_count;
};

static bbl_surface_stress_diag_t * bbl_surface_registered_ = NULL;

static int bbl_surface_parse_action(pe_t * pe, rt_t * rt);
static int bbl_surface_site_is_local_interior(cs_t * cs, int index);
static int bbl_surface_particle_copy_is_local(colloids_info_t * cinfo,
                                              colloid_t * pc);
static void bbl_surface_global_sites(cs_t * cs, int fluid_index,
                                     int solid_index, int fluid_site[3],
                                     int solid_site[3]);
static void bbl_surface_tensor_to_symm(const double a[3][3],
                                       double s[BBL_SURF_NSYM]);
static void bbl_surface_full_second_moment(lb_t * lb, int index,
                                           double pi[3][3]);
static double bbl_surface_tensor_identity_error(double rho,
                                                const double pi[3][3],
                                                const double sdev[3][3],
                                                double cs2);
static double bbl_surface_vec_diff_norm(const double a[3], const double b[3]);
static int bbl_surface_record_finite(const bbl_surface_record_t * r, int nvel);
static int bbl_surface_record_compare(const void * va, const void * vb);
static int bbl_surface_same_key(const bbl_surface_record_t * a,
                                const bbl_surface_record_t * b);
static bbl_surface_record_t * bbl_surface_find_local_record(
    bbl_surface_stress_diag_t * diag, int particle_index,
    const int fluid_site[3], int p);
static void bbl_surface_write_header(FILE * fp, int nvel);
static void bbl_surface_write_summary_header(FILE * fp);
static void bbl_surface_write_symm(FILE * fp, const double s[BBL_SURF_NSYM]);
static void bbl_surface_write_record(FILE * fp,
                                     const bbl_surface_record_t * r,
                                     int nvel);

/*****************************************************************************/
int bbl_surface_stress_diag_create(pe_t * pe, cs_t * cs, rt_t * rt,
                                   colloids_info_t * cinfo, lb_t * lb,
                                   hydro_t * hydro, int nstart,
                                   int have_viscosity_model,
                                   bbl_surface_stress_diag_t ** pdiag) {

  bbl_surface_stress_diag_t * diag = NULL;
  colloid_t * pc = NULL;
  int enabled = 0;
  int local_count = 0;
  int global_count = 0;
  int local_bbl = 0;
  int global_bbl = 0;
  int local_active = 0;
  int global_active = 0;
  int n = 0;
  physics_t * phys = NULL;
  double rho0_check = 1.0;
  double eta_shear_check = 0.0;
  double eta_bulk_check = 0.0;

  assert(pe);
  assert(cs);
  assert(rt);
  assert(cinfo);
  assert(lb);
  assert(pdiag);

  *pdiag = NULL;

  if (rt_key_present(rt, "colloid_surface_stress_output")) {
    enabled = rt_switch(rt, "colloid_surface_stress_output");
  }
  if (enabled == 0) return 0;

  if (hydro == NULL) {
    pe_fatal(pe, "colloid_surface_stress_output requires hydrodynamics\n");
  }
  if (lb->ndist != 1) {
    pe_fatal(pe,
             "Case08 surface stress diagnostic v1 supports ndist=1 only; "
             "no pressure/viscous interpretation is assumed for binary LB\n");
  }
  if (have_viscosity_model) {
    pe_fatal(pe,
             "Case08 surface stress diagnostic v1 is restricted to constant "
             "viscosity; variable-viscosity models require a local-tau "
             "extension before rheology-moment interpretation\n");
  }

  /* The only quantity labelled as a viscous estimator below deliberately
   * reproduces the existing Ludwig stats_rheology constant-viscosity factor.
   * Keep v1 scoped to the Case08B regime where that interpretation is exact. */
  physics_ref(&phys);
  physics_rho0(phys, &rho0_check);
  physics_eta_shear(phys, &eta_shear_check);
  physics_eta_bulk(phys, &eta_bulk_check);
  if (fabs(rho0_check - 1.0) > 64.0*DBL_EPSILON) {
    pe_fatal(pe,
             "Case08 surface stress diagnostic v1 requires rho0=1 to use "
             "the existing stats_rheology viscous conversion factor\n");
  }
  if (fabs(lb->model.cs2 - 1.0/3.0) > 64.0*DBL_EPSILON) {
    pe_fatal(pe,
             "Case08 surface stress diagnostic v1 requires cs2=1/3 to use "
             "the existing stats_rheology viscous conversion factor\n");
  }
  if (fabs(eta_bulk_check - eta_shear_check) >
      64.0*DBL_EPSILON*fmax(1.0, fabs(eta_shear_check))) {
    pe_fatal(pe,
             "Case08 surface stress diagnostic v1 requires bulk viscosity "
             "equal to shear viscosity; otherwise diagonal stress modes need "
             "a separate bulk conversion\n");
  }
  if (bbl_surface_registered_ != NULL) {
    pe_fatal(pe, "Only one BBL surface stress diagnostic may be registered\n");
  }

  diag = (bbl_surface_stress_diag_t *) calloc(1, sizeof(*diag));
  if (diag == NULL) pe_fatal(pe, "calloc(bbl_surface_stress_diag_t) failed\n");

  diag->pe = pe;
  diag->cs = cs;
  diag->rt = rt;
  diag->cinfo = cinfo;
  diag->lb = lb;
  diag->hydro = hydro;
  cs_cart_comm(cs, &diag->comm);
  MPI_Comm_rank(diag->comm, &diag->rank);
  MPI_Comm_size(diag->comm, &diag->size);

  diag->target_index = 1;
  diag->frequency = BBL_SURF_DEFAULT_FREQ;
  diag->output_start = nstart + 1;
  diag->output_end = nstart + 1;
  diag->current_step = nstart;
  diag->action = bbl_surface_parse_action(pe, rt);
  diag->tolerance = BBL_SURF_DEFAULT_TOL;

  snprintf(diag->filename, FILENAME_MAX, "%s",
           "bbl_surface_stress_local.csv");
  snprintf(diag->summary_filename, FILENAME_MAX, "%s",
           "bbl_surface_stress_summary.csv");

  rt_int_parameter(rt, "colloid_surface_stress_index", &diag->target_index);
  rt_int_parameter(rt, "colloid_surface_stress_output_freq", &diag->frequency);
  rt_int_parameter(rt, "colloid_surface_stress_output_start", &diag->output_start);
  rt_int_parameter(rt, "colloid_surface_stress_output_end", &diag->output_end);
  rt_double_parameter(rt, "colloid_surface_stress_integrity_tolerance",
                      &diag->tolerance);

  n = rt_string_parameter(rt, "colloid_surface_stress_filename",
                          diag->filename, FILENAME_MAX);
  if (n == 0) {
    snprintf(diag->filename, FILENAME_MAX, "%s",
             "bbl_surface_stress_local.csv");
  }
  n = rt_string_parameter(rt, "colloid_surface_stress_summary_filename",
                          diag->summary_filename, FILENAME_MAX);
  if (n == 0) {
    snprintf(diag->summary_filename, FILENAME_MAX, "%s",
             "bbl_surface_stress_summary.csv");
  }

  if (diag->target_index <= 0) {
    pe_fatal(pe, "colloid_surface_stress_index must be positive\n");
  }
  if (diag->frequency <= 0) {
    pe_fatal(pe, "colloid_surface_stress_output_freq must be positive\n");
  }
  if (diag->output_start <= nstart) {
    pe_fatal(pe,
             "colloid_surface_stress_output_start must be greater than N_start\n");
  }
  if (diag->output_end < diag->output_start) {
    pe_fatal(pe,
             "colloid_surface_stress_output_end must be >= output_start\n");
  }
  if (!(diag->tolerance > 0.0) || !isfinite(diag->tolerance)) {
    pe_fatal(pe,
             "colloid_surface_stress_integrity_tolerance must be finite and > 0\n");
  }

  colloids_info_local_head(cinfo, &pc);
  for (; pc; pc = pc->nextlocal) {
    if (pc->s.index == diag->target_index) {
      local_count += 1;
      if (pc->s.bc == COLLOID_BC_BBL) local_bbl += 1;
      if (pc->s.active) local_active += 1;
    }
  }
  MPI_Allreduce(&local_count, &global_count, 1, MPI_INT, MPI_SUM, diag->comm);
  MPI_Allreduce(&local_bbl, &global_bbl, 1, MPI_INT, MPI_SUM, diag->comm);
  MPI_Allreduce(&local_active, &global_active, 1, MPI_INT, MPI_SUM, diag->comm);

  if (global_count != 1) {
    pe_fatal(pe, "Surface stress target colloid %d owner count is %d\n",
             diag->target_index, global_count);
  }
  if (global_bbl != 1) {
    pe_fatal(pe, "Surface stress target colloid %d is not a BBL colloid\n",
             diag->target_index);
  }
  if (global_active != 0) {
    pe_fatal(pe,
             "Case08 surface stress diagnostic v1 supports passive BBL "
             "colloids only; active/squirmer pass1 changes the outside-site "
             "distribution before pass2\n");
  }

  if (diag->rank == 0) {
    diag->fp = fopen(diag->filename, "w");
    if (diag->fp == NULL) {
      pe_fatal(pe, "Cannot open surface stress file %s\n", diag->filename);
    }
    diag->summary_fp = fopen(diag->summary_filename, "w");
    if (diag->summary_fp == NULL) {
      pe_fatal(pe, "Cannot open surface stress summary %s\n",
               diag->summary_filename);
    }
    bbl_surface_write_header(diag->fp, lb->model.nvel);
    bbl_surface_write_summary_header(diag->summary_fp);
    fflush(diag->fp);
    fflush(diag->summary_fp);
  }

  bbl_surface_registered_ = diag;

  pe_info(pe, "\nCase08 BBL local distribution/momentum-flux diagnostic\n");
  pe_info(pe, "------------------------------------------------------\n");
  pe_info(pe, "Target colloid index:        %d\n", diag->target_index);
  pe_info(pe, "Output range:                %d to %d\n",
          diag->output_start, diag->output_end);
  pe_info(pe, "Output frequency:            %d\n", diag->frequency);
  pe_info(pe, "Local output:                %s\n", diag->filename);
  pe_info(pe, "Summary output:              %s\n", diag->summary_filename);
  pe_info(pe, "Integrity tolerance:         %.7e\n", diag->tolerance);
  pe_info(pe, "Owner rule:                  outside fluid site local interior\n");
  pe_info(pe, "Pre state:                   immediately before lb_collide\n");
  pe_info(pe, "Post state:                  post-collision/pre-particle-BBL\n");
  pe_info(pe, "Pressure quantity:           p_lb = cs2*rho (and gauge form)\n");
  pe_info(pe, "Viscous quantity:            Ludwig rheology moment estimator only\n");
  pe_info(pe, "Exact BBL force split:       NOT claimed by this diagnostic\n");

  *pdiag = diag;
  return 0;
}

/*****************************************************************************/
int bbl_surface_stress_diag_free(bbl_surface_stress_diag_t ** pdiag) {

  bbl_surface_stress_diag_t * diag = NULL;

  if (pdiag == NULL || *pdiag == NULL) return 0;
  diag = *pdiag;

  if (diag->rank == 0) {
    if (diag->fp) {
      fflush(diag->fp);
      fclose(diag->fp);
    }
    if (diag->summary_fp) {
      fflush(diag->summary_fp);
      fclose(diag->summary_fp);
    }
  }

  if (bbl_surface_registered_ == diag) bbl_surface_registered_ = NULL;
  free(diag->records);
  free(diag);
  *pdiag = NULL;
  return 0;
}

/*****************************************************************************/
int bbl_surface_stress_diag_step_set(bbl_surface_stress_diag_t * diag, int step) {

  if (diag == NULL) return 0;

  diag->current_step = step;
  diag->active_step = 0;
  diag->nrecord = 0;
  diag->local_post_link_count = 0;
  diag->local_unmatched_post_count = 0;

  if (step >= diag->output_start && step <= diag->output_end &&
      ((step - diag->output_start) % diag->frequency) == 0) {
    diag->active_step = 1;
  }

  return 0;
}

/*****************************************************************************/
int bbl_surface_stress_diag_pre_collision(bbl_surface_stress_diag_t * diag) {

  colloid_t * pc = NULL;
  physics_t * phys = NULL;
  double fbody[3] = {0};
  double fpulse[3] = {0};
  double fpulse_frequency = 0.0;
  double force_global[3] = {0};
  double rho0 = 1.0;
  double eta = 0.0;
  double cs2 = 0.0;
  double tau = 0.0;
  double visc_factor = 0.0;

  if (diag == NULL || diag->active_step == 0) return 0;

  assert(diag->lb);
  assert(diag->hydro);
  assert(diag->cinfo);

  /* Read-only host snapshots. On CPU these alias and do not copy. */
  lb_memcpy(diag->lb, tdpMemcpyDeviceToHost);
  hydro_memcpy(diag->hydro, tdpMemcpyDeviceToHost);

  physics_ref(&phys);
  physics_fbody(phys, fbody);
  physics_fpulse(phys, fpulse);
  physics_fpulse_frequency(phys, &fpulse_frequency);
  physics_rho0(phys, &rho0);
  physics_eta_shear(phys, &eta);

  {
    const double pi = 4.0*atan(1.0);
    const double t = 1.0*physics_control_timestep(phys);
    const double omega = 2.0*pi*fpulse_frequency;
    for (int ia = 0; ia < 3; ia++) {
      force_global[ia] = fbody[ia] + fpulse[ia]*sin(omega*t);
    }
  }

  cs2 = diag->lb->model.cs2;
  tau = 0.5 + eta/(rho0*cs2);
  /* Existing Ludwig stats_rheology factor:
   * viscous = -3*eta*2/(1 + 6*eta) = -3*eta/tau for rho0=1. */
  visc_factor = -3.0*eta*2.0/(1.0 + 6.0*eta);

  colloids_info_all_head(diag->cinfo, &pc);
  for (; pc; pc = pc->nextall) {
    colloid_link_t * link = NULL;

    if (pc->s.index != diag->target_index) continue;
    if (pc->s.bc != COLLOID_BC_BBL) continue;

    for (link = pc->lnk; link; link = link->next) {
      bbl_surface_record_t * r = NULL;
      double pi[3][3] = {{0}};
      double sdev[3][3] = {{0}};
      double g[3] = {0};
      double flocal[3] = {0};
      double u[3] = {0};
      double sneq_collision[3][3] = {{0}};
      double sneq_rheo[3][3] = {{0}};
      int fluid_site[3] = {0};
      int solid_site[3] = {0};
      int i = 0;
      int j = 0;
      int p = 0;

      if (link->status != LINK_FLUID) continue;
      i = link->i;
      j = link->j;
      p = link->p;

      if (!bbl_surface_site_is_local_interior(diag->cs, i)) continue;

      if (diag->nrecord == diag->capacity) {
        int newcap = (diag->capacity == 0) ? BBL_SURF_INITIAL_CAPACITY
                                           : 2*diag->capacity;
        bbl_surface_record_t * tmp = (bbl_surface_record_t *) realloc(
            diag->records, (size_t) newcap*sizeof(*tmp));
        if (tmp == NULL) pe_fatal(diag->pe, "realloc surface records failed\n");
        diag->records = tmp;
        diag->capacity = newcap;
      }

      r = diag->records + diag->nrecord++;
      memset(r, 0, sizeof(*r));

      bbl_surface_global_sites(diag->cs, i, j, fluid_site, solid_site);

      r->step = diag->current_step;
      r->particle_index = pc->s.index;
      r->owner_rank = cs_cart_rank(diag->cs);
      r->particle_copy_is_local =
          bbl_surface_particle_copy_is_local(diag->cinfo, pc);
      r->p = p;

      for (int ia = 0; ia < 3; ia++) {
        r->fluid_site[ia] = fluid_site[ia];
        r->solid_site[ia] = solid_site[ia];
        r->c[ia] = diag->lb->model.cv[p][ia];
        r->rb[ia] = link->rb[ia];
        r->particle_r_pre[ia] = pc->s.r[ia];
        r->particle_v_pre[ia] = pc->s.v[ia];
        r->particle_w_pre[ia] = pc->s.w[ia];
      }

      lb_0th_moment(diag->lb, i, LB_RHO, &r->rho_pre);
      lb_1st_moment(diag->lb, i, LB_RHO, g);
      hydro_f_local(diag->hydro, i, flocal);

      for (int ia = 0; ia < 3; ia++) {
        r->momentum_pre[ia] = g[ia];
        r->force_total[ia] = force_global[ia] + flocal[ia];
        if (r->rho_pre != 0.0) {
          u[ia] = (g[ia] + 0.5*r->force_total[ia])/r->rho_pre;
        }
        r->u_collision_pre[ia] = u[ia];
      }

      r->p_lb_abs = cs2*r->rho_pre;
      r->p_lb_gauge = cs2*(r->rho_pre - rho0);
      r->cs2 = cs2;
      r->rho0_reference = rho0;
      r->eta_shear = eta;
      r->tau_shear = tau;
      r->rheology_viscous_factor = visc_factor;

      bbl_surface_full_second_moment(diag->lb, i, pi);
      lb_2nd_moment(diag->lb, i, LB_RHO, sdev);

      for (int ia = 0; ia < 3; ia++) {
        for (int ib = 0; ib < 3; ib++) {
          double u_raw_a = 0.0;
          double u_raw_b = 0.0;
          if (r->rho_pre != 0.0) {
            u_raw_a = g[ia]/r->rho_pre;
            u_raw_b = g[ib]/r->rho_pre;
          }
          sneq_collision[ia][ib] =
              sdev[ia][ib] - r->rho_pre*u[ia]*u[ib];
          sneq_rheo[ia][ib] =
              sdev[ia][ib] - r->rho_pre*u_raw_a*u_raw_b;
        }
      }

      bbl_surface_tensor_to_symm(pi, r->pi_pre);
      bbl_surface_tensor_to_symm(sdev, r->sdev_pre);
      bbl_surface_tensor_to_symm(sneq_collision, r->sneq_collision_pre);
      bbl_surface_tensor_to_symm(sneq_rheo, r->sneq_rheology_pre);
      for (int is = 0; is < BBL_SURF_NSYM; is++) {
        r->viscous_rheology_pre[is] = visc_factor*r->sneq_rheology_pre[is];
      }

      r->pre_moment_identity_error = bbl_surface_tensor_identity_error(
          r->rho_pre, pi, sdev, cs2);

      for (int q = 0; q < diag->lb->model.nvel; q++) {
        lb_f(diag->lb, i, q, LB_RHO, r->f_pre + q);
      }
    }
  }

  return 0;
}

/*****************************************************************************/
int bbl_surface_stress_diag_prepare_pass2(bbl_surface_stress_diag_t * diag) {

  if (diag == NULL || diag->active_step == 0) return 0;
  assert(diag->hydro);

  /* Collision has written current rho/u on target. BBL itself synchronises
   * lb distributions to host; here we only synchronise hydro for diagnostics. */
  hydro_memcpy(diag->hydro, tdpMemcpyDeviceToHost);

  return 0;
}

/*****************************************************************************/
int bbl_surface_stress_diag_pass2_link(lb_t * lb, colloid_t * pc,
                                       colloid_link_t * link,
                                       int fluid_index, int solid_index) {

  bbl_surface_stress_diag_t * diag = bbl_surface_registered_;
  bbl_surface_record_t * r = NULL;
  double pi[3][3] = {{0}};
  double sdev[3][3] = {{0}};
  double g[3] = {0};
  double u_post[3] = {0};
  double u_hydro[3] = {0};
  int fluid_site[3] = {0};
  int solid_site[3] = {0};

  if (diag == NULL || diag->active_step == 0) return 0;
  if (pc == NULL || link == NULL || lb == NULL) return 0;
  if (pc->s.index != diag->target_index) return 0;
  if (pc->s.bc != COLLOID_BC_BBL) return 0;
  if (link->status != LINK_FLUID) return 0;
  if (!bbl_surface_site_is_local_interior(diag->cs, fluid_index)) return 0;

  diag->local_post_link_count += 1;

  bbl_surface_global_sites(diag->cs, fluid_index, solid_index,
                           fluid_site, solid_site);
  r = bbl_surface_find_local_record(diag, pc->s.index,
                                    fluid_site, link->p);
  if (r == NULL) {
    /* A missing pre-collision key indicates a time/ownership inconsistency.
     * Do not create a synthetic replacement row. Count it independently so
     * an extra pass2 physical link can never disappear silently. */
    diag->local_unmatched_post_count += 1;
    return 0;
  }

  r->post_seen += 1;
  for (int ia = 0; ia < 3; ia++) {
    r->particle_r_post[ia] = pc->s.r[ia];
    r->particle_v_post[ia] = pc->s.v[ia];
    r->particle_w_post[ia] = pc->s.w[ia];
  }

  lb_0th_moment(lb, fluid_index, LB_RHO, &r->rho_post);
  lb_1st_moment(lb, fluid_index, LB_RHO, g);
  hydro_u(diag->hydro, fluid_index, u_hydro);

  for (int ia = 0; ia < 3; ia++) {
    r->momentum_post[ia] = g[ia];
    if (r->rho_post != 0.0) {
      /* Post-collision momentum = pre momentum + total force, whereas
       * collision velocity uses pre momentum + 1/2 force. */
      u_post[ia] = (g[ia] - 0.5*r->force_total[ia])/r->rho_post;
    }
    r->u_collision_from_post[ia] = u_post[ia];
  }

  bbl_surface_full_second_moment(lb, fluid_index, pi);
  lb_2nd_moment(lb, fluid_index, LB_RHO, sdev);
  bbl_surface_tensor_to_symm(pi, r->pi_post);
  bbl_surface_tensor_to_symm(sdev, r->sdev_post);
  r->post_moment_identity_error = bbl_surface_tensor_identity_error(
      r->rho_post, pi, sdev, lb->model.cs2);

  for (int q = 0; q < lb->model.nvel; q++) {
    lb_f(lb, fluid_index, q, LB_RHO, r->f_post + q);
  }

  r->velocity_time_consistency_error = bbl_surface_vec_diff_norm(
      r->u_collision_pre, r->u_collision_from_post);

  /* Independent check against collision's hydro u. Keep the larger mismatch. */
  {
    double e = bbl_surface_vec_diff_norm(r->u_collision_pre, u_hydro);
    if (e > r->velocity_time_consistency_error) {
      r->velocity_time_consistency_error = e;
    }
  }

  return 0;
}

/*****************************************************************************/
int bbl_surface_stress_diag_end_step(void) {

  bbl_surface_stress_diag_t * diag = bbl_surface_registered_;
  bbl_surface_record_t * all = NULL;
  colloid_t * pc = NULL;
  int * counts = NULL;
  int * byte_counts = NULL;
  int * displs = NULL;
  int send_bytes = 0;
  int total = 0;
  int total_bytes = 0;
  int owner_local = 0;
  int owner_count = 0;
  int post_link_count = 0;
  int unmatched_post_count = 0;
  int duplicate_count = 0;
  int nonfinite_count = 0;
  int complete_count = 0;
  int unique_count = 0;
  int fail = 0;
  double max_pre_identity = 0.0;
  double max_post_identity = 0.0;
  double max_velocity_time = 0.0;

  if (diag == NULL || diag->active_step == 0) return 0;

  colloids_info_local_head(diag->cinfo, &pc);
  for (; pc; pc = pc->nextlocal) {
    if (pc->s.index == diag->target_index) owner_local += 1;
  }
  MPI_Reduce(&owner_local, &owner_count, 1, MPI_INT, MPI_SUM, 0, diag->comm);
  MPI_Reduce(&diag->local_post_link_count, &post_link_count,
             1, MPI_INT, MPI_SUM, 0, diag->comm);
  MPI_Reduce(&diag->local_unmatched_post_count, &unmatched_post_count,
             1, MPI_INT, MPI_SUM, 0, diag->comm);

  send_bytes = diag->nrecord*(int) sizeof(*diag->records);

  if (diag->rank == 0) {
    counts = (int *) calloc((size_t) diag->size, sizeof(int));
    byte_counts = (int *) calloc((size_t) diag->size, sizeof(int));
    displs = (int *) calloc((size_t) diag->size, sizeof(int));
    if (counts == NULL || byte_counts == NULL || displs == NULL) {
      pe_fatal(diag->pe, "calloc surface diagnostic gather metadata failed\n");
    }
  }

  MPI_Gather(&diag->nrecord, 1, MPI_INT, counts, 1, MPI_INT, 0, diag->comm);

  if (diag->rank == 0) {
    for (int ir = 0; ir < diag->size; ir++) {
      byte_counts[ir] = counts[ir]*(int) sizeof(*diag->records);
      displs[ir] = total_bytes;
      total_bytes += byte_counts[ir];
      total += counts[ir];
    }
    if (total > 0) {
      all = (bbl_surface_record_t *) malloc((size_t) total*sizeof(*all));
      if (all == NULL) pe_fatal(diag->pe, "malloc gathered surface records failed\n");
    }
  }

  MPI_Gatherv(diag->records, send_bytes, MPI_BYTE,
              all, byte_counts, displs, MPI_BYTE, 0, diag->comm);

  if (diag->rank == 0) {
    if (total > 1) {
      qsort(all, (size_t) total, sizeof(*all), bbl_surface_record_compare);
    }

    for (int ir = 0; ir < total; ir++) {
      int duplicate = 0;
      bbl_surface_record_t * r = all + ir;

      if (ir > 0) duplicate = bbl_surface_same_key(r, all + ir - 1);
      if (duplicate) {
        duplicate_count += 1;
        continue;
      }

      unique_count += 1;
      if (r->post_seen == 1) complete_count += 1;
      if (!bbl_surface_record_finite(r, diag->lb->model.nvel)) {
        nonfinite_count += 1;
      }

      if (r->pre_moment_identity_error > max_pre_identity) {
        max_pre_identity = r->pre_moment_identity_error;
      }
      if (r->post_moment_identity_error > max_post_identity) {
        max_post_identity = r->post_moment_identity_error;
      }
      if (r->velocity_time_consistency_error > max_velocity_time) {
        max_velocity_time = r->velocity_time_consistency_error;
      }

      bbl_surface_write_record(diag->fp, r, diag->lb->model.nvel);
    }

    fprintf(diag->summary_fp,
            "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%.17e,%.17e,%.17e,"
            "pre_collision_after_force_and_viscosity_update,"
            "post_collision_after_halo_open_wall_setup_and_pass1_before_particle_BBL_write,"
            "outside_fluid_site_local_interior,"
            "p_lb_isothermal_cs2_rho_not_general_thermodynamic_pressure,"
            "ludwig_rheology_second_moment_estimator_not_exact_BBL_traction,"
            "no_exact_pressure_viscous_BBL_force_split_claimed\n",
            diag->current_step, diag->target_index, owner_count,
            total, post_link_count, unmatched_post_count, unique_count,
            complete_count, duplicate_count, nonfinite_count,
            max_pre_identity, max_post_identity, max_velocity_time);

    fflush(diag->fp);
    fflush(diag->summary_fp);

    fail = (owner_count != 1 || unique_count <= 0 ||
            post_link_count != total || unmatched_post_count != 0 ||
            complete_count != unique_count || duplicate_count != 0 ||
            nonfinite_count != 0 || !isfinite(max_pre_identity) ||
            !isfinite(max_post_identity) || !isfinite(max_velocity_time) ||
            max_pre_identity > diag->tolerance ||
            max_post_identity > diag->tolerance ||
            max_velocity_time > diag->tolerance);

    if (fail && diag->action == BBL_SURF_ACTION_WARN) {
      pe_warn(diag->pe,
              "Surface stress diagnostic integrity warning step %d: "
              "owner=%d pre=%d post=%d unmatchedPost=%d unique=%d "
              "complete=%d dup=%d nonfinite=%d prePi=%.7e postPi=%.7e "
              "uTime=%.7e\n",
              diag->current_step, owner_count, total, post_link_count,
              unmatched_post_count, unique_count, complete_count,
              duplicate_count, nonfinite_count, max_pre_identity,
              max_post_identity, max_velocity_time);
    }
  }

  MPI_Bcast(&fail, 1, MPI_INT, 0, diag->comm);

  free(all);
  free(counts);
  free(byte_counts);
  free(displs);

  if (fail && diag->action == BBL_SURF_ACTION_FATAL) {
    pe_fatal(diag->pe,
             "Surface stress diagnostic integrity check failed at step %d\n",
             diag->current_step);
  }

  return 0;
}

/*****************************************************************************/
static int bbl_surface_parse_action(pe_t * pe, rt_t * rt) {

  char action[BUFSIZ] = {0};
  int n = 0;

  assert(pe);
  assert(rt);

  n = rt_string_parameter(rt, "colloid_surface_stress_action",
                          action, BUFSIZ);
  if (n == 0 || strcmp(action, "warn") == 0) return BBL_SURF_ACTION_WARN;
  if (strcmp(action, "off") == 0) return BBL_SURF_ACTION_OFF;
  if (strcmp(action, "fatal") == 0) return BBL_SURF_ACTION_FATAL;

  pe_fatal(pe, "colloid_surface_stress_action must be off, warn, or fatal\n");
  return BBL_SURF_ACTION_WARN;
}

static int bbl_surface_site_is_local_interior(cs_t * cs, int index) {

  int ijk[3] = {0};
  int nlocal[3] = {0};

  assert(cs);
  cs_index_to_ijk(cs, index, ijk);
  cs_nlocal(cs, nlocal);

  return (ijk[X] >= 1 && ijk[X] <= nlocal[X] &&
          ijk[Y] >= 1 && ijk[Y] <= nlocal[Y] &&
          ijk[Z] >= 1 && ijk[Z] <= nlocal[Z]);
}

static int bbl_surface_particle_copy_is_local(colloids_info_t * cinfo,
                                              colloid_t * pc) {

  colloid_t * local = NULL;
  assert(cinfo);
  assert(pc);

  colloids_info_local_head(cinfo, &local);
  for (; local; local = local->nextlocal) {
    if (local == pc) return 1;
  }
  return 0;
}

static void bbl_surface_global_sites(cs_t * cs, int fluid_index,
                                     int solid_index, int fluid_site[3],
                                     int solid_site[3]) {

  int fi[3] = {0};
  int si[3] = {0};
  int noffset[3] = {0};
  int ntotal[3] = {0};
  int periodic[3] = {0};

  cs_index_to_ijk(cs, fluid_index, fi);
  cs_index_to_ijk(cs, solid_index, si);
  cs_nlocal_offset(cs, noffset);
  cs_ntotal(cs, ntotal);
  cs_periodic(cs, periodic);

  for (int ia = 0; ia < 3; ia++) {
    fluid_site[ia] = noffset[ia] + fi[ia];
    solid_site[ia] = noffset[ia] + si[ia];
    if (periodic[ia]) {
      while (fluid_site[ia] < 1) fluid_site[ia] += ntotal[ia];
      while (fluid_site[ia] > ntotal[ia]) fluid_site[ia] -= ntotal[ia];
      while (solid_site[ia] < 1) solid_site[ia] += ntotal[ia];
      while (solid_site[ia] > ntotal[ia]) solid_site[ia] -= ntotal[ia];
    }
  }
}

static void bbl_surface_tensor_to_symm(const double a[3][3],
                                       double s[BBL_SURF_NSYM]) {
  s[0] = a[X][X];
  s[1] = a[X][Y];
  s[2] = a[X][Z];
  s[3] = a[Y][Y];
  s[4] = a[Y][Z];
  s[5] = a[Z][Z];
}

static void bbl_surface_full_second_moment(lb_t * lb, int index,
                                           double pi[3][3]) {

  assert(lb);
  memset(pi, 0, 9*sizeof(double));

  for (int p = 0; p < lb->model.nvel; p++) {
    double f = 0.0;
    lb_f(lb, index, p, LB_RHO, &f);
    for (int ia = 0; ia < 3; ia++) {
      for (int ib = 0; ib < 3; ib++) {
        pi[ia][ib] += f*lb->model.cv[p][ia]*lb->model.cv[p][ib];
      }
    }
  }
}

static double bbl_surface_tensor_identity_error(double rho,
                                                const double pi[3][3],
                                                const double sdev[3][3],
                                                double cs2) {
  double err = 0.0;
  for (int ia = 0; ia < 3; ia++) {
    for (int ib = 0; ib < 3; ib++) {
      double dab = (ia == ib) ? 1.0 : 0.0;
      double d = fabs(pi[ia][ib] - (sdev[ia][ib] + cs2*rho*dab));
      if (d > err) err = d;
    }
  }
  return err;
}

static double bbl_surface_vec_diff_norm(const double a[3], const double b[3]) {
  double d0 = a[0] - b[0];
  double d1 = a[1] - b[1];
  double d2 = a[2] - b[2];
  return sqrt(d0*d0 + d1*d1 + d2*d2);
}

static int bbl_surface_record_finite(const bbl_surface_record_t * r, int nvel) {

  const double * blocks[] = {
    r->rb, r->particle_r_pre, r->particle_v_pre, r->particle_w_pre,
    r->particle_r_post, r->particle_v_post, r->particle_w_post,
    r->momentum_pre, r->force_total, r->u_collision_pre,
    r->pi_pre, r->sdev_pre, r->sneq_collision_pre,
    r->sneq_rheology_pre, r->viscous_rheology_pre,
    r->momentum_post, r->u_collision_from_post, r->pi_post, r->sdev_post
  };
  const int sizes[] = {3,3,3,3,3,3,3,3,3,3,6,6,6,6,6,3,3,6,6};

  if (r->post_seen != 1) return 0;
  if (!isfinite(r->rho_pre) || !isfinite(r->rho_post) ||
      !isfinite(r->p_lb_abs) || !isfinite(r->p_lb_gauge) ||
      !isfinite(r->cs2) || !isfinite(r->rho0_reference) ||
      !isfinite(r->eta_shear) || !isfinite(r->tau_shear) ||
      !isfinite(r->rheology_viscous_factor) ||
      !isfinite(r->pre_moment_identity_error) ||
      !isfinite(r->post_moment_identity_error) ||
      !isfinite(r->velocity_time_consistency_error)) return 0;

  for (size_t ib = 0; ib < sizeof(blocks)/sizeof(blocks[0]); ib++) {
    for (int i = 0; i < sizes[ib]; i++) {
      if (!isfinite(blocks[ib][i])) return 0;
    }
  }
  for (int p = 0; p < nvel; p++) {
    if (!isfinite(r->f_pre[p]) || !isfinite(r->f_post[p])) return 0;
  }
  return 1;
}

#define BBL_SURF_CMP_INT(field) \
  do { if (a->field < b->field) return -1; if (a->field > b->field) return 1; } while (0)

static int bbl_surface_record_compare(const void * va, const void * vb) {
  const bbl_surface_record_t * a = (const bbl_surface_record_t *) va;
  const bbl_surface_record_t * b = (const bbl_surface_record_t *) vb;
  BBL_SURF_CMP_INT(step);
  BBL_SURF_CMP_INT(particle_index);
  BBL_SURF_CMP_INT(fluid_site[0]);
  BBL_SURF_CMP_INT(fluid_site[1]);
  BBL_SURF_CMP_INT(fluid_site[2]);
  BBL_SURF_CMP_INT(p);
  BBL_SURF_CMP_INT(owner_rank);
  return 0;
}
#undef BBL_SURF_CMP_INT

static int bbl_surface_same_key(const bbl_surface_record_t * a,
                                const bbl_surface_record_t * b) {
  return (a->step == b->step &&
          a->particle_index == b->particle_index &&
          a->fluid_site[0] == b->fluid_site[0] &&
          a->fluid_site[1] == b->fluid_site[1] &&
          a->fluid_site[2] == b->fluid_site[2] &&
          a->p == b->p);
}

static bbl_surface_record_t * bbl_surface_find_local_record(
    bbl_surface_stress_diag_t * diag, int particle_index,
    const int fluid_site[3], int p) {

  assert(diag);
  for (int i = 0; i < diag->nrecord; i++) {
    bbl_surface_record_t * r = diag->records + i;
    if (r->particle_index == particle_index && r->p == p &&
        r->fluid_site[0] == fluid_site[0] &&
        r->fluid_site[1] == fluid_site[1] &&
        r->fluid_site[2] == fluid_site[2]) {
      return r;
    }
  }
  return NULL;
}

static void bbl_surface_write_header(FILE * fp, int nvel) {
  assert(fp);

  fprintf(fp,
    "step,particle_index,owner_cart_rank,particle_copy_is_local,"
    "fluid_gx,fluid_gy,fluid_gz,solid_gx,solid_gy,solid_gz,p,cx,cy,cz,"
    "rbx,rby,rbz,"
    "particle_x_pre,particle_y_pre,particle_z_pre,"
    "particle_vx_pre,particle_vy_pre,particle_vz_pre,"
    "particle_wx_pre,particle_wy_pre,particle_wz_pre,"
    "particle_x_post,particle_y_post,particle_z_post,"
    "particle_vx_post,particle_vy_post,particle_vz_post,"
    "particle_wx_post,particle_wy_post,particle_wz_post,"
    "rho_pre,gx_pre,gy_pre,gz_pre,Fx_collision,Fy_collision,Fz_collision,"
    "ux_collision_pre,uy_collision_pre,uz_collision_pre,p_lb_abs,p_lb_gauge,"
    "cs2,rho0_reference,eta_shear,tau_shear,"
    "Pi_xx_pre,Pi_xy_pre,Pi_xz_pre,Pi_yy_pre,Pi_yz_pre,Pi_zz_pre,"
    "Sdev_xx_pre,Sdev_xy_pre,Sdev_xz_pre,Sdev_yy_pre,Sdev_yz_pre,Sdev_zz_pre,"
    "Sneq_collision_xx_pre,Sneq_collision_xy_pre,Sneq_collision_xz_pre,"
    "Sneq_collision_yy_pre,Sneq_collision_yz_pre,Sneq_collision_zz_pre,"
    "Sneq_rheology_xx_pre,Sneq_rheology_xy_pre,Sneq_rheology_xz_pre,"
    "Sneq_rheology_yy_pre,Sneq_rheology_yz_pre,Sneq_rheology_zz_pre,"
    "visc_rheology_xx_pre,visc_rheology_xy_pre,visc_rheology_xz_pre,"
    "visc_rheology_yy_pre,visc_rheology_yz_pre,visc_rheology_zz_pre,"
    "rheology_viscous_factor,"
    "rho_post,gx_post,gy_post,gz_post,"
    "ux_collision_from_post,uy_collision_from_post,uz_collision_from_post,"
    "Pi_xx_post,Pi_xy_post,Pi_xz_post,Pi_yy_post,Pi_yz_post,Pi_zz_post,"
    "Sdev_xx_post,Sdev_xy_post,Sdev_xz_post,Sdev_yy_post,Sdev_yz_post,Sdev_zz_post,"
    "pre_moment_identity_error,post_moment_identity_error,"
    "velocity_time_consistency_error,post_seen");

  for (int p = 0; p < nvel; p++) fprintf(fp, ",f_pre_%d", p);
  for (int p = 0; p < nvel; p++) fprintf(fp, ",f_post_%d", p);
  fprintf(fp, "\n");
}

static void bbl_surface_write_summary_header(FILE * fp) {
  assert(fp);
  fprintf(fp,
    "step,particle_index,particle_owner_count,raw_pre_record_count,"
    "raw_post_link_count,unmatched_post_count,unique_record_count,"
    "complete_record_count,duplicate_count,nonfinite_count,"
    "max_pre_moment_identity_error,max_post_moment_identity_error,"
    "max_velocity_time_consistency_error,pre_time_semantics,post_time_semantics,"
    "owner_rule,pressure_semantics,viscous_semantics,closure_semantics\n");
}

static void bbl_surface_write_symm(FILE * fp, const double s[BBL_SURF_NSYM]) {
  for (int i = 0; i < BBL_SURF_NSYM; i++) {
    fprintf(fp, "%.17e,", s[i]);
  }
}

static void bbl_surface_write_record(FILE * fp,
                                     const bbl_surface_record_t * r,
                                     int nvel) {
  assert(fp);
  assert(r);

  fprintf(fp,
          "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,",
          r->step, r->particle_index, r->owner_rank,
          r->particle_copy_is_local,
          r->fluid_site[0], r->fluid_site[1], r->fluid_site[2],
          r->solid_site[0], r->solid_site[1], r->solid_site[2],
          r->p, r->c[0], r->c[1], r->c[2]);

  for (int i = 0; i < 3; i++) fprintf(fp, "%.17e,", r->rb[i]);
  for (int i = 0; i < 3; i++) fprintf(fp, "%.17e,", r->particle_r_pre[i]);
  for (int i = 0; i < 3; i++) fprintf(fp, "%.17e,", r->particle_v_pre[i]);
  for (int i = 0; i < 3; i++) fprintf(fp, "%.17e,", r->particle_w_pre[i]);
  for (int i = 0; i < 3; i++) fprintf(fp, "%.17e,", r->particle_r_post[i]);
  for (int i = 0; i < 3; i++) fprintf(fp, "%.17e,", r->particle_v_post[i]);
  for (int i = 0; i < 3; i++) fprintf(fp, "%.17e,", r->particle_w_post[i]);

  fprintf(fp, "%.17e,", r->rho_pre);
  for (int i = 0; i < 3; i++) fprintf(fp, "%.17e,", r->momentum_pre[i]);
  for (int i = 0; i < 3; i++) fprintf(fp, "%.17e,", r->force_total[i]);
  for (int i = 0; i < 3; i++) fprintf(fp, "%.17e,", r->u_collision_pre[i]);
  fprintf(fp, "%.17e,%.17e,", r->p_lb_abs, r->p_lb_gauge);
  fprintf(fp, "%.17e,%.17e,%.17e,%.17e,",
          r->cs2, r->rho0_reference, r->eta_shear, r->tau_shear);

  bbl_surface_write_symm(fp, r->pi_pre);
  bbl_surface_write_symm(fp, r->sdev_pre);
  bbl_surface_write_symm(fp, r->sneq_collision_pre);
  bbl_surface_write_symm(fp, r->sneq_rheology_pre);
  bbl_surface_write_symm(fp, r->viscous_rheology_pre);
  fprintf(fp, "%.17e,", r->rheology_viscous_factor);

  fprintf(fp, "%.17e,", r->rho_post);
  for (int i = 0; i < 3; i++) fprintf(fp, "%.17e,", r->momentum_post[i]);
  for (int i = 0; i < 3; i++) fprintf(fp, "%.17e,", r->u_collision_from_post[i]);
  bbl_surface_write_symm(fp, r->pi_post);
  bbl_surface_write_symm(fp, r->sdev_post);

  fprintf(fp, "%.17e,%.17e,%.17e,%d",
          r->pre_moment_identity_error,
          r->post_moment_identity_error,
          r->velocity_time_consistency_error,
          r->post_seen);

  for (int p = 0; p < nvel; p++) fprintf(fp, ",%.17e", r->f_pre[p]);
  for (int p = 0; p < nvel; p++) fprintf(fp, ",%.17e", r->f_post[p]);
  fprintf(fp, "\n");
}
