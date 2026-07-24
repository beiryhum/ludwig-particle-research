/*****************************************************************************
 *
 *  stats_colloid_force_io.c
 *
 *  Optional per-colloid hydrodynamic force and torque time-series output.
 *
 *  The diagnostic values are those already assembled by the BBL update:
 *
 *    pc->diagnostic.fhydro[]
 *    pc->diagnostic.Thydro[]
 *
 *  No second force calculation is performed here. Only the physical owner
 *  list (nextlocal) is inspected, so halo copies are not double counted.
 *
 *****************************************************************************/

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "stats_colloid_force_io.h"

struct stats_colloid_force_io_s {
  pe_t * pe;
  colloids_info_t * cinfo;
  int index;
  int frequency;
  int flush_frequency;
  int output_count;
  FILE * fp;
  char filename[FILENAME_MAX];
};

static int stats_colloid_force_io_target_check(stats_colloid_force_io_t * obj);
static int stats_colloid_force_io_file_open(stats_colloid_force_io_t * obj,
                                             int nstart);
static int stats_colloid_force_io_file_state(const char * filename,
                                              int * exists,
                                              int * empty,
                                              int * last_step);

/*****************************************************************************
 *
 *  stats_colloid_force_io_create
 *
 *****************************************************************************/

__host__ int stats_colloid_force_io_create(pe_t * pe, rt_t * rt,
                                            colloids_info_t * cinfo,
                                            int nstart,
                                            int hydrodynamics,
                                            stats_colloid_force_io_t ** pobj) {

  int active = 0;
  stats_colloid_force_io_t * obj = NULL;

  assert(pe);
  assert(rt);
  assert(cinfo);
  assert(pobj);

  *pobj = NULL;
  active = rt_switch(rt, "colloid_force_output");

  /* Default off: no allocation, no file, and no collectives in time steps. */
  if (active == 0) return 0;

  if (hydrodynamics == 0) {
    pe_fatal(pe, "colloid_force_output requires hydrodynamics on\n");
  }

  obj = (stats_colloid_force_io_t *) calloc(1, sizeof(*obj));
  assert(obj);
  if (obj == NULL) pe_fatal(pe, "calloc(stats_colloid_force_io_t) failed\n");

  obj->pe = pe;
  obj->cinfo = cinfo;
  obj->index = 1;
  obj->frequency = 1;
  obj->flush_frequency = 100;
  obj->output_count = 0;
  snprintf(obj->filename, FILENAME_MAX, "%s", "colloid_force_torque.csv");

  rt_int_parameter(rt, "colloid_force_index", &obj->index);
  rt_int_parameter(rt, "colloid_force_freq", &obj->frequency);
  rt_int_parameter(rt, "colloid_force_flush_freq", &obj->flush_frequency);
  rt_string_parameter(rt, "colloid_force_filename", obj->filename,
                      FILENAME_MAX);

  if (obj->index < 1) {
    pe_fatal(pe, "colloid_force_index must be greater than zero\n");
  }
  if (obj->frequency < 1) {
    pe_fatal(pe, "colloid_force_freq must be greater than zero\n");
  }
  if (obj->flush_frequency < 1) {
    pe_fatal(pe, "colloid_force_flush_freq must be greater than zero\n");
  }
  if (obj->filename[0] == '\0') {
    pe_fatal(pe, "colloid_force_filename must not be empty\n");
  }

  stats_colloid_force_io_target_check(obj);
  stats_colloid_force_io_file_open(obj, nstart);

  pe_info(pe, "\n");
  pe_info(pe, "Colloid hydrodynamic force/torque time series\n");
  pe_info(pe, "---------------------------------------------\n");
  pe_info(pe, "Output file:              %s\n", obj->filename);
  pe_info(pe, "Colloid index:            %d\n", obj->index);
  pe_info(pe, "Output frequency:         %d\n", obj->frequency);
  pe_info(pe, "Flush frequency (rows):   %d\n", obj->flush_frequency);
  pe_info(pe, "Force convention:         fluid acting on colloid\n");
  pe_info(pe, "Torque reference point:   current colloid centre\n");
  pe_info(pe, "Time semantics:           instantaneous current BBL step\n");

  *pobj = obj;

  return 0;
}

/*****************************************************************************
 *
 *  stats_colloid_force_io_target_check
 *
 *  A global colloid index is represented by exactly one item in the local
 *  owner lists. Halo copies are intentionally excluded.
 *
 *****************************************************************************/

static int stats_colloid_force_io_target_check(stats_colloid_force_io_t * obj) {

  int count_local = 0;
  int count_global = 0;
  int bbl_local = 0;
  int bbl_global = 0;
  MPI_Comm comm = MPI_COMM_NULL;
  colloid_t * pc = NULL;

  assert(obj);
  assert(obj->pe);
  assert(obj->cinfo);

  pe_mpi_comm(obj->pe, &comm);
  colloids_info_local_head(obj->cinfo, &pc);

  for ( ; pc; pc = pc->nextlocal) {
    if (pc->s.index == obj->index) {
      count_local += 1;
      if (pc->s.bc == COLLOID_BC_BBL) bbl_local += 1;
    }
  }

  MPI_Allreduce(&count_local, &count_global, 1, MPI_INT, MPI_SUM, comm);
  MPI_Allreduce(&bbl_local, &bbl_global, 1, MPI_INT, MPI_SUM, comm);

  if (count_global != 1) {
    pe_fatal(obj->pe,
             "colloid_force_index %d has global local-owner count %d; "
             "expected exactly one\n",
             obj->index, count_global);
  }

  if (bbl_global != 1) {
    pe_fatal(obj->pe,
             "colloid_force_index %d is not a resolved BBL colloid\n",
             obj->index);
  }

  return 0;
}

/*****************************************************************************
 *
 *  stats_colloid_force_io_file_state
 *
 *****************************************************************************/

static int stats_colloid_force_io_file_state(const char * filename,
                                              int * exists,
                                              int * empty,
                                              int * last_step) {

  FILE * fp = NULL;
  char line[BUFSIZ] = {0};

  assert(filename);
  assert(exists);
  assert(empty);
  assert(last_step);

  *exists = 0;
  *empty = 1;
  *last_step = -1;

  fp = fopen(filename, "r");
  if (fp == NULL) return 0;

  *exists = 1;

  while (fgets(line, BUFSIZ, fp)) {
    int step = -1;

    *empty = 0;
    if (sscanf(line, "%d,", &step) == 1) {
      if (step > *last_step) *last_step = step;
    }
  }

  fclose(fp);

  return 0;
}

/*****************************************************************************
 *
 *  stats_colloid_force_io_file_open
 *
 *****************************************************************************/

static int stats_colloid_force_io_file_open(stats_colloid_force_io_t * obj,
                                             int nstart) {

  int exists = 0;
  int empty = 1;
  int last_step = -1;
  int rank = -1;

  assert(obj);

  rank = pe_mpi_rank(obj->pe);

  if (rank == 0) {

    stats_colloid_force_io_file_state(obj->filename, &exists, &empty,
                                      &last_step);

    if (nstart == 0) {
      obj->fp = fopen(obj->filename, "w");
      empty = 1;
    }
    else {
      if (exists && last_step > nstart) {
        pe_fatal(obj->pe,
                 "Existing %s ends at step %d, beyond N_start %d\n",
                 obj->filename, last_step, nstart);
      }

      obj->fp = fopen(obj->filename, "a");

      if (exists && last_step >= 0 && last_step < nstart) {
        pe_warn(obj->pe,
                "Appending %s: previous last step %d is below N_start %d\n",
                obj->filename, last_step, nstart);
      }
    }

    if (obj->fp == NULL) {
      pe_fatal(obj->pe, "Could not open colloid force file %s\n",
               obj->filename);
    }

    if (empty) {
      fprintf(obj->fp,
              "step,index,"
              "Fx_hydro,Fy_hydro,Fz_hydro,"
              "Tx_hydro,Ty_hydro,Tz_hydro,"
              "x,y,z,vx,vy,vz,wx,wy,wz\n");
      fflush(obj->fp);
    }
  }

  return 0;
}

/*****************************************************************************
 *
 *  stats_colloid_force_io_write
 *
 *  The six force/torque components are reduced in one collective. State is
 *  reduced separately because only the physical owner contributes non-zero
 *  values. All ranks return before collectives on non-output time steps.
 *
 *****************************************************************************/

__host__ int stats_colloid_force_io_write(stats_colloid_force_io_t * obj,
                                           int timestep) {

  int count_local = 0;
  int count_global = 0;
  int rank = -1;
  double ft_local[6] = {0.0};
  double ft_global[6] = {0.0};
  double state_local[9] = {0.0};
  double state_global[9] = {0.0};
  MPI_Comm comm = MPI_COMM_NULL;
  colloid_t * pc = NULL;

  if (obj == NULL) return 0;
  if (timestep % obj->frequency != 0) return 0;

  assert(obj->pe);
  assert(obj->cinfo);

  pe_mpi_comm(obj->pe, &comm);
  rank = pe_mpi_rank(obj->pe);

  colloids_info_local_head(obj->cinfo, &pc);

  for ( ; pc; pc = pc->nextlocal) {
    if (pc->s.index == obj->index) {
      count_local += 1;

      ft_local[0] = pc->diagnostic.fhydro[X];
      ft_local[1] = pc->diagnostic.fhydro[Y];
      ft_local[2] = pc->diagnostic.fhydro[Z];
      ft_local[3] = pc->diagnostic.Thydro[X];
      ft_local[4] = pc->diagnostic.Thydro[Y];
      ft_local[5] = pc->diagnostic.Thydro[Z];

      state_local[0] = pc->s.r[X];
      state_local[1] = pc->s.r[Y];
      state_local[2] = pc->s.r[Z];
      state_local[3] = pc->s.v[X];
      state_local[4] = pc->s.v[Y];
      state_local[5] = pc->s.v[Z];
      state_local[6] = pc->s.w[X];
      state_local[7] = pc->s.w[Y];
      state_local[8] = pc->s.w[Z];
    }
  }

  MPI_Allreduce(&count_local, &count_global, 1, MPI_INT, MPI_SUM, comm);

  if (count_global != 1) {
    pe_fatal(obj->pe,
             "colloid_force_index %d has owner count %d at step %d\n",
             obj->index, count_global, timestep);
  }

  MPI_Reduce(ft_local, ft_global, 6, MPI_DOUBLE, MPI_SUM, 0, comm);
  MPI_Reduce(state_local, state_global, 9, MPI_DOUBLE, MPI_SUM, 0, comm);

  if (rank == 0) {
    int n = 0;

    for (n = 0; n < 6; n++) {
      if (!isfinite(ft_global[n])) {
        pe_fatal(obj->pe,
                 "Non-finite colloid force/torque at step %d component %d\n",
                 timestep, n);
      }
    }
    for (n = 0; n < 9; n++) {
      if (!isfinite(state_global[n])) {
        pe_fatal(obj->pe,
                 "Non-finite colloid state at step %d component %d\n",
                 timestep, n);
      }
    }

    fprintf(obj->fp,
            "%d,%d,"
            "%.17e,%.17e,%.17e,"
            "%.17e,%.17e,%.17e,"
            "%.17e,%.17e,%.17e,"
            "%.17e,%.17e,%.17e,"
            "%.17e,%.17e,%.17e\n",
            timestep, obj->index,
            ft_global[0], ft_global[1], ft_global[2],
            ft_global[3], ft_global[4], ft_global[5],
            state_global[0], state_global[1], state_global[2],
            state_global[3], state_global[4], state_global[5],
            state_global[6], state_global[7], state_global[8]);

    obj->output_count += 1;

    if (obj->output_count % obj->flush_frequency == 0) {
      fflush(obj->fp);
    }
  }

  return 0;
}

/*****************************************************************************
 *
 *  stats_colloid_force_io_free
 *
 *****************************************************************************/

__host__ int stats_colloid_force_io_free(stats_colloid_force_io_t ** pobj) {

  stats_colloid_force_io_t * obj = NULL;

  if (pobj == NULL || *pobj == NULL) return 0;

  obj = *pobj;

  if (pe_mpi_rank(obj->pe) == 0 && obj->fp) {
    fflush(obj->fp);
    fclose(obj->fp);
    obj->fp = NULL;
  }

  free(obj);
  *pobj = NULL;

  return 0;
}
