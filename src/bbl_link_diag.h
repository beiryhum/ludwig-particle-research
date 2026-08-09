/*****************************************************************************
 *
 *  bbl_link_diag.h
 *
 *  Read-only exact moving-BBL link and implicit-force diagnostics.
 *
 *****************************************************************************/

#ifndef LUDWIG_BBL_LINK_DIAG_H
#define LUDWIG_BBL_LINK_DIAG_H

#include "coords.h"
#include "colloids.h"
#include "pe.h"
#include "runtime.h"

typedef struct bbl_link_diag_s bbl_link_diag_t;
typedef struct bbl_link_record_s bbl_link_record_t;

typedef enum bbl_link_diag_pass1_component_enum {
  BBL_LINK_DIAG_PASS1_FLUID = 0,
  BBL_LINK_DIAG_PASS1_COLLOID,
  BBL_LINK_DIAG_PASS1_BOUNDARY,
  BBL_LINK_DIAG_PASS1_ACTIVE
} bbl_link_diag_pass1_component_t;

struct bbl_link_record_s {
  int step;
  int particle_index;
  int owner_rank;
  int particle_copy_is_local;

  int fluid_site[3];
  int solid_site[3];
  int p;
  int c[3];

  double rb[3];
  double v[3];
  double w[3];
  double ub[3];

  double fdist_in;
  double dm;
  double df;
  double deltam;
  double dms;
  double sump;

  double dF[3];
  double dT[3];
};

int bbl_link_diag_create(pe_t * pe, cs_t * cs, rt_t * rt,
                         colloids_info_t * cinfo, int nstart,
                         bbl_link_diag_t ** pdiag);
int bbl_link_diag_free(bbl_link_diag_t ** pdiag);

/* Called once near the beginning of each LB time step, before colloid rebuild. */
int bbl_link_diag_step_set(bbl_link_diag_t * diag, int step);

int bbl_link_diag_active_for(const bbl_link_diag_t * diag,
                             int particle_index);

/* Exact source hook for momentum added to f0/t0 while lattice sites are
 * swallowed or exposed during particle-map rebuilding. A global registration
 * is used so build.c need not acquire or own the bbl object. */
int bbl_link_diag_build_source_add(cs_t * cs, int particle_index,
                                   int site_index, const double force[3],
                                   const double torque[3]);

/* Pass-1 decomposition. Caller must apply the unique lattice-site owner rule. */
int bbl_link_diag_pass1_component_add(
    bbl_link_diag_t * diag, int particle_index,
    bbl_link_diag_pass1_component_t component,
    const double force[3], const double torque[3]);

/* Capture Ludwig's actual post-halo pass-1 f0/t0 before the 6x6 update. */
int bbl_link_diag_pass1_total_capture(bbl_link_diag_t * diag,
                                      colloids_info_t * cinfo);

/* Independent pass-2 reference accumulator. Caller must apply the same unique
 * outside-fluid-site owner rule used by the physical-link CSV path. */
int bbl_link_diag_pass2_reference_add(bbl_link_diag_t * diag,
                                      int particle_index,
                                      const double force[3],
                                      const double torque[3]);

int bbl_link_diag_record(bbl_link_diag_t * diag,
                         const bbl_link_record_t * record);
int bbl_link_diag_end(bbl_link_diag_t * diag, colloids_info_t * cinfo);

#endif
