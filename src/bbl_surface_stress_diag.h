/*****************************************************************************
 *
 *  bbl_surface_stress_diag.h
 *
 *  Read-only local distribution / momentum-flux diagnostic on physical
 *  colloid BBL links. Diagnostic only; no BBL or colloid dynamics changes.
 *
 *****************************************************************************/

#ifndef LUDWIG_BBL_SURFACE_STRESS_DIAG_H
#define LUDWIG_BBL_SURFACE_STRESS_DIAG_H

#include "coords.h"
#include "colloid_link.h"
#include "colloids.h"
#include "hydro.h"
#include "lb_data.h"
#include "pe.h"
#include "runtime.h"

typedef struct bbl_surface_stress_diag_s bbl_surface_stress_diag_t;

int bbl_surface_stress_diag_create(pe_t * pe, cs_t * cs, rt_t * rt,
                                   colloids_info_t * cinfo, lb_t * lb,
                                   hydro_t * hydro, int nstart,
                                   int have_viscosity_model,
                                   bbl_surface_stress_diag_t ** pdiag);
int bbl_surface_stress_diag_free(bbl_surface_stress_diag_t ** pdiag);

/* Must be called at the beginning of each LB time step, before colloid rebuild. */
int bbl_surface_stress_diag_step_set(bbl_surface_stress_diag_t * diag, int step);

/* Read-only pre-collision snapshot. This is called after all fluid force and
 * local-viscosity updates, immediately before lb_collide(). */
int bbl_surface_stress_diag_pre_collision(bbl_surface_stress_diag_t * diag);

/* Synchronise current post-collision hydro fields to host before host-side BBL.
 * No physical data are modified. */
int bbl_surface_stress_diag_prepare_pass2(bbl_surface_stress_diag_t * diag);

/* Called from bbl_pass2() for LINK_FLUID before the BBL distribution write.
 * A process-global registration is used only to avoid changing the BBL API. */
int bbl_surface_stress_diag_pass2_link(lb_t * lb, colloid_t * pc,
                                       colloid_link_t * link,
                                       int fluid_index, int solid_index);

/* Called once after bbl_pass2() has visited all colloids/links. */
int bbl_surface_stress_diag_end_step(void);

#endif
