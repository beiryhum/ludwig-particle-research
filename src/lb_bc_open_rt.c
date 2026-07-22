/*****************************************************************************
 *
 *  lb_bc_open_rt.c
 *
 *  LB open boundary condition initialisation.
 *
 *  Currently acceptable positions:
 *    We must have wall in two coordinate directions (3d)
 *    We must have no periodic conditions at all "periodicity 0_0_0"
 *    The inflow/outflow direction must be in the non-wall direction.
 *
 *  Edinburgh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *
 *  (c) 2021-2022 The University of Edinburgh
 *
 *  Contibuting authors:
 *  Kevin Stratford (kevin@epcc.ed.ac.uk)
 *
 *****************************************************************************/

#include <assert.h>
#include <stdio.h>
#include <string.h>

#include "lb_bc_open_rt.h"

/* Available implementations */

#include "lb_bc_inflow_rhou.h"
#include "lb_bc_outflow_rhou.h"

/*****************************************************************************
 *
 *  lb_bc_open_rt
 *
 *  As we are not generally going to have an inflow without an outflow
 *  (and vice-versa) these are together...
 *
 *****************************************************************************/

__host__ int lb_bc_open_rt(pe_t * pe, rt_t * rt, cs_t * cs, lb_t * lb,
			   lb_bc_open_t ** inflow,
			   lb_bc_open_t ** outflow) {

  int have_bc = 0;   /* inflow/outflow required */
  int flow[3] = {0}; /* flow direction */
  
  assert(pe);
  assert(rt);
  assert(cs);
  assert(lb);

  have_bc = rt_switch(rt, "lb_bc_open");

    if (have_bc) {
    int wall[3] = {0};
    int periodic[3] = {0};
    int nflow = 0;
    int id = -1;

    rt_int_parameter_vector(rt, "boundary_walls", wall);
    cs_periodic(cs, periodic);

    /*
     * Prefer an explicitly supplied open-flow direction.
     * Retain the old wall-based inference as a fallback for existing
     * fully walled channel inputs.
     */
    nflow = rt_int_parameter_vector(rt, "lb_bc_open_flow", flow);

    if (nflow == 0) {
      flow[X] = 1 - wall[X];
      flow[Y] = 1 - wall[Y];
      flow[Z] = 1 - wall[Z];
    }

    /*
     * Exactly one coordinate direction must be the open-flow direction.
     */
    if (lb_bc_inflow_opts_flow_valid(flow) == 0) {
      pe_fatal(pe,
               "Inflow/outflow requires exactly one open-flow direction\n");
    }

    for (int ia = 0; ia < 3; ia++) {
      if (flow[ia]) id = ia;
    }

    assert(id == X || id == Y || id == Z);

    /*
     * The streamwise direction must be open:
     * it cannot simultaneously be a wall or periodic.
     */
    if (wall[id]) {
      pe_fatal(pe,
               "Open-flow direction cannot also contain boundary walls\n");
    }

    if (periodic[id]) {
      pe_fatal(pe,
               "Open-flow direction cannot be periodic\n");
    }

    /*
     * Each transverse direction must be either:
     *
     *   (a) bounded by physical walls; or
     *   (b) periodic.
     *
     * It cannot be both and cannot be neither.
     */
    for (int ia = 0; ia < 3; ia++) {

      if (ia == id) continue;

      if (wall[ia] && periodic[ia]) {
        pe_fatal(pe,
                 "A transverse direction cannot be both wall and periodic\n");
      }

      if (!wall[ia] && !periodic[ia]) {
        pe_fatal(pe,
                 "Each transverse direction must be wall or periodic\n");
      }
    }

    rt_key_required(rt, "lb_bc_inflow_type",  RT_INFO);
    rt_key_required(rt, "lb_bc_outflow_type", RT_INFO);
  }

    /* X, Y, and Z flow directions are supported in serial and MPI. */

    /* Test periodicity. Must be (0,0,0). */

    
  
  /* Inflow */

  if (have_bc) {
 
    char intype[BUFSIZ] = {0};
    double u0[3] = {0};

    rt_string_parameter(rt, "lb_bc_inflow_type", intype, BUFSIZ/2);
    rt_double_parameter_vector(rt, "lb_bc_inflow_rhou_u0", u0);

    pe_info(pe, "\n");
    pe_info(pe, "Hydrodynamic open boundary condition for inflow\n");

    if (strncmp(intype, "rhou", BUFSIZ) == 0) {
      /* Give me a rhou inflow */
      lb_bc_inflow_rhou_t * rhou = NULL;
      lb_bc_inflow_opts_t options = {.nvel = lb->model.nvel,
	                             .flow = {flow[X], flow[Y], flow[Z]},
	                             .u0   = {u0[X],u0[Y],u0[Z]}};

      /* Check options are valid */
      if (lb_bc_inflow_opts_valid(options) == 0) {
	/* Print them out */
	pe_fatal(pe, "Please check inflow options and try again\n.");
      }

      lb_bc_inflow_rhou_create(pe, cs, &options, &rhou);
      *inflow = (lb_bc_open_t *) rhou;

      /* Might be nice to delegate this elsewhere ... */
      pe_info(pe, "Inflow type:              %s\n", "rhou");
      pe_info(pe, "Inflow flow profile:      %s\n", "uniform");
      pe_info(pe, "Inflow flow direction:    %d %d %d\n",
	      flow[X], flow[Y], flow[Z]);
      pe_info(pe, "Inflow flow value u0:    %14.7e %14.7e %14.7e\n",
	      u0[X], u0[Y], u0[Z]);
    }
    else {
      /* Not recognised */
      pe_fatal(pe, "lb_bc_inflow_type not recognised\n");
    }
  }

  /* Outflow */
  if (have_bc) {

    char outtype[BUFSIZ] = {0};
    double rho0 = 1.0;

    rt_string_parameter(rt, "lb_bc_outflow_type", outtype, BUFSIZ/2);
    rt_double_parameter(rt, "rho0", &rho0);
    rt_double_parameter(rt, "lb_bc_outflow_rhou_rho0", &rho0);

    pe_info(pe, "\n");
    pe_info(pe, "Hydrodynamic open boundary condition at outflow\n");
    
    if (strncmp(outtype, "rhou", BUFSIZ) == 0) {
      lb_bc_outflow_rhou_t * rhou = NULL;
      lb_bc_outflow_opts_t options = {.nvel = lb->model.nvel,
	                              .flow = {flow[X], flow[Y], flow[Z]},
	                              .rho0 = rho0};
      /* Check options valid */
      if (lb_bc_outflow_opts_valid(options) == 0) {
	/* Print out? */
	pe_fatal(pe, "Please check outflow options and try again\n");
      }

      lb_bc_outflow_rhou_create(pe, cs, &options, &rhou);
      *outflow = (lb_bc_open_t *) rhou;

      pe_info(pe, "Outflow type:             %s\n", "rhou");
      pe_info(pe, "Outflow flow direction:   %d %d %d\n",
	      flow[X], flow[Y], flow[Z]);
      pe_info(pe, "Outflow flow rho0:       %14.7e\n", rho0);
    }
    else {
      pe_fatal(pe, "lb_bc_outflow_type not recognised\n");
    }
  }

  return 0;
}
