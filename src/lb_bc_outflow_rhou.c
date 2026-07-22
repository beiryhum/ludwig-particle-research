/*****************************************************************************
 *
 *  lb_bc_outflow_rhou.c
 *
 *  Implementation of an outflow boundary condition which assumes
 *  a channel geometry with walls in two dimensions.
 *
 *  Hydrodynamic conditions at the outflow are, .e.g.,:
 *
 *    rho_outflow(x=L,y,z) = rho0
 *    u_outflow(x=L,y,z)   = u(L-1,y,z)
 *
 *  Incoming distributions are then set f_i = f^eq_i(rho, u) in
 *  the outflow boundary region. Outgoing distributions from the
 *  domain proper merely "disappear".
 *
 *
 *  Edinburgh Soft Matter and Statistical Phsyics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *  (c) 2021-2022 The University of Edinburgh
 *
 *****************************************************************************/

#include <assert.h>
#include <stdlib.h>

#include "lb_model.h"
#include "lb_bc_outflow_rhou.h"
#include "cs_limits.h"

typedef enum {LINK_COUNT, LINK_ASSIGN} link_init_enum_t;

__host__ int lb_bc_outflow_rhou_init_internal(lb_bc_outflow_rhou_t * outflow);
__host__ int lb_bc_outflow_init_link(lb_bc_outflow_rhou_t * outflow,
				     link_init_enum_t flag, int id);

static int int_max(int a, int b) {return (a > b) ? a : b;}

static const lb_bc_open_vtable_t vt_ = {
  (lb_bc_open_free_ft)   lb_bc_outflow_rhou_free,
  (lb_bc_open_update_ft) lb_bc_outflow_rhou_update,
  (lb_bc_open_impose_ft) lb_bc_outflow_rhou_impose,
  (lb_bc_open_stats_ft)  lb_bc_outflow_rhou_stats
};

/*****************************************************************************
 *
 *  lb_bc_outflow_rhou_create
 *
 *****************************************************************************/

__host__ int lb_bc_outflow_rhou_create(pe_t * pe, cs_t * cs,
				       const lb_bc_outflow_opts_t * options,
				       lb_bc_outflow_rhou_t ** outflow) {
  lb_bc_outflow_rhou_t * bc = NULL;

  assert(pe);
  assert(cs);
  assert(options);
  assert(outflow);

  bc = (lb_bc_outflow_rhou_t *) calloc(1, sizeof(lb_bc_outflow_rhou_t));
  assert(bc);
  if (bc == NULL) pe_fatal(pe, "Failed to allocate lb_bc_outflow_rhou_t");

  /* Pointers; superclass block */
  bc->pe = pe;
  bc->cs = cs;

  bc->super.func = &vt_;
  bc->super.id   = LB_BC_OUTFLOW_RHOU;

  if (!lb_bc_outflow_opts_valid(*options)) {
    /* Internal error if we reach this point. */
    pe_fatal(pe, "lb_bc_outflow_rhou_create: options invalid\n");
  }
  bc->options = *options;

  lb_bc_outflow_rhou_init_internal(bc);

  *outflow = bc;

  return 0;
}

/*****************************************************************************
 *
 *  lb_bc_outflow_rhou_init_internal
 *
 *****************************************************************************/

__host__ int lb_bc_outflow_rhou_init_internal(lb_bc_outflow_rhou_t * outflow) {

  assert(outflow);

  if (outflow->options.flow[X]) lb_bc_outflow_init_link(outflow, LINK_COUNT, X);
  if (outflow->options.flow[Y]) lb_bc_outflow_init_link(outflow, LINK_COUNT, Y);
  if (outflow->options.flow[Z]) lb_bc_outflow_init_link(outflow, LINK_COUNT, Z);

  {
    pe_t * pe = outflow->pe;
    int nlink = int_max(1, outflow->nlink); /* No zero sized allocations. */

    outflow->linkp = (int8_t *) calloc(nlink, sizeof(int8_t));
    outflow->linki = (int *)    calloc(nlink, sizeof(int));
    outflow->linkj = (int *)    calloc(nlink, sizeof(int));

    assert(outflow->linkp);
    assert(outflow->linki);
    assert(outflow->linkj);

    if (outflow->linkp == NULL) pe_fatal(pe, "calloc(outflow->linkp) NULL\n");
    if (outflow->linki == NULL) pe_fatal(pe, "calloc(outflow->linki) NULL\n");
    if (outflow->linkj == NULL) pe_fatal(pe, "calloc(outflow->linkj) NULL\n");
  }

  if (outflow->options.flow[X]) lb_bc_outflow_init_link(outflow, LINK_ASSIGN, X);
  if (outflow->options.flow[Y]) lb_bc_outflow_init_link(outflow, LINK_ASSIGN, Y);
  if (outflow->options.flow[Z]) lb_bc_outflow_init_link(outflow, LINK_ASSIGN, Z);

  return 0;
}

/*****************************************************************************
 *
 *  lb_bc_outflow_init_link
 *
 *  Identify links representing incoming distributions at the outflow
 *  which are fluid to fluid in the coordinate direction "id".
 *
 *  This assumes that the other two directions are walls, and so no
 *  duplicates of incoming solid-fluid (bbl) links are wanted.
 *
 *****************************************************************************/

__host__ int lb_bc_outflow_init_link(lb_bc_outflow_rhou_t * outflow,
				     link_init_enum_t init, int id) {

  cs_t * cs = NULL;
    int noffset[3] = {0};
  int ntotal[3] = {0};
  int nlocal[3] = {0};
  int periodic[3] = {0};
  int nlink = 0;

  assert(outflow);
  assert(id == X || id == Y || id == Z);

  cs = outflow->cs;

  cs_ntotal(cs, ntotal);
  cs_nlocal(cs, nlocal);
  cs_nlocal_offset(cs, noffset);
  cs_periodic(cs, periodic);

  if (noffset[id] + nlocal[id] == ntotal[id]) {
    int nmin[3]   = {1, 1, 1};

    lb_model_t model = {0};
    lb_model_create(outflow->options.nvel, &model);

    nmin[id] = nlocal[id]; /* Only rightmost edge in the relevant direction */

    for (int ic = nmin[X]; ic <= nlocal[X]; ic++) {
      for (int jc = nmin[Y]; jc <= nlocal[Y]; jc++) {
	for (int kc = nmin[Z]; kc <= nlocal[Z]; kc++) {

	  for (int p = 1; p < model.nvel; p++) {

	    if (model.cv[p][id] != +1) continue;

	    int ic1 = ic + model.cv[p][X];
	    int jc1 = jc + model.cv[p][Y];
	    int kc1 = kc + model.cv[p][Z];

	    /* Some shuffling to get "other 2 directions" */
	    int id1 = (id + 1) % 3;
	    int id2 = (id + 2) % 3;

	    int ijk[3] = {ic1, jc1, kc1};

	    	    /*
	     * Retain links crossing a periodic transverse boundary.
	     * Only links leaving a non-periodic physical direction are
	     * excluded.
	     */

	    if (!periodic[id1]) {
	      if (noffset[id1] + ijk[id1] < 1)           continue;
	      if (noffset[id1] + ijk[id1] > ntotal[id1]) continue;
	    }

	    if (!periodic[id2]) {
	      if (noffset[id2] + ijk[id2] < 1)           continue;
	      if (noffset[id2] + ijk[id2] > ntotal[id2]) continue;
	    }

	    if (init == LINK_ASSIGN) {
	      outflow->linkp[nlink] = model.nvel - p;
	      outflow->linki[nlink] = cs_index(cs, ic1, jc1, kc1);
	      outflow->linkj[nlink] = cs_index(cs, ic,  jc,  kc);
	    }
	    nlink += 1;
	  }
	}
      }
    }

    lb_model_free(&model);
  }

  outflow->nlink = nlink;

  return 0;
}

/*****************************************************************************
 *
 *  lb_bc_outflow_rhou_free
 *
 *****************************************************************************/

__host__ int lb_bc_outflow_rhou_free(lb_bc_outflow_rhou_t * outflow) {

  assert(outflow);

  free(outflow->linkj);
  free(outflow->linki);
  free(outflow->linkp);
  free(outflow);

  return 0;
}

/*****************************************************************************
 *
 *  lb_bc_outflow_rhou_update
 *
 *  Set any relevant hydrodynamic conditions at the outflow:
 *    (1) rho is constant as prescribed
 *    (2) u is set to the value relevant at outflow.
 *
 *****************************************************************************/

__host__ int lb_bc_outflow_rhou_update(lb_bc_outflow_rhou_t * outflow,
				      hydro_t * hydro) {

  cs_t * cs = NULL;
  int id = -1;
  int nlocal[3] = {0};
  int ntotal[3] = {0};
  int noffset[3] = {0};

  assert(outflow);
  assert(hydro);

  cs = outflow->cs;

  cs_ntotal(cs, ntotal);
  cs_nlocal(cs, nlocal);
  cs_nlocal_offset(cs, noffset);

  if (outflow->options.flow[X]) id = X;
  if (outflow->options.flow[Y]) id = Y;
  if (outflow->options.flow[Z]) id = Z;


  if (noffset[id] + nlocal[id] == ntotal[id]) {
    int nmin[3]   = {1,1,1};
    int idx = outflow->options.flow[X];
    int jdy = outflow->options.flow[Y];
    int kdz = outflow->options.flow[Z];

    nmin[id] = nlocal[id]; /* Only rightmost edge in the relevant direction */

    for (int ic = nmin[X]; ic <= nlocal[X]; ic++) {
      for (int jc = nmin[Y]; jc <= nlocal[Y]; jc++) {
	for (int kc = nmin[Z]; kc <= nlocal[Z]; kc++) {

	  int index0 = cs_index(cs, ic + idx, jc + jdy, kc + kdz);
	  int index1 = cs_index(cs, ic, jc, kc);
	  double u[3] = {0};

	  hydro_rho_set(hydro, index0, outflow->options.rho0);
	  hydro_u(hydro, index1, u);
	  hydro_u_set(hydro, index0, u);
	}
      }
    }
  }

  return 0;
}

/*****************************************************************************
 *
 *  lb_bc_outflow_rhou_halo
 *
 *  There is a lot in common with the inflow hlow which should be
 *  abstracted away as only dependent on hydro_t.
 *
 *  Apart from being at the other end of the system.
 *
 *****************************************************************************/

__host__ int lb_bc_outflow_rhou_halo_update(lb_bc_outflow_rhou_t * outflow,
                                           hydro_t * hydro) {
  cs_t * cs = NULL;
  int id = -1;
  int id1 = -1;
  int id2 = -1;
  int noffset[3] = {0};
  int nlocal[3]  = {0};
  int ntotal[3]  = {0};

  const int tag = 12348;

  assert(outflow);
  assert(hydro);

  cs = outflow->cs;
  cs_nlocal(cs, nlocal);
  cs_nlocal_offset(cs, noffset);
  cs_ntotal(cs, ntotal);

  if (outflow->options.flow[X]) id = X;
  if (outflow->options.flow[Y]) id = Y;
  if (outflow->options.flow[Z]) id = Z;

  assert(id == X || id == Y || id == Z);

  id1 = (id + 1) % 3;
  id2 = (id + 2) % 3;

  if (noffset[id] + nlocal[id] == ntotal[id]) {
    MPI_Comm comm = MPI_COMM_NULL;
    MPI_Request req[8] = {0};
    double * buf[8] = {0};
    int ngbr[8] = {0};
    int count[8] = {0};
    const int nhalo = 1;
    const int count1 = 4*nlocal[id2]*nhalo;
    const int count2 = 4*nlocal[id1]*nhalo;

    cs_cart_comm(cs, &comm);

    /* Receive from the two transverse coordinate directions. */
    ngbr[0] = cs_cart_neighb(cs, FORWARD,  id1);
    ngbr[1] = cs_cart_neighb(cs, BACKWARD, id1);
    ngbr[2] = cs_cart_neighb(cs, FORWARD,  id2);
    ngbr[3] = cs_cart_neighb(cs, BACKWARD, id2);

    /* Send to the opposite neighbours. */
    ngbr[4] = cs_cart_neighb(cs, BACKWARD, id1);
    ngbr[5] = cs_cart_neighb(cs, FORWARD,  id1);
    ngbr[6] = cs_cart_neighb(cs, BACKWARD, id2);
    ngbr[7] = cs_cart_neighb(cs, FORWARD,  id2);

    count[0] = count1;
    count[1] = count1;
    count[2] = count2;
    count[3] = count2;
    count[4] = count1;
    count[5] = count1;
    count[6] = count2;
    count[7] = count2;

    for (int m = 0; m < 8; m++) {
      pe_t * pe = outflow->pe;

      req[m] = MPI_REQUEST_NULL;
      buf[m] = (double *) malloc(sizeof(double)*count[m]);
      assert(buf[m]);
      if (buf[m] == NULL) pe_fatal(pe, "malloc() failed\n");
    }

    for (int m = 0; m < 4; m++) {
      MPI_Irecv(buf[m], count[m], MPI_DOUBLE, ngbr[m], tag, comm, req + m);
    }

    /* Pack the open-boundary halo plane and post the four sends. */
    for (int ms = 0; ms < 4; ms++) {
      int m = 4 + ms;
      int ib = 0;
      int lo1 = 1;
      int hi1 = nlocal[id1];
      int lo2 = 1;
      int hi2 = nlocal[id2];

      if (m == 4) {
        lo1 = 1;
        hi1 = nhalo;
      }
      if (m == 5) {
        lo1 = nlocal[id1] - nhalo + 1;
        hi1 = nlocal[id1];
      }
      if (m == 6) {
        lo2 = 1;
        hi2 = nhalo;
      }
      if (m == 7) {
        lo2 = nlocal[id2] - nhalo + 1;
        hi2 = nlocal[id2];
      }

      for (int i1 = lo1; i1 <= hi1; i1++) {
        for (int i2 = lo2; i2 <= hi2; i2++) {
          int ijk[3] = {0};
          int index = -1;
          double rho = 0.0;
          double u[3] = {0};

          ijk[id]  = nlocal[id] + 1;
          ijk[id1] = i1;
          ijk[id2] = i2;
          index = cs_index(cs, ijk[X], ijk[Y], ijk[Z]);

          hydro_rho(hydro, index, &rho);
          hydro_u(hydro, index, u);

          buf[m][ib++] = rho;
          buf[m][ib++] = u[X];
          buf[m][ib++] = u[Y];
          buf[m][ib++] = u[Z];
        }
      }

      assert(ib == count[m]);
      MPI_Isend(buf[m], count[m], MPI_DOUBLE, ngbr[m], tag, comm, req + m);
    }

    for (int ms = 0; ms < 8; ms++) {
      int m = -1;
      MPI_Status status = {0};

      MPI_Waitany(8, req, &m, &status);
      if (m == MPI_UNDEFINED) continue;
      if (status.MPI_SOURCE == MPI_PROC_NULL) continue;

      if (m >= 4) {
        free(buf[m]);
        buf[m] = NULL;
      }
      else {
        int ib = 0;
        int lo1 = 1;
        int hi1 = nlocal[id1];
        int lo2 = 1;
        int hi2 = nlocal[id2];

        if (m == 0) {
          lo1 = nlocal[id1] + 1;
          hi1 = nlocal[id1] + nhalo;
        }
        if (m == 1) {
          lo1 = 1 - nhalo;
          hi1 = 0;
        }
        if (m == 2) {
          lo2 = nlocal[id2] + 1;
          hi2 = nlocal[id2] + nhalo;
        }
        if (m == 3) {
          lo2 = 1 - nhalo;
          hi2 = 0;
        }

        for (int i1 = lo1; i1 <= hi1; i1++) {
          for (int i2 = lo2; i2 <= hi2; i2++) {
            int ijk[3] = {0};
            int index = -1;
            double rho = 0.0;
            double u[3] = {0};

            ijk[id]  = nlocal[id] + 1;
            ijk[id1] = i1;
            ijk[id2] = i2;
            index = cs_index(cs, ijk[X], ijk[Y], ijk[Z]);

            rho  = buf[m][ib++];
            u[X] = buf[m][ib++];
            u[Y] = buf[m][ib++];
            u[Z] = buf[m][ib++];

            hydro_rho_set(hydro, index, rho);
            hydro_u_set(hydro, index, u);
          }
        }

        assert(ib == count[m]);
        free(buf[m]);
        buf[m] = NULL;
      }
    }

    for (int m = 0; m < 8; m++) {
      free(buf[m]);
    }
  }

  return 0;
}

/*****************************************************************************
 *
 *  lb_bc_outflow_rhou_impose
 *
 *  Very largely the same as the inflow routine; there could be
 *  some consolidation.
 *  Intent: After lattice halo swap; before propagation.
 *
 *****************************************************************************/

__host__ int lb_bc_outflow_rhou_impose(lb_bc_outflow_rhou_t * outflow,
				       hydro_t * hydro,
				       lb_t * lb) {
  assert(outflow);
  assert(hydro);
  assert(lb);

  /*
   * The boundary update sets rho and u on the normal outflow halo plane.
   * Exchange the transverse parts of that plane before reconstructing
   * incoming distributions at linki.
   */
  if (pe_mpi_size(outflow->pe) > 1) {
    lb_bc_outflow_rhou_halo_update(outflow, hydro);
  }

  for (int n = 0; n < outflow->nlink; n++) {

    int index = outflow->linki[n];
    int8_t p  = outflow->linkp[n];

    double rho = 0.0;
    double u[3] = {0.0, 0.0, 0.0};

    hydro_rho(hydro, index, &rho);
    hydro_u(hydro, index, u);

    {
      double cs2   = lb->model.cs2;
      double rcs2  = 1.0/cs2;
      double udotc = 0.0;
      double sdotq = 0.0;
      double fp    = 0.0;

      for (int ia = 0; ia < 3; ia++) {
	udotc += u[ia]*lb->model.cv[p][ia];

	for (int ib = 0; ib < 3; ib++) {
	  double d_ab = (ia == ib);
	  double s_ab =
	    lb->model.cv[p][ia]*lb->model.cv[p][ib] - cs2*d_ab;

	  sdotq += s_ab*u[ia]*u[ib];
	}
      }

      fp = rho*lb->model.wv[p]*
	(1.0 + rcs2*udotc + 0.5*rcs2*rcs2*sdotq);

      lb_f_set(lb, index, p, LB_RHO, fp);
    }
  }

  return 0;
}

/*****************************************************************************
 *
 *  lb_bc_outflow_rhou_stats
 *
 *  Placeholder. No operation at the moment.
 *
 *****************************************************************************/

__host__ int lb_bc_outflow_rhou_stats(lb_bc_outflow_rhou_t * outflow) {

  return 0;
}
