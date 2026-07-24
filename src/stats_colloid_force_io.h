/*****************************************************************************
 *
 *  stats_colloid_force_io.h
 *
 *  Optional per-colloid hydrodynamic force and torque time-series output.
 *
 *****************************************************************************/

#ifndef LUDWIG_STATS_COLLOID_FORCE_IO_H
#define LUDWIG_STATS_COLLOID_FORCE_IO_H

#include "pe.h"
#include "runtime.h"
#include "colloids.h"

typedef struct stats_colloid_force_io_s stats_colloid_force_io_t;

__host__ int stats_colloid_force_io_create(pe_t * pe, rt_t * rt,
                                            colloids_info_t * cinfo,
                                            int nstart,
                                            int hydrodynamics,
                                            stats_colloid_force_io_t ** pobj);

__host__ int stats_colloid_force_io_write(stats_colloid_force_io_t * obj,
                                           int timestep);

__host__ int stats_colloid_force_io_free(stats_colloid_force_io_t ** pobj);

#endif
