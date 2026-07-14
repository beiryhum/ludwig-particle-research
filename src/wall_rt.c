/*****************************************************************************
 *
 *  wall_rt.c
 *
 *  Edinburgh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *  (c) 2015-2022 The University of Edinburgh
 *
 *  Contributing authors:
 *  Kevin Stratford (kevin@epcc.ed.ac.uk)
 *
 *****************************************************************************/

#include <assert.h>

#include "wall_rt.h"

/*****************************************************************************
 *
 *  wall_rt_init
 *
 ****************************************************************************/

int wall_rt_init(pe_t * pe, cs_t * cs, rt_t * rt, lb_t * lb, map_t * map,
		 wall_t ** wall) {
  int i;
  wall_slip_t ws = {0};
  wall_param_t p = {0};

  assert(rt);
  assert(cs);
  assert(lb);
  assert(map);
  assert(wall);

  /* We may have porous media, or wall boundaries, but not usually
     both. Some care may be required. */

  map_pm(map, &p.isporousmedia);
  rt_int_parameter_vector(rt, "boundary_walls", &p.isboundary[X]);

  p.iswall = (p.isboundary[X] || p.isboundary[Y] || p.isboundary[Z]);

  /* Run through input parameters */

if (p.iswall) {
    /* 1. 初始化速度矢量为 0 */
    for (i = 0; i < 3; i++) {
      p.ubot[i] = 0.0;
      p.utop[i] = 0.0;
      p.osc_A_bot[i] = 0.0; /* 新增初始化 */
      p.osc_A_top[i] = 0.0;
      p.osc_T_bot[i] = 0.0;
      p.osc_T_top[i] = 0.0;
    }

    /* 2. 使用矢量读取：允许用户在 input 中设置三个方向的速度分量 */
    rt_double_parameter_vector(rt, "boundary_speed_bottom", p.ubot);
    rt_double_parameter_vector(rt, "boundary_speed_top",    p.utop);
    /* --- 新增：读取输入文件中的震荡参数 --- */
    rt_double_parameter_vector(rt, "boundary_osc_amp_bottom", p.osc_A_bot);
    rt_double_parameter_vector(rt, "boundary_osc_amp_top",    p.osc_A_top);
    rt_double_parameter_vector(rt, "boundary_osc_period_bottom", p.osc_T_bot);
    rt_double_parameter_vector(rt, "boundary_osc_period_top",    p.osc_T_top);

    /* 3. 通用物理检查：如果某方向有墙，则该方向的速度分量必须为 0 (不能穿墙) */
    for (i = 0; i < 3; i++) {
      if (p.isboundary[i]) {
        if (p.ubot[i] != 0.0 || p.utop[i] != 0.0) {
          pe_info(pe, "Error: Wall normal to direction %d cannot have a normal velocity!\n", i);
          pe_fatal(pe, "Please ensure velocity vectors are tangential to the wall.\n");
        }
      }
    }

    /* 4. 读取润滑力参数 (各向同性应用) */
    rt_double_parameter(rt, "boundary_lubrication_rcnormal", &p.lubr_rc[X]);
    p.lubr_rc[Y] = p.lubr_rc[X];
    p.lubr_rc[Z] = p.lubr_rc[X];

    rt_double_parameter(rt, "boundary_lubrication_dhnormal", &p.lubr_dh[X]);
    p.lubr_dh[Y] = p.lubr_dh[X];
    p.lubr_dh[Z] = p.lubr_dh[X];
  }
  /* Slip properties [optional] */

  if (p.iswall) {
    double sbot[3] = {0.0, 0.0, 0.0};
    double stop[3] = {0.0, 0.0, 0.0};

    rt_double_parameter_vector(rt, "boundary_walls_slip_fraction_bot", sbot);
    rt_double_parameter_vector(rt, "boundary_walls_slip_fraction_top", stop);

    ws = wall_slip(sbot, stop);
    if (wall_slip_valid(&ws) == 0) {
      /* Could move elsewhere */
      pe_info(pe,  "Wall slip parameters must be 0 <= s <= 1 everywhere\n");
      pe_fatal(pe, "Please check and try again!\n");
    }

    /* Allow user to force use of slip (even if s = 0 everywhere) */
    if (rt_switch(rt, "boundary_walls_slip_active")) ws.active = 1;

    p.slip = ws;
  }

  /* Allocate */

  wall_create(pe, cs, map, lb, wall);

  if (p.iswall || p.isporousmedia) {
    wall_commit(*wall, &p);
    wall_info(*wall);
  }

  return 0;
}
