#!/usr/bin/env python3
from __future__ import print_function

import argparse
import io
import os

CASES = {
    'B01': ('08B01_from_A02_inner', 'from the inner side of the candidate stable region'),
    'B02': ('08B02_from_A03_outer', 'from the outer side of the candidate stable region'),
    'B04': ('08B04_from_A06_near_eq', 'near the candidate stable region'),
}


def render_input(case, segment):
    if case not in CASES:
        raise ValueError('case must be one of B01, B02, B04')
    if segment not in (3, 4, 5):
        raise ValueError('segment must be 3, 4, or 5')
    directory, purpose = CASES[case]
    nstart = 600000 + segment * 200000
    nend = nstart + 200000
    return '''# Case08 defense mechanism extension: free wall-normal migration
# Case: {case}
# Directory: {directory}
# Purpose: {purpose}
# Run: defense extension segment{segment} {nstart}->{nend}
# x=streamwise, y=wall-normal, z=spanwise
# Restart config must retain isfixedrxyz=0_0_1, isfixedvxyz=0_0_1, isfixedw=0.
# Long migration run: exact per-link BBL output remains OFF.

N_start {nstart}
N_cycles 200000

size 480_192_192
grid 8_4_3
periodicity 1_0_1
reorder 1

free_energy none
hydrodynamics on

fluid_rho0 1.0
viscosity 0.024000
viscosity_bulk 0.024000

lb_relaxation_scheme M10
lb_fluctuations no
lb_ghost_modes yes

boundary_walls 0_1_0
boundary_speed_bottom 0.0_0.0_0.0
boundary_speed_top 0.0_0.0_0.0

force 1.0416666666666667e-07_0.0_0.0

colloid_init from_file
colloid_rho0 1.0
colloid_rebuild_freq 1

colloid_force_output yes
colloid_force_index 1
colloid_force_freq 10
colloid_force_flush_freq 100
colloid_force_filename colloid_force_torque.csv

colloid_bbl_link_output no

colloid_io_grid 1_1_1
colloid_io_format ASCII
colloid_io_freq 200000

default_io_mode mpiio
default_io_format binary
default_io_grid 1_1_1
lb_io_grid 1_1_1
rho_io_grid 1_1_1
vel_io_grid 1_1_1

freq_statistics 5000
freq_measure 0
freq_config 200000
rho_io_freq 200000
vel_io_freq 200000

config_at_end yes
stats_vel_print_vol_flux yes
'''.format(case=case, directory=directory, purpose=purpose, segment=segment, nstart=nstart, nend=nend)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--output-dir', required=True)
    args = p.parse_args()
    if not os.path.isdir(args.output_dir):
        os.makedirs(args.output_dir)
    for case in sorted(CASES):
        for segment in (3, 4, 5):
            path = os.path.join(args.output_dir, 'input.{0}.segment{1}'.format(case, segment))
            with io.open(path, 'w', encoding='utf-8', newline='\n') as handle:
                handle.write(render_input(case, segment))
            print(path)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
