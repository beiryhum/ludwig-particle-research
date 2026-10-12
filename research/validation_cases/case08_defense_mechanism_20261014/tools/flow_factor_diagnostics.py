#!/usr/bin/env python3
from __future__ import print_function

import argparse
import csv
import io
import math

Y_CENTER = 96.5
HALF_HEIGHT = 96.0
RADIUS = 12.0
UMAX = 0.02


def eta_from_y(y, y_center=Y_CENTER, half_height=HALF_HEIGHT):
    return (float(y) - y_center) / half_height


def u0_from_eta(eta, umax=UMAX):
    return umax * (1.0 - float(eta) ** 2)


def du0_dy(eta, umax=UMAX, half_height=HALF_HEIGHT):
    return -2.0 * umax * float(eta) / half_height


def d2u0_dy2(umax=UMAX, half_height=HALF_HEIGHT):
    return -2.0 * umax / (half_height ** 2)


def wall_gap_from_y(y, y_center=Y_CENTER, half_height=HALF_HEIGHT, radius=RADIUS):
    return half_height - abs(float(y) - y_center) - radius


def fluid_vorticity_z(eta):
    return -du0_dy(eta)


def angular_slip_z(wz, eta):
    return float(wz) - 0.5 * fluid_vorticity_z(eta)


def mean(values):
    return sum(values) / float(len(values)) if values else float('nan')


def population_std(values):
    if not values:
        return float('nan')
    m = mean(values)
    return math.sqrt(sum((v - m) ** 2 for v in values) / float(len(values)))


def parse_args():
    p = argparse.ArgumentParser(description='Compute mechanism-oriented flow factors from Ludwig colloid_force_torque.csv')
    p.add_argument('--csv', required=True)
    p.add_argument('--case', required=True)
    p.add_argument('--start', type=int, required=True)
    p.add_argument('--end', type=int, required=True)
    p.add_argument('--window', type=int, default=100000)
    p.add_argument('--output-prefix', required=True)
    return p.parse_args()


def main():
    args = parse_args()
    rows = []
    with io.open(args.csv, 'r', encoding='utf-8', newline='') as handle:
        reader = csv.DictReader(handle)
        for raw in reader:
            if raw.get('step') == 'step':
                continue
            try:
                step = int(raw['step'])
            except Exception:
                continue
            if not (args.start < step <= args.end):
                continue
            y = float(raw['y'])
            eta = eta_from_y(y)
            vx = float(raw['vx'])
            rec = {
                'case': args.case,
                'step': step,
                'eta': eta,
                'wall_gap': wall_gap_from_y(y),
                'vx': vx,
                'vy': float(raw['vy']),
                'wz': float(raw['wz']),
                'Fy_hydro': float(raw['Fy_hydro']),
                'u0_background': u0_from_eta(eta),
                'du0_dy': du0_dy(eta),
                'd2u0_dy2': d2u0_dy2(),
                'slip_x': vx - u0_from_eta(eta),
                'fluid_vorticity_z': fluid_vorticity_z(eta),
                'angular_slip_z': angular_slip_z(float(raw['wz']), eta),
            }
            rows.append(rec)

    fields = ['case','step','eta','wall_gap','vx','vy','wz','Fy_hydro','u0_background','du0_dy','d2u0_dy2','slip_x','fluid_vorticity_z','angular_slip_z']
    with io.open(args.output_prefix + '_per_step.csv', 'w', encoding='utf-8', newline='') as handle:
        w = csv.DictWriter(handle, fieldnames=fields)
        w.writeheader()
        for rec in rows:
            w.writerow(rec)

    summary_fields = ['case','window_start','window_end','row_count','mean_eta','std_eta','mean_wall_gap','mean_Fy_hydro','std_Fy_hydro','mean_vx','mean_vy','mean_wz','mean_u0_background','mean_du0_dy','mean_d2u0_dy2','mean_slip_x','mean_fluid_vorticity_z','mean_angular_slip_z']
    summaries = []
    ws = args.start
    while ws < args.end:
        we = min(ws + args.window, args.end)
        chunk = [r for r in rows if ws < r['step'] <= we]
        if chunk:
            def vals(k): return [r[k] for r in chunk]
            summaries.append({
                'case': args.case,
                'window_start': ws,
                'window_end': we,
                'row_count': len(chunk),
                'mean_eta': mean(vals('eta')),
                'std_eta': population_std(vals('eta')),
                'mean_wall_gap': mean(vals('wall_gap')),
                'mean_Fy_hydro': mean(vals('Fy_hydro')),
                'std_Fy_hydro': population_std(vals('Fy_hydro')),
                'mean_vx': mean(vals('vx')),
                'mean_vy': mean(vals('vy')),
                'mean_wz': mean(vals('wz')),
                'mean_u0_background': mean(vals('u0_background')),
                'mean_du0_dy': mean(vals('du0_dy')),
                'mean_d2u0_dy2': mean(vals('d2u0_dy2')),
                'mean_slip_x': mean(vals('slip_x')),
                'mean_fluid_vorticity_z': mean(vals('fluid_vorticity_z')),
                'mean_angular_slip_z': mean(vals('angular_slip_z')),
            })
        ws = we

    with io.open(args.output_prefix + '_window_summary.csv', 'w', encoding='utf-8', newline='') as handle:
        w = csv.DictWriter(handle, fieldnames=summary_fields)
        w.writeheader()
        for rec in summaries:
            w.writerow(rec)

    print('case={0} selected_rows={1} windows={2}'.format(args.case, len(rows), len(summaries)))
    if rows:
        print('eta_first={0:.10f} eta_last={1:.10f}'.format(rows[0]['eta'], rows[-1]['eta']))
        print('output_prefix={0}'.format(args.output_prefix))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
