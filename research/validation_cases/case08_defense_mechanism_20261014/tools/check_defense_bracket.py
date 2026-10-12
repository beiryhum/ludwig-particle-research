#!/usr/bin/env python3
from __future__ import print_function
import argparse, csv, io

CENTER=96.5
HALF=96.0

def eta(y): return (float(y)-CENTER)/HALF

def mean(v): return sum(v)/float(len(v)) if v else float('nan')

def tail(path, end, width):
    vals=[]
    with io.open(path,'r',encoding='utf-8',newline='') as h:
        for r in csv.DictReader(h):
            if r.get('step')=='step': continue
            try: s=int(r['step'])
            except Exception: continue
            if end-width < s <= end:
                vals.append((eta(r['y']), float(r['vy']), float(r['Fy_hydro'])))
    return vals

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--b01',required=True); p.add_argument('--b02',required=True); p.add_argument('--b04',required=True)
    p.add_argument('--end',type=int,required=True); p.add_argument('--tail-width',type=int,default=100000)
    p.add_argument('--separation-tolerance',type=float,default=0.02)
    a=p.parse_args()
    data={k:tail(v,a.end,a.tail_width) for k,v in [('B01',a.b01),('B02',a.b02),('B04',a.b04)]}
    for k in ('B01','B02','B04'):
        vals=data[k]
        if not vals: raise SystemExit('no tail rows for '+k)
        print('{0}: tail_mean_eta={1:.10f} tail_mean_vy={2:.6e} tail_mean_Fy={3:.6e} rows={4}'.format(k,mean([x[0] for x in vals]),mean([x[1] for x in vals]),mean([x[2] for x in vals]),len(vals)))
    e1=mean([x[0] for x in data['B01']]); e2=mean([x[0] for x in data['B02']])
    sep=abs(e2-e1)
    print('B01_B02_tail_eta_separation={0:.10f}'.format(sep))
    print('bracket_tolerance={0:.10f}'.format(a.separation_tolerance))
    print('BRACKET_TOLERANCE_PASS={0}'.format('YES' if sep <= a.separation_tolerance else 'NO'))
    print('NOTE: bracket tolerance pass is not by itself publication-grade final-equilibrium convergence.')
    return 0
if __name__=='__main__': raise SystemExit(main())
