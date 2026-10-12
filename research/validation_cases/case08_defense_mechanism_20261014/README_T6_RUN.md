# CASE08_DEFENSE_MECHANISM_20261014 — T6 run package

## Purpose
Defense-priority continuation of Case08B free-Y migration for B01/B02/B04, plus mechanism-oriented flow-factor diagnostics.

Primary question: which flow conditions co-vary with the sign and magnitude of the lateral hydrodynamic force during migration?

## Frozen physical setup
- Ludwig executable: `/public1/home/sch10041/lzh/ludwig_t6/build/prod/ludwig-0.23.1/src/Ludwig.exe`
- expected SHA256: `f82b1b56db9aad42af16975f66d92ca761a1e41e98ef80c9921535235da03ecd`
- lattice: `480 x 192 x 192`
- MPI decomposition: `8 x 4 x 3 = 96 ranks`
- `nu = 0.024`
- body force: `1.0416666666666667e-07` in +x
- wall-normal coordinate: `eta=(y-96.5)/96`
- sphere radius: `a=12`
- free X/Y translation, fixed Z; free rotation.

## Extension segments
- segment3: 1,200,000 -> 1,400,000
- segment4: 1,400,000 -> 1,600,000
- segment5: 1,600,000 -> 1,800,000

Long runs keep `colloid_bbl_link_output no`. Exact per-link BBL remains a separate short diagnostic job.

## Generate inputs
```bash
python3 tools/generate_extension_inputs.py --output-dir inputs
```

## T6 install
Copy the complete archived package `CASE08_DEFENSE_MECHANISM_20261014.tar.gz` to:
`/public1/home/sch10041/lzh/ludwig_t6/runs/validation/08_free_sphere_migration/`

Then:
```bash
chmod +x *.sh *.slurm tools/*.py
./preflight_check.sh
```

## Submit
```bash
LAST_SEGMENT=3 ./submit_defense_extension.sh   # first 1.4M checkpoint
LAST_SEGMENT=4 ./submit_defense_extension.sh   # through 1.6M
./submit_defense_extension.sh                  # through 1.8M
```

## Post-process
```bash
END=1400000 ./postprocess_defense.sh
```

Mechanism diagnostics report eta, wall gap, full Ludwig `Fy_hydro`, ideal background Poiseuille velocity/shear/curvature, streamwise particle-background slip, fluid vorticity and angular slip.

These are diagnostic associations; they do not by themselves prove a unique causal lift decomposition.

## Local-load resolution
15 degrees is the primary quantitative scale. 7.5 degrees remains discretization sensitivity/SI only.

## Semantic boundary
Exact moving-BBL local output is the current pass-2 physical momentum-exchange component. Its integral closes to the pass-2 physical reference, not directly to full Ludwig `Fy_hydro`.
