#!/bin/bash
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE="/public1/home/sch10041/lzh/ludwig_t6/runs/validation/08_free_sphere_migration"
END="${END:-1400000}"
START="${START:-600000}"
OUT="$HERE/results_${END}"
mkdir -p "$OUT"

for case_name in B01 B02 B04; do
  case "$case_name" in
    B01) d="$BASE/08B01_from_A02_inner" ;;
    B02) d="$BASE/08B02_from_A03_outer" ;;
    B04) d="$BASE/08B04_from_A06_near_eq" ;;
  esac
  python3 "$HERE/tools/flow_factor_diagnostics.py" \
    --csv "$d/colloid_force_torque.csv" --case "$case_name" \
    --start "$START" --end "$END" --window 100000 \
    --output-prefix "$OUT/${case_name}_flow_factors"
done

python3 "$HERE/tools/check_defense_bracket.py" \
  --b01 "$BASE/08B01_from_A02_inner/colloid_force_torque.csv" \
  --b02 "$BASE/08B02_from_A03_outer/colloid_force_torque.csv" \
  --b04 "$BASE/08B04_from_A06_near_eq/colloid_force_torque.csv" \
  --end "$END" --tail-width 100000 \
  | tee "$OUT/DEFENSE_BRACKET_CHECK.txt"

echo "Postprocess outputs: $OUT"
