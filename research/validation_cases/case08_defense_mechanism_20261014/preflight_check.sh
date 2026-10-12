#!/bin/bash
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE="/public1/home/sch10041/lzh/ludwig_t6/runs/validation/08_free_sphere_migration"
EXE="/public1/home/sch10041/lzh/ludwig_t6/build/prod/ludwig-0.23.1/src/Ludwig.exe"
EXPECTED="f82b1b56db9aad42af16975f66d92ca761a1e41e98ef80c9921535235da03ecd"

bash -n "$HERE/run_case08b_extension.slurm"
bash -n "$HERE/submit_defense_extension.sh"
python3 -m py_compile "$HERE"/tools/*.py

actual="$(sha256sum "$EXE" | awk '{print $1}')"
[[ "$actual" == "$EXPECTED" ]] || { echo "EXE SHA mismatch"; exit 10; }

for case_name in B01 B02 B04; do
  case "$case_name" in
    B01) d="$BASE/08B01_from_A02_inner" ;;
    B02) d="$BASE/08B02_from_A03_outer" ;;
    B04) d="$BASE/08B04_from_A06_near_eq" ;;
  esac
  echo "Checking $case_name: $d"
  test -f "$d/SEGMENT_2_COMPLETE"
  test -s "$d/config.cds01200000.001-001"
  test -s "$d/dist-001200000.001-001"
  test -s "$d/rho-001200000.001-001"
  test -s "$d/vel-001200000.001-001"
  test -r "$d/tools/check_config08b.py"
  test -r "$d/tools/check_case08b_segment.py"
  for seg in 3 4 5; do
    test -r "$HERE/inputs/input.${case_name}.segment${seg}"
  done
done

echo "PREFLIGHT: PASS"
