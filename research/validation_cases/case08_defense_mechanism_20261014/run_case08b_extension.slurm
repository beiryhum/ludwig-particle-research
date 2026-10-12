#!/bin/bash
#SBATCH --job-name=C08D_ext
#SBATCH --partition=v6_384
#SBATCH --nodes=1
#SBATCH --ntasks=96
#SBATCH --ntasks-per-node=96
#SBATCH --cpus-per-task=1
#SBATCH --mem=200G
#SBATCH --time=72:00:00
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err

set -euo pipefail
umask 027

CASE="${CASE:-}"
SEGMENT="${SEGMENT:-}"

case "$CASE" in
  B01|B02|B04) ;;
  *) echo "ERROR: CASE must be B01, B02, or B04"; exit 2 ;;
esac
case "$SEGMENT" in
  3|4|5) ;;
  *) echo "ERROR: SEGMENT must be 3, 4, or 5"; exit 2 ;;
esac

NSTART=$((600000 + SEGMENT * 200000))
NEND=$((NSTART + 200000))
STEP_START9="$(printf '%09d' "$NSTART")"
STEP_END9="$(printf '%09d' "$NEND")"
STEP_START8="$(printf '%08d' "$NSTART")"
STEP_END8="$(printf '%08d' "$NEND")"

PROJECT="/public1/home/sch10041/lzh/ludwig_t6"
BASE="$PROJECT/runs/validation/08_free_sphere_migration"
PACKAGE="$BASE/CASE08_DEFENSE_MECHANISM_20261014"
EXE="$PROJECT/build/prod/ludwig-0.23.1/src/Ludwig.exe"
EXPECTED_EXE_SHA256="f82b1b56db9aad42af16975f66d92ca761a1e41e98ef80c9921535235da03ecd"

case "$CASE" in
  B01) CASE_DIR="$BASE/08B01_from_A02_inner" ;;
  B02) CASE_DIR="$BASE/08B02_from_A03_outer" ;;
  B04) CASE_DIR="$BASE/08B04_from_A06_near_eq" ;;
esac

INPUT_TEMPLATE="$PACKAGE/inputs/input.${CASE}.segment${SEGMENT}"
CONFIG_START="$CASE_DIR/config.cds${STEP_START8}.001-001"
REPORT_PREFIX="CASE08_DEFENSE_SEGMENT${SEGMENT}_CHECK"
PREV=$((SEGMENT - 1))

source /public1/soft/modules/module.sh
module load mpi/oneAPI/2022.1

test -d "$CASE_DIR"
test -x "$EXE"
test -r "$INPUT_TEMPLATE"
test -r "$CONFIG_START"
test -f "$CASE_DIR/SEGMENT_${PREV}_COMPLETE"

actual_sha="$(sha256sum "$EXE" | awk '{print $1}')"
if [[ "$actual_sha" != "$EXPECTED_EXE_SHA256" ]]; then
  echo "ERROR: production executable SHA256 mismatch"
  echo "expected=$EXPECTED_EXE_SHA256"
  echo "actual=$actual_sha"
  exit 10
fi

for name in \
  "dist-${STEP_START9}.001-001" \
  "rho-${STEP_START9}.001-001" \
  "vel-${STEP_START9}.001-001" \
  dist-metadata.001-001 rho-metadata.001-001 vel-metadata.001-001; do
  test -s "$CASE_DIR/$name"
done

if [[ -e "$CASE_DIR/SEGMENT_${SEGMENT}_COMPLETE" ||
      -e "$CASE_DIR/dist-${STEP_END9}.001-001" ||
      -e "$CASE_DIR/rho-${STEP_END9}.001-001" ||
      -e "$CASE_DIR/vel-${STEP_END9}.001-001" ||
      -e "$CASE_DIR/config.cds${STEP_END8}.001-001" ]]; then
  echo "ERROR: segment end output already exists; refusing overwrite."
  exit 3
fi

cd "$CASE_DIR"
python3 tools/check_config08b.py --file "$(basename "$CONFIG_START")" --json "CONFIG_CHECK_START_DEFENSE_SEGMENT${SEGMENT}.json"
cp -f "$INPUT_TEMPLATE" input

export OMP_NUM_THREADS=1
export OMP_DYNAMIC=false
export OMP_PROC_BIND=close
export OMP_PLACES=cores
export I_MPI_PIN=0
export I_MPI_SHM=off

sha256sum "$INPUT_TEMPLATE" > "input.defense.segment${SEGMENT}.sha256"
sha256sum "$EXE" > executable.sha256
python3 --version > python_version.txt 2>&1
module list 2> "module_list_defense_segment${SEGMENT}.txt" || true

before_rows=0
if [[ -f colloid_force_torque.csv ]]; then before_rows="$(wc -l < colloid_force_torque.csv)"; fi
echo "$before_rows" > "force_rows_before_defense_segment${SEGMENT}.txt"

LOG="ludwig_${CASE}_defense_segment${SEGMENT}_${SLURM_JOB_ID}.log"
set +e
srun --cpu-bind=cores -n 96 "$EXE" 2>&1 | tee "$LOG"
run_status=${PIPESTATUS[0]}
set -e
echo "$run_status" > "defense_segment${SEGMENT}_srun_exit_status.txt"

if [[ "$run_status" -ne 0 ]]; then touch "DEFENSE_SEGMENT_${SEGMENT}_FAILED"; exit "$run_status"; fi
if ! grep -q "Ludwig finished normally." "$LOG"; then touch "DEFENSE_SEGMENT_${SEGMENT}_INCOMPLETE"; exit 4; fi
if grep -Ein '(^|[^[:alnum:]_])(FATAL|NaN|[-+]?Inf(inity)?)([^[:alnum:]_]|$)' "$LOG" >/dev/null; then touch "DEFENSE_SEGMENT_${SEGMENT}_NUMERICAL_ERROR"; exit 5; fi

for name in "dist-${STEP_END9}.001-001" "rho-${STEP_END9}.001-001" "vel-${STEP_END9}.001-001" "config.cds${STEP_END8}.001-001"; do test -s "$name"; done
test -s colloid_force_torque.csv

after_rows="$(wc -l < colloid_force_torque.csv)"
echo "$after_rows" > "force_rows_after_defense_segment${SEGMENT}.txt"
if [[ "$after_rows" -le "$before_rows" ]]; then echo "ERROR: force CSV did not gain rows."; exit 6; fi

set +e
python3 tools/check_case08b_segment.py --csv colloid_force_torque.csv --config "$(basename "$CONFIG_START")" --case "$CASE" --start "$NSTART" --end "$NEND" --frequency 10 --report-prefix "$REPORT_PREFIX" 2>&1 | tee "python_defense_segment${SEGMENT}_check.log"
check_status=${PIPESTATUS[0]}
set -e
echo "$check_status" > "defense_segment${SEGMENT}_python_check_status.txt"
if [[ "$check_status" -ne 0 ]]; then touch "DEFENSE_SEGMENT_${SEGMENT}_CHECK_FAIL"; exit "$check_status"; fi

touch "SEGMENT_${SEGMENT}_COMPLETE"
if [[ "$SEGMENT" -eq 5 ]]; then touch CASE08_DEFENSE_EXTENSION_COMPLETE; fi
find . -maxdepth 1 -type f -printf '%f\n' | sort > output_manifest.txt
echo "CASE08 DEFENSE ${CASE} SEGMENT ${SEGMENT}: COMPLETED (${NSTART}->${NEND})"
