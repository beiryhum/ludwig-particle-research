#!/bin/bash
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

LAST_SEGMENT="${LAST_SEGMENT:-5}"
case "$LAST_SEGMENT" in
  3|4|5) ;;
  *) echo "ERROR: LAST_SEGMENT must be 3, 4, or 5"; exit 2 ;;
esac

for case_name in B01 B02 B04; do
  echo "===== sbatch --test-only CASE=${case_name} SEGMENT=3 ====="
  sbatch --test-only --export=ALL,CASE="$case_name",SEGMENT=3 run_case08b_extension.slurm
done

echo "All Slurm test-only checks passed."

submit_chain() {
  case_name="$1"
  jid="$(sbatch --parsable --export=ALL,CASE="$case_name",SEGMENT=3 run_case08b_extension.slurm)"
  echo "$case_name segment3 job=$jid"
  seg=4
  while [[ "$seg" -le "$LAST_SEGMENT" ]]; do
    jid="$(sbatch --parsable --dependency=afterok:"$jid" --export=ALL,CASE="$case_name",SEGMENT="$seg" run_case08b_extension.slurm)"
    echo "$case_name segment${seg} job=$jid"
    seg=$((seg + 1))
  done
}

submit_chain B01
submit_chain B02
submit_chain B04

echo
squeue -u "$USER" -o "%.18i %.20j %.10T %.30R %.19S" || true
