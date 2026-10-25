#!/usr/bin/env bash
# 23_run_rexach_container.sh
# Step 3 of the Rexach container test in WSL: runs rexach-run all --compare in the image
# rexach-cepip, on the rerun folder made by step 2.
#   - /data      the rerun folder (~/Rexach_rerun)
#   - /ref       the original External_Dataset_2 folder, read-only
#   - /ref_supp  the original Manuscript/Supp_Tables folder, read-only (for Table S21)
# The steps print their progress as they run. At the end, the run summary, the comparison log, and
# the guardrail logs of R14 to R17 are copied to a results folder in Downloads, for upload.
#
# Usage, in WSL:
#   cd /mnt/c/Users/Wong/Downloads/rexach_container
#   bash 23_run_rexach_container.sh
# To resume after a failed step (the summary prints the step name):
#   FROM=R15_Rexach_diagnosis_contrasts bash 23_run_rexach_container.sh
# Settings: WORKERS (cepip workers in R14, default 6), FIGURE (Table S21 figure number, default 19).
# Log: 23_run_rexach_container_log.txt in this folder.

set -uo pipefail
SRC="$(cd "$(dirname "$0")" && pwd)"
RERUN="${RERUN:-$HOME/Rexach_rerun}"
ORIG="${ORIG:-/mnt/c/Users/Wong/Desktop/SEA-AD/External_Dataset_2}"
SUPP="${SUPP:-/mnt/c/Users/Wong/Desktop/SEA-AD/Manuscript/Supp_Tables}"
OUTC="${OUTC:-/mnt/c/Users/Wong/Downloads/rexach_results}"
WORKERS="${WORKERS:-6}"
FIGURE="${FIGURE:-19}"
IMAGE=rexach-cepip
LOG="$SRC/23_run_rexach_container_log.txt"

: > "$LOG"
exec > >(tee -a "$LOG") 2>&1
echo "== 23_run_rexach_container.sh started $(date '+%Y-%m-%d %H:%M:%S')"
if command -v podman >/dev/null 2>&1; then ENG=podman; else ENG=docker; fi
$ENG image inspect "$IMAGE" >/dev/null 2>&1 || { echo "STOPPED: image $IMAGE not found; run 21_build_rexach_image.sh"; exit 1; }
for d in "$RERUN" "$ORIG" "$SUPP"; do [ -d "$d" ] || { echo "STOPPED: $d not found"; exit 1; }; done
echo "Engine $ENG | rerun $RERUN | reference $ORIG | Supp_Tables $SUPP | workers $WORKERS | figure $FIGURE"
free -g | sed -n 1,2p

# WSL folders do not support SELinux labels, so the volumes are mounted without :Z.
ARGS=(all --compare)
[ -n "${FROM:-}" ] && ARGS+=(--from "$FROM")
t0=$(date +%s)
$ENG run --rm \
    -e REXACH_WORKERS="$WORKERS" -e REXACH_FIGURE="$FIGURE" \
    -v "$RERUN":/data \
    -v "$ORIG":/ref:ro \
    -v "$SUPP":/ref_supp:ro \
    "$IMAGE" rexach-run "${ARGS[@]}"
rc=$?
echo "== rexach-run exit code $rc | $(( ($(date +%s) - t0 + 30) / 60 )) min"

# ---- Copy the logs for upload ---------------------------------------------------------------------
mkdir -p "$OUTC"
for f in rexach-run_logs/rexach-run_summary.txt compare/compare_rexach_all_log.txt compare/compare_rexach_all_detail.csv \
         R14_readings/R14_Rexach_Guardrail_Log.txt R15_readings/R15_Rexach_Guardrail_Log.txt \
         R17_posthoc/R17_Rexach_Guardrail_Log.txt R16_Figure_Rexach/R16_Rexach_Guardrail_Log.txt \
         Supp_Tables/build_SuppTable_S21_log.txt; do
    [ -f "$RERUN/$f" ] && cp "$RERUN/$f" "$OUTC/" && echo "  copied $f"
done
cp "$LOG" "$OUTC/" 2>/dev/null
echo "Results for upload: $OUTC"
echo "== 23_run_rexach_container.sh finished $(date '+%Y-%m-%d %H:%M:%S')"
exit "$rc"
