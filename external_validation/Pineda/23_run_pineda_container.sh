#!/usr/bin/env bash
# 23_run_pineda_container.sh
# Step 3 of the Pineda container test in WSL: runs pineda-run all --compare in the image
# pineda-cepip, on the rerun folder made by step 2.
#   - /data      the rerun folder (~/Pineda_rerun)
#   - /ref       the original External_Dataset_1 folder, read-only
#   - /ref_supp  the original Manuscript/Supp_Tables folder, read-only (for Table S20)
# The steps print their progress as they run. At the end, the run summary, the comparison logs, and
# the guardrail logs are copied to a results folder in Downloads, for upload.
#
# Usage, in WSL:
#   cd /mnt/c/Users/Wong/Downloads/pineda_container
#   bash 23_run_pineda_container.sh
# To resume after a failed step (the summary prints the step name):
#   FROM=R14_Pineda_core_cepip bash 23_run_pineda_container.sh
# Settings: WORKERS (libraries at the same time in R02 and cepip workers in R14, default 4, because WSL has
# half the memory of Windows by default; R02 ran with 6 on a 48 GB laptop in Windows), FIGURE (Table S20 figure number, default 19).
# Log: 23_run_pineda_container_log.txt in this folder.

set -uo pipefail
SRC="$(cd "$(dirname "$0")" && pwd)"
RERUN="${RERUN:-$HOME/Pineda_rerun}"
ORIG="${ORIG:-/mnt/c/Users/Wong/Desktop/SEA-AD/External_Dataset_1}"
SUPP="${SUPP:-/mnt/c/Users/Wong/Desktop/SEA-AD/Manuscript/Supp_Tables}"
OUTC="${OUTC:-/mnt/c/Users/Wong/Downloads/pineda_results}"
WORKERS="${WORKERS:-4}"
FIGURE="${FIGURE:-19}"
IMAGE=pineda-cepip
LOG="$SRC/23_run_pineda_container_log.txt"

: > "$LOG"
exec > >(tee -a "$LOG") 2>&1
echo "== 23_run_pineda_container.sh started $(date '+%Y-%m-%d %H:%M:%S')"
if command -v podman >/dev/null 2>&1; then ENG=podman; else ENG=docker; fi
$ENG image inspect "$IMAGE" >/dev/null 2>&1 || { echo "STOPPED: image $IMAGE not found; run 21_build_pineda_image.sh"; exit 1; }
for d in "$RERUN" "$ORIG" "$SUPP"; do [ -d "$d" ] || { echo "STOPPED: $d not found"; exit 1; }; done
echo "Engine $ENG | rerun $RERUN | reference $ORIG | Supp_Tables $SUPP | workers $WORKERS | figure $FIGURE"
free -g | sed -n 1,2p

# WSL folders do not support SELinux labels, so the volumes are mounted without :Z.
ARGS=(all --compare)
[ -n "${FROM:-}" ] && ARGS+=(--from "$FROM")
# Out-of-memory lines already in the kernel log before this run are not reported again.
OOM0=$(dmesg 2>/dev/null | grep -c -i -E "out of memory|killed process")
t0=$(date +%s)
$ENG run --rm \
    -e PINEDA_WORKERS="$WORKERS" -e PINEDA_FIGURE="$FIGURE" \
    -v "$RERUN":/data \
    -v "$ORIG":/ref:ro \
    -v "$SUPP":/ref_supp:ro \
    "$IMAGE" pineda-run "${ARGS[@]}"
rc=$?
echo "== pineda-run exit code $rc | $(( ($(date +%s) - t0 + 30) / 60 )) min"

# ---- Copy the logs for upload ---------------------------------------------------------------------
mkdir -p "$OUTC"
for f in pineda-run_logs/pineda-run_summary.txt compare/compare_pineda_all_log.txt \
         compare/compare_pineda_all_detail.csv compare/compare_pineda_all_R02_libraries.csv \
         R02_cell_calling/R02_cell_calling_log.txt R13_count_checks/R13_Pineda_log.txt \
         R14_core/R14_Pineda_Guardrail_Log.txt R15_contrasts/R15_Pineda_Guardrail_Log.txt \
         R17_posthoc/R17_Pineda_Guardrail_Log.txt R16_Figure_Pineda/R16_Pineda_Guardrail_Log.txt \
         Supp_Tables/build_SuppTable_S20_log.txt; do
    [ -f "$RERUN/$f" ] && cp "$RERUN/$f" "$OUTC/" && echo "  copied $f"
done
# The console log of every step that ran (a failed step's error is there), and the R02 progress file.
for f in "$RERUN"/pineda-run_logs/*.txt "$RERUN/R02_cell_calling/R02_progress.txt"; do
    [ -f "$f" ] || continue
    case "$(basename "$f")" in pineda-run_summary.txt) continue ;; esac
    cp "$f" "$OUTC/step_$(basename "$f")" && echo "  copied ${f#$RERUN/}"
done
# A worker killed for lack of memory leaves no R error; the kernel log records it.
OOM=$(dmesg 2>/dev/null | grep -i -E "out of memory|killed process" | tail -n +$((OOM0 + 1)))
if [ -n "$OOM" ]; then
    echo "  The kernel log reports processes killed for lack of memory during this run:"
    echo "$OOM" | tail -n 5 | tee "$OUTC/oom_messages.txt"
fi
cp "$LOG" "$OUTC/" 2>/dev/null
echo "Results for upload: $OUTC"
echo "== 23_run_pineda_container.sh finished $(date '+%Y-%m-%d %H:%M:%S')"
exit "$rc"
