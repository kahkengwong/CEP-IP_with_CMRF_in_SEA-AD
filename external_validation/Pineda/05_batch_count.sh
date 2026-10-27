#!/usr/bin/env bash
# Script 5: cellranger count on each library downloaded by scripts 7 and 8, one library at a time.
# Runs alongside the download scripts. Waits when no downloaded library is ready.
# For each library: count, check outputs, move the 4 kept files to counts/, then delete
# the FASTQ files and the Cell Ranger working folder to free disk for the downloads.
# Waits while a download script (4, 7 or 8) is running; stops when none is running and nothing is left.
# Safe to stop and run again: libraries with a counted marker are skipped.
# Output is saved to 05_batch_count_log.txt; Cell Ranger's own output goes to cellranger_logs/.
# Run in the background:  nohup bash 05_batch_count.sh > /dev/null 2>&1 &
set -uo pipefail
BASE="/home/wong/Desktop/Ext_Validation_1"
export PATH="$BASE/CellRanger/cellranger-10.1.0:$PATH"
REF="$BASE/Reference/refdata-gex-GRCh38-2020-A"
MAN="$BASE/Pineda_BA9_210526_manifest.tsv"
PILOT="210526_FTLD_228_snRNA-C8"
CORES=31
MEM=56
MARK="$BASE/markers"
OUT="$BASE/cellranger_out"
KEEP="$BASE/counts"
LOGDIR="$BASE/cellranger_logs"
LOG="$BASE/05_batch_count_log.txt"
MET="$BASE/05_cellranger_metrics.csv"
mkdir -p "$MARK" "$OUT" "$KEEP" "$LOGDIR"

keep_outputs () {   # argument: LIB. Moves the 4 kept files to counts/ and adds the metrics row.
  local LIB="$1" O="$OUT/$1/outs" F
  for F in raw_feature_bc_matrix.h5 filtered_feature_bc_matrix.h5 metrics_summary.csv web_summary.html; do
    [ -s "$O/$F" ] || { echo "  missing output: $F"; return 1; }
  done
  [ -s "$MET" ] || printf 'library,%s\n' "$(head -n 1 "$O/metrics_summary.csv")" > "$MET"
  printf '%s,%s\n' "$LIB" "$(sed -n 2p "$O/metrics_summary.csv")" >> "$MET"
  for F in raw_feature_bc_matrix.h5 filtered_feature_bc_matrix.h5 metrics_summary.csv web_summary.html; do
    mv "$O/$F" "$KEEP/${LIB}_$F"
  done
}

n_marks () { ls "$MARK" | grep -c "\.$1\$"; }

{
  echo "Start: $(date)"
  cellranger --version
  N_LIB=$(( $(wc -l < "$MAN") - 1 ))
  echo "Libraries in manifest: $N_LIB"

  # The pilot was counted by script 2 (24 cores); keep its outputs instead of counting it again.
  if [ ! -e "$MARK/$PILOT.counted" ] && [ -s "$OUT/$PILOT/outs/raw_feature_bc_matrix.h5" ]; then
    echo "== $PILOT: pilot outputs from script 2"
    if keep_outputs "$PILOT"; then
      rm -rf "$OUT/$PILOT" "$BASE/fastq/$PILOT"
      touch "$MARK/$PILOT.counted"
      echo "DONE $PILOT (pilot); FASTQ and working folder deleted"
    else
      echo "Pilot outputs incomplete; the pilot is left as it is"
    fi
  fi

  while true; do
    NEXT=""
    while IFS=$'\t' read -r -u 3 LIB REST; do
      if [ -e "$MARK/$LIB.downloaded" ] && [ ! -e "$MARK/$LIB.counted" ] && [ ! -e "$MARK/$LIB.count_failed" ]; then
        NEXT="$LIB"; break
      fi
    done 3< <(tail -n +2 "$MAN")

    if [ -z "$NEXT" ]; then
      N_END=$(( $(n_marks counted) + $(n_marks count_failed) ))
      if [ "$N_END" -ge "$N_LIB" ]; then echo "All $N_LIB libraries finished: $(date)"; break; fi
      if pgrep -f "bash .*0[478]_(batch_download|ncbi_fallback|batch_ncbi)\.sh" > /dev/null; then sleep 300; continue; fi
      echo "No download script is running and no downloaded library is waiting; stopping: $(date)"
      break
    fi

    echo "== $NEXT: count start: $(date)"
    rm -rf "$OUT/$NEXT"
    cd "$OUT"
    CLOG="$LOGDIR/${NEXT}_cellranger.txt"
    /usr/bin/time -v cellranger count --id="$NEXT" --sample="$NEXT" --fastqs="$BASE/fastq/$NEXT" \
      --transcriptome="$REF" --chemistry=auto --include-introns=true --create-bam=false \
      --nosecondary --localcores=$CORES --localmem=$MEM > "$CLOG" 2>&1
    RC=$?
    cd "$BASE"
    WALL=$(grep "Elapsed (wall clock)" "$CLOG" | awk -F': ' '{print $2}')
    RSS=$(grep "Maximum resident set size" "$CLOG" | awk -F': ' '{printf "%.1f GB", $2/1048576}')
    if [ "$RC" -eq 0 ] && grep -q "Pipestance completed successfully" "$CLOG" && keep_outputs "$NEXT"; then
      rm -rf "$OUT/$NEXT" "$BASE/fastq/$NEXT"
      touch "$MARK/$NEXT.counted"
      echo "DONE $NEXT: wall clock $WALL; largest process $RSS; FASTQ and working folder deleted"
    else
      touch "$MARK/$NEXT.count_failed"
      echo "FAILED $NEXT: exit $RC; FASTQ and working folder kept; see $CLOG"
    fi
    df -h "$BASE" | tail -n 1
  done

  echo
  echo "Counted: $(n_marks counted); failed: $(n_marks count_failed)"
  echo "End: $(date)"
} 2>&1 | tee -a "$LOG"
