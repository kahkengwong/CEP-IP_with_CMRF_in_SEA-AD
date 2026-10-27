#!/usr/bin/env bash
# Script 2: pilot cellranger count on one library. Output is saved to 02_pilot_count_log.txt.
set -euo pipefail
BASE="/home/wong/Desktop/Ext_Validation_1"
export PATH="$BASE/CellRanger/cellranger-10.1.0:$PATH"
REF="$BASE/Reference/refdata-gex-GRCh38-2020-A"
LIB="210526_FTLD_228_snRNA-C8"
FQ="$BASE/fastq/$LIB"
OUTDIR="$BASE/cellranger_out"
mkdir -p "$OUTDIR"
cd "$OUTDIR"
LOG="$BASE/02_pilot_count_log.txt"
{
  date
  cellranger --version
  /usr/bin/time -v cellranger count \
    --id="$LIB" \
    --sample="$LIB" \
    --fastqs="$FQ" \
    --transcriptome="$REF" \
    --chemistry=auto \
    --include-introns=true \
    --create-bam=false \
    --localcores=24 \
    --localmem=56
  echo
  echo "== Outputs =="
  ls -l "$OUTDIR/$LIB/outs"
  du -sh "$OUTDIR/$LIB"
  echo
  echo "== metrics_summary.csv =="
  cat "$OUTDIR/$LIB/outs/metrics_summary.csv"
  df -h "$BASE"
  date
} 2>&1 | tee "$LOG"
