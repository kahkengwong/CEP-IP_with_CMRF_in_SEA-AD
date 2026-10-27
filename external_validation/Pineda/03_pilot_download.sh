#!/usr/bin/env bash
# Script 3: pilot FASTQ download with aria2c, gzip check, and automatic re-download.
# Output is saved to 03_pilot_download_log.txt.
# Run in the background:  nohup bash 03_pilot_download.sh > /dev/null 2>&1 &
set -uo pipefail
BASE="/home/wong/Desktop/Ext_Validation_1"
RUN="SRR27882204"
LIB="210526_FTLD_228_snRNA-C8"
FQ="$BASE/fastq/$LIB"
ENA="https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR278/004/${RUN}"
LOG="$BASE/03_pilot_download_log.txt"
mkdir -p "$FQ"
cd "$FQ"
{
  echo "Start: $(date)"
  curl -sS -o /dev/null --max-time 60 -w 'Line check (Cloudflare 30 MB): %{speed_download} bytes/s\n' \
    "https://speed.cloudflare.com/__down?bytes=31457280"
  ALL_OK=1
  for R in 1 2; do
    OUT="${LIB}_S1_L001_R${R}_001.fastq.gz"
    OK=0
    for TRY in 1 2 3; do
      if [ -s "$OUT" ] && [ ! -e "$OUT.aria2" ]; then
        echo "R${R}: checking existing file (try $TRY): $(date)"
        if gzip -t "$OUT" 2>/dev/null; then echo "R${R}: gzip OK"; OK=1; break; fi
        echo "R${R}: gzip FAILED; deleting and downloading again"
        rm -f "$OUT" "$OUT.aria2"
      fi
      echo "R${R}: download (try $TRY): $(date)"
      aria2c -c -x16 -s16 -k1M --max-tries=0 --retry-wait=30 --timeout=60 \
        --summary-interval=300 --console-log-level=warn \
        -o "$OUT" "${ENA}/${RUN}_${R}.fastq.gz"
      echo "R${R}: download finished: $(date); size $(stat -c %s "$OUT" 2>/dev/null)"
      if gzip -t "$OUT" 2>/dev/null; then echo "R${R}: gzip OK"; OK=1; break; fi
      echo "R${R}: gzip FAILED after download; deleting"
      rm -f "$OUT" "$OUT.aria2"
    done
    [ "$OK" -eq 1 ] || { echo "R${R}: still failing after 3 tries"; ALL_OK=0; }
  done
  ls -l "$FQ"
  echo "Expected sizes: R1 7715317270, R2 7210339319"
  if [ "$ALL_OK" -eq 1 ]; then echo "PILOT FASTQ READY"; else echo "PILOT FASTQ NOT READY"; fi
  echo "End: $(date)"
} 2>&1 | tee -a "$LOG"
