#!/usr/bin/env bash
# Script 1 (v2): Cell Ranger path check, GRCh38-2020-A reference, pilot FASTQ download.
# Uses curl over https with resume. Output is saved to 01_setup_log.txt.
# Run in the background:  nohup bash 01_setup.sh > /dev/null 2>&1 &
set -uo pipefail
BASE="/home/wong/Desktop/Ext_Validation_1"
CR="$BASE/CellRanger/cellranger-10.1.0"
REFDIR="$BASE/Reference"
REFNAME="refdata-gex-GRCh38-2020-A"
REFURL="https://cf.10xgenomics.com/supp/cell-exp/${REFNAME}.tar.gz"
RUN="SRR27882204"
LIB="210526_FTLD_228_snRNA-C8"
FQ="$BASE/fastq/$LIB"
ENA="https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR278/004/${RUN}"
mkdir -p "$REFDIR" "$FQ"
LOG="$BASE/01_setup_log.txt"

fetch() {  # fetch URL OUTFILE: resume until curl reports success, up to 20 tries
  local url="$1" out="$2" i
  for i in $(seq 1 20); do
    if curl -sS -L -C - --retry 5 --retry-delay 10 -o "$out" "$url"; then
      echo "  done: $out ($(du -h "$out" | cut -f1))"; return 0
    fi
    echo "  curl stopped (try $i); resuming in 30 s"; sleep 30
  done
  echo "  FAILED: $url"; return 1
}

{
  echo "Start: $(date)"
  echo "== Cell Ranger =="
  export PATH="$CR:$PATH"
  grep -q "cellranger-10.1.0" ~/.bashrc || echo "export PATH=\"$CR:\$PATH\"" >> ~/.bashrc
  cellranger --version
  echo
  echo "== Machine =="
  nproc; free -g; df -h "$BASE"
  echo
  echo "== Reference $REFNAME =="
  cd "$REFDIR"
  if [ ! -s "$REFNAME/genes/genes.gtf" ] && [ ! -s "$REFNAME/genes/genes.gtf.gz" ]; then
    rm -rf "$REFNAME"
    echo "  downloading: $(date)"
    fetch "$REFURL" "${REFNAME}.tar.gz" || exit 1
    echo "  unpacking: $(date)"
    tar -xzf "${REFNAME}.tar.gz" || { echo "  tar failed: archive incomplete; rerun to resume"; exit 1; }
  else
    echo "  already present"
  fi
  ls "$REFNAME"
  echo "  genes in reference: $(zcat -f $REFNAME/genes/genes.gtf* | awk '$3=="gene"' | wc -l) (expected 36601)"
  echo
  echo "== Pilot FASTQ: $LIB ($RUN) =="
  cd "$FQ"
  for R in 1 2; do
    OUT="${LIB}_S1_L001_R${R}_001.fastq.gz"
    echo "  R${R} start: $(date)"
    fetch "${ENA}/${RUN}_${R}.fastq.gz" "$OUT" || exit 1
    gzip -t "$OUT" && echo "  R${R} gzip integrity OK" || echo "  R${R} gzip integrity FAILED"
  done
  ls -l
  echo "  R1 read length: $(zcat ${LIB}_S1_L001_R1_001.fastq.gz | sed -n 2p | tr -d '\n' | wc -c) (expected 28)"
  echo "  R2 read length: $(zcat ${LIB}_S1_L001_R2_001.fastq.gz | sed -n 2p | tr -d '\n' | wc -c) (expected 90)"
  df -h "$BASE"
  echo "End: $(date)"
} 2>&1 | tee -a "$LOG"
