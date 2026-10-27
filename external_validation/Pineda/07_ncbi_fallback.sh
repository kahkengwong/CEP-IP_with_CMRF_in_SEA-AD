#!/usr/bin/env bash
# Script 7: get one library from NCBI SRA when its ENA FASTQ fails the gzip check.
# Arguments: library name and the SRA Lite https link from the NCBI "Data access" tab, e.g.
#   bash 07_ncbi_fallback.sh 210526_PN_301_snRNA-D7 https://sra-downloadb.be-md.ncbi.nlm.nih.gov/.../SRR27882193.lite.1
# Steps, each with a progress line every 30 s:
#   1) aria2c download of the SRA Lite file
#   2) fasterq-dump to R1 and R2 FASTQ
#   3) pigz to the file names Cell Ranger expects, then gzip -t on both
#   4) delete the SRA Lite and plain FASTQ files, write the download marker
# 05_batch_count.sh then counts the library as usual.
# Output is saved to 07_ncbi_fallback_log.txt.
# Run in the background:  nohup bash 07_ncbi_fallback.sh LIB URL > /dev/null 2>&1 &
# Watch continuously:     tail -f 07_ncbi_fallback_log.txt   (Ctrl+C stops the view only)
set -uo pipefail
BASE="/home/wong/Desktop/Ext_Validation_1"
SRA_BIN="/home/wong/miniconda3/envs/sra/bin"
THREADS=8
MAN="$BASE/Pineda_BA9_210526_manifest.tsv"
MARK="$BASE/markers"
LOG="$BASE/07_ncbi_fallback_log.txt"
LIB="${1:?library name missing}"
URL="${2:?SRA Lite link missing}"

pct_done () {   # argument: .aria2 control file. Prints done GB, total GB and %.
  python3 -I - "$1" <<'PY'
import struct, sys
try:
    b = open(sys.argv[1], 'rb').read()
    e = '>' if struct.unpack('>H', b[0:2])[0] == 1 else '<'
    p = 6
    p += 4 + struct.unpack(e + 'I', b[p:p+4])[0]
    plen = struct.unpack(e + 'I', b[p:p+4])[0]; p += 4
    total = struct.unpack(e + 'Q', b[p:p+8])[0]; p += 16
    bfl = struct.unpack(e + 'I', b[p:p+4])[0]; p += 4
    done = sum(bin(x).count('1') for x in b[p:p+bfl]) * plen; p += bfl
    n = struct.unpack(e + 'I', b[p:p+4])[0]; p += 4
    for _ in range(n):
        p += 8
        l = struct.unpack(e + 'I', b[p:p+4])[0]; p += 4
        done += sum(bin(x).count('1') for x in b[p:p+l]) * 16384; p += l
    done = min(done, total)
    print("%.1f/%.1f GB (%.0f%%)" % (done / 1e9, total / 1e9, 100.0 * done / total))
except Exception:
    print("progress file not saved yet")
PY
}

watch_pid () {   # arguments: PID LABEL FILES... Prints GB written and speed every 30 s while PID runs.
  local P="$1" L="$2" I S0 S1; shift 2
  S0=$(du -cb "$@" 2>/dev/null | tail -n 1 | cut -f1); S0=${S0:-0}
  while kill -0 "$P" 2>/dev/null; do
    for I in $(seq 1 30); do kill -0 "$P" 2>/dev/null || break; sleep 1; done
    kill -0 "$P" 2>/dev/null || break
    S1=$(du -cb "$@" 2>/dev/null | tail -n 1 | cut -f1); S1=${S1:-0}
    echo "  $(date +%H:%M:%S) $L: $(awk -v b="$S1" 'BEGIN {printf "%.1f GB", b / 1e9}') written, $(awk -v a="$S0" -v b="$S1" 'BEGIN {printf "%.1f MiB/s", (b - a) / 30 / 1048576}')"
    S0=$S1
  done
}

{
  echo "Start: $(date)"
  RUN=$(awk -F'\t' -v l="$LIB" '$1 == l {print $2}' "$MAN")
  [ -n "$RUN" ] || { echo "Library not in manifest: $LIB"; exit 1; }
  case "$URL" in *"$RUN"*lite*) ;; *) echo "Link does not contain $RUN and lite: $URL"; exit 1;; esac
  if [ -e "$MARK/$LIB.downloaded" ] || [ -e "$MARK/$LIB.counted" ]; then echo "$LIB already has a marker; nothing to do"; exit 0; fi
  echo "Library $LIB, run $RUN"
  "$SRA_BIN/fasterq-dump" --version | head -n 2
  WORK="$BASE/ncbi_work/$LIB"; FQ="$BASE/fastq/$LIB"
  SRA="$BASE/ncbi_work/$RUN.sralite"
  mkdir -p "$WORK/tmp" "$FQ"; cd "$WORK"

  # 1) Download
  echo "== 1) download SRA Lite: $(date)"
  OK=0
  for A in $(seq 1 30); do
    aria2c -c -x16 -s16 -k1M --max-tries=0 --retry-wait=30 --timeout=60 \
      --show-console-readout=false --summary-interval=0 --auto-save-interval=30 --console-log-level=warn \
      --download-result=hide -d "$BASE/ncbi_work" -o "$RUN.sralite" "$URL" < /dev/null &
    APID=$!
    W0=$(awk '/^wchar/ {print $2}' /proc/$APID/io 2>/dev/null || echo 0)
    while kill -0 "$APID" 2>/dev/null; do
      for I in $(seq 1 30); do kill -0 "$APID" 2>/dev/null || break; sleep 1; done
      kill -0 "$APID" 2>/dev/null || break
      W1=$(awk '/^wchar/ {print $2}' /proc/$APID/io 2>/dev/null || echo "$W0")
      echo "  $(date +%H:%M:%S) download $RUN.sralite: $(pct_done "$SRA.aria2"), $(awk -v a="$W0" -v b="$W1" 'BEGIN {printf "%.1f MiB/s", (b - a) / 30 / 1048576}')"
      W0=$W1
    done
    wait "$APID" && { OK=1; break; }
    echo "  aria2c exit nonzero (attempt $A); resuming in 30 s"; sleep 30
  done
  [ "$OK" -eq 1 ] && [ ! -e "$SRA.aria2" ] || { echo "Download FAILED"; exit 1; }
  echo "  downloaded $(du -b "$SRA" | cut -f1) bytes: $(date)"

  # 2) SRA Lite to plain FASTQ
  echo "== 2) fasterq-dump: $(date)"
  rm -f "${RUN}"_*.fastq
  "$SRA_BIN/fasterq-dump" "$SRA" --split-files --include-technical --outfile "$RUN.fastq" \
    -O "$WORK" -t "$WORK/tmp" -e "$THREADS" > fasterq_dump_messages.txt 2>&1 &
  DPID=$!
  watch_pid "$DPID" "fasterq-dump" "$WORK"
  wait "$DPID"; RC=$?
  cat fasterq_dump_messages.txt
  ls -l "${RUN}"_*.fastq 2>/dev/null
  N_OUT=$(ls "${RUN}"_*.fastq 2>/dev/null | wc -l)
  if [ "$RC" -ne 0 ] || [ "$N_OUT" -ne 2 ] || [ ! -s "${RUN}_1.fastq" ] || [ ! -s "${RUN}_2.fastq" ]; then
    echo "fasterq-dump FAILED or did not give exactly 2 files (exit $RC, files $N_OUT); work folder kept: $WORK"; exit 1
  fi
  echo "  R1 first read:"; head -n 4 "${RUN}_1.fastq"
  echo "  R2 first read:"; head -n 4 "${RUN}_2.fastq"

  # 3) Compress and check
  echo "== 3) pigz and gzip -t: $(date)"
  for R in 1 2; do
    OUT="$FQ/${LIB}_S1_L001_R${R}_001.fastq.gz"
    "$SRA_BIN/pigz" -p "$THREADS" -c "${RUN}_${R}.fastq" > "$OUT" &
    ZPID=$!
    watch_pid "$ZPID" "pigz R$R" "$OUT"
    wait "$ZPID" || { echo "pigz FAILED on R$R"; exit 1; }
    echo "  R$R: gzip -t check: $(date +%H:%M:%S)"
    gzip -t "$OUT" || { echo "gzip -t FAILED on R$R"; exit 1; }
    echo "  R$R: $(du -b "$OUT" | cut -f1) bytes, gzip OK"
  done

  # 4) Clean up and hand over to script 5
  cd "$BASE"
  rm -rf "$WORK" "$SRA"
  touch "$MARK/$LIB.downloaded"
  echo "DONE $LIB from NCBI SRA Lite: $(date). Marker written; 05_batch_count.sh will count it."
  df -h "$BASE" | tail -n 1
} 2>&1 | tee -a "$LOG"
