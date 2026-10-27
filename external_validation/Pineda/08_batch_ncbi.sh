#!/usr/bin/env bash
# Script 8: get the remaining BA9 libraries from NCBI SRA Lite. Replaces script 4 (ENA FASTQ files failed gzip -t).
# Order: PN, then FTLD, then ALS. Skips the pilot and any library that already has a marker.
# For each library, with a progress line every 30 s:
#   1) Look up the SRA Lite link and size (NCBI locate service). Download with aria2c. Check the size.
#   2) fastq-dump streams the read pairs (no temporary files). mawk checks that both reads of each pair
#      have the same name and writes R1 and R2. pigz compresses both at once. No plain FASTQ on disk.
#      (fasterq-dump is not used: its streamed output does not keep the two reads of a pair together,
#      and it writes very large temporary files.)
#   3) Check that pairs written = spots read by fastq-dump. Run gzip -t on R1 and R2.
# Stops if 2 libraries in a row fail, so a systematic fault cannot fill the disk overnight.
#   4) Delete the SRA Lite file and write the download marker. 05_batch_count.sh then counts the library.
# When free disk is below the floor, the script waits until script 5 frees space.
# Safe to stop and run again: libraries with a marker are skipped; an unfinished library starts again.
# Output is saved to 08_batch_ncbi_log.txt.
# Run in the background:  nohup bash 08_batch_ncbi.sh > /dev/null 2>&1 &
# Watch continuously:     tail -f 08_batch_ncbi_log.txt   (Ctrl+C stops the view only)
set -uo pipefail
BASE="/home/wong/Desktop/Ext_Validation_1"
SRA_BIN="/home/wong/miniconda3/envs/sra/bin"
MAN="$BASE/Pineda_BA9_210526_manifest.tsv"
PILOT="210526_FTLD_228_snRNA-C8"
PIGZ_THREADS=12      # per read file; R1 and R2 are compressed at the same time
FLOOR_GB=100         # free disk kept for Cell Ranger and libraries waiting to be counted
GZ_PER_LITE=1.75     # compressed R1 + R2 size per SRA Lite byte (PN_301: 24.4 GB from 14.0 GB)
MARK="$BASE/markers"
WROOT="$BASE/ncbi_work"
LOG="$BASE/08_batch_ncbi_log.txt"
mkdir -p "$MARK" "$WROOT" "$BASE/fastq"

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

lookup () {   # argument: RUN. Prints "LINK SIZE" of the SRA Lite file at NCBI, or nothing.
  curl -sS --max-time 60 "https://locate.ncbi.nlm.nih.gov/sdl/2/retrieve?acc=$1&accept-alternate-locations=yes" | \
  python3 -I -c '
import json, sys
run = sys.argv[1]
try:
    d = json.load(sys.stdin)
except Exception:
    sys.exit(0)
for r in d.get("result", []):
    for f in r.get("files", []):
        if f.get("name") == run + ".lite":
            for l in f.get("locations", []):
                if l.get("service") == "ncbi" and "link" in l:
                    print(l["link"], f["size"]); sys.exit(0)
' "$1"
}

download () {   # arguments: LINK SIZE FILE. aria2c with a progress line every 30 s, then a size check.
  local URL="$1" SIZE="$2" F="$3" TRY A I W0 W1 APID OK
  for TRY in 1 2 3; do
    if [ -s "$F" ] && [ ! -e "$F.aria2" ]; then
      [ "$(stat -c %s "$F")" = "$SIZE" ] && { echo "  size OK ($SIZE bytes)"; return 0; }
      echo "  size $(stat -c %s "$F") (expected $SIZE); deleting (try $TRY)"; rm -f "$F" "$F.aria2"
    fi
    echo "  download (try $TRY): $(date +%H:%M:%S)"
    OK=0
    for A in $(seq 1 30); do
      aria2c -c -x16 -s16 -k1M --max-tries=0 --retry-wait=30 --timeout=60 \
        --show-console-readout=false --summary-interval=0 --auto-save-interval=30 --console-log-level=warn \
        --download-result=hide -d "$(dirname "$F")" -o "$(basename "$F")" "$URL" < /dev/null &
      APID=$!
      W0=$(awk '/^wchar/ {print $2}' /proc/$APID/io 2>/dev/null || echo 0)
      while kill -0 "$APID" 2>/dev/null; do
        for I in $(seq 1 30); do kill -0 "$APID" 2>/dev/null || break; sleep 1; done
        kill -0 "$APID" 2>/dev/null || break
        W1=$(awk '/^wchar/ {print $2}' /proc/$APID/io 2>/dev/null || echo "$W0")
        echo "  $(date +%H:%M:%S) download $(basename "$F"): $(pct_done "$F.aria2"), $(awk -v a="$W0" -v b="$W1" 'BEGIN {printf "%.1f MiB/s", (b - a) / 30 / 1048576}')"
        W0=$W1
      done
      wait "$APID" && { OK=1; break; }
      echo "  aria2c exit nonzero (attempt $A); resuming in 30 s"; sleep 30
    done
    if [ "$OK" -eq 1 ] && [ ! -e "$F.aria2" ] && [ "$(stat -c %s "$F" 2>/dev/null)" = "$SIZE" ]; then
      echo "  size OK ($SIZE bytes): $(date +%H:%M:%S)"; return 0
    fi
    echo "  size $(stat -c %s "$F" 2>/dev/null) (expected $SIZE) after download; deleting"; rm -f "$F" "$F.aria2"
  done
  return 1
}

convert () {   # arguments: SRA_FILE LIB WORK SIZE. Streams to R1/R2 .fastq.gz with pair-name check.
  local SRA="$1" LIB="$2" WORK="$3" SIZE="$4" FQ="$BASE/fastq/$2" CPID Z1 Z2 RC RZ1 RZ2 I S0 S1 EXP SPOTS PAIRS
  local R1="$FQ/${LIB}_S1_L001_R1_001.fastq.gz" R2="$FQ/${LIB}_S1_L001_R2_001.fastq.gz"
  mkdir -p "$FQ"
  find "$FQ" -maxdepth 1 -type f -delete          # removes any old ENA or partial files of this library
  rm -f "$WORK/r1.fifo" "$WORK/r2.fifo"; mkfifo "$WORK/r1.fifo" "$WORK/r2.fifo"
  "$SRA_BIN/pigz" -p "$PIGZ_THREADS" -c < "$WORK/r1.fifo" > "$R1" & Z1=$!
  "$SRA_BIN/pigz" -p "$PIGZ_THREADS" -c < "$WORK/r2.fifo" > "$R2" & Z2=$!
  ( "$SRA_BIN/fastq-dump" --split-spot -Z "$SRA" 2> "$WORK/dump_messages.txt" | \
    "$SRA_BIN/mawk" -v r1="$WORK/r1.fifo" -v r2="$WORK/r2.fifo" '
      NR % 8 == 1 { n = $1 }
      NR % 8 == 5 { if ($1 != n) { print "PAIR NAME MISMATCH at line " NR > "/dev/stderr"; exit 3 } }
      { if ((NR - 1) % 8 < 4) print > r1; else print > r2 }
      END { if (NR % 8 != 0) { print "LINE COUNT NOT A MULTIPLE OF 8: " NR > "/dev/stderr"; exit 4 }
            printf "pairs written: %d\n", NR / 8 > "/dev/stderr" }' 2> "$WORK/split_messages.txt" ) &
  CPID=$!
  EXP=$(awk -v s="$SIZE" -v k="$GZ_PER_LITE" 'BEGIN {print s * k}')
  S0=0
  while kill -0 "$CPID" 2>/dev/null; do
    for I in $(seq 1 30); do kill -0 "$CPID" 2>/dev/null || break; sleep 1; done
    kill -0 "$CPID" 2>/dev/null || break
    S1=$(du -cb "$R1" "$R2" 2>/dev/null | tail -n 1 | cut -f1); S1=${S1:-0}
    echo "  $(date +%H:%M:%S) convert: $(awk -v b="$S1" -v e="$EXP" 'BEGIN {printf "%.1f GB compressed (about %.0f%% of expected)", b / 1e9, 100 * b / e}'), $(awk -v a="$S0" -v b="$S1" 'BEGIN {printf "%.1f MiB/s", (b - a) / 30 / 1048576}')"
    S0=$S1
  done
  wait "$CPID"; RC=$?
  # If the stream stopped before mawk opened a pipe, pigz is still waiting for it; open and close it to end pigz.
  exec 7<>"$WORK/r1.fifo" 8<>"$WORK/r2.fifo"; exec 7>&- 8>&-
  wait "$Z1"; RZ1=$?; wait "$Z2"; RZ2=$?
  cat "$WORK/dump_messages.txt" "$WORK/split_messages.txt"
  SPOTS=$(sed -n 's/^Read \([0-9,]*\) spots.*/\1/p' "$WORK/dump_messages.txt" | head -n 1 | tr -d ',')
  PAIRS=$(sed -n 's/^pairs written: \([0-9]*\).*/\1/p' "$WORK/split_messages.txt" | head -n 1)
  if [ "$RC" -ne 0 ] || [ "$RZ1" -ne 0 ] || [ "$RZ2" -ne 0 ] || [ -z "$SPOTS" ] || [ "$SPOTS" != "$PAIRS" ]; then
    echo "  convert FAILED (stream exit $RC, pigz exits $RZ1/$RZ2, spots read ${SPOTS:-none}, pairs written ${PAIRS:-none})"
    rm -f "$R1" "$R2"; return 1
  fi
  echo "  spots read = pairs written = $SPOTS"
  echo "  R1 first read:"; gzip -dc "$R1" 2>/dev/null | head -n 4
  echo "  R2 first read:"; gzip -dc "$R2" 2>/dev/null | head -n 4
  echo "  gzip -t on R1 and R2 (about 2 min): $(date +%H:%M:%S)"
  gzip -t "$R1" & Z1=$!; gzip -t "$R2" & Z2=$!
  wait "$Z1"; RZ1=$?; wait "$Z2"; RZ2=$?
  if [ "$RZ1" -ne 0 ] || [ "$RZ2" -ne 0 ]; then echo "  gzip -t FAILED (R1 $RZ1, R2 $RZ2)"; rm -f "$R1" "$R2"; return 1; fi
  echo "  gzip OK: R1 $(stat -c %s "$R1") bytes, R2 $(stat -c %s "$R2") bytes: $(date +%H:%M:%S)"
  return 0
}

{
  echo "Start: $(date)"
  [ -s "$MAN" ] || { echo "Manifest not found: $MAN"; exit 1; }
  for T in fastq-dump pigz mawk; do [ -x "$SRA_BIN/$T" ] || { echo "Missing $SRA_BIN/$T"; exit 1; }; done
  "$SRA_BIN/fastq-dump" --version | grep -m 1 .
  N_DONE=0; N_FAIL=0; N_RUN_FAIL=0
  while IFS=$'\t' read -r -u 3 LIB RUN REST; do
    [ "$LIB" = "$PILOT" ] && continue
    if [ -e "$MARK/$LIB.downloaded" ] || [ -e "$MARK/$LIB.counted" ] || [ -e "$MARK/$LIB.ncbi_failed" ]; then continue; fi
    echo "== $LIB ($RUN): $(date)"
    T0=$(date +%s)
    L=""; for A in 1 2 3; do L=$(lookup "$RUN"); [ -n "$L" ] && break; sleep 30; done
    if [ -z "$L" ]; then echo "FAILED $LIB: no SRA Lite link at NCBI"; touch "$MARK/$LIB.ncbi_failed"; N_FAIL=$((N_FAIL + 1)); continue; fi
    URL=${L% *}; SIZE=${L##* }
    echo "  SRA Lite: $SIZE bytes; $URL"
    NEED=$(awk -v s="$SIZE" -v k="$GZ_PER_LITE" -v f="$FLOOR_GB" 'BEGIN {printf "%.0f", s * (1 + k) + f * 1e9}')
    while [ "$(df --output=avail -B1 "$BASE" | tail -n 1)" -lt "$NEED" ]; do
      echo "  $(date +%H:%M:%S) waiting for disk: free $(df -h --output=avail "$BASE" | tail -n 1 | tr -d ' '), need $(( NEED / 1000000000 )) GB"
      sleep 300
    done
    WORK="$WROOT/$LIB"; mkdir -p "$WORK"
    SRA="$WORK/$RUN.sralite"
    if ! download "$URL" "$SIZE" "$SRA"; then
      echo "FAILED $LIB: download"; touch "$MARK/$LIB.ncbi_failed"; N_FAIL=$((N_FAIL + 1)); continue
    fi
    OK=0
    for A in 1 2; do
      echo "  convert (try $A): $(date +%H:%M:%S)"
      convert "$SRA" "$LIB" "$WORK" "$SIZE" && { OK=1; break; }
    done
    if [ "$OK" -eq 1 ]; then
      rm -rf "$WORK"
      touch "$MARK/$LIB.downloaded"; N_DONE=$((N_DONE + 1)); N_RUN_FAIL=0
      echo "DONE $LIB: $(( ($(date +%s) - T0) / 60 )) min; marker written"
    else
      echo "FAILED $LIB: convert; SRA Lite file kept in $WORK"; touch "$MARK/$LIB.ncbi_failed"; N_FAIL=$((N_FAIL + 1))
      N_RUN_FAIL=$((N_RUN_FAIL + 1))
      if [ "$N_RUN_FAIL" -ge 2 ]; then echo "STOPPED: 2 libraries in a row failed to convert: $(date)"; break; fi
    fi
    df -h "$BASE" | tail -n 1
  done 3< <(tail -n +2 "$MAN")
  echo
  echo "Done this session: $N_DONE; failed: $N_FAIL"
  echo "End: $(date)"
} 2>&1 | tee -a "$LOG"
