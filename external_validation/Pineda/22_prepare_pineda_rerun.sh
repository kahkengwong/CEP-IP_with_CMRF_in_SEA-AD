#!/usr/bin/env bash
# 22_prepare_pineda_rerun.sh
# Step 2 of the Pineda container test in WSL: makes a new data folder for the rerun that holds only
# the inputs, copied from the original External_Dataset_1 folder:
#   - counts/<LIB>_raw_feature_bc_matrix.h5 (the 63 Cell Ranger raw count matrices, about 10 GB);
#   - 05_cellranger_metrics.csv and Pineda_BA9_210526_manifest.tsv;
#   - the S01 workbook (S01_sample_information.xlsx, Synapse syn51105515);
#   - Pineda2024_metadata/ (the authors' <LIB>_col_metadata.tsv files);
#   - Stats_15H.csv (SEA-AD Figure 15H source data);
#   - the three MapMyCells folders: R13_mapmycells_output/ (12 csv files, the nuclei of the authors'
#     QC), R13_mapmycells_output_allKept/ (26, all R02 kept nuclei), and
#     R13_mapmycells_output_allKept_repeat/ (26, the same file mapped again).
# MapMyCells is a web service, and the same file mapped twice gave other labels for 0.8% of nuclei,
# so the rerun uses these mapping results.
# The rerun folder is on the Linux disk (~/Pineda_rerun), which the container reads faster than
# /mnt/c. The script stops if the folder already holds files, so an earlier rerun is never mixed in.
# Copy progress is shown as it runs.
#
# Usage, in WSL:
#   cd /mnt/c/Users/Wong/Downloads/pineda_container
#   bash 22_prepare_pineda_rerun.sh
# Other folders: ORIG=/path/to/External_Dataset_1 RERUN=/path/to/new_folder bash 22_prepare_pineda_rerun.sh
# Log: 22_prepare_pineda_rerun_log.txt in this folder.

set -uo pipefail
SRC="$(cd "$(dirname "$0")" && pwd)"
ORIG="${ORIG:-/mnt/c/Users/Wong/Desktop/SEA-AD/External_Dataset_1}"
RERUN="${RERUN:-$HOME/Pineda_rerun}"
LOG="$SRC/22_prepare_pineda_rerun_log.txt"
MAPS=("R13_mapmycells_output|12" "R13_mapmycells_output_allKept|26" "R13_mapmycells_output_allKept_repeat|26")
TOP=(05_cellranger_metrics.csv Pineda_BA9_210526_manifest.tsv Stats_15H.csv)

: > "$LOG"
exec > >(tee -a "$LOG") 2>&1
echo "== 22_prepare_pineda_rerun.sh started $(date '+%Y-%m-%d %H:%M:%S')"
echo "Original folder: $ORIG"
echo "Rerun folder:    $RERUN"

[ -d "$ORIG" ] || { echo "STOPPED: $ORIG not found"; exit 1; }
if [ -d "$RERUN" ] && [ -n "$(ls -A "$RERUN" 2>/dev/null)" ]; then
    echo "STOPPED: $RERUN already holds files. Move it away or give another folder with RERUN=..."
    exit 1
fi

# ---- Inputs in the original folder ---------------------------------------------------------------
n_h5=$(find "$ORIG/counts" -maxdepth 1 -type f -name '*_raw_feature_bc_matrix.h5' -size +0 2>/dev/null | wc -l)
[ "$n_h5" -eq 63 ] || { echo "STOPPED: expected 63 non-empty raw .h5 files in $ORIG/counts, found $n_h5"; exit 1; }
for e in "${MAPS[@]}"; do
    d=${e%%|*}; k=${e##*|}
    n=$(find "$ORIG/$d" -maxdepth 1 -type f -name '*.csv' 2>/dev/null | wc -l)
    [ "$n" -eq "$k" ] || { echo "STOPPED: expected $k csv files in $ORIG/$d, found $n"; exit 1; }
done
for f in "${TOP[@]}"; do [ -f "$ORIG/$f" ] || { echo "STOPPED: $ORIG/$f not found"; exit 1; }; done
mapfile -t s01 < <(find "$ORIG" -maxdepth 1 -type f -iname 'S01*.xlsx' | sort)
[ "${#s01[@]}" -eq 1 ] || { echo "STOPPED: expected one S01*.xlsx in $ORIG, found ${#s01[@]}"; exit 1; }
n_meta=$(find "$ORIG/Pineda2024_metadata" -type f -name '*_col_metadata.tsv' 2>/dev/null | wc -l)
[ "$n_meta" -ge 63 ] || { echo "STOPPED: expected the 63 authors' <LIB>_col_metadata.tsv files in $ORIG/Pineda2024_metadata, found $n_meta"; exit 1; }

need_kb=$( { du -sk "$ORIG"/counts/*_raw_feature_bc_matrix.h5 "$ORIG/Pineda2024_metadata" "${s01[0]}";
             for e in "${MAPS[@]}"; do du -sk "$ORIG/${e%%|*}"; done;
             for f in "${TOP[@]}"; do du -sk "$ORIG/$f"; done; } | awk '{s += $1} END {print s}')
mkdir -p "$RERUN"
free_kb=$(df --output=avail -k "$RERUN" | tail -1 | tr -dc '0-9')
awk -v n="$need_kb" -v f="$free_kb" 'BEGIN {printf "Inputs: %.1f GB | free on the Linux disk: %.1f GB\n", n / 1048576, f / 1048576}'
# The rerun writes about 15 GB more (R02 kept counts, MapMyCells input files, and outputs).
if [ "$free_kb" -lt $((need_kb + 25 * 1048576)) ]; then
    echo "STOPPED: less than the inputs plus 25 GB free"; exit 1
fi

# ---- Copy ------------------------------------------------------------------------------------------
copy() {   # copy <source> <destination folder>
    if command -v rsync >/dev/null 2>&1; then rsync -a --info=progress2 "$1" "$2"
    else cp -a "$1" "$2"; fi
}
echo
echo "== Copy"
mkdir -p "$RERUN/counts"
echo "  counts/: 63 raw .h5 files"
if command -v rsync >/dev/null 2>&1; then
    rsync -a --info=progress2 --include='*_raw_feature_bc_matrix.h5' --exclude='*' "$ORIG/counts/" "$RERUN/counts/"
else
    cp -a "$ORIG"/counts/*_raw_feature_bc_matrix.h5 "$RERUN/counts/"
fi
for f in "${TOP[@]}"; do copy "$ORIG/$f" "$RERUN/"; done
copy "${s01[0]}" "$RERUN/"
copy "$ORIG/Pineda2024_metadata" "$RERUN/"
for e in "${MAPS[@]}"; do
    d=${e%%|*}; mkdir -p "$RERUN/$d"
    for f in "$ORIG/$d"/*.csv; do copy "$f" "$RERUN/$d/"; done
done

echo
echo "== Rerun folder"
(cd "$RERUN" && find . -maxdepth 1 -type f -printf '%10s  %p\n' | sort -k2)
fail=0
chk() {   # chk <label> <found> <expected>
    if [ "$2" -eq "$3" ]; then echo "  [PASS] $1: $2"; else echo "  [FAIL] $1: $2 (expected $3)"; fail=1; fi
}
chk "raw .h5 files in counts/" "$(find "$RERUN/counts" -maxdepth 1 -name '*_raw_feature_bc_matrix.h5' -size +0 | wc -l)" 63
for e in "${MAPS[@]}"; do
    d=${e%%|*}; chk "csv files in $d/" "$(find "$RERUN/$d" -maxdepth 1 -name '*.csv' | wc -l)" "${e##*|}"
done
chk "authors' metadata files" "$(find "$RERUN/Pineda2024_metadata" -type f -name '*_col_metadata.tsv' | wc -l)" "$n_meta"
chk "S01 workbook" "$(find "$RERUN" -maxdepth 1 -iname 'S01*.xlsx' | wc -l)" 1
for f in "${TOP[@]}"; do chk "$f" "$([ -f "$RERUN/$f" ] && echo 1 || echo 0)" 1; done
# Sizes must equal the originals.
bad=$(cd "$ORIG" && for f in counts/*_raw_feature_bc_matrix.h5; do
          [ "$(stat -c %s "$f")" = "$(stat -c %s "$RERUN/$f" 2>/dev/null)" ] || echo "$f"; done)
if [ -z "$bad" ]; then echo "  [PASS] raw .h5 sizes equal the originals"; else echo "  [FAIL] sizes differ: $bad"; fail=1; fi
if [ "$fail" -eq 0 ]; then echo "PASS: the rerun folder holds all inputs"; else echo "STOPPED: the rerun folder is incomplete"; fi
echo "== 22_prepare_pineda_rerun.sh finished $(date '+%Y-%m-%d %H:%M:%S')"
exit "$fail"
