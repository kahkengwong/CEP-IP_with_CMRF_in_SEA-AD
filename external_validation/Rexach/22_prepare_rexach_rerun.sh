#!/usr/bin/env bash
# 22_prepare_rexach_rerun.sh
# Step 2 of the Rexach container test in WSL: makes a new data folder for the rerun that holds only
# the inputs, copied from the original External_Dataset_2 folder:
#   - the CELLxGENE h5ad (the only .h5ad in External_Dataset_2);
#   - R13_mapmycells_output/ (the six MapMyCells csv files);
#   - mmc1.xlsx (Table S1 of Rexach et al. 2024);
#   - Part12c_cell_metrics.csv (SEA-AD %MT per cell).
# The rerun folder is on the Linux disk (~/Rexach_rerun), which the container reads faster than
# /mnt/c. The script stops if the folder already holds files, so an earlier rerun is never mixed in.
# Copy progress is shown as it runs.
#
# Usage, in WSL:
#   cd /mnt/c/Users/Wong/Downloads/rexach_container
#   bash 22_prepare_rexach_rerun.sh
# Other folders: ORIG=/path/to/External_Dataset_2 RERUN=/path/to/new_folder bash 22_prepare_rexach_rerun.sh
# Log: 22_prepare_rexach_rerun_log.txt in this folder.

set -uo pipefail
SRC="$(cd "$(dirname "$0")" && pwd)"
ORIG="${ORIG:-/mnt/c/Users/Wong/Desktop/SEA-AD/External_Dataset_2}"
RERUN="${RERUN:-$HOME/Rexach_rerun}"
LOG="$SRC/22_prepare_rexach_rerun_log.txt"

: > "$LOG"
exec > >(tee -a "$LOG") 2>&1
echo "== 22_prepare_rexach_rerun.sh started $(date '+%Y-%m-%d %H:%M:%S')"
echo "Original folder: $ORIG"
echo "Rerun folder:    $RERUN"

[ -d "$ORIG" ] || { echo "STOPPED: $ORIG not found"; exit 1; }
if [ -d "$RERUN" ] && [ -n "$(ls -A "$RERUN" 2>/dev/null)" ]; then
    echo "STOPPED: $RERUN already holds files. Move it away or give another folder with RERUN=..."
    exit 1
fi

# ---- Inputs in the original folder ---------------------------------------------------------------
mapfile -t h5 < <(find "$ORIG" -maxdepth 1 -type f -name '*.h5ad' | sort)
[ "${#h5[@]}" -eq 1 ] || { echo "STOPPED: expected one .h5ad in $ORIG, found ${#h5[@]}"; exit 1; }
n_map=$(find "$ORIG/R13_mapmycells_output" -maxdepth 1 -type f -name '*.csv' 2>/dev/null | wc -l)
[ "$n_map" -eq 6 ] || { echo "STOPPED: expected 6 csv files in $ORIG/R13_mapmycells_output, found $n_map"; exit 1; }
for f in mmc1.xlsx Part12c_cell_metrics.csv; do
    [ -f "$ORIG/$f" ] || { echo "STOPPED: $ORIG/$f not found"; exit 1; }
done
need_kb=$(du -sk "${h5[0]}" "$ORIG/R13_mapmycells_output" "$ORIG/mmc1.xlsx" "$ORIG/Part12c_cell_metrics.csv" |
          awk '{s += $1} END {print s}')
mkdir -p "$RERUN"
free_kb=$(df --output=avail -k "$RERUN" | tail -1 | tr -dc '0-9')
awk -v n="$need_kb" -v f="$free_kb" 'BEGIN {printf "Inputs: %.1f GB | free on the Linux disk: %.1f GB\n", n / 1048576, f / 1048576}'
# The rerun writes about 5 GB more (MapMyCells input files, counts, and outputs).
if [ "$free_kb" -lt $((need_kb + 10 * 1048576)) ]; then
    echo "STOPPED: less than the inputs plus 10 GB free"; exit 1
fi

# ---- Copy ------------------------------------------------------------------------------------------
copy() {   # copy <source> <destination folder>
    if command -v rsync >/dev/null 2>&1; then rsync -a --info=progress2 "$1" "$2"
    else cp -av "$1" "$2"; fi
}
echo
echo "== Copy"
copy "${h5[0]}" "$RERUN/"
mkdir -p "$RERUN/R13_mapmycells_output"
for f in "$ORIG"/R13_mapmycells_output/*.csv; do copy "$f" "$RERUN/R13_mapmycells_output/"; done
copy "$ORIG/mmc1.xlsx" "$RERUN/"
copy "$ORIG/Part12c_cell_metrics.csv" "$RERUN/"

echo
echo "== Rerun folder"
(cd "$RERUN" && find . -type f -printf '%10s  %p\n' | sort -k2)
n_h5=$(find "$RERUN" -maxdepth 1 -name '*.h5ad' | wc -l)
n_csv=$(find "$RERUN/R13_mapmycells_output" -maxdepth 1 -name '*.csv' | wc -l)
if [ "$n_h5" -eq 1 ] && [ "$n_csv" -eq 6 ] && [ -f "$RERUN/mmc1.xlsx" ] && [ -f "$RERUN/Part12c_cell_metrics.csv" ]; then
    echo "PASS: one h5ad, 6 MapMyCells csv files, mmc1.xlsx, and Part12c_cell_metrics.csv"
else
    echo "STOPPED: the rerun folder is incomplete"; exit 1
fi
echo "== 22_prepare_rexach_rerun.sh finished $(date '+%Y-%m-%d %H:%M:%S')"
