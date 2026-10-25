#!/usr/bin/env bash
# 22b_download_rexach_inputs.sh
# Step 2 for a tester or reader without the original External_Dataset_2 folder (use it in place of
# 22_prepare_rexach_rerun.sh). Makes a new data folder for the rerun (~/Rexach_rerun) and downloads
# the inputs into it:
#   - the CELLxGENE h5ad, dataset version aa87f914-07e0-48aa-b915-8de906c95baf (5.0 GB);
#   - Part12c_cell_metrics.csv and the six MapMyCells csv files, from the Hugging Face dataset
#     kahkengwong/CEP-IP_with_CMRF_in_SEA-AD, folder external_validation/Rexach (107 MB).
# Table S1 of Rexach et al. 2024 (mmc1.xlsx) is the publisher's file and is not downloaded here.
# Download it from the paper's page (Cell 187:5753, Supplemental Information, Table S1) and give its
# path with MMC1=..., or place it in the rerun folder before running rexach-run.
# Every download shows its progress as it runs, can be resumed if it stops (run the script again),
# and is checked against the size of the original file.
#
# Usage, in WSL:
#   cd /path/to/rexach_container
#   MMC1=/mnt/c/Users/<name>/Downloads/mmc1.xlsx bash 22b_download_rexach_inputs.sh
# Other folder: RERUN=/path/to/new_folder bash 22b_download_rexach_inputs.sh
# Log: 22b_download_rexach_inputs_log.txt in this folder.

set -uo pipefail
SRC="$(cd "$(dirname "$0")" && pwd)"
RERUN="${RERUN:-$HOME/Rexach_rerun}"
MMC1="${MMC1:-}"
LOG="$SRC/22b_download_rexach_inputs_log.txt"
H5AD=aa87f914-07e0-48aa-b915-8de906c95baf.h5ad
H5AD_URL="https://datasets.cellxgene.cziscience.com/$H5AD"
HF="https://huggingface.co/datasets/kahkengwong/CEP-IP_with_CMRF_in_SEA-AD/resolve/main/external_validation/Rexach"
# file name | size in bytes (of the files used for the paper)
FILES=(
  "Part12c_cell_metrics.csv|35453578"
  "R13_mapmycells_output/R12_Rexach_MapMyCells_part1_10xHumanMTGSEA-AD(CCN20230505)_DeepGenerativeMapping_UTC_1792896410850.csv|11991951"
  "R13_mapmycells_output/R12_Rexach_MapMyCells_part1_10xHumanMTGSEA-AD(CCN20230505)_SEA-ADHierarchicalMapping_UTC_1792896410850.csv|13193219"
  "R13_mapmycells_output/R12_Rexach_MapMyCells_part2_10xHumanMTGSEA-AD(CCN20230505)_DeepGenerativeMapping_UTC_1792896470400.csv|11647013"
  "R13_mapmycells_output/R12_Rexach_MapMyCells_part2_10xHumanMTGSEA-AD(CCN20230505)_SEA-ADHierarchicalMapping_UTC_1792896470400.csv|12772135"
  "R13_mapmycells_output/R12_Rexach_MapMyCells_part3_10xHumanMTGSEA-AD(CCN20230505)_DeepGenerativeMapping_UTC_1792896477547.csv|10262359"
  "R13_mapmycells_output/R12_Rexach_MapMyCells_part3_10xHumanMTGSEA-AD(CCN20230505)_SEA-ADHierarchicalMapping_UTC_1792896477547.csv|11311396"
)
H5AD_SIZE=4993733632

: > "$LOG"
exec > >(tee -a "$LOG") 2>&1
echo "== 22b_download_rexach_inputs.sh started $(date '+%Y-%m-%d %H:%M:%S')"
echo "Rerun folder: $RERUN"
command -v curl >/dev/null 2>&1 || { echo "STOPPED: curl not found (sudo apt install curl)"; exit 1; }
mkdir -p "$RERUN/R13_mapmycells_output"
if find "$RERUN" -mindepth 1 -maxdepth 1 ! -name '*.h5ad' ! -name 'R13_mapmycells_output' ! -name '*.csv' ! -name 'mmc1.xlsx' \
        | grep -q .; then
    echo "STOPPED: $RERUN holds outputs of an earlier run. Give a new folder with RERUN=..."
    exit 1
fi
free_kb=$(df --output=avail -k "$RERUN" | tail -1 | tr -dc '0-9')
awk -v f="$free_kb" 'BEGIN {printf "Free on the Linux disk: %.1f GB (needed: about 15 GB)\n", f / 1048576}'
[ "$free_kb" -ge $((15 * 1048576)) ] || { echo "STOPPED: less than 15 GB free"; exit 1; }

size_of() { stat -c %s "$1" 2>/dev/null || echo 0; }
fetch() {   # fetch <url> <destination> <size>
    local url=$1 dst=$2 want=$3
    if [ "$(size_of "$dst")" = "$want" ]; then echo "  already complete: ${dst#$RERUN/}"; return 0; fi
    echo "  downloading ${dst#$RERUN/}"
    curl -fL --retry 5 --retry-delay 10 --retry-all-errors -C - --progress-bar -o "$dst" "$url" || return 1
    [ "$(size_of "$dst")" = "$want" ]
}

echo
echo "== Hugging Face files (107 MB)"
fail=0
for e in "${FILES[@]}"; do
    f=${e%%|*}; s=${e##*|}
    url="$HF/$(echo "$f" | sed 's/(/%28/g; s/)/%29/g')"
    if fetch "$url" "$RERUN/$f" "$s"; then echo "  [PASS] $s bytes  $f"
    else echo "  [FAIL] ${f}: $(size_of "$RERUN/$f") bytes, expected $s"; fail=$((fail + 1)); fi
done

echo
echo "== CELLxGENE h5ad (5.0 GB; about 10 to 60 min, depending on the connection)"
if fetch "$H5AD_URL" "$RERUN/$H5AD" "$H5AD_SIZE"; then echo "  [PASS] $H5AD_SIZE bytes  $H5AD"
else
    echo "  [FAIL] $H5AD: $(size_of "$RERUN/$H5AD") bytes, expected $H5AD_SIZE"
    echo "         Run the script again to resume. If it keeps failing, download the h5ad of dataset"
    echo "         aa87f914-07e0-48aa-b915-8de906c95baf from https://cellxgene.cziscience.com and place"
    echo "         it in $RERUN as $H5AD"
    fail=$((fail + 1))
fi

echo
echo "== Table S1 (mmc1.xlsx)"
if [ -f "$RERUN/mmc1.xlsx" ]; then echo "  [PASS] mmc1.xlsx is in the rerun folder"
elif [ -n "$MMC1" ] && [ -f "$MMC1" ]; then cp "$MMC1" "$RERUN/mmc1.xlsx" && echo "  [PASS] copied $MMC1"
else
    echo "  [FAIL] mmc1.xlsx not found. Download Table S1 from the paper's page (Rexach et al. 2024,"
    echo "         Cell 187:5753) and run again with MMC1=/path/to/mmc1.xlsx"
    fail=$((fail + 1))
fi

echo
echo "== Rerun folder"
(cd "$RERUN" && find . -type f -printf '%10s  %p\n' | sort -k2)
if [ "$fail" -eq 0 ]; then
    echo "PASS: one h5ad, 6 MapMyCells csv files, mmc1.xlsx, and Part12c_cell_metrics.csv"
else
    echo "STOPPED: $fail input(s) missing or incomplete (see above)"
fi
echo "== 22b_download_rexach_inputs.sh finished $(date '+%Y-%m-%d %H:%M:%S')"
exit "$fail"
