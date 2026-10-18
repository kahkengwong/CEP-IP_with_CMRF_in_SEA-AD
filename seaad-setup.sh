#!/usr/bin/env bash
# seaad-setup
#
# Prepares the mounted SEA-AD folder (SEAAD_ROOT, /data in the container) for
# the deposited scripts, and reports which inputs it finds.
#
# Files copied if missing:
#   CEP_arms.R, R_config.R, R_config_HKG8.R          into SEAAD_ROOT
#   Part 12 scripts, settings (.json) and pre-specified readings (.md),
#   and Part12_02_Extract.py                         into SEAAD_ROOT/Part12_rerun
# A file that is already there and differs from the deposited copy is left
# unchanged and reported. Nothing else in the data folder is changed.

set -uo pipefail
H="${SEAAD_ROOT:?SEAAD_ROOT is not set}"
R=/opt/seaad/R
PY=/opt/seaad/Python

if [ ! -d "$H" ]; then
    echo "Data folder $H not found. Start the container with -v /path/to/SEA-AD:/data:Z"
    exit 1
fi

copy_if_absent() {
    local src=$1 dst=$2
    mkdir -p "$(dirname "$dst")"
    if [ ! -e "$dst" ]; then
        cp "$src" "$dst" && echo "  copied     ${dst#$H/}"
    elif cmp -s "$src" "$dst"; then
        echo "  same       ${dst#$H/}"
    else
        echo "  DIFFERENT  ${dst#$H/} (left unchanged; deposited copy: $src)"
    fi
}

echo "Settings files:"
for f in CEP_arms.R R_config.R R_config_HKG8.R; do
    copy_if_absent "$R/$f" "$H/$f"
done

echo
echo "Part 12 files in Part12_rerun:"
for src in "$R"/Part12_* "$PY"/Part12_02_Extract.py; do
    [ -e "$src" ] && copy_if_absent "$src" "$H/Part12_rerun/$(basename "$src")"
done

echo
echo "Inputs:"
check() {
    if [ -e "$H/$1" ]; then echo "  found      $1"; else echo "  missing    $1  ($2)"; fi
}
check_pattern() {
    local n
    n=$(find "$H" -path "$H/Part12_rerun" -prune -o -path "$H/Part12_archive" -prune -o \
            -type f -name "$1" -print 2>/dev/null | wc -l)
    if [ "$n" -gt 0 ]; then echo "  found      $1 ($n file(s))"; else echo "  missing    $1  ($2)"; fi
}
check "R_export_CP10K" "Part 8 export, read by Part 9 onwards"
check_pattern "*_S5b_donor_clinical.csv" "Part 2 donor tables, read by Part 10h and Part 12"
check "Part7_DecontX/decontX_contamination_all.csv" "Part 7, read by Part 12"
check "SEAAD_Immune_HIP_MEC_LEC_ITG_MTG_FI_STG_DFC_AnG_V1C_RNAseq_final-nuclei.2026-06-22.h5ad" "SEA-AD h5ad, read by Parts 1 to 8 and Part12_02"
check_pattern "sea-ad-quantitative-neuropathology-063026.xlsx" "SEA-AD neuropathology, read by Part 12"
check_pattern "sea-ad_global_and_local_cps_050126.xlsx" "SEA-AD CPS workbook, read by Part 12"
check "Gerrits/Metadata/GSE148822_series_matrix.txt.gz" "GSE148822, read by Part12_06 and Battery 1"

echo
echo "Scripts: /opt/seaad/Python and /opt/seaad/R. Run order: README_container.md."
echo "Package versions: /opt/seaad/versions.txt"
