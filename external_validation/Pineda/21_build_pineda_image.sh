#!/usr/bin/env bash
# 21_build_pineda_image.sh
# Step 1 of the Pineda container test in WSL: builds the image pineda-cepip on top of the SEA-AD
# image (localhost/seaad-cepip).
#
#   - Copies this build folder to ~/pineda_build (the Linux disk builds faster than /mnt/c), and
#     removes Windows line endings from the copies, in case an editor added them.
#   - Checks that the SEA-AD image is present. If it is not, and SEAAD_TAR names the Zenodo archive
#     (seaad-cepip_image.tar.gz), it loads it first.
#   - Builds pineda-cepip. The build needs an internet connection (GitHub for cepip, Posit Package
#     Manager for brglm2 and xgboost, Bioconductor for scDblFinder and its packages) and takes about
#     20 to 40 minutes, because the Bioconductor packages are compiled. Its output is shown as it runs.
#   - Checks that the build log holds "check_lock_pineda.R: PASS".
#
# Usage, in WSL:
#   cd /mnt/c/Users/Wong/Downloads/pineda_container
#   bash 21_build_pineda_image.sh
# If the SEA-AD image is not on this computer yet:
#   SEAAD_TAR=/mnt/c/Users/Wong/Downloads/seaad-cepip_image.tar.gz bash 21_build_pineda_image.sh
# Log: 21_build_pineda_image_log.txt in this folder.

set -uo pipefail
SRC="$(cd "$(dirname "$0")" && pwd)"
BUILD="$HOME/pineda_build"
LOG="$SRC/21_build_pineda_image_log.txt"
IMAGE=pineda-cepip
FILES=(Containerfile check_lock_pineda.R install_pineda_packages.R pineda_r_packages.tsv pineda-run.sh
       compare_pineda.py build_SuppTable_S20.py
       R02_Pineda_cell_calling_QC.R R11_Pineda_nuclei_table.R R12_Pineda_mapping_input.py
       R13_Pineda_count_checks.py R13b_Pineda_mapping_comparison.py R13c_Pineda_noAuthorQC_arm.py
       R14_Pineda_core_cepip.R R15_Pineda_diagnosis_contrasts.R R16_Pineda_main_figure.R
       R17_Pineda_posthoc_beta.R README_container_Pineda.md)

: > "$LOG"
exec > >(tee -a "$LOG") 2>&1
echo "== 21_build_pineda_image.sh started $(date '+%Y-%m-%d %H:%M:%S')"

# ---- Container engine ------------------------------------------------------------------------
if command -v podman >/dev/null 2>&1; then ENG=podman
elif command -v docker >/dev/null 2>&1; then ENG=docker
else echo "STOPPED: neither podman nor docker found in WSL. Install podman: sudo apt install podman"; exit 1; fi
echo "Engine: $ENG ($($ENG --version))"

# ---- Build folder ----------------------------------------------------------------------------
missing=()
for f in "${FILES[@]}"; do [ -f "$SRC/$f" ] || missing+=("$f"); done
if [ "${#missing[@]}" -gt 0 ]; then echo "STOPPED: missing in $SRC: ${missing[*]}"; exit 1; fi
mkdir -p "$BUILD"
for f in "${FILES[@]}"; do
    cp "$SRC/$f" "$BUILD/$f"
    if grep -q $'\r' "$BUILD/$f"; then sed -i 's/\r$//' "$BUILD/$f"; echo "  Windows line endings removed: $f"; fi
done
echo "Build folder: $BUILD (${#FILES[@]} files)"

# ---- SEA-AD base image -----------------------------------------------------------------------
BASE=""
for cand in localhost/seaad-cepip:latest seaad-cepip:latest docker.io/library/seaad-cepip:latest; do
    if $ENG image inspect "$cand" >/dev/null 2>&1; then BASE="$cand"; break; fi
done
if [ -z "$BASE" ] && [ -n "${SEAAD_TAR:-}" ]; then
    echo "Loading the SEA-AD image from $SEAAD_TAR (a few minutes)"
    $ENG load -i "$SEAAD_TAR" || { echo "STOPPED: $ENG load failed"; exit 1; }
    for cand in localhost/seaad-cepip:latest seaad-cepip:latest; do
        if $ENG image inspect "$cand" >/dev/null 2>&1; then BASE="$cand"; break; fi
    done
fi
if [ -z "$BASE" ]; then
    echo "STOPPED: the SEA-AD image (seaad-cepip) is not in $ENG on this computer."
    echo "  Download seaad-cepip_image.tar.gz from https://doi.org/10.5281/zenodo.23446228, then run:"
    echo "  SEAAD_TAR=/mnt/c/Users/Wong/Downloads/seaad-cepip_image.tar.gz bash 21_build_pineda_image.sh"
    $ENG images
    exit 1
fi
echo "SEA-AD image: $BASE ($($ENG image inspect --format '{{.Id}}' "$BASE" | cut -c1-19))"
df -h "$HOME" | tail -1 | awk '{print "Free space in home: " $4}'

# ---- Build -----------------------------------------------------------------------------------
echo
echo "== Build (output follows)"
cd "$BUILD" || exit 1
t0=$(date +%s)
$ENG build --build-arg BASE="$BASE" -t "$IMAGE" -f Containerfile .
rc=$?
echo "== Build exit code $rc | $(( ($(date +%s) - t0 + 30) / 60 )) min"
if [ "$rc" -ne 0 ]; then echo "STOPPED: the build failed; see the output above"; exit 1; fi
if grep -q "check_lock_pineda.R: PASS" "$LOG"; then
    echo "check_lock_pineda.R: PASS found"
else
    echo "WARN: 'check_lock_pineda.R: PASS' not found in the build output (a cached layer hides it)."
    echo "      Rebuild with: $ENG build --no-cache --build-arg BASE=$BASE -t $IMAGE -f Containerfile ."
fi
echo
echo "== Versions recorded in the image"
$ENG run --rm "$IMAGE" bash -c 'grep -E "^R version|^Python|^cepip|^brglm2|^scDblFinder|^DropletUtils|^xgboost|^scater|^scran|^mgcv|^openxlsx|^ggplot2" /opt/pineda/versions.txt'
echo
echo "== 21_build_pineda_image.sh finished $(date '+%Y-%m-%d %H:%M:%S') | image $IMAGE"
