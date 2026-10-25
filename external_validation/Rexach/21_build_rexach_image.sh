#!/usr/bin/env bash
# 21_build_rexach_image.sh
# Step 1 of the Rexach container test in WSL: builds the image rexach-cepip on top of the SEA-AD
# image (localhost/seaad-cepip).
#
#   - Copies this build folder to ~/rexach_build (the Linux disk builds faster than /mnt/c), and
#     removes Windows line endings from the copies, in case an editor added them.
#   - Checks that the SEA-AD image is present. If it is not, and SEAAD_TAR names the Zenodo archive
#     (seaad-cepip_image.tar.gz), it loads it first.
#   - Builds rexach-cepip. The build needs an internet connection (GitHub for cepip, Posit Package
#     Manager for brglm2) and takes a few minutes. Its output is shown as it runs.
#   - Checks that the build log holds "check_lock.R: PASS".
#
# Usage, in WSL:
#   cd /mnt/c/Users/Wong/Downloads/rexach_container
#   bash 21_build_rexach_image.sh
# If the SEA-AD image is not on this laptop yet:
#   SEAAD_TAR=/mnt/c/Users/Wong/Downloads/seaad-cepip_image.tar.gz bash 21_build_rexach_image.sh
# Log: 21_build_rexach_image_log.txt in this folder.

set -uo pipefail
SRC="$(cd "$(dirname "$0")" && pwd)"
BUILD="$HOME/rexach_build"
LOG="$SRC/21_build_rexach_image_log.txt"
IMAGE=rexach-cepip
FILES=(Containerfile check_lock.R rexach-run.sh compare_rexach.py build_SuppTable_S21.py
       R11_Rexach_cellxgene_check.py R12_Rexach_mapping_input.py R13_Rexach_count_checks.py
       R14_Rexach_core_cepip.R R15_Rexach_diagnosis_contrasts.R R16_Rexach_main_figure.R
       R17_Rexach_posthoc_beta.R README_container_Rexach.md)

: > "$LOG"
exec > >(tee -a "$LOG") 2>&1
echo "== 21_build_rexach_image.sh started $(date '+%Y-%m-%d %H:%M:%S')"

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
    echo "STOPPED: the SEA-AD image (seaad-cepip) is not in $ENG on this laptop."
    echo "  Download seaad-cepip_image.tar.gz from https://doi.org/10.5281/zenodo.23446228, then run:"
    echo "  SEAAD_TAR=/mnt/c/Users/Wong/Downloads/seaad-cepip_image.tar.gz bash 21_build_rexach_image.sh"
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
if grep -q "check_lock.R: PASS" "$LOG"; then
    echo "check_lock.R: PASS found"
else
    echo "WARN: 'check_lock.R: PASS' not found in the build output (a cached layer hides it)."
    echo "      Rebuild with: $ENG build --no-cache --build-arg BASE=$BASE -t $IMAGE -f Containerfile ."
fi
echo
echo "== Versions recorded in the image"
$ENG run --rm "$IMAGE" bash -c 'grep -E "^R version|^Python|^cepip|^brglm2|^mgcv|^openxlsx|^ggplot2" /opt/rexach/versions.txt'
echo
echo "== 21_build_rexach_image.sh finished $(date '+%Y-%m-%d %H:%M:%S') | image $IMAGE"
