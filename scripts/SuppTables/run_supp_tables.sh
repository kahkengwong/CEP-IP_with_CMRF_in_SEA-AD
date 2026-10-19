#!/usr/bin/env bash
# run_supp_tables.sh
# Builds Supplementary Tables from the pipeline outputs under the data folder,
# writes them in their published form, and compares them with the published
# tables when those are available.
#
# Usage (inside the container, data folder at /data):
#   bash run_supp_tables.sh                  # S1 to S19
#   bash run_supp_tables.sh 6-12             # selected tables, for example 6-12 or 7,9
#   PUBLISHED=/ref/Supp_Tables bash run_supp_tables.sh 6-12   # and compare
#
# S18 runs Part 13c (scripts/R), which also draws Supplementary Figures S20
# and S21. S19 is finalized from Fig17_Summary.xlsx, which Fig17_Summary_code.r
# writes, so it must run after that script.
# Each builder writes a log beside its workbook in Manuscript/Supp_Tables/build/.
# A table whose builder fails is not finalized. The summary is saved to
# Manuscript/Supp_Tables/run_supp_tables_log.txt.

set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="${SEAAD_ROOT:-/data}"
export SEAAD_ROOT="$ROOT"
BUILD="$ROOT/Manuscript/Supp_Tables/build"
OUT="$ROOT/Manuscript/Supp_Tables"
LOG="$OUT/run_supp_tables_log.txt"
mkdir -p "$BUILD"
: > "$LOG"
say() { echo "$*" | tee -a "$LOG"; }

# Table numbers from the argument (default 1-19).
expand() {
    local out=() part a b
    IFS=',' read -ra parts <<< "${1:-1-19}"
    for part in "${parts[@]}"; do
        a="${part%-*}"; b="${part#*-}"
        for ((i = a; i <= b; i++)); do out+=("$i"); done
    done
    echo "${out[@]}"
}
TABLES=($(expand "${1:-1-19}"))

py() { python "$HERE/$1" "${@:2}"; }
rs() { Rscript -e "source('$HERE/$1')"; }

# Builder steps of each table, in order.
steps() {
    case "$1" in
        1|2|3|4|5|6) py "build_SuppTable_S$1.py" ;;
        7)  py build_SuppTable_S7.py ;;
        8)  py build_SuppTable_S8_S9_step0_B1_counts.py && py build_SuppTable_S8.py ;;
        9)  { [ -f "$BUILD/B1_group_donor_counts.csv" ] || py build_SuppTable_S8_S9_step0_B1_counts.py; } \
                && py build_SuppTable_S9.py ;;
        10) rs build_SuppTable_S10_step0_Fig8_cells.r && py build_SuppTable_S10.py --allow-figure-drift ;;
        11) rs build_SuppTable_S11.r ;;
        12) py build_SuppTable_S12.py ;;
        13|14|15|16|17) py "build_SuppTable_S$1.py" ;;
        18) rs ../R/Part13c_CP10K_M4_SuppConsolidate_code.r ;;
        19) return 0 ;;
        *)  echo "no builder for S$1"; return 1 ;;
    esac
}

say "run_supp_tables.sh started $(date '+%Y-%m-%d %H:%M:%S') | root $ROOT | tables ${TABLES[*]}"
OK=()
FAILED=()
FIN_FAIL=0
for t in "${TABLES[@]}"; do
    # Remove a workbook from an earlier run, so a failed build cannot leave it behind.
    rm -f "$BUILD"/SuppTable_S"$t"_*.xlsx
    say ""
    say "== S$t"
    steps "$t" 2>&1 | tee -a "$LOG" | tail -3
    rc=${PIPESTATUS[0]}
    # An R builder can stop without a non-zero exit code; S1 to S12 must also leave their workbook.
    if [ "$t" -le 12 ] && ! ls "$BUILD"/SuppTable_S"$t"_*.xlsx > /dev/null 2>&1; then rc=1; fi
    if [ "$rc" -eq 0 ]; then
        OK+=("$t")
        say "   S$t builder: done"
    else
        FAILED+=("$t")
        say "   S$t builder: FAILED (see the build log in $BUILD)"
    fi
done

if [ "${#OK[@]}" -gt 0 ]; then
    LIST="$(IFS=,; echo "${OK[*]}")"
    say ""
    say "== finalize S$LIST"
    python "$HERE/finalize_supp_tables.py" --root "$ROOT" --tables "$LIST" 2>&1 | tee -a "$LOG" | grep -E "^S[0-9]|FAIL|Written|Not written|VERDICT"
    [ "${PIPESTATUS[0]}" -eq 0 ] || FIN_FAIL=1
    if [ -n "${PUBLISHED:-}" ]; then
        say ""
        say "== compare with $PUBLISHED"
        python "$HERE/compare_supp_tables.py" --rebuilt "$OUT" --published "$PUBLISHED" --tables "$LIST" 2>&1 | tee -a "$LOG" | tail -12
    fi
fi

say ""
say "Builders done: ${OK[*]:-none}"
say "Builders failed: ${FAILED[*]:-none}"
say "Summary saved to $LOG"
# Exit code 1 when a builder or the finalize step failed (used by seaad-run).
[ "${#FAILED[@]}" -eq 0 ] && [ "$FIN_FAIL" -eq 0 ] || exit 1
