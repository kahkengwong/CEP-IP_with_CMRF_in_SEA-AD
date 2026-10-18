#!/usr/bin/env bash
# =============================================================================
# run_rest.sh
# Runs the remaining scripts in order inside the container, from Part 10h
# (Prognosis_2) to the Figure 17 summary, then the full comparison with the
# reference results.
#
# Usage (inside the container, from /data):
#   bash /data/run_rest.sh                      # start at the first script
#   bash /data/run_rest.sh Part12_04_Labels     # start at a named script
#
# Each script's console output is saved to console_logs/<script>_console_log.txt.
# A failed script is recorded and the run continues with the next one. A script
# that exits normally but reports VERDICT: FAIL or a stop in its log is also
# recorded as failed.
# Part12_06_GerritsCells runs with P12_WORKERS_06 workers (default 6), since
# each worker holds one raw Gerrits matrix in memory. Each sample has its own
# seed, so the number of workers does not change the results.
# The run summary (start, end, minutes and exit code per script) is saved to
# console_logs/run_rest_summary.txt.
#
# Two scripts read supplementary tables that the deposit does not rebuild:
#   the S11 follow-up reads Manuscript/Supp_Tables/Supplementary_Table_S11_DEGs_ModvsNonMod.xlsx
#   Fig17_Summary reads the S12 workbook and Part12_rerun/07..11_FigNN/SuppTable_FigNN.xlsx
# These are copied from the reference (/win) just before the script and removed
# right after it, so the comparison does not count them as rerun results. A
# workbook not found at its usual path is searched for by name in /win.
# =============================================================================

cd /data || exit 1
mkdir -p console_logs
SUM=console_logs/run_rest_summary.txt
START=${1:-}
P12_WORKERS_06=${P12_WORKERS_06:-6}

STEPS=(
  "R  Part10h_CP10K_Prognosis_2_CogVsPath_code"
  "R  Part10h_CP10K_Prognosis_4_Multivariate_code"
  "R  Part10h_CP10K_Comparator_Fig13_code"
  "R  Part11_CP10K_CogVsPath_code"
  "R  Part11_CP10K_CogVsPath_HighLevel_code"
  "R  Part11_Identity_BandAblation_code"
  "R  Part11_BandAblation_Fig11Export_code"
  "PY Part12_02_Extract"
  "P12 Part12_04_Labels"
  "P12 Part12_05_Units"
  "P12 Part12_06_GerritsCells"
  "P12 Part12_Fig15_Analysis"
  "P12 Part12_Fig15_Plots"
  "P12 Part12_Fig16_Analysis"
  "P12 Part12_Fig16_Plots"
  "P12 Part12_Fig17_Analysis"
  "P12 Part12_Fig17_Plots"
  "P12 Part12_Fig18_Analysis"
  "P12 Part12_Fig18_Plots"
  "P12 Part12_Fig19_Analysis"
  "P12 Part12_Fig19_Plots"
  "P12 Part12_ModulePostIP_Overview"
  "R  Part10g_CP10K_DEGs_Part1_code_S11_followup_analyses"
  "R  Part13a_CP10K_Find14DpctMT_code"
  "R  Part13_CP10K_M4Share_code"
  "R  Part13b_CP10K_M4RandomSubset_code"
  "R  Fig17_Summary_code"
)

S11_IN=("Manuscript/Supp_Tables/Supplementary_Table_S11_DEGs_ModvsNonMod.xlsx")
F17_IN=("Manuscript/Supp_Tables/Supplementary_Table_S12_Prognosis_&_Ablation.xlsx"
        "Part12_rerun/07_Fig15/SuppTable_Fig15.xlsx"
        "Part12_rerun/08_Fig16/SuppTable_Fig16.xlsx"
        "Part12_rerun/09_Fig17/SuppTable_Fig17.xlsx"
        "Part12_rerun/10_Fig18/SuppTable_Fig18.xlsx"
        "Part12_rerun/11_Fig19/SuppTable_Fig19.xlsx")
COPIED=()

say() { echo "$*" | tee -a "$SUM"; }

copy_in() {            # copies reference inputs that are not already in /data
    COPIED=()
    for f in "$@"; do
        if [ -e "/data/$f" ]; then
            say "  input already present: $f"
        elif [ -e "/win/$f" ]; then
            mkdir -p "$(dirname "/data/$f")"
            cp "/win/$f" "/data/$f" && COPIED+=("$f") && say "  input copied from reference: $f"
        else
            src=$(find /win -type f -name "$(basename "$f")" -printf '%T@ %p\n' 2>/dev/null | sort -rn | head -1 | cut -d' ' -f2-)
            if [ -n "$src" ]; then
                mkdir -p "$(dirname "/data/$f")"
                cp "$src" "/data/$f" && COPIED+=("$f") && say "  input copied from reference: $f (found at ${src#/win/})"
            else
                say "  input MISSING in reference: $f"
            fi
        fi
    done
}

remove_copied() {      # removes only the inputs copied by copy_in
    for f in "${COPIED[@]}"; do
        rm -f "/data/$f" && say "  input removed after use: $f"
    done
    rmdir -p --ignore-fail-on-non-empty /data/Manuscript/Supp_Tables 2>/dev/null
    COPIED=()
}

say "run_rest.sh started $(date '+%Y-%m-%d %H:%M:%S')"
[ -n "$START" ] && say "Starting at: $START"
say ""

started=0
[ -z "$START" ] && started=1
n_ok=0; n_fail=0; failed=()

for step in "${STEPS[@]}"; do
    kind=${step%% *}
    s=$(echo "${step#* }" | xargs)
    if [ "$started" -eq 0 ]; then
        [ "$s" = "$START" ] && started=1 || continue
    fi

    [ "$s" = "Part10g_CP10K_DEGs_Part1_code_S11_followup_analyses" ] && copy_in "${S11_IN[@]}"
    [ "$s" = "Fig17_Summary_code" ] && copy_in "${F17_IN[@]}"

    log="console_logs/${s}_console_log.txt"
    t0=$(date +%s)
    say "== $s   start $(date '+%H:%M:%S')"
    case $kind in
        R)   Rscript -e "source('/opt/seaad/R/$s.r')" 2>&1 | tee "$log"; rc=${PIPESTATUS[0]} ;;
        P12) pre=""
             [ "$s" = "Part12_06_GerritsCells" ] && pre="P12_WORKERS <- ${P12_WORKERS_06}L; "
             Rscript -e "${pre}source('Part12_rerun/$s.r')" 2>&1 | tee "$log"; rc=${PIPESTATUS[0]} ;;
        PY)  python "Part12_rerun/$s.py" 2>&1 | tee "$log"; rc=${PIPESTATUS[0]} ;;
    esac
    mins=$(( ($(date +%s) - t0 + 30) / 60 ))
    if [ "$rc" -eq 0 ] && grep -qE '^VERDICT: FAIL|Run stopped|STOPPED|^ERROR:|Execution halted' "$log"; then
        rc=2      # finished, but its own checks stopped or failed it
    fi
    if [ "$rc" -eq 0 ]; then
        n_ok=$((n_ok + 1)); say "   finished $(date '+%H:%M:%S') | ${mins} min | exit 0"
    else
        n_fail=$((n_fail + 1)); failed+=("$s"); say "   FAILED $(date '+%H:%M:%S') | ${mins} min | exit $rc"
    fi

    [ "${#COPIED[@]}" -gt 0 ] && remove_copied
done

if [ "$started" -eq 0 ]; then
    say "Start script not found: $START"
    exit 1
fi

say ""
say "Scripts finished: $n_ok ok, $n_fail failed"
for f in "${failed[@]}"; do say "   failed: $f"; done
say ""
say "== Full comparison   start $(date '+%H:%M:%S')"
COMPARE_CELLS=1 python /data/compare_rerun.py 2>&1 | tee console_logs/compare_all_console_log.txt
say "run_rest.sh ended $(date '+%Y-%m-%d %H:%M:%S')"
