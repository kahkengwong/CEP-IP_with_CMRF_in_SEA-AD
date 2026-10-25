#!/usr/bin/env bash
# rexach-run
#
# Runs external validation 2 (Rexach et al. 2024, insula) in order inside the container, from the
# CELLxGENE h5ad to Supplementary Table S21, with one command.
#
# Usage (the Rexach data folder is mounted at /data):
#   rexach-run all            R11 to R13 (Python), R14 to R17 (R), and Table S21
#   rexach-run r              R14 to R17 and Table S21, from the R12 and R13 outputs
#   rexach-run list [mode]    print the steps of a mode (default all)
# Options:
#   --from <step>             start at a named step (as printed by list), to resume
#   --compare                 compare the results with the original results mounted read-only
#                             at /ref (the original External_Dataset_2 folder), and Table S21
#                             with /ref_supp (the original Manuscript/Supp_Tables), if mounted
#
# Settings (set with -e on the run line):
#   REXACH_FIGURE   main figure number in the Table S21 title and file name (default 19)
#   REXACH_WORKERS  workers for cepip() in R14 (default 6; the results do not depend on it)
#
# Each step's console output is saved to rexach-run_logs/<step>.txt, and the start, end,
# minutes, and result of every step to rexach-run_logs/rexach-run_summary.txt. The run stops at
# the first step that fails, because later steps read earlier outputs. A step fails when it exits
# with an error, or when its log shows a Python traceback, a failed check ([FAIL], VERDICT: FAIL)
# or a stop (STOPPED, ERROR:, Execution halted).

set -uo pipefail
H="${REXACH_ROOT:-/data}"
SDIR="${REXACH_SCRIPTS:-/opt/rexach}"
REF="${REF_ROOT:-/ref}"
REF_SUPP="${REF_SUPP:-/ref_supp}"
LOGDIR="$H/rexach-run_logs"
FIGURE="${REXACH_FIGURE:-19}"
export REXACH_ROOT="$H"

##### Steps #####
# Each line: modes | kind | step name
#   PY   python /opt/rexach/<name>.py (R11 and R12 are given the h5ad path)
#   R    Rscript, source("/opt/rexach/<name>.R")
#   S21  python /opt/rexach/build_SuppTable_S21.py
STEPS=(
  "all     | PY  | R11_Rexach_cellxgene_check"
  "all     | PY  | R12_Rexach_mapping_input"
  "all     | PY  | R13_Rexach_count_checks"
  "all r   | R   | R14_Rexach_core_cepip"
  "all r   | R   | R15_Rexach_diagnosis_contrasts"
  "all r   | R   | R17_Rexach_posthoc_beta"
  "all r   | R   | R16_Rexach_main_figure"
  "all r   | S21 | build_SuppTable_S21"
)

field() { echo "$1" | cut -d'|' -f"$2" | xargs; }
usage() { sed -n '2,27p' "$0" | sed 's/^# \{0,1\}//'; }

##### Arguments #####
MODE=""; FROM=""; COMPARE=0; LIST=0
while [ $# -gt 0 ]; do
    case "$1" in
        all|r) MODE="$1" ;;
        list) LIST=1 ;;
        --from) shift; FROM="${1:-}" ;;
        --compare) COMPARE=1 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown argument: $1"; usage; exit 1 ;;
    esac
    shift
done
if [ "$LIST" = 1 ]; then
    m="${MODE:-all}"
    for s in "${STEPS[@]}"; do
        case " $(field "$s" 1) " in *" $m "*) echo "  $(field "$s" 3)" ;; esac
    done
    exit 0
fi
[ -n "$MODE" ] || { usage; exit 1; }

if [ ! -d "$H" ]; then
    echo "Data folder $H not found. Start the container with -v /path/to/Rexach_rerun:/data:Z"
    exit 1
fi
mkdir -p "$LOGDIR"
SUMMARY="$LOGDIR/rexach-run_summary.txt"
say() { echo "$*" | tee -a "$SUMMARY"; }
if [ -z "$FROM" ]; then : > "$SUMMARY"; fi
say "rexach-run $MODE started $(date '+%Y-%m-%d %H:%M:%S') | data $H"
say "R: $(R --version | head -1) | Python: $(python --version 2>&1)"
say "cepip: $(Rscript -e 'd <- packageDescription("cepip"); cat(d$Version, d$RemoteSha)' 2>/dev/null)"

##### Inputs #####
say ""
say "== Inputs"
missing=()
need() {   # need <path relative to /data> <description>
    if [ -e "$H/$1" ]; then say "  found      $1"; else say "  missing    $1  ($2)"; missing+=("$1"); fi
}
H5AD=""
if [ "$MODE" = all ]; then
    mapfile -t h5 < <(find "$H" -maxdepth 1 -type f -name '*.h5ad' | sort)
    if [ "${#h5[@]}" -eq 1 ]; then
        H5AD="${h5[0]}"; say "  found      $(basename "$H5AD")  (CELLxGENE h5ad)"
    else
        say "  missing    one .h5ad in the data folder (found ${#h5[@]})  (CELLxGENE aa87f914-07e0-48aa-b915-8de906c95baf)"
        missing+=("h5ad")
    fi
    n_map=$(find "$H/R13_mapmycells_output" -maxdepth 1 -type f -name '*.csv' 2>/dev/null | wc -l)
    if [ "$n_map" -eq 6 ]; then say "  found      R13_mapmycells_output/ (6 MapMyCells csv files)"
    else say "  missing    R13_mapmycells_output/ with 6 MapMyCells csv files (found $n_map)"; missing+=("mapmycells"); fi
fi
if [ "$MODE" = r ]; then
    need "R12_mapmycells_input/R12_Rexach_nuclei.csv.gz" "R12 nuclei table"
    need "R13_count_checks/R13_Rexach_MicroPVM3_counts.mtx.gz" "R13 counts"
    need "R13_count_checks/R13_Rexach_MicroPVM3_genes.csv.gz" "R13 genes"
    need "R13_count_checks/R13_Rexach_MicroPVM3_nuclei.csv.gz" "R13 nuclei"
fi
need "mmc1.xlsx" "Table S1 of Rexach et al. 2024, read by R15"
need "Part12c_cell_metrics.csv" "SEA-AD %MT per cell, read by R15"
if [ "${#missing[@]}" -gt 0 ]; then
    say ""
    say "STOPPED: ${#missing[@]} input(s) missing. README_container_Rexach.md section 3 lists the inputs."
    exit 1
fi

##### Run #####
fails_in_log() {
    grep -E -q 'Traceback \(most recent call last\)|\[FAIL\]|VERDICT: FAIL|STOPPED|ERROR:|Execution halted' "$1"
}
run_step() {   # run_step <kind> <name> <log>
    case $1 in
        PY)  if [ "$2" = R13_Rexach_count_checks ]; then python "$SDIR/$2.py"
             else python "$SDIR/$2.py" "$H5AD"; fi ;;
        R)   Rscript -e "source('$SDIR/$2.R')" ;;
        S21) mkdir -p "$H/Supp_Tables"
             python "$SDIR/build_SuppTable_S21.py" --root "$H" --supp "$H/Supp_Tables" --figure "$FIGURE" ;;
    esac > "$3" 2>&1
}

say ""
started=0; [ -z "$FROM" ] && started=1
n_done=0; T_RUN=$(date +%s)
cd "$H" || exit 1
for line in "${STEPS[@]}"; do
    case " $(field "$line" 1) " in *" $MODE "*) ;; *) continue ;; esac
    kind=$(field "$line" 2); s=$(field "$line" 3)
    if [ "$started" = 0 ]; then [ "$s" = "$FROM" ] && started=1 || continue; fi
    t0=$(date +%s)
    say "== $s   start $(date '+%H:%M:%S')"
    run_step "$kind" "$s" "$LOGDIR/$s.txt"; rc=$?
    mins=$(( ($(date +%s) - t0 + 30) / 60 ))
    if [ "$rc" -ne 0 ] || fails_in_log "$LOGDIR/$s.txt"; then
        say "   FAILED $(date '+%H:%M:%S') | $mins min | exit code $rc | log rexach-run_logs/$s.txt"
        say ""
        say "Resume after fixing the cause with: rexach-run $MODE --from $s"
        exit 1
    fi
    say "   finished $(date '+%H:%M:%S') | $mins min"
    n_done=$((n_done + 1))
done
if [ "$started" = 0 ]; then
    say "Step not found in mode $MODE: $FROM (rexach-run list $MODE prints the steps)"; exit 1
fi
say ""
say "== $n_done step(s) finished in $(( ($(date +%s) - T_RUN + 30) / 60 )) min"

##### Compare #####
if [ "$COMPARE" = 1 ]; then
    say ""
    if [ ! -d "$REF" ]; then
        say "== Compare: $REF not mounted; add -v /path/to/original/External_Dataset_2:/ref:ro,Z"
        exit 1
    fi
    say "== Compare with $REF"
    python "$SDIR/compare_rexach.py" "$MODE" > "$LOGDIR/compare_rexach.txt" 2>&1; rc=$?
    tail -n 25 "$LOGDIR/compare_rexach.txt" | tee -a "$SUMMARY"
    [ "$rc" -eq 0 ] || { say "Compare stopped with exit code $rc (log rexach-run_logs/compare_rexach.txt)"; exit 1; }
fi
say ""
say "rexach-run $MODE finished $(date '+%Y-%m-%d %H:%M:%S')"
