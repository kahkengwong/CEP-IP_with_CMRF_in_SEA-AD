#!/usr/bin/env bash
# pineda-run
#
# Runs external validation 1 (Pineda et al. 2024, BA9 batch 210526) in order inside the container,
# from the Cell Ranger raw count matrices to Supplementary Table S20, with one command.
#
# Usage (the Pineda data folder is mounted at /data):
#   pineda-run all            R02 and R11 (R), R12 to R13c (Python), R14 to R17 (R), and Table S20
#   pineda-run r              R14 to R17 and Table S20, from the R13, R13b, and R13c outputs
#   pineda-run list [mode]    print the steps of a mode (default all)
# Options:
#   --from <step>             start at a named step (as printed by list), to resume
#   --compare                 compare the results with the original results mounted read-only
#                             at /ref (the original External_Dataset_1 folder), and Table S20
#                             with /ref_supp (the original Manuscript/Supp_Tables), if mounted
#
# Settings (set with -e on the run line):
#   PINEDA_FIGURE   main figure number in the Table S20 title and file name (default 19)
#   PINEDA_WORKERS  libraries processed at the same time in R02, and workers for cepip() in R14
#                   (default 6; the results do not depend on it; R02 ran with 6 on a 48 GB laptop)
#
# Each step's console output is saved to pineda-run_logs/<step>.txt, and the start, end, minutes,
# and result of every step to pineda-run_logs/pineda-run_summary.txt. The run stops at the first
# step that fails, because later steps read earlier outputs. A step fails when it exits with an
# error, or when its log shows a Python traceback, a failed check ([FAIL], VERDICT: FAIL) or a stop
# (STOPPED, ERROR:, Execution halted).

set -uo pipefail
H="${PINEDA_ROOT:-/data}"
SDIR="${PINEDA_SCRIPTS:-/opt/pineda}"
REF="${REF_ROOT:-/ref}"
REF_SUPP="${REF_SUPP:-/ref_supp}"
LOGDIR="$H/pineda-run_logs"
FIGURE="${PINEDA_FIGURE:-19}"
export PINEDA_ROOT="$H"

##### Steps #####
# Each line: modes | kind | step name
#   R     Rscript, source("/opt/pineda/<name>.R")
#   PY    python /opt/pineda/<name>.py
#   R13B  python /opt/pineda/R13b_Pineda_mapping_comparison.py with the two folders and a label
#   S20   python /opt/pineda/build_SuppTable_S20.py
STEPS=(
  "all     | R    | R02_Pineda_cell_calling_QC"
  "all     | R    | R11_Pineda_nuclei_table"
  "all     | PY   | R12_Pineda_mapping_input"
  "all     | PY   | R13_Pineda_count_checks"
  "all     | R13B | R13b_repeat"
  "all     | R13B | R13b_QC"
  "all     | PY   | R13c_Pineda_noAuthorQC_arm"
  "all r   | R    | R14_Pineda_core_cepip"
  "all r   | R    | R15_Pineda_diagnosis_contrasts"
  "all r   | R    | R17_Pineda_posthoc_beta"
  "all r   | R    | R16_Pineda_main_figure"
  "all r   | S20  | build_SuppTable_S20"
)

field() { echo "$1" | cut -d'|' -f"$2" | xargs; }
usage() { sed -n '2,28p' "$0" | sed 's/^# \{0,1\}//'; }

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
    echo "Data folder $H not found. Start the container with -v /path/to/Pineda_rerun:/data:Z"
    exit 1
fi
mkdir -p "$LOGDIR"
SUMMARY="$LOGDIR/pineda-run_summary.txt"
say() { echo "$*" | tee -a "$SUMMARY"; }
if [ -z "$FROM" ]; then : > "$SUMMARY"; fi
say "pineda-run $MODE started $(date '+%Y-%m-%d %H:%M:%S') | data $H | workers ${PINEDA_WORKERS:-6}"
say "R: $(R --version | head -1) | Python: $(python --version 2>&1)"
say "cepip: $(Rscript -e 'd <- packageDescription("cepip"); cat(d$Version, d$RemoteSha)' 2>/dev/null)"
say "scDblFinder: $(Rscript -e 'cat(as.character(packageVersion("scDblFinder")))' 2>/dev/null)"

##### Inputs #####
say ""
say "== Inputs"
missing=()
need() {   # need <path relative to /data> <description>
    if [ -e "$H/$1" ]; then say "  found      $1"; else say "  missing    $1  ($2)"; missing+=("$1"); fi
}
need_n() {   # need_n <folder> <find pattern> <number> <description>
    local n
    n=$(find "$H/$1" -maxdepth 1 -type f -name "$2" 2>/dev/null | wc -l)
    if [ "$n" -eq "$3" ]; then say "  found      $1/ ($3 $4)"
    else say "  missing    $1/ with $3 $4 (found $n)"; missing+=("$1"); fi
}
if [ "$MODE" = all ]; then
    need_n counts '*_raw_feature_bc_matrix.h5' 63 "Cell Ranger raw count matrices"
    need "05_cellranger_metrics.csv" "Cell Ranger metrics of the 63 libraries, read by R02"
    n_s01=$(find "$H" -maxdepth 1 -type f -iname 'S01*.xlsx' | wc -l)
    if [ "$n_s01" -eq 1 ]; then say "  found      S01 workbook (Synapse syn51105515)"
    else say "  missing    one S01*.xlsx (Synapse syn51105515), read by R11 and R15 (found $n_s01)"; missing+=("S01"); fi
    n_meta=$(find "$H/Pineda2024_metadata" -type f -name '*_col_metadata.tsv' 2>/dev/null | wc -l)
    if [ "$n_meta" -ge 63 ]; then say "  found      Pineda2024_metadata/ ($n_meta authors' <LIB>_col_metadata.tsv files)"
    else say "  missing    Pineda2024_metadata/ with the 63 authors' <LIB>_col_metadata.tsv files (found $n_meta)"; missing+=("metadata"); fi
    need "Stats_15H.csv" "SEA-AD Figure 15H source data, read by R13"
    need_n R13_mapmycells_output '*.csv' 12 "MapMyCells csv files of the authors' QC nuclei"
    need_n R13_mapmycells_output_allKept '*.csv' 26 "MapMyCells csv files of all R02 kept nuclei"
    need_n R13_mapmycells_output_allKept_repeat '*.csv' 26 "MapMyCells csv files of the repeat mapping"
fi
if [ "$MODE" = r ]; then
    need "R13_count_checks/R13_Pineda_MicroPVM3_counts.mtx.gz" "R13 counts"
    need "R13_count_checks/R13_Pineda_MicroPVM3_genes.csv.gz" "R13 genes"
    need "R13_count_checks/R13_Pineda_MicroPVM3_nuclei.csv.gz" "R13 nuclei"
    need "R13_count_checks/R13_Pineda_count_checks.xlsx" "R13 workbook, read by R16"
    need "R13b_mapping_comparison_repeat/R13b_Pineda_mapping_comparison.xlsx" "R13b repeat run, read by Table S20"
    need "R13b_mapping_comparison_QC/R13b_Pineda_mapping_comparison.xlsx" "R13b QC before mapping, read by Table S20"
    need "R13c_noAuthorQC_arm/R13c_Pineda_noAuthorQC_arm.xlsx" "R13c, read by Table S20"
    n_s01=$(find "$H" -maxdepth 1 -type f -iname 'S01*.xlsx' | wc -l)
    [ "$n_s01" -eq 1 ] || { say "  missing    one S01*.xlsx, read by R15 (found $n_s01)"; missing+=("S01"); }
fi
need "Pineda_BA9_210526_manifest.tsv" "the 63 libraries"
if [ "${#missing[@]}" -gt 0 ]; then
    say ""
    say "STOPPED: ${#missing[@]} input(s) missing. README_container_Pineda.md section 3 lists the inputs."
    exit 1
fi

##### Run #####
fails_in_log() {   # fails_in_log <log> <step>
    grep -E -q 'Traceback \(most recent call last\)|\[FAIL\]|VERDICT: FAIL|STOPPED|ERROR:|Execution halted' "$1" && return 0
    # R02 records a failed library and goes on with the others; the step fails unless all 63 finished.
    if [ "$2" = R02_Pineda_cell_calling_QC ]; then
        grep -E -q '^  FAILED ' "$1" && return 0
        grep -q 'Processed 63 of 63 libraries' "$1" || return 0
    fi
    return 1
}
run_step() {   # run_step <kind> <name> <log>
    case $1 in
        R)    Rscript -e "source('$SDIR/$2.R')" ;;
        PY)   python "$SDIR/$2.py" ;;
        R13B) if [ "$2" = R13b_repeat ]; then
                  python "$SDIR/R13b_Pineda_mapping_comparison.py" R13_mapmycells_output_allKept \
                      R13_mapmycells_output_allKept_repeat repeat
              else
                  python "$SDIR/R13b_Pineda_mapping_comparison.py" R13_mapmycells_output_allKept \
                      R13_mapmycells_output QC
              fi ;;
        S20)  mkdir -p "$H/Supp_Tables"
              python "$SDIR/build_SuppTable_S20.py" --root "$H" --supp "$H/Supp_Tables" --figure "$FIGURE" ;;
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
    if [ "$rc" -ne 0 ] || fails_in_log "$LOGDIR/$s.txt" "$s"; then
        say "   FAILED $(date '+%H:%M:%S') | $mins min | exit code $rc | log pineda-run_logs/$s.txt"
        say ""
        say "Resume after fixing the cause with: pineda-run $MODE --from $s"
        exit 1
    fi
    say "   finished $(date '+%H:%M:%S') | $mins min"
    n_done=$((n_done + 1))
done
if [ "$started" = 0 ]; then
    say "Step not found in mode $MODE: $FROM (pineda-run list $MODE prints the steps)"; exit 1
fi
say ""
say "== $n_done step(s) finished in $(( ($(date +%s) - T_RUN + 30) / 60 )) min"

##### Compare #####
if [ "$COMPARE" = 1 ]; then
    say ""
    if [ ! -d "$REF" ]; then
        say "== Compare: $REF not mounted; add -v /path/to/original/External_Dataset_1:/ref:ro,Z"
        exit 1
    fi
    say "== Compare with $REF"
    python "$SDIR/compare_pineda.py" "$MODE" > "$LOGDIR/compare_pineda.txt" 2>&1; rc=$?
    tail -n 40 "$LOGDIR/compare_pineda.txt" | tee -a "$SUMMARY"
    [ "$rc" -eq 0 ] || { say "Compare stopped with exit code $rc (log pineda-run_logs/compare_pineda.txt)"; exit 1; }
fi
say ""
say "pineda-run $MODE finished $(date '+%Y-%m-%d %H:%M:%S')"
