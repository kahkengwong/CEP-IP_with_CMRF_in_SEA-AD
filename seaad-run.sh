#!/usr/bin/env bash
# seaad-run
#
# Runs the SEA-AD CEP-IP/CMRF analysis in order inside the container, from the
# SEA-AD h5ad to the 19 supplementary tables, with one command.
#
# Usage (the data folder is mounted at /data):
#   seaad-run all            Python Parts 1 to 8, R Parts 9 to 13, Figure 17 and
#                            Supplementary Tables S1 to S19 (about 6 hours)
#   seaad-run python         Python Parts 1 to 8 and Tables S1 to S6
#   seaad-run r              R Parts 9 to 13, Figure 17 and Tables S7 to S19
#                            (needs the Part 8 export, the Part 2 donor tables and the
#                            Part 7 contamination table, for example from Hugging Face)
#   seaad-run hf             R Parts 9 to 11 and Tables S7 to S10 and S12, from the
#                            Hugging Face files only (no h5ad, no Part 12 inputs)
#   seaad-run tables         Tables S1 to S19 from outputs already in the data folder
#   seaad-run list [mode]    print the steps of a mode (default all)
# Options:
#   --from <step>            start at a named step (as printed by list), to resume
#   --compare                compare the results with the original results mounted
#                            read-only at /ref (the folder holding the original Parts and
#                            Manuscript/Supp_Tables)
#
# Each step's console output is saved to seaad-run_logs/<step>.txt, and the
# run summary (start, end, minutes and result per step) to
# seaad-run_logs/seaad-run_summary.txt. The run stops at the first step that
# fails, because later steps read earlier outputs, and prints the command to
# resume. A step fails when it exits with an error, or when its log shows a
# Python traceback, a failed check ([FAIL], GUARDRAIL FAILED, VERDICT: FAIL) or a
# stop (Run stopped, STOPPED, ERROR:, Execution halted).
#
# Part 7 (DecontX) needs the container started with --shm-size=16g.
# Part12_06_GerritsCells uses P12_WORKERS_06 workers (default 6), since each
# worker holds one raw Gerrits matrix in memory; the results do not depend on it.

set -uo pipefail
H="${SEAAD_ROOT:-/data}"
PYDIR=/opt/seaad/Python
RDIR=/opt/seaad/R
TDIR=/opt/seaad/SuppTables
VDIR=/opt/seaad/verification
REF="${REF_ROOT:-/ref}"
LOGDIR="$H/seaad-run_logs"
P8=Part8_Data-export-CP10K_MT-CO3_8DFG
P12_WORKERS_06=${P12_WORKERS_06:-6}

##### Steps #####
# Each line: modes | kind | step name
#   PY   python /opt/seaad/Python/<name>.py
#   P8   Part 8 for one arm (the deposited script exports HKG8_matched; the other
#        arms use a copy with that arm's settings block uncommented)
#   R    Rscript, source("/opt/seaad/R/<name>.r")
#   P12  Rscript, source("Part12_rerun/<name>.r");  P12PY  python Part12_rerun/<name>.py
#   TAB  one supplementary table (scripts in /opt/seaad/SuppTables), named SuppTable_S<n>
STEPS=(
  "all python      | PY  | Part1a_NormForensics_code"
  "all python      | PY  | Part1b_MT-audit_code"
  "all python      | PY  | Part1c_WithinDonor_code"
  "all python      | PY  | Part1d_Errata_code"
  "all python      | PY  | Part1e_TechQC_code"
  "all python      | PY  | Part1f_Ordinals_code"
  "all python      | PY  | Part1g_TrendReconcile_code"
  "all python      | PY  | Part2_NeoC_CP10K_code"
  "all python      | PY  | Part2_NonNeoC_CP10K_code"
  "all python      | PY  | Part3_M4Composition_code"
  "all python      | PY  | Part4_NeoC_ScreenScatter"
  "all python      | PY  | Part4_NonNeoC_ScreenScatter"
  "all python      | PY  | Part5_NeoC_DEGs"
  "all python      | PY  | Part5_NonNeoC_DEGs"
  "all python      | PY  | Part5b_NeoC_Validation"
  "all python      | PY  | Part5b_NonNeoC_Validation"
  "all python      | PY  | Part6_NeoC_GOI-DFG_code"
  "all python      | PY  | Part6_NonNeoC_GOI-DFG_code"
  "all python      | PY  | Part6b_NeoC_Pooled_GOI-DFG"
  "all python      | PY  | Part6b_NonNeoC_Pooled_GOI-DFG"
  "all python      | PY  | Part7_DecontX"
  "all python      | P8  | Part8_PVM3_all"
  "all python      | P8  | Part8_PVM1_all"
  "all python      | P8  | Part8_Mic13_PVM3"
  "all python      | P8  | Part8_HKG8_matched"
  "all python tables | TAB | SuppTable_S1"
  "all python tables | TAB | SuppTable_S2"
  "all python tables | TAB | SuppTable_S3"
  "all python tables | TAB | SuppTable_S4"
  "all python tables | TAB | SuppTable_S5"
  "all python tables | TAB | SuppTable_S6"
  "all r hf        | R   | Part9_CP10K_QC_GAM_code"
  "all r hf        | R   | Part9b_CP10K_SiKi_code"
  "all r hf        | R   | Part10_CP10K_SKCS_code"
  "all r hf        | R   | Part10b_CP10K_IPRS_code"
  "all r hf        | R   | Part10b_CP10K_SKCS_TREP_code"
  "all r hf        | R   | Part10c_CP10K_Copula1_code"
  "all r hf        | R   | Part10d_CP10K_Copula2_code"
  "all r hf        | R   | Part10e_CP10K_CMRF_code"
  "all r hf        | R   | Part10e_CP10K_Validation1_code"
  "all r hf        | R   | Part10e_CP10K_Validation2_UMAP_code"
  "all r hf        | R   | Part10f_CP10K_Band_Formation_code"
  "all r hf        | R   | Part10f_CP10K_CMRF_RSkfCV_Mechanism_code"
  "all r hf        | R   | Part10g_CP10K_DEGs_Part1_code"
  "all r hf        | R   | Part10g_Step7to9_rebuild_code"
  "all r hf        | R   | Part10h_CP10K_Prognosis_code"
  "all r hf        | R   | Part10h_CP10K_Prognosis_1_Comparator_code"
  "all r hf        | R   | Part10h_CP10K_Prognosis_2_CogVsPath_code"
  "all r hf        | R   | Part10h_CP10K_Prognosis_4_Multivariate_code"
  "all r hf        | R   | Part10h_CP10K_Comparator_Fig13_code"
  "all r hf        | R   | Part11_CP10K_CogVsPath_code"
  "all r hf        | R   | Part11_CP10K_CogVsPath_HighLevel_code"
  "all r hf        | R   | Part11_Identity_BandAblation_code"
  "all r hf        | R   | Part11_BandAblation_Fig11Export_code"
  "all r hf tables | TAB | SuppTable_S7"
  "all r hf tables | TAB | SuppTable_S8"
  "all r hf tables | TAB | SuppTable_S9"
  "all r hf tables | TAB | SuppTable_S10"
  "all r hf tables | TAB | SuppTable_S12"
  "all r           | P12PY | Part12_02_Extract"
  "all r           | P12 | Part12_04_Labels"
  "all r           | P12 | Part12_05_Units"
  "all r           | P12 | Part12_06_GerritsCells"
  "all r           | P12 | Part12_Fig15_Analysis"
  "all r           | P12 | Part12_Fig15_Plots"
  "all r           | P12 | Part12_Fig16_Analysis"
  "all r           | P12 | Part12_Fig16_Plots"
  "all r           | P12 | Part12_Fig17_Analysis"
  "all r           | P12 | Part12_Fig17_Plots"
  "all r           | P12 | Part12_Fig18_Analysis"
  "all r           | P12 | Part12_Fig18_Plots"
  "all r           | P12 | Part12_Fig19_Analysis"
  "all r           | P12 | Part12_Fig19_Plots"
  "all r           | P12 | Part12_ModulePostIP_Overview"
  "all r tables    | TAB | SuppTable_S11"
  "all r           | R   | Part10g_CP10K_DEGs_Part1_code_S11_followup_analyses"
  "all r           | R   | Part13a_CP10K_Find14DpctMT_code"
  "all r           | R   | Part13_CP10K_M4Share_code"
  "all r           | R   | Part13b_CP10K_M4RandomSubset_code"
  "all r tables    | TAB | SuppTable_S18"
  "all r tables    | TAB | SuppTable_S13"
  "all r tables    | TAB | SuppTable_S14"
  "all r tables    | TAB | SuppTable_S15"
  "all r tables    | TAB | SuppTable_S16"
  "all r tables    | TAB | SuppTable_S17"
  "all r           | R   | Fig17_Summary_code"
  "all r tables    | TAB | SuppTable_S19"
)
# SuppTable_S18 runs Part13c_CP10K_M4_SuppConsolidate_code.r, which also draws
# Supplementary Figures S20 and S21.

field() { echo "$1" | cut -d'|' -f"$2" | xargs; }

usage() { sed -n '2,33p' "$0" | sed 's/^# \{0,1\}//'; }

##### Arguments #####
MODE=""
FROM=""
COMPARE=0
LIST=0
while [ $# -gt 0 ]; do
    case "$1" in
        all|python|r|hf|tables) MODE="$1" ;;
        list) LIST=1 ;;
        --from) shift; FROM="${1:-}" ;;
        --compare) COMPARE=1 ;;
        -h|--help|help) usage; exit 0 ;;
        *) echo "Unknown argument: $1"; echo; usage; exit 1 ;;
    esac
    shift
done
if [ "$LIST" -eq 1 ]; then
    MODE="${MODE:-all}"
    n=0
    for s in "${STEPS[@]}"; do
        case " $(field "$s" 1) " in *" $MODE "*)
            n=$((n + 1)); printf "%3d  %-6s %s\n" "$n" "$(field "$s" 2)" "$(field "$s" 3)" ;;
        esac
    done
    exit 0
fi
if [ -z "$MODE" ]; then usage; exit 1; fi

if [ ! -d "$H" ]; then
    echo "Data folder $H not found. Start the container with -v /path/to/SEA-AD:/data:Z"
    exit 1
fi
cd "$H" || exit 1
mkdir -p "$LOGDIR"
SUM="$LOGDIR/seaad-run_summary.txt"
say() { echo "$*" | tee -a "$SUM"; }

##### Inputs #####
H5AD=SEAAD_Immune_HIP_MEC_LEC_ITG_MTG_FI_STG_DFC_AnG_V1C_RNAseq_final-nuclei.2026-06-22.h5ad
missing=()
need() { [ -e "$H/$1" ] || missing+=("$1  ($2)"); }
need_pattern() {
    [ -n "$(find "$H" -maxdepth 4 -type f -name "$1" -print -quit 2>/dev/null)" ] || missing+=("$1  ($2)")
}
case "$MODE" in
    all|python)
        need "$H5AD" "SEA-AD immune h5ad"
        need "MetaData/sea-ad_cohort_donor_metadata_072524.xlsx" "SEA-AD donor metadata, Parts 1d and 1e" ;;
esac
case "$MODE" in
    all|r)
        need_pattern "sea-ad-quantitative-neuropathology-063026.xlsx" "SEA-AD neuropathology workbook, Part 12"
        need_pattern "sea-ad_global_and_local_cps_050126.xlsx" "SEA-AD CPS workbook, Part 12"
        need "Gerrits/Metadata/GSE148822_series_matrix.txt.gz" "Gerrits et al. GSE148822, Part 12"
        [ "$MODE" = r ] && need "$H5AD" "SEA-AD immune h5ad, read by Part12_02" ;;
esac
case "$MODE" in
    r|hf)
        need "R_export_CP10K" "Part 8 export"
        need_pattern "*_S5b_donor_clinical.csv" "Part 2 donor tables"
        need "Part7_DecontX/decontX_contamination_all.csv" "Part 7 contamination table" ;;
esac
if [ "${#missing[@]}" -gt 0 ]; then
    echo "Inputs missing from $H:"
    for m in "${missing[@]}"; do echo "  $m"; done
    echo "See README_container.md, section 3."
    exit 1
fi
if [ "$COMPARE" -eq 1 ] && [ ! -d "$REF" ]; then
    echo "--compare needs the original results mounted at $REF (for example -v /path/to/original/SEA-AD:/ref:ro,Z)."
    exit 1
fi

##### Part 8 arm copies #####
make_part8() {   # writes seaad-run_logs/part8_scripts/Part8_<arm>.py with that arm's settings
mkdir -p "$LOGDIR/part8_scripts"
python - "$1" "$PYDIR/$P8.py" "$LOGDIR/part8_scripts/Part8_$1.py" <<'EOF'
import re, sys
arm, src, dst = sys.argv[1:4]
lines = open(src, encoding="utf-8").read().split("\n")
names = ["ARM_TAG", "TARGET_STYPE", "USE_M4_TOP", "SELECT_GENES", "PANEL_GOI",
         "PANEL_DFG", "PANEL_LABEL", "RUN_SUFFIX", "NEST_SUBSET_OF"]
i = [k for k, l in enumerate(lines) if re.match(r'^# ARM_TAG\s*=\s*"%s"\s*$' % re.escape(arm), l)]
if len(i) != 1:
    sys.exit("arm block for %s not found once in %s" % (arm, src))
block = lines[i[0]:i[0] + len(names)]
for l, n in zip(block, names):
    if not re.match(r"^# %s\s*=" % n, l):
        sys.exit("unexpected line in the %s block: %r" % (arm, l))
block = [l[2:] for l in block]
p = [k for k, l in enumerate(lines) if l.startswith("# == Paths")]
if len(p) != 1:
    sys.exit("'# == Paths' not found once in %s" % src)
ins = ["# seaad-run: arm %s, settings copied from its commented block above" % arm] + block + [""]
open(dst, "w", encoding="utf-8").write("\n".join(lines[:p[0]] + ins + lines[p[0]:]))
print("  arm block used:")
for l in block:
    print("    " + l)
EOF
}

run_step() {     # kind name log
    local kind=$1 s=$2 log=$3
    case $kind in
        PY)    python -u "$PYDIR/$s.py" ;;
        P8)    local arm=${s#Part8_}
               if [ "$arm" = HKG8_matched ]; then
                   python -u "$PYDIR/$P8.py"
               else
                   make_part8 "$arm" && python -u "$LOGDIR/part8_scripts/Part8_$arm.py"
               fi ;;
        R)     Rscript -e "source('$RDIR/$s.r')" ;;
        P12)   local pre=""
               [ "$s" = Part12_06_GerritsCells ] && pre="P12_WORKERS <- ${P12_WORKERS_06}L; "
               Rscript -e "${pre}source('Part12_rerun/$s.r')" ;;
        P12PY) python -u "Part12_rerun/$s.py" ;;
        TAB)   bash "$TDIR/run_supp_tables.sh" "${s#SuppTable_S}" ;;
    esac > "$log" 2>&1
}

##### Run #####
say "seaad-run $MODE started $(date '+%Y-%m-%d %H:%M:%S') | data $H"
[ -n "$FROM" ] && say "Starting at: $FROM"
say "R: $(R --version | head -1) | Python: $(python --version 2>&1)"
say ""
say "== seaad-setup"
seaad-setup > "$LOGDIR/seaad-setup.txt" 2>&1
grep -E "DIFFERENT|missing" "$LOGDIR/seaad-setup.txt" | tee -a "$SUM"

started=1
[ -n "$FROM" ] && started=0
n_done=0
T_RUN=$(date +%s)
for line in "${STEPS[@]}"; do
    case " $(field "$line" 1) " in *" $MODE "*) ;; *) continue ;; esac
    kind=$(field "$line" 2)
    s=$(field "$line" 3)
    if [ "$started" -eq 0 ]; then
        [ "$s" = "$FROM" ] && started=1 || continue
    fi
    log="$LOGDIR/$s.txt"
    t0=$(date +%s)
    say "== $s   start $(date '+%H:%M:%S')"
    run_step "$kind" "$s" "$log"
    rc=$?
    if [ "$rc" -eq 0 ] && grep -qE 'Traceback \(most recent call last\)|GUARDRAIL FAILED|\[FAIL\]|^VERDICT: FAIL|Run stopped|STOPPED|^ERROR:|Execution halted' "$log"; then
        rc=2      # finished, but its own checks failed or it stopped
    fi
    mins=$(( ($(date +%s) - t0 + 30) / 60 ))
    if [ "$rc" -ne 0 ]; then
        say "   FAILED $(date '+%H:%M:%S') | ${mins} min | exit $rc | log: ${log#$H/}"
        tail -15 "$log" | sed 's/^/   | /' | tee -a "$SUM"
        say ""
        say "Run stopped at $s. After fixing the cause, resume with:"
        say "   seaad-run $MODE --from $s$([ "$COMPARE" -eq 1 ] && echo " --compare")"
        exit 1
    fi
    n_done=$((n_done + 1))
    say "   finished $(date '+%H:%M:%S') | ${mins} min"
done
if [ "$started" -eq 0 ]; then
    say "Step not found in mode $MODE: $FROM (seaad-run list $MODE prints the steps)"
    exit 1
fi
say ""
say "$n_done steps finished in $(( ($(date +%s) - T_RUN + 30) / 60 )) min"

##### Comparison with the original results #####
if [ "$COMPARE" -eq 1 ]; then
    say ""
    say "== compare with the original results in $REF   start $(date '+%H:%M:%S')"
    case "$MODE" in
        all|python)
            REF_ROOT="$REF" python "$VDIR/compare_python.py" > "$LOGDIR/compare_python.txt" 2>&1
            say "   Python Parts: compare/compare_python_tables_log.txt and compare_python_matrices_log.txt" ;;
    esac
    case "$MODE" in
        all|r|hf)
            parts=$(cd "$H" && ls -d Part9* Part10* Part11* Part12_rerun Part13* 2>/dev/null | tr '\n' ' ')
            WIN_ROOT="$REF" COMPARE_CELLS=1 python "$VDIR/compare_rerun.py" $parts > "$LOGDIR/compare_rerun.txt" 2>&1
            say "   R Parts: compare/compare_<parts>_log.txt (summary in seaad-run_logs/compare_rerun.txt)" ;;
    esac
    if [ -d "$REF/Manuscript/Supp_Tables" ]; then
        case "$MODE" in python) tl=1-6 ;; hf) tl=7-10,12 ;; r) tl=7-19 ;; *) tl=1-19 ;; esac
        python "$TDIR/compare_supp_tables.py" --rebuilt "$H/Manuscript/Supp_Tables" \
            --published "$REF/Manuscript/Supp_Tables" --tables "$tl" > "$LOGDIR/compare_supp_tables.txt" 2>&1
        say "   Supplementary Tables: Manuscript/Supp_Tables/compare_supp_tables_log.txt"
        sed -n '/^Sheets per result:/,$p' "$LOGDIR/compare_supp_tables.txt" | head -10 | tee -a "$SUM"
    else
        say "   $REF/Manuscript/Supp_Tables not found; the supplementary tables were not compared."
    fi
fi
say ""
say "seaad-run $MODE ended $(date '+%Y-%m-%d %H:%M:%S')"
