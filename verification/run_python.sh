#!/usr/bin/env bash
# =============================================================================
# run_python.sh
# Reruns the Python scripts of Parts 1 to 8 in order inside the container, then
# compares every output with the original results mounted read-only at /ref.
#
# Usage (inside the container, from /data):
#   bash /data/run_python.sh                              # start at the first script
#   bash /data/run_python.sh Part6_NeoC_GOI-DFG_code      # resume at a named step
#   bash /data/run_python.sh compare                      # comparison only
#
# Each script's console output is saved to console_logs_py/<step>_console_log.txt,
# and the run summary (start, end, minutes and exit code per step) to
# console_logs_py/run_python_summary.txt. The run stops at the first step that
# exits with an error, prints a Python traceback, or reports a failed check
# ([FAIL] or GUARDRAIL FAILED), because every later part reads earlier outputs.
#
# Part 8 exports one arm per run, and the arm is chosen by the block left
# uncommented at the top of the script (HKG8_matched in the deposited copy).
# For the other three arms, a copy is written to pyrerun_scripts/ with that
# arm's block, taken from the script's own commented lines, placed after the
# arm settings. The arms run in nesting order (PVM3_all, PVM1_all, Mic13_PVM3,
# HKG8_matched), so the subset checks of Part 8 find their parent arm.
# =============================================================================

cd /data || exit 1
mkdir -p console_logs_py pyrerun_scripts
SUM=console_logs_py/run_python_summary.txt
START=${1:-}
PYDIR=/opt/seaad/Python
P8=Part8_Data-export-CP10K_MT-CO3_8DFG

STEPS=(
  Part1a_NormForensics_code
  Part1b_MT-audit_code
  Part1c_WithinDonor_code
  Part1d_Errata_code
  Part1e_TechQC_code
  Part1f_Ordinals_code
  Part1g_TrendReconcile_code
  Part2_NeoC_CP10K_code
  Part2_NonNeoC_CP10K_code
  Part3_M4Composition_code
  Part4_NeoC_ScreenScatter
  Part4_NonNeoC_ScreenScatter
  Part5_NeoC_DEGs
  Part5_NonNeoC_DEGs
  Part5b_NeoC_Validation
  Part5b_NonNeoC_Validation
  Part6_NeoC_GOI-DFG_code
  Part6_NonNeoC_GOI-DFG_code
  Part6b_NeoC_Pooled_GOI-DFG
  Part6b_NonNeoC_Pooled_GOI-DFG
  Part7_DecontX
  Part8_PVM3_all
  Part8_PVM1_all
  Part8_Mic13_PVM3
  Part8_HKG8_matched
  compare
)

say() { echo "$*" | tee -a "$SUM"; }

make_part8() {   # writes pyrerun_scripts/Part8_<arm>.py and prints the arm block used
python - "$1" "$PYDIR/$P8.py" "pyrerun_scripts/Part8_$1.py" <<'EOF'
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
ins = ["# Rerun: arm %s, settings copied from its commented block above" % arm] + block + [""]
open(dst, "w", encoding="utf-8").write("\n".join(lines[:p[0]] + ins + lines[p[0]:]))
print("  arm block used:")
for l in block:
    print("    " + l)
EOF
}

say "run_python.sh started $(date '+%Y-%m-%d %H:%M:%S')"
[ -n "$START" ] && say "Starting at: $START"
say "Python: $(python --version 2>&1)"
say ""

started=0
[ -z "$START" ] && started=1

for s in "${STEPS[@]}"; do
    if [ "$started" -eq 0 ]; then
        [ "$s" = "$START" ] && started=1 || continue
    fi
    log="console_logs_py/${s}_console_log.txt"
    t0=$(date +%s)
    say "== $s   start $(date '+%H:%M:%S')"
    case $s in
        compare)
            python -u /data/compare_python.py 2>&1 | tee "$log"; rc=${PIPESTATUS[0]} ;;
        Part8_HKG8_matched)
            python -u "$PYDIR/$P8.py" 2>&1 | tee "$log"; rc=${PIPESTATUS[0]} ;;
        Part8_*)
            arm=${s#Part8_}
            make_part8 "$arm" 2>&1 | tee -a "$SUM"; rc=${PIPESTATUS[0]}
            if [ "$rc" -eq 0 ]; then
                python -u "pyrerun_scripts/Part8_$arm.py" 2>&1 | tee "$log"; rc=${PIPESTATUS[0]}
            fi ;;
        *)
            python -u "$PYDIR/$s.py" 2>&1 | tee "$log"; rc=${PIPESTATUS[0]} ;;
    esac
    mins=$(( ($(date +%s) - t0 + 30) / 60 ))
    if [ "$rc" -eq 0 ] && [ -f "$log" ] && \
       grep -qE 'Traceback \(most recent call last\)|GUARDRAIL FAILED|\[FAIL\]' "$log"; then
        rc=2      # finished, but a traceback or a failed check was printed
    fi
    if [ "$rc" -eq 0 ]; then
        say "   finished $(date '+%H:%M:%S') | ${mins} min | exit 0"
    else
        say "   FAILED $(date '+%H:%M:%S') | ${mins} min | exit $rc"
        say ""
        say "Run stopped at $s. After fixing the cause, resume with:"
        say "   bash /data/run_python.sh $s"
        exit 1
    fi
done

if [ "$started" -eq 0 ]; then
    say "Start step not found: $START"
    exit 1
fi
say ""
say "run_python.sh ended $(date '+%Y-%m-%d %H:%M:%S')"
