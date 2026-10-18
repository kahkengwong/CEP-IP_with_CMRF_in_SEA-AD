#!/usr/bin/env python3
# =============================================================================
# compare_rerun.py
# Compares a rerun with the reference results, file by file.
#   For every CSV (also .csv.gz) and Excel workbook in the chosen part folders,
#   the matching reference file is read and compared cell by cell:
#   numbers within a relative tolerance, text exactly (cells holding a path,
#   a date or a time, or a settings version, are counted separately).
#   Figures, RDS files and logs are listed by presence only.
#
# Each table is given one of these results:
#   SAME                      all numbers within 1e-6 (relative), text identical
#   SAME UP TO SIGN           the only differences are numbers with the opposite
#                             sign and the same size (eigenvector signs)
#   SAME VALUES, OTHER ORDER  the same rows in a different order
#   TEXT ONLY                 numbers the same, some text differs
#   WITHIN 0.1%               some numbers differ by more than 1e-6 but all by
#                             at most 0.1% (relative)
#   OUTSIDE 0.1%              some numbers differ by more than 0.1%, or a cell
#                             is empty in one table only
#
# Usage: copy this file into the rerun data folder, start the container with
# the reference results mounted read-only at /win, then:
#   python /data/compare_rerun.py Part9 Part9b Part10 ...
#   python /data/compare_rerun.py            # every folder found in both
#
# Outputs (/data/compare/): compare_<parts>_log.txt (summary) and
#   compare_<parts>_detail.csv (one row per file and sheet). With
#   COMPARE_CELLS=1, also compare_<parts>_cells.csv (each differing cell, up
#   to 500 per table). Rows are matched on the first column when it is a
#   unique key in both tables, since tied values can sort differently.
# =============================================================================
import math
import os
import re
import sys
import time

import pandas as pd

NEW = os.environ.get("SEAAD_ROOT", "/data")
WIN = os.environ.get("WIN_ROOT", "/win")
OUT = os.path.join(NEW, "compare")
REL_TOL = 1e-6      # numbers agree when |a - b| <= REL_TOL * max(|a|, |b|) or <= ABS_TOL
ABS_TOL = 1e-9
SKIP_TOP = {"R_export_CP10K", "Gerrits", "compare", "console_logs", "Manuscript",
            "[1]_NeoC", "[2]_Non-NeoC", "Part7_DecontX"}
LARGE_TOL = 1e-3    # numbers within 0.1% (relative) are reported as WITHIN 0.1%
TABLE_EXT = (".csv", ".csv.gz", ".xlsx")
# Columns that record which settings file version wrote the table, not a result.
# Their differences are counted with the path and date cells.
PROVENANCE_COLS = {"config_version"}
ORDER = ["SAME", "SAME UP TO SIGN", "SAME VALUES, OTHER ORDER", "TEXT ONLY",
         "WITHIN 0.1%", "OUTSIDE 0.1%", "DIFFERENT SIZE", "SHEET ONLY IN RERUN",
         "SHEET ONLY IN REFERENCE", "ONLY IN RERUN", "ONLY IN REFERENCE", "READ ERROR"]
# Known causes, matched on "part/file :: sheet". Printed under the matching
# entries of the log. Extend as needed.
EXPLAIN = [
    (r"Validation2.*::\s*(NeoC|NonNeoC)_(TOST_Equivalence|PERMANOVA_PERMDISP|kBET|iLISI|Neighbour_Proportions)",
     "Secondary statistic computed from UMAP neighbourhoods, subsamples or permutations. "
     "The UMAP coordinates differ slightly between systems, which changes these values."),
    (r"Band_Formation.*::\s*Donor_Spread",
     "TVD_null_p95 is a permutation null threshold; it varies with the order of the permuted cells."),
    (r"15H_DepthThinning",
     "Values from random thinning of the counts."),
    (r"04_module_cuts\.csv|Part13b?_M4(Share|RandomSubset).*::\s*(K_Inputs|README)$",
     "The checksum column and its description were removed from the deposited scripts."),
    (r"Part13a_CP10K_Find14DpctMT_candidates\.csv",
     "Candidate %MT columns are found by scanning the data folder, and the reference folder holds more files."),
    (r"::\s*Formulas_Legend$",
     "Explanatory sheet removed from the deposited script; it holds no results."),
]
VOLATILE = re.compile(r"(/data|/win|C:/|C:\\|/home/|\d{4}-\d{2}-\d{2}|\d{8}[_-]\d{4,6}|"
                      r"\d{1,2}:\d{2}(:\d{2})?|seconds|mins?\b|elapsed)", re.I)

os.makedirs(OUT, exist_ok=True)
parts = sys.argv[1:] or sorted(d for d in os.listdir(NEW)
                               if os.path.isdir(os.path.join(NEW, d))
                               and d not in SKIP_TOP
                               and os.path.isdir(os.path.join(WIN, d)))
# Figure 17 summary outputs sit in the Manuscript folder, which is otherwise skipped.
EXTRA = ["Manuscript/Fig17"]
if not sys.argv[1:]:
    parts += [d for d in EXTRA if os.path.isdir(os.path.join(NEW, d)) and os.path.isdir(os.path.join(WIN, d))]
tag = "all" if not sys.argv[1:] else re.sub(r"[^A-Za-z0-9_.-]", "_", "_".join(parts))[:80]
LOG = os.path.join(OUT, "compare_%s_log.txt" % tag)
DET = os.path.join(OUT, "compare_%s_detail.csv" % tag)
_log = open(LOG, "w", encoding="utf-8")


def say(s=""):
    print(s)
    _log.write(s + "\n")
    _log.flush()


def num(x):
    if x is None:
        return None
    if isinstance(x, (int, float)) and not isinstance(x, bool):
        return float(x)
    s = str(x).strip()
    if s == "" or s.lower() in ("nan", "na", "none", "null"):
        return None
    try:
        return float(s)
    except ValueError:
        return None


def is_blank(x):
    return x is None or (isinstance(x, float) and math.isnan(x)) or str(x).strip() in ("", "nan", "NA", "None")


CELLS = []          # cell-level differences, written when COMPARE_CELLS=1
DUMP = os.environ.get("COMPARE_CELLS", "0") == "1"
DUMP_MAX = 500      # rows of cell differences kept per table


def align_rows(a, b):
    """Rows matched on the first column when it is a unique key in both tables
    (row order can differ between systems when values tie)."""
    if a.shape != b.shape or a.shape[0] < 3:
        return a, b, False
    ka, kb = a.iloc[1:, 0].astype(str), b.iloc[1:, 0].astype(str)
    if ka.is_unique and kb.is_unique and set(ka) == set(kb) and not ka.equals(kb):
        a2 = pd.concat([a.iloc[:1], a.iloc[1:].assign(_k=ka.values).sort_values("_k").drop(columns="_k")])
        b2 = pd.concat([b.iloc[:1], b.iloc[1:].assign(_k=kb.values).sort_values("_k").drop(columns="_k")])
        return a2.reset_index(drop=True), b2.reset_index(drop=True), True
    return a, b, False


def same_rows_any_order(a, b):
    """True when both tables hold the same rows (header kept first)."""
    if a.shape != b.shape or a.shape[0] < 2:
        return False
    ra = sorted(map(tuple, a.iloc[1:].astype(str).values))
    rb = sorted(map(tuple, b.iloc[1:].astype(str).values))
    return ra == rb and list(a.iloc[0].astype(str)) == list(b.iloc[0].astype(str))


def explain(where):
    txt = "%s/%s :: %s" % where
    return " ".join(note for pat, note in EXPLAIN if re.search(pat, txt))


def compare_grids(a, b, where=("", "", "")):
    """Cell-by-cell comparison of two raw grids (DataFrames, no header)."""
    r = dict(shape_new="%dx%d" % a.shape, shape_win="%dx%d" % b.shape,
             n_num=0, n_num_diff=0, n_sign=0, n_within=0, n_outside=0, max_abs=0.0, max_rel=0.0,
             n_text=0, n_text_diff=0, n_volatile_diff=0, n_blank_diff=0, example="",
             rows_aligned_on_key=False)
    if a.shape != b.shape:
        r["status"] = "DIFFERENT SIZE"
        return r
    a0, b0 = a, b
    a, b, r["rows_aligned_on_key"] = align_rows(a, b)
    av, bv = a.values, b.values
    head = [str(x) for x in av[0]] if av.shape[0] else []
    kept = [0]

    def keep(i, j, x, y, kind):
        if DUMP and kept[0] < DUMP_MAX:
            kept[0] += 1
            CELLS.append(dict(part=where[0], file=where[1], sheet=where[2], row=i + 1, col=j + 1,
                              column=head[j] if j < len(head) else "", key=str(av[i, 0]),
                              rerun=x, reference=y, kind=kind))
    for i in range(av.shape[0]):
        for j in range(av.shape[1]):
            x, y = av[i, j], bv[i, j]
            if i > 0 and j < len(head) and head[j] in PROVENANCE_COLS:
                if str(x).strip() != str(y).strip():
                    r["n_volatile_diff"] += 1
                continue
            bx, by = is_blank(x), is_blank(y)
            if bx and by:
                continue
            if bx != by:
                r["n_blank_diff"] += 1
                keep(i, j, x, y, "blank")
                if not r["example"]:
                    r["example"] = "r%d c%d: %r vs %r" % (i + 1, j + 1, x, y)
                continue
            nx, ny = num(x), num(y)
            if nx is not None and ny is not None:
                r["n_num"] += 1
                if math.isnan(nx) and math.isnan(ny):
                    continue
                if math.isinf(nx) or math.isinf(ny):
                    same = nx == ny
                    d = 0.0 if same else float("inf")
                    rel = d
                else:
                    d = abs(nx - ny)
                    rel = d / max(abs(nx), abs(ny)) if max(abs(nx), abs(ny)) > 0 else 0.0
                    same = d <= ABS_TOL or rel <= REL_TOL
                if not same:
                    s_ = abs(nx + ny)
                    if not math.isinf(nx) and (s_ <= ABS_TOL or s_ <= REL_TOL * max(abs(nx), abs(ny))):
                        r["n_sign"] += 1
                        keep(i, j, x, y, "sign")
                        continue
                    r["max_abs"] = max(r["max_abs"], d)
                    r["max_rel"] = max(r["max_rel"], rel)
                    r["n_num_diff"] += 1
                    if rel <= LARGE_TOL:
                        r["n_within"] += 1
                        keep(i, j, x, y, "within 0.1%")
                    else:
                        r["n_outside"] += 1
                        keep(i, j, x, y, "outside 0.1%")
                    if not r["example"]:
                        r["example"] = "r%d c%d: %r vs %r" % (i + 1, j + 1, x, y)
            else:
                r["n_text"] += 1
                if str(x).strip() != str(y).strip():
                    if VOLATILE.search(str(x)) or VOLATILE.search(str(y)):
                        r["n_volatile_diff"] += 1
                    else:
                        r["n_text_diff"] += 1
                        keep(i, j, x, y, "text")
                        if not r["example"]:
                            r["example"] = "r%d c%d: %r vs %r" % (i + 1, j + 1, str(x)[:80], str(y)[:80])
    if r["n_outside"] or r["n_blank_diff"]:
        r["status"] = "OUTSIDE 0.1%"
    elif r["n_within"]:
        r["status"] = "WITHIN 0.1%"
    elif r["n_text_diff"]:
        r["status"] = "TEXT ONLY"
    elif r["n_sign"]:
        r["status"] = "SAME UP TO SIGN"
    else:
        r["status"] = "SAME"
    if r["status"] != "SAME" and same_rows_any_order(a0, b0):
        r["status"] = "SAME VALUES, OTHER ORDER"
    return r


def read_tables(path):
    if path.endswith(".xlsx"):
        return pd.read_excel(path, sheet_name=None, header=None, dtype=object)
    return {"": pd.read_csv(path, header=None, dtype=str, keep_default_na=False,
                            na_values=[], low_memory=False)}


t0 = time.time()
say("compare_rerun.py   %s" % time.strftime("%Y-%m-%d %H:%M:%S"))
say("Rerun   : %s" % NEW)
say("Reference: %s" % WIN)
say("Tolerance: relative %.0e, absolute %.0e" % (REL_TOL, ABS_TOL))
say("Parts   : %s" % ", ".join(parts))
say("")
say("How to read this log")
say("  The reference results were produced on Windows; the rerun used the container on")
say("  Linux. Results are not expected to be identical to the last digit across operating")
say("  systems and linear algebra libraries, and statistics based on random subsamples,")
say("  permutations or nearest neighbours can change more. Each table is given one result:")
say("    SAME                      all numbers within 1e-6 (relative), text identical")
say("    SAME UP TO SIGN           only opposite signs of equal size (eigenvector signs)")
say("    SAME VALUES, OTHER ORDER  the same rows in a different order")
say("    TEXT ONLY                 numbers the same, some text differs")
say("    WITHIN 0.1%               numbers differ by more than 1e-6 but at most 0.1%")
say("    OUTSIDE 0.1%              some numbers differ by more than 0.1%")
say("  Tables only in the reference come from scripts or steps that are not part of the")
say("  deposit. Known causes are given as notes under the entries.")
say("")

rows = []
for part in parts:
    pn, pw = os.path.join(NEW, part), os.path.join(WIN, part)
    if not os.path.isdir(pw):
        say("== %s: not in the reference" % part)
        continue
    fn = {os.path.relpath(os.path.join(d, f), pn) for d, _, fs in os.walk(pn) for f in fs}
    fw = {os.path.relpath(os.path.join(d, f), pw) for d, _, fs in os.walk(pw) for f in fs}
    tab_n = sorted(f for f in fn if f.endswith(TABLE_EXT) and not os.path.basename(f).startswith("~$"))
    tab_w = {f for f in fw if f.endswith(TABLE_EXT)}
    other_n = {f for f in fn if not f.endswith(TABLE_EXT)}
    other_w = {f for f in fw if not f.endswith(TABLE_EXT)}
    cnt = {k: 0 for k in ORDER}
    for f in tab_n:
        if f not in tab_w:
            cnt["ONLY IN RERUN"] += 1
            rows.append(dict(part=part, file=f, sheet="", status="ONLY IN RERUN"))
            continue
        try:
            A, B = read_tables(os.path.join(pn, f)), read_tables(os.path.join(pw, f))
        except Exception as e:  # noqa: BLE001
            cnt["READ ERROR"] += 1
            rows.append(dict(part=part, file=f, sheet="", status="READ ERROR", example=str(e)[:200]))
            continue
        for sh in sorted(set(A) | set(B)):
            if sh not in A or sh not in B:
                st = "SHEET ONLY IN " + ("RERUN" if sh in A else "REFERENCE")
                rows.append(dict(part=part, file=f, sheet=sh, status=st))
                cnt[st] += 1
                continue
            r = compare_grids(A[sh].fillna(""), B[sh].fillna(""), (part, f, sh))
            r.update(part=part, file=f, sheet=sh)
            rows.append(r)
            cnt[r["status"]] += 1
    only_win = sorted(tab_w - set(tab_n))
    for f in only_win:
        rows.append(dict(part=part, file=f, sheet="", status="ONLY IN REFERENCE"))
        cnt["ONLY IN REFERENCE"] += 1
    say("== %s" % part)
    say("   tables (file or sheet): " + " | ".join("%s %d" % (k, v) for k, v in cnt.items() if v))
    say("   other files (figures, RDS, logs): rerun %d | reference %d | same name in both %d"
        % (len(other_n), len(other_w), len(other_n & other_w)))
    for r in rows:
        if r["part"] == part and r["status"] != "SAME":
            ex = r.get("example", "")
            extra = ""
            if "n_num" in r and r["status"] not in ("DIFFERENT SIZE", "SAME VALUES, OTHER ORDER"):
                bits = []
                if r["n_sign"]:
                    bits.append("%d of %d numbers differ in sign only" % (r["n_sign"], r["n_num"]))
                if r["n_num_diff"]:
                    bits.append("%d of %d numbers beyond 1e-6, %d beyond 0.1%%, max rel %.2e"
                                % (r["n_num_diff"], r["n_num"], r["n_outside"], r["max_rel"]))
                if r["n_blank_diff"]:
                    bits.append("%d cells empty in one table only" % r["n_blank_diff"])
                if r["n_text_diff"]:
                    bits.append("%d text cells differ" % r["n_text_diff"])
                extra = "".join(" | " + x for x in bits)
            elif r["status"] == "DIFFERENT SIZE":
                extra = " | rerun %s, reference %s" % (r["shape_new"], r["shape_win"])
            note = explain((part, r["file"], r.get("sheet", "")))
            if r["status"] == "SAME VALUES, OTHER ORDER":
                ex = ""
            say("   [%s] %s%s%s%s%s" % (r["status"], r["file"], (" :: " + r["sheet"]) if r.get("sheet") else "",
                                        extra, ("\n        e.g. " + ex) if ex else "",
                                        ("\n        note: " + note) if note else ""))
    say("")

pd.DataFrame(rows).to_csv(DET, index=False)
if DUMP:
    CEL = DET.replace("_detail.csv", "_cells.csv")
    pd.DataFrame(CELLS).to_csv(CEL, index=False)
tot = pd.Series([r["status"] for r in rows]).value_counts()
say("== Overall")
for k in ORDER:
    if k in tot:
        say("   %-26s %d" % (k, tot[k]))
say("")
say("Detail: %s" % DET)
if DUMP:
    say("Cells : %s" % CEL)
say("Done in %.1f min" % ((time.time() - t0) / 60))
