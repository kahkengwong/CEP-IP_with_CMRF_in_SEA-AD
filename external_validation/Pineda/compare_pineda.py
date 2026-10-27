#!/usr/bin/env python3
# =============================================================================
# compare_pineda.py
# Compares a container rerun of external validation 1 (Pineda et al. 2024) with the original
# results, file by file. The table comparison is that of compare_rexach.py and of
# verification/compare_rerun.py of the SEA-AD analysis (the same tolerances, results, and rules);
# the list of files differs, and two parts are added:
#   - the R02 cell calls, nucleus by nucleus and library by library (from the R11 nuclei table):
#     the called nuclei, the emptyDrops FDR, the scDblFinder score and class, and the kept nuclei;
#   - tables of more than 1 million cells are compared column by column with the same tolerances.
#
#   Tables (CSV, .csv.gz, and every sheet of each Excel workbook): numbers within a relative
#   tolerance, text exactly (cells holding a path, a date, a time, or a BLAS library are
#   counted separately).
#   R12 MapMyCells input (.h5ad): the nuclei, the genes, and every count must be identical.
#   R02 per-library .rds files are listed by presence; their nuclei are compared through R11.
#   R13 count matrix (.mtx.gz): the decompressed file must be identical.
#   Figures, RDS files, and logs are listed by presence only.
#
# Each table is given one of these results:
#   SAME                      all numbers within 1e-6 (relative), text identical
#   SAME UP TO SIGN           the only differences are numbers with the opposite sign and the
#                             same size
#   SAME VALUES, OTHER ORDER  the same rows in a different order
#   TEXT ONLY                 numbers the same, some text differs
#   WITHIN 0.1%               some numbers differ by more than 1e-6 but all by at most 0.1%
#   OUTSIDE 0.1%              some numbers differ by more than 0.1%, or a cell is empty in one
#                             table only
#   PROVENANCE ONLY           a Session sheet (versions, paths, times, BLAS, workers) differs
#
# Usage (inside the container, after pineda-run; pineda-run --compare calls it):
#   python /opt/pineda/compare_pineda.py all     # R02 to R17 and Table S20
#   python /opt/pineda/compare_pineda.py r       # R14 to R17 and Table S20
# The rerun is PINEDA_ROOT (/data), the original External_Dataset_1 folder is REF_ROOT (/ref),
# and the original Manuscript/Supp_Tables folder is REF_SUPP (/ref_supp, optional).
#
# Outputs (/data/compare/): compare_pineda_<mode>_log.txt (summary),
#   compare_pineda_<mode>_detail.csv (one row per file and sheet), and, for mode all,
#   compare_pineda_all_R02_libraries.csv (the cell calls of each library). With COMPARE_CELLS=1,
#   also compare_pineda_<mode>_cells.csv (each differing cell, up to 500 per table).
# =============================================================================
import glob
import gzip
import math
import os
import re
import sys
import time

import warnings

import pandas as pd

warnings.filterwarnings("ignore", message="This pattern is interpreted as a regular expression")

NEW = os.environ.get("PINEDA_ROOT", "/data")
REF = os.environ.get("REF_ROOT", "/ref")
REF_SUPP = os.environ.get("REF_SUPP", "/ref_supp")
OUT = os.path.join(NEW, "compare")
REL_TOL = 1e-6      # numbers agree when |a - b| <= REL_TOL * max(|a|, |b|) or <= ABS_TOL
ABS_TOL = 1e-9
LARGE_TOL = 1e-3    # numbers within 0.1% (relative) are reported as WITHIN 0.1%
TABLE_EXT = (".csv", ".csv.gz", ".xlsx")
# Columns that record the run, not results: R02 "minutes" is the time each library took.
PROVENANCE_COLS = {"config_version", "minutes"}
BIG_CELLS = 1_000_000   # tables with more cells are compared column by column
ORDER = ["SAME", "SAME UP TO SIGN", "SAME VALUES, OTHER ORDER", "TEXT ONLY",
         "WITHIN 0.1%", "OUTSIDE 0.1%", "DIFFERENT SIZE", "SHEET ONLY IN RERUN",
         "SHEET ONLY IN REFERENCE", "ONLY IN RERUN", "ONLY IN REFERENCE", "READ ERROR",
         "IDENTICAL", "DIFFERENT", "PROVENANCE ONLY"]
# Sheets that record the run, not results: differences confined to their text are PROVENANCE ONLY.
PROVENANCE_SHEET = re.compile(r"^(\d+_)?Session$")
# What is compared, per mode: the folders, each with the same name in the rerun and the reference
R_PART = ["R14_core", "R14_PVM_3", "R14_HKG8", "R14_HKG_control", "R15_contrasts", "R17_posthoc",
          "R16_Figure_Pineda"]
FOLDERS = {
    "all": ["R02_cell_calling", "R11_nuclei", "R12_mapmycells_input", "R13_count_checks",
            "R13b_mapping_comparison_repeat", "R13b_mapping_comparison_QC", "R13c_noAuthorQC_arm"] + R_PART,
    "r": R_PART,
}
# Known causes, matched on "part/file :: sheet". Printed under the matching entries.
EXPLAIN = [
    (r"R17_posthoc/.*::\s*(01_Summary|03_Permutation|06_Bootstrap)$",
     "Monte Carlo statistic (Freedman-Lane permutations or the donor bootstrap) with fixed seeds."),
    (r"::\s*(\d+_)?Session$",
     "Run record: start and end times, paths, and the BLAS and LAPACK libraries of the system."),
    (r"::\s*(\d+_)?Settings$",
     "Settings record: holds the paths of the input files."),
    (r"R02_cell_calling/.*::\s*Libraries$",
     "The column minutes (time per library) is not compared."),
    (r"(R11_nuclei/R11_Pineda_nuclei|R12_mapmycells_input/R12_Pineda_nuclei)\.csv\.gz",
     "The differences are in the column scDblFinder_score, a continuous score computed after a PCA, "
     "which depends on the linear algebra library. The R02 cell calls above show whether any "
     "doublet class changed."),
]
VOLATILE = re.compile(r"(/data|/ref|C:/|C:\\|/home/|/mnt/|/usr/|/opt/|\.so(\.\d+)*\b|blas|lapack|"
                      r"\d{4}-\d{2}-\d{2}|\d{8}[_-]\d{4,6}|\d{1,2}:\d{2}(:\d{2})?|seconds|mins?\b|elapsed)",
                      re.I)

mode = sys.argv[1] if len(sys.argv) > 1 else "all"
if mode not in FOLDERS:
    sys.exit("Usage: compare_pineda.py all|r")
os.makedirs(OUT, exist_ok=True)
LOG = os.path.join(OUT, "compare_pineda_%s_log.txt" % mode)
DET = os.path.join(OUT, "compare_pineda_%s_detail.csv" % mode)
R02_LIB = os.path.join(OUT, "compare_pineda_all_R02_libraries.csv")
_log = open(LOG, "w", encoding="utf-8")


def say(s=""):
    print(s)
    _log.write(s + "\n")
    _log.flush()


# ---- table comparison: as verification/compare_rerun.py ------------------------------------------
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


CELLS = []
DUMP = os.environ.get("COMPARE_CELLS", "0") == "1"
DUMP_MAX = 500


def align_rows(a, b):
    """Rows matched on the first column when it is a unique key in both tables."""
    if a.shape != b.shape or a.shape[0] < 3:
        return a, b, False
    ka, kb = a.iloc[1:, 0].astype(str), b.iloc[1:, 0].astype(str)
    if ka.is_unique and kb.is_unique and set(ka) == set(kb) and not ka.equals(kb):
        a2 = pd.concat([a.iloc[:1], a.iloc[1:].assign(_k=ka.values).sort_values("_k").drop(columns="_k")])
        b2 = pd.concat([b.iloc[:1], b.iloc[1:].assign(_k=kb.values).sort_values("_k").drop(columns="_k")])
        return a2.reset_index(drop=True), b2.reset_index(drop=True), True
    return a, b, False


def same_rows_any_order(a, b):
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
    r = dict(shape_new="%dx%d" % a.shape, shape_ref="%dx%d" % b.shape,
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


def compare_big(a, b):
    """Column-by-column comparison of two large raw grids, with the rules of compare_grids()."""
    import numpy as np
    r = dict(shape_new="%dx%d" % a.shape, shape_ref="%dx%d" % b.shape,
             n_num=0, n_num_diff=0, n_sign=0, n_within=0, n_outside=0, max_abs=0.0, max_rel=0.0,
             n_text=0, n_text_diff=0, n_volatile_diff=0, n_blank_diff=0, example="",
             rows_aligned_on_key=False)
    if a.shape != b.shape:
        r["status"] = "DIFFERENT SIZE"
        return r
    head = [str(x) for x in a.iloc[0]]
    if head != [str(x) for x in b.iloc[0]]:
        r["n_text_diff"] += 1
        r["example"] = "header row differs"
    A, B = a.iloc[1:].reset_index(drop=True), b.iloc[1:].reset_index(drop=True)
    ka, kb = A.iloc[:, 0].astype(str), B.iloc[:, 0].astype(str)
    if not ka.equals(kb) and ka.is_unique and kb.is_unique and set(ka) == set(kb):
        B = B.set_index(kb).loc[ka.values].reset_index(drop=True)
        r["rows_aligned_on_key"] = True
    blanks = ["", "nan", "NA", "None", "NaN"]
    for j in range(A.shape[1]):
        if head[j] in PROVENANCE_COLS:
            continue
        x, y = A.iloc[:, j].astype(str).str.strip(), B.iloc[:, j].astype(str).str.strip()
        bx, by = x.isin(blanks), y.isin(blanks)
        r["n_blank_diff"] += int((bx != by).sum())
        both = ~bx & ~by
        nx = pd.to_numeric(x.where(both), errors="coerce")
        ny = pd.to_numeric(y.where(both), errors="coerce")
        isnum = both & nx.notna() & ny.notna()
        r["n_num"] += int(isnum.sum())
        d = (nx - ny).abs()
        mx = np.maximum(nx.abs(), ny.abs())
        rel = (d / mx.where(mx > 0)).fillna(0.0)
        close = (d <= ABS_TOL) | (rel <= REL_TOL) | (x == y)
        diff = isnum & ~close
        s_ = (nx + ny).abs()
        sign = diff & ((s_ <= ABS_TOL) | (s_ <= REL_TOL * mx))
        real = diff & ~sign
        within = real & (rel <= LARGE_TOL)
        r["n_sign"] += int(sign.sum())
        r["n_num_diff"] += int(real.sum())
        r["n_within"] += int(within.sum())
        r["n_outside"] += int((real & ~within).sum())
        if real.any():
            r["max_abs"] = max(r["max_abs"], float(d[real].max()))
            r["max_rel"] = max(r["max_rel"], float(rel[real].max()))
        txt = both & ~isnum
        r["n_text"] += int(txt.sum())
        td = txt & (x != y)
        if td.any():
            vol = td & (x.str.contains(VOLATILE) | y.str.contains(VOLATILE))
            r["n_volatile_diff"] += int(vol.sum())
            r["n_text_diff"] += int((td & ~vol).sum())
        bad = (bx != by) | real | (td if td.any() else td)
        if not r["example"] and bad.any():
            i = int(np.flatnonzero(bad.values)[0])
            r["example"] = "column %s, row %d: %r vs %r" % (head[j], i + 2, x.iloc[i][:80], y.iloc[i][:80])
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
    return r


def r02_calls():
    """R02 cell calls of each library, nucleus by nucleus, from the R11 nuclei table."""
    import numpy as np
    f = os.path.join("R11_nuclei", "R11_Pineda_nuclei.csv.gz")
    pn, pw = os.path.join(NEW, f), os.path.join(REF, f)
    say("== R02 cell calls (from %s)" % f)
    if not (os.path.exists(pn) and os.path.exists(pw)):
        say("   not compared: %s missing in the %s" % (f, "rerun" if not os.path.exists(pn) else "reference"))
        say("")
        return dict(part="R02 cell calls", file=f, sheet="", status="ONLY IN " + ("REFERENCE" if os.path.exists(pw) else "RERUN"))
    cols = ["cell_id", "library_id", "emptydrops_FDR", "pass_gene_floor", "scDblFinder_score",
            "scDblFinder_class", "keep"]
    rd = lambda p: pd.read_csv(p, usecols=cols, dtype={"cell_id": str, "library_id": str, "scDblFinder_class": str})
    a, b = rd(pn), rd(pw)
    for d in (a, b):
        for c in ("pass_gene_floor", "keep"):
            d[c] = d[c].astype(str).str.upper().isin(["TRUE", "T", "1"])
    m = a.merge(b, on="cell_id", how="outer", suffixes=("_new", "_ref"), indicator=True)
    m["library"] = m["library_id_new"].fillna(m["library_id_ref"])
    inb = m["_merge"] == "both"
    m["fdr_rel"] = ((m["emptydrops_FDR_new"] - m["emptydrops_FDR_ref"]).abs() /
                    np.maximum(m["emptydrops_FDR_new"].abs(), m["emptydrops_FDR_ref"].abs()).where(lambda v: v > 0)).fillna(0.0)
    m["score_abs"] = (m["scDblFinder_score_new"] - m["scDblFinder_score_ref"]).abs()
    out = []
    for lib, g in m.groupby("library", sort=True):
        gb = g[g["_merge"] == "both"]
        fl = gb[gb["pass_gene_floor_new"] & gb["pass_gene_floor_ref"]]
        out.append(dict(
            library=lib, called_rerun=int((g["_merge"] != "right_only").sum()),
            called_reference=int((g["_merge"] != "left_only").sum()), called_both=len(gb),
            gene_floor_differs=int((gb["pass_gene_floor_new"] != gb["pass_gene_floor_ref"]).sum()),
            doublet_class_differs=int((fl["scDblFinder_class_new"] != fl["scDblFinder_class_ref"]).sum()),
            kept_rerun=int(g["keep_new"].fillna(False).astype(bool).sum()),
            kept_reference=int(g["keep_ref"].fillna(False).astype(bool).sum()),
            kept_both=int((gb["keep_new"] & gb["keep_ref"]).sum()),
            max_rel_diff_FDR=float(gb["fdr_rel"].max()) if len(gb) else float("nan"),
            max_abs_diff_score=float(gb["score_abs"].max()) if len(fl) else float("nan")))
    t = pd.DataFrame(out)
    t["same_calls"] = ((t["called_rerun"] == t["called_both"]) & (t["called_reference"] == t["called_both"]) &
                       (t["gene_floor_differs"] == 0) & (t["doublet_class_differs"] == 0) &
                       (t["kept_rerun"] == t["kept_both"]) & (t["kept_reference"] == t["kept_both"]))
    t.to_csv(R02_LIB, index=False)
    n_same = int(t["same_calls"].sum())
    say("   libraries with the same calls: %d of %d" % (n_same, len(t)))
    say("   called nuclei: rerun %d | reference %d | in both %d" %
        (t["called_rerun"].sum(), t["called_reference"].sum(), t["called_both"].sum()))
    say("   kept nuclei:   rerun %d | reference %d | in both %d" %
        (t["kept_rerun"].sum(), t["kept_reference"].sum(), t["kept_both"].sum()))
    say("   doublet class differs for %d nuclei | largest emptyDrops FDR difference %.2e (relative) | "
        "largest scDblFinder score difference %.2e" %
        (t["doublet_class_differs"].sum(), t["max_rel_diff_FDR"].max(), t["max_abs_diff_score"].max()))
    for _, x in t[~t["same_calls"]].iterrows():
        say("   [DIFFERENT] %s: called %d / %d (both %d) | doublet class differs %d | kept %d / %d (both %d)" %
            (x["library"], x["called_rerun"], x["called_reference"], x["called_both"],
             x["doublet_class_differs"], x["kept_rerun"], x["kept_reference"], x["kept_both"]))
    say("   per library: %s" % R02_LIB)
    say("")
    return dict(part="R02 cell calls", file=f, sheet="", status="IDENTICAL" if n_same == len(t) else "DIFFERENT",
                example="%d of %d libraries with the same calls" % (n_same, len(t)))


def read_tables(path):
    if path.endswith(".xlsx"):
        return pd.read_excel(path, sheet_name=None, header=None, dtype=object)
    return {"": pd.read_csv(path, header=None, dtype=str, keep_default_na=False,
                            na_values=[], low_memory=False)}


# ---- exact comparisons: the R12 h5ad parts and the R13 count matrix ---------------------------------
def same_h5ad(p, q):
    import anndata as ad
    import numpy as np
    a, b = ad.read_h5ad(p), ad.read_h5ad(q)
    if a.shape != b.shape:
        return False, "shape %s vs %s" % (a.shape, b.shape)
    if not (list(a.obs_names) == list(b.obs_names) and list(a.var_names) == list(b.var_names)):
        return False, "nuclei or genes differ"
    xa, xb = a.X, b.X
    diff = (xa != xb)
    n = int(diff.nnz) if hasattr(diff, "nnz") else int(np.sum(diff))
    return n == 0, "%d nuclei x %d genes, %d counts differ" % (a.shape[0], a.shape[1], n)


def same_gz(p, q):
    with gzip.open(p, "rb") as f, gzip.open(q, "rb") as g:
        while True:
            x, y = f.read(1 << 20), g.read(1 << 20)
            if x != y:
                return False, "decompressed content differs"
            if not x:
                return True, "decompressed content identical"


# ---- run -----------------------------------------------------------------------------------------
t0 = time.time()
say("compare_pineda.py   %s | mode %s" % (time.strftime("%Y-%m-%d %H:%M:%S"), mode))
say("Rerun    : %s" % NEW)
say("Reference: %s (Table S20: %s)" % (REF, REF_SUPP if os.path.isdir(REF_SUPP) else "not mounted"))
say("Tolerance: relative %.0e, absolute %.0e" % (REL_TOL, ABS_TOL))
say("")
say("How to read this log")
say("  The reference results were produced on Windows (R) and in WSL (Python); the rerun used the")
say("  container on Linux. Numbers can differ in the last digits across operating systems and")
say("  linear algebra libraries. Each table is given one result:")
say("    SAME                      all numbers within 1e-6 (relative), text identical")
say("    SAME UP TO SIGN           only opposite signs of equal size")
say("    SAME VALUES, OTHER ORDER  the same rows in a different order")
say("    TEXT ONLY                 numbers the same, some text differs")
say("    WITHIN 0.1%               numbers differ by more than 1e-6 but at most 0.1%")
say("    OUTSIDE 0.1%              some numbers differ by more than 0.1%")
say("    PROVENANCE ONLY           a Session sheet (versions, paths, times, BLAS, workers) differs")
say("  The R12 MapMyCells input files and the R13 count matrix are compared exactly")
say("  (IDENTICAL or DIFFERENT). Known causes are given as notes under the entries.")
say("")

targets = [(d, os.path.join(NEW, d), os.path.join(REF, d), None) for d in FOLDERS[mode]]
targets.append(("Table S20", os.path.join(NEW, "Supp_Tables"), REF_SUPP, ["Supplementary_Table_S20_*.xlsx"]))


def list_files(base, pats):
    if not os.path.isdir(base):
        return set()
    if pats:
        return {os.path.basename(p) for pat in pats for p in glob.glob(os.path.join(base, pat))
                if "_vOri" not in p}
    return {os.path.relpath(os.path.join(d, f), base) for d, _, fs in os.walk(base) for f in fs}


rows = []
if mode == "all":
    rows.append(r02_calls())
for label, pn, pw, pats in targets:
    if not os.path.isdir(pw):
        say("== %s: reference folder %s not found" % (label, pw))
        say("")
        continue
    fn, fw = list_files(pn, pats), list_files(pw, pats)
    tab_n = sorted(f for f in fn if f.endswith(TABLE_EXT) and not os.path.basename(f).startswith("~$"))
    tab_w = {f for f in fw if f.endswith(TABLE_EXT)}
    exact_n = sorted(f for f in fn if f.endswith((".h5ad", ".mtx.gz")))
    other_n = {f for f in fn if not f.endswith(TABLE_EXT) and f not in exact_n}
    other_w = {f for f in fw if not f.endswith(TABLE_EXT) and not f.endswith((".h5ad", ".mtx.gz"))}
    cnt = {k: 0 for k in ORDER}
    for f in tab_n:
        if f not in tab_w:
            cnt["ONLY IN RERUN"] += 1
            rows.append(dict(part=label, file=f, sheet="", status="ONLY IN RERUN"))
            continue
        try:
            A, B = read_tables(os.path.join(pn, f)), read_tables(os.path.join(pw, f))
        except Exception as e:  # noqa: BLE001
            cnt["READ ERROR"] += 1
            rows.append(dict(part=label, file=f, sheet="", status="READ ERROR", example=str(e)[:200]))
            continue
        for sh in sorted(set(A) | set(B)):
            if sh not in A or sh not in B:
                st = "SHEET ONLY IN " + ("RERUN" if sh in A else "REFERENCE")
                rows.append(dict(part=label, file=f, sheet=sh, status=st))
                cnt[st] += 1
                continue
            ga, gb = A[sh].fillna(""), B[sh].fillna("")
            r = (compare_big(ga, gb) if ga.size > BIG_CELLS else compare_grids(ga, gb, (label, f, sh)))
            # A Session sheet records the run (versions, paths, times, BLAS, workers, and the lines of
            # sessionInfo()), so any difference in it is PROVENANCE ONLY. The result it would have had
            # is kept in status_before.
            if PROVENANCE_SHEET.search(str(sh)) and r["status"] != "SAME":
                r["status_before"] = r["status"]
                r["status"] = "PROVENANCE ONLY"
            r.update(part=label, file=f, sheet=sh)
            rows.append(r)
            cnt[r["status"]] += 1
    for f in sorted(tab_w - set(tab_n)):
        rows.append(dict(part=label, file=f, sheet="", status="ONLY IN REFERENCE"))
        cnt["ONLY IN REFERENCE"] += 1
    for f in exact_n:
        q = os.path.join(pw, f)
        if not os.path.exists(q):
            rows.append(dict(part=label, file=f, sheet="", status="ONLY IN RERUN"))
            cnt["ONLY IN RERUN"] += 1
            continue
        try:
            ok, what = (same_h5ad if f.endswith(".h5ad") else same_gz)(os.path.join(pn, f), q)
        except Exception as e:  # noqa: BLE001
            ok, what = False, "read error: %s" % str(e)[:150]
        st = "IDENTICAL" if ok else "DIFFERENT"
        rows.append(dict(part=label, file=f, sheet="", status=st, example=what))
        cnt[st] += 1
    say("== %s" % label)
    say("   files compared: " + " | ".join("%s %d" % (k, v) for k, v in cnt.items() if v))
    say("   other files (figures, RDS, logs): rerun %d | reference %d | same name in both %d"
        % (len(other_n), len(other_w), len(other_n & other_w)))
    for r in rows:
        if r["part"] != label or r["status"] in ("SAME", "IDENTICAL"):
            continue
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
            extra = " | rerun %s, reference %s" % (r["shape_new"], r["shape_ref"])
        note = explain((label, r["file"], r.get("sheet", "")))
        if r["status"] == "SAME VALUES, OTHER ORDER":
            ex = ""
        say("   [%s] %s%s%s%s%s" % (r["status"], r["file"], (" :: " + r["sheet"]) if r.get("sheet") else "",
                                    extra, ("\n        e.g. " + ex) if ex else "",
                                    ("\n        note: " + note) if note else ""))
    say("")

pd.DataFrame(rows).to_csv(DET, index=False)
if DUMP:
    pd.DataFrame(CELLS).to_csv(DET.replace("_detail.csv", "_cells.csv"), index=False)
tot = pd.Series([r["status"] for r in rows]).value_counts()
say("== Overall")
for k in ORDER:
    if k in tot:
        say("   %-26s %d" % (k, tot[k]))
agree = sum(tot.get(k, 0) for k in ("SAME", "SAME UP TO SIGN", "SAME VALUES, OTHER ORDER", "TEXT ONLY",
                                     "WITHIN 0.1%", "IDENTICAL", "PROVENANCE ONLY"))
say("   agree (SAME to WITHIN 0.1%%, IDENTICAL, or PROVENANCE ONLY): %d of %d" % (agree, len(rows)))
say("")
say("Detail: %s" % DET)
say("Done in %.1f min" % ((time.time() - t0) / 60))
