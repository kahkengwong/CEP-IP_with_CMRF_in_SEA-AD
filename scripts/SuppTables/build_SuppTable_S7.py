#!/usr/bin/env python3
"""
build_SuppTable_S7.py

Builds SuppTable_S7_Fig2_SKCS.xlsx, the workbook behind Supplementary Table S7
(Figures 1H and 2 to 4, Supplementary Figures S2 to S10): one row per
GOI-positive cell of the four arms, the GAM fits, the S_i floor diagnostic of
Part 9b and the pitcher-plot summary. finalize_supp_tables.py turns this
workbook into Supplementary Table S7; it keeps the data sheets A_ to H_ and
replaces the Contents sheet written here.

Usage:
    python3 build_SuppTable_S7.py
    python3 build_SuppTable_S7.py --root /data
--root defaults to $SEAAD_ROOT, else the current folder.

Inputs, read from the root (nothing is modified):
  Part9/Part9_<arm>/Part9_CP10K_<tag>_<ds>_GAM_Results_1.xlsx
      (Part9_CP10K_QC_GAM_code.r) Cell_Count_Audit, GAM_Performance_DFG,
      GAM_Performance_HKG, Observed_Predicted (population-level prediction of
      every cell, used to label TREP among all cells)
  Part10/Part10_<arm>/Part10_<arm>_<ds>/Part10_CP10K_<tag>_<stem>_<ds>_GAM_Results.xlsx
      (Part10_CP10K_SKCS_code.r) GAM_Performance_SKCS, DE_Comparison
  Part10/Part10_<arm>/Part10_<arm>_<ds>/Part10_CP10K_<tag>_<stem>_<ds>_Pooled_SKCS_scores.csv
      (Part10_CP10K_SKCS_code.r) every GOI-positive cell: expression of the
      GOI, the composite and the 17 module and housekeeping genes, S_i, K_i,
      SKCS-core, SKCS-core TREP
  Part10/Part10_<arm>/Part10_CP10K_<tag>_GAM_Log.txt
      (Part10_CP10K_SKCS_code.r) pooled r_s, tau_b and SKCS-core counts,
      including the housekeeping genes used as GOI in the HKG8 arm
  Part10b_SKCS/Part10b_SKCS_<arm>/Part10b_CP10K_<tag>_SKCS_<ds>_EP_TREP_values.csv
      (Part10b_CP10K_SKCS_TREP_code.r) EP, TREP label and IP group of every
      SKCS-core cell
  R_export_CP10K/<export arm>/<neocortex|nonneocortex>/gam_data/*_selected.csv
      (Part8_Data-export-CP10K_MT-CO3_8DFG.py) donor of every cell
  Part9b/Part9b_Mic13_PVM3/Part9b_CP10K_M4_SiKi_threshold_comparison.xlsx
  Part9b/Part9b_Mic13_PVM3/Part9b_CP10K_M4_SiKi_Log.txt
      (Part9b_CP10K_SiKi_code.r)
<arm> is Mic13_PVM3, PVM_3, PVM_1 or HKG8; <tag> is M4, PVM_3, PVM_1 or HKG8;
<ds> is NeoC or NonNeoC.

Matched from Part 9 by cell: fit_all_cells, the population-level prediction of
the all-cells GAM (Observed_Predicted, Predicted_pop_value).

Recomputed here, because no part saves them per cell (every recomputed value
is printed next to the value in the logs or the figures, and the script stops
if a check fails):
  TREP among all GOI-positive cells: EP and the top round(n x DE) rule of
      Part 10, applied to the Part 9 population-level predictions. Part 9 has
      no barcodes, so its cells are matched to the score CSV by their GOI and
      composite values (tolerance 1e-4).
  S_i and K_i with YWHAZ and RAB7A as GOI in HKG8: the Part 10 formula
      (compute_skcs) on the HKG8 expression columns, as in Part 10 section D5b.
  Copula coordinates u, v, u_SK, v_SK: rank / (n + 1), ties averaged, within
      the arm and region set.

The workbooks are read with a small XML reader rather than openpyxl, because
the Part 9b workbook carries a drawing reference without the drawing and
openpyxl refuses to open it. The output is written in openpyxl write-only
mode, because the per-cell sheet holds about 37,000 rows.

Output:
  <root>/Manuscript/Supp_Tables/build/SuppTable_S7_Fig2_SKCS.xlsx
      sheets Contents, A_cells, B_HKG_pitcher_cells, C_GAM, D_Si_mixture,
      E_Si_boundaries, F_Si_posterior, G_Si_rank_tests, H_pitcher_summary
  <root>/Manuscript/Supp_Tables/build/SuppTable_S7_Fig2_SKCS_build_log.txt
      the console output, ending with VERDICT: PASS or VERDICT: FAIL
The workbook is written only if every check passes.

Formatting: Arial 10, bold header row, frozen top pane, no fills and no
borders.
"""
import argparse
import atexit
import bisect
import glob
import io
import os
import re
import sys
import zipfile
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd
from openpyxl import Workbook
from openpyxl.cell import WriteOnlyCell
from openpyxl.styles import Font, Alignment
from openpyxl.utils import get_column_letter

try:
    from scipy.stats import spearmanr, kendalltau
    HAVE_SCIPY = True
except ImportError:
    HAVE_SCIPY = False

# ---------------------------------------------------------------- settings
_ap = argparse.ArgumentParser(description="Build SuppTable_S7_Fig2_SKCS.xlsx")
_ap.add_argument("--root", default=os.environ.get("SEAAD_ROOT") or os.getcwd(),
                 help="SEA-AD folder (default: $SEAAD_ROOT, else the current folder)")
ROOT = os.path.abspath(_ap.parse_args().root)
if not os.path.isdir(ROOT):
    raise SystemExit("SEA-AD folder not found: %s" % ROOT)
OUT_DIR = os.path.join(ROOT, "Manuscript", "Supp_Tables", "build")
OUT = os.path.join(OUT_DIR, "SuppTable_S7_Fig2_SKCS.xlsx")
LOG = os.path.join(OUT_DIR, "SuppTable_S7_Fig2_SKCS_build_log.txt")


# ---------------------------------------------------------------- build log
class _Tee(io.TextIOBase):
    """Writes to the console and keeps a copy for the build log."""

    def __init__(self, stream, buf):
        self.stream, self.buf = stream, buf

    def write(self, s):
        self.stream.write(s)
        self.buf.write(s)
        return len(s)

    def flush(self):
        self.stream.flush()


_LOG_BUF = io.StringIO()
sys.stdout = _Tee(sys.stdout, _LOG_BUF)
sys.stderr = _Tee(sys.stderr, _LOG_BUF)
WRITTEN = False


def _write_build_log():
    """Saves the console output beside the workbook, with the verdict last."""
    sys.stdout.write("VERDICT: %s\n" % ("PASS" if WRITTEN else "FAIL"))
    os.makedirs(OUT_DIR, exist_ok=True)
    with open(LOG, "w", encoding="utf-8") as fh:
        fh.write(_LOG_BUF.getvalue())


atexit.register(_write_build_log)
print("build_SuppTable_S7.py")
print("Root: %s" % ROOT)
print("Output: %s\n" % OUT)

SETS = [
    dict(set="M4", p9="Part9_Mic13_PVM3", p10="Part10_Mic13_PVM3", tag="M4",
         stem="Module1_MTCO3_cDFG", goi="MT-CO3", comp="cDFG"),
    dict(set="PVM_3", p9="Part9_PVM_3", p10="Part10_PVM_3", tag="PVM_3",
         stem="Module1_MTCO3_cDFG", goi="MT-CO3", comp="cDFG"),
    dict(set="PVM_1", p9="Part9_PVM_1", p10="Part10_PVM_1", tag="PVM_1",
         stem="Module1_MTCO3_cDFG", goi="MT-CO3", comp="cDFG"),
    dict(set="HKG8", p9="Part9_HKG8", p10="Part10_HKG8", tag="HKG8",
         stem="Module1_ACTB_cHKG8", goi="ACTB", comp="cHKG8"),
]
DATASETS = ["NeoC", "NonNeoC"]
REGIONS = {"NeoC": ["DFC", "MTG", "AnG", "STG", "V1C", "ITG", "Pooled"],
           "NonNeoC": ["HIP", "MEC", "LEC", "FI", "Pooled"]}

P9B_DIR = os.path.join(ROOT, "Part9b", "Part9b_Mic13_PVM3")
P9B_XLSX = os.path.join(P9B_DIR, "Part9b_CP10K_M4_SiKi_threshold_comparison.xlsx")
P9B_LOG = os.path.join(P9B_DIR, "Part9b_CP10K_M4_SiKi_Log.txt")

# Reference values for the checks at the end. n, SKCS-core, SKCS-core TREP
# and SKCS-core donors are from the Part 10 guardrail and GAM logs of
# 2026-08-17. trep_all is read from the Figure 2 and Supplementary Figure 5
# annotations, because no part writes it to a file.
REF = {
    ("M4", "NeoC"): dict(n=2130, core=1419, trep_core=1308, don_core=74, trep_all=1565),
    ("M4", "NonNeoC"): dict(n=1987, core=1303, trep_core=1202, don_core=79, trep_all=1499),
    ("PVM_3", "NeoC"): dict(n=13331, core=8407, trep_core=7617, don_core=84, trep_all=9072),
    ("PVM_3", "NonNeoC"): dict(n=11132, core=6941, trep_core=6317, don_core=83, trep_all=7904),
    ("PVM_1", "NeoC"): dict(n=3218, core=1862, trep_core=1708, don_core=82, trep_all=2199),
    ("PVM_1", "NonNeoC"): dict(n=2658, core=1677, trep_core=1547, don_core=82, trep_all=1990),
    ("HKG8", "NeoC"): dict(n=1432, core=250, trep_core=218, don_core=55, trep_all=160),
    ("HKG8", "NonNeoC"): dict(n=1079, core=199, trep_core=159, don_core=56, trep_all=208),
}
# Part 9 pooled DE (per cent) for M4, from the Part 10 log baseline lines.
REF_DE_ALL_M4 = {"NeoC": 67.37, "NonNeoC": 65.2}


# Part 8 export folder of each arm (donor per barcode) and the region-set
# folder names used by the export.
EXPORT_DIR = {"M4": "Mic13_PVM3", "PVM_3": "PVM3_all", "PVM_1": "PVM1_all",
              "HKG8": "HKG8_matched"}
EXPORT_DS = {"NeoC": "neocortex", "NonNeoC": "nonneocortex"}

DFG_GENES = ["MT-ATP6", "MT-CO2", "MT-ND3", "MT-ND4", "MT-CYB", "MT-ND1",
             "MT-ND2", "MT-CO1"]
HKG_GENES = ["ACTB", "EEF1A1", "RPL13A", "B2M", "YWHAZ", "FTL", "PSAP",
             "RAB7A", "HNRNPA2B1"]
ALT_GOIS = ["YWHAZ", "RAB7A"]          # EEF1A1 is computed by Part 10 but not shown
MATCH_TOL = 1e-4                       # Part 9 to score CSV, log1p CP10K units

# Group counts of the SKCS-core (Figure 3 and Supplementary Figure S7):
# pre-IP TREP, pre-IP non-TREP, post-IP TREP, post-IP non-TREP.
REF_GROUPS = {
    ("M4", "NeoC"): (640, 31, 668, 80), ("M4", "NonNeoC"): (582, 44, 620, 57),
    ("PVM_3", "NeoC"): (3781, 101, 3836, 689), ("PVM_3", "NonNeoC"): (3103, 116, 3214, 508),
    ("PVM_1", "NeoC"): (787, 1, 921, 153), ("PVM_1", "NonNeoC"): (743, 6, 804, 124),
    ("HKG8", "NeoC"): (71, 0, 147, 32), ("HKG8", "NonNeoC"): (43, 0, 116, 40),
}
GROUPS = ["Pre-IP TREP", "Pre-IP non-TREP", "Post-IP TREP", "Post-IP non-TREP"]

# Worked cells printed in Figure 4 and Supplementary Figures 8 to 10, used to
# check the copula columns. Expression space: (arm, region set, GOI, rank,
# composite, rank). Score space: (arm, region set, K_i, rank S_i, rank K_i).
EXPR_CELLS = {
    "A": ("M4", "NeoC", 1.8209, 725, 0.9344, 459),
    "B": ("HKG8", "NeoC", 2.7345, 1346, 1.5275, 1032),
    "E": ("PVM_3", "NeoC", 3.2135, 11232, 2.4900, 10686),
    "F": ("PVM_1", "NeoC", 2.7748, 2499, 1.5316, 1834),
    "I": ("M4", "NonNeoC", 2.2005, 782, 1.1033, 419),
    "J": ("PVM_3", "NonNeoC", 2.9878, 7922, 2.6856, 8766),
    "K": ("PVM_1", "NonNeoC", 3.2983, 2061, 1.8387, 1460),
    "L": ("HKG8", "NonNeoC", 0.6701, 45, 0.7885, 211),
}
SCORE_CELLS = {
    "C": ("M4", "NeoC", 0.9079, 1705, 1939),
    "D": ("HKG8", "NeoC", 0.3627, 1146, 1022),
    "G": ("PVM_3", "NeoC", 0.6663, 9202, 7585),
    "H": ("PVM_1", "NeoC", 0.3031, 644, 311),
    "M": ("M4", "NonNeoC", 0.6596, 1382, 860.5),
    "N": ("PVM_3", "NonNeoC", 0.5064, 7799, 2323),
    "O": ("PVM_1", "NonNeoC", 0.8668, 2127, 2319.5),
    "P": ("HKG8", "NonNeoC", 0.3933, 864, 781),
}

# ---------------------------------------------------------------- xlsx reader
_NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
       "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships"}
_REL = "http://schemas.openxmlformats.org/package/2006/relationships"


def _col_index(ref):
    letters = re.match(r"[A-Z]+", ref).group(0)
    n = 0
    for ch in letters:
        n = n * 26 + (ord(ch) - 64)
    return n - 1


def read_xlsx(path, sheet):
    """One sheet as a DataFrame, header from the first row.

    Reads the XML directly: shared strings, inline strings, booleans and
    numbers. Empty cells are None. Numeric-looking columns are converted.
    """
    with zipfile.ZipFile(path) as z:
        ss = []
        if "xl/sharedStrings.xml" in z.namelist():
            root = ET.fromstring(z.read("xl/sharedStrings.xml"))
            for si in root.findall("m:si", _NS):
                ss.append("".join(t.text or "" for t in
                                  si.iter("{%s}t" % _NS["m"])))
        wb = ET.fromstring(z.read("xl/workbook.xml"))
        rels = ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
        targets = {r.get("Id"): r.get("Target")
                   for r in rels.findall("{%s}Relationship" % _REL)}
        rid = None
        for s in wb.iter("{%s}sheet" % _NS["m"]):
            if s.get("name") == sheet:
                rid = s.get("{%s}id" % _NS["r"])
        if rid is None:
            return None
        target = targets[rid].lstrip("/")
        if not target.startswith("xl/"):
            target = "xl/" + target
        root = ET.fromstring(z.read(target))
    rows = []
    for row in root.iter("{%s}row" % _NS["m"]):
        vals = {}
        for c in row.findall("m:c", _NS):
            t = c.get("t")
            v = c.find("m:v", _NS)
            if t == "s":
                val = ss[int(v.text)]
            elif t == "inlineStr":
                val = "".join(x.text or "" for x in c.iter("{%s}t" % _NS["m"]))
            elif t == "b":
                val = (v is not None and v.text == "1")
            elif v is None:
                val = None
            elif t == "str":
                val = v.text
            else:
                val = float(v.text)
            vals[_col_index(c.get("r"))] = val
        rows.append(vals)
    if not rows:
        return pd.DataFrame()
    width = max(max(r.keys(), default=-1) for r in rows) + 1
    grid = [[r.get(i) for i in range(width)] for r in rows]
    df = pd.DataFrame(grid[1:], columns=[str(x) for x in grid[0]])
    for col in df.columns:
        conv = pd.to_numeric(df[col], errors="coerce")
        if conv.notna().sum() == df[col].notna().sum() and df[col].notna().any():
            if not df[col].map(lambda x: isinstance(x, bool)).any():
                df[col] = conv
    return df


def need(path):
    if not os.path.isfile(path):
        raise SystemExit("input not found: %s" % path)
    return path


def p9_book(s, ds):
    return need(os.path.join(ROOT, "Part9", s["p9"],
                             "Part9_CP10K_%s_%s_GAM_Results_1.xlsx" % (s["tag"], ds)))


def p10_dir(s, ds):
    return os.path.join(ROOT, "Part10", s["p10"], "%s_%s" % (s["p10"], ds))


def p10_book(s, ds):
    return need(os.path.join(p10_dir(s, ds), "Part10_CP10K_%s_%s_%s_GAM_Results.xlsx"
                             % (s["tag"], s["stem"], ds)))


def p10_scores(s, ds):
    return need(os.path.join(p10_dir(s, ds), "Part10_CP10K_%s_%s_%s_Pooled_SKCS_scores.csv"
                             % (s["tag"], s["stem"], ds)))


def p10_log(s):
    return need(os.path.join(ROOT, "Part10", s["p10"],
                             "Part10_CP10K_%s_GAM_Log.txt" % s["tag"]))


def gam_data_files(s, ds):
    d = os.path.join(ROOT, "R_export_CP10K", EXPORT_DIR[s["set"]], EXPORT_DS[ds], "gam_data")
    files = sorted(glob.glob(os.path.join(d, "*_selected.csv")))
    if not files:
        raise SystemExit("no *_selected.csv files in %s" % d)
    return files


def p10b_ep(s, ds):
    fn = "Part10b_CP10K_%s_SKCS_%s_EP_TREP_values.csv" % (s["tag"], ds)
    hits = sorted(glob.glob(os.path.join(ROOT, "Part10b_SKCS", "*", fn)))
    if len(hits) != 1:
        raise SystemExit("expected one %s under Part10b_SKCS, found %d" % (fn, len(hits)))
    return hits[0]


def norm_region(r):
    r = str(r)
    return "Pooled" if r.endswith("Pooled") else r


def region_order(df):
    key = {(ds, r): i for ds in DATASETS for i, r in enumerate(REGIONS[ds])}
    sk = {s["set"]: i for i, s in enumerate(SETS)}
    return df.assign(_k1=df["set"].map(sk), _k2=df["region_set"].map({"NeoC": 0, "NonNeoC": 1}),
                     _k3=[key.get((a, b), 99) for a, b in zip(df["region_set"], df["region"])]
                     ).sort_values(["_k1", "_k2", "_k3"]).drop(columns=["_k1", "_k2", "_k3"])


# ---------------------------------------------------------------- shared helpers
checks = []


def show(label, got, ref, tol=0):
    ok = (got is not None and ref is not None and
          (abs(got - ref) <= tol if tol else got == ref))
    print("  %-52s %12s %12s  %s" % (label, got, ref, "ok" if ok else "CHECK"))
    checks.append(ok)
    return ok


def rank_avg(v):
    return pd.Series(v).rank(method="average").to_numpy(float)


def compute_skcs_py(x, y, chunk=500):
    """S_i (spear_norm) and K_i (kend_contrib) as in Part 10 compute_skcs()."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    n = len(x)
    rx, ry = rank_avg(x), rank_avg(y)
    spear_norm = (rx - rx.mean()) * (ry - ry.mean()) / (n - 1) / ((n - 1) / 4)
    cmd = np.empty(n)
    tx = np.empty(n)
    ty = np.empty(n)
    for a in range(0, n, chunk):
        b = min(a + chunk, n)
        sg = np.sign(x[a:b, None] - x[None, :])
        sc = np.sign(y[a:b, None] - y[None, :])
        prod = sg * sc                      # self-pairs are 0 because sg = 0 there
        cmd[a:b] = (prod > 0).sum(1) - (prod < 0).sum(1)
        tx[a:b] = (sg == 0).sum(1) - 1      # ties, self excluded
        ty[a:b] = (sc == 0).sum(1) - 1
    denom = np.sqrt((n - 1 - tx) * (n - 1 - ty))
    denom[denom == 0] = 1
    return spear_norm, cmd / denom


def match_rows(op_g, op_y, sc_g, sc_y, tol=MATCH_TOL):
    """Index into the score CSV for each Part 9 row, matched on (GOI, composite).

    Rows are taken in Part 9 order; each score row is used once. A Part 9 row
    with no unused score row within tol on both values is left at -1.
    """
    order = np.argsort(sc_g, kind="stable")
    gs = sc_g[order]
    used = np.zeros(len(sc_g), bool)
    out = np.full(len(op_g), -1, int)
    for i, (g, yv) in enumerate(zip(op_g, op_y)):
        lo = bisect.bisect_left(gs, g - tol)
        hi = bisect.bisect_right(gs, g + tol)
        cand = [j for j in order[lo:hi] if not used[j] and abs(sc_y[j] - yv) <= tol]
        if not cand:
            continue
        j = min(cand, key=lambda j: abs(sc_y[j] - yv) + abs(sc_g[j] - g))
        used[j] = True
        out[i] = j
    return out


def trep_all_cells(op):
    """TREP labels among all GOI-positive cells, Part 10 rule.

    EP = 1 - (y - yhat)^2 / (y - ybar)^2 with yhat the population-level
    prediction and ybar the unweighted mean (0 where y = ybar); the top
    round(n x DE) cells by EP are TREP, with DE the donor-weighted deviance
    explained of the same prediction, clamped to [0, 1]. Ties keep their
    original order, as R's order() does.
    """
    y = op["Observed_value"].to_numpy(float)
    yhat = op["Predicted_pop_value"].to_numpy(float)
    donor = op["Donor"].astype(str)
    w = 1.0 / donor.map(donor.value_counts()).to_numpy(float)
    ybar_w = np.sum(w * y) / np.sum(w)
    de = 1 - np.sum(w * (y - yhat) ** 2) / np.sum(w * (y - ybar_w) ** 2)
    null_sq = (y - y.mean()) ** 2
    with np.errstate(divide="ignore", invalid="ignore"):
        ep = np.where(null_sq > 0, 1 - (y - yhat) ** 2 / null_sq, 0.0)
    n = len(y)
    n_trep = max(1, min(int(round(n * max(0.0, min(1.0, de)))), n))
    idx = np.argsort(-ep, kind="stable")[:n_trep]
    lab = np.zeros(n, bool)
    lab[idx] = True
    return lab, de


# ---------------------------------------------------------------- A, B, C (per arm and region set)
A_parts, B_parts, C_rows = [], [], []
print("Per-cell checks (this run against logs, figures or the other source)")
print("  %-52s %12s %12s" % ("check", "this run", "reference"))

for s in SETS:
    is_hkg8 = s["set"] == "HKG8"
    for ds in DATASETS:
        tag = "%s %s" % (s["set"], ds)
        ref = REF[(s["set"], ds)]

        # Score CSV: one row per GOI-positive cell, pooled region set
        sc = pd.read_csv(p10_scores(s, ds))
        sc = sc[~sc["is_GOI_neg"].astype(bool)].reset_index(drop=True)
        n = len(sc)
        show(tag + " cells", n, ref["n"])

        # Donor of every barcode, Part 8 export
        gd = pd.concat([pd.read_csv(f, usecols=["barcode", "donor"]) for f in gam_data_files(s, ds)],
                       ignore_index=True).drop_duplicates("barcode")
        donor = sc["barcode"].map(gd.set_index("barcode")["donor"])
        show(tag + " cells without a donor", int(donor.isna().sum()), 0)

        # Composite and cHKG against the gene columns
        if is_hkg8:
            hk = [g for g in HKG_GENES if g != "ACTB"]
            comp_chk = sc[hk].mean(axis=1)
            chkg = pd.Series(np.nan, index=sc.index)
        else:
            comp_chk = sc[DFG_GENES].mean(axis=1)
            chkg = sc[HKG_GENES].mean(axis=1)
        show(tag + " composite = mean of its genes (max diff < 1e-6)",
             bool((comp_chk - sc["cDFG_expr"]).abs().max() < 1e-6), True)

        # SKCS-core EP, TREP and IP group, Part 10b
        e = pd.read_csv(p10b_ep(s, ds)).set_index("barcode")
        core = sc["is_skcs"].astype(bool).to_numpy()
        trep_core = sc["is_trep"].astype(bool).to_numpy()
        show(tag + " SKCS-core", int(core.sum()), ref["core"])
        show(tag + " TREP in SKCS-core", int(trep_core.sum()), ref["trep_core"])
        show(tag + " SKCS-core cells in Part 10b", int(sc.loc[core, "barcode"].isin(e.index).sum()),
             int(core.sum()))
        show(tag + " Part 10b cells outside the SKCS-core", int((~e.index.isin(sc.loc[core, "barcode"])).sum()), 0)
        t10b = sc["barcode"].map(e["is_TREP"]).fillna(False).astype(bool).to_numpy()
        show(tag + " TREP label, Part 10 = Part 10b", int((t10b != trep_core).sum()), 0)
        group = sc["barcode"].map(e["group_4"]).where(core, "Non-SKCS-core")
        for gname, gref in zip(GROUPS, REF_GROUPS[(s["set"], ds)]):
            show("%s %s" % (tag, gname), int((group == gname).sum()), gref)
        show(tag + " SKCS-core donors", int(donor[core].nunique()), ref["don_core"])

        # TREP among all GOI-positive cells, Part 9 predictions
        op = read_xlsx(p9_book(s, ds), "Observed_Predicted")
        op_reg = op["Region"].map(norm_region)
        opd = op[(op["Gene_Set"] == "DFG") & (op_reg == "Pooled")].reset_index(drop=True)
        lab, de_pop = trep_all_cells(opd)
        idx = match_rows(opd["GOI"].to_numpy(float), opd["Observed_value"].to_numpy(float),
                         sc["GOI_expr"].to_numpy(float), sc["cDFG_expr"].to_numpy(float))
        show(tag + " Part 9 cells matched to barcodes", int((idx >= 0).sum()), n)
        trep_all = np.zeros(n, bool)
        trep_all[idx[idx >= 0]] = lab[idx >= 0]
        fit_all = np.full(n, np.nan)
        fit_all[idx[idx >= 0]] = opd["Predicted_pop_value"].to_numpy(float)[idx >= 0]
        show(tag + " cells with a fitted value (all-cells GAM)", int(np.isfinite(fit_all).sum()), n)
        show(tag + " TREP among all cells (vs figure)", int(trep_all.sum()), ref["trep_all"])
        if not is_hkg8:
            oph = op[(op["Gene_Set"] == "HKG") & (op_reg == "Pooled")].reset_index(drop=True)
            if len(oph) == len(opd) and np.allclose(oph["GOI"], opd["GOI"]):
                d = np.abs(oph["Observed_value"].to_numpy(float)[idx >= 0] -
                           chkg.to_numpy(float)[idx[idx >= 0]])
                show(tag + " cHKG = Part 9 cHKG (max diff < 1e-4)", bool(d.max() < 1e-4), True)

        # Copula coordinates
        u = rank_avg(sc["GOI_expr"]) / (n + 1)
        v = rank_avg(sc["cDFG_expr"]) / (n + 1)
        u_sk = rank_avg(sc["spear_norm"]) / (n + 1)
        v_sk = rank_avg(sc["kend_contrib"]) / (n + 1)
        for c, (arm, rs_, x, rx, yv, ry) in EXPR_CELLS.items():
            if (arm, rs_) != (s["set"], ds):
                continue
            hit = np.where(((sc["GOI_expr"] - x).abs() < 5e-5).to_numpy() &
                           ((sc["cDFG_expr"] - yv).abs() < 5e-5).to_numpy())[0]
            ok = len(hit) == 1 and abs(u[hit[0]] * (n + 1) - rx) < 1e-6 and abs(v[hit[0]] * (n + 1) - ry) < 1e-6
            show("%s worked cell %s ranks as in the figure" % (tag, c), bool(ok), True)
        for c, (arm, rs_, ki, rs_rank, rk_rank) in SCORE_CELLS.items():
            if (arm, rs_) != (s["set"], ds):
                continue
            hit = np.where((sc["kend_contrib"] - ki).abs().to_numpy() < 5e-5)[0]
            ok = any(abs(u_sk[h] * (n + 1) - rs_rank) < 1e-6 and abs(v_sk[h] * (n + 1) - rk_rank) < 1e-6
                     for h in hit)
            show("%s worked cell %s ranks as in the figure" % (tag, c), bool(ok), True)

        a = pd.DataFrame({
            "arm": s["set"], "region_set": ds, "region": sc["region"], "barcode": sc["barcode"],
            "donor": donor, "GOI": sc["GOI_expr"], "composite": sc["cDFG_expr"],
            "fit_all_cells": fit_all, "cHKG": chkg})
        for g in DFG_GENES:
            a[g] = np.nan if is_hkg8 else sc[g]
        for g in HKG_GENES:
            a[g] = sc["GOI_expr"] if (is_hkg8 and g == "ACTB") else sc[g]
        a["S_i"] = sc["spear_norm"]
        a["K_i"] = sc["kend_contrib"]
        a["SKCS_core"] = core
        a["TREP_all_cells"] = trep_all
        a["TREP_SKCS_core"] = trep_core
        a["EP_SKCS_core"] = sc["barcode"].map(e["EP_value"]).where(core)
        a["group"] = group
        a["u"], a["v"], a["u_SK"], a["v_SK"] = u, v, u_sk, v_sk
        A_parts.append(a)

        # Housekeeping genes as GOI (HKG8), Part 10 section D5b
        if is_hkg8:
            hk = [g for g in HKG_GENES if g != "ACTB"]
            for g_alt in ALT_GOIS:
                comp = pd.concat([sc["GOI_expr"]] + [sc[c] for c in hk if c != g_alt], axis=1).mean(axis=1)
                goi = sc[g_alt]
                keep = (np.isfinite(goi) & (goi > 0) & np.isfinite(comp)).to_numpy()
                si, ki = compute_skcs_py(goi[keep], comp[keep])
                B_parts.append(pd.DataFrame({
                    "arm": "HKG8", "region_set": ds, "region": sc.loc[keep, "region"].to_numpy(),
                    "barcode": sc.loc[keep, "barcode"].to_numpy(), "GOI_gene": g_alt,
                    "S_i": si, "K_i": ki, "SKCS_core": (si > 0.05) & (ki > 0.5)}))

        # C: GAM fits, all cells (Part 9) and SKCS-core (Part 10)
        cc = read_xlsx(p9_book(s, ds), "Cell_Count_Audit")
        n_all = dict(zip(cc["Region"].map(norm_region), cc["N_Used_for_GAM"].astype(int)))
        dc = read_xlsx(p10_book(s, ds), "DE_Comparison")
        n_core = dict(zip(dc["region"].map(norm_region), dc["n_skcs"].astype(int)))
        de_cmp = {norm_region(r["region"]): (r["DE_full_pct"], r["DE_skcs_pct"]) for _, r in dc.iterrows()}
        for sheet in ("GAM_Performance_DFG", "GAM_Performance_HKG"):
            g = read_xlsx(p9_book(s, ds), sheet)
            # HKG8 writes GAM_Performance_HKG as a one-column note: no cHKG fit.
            if g is None or g.empty or "Region" not in g.columns:
                continue
            for _, r in g.iterrows():
                reg = norm_region(r["Region"])
                C_rows.append(dict(set=s["set"], region_set=ds, region=reg, population="All GOI-positive cells",
                                   composite=s["comp"] if r["Gene_Set"] == "DFG" else "cHKG",
                                   n_cells=n_all.get(reg), n_donors=int(r["N_Donors"]), k=int(r["k"]),
                                   lambda_REML=r["Lambda_REML"], smooth_EDF=r["Smooth_EDF"],
                                   DE_pct=100 * r["DE"], smooth_p=r["Smooth_p_value"],
                                   smooth_p_BH_FDR=r["Smooth_p_BH_FDR"]))
        gs = read_xlsx(p10_book(s, ds), "GAM_Performance_SKCS")
        for _, r in gs.iterrows():
            reg = norm_region(r["Region"])
            C_rows.append(dict(set=s["set"], region_set=ds, region=reg, population="SKCS-core",
                               composite=s["comp"] if r["Gene_Set"] == "DFG_SKCS" else "cHKG",
                               n_cells=n_core.get(reg), n_donors=int(r["N_Donors"]), k=int(r["k"]),
                               lambda_REML=r["Lambda_REML"], smooth_EDF=r["Smooth_EDF"],
                               DE_pct=100 * r["DE"], smooth_p=r["Smooth_p_value"],
                               smooth_p_BH_FDR=r["Smooth_p_BH_FDR"]))
        de_all = [r["DE_pct"] for r in C_rows if r["set"] == s["set"] and r["region_set"] == ds and
                  r["region"] == "Pooled" and r["composite"] == s["comp"]]
        show(tag + " DE all cells, Part 9 = Part 10 comparison", round(de_all[0], 2),
             round(de_cmp["Pooled"][0], 2), tol=0.011)
        show(tag + " DE SKCS-core, GAM sheet = Part 10 comparison", round(de_all[1], 2),
             round(de_cmp["Pooled"][1], 2), tol=0.011)

A_cells = pd.concat(A_parts, ignore_index=True)
B_cells = pd.concat(B_parts, ignore_index=True)
C_gam = region_order(pd.DataFrame(C_rows))
C_gam["_p"] = C_gam["population"].map({"All GOI-positive cells": 0, "SKCS-core": 1})
C_gam["_c"] = C_gam["composite"].map(lambda x: 1 if x == "cHKG" else 0)
C_gam = region_order(C_gam.sort_values(["_p", "_c"], kind="stable")).drop(columns=["_p", "_c"])

# ---------------------------------------------------------------- D to G (Part 9b, S_i floor)
def p9b_blocks(log_lines):
    """Split the Part 9b log into per-dataset blocks."""
    out, cur = {}, None
    for ln in log_lines:
        m = re.match(r"\s*Dataset:\s*(NeoC|NonNeoC)\s*$", ln)
        if m:
            cur = m.group(1)
            out[cur] = []
            continue
        if re.match(r"\s*-- Section 6", ln):
            cur = None
        if cur:
            out[cur].append(ln)
    return out


def p9b_log_tables(blocks):
    """Sheets D (mixture), F (posterior) and the S_i bins, read from the Part 9b log blocks."""
    F, H, I = [], [], []
    kde = {}
    for ds in DATASETS:
        lines = blocks.get(ds, [])
        waist_c = None
        for ln in lines:
            m = re.search(r"Waist spike: C(\d)", ln)
            if m:
                waist_c = int(m.group(1))
            m = re.search(r"Waist peak: \|Si_norm\| = ([\d.]+) \(density = ([\d.]+)\)", ln)
            if m:
                kde.setdefault(ds, {})["peak"] = (float(m.group(1)), float(m.group(2)))
            m = re.search(r"Density trough at \|Si_norm\| = ([\d.]+) \(density = ([\d.]+)\)", ln)
            if m:
                kde.setdefault(ds, {})["trough"] = (float(m.group(1)), float(m.group(2)))
        for ln in lines:
            m = re.search(r"C(\d): mean=([-\d.]+)\s+sd=([\d.]+)\s+prop=([\d.]+)%", ln)
            if m:
                F.append(dict(region_set=ds, component="C" + m.group(1),
                              mean_S_i=float(m.group(2)), SD_S_i=float(m.group(3)),
                              pct_of_cells=float(m.group(4)),
                              waist_component=(int(m.group(1)) == waist_c)))
            m = re.search(r"Si_norm = ([\d.]+): P\(waist\) = ([\d.]+)%", ln)
            if m:
                H.append(dict(region_set=ds, S_i=float(m.group(1)),
                              P_waist_pct=float(m.group(2))))
        bins = {}
        for ln in lines:
            m = re.match(r"\s{4}(Waist: \|Si\|<[\d.]+|\[[\d.]+, [\d.]+\))\s+(\d+)\s+([\d.]+)\s+([\d.]+)\s*$", ln)
            if m:
                bins[m.group(1)] = dict(region_set=ds, S_i_bin=m.group(1), n_cells=int(m.group(2)),
                                        median_rank_deviation=float(m.group(3)),
                                        IQR_rank_deviation=float(m.group(4)))
            m = re.match(r"\s{6}(Waist: \|Si\|<[\d.]+|\[[\d.]+, [\d.]+\))\s+W=([\d.]+), p=([\d.eE+-]+)", ln)
            if m and m.group(1) in bins:
                bins[m.group(1)].update(Shapiro_W=float(m.group(2)), Shapiro_p=float(m.group(3)))
        I.extend(bins.values())
    return pd.DataFrame(F), pd.DataFrame(H), pd.DataFrame(I)


p9b_log = open(need(P9B_LOG), encoding="utf-8", errors="replace").read().splitlines()
need(P9B_XLSX)
blocks = p9b_blocks(p9b_log)
F, H, I = p9b_log_tables(blocks)

# G: waist boundaries, from the workbook (rank-based, KDE) and Thresholds
wr = read_xlsx(P9B_XLSX, "Waist_RankBased")
th = read_xlsx(P9B_XLSX, "Thresholds")
G = []
for ds in DATASETS:
    for _, r in wr[wr["Dataset"] == ds].iterrows():
        G.append(dict(region_set=ds, boundary=r["Option"], K_ranks=r["K_ranks"],
                      theory_bound=r["Theory_bound"], empirical_bound=r["Empirical_bound"],
                      n_cells_below=r["N_waist_cells"], pct_cells_below=r["Pct_of_total"],
                      n_between_bound_and_0_05=r["N_recovery_zone"]))
    for _, r in th[th["Dataset"] == ds].iterrows():
        G.append(dict(region_set=ds, boundary=r["Approach"], empirical_bound=r["Si_thresh"]))
G = pd.DataFrame(G)
G["boundary"] = G["boundary"].str.replace("SKCS (adopted)", "S_i floor adopted", regex=False)

# J: rank-dispersion tests, from the workbook, KW df from the log
rv = read_xlsx(P9B_XLSX, "Rank_Verification")
kw_df = {}
for ds in DATASETS:
    for ln in blocks.get(ds, []):
        m = re.search(r"Kruskal-Wallis: chi2 = [\d.]+, df = (\d+)", ln)
        if m:
            kw_df[ds] = int(m.group(1))
J = pd.DataFrame([dict(region_set=r["Dataset"],
                       waist_bin_median_rank_deviation=r["Waist_median_GOI_dev"],
                       bin_0_04_0_06_median_rank_deviation=r["Thresh_median_GOI_dev"],
                       dispersion_ratio=r["Dispersion_ratio"],
                       KW_chi2=r["KW_chi2"], KW_df=kw_df.get(r["Dataset"]),
                       KW_p=r["KW_p"], KW_FDR=r["KW_FDR"],
                       Dunn_z_waist_vs_0_04_0_06=r["Dunn_z"],
                       Dunn_p_BH=r["Dunn_p_adj"],
                       Kendall_tau_S_i_vs_rank_deviation=r["Kendall_tau"],
                       Kendall_p=r["Kendall_p"], Kendall_FDR=r["Kendall_FDR"])
                  for _, r in rv.iterrows()])
J["region_set"] = pd.Categorical(J["region_set"], DATASETS)
J = J.sort_values("region_set")
gp = read_xlsx(P9B_XLSX, "GMM_Posterior")

# ---------------------------------------------------------------- H (pitcher summary, Part 10 logs)
pat_head = re.compile(r"\[(NeoC|NonNeoC)_Pooled(?:_([A-Za-z0-9-]+))?\] Global Spearman rho = ([-\d.]+)")
pat_tau = re.compile(r"\[(NeoC|NonNeoC)_Pooled(?:_([A-Za-z0-9-]+))?\] Global Kendall tau\s+= ([-\d.]+)")
pat_core = re.compile(r"\[(NeoC|NonNeoC)_Pooled(?:_([A-Za-z0-9-]+))?\] SKCS-core: (\d+) \(([\d.]+)%\)\s+\|\s+Non-SKCS: (\d+)")
SHOWN = {("M4", "MT-CO3"): "Figure 2B, 2D",
         ("PVM_3", "MT-CO3"): "Supplementary Figure 5",
         ("PVM_1", "MT-CO3"): "Supplementary Figure 5",
         ("HKG8", "ACTB"): "Figure 2B, 2D; Supplementary Figure 5",
         ("HKG8", "YWHAZ"): "Figure 2B, 2D",
         ("HKG8", "RAB7A"): "Figure 2B, 2D"}
pat_alt = re.compile(r"Alternate GOI ([A-Za-z0-9-]+): (\d+) cells \((\d+) dropped")


def pitcher_records(log_path, s):
    """Pooled r_s, tau_b and SKCS-core counts of one arm, from its Part 10 log."""
    rec, pending = {}, {}
    for ln in open(log_path, encoding="utf-8", errors="replace"):
        m = pat_alt.search(ln)
        if m:
            # the "Alternate GOI" line precedes the dataset-tagged lines
            pending[m.group(1)] = int(m.group(3))
        for pat, key in ((pat_head, "r_s"), (pat_tau, "tau_b"), (pat_core, "core")):
            m = pat.search(ln)
            if m:
                k = (m.group(1), m.group(2) or s["goi"])
                d = rec.setdefault(k, {})
                if key == "r_s" and m.group(2) in pending:
                    d["dropped"] = pending.pop(m.group(2))
                if key == "core":
                    d["n_SKCS_core"] = int(m.group(3))
                    d["n_cells"] = int(m.group(3)) + int(m.group(5))
                else:
                    d[key] = float(m.group(3))
    return rec


K = []
for s in SETS:
    rec = pitcher_records(p10_log(s), s)
    for (ds, goi), d in rec.items():
        alt = goi != s["goi"]
        K.append(dict(set=s["set"], region_set=ds, GOI=goi,
                      composite=("other eight housekeeping genes" if alt else s["comp"]),
                      shown_in=SHOWN.get((s["set"], goi), "not shown"),
                      n_cells=d.get("n_cells"), n_GOI_zero_excluded=d.get("dropped"),
                      r_s=d.get("r_s"), tau_b=d.get("tau_b"),
                      n_SKCS_core=d.get("n_SKCS_core"),
                      pct_SKCS_core=100 * d["n_SKCS_core"] / d["n_cells"] if d.get("n_cells") else None))
K = pd.DataFrame(K)
K["_k1"] = K["set"].map({s["set"]: i for i, s in enumerate(SETS)})
K["_k2"] = K["region_set"].map({"NeoC": 0, "NonNeoC": 1})
K["_k3"] = K["GOI"].map({"MT-CO3": 0, "ACTB": 0, "YWHAZ": 1, "RAB7A": 2, "EEF1A1": 3})
K = K.sort_values(["_k1", "_k2", "_k3"]).drop(columns=["_k1", "_k2", "_k3"])

# ---------------------------------------------------------------- remaining checks
print("\nSheet B and sheet A against the Part 10 log (sheet H)")
for ds in DATASETS:
    for g_alt in ALT_GOIS:
        b = B_cells[(B_cells["region_set"] == ds) & (B_cells["GOI_gene"] == g_alt)]
        k = K[(K["set"] == "HKG8") & (K["region_set"] == ds) & (K["GOI"] == g_alt)]
        tag = "HKG8 %s %s as GOI" % (ds, g_alt)
        show(tag + " cells", len(b), int(k["n_cells"].iloc[0]))
        show(tag + " SKCS-core", int(b["SKCS_core"].sum()), int(k["n_SKCS_core"].iloc[0]))
        if HAVE_SCIPY:
            a = A_cells[(A_cells["arm"] == "HKG8") & (A_cells["region_set"] == ds)].set_index("barcode")
            comp = a.loc[b["barcode"], [c for c in HKG_GENES if c != g_alt]].mean(axis=1).to_numpy()
            x = a.loc[b["barcode"], g_alt].to_numpy()
            show(tag + " r_s", round(spearmanr(x, comp).correlation, 4), float(k["r_s"].iloc[0]), tol=1e-4)
            show(tag + " tau_b", round(kendalltau(x, comp).correlation, 4), float(k["tau_b"].iloc[0]), tol=1e-4)
if HAVE_SCIPY:
    for s in SETS:
        for ds in DATASETS:
            a = A_cells[(A_cells["arm"] == s["set"]) & (A_cells["region_set"] == ds)]
            k = K[(K["set"] == s["set"]) & (K["region_set"] == ds) & (K["GOI"] == s["goi"])]
            tag = "%s %s" % (s["set"], ds)
            show(tag + " r_s from sheet A", round(spearmanr(a["GOI"], a["composite"]).correlation, 4),
                 float(k["r_s"].iloc[0]), tol=1e-4)
            show(tag + " tau_b from sheet A", round(kendalltau(a["GOI"], a["composite"]).correlation, 4),
                 float(k["tau_b"].iloc[0]), tol=1e-4)
for ds in DATASETS:
    b = C_gam[(C_gam["set"] == "M4") & (C_gam["region_set"] == ds) & (C_gam["region"] == "Pooled") &
              (C_gam["composite"] == "cDFG") & (C_gam["population"] == "All GOI-positive cells")]
    show("M4 %s DE all cells, Part 9 (%%)" % ds, round(float(b["DE_pct"].iloc[0]), 2),
         REF_DE_ALL_M4[ds], tol=0.06)

print("\nPart 9b: dispersion ratio = bin median / waist median")
for _, r in J.iterrows():
    show("%s ratio recomputed" % r["region_set"],
         round(r["bin_0_04_0_06_median_rank_deviation"] / r["waist_bin_median_rank_deviation"], 6),
         round(r["dispersion_ratio"], 6), tol=1e-6)
print("Part 9b: P(waist | S_i = 0.05), sheet F against GMM_Posterior")
for ds in DATASETS:
    h = H[(H["region_set"] == ds) & (np.isclose(H["S_i"], 0.05))]
    g = gp[gp["Dataset"] == ds]
    show("%s P(waist) at 0.05 (%%)" % ds, float(h["P_waist_pct"].iloc[0]),
         round(float(g["P_waist_at_005_pct"].iloc[0]), 4), tol=1e-4)

print("\n%d checks, %d need attention" % (len(checks), checks.count(False)))
if not all(checks):
    raise SystemExit("Checks failed; workbook not written.")

# ---------------------------------------------------------------- Contents
TITLE = "Supplementary Table [Fig2]. Cells of the four arms, GAM fits, the S_i floor and the pitcher plots."
DESC = ("Values behind Figures 2 to 4 and Supplementary Figures 2 to 10. Sheet A lists every "
        "GOI-positive cell once per arm and region set; Supplementary Tables [Fig3] and [Fig4] "
        "refer to it for cell-level values.")
SRC9 = "Part9_CP10K_{tag}_{ds}_GAM_Results_1.xlsx"
SRC10 = "Part10_CP10K_{tag}_{stem}_{ds}_GAM_Results.xlsx"
SRC10S = "Part10_CP10K_{tag}_{stem}_{ds}_Pooled_SKCS_scores.csv"
SRC9B = "Part9b_CP10K_M4_SiKi_threshold_comparison.xlsx"
SRC9BL = "Part9b_CP10K_M4_SiKi_Log.txt"

ren = {"set": "arm"}
sheets = [
    ("A_cells", A_cells,
     SRC10S + "; Part10b_CP10K_{tag}_SKCS_{ds}_EP_TREP_values.csv; " + SRC9 +
     "  [Observed_Predicted]; gam_data/*_selected.csv",
     "Part10_CP10K_SKCS_code.r; Part10b_CP10K_SKCS_TREP_code.r; this script",
     "One row per GOI-positive cell: donor, expression of the GOI, the composite, the fitted "
     "composite of the all-cells GAM, cHKG and the 17 genes, S_i, K_i, SKCS-core, TREP labels, EP, IP group and copula coordinates. "
     "Figures 2 to 4; Supplementary Figures 2, 3, 5 and 6 (cells) to 10. See notes 1 to 3."),
    ("B_HKG_pitcher_cells", B_cells, SRC10S + " (HKG8)", "this script (Part 10 compute_skcs)",
     "HKG8 cells with YWHAZ or RAB7A as GOI: S_i, K_i and SKCS-core. Expression is in sheet A "
     "under the same barcode. Figure 2B and 2D. See note 4."),
    ("C_GAM", C_gam.rename(columns=ren),
     SRC9 + "  [GAM_Performance_DFG, _HKG, Cell_Count_Audit]; " + SRC10 +
     "  [GAM_Performance_SKCS, DE_Comparison]",
     "Part9_CP10K_QC_GAM_code.r; Part10_CP10K_SKCS_code.r",
     "GAM of the composite (and of cHKG) on the GOI, over all GOI-positive cells and over the "
     "SKCS-core, per region and pooled: cells, donors, basis size, smoothing parameter, effective "
     "degrees of freedom, deviance explained and smooth-term test. Figure 1H, Figure 2 and "
     "Supplementary Figure 5. See note 5."),
    ("D_Si_mixture", F, SRC9BL, "Part9b_CP10K_SiKi_code.r",
     "Gaussian mixture of S_i in M4, by region set; the component centred on zero is the peak "
     "(waist). Supplementary Figure 4A and 4E."),
    ("E_Si_boundaries", G, SRC9B + "  [Waist_RankBased, Thresholds]", "Part9b_CP10K_SiKi_code.r",
     "Candidate upper ends of the peak of S_i at zero: density trough, three rank-based bounds, "
     "the mixture intersection and the S_i values where P(waist) falls to 5% and 1%. "
     "Supplementary Figure 4B, 4C, 4F and 4G."),
    ("F_Si_posterior", H, SRC9BL, "Part9b_CP10K_SiKi_code.r",
     "Posterior probability (per cent) that a cell with the given S_i belongs to the peak "
     "component. Supplementary Figure 4D and 4H."),
    ("G_Si_rank_tests", J, SRC9B + "  [Rank_Verification]", "Part9b_CP10K_SiKi_code.r",
     "Kruskal-Wallis test of GOI rank deviation over bins of S_i, Dunn's test of the peak bin "
     "against [0.04, 0.06), the dispersion ratio, and Kendall tau between S_i and rank deviation. "
     "Supplementary Figure 4D and 4H. See note 6."),
    ("H_pitcher_summary", K.rename(columns=ren), "Part10_CP10K_{tag}_GAM_Log.txt",
     "Part10_CP10K_SKCS_code.r",
     "Pooled r_s, tau_b and SKCS-core of each arm and region set, and of the HKG8 cells with "
     "EEF1A1, YWHAZ or RAB7A as GOI. Figure 2B and 2D and Supplementary Figure 5. See note 7."),
]

NOTES = [
    "1) Arms. M4: M4-selected PVM_3 cells with MT-CO3 >0, a subset of PVM_3. PVM_3 and PVM_1: "
    "MT-CO3-positive cells of each supertype. HKG8: M4 cells that are also ACTB-positive. GOI and "
    "composite: MT-CO3 and cDFG, the mean of the eight DFGs (M4, PVM_3, PVM_1); ACTB and cHKG8, the "
    "mean of the other eight housekeeping genes (HKG8). cHKG, the mean of the nine housekeeping "
    "genes, is given for the main arms only. Region sets: NeoC (DFC, MTG, AnG, STG, V1C, ITG) and "
    "NonNeoC (HIP, MEC, LEC, FI). Expression values are log1p CP10K. Cells are single-nucleus "
    "profiles.",
    "2) Sheet A. One row per GOI-positive cell of each arm and region set. The arms overlap, so a "
    "barcode can appear under more than one arm. S_i, K_i, the ranks and the copula coordinates "
    "are computed within the arm and region set, pooled over its regions: u = rank(GOI)/(n + 1), "
    "v = rank(composite)/(n + 1), u_SK = rank(S_i)/(n + 1), v_SK = rank(K_i)/(n + 1), ties "
    "averaged. SKCS-core: S_i >0.05 and K_i >0.5. Gene columns that are not part of an arm are blank.",
    "3) TREP and groups. EP = 1 - (y - yhat)^2 / (y - ybar)^2, with yhat the population-level GAM "
    "prediction (donor term excluded) and ybar the unweighted mean; the top round(n x DE) cells by "
    "EP are TREP, where DE is the donor-weighted deviance explained of the same prediction. "
    "TREP_all_cells uses the GAM on all GOI-positive cells (Part 9 predictions; labels recomputed "
    "here, because no part saves them per cell). fit_all_cells is that population-level prediction "
    "(donor term excluded), so composite - fit_all_cells is the residual behind R(x) in "
    "Supplementary Figure 6; TREP_SKCS_core and EP_SKCS_core use the GAM on "
    "the SKCS-core (Part 10). The two sets are defined in different populations, so one is not a "
    "subset of the other. group: pre-IP (composite below IP_y) or post-IP, TREP or non-TREP, for "
    "SKCS-core cells (Part 10b); Non-SKCS-core otherwise.",
    "4) Sheet B. For YWHAZ and RAB7A as GOI, HKG8 cells with zero expression of that gene were "
    "excluded, and the composite is the mean of the other eight housekeeping genes, including ACTB. "
    "S_i and K_i were recomputed here with the formula of Part 10, because Part 10 saves only the "
    "plots; the counts, r_s and tau_b agree with sheet H.",
    "5) Sheet C. Composite ~ GOI + s(GOI) + donor random intercept, REML, gamma 1.5, weights "
    "1/(cells of the donor) recomputed for each population, basis size chosen by penalized "
    "residual sum of squares between 3 and 10. DE_pct is measured against a model with the donor "
    "random intercept only. HKG8 has no cHKG fit, because its GOI is a housekeeping gene. Regions "
    "with fewer than 50 SKCS-core cells were not fitted. p-values below 2.2 \u00d7 10\u207b\u00b9\u2076, "
    "which R stores as 0, are written as < 2.2 \u00d7 10\u207b\u00b9\u2076.",
    "6) S_i floor (sheets D to G). M4 only, pooled per region set. Sheet G uses the mixture "
    "intersection as the upper end of the peak (0.014 NeoC, 0.019 NonNeoC). GOI rank deviation is "
    "|rank(GOI)/n - 0.5|, the distance from the middle GOI rank as a fraction of n; the dispersion "
    "ratio is the median deviation of cells with S_i in [0.04, 0.06) divided by that of cells in the "
    "peak. P(waist) is in per cent, rounded to four decimals as logged (NeoC at S_i = 0.05 is "
    "2.95e-5 per cent unrounded).",
    "7) Sheet H. Read from the Part 10 logs. EEF1A1 as GOI is computed by Part 10 but not shown in "
    "Figure 2.",
    "8) build_SuppTable_Fig2.py prints every recomputed value next to the value in the logs, the "
    "figures or another source file, and writes the workbook only if all checks pass.",
]

# ---------------------------------------------------------------- write (write-only mode)
P_FLOOR = 2.2e-16
P_FLOOR_TEXT = "< 2.2 \u00d7 10\u207b\u00b9\u2076"
P_FLOOR_COLS = {"smooth_p", "smooth_p_BH_FDR", "Kendall_p", "Kendall_FDR"}
FLOAT_KEEP = re.compile(r"(^|_)(p|FDR|p_BH|smooth_p|KW_p|Kendall_p|Shapiro_p)($|_)")
RAW_SHEETS = {"A_cells", "B_HKG_pitcher_cells"}

ARIAL = Font(name="Arial", size=10)
HEAD = Font(name="Arial", size=10, bold=True)


def tidy(df):
    df = df.copy()
    for c in df.columns:
        if df[c].dtype.kind == "f" and not FLOAT_KEEP.search(c):
            df[c] = df[c].round(6)
    return df


def cell(ws, v, font=ARIAL, wrap=False):
    c = WriteOnlyCell(ws, value=v)
    c.font = font
    if wrap:
        c.alignment = Alignment(vertical="top", wrap_text=True)
    return c


def py(v):
    if isinstance(v, (np.bool_, bool)):
        return bool(v)
    if isinstance(v, np.integer):
        return int(v)
    if isinstance(v, (float, np.floating)):
        return None if not np.isfinite(v) else float(v)
    if v is pd.NA or v is None:
        return None
    return v


wb = Workbook(write_only=True)
ws = wb.create_sheet("Contents")
for col, wdt in zip("ABCDE", (22, 46, 30, 8, 82)):
    ws.column_dimensions[col].width = wdt
ws.freeze_panes = "A5"
ws.append([cell(ws, TITLE, HEAD)])
ws.append([cell(ws, DESC)])
ws.append([])
ws.append([cell(ws, h, HEAD) for h in ("Sheet", "Source file", "Source script", "Rows", "Description")])
for name, df, src, scr, desc in sheets:
    ws.append([cell(ws, name, HEAD), cell(ws, src), cell(ws, scr),
               cell(ws, "=COUNTA(%s!A:A)-1" % name), cell(ws, desc, wrap=True)])
ws.append([])
ws.append([cell(ws, "Notes", HEAD)])
for n_ in NOTES:
    ws.append([cell(ws, n_)])

for name, df, *_ in sheets:
    if name not in RAW_SHEETS:
        df = tidy(df)
    sh = wb.create_sheet(name[:31])
    for i, col in enumerate(df.columns, start=1):
        vals = [len(str(v)) for v in df[col].head(200)]
        sh.column_dimensions[get_column_letter(i)].width = min(max([len(str(col))] + vals) + 2, 46)
    sh.freeze_panes = "A2"
    sh.append([cell(sh, str(c), HEAD) for c in df.columns])
    p_idx = [j for j, c in enumerate(df.columns) if c in P_FLOOR_COLS]
    for row in df.itertuples(index=False):
        out = [py(v) for v in row]
        for j in p_idx:
            if isinstance(out[j], float) and out[j] < P_FLOOR:
                out[j] = P_FLOOR_TEXT
        sh.append([cell(sh, v) for v in out])

os.makedirs(OUT_DIR, exist_ok=True)
wb.save(OUT)
WRITTEN = True
print("\nsaved", OUT)
for name, df, *_ in sheets:
    print("  %-20s %d rows" % (name, len(df)))
