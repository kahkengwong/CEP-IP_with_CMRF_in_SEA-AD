#!/usr/bin/env python3
"""
build_SuppTable_S12.py

Builds the workbook behind Supplementary Table S12 (Figures 9 to 11 and
Supplementary Figure S14: clinical association, benchmark against
competitors, boundary ablation). It reads the workbooks of Part 10h and its
stages 1 (comparator), 2 (CogVsPath) and 4 (multivariate) for the four arms,
the Part 11 CogVsPath workbook, the Part 11 BandAblation workbook and its
Figure 11 export, and the plotted-data CSVs of the Figure 13 comparator
script. It writes one plain workbook (Arial 10, bold header, frozen top row)
with the 19 data sheets A_donor_frame to S_SuppFig14_tissue_disease, and a
build log with a PASS/WARN/FAIL check list ending in "VERDICT: PASS" or
"VERDICT: FAIL". Each data sheet is checked against a second source (another
workbook or the plotted-data CSVs of the figure script) or against values
printed in the figure. Leave-one-region-out results are not included.
finalize_supp_tables.py then keeps the data sheets, applies the published
text edits, writes the Contents sheet and saves Supplementary Table S12.

Usage:
    python build_SuppTable_S12.py
    python build_SuppTable_S12.py --root /data
Options:
    --root  pipeline output folder (default: $SEAAD_ROOT if set, otherwise
            the current folder)
    --out   output folder (default: <root>/Manuscript/Supp_Tables/build)

Inputs (under --root; <fold> is Mic13_PVM3, PVM_3, PVM_1 or HKG8 and <tag>
is M4, PVM_3, PVM_1 or HKG8):
    Part10h_Prognosis/Part10h_Prognosis_<fold>/Part10h_CP10K_<tag>_Prognosis.xlsx
        (Part10h_CP10K_Prognosis_code.r)
    Part10h_Prognosis_1_Comparator/Part10h_Prognosis_1_Comparator_<fold>/
        Part10h_CP10K_<tag>_Comparator.xlsx
        (Part10h_CP10K_Prognosis_1_Comparator_code.r)
    Part10h_Prognosis_2_CogVsPath/Part10h_Prognosis_2_CogVsPath_<fold>/
        Part10h_CP10K_<tag>_CogVsPath.xlsx
        (Part10h_CP10K_Prognosis_2_CogVsPath_code.r)
    Part10h_Prognosis_4_Multivariate/Part10h_Prognosis_4_Multivariate_<fold>/
        Part10h_CP10K_<tag>_Multivariate.xlsx
        (Part10h_CP10K_Prognosis_4_Multivariate_code.r)
    Part10h_Prognosis_1_Comparator/Fig13/Fig13A_DeltaR2_grid_plotted.csv,
        Fig13B_Duel_verdicts_plotted.csv, Fig13C_Pub_Decomposition_plotted.csv,
        Fig13C2_Pub_Contrast_plotted.csv (Part10h_CP10K_Comparator_Fig13_code.r)
    Part11_Identity/Part11_CogVsPath/Part11_CP10K_CogVsPath.xlsx
        (Part11_CP10K_CogVsPath_code.r)
    Part11_Identity/Part11_Identity_BandAblation_P10hDenom/
        Part11_Identity_BandAblation.xlsx (Part11_Identity_BandAblation_code.r)
    Part11_Identity/Part11_BandAblation_Fig11Export/BandAblation_cells.csv,
        BandAblation_Fig11B_curves.csv, BandAblation_set_summary.csv
        (Part11_BandAblation_Fig11Export_code.r, run after the BandAblation
        script)

Output:
    <out>/SuppTable_S12_Prognosis_&_Ablation.xlsx
    <out>/SuppTable_S12_build_log.txt
"""
import argparse
import csv
import datetime
import math
import os
import statistics
import sys
import traceback

import numpy as np
import openpyxl
from openpyxl.styles import Font
from openpyxl.utils import get_column_letter

# ---- Settings ----------------------------------------------------------------
DEFAULT_ROOT = os.environ.get("SEAAD_ROOT", os.getcwd())
OUT_NAME = "SuppTable_S12_Prognosis_&_Ablation.xlsx"
LOG_NAME = "SuppTable_S12_build_log.txt"

# (manuscript label, folder suffix, file tag)
ARMS = [("M4", "Mic13_PVM3", "M4"), ("PVM_3", "PVM_3", "PVM_3"),
        ("PVM_1", "PVM_1", "PVM_1"), ("HKG8", "HKG8", "HKG8")]
ARM_ORDER = [a[0] for a in ARMS]
DS_ORDER = ["NeoC", "NonNeoC"]
N_DONORS = {"NeoC": 84, "NonNeoC": 83}

OUT9 = ["Cognitive status", "MMSE score", "CASI score", "CERAD score",
        "Braak stage", "Thal phase", "ADNC rank", "CAA score", "APOE e4 dose"]
COG3 = ["Cognitive status", "MMSE score", "CASI score"]
PATH5 = ["CERAD score", "Braak stage", "Thal phase", "ADNC rank", "CAA score"]
PATH4 = ["ADNC rank", "CERAD score", "Braak stage", "Thal phase"]
COG2 = ["MMSE score", "CASI score"]

POP9A = {"pct_pre": "Module pre-IP", "pct_post": "Module post-IP",
         "pct_excl_pre": "Non-module pre-IP",
         "pct_excl_post": "Non-module post-IP"}
POP_ORDER = list(POP9A.values())

BLOCK9B = {"cov4": "none", "+ADNC": "ADNC",
           "+BTC": "Braak + Thal + CERAD",
           "+ALL": "Braak + Thal + CERAD + ADNC"}
BLOCK9C = {"cov4": ("cov4", "none"), "+COG": ("+MMSE + CASI", "MMSE + CASI")}

PRED10 = [("CMRF module Post-IP", "Mod"),
          ("Cluster proportion (supertype)", "CP"),
          ("Module score (AUCell)", "AUC"),
          ("Changepoint split (segmented)", "CS"),
          ("Module score (UCell)", "UC"),
          ("Module score (AddModuleScore)", "AMS")]
PRED10_ABBR = dict(PRED10)
PRED10_ORDER = [p[0] for p in PRED10]
SCORE_ORDER = ["AUCell", "UCell", "AddModuleScore"]

RUNG = {"A0_OLS": "OLS", "A1_spline_fixed": "Fixed spline",
        "A2_spline_local": "Local spline", "A3_CMRF": "CMRF"}
RUNG_ORDER = ["OLS", "Fixed spline", "Local spline", "CMRF"]
GROUP11 = {"M4": "M4", "PVM_3_all": "PVM_3", "PVM_1_all": "PVM_1",
           "HKG8": "HKG8"}
ARM11D = {"Mic13_PVM3": "M4", "PVM3_all": "PVM_3", "PVM1_all": "PVM_1",
          "HKG8_matched": "HKG8"}
OUT11D = {"MMSE": "MMSE score", "CASI": "CASI score", "ADNC": "ADNC rank",
          "CERAD": "CERAD score", "Braak": "Braak stage", "Thal": "Thal phase"}
GROUP13 = {"M4 (masked, primary)": "M4", "PVM_3, all cells": "PVM_3",
           "PVM_1, all cells": "PVM_1", "HKG8 (housekeeping control)": "HKG8",
           "M4": "M4"}

# Clinical columns of the comparator Donor_Frame -> S12 column names.
CLIN_COLS = [("CogStatus_num", "cognitive_status"), ("MMSE_num", "MMSE"),
             ("CASI_num", "CASI"), ("CERAD_num", "CERAD"),
             ("Braak_num", "Braak"), ("Thal_num", "Thal"),
             ("ADNC_rank", "ADNC"), ("CAA_num", "CAA"),
             ("APOE_e4_dose", "APOE_e4_dose"), ("Age_num", "age_at_death"),
             ("Sex_num", "sex"), ("PMI_num", "PMI"), ("RIN_num", "RIN"),
             ("BrainPH_num", "brain_pH"),
             ("BrainWeight_num", "fresh_brain_weight")]

# ---- Values printed in the figures (used as checks) --------------------------
# Figure 10C: (FDR, max VIF, verdict as shown), order: cov4 MMSE, cov4 CASI,
# +BTC MMSE, +BTC CASI. PVM_3 NonNeoC VIFs and PVM_1 NonNeoC +BTC CASI FDR are
# the source values.
DH, DN, UR = "HOLDS", "doesn't hold", "unreadable"
EXP_10C = {
    ("NeoC", "M4"): [(0.011, 3.36, DH), (0.030, 3.16, DH), (0.014, 3.57, DH),
                     (0.053, 3.36, DN)],
    ("NeoC", "PVM_3"): [(0.098, 16.37, UR), (0.052, 14.26, UR),
                        (0.137, 17.17, UR), (0.105, 14.79, UR)],
    ("NeoC", "PVM_1"): [(0.416, 13.58, UR), (0.416, 12.16, UR),
                        (0.458, 14.22, UR), (0.458, 12.66, UR)],
    ("NeoC", "HKG8"): [(0.589, 1.38, DN), (0.449, 1.51, DN),
                       (0.845, 1.46, DN), (0.845, 1.69, DN)],
    ("NonNeoC", "M4"): [(0.254, 2.32, DN), (0.284, 2.38, DN),
                        (0.225, 2.33, DN), (0.225, 2.39, DN)],
    ("NonNeoC", "PVM_3"): [(0.128, 8.90, UR), (0.052, 8.33, UR),
                           (0.166, 9.48, UR), (0.106, 8.85, UR)],
    ("NonNeoC", "PVM_1"): [(0.973, 10.86, UR), (0.239, 11.09, UR),
                           (0.733, 11.33, UR), (0.086, 11.43, UR)],
    ("NonNeoC", "HKG8"): [(0.139, 1.58, DN), (0.245, 1.76, DN),
                          (0.219, 1.61, DN), (0.309, 1.76, DN)],
}
EXP_RS = {("M4", "NeoC"): 0.815, ("PVM_3", "NeoC"): 0.972,
          ("PVM_1", "NeoC"): 0.928, ("HKG8", "NeoC"): 0.426,
          ("M4", "NonNeoC"): 0.729, ("PVM_3", "NonNeoC"): 0.975,
          ("PVM_1", "NonNeoC"): 0.954, ("HKG8", "NonNeoC"): 0.513}
# Figure 10D (M4): (beta module, beta score, FDR), score order AUCell, UCell,
# AddModuleScore.
EXP_10D = {
    ("NeoC", "cov4", "MMSE score"): [(-0.427, -0.073, 0.006),
                                     (-0.423, -0.077, 0.008),
                                     (-0.389, -0.128, 0.011)],
    ("NeoC", "cov4", "CASI score"): [(-0.362, -0.163, 0.015),
                                     (-0.350, -0.179, 0.022),
                                     (-0.312, -0.241, 0.028)],
    ("NeoC", "+BTC", "MMSE score"): [(-0.420, -0.025, 0.010),
                                     (-0.418, -0.027, 0.012),
                                     (-0.383, -0.083, 0.017)],
    ("NeoC", "+BTC", "CASI score"): [(-0.320, -0.098, 0.039),
                                     (-0.312, -0.109, 0.048),
                                     (-0.272, -0.178, 0.070)],
    ("NonNeoC", "cov4", "MMSE score"): [(-0.212, -0.270, 0.052),
                                        (-0.201, -0.286, 0.065),
                                        (-0.213, -0.297, 0.081)],
    ("NonNeoC", "cov4", "CASI score"): [(-0.177, -0.256, 0.052),
                                        (-0.163, -0.280, 0.065),
                                        (-0.134, -0.380, 0.081)],
    ("NonNeoC", "+BTC", "MMSE score"): [(-0.213, -0.196, 0.086),
                                        (-0.204, -0.211, 0.101),
                                        (-0.212, -0.224, 0.150)],
    ("NonNeoC", "+BTC", "CASI score"): [(-0.191, -0.137, 0.086),
                                        (-0.181, -0.157, 0.101),
                                        (-0.144, -0.262, 0.150)],
}
# Figure 11C: overlap median, min, max (None for CMRF), C1, C2, C3, C4,
# median |beta| of discarded cells.
EXP_11C = {
    ("three real groups", "OLS"): (0.821, 0.792, 0.875, 44, 15, 4, 17, 0.129),
    ("three real groups", "Fixed spline"): (0.798, 0.754, 0.857, 44, 16, 3,
                                            17, 0.123),
    ("three real groups", "Local spline"): (0.748, 0.688, 0.813, 43, 13, 10,
                                            17, 0.197),
    ("three real groups", "CMRF"): (None, None, None, 44, 13, 0, 17, 0.106),
    ("HKG8 control", "OLS"): (0.197, 0.145, 0.249, 0, 0, 3, 0, 0.120),
    ("HKG8 control", "Fixed spline"): (0.182, 0.147, 0.216, 3, 0, 0, 0,
                                       0.134),
    ("HKG8 control", "Local spline"): (0.188, 0.154, 0.221, 6, 1, 0, 0,
                                       0.110),
    ("HKG8 control", "CMRF"): (None, None, None, 7, 0, 0, 1, 0.103),
}
EXP_11C_OF = {"three real groups": (54, 36, 54, 18),
              "HKG8 control": (18, 12, 18, 6)}
# Figure 9C: rows significant after adding MMSE + CASI.
EXP_9C_SIG = {("M4", "NeoC", o) for o in ("ADNC rank", "CERAD score",
                                           "Thal phase")} | \
             {("PVM_3", "NeoC", o) for o in ("ADNC rank", "CERAD score",
                                              "Thal phase")}
TOL3 = 0.0005 + 1e-9
TOL2 = 0.005 + 1e-9
TOL_ROUND4 = 6e-5          # values rounded to 4 dp in one source


# ---- Utilities ---------------------------------------------------------------
class Tee:
    def __init__(self, path):
        self.f = open(path, "w", encoding="utf-8")
        self.out = sys.stdout

    def write(self, s):
        self.out.write(s)
        self.f.write(s)

    def flush(self):
        self.out.flush()
        self.f.flush()


def num(v):
    if v is None or isinstance(v, bool):
        return v
    if isinstance(v, (int, float)):
        return float(v)
    s = str(v).strip()
    if s in ("", "NA", "NaN", "nan", "Inf", "-Inf"):
        return None
    try:
        return float(s)
    except ValueError:
        return None


def tb(v):
    if isinstance(v, bool):
        return v
    return str(v).strip().upper() in ("TRUE", "T", "1")


def read_xlsx(path, sheets):
    """Returns {sheet: [row dict, ...]} for the requested sheets."""
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    try:
        out = {}
        for sh in sheets:
            if sh not in wb.sheetnames:
                raise KeyError(f"sheet '{sh}' not found in {path}")
            it = wb[sh].iter_rows(values_only=True)
            hdr = [str(h) for h in next(it)]
            rows = []
            for r in it:
                if r is None or all(v is None for v in r):
                    continue
                rows.append(dict(zip(hdr, r)))
            out[sh] = rows
        return out
    finally:
        wb.close()


def read_csv(path):
    with open(path, "r", encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh))


def near(a, b, tol):
    a, b = num(a), num(b)
    if a is None or b is None:
        return a is None and b is None
    return abs(a - b) <= tol


def _rank(x):
    """Average ranks (1-based) with ties averaged."""
    order = sorted(range(len(x)), key=lambda i: x[i])
    r = [0.0] * len(x)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and x[order[j + 1]] == x[order[i]]:
            j += 1
        for k in range(i, j + 1):
            r[order[k]] = (i + j) / 2.0 + 1.0
        i = j + 1
    return r


def _betacf(a, b, x):
    """Continued fraction for the regularized incomplete beta."""
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c, d = 1.0, 1.0 - qab * x / qap
    d = 1.0 / (d if abs(d) > 1e-300 else 1e-300)
    h = d
    for m in range(1, 300):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        d = 1.0 / (d if abs(d) > 1e-300 else 1e-300)
        c = 1.0 + aa / c if abs(c) > 1e-300 else 1e300
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        d = 1.0 / (d if abs(d) > 1e-300 else 1e-300)
        c = 1.0 + aa / c if abs(c) > 1e-300 else 1e300
        de = d * c
        h *= de
        if abs(de - 1.0) < 1e-12:
            break
    return h


def _betai(a, b, x):
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    lb = (math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
          + a * math.log(x) + b * math.log(1.0 - x))
    if x < (a + 1.0) / (a + b + 2.0):
        return math.exp(lb) * _betacf(a, b, x) / a
    return 1.0 - math.exp(lb) * _betacf(b, a, 1.0 - x) / b


def spearman(x, y):
    """Spearman rho with a two-sided t-approximation p-value."""
    pairs = [(a, b) for a, b in zip(x, y) if a is not None and b is not None]
    n = len(pairs)
    if n < 4:
        return n, None, None
    rx = _rank([a for a, _ in pairs])
    ry = _rank([b for _, b in pairs])
    mx, my = sum(rx) / n, sum(ry) / n
    sxy = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    sxx = sum((a - mx) ** 2 for a in rx)
    syy = sum((b - my) ** 2 for b in ry)
    if sxx == 0 or syy == 0:
        return n, None, None
    rho = sxy / math.sqrt(sxx * syy)
    if abs(rho) >= 1.0:
        return n, rho, 0.0
    df = n - 2
    t2 = rho * rho * df / (1.0 - rho * rho)
    return n, rho, _betai(df / 2.0, 0.5, df / (df + t2))


def bh(pvals):
    """Benjamini-Hochberg adjusted p-values (None kept as None)."""
    idx = [i for i, p in enumerate(pvals) if p is not None]
    m = len(idx)
    out = [None] * len(pvals)
    prev = 1.0
    for rank, i in enumerate(sorted(idx, key=lambda k: pvals[k],
                                    reverse=True)):
        k = m - rank
        prev = min(prev, pvals[i] * m / k)
        out[i] = prev
    return out


class Checks:
    def __init__(self):
        self.items = []

    def add(self, cid, ok, msg, warn_only=False):
        status = "PASS" if ok else ("WARN" if warn_only else "FAIL")
        self.items.append((cid, status, msg))
        print(f"  [{status}] {cid}: {msg}")

    def info(self, msg):
        print(f"         {msg}")

    def counts(self):
        c = {"PASS": 0, "WARN": 0, "FAIL": 0}
        for _, s, _ in self.items:
            c[s] += 1
        return c



# ---- Paths -------------------------------------------------------------------
def paths(root):
    P = {}
    for lab, fold, tag in ARMS:
        P[("p10h", lab)] = os.path.join(
            root, "Part10h_Prognosis", f"Part10h_Prognosis_{fold}",
            f"Part10h_CP10K_{tag}_Prognosis.xlsx")
        P[("comp", lab)] = os.path.join(
            root, "Part10h_Prognosis_1_Comparator",
            f"Part10h_Prognosis_1_Comparator_{fold}",
            f"Part10h_CP10K_{tag}_Comparator.xlsx")
        P[("cvp", lab)] = os.path.join(
            root, "Part10h_Prognosis_2_CogVsPath",
            f"Part10h_Prognosis_2_CogVsPath_{fold}",
            f"Part10h_CP10K_{tag}_CogVsPath.xlsx")
        P[("mv", lab)] = os.path.join(
            root, "Part10h_Prognosis_4_Multivariate",
            f"Part10h_Prognosis_4_Multivariate_{fold}",
            f"Part10h_CP10K_{tag}_Multivariate.xlsx")
    P["p11cvp"] = os.path.join(root, "Part11_Identity", "Part11_CogVsPath",
                               "Part11_CP10K_CogVsPath.xlsx")
    P["p11ba"] = os.path.join(root, "Part11_Identity",
                              "Part11_Identity_BandAblation_P10hDenom",
                              "Part11_Identity_BandAblation.xlsx")
    exp = os.path.join(root, "Part11_Identity",
                       "Part11_BandAblation_Fig11Export")
    P["ba_cells"] = os.path.join(exp, "BandAblation_cells.csv")
    P["ba_curves"] = os.path.join(exp, "BandAblation_Fig11B_curves.csv")
    P["ba_sum"] = os.path.join(exp, "BandAblation_set_summary.csv")
    f13 = os.path.join(root, "Part10h_Prognosis_1_Comparator", "Fig13")
    P["f13A"] = os.path.join(f13, "Fig13A_DeltaR2_grid_plotted.csv")
    P["f13B"] = os.path.join(f13, "Fig13B_Duel_verdicts_plotted.csv")
    P["f13C"] = os.path.join(f13, "Fig13C_Pub_Decomposition_plotted.csv")
    P["f13C2"] = os.path.join(f13, "Fig13C2_Pub_Contrast_plotted.csv")
    return P


# ---- Main --------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=DEFAULT_ROOT)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    root = args.root
    out_dir = args.out or os.path.join(root, "Manuscript", "Supp_Tables",
                                       "build")
    os.makedirs(out_dir, exist_ok=True)
    log_path = os.path.join(out_dir, LOG_NAME)
    tee = Tee(log_path)
    sys.stdout = tee
    failed = False
    try:
        failed = run(root, out_dir, log_path) == "FAIL"
    except Exception:
        # An unexpected error (for example a missing sheet or column) still
        # ends the log with a verdict.
        failed = True
        traceback.print_exc(file=sys.stdout)
        print("VERDICT: FAIL")
    finally:
        sys.stdout = tee.out
        tee.f.close()
    print(f"Log written: {log_path}")
    if failed:
        sys.exit(1)


def run(root, out_dir, log_path):
    ck = Checks()
    print("build_SuppTable_S12")
    print(f"Started : {datetime.datetime.now():%Y-%m-%d %H:%M:%S}")
    print(f"Builder : {os.path.abspath(__file__)}")
    print(f"Root    : {root}")
    print(f"Output  : {os.path.join(out_dir, OUT_NAME)}")
    P = paths(root)

    # ---- 0. Sources --------------------------------------------------------
    print("\n0. SOURCES")
    missing = [p for p in P.values() if not os.path.isfile(p)]
    for k, p in P.items():
        if os.path.isfile(p):
            ts = datetime.datetime.fromtimestamp(os.path.getmtime(p))
            print(f"  {os.path.relpath(p, root)}  [{ts:%Y-%m-%d %H:%M}]")
    ck.add("S0", not missing, f"{len(P) - len(missing)} of {len(P)} source "
           "files found" + ("" if not missing else "; missing: " +
                            "; ".join(os.path.relpath(m, root)
                                      for m in missing)))
    if missing:
        return finish(ck, None, out_dir)

    p10h, comp, cvp, mv = {}, {}, {}, {}
    for lab in ARM_ORDER:
        p10h[lab] = read_xlsx(P[("p10h", lab)],
                              ["Multivariate", "Donor_pct_NeoC",
                               "Donor_pct_NonNeoC", "Cov4_vs_Cov7",
                               "LeaveOneIn_Scan"])
        comp[lab] = read_xlsx(P[("comp", lab)],
                              ["Donor_Frame", "Marginal", "Duels"])
        cvp[lab] = read_xlsx(P[("cvp", lab)], ["Models", "Reverse"])
        mv[lab] = read_xlsx(P[("mv", lab)], ["Decomposition", "Correlations"])
    p11 = read_xlsx(P["p11cvp"], ["Contrast"])
    ba = read_xlsx(P["p11ba"], ["Overlap", "Donor_Assoc", "Full_Panel",
                                "Scorecard", "Coverage", "Donor_Pct"])
    sheets = []          # (name, header, rows)

    # ---- A. Donor frame ----------------------------------------------------
    print("\nA. DONOR FRAME")
    ba_pct = {}
    for r in ba["Donor_Pct"]:
        arm = GROUP11.get(r["group"])
        ba_pct[(arm, r["dataset"], RUNG[r["rung"]], str(r["donor"]))] = r
    hdr_a = (["arm", "region_set", "donor_id", "n_supertype_cells",
              "pct_module_preIP", "pct_module_postIP",
              "pct_nonmodule_preIP", "pct_nonmodule_postIP",
              "n_microglia_PVM_cells", "B1_cluster_proportion_pct",
              "B2_AUCell", "B3_changepoint_split_pct", "B4_UCell",
              "B5_AddModuleScore",
              "pct_module_postIP_OLS", "pct_module_postIP_fixed_spline",
              "pct_module_postIP_local_spline",
              "pct_nonmodule_postIP_OLS", "pct_nonmodule_postIP_fixed_spline",
              "pct_nonmodule_postIP_local_spline",
              "pct_nonmodule_postIP_CMRF_boundary"]
             + [c[1] for c in CLIN_COLS])
    rows_a = []
    ok_n = ok_set = True
    d_mod = d_ba = d_nm = 0.0
    nm_diff = []        # (arm, ds, donor, diff in pp, diff in cells)
    den_ok = True
    donor_frame = {}
    for lab in ARM_ORDER:
        for ds in DS_ORDER:
            pd = {str(r["donor_id"]): r for r in p10h[lab][f"Donor_pct_{ds}"]}
            cf = {str(r[".donor"]): r for r in comp[lab]["Donor_Frame"]
                  if r["dataset"] == ds}
            bad = {d for (a, s, rg, d) in ba_pct
                   if a == lab and s == ds and rg == "CMRF"}
            if len(pd) != N_DONORS[ds]:
                ok_n = False
            if not (set(pd) == set(cf) == bad):
                ok_set = False
                ck.info(f"donor sets differ for {lab} {ds}: Part10h "
                        f"{len(pd)}, comparator {len(cf)}, BandAblation "
                        f"{len(bad)}")
            for d in sorted(pd):
                a, c = pd[d], cf.get(d, {})
                b = {rg: ba_pct.get((lab, ds, rg, d), {}) for rg in RUNG_ORDER}
                mp = num(a["pct_post"])
                if c:
                    d_mod = max(d_mod, abs(num(c["module"]) - mp))
                if b["CMRF"]:
                    d_ba = max(d_ba, abs(num(b["CMRF"]["pct_mod_post"]) - mp))
                    dd = (num(b["CMRF"]["pct_nm_post"]) -
                          num(a["pct_excl_post"]))
                    d_nm = max(d_nm, abs(dd))
                    if abs(dd) > 1e-9:
                        nm_diff.append((lab, ds, d, dd,
                                        dd * num(a["n_total"]) / 100.0))
                    if num(b["CMRF"]["n_denominator"]) != num(a["n_total"]):
                        den_ok = False
                row = [lab, ds, d, num(a["n_total"]), num(a["pct_pre"]), mp,
                       num(a["pct_excl_pre"]), num(a["pct_excl_post"]),
                       num(c.get("n_micro")), num(c.get("B1")),
                       num(c.get("B2")), num(c.get("B3")), num(c.get("B4")),
                       num(c.get("B5"))]
                row += [num(b[rg].get("pct_mod_post"))
                        for rg in RUNG_ORDER[:3]]
                row += [num(b[rg].get("pct_nm_post"))
                        for rg in RUNG_ORDER]
                row += [num(c.get(src)) for src, _ in CLIN_COLS]
                rows_a.append(row)
                donor_frame[(lab, ds, d)] = row
    sheets.append(("A_donor_frame", hdr_a, rows_a))
    ck.add("A1", ok_n and len(rows_a) == 4 * (84 + 83),
           f"{len(rows_a)} donor rows (expected 668; 84 NeoC and 83 NonNeoC "
           "per arm)")
    ck.add("A2", ok_set, "donor sets identical in Part 10h, comparator and "
           "BandAblation for every arm and region set")
    ck.add("A3", d_mod < 1e-6, f"comparator module % equals Part 10h module "
           f"post-IP % (max |diff| {d_mod:.2e})")
    ck.add("A4", d_ba < 1e-6 and den_ok, f"BandAblation CMRF boundary module "
           f"post-IP % equals Part 10h (max |diff| {d_ba:.2e}); denominators "
           f"{'equal' if den_ok else 'DIFFER'}")
    ck.add("A5", d_nm < 1e-6, f"BandAblation CMRF boundary non-module "
           f"post-IP % vs Part 10h non-module post-IP % (max |diff| "
           f"{d_nm:.2e} pp; {len(nm_diff)} donor row(s) differ)",
           warn_only=True)
    for lab, ds, d, dd, dc in nm_diff:
        ck.info(f"{lab} {ds} {d}: {dd:+.4f} pp = {dc:+.2f} cell(s)")

    # ---- B. Figure 9A ------------------------------------------------------
    print("\nB. FIGURE 9A")
    hdr_b = ["arm", "region_set", "population", "outcome", "n_donors",
             "beta_std", "se", "p_value", "FDR", "delta_R2",
             "R2_covariates_only", "R2_full", "VIF_population", "max_VIF"]
    rows_b = []
    for lab in ARM_ORDER:
        for r in p10h[lab]["Multivariate"]:
            if (r["tier"] == "all" and r["covariate_set"] == "cov4" and
                    r["predictor"] in POP9A and r["outcome"] in OUT9):
                rows_b.append([lab, r["dataset"], POP9A[r["predictor"]],
                               r["outcome"], num(r["n_donors"]),
                               num(r["beta_std"]), num(r["se"]),
                               num(r["p_value"]), num(r["p_adj"]),
                               num(r["delta_R2"]),
                               num(r["R2_covariates_only"]),
                               num(r["R2_full"]), num(r["arm_VIF"]),
                               num(r["max_VIF"])])
    rows_b.sort(key=lambda x: (ARM_ORDER.index(x[0]), DS_ORDER.index(x[1]),
                               POP_ORDER.index(x[2]), OUT9.index(x[3])))
    sheets.append(("B_Fig9A_variance", hdr_b, rows_b))
    ck.add("B1", len(rows_b) == 288, f"{len(rows_b)} rows (expected 288 = "
           "4 arms x 2 region sets x 4 populations x 9 outcomes)")
    mism = 0
    for lab in ARM_ORDER:
        mg = {(r["dataset"], r["outcome"]): r for r in comp[lab]["Marginal"]
              if r["block"] == "cov4" and r["key"] == "module"}
        for x in rows_b:
            if x[0] != lab or x[2] != "Module post-IP":
                continue
            m = mg.get((x[1], x[3]))
            if (m is None or num(m["n"]) != x[4] or
                    not near(m["beta"], x[5], TOL_ROUND4) or
                    not near(m["delta_r2"], x[9], TOL_ROUND4)):
                mism += 1
    ck.add("B2", mism == 0, f"Module post-IP rows agree with the comparator "
           f"cov4 module rows on n, beta and delta R2 ({mism} mismatch(es))")
    sig = sum(1 for x in rows_b if x[8] is not None and x[8] < 0.05)
    ck.info(f"rows with FDR < 0.05: {sig} of {len(rows_b)}")

    # ---- C. Figure 9B ------------------------------------------------------
    print("\nC. FIGURE 9B")
    hdr_c = ["arm", "region_set", "outcome", "block", "pathology_added",
             "n_donors", "beta_std", "se", "CI_low", "CI_high", "p_value",
             "FDR", "delta_R2", "VIF_module", "max_VIF"]
    rows_c = []
    blocks = list(BLOCK9B)
    for lab in ARM_ORDER:
        for r in cvp[lab]["Models"]:
            if r["key"] == "module" and r["outcome"] in COG2:
                b, se = num(r["beta"]), num(r["se"])
                rows_c.append([lab, r["dataset"], r["outcome"], r["block"],
                               BLOCK9B[r["block"]], num(r["n"]), b, se,
                               b - 1.96 * se, b + 1.96 * se,
                               num(r["p_value"]), num(r["p_adj"]),
                               num(r["delta_r2"]), num(r["vif_module"]),
                               num(r["max_vif"])])
    rows_c.sort(key=lambda x: (ARM_ORDER.index(x[0]), DS_ORDER.index(x[1]),
                               COG2.index(x[2]), blocks.index(x[3])))
    sheets.append(("C_Fig9B_cognition", hdr_c, rows_c))
    ck.add("C1", len(rows_c) == 64, f"{len(rows_c)} rows (expected 64)")
    mism = 0
    for lab in ARM_ORDER:
        mvr = {(r["dataset"], r["outcome"]): r
               for r in p10h[lab]["Multivariate"]
               if r["tier"] == "all" and r["covariate_set"] == "cov4"
               and r["predictor"] == "pct_post"}
        for x in rows_c:
            if x[0] == lab and x[3] == "cov4":
                m = mvr.get((x[1], x[2]))
                if (m is None or num(m["n_donors"]) != x[5] or
                        not near(m["beta_std"], x[6], TOL_ROUND4)):
                    mism += 1
    ck.add("C2", mism == 0, "cov4 rows agree with the Part 10h Multivariate "
           f"sheet on n and beta ({mism} mismatch(es))")
    for lab in ARM_ORDER:
        s = sum(1 for x in rows_c if x[0] == lab and x[11] < 0.05)
        ck.info(f"{lab}: {s} of 16 rows with FDR < 0.05")

    # ---- D. Figure 9C ------------------------------------------------------
    print("\nD. FIGURE 9C")
    hdr_d = ["arm", "region_set", "outcome", "block", "cognition_added",
             "n_donors", "beta_std", "se", "CI_low", "CI_high", "p_value",
             "FDR", "delta_R2", "VIF_module", "max_VIF"]
    rows_d = []
    for lab in ARM_ORDER:
        for r in cvp[lab]["Reverse"]:
            if r["key"] == "module" and r["outcome"] in PATH4:
                b, se = num(r["beta"]), num(r["se"])
                blab, added = BLOCK9C[r["block"]]
                rows_d.append([lab, r["dataset"], r["outcome"], blab, added,
                               num(r["n"]), b, se, b - 1.96 * se,
                               b + 1.96 * se, num(r["p_value"]),
                               num(r["p_adj"]), num(r["delta_r2"]),
                               num(r["vif_module"]), num(r["max_vif"])])
    bl9c = [v[0] for v in BLOCK9C.values()]
    rows_d.sort(key=lambda x: (ARM_ORDER.index(x[0]), DS_ORDER.index(x[1]),
                               PATH4.index(x[2]), bl9c.index(x[3])))
    sheets.append(("D_Fig9C_pathology", hdr_d, rows_d))
    ck.add("D1", len(rows_d) == 64, f"{len(rows_d)} rows (expected 64)")
    got = {(x[0], x[1], x[2]) for x in rows_d
           if x[3] == "+MMSE + CASI" and x[11] < 0.05}
    ck.add("D2", got == EXP_9C_SIG, "rows with FDR < 0.05 after adding MMSE "
           "+ CASI are ADNC, CERAD and Thal in M4 and PVM_3 NeoC only (found "
           f"{len(got)}: " + "; ".join(" ".join(g) for g in sorted(got)) + ")")
    mism = 0
    for lab in ARM_ORDER:
        mvr = {(r["dataset"], r["outcome"]): r
               for r in p10h[lab]["Multivariate"]
               if r["tier"] == "all" and r["covariate_set"] == "cov4"
               and r["predictor"] == "pct_post"}
        for x in rows_d:
            if x[0] == lab and x[3] == "cov4":
                m = mvr.get((x[1], x[2]))
                if (m is None or num(m["n_donors"]) != x[5] or
                        not near(m["beta_std"], x[6], TOL_ROUND4)):
                    mism += 1
    ck.add("D3", mism == 0, "cov4 rows agree with the Part 10h Multivariate "
           f"sheet on n and beta ({mism} mismatch(es))")

    # ---- E. Figure 9D ------------------------------------------------------
    print("\nE. FIGURE 9D")
    hdr_e = ["arm", "region_set", "outcome", "domain", "n_donors",
             "beta_module", "beta_nonmodule", "difference", "se", "CI_low",
             "CI_high", "p_value", "FDR", "boot_CI_low", "boot_CI_high",
             "MDES_80pct", "verdict", "p_module", "FDR_module",
             "p_nonmodule", "FDR_nonmodule"]
    rows_e = []
    o6 = list(OUT11D.values())
    for r in p11["Contrast"]:
        rows_e.append([ARM11D[r["arm"]], r["dataset"], OUT11D[r["outcome"]],
                       r["domain"], num(r["n_donors"]), num(r["b_module"]),
                       num(r["b_nonmodule"]), num(r["delta"]), num(r["se"]),
                       num(r["exact_lo"]), num(r["exact_hi"]),
                       num(r["p_exact"]), num(r["p_bh"]), num(r["boot_lo"]),
                       num(r["boot_hi"]), num(r["mdes80"]), r["verdict"]])
    rows_e.sort(key=lambda x: (ARM_ORDER.index(x[0]), DS_ORDER.index(x[1]),
                               o6.index(x[2])))
    # Each coefficient of the joint model tested on its own, refitted from
    # sheet A: standardized outcome (MMSE and CASI x -1) on standardized
    # cov4 + Module post-IP % + non-module post-IP %, donors complete on all
    # six outcomes; BH within arm, region set and predictor over 6 outcomes.
    col = {h: i for i, h in enumerate(hdr_a)}
    y_col = {"MMSE score": "MMSE", "CASI score": "CASI", "ADNC rank": "ADNC",
             "CERAD score": "CERAD", "Braak stage": "Braak",
             "Thal phase": "Thal"}
    sign = {"MMSE score": -1.0, "CASI score": -1.0}
    x_cols = ["age_at_death", "sex", "APOE_e4_dose", "PMI",
              "pct_module_postIP", "pct_nonmodule_postIP"]
    need = list(y_col.values()) + x_cols
    fit = {}
    for lab in ARM_ORDER:
        for ds in DS_ORDER:
            dn = [v for (a, s, _), v in sorted(donor_frame.items())
                  if a == lab and s == ds
                  and all(v[col[c]] is not None for c in need)]
            X = np.array([[v[col[c]] for c in x_cols] for v in dn], float)
            X = (X - X.mean(0)) / X.std(0, ddof=1)
            X = np.column_stack([np.ones(len(dn)), X])
            xtx_inv = np.linalg.inv(X.T @ X)
            dfr = len(dn) - X.shape[1]
            for o, yc in y_col.items():
                y = np.array([v[col[yc]] for v in dn], float)
                y = sign.get(o, 1.0) * (y - y.mean()) / y.std(ddof=1)
                b = xtx_inv @ X.T @ y
                res = y - X @ b
                se = np.sqrt(np.diag(xtx_inv) * (res @ res) / dfr)
                out = []
                for k in (5, 6):                      # module, non-module
                    t = b[k] / se[k]
                    out += [float(b[k]),
                            _betai(dfr / 2.0, 0.5, dfr / (dfr + t * t))]
                fit[(lab, ds, o)] = (len(dn), out)
    for lab in ARM_ORDER:
        for ds in DS_ORDER:
            for j in (1, 3):                          # p_module, p_nonmodule
                keys = [(lab, ds, o) for o in o6]
                q = bh([fit[k][1][j] for k in keys])
                for k, v in zip(keys, q):
                    fit[k][1].append(v)
    mism = 0
    for x in rows_e:
        n, (bm, pm, bn, pn, qm, qn) = fit[(x[0], x[1], x[2])]
        if (n != x[4] or not near(bm, x[5], TOL_ROUND4)
                or not near(bn, x[6], TOL_ROUND4)):
            mism += 1
        x += [pm, qm, pn, qn]
    sheets.append(("E_Fig9D_mod_vs_nonmod", hdr_e, rows_e))
    ck.add("E1", len(rows_e) == 48, f"{len(rows_e)} rows (expected 48)")
    sig = [x for x in rows_e if x[12] < 0.05]
    ok = (len(sig) == 6 and all(x[3] == "cognition" for x in sig)
          and all(x[0] != "HKG8" for x in sig))
    ck.add("E2", ok, f"{len(sig)} rows with FDR < 0.05 (expected 6, all "
           "cognition, none HKG8): " +
           "; ".join(f"{x[0]} {x[1]} {x[2]}" for x in sig))
    ck.add("E3", mism == 0, "refitted module and non-module coefficients "
           f"agree with the Contrast sheet on n and beta ({mism} "
           "mismatch(es))")
    main3 = [x for x in rows_e if x[0] != "HKG8"]
    m_cog = sum(1 for x in main3 if x[3] == "cognition" and x[18] < 0.05)
    m_pat = sum(1 for x in main3 if x[3] == "pathology" and x[18] < 0.05)
    n_any = sum(1 for x in rows_e if x[20] < 0.05)
    ck.add("E4", (m_cog, m_pat, n_any) == (12, 22, 0),
           f"main arms, FDR_module < 0.05: cognition {m_cog} of 12 "
           f"(expected 12), pathology {m_pat} of 24 (expected 22); "
           f"FDR_nonmodule < 0.05 in {n_any} of 48 rows (expected 0)")
    ck.info("smallest FDR_nonmodule: "
            f"{min(x[20] for x in rows_e):.3f}")

    # ---- F. Figure 10A -----------------------------------------------------
    print("\nF. FIGURE 10A")
    hdr_f = ["arm", "region_set", "predictor", "predictor_abbrev", "outcome",
             "n_donors", "beta_std", "se", "p_value", "FDR", "delta_R2",
             "R2_covariates_only", "VIF_predictor"]
    rows_f = []
    for lab in ARM_ORDER:
        for r in comp[lab]["Marginal"]:
            if r["block"] == "cov4" and r["outcome"] in OUT9:
                rows_f.append([lab, r["dataset"], r["predictor"],
                               PRED10_ABBR[r["predictor"]], r["outcome"],
                               num(r["n"]), num(r["beta"]), num(r["se"]),
                               num(r["p_value"]), num(r["p_adj"]),
                               num(r["delta_r2"]),
                               num(r["r2_covariates_only"]),
                               num(r["vif_pred"])])
    rows_f.sort(key=lambda x: (ARM_ORDER.index(x[0]), DS_ORDER.index(x[1]),
                               PRED10_ORDER.index(x[2]), OUT9.index(x[4])))
    sheets.append(("F_Fig10A_benchmark", hdr_f, rows_f))
    ck.add("F1", len(rows_f) == 432, f"{len(rows_f)} rows (expected 432)")
    pl = {(GROUP13[r["group"]], r["dataset"], r["predictor"], r["outcome"]): r
          for r in read_csv(P["f13A"])}
    mism = sum(1 for x in rows_f
               if (x[0], x[1], x[2], x[4]) not in pl
               or num(pl[(x[0], x[1], x[2], x[4])]["n"]) != x[5]
               or not near(pl[(x[0], x[1], x[2], x[4])]["beta"], x[6], 1e-9)
               or not near(pl[(x[0], x[1], x[2], x[4])]["delta_r2"], x[10],
                           1e-9)
               or not near(pl[(x[0], x[1], x[2], x[4])]["p_adj"], x[9],
                           1e-12))
    ck.add("F2", mism == 0 and len(pl) == len(rows_f),
           f"all rows agree with Fig13A_DeltaR2_grid_plotted.csv "
           f"({len(pl)} plotted rows; {mism} mismatch(es))")

    # ---- G. Figure 10B -----------------------------------------------------
    print("\nG. FIGURE 10B")
    hdr_g = ["arm", "region_set", "competitor", "competitor_abbrev",
             "domain", "outcome", "n_donors", "beta_module",
             "beta_competitor", "p_module", "p_competitor", "FDR_module",
             "FDR_competitor", "R2_module_only", "R2_competitor_only",
             "R2_joint", "VIF_joint", "readable", "verdict", "verdict_shown"]
    rows_g = []
    out8 = COG3 + PATH5
    for lab in ARM_ORDER:
        for r in comp[lab]["Duels"]:
            if r["outcome"] not in out8:
                continue
            vif = num(r["vif_joint"])
            readable = vif is not None and vif < 5
            rows_g.append([lab, r["dataset"], r["rival_label"],
                           PRED10_ABBR[r["rival_label"]],
                           "Cognition" if r["outcome"] in COG3
                           else "Pathology", r["outcome"], num(r["n"]),
                           num(r["beta_module"]), num(r["beta_rival"]),
                           num(r["p_module"]), num(r["p_rival"]),
                           num(r["p_module_adj"]), num(r["p_rival_adj"]),
                           num(r["r2_module_only"]),
                           num(r["r2_rival_only"]), num(r["r2_joint"]), vif,
                           readable, r["verdict"],
                           r["verdict"] if readable else "not readable"])
    rows_g.sort(key=lambda x: (ARM_ORDER.index(x[0]), DS_ORDER.index(x[1]),
                               PRED10_ORDER.index(x[2]), out8.index(x[5])))
    sheets.append(("G_Fig10B_competitors", hdr_g, rows_g))
    ck.add("G1", len(rows_g) == 320, f"{len(rows_g)} rows (expected 320)")
    pl = {(GROUP13[r["group"]], r["dataset"], r["rival_label"],
           r["outcome"]): r for r in read_csv(P["f13B"])}
    mism = 0
    for x in rows_g:
        q = pl.get((x[0], x[1], x[2], x[5]))
        if (q is None or num(q["n"]) != x[6]
                or not near(q["beta_module"], x[7], 1e-9)
                or not near(q["beta_rival"], x[8], 1e-9)
                or not near(q["p_module_adj"], x[11], 1e-12)
                or not near(q["p_rival_adj"], x[12], 1e-12)
                or not near(q["vif_joint"], x[16], 1e-9)
                or q["verdict"] != x[18]
                or tb(q["unreadable"]) == x[17]
                or q["verdict_disp"] != x[19]):
            mism += 1
    ck.add("G2", mism == 0 and len(pl) == len(rows_g),
           f"all rows agree with Fig13B_Duel_verdicts_plotted.csv on "
           f"coefficients, FDRs, VIF and verdict ({mism} mismatch(es))")
    ck.info(f"not readable (joint VIF >= 5): "
            f"{sum(1 for x in rows_g if not x[17])} of {len(rows_g)}")

    # ---- H. Figure 10C and 10D ----------------------------------------------
    print("\nH. FIGURES 10C AND 10D")
    hdr_h = ["arm", "region_set", "block", "outcome", "score",
             "figure_panel", "n_donors", "rs_module_vs_score",
             "beta_module", "p_module", "FDR_module", "beta_score",
             "p_score", "beta_cluster", "p_cluster", "max_VIF",
             "delta_R2_module", "F_1df", "verdict",
             "beta_module_without_CP", "p_module_without_CP",
             "FDR_module_without_CP", "VIF_module_without_CP"]
    rows_h = []
    vmap = {"HOLDS": "HOLDS", "does not hold": "doesn't hold",
            "not readable": "unreadable"}
    for lab in ARM_ORDER:
        rs = {}
        for r in mv[lab]["Correlations"]:
            if str(r["x"]).startswith("CMRF module"):
                for s in SCORE_ORDER:
                    if str(r["y"]).startswith(f"Module score ({s})"):
                        rs[(r["dataset"], s)] = num(r["spearman"])
        for r in mv[lab]["Decomposition"]:
            s = r["score_arm"]
            panels = []
            if s == "AUCell":
                panels.append("10C")
            if lab == "M4":
                panels.append("10D")
            if not panels:
                continue
            rows_h.append([lab, r["dataset"], r["block"], r["outcome"], s,
                           ", ".join(panels), num(r["n"]),
                           rs.get((r["dataset"], s)), num(r["beta_module"]),
                           num(r["p_module_t"]), num(r["q_module_t"]),
                           num(r["beta_score"]), num(r["p_score"]),
                           num(r["beta_cluster"]), num(r["p_cluster"]),
                           num(r["max_vif"]),
                           num(r["dR2_module_over_rival"]), num(r["F_1df"]),
                           vmap.get(r["verdict"], r["verdict"]),
                           num(r["beta_module_vs_score"]),
                           num(r["p_module_vs_score"]),
                           num(r["q_module_vs_score"]),
                           num(r["vif_module_vs_score"])])
    bl2 = ["cov4", "+BTC"]
    rows_h.sort(key=lambda x: (ARM_ORDER.index(x[0]), DS_ORDER.index(x[1]),
                               bl2.index(x[2]), COG2.index(x[3]),
                               SCORE_ORDER.index(x[4])))
    sheets.append(("H_Fig10CD_decomposition", hdr_h, rows_h))
    ck.add("H1", len(rows_h) == 48, f"{len(rows_h)} rows (expected 48: 32 "
           "for 10C, 24 for 10D, 8 shared)")
    idx = {(x[0], x[1], x[2], x[3], x[4]): x for x in rows_h}
    order4 = [("cov4", "MMSE score"), ("cov4", "CASI score"),
              ("+BTC", "MMSE score"), ("+BTC", "CASI score")]
    bad = []
    for (ds, lab), vals in EXP_10C.items():
        for (bk, oc), (q, v, vd) in zip(order4, vals):
            x = idx.get((lab, ds, bk, oc, "AUCell"))
            if (x is None or not near(x[10], q, TOL3) or
                    not near(x[15], v, TOL2) or x[18] != vd):
                bad.append(f"{lab} {ds} {bk} {oc}")
    ck.add("H2", not bad, "10C FDR, VIF and verdict match the figure for "
           f"32 cells ({len(bad)} differ" + (": " + "; ".join(bad)
                                             if bad else "") + ")")
    bad = [f"{a} {d}" for (a, d), v in EXP_RS.items()
           if not near(idx.get((a, d, "cov4", "MMSE score", "AUCell"),
                               [None] * 8)[7], v, TOL3)]
    ck.add("H3", not bad, "10C module vs AUCell rs match the figure for 8 "
           f"arm and region-set pairs ({len(bad)} differ)")
    bad = []
    for (ds, bk, oc), vals in EXP_10D.items():
        for s, (bm, bs, q) in zip(SCORE_ORDER, vals):
            x = idx.get(("M4", ds, bk, oc, s))
            if (x is None or not near(x[8], bm, TOL3) or
                    not near(x[11], bs, TOL3) or not near(x[21], q, TOL3)):
                bad.append(f"{ds} {bk} {oc} {s}")
    ck.add("H4", not bad, "10D betas (module, score) and FDR match the "
           f"figure for 24 cells ({len(bad)} differ)")
    mism = 0
    for r in read_csv(P["f13C"]):
        x = idx.get((GROUP13[r["group"]], r["dataset"], r["block"],
                     r["outcome"], r["score_arm"]))
        if (x is None or not near(r["beta_module"], x[8], 1e-9)
                or not near(r["q_module_t"], x[10], 1e-9)
                or not near(r["max_vif"], x[15], 1e-9)):
            mism += 1
    for r in read_csv(P["f13C2"]):
        x = idx.get((GROUP13[r["group"]], r["dataset"], r["block"],
                     r["outcome"], r["score_arm"]))
        if (x is None or not near(r["beta_score"], x[11], 1e-9)
                or not near(r["q_module_vs_score"], x[21], 1e-9)
                or not near(r["beta_module_vs_score"], x[19], 1e-9)):
            mism += 1
    ck.add("H5", mism == 0, "rows agree with Fig13C_Pub_Decomposition and "
           f"Fig13C2_Pub_Contrast plotted CSVs ({mism} mismatch(es))")

    # ---- I. Figure 11A -----------------------------------------------------
    print("\nI. FIGURE 11A")
    hdr_i = ["arm", "region_set", "boundary", "outcome", "block", "plotted",
             "n_donors", "beta_std", "p_value", "FDR", "delta_R2"]
    rows_i = []
    for r in ba["Donor_Assoc"]:
        blk = "+BTC" if r["block"] == "BTC" else r["block"]
        oc = "MMSE score" if r["outcome"] == "MMSE" else "CASI score"
        rows_i.append([GROUP11[r["group"]], r["dataset"], RUNG[r["rung"]],
                       oc, blk, blk == "+BTC", num(r["n"]), num(r["beta"]),
                       num(r["p_value"]), num(r["FDR"]), num(r["delta_r2"])])
    rows_i.sort(key=lambda x: (ARM_ORDER.index(x[0]), DS_ORDER.index(x[1]),
                               ["CASI score", "MMSE score"].index(x[3]),
                               RUNG_ORDER.index(x[2]), bl2.index(x[4])))
    sheets.append(("I_Fig11A_boundary_assoc", hdr_i, rows_i))
    ck.add("I1", len(rows_i) == 128, f"{len(rows_i)} rows (expected 128)")
    btc = [x for x in rows_i if x[4] == "+BTC"]
    s_main = sum(1 for x in btc if x[0] != "HKG8" and x[9] < 0.05)
    s_hkg = sum(1 for x in btc if x[0] == "HKG8" and x[9] < 0.05)
    ck.add("I2", s_main == 48 and s_hkg == 0, f"+BTC rows with FDR < 0.05: "
           f"{s_main} of 48 on M4, PVM_3, PVM_1 and {s_hkg} of 16 on HKG8 "
           "(figure: 48 and 0)")
    mism = 0
    for x in rows_i:
        if x[2] != "CMRF":
            continue
        m = [r for r in cvp[x[0]]["Models"] if r["key"] == "module"
             and r["dataset"] == x[1] and r["outcome"] == x[3]
             and r["block"] == x[4]]
        if (len(m) != 1 or num(m[0]["n"]) != x[6] or
                not near(m[0]["beta"], x[7], 1e-4)):
            mism += 1
    ck.add("I3", mism == 0, "CMRF boundary rows reproduce the stage 2 module "
           f"rows (cov4 and +BTC) on n and beta ({mism} mismatch(es))")

    # ---- J, K. Figure 11B and cell-level boundary membership ---------------
    print("\nJ-K. FIGURE 11B AND CELL-LEVEL BOUNDARIES")
    cov_run = {(GROUP11[r["group"]], r["dataset"]): r for r in ba["Coverage"]}
    hdr_jc = ["arm", "region_set", "boundary", "u", "center_line", "lower",
              "upper"]
    rows_jc = []
    for r in read_csv(P["ba_curves"]):
        rows_jc.append([GROUP11[r["group"]], r["dataset"], RUNG[r["rung"]],
                        num(r["u"]), num(r["fit"]), num(r["lo"]),
                        num(r["hi"])])
    rows_jc.sort(key=lambda x: (RUNG_ORDER.index(x[2]), x[3]))
    sheets.append(("J_Fig11B_boundary_curves", hdr_jc, rows_jc))
    got = sorted({x[2] for x in rows_jc}, key=RUNG_ORDER.index)
    sets_jc = {(x[0], x[1]) for x in rows_jc}
    ck.add("J1", got == RUNG_ORDER and sets_jc == {("M4", "NeoC")},
           f"{len(rows_jc)} curve rows for M4 NeoC, boundaries: "
           + ", ".join(got))
    cv = cov_run[("M4", "NeoC")]
    dw = 0.0
    for x in rows_jc:
        w = {"OLS": num(cv["width_A0"]),
             "Fixed spline": num(cv["width_A1"])}.get(x[2])
        if w is not None:
            dw = max(dw, abs((x[6] - x[4]) - w), abs((x[4] - x[5]) - w))
    ck.add("J2", dw < 1e-9, "OLS and fixed-spline half-widths equal the "
           f"run's boundary widths (max |diff| {dw:.2e})")

    hdr_kc = ["arm", "region_set", "barcode", "donor_id", "cell_class",
              "cDFG", "side", "u", "v", "fit_OLS", "fit_spline",
              "half_width_local_spline", "in_OLS", "in_fixed_spline",
              "in_local_spline", "in_CMRF", "in_Fig11B"]
    rows_kc = []
    for r in read_csv(P["ba_cells"]):
        arm = GROUP11[r["group"]]
        dn = r["donor"] if r["donor"] not in ("", "NA") else None
        rows_kc.append([arm, r["dataset"], r["barcode"], dn,
                        r["cell_class"], num(r["cDFG"]), r["side"],
                        num(r["u"]), num(r["v"]), num(r["fit_OLS"]),
                        num(r["fit_spline"]), num(r["width_local"]),
                        tb(r["in_OLS"]), tb(r["in_fixed_spline"]),
                        tb(r["in_local_spline"]), tb(r["in_CMRF"]),
                        arm == "M4" and r["dataset"] == "NeoC"])
    rows_kc.sort(key=lambda x: (ARM_ORDER.index(x[0]), DS_ORDER.index(x[1]),
                                x[5]))
    sheets.append(("K_Fig11BC_cells", hdr_kc, rows_kc))
    mcol = {"OLS": 12, "Fixed spline": 13, "Local spline": 14, "CMRF": 15}
    # K1: counts per set against the run's Coverage sheet.
    bad = []
    by_set = {}
    for x in rows_kc:
        by_set.setdefault((x[0], x[1]), []).append(x)
    for key, cv in cov_run.items():
        z = by_set.get(key, [])
        have = [len(z)] + [sum(1 for x in z if x[mcol[rg]])
                           for rg in RUNG_ORDER]
        want = [num(cv[k]) for k in ("n_cells", "n_A0", "n_A1", "n_A2",
                                     "n_A3")]
        if [float(h) for h in have] != want:
            bad.append(f"{key[0]} {key[1]} {have} vs {want}")
    ck.add("K1", not bad and len(by_set) == 8, "cell and boundary counts "
           f"equal the run's Coverage sheet in all 8 sets ({len(bad)} differ"
           + (": " + "; ".join(bad) if bad else "") + ")")
    # K2: Jaccard recomputed from membership.
    dj = 0.0
    for r in ba["Overlap"]:
        z = by_set.get((GROUP11[r["group"]], r["dataset"]), [])
        a = {x[2] for x in z if x[mcol[RUNG[r["a"]]]]}
        b = {x[2] for x in z if x[mcol[RUNG[r["b"]]]]}
        j = len(a & b) / len(a | b) if (a | b) else None
        dj = max(dj, abs(j - num(r["jaccard"])))
    ck.add("K2", dj < 1e-9, "Jaccard overlaps recomputed from the membership "
           f"columns equal the run's Overlap sheet (max |diff| {dj:.2e})")
    # K3: donor percentages recomputed from the cells.
    ci = {c: hdr_a.index(c) for c in hdr_a}
    pcols = {("OLS", True): "pct_module_postIP_OLS",
             ("Fixed spline", True): "pct_module_postIP_fixed_spline",
             ("Local spline", True): "pct_module_postIP_local_spline",
             ("CMRF", True): "pct_module_postIP",
             ("OLS", False): "pct_nonmodule_postIP_OLS",
             ("Fixed spline", False): "pct_nonmodule_postIP_fixed_spline",
             ("Local spline", False): "pct_nonmodule_postIP_local_spline",
             ("CMRF", False): "pct_nonmodule_postIP_CMRF_boundary"}
    cnt = {}
    for x in rows_kc:
        if x[3] is None or x[6] != "Post":
            continue
        for rg in RUNG_ORDER:
            k = (x[0], x[1], x[3], rg, bool(x[mcol[rg]]))
            cnt[k] = cnt.get(k, 0) + 1
    dp = 0.0
    for (lab, ds, d), row in donor_frame.items():
        nt = row[ci["n_supertype_cells"]]
        for (rg, mod), col in pcols.items():
            pct = 100.0 * cnt.get((lab, ds, d, rg, mod), 0) / nt
            dp = max(dp, abs(pct - row[ci[col]]))
    ck.add("K3", dp < 1e-9, "donor module and non-module post-IP % for all "
           "four boundaries recomputed from the cells equal sheet A (max "
           f"|diff| {dp:.2e} pp)")
    summ = {(GROUP11[r["group"]], r["dataset"]): r
            for r in read_csv(P["ba_sum"])}
    di = max(abs(num(summ[k]["ip_y"]) - num(cov_run[k]["ip_y"]))
             for k in cov_run)
    ck.add("K4", di < 1e-9, "export IP_y equals the run's Coverage IP_y "
           f"(max |diff| {di:.2e})")
    z = by_set.get(("M4", "NeoC"), [])
    ck.info(f"Figure 11B cells (M4 NeoC): {len(z)}, CMRF module "
            f"{sum(1 for x in z if x[15])}; cells without a donor: "
            f"{sum(1 for x in rows_kc if x[3] is None)}")

    # ---- L, M, N, O. Figure 11C --------------------------------------------
    print("\nL-O. FIGURE 11C")
    sets = {"three real groups": ["M4", "PVM_3", "PVM_1"],
            "HKG8 control": ["HKG8"]}
    ov = {}
    for r in ba["Overlap"]:
        if r["b"] == "A3_CMRF":
            ov.setdefault(RUNG[r["a"]], []).append(
                (GROUP11[r["group"]], num(r["jaccard"])))
    hdr_j = ["set", "boundary", "overlap_with_CMRF_median",
             "overlap_with_CMRF_min", "overlap_with_CMRF_max",
             "clinical_tests_hit", "clinical_tests_of",
             "negative_control_tests_hit", "negative_control_tests_of",
             "discarded_cell_tests_hit", "discarded_cell_tests_of",
             "cognition_tests_surviving_BTC", "cognition_tests_of",
             "median_abs_beta_discarded_cells"]
    rows_j = []
    bad = []
    for r in ba["Scorecard"]:
        st, rg = r["set"], RUNG[r["rung"]]
        js = [j for g, j in ov.get(rg, []) if g in sets[st]]
        mean_ = statistics.median(js) if js else None
        mn, mx = (min(js), max(js)) if js else (None, None)
        row = [st, rg, mean_, mn, mx, num(r["C1_clinical_hits"]),
               num(r["C1_of"]), num(r["C2_control_hits"]), num(r["C2_of"]),
               num(r["C3_nonmodule_hits"]), num(r["C3_of"]),
               num(r["C4_cog_BTC"]), num(r["C4_of"]),
               num(r["median_abs_beta_nonmodule"])]
        rows_j.append(row)
        e = EXP_11C.get((st, rg))
        of = EXP_11C_OF[st]
        okr = (e is not None and
               all(near(a, b, TOL3) for a, b in zip(row[2:5], e[:3])) and
               [row[5], row[7], row[9], row[11]] == list(map(float, e[3:7]))
               and [row[6], row[8], row[10], row[12]] == list(map(float, of))
               and near(row[13], e[7], TOL3))
        if not okr:
            bad.append(f"{st} {rg}")
    rows_j.sort(key=lambda x: (list(sets).index(x[0]),
                               RUNG_ORDER.index(x[1])))
    sheets.append(("L_Fig11C_scorecard", hdr_j, rows_j))
    ck.add("L1", len(rows_j) == 8 and not bad, "8 scorecard rows match the "
           f"Figure 11C table ({len(bad)} differ" +
           (": " + "; ".join(bad) if bad else "") + ")")

    hdr_k = ["arm", "region_set", "boundary", "model", "predictor", "panel",
             "outcome", "block", "n_donors", "beta_std", "p_value", "FDR",
             "significant", "delta_R2", "scorecard_count",
             "in_median_abs_beta"]
    pred_lab = {"pct_mod_post": "Module post-IP %",
                "pct_nm_post": "Non-module post-IP %"}
    rows_k = []
    for r in ba["Full_Panel"]:
        rows_k.append([GROUP11[r["group"]], r["dataset"], RUNG[r["rung"]],
                       r["model"], pred_lab.get(r["predictor"],
                                                r["predictor"]),
                       r["panel"], r["outcome"],
                       "+BTC" if r["block"] == "BTC" else r["block"],
                       num(r["n"]), num(r["beta"]), num(r["p_value"]),
                       num(r["FDR"]), tb(r["significant"]),
                       num(r["delta_r2"]), "", False])
    # Scorecard columns recomputed from the test rows.
    def col_of(x):
        clin = x[5] in ("cognition", "pathology")
        if x[3] == "module_only" and x[7] == "cov4" and clin:
            return "C1"
        if x[3] == "module_only" and x[7] == "cov4" and \
                x[5] == "negative control":
            return "C2"
        if x[3] == "encompassing" and x[4].startswith("Non-module") and \
                x[7] == "cov4" and clin:
            return "C3"
        if x[3] == "module_only" and x[7] == "+BTC" and x[5] == "cognition":
            return "C4"
        return ""
    def in_med(x):
        return (x[3] == "encompassing" and x[4].startswith("Non-module")
                and x[5] in ("cognition", "pathology"))
    rec_ok = True
    for st, arms in sets.items():
        for rg in RUNG_ORDER:
            sub = [x for x in rows_k if x[0] in arms and x[2] == rg]
            cnt = {c: [0, 0] for c in ("C1", "C2", "C3", "C4")}
            nm_beta = []
            for x in sub:
                c = col_of(x)
                if c:
                    cnt[c][1] += 1
                    cnt[c][0] += int(bool(x[12]))
                if in_med(x):
                    nm_beta.append(abs(x[9]))
            j = [y for y in rows_j if y[0] == st and y[1] == rg][0]
            want = [(j[5], j[6]), (j[7], j[8]), (j[9], j[10]),
                    (j[11], j[12])]
            have = [tuple(map(float, cnt[c])) for c in ("C1", "C2", "C3",
                                                        "C4")]
            med = statistics.median(nm_beta) if nm_beta else None
            if have != want or not near(med, j[13], 1e-9):
                rec_ok = False
                ck.info(f"recompute {st} {rg}: counts {have} vs {want}; "
                        f"median {med} vs {j[13]}")
    if rec_ok:
        for x in rows_k:
            x[14] = col_of(x)
            x[15] = in_med(x)
    ck.add("M1", rec_ok, "scorecard counts and median |beta| recomputed from "
           "the test rows (C1 module-only model, cov4, clinical; C2 "
           "module-only, cov4, negative control; C3 non-module term of the "
           "encompassing model, cov4, clinical; C4 module-only, +BTC, "
           "cognition; median |beta| over the non-module term of the "
           "encompassing model, cov4 and +BTC, clinical)")
    rows_k.sort(key=lambda x: (ARM_ORDER.index(x[0]), DS_ORDER.index(x[1]),
                               RUNG_ORDER.index(x[2]), x[3], x[5], x[6],
                               x[7], x[4]))
    sheets.append(("M_Fig11C_tests", hdr_k, rows_k))
    ck.add("M2", len(rows_k) == 2304, f"{len(rows_k)} rows (expected 2304)")

    hdr_l = ["arm", "region_set", "boundary_a", "boundary_b", "n_cells_a",
             "n_cells_b", "jaccard"]
    rows_l = [[GROUP11[r["group"]], r["dataset"], RUNG[r["a"]],
               RUNG[r["b"]], num(r["n_a"]), num(r["n_b"]),
               num(r["jaccard"])] for r in ba["Overlap"]]
    rows_l.sort(key=lambda x: (ARM_ORDER.index(x[0]), DS_ORDER.index(x[1]),
                               RUNG_ORDER.index(x[2]),
                               RUNG_ORDER.index(x[3])))
    sheets.append(("N_Fig11C_overlap", hdr_l, rows_l))
    ck.add("N1", len(rows_l) == 48, f"{len(rows_l)} rows (expected 48)")

    hdr_m = ["arm", "region_set", "n_cells", "IP_y", "IP_y_rule",
             "n_cells_without_donor", "target_coverage",
             "width_OLS", "width_fixed_spline", "width_local_spline_min",
             "width_local_spline_max", "n_module_cells_OLS",
             "n_module_cells_fixed_spline", "n_module_cells_local_spline",
             "n_module_cells_CMRF"]
    rows_m = [[GROUP11[r["group"]], r["dataset"], num(r["n_cells"]),
               num(r["ip_y"]),
               summ[(GROUP11[r["group"]], r["dataset"])]["ip_y_rule"],
               num(summ[(GROUP11[r["group"]], r["dataset"])]["n_no_donor"]),
               num(r["target_coverage"]),
               num(r["width_A0"]), num(r["width_A1"]),
               num(r["width_A2_min"]), num(r["width_A2_max"]),
               num(r["n_A0"]), num(r["n_A1"]), num(r["n_A2"]),
               num(r["n_A3"])] for r in ba["Coverage"]]
    rows_m.sort(key=lambda x: (ARM_ORDER.index(x[0]), DS_ORDER.index(x[1])))
    sheets.append(("O_Fig11BC_boundary_summary", hdr_m, rows_m))
    ck.add("O1", len(rows_m) == 8, f"{len(rows_m)} rows (expected 8)")

    # ---- P. Supplementary Figure 14: cov4 against cov7 ----------------------
    print("\nP. SUPPLEMENTARY FIGURE 14 (cov4 AND cov7)")
    BLOCKS = ("cov4", "cov7")
    CTRL6 = ["Age at death", "Sex", "Brain pH", "PMI", "RIN",
             "Fresh brain weight"]
    OUT15 = OUT9 + CTRL6
    hdr_p = ["arm", "region_set", "outcome", "block", "covariates",
             "n_donors", "beta_std", "se", "p_value", "FDR", "delta_R2",
             "VIF_module", "max_VIF"]
    rows_p = []
    for lab in ARM_ORDER:
        for r in p10h[lab]["Multivariate"]:
            if (r["tier"] == "all" and r["predictor"] == "pct_post" and
                    r["covariate_set"] in BLOCKS and r["outcome"] in OUT15):
                rows_p.append([lab, r["dataset"], r["outcome"],
                               r["covariate_set"], r.get("covariates"),
                               num(r["n_donors"]), num(r["beta_std"]),
                               num(r["se"]), num(r["p_value"]),
                               num(r["p_adj"]), num(r["delta_R2"]),
                               num(r["arm_VIF"]), num(r["max_VIF"])])
    rows_p.sort(key=lambda x: (ARM_ORDER.index(x[0]), DS_ORDER.index(x[1]),
                               OUT15.index(x[2]), x[3]))
    sheets.append(("P_SuppFig14_cov4_cov7", hdr_p, rows_p))
    ck.add("P1", len(rows_p) == 240, f"{len(rows_p)} rows (expected 240 = "
           "4 arms x 2 region sets x 15 outcomes x 2 blocks)")
    idx_b = {(x[0], x[1], x[3]): x for x in rows_b
             if x[2] == "Module post-IP"}
    mism = sum(1 for x in rows_p if x[3] == "cov4" and x[2] in OUT9 and (
        not near(idx_b[(x[0], x[1], x[2])][9], x[10], 1e-12) or
        not near(idx_b[(x[0], x[1], x[2])][8], x[9], 1e-12)))
    ck.add("P2", mism == 0, "cov4 rows equal the Module post-IP rows of "
           f"sheet B ({mism} mismatch(es))")
    mism = 0
    for lab in ARM_ORDER:
        cc = {(r["dataset"], r["outcome"]): r
              for r in p10h[lab]["Cov4_vs_Cov7"]
              if r["predictor"] == "pct_post"}
        for x in rows_p:
            if x[0] != lab:
                continue
            r = cc.get((x[1], x[2]))
            k = x[3]
            if (r is None or num(r[f"n_{k}"]) != x[5]
                    or not near(r[f"dR2_{k}"], x[10], TOL_ROUND4)
                    or not near(r[f"beta_{k}"], x[6], TOL_ROUND4)
                    or not near(r[f"padj_{k}"], x[9], 1e-4)):
                mism += 1
    ck.add("P3", mism == 0, "cov4 and cov7 rows agree with the Part 10h "
           f"Cov4_vs_Cov7 sheet on n, beta, delta R2 and FDR ({mism} "
           "mismatch(es))")
    star = lambda q: ("***" if q < 0.001 else "**" if q < 0.01 else
                      "*" if q < 0.05 else "ns")
    exp_st = {"M4": ["***/***", "***/**", "***/**"],
              "PVM_3": ["***/**", "***/**", "***/**"],
              "PVM_1": ["***/*", "***/*", "***/*"],
              "HKG8": ["ns/ns", "*/ns", "ns/ns"]}
    ip = {(x[0], x[1], x[2], x[3]): x for x in rows_p}
    bad = []
    for lab, vals in exp_st.items():
        for oc, ev in zip(["MMSE score", "CASI score", "Cognitive status"],
                          vals):
            got = (star(ip[(lab, "NeoC", oc, "cov4")][9]) + "/" +
                   star(ip[(lab, "NeoC", oc, "cov7")][9]))
            if got != ev:
                bad.append(f"{lab} {oc} {got} vs {ev}")
    ck.add("P4", not bad, "NeoC cognition stars (cov4/cov7) match the "
           f"earlier Supplementary Figure 14 ({len(bad)} differ" +
           (": " + "; ".join(bad) if bad else "") + ")")
    for blk in BLOCKS:
        n = sum(1 for x in rows_p if x[3] == blk and x[0] != "HKG8"
                and x[2] in OUT9 and x[9] < 0.05)
        m = sum(1 for x in rows_p if x[3] == blk and x[0] != "HKG8"
                and x[2] in CTRL6 and x[9] < 0.05)
        ck.info(f"{blk}: {n} of 54 main-arm clinical rows and {m} of 36 "
                "main-arm negative-control rows with FDR < 0.05")

    # ---- Q. Supplementary Figure 14: leave-one-in scan ----------------------
    print("\nQ. SUPPLEMENTARY FIGURE 14 (LEAVE-ONE-IN SCAN)")
    TECH = ("mito", "numi", "ngene", "doublet", "intronexon", "satur")
    hdr_q = ["arm", "region_set", "outcome", "added_covariate",
             "covariate_type", "n_core", "beta_core", "p_core", "n_added",
             "beta_added", "p_added", "attenuation_pct", "max_VIF",
             "sign_flipped", "lost_significance"]
    rows_q = []
    for lab in ARM_ORDER:
        for r in p10h[lab]["LeaveOneIn_Scan"]:
            if r["arm"] != "Module Post-IP":
                continue
            bc, ba = num(r["beta_core"]), num(r["beta_added"])
            att = (None if bc is None or ba is None or abs(bc) < 0.10
                   else (1.0 - ba / bc) * 100.0)
            rows_q.append([lab, r["dataset"], r["outcome"],
                           r["added_covariate"],
                           "technical" if r["covariate_key"] in TECH
                           else "biological or tissue",
                           num(r["n_core"]), bc, num(r["p_core"]),
                           num(r["n_added"]), ba, num(r["p_added"]), att,
                           num(r["max_VIF"]), tb(r["sign_flipped"]),
                           tb(r["lost_significance"])])
    rows_q.sort(key=lambda x: (ARM_ORDER.index(x[0]), DS_ORDER.index(x[1]),
                               x[2], x[4], x[3]))
    sheets.append(("Q_SuppFig14_leave_one_in", hdr_q, rows_q))
    ck.add("Q1", len(rows_q) == 1440, f"{len(rows_q)} rows (expected 1440 = "
           "4 arms x 2 region sets x 10 outcomes x 18 covariates)")
    m4 = [x for x in rows_q if x[0] == "M4" and x[14]]
    m4_cog = [x for x in m4 if x[2] in ("MMSE score", "CASI score")]
    m4_bio = [x for x in m4 if x[4] != "technical"]
    ck.add("Q2", len(m4) == 13 and not m4_cog and
           all(x[3] == "RIN" for x in m4_bio) and len(m4_bio) == 2,
           f"M4 rows losing significance: {len(m4)} (record 13), on MMSE "
           f"or CASI: {len(m4_cog)} (record 0), non-technical: "
           f"{len(m4_bio)}, all RIN (record 2)", warn_only=True)
    for lab in ARM_ORDER:
        lost = [x for x in rows_q if x[0] == lab and x[14]]
        cog = sum(1 for x in lost if x[2] in ("MMSE score", "CASI score"))
        bio = sum(1 for x in lost if x[4] != "technical")
        ck.info(f"{lab}: {len(lost)} rows lose significance, {cog} on MMSE "
                f"or CASI, {bio} to a biological or tissue covariate")

    # ---- R. Supplementary Figure 14: without the eleven SA donors -----------
    print("\nR. SUPPLEMENTARY FIGURE 14 (WITHOUT THE 11 SA DONORS)")
    TIERS = {"all": "all donors", "drop_severe": "without the 11 SA donors"}
    hdr_r = ["arm", "region_set", "outcome", "donor_set", "n_donors",
             "beta_std", "se", "p_value", "FDR", "delta_R2"]
    rows_r = []
    for lab in ARM_ORDER:
        for r in p10h[lab]["Multivariate"]:
            if (r["tier"] in TIERS and r["predictor"] == "pct_post" and
                    r["covariate_set"] == "cov4" and r["outcome"] in CTRL6):
                rows_r.append([lab, r["dataset"], r["outcome"],
                               TIERS[r["tier"]], num(r["n_donors"]),
                               num(r["beta_std"]), num(r["se"]),
                               num(r["p_value"]), num(r["p_adj"]),
                               num(r["delta_R2"])])
    tord = list(TIERS.values())
    rows_r.sort(key=lambda x: (ARM_ORDER.index(x[0]), DS_ORDER.index(x[1]),
                               CTRL6.index(x[2]), tord.index(x[3])))
    sheets.append(("R_SuppFig14_without_SA", hdr_r, rows_r))
    ck.add("R1", len(rows_r) == 96, f"{len(rows_r)} rows (expected 96 = "
           "4 arms x 2 region sets x 6 negative-control outcomes x 2 donor "
           "sets; the clinical outcomes are in Battery 5)")
    ip4 = {(x[0], x[1], x[2]): x for x in rows_p if x[3] == "cov4"}
    mism = sum(1 for x in rows_r if x[3] == "all donors" and (
        not near(ip4[(x[0], x[1], x[2])][10], x[9], 1e-12) or
        not near(ip4[(x[0], x[1], x[2])][9], x[8], 1e-12)))
    ck.add("R2", mism == 0, "all-donor rows equal the cov4 rows of sheet P "
           f"({mism} mismatch(es))")
    ir = {(x[0], x[1], x[2], x[3]): x for x in rows_r}
    k = ("M4", "NeoC")
    rin = ir[k + ("RIN", TIERS["drop_severe"])]
    ck.add("R3", near(rin[9], 0.033, 0.002) and rin[8] >= 0.05,
           f"M4 NeoC without the SA donors: RIN delta R2 {rin[9]:.3f}, FDR "
           f"{rin[8]:.3f} (record 0.033, not significant)", warn_only=True)

    # ---- S. Supplementary Figure 14: tissue covariates vs disease ----------
    print("\nS. SUPPLEMENTARY FIGURE 14 (TISSUE COVARIATES AND DISEASE)")
    col_a = {h: i for i, h in enumerate(hdr_a)}
    TISSUE = [("RIN", "RIN"), ("brain_pH", "Brain pH"),
              ("fresh_brain_weight", "Fresh brain weight")]
    DIS = [("MMSE", "MMSE score"), ("CASI", "CASI score"),
           ("cognitive_status", "Cognitive status"), ("CERAD", "CERAD score"),
           ("Braak", "Braak stage"), ("Thal", "Thal phase"),
           ("ADNC", "ADNC rank"), ("CAA", "CAA score")]
    donors = {}
    for x in rows_a:
        donors.setdefault(str(x[col_a["donor_id"]]), x)
    raw = []
    for tc, tl in TISSUE:
        for dc, dl in DIS:
            xs = [num(r[col_a[tc]]) for r in donors.values()]
            ys = [num(r[col_a[dc]]) for r in donors.values()]
            n, rho, pv = spearman(xs, ys)
            raw.append([tl, dl, n, rho, pv])
    adj = bh([r[4] for r in raw])
    rows_s = [r + [q] for r, q in zip(raw, adj)]
    hdr_s = ["tissue_covariate", "outcome", "n_donors", "spearman_rho",
             "p_value", "FDR"]
    sheets.append(("S_SuppFig14_tissue_disease", hdr_s, rows_s))
    ck.add("S1", len(rows_s) == 24 and len(donors) == 84,
           f"{len(rows_s)} rows (expected 24 = 3 tissue covariates x 8 "
           f"outcomes) over {len(donors)} donors (expected 84)")
    for r in rows_s:
        if r[0] == "RIN" and r[1] == "ADNC rank" and r[3] is not None:
            ck.info(f"RIN vs ADNC: rho {r[3]:+.3f}, FDR {r[5]:.3g}, n {r[2]}"
                    " (Part 10h record about -0.23)")
    ck.info(f"tissue-disease pairs with FDR < 0.05: "
            f"{sum(1 for r in rows_s if r[5] is not None and r[5] < 0.05)}"
            " of 24")

    # ---- Contents ------------------------------------------------------------
    rel = lambda k: os.path.relpath(P[k], root)
    arm_pat = lambda k: os.path.relpath(P[(k, "M4")], root).replace(
        "Mic13_PVM3", "<arm folder>").replace("_M4_", "_<arm tag>_")
    contents = [
        ("A_donor_frame", "9-11", "Donor-level predictors, outcomes and "
         "covariates for every arm and region set",
         f"{arm_pat('p10h')} (Donor_pct_*); {arm_pat('comp')} (Donor_Frame);"
         f" {rel('p11ba')} (Donor_Pct)"),
        ("B_Fig9A_variance", "9A", "cov4 models, four cell populations",
         f"{arm_pat('p10h')} (Multivariate, tier all, cov4)"),
        ("C_Fig9B_cognition", "9B", "Module post-IP on MMSE and CASI with "
         "pathology added", f"{arm_pat('cvp')} (Models, module)"),
        ("D_Fig9C_pathology", "9C", "Module post-IP on pathology with MMSE "
         "+ CASI added", f"{arm_pat('cvp')} (Reverse, module)"),
        ("E_Fig9D_mod_vs_nonmod", "9D", "Module minus non-module "
         "coefficient, one model; each coefficient also tested on its own "
         "(p_module, FDR_module, p_nonmodule, FDR_nonmodule; BH within arm, "
         "region set and predictor over the six outcomes)",
         f"{rel('p11cvp')} (Contrast); single-coefficient tests computed "
         "here from sheet A"),
        ("F_Fig10A_benchmark", "10A", "Module post-IP and five competitors, "
         "separate cov4 models", f"{arm_pat('comp')} (Marginal, cov4)"),
        ("G_Fig10B_competitors", "10B", "Module post-IP and one competitor "
         "in one model", f"{arm_pat('comp')} (Duels)"),
        ("H_Fig10CD_decomposition", "10C, 10D", "Module post-IP + cluster "
         "proportion + score", f"{arm_pat('mv')} (Decomposition, "
         "Correlations)"),
        ("I_Fig11A_boundary_assoc", "11A", "Module post-IP per boundary on "
         "MMSE and CASI", f"{rel('p11ba')} (Donor_Assoc)"),
        ("J_Fig11B_boundary_curves", "11B", "Center line and boundary lines "
         "of each boundary, M4 NeoC", f"{rel('ba_curves')}"),
        ("K_Fig11BC_cells", "11B, 11C", "Every cell of the 8 sets with "
         "u, v, side of the IP and membership of each boundary",
         f"{rel('ba_cells')}"),
        ("L_Fig11C_scorecard", "11C", "Boundary scorecard and overlap with "
         "CMRF", f"{rel('p11ba')} (Scorecard, Overlap)"),
        ("M_Fig11C_tests", "11C", "Every test behind the scorecard",
         f"{rel('p11ba')} (Full_Panel)"),
        ("N_Fig11C_overlap", "11C", "Jaccard overlap between boundaries",
         f"{rel('p11ba')} (Overlap)"),
        ("O_Fig11BC_boundary_summary", "11B, 11C", "IP_y, coverage, "
         "widths and cell counts per set", f"{rel('p11ba')} (Coverage); "
         f"{rel('ba_sum')}"),
        ("P_SuppFig14_cov4_cov7", "SF14", "Module post-IP% at cov4 and "
         "cov7, nine clinical and six negative-control outcomes",
         f"{arm_pat('p10h')} (Multivariate, "
         "tier all, cov4 and cov7; checked against Cov4_vs_Cov7)"),
        ("Q_SuppFig14_leave_one_in", "SF14", "Module post-IP% at cov4 "
         "with one extra covariate at a time (18 covariates)",
         f"{arm_pat('p10h')} (LeaveOneIn_Scan, Module Post-IP)"),
        ("R_SuppFig14_without_SA", "SF14", "Module post-IP% at cov4 on "
         "the negative-control outcomes, all donors and without the 11 SA "
         "donors",
         f"{arm_pat('p10h')} (Multivariate, tiers all and drop_severe)"),
        ("S_SuppFig14_tissue_disease", "SF14", "Spearman correlation of "
         "RIN, brain pH and brain weight with eight disease outcomes",
         "computed here from sheet A, one row per donor"),
    ]
    return finish(ck, (sheets, contents), out_dir)


def finish(ck, payload, out_dir):
    out_path = os.path.join(out_dir, OUT_NAME)
    if payload is not None:
        sheets, contents = payload
        try:
            write_workbook(out_path, sheets, contents)
            print(f"\nWorkbook written: {out_path}")
            for name, hdr, rows in sheets:
                print(f"  {name}: {len(rows)} rows x {len(hdr)} columns")
        except PermissionError:
            ck.add("W1", False, "workbook could not be written (file open in "
                   "Excel?)")
    c = ck.counts()
    verdict = "FAIL" if c["FAIL"] else "PASS"
    print(f"\nChecks: PASS {c['PASS']} | WARN {c['WARN']} | FAIL {c['FAIL']}")
    print(f"Finished: {datetime.datetime.now():%Y-%m-%d %H:%M:%S}")
    print(f"VERDICT: {verdict}")
    return verdict


def write_workbook(path, sheets, contents):
    font = Font(name="Arial", size=10)
    bold = Font(name="Arial", size=10, bold=True)
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "Contents"
    put(ws, ["sheet", "figure_panel", "content", "source"],
        [list(c) for c in contents], font, bold)
    for name, hdr, rows in sheets:
        put(wb.create_sheet(name), hdr, rows, font, bold)
    tmp = path + ".tmp.xlsx"
    wb.save(tmp)
    os.replace(tmp, path)


def put(ws, hdr, rows, font, bold):
    ws.append(hdr)
    for c in ws[1]:
        c.font = bold
    width = [len(str(h)) for h in hdr]
    for r in rows:
        ws.append(r)
    for row in ws.iter_rows(min_row=2):
        for c in row:
            c.font = font
    for i, r in enumerate(rows[:200]):
        for j, v in enumerate(r):
            width[j] = max(width[j], min(len(str(v)), 60))
    for j, w in enumerate(width):
        ws.column_dimensions[get_column_letter(j + 1)].width = min(w + 2, 62)
    ws.freeze_panes = "A2"


if __name__ == "__main__":
    main()
