#!/usr/bin/env python3
"""
build_SuppTable_S10.py

Builds SuppTable_S10_Fig5_CMRF.xlsx, the copula monotone recovery filter
(CMRF) table behind Figures 5 to 8 and Supplementary Figures S11 and S12.
finalize_supp_tables.py turns this workbook into Supplementary Table S10.

Sheets
  A_cells            one row per GOI-positive cell of each arm and region set:
                     copula coordinates, the four SKCS-core groups, spline
                     training cells, g(u), |v - g(u)|, delta_adapt(u), Stage 1
                     group, Stage 2 removal step, CMRF class, final layer, and
                     the Figure 8 columns (u_SK, v_SK, refit group, 2D UMAP)
  B_band             band and spline parameters and class counts, per dataset
  C_gate             RSk-fCV gate summary, per dataset
  D_gate_trajectory  p and dRMSE at every gate iteration, M4 NeoC
  E_fold_RMSE        the 200 paired fold RMSEs at iterations 1 and 55, M4 NeoC
  F_candidates       the 352 candidates ranked by |v - PAVA(u)|, M4 NeoC
  G_tests            recovered-cell and class tests, per dataset
  H_resid_by_class   |v - g(u)| by CMRF class, per dataset
  I_Fig8_summary     the Figure 8C table, per dataset
  J_Fig8_tests       omnibus, trend and pairwise tests behind Figure 8C

Usage
  Rscript build_SuppTable_S10_step0_Fig8_cells.r     (first, same root)
  python3 build_SuppTable_S10.py                     root: $SEAAD_ROOT, else the current folder
  python3 build_SuppTable_S10.py --root /data
  python3 build_SuppTable_S10.py --root /data --allow-figure-drift

Inputs under the root (read only). <tag>/<folder> is M4/Mic13_PVM3,
PVM_3/PVM_3, PVM_1/PVM_1 or HKG8/HKG8; <ds> is NeoC or NonNeoC.
  Written by Part10e_CP10K_CMRF_code.r:
  Part10e_CMRF/Part10e_CMRF_<folder>/Part10e_CP10K_<tag>_CMRF_<ds>_SK_interaction_analysis.xlsx
      Step4_Refined_Cells, Step4_Boundary_Params, Step4_CMRF_Stage2_Gate,
      Step5_Recovery_Summary, Step4_Spine_PRSS
  Part10e_CMRF/Part10e_CMRF_<folder>/Part10e_CP10K_<tag>_CMRF_<ds>_RSk_fCV_Audit.xlsx
      Fold_RMSE, Wilcoxon_Summary, Gate_Iteration, Gate_PValue_Trajectory
  Part10e_CMRF/Part10e_CMRF_<folder>/Part10e_CP10K_<tag>_CMRF_<ds>_VP_PlotData_T1_to_T7A.xlsx
      T3_JT_ordered_trend
  Part10e_CMRF/Part10e_CMRF_<folder>/Part10e_CP10K_<tag>_CMRF_Verification.xlsx
      <ds>_StatTests
  Written by Part10f_CP10K_CMRF_RSkfCV_Mechanism_code.r:
  Part10f_RSkfCV_Mechanism/Part10f_RSkfCV_Mechanism_Mic13_PVM3/Part10e_CP10K_M4_RSkfCV_NeoC_RSkfCV_Mechanism.xlsx
      Gate_Trajectory, Fold_RMSE, Candidates
  Written by Part10e_CP10K_Validation1_code.r:
  Part10e_Validation1/Part10e_Validation1_<folder>/Part10e_CP10K_<tag>_V1_<ds>_SK_interaction_analysis.xlsx
      Step4_Refined_Cells (refit group labels)
  Written by Part10e_CP10K_Validation2_UMAP_code.r:
  Part10e_Validation2/Part10e_Validation2_<folder>/Part10e_CP10K_<tag>_V2_Validation2_Summary.xlsx
      Figure 8C statistics
  Written by Part10_CP10K_SKCS_code.r:
  Part10/Part10_<folder>/Part10_<folder>_<ds>/Part10_CP10K_<tag>_Module1_<GOI>_<composite>_<ds>_Pooled_SKCS_scores.csv
      S_i and K_i, to check u_SK and v_SK
  Written by build_SuppTable_S10_step0_Fig8_cells.r (from the Validation 1 and 2 RDS):
  Manuscript/Supp_Tables/build/Part10e_CP10K_<tag>_V2_<ds>_Fig8_cells.csv
      2D UMAP coordinates, refit group, u_SK and v_SK
  Manuscript/Supp_Tables/build/build_SuppTable_S10_step0_Fig8_cells_log.txt
      must end with VERDICT: PASS

Every value is copied from these files except: the four-group label, Stage 1
group and final-layer flag (derived from the per-cell columns); the absolute
dRMSE and BV (from the 200 fold pairs; the workbooks store them as ratios to
the median RMSE_A); and the candidate barcodes (matched on u and v). The final
Wilcoxon p is read from Wilcoxon_Summary, because the gate summary sheets store
it rounded to 0.

Checks. The script prints every check. Two kinds are kept apart:
  - consistency checks between the input files (for example the CMRF class
    against the Stage 1 and Stage 2 columns); any failure stops the build;
  - figure checks against values typed from the published figures (REF_N,
    REF_KEYS, REF_STAGE, REF_SPLINE_N, REF_FIG7, FIG8_N, FIG8_TABLE). By default
    a failure also stops the build. With --allow-figure-drift these failures
    are reported as DRIFT and the workbook is still written. The Figure 8C
    neighbor-enrichment, kBET and position values come from the UMAP and PCA of
    Validation 2, which can differ slightly between computers.
The Figure 8C source of each column (cell or donor unit, mean or median
enrichment) is the one that reproduces the figure in all 8 datasets; if none
does, the closest source is used and the figure check fails.

Output, in <root>/Manuscript/Supp_Tables/build/ (created if missing):
  SuppTable_S10_Fig5_CMRF.xlsx
  build_SuppTable_S10_log.txt, ending with VERDICT: PASS or VERDICT: FAIL
"""
import argparse
import atexit
import datetime as dt
import glob
import os
import sys
import traceback

import numpy as np
import pandas as pd
from openpyxl import Workbook
from openpyxl.cell import WriteOnlyCell
from openpyxl.styles import Font, Alignment
from openpyxl.utils import get_column_letter

try:
    from scipy.stats import wilcoxon, chi2
    HAVE_SCIPY = True
except ImportError:
    HAVE_SCIPY = False

#####
# Settings
#####
ap = argparse.ArgumentParser(description="Build SuppTable_S10_Fig5_CMRF.xlsx")
ap.add_argument("--root", default=os.environ.get("SEAAD_ROOT") or os.getcwd(),
                help="folder holding the pipeline outputs (default: $SEAAD_ROOT, else the current folder)")
ap.add_argument("--allow-figure-drift", action="store_true",
                help="report failed figure checks as DRIFT and still write the workbook")
ARGS = ap.parse_args()
ROOT = os.path.abspath(ARGS.root)
if not os.path.isdir(ROOT):
    raise SystemExit("root folder not found: %s" % ROOT)
OUT_DIR = os.path.join(ROOT, "Manuscript", "Supp_Tables", "build")
OUT = os.path.join(OUT_DIR, "SuppTable_S10_Fig5_CMRF.xlsx")
LOG = os.path.join(OUT_DIR, "build_SuppTable_S10_log.txt")
STEP0_LOG = os.path.join(OUT_DIR, "build_SuppTable_S10_step0_Fig8_cells_log.txt")
os.makedirs(OUT_DIR, exist_ok=True)
if os.path.isfile(OUT):
    os.remove(OUT)   # a failed run must not leave the workbook of an earlier run


class Tee:
    """Copies everything printed to the console into the build log."""

    def __init__(self, stream, fh):
        self.stream, self.fh = stream, fh

    def write(self, s):
        self.stream.write(s)
        self.fh.write(s)

    def flush(self):
        self.stream.flush()
        self.fh.flush()


_LOG_FH = open(LOG, "w", encoding="utf-8")
sys.stdout = Tee(sys.__stdout__, _LOG_FH)
SAVED = False


def _excepthook(etype, value, tb):
    print("".join(traceback.format_exception(etype, value, tb)))


def _verdict():
    print("VERDICT: PASS" if SAVED else "VERDICT: FAIL")
    sys.stdout.flush()
    sys.stdout = sys.__stdout__
    _LOG_FH.close()


sys.excepthook = _excepthook
atexit.register(_verdict)


def die(msg):
    print("STOP: %s" % msg)
    sys.exit(1)


print("build_SuppTable_S10.py")
print("Date  : %s" % dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
print("Root  : %s" % ROOT)
print("Output: %s" % OUT)
print("Figure checks: %s\n" % ("DRIFT allowed" if ARGS.allow_figure_drift else "must pass"))

ARMS = [("M4", "Mic13_PVM3"), ("PVM_3", "PVM_3"), ("PVM_1", "PVM_1"), ("HKG8", "HKG8")]
DATASETS = ["NeoC", "NonNeoC"]
GROUPS = ["Pre-IP TREP", "Pre-IP non-TREP", "Post-IP TREP", "Post-IP non-TREP"]
CLASSES = ["Original + Refined", "Recovered", "Original only", "Excluded"]
TOL_BAND = 1e-8

# Values printed in the figures, checked against the rebuilt cells.
# Final-module keys (Figure 6 panel 8, Supplementary Figure 12): OriRef by the
# four groups, Recovered, non-module cells.
REF_KEYS = {
    ("M4", "NeoC"): (615, 3, 668, 72, 298, 474),
    ("M4", "NonNeoC"): (560, 3, 619, 50, 242, 513),
    ("PVM_3", "NeoC"): (3633, 6, 3836, 628, 1351, 3877),
    ("PVM_3", "NonNeoC"): (3001, 31, 3210, 404, 1089, 3397),
    ("PVM_1", "NeoC"): (755, 1, 921, 144, 381, 1016),
    ("PVM_1", "NonNeoC"): (701, 0, 803, 114, 318, 722),
    ("HKG8", "NeoC"): (71, 0, 147, 21, 287, 906),
    ("HKG8", "NonNeoC"): (43, 0, 116, 32, 240, 648),
}
REF_N = {("M4", "NeoC"): 2130, ("M4", "NonNeoC"): 1987, ("PVM_3", "NeoC"): 13331,
         ("PVM_3", "NonNeoC"): 11132, ("PVM_1", "NeoC"): 3218, ("PVM_1", "NonNeoC"): 2658,
         ("HKG8", "NeoC"): 1432, ("HKG8", "NonNeoC"): 1079}
# Stage 1 candidates within delta_base and within delta_adapt only; Original
# only; removed by the gate (Figure 6 panels 5-7, Supplementary Figure 12A)
REF_STAGE = {("M4", "NeoC"): (292, 60, 61, 54), ("PVM_3", "NeoC"): (1583, 411, 304, 643)}
REF_SPLINE_N = {("M4", "NeoC"): 1308}
REF_FIG7 = dict(p_iter1=1.40e-34, p_final=1.39e-25, d_iter1=0.0099, d_final=0.0044,
                bv=0.0046, n_iter=55, dcor=0.557)

# Figure 8 (B) group sizes: pre-IP TREP, Recovered, post-IP TREP (TREP of the
# Validation 1 refit) and (C) the summary table as printed in the figure.
FIG8_N = {("M4", "NeoC"): (618, 298, 702), ("M4", "NonNeoC"): (563, 242, 644),
          ("PVM_3", "NeoC"): (3635, 1351, 4039), ("PVM_3", "NonNeoC"): (3032, 1089, 3327),
          ("PVM_1", "NeoC"): (755, 381, 934), ("PVM_1", "NonNeoC"): (701, 318, 853),
          ("HKG8", "NeoC"): (71, 287, 156), ("HKG8", "NonNeoC"): (43, 240, 132)}
F8_COLS = ["br_pre", "br_post", "br_rec", "br_exc", "fdr_br", "ratio", "pct_si", "pct_ki",
           "pos_med", "pos_q1", "pos_q3", "fdr_pos", "ne_pre", "ne_post", "kb_pre", "kb_post"]
FIG8_TABLE = {
    ("M4", "NeoC"): (0.033, 0.138, 0.166, 0.204, 3.70e-6, 1.9, 99.7, 0.3, 0.28, 0.08, 0.58, 2.43e-20, 0.94, 0.84, 0.28, 0.57),
    ("M4", "NonNeoC"): (0.025, 0.050, 0.104, 0.211, 1.38e-15, 2.8, 95.9, 4.1, 0.31, 0.05, 0.58, 1.87e-12, 0.98, 0.92, 0.17, 0.48),
    ("PVM_3", "NeoC"): (0.012, 0.112, 0.137, 0.183, 2.54e-13, 2.2, 98.7, 1.3, 0.23, -0.03, 0.59, 4.63e-29, 0.93, 0.83, 0.44, 0.59),
    ("PVM_3", "NonNeoC"): (0.010, 0.082, 0.116, 0.179, 1.96e-13, 2.5, 98.6, 1.4, 0.28, 0.01, 0.61, 5.55e-29, 0.93, 0.87, 0.37, 0.56),
    ("PVM_1", "NeoC"): (0.035, 0.138, 0.171, 0.212, 4.09e-9, 2.0, 95.0, 5.8, 0.35, 0.07, 0.65, 6.11e-22, 0.97, 0.89, 0.35, 0.51),
    ("PVM_1", "NonNeoC"): (0.030, 0.141, 0.168, 0.206, 1.90e-9, 2.0, 98.1, 1.9, 0.27, 0.06, 0.54, 8.29e-20, 0.97, 0.90, 0.29, 0.51),
    ("HKG8", "NeoC"): (0.010, 0.020, 0.318, 0.735, 1.37e-14, 21.0, 29.3, 100.0, 0.22, 0.07, 0.44, 2.34e-11, 0.94, 0.54, 0.46, 0.87),
    ("HKG8", "NonNeoC"): (0.011, 0.026, 0.319, 0.751, 0.002, 17.2, 28.3, 100.0, 0.27, 0.09, 0.49, 3.41e-7, 0.94, 0.68, 0.42, 0.89),
}
FOUR_GROUPS = ["Pre-IP TREP", "Pre-IP non-TREP", "Post-IP TREP", "Post-IP non-TREP"]

#####
# Helpers
#####
CHECKS = []       # (passed, is_figure_check, label)


def show(label, got, ref, tol=0.0, rel=False, fig=False):
    """fig=True: the reference is a value typed from a published figure."""
    if got is None or ref is None or (isinstance(got, float) and not np.isfinite(got)):
        ok = False
    elif isinstance(ref, (bool, str)):
        ok = got == ref
    elif rel:
        ok = abs(got - ref) <= tol * abs(ref)
    else:
        ok = abs(got - ref) <= tol
    flag = "ok" if ok else ("DRIFT" if fig and ARGS.allow_figure_drift else "CHECK")
    print("  %-62s %14s %14s  %s" % (label, got, ref, flag))
    CHECKS.append((ok, fig, label))
    return ok


def need(path):
    if not os.path.isfile(path):
        die("input not found: %s" % path)
    return path


def sheet(path, name):
    return pd.read_excel(need(path), sheet_name=name)


def as_bool(s):
    """Excel booleans come back as bool, 0/1 or TRUE/FALSE strings."""
    if s.dtype == bool:
        return s
    m = {"TRUE": True, "FALSE": False, "1": True, "0": False, "1.0": True, "0.0": False}
    out = s.map(lambda x: x if isinstance(x, (bool, np.bool_)) else m.get(str(x).strip().upper()))
    if out.isna().any():
        die("column %s has values that are not TRUE/FALSE" % s.name)
    return out.astype(bool)


def p10e_file(tag, folder, ds, kind):
    return os.path.join(ROOT, "Part10e_CMRF", "Part10e_CMRF_%s" % folder,
                        "Part10e_CP10K_%s_CMRF_%s_%s.xlsx" % (tag, ds, kind))


def verification_file(tag, folder):
    return os.path.join(ROOT, "Part10e_CMRF", "Part10e_CMRF_%s" % folder,
                        "Part10e_CP10K_%s_CMRF_Verification.xlsx" % tag)


def v1_file(tag, folder, ds):
    return os.path.join(ROOT, "Part10e_Validation1", "Part10e_Validation1_%s" % folder,
                        "Part10e_CP10K_%s_V1_%s_SK_interaction_analysis.xlsx" % (tag, ds))


def v2_dir(folder):
    return os.path.join(ROOT, "Part10e_Validation2", "Part10e_Validation2_%s" % folder)


def iqr7(x):
    q = np.percentile(x, [25, 75])   # linear interpolation = R type 7
    return float(q[1] - q[0])


#####
# Per dataset
#####
A_parts, B_rows, C_rows, G_parts, H_parts = [], [], [], [], []
# The Figure 8 cell files must come from a passing step0 run.
_last = ""
if os.path.isfile(STEP0_LOG):
    with open(STEP0_LOG, encoding="utf-8") as fh:
        _lines = [x.strip() for x in fh if x.strip()]
    _last = _lines[-1] if _lines else ""
print("step0 log %s: %s" % (STEP0_LOG, _last or "MISSING"))
if _last != "VERDICT: PASS":
    die("run build_SuppTable_S10_step0_Fig8_cells.r with the same root first; its log must end with VERDICT: PASS")
print("Checks (this run against the source workbooks and the figures)")
print("  %-62s %14s %14s" % ("check", "this run", "reference"))

for tag, folder in ARMS:
    ver = verification_file(tag, folder)
    for ds in DATASETS:
        key = (tag, ds)
        lab = "%s %s" % key
        ska = p10e_file(tag, folder, ds, "SK_interaction_analysis")
        aud = p10e_file(tag, folder, ds, "RSk_fCV_Audit")
        vp = p10e_file(tag, folder, ds, "VP_PlotData_T1_to_T7A")

        rc = sheet(ska, "Step4_Refined_Cells")
        bp = sheet(ska, "Step4_Boundary_Params").iloc[0]
        gate = sheet(ska, "Step4_CMRF_Stage2_Gate").iloc[0]
        rec = sheet(ska, "Step5_Recovery_Summary")
        prss = sheet(ska, "Step4_Spine_PRSS")
        folds = sheet(aud, "Fold_RMSE")
        wsum = sheet(aud, "Wilcoxon_Summary").iloc[0]
        gi = sheet(aud, "Gate_Iteration")
        gtraj = sheet(aud, "Gate_PValue_Trajectory")

        n = len(rc)
        show(lab + " cells", n, REF_N[key], fig=True)
        if rc["barcode"].duplicated().any():
            die("%s: duplicated barcodes in Step4_Refined_Cells" % lab)

        core = as_bool(rc["is_orig_skcs"])
        trep = as_bool(rc["is_TREP"].fillna(False)) & core
        s1 = as_bool(rc["pass_stage1"])
        removed = as_bool(rc["removed_by_stage2"])
        cls = rc["cell_class"].astype(str)
        u = rc["u_all"].to_numpy(float)
        v = rc["v_all"].to_numpy(float)
        g = rc["g_u"].to_numpy(float)
        resid = rc["resid_abs"].to_numpy(float)
        dadapt = rc["delta_adapt"].to_numpy(float)
        # delta_base and K(u) are taken from the per-cell columns, which keep
        # full precision; the boundary sheet may round them.
        dbase_cells = rc["delta_base"].to_numpy(float)
        dbase = float(dbase_cells[0])
        show(lab + " delta_base constant over cells", float(np.ptp(dbase_cells)), 0.0, tol=1e-12)
        show(lab + " delta_base = boundary sheet", dbase, float(bp["delta_base"]), tol=1e-6)
        kappa, bw, ipu = float(bp["kappa"]), float(bp["bw_gaussian"]), float(bp["IP_u"])
        kden = rc["K_density"].to_numpy(float)

        # Band arithmetic, recomputed cell by cell
        k_re = np.exp(-(u - ipu) ** 2 / (2 * bw ** 2))
        show(lab + " K(u) = Gaussian kernel at IP_u (max |dev|)", float(np.max(np.abs(k_re - kden))), 0.0, tol=1e-5)
        dad_re = dbase * (1 + kappa * kden)
        show(lab + " delta_adapt = delta_base(1 + kappa K) (max |dev|)", float(np.max(np.abs(dad_re - dadapt))), 0.0, tol=TOL_BAND)
        show(lab + " |v - g(u)| = column (max |dev|)", float(np.max(np.abs(np.abs(v - g) - resid))), 0.0, tol=1e-9)
        show(lab + " Stage 1 = |v - g(u)| <= delta_adapt", int(((resid <= dadapt) != s1.to_numpy()).sum()), 0)
        show(lab + " delta_max = delta_base(1 + kappa)", float(bp["delta_max_waist"]), dbase * (1 + kappa), tol=1e-6)

        # Four SKCS-core groups from IP status and TREP
        ip = rc["IP_status"].astype(str)
        pre = ip.str.contains("Pre", case=False)
        post = ip.str.contains("Post", case=False)
        if (core & ~(pre ^ post)).any():
            die("%s: SKCS-core cells without a single pre-IP or post-IP label" % lab)
        grp = np.where(~core, "Non-SKCS-core",
                       np.where(pre, "Pre-IP", "Post-IP") + np.where(trep, " TREP", " non-TREP"))

        # Stage 1 groups and CMRF classes
        stage1 = np.where(core & s1, "OriRef", np.where(core, "Original only",
                          np.where(s1, "Candidate", "Not admitted")))
        exp_cls = np.where(core & s1, "Original + Refined", np.where(core, "Original only",
                           np.where(s1 & ~removed, "Recovered", "Excluded")))
        show(lab + " CMRF class = Stage 1 and Stage 2 columns", int((exp_cls != cls.to_numpy()).sum()), 0)
        in_base = resid <= dbase
        cand = (~core) & s1
        show(lab + " removed cells are all candidates", int((removed & ~cand).sum()), 0)

        # Removal step from Gate_Iteration
        step = pd.Series(np.nan, index=rc.index)
        pos = dict(zip(rc["barcode"], rc.index))
        miss = [b for b in gi["barcode"] if b not in pos]
        show(lab + " Gate_Iteration barcodes found", len(gi) - len(miss), len(gi))
        for b, st in zip(gi["barcode"], gi["step"]):
            if b in pos:
                step[pos[b]] = st
        show(lab + " removed = Gate_Iteration cells", int((removed != step.notna()).sum()), 0)

        # Recovered cells: pre-IP or post-IP on the composite at IP_y
        gref = rc["group_refined"].astype(str)
        is_rec = cls == "Recovered"
        rec_side = np.where(is_rec & gref.str.contains("Pre-IP"), "Pre-IP",
                            np.where(is_rec & gref.str.contains("Post-IP"), "Post-IP", ""))
        show(lab + " every Recovered cell has a side", int((is_rec & (pd.Series(rec_side) == "")).sum()), 0)
        if "cDFG_expr" in rc.columns and np.isfinite(float(bp["IP_y"])):
            side_re = np.where(rc["cDFG_expr"].to_numpy(float) < float(bp["IP_y"]), "Pre-IP", "Post-IP")
            show(lab + " Recovered side = composite vs IP_y", int(((side_re != rec_side) & is_rec).sum()), 0)
        recd = dict(zip(rec["category"].astype(str), rec["n_cells"]))
        show(lab + " Recovered pre-IP vs Step5", int((rec_side == "Pre-IP").sum()), int(recd.get("Recovered Pre-IP", -1)))
        show(lab + " Recovered post-IP vs Step5", int((rec_side == "Post-IP").sum()), int(recd.get("Recovered Post-IP", -1)))

        # Counts against the workbooks and the figure keys
        n_ori = int((cls == "Original + Refined").sum())
        n_rec = int(is_rec.sum())
        n_oo = int((cls == "Original only").sum())
        n_exc = int((cls == "Excluded").sum())
        n_rem = int(removed.sum())
        n_cand = int(cand.sum())
        show(lab + " OriRef vs gate sheet", n_ori, int(gate["n_grp1a"]))
        show(lab + " candidates vs gate sheet", n_cand, int(gate["n_stage1_candidates"]))
        show(lab + " removed vs gate sheet", n_rem, int(gate["n_removed_by_gate"]))
        show(lab + " recovered vs gate sheet", n_rec, int(gate["n_accepted"]))
        show(lab + " Original only vs boundary sheet", n_oo, int(bp["n_orig_only"]))
        ori_grp = [int(((cls == "Original + Refined") & (grp == gname)).sum()) for gname in GROUPS]
        got_keys = tuple(ori_grp) + (n_rec, n_oo + n_exc)
        show(lab + " final-module key (4 groups, Rec, non-module)", str(got_keys), str(REF_KEYS[key]), fig=True)
        if key in REF_STAGE:
            nb = int((cand & in_base).sum())
            show(lab + " candidates within delta_base / adapt only / Orig only / removed",
                 str((nb, n_cand - nb, n_oo, n_rem)), str(REF_STAGE[key]), fig=True)
        n_spline = int(trep.sum())
        show(lab + " spline training cells vs boundary sheet", n_spline, int(bp["spine_n_cells"]))
        if key in REF_SPLINE_N:
            show(lab + " spline training cells vs figure", n_spline, REF_SPLINE_N[key], fig=True)

        # Gate: absolute dRMSE and BV from the 200 final fold pairs
        a = folds["RMSE_A"].to_numpy(float)
        b = folds["RMSE_B"].to_numpy(float)
        ok_pair = np.isfinite(a) & np.isfinite(b)
        med_a = float(np.median(a[np.isfinite(a)]))
        iqr_a = iqr7(a[np.isfinite(a)])
        d_abs = float(np.median(b[ok_pair] - a[ok_pair]))
        show(lab + " dRMSE ratio recomputed vs gate sheet", round(d_abs / med_a, 6), round(float(gate["dRMSE_final"]), 6), tol=2e-6)
        show(lab + " BV ratio recomputed vs gate sheet", round(iqr_a / med_a, 6), round(float(gate["BV"]), 6), tol=2e-6)
        show(lab + " closed on the effect size (dRMSE <= BV)", bool(d_abs <= iqr_a), True)
        p_final = float(wsum["p_value"])
        show(lab + " final p < 0.05", bool(p_final < 0.05), True)
        p_last = float(gtraj["p_value"].iloc[-1])
        if p_last > 0:
            show(lab + " final p = last trajectory p", p_last, p_final, tol=1e-4, rel=True)
        else:
            print("  %-62s %14s %14s  (trajectory stores p rounded to 0; not compared)"
                  % (lab + " final p = last trajectory p", p_last, p_final))
        if HAVE_SCIPY:
            p_re = float(wilcoxon(b[ok_pair], a[ok_pair], alternative="greater",
                                  correction=True, method="approx").pvalue)
            show(lab + " final p recomputed (log10)", round(np.log10(p_re), 2), round(np.log10(p_final), 2), tol=0.05)
        iters = len(gtraj)
        show(lab + " iterations = removals + 1", iters, n_rem + 1)

        # Spline basis size
        best_k = None
        if "best_k" in prss.columns and prss["best_k"].notna().any():
            best_k = int(prss["best_k"].dropna().iloc[-1])

        A_parts.append(pd.DataFrame({
            "arm": tag, "region_set": ds, "region": rc["region"], "barcode": rc["barcode"],
            "u": u, "v": v, "SKCS_group": grp, "spline_training": trep.to_numpy(),
            "g_u": g, "abs_v_minus_g_u": resid, "delta_adapt_u": dadapt,
            "within_delta_base": in_base, "within_delta_adapt": s1.to_numpy(),
            "stage1_group": stage1, "removal_step": step.to_numpy(),
            "CMRF_class": cls.to_numpy(), "final_layer": cls.isin(["Original + Refined", "Recovered"]).to_numpy(),
            "recovered_IP_side": np.where(rec_side == "", None, rec_side)}))

        # Figure 8 per-cell columns: Validation 1 refit labels, Validation 2 UMAP,
        # copula score coordinates (written from the RDS by build_SuppTable_S10_step0_Fig8_cells.r)
        v1 = sheet(v1_file(tag, folder, ds), "Step4_Refined_Cells")
        f8p = need(os.path.join(OUT_DIR, "Part10e_CP10K_%s_V2_%s_Fig8_cells.csv" % (tag, ds)))
        f8 = pd.read_csv(f8p)
        show(lab + " Figure 8 cells (UMAP export)", len(f8), n)
        f8 = f8.set_index("barcode")
        v1 = v1.set_index("barcode")
        miss8 = int((~rc["barcode"].isin(f8.index)).sum())
        show(lab + " cells without a Figure 8 row", miss8, 0)
        f8 = f8.reindex(rc["barcode"])
        v1 = v1.reindex(rc["barcode"])
        rec8 = (f8["cell_class"].astype(str) == "Recovered").to_numpy()
        gref = f8["group_ref"].astype(str)
        ori8 = (f8["cell_class"].astype(str) == "Original + Refined").to_numpy()
        fig8_group = np.where(rec8, "Recovered",
                              np.where(ori8 & gref.isin(FOUR_GROUPS).to_numpy(), gref, "Non-Module"))
        show(lab + " Figure 8 class = CMRF class", int(((f8["cell_class"].astype(str).to_numpy()) != cls.to_numpy()).sum()), 0)
        n8 = (int((fig8_group == "Pre-IP TREP").sum()), int(rec8.sum()), int((fig8_group == "Post-IP TREP").sum()))
        show(lab + " Figure 8B groups (pre-IP TREP, Rec, post-IP TREP)", str(n8), str(FIG8_N[key]), fig=True)
        ori = (cls == "Original + Refined").to_numpy()
        show(lab + " refit group: UMAP export = Validation 1 sheet",
             int((gref.to_numpy()[ori] != v1["group_refined"].astype(str).to_numpy()[ori]).sum()), 0)
        sc_hits = glob.glob(os.path.join(ROOT, "Part10", "*", "*",
                                         "Part10_CP10K_%s_Module1_*_%s_Pooled_SKCS_scores.csv" % (tag, ds)))
        if len(sc_hits) != 1:
            die("%s: expected one Part 10 pooled score CSV, found %d" % (lab, len(sc_hits)))
        sc10 = pd.read_csv(sc_hits[0])
        sc10 = sc10[~sc10["is_GOI_neg"].astype(bool)]
        n10 = len(sc10)
        sc10 = sc10.set_index("barcode")
        us_re = (sc10["spear_norm"].rank(method="average") / (n10 + 1)).reindex(rc["barcode"]).to_numpy()
        vs_re = (sc10["kend_contrib"].rank(method="average") / (n10 + 1)).reindex(rc["barcode"]).to_numpy()
        show(lab + " u_SK = rank(S_i)/(n+1), Part 10 scores (max |dev|)", float(np.nanmax(np.abs(us_re - f8["u_SK"].to_numpy()))), 0.0, tol=1e-9)
        show(lab + " v_SK = rank(K_i)/(n+1), Part 10 scores (max |dev|)", float(np.nanmax(np.abs(vs_re - f8["v_SK"].to_numpy()))), 0.0, tol=1e-9)
        A_parts[-1]["u_SK"] = f8["u_SK"].to_numpy()
        A_parts[-1]["v_SK"] = f8["v_SK"].to_numpy()
        A_parts[-1]["refit_group_Fig8"] = fig8_group
        A_parts[-1]["UMAP_1"] = f8["UMAP_1"].to_numpy()
        A_parts[-1]["UMAP_2"] = f8["UMAP_2"].to_numpy()

        B_rows.append(dict(
            arm=tag, region_set=ds, n_cells=n, IP_x=float(bp["IP_x"]), IP_u=ipu, IP_y=float(bp["IP_y"]),
            IP_v=float(bp["IP_v"]), delta_base=dbase, delta_max=float(bp["delta_max_waist"]),
            kappa=kappa, bw=bw, coverage=float(bp["coverage"]), spline_cells=n_spline,
            spline_donors=bp.get("spine_n_donors"), spline_k=best_k,
            spline_non_decreasing=bp.get("spine_monotone"), spline_min_step=bp.get("spine_min_diff"),
            n_SKCS_core=int(core.sum()), n_OriRef=n_ori, n_Original_only=n_oo,
            n_candidates_delta_base=int((cand & in_base).sum()),
            n_candidates_delta_adapt_only=int((cand & ~in_base).sum()),
            n_removed=n_rem, n_Recovered=n_rec, n_Recovered_pre_IP=int((rec_side == "Pre-IP").sum()),
            n_Recovered_post_IP=int((rec_side == "Post-IP").sum()), n_Excluded=n_exc,
            n_final_layer=n_ori + n_rec))

        C_rows.append(dict(
            arm=tag, region_set=ds, n_OriRef=n_ori, n_candidates=n_cand, iterations=iters,
            n_removed=n_rem, n_Recovered=n_rec, median_RMSE_A=med_a, BV=iqr_a, dRMSE_final=d_abs,
            BV_ratio=iqr_a / med_a, dRMSE_final_ratio=d_abs / med_a,
            Wilcoxon_V_final=wsum.get("statistic"), p_final=p_final,
            closed_on="dRMSE <= BV" if d_abs <= iqr_a else ("p >= 0.05" if p_final >= 0.05 else "not closed")))

        # Tests and residuals by class
        st = sheet(ver, "%s_StatTests" % ds)
        st.insert(0, "region_set", ds)
        st.insert(0, "arm", tag)
        G_parts.append(st)
        t3 = sheet(vp, "T3_JT_ordered_trend")
        t3 = t3.drop(columns=[c for c in t3.columns if c.endswith("_f")])
        t3.insert(0, "region_set", ds)
        t3.insert(0, "arm", tag)
        H_parts.append(t3)
        for cname in CLASSES:
            r = t3[t3["class"] == cname]
            if len(r):
                show("%s median |v - g(u)| %s" % (lab, cname), round(float(r["med_resid"].iloc[0]), 6),
                     round(float(np.median(resid[cls.to_numpy() == cname])), 6), tol=1e-6)

A_cells = pd.concat(A_parts, ignore_index=True)
B_band = pd.DataFrame(B_rows)
C_gate = pd.DataFrame(C_rows)

# Tests: keep the columns that carry a value in at least one row
G_tests = pd.concat(G_parts, ignore_index=True)
DROP = {"dataset", "H0", "normality_SW", "all_normal", "posthoc_summary", "interpretation",
        "SW_p_g1", "SW_p_g2", "normal_g1", "normal_g2", "cv_method", "alpha_used",
        "n_total_folds", "mean_grp1a", "mean_union"}
G_tests = G_tests[[c for c in G_tests.columns if c not in DROP and G_tests[c].notna().any()]]
H_resid = pd.concat(H_parts, ignore_index=True)

# The Verification workbooks round p to six decimals (and R returns 0 below
# the floating-point range), so several p-values are stored as 0. Restore what
# can be restored exactly and label the rest.
#   Wilcoxon signed-rank (the RSk-fCV gate): exact p from sheet C.
#   Kruskal-Wallis: chi-squared tail of the stored statistic, df = classes - 1.
#   Other p stored as 0 (two-sided Mann-Whitney on the large arms): text.
#   BH-adjusted p stored as 0: text, since the source rounds them to 0.
P_TINY_RAW = "< 1 \u00d7 10\u207b\u00b3\u2070\u2070 (stored as 0)"
P_TINY_BH = "< 1 \u00d7 10\u207b\u2076 (rounded to 0 in the source)"
print("\nTest p-values restored from rounded source values")
for col in ("p_value", "p_value_BH_FDR", "overall_stat"):
    if col not in G_tests.columns:
        G_tests[col] = np.nan
G_tests["p_value"] = G_tests["p_value"].astype(object)
G_tests["p_value_BH_FDR"] = G_tests["p_value_BH_FDR"].astype(object)
for i, r in G_tests.iterrows():
    tname = str(r["test"])
    lab = "%s %s %s" % (r["arm"], r["region_set"], tname.split(" (")[0])
    p = r["p_value"]
    if "Wilcoxon signed-rank" in tname:
        pc = float(C_gate[(C_gate["arm"] == r["arm"]) & (C_gate["region_set"] == r["region_set"])]["p_final"].iloc[0])
        show(lab + ": stored p = exact p rounded", float(p), round(pc, 6), tol=1e-12)
        G_tests.at[i, "p_value"] = pc
    elif "Kruskal-Wallis" in tname and HAVE_SCIPY:
        hk = H_resid[(H_resid["arm"] == r["arm"]) & (H_resid["region_set"] == r["region_set"])]
        k = int((hk["n"] > 0).sum())
        pk = float(chi2.sf(float(r["overall_stat"]), k - 1))
        show(lab + ": stored p = chi-squared tail rounded", float(p), round(pk, 6), tol=1e-12)
        G_tests.at[i, "p_value"] = pk if pk > 0 else P_TINY_RAW
    elif isinstance(p, (int, float)) and float(p) == 0:
        G_tests.at[i, "p_value"] = P_TINY_RAW
    pb = r["p_value_BH_FDR"]
    if isinstance(pb, (int, float)) and np.isfinite(pb) and float(pb) == 0:
        G_tests.at[i, "p_value_BH_FDR"] = P_TINY_BH

#####
# Figure 8C: the summary table, from each arm's Validation 2 workbook
#####
print("\nFigure 8C (Validation 2 summary workbooks)")


def near(got, ref, kind):
    if got is None or not np.isfinite(got):
        return False
    if kind == "p":
        if abs(ref) >= 1e-3:
            return abs(got - ref) <= 5e-4 + 1e-12
        return abs(got - ref) <= 0.006 * abs(ref)
    return abs(got - ref) <= {"d3": 6e-4, "d2": 6e-3, "d1": 6e-2}[kind] + 1e-12


def miss_score(got, ref, kind):
    """Distance of a value from the figure value, in units of the rounding tolerance."""
    if got is None or not np.isfinite(got):
        return np.inf
    if kind == "p":
        tol = 5e-4 if abs(ref) >= 1e-3 else 0.006 * abs(ref)
    else:
        tol = {"d3": 6e-4, "d2": 6e-3, "d1": 6e-2}[kind]
    return abs(got - ref) / tol


def pick(alts, kind, col):
    """Choose the source that reproduces the figure value in all 8 datasets.
    If none does, return the closest source and ok = False."""
    for name, vals in alts.items():
        if all(near(vals.get(k), FIG8_TABLE[k][F8_COLS.index(col)], kind) for k in FIG8_TABLE):
            return name, vals, True
    print("  no source reproduces %s; candidates:" % col)
    for name, vals in alts.items():
        print("    %-26s %s" % (name, ", ".join("%s %s" % ("%s/%s" % k, None if vals.get(k) is None else "%.4g" % vals[k]) for k in FIG8_TABLE)))
    score = {name: sum(miss_score(vals.get(k), FIG8_TABLE[k][F8_COLS.index(col)], kind) for k in FIG8_TABLE)
             for name, vals in alts.items()}
    best = min(score, key=score.get)
    if not np.isfinite(score[best]):
        return None, {}, False
    print("    closest: %s" % best)
    return best, alts[best], False


raw = {}
for tag, folder in ARMS:
    wbp = need(os.path.join(v2_dir(folder), "Part10e_CP10K_%s_V2_Validation2_Summary.xlsx" % tag))
    omn = sheet(wbp, "Fig11_Omnibus_BH")
    for ds in DATASETS:
        raw[(tag, ds)] = dict(
            desc=sheet(wbp, "%s_Fig11_Descriptives" % ds), post=sheet(wbp, "%s_Fig11_PostHoc" % ds),
            nb=sheet(wbp, "%s_Neighbour_Proportions" % ds), kb=sheet(wbp, "%s_kBET" % ds),
            cut=sheet(wbp, "%s_Fig11_SKCS_Cuts" % ds).iloc[0], omn=omn[omn["dataset"] == ds])


def dval(d, stat, grp, unit, col):
    r = d[(d["statistic"] == stat) & (d["group"] == grp) & (d["unit"] == unit)]
    return float(r[col].iloc[0]) if len(r) else None


def post_p(d, stat, unit, ga, gb):
    r = d[(d["statistic"] == stat) & (d["unit"] == unit) &
          (((d["group_a"] == ga) & (d["group_b"] == gb)) | ((d["group_a"] == gb) & (d["group_b"] == ga)))]
    return float(r["p_BH_FDR"].iloc[0]) if len(r) else None


def omn_p(d, stat, unit):
    r = d[(d["statistic"] == stat) & (d["unit"] == unit)]
    return float(r["trend_p_BH"].iloc[0]) if len(r) else None


def nb_val(d, ref, col):
    r = d[(d["focal"] == "Recovered") & (d["reference"] == ref)]
    if "unit" in d.columns:
        r = r[r["unit"] == "cell"] if (r["unit"] == "cell").any() else r
    return float(r[col].iloc[0]) if len(r) else None


BR, PP = "Bridge_distance_copula", "PrePost_projection"
choice, F8 = {}, {k: {} for k in FIG8_TABLE}
specs = [
    ("br_pre", "d3", {u: {k: dval(raw[k]["desc"], BR, "Pre-IP TREP", u, "median") for k in raw} for u in ("cell", "donor")}),
    ("br_post", "d3", {u: {k: dval(raw[k]["desc"], BR, "Post-IP TREP", u, "median") for k in raw} for u in ("cell", "donor")}),
    ("br_rec", "d3", {u: {k: dval(raw[k]["desc"], BR, "Recovered", u, "median") for k in raw} for u in ("cell", "donor")}),
    ("br_exc", "d3", {u: {k: dval(raw[k]["desc"], BR, "Excluded", u, "median") for k in raw} for u in ("cell", "donor")}),
    ("fdr_br", "p", {u: {k: post_p(raw[k]["post"], BR, u, "Recovered", "Excluded") for k in raw} for u in ("cell", "donor")}),
    ("pos_med", "d2", {u: {k: dval(raw[k]["desc"], PP, "Recovered", u, "median") for k in raw} for u in ("cell", "donor")}),
    ("pos_q1", "d2", {u: {k: dval(raw[k]["desc"], PP, "Recovered", u, "q25") for k in raw} for u in ("cell", "donor")}),
    ("pos_q3", "d2", {u: {k: dval(raw[k]["desc"], PP, "Recovered", u, "q75") for k in raw} for u in ("cell", "donor")}),
    ("fdr_pos", "p", {u: {k: omn_p(raw[k]["omn"], PP, u) for k in raw} for u in ("cell", "donor")}),
    ("ne_pre", "d2", {c: {k: nb_val(raw[k]["nb"], "Pre-IP TREP", c) for k in raw} for c in ("mean_enrichment", "median_enrichment")}),
    ("ne_post", "d2", {c: {k: nb_val(raw[k]["nb"], "Post-IP TREP", c) for k in raw} for c in ("mean_enrichment", "median_enrichment")}),
    ("kb_pre", "d2", {"rejection_rate": {k: nb_val(raw[k]["kb"], "Pre-IP TREP", "rejection_rate") for k in raw}}),
    ("kb_post", "d2", {"rejection_rate": {k: nb_val(raw[k]["kb"], "Post-IP TREP", "rejection_rate") for k in raw}}),
    ("pct_si", "d1", {"pct_fail_Si_any": {k: float(raw[k]["cut"]["pct_fail_Si_any"]) for k in raw}}),
    ("pct_ki", "d1", {"pct_fail_Ki_any": {k: float(raw[k]["cut"]["pct_fail_Ki_any"]) for k in raw}}),
]
for col, kind, alts in specs:
    name, vals, ok_src = pick(alts, kind, col)
    choice[col] = name
    show("Figure 8C %-8s reproduced in all 8 datasets (source: %s)" % (col, name), ok_src, True, fig=True)
    if name is None:
        die("Figure 8C %s: no source has a value in every dataset" % col)
    for k in FIG8_TABLE:
        F8[k][col] = vals.get(k)
for k in FIG8_TABLE:
    if all(F8[k].get(c) is not None for c in ("br_pre", "br_post", "br_rec")):
        F8[k]["ratio"] = F8[k]["br_rec"] / ((F8[k]["br_pre"] + F8[k]["br_post"]) / 2)
        show("Figure 8C ratio %s %s" % k, round(F8[k]["ratio"], 1), FIG8_TABLE[k][F8_COLS.index("ratio")], tol=0.051, fig=True)

J_rows = []
for (tag, ds), v in F8.items():
    J_rows.append(dict(
        arm=tag, region_set=ds,
        bridge_median_PreIP_TREP=v.get("br_pre"), bridge_median_PostIP_TREP=v.get("br_post"),
        bridge_median_Recovered=v.get("br_rec"), bridge_median_Excluded=v.get("br_exc"),
        bridge_FDR_Recovered_vs_Excluded=v.get("fdr_br"),
        bridge_Recovered_over_mean_TREP=v.get("ratio"),
        pct_Recovered_fail_Si=v.get("pct_si"), pct_Recovered_fail_Ki=v.get("pct_ki"),
        position_Recovered_median=v.get("pos_med"), position_Recovered_q25=v.get("pos_q1"),
        position_Recovered_q75=v.get("pos_q3"), position_trend_FDR=v.get("fdr_pos"),
        neighbor_enrichment_vs_PreIP_TREP=v.get("ne_pre"), neighbor_enrichment_vs_PostIP_TREP=v.get("ne_post"),
        kBET_rejection_vs_PreIP_TREP=v.get("kb_pre"), kBET_rejection_vs_PostIP_TREP=v.get("kb_post")))
J_fig8 = pd.DataFrame(J_rows)
F8_SOURCE = "; ".join("%s: %s" % (c, choice[c]) for c in choice)

# Full test rows behind Figure 8C (both units), for the two statistics in the table
K_parts = []
for (tag, ds), r in raw.items():
    o = r["omn"][r["omn"]["statistic"].isin([BR, PP])].copy()
    o.insert(0, "arm", tag)
    K_parts.append(o.assign(test_level="omnibus and trend"))
    p = r["post"][r["post"]["statistic"].isin([BR, PP])].copy()
    p.insert(0, "dataset", ds)
    p.insert(0, "arm", tag)
    K_parts.append(p.assign(test_level="pairwise"))
K_tests8 = pd.concat(K_parts, ignore_index=True)

#####
# Figure 7: the mechanism workbook (M4 NeoC)
#####
mech = sorted(glob.glob(os.path.join(ROOT, "Part10f_RSkfCV_Mechanism", "*",
                                     "*M4*NeoC*RSkfCV_Mechanism.xlsx")))
mech = [m for m in mech if "NonNeoC" not in os.path.basename(m)]
if len(mech) != 1:
    die("expected one M4 NeoC mechanism workbook, found %d" % len(mech))
mech = mech[0]
mt = sheet(mech, "Gate_Trajectory")
mf = sheet(mech, "Fold_RMSE")
mc = sheet(mech, "Candidates")
c4 = C_gate[(C_gate["arm"] == "M4") & (C_gate["region_set"] == "NeoC")].iloc[0]
a4 = A_cells[(A_cells["arm"] == "M4") & (A_cells["region_set"] == "NeoC")].reset_index(drop=True)
aud4 = sheet(p10e_file("M4", "Mic13_PVM3", "NeoC", "RSk_fCV_Audit"), "Fold_RMSE")

print("\nFigure 7 (mechanism workbook, M4 NeoC)")
med_a4 = float(np.median(mf["RMSE_A"]))
mm = mf.merge(aud4[["repeat_id", "fold_id", "RMSE_A", "RMSE_B"]], on=["repeat_id", "fold_id"],
              suffixes=("", "_audit"), validate="one_to_one")
show("folds matched to the Part 10e audit", len(mm), 200)
show("RMSE_A identical to the Part 10e audit", float(np.max(np.abs(mm["RMSE_A"] - mm["RMSE_A_audit"]))), 0.0, tol=1e-9)
show("RMSE_B iteration 55 identical to the audit", float(np.max(np.abs(mm["RMSE_B_final"] - mm["RMSE_B"]))), 0.0, tol=1e-9)
d1 = float(np.median(mf["RMSE_B_iter1"] - mf["RMSE_A"]))
dF = float(np.median(mf["RMSE_B_final"] - mf["RMSE_A"]))
bv4 = iqr7(mf["RMSE_A"].to_numpy(float))
show("dRMSE iteration 1 (figure 0.0099)", round(d1, 4), REF_FIG7["d_iter1"], tol=1e-9, fig=True)
show("dRMSE iteration 55 (figure 0.0044)", round(dF, 4), REF_FIG7["d_final"], tol=1e-9, fig=True)
show("BV (figure 0.0046)", round(bv4, 4), REF_FIG7["bv"], tol=1e-9, fig=True)
show("trajectory iterations", len(mt), REF_FIG7["n_iter"], fig=True)
show("trajectory p, iteration 1 (figure 1.40e-34)", float("%.2e" % mt["wilcox_p"].iloc[0]), REF_FIG7["p_iter1"], tol=1e-9, rel=True, fig=True)
show("trajectory p, iteration 55 (figure 1.39e-25)", float("%.2e" % mt["wilcox_p"].iloc[-1]), REF_FIG7["p_final"], tol=1e-9, rel=True, fig=True)
show("trajectory dRMSE ratio, iteration 1 = fold pairs", round(float(mt["dRMSE"].iloc[0]), 4), round(d1 / med_a4, 4), tol=1e-4)
show("trajectory last dRMSE ratio = gate sheet", round(float(mt["dRMSE"].iloc[-1]), 6), round(float(c4["dRMSE_final_ratio"]), 6), tol=2e-6)

# Candidate barcodes, matched on (u, v) among the M4 NeoC candidates
cand4 = a4[a4["stage1_group"] == "Candidate"]
show("candidates in the mechanism workbook", len(mc), len(cand4))
bc, dist = [], []
cu, cv_ = cand4["u"].to_numpy(), cand4["v"].to_numpy()
for uu, vv in zip(mc["u"], mc["v"]):
    d = np.abs(cu - uu) + np.abs(cv_ - vv)
    j = int(np.argmin(d))
    bc.append(cand4["barcode"].iloc[j])
    dist.append(float(d[j]))
show("candidates matched on (u, v) (max distance)", float(max(dist)), 0.0, tol=1e-6)
show("candidate barcodes unique", int(pd.Series(bc).duplicated().sum()), 0)
mc = mc.assign(barcode=bc)
stepmap = dict(zip(a4["barcode"], a4["removal_step"]))
st_mech = mc["removal_step"].to_numpy(float)
st_cells = np.array([stepmap[b] for b in bc], float)
same = (np.isnan(st_mech) & np.isnan(st_cells)) | (st_mech == st_cells)
show("candidate removal steps = A_cells", int((~same).sum()), 0)
show("candidates accepted (figure 298)", int(as_bool(mc["accepted"]).sum()), 298, fig=True)

D_traj = mt.rename(columns={"wilcox_p": "p", "dRMSE": "dRMSE_ratio", "BV": "BV_ratio"})
D_traj = D_traj.drop(columns=[c for c in ["alpha"] if c in D_traj.columns])
D_traj.insert(0, "region_set", "NeoC")
D_traj.insert(0, "arm", "M4")
D_traj["dRMSE"] = D_traj["dRMSE_ratio"] * med_a4
D_traj["BV"] = bv4
E_folds = mf[["repeat_id", "fold_id", "RMSE_A", "RMSE_B_iter1", "RMSE_B_final"]].rename(
    columns={"RMSE_B_iter1": "RMSE_B_iteration_1", "RMSE_B_final": "RMSE_B_iteration_55"})
E_folds.insert(0, "region_set", "NeoC")
E_folds.insert(0, "arm", "M4")
F_cand = mc[["rank_by_pava", "barcode", "u", "v", "resid_pava", "removal_step"]].rename(
    columns={"rank_by_pava": "rank", "resid_pava": "abs_v_minus_PAVA_u"})
F_cand["outcome"] = np.where(as_bool(mc["accepted"]), "Recovered", "Removed")
F_cand.insert(0, "region_set", "NeoC")
F_cand.insert(0, "arm", "M4")

g4 = G_tests[(G_tests["arm"] == "M4") & (G_tests["region_set"] == "NeoC")]
if "dCor" in g4.columns and g4["dCor"].notna().any():
    show("distance correlation, recovered cells (figure 0.557)", round(float(g4["dCor"].dropna().iloc[0]), 3), REF_FIG7["dcor"], tol=1e-9, fig=True)

n_hard = sum(1 for ok, fig, _ in CHECKS if not ok and not fig)
n_fig = sum(1 for ok, fig, _ in CHECKS if not ok and fig)
print("\n%d checks: %d consistency checks failed, %d figure checks failed" % (len(CHECKS), n_hard, n_fig))
for ok, fig, label in CHECKS:
    if not ok:
        print("  %s %s" % ("figure     " if fig else "consistency", label))
if n_hard or (n_fig and not ARGS.allow_figure_drift):
    die("checks failed; workbook not written")
if n_fig:
    print("Figure checks failed but --allow-figure-drift is set; the workbook is written "
          "with the values of this run.")

#####
# Contents and notes
#####
TITLE = "Supplementary Table S10. Copula monotone recovery filter (CMRF)."
DESC = ("Values behind Figures 5 to 8 and Supplementary Figures 11 and 12. Sheet A lists every "
        "GOI-positive cell once per arm and region set; sheets D to F hold the RSk-fCV gate of "
        "M4 NeoC shown in Figure 7.")
SRC_SK = "Part10e_CP10K_<tag>_CMRF_<ds>_SK_interaction_analysis.xlsx"
SRC_AU = "Part10e_CP10K_<tag>_CMRF_<ds>_RSk_fCV_Audit.xlsx"
SRC_ME = "Part10e_CP10K_M4_RSkfCV_NeoC_RSkfCV_Mechanism.xlsx"
SHEETS = [
    ("A_cells", A_cells, SRC_SK + " [Step4_Refined_Cells]; " + SRC_AU + " [Gate_Iteration]",
     "One row per GOI-positive cell: copula coordinates, SKCS-core group, spline training cells, "
     "g(u), |v - g(u)|, delta_adapt(u), Stage 1 group, Stage 2 removal step, CMRF class, final "
     "layer and the IP side of Recovered cells; for Figure 8, the copula score coordinates, the "
     "refit group and the 2D UMAP coordinates. Figures 6 and 8; Figure 7Q; Supplementary Figure 12."),
    ("B_band", B_band, SRC_SK + " [Step4_Boundary_Params, Step4_Spine_PRSS]",
     "Band and spline parameters and class counts per arm and region set. Figure 6; "
     "Supplementary Figures 11 and 12."),
    ("C_gate", C_gate, SRC_SK + " [Step4_CMRF_Stage2_Gate]; " + SRC_AU + " [Fold_RMSE, Wilcoxon_Summary]",
     "RSk-fCV gate summary per arm and region set: candidates, iterations, removed, recovered, "
     "median RMSE_A, BV, final dRMSE (absolute and as ratios), final Wilcoxon p and the closing "
     "criterion. Figure 5; Figure 7; Supplementary Figure 12."),
    ("D_gate_trajectory", D_traj, SRC_ME + " [Gate_Trajectory]",
     "Wilcoxon p, dRMSE and active candidates at each of the 55 iterations, M4 NeoC. Figure 7I."),
    ("E_fold_RMSE", E_folds, SRC_ME + " [Fold_RMSE]",
     "RMSE_A and RMSE_B of the 200 folds (40 repeats x 5) at iterations 1 and 55, M4 NeoC. "
     "Figure 7D-F and K-M."),
    ("F_candidates", F_cand, SRC_ME + " [Candidates]",
     "The 352 Stage 1 candidates ranked by |v - PAVA(u)|, with removal step and outcome, M4 NeoC. "
     "Figure 7G, H and N."),
    ("G_tests", G_tests, "Part10e_CP10K_<tag>_CMRF_Verification.xlsx [<ds>_StatTests]",
     "Tests on the CMRF classes and the recovered cells, per arm and region set: the RSk-fCV "
     "comparison, |v - g(u)| across classes (Kruskal-Wallis and Jonckheere-Terpstra), Recovered "
     "against Excluded, and the dependence of u and v in Recovered cells (Hoeffding's D, distance "
     "correlation). Figure 7O and P."),
    ("H_resid_by_class", H_resid, "Part10e_CP10K_<tag>_CMRF_<ds>_VP_PlotData_T1_to_T7A.xlsx [T3_JT_ordered_trend]",
     "Median and quartiles of |v - g(u)| by CMRF class, per arm and region set. Figure 7P."),
    ("I_Fig8_summary", J_fig8,
     "Part10e_CP10K_<tag>_V2_Validation2_Summary.xlsx [<ds>_Fig11_Descriptives, _Fig11_PostHoc, "
     "Fig11_Omnibus_BH, _Neighbour_Proportions, _kBET, _Fig11_SKCS_Cuts]",
     "The Figure 8C table per arm and region set: bridge distances, the Recovered/TREP ratio, "
     "gates failed, position on the pre-IP to post-IP axis, neighbor enrichment and kBET. "
     "Sources chosen by matching the figure: " + F8_SOURCE + "."),
    ("J_Fig8_tests", K_tests8, "Part10e_CP10K_<tag>_V2_Validation2_Summary.xlsx [Fig11_Omnibus_BH, <ds>_Fig11_PostHoc]",
     "Omnibus, ordered-trend and pairwise tests of the bridge distance and the position, at the "
     "cell and donor units, with BH-adjusted p. Figure 8C."),
]
NOTES = [
    "1) Copula coordinates and arms as in Supplementary Table S7: u and v are the ranks of the "
    "GOI and the composite divided by n + 1 within each arm and region set. SKCS_group gives the "
    "SKCS-core cells their pre-IP or post-IP and TREP or non-TREP label; other cells are "
    "Non-SKCS-core. spline_training marks the TREP cells of the SKCS-core, on which g(u) was fitted.",
    "2) Band. delta_adapt(u) = delta_base[1 + kappa exp(-(u - IP_u)^2/(2 bw^2))]; delta_max = "
    "delta_base(1 + kappa). within_delta_adapt is the Stage 1 admission; within_delta_base marks the "
    "flat band. stage1_group: OriRef (SKCS-core, admitted), Original only (SKCS-core, not admitted), "
    "Candidate (non-SKCS, admitted), Not admitted (non-SKCS).",
    "3) Gate. removal_step is the iteration at which a candidate was removed (blank if never). "
    "CMRF_class: Original + Refined, Recovered (candidates not removed), Original only, Excluded "
    "(non-SKCS cells not admitted or removed). final_layer = Original + Refined or Recovered. "
    "recovered_IP_side is pre-IP or post-IP on the composite at IP_y.",
    "4) dRMSE = median(RMSE_B - RMSE_A) and BV = IQR(RMSE_A) over the 200 fold pairs, in RMSE units; "
    "the _ratio columns divide both by the median RMSE_A, which is how the source workbooks store "
    "them. The gate removes a candidate while p < 0.05 and dRMSE > BV.",
    "5) The final Wilcoxon p is read from Wilcoxon_Summary, because the gate summary sheets store it "
    "rounded to 0. RMSE_A is computed once from folds over OriRef alone and is the same at every "
    "iteration; RMSE_B is recomputed at each iteration on folds over OriRef plus the remaining "
    "candidates. The pairs share repeat and fold index.",
    "6) Candidate barcodes in sheet F were matched to sheet A on u and v (the mechanism workbook "
    "records candidate positions only); every match is exact and unique.",
    "7) Figure 7C and J draw isotonic (PAVA) fits of v on u for Curve A (OriRef) and Curve B (OriRef "
    "plus candidates), with each cell's residual; they can be refitted from sheet A. The right-hand "
    "panels draw both curves on one fold for comparison.",
    "8) Sheet G. The source workbooks round p to six decimals. The Wilcoxon p of the gate is replaced "
    "by its exact value from sheet C, and the Kruskal-Wallis p by the chi-squared tail of the stored "
    "statistic (df = number of classes - 1); other p-values stored as 0 are labelled. "
    "Jonckheere-Terpstra (10,000 permutations), Hoeffding's D (9,999 permutations) and distance "
    "correlation (999 permutations) give permutation p-values whose smallest attainable values are "
    "1x10^-4, 1x10^-4 and 1x10^-3.",
    "9) Figure 8 columns of sheet A. u_SK = rank(S_i)/(n + 1) and v_SK = rank(K_i)/(n + 1); "
    "refit_group_Fig8 is the four-group label of the GAM refitted on the final Module (Validation 1), "
    "Recovered, or Non-Module; UMAP_1 and UMAP_2 are the 2D Monocle3 UMAP of Validation 2, built from "
    "raw counts. Neighbor enrichment and kBET were computed in the 3D UMAP. The UMAP coordinates and "
    "u_SK, v_SK are copied from the Validation RDS files by build_SuppTable_S10_step0_Fig8_cells.r.",
    "10) build_SuppTable_S10.py prints every check against the source workbooks and the values in "
    "the figures, and writes the workbook only if all checks pass.",
]

#####
# Write (write-only mode)
#####
ARIAL = Font(name="Arial", size=10)
HEAD = Font(name="Arial", size=10, bold=True)
RAW = {"A_cells", "E_fold_RMSE", "F_candidates", "D_gate_trajectory", "J_Fig8_tests"}


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
    if v is None or v is pd.NA or (isinstance(v, float) and np.isnan(v)):
        return None
    return v


def tidy(df):
    df = df.copy()
    for c in df.columns:
        if df[c].dtype.kind == "f" and not any(k in c for k in ("p_value", "p_final", "FDR", "_p", "p_BH")) \
                and c not in ("p",):
            df[c] = df[c].round(6)
    return df


wb = Workbook(write_only=True)
ws = wb.create_sheet("Contents")
for col, wdt in zip("ABCD", (22, 60, 8, 90)):
    ws.column_dimensions[col].width = wdt
ws.freeze_panes = "A5"
ws.append([cell(ws, TITLE, HEAD)])
ws.append([cell(ws, DESC)])
ws.append([])
ws.append([cell(ws, h, HEAD) for h in ("Sheet", "Source", "Rows", "Description")])
for name, df, src, desc in SHEETS:
    ws.append([cell(ws, name, HEAD), cell(ws, src), cell(ws, "=COUNTA(%s!A:A)-1" % name),
               cell(ws, desc, wrap=True)])
ws.append([])
ws.append([cell(ws, "Notes", HEAD)])
for n_ in NOTES:
    ws.append([cell(ws, n_)])

for name, df, *_ in SHEETS:
    if name not in RAW:
        df = tidy(df)
    sh = wb.create_sheet(name)
    for i, col in enumerate(df.columns, start=1):
        vals = [len(str(x)) for x in df[col].head(200)]
        sh.column_dimensions[get_column_letter(i)].width = min(max([len(str(col))] + vals) + 2, 46)
    sh.freeze_panes = "A2"
    sh.append([cell(sh, str(c), HEAD) for c in df.columns])
    for row in df.itertuples(index=False):
        sh.append([cell(sh, py(x)) for x in row])

wb.save(OUT)
print("\nsaved", OUT)
for name, df, *_ in SHEETS:
    print("  %-20s %d rows" % (name, len(df)))
print("Log:", LOG)
SAVED = True
