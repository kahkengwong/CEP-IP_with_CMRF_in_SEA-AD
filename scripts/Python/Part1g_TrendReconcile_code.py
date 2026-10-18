# =============================================================================
# Part1g_TrendReconcile_code.py
# SEA-AD CEP-IP/CMRF, Part 1g, detectability trend on two abundance axes
#   (Figure 1G, Table S1)
#   Fits rho(MT-CO3, gene) against log10 mean counts over the mitochondrial
#   genes, with mean counts taken over all Micro-PVM_3-SEAAD nuclei and over
#   the analysed nuclei (MT-CO3 > 0), and reports the ACTB residual for both.
#
# Usage (from the SEA-AD folder, or set SEAAD_ROOT)
#   python3 Part1g_TrendReconcile_code.py
#
# Outputs (Part1_PreChecks/): S1g_trend_reconcile.csv, S1g_trend_fits.csv and
#   report.
# =============================================================================

import os

import numpy as np
import pandas as pd
import anndata as ad
from scipy.stats import spearmanr

# == Config ==================================================================
H5AD = os.environ.get("SEAAD_ROOT", os.getcwd()) + "/SEAAD_Immune_HIP_MEC_LEC_ITG_MTG_FI_STG_DFC_AnG_V1C_RNAseq_final-nuclei.2026-06-22.h5ad"
OUT_DIR = os.environ.get("SEAAD_ROOT", os.getcwd()) + "/Part1_PreChecks"

TARGET_STYPE = "Micro-PVM_3-SEAAD"
GOI = "MT-CO3"
REF_GENE = "ACTB"
COUNT_LAYER = "UMIs"          # raw counts; adata.X is ln(UP10K+1)
STYPE_COL_CANDIDATES = ["Supertype", "supertype"]

# Recorded values, checked against this run
RECORDED = {
    "all_PVM3":  dict(slope=+0.3542, intercept=+0.4573, r2=0.9877,
                      resid=-0.5064, ratio=14.4),
    "analysed":  dict(slope=+0.352, intercept=+0.432, r2=0.9876,
                      resid=-0.495, ratio=13.9),
}

_report = []


def say(s=""):
    print(s)
    _report.append(s)


def hdr(s):
    say("")
    say("=" * 74)
    say(f"  {s}")
    say("=" * 74)


def find_col(obs, candidates):
    for c in candidates:
        if c in obs.columns:
            return c
    raise KeyError(f"none of {candidates} found in obs; columns are {list(obs.columns)}")


def fit(x, y):
    """OLS of y on x. Returns slope, intercept, R2, predictions."""
    slope, intercept = np.polyfit(x, y, 1)
    pred = intercept + slope * x
    ss_res = ((y - pred) ** 2).sum()
    ss_tot = ((y - y.mean()) ** 2).sum()
    return slope, intercept, 1.0 - ss_res / ss_tot, pred


# == Load ====================================================================
hdr("Loading")

adata = ad.read_h5ad(H5AD, backed="r")
say(f"  file  : {os.path.basename(H5AD)}")
say(f"  shape : {adata.shape}")

stype_col = find_col(adata.obs, STYPE_COL_CANDIDATES)
is_target = (adata.obs[stype_col].astype(str).values == TARGET_STYPE)
say(f"  {TARGET_STYPE} : {int(is_target.sum()):,} nuclei")

mt_all = [g for g in adata.var_names if str(g).startswith("MT-")]
genes = [GOI, REF_GENE] + [g for g in mt_all if g != GOI]
genes = [g for g in genes if g in adata.var_names]
say(f"  genes extracted : {len(genes)}  ({len(mt_all)} mitochondrial)")

idx = np.where(is_target)[0]
sub = adata[idx, genes].to_memory()

Cn, E = {}, {}
for j, g in enumerate(genes):
    counts = sub.layers[COUNT_LAYER][:, j]
    norm = sub.X[:, j]
    Cn[g] = np.asarray(counts.todense()).ravel() if hasattr(counts, "todense") else np.asarray(counts).ravel()
    E[g] = np.asarray(norm.todense()).ravel() if hasattr(norm, "todense") else np.asarray(norm).ravel()

adata.file.close()

ANALYSED = Cn[GOI] > 0                      # index into the target subset
say(f"  analysed ({GOI} > 0)  : {int(ANALYSED.sum()):,} nuclei")
say(f"  excluded ({GOI} == 0) : {int((~ANALYSED).sum()):,} nuclei")

# == Per gene ================================================================
hdr("Per-gene table, both abundance axes")

rows = []
for g in genes:
    if g == GOI:
        rho = np.nan
    else:
        rho, _ = spearmanr(E[GOI][ANALYSED], E[g][ANALYSED])
    rows.append(dict(
        gene=g,
        is_mito=g.startswith("MT-"),
        rho_analysed=rho,
        mean_counts_all=float(Cn[g].mean()),
        mean_counts_analysed=float(Cn[g][ANALYSED].mean()),
        detect_pct_all=100.0 * float((Cn[g] > 0).mean()),
        detect_pct_analysed=100.0 * float((Cn[g][ANALYSED] > 0).mean()),
    ))
inv = pd.DataFrame(rows)
inv["ratio_analysed_over_all"] = inv.mean_counts_analysed / inv.mean_counts_all

say(inv[["gene", "rho_analysed", "mean_counts_all", "mean_counts_analysed",
         "ratio_analysed_over_all", "detect_pct_all",
         "detect_pct_analysed"]].to_string(index=False,
                                           float_format=lambda v: f"{v:9.4f}"))

# == Both fits ===============================================================
hdr("The two fits")

mito = inv[inv.is_mito & inv.rho_analysed.notna() & (inv.mean_counts_all > 0)]
ref = inv[inv.gene == REF_GENE]
if len(ref) != 1:
    raise RuntimeError(f"{REF_GENE} not found in the object")
ref = ref.iloc[0]

results = {}
for label, col in [("all_PVM3", "mean_counts_all"),
                   ("analysed", "mean_counts_analysed")]:
    x = np.log10(mito[col].values)
    y = mito.rho_analysed.values
    slope, intercept, r2, pred = fit(x, y)

    xr = float(np.log10(ref[col]))
    yr = float(ref.rho_analysed)
    pred_ref = intercept + slope * xr
    resid = yr - pred_ref
    worst = float(np.max(np.abs(y - pred)))
    worst_gene = mito.gene.values[int(np.argmax(np.abs(y - pred)))]

    results[label] = dict(slope=slope, intercept=intercept, r2=r2,
                          xr=xr, yr=yr, pred_ref=pred_ref, resid=resid,
                          worst=worst, worst_gene=worst_gene,
                          ratio=abs(resid) / worst, n_genes=len(x))

say(f"  fitted on {results['all_PVM3']['n_genes']} non-GOI mitochondrial genes;")
say("")
say(f"  {'':<34}{'all PVM_3':>14}{'analysed':>14}{'difference':>14}")
say("  " + "-" * 74)
for key, fmt in [("slope", "{:+.4f}"), ("intercept", "{:+.4f}"),
                 ("r2", "{:.4f}"), ("xr", "{:+.4f}"),
                 ("pred_ref", "{:+.4f}"), ("resid", "{:+.4f}"),
                 ("worst", "{:+.4f}"), ("ratio", "{:.2f}")]:
    a, b = results["all_PVM3"][key], results["analysed"][key]
    labels = {"xr": f"log10 mean counts, {REF_GENE}",
              "pred_ref": f"{REF_GENE} predicted rho",
              "resid": f"{REF_GENE} residual",
              "worst": "worst mito residual",
              "ratio": "residual / worst mito residual"}
    say(f"  {labels.get(key, key):<34}{fmt.format(a):>14}{fmt.format(b):>14}"
        f"{fmt.format(b - a):>14}")
say("")
say(f"  {REF_GENE} observed rho : {results['all_PVM3']['yr']:+.4f}  (same in both)")
say(f"  worst mito residual carried by: "
    f"{results['all_PVM3']['worst_gene']} / {results['analysed']['worst_gene']}")

# Detection-rate variant
say("")
for label, col in [("all_PVM3", "detect_pct_all"),
                   ("analysed", "detect_pct_analysed")]:
    xd = mito[col].values
    if np.ptp(xd) == 0:
        say(f"  vs detection% ({label:<9}): not defined "
            f"(every gene is 100% by construction)")
        continue
    _, _, r2d, _ = fit(xd, mito.rho_analysed.values)
    say(f"  vs detection% ({label:<9}): R2 {r2d:.4f}")

# == Check against the recorded values =======================================
hdr("Against the recorded values")

for label in ("all_PVM3", "analysed"):
    say(f"  {label}")
    for key in ("slope", "intercept", "r2", "resid", "ratio"):
        got, exp = results[label][key], RECORDED[label][key]
        # Tolerances follow the rounding of the recorded values
        tol = {"slope": 5e-4, "intercept": 5e-4, "r2": 5e-4,
               "resid": 5e-4, "ratio": 5e-2}[key]
        flag = "PASS" if abs(got - exp) <= tol else "DIFFERS"
        say(f"    {key:<12} run {got:+.4f}   recorded {exp:+.4f}   {flag}")
    say("")

# == Save ====================================================================
hdr("Saving")
os.makedirs(OUT_DIR, exist_ok=True)
path = os.path.join(OUT_DIR, "S1g_trend_reconcile.csv")
inv.to_csv(path, index=False)
say(f"  saved: {path}")

fits = pd.DataFrame(results).T.reset_index().rename(columns={"index": "abundance_axis"})
path2 = os.path.join(OUT_DIR, "S1g_trend_fits.csv")
fits.to_csv(path2, index=False)
say(f"  saved: {path2}")

with open(os.path.join(OUT_DIR, "S1g_trend_reconcile_report.txt"), "w",
          encoding="utf-8") as f:
    f.write("\n".join(_report))
say(f"  saved: {OUT_DIR}/S1g_trend_reconcile_report.txt")

say("")
say("=" * 74)
say("  PART 1g COMPLETE")
say("=" * 74)
