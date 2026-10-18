# =============================================================================
# Part1c_WithinDonor_code.py
# SEA-AD CEP-IP/CMRF, Part 1c, within-donor decomposition of the module
#   Splits the pooled rho(MT-CO3, cDFG) in GOI-positive PVM_3 cells into a
#   within-donor rho (median and IQR over donors with at least 30 cells) and a
#   between-donor rho (donor means), each with and without depth partialled
#   out. A simulation gives the IQR expected from estimation noise alone.
#
# Usage (from the SEA-AD folder, or set SEAAD_ROOT)
#   python3 Part1c_WithinDonor_code.py
#
# Outputs (Part1_PreChecks/): per-donor rho (CSV) and report.
# =============================================================================

import os, gc, warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import spearmanr, rankdata
import scanpy as sc

np.random.seed(42)
sc.settings.verbosity = 1

os.chdir(os.environ.get("SEAAD_ROOT", os.getcwd()))

FULL_H5AD = ("SEAAD_Immune_HIP_MEC_LEC_ITG_MTG_FI_STG_DFC_AnG_V1C_"
             "RNAseq_final-nuclei.2026-06-22.h5ad")
OUT_DIR = "Part1_PreChecks"
os.makedirs(OUT_DIR, exist_ok=True)

# == Module genes (as in R_config.R) =========================================
GOI = "MT-CO3"
DFG = ["MT-ATP6", "MT-CO2", "MT-ND3", "MT-ND4",
       "MT-CYB", "MT-ND1", "MT-ND2", "MT-CO1"]
MODULE1 = [GOI] + DFG

# == obs column defaults (resolved against the file below) ==================
stype_col  = "Supertype"
region_col = "Brain Region"
donor_col  = "Donor ID"
adnc_col   = "Overall AD neuropathological Change"
TARGET_STYPE = "Micro-PVM_3-SEAAD"

# == Parameters ==============================================================
# Minimum GOI-positive PVM_3 cells per donor, pooled over libraries and regions.
MIN_CELLS_DONOR = 30

# Simulation replicates for the noise floor
N_SIM = 500

# ADNC rank map
ADNC_RANK = {"Not AD": 0, "Low": 1, "Intermediate": 2, "High": 3}


# == Utility functions =======================================================
def find_col(cols, keywords):
    """Resolve one obs column name."""
    cols = list(cols)
    for kw in keywords:
        for c in cols:
            if c == kw:
                return c
    for kw in keywords:
        hits = [c for c in cols if kw.lower() in str(c).lower()]
        if hits:
            return sorted(hits, key=len)[0]
    return None


def partial_spearman(x, y, z):
    """Spearman partial correlation of x,y given covariate(s) z."""
    z = np.atleast_2d(np.asarray(z, dtype=float))
    if z.shape[0] > z.shape[1]:
        z = z.T
    R = np.column_stack([rankdata(x), rankdata(y)] + [rankdata(c) for c in z])
    sd = R.std(0); sd[sd == 0] = 1.0
    R = (R - R.mean(0)) / sd
    Z = np.column_stack([np.ones(len(R)), R[:, 2:]])
    ex = R[:, 0] - Z @ np.linalg.lstsq(Z, R[:, 0], rcond=None)[0]
    ey = R[:, 1] - Z @ np.linalg.lstsq(Z, R[:, 1], rcond=None)[0]
    return spearmanr(ex, ey)


_report_1c = []
def say_c(*a):
    line = " ".join(str(x) for x in a)
    print(line)
    _report_1c.append(line)

def hdr_c(t):
    say_c("")
    say_c("=" * 74)
    say_c("  " + t)
    say_c("=" * 74)


# =============================================================================
# SECTION C0 - Load the h5ad and extract the module genes
# =============================================================================
hdr_c("Part 1c - Within-donor decomposition of rho(MT-CO3, cDFG)")

say_c(f"  MIN_CELLS_DONOR : {MIN_CELLS_DONOR}")
say_c(f"  N_SIM           : {N_SIM}")
say_c(f"  GOI             : {GOI}")
say_c(f"  DFG             : {DFG}")
say_c(f"  TARGET_STYPE    : {TARGET_STYPE}")

hdr_c("C0 - Loading h5ad")

ad = sc.read_h5ad(FULL_H5AD)
say_c(f"  shape : {ad.shape}")

# Resolve obs columns
_oc = list(ad.obs.columns)
stype_col  = find_col(_oc, [stype_col, "Supertype", "supertype"]) or stype_col
region_col = find_col(_oc, [region_col, "Brain Region", "brain_region",
                             "Region", "region"]) or region_col
donor_col  = find_col(_oc, [donor_col, "Donor ID", "donor_id", "Donor"]) or donor_col
adnc_col   = find_col(_oc, [adnc_col, "AD neuropath", "ADNC"]) or adnc_col

say_c(f"  stype_col  : {stype_col}")
say_c(f"  donor_col  : {donor_col}")
say_c(f"  adnc_col   : {adnc_col}")

for col_name in [stype_col, donor_col]:
    if col_name not in _oc:
        raise RuntimeError(f"Critical column '{col_name}' not in obs. "
                           f"Check column names.")

# Identify count layer
COUNT_LAYER = None
for cand in ["UMIs", "counts", "raw_counts"]:
    if ad.layers and cand in ad.layers:
        COUNT_LAYER = cand
        break
say_c(f"  count layer: {COUNT_LAYER if COUNT_LAYER else 'NONE - using X'}")

# Extract only the genes this script needs
var_names = list(ad.var_names)
need = [g for g in dict.fromkeys(MODULE1) if g in var_names]
need_idx = [var_names.index(g) for g in need]
say_c(f"  extracting {len(need)}/{len(MODULE1)} module genes ...")

Craw = ad.layers[COUNT_LAYER] if COUNT_LAYER else ad.X
Xn   = ad.X

total_counts = np.asarray(Craw.sum(axis=1)).ravel().astype(float)

Cblk = Craw[:, need_idx]
Cblk = np.asarray(Cblk.toarray() if sp.issparse(Cblk) else Cblk, dtype=float)
Eblk = Xn[:, need_idx]
Eblk = np.asarray(Eblk.toarray() if sp.issparse(Eblk) else Eblk, dtype=float)
Cn = {g: Cblk[:, i] for i, g in enumerate(need)}
E  = {g: Eblk[:, i] for i, g in enumerate(need)}

obs = ad.obs.copy()
del ad, Cblk, Eblk, Craw, Xn; gc.collect()
say_c("  done")

# Build masks
is_target = (obs[stype_col].astype(str).values == TARGET_STYPE)
sel = is_target & (Cn[GOI] > 0)
say_c(f"  {TARGET_STYPE} : {int(is_target.sum()):,} cells")
say_c(f"  sel (GOI > 0)        : {int(sel.sum()):,} cells")


# =============================================================================
# SECTION C1 - Per-donor cell table
# =============================================================================
hdr_c("C1 - Donor table: GOI-positive PVM_3 cells per donor")

donor_arr = obs[donor_col].astype(str).values
donors_sel = donor_arr[sel]

# Count GOI+ PVM_3 cells per donor (pooled across libraries/regions)
donor_cell_counts = pd.Series(donors_sel).value_counts().sort_values(ascending=False)
n_donors_total = len(donor_cell_counts)

qualifying = donor_cell_counts[donor_cell_counts >= MIN_CELLS_DONOR]
excluded   = donor_cell_counts[donor_cell_counts < MIN_CELLS_DONOR]
n_qual = len(qualifying)
n_excl = len(excluded)

say_c(f"  Total donors with >=1 GOI+ PVM_3 cell : {n_donors_total}")
say_c(f"  Qualifying (>= {MIN_CELLS_DONOR} cells)              : {n_qual}")
say_c(f"  Excluded   (<  {MIN_CELLS_DONOR} cells)              : {n_excl}")
say_c(f"  Cells in qualifying donors             : "
      f"{int(qualifying.sum()):,} / {int(sel.sum()):,}")
say_c("")
say_c(f"  Cell-count distribution (qualifying donors):")
say_c(f"    min={qualifying.min()}  q25={int(np.percentile(qualifying, 25))}  "
      f"median={int(np.median(qualifying))}  "
      f"q75={int(np.percentile(qualifying, 75))}  max={qualifying.max()}")

if n_excl > 0:
    say_c(f"  Cell-count distribution (excluded donors):")
    say_c(f"    min={excluded.min()}  max={excluded.max()}")


# =============================================================================
# SECTION C2 - ADNC distribution of included and excluded donors
# =============================================================================
hdr_c("C2 - ADNC distribution of included vs excluded donors")

_has_adnc = False
try:
    if adnc_col is not None and adnc_col in obs.columns:
        _adnc_raw = obs[adnc_col].astype(str)
        _adnc_mapped = _adnc_raw.map(ADNC_RANK)
        _donor_adnc = pd.DataFrame({
            "donor": donor_arr,
            "adnc": _adnc_mapped.values
        }).groupby("donor")["adnc"].first().dropna()
        _has_adnc = len(_donor_adnc) > 0
except Exception:
    pass

if _has_adnc:
    for label, donor_set in [("Included", qualifying.index),
                              ("Excluded", excluded.index)]:
        _matched = _donor_adnc.reindex(donor_set).dropna()
        if len(_matched) == 0:
            say_c(f"  {label:<10}: no ADNC data")
            continue
        _vc = _matched.value_counts().sort_index()
        parts = []
        for rank_val, rank_label in enumerate(["NotAD", "Low", "Interm", "High"]):
            parts.append(f"{rank_label}={int(_vc.get(rank_val, 0))}")
        say_c(f"  {label:<10} (n={len(_matched)}): {', '.join(parts)}")
else:
    say_c("  ADNC column not available.")


# =============================================================================
# SECTION C3 - Within-donor rho for each qualifying donor
# =============================================================================
hdr_c("C3 - Within-donor rho(MT-CO3, cDFG)")

dfg_have = [g for g in DFG if g in E]
say_c(f"  cDFG built from {len(dfg_have)}/{len(DFG)} DFGs: {dfg_have}")

# Full arrays
_goi_all = E[GOI]
_dfg_all = np.vstack([E[g] for g in dfg_have])   # (n_dfg, n_cells)
_tc_all  = total_counts

rows_1c = []
for donor_id in qualifying.index:
    d_mask = sel & (donor_arr == donor_id)
    n_i = int(d_mask.sum())

    goi_d  = _goi_all[d_mask]
    cdfg_d = _dfg_all[:, d_mask].mean(axis=0)
    tc_d   = _tc_all[d_mask]

    # Within-donor rho (unconditional)
    if goi_d.std() == 0 or cdfg_d.std() == 0:
        rho_w, rho_w_depth = np.nan, np.nan
    else:
        rho_w, _ = spearmanr(goi_d, cdfg_d)
        # Within-donor partial rho (conditioned on depth)
        try:
            rho_w_depth, _ = partial_spearman(goi_d, cdfg_d, tc_d)
        except Exception:
            rho_w_depth = np.nan

    # Donor-level means (for between-donor rho later)
    mean_goi  = float(goi_d.mean())
    mean_cdfg = float(cdfg_d.mean())
    mean_tc   = float(tc_d.mean())

    # ADNC
    adnc_val = np.nan
    if _has_adnc and donor_id in _donor_adnc.index:
        adnc_val = float(_donor_adnc[donor_id])

    rows_1c.append({
        "donor": donor_id,
        "n_cells": n_i,
        "rho_within": round(rho_w, 4) if np.isfinite(rho_w) else np.nan,
        "rho_within_depth": round(rho_w_depth, 4) if np.isfinite(rho_w_depth) else np.nan,
        "mean_GOI": round(mean_goi, 4),
        "mean_cDFG": round(mean_cdfg, 4),
        "mean_depth": round(mean_tc, 1),
        "adnc": adnc_val,
    })

df_1c = pd.DataFrame(rows_1c)
say_c(f"  Computed rho for {len(df_1c)} donors")


# =============================================================================
# SECTION C4 - Pooled, within-donor and between-donor rho
# =============================================================================
hdr_c("C4 - Three-number decomposition: pooled / within-donor / between-donor")

# 1. Pooled rho
goi_pooled  = _goi_all[sel]
cdfg_pooled = _dfg_all[:, sel].mean(axis=0)
tc_pooled   = _tc_all[sel]

rho_pooled, _ = spearmanr(goi_pooled, cdfg_pooled)
rho_pooled_depth, _ = partial_spearman(goi_pooled, cdfg_pooled, tc_pooled)

say_c(f"  1. POOLED rho(MT-CO3, cDFG)")
say_c(f"       unconditional : {rho_pooled:+.4f}  (n={int(sel.sum()):,} cells)")
say_c(f"       | depth       : {rho_pooled_depth:+.4f}")

# 2. Within-donor rho (median and IQR)
_rho_w = df_1c["rho_within"].dropna()
_rho_wd = df_1c["rho_within_depth"].dropna()

say_c("")
say_c(f"  2. WITHIN-DONOR rho  (n={len(_rho_w)} donors, "
      f">={MIN_CELLS_DONOR} GOI+ PVM_3 cells each)")
say_c(f"       unconditional : median={_rho_w.median():+.4f}  "
      f"IQR=[{_rho_w.quantile(0.25):+.4f}, {_rho_w.quantile(0.75):+.4f}]  "
      f"range=[{_rho_w.min():+.4f}, {_rho_w.max():+.4f}]")
say_c(f"       | depth       : median={_rho_wd.median():+.4f}  "
      f"IQR=[{_rho_wd.quantile(0.25):+.4f}, {_rho_wd.quantile(0.75):+.4f}]")
say_c(f"       fraction > 0  : {(_rho_w > 0).mean():.3f}  "
      f"({int((_rho_w > 0).sum())}/{len(_rho_w)})")

# 3. Between-donor rho
_ok_between = df_1c.dropna(subset=["mean_GOI", "mean_cDFG"])
rho_between, p_between = spearmanr(_ok_between["mean_GOI"],
                                    _ok_between["mean_cDFG"])
rho_between_depth, p_between_depth = partial_spearman(
    _ok_between["mean_GOI"].values,
    _ok_between["mean_cDFG"].values,
    _ok_between["mean_depth"].values)

say_c("")
say_c(f"  3. BETWEEN-DONOR rho  (donor-mean GOI vs donor-mean cDFG)")
say_c(f"       unconditional : {rho_between:+.4f}  p={p_between:.4g}  "
      f"(n={len(_ok_between)} donors)")
say_c(f"       | depth       : {rho_between_depth:+.4f}  "
      f"p={p_between_depth:.4g}")

# Summary table
say_c("")
say_c("  +--------------------------------------------------------------+")
say_c(f"  |  Pooled rho        : {rho_pooled:+.4f}  "
      f"(depth-partialled: {rho_pooled_depth:+.4f})   |")
say_c(f"  |  Within-donor med  : {_rho_w.median():+.4f}  "
      f"(depth-partialled: {_rho_wd.median():+.4f})   |")
say_c(f"  |  Between-donor     : {rho_between:+.4f}  "
      f"(depth-partialled: {rho_between_depth:+.4f})   |")
say_c("  +--------------------------------------------------------------+")


# =============================================================================
# SECTION C5 - Noise-floor simulation
# =============================================================================
hdr_c("C5 - Noise-floor simulation: expected IQR under true rho = 0.80")

say_c(f"  {N_SIM} replicates per donor at each donor's n_i, true rho = 0.80")

TRUE_RHO = 0.80
np.random.seed(42)

sim_rhos = []
for _, row_sim in df_1c.iterrows():
    n_i = int(row_sim["n_cells"])
    mean = [0, 0]
    cov = [[1, TRUE_RHO], [TRUE_RHO, 1]]
    donor_sim_rhos = []
    for _ in range(N_SIM):
        xy = np.random.multivariate_normal(mean, cov, size=n_i)
        r_sim, _ = spearmanr(xy[:, 0], xy[:, 1])
        donor_sim_rhos.append(r_sim)
    sim_rhos.append(np.median(donor_sim_rhos))

sim_rhos = np.array(sim_rhos)
sim_q25 = np.percentile(sim_rhos, 25)
sim_q75 = np.percentile(sim_rhos, 75)
sim_iqr = sim_q75 - sim_q25

obs_iqr = float(_rho_w.quantile(0.75) - _rho_w.quantile(0.25))

say_c(f"  Simulated  : median={np.median(sim_rhos):+.4f}  "
      f"IQR=[{sim_q25:+.4f}, {sim_q75:+.4f}]  width={sim_iqr:.4f}")
say_c(f"  Observed   : median={_rho_w.median():+.4f}  "
      f"IQR=[{_rho_w.quantile(0.25):+.4f}, {_rho_w.quantile(0.75):+.4f}]  "
      f"width={obs_iqr:.4f}")
say_c("")
if obs_iqr > sim_iqr * 1.5:
    say_c("  Observed IQR substantially wider than simulated.")
elif obs_iqr > sim_iqr * 1.1:
    say_c("  Observed IQR moderately wider than simulated.")
else:
    say_c("  Observed IQR close to or narrower than simulated.")


# =============================================================================
# SECTION C6 - Within-donor rho for each DFG
# =============================================================================
hdr_c("C6 - Per-DFG within-donor rho (median across donors)")

say_c(f"  {'gene':<10} {'med_rho':>8} {'q25':>8} {'q75':>8} {'n_donors':>9}")
say_c("  " + "-" * 48)

for dfg_gene in dfg_have:
    _per_donor_rhos = []
    for donor_id in qualifying.index:
        d_mask = sel & (donor_arr == donor_id)
        goi_d = _goi_all[d_mask]
        dfg_d = E[dfg_gene][d_mask]
        if goi_d.std() == 0 or dfg_d.std() == 0:
            continue
        r_d, _ = spearmanr(goi_d, dfg_d)
        if np.isfinite(r_d):
            _per_donor_rhos.append(r_d)
    _arr = np.array(_per_donor_rhos)
    if len(_arr) > 0:
        say_c(f"  {dfg_gene:<10} {np.median(_arr):>+8.4f} "
              f"{np.percentile(_arr, 25):>+8.4f} "
              f"{np.percentile(_arr, 75):>+8.4f} "
              f"{len(_arr):>9}")


# =============================================================================
# SECTION C7 - Verdict
# =============================================================================
hdr_c("C7 - VERDICT")

med_w = float(_rho_w.median())

if med_w >= 0.75:
    verdict = "PASS"
    say_c(f"  Within-donor median = {med_w:+.4f} (>= 0.75)")
elif med_w >= 0.40:
    verdict = "PARTIAL"
    say_c(f"  Within-donor median = {med_w:+.4f} (0.40 - 0.75)")
else:
    verdict = "FAIL"
    say_c(f"  Within-donor median = {med_w:+.4f} (< 0.40)")

say_c("")
say_c(f"  VERDICT: {verdict}")
say_c(f"    pooled rho           = {rho_pooled:+.4f}")
say_c(f"    within-donor median  = {med_w:+.4f}")
say_c(f"    between-donor rho    = {rho_between:+.4f}")
say_c(f"    depth-partialled:")
say_c(f"      pooled             = {rho_pooled_depth:+.4f}")
say_c(f"      within-donor med   = {_rho_wd.median():+.4f}")
say_c(f"      between-donor      = {rho_between_depth:+.4f}")


# =============================================================================
# SECTION C8 - Save
# =============================================================================
hdr_c("Saving")

df_1c.to_csv(os.path.join(OUT_DIR, "S1c_within_donor_rho.csv"), index=False)
say_c(f"  saved: {OUT_DIR}/S1c_within_donor_rho.csv  ({len(df_1c)} donors)")

with open(os.path.join(OUT_DIR, "S1c_within_donor_report.txt"), "w",
          encoding="utf-8") as f:
    f.write("\n".join(_report_1c))
say_c(f"  saved: {OUT_DIR}/S1c_within_donor_report.txt")

say_c("")
say_c("=" * 74)
say_c("  PART 1c COMPLETE")
say_c("=" * 74)
