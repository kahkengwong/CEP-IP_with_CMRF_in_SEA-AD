# =============================================================================
# Part8_Data-export-CP10K_MT-CO3_8DFG.py
# SEA-AD CEP-IP/CMRF, Part 8, export of one arm for the R parts
#   Selects the cells of one arm and writes the inputs R Parts 9 to 13 read:
#   per-region GOI and composite tables for the GAMs, raw and CP10K count
#   matrices, cell and donor metadata, embeddings, a manifest, an export record
#   and a guardrail log. The arm is set by leaving one ARM CONFIGURATION block
#   uncommented (Mic13_PVM3, PVM3_all, PVM1_all or HKG8_matched); the M4-top
#   arms read the Part 2 joblibs.
#
# Usage (from the SEA-AD folder, or set SEAAD_ROOT), after Part 2
#   python3 Part8_Data-export-CP10K_MT-CO3_8DFG.py
#
# Outputs (R_export_CP10K/<arm>/ and R_export_CP10K/shared/): GAM inputs,
#   count matrices, metadata, embeddings, manifest, export record and log.
# =============================================================================

import os
import gc
import time
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.io
import scanpy as sc
import joblib


# == Module 1 panel ===========================================================
MOD1_GOI = "MT-CO3"
MOD1_DFG = [
    "MT-ATP6", "MT-CO2", "MT-ND3", "MT-ND4",
    "MT-CYB",  "MT-ND1", "MT-ND2", "MT-CO1",
]

# == HKG panel: ACTB as GOI, the other eight as the composite =================
HKG_GENES = [
    "ACTB", "EEF1A1", "RPL13A", "B2M", "YWHAZ",
    "FTL",  "PSAP",   "RAB7A",  "HNRNPA2B1",
]
HKG_GOI = "ACTB"
HKG_DFG = [g for g in HKG_GENES if g != HKG_GOI]


# =============================================================================
# ARM CONFIGURATION: leave one block uncommented
# =============================================================================

# (a) Mic13-PVM_3 - primary result (M4-top PVM_3, MT-CO3>0)
# ARM_TAG        = "Mic13_PVM3"
# TARGET_STYPE   = "Micro-PVM_3-SEAAD"
# USE_M4_TOP     = True
# SELECT_GENES   = [MOD1_GOI]
# PANEL_GOI      = MOD1_GOI
# PANEL_DFG      = list(MOD1_DFG)
# PANEL_LABEL    = "cDFG"
# RUN_SUFFIX     = ""
# NEST_SUBSET_OF = "PVM3_all"

# # (b) PVM_3 all cells - mask specificity (all PVM_3, MT-CO3>0, no M4 mask)
# ARM_TAG        = "PVM3_all"
# TARGET_STYPE   = "Micro-PVM_3-SEAAD"
# USE_M4_TOP     = False
# SELECT_GENES   = [MOD1_GOI]
# PANEL_GOI      = MOD1_GOI
# PANEL_DFG      = list(MOD1_DFG)
# PANEL_LABEL    = "cDFG"
# RUN_SUFFIX     = "_PVM3all"
# NEST_SUBSET_OF = None

# # (c) PVM_1 all cells - supertype specificity (all PVM_1, MT-CO3>0, no M4 mask)
# #     PVM_1 is labelled "Micro-PVM_1" (no -SEAAD suffix).
# ARM_TAG        = "PVM1_all"
# TARGET_STYPE   = "Micro-PVM_1"
# USE_M4_TOP     = False
# SELECT_GENES   = [MOD1_GOI]
# PANEL_GOI      = MOD1_GOI
# PANEL_DFG      = list(MOD1_DFG)
# PANEL_LABEL    = "cDFG"
# RUN_SUFFIX     = "_PVM1all"
# NEST_SUBSET_OF = None

# # (d) HKG8 cell-matched - panel specificity
# #     Same cells as arm (a) with ACTB > 0; HKG panel exported.
ARM_TAG        = "HKG8_matched"
TARGET_STYPE   = "Micro-PVM_3-SEAAD"
USE_M4_TOP     = True
SELECT_GENES   = [MOD1_GOI, HKG_GOI]
PANEL_GOI      = HKG_GOI
PANEL_DFG      = list(HKG_DFG)
PANEL_LABEL    = "cHKG8"
RUN_SUFFIX     = "_HKG8"
NEST_SUBSET_OF = "Mic13_PVM3"


# == Paths ====================================================================
BASE_DIR     = os.environ.get("SEAAD_ROOT", os.getcwd())

H5AD_PATH    = os.path.join(BASE_DIR,
               "SEAAD_Immune_HIP_MEC_LEC_ITG_MTG_FI_STG_DFC_AnG_V1C_"
               "RNAseq_final-nuclei.2026-06-22.h5ad")

# Part 2 env files are searched in <arm folder>/Part2_<Arm>_M4/, the arm folder,
# then BASE_DIR (os.listdir, because the folder names contain [ and ]).
JOBLIB_SEARCH_DIRS = {
    "neocortex": [
        os.path.join(BASE_DIR, "[1]_NeoC", "Part2_NeoC_M4"),
        os.path.join(BASE_DIR, "[1]_NeoC"),
        BASE_DIR,
    ],
    "nonneocortex": [
        os.path.join(BASE_DIR, "[2]_Non-NeoC", "Part2_NonNeoC_M4"),
        os.path.join(BASE_DIR, "[2]_Non-NeoC"),
        BASE_DIR,
    ],
}
JOBLIB_BASENAME = {
    "neocortex":    "Part2_neocortex_env.joblib",
    "nonneocortex": "Part2_nonneocortex_env.joblib",
}

# Paths from the deprecated pre-release trees are refused.
DEPRECATED_PATH_MARKERS = ["Allen_Institute_2026_apps", "R_export_scran",
                           "R_export_vPreRelease"]

OUT_ROOT     = os.path.join(BASE_DIR, "R_export_CP10K")
ARM_ROOT     = os.path.join(OUT_ROOT, ARM_TAG)

# Rewrite shared/ files even when they already exist.
FORCE_SHARED_REWRITE = False


def _reject_deprecated(path, what):
    for marker in DEPRECATED_PATH_MARKERS:
        if marker in path:
            raise RuntimeError(
                f"{what} resolves to '{path}', which contains the deprecated "
                f"path marker '{marker}'.")
    return path


def resolve_joblib(dataset_tag):
    """Locate the Part 2 env for one arm, searching Part2_<arm>_M4 then root."""
    basename = JOBLIB_BASENAME[dataset_tag]
    for d in JOBLIB_SEARCH_DIRS[dataset_tag]:
        if not os.path.isdir(d):
            continue
        if basename in os.listdir(d):
            return _reject_deprecated(os.path.join(d, basename), "Part 2 env")
    raise FileNotFoundError(
        f"'{basename}' not found for arm '{dataset_tag}' in "
        f"{JOBLIB_SEARCH_DIRS[dataset_tag]}.")


_reject_deprecated(H5AD_PATH, "H5AD_PATH")
_reject_deprecated(OUT_ROOT, "OUT_ROOT")

# == Normalisation constants (Part 1a) ========================================
NORM_TARGET_SUM = 1e4
NORM_CEILING    = np.log1p(NORM_TARGET_SUM)   # 9.21034...
NORM_CEILING_TOL = 0.01                        # float tolerance on ceiling check
NORM_TAG        = "CP10K"
NORM_UNS_KEY    = "X_normalization"
NORM_UNS_VALUE  = "ln(UP10K+1)"

# == Panel actually exported by this arm ======================================
GOI = PANEL_GOI
DFG = list(PANEL_DFG)
# Same DFG set in NeoC and NonNeoC.
DFG_NEOC  = list(DFG)
DFG_NONNC = list(DFG)

# == Region definitions =======================================================
NEOC_REGIONS  = ["DFC", "MTG", "AnG", "STG", "V1C", "ITG"]
NONNC_REGIONS = ["HIP", "MEC", "LEC", "FI"]

# Other names for DFC, remapped onto DFC.
REGION_ALIASES = {"PFC": "DFC", "DLPFC": "DFC", "A9": "DFC"}

print("=" * 72)
print(f"Part 8 - CEP-IP Data Extraction for R  [{NORM_TAG}]")
print(f"  ARM         : {ARM_TAG}")
print(f"  Supertype   : {TARGET_STYPE}")
print(f"  M4-top mask : {USE_M4_TOP}")
print(f"  Select on   : {' > 0 AND '.join(SELECT_GENES)} > 0")
print(f"  Panel GOI   : {GOI}")
print(f"  Panel DFG ({len(DFG)}): {DFG}   [{PANEL_LABEL}]")
print(f"  Norm        : {NORM_TAG} (log1p(UMIs/total × {NORM_TARGET_SUM:.0f}))")
print(f"  Output root : {ARM_ROOT}")
print("=" * 72)


# =============================================================================
# SECTION 0: Directory setup
# =============================================================================
def mkd(*parts):
    p = os.path.join(*parts)
    os.makedirs(p, exist_ok=True)
    return p

for dataset in ("neocortex", "nonneocortex"):
    for subdir in ("gam_data", "counts", "metadata", "embeddings", "deg_support"):
        mkd(ARM_ROOT, dataset, subdir)
    # Dataset-wide raw counts do not vary by arm, so they live under shared/.
    mkd(OUT_ROOT, "shared", dataset, "counts")
mkd(OUT_ROOT, "shared")
print("\nSECTION 0: Output directories created.")


# =============================================================================
# SECTION 1: Load h5ad once - extract everything needed
# =============================================================================
print("\nSECTION 1: Loading h5ad ...")
t0 = time.time()

adata_full = sc.read_h5ad(H5AD_PATH)
print(f"  Loaded: {adata_full.shape}  [{time.time()-t0:.1f}s]")

obs_all   = adata_full.obs.copy()
var_names = list(adata_full.var_names)
n_total   = adata_full.n_obs


# =============================================================================
# SECTION 1b: Normalisation guardrails (CP10K)
# =============================================================================
# Checks that adata.X = log1p(UMIs / total x 10000), natural log: provenance
# string, ceiling, recomputation on a sample, row sums and no double
# normalisation.

_gr_log = []
_gr_pass = 0
_gr_fail = 0

def _gr(check_id, condition, detail_pass, detail_fail):
    global _gr_pass, _gr_fail
    status = "PASS" if condition else "FAIL"
    detail = detail_pass if condition else detail_fail
    _gr_log.append((check_id, status, detail))
    tag = "[PASS]" if condition else "[FAIL]"
    print(f"  {tag}  {check_id}  {detail}")
    if condition:
        _gr_pass += 1
    else:
        _gr_fail += 1
    return condition

def _gr_write_log():
    log_path = os.path.join(ARM_ROOT, "Part8_Guardrail_Log.txt")
    with open(log_path, "w") as f:
        f.write("=" * 72 + "\n")
        f.write(f"Part 8 - Normalisation Guardrail Log ({NORM_TAG})\n")
        f.write(f"Date: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"ARM: {ARM_TAG}\n")
        f.write(f"h5ad: {H5AD_PATH}\n")
        f.write(f"Normalisation: adata.X ({NORM_TAG} - "
                f"log1p(UMIs/total × {NORM_TARGET_SUM:.0f}), natural log)\n")
        f.write("=" * 72 + "\n\n")
        for cid, status, detail in _gr_log:
            f.write(f"  [{status}]  {cid:<16}  {detail}\n")
        f.write("\n" + "=" * 72 + "\n")
        f.write(f"  Total : {_gr_pass + _gr_fail}\n")
        f.write(f"  PASS  : {_gr_pass}\n")
        f.write(f"  FAIL  : {_gr_fail}\n")
        if _gr_fail == 0:
            f.write("\n  All guardrails passed\n")
        else:
            f.write("\n  Guardrail failures\n")
        f.write("=" * 72 + "\n")
    print(f"\n  Guardrail log written: {log_path}")
    return log_path

print(f"\n== Normalisation guardrails ({NORM_TAG}) "
      "==================================")

# == Provenance string (advisory; the recomputation decides) ==================
_uns_val = adata_full.uns.get(NORM_UNS_KEY, None)
_gr("GR-00",
    _uns_val == NORM_UNS_VALUE,
    f"uns['{NORM_UNS_KEY}'] = '{_uns_val}'  (expected '{NORM_UNS_VALUE}')",
    f"uns['{NORM_UNS_KEY}'] = '{_uns_val}'  (expected '{NORM_UNS_VALUE}')")

# == obs index is unique ======================================================
# Barcode-to-row maps assume unique barcodes.
_idx_unique = obs_all.index.is_unique
_gr("GR-00b",
    _idx_unique,
    f"obs index unique across {n_total:,} cells",
    f"obs index has {n_total - obs_all.index.nunique():,} duplicate barcodes")
if not _idx_unique:
    _gr_write_log()
    raise RuntimeError("GR-00b FAIL: obs index is not unique.")

# == layers['UMIs'] exists ====================================================
_gr("GR-01",
    "UMIs" in adata_full.layers,
    f"layers['UMIs'] present  (available: {list(adata_full.layers.keys())})",
    f"layers['UMIs'] MISSING  (available: {list(adata_full.layers.keys())})")
if "UMIs" not in adata_full.layers:
    _gr_write_log()
    raise RuntimeError("GR-01 FAIL: layers['UMIs'] not found.")

# == Load X_raw from layers['UMIs'] ===========================================
X_raw = adata_full.layers["UMIs"]
if not sp.issparse(X_raw):
    X_raw = sp.csr_matrix(X_raw)
else:
    X_raw = X_raw.tocsr()

# == Raw counts are integer-like ==============================================
_raw_max = float(X_raw.max())
_raw_sample = X_raw.data[:min(10000, len(X_raw.data))]
_raw_int = np.allclose(_raw_sample, np.round(_raw_sample), atol=0.01)
_gr("GR-02",
    _raw_int,
    f"layers['UMIs'] values are integer-like (sampled {len(_raw_sample):,} nz)",
    f"layers['UMIs'] not integer-like")

# == Raw count maximum above 20 ===============================================
_gr("GR-03",
    _raw_max > 20,
    f"layers['UMIs'] max = {_raw_max:.0f} (>> 20, consistent with raw UMIs)",
    f"layers['UMIs'] max = {_raw_max:.2f} (<= 20, too low for raw UMIs)")

# == adata.raw status (not counted) ===========================================
_raw_is_none = adata_full.raw is None
_gr_log.append(("GR-04-info", "INFO",
                f"adata.raw is None = {_raw_is_none}"))
print(f"  [INFO]  GR-04-info  adata.raw is None = {_raw_is_none}")

# == Load X_norm from adata.X =================================================
print(f"  Loading X_norm = adata.X ({NORM_TAG}, no further transform) ...")
X_norm = adata_full.X
if not sp.issparse(X_norm):
    X_norm = sp.csr_matrix(X_norm)
else:
    X_norm = X_norm.tocsr()
print(f"  X_norm loaded: shape={X_norm.shape}  nnz={X_norm.nnz:,}")

# == adata.X is float =========================================================
_gr("GR-05",
    X_norm.dtype in (np.float32, np.float64),
    f"adata.X dtype = {X_norm.dtype}  (float)",
    f"adata.X dtype = {X_norm.dtype}  (not float)")

# == adata.X maximum within log1p(1e4) ========================================
_norm_max = float(X_norm.max())
_gr("GR-06",
    _norm_max <= NORM_CEILING + NORM_CEILING_TOL,
    f"adata.X max = {_norm_max:.4f}  (<= {NORM_CEILING:.4f} ceiling, "
    f"consistent with {NORM_TAG})",
    f"adata.X max = {_norm_max:.4f}  (> {NORM_CEILING:.4f} ceiling)")

# == adata.X has no negative values ===========================================
_norm_min = float(X_norm.data.min()) if X_norm.nnz > 0 else 0.0
_gr("GR-07",
    _norm_min >= 0.0,
    f"adata.X min(nz) = {_norm_min:.6f}  (>= 0, consistent with log1p(...))",
    f"adata.X min(nz) = {_norm_min:.6f}  (negative)")

# == Raw maximum far above the normalised maximum =============================
_gr("GR-08",
    _raw_max > _norm_max * 10,
    f"max(X_raw)/max(X_norm) = {_raw_max/_norm_max:.1f}x  "
    f"({_raw_max:.0f} vs {_norm_max:.4f} - raw counts >> {NORM_TAG})",
    f"max(X_raw)/max(X_norm) = {_raw_max/_norm_max:.2f}x  "
    f"({_raw_max:.2f} vs {_norm_max:.4f})")

# == Raw and normalised matrices are distinct objects =========================
_gr("GR-09",
    X_raw is not X_norm,
    "X_raw and X_norm are distinct objects",
    "X_raw is X_norm (same object)")

# == CP10K recomputed on a sample of cells and genes ==========================
_n_check = min(200, X_raw.shape[0])
np.random.seed(42)
_check_idx = np.random.choice(X_raw.shape[0], _n_check, replace=False)
_check_genes = np.random.choice(X_raw.shape[1], min(500, X_raw.shape[1]),
                                replace=False)

_Xnorm_sample = np.asarray(X_norm[_check_idx][:, _check_genes].todense())

# 10a: X_norm matches adata.X
_adataX_sample = np.asarray(
    sp.csr_matrix(adata_full.X)[_check_idx][:, _check_genes].todense())
_match_adataX = np.allclose(_Xnorm_sample, _adataX_sample, atol=1e-6)
_max_diff_adataX = float(np.max(np.abs(_Xnorm_sample - _adataX_sample)))
_gr("GR-10a",
    _match_adataX,
    f"X_norm == adata.X spot-check MATCH  "
    f"(max_diff={_max_diff_adataX:.2e}, n={_n_check} cells × "
    f"{len(_check_genes)} genes)",
    f"X_norm != adata.X  (max_diff={_max_diff_adataX:.2e})")

# 10b: X_norm equals log1p(CP10K(UMIs))
_Xraw_check = np.asarray(X_raw[_check_idx][:, _check_genes].todense())
_full_rowsums = np.asarray(X_raw[_check_idx].sum(axis=1)).ravel()[:, None]
_full_rowsums[_full_rowsums == 0] = 1.0
_expected_cp10k = np.log1p(_Xraw_check / _full_rowsums * NORM_TARGET_SUM)
_match_cp10k = np.allclose(_Xnorm_sample, _expected_cp10k, atol=1e-4)
_max_diff_cp10k = float(np.max(np.abs(_Xnorm_sample - _expected_cp10k)))
_gr("GR-10b",
    _match_cp10k,
    f"X_norm == log1p(UMIs/total × {NORM_TARGET_SUM:.0f})  "
    f"(max_diff={_max_diff_cp10k:.2e})",
    f"X_norm != log1p(UMIs/total × {NORM_TARGET_SUM:.0f})  "
    f"(max_diff={_max_diff_cp10k:.4f})")

# 10c: no double normalisation
_adataX_full = np.asarray(
    sp.csr_matrix(adata_full.X)[_check_idx].todense())
_xrs = _adataX_full.sum(axis=1, keepdims=True)
_xrs[_xrs == 0] = 1.0
_double_norm = np.log2(_adataX_full[:, _check_genes] / _xrs * NORM_TARGET_SUM + 1)
_match_double = np.allclose(_Xnorm_sample, _double_norm, atol=1e-4)
_max_diff_double = float(np.max(np.abs(_Xnorm_sample - _double_norm)))
_gr("GR-10c",
    not _match_double,
    f"X_norm != log2(CP10K(adata.X)+1)  (max_diff={_max_diff_double:.4f}) "
    f"- not double-normalised",
    f"X_norm ~ log2(CP10K(adata.X)+1)  (max_diff={_max_diff_double:.2e}) "
    f"- double normalisation")

# 10d: expm1(X_norm) row sums pin NORM_TARGET_SUM (all genes, a few cells)
_n_expm1 = min(50, _n_check)
_expm1_full = np.asarray(
    sp.csr_matrix(adata_full.X)[_check_idx[:_n_expm1]].todense())
_expm1_full_rs = np.expm1(_expm1_full).sum(axis=1)
_expm1_cv = np.std(_expm1_full_rs) / np.mean(_expm1_full_rs)
_expm1_mean = np.mean(_expm1_full_rs)
_gr("GR-10d",
    abs(_expm1_mean - NORM_TARGET_SUM) < 1.0 and _expm1_cv < 1e-6,
    f"expm1(X_norm) row sums pin {NORM_TARGET_SUM:.0f}  "
    f"(mean={_expm1_mean:.6f}, CV={_expm1_cv:.2e}, n={_n_expm1})",
    f"expm1(X_norm) row sums do not pin {NORM_TARGET_SUM:.0f}  "
    f"(mean={_expm1_mean:.4f}, CV={_expm1_cv:.2e})")

# == Same sparsity pattern ====================================================
_gr("GR-11",
    X_raw.nnz == X_norm.nnz,
    f"nnz(X_raw) == nnz(X_norm) = {X_raw.nnz:,}  (same sparsity pattern)",
    f"nnz(X_raw)={X_raw.nnz:,} != nnz(X_norm)={X_norm.nnz:,}")

del _Xnorm_sample, _adataX_sample, _Xraw_check, _adataX_full
del _expected_cp10k, _double_norm, _full_rowsums, _xrs
del _expm1_full, _expm1_full_rs
gc.collect()

# == Guardrail checkpoint =====================================================
print(f"\n  Guardrail summary: {_gr_pass} PASS / {_gr_fail} FAIL "
      f"out of {_gr_pass + _gr_fail} checks")
if _gr_fail > 0:
    _gr_write_log()
    raise RuntimeError(
        f"Part 8 ABORTED: {_gr_fail} guardrail(s) FAILED.  "
        f"See {os.path.join(ARM_ROOT, 'Part8_Guardrail_Log.txt')}")

print(f"\n  X_raw  shape: {X_raw.shape}  nnz={X_raw.nnz:,}  "
      f"(layers['UMIs'], raw integer counts)")
print(f"  X_norm shape: {X_norm.shape}  nnz={X_norm.nnz:,}  "
      f"(adata.X, {NORM_TAG} - no further transform)")


# =============================================================================
# Column indices for GOI and DFGs
# =============================================================================
def gene_idx(gene, vnames):
    try:
        return vnames.index(gene)
    except ValueError:
        print(f"  WARNING: {gene} not found in var_names")
        return None

goi_idx = gene_idx(GOI, var_names)

# Build paired (gene_name, var_index) lists, dropping any that don't resolve
_nc_pairs = [(g, gene_idx(g, var_names)) for g in DFG_NEOC]
_nn_pairs = [(g, gene_idx(g, var_names)) for g in DFG_NONNC]
dfg_nc_genes = [g for g, i in _nc_pairs if i is not None]
dfg_nc_idx   = [i for g, i in _nc_pairs if i is not None]
dfg_nn_genes = [g for g, i in _nn_pairs if i is not None]
dfg_nn_idx   = [i for g, i in _nn_pairs if i is not None]

# Cell-selection genes - every one must be > 0 for a cell to be selected.
_sel_pairs   = [(g, gene_idx(g, var_names)) for g in SELECT_GENES]
select_genes = [g for g, i in _sel_pairs if i is not None]
select_idx   = [i for g, i in _sel_pairs if i is not None]

print(f"\n  Panel GOI '{GOI}' index: {goi_idx}")
print(f"  DFG neocortex present    : {len(dfg_nc_idx)}/{len(DFG_NEOC)}"
      f"  {dfg_nc_genes}")
print(f"  DFG non-neocortex present: {len(dfg_nn_idx)}/{len(DFG_NONNC)}"
      f"  {dfg_nn_genes}")
print(f"  Selection genes resolved : {len(select_idx)}/{len(SELECT_GENES)}"
      f"  {select_genes}")

# == Panel GOI, DFGs and selection genes must resolve =========================
_gr("GR-15-goi",
    goi_idx is not None,
    f"panel GOI '{GOI}' resolved at var_names index {goi_idx}",
    f"panel GOI '{GOI}' not found in var_names")
_gr("GR-15-neoc",
    len(dfg_nc_idx) == len(DFG_NEOC),
    f"all {len(DFG_NEOC)} NeoC DFGs resolved",
    f"only {len(dfg_nc_idx)}/{len(DFG_NEOC)} NeoC DFGs resolved; "
    f"missing {sorted(set(DFG_NEOC) - set(dfg_nc_genes))}")
_gr("GR-15-nonnc",
    len(dfg_nn_idx) == len(DFG_NONNC),
    f"all {len(DFG_NONNC)} NonNeoC DFGs resolved",
    f"only {len(dfg_nn_idx)}/{len(DFG_NONNC)} NonNeoC DFGs resolved; "
    f"missing {sorted(set(DFG_NONNC) - set(dfg_nn_genes))}")
_gr("GR-15-select",
    len(select_idx) == len(SELECT_GENES),
    f"all {len(SELECT_GENES)} selection genes resolved: {select_genes}",
    f"only {len(select_idx)}/{len(SELECT_GENES)} selection genes resolved; "
    f"missing {sorted(set(SELECT_GENES) - set(select_genes))}")

if _gr_fail:
    _gr_write_log()
    raise RuntimeError(
        f"Part 8 ABORTED at GR-15: {_gr_fail} guardrail(s) FAILED.")


# =============================================================================
# Identify obs column names
# =============================================================================
obsm_keys = list(adata_full.obsm.keys())
print(f"\n  obsm keys available: {obsm_keys}")

def find_col(obs, keywords):
    for kw in keywords:
        # Exact match first, then shortest substring
        exact = [c for c in obs.columns if c.lower() == kw.lower()]
        if exact:
            return exact[0]
        hits = [c for c in obs.columns if kw.lower() in c.lower()]
        if hits:
            return min(hits, key=len)
    return None

region_col = find_col(obs_all, ["Brain Region", "region"])   or "Brain Region"
stype_col  = find_col(obs_all, ["Supertype", "supertype"])   or "Supertype"
donor_col  = find_col(obs_all, ["Donor", "donor"])           or "Donor ID"
adnc_col   = find_col(obs_all, ["AD neuropath", "ADNC"])     or "Overall AD neuropathological Change"
braak_col  = find_col(obs_all, ["Braak", "braak"])           or "Braak"
cerad_col  = find_col(obs_all, ["CERAD", "cerad"])           or "CERAD score"
thal_col   = find_col(obs_all, ["Thal", "thal"])             or "Thal"
cog_col    = find_col(obs_all, ["Cognitive Status"])          or "Cognitive Status"
apoe_col   = find_col(obs_all, ["APOE", "apoe"])             or "APOE Genotype"
mmse_col   = find_col(obs_all, ["MMSE", "mmse"])             or "Last MMSE Score"
age_col    = find_col(obs_all, ["Age at Death"])              or "Age at Death"
caa_col    = find_col(obs_all, ["CAA", "caa"])               or "Overall CAA Score"

# Donor-level clinical and laboratory columns
pmi_col    = find_col(obs_all, ["PMI"])
rin_col    = find_col(obs_all, ["RIN"])
ph_col     = find_col(obs_all, ["Brain pH", "pH"])
weight_col = find_col(obs_all, ["Fresh Brain Weight", "Brain Weight"])
severe_col = find_col(obs_all, ["Severely Affected"])

# These five must resolve.
CRITICAL_COLS = [("Region", region_col), ("Supertype", stype_col),
                 ("Donor", donor_col), ("ADNC", adnc_col),
                 ("Braak", braak_col)]
OPTIONAL_COLS = [("CERAD", cerad_col), ("Thal", thal_col),
                 ("CogStatus", cog_col), ("APOE", apoe_col),
                 ("MMSE", mmse_col), ("Age", age_col), ("CAA", caa_col),
                 ("PMI", pmi_col), ("RIN", rin_col), ("BrainPH", ph_col),
                 ("BrainWeight", weight_col), ("SevereFlag", severe_col)]

print(f"\n  Column mapping:")
_missing_critical = []
for nm, col in CRITICAL_COLS:
    ok = col is not None and col in obs_all.columns
    print(f"    {nm:<12}: {col}  [{'OK' if ok else 'missing (critical)'}]")
    if not ok:
        _missing_critical.append(nm)
for nm, col in OPTIONAL_COLS:
    ok = col is not None and col in obs_all.columns
    print(f"    {nm:<12}: {col}  "
          f"[{'OK' if ok else 'missing, dropped from donor_clinical'}]")

if _missing_critical:
    raise RuntimeError(
        f"Critical obs columns did not resolve: {_missing_critical}.")

# Region aliases remapped onto DFC; an alias present together with DFC aborts.
obs_all[region_col] = obs_all[region_col].astype(str)
_alias_present = sorted({a for a in REGION_ALIASES
                         if (obs_all[region_col] == a).any()})
if _alias_present:
    _collisions = [a for a in _alias_present
                   if (obs_all[region_col] == REGION_ALIASES[a]).any()]
    if _collisions:
        raise RuntimeError(
            f"Region alias collision: obs contains {_collisions} and their "
            f"target code(s) {[REGION_ALIASES[a] for a in _collisions]} "
            f"simultaneously.")
    obs_all[region_col] = obs_all[region_col].replace(REGION_ALIASES)
    print(f"\n  Region aliases remapped: "
          f"{ {a: REGION_ALIASES[a] for a in _alias_present} }")

# Every declared region must be present.
_region_counts = obs_all[region_col].value_counts()
_absent = [r for r in NEOC_REGIONS + NONNC_REGIONS
           if int(_region_counts.get(r, 0)) == 0]
if _absent:
    raise RuntimeError(
        f"Declared regions absent from obs['{region_col}']: {_absent}. "
        f"Present: {sorted(_region_counts.index.tolist())}")
print(f"  All {len(NEOC_REGIONS) + len(NONNC_REGIONS)} declared regions "
      f"present in obs")

# TARGET_STYPE must exist in obs
_stype_vals = obs_all[stype_col].astype(str).unique().tolist()
if TARGET_STYPE not in _stype_vals:
    raise RuntimeError(
        f"TARGET_STYPE '{TARGET_STYPE}' not found in obs['{stype_col}']. "
        f"Available: {sorted(_stype_vals)}")
print(f"  TARGET_STYPE '{TARGET_STYPE}' found in obs ({stype_col})")


# == Ordinal maps from the Part 2 NeoC env (the same for every arm) ===========
ORDINAL_SRC = joblib.load(resolve_joblib("neocortex"))
_missing_maps = [k for k in ("BRAAK_MAP", "CERAD_MAP", "THAL_MAP", "CAA_MAP")
                 if k not in ORDINAL_SRC]
if _missing_maps:
    raise RuntimeError(
        f"Part 2 env at '{resolve_joblib('neocortex')}' does not contain "
        f"{_missing_maps}.")
BRAAK_MAP = ORDINAL_SRC["BRAAK_MAP"]
CERAD_MAP = ORDINAL_SRC["CERAD_MAP"]
THAL_MAP  = ORDINAL_SRC["THAL_MAP"]
CAA_MAP   = ORDINAL_SRC["CAA_MAP"]

# Every observed level must be a key in its map.
for _name, _col, _map in (("braak", braak_col, BRAAK_MAP),
                          ("cerad", cerad_col, CERAD_MAP),
                          ("thal",  thal_col,  THAL_MAP),
                          ("caa",   caa_col,   CAA_MAP)):
    if _col is None or _col not in obs_all.columns:
        _gr_log.append((f"GR-20-{_name}", "SKIP",
                        f"column for {_name} not present in obs"))
        continue
    _levels = sorted(str(v) for v in obs_all[_col].dropna().unique())
    _uncov = [v for v in _levels if v not in _map]
    _gr(f"GR-20-{_name}",
        len(_uncov) == 0,
        f"{_col}: all {len(_levels)} observed levels covered by the map",
        f"{_col}: {len(_uncov)} observed level(s) absent from the Part 2 map: "
        f"{_uncov}")

if _gr_fail:
    _gr_write_log()
    raise RuntimeError(
        f"Part 8 ABORTED at GR-20: {_gr_fail} ordinal map(s) do not cover the "
        f"observed levels.")

# Per-dataset selected counts, filled by extract_dataset for the export record.
N_SELECTED = {}


# =============================================================================
# SECTION 2: Shared exports - gene metadata and donor clinical table
# =============================================================================
print("\nSECTION 2: Shared exports ...")

shared_dir = os.path.join(OUT_ROOT, "shared")


def _write_shared(path, writer, label):
    """Write a shared file unless it exists and FORCE_SHARED_REWRITE is False."""
    if os.path.exists(path) and not FORCE_SHARED_REWRITE:
        print(f"  {label}: already exists, skipping "
              f"(set FORCE_SHARED_REWRITE = True to rebuild)")
        return False
    writer(path)
    return True


def _obs_to_frame(df):
    """obs -> export frame keyed on 'barcode'.

    The index carries the barcode. Any pre-existing obs column is left under
    its own name so the export never contains two columns called 'barcode'.
    """
    out_df = df.copy()
    out_df.index.name = "barcode"
    out_df = out_df.reset_index()
    return out_df


# Gene metadata
def _w_gene_meta(path):
    gm = adata_full.var.copy()
    gm.index.name = "gene_symbol"
    gm.reset_index(inplace=True)
    gm.to_csv(path, index=False)
    print(f"  gene_metadata.csv: {len(gm)} genes")


_write_shared(os.path.join(shared_dir, "gene_metadata.csv"),
              _w_gene_meta, "gene_metadata.csv")


# Donor-level clinical table
def _w_donor_clin(path):
    # region_col is cell-level, so it is left out.
    clin_cols = [donor_col, adnc_col, braak_col, cerad_col, thal_col,
                 cog_col, apoe_col, mmse_col, age_col, caa_col,
                 pmi_col, rin_col, ph_col, weight_col, severe_col]
    clin_present = [c for c in clin_cols
                    if c is not None and c in obs_all.columns]
    clin_present = list(dict.fromkeys(clin_present))   # de-duplicate

    # Columns that vary within a donor are dropped.
    nonconstant = []
    for c in clin_present:
        if c == donor_col:
            continue
        if obs_all.groupby(donor_col, observed=True)[c].nunique(
                dropna=False).max() > 1:
            nonconstant.append(c)
    if nonconstant:
        _msg = f"not constant within donor, dropped: {nonconstant}"
        print(f"    WARNING: {_msg}")
        _gr_log.append(("GR-18-donorconst", "WARN", _msg))
        clin_present = [c for c in clin_present if c not in nonconstant]
        # ADNC and Braak are required.
        _required = [c for c in (adnc_col, braak_col) if c in nonconstant]
        if _required:
            raise RuntimeError(
                f"donor_clinical: {_required} vary within donor and were "
                f"dropped.")
    else:
        _gr_log.append(("GR-18-donorconst", "INFO",
                        f"all {len(clin_present) - 1} clinical columns are "
                        f"constant within donor"))

    donor_clin = (obs_all[clin_present]
                  .drop_duplicates(subset=[donor_col])
                  .sort_values(donor_col)
                  .reset_index(drop=True))
    donor_clin.to_csv(path, index=False)
    print(f"  donor_clinical.csv: {len(donor_clin)} donors × "
          f"{len(donor_clin.columns)} columns  {clin_present}")


_write_shared(os.path.join(shared_dir, "donor_clinical.csv"),
              _w_donor_clin, "donor_clinical.csv")

print(f"\n  h5ad obs index name: '{obs_all.index.name or '(unnamed)'}' "
      f"-> exported as 'barcode'")


# Full obs metadata (all cells)
def _w_obs_all(path):
    df = _obs_to_frame(obs_all)
    df.to_csv(path, index=False)
    print(f"  obs_all_cells.csv: {len(df)} cells × {len(df.columns)} columns")


_write_shared(os.path.join(shared_dir, "obs_all_cells.csv"),
              _w_obs_all, "obs_all_cells.csv")

# All embeddings
for key in obsm_keys:
    _emb_name = f"embedding_{key.replace('X_','').replace('/','_')}.csv"

    def _w_emb(path, key=key, name=_emb_name):
        arr = adata_full.obsm[key]
        cols = [f"{key}_{i+1}" for i in range(arr.shape[1])]
        emb_df = pd.DataFrame(arr, columns=cols, index=obs_all.index)
        emb_df.index.name = "barcode"
        emb_df.reset_index().to_csv(path, index=False)
        print(f"  {name}: {arr.shape}")

    _write_shared(os.path.join(shared_dir, _emb_name), _w_emb, _emb_name)

print("  Section 2 complete.")


# =============================================================================
# Helper: sparse matrix export in 10x CellRanger format
# =============================================================================
def _export_sparse_matrix(X_sparse, barcodes, features, out_dir, var_df,
                          tag="", prefix="matrix"):
    os.makedirs(out_dir, exist_ok=True)
    X_csc = X_sparse.tocsc()

    matrix_path   = os.path.join(out_dir, f"{prefix}.mtx")
    barcodes_path = os.path.join(out_dir, f"{prefix}_barcodes.tsv")
    features_path = os.path.join(out_dir, f"{prefix}_features.tsv")

    print(f"    Writing {tag}: {X_csc.shape}  nnz={X_csc.nnz:,}")

    scipy.io.mmwrite(matrix_path, X_csc)

    with open(barcodes_path, "w") as f:
        f.write("\n".join(barcodes) + "\n")

    # Column 1 = Ensembl ID where var carries one, column 2 = gene symbol.
    gene_ids = None
    if var_df is not None and len(var_df) == len(features):
        for cand in ("gene_ids", "gene_id", "ensembl_id", "Ensembl"):
            if cand in var_df.columns:
                gene_ids = var_df[cand].astype(str).tolist()
                break
    if gene_ids is None:
        gene_ids = list(features)

    with open(features_path, "w") as f:
        f.write("\n".join(f"{gid}\t{sym}\tGene Expression"
                          for gid, sym in zip(gene_ids, features)) + "\n")

    print(f"      -> {matrix_path}")
    print(f"      -> {barcodes_path}")
    print(f"      -> {features_path}")


# =============================================================================
# SECTION 3: Cell selection - selects cells per arm configuration
# =============================================================================
def select_cells(dataset_tag, regions, dfg_idx_list, dfg_genes_list):
    """Select the cells of the current arm in one dataset (M4-top mask from the
    Part 2 env, or all cells of the supertype), keeping cells with every
    selection gene > 0."""
    print(f"\n{'='*72}")
    print(f"  SECTION 3 [{dataset_tag.upper()} / {ARM_TAG}] - selecting cells")
    print(f"{'='*72}")

    obs_all_idx_list = obs_all.index.tolist()
    obs_all_idx_map = {b: i for i, b in enumerate(obs_all_idx_list)}

    if USE_M4_TOP:
        # == M4-top mode: load joblib for the mask ============================
        jl_path = resolve_joblib(dataset_tag)
        print(f"  Loading joblib: {jl_path}")
        env = joblib.load(jl_path)
        obs_jl       = env["obs_full"]
        M4_TOP_MASK  = env["M4_TOP_MASK"]
        score_cols_d = env["score_cols"]
        umap_jl      = env.get("umap_coords")
        TOP_PCT      = env.get("TOP_PCT", 0.030)

        # == The Part 2 env must come from this h5ad and normalisation ========
        _jl_h5ad = env.get("FULL_H5AD")
        _h5ad_ok = (_jl_h5ad is None or
                    os.path.basename(str(_jl_h5ad)) ==
                    os.path.basename(H5AD_PATH))
        _gr(f"GR-19-h5ad-{dataset_tag[:4]}",
            _h5ad_ok,
            f"Part 2 env scored '{os.path.basename(str(_jl_h5ad))}', "
            f"matching this run" if _jl_h5ad else
            "Part 2 env records no FULL_H5AD",
            f"Part 2 env scored '{os.path.basename(str(_jl_h5ad))}' but this "
            f"run reads '{os.path.basename(H5AD_PATH)}'")

        _jl_target = env.get("NORM_TARGET_SUM")
        _target_ok = (_jl_target is None or
                      abs(float(_jl_target) - NORM_TARGET_SUM) < 1e-6)
        _gr(f"GR-19-norm-{dataset_tag[:4]}",
            _target_ok,
            f"Part 2 NORM_TARGET_SUM = {_jl_target} matches this run"
            if _jl_target else
            "Part 2 env records no NORM_TARGET_SUM",
            f"Part 2 NORM_TARGET_SUM = {_jl_target}, this run uses "
            f"{NORM_TARGET_SUM:.0f}")

        if not (_h5ad_ok and _target_ok):
            _gr_write_log()
            raise RuntimeError(
                f"Part 8 ABORTED at GR-19: the Part 2 env at '{jl_path}' was "
                f"not produced from this h5ad or this normalisation.")

        print(f"  Inherited contract: RELEASE={env.get('RELEASE')}  "
              f"TOP_PCT={TOP_PCT}  N_M4_TOP={env.get('N_M4_TOP')}  "
              f"KEEP_FLAGGED_NUCLEI={env.get('KEEP_FLAGGED_NUCLEI')}")
        print(f"  Scorers in env: {sorted(score_cols_d.keys())}")

        n_ds = len(obs_jl)
        jl_barcodes = obs_jl.index.tolist()

        # The env must cover only this dataset's regions.
        _jl_regions = set(obs_jl[region_col].astype(str).unique())
        _jl_regions = {REGION_ALIASES.get(r, r) for r in _jl_regions}
        if not _jl_regions.issubset(set(regions)):
            raise RuntimeError(
                f"Arm mismatch: '{jl_path}' covers regions "
                f"{sorted(_jl_regions)}, which is not a subset of the "
                f"'{dataset_tag}' regions {regions}.")
        print(f"  Arm check: joblib regions {sorted(_jl_regions)} "
              f"within {dataset_tag}")

        # Verify barcodes match h5ad
        obs_all_idx_set = set(obs_all_idx_list)
        n_matched = sum(1 for b in jl_barcodes if b in obs_all_idx_set)
        print(f"  Barcode match: {n_matched:,}/{len(jl_barcodes):,} "
              f"found in h5ad")

        if n_matched == 0:
            # Positional fallback
            print("  WARNING: barcodes not matched - "
                  "attempting positional fallback ...")
            _region_mask = obs_all[region_col].astype(str).isin(regions)
            _region_barcodes = obs_all.index[_region_mask].tolist()
            if len(_region_barcodes) == n_ds:
                jl_barcodes = _region_barcodes
                print(f"  Positional fallback succeeded: "
                      f"{len(jl_barcodes):,} barcodes")
            else:
                raise ValueError(
                    f"Barcode mismatch: joblib has {n_ds} cells but region "
                    f"filter gives {len(_region_barcodes)} for {regions}.")

        # Global h5ad row indices for ALL dataset cells
        global_row_idx = np.array(
            [obs_all_idx_map[b] for b in jl_barcodes], dtype=np.int64)

        # Build combined mask: M4-top AND target supertype AND GOI>0
        mask = M4_TOP_MASK.copy()
        print(f"  M4-top cells: {int(mask.sum()):,} / {n_ds:,}")

        # Filter to target supertype
        stype_mask = (obs_jl[stype_col].astype(str) == TARGET_STYPE).values
        mask &= stype_mask
        print(f"  After supertype filter ({TARGET_STYPE}): "
              f"{int(mask.sum()):,}")

        # Selection genes - every one must be > 0
        for _sg_name, _sg_idx in zip(select_genes, select_idx):
            _sg_vals = np.asarray(
                X_norm[global_row_idx, _sg_idx].todense()).ravel()
            mask &= (_sg_vals > 0)
            print(f"  After {_sg_name}>0 filter: {int(mask.sum()):,}")

        # Apply mask
        sel_local = np.where(mask)[0]
        sel_barcodes = [jl_barcodes[i] for i in sel_local]
        sel_global = global_row_idx[sel_local]
        sel_obs = obs_jl.iloc[sel_local].copy()

        # Score columns: subset to selected cells
        sel_scores = {}
        for mkey, sc_arr in score_cols_d.items():
            sel_scores[mkey] = sc_arr[sel_local]

        # UMAP for selected
        sel_umap = umap_jl[sel_local] if umap_jl is not None else None

        # Dataset-wide data for raw_dataset export
        ds_barcodes = jl_barcodes
        ds_global = global_row_idx

        del env
        gc.collect()

        return {
            "barcodes":           sel_barcodes,
            "global_idx":         sel_global,
            "obs_selected":       sel_obs,
            "n_dataset":          n_ds,
            "region_list":        regions,
            "score_cols_d":       sel_scores,
            "umap_coords":        sel_umap,
            "dataset_barcodes":   ds_barcodes,
            "dataset_global_idx": ds_global,
            "obs_dataset":        obs_jl,
        }

    else:
        # == All-cells mode: filter obs_all directly ===========================
        print("  All-cells mode - filtering obs_all directly (no joblib)")

        # Dataset mask: cells in the target regions
        ds_mask = obs_all[region_col].astype(str).isin(regions)
        n_ds = int(ds_mask.sum())
        print(f"  Dataset cells (all regions): {n_ds:,}")

        ds_barcodes = obs_all.index[ds_mask].tolist()
        ds_global = np.array(
            [obs_all_idx_map[b] for b in ds_barcodes], dtype=np.int64)

        # Selection mask: target supertype within the dataset regions
        stype_mask = ((obs_all[stype_col].astype(str) == TARGET_STYPE).values
                      & ds_mask.values)
        print(f"  {TARGET_STYPE} cells in dataset: {int(stype_mask.sum()):,}")

        # Selection genes - every one must be > 0
        sel_global = np.where(stype_mask)[0]
        for _sg_name, _sg_idx in zip(select_genes, select_idx):
            _sg_vals = np.asarray(
                X_norm[sel_global, _sg_idx].todense()).ravel()
            sel_global = sel_global[_sg_vals > 0]
            print(f"  After {_sg_name}>0 filter: {len(sel_global):,}")

        sel_barcodes = [obs_all_idx_list[i] for i in sel_global]
        sel_obs = obs_all.iloc[sel_global].copy()

        return {
            "barcodes":           sel_barcodes,
            "global_idx":         sel_global,
            "obs_selected":       sel_obs,
            "n_dataset":          n_ds,
            "region_list":        regions,
            "score_cols_d":       None,
            "umap_coords":        None,
            "dataset_barcodes":   ds_barcodes,
            "dataset_global_idx": ds_global,
            "obs_dataset":        obs_all[ds_mask].copy(),
        }


# =============================================================================
# SECTION 4: Dataset extraction function
# =============================================================================
def extract_dataset(dataset_tag, sel, dfg_idx_list, dfg_genes_list):
    """
    Extract all R-required data for one dataset (neocortex or nonneocortex)
    using the pre-selected cells from select_cells().

    sel: dict returned by select_cells()
    """
    print(f"\n{'='*72}")
    print(f"  SECTION 4 [{dataset_tag.upper()} / {ARM_TAG}] - extracting")
    print(f"{'='*72}")
    t_ds = time.time()

    sel_barcodes  = sel["barcodes"]
    sel_global    = sel["global_idx"]
    sel_obs       = sel["obs_selected"]
    n_ds          = sel["n_dataset"]
    score_cols_d  = sel["score_cols_d"]
    ds_barcodes   = sel["dataset_barcodes"]
    ds_global     = sel["dataset_global_idx"]
    obs_ds        = sel["obs_dataset"]

    n_sel = len(sel_barcodes)
    N_SELECTED[dataset_tag] = n_sel
    regions = sel["region_list"]

    out = lambda *p: os.path.join(ARM_ROOT, dataset_tag, *p)

    print(f"  Dataset cells: {n_ds:,}  |  Selected: {n_sel:,}  "
          f"({100*n_sel/n_ds:.2f}%)")

    # == METADATA EXPORT ======================================================
    print("  Exporting metadata ...")

    # Selected cells obs
    sel_obs_out = _obs_to_frame(sel_obs)
    sel_obs_out.to_csv(out("metadata", "obs_selected.csv"), index=False)
    print(f"    obs_selected.csv: {len(sel_obs_out)} cells × "
          f"{len(sel_obs_out.columns)} cols")

    # Full dataset obs
    obs_ds_out = _obs_to_frame(obs_ds)
    obs_ds_out.to_csv(out("metadata", "obs_dataset.csv"), index=False)
    print(f"    obs_dataset.csv: {len(obs_ds_out)} cells")

    # Cell index audit table
    idx_df = pd.DataFrame({
        "barcode":    sel_barcodes,
        "global_idx": sel_global,
        "region":     sel_obs[region_col].values,
        "supertype":  sel_obs[stype_col].values,
        "donor":      sel_obs[donor_col].values,
    })
    idx_df.to_csv(out("metadata", "cell_index_selected.csv"), index=False)
    print(f"    cell_index_selected.csv: {len(idx_df)} cells")

    # Score table (M4-top arms only)
    if score_cols_d is not None:
        score_df = pd.DataFrame({"barcode": sel_barcodes})
        for mkey, sc_arr in score_cols_d.items():
            score_df[f"score_{mkey}"] = sc_arr
        score_df.to_csv(
            out("metadata", "scores_selected_all_methods.csv"), index=False)
        print(f"    scores_selected_all_methods.csv: "
              f"{len(score_df)} cells × {len(score_df.columns)} cols")

    # == EMBEDDINGS EXPORT ====================================================
    print("  Exporting embeddings ...")

    umap_coords = sel["umap_coords"]
    if umap_coords is not None:
        umap_df = pd.DataFrame(
            umap_coords, columns=["UMAP_1", "UMAP_2"],
            index=sel_barcodes)
        umap_df.index.name = "barcode"
        umap_df.reset_index().to_csv(
            out("embeddings", "umap_selected.csv"), index=False)
        print(f"    umap_selected.csv exported")

    # h5ad obsm embeddings for selected cells
    for key in obsm_keys:
        arr = adata_full.obsm[key][sel_global]
        n_d = arr.shape[1]
        cols = [f"{key}_{i+1}" for i in range(n_d)]
        emb_df = pd.DataFrame(arr, columns=cols, index=sel_barcodes)
        emb_df.index.name = "barcode"
        fname = f"embedding_{key.replace('X_','').replace('/','_')}_selected.csv"
        emb_df.reset_index().to_csv(out("embeddings", fname), index=False)
    print(f"    Exported {len(obsm_keys)} embedding keys for selected cells")

    # == GAM DATA EXPORT, per region ==========================================
    print("  Exporting GAM input data ...")

    # Extract GOI and cDFG expression for selected cells
    print(f"    Extracting GOI ({GOI}) expression ...")
    goi_expr = np.asarray(
        X_norm[sel_global, goi_idx].todense()).ravel().astype(np.float32)

    print(f"    Extracting cDFG expression ({len(dfg_idx_list)} genes) ...")
    dfg_mat = np.asarray(
        X_norm[sel_global, :][:, dfg_idx_list].todense()).astype(np.float32)
    cdfg_expr = dfg_mat.mean(axis=1).ravel()

    # Per-gene DFG expressions (for LOO analysis in R Part 9)
    dfg_indiv = {}
    for gname, gidx in zip(dfg_genes_list, dfg_idx_list):
        dfg_indiv[gname] = np.asarray(
            X_norm[sel_global, gidx].todense()).ravel().astype(np.float32)

    # Build master GAM dataframe
    gam_master = pd.DataFrame({
        "barcode":    sel_barcodes,
        "region":     sel_obs[region_col].values,
        "supertype":  sel_obs[stype_col].values,
        "donor":      sel_obs[donor_col].values,
        "GOI":        goi_expr,
        "Expression": cdfg_expr,
    })
    # Individual DFG columns
    for gname, arr in dfg_indiv.items():
        gam_master[gname] = arr

    # Score columns (M4-top arms only)
    if score_cols_d is not None:
        for mkey, sc_arr in score_cols_d.items():
            gam_master[f"score_{mkey}"] = sc_arr

    gam_master.to_csv(
        out("gam_data", "gam_master_all_cells.csv"), index=False)
    print(f"    gam_master_all_cells.csv: {len(gam_master)} cells × "
          f"{len(gam_master.columns)} cols")

    # Per-region CSVs
    summary_rows = []
    skipped_regions = []
    for region in regions:
        # Selected cells in this region × target supertype
        mask_sel = ((sel_obs[region_col].astype(str) == region).values &
                    (sel_obs[stype_col].astype(str) == TARGET_STYPE).values)
        n_sel_r = int(mask_sel.sum())

        if n_sel_r < 5:
            skipped_regions.append((region, n_sel_r))
            print(f"      {region}: n_selected={n_sel_r} (<5) - no CSV written")
            continue

        sub_df = gam_master[mask_sel].copy()
        fname = (f"{region}_{TARGET_STYPE.replace('-','_').replace('/','_')}"
                 f"_selected.csv")
        sub_df.to_csv(out("gam_data", fname), index=False)

        n_unique_goi = sub_df["GOI"].nunique()
        summary_rows.append({
            "region":       region,
            "supertype":    TARGET_STYPE,
            "n_selected":   n_sel_r,
            "n_unique_GOI": n_unique_goi,
            "GOI_mean":     round(sub_df["GOI"].mean(), 4),
            "cDFG_mean":    round(sub_df["Expression"].mean(), 4),
            "file":         fname,
        })
        print(f"      {region} × {TARGET_STYPE}: n_selected={n_sel_r:,}")

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(out("gam_data", "summary_gam_inputs.csv"), index=False)
    print(f"    summary_gam_inputs.csv: {len(summary_rows)} "
          f"region×supertype combinations")

    # Regions below the 5-cell minimum are recorded, not counted.
    if skipped_regions:
        _msg = (f"{len(skipped_regions)} {dataset_tag} region(s) below the "
                f"5-cell minimum, no GAM CSV written: {skipped_regions}")
        _gr_log.append((f"GR-16-{dataset_tag[:4]}", "WARN", _msg))
        print(f"    [WARN]  GR-16-{dataset_tag[:4]}  {_msg}")
    else:
        _gr_log.append((f"GR-16-{dataset_tag[:4]}", "INFO",
                        f"all {len(regions)} {dataset_tag} regions produced "
                        f"a GAM CSV"))

    # == RAW COUNT MATRIX EXPORT ==============================================
    print("  Exporting raw count matrices ...")

    # Selected cells
    X_raw_sel = X_raw[sel_global, :]
    _export_sparse_matrix(
        X_raw_sel,
        barcodes=sel_barcodes,
        features=var_names,
        out_dir=out("counts", "raw_selected"),
        var_df=adata_full.var,
        tag=f"{dataset_tag} selected raw counts"
    )

    # All cells in the dataset (for Monocle3), written once under shared/.
    _ds_counts_dir = os.path.join(OUT_ROOT, "shared", dataset_tag,
                                  "counts", "raw_dataset")
    if os.path.exists(os.path.join(_ds_counts_dir, "matrix.mtx")) \
            and not FORCE_SHARED_REWRITE:
        print(f"    shared/{dataset_tag}/counts/raw_dataset: already exists, "
              f"skipping")
    else:
        _export_sparse_matrix(
            X_raw[ds_global, :],
            barcodes=ds_barcodes,
            features=var_names,
            out_dir=_ds_counts_dir,
            var_df=adata_full.var,
            tag=f"{dataset_tag} all-dataset raw counts (shared)"
        )
        gc.collect()

    # == NORMALISED MATRIX EXPORT (DEG support) ===============================
    print("  Exporting normalised matrix for DEG support ...")

    X_norm_sel = X_norm[sel_global, :]
    _export_sparse_matrix(
        X_norm_sel,
        barcodes=sel_barcodes,
        features=var_names,
        out_dir=out("deg_support"),
        var_df=adata_full.var,
        tag=f"{dataset_tag} selected {NORM_TAG}",
        prefix="norm"
    )

    # Dense CSV for GOI+DFG genes only (fast load in R)
    goi_dfg_genes = [GOI] + dfg_genes_list
    goi_dfg_idx   = [goi_idx] + dfg_idx_list
    goi_dfg_mat = np.asarray(
        X_norm[sel_global, :][:, goi_dfg_idx].todense()).astype(np.float32)
    goi_dfg_df = pd.DataFrame(
        goi_dfg_mat, columns=goi_dfg_genes, index=sel_barcodes)
    goi_dfg_df.index.name = "barcode"
    goi_dfg_df.reset_index().to_csv(
        out("deg_support", "norm_GOI_DFG_selected.csv"), index=False)
    print(f"    norm_GOI_DFG_selected.csv: {goi_dfg_df.shape}")

    # == MANIFEST =============================================================
    _write_manifest(dataset_tag, out, regions, dfg_genes_list,
                    summary_rows, n_ds, n_sel)

    print(f"\n  [{dataset_tag.upper()}] complete in {time.time()-t_ds:.1f}s")
    gc.collect()


# =============================================================================
# Helper: manifest file
# =============================================================================
def _write_manifest(tag, out_fn, regions, dfg_genes, summary_rows,
                    n_ds, n_sel):
    sel_label = (("M4-top" if USE_M4_TOP else "all") + f" {TARGET_STYPE}, "
                 + " AND ".join(f"{g}>0" for g in select_genes))

    lines = [
        f"CEP-IP R Export Manifest - {tag.upper()} / {ARM_TAG}",
        f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}",
        f"",
        f"ARM         : {ARM_TAG}",
        f"DATASET     : {tag}",
        f"Selection   : {sel_label}",
        f"Regions     : {regions}",
        f"Total cells : {n_ds:,}",
        f"Selected    : {n_sel:,}  ({100*n_sel/n_ds:.2f}%)",
        f"Select on   : {' AND '.join(f'{g}>0' for g in select_genes)}",
        f"Panel GOI   : {GOI}",
        f"Panel DFGs  : {dfg_genes}  [{PANEL_LABEL}]",
        f"Norm        : {NORM_TAG} (log1p(UMIs/total x {NORM_TARGET_SUM:.0f}))",
        f"",
        f"FILES:",
        f"",
        f"  metadata/",
        f"    obs_selected.csv              - obs for selected cells",
        f"    obs_dataset.csv               - obs for ALL cells in dataset",
        f"    cell_index_selected.csv       - barcode <-> global_h5ad_idx",
        f"    scores_selected_all_methods.csv - M1-M6 scores (M4-top arms only)",
        f"",
        f"  embeddings/",
        f"    umap_selected.csv             - UMAP for selected cells (if in joblib)",
        f"    embedding_*_selected.csv      - all h5ad obsm keys, selected cells",
        f"",
        f"  gam_data/",
        f"    gam_master_all_cells.csv      - master table (all selected cells):",
        f"                                    barcode, region, supertype, donor,",
        f"                                    GOI ({NORM_TAG}, {GOI})",
        f"                                    Expression ({NORM_TAG}, cDFG mean),",
        f"                                    individual DFG columns, score columns",
        f"    {{region}}_*_selected.csv       - GAM input per region",
        f"    summary_gam_inputs.csv         - n per region",
        f"",
        f"  counts/",
        f"    raw_selected/                 - 10x-format raw count matrix (selected)",
        f"      matrix.mtx, matrix_barcodes.tsv, matrix_features.tsv",
        f"",
        f"  ../../shared/{tag}/counts/raw_dataset/",
        f"                                  - 10x-format (ALL dataset cells);",
        f"                                    arm-invariant, written once",
        f"",
        f"  deg_support/",
        f"    norm.mtx                      - {NORM_TAG} (selected x all genes)",
        f"    norm_barcodes.tsv, norm_features.tsv",
        f"    norm_GOI_DFG_selected.csv     - GOI+DFG dense CSV",
        f"",
        f"GAM DATA SUMMARY:",
    ]
    for row in summary_rows:
        lines.append(
            f"  {row['region']:<6} x {row['supertype']:<40} "
            f"n_sel={row['n_selected']:>6,}")
    lines.append("")

    with open(out_fn("manifest.txt"), "w") as f:
        f.write("\n".join(lines))
    print(f"    manifest.txt written")


# =============================================================================
# SECTION 5: Run extractions for both datasets
# =============================================================================
print("\n" + "=" * 72)
print("SECTION 5: Running dataset extractions")
print("=" * 72)

# Neocortex
sel_neoc = select_cells("neocortex", NEOC_REGIONS, dfg_nc_idx, dfg_nc_genes)
extract_dataset("neocortex", sel_neoc, dfg_nc_idx, dfg_nc_genes)

# Non-neocortex
sel_nonnc = select_cells("nonneocortex", NONNC_REGIONS, dfg_nn_idx, dfg_nn_genes)
extract_dataset("nonneocortex", sel_nonnc, dfg_nn_idx, dfg_nn_genes)


# =============================================================================
# SECTION 6: Export record
# =============================================================================
# Part8_export_record.csv records what this arm exported; R_config.R and
# R_config_HKG8.R check their own definitions against it.
print("\nSECTION 6: Writing export record ...")


def _fmt_map(d):
    # int() so that values are written as "4", not "4.0".
    return ";".join(f"{k}={int(v)}" for k, v in d.items())


record_rows = [
    ("arm_tag",             ARM_TAG),
    ("target_stype",        TARGET_STYPE),
    ("use_m4_top",          str(USE_M4_TOP)),
    ("select_genes",        ";".join(select_genes)),
    ("goi",                 GOI),
    ("dfg",                 ";".join(DFG)),
    ("panel_label",         PANEL_LABEL),
    ("run_suffix",          RUN_SUFFIX),
    ("norm_tag",            NORM_TAG),
    ("norm_target_sum",     f"{NORM_TARGET_SUM:.0f}"),
    ("h5ad",                os.path.basename(H5AD_PATH)),
    ("release",             str(ORDINAL_SRC.get("RELEASE"))),
    ("neoc_regions",        ";".join(NEOC_REGIONS)),
    ("nonnc_regions",       ";".join(NONNC_REGIONS)),
    ("n_selected_neoc",     str(N_SELECTED.get("neocortex", 0))),
    ("n_selected_nonneoc",  str(N_SELECTED.get("nonneocortex", 0))),
    ("braak_map",           _fmt_map(BRAAK_MAP)),
    ("cerad_map",           _fmt_map(CERAD_MAP)),
    ("thal_map",            _fmt_map(THAL_MAP)),
    ("caa_map",             _fmt_map(CAA_MAP)),
    ("exported_utc",        time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())),
]

record_path = os.path.join(ARM_ROOT, "Part8_export_record.csv")
pd.DataFrame(record_rows, columns=["field", "value"]).to_csv(
    record_path, index=False)
print(f"  Saved: {record_path}")


# SECTION 7: Final summary
# =============================================================================
print("\n" + "=" * 72)
print("SECTION 7: Extraction complete - final inventory")
print("=" * 72)

total_files = 0
total_size  = 0
for root, dirs, files in os.walk(ARM_ROOT):
    for fname in files:
        fpath = os.path.join(root, fname)
        sz    = os.path.getsize(fpath)
        total_files += 1
        total_size  += sz

print(f"\n  Total files  : {total_files:,}  (arm tree only)")
print(f"  Total size   : {total_size/1e9:.2f} GB  (arm tree only)")
print(f"  ARM          : {ARM_TAG}")
print(f"  Supertype    : {TARGET_STYPE}")
print(f"  M4-top       : {USE_M4_TOP}")
print(f"  Select on    : {' AND '.join(f'{g}>0' for g in select_genes)}")
print(f"  Panel        : {GOI} + {len(DFG)} DFG  [{PANEL_LABEL}]")
print(f"  Output root  : {ARM_ROOT}")

print(f"\n  Key files for R:")
print(f"    {record_path}  (export record)")
print(f"    {OUT_ROOT}/shared/  (gene metadata, donor clinical, obs, embeddings)")
for ds in ("neocortex", "nonneocortex"):
    print(f"    {ARM_ROOT}/{ds}/gam_data/  (per-region CSVs)")
    print(f"    {ARM_ROOT}/{ds}/counts/raw_selected/  (matrix.mtx)")
    print(f"    {ARM_ROOT}/{ds}/deg_support/  (norm.mtx)")
    print(f"    {OUT_ROOT}/shared/{ds}/counts/raw_dataset/  (matrix.mtx)")



# =============================================================================
# SECTION 8: Post-export checks of the exported files
# =============================================================================
print("\n" + "=" * 72)
print(f"SECTION 8: Post-export normalisation cross-validation ({NORM_TAG})")
print("=" * 72)

for _ds_tag, _ds_dir in [("neocortex", "neocortex"),
                          ("nonneocortex", "nonneocortex")]:
    _norm_mtx_path = os.path.join(ARM_ROOT, _ds_dir, "deg_support", "norm.mtx")
    _norm_bar_path = os.path.join(
        ARM_ROOT, _ds_dir, "deg_support", "norm_barcodes.tsv")
    _norm_feat_path = os.path.join(
        ARM_ROOT, _ds_dir, "deg_support", "norm_features.tsv")
    _raw_mtx_path = os.path.join(
        ARM_ROOT, _ds_dir, "counts", "raw_selected", "matrix.mtx")

    # norm.mtx exists and its maximum is within the CP10K ceiling
    if os.path.exists(_norm_mtx_path):
        _Xn = scipy.io.mmread(_norm_mtx_path).tocsr()
        _xn_max = float(_Xn.max())
        _gr(f"GR-12-{_ds_tag[:4]}",
            _xn_max <= NORM_CEILING + NORM_CEILING_TOL,
            f"norm.mtx max = {_xn_max:.4f}  "
            f"(<= {NORM_CEILING:.4f}, consistent with {NORM_TAG})",
            f"norm.mtx max = {_xn_max:.2f}  "
            f"(> {NORM_CEILING:.4f})")

        # norm.mtx on disk matches adata.X for the same cells and genes
        with open(_norm_bar_path) as _bf:
            _nbarcodes = [l.strip() for l in _bf.readlines()]
        _nbar_pos = {b: i for i, b in enumerate(_nbarcodes)}

        _disk_features = []
        if os.path.exists(_norm_feat_path):
            with open(_norm_feat_path) as _ff:
                for _line in _ff:
                    _parts = _line.rstrip("\n").split("\t")
                    _disk_features.append(
                        _parts[1] if len(_parts) > 1 else _parts[0])
        _disk_gene_map = {g: i for i, g in enumerate(_disk_features)}
        _h5_gene_map = {g: i for i, g in enumerate(var_names)}

        np.random.seed(99)
        _dfg_set = DFG_NEOC if _ds_tag.startswith("neo") else DFG_NONNC
        _gr13_genes = [GOI] + _dfg_set
        _random_pool = [g for g in var_names
                        if g not in _gr13_genes and g in _disk_gene_map]
        if len(_random_pool) >= 5:
            _gr13_genes += list(
                np.random.choice(_random_pool, 5, replace=False))
        _gr13_genes = [g for g in _gr13_genes
                       if g in _disk_gene_map and g in _h5_gene_map]

        _obs_idx_map = {b: i for i, b in enumerate(obs_all.index.tolist())}
        _sample_bc = [b for b in _nbarcodes[:50]
                      if b in _obs_idx_map and b in _nbar_pos]
        _global_idx = np.array([_obs_idx_map[b] for b in _sample_bc])
        _local_idx = np.array([_nbar_pos[b] for b in _sample_bc])

        if len(_gr13_genes) < 5 or len(_global_idx) < 10:
            _gr(f"GR-13-{_ds_tag[:4]}", False, "",
                f"Insufficient genes ({len(_gr13_genes)}) or barcodes "
                f"({len(_global_idx)}) for cross-validation")
        else:
            _all_diffs = []
            for _g in _gr13_genes:
                _di = _disk_gene_map[_g]
                _hi = _h5_gene_map[_g]
                _disk_vals = np.asarray(
                    _Xn[_local_idx, _di].todense()).ravel()
                _adataX_vals = np.asarray(
                    X_norm[_global_idx, _hi].todense()).ravel()
                _all_diffs.append(
                    float(np.max(np.abs(_disk_vals - _adataX_vals))))

            _max_diff = max(_all_diffs)
            _gr(f"GR-13-{_ds_tag[:4]}",
                _max_diff < 1e-4,
                f"norm.mtx on disk vs adata.X ({NORM_TAG}): "
                f"max_diff = {_max_diff:.2e}  "
                f"({len(_gr13_genes)} genes x {len(_global_idx)} cells)",
                f"norm.mtx on disk vs adata.X ({NORM_TAG}): "
                f"max_diff = {_max_diff:.2e}")
        del _Xn
    else:
        _gr(f"GR-12-{_ds_tag[:4]}", False,
            "", f"norm.mtx not found at {_norm_mtx_path}")

    # Raw count matrix is integer-valued
    if os.path.exists(_raw_mtx_path):
        _Xr = scipy.io.mmread(_raw_mtx_path).tocsr()
        _xr_max = float(_Xr.max())
        _xr_sample = _Xr.data[:min(10000, len(_Xr.data))]
        _xr_int = np.allclose(_xr_sample, np.round(_xr_sample), atol=0.01)
        _gr(f"GR-14-{_ds_tag[:4]}",
            _xr_int,
            f"raw_selected/matrix.mtx: integer-valued, max={_xr_max:.0f}",
            f"raw_selected/matrix.mtx: not integer-valued, max={_xr_max:.2f}")
        del _Xr
    else:
        _gr(f"GR-14-{_ds_tag[:4]}", False,
            "", f"raw_selected/matrix.mtx not found at {_raw_mtx_path}")

gc.collect()

# == Arms nest: HKG8_matched in Mic13_PVM3 in PVM3_all (if the parent exists) ====
if NEST_SUBSET_OF:
    _parent_root = os.path.join(OUT_ROOT, NEST_SUBSET_OF)
    for _ds_tag in ("neocortex", "nonneocortex"):
        _child = os.path.join(ARM_ROOT, _ds_tag, "metadata",
                              "cell_index_selected.csv")
        _parent = os.path.join(_parent_root, _ds_tag, "metadata",
                               "cell_index_selected.csv")
        if not os.path.exists(_parent):
            _gr_log.append((f"GR-17-{_ds_tag[:4]}", "SKIP",
                            f"parent arm '{NEST_SUBSET_OF}' not exported yet"))
            print(f"  [SKIP]  GR-17-{_ds_tag[:4]}  parent arm "
                  f"'{NEST_SUBSET_OF}' not exported yet")
            continue
        _cb = set(pd.read_csv(_child)["barcode"])
        _pb = set(pd.read_csv(_parent)["barcode"])
        _extra = len(_cb - _pb)
        _gr(f"GR-17-{_ds_tag[:4]}",
            _extra == 0,
            f"{ARM_TAG} ({len(_cb):,}) is a subset of {NEST_SUBSET_OF} "
            f"({len(_pb):,})",
            f"{_extra:,} of {len(_cb):,} {ARM_TAG} cells are absent from "
            f"{NEST_SUBSET_OF}")

# == Write final guardrail log ================================================
_gr_log_path = _gr_write_log()

print(f"\n  Guardrail final tally: {_gr_pass} PASS / {_gr_fail} FAIL")
if _gr_fail > 0:
    print(f"\n  {_gr_fail} guardrail(s) failed")
    print(f"  Review: {_gr_log_path}")
    raise RuntimeError(
        f"Part 8 [{ARM_TAG}]: {_gr_fail} post-export guardrail(s) FAILED. "
        f"See {_gr_log_path}")

print(f"\n  All {_gr_pass} guardrails passed")
print("\n" + "=" * 72)
print(f"Part 8 + guardrails complete  [{ARM_TAG}].")
print("=" * 72)
