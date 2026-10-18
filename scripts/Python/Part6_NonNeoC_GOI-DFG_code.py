# =============================================================================
# Part6_NonNeoC_GOI-DFG_code.py
# SEA-AD CEP-IP/CMRF, Part 6, GOI-DFG screen by region and supertype,
# non-neocortex
#   Within the M4-top cells from Part 2, screens every gene as a gene of
#   interest (GOI) against all expressed partner genes, in six populations:
#   PVM_3 in each non-neocortical region (HIP, MEC, LEC, FI), and PVM_2_3 and
#   PVM_4 pooled over the regions. A partner is a dual-filter gene (DFG) when
#   Spearman |rho| >= 0.6 and Kendall |tau| >= 0.5 with the same sign.
#
# Usage (from the SEA-AD folder, or set SEAAD_ROOT), after Part 2
#   python3 Part6_NonNeoC_GOI-DFG_code.py
#   An interrupted run resumes from the Part<N>_checkpoint.pkl files.
#
# Outputs ([2]_Non-NeoC/Part6_NonNeoC_GOI-DFG/): per population, the Spearman,
#   Kendall and DualFilter workbooks per GOI, the MasterDFG ranking and logs;
#   the population tally (CSV) and the cross-population summary workbook.
# =============================================================================

# == Imports ==================================================================
import os
import sys
import gc
import time
import pickle
import logging
import warnings
warnings.filterwarnings("ignore")

# One BLAS thread per worker process
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"

import numpy as np
import pandas as pd
import scipy.sparse as sp_sparse
from scipy.stats import spearmanr, kendalltau
import h5py
import openpyxl
from openpyxl import Workbook
from openpyxl.styles import PatternFill, Font, Alignment, Border, Side
from openpyxl.utils import get_column_letter

import concurrent.futures
import random as _random_stdlib

np.random.seed(42)
_random_stdlib.seed(42)
os.environ.setdefault("PYTHONHASHSEED", "42")

# =============================================================================
# Configuration
# =============================================================================

BASE_DIR     = os.environ.get("SEAAD_ROOT", os.getcwd())
ARM_FOLDER   = r"[2]_Non-NeoC"
ARM_TAG      = "NonNeoC"
ARM_NAME     = "non-neocortical"
ARM_ROOT     = os.path.join(BASE_DIR, ARM_FOLDER)

FULL_H5AD    = os.path.join(BASE_DIR,
               "SEAAD_Immune_HIP_MEC_LEC_ITG_MTG_FI_STG_DFC_AnG_V1C_"
               "RNAseq_final-nuclei.2026-06-22.h5ad")

PART2_DIR_TAG  = "M4"
PART2_DIR      = os.path.join(ARM_ROOT, f"Part2_{ARM_TAG}_{PART2_DIR_TAG}")
PART2_ENV_NAME = "Part2_nonneocortex_env.joblib"
_JOBLIB_PATH   = os.path.join(PART2_DIR, PART2_ENV_NAME)

BASE_OUT = os.path.join(ARM_ROOT, f"Part6_{ARM_TAG}_GOI-DFG")

# == Target supertypes ========================================================
CANDIDATE_SUPERTYPES = [
    "Micro-PVM_3-SEAAD",
    "Micro-PVM_2_3-SEAAD",
    "Micro-PVM_4-SEAAD",
]
SHORT_ST = {
    "Micro-PVM_3-SEAAD"  : "PVM_3",
    "Micro-PVM_2_3-SEAAD": "PVM_2_3",
    "Micro-PVM_4-SEAAD"  : "PVM_4",
}

# == Regions ==================================================================
NON_NEOCORTEX_REGIONS = ["HIP", "MEC", "LEC", "FI"]

# == Sentinel for pooled (all-region) populations =============================
_POOLED_SENTINEL = "ALL"   # used in POPULATION_MAP region field

# == Populations ==============================================================
# PVM_3 by region; PVM_2_3 and PVM_4 pooled (too few cells per region).
ALL_POPULATIONS = (
    # PVM_3: by region
    [f"PVM_3_{reg}" for reg in NON_NEOCORTEX_REGIONS]
    +
    # PVM_2_3 and PVM_4: pooled
    ["PVM_2_3_pooled", "PVM_4_pooled"]
)

# == Thresholds ===============================================================
MIN_DONOR_POP_CELLS = 5      # min cells per donor within this population
MIN_POP_CELLS       = 10     # skip entire population if n_M4_cells < this
GOI_PCT_THRESHOLD   = 0.25   # GOI expressed in >= 25% of pooled cells
GOI_EXPRESSOR_MEAN  = 0.30   # mean expr among expressors >= 0.3
MIN_CELLS           = 10     # partner gene expressed in >= 10 cells
SP_THRESHOLD        = 0.6    # positive Spearman DFG threshold
KD_THRESHOLD        = 0.5    # positive Kendall DFG threshold
SP_NEG_THRESHOLD    = 0.6    # negative Spearman threshold (symmetric with positive)
KD_NEG_THRESHOLD    = 0.5    # negative Kendall threshold (symmetric with positive)
N_PARTS             = 2      # number of GOI list parts per population
N_WORKERS           = 28     # parallel processes

GLOBAL_SEED = 42

print("=" * 70)
print(f"Part 6  -  {ARM_NAME} M4 GOI-DFG Screening  [{ARM_TAG} arm]")
print(f"  Populations : {len(ALL_POPULATIONS)}")
print(f"  Base output : {BASE_OUT}/")
print(f"  Joblib      : {_JOBLIB_PATH}")
print("=" * 70)


# =============================================================================
# SECTION 1: Auto-loader, guardrails, adata load
# =============================================================================
print("\n" + "=" * 70)
print("SECTION 1: Auto-loader, guardrails, adata load")
print("=" * 70)

_REQUIRED_VARS = [
    "obs_full", "M4_TOP_MASK", "N_M4_TOP",
    "M4_THR_M1", "M4_THR_M2", "TOP_PCT",
    "score_cols", "donor_col", "stype_col", "region_col",
]
_missing_vars = [v for v in _REQUIRED_VARS if v not in dir()]

if _missing_vars:
    print(f"  [AUTO-LOADER] {len(_missing_vars)} variables missing "
          f" -  loading from joblib ...")
    import joblib as _jl6
    if not os.path.exists(_JOBLIB_PATH):
        raise FileNotFoundError(
            f"\n  FATAL: joblib not found: {_JOBLIB_PATH}\n"
            f"  Run Part 2 first.")
    _env6 = _jl6.load(_JOBLIB_PATH)
    obs_full    = _env6["obs_full"]
    M4_TOP_MASK = _env6["M4_TOP_MASK"]
    N_M4_TOP    = _env6["N_M4_TOP"]
    M4_THR_M1   = _env6["M4_THR_M1"]
    M4_THR_M2   = _env6["M4_THR_M2"]
    TOP_PCT     = _env6["TOP_PCT"]
    score_cols  = _env6["score_cols"]
    donor_col   = _env6["donor_col"]
    stype_col   = _env6["stype_col"]
    region_col  = _env6["region_col"]
    del _env6; gc.collect()
    print("  [AUTO-LOADER] All variables loaded successfully.")
else:
    print(f"  [AUTO-LOADER] All {len(_REQUIRED_VARS)} variables already "
          f"in namespace  -  skipping load.")

TOP_PCT_STR = f"{TOP_PCT * 100:.1f}%"

# == Guardrail 1: M4_TOP_MASK integrity =======================================
print("\n  [GUARDRAIL 1] M4_TOP_MASK integrity check")
assert isinstance(M4_TOP_MASK, np.ndarray) and M4_TOP_MASK.dtype == bool, (
    "M4_TOP_MASK must be a boolean numpy array.")
assert len(M4_TOP_MASK) == len(obs_full), (
    f"M4_TOP_MASK length {len(M4_TOP_MASK)} != obs_full {len(obs_full)} "
    f" -  row alignment broken.")
assert int(M4_TOP_MASK.sum()) == N_M4_TOP, (
    f"M4_TOP_MASK.sum()={M4_TOP_MASK.sum()} != N_M4_TOP={N_M4_TOP}")
print(f"    PASS  n={N_M4_TOP:,} cells | M1_thr={M4_THR_M1:.6f} | "
      f"M2_thr={M4_THR_M2:.6f} | TOP_PCT={TOP_PCT_STR}")

# == Guardrail 2: regions =====================================================
print("\n  [GUARDRAIL 2] Non-neocortical region check")
_expected_regions = set(NON_NEOCORTEX_REGIONS)
_obs_regions      = set(obs_full[region_col].astype(str).unique())
_unexpected       = _obs_regions - _expected_regions
if _unexpected:
    raise AssertionError(
        f"\n  GUARDRAIL FAILED: obs_full contains unexpected regions: "
        f"{_unexpected}\n"
        f"  obs_full must be the non-neocortical subset "
        f"(HIP/MEC/LEC/FI) only.")
_missing_regs = _expected_regions - _obs_regions
if _missing_regs:
    print(f"    WARNING: expected regions absent in obs_full: {_missing_regs}")
print(f"    PASS  obs_full regions = {sorted(_obs_regions)} "
      f"(n={len(obs_full):,} cells)")

# == Load adata (the h5ad; not in the joblib) =================================
print(f"\n  Loading adata: {FULL_H5AD}")
import scanpy as sc
adata = sc.read_h5ad(FULL_H5AD, backed="r")
print(f"  adata: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
if adata.X is None:
    raise AssertionError(
        "adata.X is empty.")
print(f"  adata.X confirmed present  -  shape: {adata.shape}, dtype: {adata.X.dtype}")

# == Map obs_full barcodes to adata row indices ===============================
print("\n  Building obs_full → adata global index map ...")
_all_barcodes_adata = np.array(adata.obs_names.tolist())
_barcode_to_global  = {b: i for i, b in enumerate(_all_barcodes_adata)}
_obs_full_barcodes  = obs_full.index.astype(str).tolist()
GLOBAL_IDX_OBSF     = np.array(
    [_barcode_to_global[b] for b in _obs_full_barcodes], dtype=np.int64)
_n_found = len(GLOBAL_IDX_OBSF)
_n_missing_bc = len(_obs_full_barcodes) - _n_found
if _n_missing_bc > 0:
    raise AssertionError(
        f"{_n_missing_bc} obs_full barcodes not found in adata.obs_names.")
print(f"  Mapped {_n_found:,} obs_full cells → adata global indices. "
      f"No missing barcodes.")

VAR_NAMES  = np.array(adata.var_names.tolist())   # gene names
N_GENES    = adata.n_vars
N_OBS_FULL = len(obs_full)
os.makedirs(BASE_OUT, exist_ok=True)


# =============================================================================
# SECTION 2: Protein-coding gene list (from h5ad var metadata, cached)
# =============================================================================
print("\n" + "=" * 70)
print("SECTION 2: Protein-coding gene list")
print("=" * 70)

def get_all_genes_with_biotype(
        h5ad_path, cache_path=os.path.join(BASE_OUT, "gene_biotype_map.pkl")):
    """Biotype label for every gene (annotation only; every gene is screened).
    Taken from a var biotype column if present, otherwise inferred from the
    gene name. Cached to pkl."""
    import re as _re

    os.makedirs(os.path.dirname(cache_path), exist_ok=True)
    if os.path.exists(cache_path):
        with open(cache_path, "rb") as f:
            result = pickle.load(f)
        if ("genes" in result and "biotype_map" in result
                and len(result["genes"]) > 5000):
            print(f"  Gene biotype map loaded from cache  "
                  f"({len(result['genes']):,} genes, "
                  f"method='{result.get('method', 'unknown')}'): "
                  f"{cache_path}")
            return result["genes"], result["biotype_map"]
        else:
            print("  Cache invalid  -  re-extracting ...")

    print("  Extracting gene list and biotype annotations from h5ad ...")
    import anndata as _ad
    _adata_tmp = _ad.read_h5ad(h5ad_path, backed="r")
    var        = _adata_tmp.var.copy()
    all_genes  = var.index.tolist()
    _adata_tmp.file.close()
    print(f"  Total genes in matrix : {len(all_genes):,}")
    print(f"  var columns           : {list(var.columns)}")

    # == Biotype column in var ================================================
    biotype_col = None
    for candidate in ["feature_biotype", "gene_biotype", "biotype",
                       "feature_type", "gene_type"]:
        if candidate in var.columns:
            biotype_col = candidate
            print(f"  Biotype column found: '{biotype_col}'")
            print(f"  Biotype value counts (top 12):\n"
                  f"{var[biotype_col].value_counts().head(12).to_string()}")
            break

    if biotype_col is not None:
        biotype_map = {g: str(var.loc[g, biotype_col]) for g in all_genes}
        method      = "var_biotype_column"
        _n_pc = sum(1 for v in biotype_map.values()
                    if "protein_coding" in v.lower())
        print(f"  Protein-coding genes : {_n_pc:,}")
        print(f"  All genes retained   : {len(all_genes):,} (lncRNAs included)")

    else:
        # == Otherwise, infer the biotype from the gene name ==================
        print(f"  No biotype column in var  -  inferring biotype from gene name ...")

        # Pattern -> biotype label (first match wins)
        _PATTERNS = [
            # lncRNA families
            (_re.compile(r'^LINC\d+'),                          "lncRNA_inferred"),
            (_re.compile(r'^(AL|AC|AP|BX|CR|AJ)\d{5,}'),       "lncRNA_inferred"),
            (_re.compile(r'^SNHG\d+'),                          "lncRNA_inferred"),
            (_re.compile(r'^PURPL|NORAD|NEAT1|MALAT1|XIST|KCNQ1OT1|MEG3|H19'),
                                                                "lncRNA_inferred"),
            # miRNA: pure digit suffix (MIR123, MIRLET7A1), not host genes
            (_re.compile(r'^MIR(LET)?\d+[A-Z]?\d*$'),          "miRNA_inferred"),
            # snoRNA / scaRNA
            (_re.compile(r'^(SNORD|SNORA|SCARNA)\w+'),          "snoRNA_inferred"),
            # pseudogenes (Ensembl-style uncharacterised patches)
            (_re.compile(r'^LOC\d+'),                           "pseudogene_inferred"),
            (_re.compile(r'^RP\d+[-C]\d+'),                     "pseudogene_inferred"),
            # mitochondrial rRNA
            (_re.compile(r'^MTRNR'),                            "rRNA_inferred"),
        ]

        def _infer_biotype(gene_name):
            for pat, label in _PATTERNS:
                if pat.match(gene_name):
                    return label
            return "protein_coding_inferred"

        biotype_map = {g: _infer_biotype(g) for g in all_genes}
        method      = "name_pattern_inference"

        # Summarise the inferred distribution
        from collections import Counter as _Counter
        _bt_counts = _Counter(biotype_map.values())
        print(f"  Inferred biotype distribution:")
        for _bt, _n in sorted(_bt_counts.items(), key=lambda x: -x[1]):
            print(f"    {_bt:<35}: {_n:>6,}")
        print(f"  All {len(all_genes):,} genes retained as GOI candidates "
              f"(lncRNAs included).")

    result = {"genes": all_genes, "biotype_map": biotype_map, "method": method}
    with open(cache_path, "wb") as f:
        pickle.dump(result, f, protocol=pickle.HIGHEST_PROTOCOL)
    print(f"  Gene biotype map cached ({method}): {cache_path}")
    return all_genes, biotype_map


ALL_GENES, BIOTYPE_MAP = get_all_genes_with_biotype(FULL_H5AD)

# ALL_GENES: every gene in the matrix; BIOTYPE_MAP: annotation only.
print(f"\n  GOI candidate list: {len(ALL_GENES):,} genes "
      f"(all biotypes screened; biotype used for annotation only)")


# =============================================================================
# SECTION 3: Expression extraction helpers
# =============================================================================
print("\n" + "=" * 70)
print("SECTION 3: Expression extraction helpers")
print("=" * 70)

def _read_rows_h5py(h5ad_path, row_indices, n_total_rows, n_genes_total):
    """Rows of adata.X (CP10K) from the h5ad CSR matrix, as dense float32."""
    row_indices = np.sort(np.asarray(row_indices, dtype=np.int64))
    n_rows  = len(row_indices)
    result  = np.zeros((n_rows, n_genes_total), dtype=np.float32)
    CHUNK   = 500

    with h5py.File(h5ad_path, "r") as hf:
        # Locate the CSR group
        layer_priority = ["X"]
        csr_group = None
        for lp in layer_priority:
            parts = lp.split("/")
            node  = hf
            try:
                for p in parts:
                    node = node[p]
                if all(k in node for k in ("data", "indices", "indptr")):
                    csr_group = lp
                    break
            except (KeyError, AttributeError):
                continue

        if csr_group is None:
            raise RuntimeError(
                "Cannot locate a valid CSR matrix in adata.X.")

        # Navigate to chosen group
        node = hf
        for p in csr_group.split("/"):
            node = node[p]
        data_ds    = node["data"]
        indices_ds = node["indices"]
        indptr_ds  = node["indptr"]

        for start in range(0, n_rows, CHUNK):
            end       = min(start + CHUNK, n_rows)
            rows_sub  = row_indices[start:end]
            row_min   = int(rows_sub.min())
            row_max   = int(rows_sub.max())
            indptr    = indptr_ds[row_min: row_max + 2]
            d_s       = int(indptr[0])
            d_e       = int(indptr[-1])
            if d_e > d_s:
                dv  = data_ds[d_s:d_e].astype(np.float32)
                ci  = indices_ds[d_s:d_e].astype(np.int32)
                for k in range(end - start):
                    local  = int(rows_sub[k]) - row_min
                    ps     = int(indptr[local])     - d_s
                    pe     = int(indptr[local + 1]) - d_s
                    if pe > ps:
                        result[start + k, ci[ps:pe]] = dv[ps:pe]

    return result


def _normalise_log1p(X_expr):
    """Cast to float32 (adata.X is already CP10K-normalised)."""
    return X_expr.astype(np.float32)


def _runtime_dtype_check(X_sample, label=""):
    """Guardrail 3: stop if adata.X looks like raw integer counts."""
    mx     = float(X_sample.max())
    mn_nz  = float(X_sample[X_sample > 0].min()) if (X_sample > 0).any() else 0.0
    is_int = np.allclose(X_sample, X_sample.astype(int), atol=0.01)
    if is_int and mx > 50:
        raise RuntimeError(
            f"\n  [GUARDRAIL 3] FATAL: adata.X looks like raw integer counts, "
            f"not pre-normalised expression [{label}].\n"
            f"  max value = {mx:.4f}  |  integer-like = {is_int}\n"
            f"  Expected CP10K-normalised floats (max ~10, min_nz ~0.5-1.0).")
    print(f"  [GUARDRAIL 3 PASS] adata.X dtype: max={mx:.2f}, "
          f"min_nz={mn_nz:.4f}, integer-like={is_int} [{label}]")


_DTYPE_CHECKED = False   # set True after first successful check


print("  Expression helpers ready.")


# =============================================================================
# SECTION 4: GOI gate check
# =============================================================================

def check_goi_gate(goi_vec, pop_label, logger=None):
    """GOI expression gate; the result is logged to file only."""
    n_cells  = len(goi_vec)
    n_expr   = int((goi_vec > 0).sum())
    pct_expr = n_expr / n_cells if n_cells > 0 else 0.0
    mean_pos = (float(goi_vec[goi_vec > 0].mean())
                if n_expr > 0 else 0.0)
    passes   = (pct_expr >= GOI_PCT_THRESHOLD) and (mean_pos >= GOI_EXPRESSOR_MEAN)
    msg      = (f"    Gate [{pop_label}]: cells={n_cells} | "
                f"expr={n_expr} ({100*pct_expr:.1f}%) | "
                f"mean(>0)={mean_pos:.3f} | "
                f"{'PASS' if passes else 'SKIP'}")
    if logger:
        logger.debug(msg)   # file only
    return passes, pct_expr, mean_pos


# =============================================================================
# SECTION 5: Vectorised Spearman
# =============================================================================

def run_spearman(expr_partner, goi_vec, goi, pop_label, partner_genes, logger):
    """Vectorised Spearman (Pearson of ranks), GOI against all partner genes.
    Returns the passers: Gene, Spearman_Correlation, Direction."""
    from scipy.stats import rankdata

    n_cells, n_partner = expr_partner.shape
    logger.debug(f"    [Spearman] {pop_label}: {n_partner:,} genes × "
                 f"{n_cells} cells ...")
    t0 = time.time()

    goi_ranks  = rankdata(goi_vec).astype(np.float32)
    goi_ranks -= goi_ranks.mean()
    goi_norm   = float(np.linalg.norm(goi_ranks))
    if goi_norm == 0:
        logger.debug("  GOI has zero variance  -  skipping")
        return pd.DataFrame(columns=["Gene", "Spearman_Correlation", "Direction"])

    CHUNK    = 3000
    sp_vals  = np.empty(n_partner, dtype=np.float32)

    for start in range(0, n_partner, CHUNK):
        end   = min(start + CHUNK, n_partner)
        chunk = expr_partner[:, start:end].astype(np.float64)

        from scipy.stats import rankdata as _rd
        chunk_ranked = np.apply_along_axis(_rd, 0, chunk).astype(np.float32)
        chunk_ranked -= chunk_ranked.mean(axis=0)
        chunk_norms   = np.linalg.norm(chunk_ranked, axis=0)
        chunk_norms[chunk_norms == 0] = 1.0

        sp_vals[start:end] = (
            (goi_ranks @ chunk_ranked) / (goi_norm * chunk_norms)
        ).astype(np.float32)

    gene_arr = np.array(partner_genes)
    mask_pos = (sp_vals >=  SP_THRESHOLD)      & (gene_arr != goi)
    mask_neg = (sp_vals <= -SP_NEG_THRESHOLD)  & (gene_arr != goi)
    mask     = mask_pos | mask_neg

    sp_df = (pd.DataFrame({
                 "Gene":                 gene_arr[mask],
                 "Spearman_Correlation": sp_vals[mask].astype(float),
                 "Direction": np.where(mask_pos[mask], "POS", "NEG"),
             })
             .sort_values("Spearman_Correlation", ascending=False)
             .reset_index(drop=True))

    n_pos = int(mask_pos.sum()); n_neg = int(mask_neg.sum())
    logger.debug(f"    [Spearman done] pos(>={SP_THRESHOLD}): {n_pos} | "
                 f"neg(<=-{SP_NEG_THRESHOLD}): {n_neg} | "
                 f"[{time.time()-t0:.1f}s]")
    return sp_df


# =============================================================================
# SECTION 6: Sequential Kendall (Spearman passers only)
# =============================================================================

def run_kendall(expr_kd, goi_vec, pop_label, sp_genes, logger):
    """Kendall tau for the Spearman passers: Gene, Kendall_Correlation,
    Direction."""
    from scipy.stats import kendalltau as _kt
    n_kd = expr_kd.shape[1]
    logger.debug(f"    [Kendall] {pop_label}: {n_kd} genes × "
                 f"{len(goi_vec)} cells ...")
    t0 = time.time()

    goi_f64  = goi_vec.astype(np.float64)
    kd_vals  = np.empty(n_kd, dtype=np.float32)
    for i in range(n_kd):
        vec = expr_kd[:, i].astype(np.float64)
        if vec.std() == 0:
            kd_vals[i] = 0.0
        else:
            kd_vals[i] = float(_kt(goi_f64, vec, method="auto").statistic)

    sp_genes_arr = np.array(sp_genes)
    mask_pos = kd_vals >=  KD_THRESHOLD
    mask_neg = kd_vals <= -KD_NEG_THRESHOLD
    mask     = mask_pos | mask_neg

    kd_df = (pd.DataFrame({
                 "Gene":                sp_genes_arr[mask],
                 "Kendall_Correlation": kd_vals[mask].astype(float),
                 "Direction": np.where(mask_pos[mask], "POS", "NEG"),
             })
             .sort_values("Kendall_Correlation", ascending=False)
             .reset_index(drop=True))

    n_pos = int(mask_pos.sum()); n_neg = int(mask_neg.sum())
    logger.debug(f"    [Kendall done] pos(>={KD_THRESHOLD}): {n_pos} | "
                 f"neg(<=-{KD_NEG_THRESHOLD}): {n_neg} | "
                 f"[{time.time()-t0:.1f}s]")
    return kd_df


# =============================================================================
# SECTION 7: Dual filter (inner join + sign consistency)
# =============================================================================

def dual_filter(sp_df, kd_df, pop_label, logger):
    """Genes passing both Spearman and Kendall with the same sign."""
    if sp_df.empty or kd_df.empty:
        logger.debug(f"    DualFilter [{pop_label}]: 0 DFGs")
        return pd.DataFrame(columns=["Gene", "Spearman_Correlation",
                                      "Kendall_Correlation", "Direction", "State"])
    hits = (sp_df[["Gene", "Spearman_Correlation", "Direction"]]
            .merge(kd_df[["Gene", "Kendall_Correlation"]], on="Gene")
            .assign(State=pop_label))
    sign_ok = (
        ((hits["Direction"] == "POS") & (hits["Kendall_Correlation"] >= 0)) |
        ((hits["Direction"] == "NEG") & (hits["Kendall_Correlation"] <= 0))
    )
    hits = (hits[sign_ok]
            .sort_values("Spearman_Correlation", ascending=False)
            .reset_index(drop=True))
    n_pos = int((hits["Direction"] == "POS").sum())
    n_neg = int((hits["Direction"] == "NEG").sum())
    logger.debug(f"    DualFilter [{pop_label}]: {len(hits)} DFGs "
                f"(POS={n_pos}, NEG={n_neg})")
    return hits


# =============================================================================
# SECTION 8: GOI expression statistics (MasterDFG metadata)
# =============================================================================

def goi_expr_stats(goi_vec):
    """GOI expression summary for the MasterDFG table."""
    pos = goi_vec[goi_vec > 0]
    return {
        "Pct_expressing": round(float((goi_vec > 0).mean()), 4),
        "Mean_all"      : round(float(goi_vec.mean()), 4),
        "Median_all"    : round(float(np.median(goi_vec)), 4),
        "Q1_all"        : round(float(np.percentile(goi_vec, 25)), 4),
        "Q3_all"        : round(float(np.percentile(goi_vec, 75)), 4),
        "IQR_all"       : round(float(np.percentile(goi_vec, 75) -
                                       np.percentile(goi_vec, 25)), 4),
        "Mean_pos"      : round(float(pos.mean())     if len(pos) > 0 else 0.0, 4),
        "Median_pos"    : round(float(np.median(pos)) if len(pos) > 0 else 0.0, 4),
    }


# =============================================================================
# SECTION 9: Workbook writers
# =============================================================================

def _df_to_sheet(ws, df):
    """Write a DataFrame to an openpyxl worksheet."""
    ws.append(list(df.columns))
    for row in df.itertuples(index=False):
        ws.append(list(row))


def write_goi_excels(goi, sp_df, kd_df, hits_df, pop_label, out_dir):
    """Write the Spearman, Kendall and DualFilter workbooks for one GOI."""
    # (A) Spearman
    wb_sp = Workbook(); ws = wb_sp.active; ws.title = pop_label[:31]
    sp_cols_out = [c for c in ["Gene","Spearman_Correlation","Direction"]
                   if c in sp_df.columns]
    _df_to_sheet(ws, sp_df.assign(State=pop_label)[sp_cols_out + ["State"]])
    wb_sp.save(os.path.join(out_dir, f"{goi}_Spearman.xlsx"))

    # (B) Kendall
    wb_kd = Workbook(); ws = wb_kd.active; ws.title = pop_label[:31]
    if not kd_df.empty:
        kd_cols_out = [c for c in ["Gene","Kendall_Correlation","Direction"]
                       if c in kd_df.columns]
        _df_to_sheet(ws, kd_df.assign(State=pop_label)[kd_cols_out + ["State"]])
    else:
        ws.append(["Gene", "Kendall_Correlation", "State"])
    wb_kd.save(os.path.join(out_dir, f"{goi}_Kendall.xlsx"))

    # (C) DualFilter
    wb_df = Workbook(); ws_all = wb_df.active; ws_all.title = "All_States_Combined"
    cols = ["Gene","Spearman_Correlation","Kendall_Correlation","Direction","State"]
    if not hits_df.empty:
        _df_to_sheet(ws_all, hits_df[cols])
        ws_pop = wb_df.create_sheet(pop_label[:31])
        _df_to_sheet(ws_pop, hits_df[cols])
    else:
        ws_all.append(cols)
    wb_df.save(os.path.join(out_dir, f"{goi}_DualFilter.xlsx"))


# =============================================================================
# SECTION 10: MasterDFG workbook (cumulative per population)
# =============================================================================

# =============================================================================
# GOI screener for the fork workers
# =============================================================================
# The arrays are set in _MP before the pool starts; fork workers inherit them.

_MP = {}

def _mp_screen_one_goi(goi):
    """Screen a single GOI using module-level shared state. Process-safe."""
    expr_mat            = _MP["expr_mat"]
    expr_partner        = _MP["expr_partner"]
    partner_genes       = _MP["partner_genes"]
    var_name_to_idx     = _MP["var_name_to_idx"]
    partner_gene_to_pos = _MP["partner_gene_to_pos"]
    pop_name            = _MP["pop_name"]
    out_dir             = _MP["out_dir"]

    _null_log = logging.getLogger("_mp_null")
    if not _null_log.handlers:
        _null_log.addHandler(logging.NullHandler())
        _null_log.setLevel(logging.WARNING)

    goi_global_idx = var_name_to_idx[goi]
    goi_vec = expr_mat[:, goi_global_idx].astype(np.float32)

    passes, pct_exp, mean_pos_val = check_goi_gate(goi_vec, pop_name, _null_log)
    e_stats   = goi_expr_stats(goi_vec)
    dfg_total = 0
    max_sp    = 0.0
    max_kd    = 0.0

    if passes:
        sp_df = run_spearman(
            expr_partner, goi_vec, goi, pop_name,
            partner_genes, _null_log)

        if len(sp_df) > 0:
            sp_genes    = sp_df["Gene"].tolist()
            kd_col_idx  = np.array([partner_gene_to_pos[g]
                                    for g in sp_genes
                                    if g in partner_gene_to_pos],
                                   dtype=np.int64)
            sp_genes_kd = [partner_genes[i] for i in kd_col_idx]
            expr_kd     = expr_partner[:, kd_col_idx]

            kd_df = run_kendall(expr_kd, goi_vec, pop_name,
                                sp_genes_kd, _null_log)
            hits  = dual_filter(sp_df, kd_df, pop_name, _null_log)

            dfg_total = len(hits)
            max_sp    = float(sp_df["Spearman_Correlation"].max()) \
                        if len(sp_df) > 0 else 0.0
            max_kd    = float(kd_df["Kendall_Correlation"].max()) \
                        if len(kd_df) > 0 else 0.0

            if dfg_total > 0:
                write_goi_excels(goi, sp_df, kd_df, hits,
                                 pop_name, out_dir)

    master_row = {
        "GOI"         : goi,
        "DFG_Total"   : dfg_total,
        "max_Spearman": round(max_sp, 4),
        "max_Kendall" : round(max_kd, 4),
    }
    master_row.update(e_stats)
    return master_row, dfg_total, max_sp, max_kd, passes, goi


def update_master_excel(out_dir, n_parts, pop_label):
    """Read all checkpoint pkls, rebuild and save MasterDFG Excel."""
    all_rows = []
    for p in range(1, n_parts + 1):
        pkl = os.path.join(out_dir, f"Part{p:02d}_checkpoint.pkl")
        if os.path.exists(pkl):
            with open(pkl, "rb") as f:
                all_rows.extend(pickle.load(f)["master_rows"])

    if not all_rows:
        return

    master_df = (pd.DataFrame(all_rows)
                 .query("DFG_Total >= 1")
                 .sort_values(["DFG_Total", "max_Spearman"],
                               ascending=[False, False])
                 .reset_index(drop=True))
    if master_df.empty:
        print(f"  MasterDFG [{pop_label}]: 0 GOIs with DFG hits.")
        return

    master_df.insert(0, "Rank", master_df.index + 1)
    master_df["Biotype"] = master_df["GOI"].map(BIOTYPE_MAP).fillna("unknown")

    col_order = ["Rank", "GOI", "Biotype", "DFG_Total",
                 "max_Spearman", "max_Kendall",
                 "Pct_expressing", "Mean_all", "Median_all",
                 "Q1_all", "Q3_all", "IQR_all", "Mean_pos", "Median_pos"]
    master_df = master_df[[c for c in col_order if c in master_df.columns]]

    master_path = os.path.join(out_dir, f"{pop_label}_MasterDFG_Ranking.xlsx")
    wb = Workbook(); ws = wb.active; ws.title = "Master_DFG_Ranking"
    _df_to_sheet(ws, master_df)
    wb.save(master_path)
    print(f"  MasterDFG [{pop_label}]: {len(master_df)} GOIs with DFG hits  -  "
          f"{os.path.basename(master_path)}")


# =============================================================================
# SECTION 11: Logger setup
# =============================================================================

def setup_logger(log_path):
    """Logger: DEBUG and above to the log file, INFO and above to stdout."""
    logger = logging.getLogger(log_path)
    logger.setLevel(logging.DEBUG)   # capture everything at logger level
    if logger.handlers:
        logger.handlers.clear()

    fmt = logging.Formatter("%(message)s")

    # File: DEBUG and above
    fh = logging.FileHandler(log_path, mode="a", encoding="utf-8")
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(fmt)

    # Stdout: INFO and above
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    ch.setFormatter(fmt)

    logger.addHandler(fh)
    logger.addHandler(ch)
    return logger


# =============================================================================
# SECTION 12: Population definition, cell count tally, and independence check
# =============================================================================
print("\n" + "=" * 70)
print("SECTION 12: Population catalogue and cell count tally")
print("=" * 70)

# == Population map: name -> (supertype, region or the pooled sentinel) =======
POPULATION_MAP = {}
for reg in NON_NEOCORTEX_REGIONS:
    POPULATION_MAP[f"PVM_3_{reg}"] = ("Micro-PVM_3-SEAAD", reg)
POPULATION_MAP["PVM_2_3_pooled"] = ("Micro-PVM_2_3-SEAAD", _POOLED_SENTINEL)
POPULATION_MAP["PVM_4_pooled"]   = ("Micro-PVM_4-SEAAD",   _POOLED_SENTINEL)

# Label arrays from obs_full
cell_labels_obs       = obs_full[stype_col].astype(str).values
region_vals_obs       = obs_full[region_col].astype(str).values
obs_full_barcodes_arr = obs_full.index.astype(str).values

# =============================================================================
# Guardrail 4: M4-top cells per supertype and region
# =============================================================================
print("\n  [GUARDRAIL 4] M4-TOP CELL COUNT TALLY")
print("  PVM_3: by region | PVM_2_3, PVM_4: pooled")
print("  Cross-tabulation: canonical M4_TOP_MASK cells only")

# Build the full tally table (per region, for all supertypes)
_tally_regions  = NON_NEOCORTEX_REGIONS
_tally_st_short = [SHORT_ST[st] for st in CANDIDATE_SUPERTYPES]
_tally_st_full  = CANDIDATE_SUPERTYPES

_col_w  = 10; _lbl_w = 12
_hdr_r  = "  ".join(f"{r:>{_col_w}}" for r in _tally_regions)
print(f"  {'Supertype':<{_lbl_w}}  {_hdr_r}  {'Pooled':>{_col_w}}  {'Screened_as'}")
print("  " + "-" * (_lbl_w + 2 + (_col_w + 2)*len(_tally_regions) + _col_w + 18))

_col_totals = {r: 0 for r in _tally_regions}; _col_totals["Pooled"] = 0
_tally_data = {}   # (st_full, reg) -> n_m4
_pooled_tally = {}   # st_full -> M4 cells over all regions

for st_full, short in zip(_tally_st_full, _tally_st_short):
    st_mask_t = (cell_labels_obs == st_full)
    row_counts = []; row_total = 0
    for reg in _tally_regions:
        n_t = int((M4_TOP_MASK & st_mask_t & (region_vals_obs == reg)).sum())
        _tally_data[(st_full, reg)] = n_t
        _col_totals[reg] += n_t
        row_counts.append(n_t); row_total += n_t
    _col_totals["Pooled"] += row_total
    _pooled_tally[st_full] = row_total
    _tally_data[(st_full, _POOLED_SENTINEL)] = row_total   # sentinel key for pooled pops
    _row_str = "  ".join(f"{n:>{_col_w},}" for n in row_counts)
    # Annotate screening strategy
    _strat = "by region" if short == "PVM_3" else "POOLED"
    print(f"  {short:<{_lbl_w}}  {_row_str}  {row_total:>{_col_w},}  [{_strat}]")

print("  " + "-" * (_lbl_w + 2 + (_col_w + 2)*len(_tally_regions) + _col_w + 18))
_tot_str = "  ".join(f"{_col_totals[r]:>{_col_w},}" for r in _tally_regions)
print(f"  {'ALL ST':<{_lbl_w}}  {_tot_str}  {_col_totals['Pooled']:>{_col_w},}")

# M4 total cross-check
print(f"\n  M4_TOP_MASK total (from Part 2 joblib): {N_M4_TOP:,} cells")
_sum_3st = sum(_tally_data.get((st, _POOLED_SENTINEL), 0) for st in CANDIDATE_SUPERTYPES)
print(f"  Sum across 3 supertypes (all regions) : {_sum_3st:,} cells")
_diff = N_M4_TOP - _sum_3st
if _diff == 0:
    print(f"  [TALLY CHECK PASS] All {N_M4_TOP:,} M4-top cells accounted for.")
else:
    print(f"  [TALLY NOTE] {_diff:,} M4-top cells belong to other supertypes "
          f"not in CANDIDATE_SUPERTYPES.")

# Population-level breakdown matching the actual screening jobs
print(f"\n  POPULATION BREAKDOWN ({len(ALL_POPULATIONS)} screening jobs):")
print(f"  {'Population':<22}  {'Strategy':<12}  {'Cells':>7}  Notes")
print("  " + "-" * 65)
for pop in ALL_POPULATIONS:
    st_f, rg = POPULATION_MAP[pop]
    n_pop = _tally_data.get((st_f, rg), 0)
    if rg == _POOLED_SENTINEL:
        strategy = "POOLED"
        note = f"all 4 regions combined, {SHORT_ST[st_f]} only"
    else:
        strategy = f"region={rg}"
        note = f"only {rg} cells, {SHORT_ST[st_f]} only"
    print(f"    {pop:<22}  {strategy:<12}  {n_pop:>7,}  {note}")

# =============================================================================
# Per-population eligibility check (donor filter + min cells)
# =============================================================================
print(f"\n  {'─'*70}")
print(f"  ELIGIBILITY CHECK  "
      f"(min_pop_cells={MIN_POP_CELLS}, min_donor_cells={MIN_DONOR_POP_CELLS})")
print(f"  {'─'*70}")
print(f"\n  {'Population':<22} {'Strategy':<12} {'ST_short':<10} "
      f"{'M4_cells':>9} {'after_donor_filt':>17} {'elig_donors':>13}  status")
print("  " + "-" * 100)

POPULATION_VALID = {}
_tally_rows_csv  = []

for pop_name in ALL_POPULATIONS:
    st_full, reg = POPULATION_MAP[pop_name]
    n_m4 = _tally_data.get((st_full, reg), 0)
    is_pooled = (reg == _POOLED_SENTINEL)
    strategy  = "POOLED" if is_pooled else f"region={reg}"

    if n_m4 < MIN_POP_CELLS:
        status = f"SKIP (n_M4={n_m4} < {MIN_POP_CELLS})"
        POPULATION_VALID[pop_name] = False
        n_after_donor = 0; n_elig = 0
    else:
        # Pooled: all regions; regional: one region
        if is_pooled:
            _pm_prev = (M4_TOP_MASK & (cell_labels_obs == st_full))
        else:
            _pm_prev = (M4_TOP_MASK & (cell_labels_obs == st_full)
                        & (region_vals_obs == reg))
        _dc_prev   = obs_full.loc[_pm_prev, donor_col].value_counts()
        _elig_prev = _dc_prev[_dc_prev >= MIN_DONOR_POP_CELLS]
        n_elig        = len(_elig_prev)
        n_after_donor = int(_elig_prev.sum())
        if n_elig == 0:
            status = f"SKIP (0 donors >={MIN_DONOR_POP_CELLS} cells)"
            POPULATION_VALID[pop_name] = False
        else:
            status = "OK"
            POPULATION_VALID[pop_name] = True

    short_st = SHORT_ST.get(st_full, st_full[:10])
    _flag = "✓" if POPULATION_VALID[pop_name] else "✗"
    print(f"  {_flag} {pop_name:<21} {strategy:<12} {short_st:<10} "
          f"{n_m4:>9,} {n_after_donor:>17,} {n_elig:>13}  {status}")

    _tally_rows_csv.append({
        "population"                : pop_name,
        "supertype_full"            : st_full,
        "supertype_short"           : SHORT_ST.get(st_full, ""),
        "region"                    : "ALL_pooled" if is_pooled else reg,
        "screening_strategy"        : strategy,
        "n_M4_cells_raw"            : n_m4,
        "n_cells_after_donor_filter": n_after_donor,
        "n_eligible_donors"         : n_elig,
        "min_pop_cells"             : MIN_POP_CELLS,
        "min_donor_cells"           : MIN_DONOR_POP_CELLS,
        "status"                    : "OK" if POPULATION_VALID[pop_name] else "SKIPPED",
    })

_tally_df = pd.DataFrame(_tally_rows_csv)
_tally_csv = os.path.join(BASE_OUT, "Part6_population_tally.csv")
os.makedirs(BASE_OUT, exist_ok=True)
_tally_df.to_csv(_tally_csv, index=False)
print(f"\n  Tally saved: {_tally_csv}")

n_valid_pops   = sum(POPULATION_VALID.values())
n_skipped_pops = len(ALL_POPULATIONS) - n_valid_pops
print(f"\n  Valid populations  : {n_valid_pops} / {len(ALL_POPULATIONS)}")
print(f"  Skipped populations: {n_skipped_pops} / {len(ALL_POPULATIONS)}")
print("\n  [GUARDRAIL 4 COMPLETE]")


# =============================================================================
# SECTION 13: Main screening loop  -  one population at a time
# =============================================================================

def run_population(pop_name):
    """GOI-DFG screen for one population: M4-top cells of one supertype in one
    region (or all regions when pooled), donors with enough cells."""
    global _DTYPE_CHECKED

    np.random.seed(GLOBAL_SEED)
    _random_stdlib.seed(GLOBAL_SEED)

    if not POPULATION_VALID.get(pop_name, False):
        print(f"\n  [{pop_name}] SKIPPED (see catalogue above).")
        return

    st_full, reg = POPULATION_MAP[pop_name]
    is_pooled    = (reg == _POOLED_SENTINEL)
    out_dir = os.path.join(BASE_OUT, pop_name)
    os.makedirs(out_dir, exist_ok=True)

    _reg_label = "ALL non-neocortical (HIP+MEC+LEC+FI)" if is_pooled else reg
    print("\n" + "=" * 70)
    print(f"POPULATION: {pop_name}  ({st_full} | {_reg_label})")
    print("=" * 70)

    # == Step A: cells from M4_TOP_MASK =======================================
    st_mask  = (cell_labels_obs == st_full)
    if is_pooled:
        pop_mask = M4_TOP_MASK & st_mask          # all regions
    else:
        reg_mask  = (region_vals_obs == reg)
        pop_mask  = M4_TOP_MASK & st_mask & reg_mask

    obs_pop        = obs_full[pop_mask].copy()
    global_idx_pop = GLOBAL_IDX_OBSF[pop_mask]     # adata row indices
    n_cells_raw    = len(obs_pop)

    print(f"  M4-top cells (before donor filter): {n_cells_raw:,}  "
          f"[{_reg_label}]")
    print(f"  Canonical M4_TOP_MASK thresholds: "
          f"M1>={M4_THR_M1:.6f}, M2>={M4_THR_M2:.6f}")

    # == Step B: donor filter =================================================
    donor_counts_pop = obs_pop[donor_col].value_counts()
    eligible_donors  = donor_counts_pop[donor_counts_pop >= MIN_DONOR_POP_CELLS].index
    donor_mask_pop   = obs_pop[donor_col].isin(eligible_donors).values
    obs_pop          = obs_pop[donor_mask_pop].copy()
    global_idx_pop   = global_idx_pop[donor_mask_pop]
    n_cells          = len(obs_pop)
    n_excl_donor     = n_cells_raw - n_cells

    print(f"  Eligible donors (>={MIN_DONOR_POP_CELLS} cells): "
          f"{len(eligible_donors)} | "
          f"Cells after donor filter: {n_cells:,} "
          f"({n_excl_donor:,} excluded)")

    if n_cells < MIN_POP_CELLS:
        print(f"  SKIP: fewer than {MIN_POP_CELLS} cells after donor filter.")
        return

    # == Step C: load expression ==============================================
    print(f"\n  Loading expression for {n_cells:,} cells ...")
    t0_load = time.time()
    expr_mat = _read_rows_h5py(FULL_H5AD, global_idx_pop, adata.n_obs, N_GENES)
    print(f"  Expression loaded: {expr_mat.shape}  [{time.time()-t0_load:.1f}s]")

    # Guardrail 3: CP10K check (once)
    if not _DTYPE_CHECKED:
        _runtime_dtype_check(expr_mat[:min(50, n_cells), :], label=pop_name)
        _DTYPE_CHECKED = True

    # == Step D: cast to float32 ==============================================
    expr_mat = _normalise_log1p(expr_mat)   # passthrough cast to float32

    # == Step E: GOI candidates and partner genes =============================
    _var_name_set     = set(VAR_NAMES.tolist())
    _var_name_col_idx = {g: i for i, g in enumerate(VAR_NAMES.tolist())}

    all_in_mat  = [g for g in ALL_GENES if g in _var_name_set]
    all_idx     = np.array([_var_name_col_idx[g] for g in all_in_mat], dtype=np.int64)
    n_all_genes = len(all_in_mat)
    print(f"  GOI candidates in matrix : {n_all_genes:,} "
          f"(all biotypes; protein-coding + lncRNA + other)")

    # Partner gene filter: expressed in >= MIN_CELLS cells
    gene_expr_counts = (expr_mat[:, all_idx] > 0).sum(axis=0)
    keep_partner     = gene_expr_counts >= MIN_CELLS
    partner_genes    = [all_in_mat[i] for i in range(n_all_genes)
                        if keep_partner[i]]
    partner_idx      = all_idx[keep_partner]
    n_partner        = len(partner_genes)
    print(f"  Partner genes (>={MIN_CELLS} cells expressing): "
          f"{n_partner:,} / {n_all_genes:,}")

    # Partner expression matrix (used for Spearman loop)
    expr_partner = expr_mat[:, partner_idx]   # (n_cells, n_partner)
    gc.collect()

    # GOI list = all genes in matrix (same as all_in_mat)
    goi_list   = all_in_mat
    total_goi  = len(goi_list)
    part_size  = (total_goi + N_PARTS - 1) // N_PARTS
    gene_parts = [goi_list[i * part_size: (i + 1) * part_size]
                  for i in range(N_PARTS)]

    print(f"\n  GOI list: {total_goi:,} | split into {N_PARTS} parts")
    for p_i, gp in enumerate(gene_parts):
        print(f"    Part {p_i+1:02d}: {len(gp):,} genes  "
              f"({gp[0]} … {gp[-1]})")

    # == Step F: parts, with resume ===========================================
    # Lookups: partner gene -> column of expr_partner; gene -> adata column.
    partner_gene_to_pos = {g: i for i, g in enumerate(partner_genes)}
    var_name_to_idx     = {g: int(_var_name_col_idx[g]) for g in all_in_mat}

    for part_idx in range(1, N_PARTS + 1):
        pkl_file = os.path.join(out_dir, f"Part{part_idx:02d}_checkpoint.pkl")
        log_file = os.path.join(out_dir, f"Part{part_idx:02d}_log.txt")

        if os.path.exists(pkl_file):
            print(f">>> Part {part_idx:02d}/{N_PARTS}  -  already done "
                  f"(checkpoint found), skipping.")
            continue

        part_genes = gene_parts[part_idx - 1]
        n_part     = len(part_genes)
        logger     = setup_logger(log_file)

        logger.info(f"=== Part {part_idx:02d}/{N_PARTS} | {n_part:,} GOIs | "
                    f"{pop_name} | {time.strftime('%Y-%m-%d %H:%M:%S')} ===")
        logger.info(f"SP_THRESHOLD={SP_THRESHOLD}  KD_THRESHOLD={KD_THRESHOLD}  "
                    f"SEED={GLOBAL_SEED}\n")

        part_master_rows = []
        t_part    = time.time()
        n_passed  = 0   # GOIs that passed the expression gate
        n_dfg_hit = 0   # GOIs with >= 1 DFG found so far in this part
        _PROGRESS_INTERVAL = 500

        # Shared state for the fork workers
        global _MP
        _MP["expr_mat"]            = expr_mat
        _MP["expr_partner"]        = expr_partner
        _MP["partner_genes"]       = partner_genes
        _MP["var_name_to_idx"]     = var_name_to_idx
        _MP["partner_gene_to_pos"] = partner_gene_to_pos
        _MP["pop_name"]            = pop_name
        _MP["out_dir"]             = out_dir

        import multiprocessing as _mp
        _ctx = _mp.get_context("fork")
        print(f"  Using {N_WORKERS} processes for parallel GOI screening (fork)")

        with concurrent.futures.ProcessPoolExecutor(
                max_workers=N_WORKERS, mp_context=_ctx) as executor:
            results_iter = executor.map(_mp_screen_one_goi, part_genes,
                                        chunksize=100)

            for g_idx, result in enumerate(results_iter):
                master_row, dfg_total, max_sp, max_kd, passed, goi = result
                part_master_rows.append(master_row)
                if passed:
                    n_passed += 1
                if dfg_total > 0:
                    n_dfg_hit += 1
                    print(f"  [DFG HIT] {pop_name} | GOI={goi} | "
                          f"DFGs={dfg_total} | "
                          f"max_Sp={max_sp:.3f} | max_Kd={max_kd:.3f}",
                          flush=True)

                _completed = g_idx + 1
                if _completed % _PROGRESS_INTERVAL == 0 or _completed == n_part:
                    elapsed = time.time() - t_part
                    rate = elapsed / _completed
                    eta_min = rate * (n_part - _completed) / 60
                    print(f"  [{pop_name} | Part {part_idx:02d}] "
                          f"{_completed}/{n_part} GOIs | "
                          f"{elapsed/60:.1f} min | "
                          f"gate_passed={n_passed} | DFG_hits={n_dfg_hit} | "
                          f"ETA ~{eta_min:.0f} min",
                          flush=True)

        # Save checkpoint
        with open(pkl_file, "wb") as f:
            pickle.dump({"master_rows": part_master_rows}, f,
                        protocol=pickle.HIGHEST_PROTOCOL)
        logger.info(
            f"\n=== Part {part_idx:02d} complete | {n_part:,} GOIs | "
            f"saved: {os.path.basename(pkl_file)} ===")
        print(f">>> Part {part_idx:02d}/{N_PARTS} done  -  "
              f"checkpoint: {os.path.basename(pkl_file)}")

        # == Step G: update the MasterDFG =====================================
        update_master_excel(out_dir, N_PARTS, pop_name)

    # == Final MasterDFG ======================================================
    print(f"\n  Final MasterDFG consolidation for {pop_name} ...")
    update_master_excel(out_dir, N_PARTS, pop_name)

    # Population-level summary printout
    all_rows_pop = []
    for p in range(1, N_PARTS + 1):
        pkl = os.path.join(out_dir, f"Part{p:02d}_checkpoint.pkl")
        if os.path.exists(pkl):
            with open(pkl, "rb") as f:
                all_rows_pop.extend(pickle.load(f)["master_rows"])
    if all_rows_pop:
        mdf_pop          = pd.DataFrame(all_rows_pop)
        total_screened   = len(mdf_pop)
        total_hits       = int((mdf_pop["DFG_Total"] >= 1).sum())
        print(f"\n  GOIs screened      : {total_screened:,}")
        print(f"  GOIs with >= 1 DFG : {total_hits:,}")
        top10 = (mdf_pop.query("DFG_Total >= 1")
                 .sort_values(["DFG_Total", "max_Spearman"],
                               ascending=[False, False])
                 .head(10))
        if len(top10) > 0:
            print(f"\n  Top 10 GOIs by DFG_Total:")
            print(top10[["GOI","DFG_Total","max_Spearman","max_Kendall"]
                        ].to_string(index=False))

    print("\n" + "=" * 70)
    print(f"POPULATION {pop_name} COMPLETE  → {out_dir}/")
    print("=" * 70)


# =============================================================================
# SECTION 14: Cross-population summary workbook
# =============================================================================

def write_cross_population_summary():
    """Cross-population summary workbook: cell counts, population catalogue,
    DFG counts, top 50 GOIs, and the top 100 GOIs of each population."""
    fp_summary = os.path.join(BASE_OUT, "Part6_cross_population_summary.xlsx")
    wb          = Workbook(); wb.remove(wb.active)

    # Style helpers
    _hf  = Font(name="Calibri", bold=True, size=10, color="FFFFFF")
    _hfl = PatternFill("solid", fgColor="1F3864")
    _bf  = Font(name="Calibri", size=9)
    _af  = PatternFill("solid", fgColor="EFF3FB")
    _tn  = Side(style="thin", color="BBBBBB")
    _bdr = Border(left=_tn, right=_tn, top=_tn, bottom=_tn)
    _ctr = Alignment(horizontal="center", vertical="center")
    _lft = Alignment(horizontal="left",   vertical="center")

    def _write_header(ws, headers, widths, row=1):
        for ci, (h, w) in enumerate(zip(headers, widths), 1):
            c = ws.cell(row=row, column=ci, value=h)
            c.font = _hf; c.fill = _hfl; c.border = _bdr; c.alignment = _ctr
            ws.column_dimensions[get_column_letter(ci)].width = w
        ws.row_dimensions[row].height = 24

    def _safe_cell(ws, ri, ci, val):
        c = ws.cell(ri, ci)
        try:
            if isinstance(val, (bool, np.bool_)):         c.value = bool(val)
            elif isinstance(val, (float, np.floating)):
                c.value = None if (np.isnan(val) or np.isinf(val)) \
                          else round(float(val), 6)
            elif isinstance(val, (int, np.integer)):      c.value = int(val)
            else:                                          c.value = str(val)
        except Exception:                                  c.value = str(val)
        c.font = _bf; c.border = _bdr
        return c

    # == Sheet 00a: cells per supertype and region ============================
    ws00a = wb.create_sheet("00a_cell_count_matrix")
    _ct_regions = NON_NEOCORTEX_REGIONS + ["Pooled"]
    _ct_hdrs    = ["Supertype_short", "Supertype_full"] + _ct_regions + \
                  ["Screened_independently_as"]
    _ct_widths  = [14, 28] + [12]*len(_ct_regions) + [40]
    _write_header(ws00a, _ct_hdrs, _ct_widths)

    # Note row
    _note_row = 2
    ws00a.cell(_note_row, 1,
               "Each region column is a separate screening run.").font = \
        Font(name="Calibri", size=8, italic=True, color="666666")
    ws00a.merge_cells(start_row=_note_row, start_column=1,
                      end_row=_note_row, end_column=len(_ct_hdrs))

    for ri_ct, (st_full, short) in enumerate(
            zip(CANDIDATE_SUPERTYPES, [SHORT_ST[s] for s in CANDIDATE_SUPERTYPES]), 3):
        row_vals = [short, st_full]
        row_total = 0
        pop_names_list = []
        for reg in NON_NEOCORTEX_REGIONS:
            n_ct = _tally_data.get((st_full, reg), 0)
            row_vals.append(n_ct)
            row_total += n_ct
            pop_names_list.append(f"{short}_{reg}")
        row_vals.append(row_total)
        row_vals.append(" | ".join(pop_names_list))
        for ci_ct, val_ct in enumerate(row_vals, 1):
            c_ct = _safe_cell(ws00a, ri_ct, ci_ct, val_ct)
            c_ct.alignment = _lft if ci_ct in (1, 2, len(_ct_hdrs)) else _ctr
            if ri_ct % 2 == 0 and ci_ct <= len(_ct_regions) + 2:
                c_ct.fill = _af
            # Highlight Pooled total column
            if ci_ct == len(_ct_regions) + 2:
                c_ct.font = Font(name="Calibri", size=9, bold=True)

    # Column-total row (all supertypes combined)
    ri_ct_tot = 3 + len(CANDIDATE_SUPERTYPES)
    ws00a.cell(ri_ct_tot, 1, "ALL supertypes").font = \
        Font(name="Calibri", size=9, bold=True)
    ws00a.cell(ri_ct_tot, 2, "(all 3 combined)").font = \
        Font(name="Calibri", size=9, italic=True, color="555555")
    col_sums = [sum(_tally_data.get((st, r), 0) for st in CANDIDATE_SUPERTYPES)
                for r in NON_NEOCORTEX_REGIONS]
    grand_total = sum(col_sums)
    for ci_tot, val_tot in enumerate(col_sums + [grand_total], 3):
        c_tot = _safe_cell(ws00a, ri_ct_tot, ci_tot, val_tot)
        c_tot.alignment = _ctr
        c_tot.font = Font(name="Calibri", size=9, bold=True)
        c_tot.fill = PatternFill("solid", fgColor="D9E1F2")

    # M4 total check row
    _pooled_col = 2 + len(NON_NEOCORTEX_REGIONS) + 1
    ri_check = ri_ct_tot + 1
    ws00a.cell(ri_check, 1, "M4_TOP_MASK total (Part 2 joblib)").font = \
        Font(name="Calibri", size=8, italic=True, color="444444")
    _m4_cell = ws00a.cell(ri_check, _pooled_col)
    _m4_cell.value = N_M4_TOP
    _m4_cell.font  = Font(name="Calibri", size=8, bold=True, color="444444")
    ws00a.cell(ri_check + 1, 1,
               f"{N_M4_TOP - grand_total:,} M4-top cells in other "
               f"supertypes.").font = \
        Font(name="Calibri", size=8, italic=True, color="888888")

    # == Sheet 00b: population catalogue ======================================
    ws00 = wb.create_sheet("00b_population_catalogue")
    _write_header(ws00,
                  ["Population", "Supertype_full", "Supertype_short", "Region",
                   "M4_cells_raw", "M4_cells_after_donor_filter",
                   "N_eligible_donors",
                   f"Min_pop_cells (={MIN_POP_CELLS})",
                   f"Min_donor_cells (={MIN_DONOR_POP_CELLS})",
                   "Status", "Independence"],
                  [22, 28, 14, 10, 14, 26, 18, 20, 22, 12, 38])
    for ri00, row_tally in _tally_df.iterrows():
        pop_name = row_tally["population"]
        st_full  = row_tally["supertype_full"]
        reg      = row_tally["region"]
        status_c = row_tally["status"]
        indep_note = (f"Only {reg} cells, only {SHORT_ST.get(st_full, '')} "
                      f"supertype.")
        vals00 = [
            pop_name, st_full, row_tally["supertype_short"], reg,
            row_tally["n_M4_cells_raw"],
            row_tally["n_cells_after_donor_filter"],
            row_tally["n_eligible_donors"],
            MIN_POP_CELLS, MIN_DONOR_POP_CELLS,
            status_c, indep_note,
        ]
        for ci00, val00 in enumerate(vals00, 1):
            c00 = _safe_cell(ws00, ri00 + 2, ci00, val00)
            c00.alignment = _lft if ci00 in (1, 2, 3, 10, 11) else _ctr
            if (ri00 + 2) % 2 == 0: c00.fill = _af
            if status_c == "OK" and ci00 == 10:
                c00.fill = PatternFill("solid", fgColor="C6EFCE")
            elif status_c == "SKIPPED" and ci00 == 10:
                c00.fill = PatternFill("solid", fgColor="FFC7CE")

    # == Collect the MasterDFG data ===========================================
    all_pop_data = {}   # pop_name -> DataFrame
    for pop_name in ALL_POPULATIONS:
        if not POPULATION_VALID.get(pop_name, False):
            all_pop_data[pop_name] = pd.DataFrame()
            continue
        out_dir_p = os.path.join(BASE_OUT, pop_name)
        all_rows_p = []
        for p in range(1, N_PARTS + 1):
            pkl = os.path.join(out_dir_p, f"Part{p:02d}_checkpoint.pkl")
            if os.path.exists(pkl):
                with open(pkl, "rb") as f:
                    all_rows_p.extend(pickle.load(f)["master_rows"])
        if all_rows_p:
            all_pop_data[pop_name] = pd.DataFrame(all_rows_p)
        else:
            all_pop_data[pop_name] = pd.DataFrame()

    # == Sheet 01: DFG counts per population ==================================
    ws01 = wb.create_sheet("01_dfg_counts")
    _write_header(ws01,
                  ["Population", "Supertype", "Region",
                   "GOIs_screened", "GOIs_with_DFG",
                   "Max_DFG_Total", "Top_GOI_by_DFG"],
                  [22, 26, 10, 14, 14, 14, 20])
    for ri01, pop_name in enumerate(ALL_POPULATIONS, 2):
        st_full, reg = POPULATION_MAP[pop_name]
        reg_disp = "ALL_pooled" if reg == _POOLED_SENTINEL else reg
        df_p = all_pop_data.get(pop_name, pd.DataFrame())
        if df_p.empty:
            vals01 = [pop_name, st_full, reg_disp, 0, 0, 0, "SKIPPED/NOT_RUN"]
        else:
            n_screened01 = len(df_p)
            n_hits01     = int((df_p["DFG_Total"] >= 1).sum())
            max_dfg01    = int(df_p["DFG_Total"].max()) if n_hits01 > 0 else 0
            top_goi01    = (df_p.sort_values("DFG_Total", ascending=False)
                            ["GOI"].iloc[0] if n_hits01 > 0 else " - ")
            vals01 = [pop_name, st_full, reg_disp, n_screened01,
                      n_hits01, max_dfg01, top_goi01]
        for ci01, val01 in enumerate(vals01, 1):
            c01 = _safe_cell(ws01, ri01, ci01, val01)
            c01.alignment = _lft if ci01 in (1, 2, 7) else _ctr
            if ri01 % 2 == 0: c01.fill = _af
            if ci01 == 5 and isinstance(val01, int) and val01 > 0:
                c01.fill = PatternFill("solid", fgColor="C6EFCE")

    # == Sheet 02: top 50 GOIs over all populations ===========================
    ws02 = wb.create_sheet("02_top_GOIs_all_populations")
    _write_header(ws02,
                  ["GOI", "Biotype", "Population", "Supertype", "Region",
                   "DFG_Total", "max_Spearman", "max_Kendall",
                   "Pct_expressing", "Mean_all"],
                  [16, 16, 22, 26, 10, 10, 14, 14, 14, 12])
    combined_rows = []
    for pop_name in ALL_POPULATIONS:
        df_p = all_pop_data.get(pop_name, pd.DataFrame())
        if df_p.empty: continue
        st_full, reg = POPULATION_MAP[pop_name]
        df_hits = df_p[df_p["DFG_Total"] >= 1].copy()
        if df_hits.empty: continue
        df_hits["Population"] = pop_name
        df_hits["Supertype"]  = st_full
        df_hits["Region"]     = "ALL_pooled" if reg == _POOLED_SENTINEL else reg
        df_hits["Biotype"]    = df_hits["GOI"].map(BIOTYPE_MAP).fillna("unknown")
        combined_rows.append(df_hits)

    if combined_rows:
        combined_df = (pd.concat(combined_rows, ignore_index=True)
                       .sort_values(["DFG_Total", "max_Spearman"],
                                     ascending=[False, False])
                       .head(50))
        cols02 = ["GOI", "Biotype", "Population", "Supertype", "Region",
                  "DFG_Total", "max_Spearman", "max_Kendall",
                  "Pct_expressing", "Mean_all"]
        for ri02, (_, row02) in enumerate(combined_df.iterrows(), 2):
            for ci02, col02 in enumerate(cols02, 1):
                c02 = _safe_cell(ws02, ri02, ci02,
                                 row02.get(col02, None))
                c02.alignment = _lft if ci02 in (1, 2, 3, 4) else _ctr
                if ri02 % 2 == 0: c02.fill = _af

    # == Per-population sheets (top 100) ======================================
    for pop_name in ALL_POPULATIONS:
        df_p = all_pop_data.get(pop_name, pd.DataFrame())
        ws_p = wb.create_sheet(pop_name[:31])
        _write_header(ws_p,
                      ["Rank", "GOI", "Biotype", "DFG_Total",
                       "max_Spearman", "max_Kendall",
                       "Pct_expressing", "Mean_all", "Median_all",
                       "Q1_all", "Q3_all", "IQR_all", "Mean_pos", "Median_pos"],
                      [8, 16, 16, 10, 14, 14, 14, 12, 12, 10, 10, 10, 12, 12])
        if df_p.empty:
            ws_p.cell(2, 1, "SKIPPED or no data").font = Font(
                name="Calibri", size=9, color="999999", italic=True)
            continue
        top100 = (df_p[df_p["DFG_Total"] >= 1]
                  .sort_values(["DFG_Total","max_Spearman"],
                                ascending=[False, False])
                  .head(100)
                  .reset_index(drop=True))
        if top100.empty:
            ws_p.cell(2, 1, "No DFG hits").font = Font(
                name="Calibri", size=9, color="999999", italic=True)
            continue
        top100.insert(0, "Rank", top100.index + 1)
        top100["Biotype"] = top100["GOI"].map(BIOTYPE_MAP).fillna("unknown")
        cols_p = ["Rank","GOI","Biotype","DFG_Total",
                  "max_Spearman","max_Kendall",
                  "Pct_expressing","Mean_all","Median_all",
                  "Q1_all","Q3_all","IQR_all","Mean_pos","Median_pos"]
        for ri_p, (_, row_p) in enumerate(top100.iterrows(), 2):
            for ci_p, col_p in enumerate(cols_p, 1):
                c_p = _safe_cell(ws_p, ri_p, ci_p, row_p.get(col_p, None))
                c_p.alignment = _lft if ci_p in (2, 3) else _ctr
                if ri_p % 2 == 0: c_p.fill = _af
                if col_p == "DFG_Total" and isinstance(
                        row_p.get(col_p), (int, np.integer)) and row_p[col_p] >= 5:
                    c_p.fill = PatternFill("solid", fgColor="C6EFCE")
                    c_p.font = Font(name="Calibri", size=9, bold=True)

    wb.save(fp_summary)
    print(f"\n  Cross-population summary saved: {fp_summary}")
    print(f"  Sheets: {[ws.title for ws in wb.worksheets]}")


# =============================================================================
# SECTION 15: Entry point
# =============================================================================

os.makedirs(BASE_OUT, exist_ok=True)

print("\n" + "=" * 70)
print("SECTION 15: Running all valid populations")
print("=" * 70)
print(f"  Populations to run: "
      f"{[p for p in ALL_POPULATIONS if POPULATION_VALID.get(p, False)]}")
print(f"  Populations to skip: "
      f"{[p for p in ALL_POPULATIONS if not POPULATION_VALID.get(p, False)]}")
print()

for _pop in ALL_POPULATIONS:
    if POPULATION_VALID.get(_pop, False):
        run_population(_pop)
    else:
        print(f"  Skipping {_pop} (insufficient cells/donors).")

print("\n" + "=" * 70)
print("SECTION 15: Writing cross-population summary Excel ...")
print("=" * 70)
write_cross_population_summary()

print("\n" + "=" * 70)
print("Part 6  -  ALL POPULATIONS COMPLETE")
print(f"Output: {BASE_OUT}/")
print("=" * 70)

print(f"\nDirectory listing ({BASE_OUT}/):")
for _fn in sorted(os.listdir(BASE_OUT)):
    _fp = os.path.join(BASE_OUT, _fn)
    if os.path.isdir(_fp):
        _n_files = len([f for f in os.listdir(_fp)
                        if not f.startswith(".")])
        print(f"  {_fn}/   ({_n_files} files)")
    elif not _fn.startswith("."):
        _kb = os.path.getsize(_fp) / 1024
        print(f"  {_fn:<60} {_kb:>8.1f} KB")
