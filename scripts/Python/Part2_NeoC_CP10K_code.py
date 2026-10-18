# =============================================================================
# Part2_NeoC_CP10K_code.py
# SEA-AD CEP-IP/CMRF, Part 2, M4 selection, neocortex (Figure S1A, Table S1)
#   Scores the nuclei of the neocortical regions (DFC, MTG, AnG, STG, V1C, ITG)
#   with the Mic.13 signature using three gene weightings (M1 fold change, M2
#   fold change x -log10 padj, M3 padj rank) and their intersections (M4 = M1
#   and M2, M5, M6). Checks the cutoff (Sections 2b to 2g), then for each cutoff
#   of a sweep writes the M4-top mask, supertype enrichment, ADNC and clinical
#   correlations, figures and the namespace that later parts load. M4_TOP_MASK
#   (top x% of M1 and top x% of M2) is the mask all later parts use. The
#   primary cutoff (top 3.0%) runs last.
#
# Usage (from the SEA-AD folder, or set SEAAD_ROOT), after Part 1
#   python3 Part2_NeoC_CP10K_code.py
#
# Outputs ([1]_NeoC/Part2_NeoC_M4/ and its top_<pct>pct/ subfolders): cutoff
#   checks (S2b to S2g), scores, enrichment and clinical tables (CSV), figures
#   (PNG), S5b_donor_clinical.csv, Part2_neocortex_env.joblib, sweep summary.
# =============================================================================

# == Neocortical regions ====================================================
NEOCORTEX_REGIONS = ["DFC", "MTG", "AnG", "STG", "V1C", "ITG"]

# %% SECTION 0: Imports and config
import os, gc, time, warnings, zlib
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import spearmanr, mannwhitneyu, kruskal
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use("Agg")  # non-interactive backend
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as _mcolors
import scanpy as sc

np.random.seed(42)
import random
random.seed(42)
sc.settings.n_jobs    = 1
os.environ["PYTHONHASHSEED"] = "42"   # affects subprocesses only
sc.settings.verbosity = 1

# == Project folder and input file ============================================
BASE_DIR = os.environ.get("SEAAD_ROOT", os.getcwd())
os.chdir(BASE_DIR)

FULL_H5AD = ("SEAAD_Immune_HIP_MEC_LEC_ITG_MTG_FI_STG_DFC_AnG_V1C_"
             "RNAseq_final-nuclei.2026-06-22.h5ad")
# == Output folder: <arm folder>/Part2_<arm>_M4/ ================================
ARM_FOLDER  = r"[1]_NeoC"     # the arm's top-level folder under BASE_DIR
ARM_TAG     = "NeoC"             # short arm label used in paths and tables
PART_TAG    = "Part2"             # part label used in folder names
MASK_TAG    = "M4"                # the exported mask, kept in the folder name
ARM_ROOT    = os.path.join(BASE_DIR, ARM_FOLDER)
ARM_DIR     = os.path.join(ARM_ROOT, f"{PART_TAG}_{ARM_TAG}_{MASK_TAG}")
OUT_DIR     = ARM_DIR             # rebound per cutoff by the sweep driver
os.makedirs(ARM_DIR, exist_ok=True)
print(f"[OUT] arm folder   -> {ARM_ROOT}")
print(f"[OUT] this part    -> {ARM_DIR}")
print(f"[OUT] per-cutoff   -> {os.path.join(ARM_DIR, 'top_<pct>pct')}")

# == Expected normalisation of X: ln(CP10K + 1), checked in Section 2 =========
NORM_TARGET_SUM = 1e4
NORM_CEILING    = float(np.log1p(NORM_TARGET_SUM))   # 9.210340
NORM_TOL        = 1e-6
NORM_UNS_KEY    = "X_normalization"
NORM_UNS_EXPECT = "ln(UP10K+1)"
NORM_CHECK_N    = 3000        # rows sampled for the numeric reconstruction test
NORM_ABORT_ON_FAIL = True     # stop if X fails the check

# == Expected dimensions, from Parts 1b and 1e ================================
EXPECT_SHAPE     = (240651, 36601)
EXPECT_ARM_CELLS = 134088
EXPECT_REGION_N  = {"DFC": 41553, "MTG": 36686, "AnG": 10279,
                   "STG": 15214, "V1C": 14245, "ITG": 16111}
EXPECT_N_DONORS  = 84
EXPECT_N_SUPERTYPES = 8       # distinct immune supertypes in this object

# == Nuclei flagged by upstream pass/fail fields ==============================
# Kept: SEA-AD marks them 'Used in analysis' after its own QC.
KEEP_FLAGGED_NUCLEI = True
FLAG_FAIL_TOKENS    = {"fail", "failed", "true", "1", "yes", "y"}

# == Region aliases (PFC, DLPFC and A9 are DFC) ================================
REGION_ALIASES = {"PFC": "DFC", "DLPFC": "DFC", "A9": "DFC"}

# == SEA-AD supertype renames (nomenclature 2023-05-05 to 2024-02-13) =========
KNOWN_RENAMES = {
    "VLMC_2":              "Pericyte_1",
    "VLMC_2_1-SEAAD":      "SMC-SEAAD",
    "VLMC_2_2-SEAAD":      "Pericyte_2-SEAAD",
    "Micro-PVM_2_2-SEAAD": "Lymphocyte",
    "Micro-PVM_1_1-SEAAD": "Monocyte",
}
TARGET_STYPE = "Micro-PVM_3-SEAAD"   # asserted to exist before any scoring

# == obs column names (defaults; resolved against the file in Section 1) ======
adnc_col   = "Overall AD neuropathological Change"
braak_col  = "Braak"
donor_col  = "Donor ID"
stype_col  = "Supertype"
region_col = "Brain Region"
cerad_col  = "CERAD score"
thal_col   = "Thal"
cog_col    = "Cognitive Status"
apoe_col   = "APOE Genotype"
mmse_col   = "Last MMSE Score"
age_col    = "Age at Death"
caa_col    = "Overall CAA Score"

# Tissue-quality covariates for donor_clinical (used by Part 10h).
pmi_col    = "PMI"
rin_col    = "RIN"
ph_col     = "Brain pH"
weight_col = "Fresh Brain Weight"
severe_col = "Severely Affected Donor"

# A missing critical column stops the run; a missing optional one drops its test.
CRITICAL_COLS = ["adnc_col", "braak_col", "donor_col", "stype_col", "region_col"]

# ADNC order and palette are set in Section 1 from the column values.
_ADNC_CANONICAL = ["Not AD", "Low", "Intermediate", "High"]
_ADNC_PAL_CANONICAL = {
    "Not AD":       "#2166AC",
    "Low":          "#74ADD1",
    "Intermediate": "#F46D43",
    "High":         "#A50026",
}
# Placeholders, set in Section 1
adnc_order    = _ADNC_CANONICAL[:]
adnc_pal      = _ADNC_PAL_CANONICAL.copy()
adnc_rank_map = {a: i for i, a in enumerate(adnc_order)}

# == Clinical ordinal maps (levels checked against this file in Section 1) ====
BRAAK_MAP = {"Braak 0": 0, "Braak I": 1, "Braak II": 2, "Braak III": 3,
             "Braak IV": 4, "Braak V": 5, "Braak VI": 6}
CERAD_MAP = {"Absent": 0, "Sparse": 1, "Moderate": 2, "Frequent": 3}
THAL_MAP  = {"Thal 0": 0, "Thal 1": 1, "Thal 2": 2,
             "Thal 3": 3, "Thal 4": 4, "Thal 5": 5}
CAA_MAP   = {"Not identified": 0, "Mild": 1, "Moderate": 2, "Severe": 3}

ORDINAL_UNMAPPED_ABORT = 0.05   # abort above 5% unmapped on any ordinal

PADJ_FLOOR = 1e-320   # replace exact 0 p_val_adj (machine underflow in original table)

# == Cutoff sweep ============================================================
# Sections 3-11 run once per cutoff into top_<pct>pct/; the primary runs last.
PRIMARY_TOP_PCT      = 0.030     # 3.0% in both arms
SWEEP_PCTS_REQUESTED = [0.025, 0.029, 0.030, 0.031, 0.050]
SWEEP_PCTS = sorted(set(SWEEP_PCTS_REQUESTED) | {PRIMARY_TOP_PCT})
SWEEP_PCTS = ([p for p in SWEEP_PCTS if abs(p - PRIMARY_TOP_PCT) > 1e-12]
              + [PRIMARY_TOP_PCT])
SWEEP_SUBDIR_FMT = "top_{:.1f}pct"
TOP_PCT = PRIMARY_TOP_PCT            # rebound by the sweep driver each pass

# Figures: which sweep passes make and show them.
SWEEP_MAKE_FIGS = "all"            # "all" | "primary_only" | "none"
SWEEP_SHOW_FIGS = "primary_only"   # "all" | "primary_only" | "none"

# == Cutoff-selection rule (Section 2b) =======================================
# Other supertypes must sit at or below 1% of their own cells in the top set.
RULE_ENRICHED_STYPES = ["Micro-PVM_3-SEAAD", "Micro-PVM_4-SEAAD"]
RULE_MAX_PCT = 1.0     # ceiling, in % of the supertype's own cells
RULE_MIN_N   = 100     # supertypes smaller than this cannot express the rule
SCAN_GRID    = np.round(np.arange(0.010, 0.0601, 0.0005), 5)
SCAN_METHOD  = "M4_M1andM2"   # the mask that is actually exported downstream

# == Cutoff diagnostics (Sections 2d to 2g) ====================================
# 2d mixture model, base scorers only (M4 to M6 are intersections)
MIX_METHODS    = ["M1_FC", "M2_Hybrid", "M3_pRank"]
MIX_SUBSAMPLE  = 50_000   # BIC scales with n, so the null must match n exactly
MIX_NULL_B     = 200      # skew-normal draws for the unimodal null
MIX_N_INIT     = 10       # EM restarts for the observed fits
MIX_SEP_MIN    = 2.0      # component separation below this is not a boundary
MIX_SEED       = 42

# 2e donor bootstrap
BOOT_METHODS   = ["M1_FC", "M2_Hybrid", "M4_M1andM2"]
BOOT_GRID      = np.round(np.arange(0.010, 0.0601, 0.0025), 5)   # 21 points
BOOT_B         = 200      # donor resamples
BOOT_CORE_FRAC = 0.95     # "stable core" = selected in >= this fraction
BOOT_FLAT_TOL  = 0.02     # Jaccard range below this counts as flat
BOOT_SEED      = 42

# 2f arm matching; each arm writes its half here and the second run joins them.
MATCH_SHARED_DIR       = os.path.join(BASE_DIR, "_shared_cutoff")
OTHER_ARM_PRIMARY_PCT  = 3.0   # the other arm's primary cutoff (%)

# 2g intersection redundancy
OVERLAP_REDUNDANT_AT   = 0.95

ARM_COLOR = "#7C6BAF"   # arm accent, matching the Part1e figure palette

# == Marker-dilution curve (Section 2c) =======================================
# heldout_up: Mic.13 markers not in the scoring panel (the informative curve);
# inpanel_up: in the panel (reference only); homeostatic: should rise.
MIC13_MARKER_CANDIDATES = {
    "heldout_up":  ["TREM2", "ADAM10", "SMAD3"],
    "inpanel_up":  ["PTPRG", "TPRG1", "SPP1", "GPNMB", "APOE"],
    "homeostatic": ["P2RY12", "TMEM119", "SALL1"],
}
# Marginal band = top-x% minus top-(x - MARKER_BAND_PP)%.
MARKER_BAND_PP      = 0.25  # width of the marginal band, percentage points
MARKER_MARGINAL_MIN = 1.5   # marginal fold floor
MARKER_MIN_PANEL    = 2     # minimum usable held-out genes


# == Helpers ===================================================================
def find_col(cols, keywords):
    """Resolve one obs column name: exact match first, then the shortest
    substring match. Returns None when nothing matches."""
    cols = list(cols)
    for kw in keywords:
        for c in cols:
            if c == kw:
                return c
    for kw in keywords:
        hits = [c for c in cols if str(kw).lower() in str(c).lower()]
        if hits:
            return sorted(hits, key=len)[0]
    return None


def to_num(s, mapping=None):
    """Map through `mapping` (if given) and coerce to float (also for Categoricals)."""
    ser = pd.Series(s)
    if mapping is not None:
        ser = ser.astype(str).str.strip().map(mapping)
    return pd.to_numeric(ser, errors="coerce")


def audit_ordinal(name, series, mapping):
    """Print the levels of one ordinal; stop if too many are unmapped."""
    raw = pd.Series(series).astype(str).str.strip()
    num = to_num(raw, mapping)
    n   = len(raw)
    n_unmapped = int(num.isna().sum())
    frac = n_unmapped / max(n, 1)
    lev  = raw.value_counts()
    print(f"    {name:<7} n={n:>5,}  unmapped={n_unmapped:>4,} ({100*frac:5.1f}%)"
          f"  dtype_out={num.dtype}")
    print(f"             levels: "
          + ", ".join(f"{k} ({v})" for k, v in lev.items()))
    missing = sorted(set(raw[num.isna()]))
    if missing:
        print(f"             [!] UNMAPPED LEVELS: {missing}")
    if frac > ORDINAL_UNMAPPED_ABORT:
        raise SystemExit(
            f"ABORT - {name} is {100*frac:.1f}% unmapped (limit "
            f"{100*ORDINAL_UNMAPPED_ABORT:.0f}%). Fix the map before scoring.")
    return num


def stable_seed(text):
    """Deterministic per-string seed (crc32)."""
    return zlib.crc32(str(text).encode("utf-8")) % (2**31)


def _maybe_show():
    """Show the figure only when the current sweep pass allows it."""
    if globals().get("_show_figs", True):
        plt.show()


def get_cmap_compat(name):
    """Colormap lookup for old and new matplotlib."""
    try:
        return matplotlib.colormaps[name]
    except Exception:
        return plt.get_cmap(name)


MIC13_UP = {
    "PTPRG":   (2.543639, 0.0),        "SPP1":    (2.083305, 0.0),
    "GLDN":    (1.719563, 0.0),        "TPRG1":   (1.555551, 0.0),
    "MYO1E":   (1.433911, 0.0),        "DPYD":    (1.015599, 3.4452e-246),
    "ALCAM":   (0.803272, 7.6694e-216),"PPARG":   (0.756824, 8.5687e-222),
    "DIRC3":   (0.707277, 1.73442e-85),"CPM":     (0.703485, 3.25e-139),
    "ASAP1":   (0.675848, 1.183e-259), "ARID5B":  (0.664147, 1.6153e-215),
    "MITF":    (0.660149, 1.3225e-179),"FMNL2":   (0.632292, 7.618e-203),
    "SOCS6":   (0.596924, 1.3174e-147),"DSCAM":   (0.583210, 1.0136e-104),
    "MSR1":    (0.567882, 8.01e-125),  "XYLT1":   (0.550444, 9.1264e-229),
    "CADM1":   (0.547246, 3.4144e-228),"CTSD":    (0.529815, 2.2373e-190),
    "GPNMB":   (0.519895, 8.3428e-132),"OLR1":    (0.519176, 2.4084e-131),
    "SCIN":    (0.513077, 1.9023e-103),"SAMD4A":  (0.503272, 1.3275e-109),
    "HDAC9":   (0.496561, 4.2897e-107),"PADI2":   (0.488738, 1.34447e-74),
    "LPL":     (0.486333, 2.752e-280), "GAS7":    (0.480107, 4.609e-157),
    "ADAMTS17":(0.475107, 1.18991e-70),"RTTN":    (0.474131, 1.78636e-65),
    "LIMS1":   (0.464410, 7.0679e-119),"APOE":    (0.464222, 1.35365e-95),
    "FLT1":    (0.461284, 2.36273e-98),"PLXDC2":  (0.459186, 4.3499e-264),
    "ITGAX":   (0.453712, 4.6365e-155),"ATG7":    (0.448784, 2.4888e-101),
    "CACNA1A": (0.436573, 2.32297e-86),"BCL2":    (0.430362, 7.56819e-77),
    "FNIP2":   (0.430191, 9.83893e-97),"STARD13": (0.428650, 2.12511e-45),
    "CHSY3":   (0.421250, 9.99458e-69),"NPL":     (0.418654, 1.07056e-97),
    "RASGRP3": (0.412015, 6.5959e-161),"ZNF804A": (0.406739, 2.26996e-51),
    "LNCAROD": (0.402851, 2.08297e-63),"SDK1":    (0.401401, 1.01717e-52),
    "CD83":    (0.385531, 2.97156e-57),"RGCC":    (0.379567, 2.2541e-274),
    "ELL2":    (0.375538, 6.37489e-64),"MTSS1":   (0.373229, 1.3745e-56),
    "ADARB1":  (0.364487, 2.7213e-195),"ADK":     (0.363451, 6.71523e-80),
    "PRKCE":   (0.360454, 3.35701e-75),"FRMD4B":  (0.354601, 3.03763e-60),
    "SLC11A1": (0.353130, 2.03157e-63),"DTNA":    (0.352197, 1.21504e-56),
    "FNDC3B":  (0.338784, 3.32302e-60),"LIPA":    (0.333632, 7.7832e-103),
    "MALAT1":  (0.331582, 6.8469e-174),"CDK6":    (0.330770, 4.76582e-52),
    "WIPF3":   (0.326277, 1.6209e-135),"CYTH3":   (0.322623, 1.45987e-83),
    "IPCEF1":  (0.317401, 2.75625e-44),"ZNF618":  (0.315610, 3.7878e-82),
    "PDE3B":   (0.306725, 8.35622e-37),"DENND4C": (0.299626, 2.83266e-67),
    "RALGAPA2":(0.298356, 4.62747e-52),"NCK2":    (0.292940, 7.46584e-39),
    "SYNDIG1": (0.292679, 7.96235e-26),"NHSL1":   (0.288719, 7.34929e-17),
    "SAMSN1":  (0.286930, 7.11031e-42),"PTPRM":   (0.285253, 1.10507e-40),
    "MAFB":    (0.284475, 2.48643e-73),"TANC2":   (0.283339, 8.63317e-37),
    "LPP":     (0.281859, 6.42334e-56),"BICD1":   (0.275823, 9.88825e-56),
    "TGFBR1":  (0.272279, 1.86421e-28),"RASGEF1B":(0.267955, 1.71847e-75),
    "CLEC5A":  (0.266228, 9.35161e-67),"MAML2":   (0.263986, 5.11102e-40),
    "ASAH1":   (0.263746, 6.02525e-44),"PDE8A":   (0.256056, 2.10241e-40),
    "SLC1A3":  (0.252802, 1.47788e-39),"TMEM117": (0.251988, 1.49719e-29),
    "TRHDE":   (0.251163, 3.92582e-96),
}
MIC13_DOWN = {
    "PICALM":     (0.250272, 2.57161e-42), "MAF":       (0.252086, 5.15499e-31),
    "RASGEF1C":   (0.256926, 6.49307e-26), "INPP5D":    (0.257635, 4.18499e-51),
    "NAV3":       (0.264112, 2.13387e-26), "TBXAS1":    (0.267437, 1.0915e-57),
    "RPS27":      (0.267872, 1.29204e-18), "FCHSD2":    (0.268911, 1.94018e-36),
    "SRGAP2B":    (0.278794, 1.50584e-52), "PREX1":     (0.281410, 1.72959e-56),
    "MIR646HG":   (0.287985, 7.79005e-33), "HSPA1A":    (0.288355, 9.46674e-23),
    "SERPINB9":   (0.293819, 5.65094e-45), "HAMP":      (0.294443, 2.57688e-48),
    "MS4A6A":     (0.302231, 4.05155e-45), "SRGAP2":    (0.302969, 1.09111e-77),
    "MIR99AHG":   (0.307810, 9.54668e-49), "SRGN":      (0.308019, 5.1395e-14),
    "ELMO1":      (0.308191, 2.17401e-84), "NAMPT":     (0.309445, 1.00e-7),
    "EPB41L3":    (0.314441, 7.97693e-35), "CX3CR1":    (0.316787, 5.02488e-43),
    "FYB1":       (0.321048, 1.02187e-50), "CTTNBP2":   (0.321336, 4.09724e-53),
    "FAM13A":     (0.329585, 2.13396e-57), "FRMD4A":    (0.330631, 8.81325e-79),
    "TLN2":       (0.330689, 1.16691e-56), "SSBP2":     (0.343290, 1.55589e-53),
    "ADAM28":     (0.367841, 3.43841e-80), "IL6ST":     (0.374294, 2.81998e-75),
    "ABCC4":      (0.374460, 1.59647e-71), "MVB12B":    (0.381065, 2.29474e-93),
    "FAM177B":    (0.386167, 4.79934e-74), "CPED1":     (0.392684, 1.03314e-74),
    "CSGALNACT1": (0.397928, 1.48217e-65), "GRID2":     (0.408088, 1.18795e-11),
    "SLC9A9":     (0.417432, 9.8369e-119), "IGSF21":    (0.418851, 1.7458e-91),
    "MRC1":       (0.437510, 2.30913e-89), "ST6GALNAC3":(0.438452, 1.36984e-63),
    "AOAH":       (0.470385, 5.0055e-113), "CCDC26":    (0.479074, 4.45678e-32),
    "PRKCA":      (0.527630, 2.1296e-109), "SRGAP1":    (0.593008, 3.7886e-160),
    "F13A1":      (0.596145, 5.80837e-79), "ITPR2":     (0.610236, 5.3684e-255),
    "FOXP2":      (0.667381, 3.311e-174),  "TMEM163":   (0.832115, 4.8286e-131),
}

# FC-only dicts
MIC13_UP_FC   = {g: v[0] for g, v in MIC13_UP.items()}
MIC13_DOWN_FC = {g: v[0] for g, v in MIC13_DOWN.items()}

# == Gene weights for M1, M2 and M3 ===========================================
def _build_weights(gene_dict, method="M1"):
    """
    gene_dict : {gene: (fc, padj)}  FC is positive magnitude for both UP and DOWN.
    method    : "M1" = pure FC weight
                "M2" = FC x -log10(padj)  [hybrid: effect size x confidence]
                "M3" = pure p_adj rank     [NO FC; rank alone drives weight]
    Returns   : {gene: weight} normalised to sum=1
    """
    genes  = list(gene_dict.keys())
    fc_arr = np.array([gene_dict[g][0] for g in genes], dtype=np.float64)
    pa_arr = np.array([max(gene_dict[g][1], PADJ_FLOOR) for g in genes],
                      dtype=np.float64)

    if method == "M1":
        w = fc_arr.copy()

    elif method == "M2":
        w = fc_arr * (-np.log10(pa_arr))

    elif method == "M3":
        order    = np.argsort(pa_arr)
        rank_arr = np.empty_like(order)
        rank_arr[order] = np.arange(len(order))
        rank_norm = 1.0 - rank_arr / (len(rank_arr) - 1) + 1e-6
        w = rank_norm

    else:
        raise ValueError(f"Unknown method: {method}")

    total = w.sum()
    w = w / total if total > 0 else np.ones(len(genes)) / len(genes)
    return dict(zip(genes, w))

# Build all six weight dicts (UP and DOWN for each of three methods)
SCORERS = {
    "M1_FC":      (_build_weights(MIC13_UP, "M1"), _build_weights(MIC13_DOWN, "M1")),
    "M2_Hybrid":  (_build_weights(MIC13_UP, "M2"), _build_weights(MIC13_DOWN, "M2")),
    "M3_pRank":   (_build_weights(MIC13_UP, "M3"), _build_weights(MIC13_DOWN, "M3")),
}
SCORER_LABELS = {
    "M1_FC":     "M1 — FC-weighted only",
    "M2_Hybrid": "M2 — Hybrid FC x -log10(padj)",
    "M3_pRank":  "M3 — Pure p_adj rank (no FC)",
    "M4_M1andM2": "M4 — Intersection M1 ∩ M2",
    "M5_M1andM3": "M5 — Intersection M1 ∩ M3",
    "M6_M2andM3": "M6 — Intersection M2 ∩ M3",
}
INTERSECTION_PAIRS = {
    "M4_M1andM2": ("M1_FC",    "M2_Hybrid"),
    "M5_M1andM3": ("M1_FC",    "M3_pRank"),
    "M6_M2andM3": ("M2_Hybrid","M3_pRank"),
}

# == Clinical tests ============================================================
CLINICAL_TESTS = [
    ("ADNC rank",    "adnc_rank",    "spearman",    "pos",  "Composite Abeta+tau"),
    ("Braak stage",  "braak_num",    "spearman",    "pos",  "Tau spread"),
    ("CERAD score",  "cerad_num",    "spearman",    "pos",  "Abeta plaque burden"),
    ("Thal phase",   "thal_num",     "spearman",    "pos",  "Abeta staging"),
    ("APOE e4 dose", "apoe_e4_dose", "spearman",    "pos",  "APOE e4 genotype"),
    ("MMSE score",   "mmse_num",     "spearman",    "neg",  "Cognitive score (neg ctrl)"),
    ("Age at death", "age_num",      "spearman",    "null", "Age (neg ctrl)"),
    ("CAA score",    "caa_num",      "spearman",    "null", "Vascular amyloid"),
]

print("=" * 74)
print("Part 2 - Mic.13 sweep + 6-way scoring comparison [NEOCORTICAL]")
print(f"  Release : 2026-06-22 Multiregion (CP10K, ln(UP10K+1))")
print(f"  File    : {FULL_H5AD}")
print(f"  Top-%   : {TOP_PCT*100:.1f}%")
print(f"  Regions : {NEOCORTEX_REGIONS}")
print(f"  Scorers : {list(SCORERS.keys())}")
print(f"  Flagged nuclei: {'KEEP' if KEEP_FLAGGED_NUCLEI else 'DROP'}")
print(f"  Output  : {os.path.abspath(OUT_DIR)}")
print("=" * 74)


# %% SECTION 1: Load metadata + neocortical guardrail
print("\nSECTION 1: Load metadata + neocortical guardrail")
t0 = time.time()

adata = sc.read_h5ad(FULL_H5AD, backed="r")
obs   = adata.obs.copy()
print(f"  Full dataset: {adata.shape}")

# == Dimensions check =========================================================
if tuple(adata.shape) != EXPECT_SHAPE:
    print(f"  [!] shape mismatch - expected {EXPECT_SHAPE}, got {tuple(adata.shape)}")
else:
    print(f"  [OK] shape matches Part1b: {EXPECT_SHAPE[0]:,} x {EXPECT_SHAPE[1]:,}")

# == Declared normalisation (uns) =============================================
_uns_norm = adata.uns.get(NORM_UNS_KEY, None) if adata.uns is not None else None
print(f"  uns['{NORM_UNS_KEY}'] = {_uns_norm!r}  (expected {NORM_UNS_EXPECT!r})")
if _uns_norm is None:
    print("  [!] no X_normalization key")
elif str(_uns_norm).strip() != NORM_UNS_EXPECT:
    _msg = (f"ABORT - uns['{NORM_UNS_KEY}'] is {_uns_norm!r}, not "
            f"{NORM_UNS_EXPECT!r}.")
    if NORM_ABORT_ON_FAIL:
        raise SystemExit(_msg)
    print("  [!] " + _msg)
else:
    print("  [OK] declared normalisation is CP10K, natural log.")
print(f"  uns['title'] = {adata.uns.get('title', '(absent)')!r}")
print(f"  uns['batch_condition'] = {adata.uns.get('batch_condition', '(absent)')}")

# == Resolve obs column names =================================================
_oc = list(obs.columns)
print(f"\n  obs has {len(_oc)} columns; resolving names ...")

_resolved = {
    "adnc_col":   find_col(_oc, [adnc_col, "AD neuropath", "ADNC"]),
    "braak_col":  find_col(_oc, [braak_col, "Braak stage", "Braak", "braak"]),
    "donor_col":  find_col(_oc, [donor_col, "Donor ID", "donor_id", "Donor"]),
    "stype_col":  find_col(_oc, [stype_col, "Supertype", "supertype"]),
    "region_col": find_col(_oc, [region_col, "Brain Region", "brain_region",
                                 "Region", "region", "tissue"]),
    "cerad_col":  find_col(_oc, [cerad_col, "CERAD score", "CERAD", "cerad"]),
    "thal_col":   find_col(_oc, [thal_col, "Thal phase", "Thal", "thal"]),
    "cog_col":    find_col(_oc, [cog_col, "Cognitive Status", "Cognitive"]),
    "apoe_col":   find_col(_oc, [apoe_col, "APOE Genotype", "APOE"]),
    "mmse_col":   find_col(_oc, [mmse_col, "Last MMSE Score", "MMSE"]),
    "age_col":    find_col(_oc, [age_col, "Age at Death", "age at death"]),
    "caa_col":    find_col(_oc, [caa_col, "Overall CAA Score", "CAA"]),
    "pmi_col":    find_col(_oc, [pmi_col, "PMI"]),
    "rin_col":    find_col(_oc, [rin_col, "RIN"]),
    "ph_col":     find_col(_oc, [ph_col, "Brain pH", "pH"]),
    "weight_col": find_col(_oc, [weight_col, "Fresh Brain Weight",
                                 "Brain Weight"]),
    "severe_col": find_col(_oc, [severe_col, "Severely Affected Donor",
                                 "Severely Affected"]),
}

_missing_critical = [k for k in CRITICAL_COLS if _resolved.get(k) is None]
print(f"\n  {'name':<12} {'default':<38} {'resolved':<38} status")
print("  " + "-" * 96)
for _k, _v in _resolved.items():
    _default = globals()[_k]
    if _v is None:
        _status = "CRITICAL-MISSING" if _k in CRITICAL_COLS else "optional-absent"
    elif _v == _default:
        _status = "OK"
    else:
        _status = "REMAPPED"
    print(f"  {_k:<12} {str(_default):<38} {str(_v):<38} {_status}")
    globals()[_k] = _v

if _missing_critical:
    print("\n  Full obs column list:")
    for _c in _oc:
        print(f"    {_c}")
    raise SystemExit(
        f"ABORT - critical obs columns unresolved: {_missing_critical}.")

# == Supertype renames ========================================================
_st_raw = set(obs[stype_col].astype(str).unique())
_to_rename = {k: v for k, v in KNOWN_RENAMES.items() if k in _st_raw}
if _to_rename:
    print(f"\n  applying documented supertype renames: {_to_rename}")
    obs[stype_col] = obs[stype_col].astype(str).replace(_to_rename)
else:
    print(f"\n  no legacy supertype labels present")

_st_now = sorted(obs[stype_col].astype(str).unique())
print(f"  Supertypes in file ({len(_st_now)}): {_st_now}")
if len(_st_now) != EXPECT_N_SUPERTYPES:
    print(f"  [!] expected {EXPECT_N_SUPERTYPES} immune supertypes "
          f"(Part1b), found {len(_st_now)}")
if TARGET_STYPE not in _st_now:
    raise SystemExit(
        f"ABORT - TARGET_STYPE {TARGET_STYPE!r} is absent. Present: {_st_now}")
print(f"  [OK] TARGET_STYPE {TARGET_STYPE!r} present "
      f"({int((obs[stype_col].astype(str)==TARGET_STYPE).sum()):,} nuclei, "
      f"all regions).")

# == Nuclei on failed library components ======================================
print(f"\n  Upstream pass/fail flags ("
      f"{'KEEP' if KEEP_FLAGGED_NUCLEI else 'DROP'})")
_flag_cols = [c for c in _oc
              if c in ("experiment_component_failed",
                       "rna_amplification_pass_fail",
                       "library_prep_pass_fail")]
_flag_mask_total = pd.Series(False, index=obs.index)
if not _flag_cols:
    print("    (none of the three pass/fail columns present)")
for _fc in _flag_cols:
    _vals = obs[_fc].astype(str).str.strip().str.lower()
    _bad  = _vals.isin(FLAG_FAIL_TOKENS)
    _flag_mask_total |= _bad
    print(f"    {_fc:<32}: {int(_bad.sum()):>6,} flagged "
          f"({100*_bad.mean():.2f}%)  levels={sorted(set(_vals))}")
print(f"    union flagged (all regions)     : {int(_flag_mask_total.sum()):>6,}")
if not KEEP_FLAGGED_NUCLEI:
    obs = obs[~_flag_mask_total].copy()
    print(f"    DROPPED - obs now {len(obs):,} rows.")
else:
    print("    RETAINED.")

# == Region aliases and arm filter ============================================
print(f"\n  {'-'*66}")
print("  NEOCORTICAL GUARDRAIL")
_regions_raw = sorted(obs[region_col].astype(str).unique())
print(f"  Region codes in file (raw): {_regions_raw}")

_alias_hits = {k: v for k, v in REGION_ALIASES.items() if k in _regions_raw}
if _alias_hits:
    _collide = [v for v in _alias_hits.values() if v in _regions_raw]
    if _collide:
        raise SystemExit(
            f"ABORT - region alias collision: {_alias_hits} but {_collide} "
            f"already present.")
    print(f"  remapping region aliases: {_alias_hits}")
    obs[region_col] = obs[region_col].astype(str).replace(_alias_hits)

all_regions  = sorted(obs[region_col].astype(str).unique())
present_nc   = [r for r in NEOCORTEX_REGIONS if r in all_regions]
excluded     = [r for r in all_regions if r not in NEOCORTEX_REGIONS]
absent_nc    = [r for r in NEOCORTEX_REGIONS if r not in all_regions]
print(f"  All regions after alias remap: {all_regions}")
if len(all_regions) != 10:
    print(f"  [!] expected 10 regions (Part1d), found {len(all_regions)}")
if absent_nc:
    print(f"  WARNING - expected neocortical regions not found: {absent_nc}")

n_before = len(obs)
nc_mask  = obs[region_col].astype(str).isin(present_nc)
obs      = obs[nc_mask].copy()
n_after  = len(obs)
nc_cell_indices = np.where(nc_mask.values)[0]

print(f"  Excluded (non-neocortical + other): {excluded}")
print(f"  Retained (neocortical):        observed   Part1e expected")
for r in present_nc:
    _obsn = int((obs[region_col] == r).sum())
    _expn = EXPECT_REGION_N.get(r)
    _mark = "" if (_expn is None or _obsn == _expn) else "   <- DIFFERS"
    print(f"    {r:<8}: {_obsn:>9,}   {_expn if _expn is not None else '-':>9}{_mark}")
print(f"  Total: {n_before:,} -> {n_after:,}  ({n_before-n_after:,} removed)")
if n_after != EXPECT_ARM_CELLS:
    print(f"  [!] arm size {n_after:,} != Part1e expected {EXPECT_ARM_CELLS:,}")
    if KEEP_FLAGGED_NUCLEI:
        print(f"      flagged nuclei were kept")
    else:
        print(f"      flagged nuclei were dropped")
else:
    print(f"  [OK] arm size matches Part1e: {EXPECT_ARM_CELLS:,}")
print(f"  [GUARDRAIL PASS]")
print(f"  {'-'*66}\n")

_n_donors_arm = int(obs[donor_col].astype(str).nunique())
print(f"  Donors contributing to this arm: {_n_donors_arm} "
      f"(cohort = {EXPECT_N_DONORS})")

stype_counts = obs[stype_col].value_counts()
print(f"\n  Supertype breakdown ({n_after:,} neocortical cells):")
for st, n in stype_counts.items():
    print(f"    {st:<40}: {n:>7,}")

has_umap = "X_umap" in adata.obsm
if has_umap:
    umap_coords = adata.obsm["X_umap"][nc_cell_indices]
    print(f"\n  X_umap present -> {umap_coords.shape}")
else:
    umap_coords = None
    print("\n  X_umap absent - UMAP figures will be skipped.")

# == Mic.13 genes present in the matrix =======================================
var_names = list(adata.var_names)
up_avail  = [g for g in MIC13_UP_FC   if g in var_names]
dn_avail  = [g for g in MIC13_DOWN_FC if g in var_names]
_up_miss  = sorted(set(MIC13_UP_FC)   - set(up_avail))
_dn_miss  = sorted(set(MIC13_DOWN_FC) - set(dn_avail))
print(f"\n  Mic.13 UP   genes in matrix: {len(up_avail)} / {len(MIC13_UP_FC)}")
if _up_miss:   print(f"    missing UP  : {_up_miss}")
print(f"  Mic.13 DOWN genes in matrix: {len(dn_avail)} / {len(MIC13_DOWN_FC)}")
if _dn_miss:   print(f"    missing DOWN: {_dn_miss}")
if not up_avail or not dn_avail:
    raise SystemExit("ABORT - one of the Mic.13 gene sets is empty.")
adata.file.close()

# == ADNC levels ==============================================================
_adnc_raw = sorted(obs[adnc_col].astype(str).dropna().unique().tolist())
print(f"\n  ADNC column unique values (raw): {_adnc_raw}")

_canonical_lower = {c.lower(): c for c in _ADNC_CANONICAL}
_raw_to_canonical = {}
for _rv in _adnc_raw:
    _cl = _rv.lower().strip()
    if _cl in _canonical_lower:
        _raw_to_canonical[_rv] = _canonical_lower[_cl]
    elif _cl in ('not ad', 'reference', 'no ad', 'no_ad'):
        _raw_to_canonical[_rv] = 'Not AD'
    else:
        _raw_to_canonical[_rv] = _rv

_matched_canonical = [c for c in _ADNC_CANONICAL
                      if c in _raw_to_canonical.values()]
_unmatched         = [v for v in _raw_to_canonical.values()
                      if v not in _ADNC_CANONICAL]
adnc_order    = _matched_canonical + sorted(set(_unmatched))
adnc_rank_map = {a: i for i, a in enumerate(adnc_order)}

_set2 = get_cmap_compat("Set2")
_fallback_cols = [_set2(i / 8.0) for i in range(max(len(adnc_order), 1))]
adnc_pal = {}
for _ai, _av in enumerate(adnc_order):
    adnc_pal[_av] = _ADNC_PAL_CANONICAL.get(
        _av, _mcolors.to_hex(_fallback_cols[_ai]))

# Map raw ADNC labels to the canonical ones.
if _raw_to_canonical:
    obs[adnc_col] = obs[adnc_col].astype(str).map(
        lambda x: _raw_to_canonical.get(x, x))
    print(f"  ADNC remap applied: {_raw_to_canonical}")

print(f"  adnc_order resolved: {adnc_order}")
print(f"  adnc_rank_map      : {adnc_rank_map}")

# == Ordinal parsing audit, donor level ========================================
print(f"\n  Ordinal parsing audit, donor level, neocortical arm")
_don1 = obs.drop_duplicates(subset=[donor_col]).set_index(donor_col)
_ord_specs = [("ADNC", adnc_col, adnc_rank_map),
              ("Braak", braak_col, BRAAK_MAP)]
for _nm, _cl, _mp in [("CERAD", cerad_col, CERAD_MAP),
                      ("Thal",  thal_col,  THAL_MAP),
                      ("CAA",   caa_col,   CAA_MAP)]:
    if _cl is not None:
        _ord_specs.append((_nm, _cl, _mp))
    else:
        print(f"    {_nm:<7} column absent - its clinical test will be dropped.")
for _nm, _cl, _mp in _ord_specs:
    print(f"    -- {_nm}  <- obs[{_cl!r}]  dtype_in={_don1[_cl].dtype}")
    audit_ordinal(_nm, _don1[_cl], _mp)

# Drop clinical tests whose column did not resolve
_ck_to_col = {"adnc_rank": adnc_col, "braak_num": braak_col,
              "cerad_num": cerad_col, "thal_num": thal_col,
              "caa_num": caa_col, "apoe_e4_dose": apoe_col,
              "mmse_num": mmse_col, "age_num": age_col}
_dropped_tests = [t for t in CLINICAL_TESTS if _ck_to_col.get(t[1]) is None]
if _dropped_tests:
    print(f"\n  dropping clinical tests with no source column: "
          f"{[t[0] for t in _dropped_tests]}")
    CLINICAL_TESTS = [t for t in CLINICAL_TESTS
                      if _ck_to_col.get(t[1]) is not None]
print(f"  Clinical tests active: {[t[0] for t in CLINICAL_TESTS]}")
print(f"\n  Section 1 done [{time.time()-t0:.1f}s]")


# %% SECTION 2: Score computation - all three base methods simultaneously
print("\nSECTION 2: Score computation - M1, M2, M3")
t1 = time.time()

np.random.seed(42)
random.seed(42)
adata_full = sc.read_h5ad(FULL_H5AD)
X = adata_full.X
if not sp.issparse(X): X = sp.csr_matrix(X)
X = X.tocsr()
print(f"  Full matrix: {X.shape}  [{time.time()-t1:.1f}s]")

# == Numeric normalisation check: expm1(X) row sums, max and min ===============
print(f"\n  normalisation check (target_sum={NORM_TARGET_SUM:,.0f}, "
      f"ceiling={NORM_CEILING:.4f})")
_rng_norm = np.random.default_rng(42)
_ridx = _rng_norm.choice(X.shape[0], size=min(NORM_CHECK_N, X.shape[0]),
                         replace=False)
_Xs   = X[np.sort(_ridx), :]
# expm1 on the sparse data only (expm1(0) = 0).
_Xd   = _Xs.copy().astype(np.float64)
_Xd.data = np.expm1(_Xd.data)
_rowsum = np.asarray(_Xd.sum(axis=1)).ravel()
_rs_med = float(np.median(_rowsum))
_rs_cv  = float(np.std(_rowsum) / max(abs(np.mean(_rowsum)), 1e-30))
_xmax_s = float(_Xs.max())
_xmin_s = float(_Xs.min())
_xmax_g = float(X.max())     # exact global max
print(f"    T1 expm1 row-sum : median={_rs_med:,.3f}  CV={_rs_cv:.2e}  "
      f"(n={len(_rowsum):,} sampled)")
print(f"    T2 max(X) global : {_xmax_g:.4f}   ceiling {NORM_CEILING:.4f}   "
      f"headroom {NORM_CEILING-_xmax_g:+.4f}")
print(f"       max(X) sampled: {_xmax_s:.4f}")
print(f"    T3 min(X)        : {_xmin_s:.6f}")
_norm_fail = []
if abs(_rs_med - NORM_TARGET_SUM) / NORM_TARGET_SUM > 1e-4:
    _norm_fail.append(f"row-sum median {_rs_med:,.3f} != {NORM_TARGET_SUM:,.0f}")
if _rs_cv > 1e-6:
    _norm_fail.append(f"row-sum CV {_rs_cv:.2e} too high for a fixed target_sum")
if _xmax_g > NORM_CEILING + NORM_TOL:
    _norm_fail.append(f"max(X) {_xmax_g:.4f} exceeds ceiling {NORM_CEILING:.4f}")
if _xmin_s < -NORM_TOL:
    _norm_fail.append(f"negative entries present (min {_xmin_s:.6f})")
if _norm_fail:
    _msg = ("NORMALISATION GUARDRAIL FAILED:\n    - "
            + "\n    - ".join(_norm_fail)
            + "\n  X is not ln(UP10K+1).")
    if NORM_ABORT_ON_FAIL:
        raise SystemExit("ABORT - " + _msg)
    print("  [!] " + _msg)
else:
    print(f"    [OK] X is ln(UP10K+1) on CP10K.")
del _Xs, _Xd; gc.collect()

X_nc = X[nc_cell_indices, :]
print(f"\n  Neocortical submatrix: {X_nc.shape}")

up_idx = np.array([var_names.index(g) for g in up_avail])
dn_idx = np.array([var_names.index(g) for g in dn_avail])

# == Marker genes for Section 2c, taken while X is loaded ======================
_mk_want = sorted({g for _v in MIC13_MARKER_CANDIDATES.values() for g in _v})
MARKER_GENES = [g for g in _mk_want if g in var_names]
_mk_missing  = [g for g in _mk_want if g not in var_names]
if _mk_missing:
    print(f"  [2c] marker genes absent from this gene space: {_mk_missing}")
if MARKER_GENES:
    _mk_idx = np.array([var_names.index(g) for g in MARKER_GENES])
    MARKER_X = np.asarray(X_nc[:, _mk_idx].todense(), dtype=np.float32)
    print(f"  [2c] marker matrix for the dilution curve: {MARKER_X.shape} "
          f"({len(MARKER_GENES)} genes, ln(UP10K+1) as stored)")
else:
    MARKER_X = None
    print("  [2c] no marker genes found - Section 2c will be skipped.")

print("  Extracting UP/DOWN gene matrices ...")
X_up = np.asarray(X_nc[:, up_idx].todense(), dtype=np.float32)
X_dn = np.asarray(X_nc[:, dn_idx].todense(), dtype=np.float32)
del X, X_nc; gc.collect()
print(f"  [{time.time()-t1:.1f}s]")

print("  Z-scoring ...")
Z_up = (X_up - X_up.mean(0)) / (X_up.std(0) + 1e-8)
Z_dn = (X_dn - X_dn.mean(0)) / (X_dn.std(0) + 1e-8)
del X_up, X_dn; gc.collect()

# Compute net score for each base method
score_cols = {}
for mkey, (wts_up_dict, wts_dn_dict) in SCORERS.items():
    wts_up = np.array([wts_up_dict.get(g, 0.0) for g in up_avail], dtype=np.float32)
    wts_dn = np.array([wts_dn_dict.get(g, 0.0) for g in dn_avail], dtype=np.float32)
    sc_arr = (Z_up @ wts_up) - (Z_dn @ wts_dn)
    score_cols[mkey] = sc_arr
    print(f"  {mkey:15s}: mean={sc_arr.mean():.4f}  std={sc_arr.std():.4f}"
          f"  95th={np.percentile(sc_arr,95):.4f}")

del Z_up, Z_dn; gc.collect()
del adata_full; gc.collect()

# Attach to obs_full
obs_full = obs.copy()
for mkey, sc_arr in score_cols.items():
    obs_full[f"score_{mkey}"] = sc_arr

# == Intersection scores for M4, M5 and M6 =====================================
n_cells = len(obs_full)
_prank = {}
from scipy.stats import rankdata as _rankdata
for _mk in SCORERS:
    _prank[_mk] = _rankdata(score_cols[_mk], method="average") / n_cells

for _int_key, (_pa, _pb) in INTERSECTION_PAIRS.items():
    _int_score = np.minimum(_prank[_pa], _prank[_pb])
    score_cols[_int_key] = _int_score
    obs_full[f"score_{_int_key}"] = _int_score
    _mask_a = score_cols[_pa] >= np.percentile(score_cols[_pa], (1-TOP_PCT)*100)
    _mask_b = score_cols[_pb] >= np.percentile(score_cols[_pb], (1-TOP_PCT)*100)
    _n_intersect = int((_mask_a & _mask_b).sum())
    print(f"  {_int_key}: n cells in intersection = {_n_intersect:,} "
          f"({100*_n_intersect/n_cells:.2f}% of all neocortical cells)")

ALL_METHODS = list(SCORERS.keys()) + list(INTERSECTION_PAIRS.keys())
print(f"  Total methods to evaluate: {len(ALL_METHODS)}")



# %% SECTION 2b: Cutoff-selection criterion, fine scan
# Scans the 1% rule over a fine grid of cutoffs. Raw verdict: all supertypes;
# guarded verdict: supertypes with n >= RULE_MIN_N. margin_cells is the number
# of nuclei that would flip the verdict.
print("\n" + "=" * 74)
print("SECTION 2b: Cutoff-selection criterion - fine scan")
print("=" * 74)
_t2b = time.time()

_stypes_scan = sorted(obs_full[stype_col].astype(str).unique())
_st_masks    = {st: (obs_full[stype_col].astype(str) == st).values
                for st in _stypes_scan}
_n_tot_st    = {st: int(_st_masks[st].sum()) for st in _stypes_scan}
_others      = [st for st in _stypes_scan if st not in RULE_ENRICHED_STYPES]

_missing_enr = [st for st in RULE_ENRICHED_STYPES if st not in _stypes_scan]
if _missing_enr:
    print(f"  [!] rule references supertypes absent from this arm: {_missing_enr}")

# == Smallest %top each supertype can report ==================================
print(f"\n  GRANULARITY - smallest non-zero %top each supertype can report")
print(f"  rule: every supertype except {RULE_ENRICHED_STYPES} must sit "
      f"<= {RULE_MAX_PCT:.1f}%")
print(f"\n  {'supertype':<26} {'n_cells':>8} {'min %>0':>9} {'cells allowed':>14}  status")
print("  " + "-" * 78)
_unsatisfiable = []
for st in sorted(_stypes_scan, key=lambda s: -_n_tot_st[s]):
    n   = _n_tot_st[st]
    mnp = 100.0 / max(n, 1)
    allowed = int(np.floor(RULE_MAX_PCT / 100.0 * n))
    if st in RULE_ENRICHED_STYPES:
        status = "exempt (expected enriched)"
    elif n < RULE_MIN_N:
        status = f"UNDER-POWERED (n < {RULE_MIN_N})"
        _unsatisfiable.append(st)
    elif mnp > RULE_MAX_PCT:
        status = "rule = 'exactly zero cells'"
        _unsatisfiable.append(st)
    else:
        status = "ok"
    print(f"  {st:<26} {n:>8,} {mnp:>8.2f}% {allowed:>14,}  {status}")
if _unsatisfiable:
    print(f"\n  [!] {len(_unsatisfiable)} supertype(s) cannot express the rule "
          f"meaningfully: {_unsatisfiable}")


def _mask_at(mkey, pct):
    """Top-pct mask for any scorer; intersections are the AND of the two parent
    masks (the canonical M4 definition)."""
    if mkey in INTERSECTION_PAIRS:
        _a, _b = INTERSECTION_PAIRS[mkey]
        return ((score_cols[_a] >= np.percentile(score_cols[_a], (1 - pct) * 100)) &
                (score_cols[_b] >= np.percentile(score_cols[_b], (1 - pct) * 100)))
    _s = score_cols[mkey]
    return _s >= np.percentile(_s, (1 - pct) * 100)


# == Scan ======================================================================
print(f"\n  scanning {len(SCAN_GRID)} cutoffs x {len(ALL_METHODS)} scorers "
      f"x {len(_stypes_scan)} supertypes ...")
_scan_rows, _verdict_rows = [], []
for _p in SCAN_GRID:
    for _mk in ALL_METHODS:
        _m  = _mask_at(_mk, float(_p))
        _nm = int(_m.sum())
        _pct_by_st, _margin_raw, _margin_gd = {}, [], []
        _viol_raw, _viol_gd = [], []
        for st in _stypes_scan:
            _n   = _n_tot_st[st]
            _nt  = int((_m & _st_masks[st]).sum())
            _pct = 100.0 * _nt / max(_n, 1)
            _pct_by_st[st] = _pct
            _allow = int(np.floor(RULE_MAX_PCT / 100.0 * _n))
            _scan_rows.append({
                "top_pct": round(float(_p) * 100, 3), "method": _mk,
                "supertype": st, "n_total": _n, "n_top": _nt,
                "pct_top": round(_pct, 4),
                "fold": round(_pct / (100.0 * float(_p)), 4) if _p > 0 else np.nan,
                "min_nonzero_pct": round(100.0 / max(_n, 1), 4),
                "cells_allowed_under_rule": _allow,
                "margin_cells": _allow - _nt,
                "rule_exempt": st in RULE_ENRICHED_STYPES,
                "rule_guarded_in": (st not in RULE_ENRICHED_STYPES
                                    and _n >= RULE_MIN_N),
                "violates": (st not in RULE_ENRICHED_STYPES
                             and _pct > RULE_MAX_PCT),
            })
            if st in RULE_ENRICHED_STYPES:
                continue
            _margin_raw.append(_allow - _nt)
            if _pct > RULE_MAX_PCT:
                _viol_raw.append(st)
            if _n >= RULE_MIN_N:
                _margin_gd.append(_allow - _nt)
                if _pct > RULE_MAX_PCT:
                    _viol_gd.append(st)
        _rank = sorted(_pct_by_st, key=lambda s: -_pct_by_st[s])
        _verdict_rows.append({
            "top_pct": round(float(_p) * 100, 3), "method": _mk,
            "n_top_total": _nm,
            "pass_raw": len(_viol_raw) == 0,
            "pass_guarded": len(_viol_gd) == 0,
            "n_violators_raw": len(_viol_raw),
            "n_violators_guarded": len(_viol_gd),
            "violators_raw": ";".join(_viol_raw),
            "violators_guarded": ";".join(_viol_gd),
            "margin_cells_raw": min(_margin_raw) if _margin_raw else np.nan,
            "margin_cells_guarded": min(_margin_gd) if _margin_gd else np.nan,
            "enriched_are_top2": set(_rank[:2]) == set(RULE_ENRICHED_STYPES),
            **{f"pct_{s}": round(_pct_by_st[s], 3) for s in _stypes_scan},
        })

scan_df    = pd.DataFrame(_scan_rows)
verdict_df = pd.DataFrame(_verdict_rows)
scan_df.to_csv(os.path.join(ARM_DIR, "S2b_cutoff_scan.csv"), index=False)
verdict_df.to_csv(os.path.join(ARM_DIR, "S2b_cutoff_rule_verdict.csv"), index=False)
print(f"  Saved S2b_cutoff_scan.csv ({len(scan_df):,} rows) and "
      f"S2b_cutoff_rule_verdict.csv ({len(verdict_df):,} rows)  "
      f"[{time.time()-_t2b:.1f}s]")


def _window(_flag_col, _mk=None):
    """Contiguous admissible cutoff ranges for one scorer."""
    _v = verdict_df[verdict_df["method"] == (_mk or SCAN_METHOD)].sort_values("top_pct")
    _ok = _v[_v[_flag_col]]["top_pct"].tolist()
    if not _ok:
        return []
    _step = round(float(np.diff(np.sort(SCAN_GRID)).min()) * 100, 3)
    _runs, _cur = [], [_ok[0], _ok[0]]
    for _x in _ok[1:]:
        if round(_x - _cur[1], 3) <= _step + 1e-9:
            _cur[1] = _x
        else:
            _runs.append(tuple(_cur)); _cur = [_x, _x]
    _runs.append(tuple(_cur))
    return _runs


print(f"\n  ADMISSIBLE WINDOWS for {SCAN_METHOD} "
      f"(grid {SCAN_GRID.min()*100:.1f}% to {SCAN_GRID.max()*100:.1f}%)")
_wr, _wg = _window("pass_raw"), _window("pass_guarded")
print(f"    raw rule     : "
      + (", ".join(f"{a:.1f}-{b:.1f}%" for a, b in _wr) if _wr else "NONE"))
print(f"    guarded rule : "
      + (", ".join(f"{a:.1f}-{b:.1f}%" for a, b in _wg) if _wg else "NONE"))
if _wg and _wr and (sum(b - a for a, b in _wg) > sum(b - a for a, b in _wr) + 1e-9):
    print(f"    -> the guarded window is wider than the raw window.")

# == Verdict at each sweep cutoff ===============================================
print(f"\n  VERDICT AT EACH SWEEP CUTOFF  [{SCAN_METHOD}]")
print(f"  {'cutoff':>7} {'n_top':>8} {'raw':>5} {'guard':>6} "
      f"{'margin_raw':>11} {'margin_gd':>10}  violators (raw)")
print("  " + "-" * 88)
for _p in sorted(set(list(SWEEP_PCTS) + [PRIMARY_TOP_PCT])):
    _row = verdict_df[(verdict_df["method"] == SCAN_METHOD) &
                      (np.isclose(verdict_df["top_pct"], _p * 100, atol=1e-6))]
    if _row.empty:
        _m  = _mask_at(SCAN_METHOD, float(_p))
        _pc = {st: 100.0 * int((_m & _st_masks[st]).sum()) / max(_n_tot_st[st], 1)
               for st in _stypes_scan}
        _vr = [st for st in _others if _pc[st] > RULE_MAX_PCT]
        print(f"  {_p*100:>6.1f}% {int(_m.sum()):>8,} "
              f"{str(len(_vr)==0):>5} {'-':>6} {'-':>11} {'-':>10}  "
              f"{';'.join(_vr) if _vr else '(none)'}   [off-grid]")
        continue
    _r = _row.iloc[0]
    _star = "  <- PRIMARY" if abs(_p - PRIMARY_TOP_PCT) < 1e-9 else ""
    print(f"  {_p*100:>6.1f}% {int(_r['n_top_total']):>8,} "
          f"{str(bool(_r['pass_raw'])):>5} {str(bool(_r['pass_guarded'])):>6} "
          f"{int(_r['margin_cells_raw']):>11,} {int(_r['margin_cells_guarded']):>10,}  "
          f"{_r['violators_raw'] if _r['violators_raw'] else '(none)'}{_star}")
print(f"\n  margin_cells = nuclei one supertype could add before the rule breaks.")

# == Figure ====================================================================
_vm = verdict_df[verdict_df["method"] == SCAN_METHOD].sort_values("top_pct")
_sm = scan_df[scan_df["method"] == SCAN_METHOD]
fig2b, (axA, axB) = plt.subplots(1, 2, figsize=(16, 6))
for st in sorted(_stypes_scan, key=lambda s: -_n_tot_st[s]):
    _d = _sm[_sm["supertype"] == st].sort_values("top_pct")
    _small = _n_tot_st[st] < RULE_MIN_N
    if st in RULE_ENRICHED_STYPES:
        _c, _lw, _ls, _z = ARM_COLOR if st == TARGET_STYPE else "#B08A3E", 2.4, "-", 5
    else:
        _c, _lw, _ls, _z = "#9E9E9E", 1.3, ("--" if _small else "-"), 3
    axA.plot(_d["top_pct"], _d["pct_top"], color=_c, lw=_lw, ls=_ls, zorder=_z,
             label=f"{st} (n={_n_tot_st[st]:,})" + ("  [under-powered]" if _small else ""))
axA.axhline(RULE_MAX_PCT, color="#1A9850", lw=1.8, ls="--", zorder=4,
            label=f"rule: {RULE_MAX_PCT:.1f}% ceiling for non-exempt types")
for _p in SWEEP_PCTS:
    axA.axvline(_p * 100, color="black", lw=1.0,
                ls=("-" if abs(_p - PRIMARY_TOP_PCT) < 1e-9 else ":"), alpha=0.55)
axA.set_yscale("symlog", linthresh=1.0)
axA.set_xlabel("top-% cutoff")
axA.set_ylabel("% of that supertype's cells in the top set")
axA.set_title(f"Per-supertype %top vs cutoff  [{SCAN_METHOD}]\n"
              f"neocortical (DFC/MTG/AnG/STG/V1C/ITG) - symlog y, linear below 1%", fontsize=10)
axA.legend(fontsize=7, loc="upper left")
axA.spines[["top", "right"]].set_visible(False)

axB.plot(_vm["top_pct"], _vm["margin_cells_raw"], color="#9E9E9E", lw=1.8,
         label="margin (raw rule, all supertypes)")
axB.plot(_vm["top_pct"], _vm["margin_cells_guarded"], color=ARM_COLOR, lw=2.2,
         label=f"margin (guarded, n >= {RULE_MIN_N})")
axB.axhline(0, color="#1A9850", lw=1.8, ls="--", label="rule breaks below this")
for _flag, _c, _lbl in [("pass_raw", "#9E9E9E", "raw pass"),
                        ("pass_guarded", ARM_COLOR, "guarded pass")]:
    for _a, _b in _window(_flag):
        axB.axvspan(_a, _b, color=_c, alpha=0.12, lw=0)
for _p in SWEEP_PCTS:
    axB.axvline(_p * 100, color="black", lw=1.0,
                ls=("-" if abs(_p - PRIMARY_TOP_PCT) < 1e-9 else ":"), alpha=0.55)
axB.set_xlabel("top-% cutoff")
axB.set_ylabel("margin, in nuclei, before the rule breaks")
axB.set_title("Fragility of the selection rule\n"
              "shaded = admissible; solid vline = PRIMARY, dotted = other sweep points",
              fontsize=10)
axB.legend(fontsize=8)
axB.spines[["top", "right"]].set_visible(False)
plt.suptitle(f"Cutoff-selection criterion - neocortical (DFC/MTG/AnG/STG/V1C/ITG)", fontsize=12)
plt.tight_layout()
fig2b.savefig(os.path.join(ARM_DIR, "S2b_cutoff_scan.png"), dpi=200,
              bbox_inches="tight")
plt.show()
plt.close(fig2b)
print(f"  Saved S2b_cutoff_scan.png")


# %% SECTION 2c: Marker-dilution curve
# Fold enrichment of each marker panel (vs the arm mean) among the top-x% cells
# (cumulative) and among the cells added at each step (marginal).
print("\n" + "=" * 74)
print("SECTION 2c: Marker-dilution curve")
print("=" * 74)

if MARKER_X is None or not MARKER_GENES:
    print("  SKIPPED - no marker genes available.")
    marker_df = pd.DataFrame()
else:
    _t2c = time.time()
    _gi = {g: i for i, g in enumerate(MARKER_GENES)}

    # == Held-out panels must not be in the scoring panel ========================
    _scored = set(MIC13_UP) | set(MIC13_DOWN)
    _panels, _demoted = {}, []
    for _pname, _glist in MIC13_MARKER_CANDIDATES.items():
        _keep = []
        for _g in _glist:
            if _g not in _gi:
                continue
            if _pname != "inpanel_up" and _g in _scored:
                _demoted.append((_pname, _g))
                continue
            _keep.append(_g)
        _panels[_pname] = _keep
    if _demoted:
        print(f"  [!] demoted from a held-out panel (in the scoring panel):")
        for _pn, _g in _demoted:
            print(f"        {_g}  (declared under {_pn})")
    print(f"\n  Panels resolved:")
    for _pname, _keep in _panels.items():
        _tag = {"heldout_up":  "evidence - not used in scoring",
                "inpanel_up":  "reference only - circular by construction",
                "homeostatic": "negative control - should rise with dilution"}[_pname]
        print(f"    {_pname:<12} n={len(_keep)}  {_keep}")
        print(f"      -> {_tag}")
    if len(_panels.get("heldout_up", [])) < MARKER_MIN_PANEL:
        print(f"\n  [!] only {len(_panels.get('heldout_up', []))} usable held-out "
              f"gene(s), below MARKER_MIN_PANEL={MARKER_MIN_PANEL}.")

    # background = arm mean per gene; every fold is relative to this
    _bg = MARKER_X.mean(axis=0)                      # (n_marker_genes,)
    _bg_safe = np.where(_bg > 1e-9, _bg, np.nan)

    def _panel_fold(mask, panel):
        """Mean over panel genes of (mean expression in `mask`) / (arm mean)."""
        _g = [_gi[x] for x in panel if x in _gi]
        if not _g or int(mask.sum()) == 0:
            return np.nan, np.nan
        _sub = MARKER_X[mask][:, _g]
        _fold = np.nanmean(_sub.mean(axis=0) / _bg_safe[_g])
        _det  = float((_sub > 0).mean()) * 100.0
        return float(_fold), _det

    # == Masks across the grid for the exported mask (M4) ========================
    _grid = np.sort(np.asarray(SCAN_GRID, dtype=float))
    _masks = {float(_p): _mask_at(SCAN_METHOD, float(_p)) for _p in _grid}

    _step_pp = float(np.diff(_grid).min()) * 100 if len(_grid) > 1 else 0.05
    _back = max(1, int(round(MARKER_BAND_PP / _step_pp)))
    print(f"  marginal band = {MARKER_BAND_PP:.2f} pp wide "
          f"({_back} grid steps of {_step_pp:.2f} pp), moving window")

    _rows = []
    for _k, _p in enumerate(_grid):
        _m = _masks[float(_p)]
        _band = _m if _k < _back else (_m & ~_masks[float(_grid[_k - _back])])
        _r = {"top_pct": round(float(_p) * 100, 3),
              "method": SCAN_METHOD,
              "n_cum": int(_m.sum()),
              "n_band": int(_band.sum())}
        for _pname, _keep in _panels.items():
            _fc, _dc = _panel_fold(_m, _keep)
            _fb, _db = _panel_fold(_band, _keep)
            _r[f"cum_fold_{_pname}"]  = round(_fc, 4) if _fc == _fc else np.nan
            _r[f"marg_fold_{_pname}"] = round(_fb, 4) if _fb == _fb else np.nan
            _r[f"cum_det_{_pname}"]   = round(_dc, 2) if _dc == _dc else np.nan
        for _g in MARKER_GENES:      # per-gene detail, for auditing the panels
            _r[f"cum_mean_{_g}"] = round(
                float(MARKER_X[_m][:, _gi[_g]].mean()) if int(_m.sum()) else np.nan, 5)
        _rows.append(_r)
    marker_df = pd.DataFrame(_rows)
    marker_df.to_csv(os.path.join(ARM_DIR, "S2c_marker_dilution.csv"), index=False)
    print(f"\n  Saved S2c_marker_dilution.csv ({len(marker_df):,} rows)  "
          f"[{time.time()-_t2c:.1f}s]")

    # == Knee of the cumulative curve (largest distance from the chord) ===========
    def _chord_knee(x, y):
        """Kneedle-style: the point furthest from the straight line joining the
        first and last points. Parameter-free and robust to the y-scale."""
        _ok = np.isfinite(y)
        if _ok.sum() < 4:
            return np.nan
        x2, y2 = np.asarray(x)[_ok], np.asarray(y)[_ok]
        _xn = (x2 - x2[0]) / max(x2[-1] - x2[0], 1e-12)
        # Signed normalisation, so decreasing curves also map 0 -> 1.
        _den = y2[-1] - y2[0]
        if abs(_den) < 1e-12:
            return np.nan
        _yn = (y2 - y2[0]) / _den
        return float(x2[int(np.argmax(np.abs(_yn - _xn)))])

    # == Crossings of the marginal curve ==========================================
    def _last_above(x, y, thr):
        _ok = np.isfinite(y)
        if _ok.sum() == 0:
            return np.nan
        x2, y2 = np.asarray(x)[_ok], np.asarray(y)[_ok]
        _hit = np.where(y2 >= thr)[0]
        return float(x2[_hit[-1]]) if len(_hit) else np.nan

    _xs = marker_df["top_pct"].values
    print(f"\n  CUTOFF ANCHORS from the held-out panel "
          f"({_panels.get('heldout_up', [])})")
    _anchors = {}
    for _pname in _panels:
        _cum, _marg = (marker_df[f"cum_fold_{_pname}"].values,
                       marker_df[f"marg_fold_{_pname}"].values)
        _knee = _chord_knee(_xs, _cum)
        _halfmax = np.nan
        _fin = _marg[np.isfinite(_marg)]
        if len(_fin):
            _half = 1.0 + (_fin[0] - 1.0) / 2.0        # halfway back to arm mean
            _halfmax = _last_above(_xs, _marg, _half)
        _cross = _last_above(_xs, _marg, MARKER_MARGINAL_MIN)
        _anchors[_pname] = dict(knee=_knee, halfmax=_halfmax, cross=_cross)
    print(f"  {'panel':<12} {'knee(cum)':>10} {'half-max(marg)':>15} "
          f"{'last >= ' + format(MARKER_MARGINAL_MIN, '.1f') + 'x':>14}")
    print("  " + "-" * 56)
    for _pname, _a in _anchors.items():
        _f = lambda v: f"{v:.2f}%" if v == v else "  n/a"
        print(f"  {_pname:<12} {_f(_a['knee']):>10} {_f(_a['halfmax']):>15} "
              f"{_f(_a['cross']):>14}")
    _ho = _anchors.get("heldout_up", {})
    if _panels.get("homeostatic") and _panels.get("heldout_up"):
        # Spearman trend of each marginal curve against the cutoff.
        def _trend(col):
            _v = marker_df[col].values
            _ok = np.isfinite(_v)
            if _ok.sum() < 5:
                return np.nan, np.nan
            _r, _pv = spearmanr(marker_df["top_pct"].values[_ok], _v[_ok])
            return float(_r), float(_pv)
        _rh, _ph = _trend("marg_fold_homeostatic")
        _ru, _pu = _trend("marg_fold_heldout_up")
        print(f"\n  DIRECTION CHECK (Spearman of the marginal curve vs cutoff):")
        print(f"    held-out Mic.13 markers : rho = {_ru:+.3f}  p = {_pu:.3g}")
        print(f"    homeostatic markers     : rho = {_rh:+.3f}  p = {_ph:.3g}")
        if _ru < 0 and _rh > 0:
            print(f"  Both as expected (opposite directions).")
        elif _ru < 0 and _rh <= 0:
            print(f"  [!] Held-out markers dilute but homeostatic markers do not rise.")
        else:
            print(f"  [!] The held-out panel does not dilute. Either the score is not")

    # == Figure ==================================================================
    _sty = {"heldout_up":  (ARM_COLOR, 2.6, "-",  "held-out Mic.13 markers - EVIDENCE"),
            "inpanel_up":  ("#9E9E9E", 1.6, "--", "in scoring panel - circular, reference only"),
            "homeostatic": ("#1A9850", 2.0, "-.", "homeostatic - negative control")}
    fig2c, (axC, axM) = plt.subplots(1, 2, figsize=(16, 6))
    for _pname, _keep in _panels.items():
        if not _keep:
            continue
        _c, _lw, _ls, _lbl = _sty[_pname]
        axC.plot(_xs, marker_df[f"cum_fold_{_pname}"], color=_c, lw=_lw, ls=_ls,
                 label=f"{_lbl}  ({', '.join(_keep)})")
        axM.plot(_xs, marker_df[f"marg_fold_{_pname}"], color=_c, lw=_lw, ls=_ls,
                 label=_lbl)
    for _ax in (axC, axM):
        _ax.axhline(1.0, color="black", lw=1.2, ls=":", label="arm mean (no enrichment)")
        for _p in SWEEP_PCTS:
            _ax.axvline(_p * 100, color="black", lw=1.0,
                        ls=("-" if abs(_p - PRIMARY_TOP_PCT) < 1e-9 else ":"), alpha=0.5)
        _ax.set_xlabel("top-% cutoff")
        _ax.spines[["top", "right"]].set_visible(False)
    if _ho.get("knee") == _ho.get("knee"):
        axC.axvline(_ho["knee"], color=ARM_COLOR, lw=1.6, ls="--", alpha=0.8)
        axC.annotate(f"knee {_ho['knee']:.2f}%", xy=(_ho["knee"], axC.get_ylim()[1]),
                     xytext=(4, -12), textcoords="offset points",
                     color=ARM_COLOR, fontsize=9)
    axM.axhline(MARKER_MARGINAL_MIN, color=ARM_COLOR, lw=1.4, ls="--", alpha=0.7,
                label=f"marginal floor {MARKER_MARGINAL_MIN:.1f}x")
    if _ho.get("cross") == _ho.get("cross"):
        axM.axvline(_ho["cross"], color=ARM_COLOR, lw=1.6, ls="--", alpha=0.8)
        axM.annotate(f"last >= {MARKER_MARGINAL_MIN:.1f}x: {_ho['cross']:.2f}%",
                     xy=(_ho["cross"], axM.get_ylim()[1]), xytext=(4, -12),
                     textcoords="offset points", color=ARM_COLOR, fontsize=9)
    axC.set_ylabel("fold-enrichment vs arm mean, cells KEPT")
    axC.set_title("Cumulative - what you get if you pick x\n"
                  "solid vline = PRIMARY, dotted = other sweep points", fontsize=10)
    axM.set_ylabel("fold-enrichment vs arm mean, cells ADMITTED at the margin")
    axM.set_title("Marginal - is the NEXT cell still Mic.13-like?\n"
                  "falls to 1.0 = the score has stopped discriminating", fontsize=10)
    axC.legend(fontsize=7.5, loc="best")
    axM.legend(fontsize=7.5, loc="best")
    plt.suptitle(f"Marker dilution vs cutoff - neocortical (DFC/MTG/AnG/STG/V1C/ITG)  [{SCAN_METHOD}]\n"
                 f"outcome-blind: no pathology variable enters this figure", fontsize=12)
    plt.tight_layout()
    fig2c.savefig(os.path.join(ARM_DIR, "S2c_marker_dilution.png"), dpi=200,
                  bbox_inches="tight")
    plt.show()
    plt.close(fig2c)
    print(f"  Saved S2c_marker_dilution.png")
    del _masks; gc.collect()




# %% SECTION 2d: Two-component mixture model on the score distribution
# 1- vs 2-component Gaussian mixture on each base score, against a skew-normal
# null. BIMODAL needs dBIC and separation both beyond the null and separation
# >= MIX_SEP_MIN. Detection only: the implied boundary is not a cutoff.
print("\n" + "=" * 74)
print("SECTION 2d: Mixture model - is the score bimodal or continuous?")
print("=" * 74)
_t2d = time.time()

try:
    from sklearn.mixture import GaussianMixture as _GMM
    _HAVE_SK = True
except Exception as _e_sk:                                   # pragma: no cover
    _HAVE_SK = False
    print(f"  [SKIP] scikit-learn unavailable ({_e_sk}). Section 2d needs it.")

if _HAVE_SK:
    from scipy import stats as _sps

    _mix_rows = []
    _mix_rng  = np.random.default_rng(MIX_SEED)

    def _fit_bic(x, k, n_init):
        g = _GMM(n_components=k, covariance_type="full", n_init=n_init,
                 random_state=MIX_SEED, reg_covar=1e-6)
        g.fit(x.reshape(-1, 1))
        return g, float(g.bic(x.reshape(-1, 1)))

    print(f"  subsample = {MIX_SUBSAMPLE:,} cells")
    print(f"  null      = {MIX_NULL_B} skew-normal draws, seed {MIX_SEED}")

    for _mk in MIX_METHODS:
        if _mk not in score_cols:
            print(f"  [!] {_mk} not a base scorer, skipped")
            continue
        _s_full = np.asarray(score_cols[_mk], dtype=float)
        _s_full = _s_full[np.isfinite(_s_full)]
        _n_use  = int(min(MIX_SUBSAMPLE, _s_full.size))
        _idx    = _mix_rng.choice(_s_full.size, size=_n_use, replace=False)
        _s      = np.sort(_s_full[_idx])
        _s      = (_s - _s.mean()) / (_s.std(ddof=1) + 1e-12)

        _g1, _b1 = _fit_bic(_s, 1, 1)
        _g2, _b2 = _fit_bic(_s, 2, MIX_N_INIT)
        _g3, _b3 = _fit_bic(_s, 3, MIX_N_INIT)
        _dbic    = _b1 - _b2                    # positive favours 2 components

        _mu = _g2.means_.ravel()
        _sd = np.sqrt(_g2.covariances_.ravel())
        _w  = _g2.weights_.ravel()
        _hi = int(np.argmax(_mu))
        _lo = 1 - _hi
        _sep = float(abs(_mu[_hi] - _mu[_lo]) /
                     np.sqrt((_sd[_hi] ** 2 + _sd[_lo] ** 2) / 2.0))
        _w_hi = float(_w[_hi])

        # Boundary: last point where the upper-component posterior is below 0.5,
        # scanning down from the top (the posterior can be U-shaped).
        _post   = _g2.predict_proba(_s.reshape(-1, 1))[:, _hi]
        _below  = np.where(_post < 0.5)[0]
        _above  = np.where(_post >= 0.5)[0]
        _u_tail = bool(_post[0] >= 0.5)      # upper component owns the low tail
        if len(_below) and len(_above):
            _b_idx = int(_below[-1]) + 1
            _bound_pct = 100.0 * (1.0 - _b_idx / max(_n_use - 1, 1))
            _b_note = ("non-monotone posterior (U-shaped); last crossing used"
                       if _u_tail else "monotone posterior; single crossing")
        else:
            _bound_pct = np.nan
            _b_note = ("upper component never drops below 0.5 — undefined"
                       if len(_below) == 0 else
                       "upper component never reaches 0.5 — undefined")

        # Descriptive bimodality coefficient: > 0.555 is the usual flag.
        _g  = float(_sps.skew(_s, bias=False))
        _k  = float(_sps.kurtosis(_s, fisher=True, bias=False))
        _bc = (_g ** 2 + 1.0) / (_k + 3.0 * (_n_use - 1) ** 2 /
                                 ((_n_use - 2) * (_n_use - 3)))

        # == Null: skew-normal fitted to these scores =============================
        _a, _loc, _scl = _sps.skewnorm.fit(_s)
        _null_dbic, _null_sep = [], []
        for _b in range(MIX_NULL_B):
            _x = _sps.skewnorm.rvs(_a, loc=_loc, scale=_scl, size=_n_use,
                                   random_state=int(MIX_SEED + _b))
            _x = np.sort((_x - _x.mean()) / (_x.std(ddof=1) + 1e-12))
            _, _nb1 = _fit_bic(_x, 1, 1)
            _ng2, _nb2 = _fit_bic(_x, 2, 1)
            _null_dbic.append(_nb1 - _nb2)
            _nmu = _ng2.means_.ravel(); _nsd = np.sqrt(_ng2.covariances_.ravel())
            _null_sep.append(abs(_nmu[1] - _nmu[0]) /
                             np.sqrt((_nsd[0] ** 2 + _nsd[1] ** 2) / 2.0))
        _null_dbic = np.asarray(_null_dbic, dtype=float)
        _null_sep  = np.asarray(_null_sep, dtype=float)
        _p_dbic = float((np.sum(_null_dbic >= _dbic) + 1) / (MIX_NULL_B + 1))
        _p_sep  = float((np.sum(_null_sep >= _sep) + 1) / (MIX_NULL_B + 1))

        _verdict = ("BIMODAL — a boundary exists"
                    if (_p_dbic <= 0.05 and _p_sep <= 0.05
                        and _sep >= MIX_SEP_MIN)
                    else "CONTINUOUS — no boundary the data supports")

        _mix_rows.append(dict(
            method=_mk, n_used=_n_use,
            bic_k1=round(_b1, 1), bic_k2=round(_b2, 1), bic_k3=round(_b3, 1),
            dbic_1_minus_2=round(_dbic, 1),
            best_k=int(np.argmin([_b1, _b2, _b3]) + 1),
            w_upper=round(_w_hi, 4),
            implied_top_pct_from_weight=round(100.0 * _w_hi, 3),
            boundary_top_pct=(round(_bound_pct, 3)
                              if np.isfinite(_bound_pct) else np.nan),
            boundary_note=_b_note,
            separation=round(_sep, 4),
            bimodality_coef=round(float(_bc), 4),
            null_dbic_median=round(float(np.median(_null_dbic)), 1),
            null_dbic_max=round(float(_null_dbic.max()), 1),
            null_sep_median=round(float(np.median(_null_sep)), 4),
            p_dbic_vs_null=round(_p_dbic, 4),
            p_sep_vs_null=round(_p_sep, 4),
            verdict=_verdict))

        print(f"\n  -- {_mk}")
        print(f"     BIC  k=1 {_b1:,.0f} | k=2 {_b2:,.0f} | k=3 {_b3:,.0f}"
              f"   -> best k = {int(np.argmin([_b1,_b2,_b3])+1)}")
        _p_txt = (f"p < {1.0/(MIX_NULL_B+1):.4f}"
                  if _dbic > _null_dbic.max() else f"p = {_p_dbic:.4f}")
        print(f"     dBIC(1-2) {_dbic:,.0f}   null median {np.median(_null_dbic):,.0f} "
              f"max {_null_dbic.max():,.0f}   {_p_txt}")
        print(f"     separation {_sep:.3f}      null median "
              f"{np.median(_null_sep):.3f}          p = {_p_sep:.4f}")
        print(f"     upper component weight {100*_w_hi:.2f}%  |  posterior-0.5 "
              f"boundary at top {_bound_pct:.2f}%")
        print(f"     boundary shape: {_b_note}")
        print(f"     bimodality coefficient {_bc:.4f} (flag > 0.555)")
        print(f"     VERDICT: {_verdict}")

    mix_df = pd.DataFrame(_mix_rows)
    if len(mix_df):
        mix_df.to_csv(os.path.join(ARM_DIR, "S2d_mixture.csv"), index=False)
        print(f"\n  Saved S2d_mixture.csv  [{time.time()-_t2d:.1f}s]")

        _n_bim = int((mix_df["verdict"].str.startswith("BIMODAL")).sum())
        print(f"\n  {_n_bim} of {len(mix_df)} base scorers show a boundary the "
              f"null cannot produce.")

        # == Figure ==============================================================
        _nm = len(mix_df)
        fig2d, ax2d = plt.subplots(2, _nm, figsize=(5.6 * _nm, 8.4),
                                   squeeze=False)
        for _j, _mk in enumerate(mix_df["method"]):
            _s_full = np.asarray(score_cols[_mk], dtype=float)
            _s_full = _s_full[np.isfinite(_s_full)]
            _sv = (_s_full - _s_full.mean()) / (_s_full.std(ddof=1) + 1e-12)
            _row = mix_df[mix_df.method == _mk].iloc[0]
            axT = ax2d[0][_j]
            axT.hist(_sv, bins=200, density=True, color="#9E9E9E",
                     alpha=0.55, lw=0)
            _xs = np.linspace(_sv.min(), _sv.max(), 800)
            _gg = _GMM(n_components=2, covariance_type="full",
                       n_init=MIX_N_INIT, random_state=MIX_SEED)
            _gg.fit(_sv.reshape(-1, 1))
            for _c in range(2):
                axT.plot(_xs, _gg.weights_[_c] * _sps.norm.pdf(
                    _xs, _gg.means_.ravel()[_c],
                    np.sqrt(_gg.covariances_.ravel()[_c])),
                    lw=2.0, color=ARM_COLOR if _c == int(
                        np.argmax(_gg.means_.ravel())) else "#B08A3E")
            axT.axvline(np.percentile(_sv, (1 - PRIMARY_TOP_PCT) * 100),
                        color="#1A9850", lw=1.8, ls="--",
                        label=f"primary {PRIMARY_TOP_PCT*100:.1f}%")
            axT.set_yscale("log")
            axT.set_title(f"{_mk} - score density + 2-component fit\n"
                          f"separation {_row['separation']:.2f}, "
                          f"upper weight {_row['implied_top_pct_from_weight']:.1f}%",
                          fontsize=9)
            axT.set_xlabel("standardised score"); axT.set_ylabel("density (log)")
            axT.legend(fontsize=7); axT.spines[["top", "right"]].set_visible(False)

            axB = ax2d[1][_j]
            _a, _loc, _scl = _sps.skewnorm.fit(
                _sv[_mix_rng.choice(_sv.size, min(MIX_SUBSAMPLE, _sv.size),
                                    replace=False)])
            axB.axvline(_row["dbic_1_minus_2"], color=ARM_COLOR, lw=2.4,
                        label=f"observed dBIC {_row['dbic_1_minus_2']:,.0f}")
            axB.axvline(_row["null_dbic_median"], color="#9E9E9E", lw=1.8,
                        ls="--", label=f"null median {_row['null_dbic_median']:,.0f}")
            axB.axvline(_row["null_dbic_max"], color="#9E9E9E", lw=1.2, ls=":",
                        label=f"null max {_row['null_dbic_max']:,.0f}")
            axB.set_title(f"dBIC against a skew-matched unimodal null\n"
                          f"p = {_row['p_dbic_vs_null']:.4f}  -  "
                          f"{_row['verdict'].split(' — ')[0]}", fontsize=9)
            axB.set_xlabel("dBIC (1 component - 2 components)")
            axB.set_yticks([]); axB.legend(fontsize=7)
            axB.spines[["top", "right", "left"]].set_visible(False)
        plt.suptitle(f"Is the signature score bimodal? - {ARM_FOLDER}", fontsize=12)
        plt.tight_layout()
        fig2d.savefig(os.path.join(ARM_DIR, "S2d_mixture.png"), dpi=200,
                      bbox_inches="tight")
        plt.show(); plt.close(fig2d)
        print(f"  Saved S2d_mixture.png")


# %% SECTION 2e: Donor-bootstrap stability of the threshold and selected set
# Donors resampled with replacement; the threshold is recomputed on each
# resample and applied to the full arm. No clinical variable is used.
print("\n" + "=" * 74)
print("SECTION 2e: Donor-bootstrap stability of the cutoff")
print("=" * 74)
_t2e = time.time()

_donors_all = obs_full[donor_col].astype(str).values
_uniq_don   = np.unique(_donors_all)
_don_idx    = {d: np.where(_donors_all == d)[0] for d in _uniq_don}
print(f"  {len(_uniq_don)} donors, {len(_donors_all):,} cells, "
      f"{BOOT_B} bootstrap resamples")
print(f"  grid: {BOOT_GRID.min()*100:.1f}% to {BOOT_GRID.max()*100:.1f}% "
      f"in {len(BOOT_GRID)} steps")
print(f"  methods: {BOOT_METHODS}")

_boot_rng  = np.random.default_rng(BOOT_SEED)
_boot_rows = []

# Full-data reference sets and thresholds, once.
_ref_sets = {(_mk, float(_p)): _mask_at(_mk, float(_p))
             for _mk in BOOT_METHODS for _p in BOOT_GRID}

# Per-cell selection counters, one per (method, cutoff).
_sel_count = {k: np.zeros(len(_donors_all), dtype=np.int32) for k in _ref_sets}
_jac_sum   = {k: 0.0 for k in _ref_sets}
_npc_vals  = {k: [] for k in _ref_sets}     # implied top-% on the full arm

_base_needed = set()
for _mk in BOOT_METHODS:
    if _mk in INTERSECTION_PAIRS:
        _base_needed.update(INTERSECTION_PAIRS[_mk])
    else:
        _base_needed.add(_mk)
_base_needed = sorted(_base_needed)

for _b in range(BOOT_B):
    _pick = _boot_rng.choice(_uniq_don, size=len(_uniq_don), replace=True)
    _cells = np.concatenate([_don_idx[d] for d in _pick])
    _thr_b = {}
    for _bm in _base_needed:
        _sv = np.sort(np.asarray(score_cols[_bm], dtype=float)[_cells])
        _thr_b[_bm] = _sv          # sorted once; every cutoff is a lookup
    for _mk in BOOT_METHODS:
        for _p in BOOT_GRID:
            _p = float(_p)
            if _mk in INTERSECTION_PAIRS:
                _a, _bb = INTERSECTION_PAIRS[_mk]
                _ta = _thr_b[_a][int(np.floor((1 - _p) * (len(_thr_b[_a]) - 1)))]
                _tb = _thr_b[_bb][int(np.floor((1 - _p) * (len(_thr_b[_bb]) - 1)))]
                _sel = ((np.asarray(score_cols[_a], dtype=float) >= _ta) &
                        (np.asarray(score_cols[_bb], dtype=float) >= _tb))
            else:
                _t = _thr_b[_mk][int(np.floor((1 - _p) * (len(_thr_b[_mk]) - 1)))]
                _sel = np.asarray(score_cols[_mk], dtype=float) >= _t
            _key = (_mk, _p)
            _sel_count[_key] += _sel
            _ref = _ref_sets[_key]
            _inter = int(np.count_nonzero(_sel & _ref))
            _union = int(np.count_nonzero(_sel | _ref))
            _jac_sum[_key] += _inter / max(_union, 1)
            _npc_vals[_key].append(100.0 * _sel.sum() / len(_sel))
    if (_b + 1) % max(BOOT_B // 5, 1) == 0:
        print(f"    ... {_b+1}/{BOOT_B} resamples  [{time.time()-_t2e:.0f}s]")

for (_mk, _p), _cnt in _sel_count.items():
    _ref = _ref_sets[(_mk, _p)]
    _n_ref = int(_ref.sum())
    _core = int(np.count_nonzero(_cnt[_ref] >= BOOT_CORE_FRAC * BOOT_B))
    _npc = np.asarray(_npc_vals[(_mk, _p)], dtype=float)
    _boot_rows.append(dict(
        method=_mk, top_pct=round(_p * 100, 3),
        n_ref=_n_ref,
        implied_pct_mean=round(float(_npc.mean()), 4),
        implied_pct_sd=round(float(_npc.std(ddof=1)), 4),
        implied_pct_cv=round(float(_npc.std(ddof=1) / max(_npc.mean(), 1e-9)), 5),
        mean_jaccard_vs_full=round(_jac_sum[(_mk, _p)] / BOOT_B, 4),
        core_n=_core,
        core_frac=round(_core / max(_n_ref, 1), 4)))

boot_df = pd.DataFrame(_boot_rows).sort_values(["method", "top_pct"])
boot_df.to_csv(os.path.join(ARM_DIR, "S2e_bootstrap_stability.csv"), index=False)
print(f"\n  Saved S2e_bootstrap_stability.csv  [{time.time()-_t2e:.0f}s]")

print(f"\n  STABILITY vs CUTOFF  [{SCAN_METHOD}]")
_bm = boot_df[boot_df.method == SCAN_METHOD]
print(f"  {'cutoff':>7} {'n_ref':>8} {'mean Jaccard':>13} {'core frac':>10} "
      f"{'implied % SD':>13}")
print("  " + "-" * 58)
for _, _r in _bm.iterrows():
    _star = "  <- PRIMARY" if abs(_r.top_pct - PRIMARY_TOP_PCT * 100) < 1e-6 else ""
    print(f"  {_r.top_pct:>6.2f}% {int(_r.n_ref):>8,} {_r.mean_jaccard_vs_full:>13.4f} "
          f"{_r.core_frac:>10.4f} {_r.implied_pct_sd:>13.4f}{_star}")

if len(_bm):
    _best_j = _bm.loc[_bm.mean_jaccard_vs_full.idxmax()]
    _rng_j  = float(_bm.mean_jaccard_vs_full.max() - _bm.mean_jaccard_vs_full.min())
    print(f"\n  Jaccard is maximised at {_best_j.top_pct:.2f}% "
          f"({_best_j.mean_jaccard_vs_full:.4f}); total range across the grid "
          f"is {_rng_j:.4f}.")
    if _rng_j < BOOT_FLAT_TOL:
        print(f"  Range < {BOOT_FLAT_TOL}: stability is flat across the grid.")
    else:
        print(f"  Range exceeds {BOOT_FLAT_TOL}.")

fig2e, (axE1, axE2) = plt.subplots(1, 2, figsize=(15, 5.6))
for _mk in BOOT_METHODS:
    _d = boot_df[boot_df.method == _mk]
    _c = ARM_COLOR if _mk == SCAN_METHOD else "#9E9E9E"
    _lw = 2.4 if _mk == SCAN_METHOD else 1.4
    axE1.plot(_d.top_pct, _d.mean_jaccard_vs_full, color=_c, lw=_lw, label=_mk)
    axE2.plot(_d.top_pct, _d.core_frac, color=_c, lw=_lw, label=_mk)
for _ax, _lb in [(axE1, "mean Jaccard vs the full-data set"),
                 (axE2, f"fraction of selected cells kept in >= "
                        f"{BOOT_CORE_FRAC:.0%} of resamples")]:
    for _p in SWEEP_PCTS:
        _ax.axvline(_p * 100, color="black", lw=1.0, alpha=0.5,
                    ls=("-" if abs(_p - PRIMARY_TOP_PCT) < 1e-9 else ":"))
    _ax.set_xlabel("top-% cutoff"); _ax.set_ylabel(_lb)
    _ax.legend(fontsize=8); _ax.spines[["top", "right"]].set_visible(False)
axE1.set_title("Donor-bootstrap stability of the selected set\n"
               "a flat curve means no cutoff is more reproducible than another",
               fontsize=10)
axE2.set_title(f"Stable core, {BOOT_B} donor resamples\n"
               "solid vline = PRIMARY, dotted = other sweep points", fontsize=10)
plt.suptitle(f"Cutoff stability under donor resampling - {ARM_FOLDER}", fontsize=12)
plt.tight_layout()
fig2e.savefig(os.path.join(ARM_DIR, "S2e_bootstrap_stability.png"), dpi=200,
              bbox_inches="tight")
plt.show(); plt.close(fig2e)
print(f"  Saved S2e_bootstrap_stability.png")
del _sel_count; gc.collect()


# %% SECTION 2f: Arm matching - equal N, equal purity, or equal percentage
# Selected N and purity of each arm across the grid. Each arm writes its half
# to MATCH_SHARED_DIR; the second arm to run writes the joint table.
print("\n" + "=" * 74)
print("SECTION 2f: Arm matching - N, purity, or percentage")
print("=" * 74)
_t2f = time.time()

os.makedirs(MATCH_SHARED_DIR, exist_ok=True)
_enr_mask = np.zeros(len(obs_full), dtype=bool)
for _st in RULE_ENRICHED_STYPES:
    if _st in _st_masks:
        _enr_mask |= _st_masks[_st]
_tgt_mask = _st_masks.get(TARGET_STYPE, np.zeros(len(obs_full), dtype=bool))

_match_rows = []
for _p in SCAN_GRID:
    _m = _mask_at(SCAN_METHOD, float(_p))
    _n = int(_m.sum())
    _match_rows.append(dict(
        arm=ARM_TAG, top_pct=round(float(_p) * 100, 3), n_selected=_n,
        pct_of_arm=round(100.0 * _n / len(obs_full), 4),
        purity_enriched=round(100.0 * int((_m & _enr_mask).sum()) / max(_n, 1), 3),
        purity_target=round(100.0 * int((_m & _tgt_mask).sum()) / max(_n, 1), 3)))
match_half = pd.DataFrame(_match_rows)
_half_path = os.path.join(MATCH_SHARED_DIR, f"S2f_half_{ARM_TAG}.csv")
match_half.to_csv(_half_path, index=False)
print(f"  Wrote this arm's half: {_half_path}")

_pri = match_half.iloc[(match_half.top_pct - PRIMARY_TOP_PCT * 100).abs().idxmin()]
print(f"\n  THIS ARM at the primary {PRIMARY_TOP_PCT*100:.1f}%: "
      f"n = {int(_pri.n_selected):,}, purity(PVM_3+PVM_4) = "
      f"{_pri.purity_enriched:.1f}%, purity(target) = {_pri.purity_target:.1f}%")

_others_half = [f for f in os.listdir(MATCH_SHARED_DIR)
                if f.startswith("S2f_half_") and f != os.path.basename(_half_path)]
if not _others_half:
    print(f"\n  The other arm's half is not present yet.")
else:
    _oth = pd.read_csv(os.path.join(MATCH_SHARED_DIR, sorted(_others_half)[-1]))
    _oth_tag = str(_oth["arm"].iloc[0])
    print(f"\n  Found the other arm: {_oth_tag} ({len(_oth)} grid points)")

    def _closest(df, col, target):
        _i = (df[col] - target).abs().idxmin()
        return df.loc[_i]

    _oth_pri_pct = OTHER_ARM_PRIMARY_PCT
    _oth_pri = _closest(_oth, "top_pct", _oth_pri_pct)
    print(f"  {_oth_tag} at its primary {_oth_pri_pct:.1f}%: "
          f"n = {int(_oth_pri.n_selected):,}, "
          f"purity = {_oth_pri.purity_enriched:.1f}%")

    _joint = []
    for _lbl, _tgt_col, _tgt_val in [
            ("equal percentage (status quo)", None, None),
            ("match N to the other arm", "n_selected", float(_oth_pri.n_selected)),
            ("match purity to the other arm", "purity_enriched",
             float(_oth_pri.purity_enriched))]:
        if _tgt_col is None:
            _this = _pri; _that = _oth_pri
        else:
            _this = _closest(match_half, _tgt_col, _tgt_val); _that = _oth_pri
        _joint.append(dict(
            criterion=_lbl,
            this_arm=ARM_TAG, this_top_pct=float(_this.top_pct),
            this_n=int(_this.n_selected),
            this_purity=float(_this.purity_enriched),
            other_arm=_oth_tag, other_top_pct=float(_that.top_pct),
            other_n=int(_that.n_selected),
            other_purity=float(_that.purity_enriched)))
    joint_df = pd.DataFrame(_joint)
    joint_df.to_csv(os.path.join(ARM_DIR, "S2f_arm_matching.csv"), index=False)
    print(f"\n{joint_df.to_string(index=False)}")
    print(f"\n  Saved S2f_arm_matching.csv")
print(f"  [{time.time()-_t2f:.1f}s]")


# %% SECTION 2g: Is the M4 intersection doing anything?
# Spearman of the two parent scores, and how much of the intersection the top-k
# of one parent alone recovers (k matched to the intersection).
print("\n" + "=" * 74)
print("SECTION 2g: Scorer overlap - what the intersection buys")
print("=" * 74)
_t2g = time.time()

try:
    from scipy.stats import spearmanr as _spearman2g
except Exception:
    _spearman2g = None

_ovl_rows = []
for _mk in [m for m in ALL_METHODS if m in INTERSECTION_PAIRS]:
    _a, _b = INTERSECTION_PAIRS[_mk]
    _sa = np.asarray(score_cols[_a], dtype=float)
    _sb = np.asarray(score_cols[_b], dtype=float)
    _rho = (float(_spearman2g(_sa, _sb).statistic)
            if _spearman2g is not None else np.nan)
    for _p in sorted(set(list(SWEEP_PCTS) + [PRIMARY_TOP_PCT])):
        _p = float(_p)
        _mI = _mask_at(_mk, _p)
        _k  = int(_mI.sum())
        if _k == 0:
            continue
        # Top-k of each parent alone, k matched to the intersection's size.
        _topk_a = _sa >= np.sort(_sa)[-_k]
        _topk_b = _sb >= np.sort(_sb)[-_k]
        _cont_a = int(np.count_nonzero(_mI & _topk_a)) / _k
        _cont_b = int(np.count_nonzero(_mI & _topk_b)) / _k
        _nA = int((_sa >= np.percentile(_sa, (1 - _p) * 100)).sum())
        _ovl_rows.append(dict(
            intersection=_mk, parent_a=_a, parent_b=_b,
            spearman_a_b=round(_rho, 5),
            top_pct=round(_p * 100, 3),
            n_parent_a_at_pct=_nA, n_intersection=_k,
            retention_vs_parent_a=round(_k / max(_nA, 1), 4),
            effective_pct_of_arm=round(100.0 * _k / len(_sa), 4),
            contained_in_topk_A=round(_cont_a, 4),
            contained_in_topk_B=round(_cont_b, 4)))

ovl_df = pd.DataFrame(_ovl_rows)
if len(ovl_df):
    ovl_df.to_csv(os.path.join(ARM_DIR, "S2g_scorer_overlap.csv"), index=False)
    _foc = ovl_df[(ovl_df.intersection == SCAN_METHOD) &
                  (np.isclose(ovl_df.top_pct, PRIMARY_TOP_PCT * 100, atol=1e-6))]
    if len(_foc):
        _r = _foc.iloc[0]
        print(f"\n  {SCAN_METHOD} at the primary {PRIMARY_TOP_PCT*100:.1f}%")
        print(f"     spearman({_r.parent_a}, {_r.parent_b}) = {_r.spearman_a_b:+.4f}")
        print(f"     |{_r.parent_a} top| = {int(_r.n_parent_a_at_pct):,}  ->  "
              f"|intersection| = {int(_r.n_intersection):,}  "
              f"(retention {_r.retention_vs_parent_a:.1%})")
        print(f"     effective fraction of the arm = "
              f"{_r.effective_pct_of_arm:.3f}%  (nominal "
              f"{PRIMARY_TOP_PCT*100:.1f}%)")
        print(f"     contained in the top-k of {_r.parent_a} alone: "
              f"{_r.contained_in_topk_A:.1%}")
        print(f"     contained in the top-k of {_r.parent_b} alone: "
              f"{_r.contained_in_topk_B:.1%}")
        _mx = max(_r.contained_in_topk_A, _r.contained_in_topk_B)
        if _mx >= OVERLAP_REDUNDANT_AT:
            print(f"\n     VERDICT: >= {OVERLAP_REDUNDANT_AT:.0%} of the "
                  f"intersection is recoverable from one parent at matched k.")
        else:
            print(f"\n     VERDICT: the intersection is not recoverable from "
                  f"either parent alone at matched k.")
    print(f"\n  Saved S2g_scorer_overlap.csv  [{time.time()-_t2g:.1f}s]")

# %% SECTIONS 3-11: Cutoff sweep driver
# Each pass rebinds TOP_PCT and OUT_DIR ({ARM_DIR}/top_<pct>pct/). The primary
# cutoff runs last, so the namespace left at the end is the primary's.
print("\n" + "=" * 74)
print(f"CUTOFF SWEEP - {len(SWEEP_PCTS)} cutoffs, each into its own subfolder")
print(f"  order: " + ", ".join(f"{p*100:.1f}%" for p in SWEEP_PCTS)
      + f"   (PRIMARY = {PRIMARY_TOP_PCT*100:.1f}%, run last)")
print(f"  figures: make={SWEEP_MAKE_FIGS}  show={SWEEP_SHOW_FIGS}")
print("=" * 74)

_SWEEP_RESULTS = {}
_t_sweep = time.time()

for _i_sw, _sw in enumerate(SWEEP_PCTS, 1):
    TOP_PCT     = float(_sw)
    _IS_PRIMARY = abs(TOP_PCT - PRIMARY_TOP_PCT) < 1e-12
    OUT_DIR     = os.path.join(ARM_DIR, SWEEP_SUBDIR_FMT.format(TOP_PCT * 100))
    os.makedirs(OUT_DIR, exist_ok=True)
    _make_figs  = (SWEEP_MAKE_FIGS == "all"
                   or (SWEEP_MAKE_FIGS == "primary_only" and _IS_PRIMARY))
    _show_figs  = (SWEEP_SHOW_FIGS == "all"
                   or (SWEEP_SHOW_FIGS == "primary_only" and _IS_PRIMARY))
    print("\n" + "#" * 74)
    print(f"# SWEEP {_i_sw}/{len(SWEEP_PCTS)} - top {TOP_PCT*100:.1f}%"
          f"{'   [PRIMARY]' if _IS_PRIMARY else ''}")
    print(f"#   -> {os.path.abspath(OUT_DIR)}")
    print(f"#   figures: {'make+show' if (_make_figs and _show_figs) else ('make only' if _make_figs else 'skipped')}"
          f"   elapsed {time.time()-_t_sweep:.0f}s")
    print("#" * 74)

    # == Canonical M4 mask: top TOP_PCT of M1_FC AND of M2_Hybrid ===============
    # Later parts load this mask; they do not recompute it.
    _thr_M1_canon = np.percentile(score_cols["M1_FC"],     (1 - TOP_PCT) * 100)
    _thr_M2_canon = np.percentile(score_cols["M2_Hybrid"], (1 - TOP_PCT) * 100)
    M4_TOP_MASK   = (score_cols["M1_FC"]     >= _thr_M1_canon) & \
                    (score_cols["M2_Hybrid"] >= _thr_M2_canon)
    # Attach as a column so it survives DataFrame operations downstream
    obs_full["is_M4_top"] = M4_TOP_MASK
    N_M4_TOP      = int(M4_TOP_MASK.sum())
    M4_THR_M1     = float(_thr_M1_canon)   # expose for audit logging
    M4_THR_M2     = float(_thr_M2_canon)
    print(f"\n  [CANONICAL M4 MASK]")
    print(f"    M1 threshold : {M4_THR_M1:.6f}  (top {TOP_PCT*100:.1f}%)")
    print(f"    M2 threshold : {M4_THR_M2:.6f}  (top {TOP_PCT*100:.1f}%)")
    print(f"    N M4-top cells (M1 AND M2): {N_M4_TOP:,}  "
          f"({100*N_M4_TOP/n_cells:.2f}% of neocortical cells)")
    print(f"    Stored in namespace as: M4_TOP_MASK (bool array), "
          f"obs_full['is_M4_top'] (column)")

    # == Inter-scorer Pearson correlations (base methods) ========================
    print("\n  Inter-scorer Pearson correlations (base methods):")
    mkeys = list(SCORERS.keys())
    for i in range(len(mkeys)):
        for j in range(i+1, len(mkeys)):
            r = np.corrcoef(score_cols[mkeys[i]], score_cols[mkeys[j]])[0,1]
            print(f"    {mkeys[i]} vs {mkeys[j]} : r = {r:.4f}")

    # Save scores
    save_df = obs_full[[f"score_{m}" for m in ALL_METHODS]].copy()
    save_df.to_csv(os.path.join(OUT_DIR,"S0_cell_scores_all_methods.csv"))
    obs_full[["score_M1_FC"]].rename(
        columns={"score_M1_FC":"score_net"}).to_csv(
        os.path.join(OUT_DIR,"S0_cell_scores.csv"))
    print(f"\n  Saved S0_cell_scores_all_methods.csv + S0_cell_scores.csv")


    # ===== SECTION 3: Primary analysis [M1], enrichment
    np.random.seed(42); random.seed(42)
    PRIMARY      = "M1_FC"
    score_net    = score_cols[PRIMARY]
    top_thresh   = np.percentile(score_net, (1-TOP_PCT)*100)
    top_mask     = score_net >= top_thresh
    n_top        = int(top_mask.sum())
    n_total      = len(obs_full)
    obs_full["is_top"] = top_mask

    print(f"\nSECTION 3: Top-{TOP_PCT*100:.1f}% enrichment [{PRIMARY}]")
    print(f"  Threshold: {top_thresh:.4f}  n_top: {n_top:,}")

    stype_all = obs_full[stype_col].value_counts()
    stype_top = obs_full[obs_full["is_top"]][stype_col].value_counts()

    enrich_rows = []
    N_BOOT_ENR = 2000
    for st in stype_all.index:
        n_st_all = int(stype_all.get(st,0))
        n_st_top = int(stype_top.get(st,0))
        enr = (n_st_top/max(n_st_all,1))/TOP_PCT

        # Donor-level n
        _st_top_msk = top_mask & (obs_full[stype_col]==st).values
        n_st_don    = int(obs_full.loc[_st_top_msk, donor_col].astype(str).nunique())

        # Donor-bootstrap 95% CI on the enrichment fold.
        _st_mask_all = (obs_full[stype_col]==st).values
        _per_d_all = obs_full[_st_mask_all].groupby(donor_col).size()
        _per_d_top = obs_full[_st_mask_all & top_mask].groupby(donor_col).size()
        _per_d = pd.DataFrame({"n_all": _per_d_all, "n_top": _per_d_top}).fillna(0)
        _boot_folds = np.empty(N_BOOT_ENR)
        for _bi in range(N_BOOT_ENR):
            _bidx = np.random.choice(len(_per_d), len(_per_d), replace=True)
            _b_all = _per_d["n_all"].iloc[_bidx].sum()
            _b_top = _per_d["n_top"].iloc[_bidx].sum()
            _boot_folds[_bi] = (_b_top / max(_b_all, 1)) / TOP_PCT
        ci_lo, ci_hi = np.percentile(_boot_folds, [2.5, 97.5])

        enrich_rows.append({"supertype":st,"n_total":n_st_all,"n_top":n_st_top,
                            "n_donors_top":n_st_don,
                            "pct_top":round(100*n_st_top/max(n_st_all,1),2),
                            "enrichment_fold":round(enr,3),
                            "boot_ci_lo":round(ci_lo,3),
                            "boot_ci_hi":round(ci_hi,3)})

    enrich_df = pd.DataFrame(enrich_rows).sort_values("enrichment_fold",ascending=False)
    enrich_df.to_csv(os.path.join(OUT_DIR,"S1_supertype_enrichment.csv"),index=False)

    print(f"\n  {'Supertype':<40} {'n_all':>7} {'n_top':>6} {'n_don':>6} "
          f"{'%top':>6} {'fold':>6} {'95% CI':>16}")
    print("  "+"-"*95)
    for _,r in enrich_df.iterrows():
        print(f"  {r.supertype:<40} {r.n_total:>7,} {r.n_top:>6,} "
              f"{r.n_donors_top:>6,} "
              f"{r.pct_top:>5.1f}% {r.enrichment_fold:>6.2f}x "
              f"[{r.boot_ci_lo:.2f}, {r.boot_ci_hi:.2f}]")
    print("  " + "-"*95)
    print("  Enrichment fold = pct_top / arm top-%; 95% CI from 2,000 donor resamples.")


    # ===== SECTION 4: ADNC-stratified proportions [M1]
    print("\nSECTION 4: ADNC-stratified proportions [M1]")
    np.random.seed(42); random.seed(42)

    all_donor_totals  = obs_full.groupby(donor_col).size().rename("n_all")
    donor_adnc_lookup = (obs_full[[donor_col,adnc_col]]
                         .drop_duplicates().set_index(donor_col))

    adnc_results = {}
    for st in stype_all.index:
        st_top = top_mask & (obs_full[stype_col]==st).values
        if int(st_top.sum()) < 20: continue
        per_d = obs_full[st_top].groupby(donor_col).size().rename("n_top_st")
        prop  = pd.concat([all_donor_totals,per_d],axis=1).fillna(0)
        prop["pct"]       = 100*prop["n_top_st"]/prop["n_all"]
        prop[adnc_col]    = donor_adnc_lookup[adnc_col]
        prop["adnc_rank"] = to_num(prop[adnc_col], adnc_rank_map)
        _s4_valid = prop["adnc_rank"].notna() & prop["pct"].notna()
        if _s4_valid.sum() >= 3:
            rho,p_rho = spearmanr(prop.loc[_s4_valid,"adnc_rank"],
                                   prop.loc[_s4_valid,"pct"])
        else:
            rho,p_rho = np.nan,np.nan
        groups = [prop[prop[adnc_col]==a]["pct"].values for a in adnc_order]
        groups_valid = [g for g in groups if len(g)>0]
        _,kw_p = kruskal(*groups_valid) if len(groups_valid)>=2 else (np.nan,np.nan)
        # Use last element of adnc_order as "High" and first as "Not AD"
        hi = prop[prop[adnc_col]==adnc_order[-1]]["pct"].values
        no = prop[prop[adnc_col]==adnc_order[0]]["pct"].values
        _,mw_p = (mannwhitneyu(hi,no,alternative="greater")
                  if len(hi)>0 and len(no)>0 else (np.nan,1.0))
        adnc_results[st] = {
            "adnc_rho":round(rho,4),"adnc_p":round(p_rho,6),
            "kw_p":round(kw_p,6),"mw_p":round(mw_p,6),
            "n_top":int(st_top.sum()),
            "means":{a:round(prop[prop[adnc_col]==a]["pct"].mean(),3)
                     for a in adnc_order},
            "prop_df":prop,
        }

    adnc_summary = []
    for st,r in sorted(adnc_results.items(),key=lambda x:-x[1]["adnc_rho"]):
        fold = float(enrich_df.loc[enrich_df.supertype==st,"enrichment_fold"].values[0]) \
               if st in enrich_df.supertype.values else np.nan
        adnc_summary.append({
            "supertype":st,"n_top":r["n_top"],"enrichment_fold":fold,
            "ADNC_rho":r["adnc_rho"],"ADNC_p":r["adnc_p"],
            "KW_p":r["kw_p"],"MW_High_vs_NotAD":r["mw_p"],
            **{f"mean_{a.replace(' ','_')}":r["means"].get(a,np.nan) for a in adnc_order},
        })
    adnc_df = pd.DataFrame(adnc_summary)
    adnc_df.to_csv(os.path.join(OUT_DIR,"S2_ADNC_proportions_per_supertype.csv"),index=False)
    print("  Saved S2_ADNC_proportions_per_supertype.csv")


    # ===== SECTION 5: Donor-level staging [M1]
    print("\nSECTION 5: Donor-level disease staging [M1]")

    per_donor_top   = obs_full[top_mask].groupby(donor_col).size().rename("n_top")
    per_donor_total = obs_full.groupby(donor_col).size().rename("n_total")
    donor_prop = pd.concat([per_donor_total,per_donor_top],axis=1).fillna(0)
    donor_prop["pct_top"]   = 100*donor_prop["n_top"]/donor_prop["n_total"]
    donor_prop[adnc_col]    = donor_adnc_lookup[adnc_col]
    donor_prop["adnc_rank"] = to_num(donor_prop[adnc_col], adnc_rank_map)

    if braak_col in obs_full.columns:
        dbr = obs_full[[donor_col,braak_col]].drop_duplicates().set_index(donor_col)
        donor_prop[braak_col]   = dbr[braak_col]
        donor_prop["braak_num"] = to_num(donor_prop[braak_col], BRAAK_MAP)
        unmapped = donor_prop["braak_num"].isna().sum()
        if unmapped>0:
            print(f"  WARNING: {unmapped} donors unmapped Braak: "
                  f"{donor_prop.loc[donor_prop.braak_num.isna(),braak_col].unique()}")
        valid_b = ~donor_prop["braak_num"].isna()
        n_bv    = int(valid_b.sum())
        rho_braak,p_braak = (spearmanr(donor_prop.loc[valid_b,"braak_num"],
                                        donor_prop.loc[valid_b,"pct_top"])
                             if n_bv>=3 else (np.nan,np.nan))
    else:
        donor_prop["braak_num"] = np.nan
        rho_braak,p_braak,n_bv = np.nan,np.nan,0

    _adnc_valid = donor_prop["adnc_rank"].notna() & donor_prop["pct_top"].notna()
    _n_adnc_valid = int(_adnc_valid.sum())
    if _n_adnc_valid < 3:
        print(f"  WARNING: only {_n_adnc_valid} donors have valid adnc_rank.")
        print(f"  Unique ADNC values in donor_prop: "
              f"{donor_prop[adnc_col].unique().tolist()}")
        rho_adnc, p_adnc = np.nan, np.nan
    else:
        rho_adnc,p_adnc = spearmanr(
            donor_prop.loc[_adnc_valid,"adnc_rank"],
            donor_prop.loc[_adnc_valid,"pct_top"])
    print(f"  ADNC  rho = {rho_adnc:+.4f}  p = {p_adnc:.6f}  (n={_n_adnc_valid})")
    print(f"  Braak rho = {rho_braak:+.4f}  p = {p_braak:.6f}  (n={n_bv})")
    for a in adnc_order:
        v = donor_prop[donor_prop[adnc_col]==a]["pct_top"].values
        print(f"  {a:<15}: n={len(v):2d}  mean={v.mean():.3f}%")
    donor_prop.to_csv(os.path.join(OUT_DIR,"S3_donor_level_top_pct.csv"))
    print("  Saved S3_donor_level_top_pct.csv")


    # ===== SECTION 6: Figures - all six scorers (M1-M6)
    print("\nSECTION 6: Figures [M1-M6]")
    # Seed reset for reproducible jitter
    np.random.seed(42)
    random.seed(42)

    supertypes_all = sorted(obs_full[stype_col].unique())
    n_st = len(supertypes_all)
    stype_cmap = {st:plt.cm.tab20(i/max(n_st,1)) for i,st in enumerate(supertypes_all)}

    def _scorer_top(mkey):
        """Return (top_mask, threshold, score_array). M4_M1andM2 uses canonical mask."""
        if mkey == "M4_M1andM2":
            return M4_TOP_MASK, min(M4_THR_M1, M4_THR_M2), score_cols[mkey]
        if mkey in INTERSECTION_PAIRS:
            _pa, _pb = INTERSECTION_PAIRS[mkey]
            _sc_a = score_cols[_pa];  _thr_a = np.percentile(_sc_a, (1-TOP_PCT)*100)
            _sc_b = score_cols[_pb];  _thr_b = np.percentile(_sc_b, (1-TOP_PCT)*100)
            _mask = (_sc_a >= _thr_a) & (_sc_b >= _thr_b)
            return _mask, min(_thr_a, _thr_b), score_cols[mkey]
        sc  = score_cols[mkey]
        thr = np.percentile(sc, (1-TOP_PCT)*100)
        return sc >= thr, thr, sc

    SCORER_FIG_TAG = {"M1_FC": "M1", "M2_Hybrid": "M2", "M3_pRank": "M3",
                       "M4_M1andM2": "M4", "M5_M1andM3": "M5", "M6_M2andM3": "M6"}

    # Region label for figure titles / axis labels
    REGION_LABEL = "neocortical (DFC/MTG/AnG/STG/V1C/ITG)"

    for _mkey in (ALL_METHODS if _make_figs else []):   # sweep gate
        np.random.seed(42); random.seed(42)  # reset per scorer iteration
        _tag = SCORER_FIG_TAG[_mkey]
        _sc, _thr, _sc_arr = _scorer_top(_mkey)
        _n_top_m = int(_sc.sum())

        # Enrichment recomputed for this scorer (Section 3's enrich_df is M1).
        _enr_rows_m = []
        for _stm in stype_all.index:
            _na = int(stype_all.get(_stm, 0))
            _nt = int((_sc & (obs_full[stype_col]==_stm).values).sum())
            _tabm = np.array([[_nt, _n_top_m-_nt],
                              [_na-_nt, len(obs_full)-_na-(_n_top_m-_nt)]])
            try:    _, _pm, *_ = chi2_contingency(_tabm, correction=False)
            except Exception: _pm = 1.0
            _enr_rows_m.append({"supertype":_stm, "n_total":_na, "n_top":_nt,
                                "pct_top":round(100*_nt/max(_na,1),2),
                                "enrichment_fold":round((_nt/max(_na,1))/TOP_PCT,3),
                                "chi2_p":round(_pm,8)})
        _enr_df_m = (pd.DataFrame(_enr_rows_m)
                       .sort_values("enrichment_fold", ascending=False))

        _per_top   = obs_full[_sc].groupby(donor_col).size().rename("n_top")
        _per_total = obs_full.groupby(donor_col).size().rename("n_total")
        _dp = pd.concat([_per_total, _per_top], axis=1).fillna(0)
        _dp["pct_top"]   = 100 * _dp["n_top"] / _dp["n_total"]
        _dp[adnc_col]    = donor_adnc_lookup[adnc_col]
        _dp["adnc_rank"] = to_num(_dp[adnc_col], adnc_rank_map)
        _dp[braak_col]   = donor_prop[braak_col]
        _dp["braak_num"] = donor_prop["braak_num"]
        _rho_a, _p_a = spearmanr(_dp["adnc_rank"], _dp["pct_top"])
        _vb = ~_dp["braak_num"].isna()
        _rho_b, _p_b = (spearmanr(_dp.loc[_vb,"braak_num"], _dp.loc[_vb,"pct_top"])
                        if _vb.sum()>=3 else (np.nan, np.nan))
        _nvb = int(_vb.sum())

        # Fig 1: UMAP
        if has_umap:
            np.random.seed(42); random.seed(42)
            fig,axes = plt.subplots(1,3,figsize=(24,8))
            vals=_sc_arr; vmin,vmax=np.percentile(vals,[2,98])
            sc_p=axes[0].scatter(umap_coords[:,0],umap_coords[:,1],c=vals,cmap="RdYlBu_r",
                                 vmin=vmin,vmax=vmax,s=0.8,alpha=0.6,rasterized=True)
            plt.colorbar(sc_p,ax=axes[0],fraction=0.046,pad=0.04)
            axes[0].set_title(f"Mic.13 net score [{_tag}]\n({n_total:,} {REGION_LABEL})")
            axes[1].scatter(umap_coords[~_sc,0],umap_coords[~_sc,1],
                            c="#DDDDDD",s=0.6,alpha=0.4,rasterized=True)
            axes[1].scatter(umap_coords[_sc,0],umap_coords[_sc,1],
                            c="#D73027",s=1.5,alpha=0.8,rasterized=True)
            axes[1].set_title(f"Top {TOP_PCT*100:.1f}% [{_tag}]  n={_n_top_m:,}")
            for st in supertypes_all:
                msk=(obs_full[stype_col]==st).values
                axes[2].scatter(umap_coords[msk,0],umap_coords[msk,1],
                                c=[stype_cmap[st]],s=0.6,alpha=0.5,rasterized=True)
            axes[2].legend(handles=[mpatches.Patch(color=stype_cmap[st],
                           label=f"{st} (n={stype_all[st]:,})")
                           for st in supertypes_all],fontsize=6,markerscale=4,loc="best")
            axes[2].set_title("Supertype overlay")
            plt.suptitle(f"SEA-AD {REGION_LABEL} microglia - Mic.13 [{_tag}]",fontsize=11)
            plt.tight_layout()
            _fn = f"Fig1_UMAP_net_score_{_tag}.png"
            fig.savefig(os.path.join(OUT_DIR,_fn),dpi=150,bbox_inches="tight")
            _maybe_show()
            plt.close(fig); print(f"  Saved {_fn}")

        # Fig 2: Quantitative supertype breakdown
        _st_order_f2 = _enr_df_m.sort_values("enrichment_fold", ascending=False)["supertype"].tolist()
        _n_cells_f2   = [int((_sc & (obs_full[stype_col]==st).values).sum()) for st in _st_order_f2]
        _n_total_f2   = [int((obs_full[stype_col]==st).sum()) for st in _st_order_f2]
        _pct_cells_f2 = [100*nc/max(nt,1) for nc,nt in zip(_n_cells_f2,_n_total_f2)]
        _st_lbl_f2 = [s.replace("Micro-PVM_","PVM_").replace("-SEAAD","*") for s in _st_order_f2]
        _ef_map_f2 = dict(zip(_enr_df_m["supertype"], _enr_df_m["enrichment_fold"]))
        _bar_col_f2 = [
            "#D73027" if _ef_map_f2.get(st,0)>2 else
            "#FC8D59" if _ef_map_f2.get(st,0)>1.5 else
            "#FDAE61" if _ef_map_f2.get(st,0)>1 else "#AAAAAA"
            for st in _st_order_f2
        ]

        np.random.seed(42); random.seed(42)
        fig_f2, (ax_n_f2, ax_p_f2) = plt.subplots(
            1, 2, figsize=(max(14, len(_st_order_f2)*1.9), 6))

        _ymax_n = max(_n_cells_f2) if _n_cells_f2 else 1
        bars_n_f2 = ax_n_f2.bar(range(len(_st_order_f2)), _n_cells_f2,
                                 color=_bar_col_f2, edgecolor="black", lw=0.5)
        for _bv, _nc in zip(bars_n_f2, _n_cells_f2):
            ax_n_f2.text(_bv.get_x()+_bv.get_width()/2,
                         _bv.get_height()+_ymax_n*0.015,
                         f"{_nc:,}", ha="center", va="bottom", fontsize=7.5)
        ax_n_f2.set_xticks(range(len(_st_order_f2)))
        ax_n_f2.set_xticklabels(_st_lbl_f2, rotation=35, ha="right", fontsize=9)
        ax_n_f2.set_ylabel(f"Number of top-{TOP_PCT*100:.1f}% scoring cells", fontsize=10)
        ax_n_f2.set_title(
            f"Cell counts per supertype  [{_tag}]\n"
            f"top-{TOP_PCT*100:.1f}% Mic.13 scorers - {REGION_LABEL}; sorted by fold-enrichment",
            fontsize=10)
        ax_n_f2.set_xlim(-0.6, len(_st_order_f2)-0.4)
        ax_n_f2.spines[["top","right"]].set_visible(False)
        ax_n_f2.tick_params(axis="x", length=3)

        _ymax_p = max(_pct_cells_f2) if _pct_cells_f2 else 1
        bars_p_f2 = ax_p_f2.bar(range(len(_st_order_f2)), _pct_cells_f2,
                                 color=_bar_col_f2, edgecolor="black", lw=0.5)
        ax_p_f2.axhline(TOP_PCT*100, color="black", lw=1.4, ls="--",
                        label=f"Expected {TOP_PCT*100:.1f}% (uniform null)")
        for _bv, _pc in zip(bars_p_f2, _pct_cells_f2):
            ax_p_f2.text(_bv.get_x()+_bv.get_width()/2,
                         _bv.get_height()+_ymax_p*0.015,
                         f"{_pc:.1f}%", ha="center", va="bottom", fontsize=7.5)
        ax_p_f2.set_xticks(range(len(_st_order_f2)))
        ax_p_f2.set_xticklabels(_st_lbl_f2, rotation=35, ha="right", fontsize=9)
        ax_p_f2.set_ylabel(
            f"% of supertype cells in top-{TOP_PCT*100:.1f}% Mic.13 scorers",
            fontsize=10)
        ax_p_f2.set_title(
            f"Percentage per supertype  [{_tag}]\n"
            f"top-{TOP_PCT*100:.1f}% Mic.13 scorers - {REGION_LABEL}; sorted by fold-enrichment",
            fontsize=10)
        ax_p_f2.set_xlim(-0.6, len(_st_order_f2)-0.4)
        ax_p_f2.spines[["top","right"]].set_visible(False)
        ax_p_f2.tick_params(axis="x", length=3)
        ax_p_f2.legend(fontsize=9)

        from matplotlib.patches import Patch as _MPatch
        _leg_els = [
            _MPatch(facecolor="#D73027",edgecolor="black",label=">2x enriched"),
            _MPatch(facecolor="#FC8D59",edgecolor="black",label="1.5-2x enriched"),
            _MPatch(facecolor="#FDAE61",edgecolor="black",label="1-1.5x enriched"),
            _MPatch(facecolor="#AAAAAA",edgecolor="black",label="<1x (depleted)"),
        ]
        ax_n_f2.legend(handles=_leg_els, fontsize=8, loc="upper right",
                       framealpha=0.85, title="Fold-enrichment tier",
                       title_fontsize=8)

        plt.suptitle(
            f"Quantitative supertype breakdown - Mic.13 top-{TOP_PCT*100:.1f}% scorers  [{_tag}]\n"
            f"{REGION_LABEL}  |  Colour = fold-enrichment tier  |  *SEAAD-specific supertypes",
            fontsize=11)
        plt.tight_layout()
        _fn_f2 = f"Fig2_supertype_counts_pct_{_tag}.png"
        fig_f2.savefig(os.path.join(OUT_DIR,_fn_f2),dpi=150,bbox_inches="tight")
        _maybe_show()
        plt.close(fig_f2); print(f"  Saved {_fn_f2}")

        # Fig 3: enrichment fold bar
        np.random.seed(42); random.seed(42)
        fig,ax=plt.subplots(figsize=(14,6))
        se=_enr_df_m.sort_values("enrichment_fold",ascending=True)
        cb=["#D73027" if f>2 else "#FC8D59" if f>1.5 else "#FDAE61" if f>1 else "#AAAAAA"
            for f in se.enrichment_fold]
        bars=ax.barh(range(len(se)),se.enrichment_fold,color=cb,edgecolor="black",lw=0.4)
        ax.axvline(1,color="black",lw=1.5,ls="--",label="Expected 1x")
        ax.axvline(2,color="red",lw=1,ls="--",alpha=0.6,label="2x enrichment")
        ax.set_yticks(range(len(se)))
        ax.set_yticklabels([f"{r.supertype} (n={r.n_total:,})" for _,r in se.iterrows()],fontsize=8)
        ax.set_xlabel(f"Fold-enrichment [{_tag}, top-{TOP_PCT*100:.1f}%] - {REGION_LABEL}")
        ax.legend(fontsize=9)
        for bar,(_2,r) in zip(bars,se.iterrows()):
            sig="***" if r.chi2_p<0.001 else "**" if r.chi2_p<0.01 else "*" if r.chi2_p<0.05 else ""
            ax.text(bar.get_width()+0.02,bar.get_y()+bar.get_height()/2,
                    f"{r.enrichment_fold:.2f}x {sig}",va="center",fontsize=7.5)
        plt.tight_layout()
        _fn = f"Fig3_supertype_enrichment_{_tag}.png"
        fig.savefig(os.path.join(OUT_DIR,_fn),dpi=150,bbox_inches="tight")
        _maybe_show()
        plt.close(fig); print(f"  Saved {_fn}")

        # Fig 4: ADNC top-4
        np.random.seed(42); random.seed(42)
        top4_st=adnc_df.nlargest(4,"ADNC_rho").supertype.tolist()
        fig,axes_f3=plt.subplots(1,len(top4_st),figsize=(6*len(top4_st),6))
        if len(top4_st)==1: axes_f3=[axes_f3]
        for ax3,st in zip(axes_f3,top4_st):
            st_top_m = _sc & (obs_full[stype_col]==st).values
            per_d_m = obs_full[st_top_m].groupby(donor_col).size().rename("n_top_st")
            prop_m  = pd.concat([all_donor_totals,per_d_m],axis=1).fillna(0)
            prop_m["pct"]       = 100*prop_m["n_top_st"]/prop_m["n_all"]
            prop_m[adnc_col]    = donor_adnc_lookup[adnc_col]
            prop_m["adnc_rank"] = to_num(prop_m[adnc_col], adnc_rank_map)
            _f4_valid = prop_m["adnc_rank"].notna() & prop_m["pct"].notna()
            rho_m,p_m = (spearmanr(prop_m.loc[_f4_valid,"adnc_rank"],
                                    prop_m.loc[_f4_valid,"pct"])
                         if _f4_valid.sum()>=3 else (np.nan,np.nan))
            # Per-supertype jitter seed (crc32, stable across sessions).
            _jit_seed_f4 = stable_seed(st)
            for _ai, a in enumerate(adnc_order):
                sub=prop_m[prop_m[adnc_col]==a]
                xi=adnc_rank_map[a]
                jit=np.random.default_rng(_jit_seed_f4 + _ai).uniform(-0.15,0.15,len(sub))
                ax3.scatter([xi+j for j in jit],sub["pct"],c=adnc_pal[a],s=30,alpha=0.7,
                           edgecolors="black",lw=0.3,zorder=5)
                m_=sub["pct"].mean(); s_=sub["pct"].std()/max(np.sqrt(len(sub)),1)
                ax3.errorbar(xi,m_,yerr=s_,fmt="D",color=adnc_pal[a],markersize=9,
                            capsize=4,capthick=2,zorder=10,markeredgecolor="black",markeredgewidth=1)
            # NaN-guard: drop donors without a valid ADNC mapping before polyfit
            _pm_valid = prop_m[prop_m[adnc_col].isin(adnc_rank_map)].copy()
            x_r=to_num(_pm_valid[adnc_col], adnc_rank_map).values
            if len(x_r) >= 2:
                m_r,b_r=np.polyfit(x_r,_pm_valid["pct"].values,1)
            else:
                m_r,b_r = 0.0, _pm_valid["pct"].mean() if len(_pm_valid)>0 else 0.0
            xr=np.linspace(-0.3, len(adnc_order)-0.7, 50)
            ax3.plot(xr,m_r*xr+b_r,"k--",lw=1.5,alpha=0.6)
            ax3.set_xticks(range(len(adnc_order)))
            ax3.set_xticklabels(adnc_order,fontsize=9)
            ax3.set_ylabel(f"% {REGION_LABEL} microglia per donor")
            fold_v=float(_enr_df_m.loc[_enr_df_m.supertype==st,"enrichment_fold"].values[0]) \
                   if st in _enr_df_m.supertype.values else np.nan
            ax3.set_title(f"{st}\nADNC rho={rho_m:+.3f}  p={p_m:.4f}  fold={fold_v:.2f}x",fontsize=9)
            ax3.legend(handles=[mpatches.Patch(color=adnc_pal[a],label=a) for a in adnc_order],fontsize=8)
        plt.suptitle(f"Top-{TOP_PCT*100:.1f}% Mic.13 [{_tag}] - ADNC proportions ({REGION_LABEL})",fontsize=11)
        plt.tight_layout()
        _fn = f"Fig4_ADNC_top4_supertypes_{_tag}.png"
        fig.savefig(os.path.join(OUT_DIR,_fn),dpi=150,bbox_inches="tight")
        _maybe_show()
        plt.close(fig); print(f"  Saved {_fn}")




    # ===== SECTION 7: Summary CSV
    summary = enrich_df.merge(adnc_df[["supertype","ADNC_rho","ADNC_p","MW_High_vs_NotAD"]],
                               on="supertype",how="left")
    summary["combined"] = summary.enrichment_fold.fillna(0)+summary.ADNC_rho.fillna(-1)
    summary.sort_values("combined",ascending=False).to_csv(
        os.path.join(OUT_DIR,"S4_combined_summary.csv"),index=False)
    print("\nSECTION 7: Saved S4_combined_summary.csv")


    # ===== SECTION 8: Clinical correlation suite - all methods x all supertypes
    print("\n" + "="*70)
    print("SECTION 8: Clinical correlation suite - all 6 methods (M1-M6)")
    print("="*70)
    np.random.seed(42); random.seed(42)

    # == ADNC check ==============================================================
    _obs_adnc_vals = obs_full[adnc_col].astype(str).unique().tolist()
    _mapped_count  = sum(obs_full[adnc_col].isin(adnc_order))
    _unmapped_adnc = [v for v in _obs_adnc_vals if v not in adnc_order]
    print(f"  ADNC values in obs_full: {sorted(_obs_adnc_vals)}")
    print(f"  Values mapped to adnc_order: {_mapped_count:,} / {len(obs_full):,} cells")
    if _unmapped_adnc:
        print(f"  WARNING: {len(_unmapped_adnc)} unmapped ADNC values: {_unmapped_adnc}")
    else:
        print("  All ADNC values successfully mapped to adnc_order ✓")

    # == donor_clinical, including tissue-quality covariates for Part 10h ========
    _dc_cols = [c for c in [adnc_col, braak_col, cerad_col, thal_col, cog_col,
                            apoe_col, mmse_col, age_col, caa_col,
                            pmi_col, rin_col, ph_col, weight_col, severe_col]
                if c is not None and c in obs_full.columns]

    # Donor-constancy check: drop_duplicates keeps one row per donor.
    _nonconst = []
    for _c in _dc_cols:
        try:
            _nu = obs_full.groupby(donor_col, observed=True)[_c].nunique()
            _nv = int((_nu > 1).sum())
        except Exception:
            _nv = -1
        if _nv != 0:
            _nonconst.append((_c, _nv))
    if _nonconst:
        print("\n  WARNING - not constant within donor:")
        for _c, _n in _nonconst:
            print(f"          {_c:<42}: varies in {_n} donors")
    else:
        print("\n  all donor_clinical columns verified constant within donor.")

    donor_clinical = (obs_full[[donor_col] + _dc_cols]
                      .drop_duplicates(subset=[donor_col])
                      .set_index(donor_col))

    # Ordinals through to_num().
    donor_clinical["adnc_rank"]    = to_num(donor_clinical[adnc_col], adnc_rank_map)
    donor_clinical["braak_num"]    = to_num(donor_clinical[braak_col], BRAAK_MAP)
    donor_clinical["cerad_num"]    = (to_num(donor_clinical[cerad_col], CERAD_MAP)
                                      if cerad_col in donor_clinical.columns else np.nan)
    donor_clinical["thal_num"]     = (to_num(donor_clinical[thal_col], THAL_MAP)
                                      if thal_col in donor_clinical.columns else np.nan)
    donor_clinical["caa_num"]      = (to_num(donor_clinical[caa_col], CAA_MAP)
                                      if caa_col in donor_clinical.columns else np.nan)
    donor_clinical["apoe_e4_dose"] = (donor_clinical[apoe_col].astype(str).apply(
                                          lambda g: g.replace("/","").count("4"))
                                      if apoe_col in donor_clinical.columns else np.nan)
    donor_clinical["mmse_num"]     = (to_num(donor_clinical[mmse_col])
                                      if mmse_col in donor_clinical.columns else np.nan)
    donor_clinical["age_num"]      = (to_num(donor_clinical[age_col])
                                      if age_col in donor_clinical.columns else np.nan)

    # Tier-2 tissue-quality covariates, numeric
    for _tgt, _src in [("pmi_num", pmi_col), ("rin_num", rin_col),
                       ("ph_num", ph_col),   ("weight_num", weight_col)]:
        donor_clinical[_tgt] = (to_num(donor_clinical[_src])
                                if _src is not None and _src in donor_clinical.columns
                                else np.nan)
    if severe_col is not None and severe_col in donor_clinical.columns:
        donor_clinical["severe_flag"] = (
            donor_clinical[severe_col].astype(str).str.strip().str.lower()
            .isin({"y", "yes", "true", "1"}).astype(int))
    else:
        donor_clinical["severe_flag"] = np.nan

    print(f"\n  Donor clinical availability ({len(donor_clinical)} donors):")
    for nm,ck in [("ADNC rank","adnc_rank"),("Braak","braak_num"),("CERAD","cerad_num"),
                  ("Thal","thal_num"),("CAA","caa_num"),("APOE e4","apoe_e4_dose"),
                  ("MMSE","mmse_num"),("Age","age_num"),
                  ("PMI","pmi_num"),("RIN","rin_num"),("Brain pH","ph_num"),
                  ("Fresh Brain Wt","weight_num"),("Severe flag","severe_flag")]:
        n = int(pd.Series(donor_clinical[ck]).notna().sum())
        print(f"    {nm:<15}: {n:>3}/{len(donor_clinical)}")
    _t1 = ["age_num", "pmi_num"]
    _t2 = _t1 + ["rin_num", "ph_num", "weight_num"]
    print(f"    complete cases  Tier1 (Age,PMI)           : "
          f"{int(donor_clinical[_t1].notna().all(axis=1).sum())}/{len(donor_clinical)}")
    print(f"    complete cases  Tier1+RIN/pH/Weight       : "
          f"{int(donor_clinical[_t2].notna().all(axis=1).sum())}/{len(donor_clinical)}")
    donor_clinical.to_csv(os.path.join(OUT_DIR, "S5b_donor_clinical.csv"))
    print(f"    Saved S5b_donor_clinical.csv "
          f"({len(donor_clinical)} donors x {donor_clinical.shape[1]} cols)")

    def _run_clinical_suite(obs_df, mkey, eligible_st, all_donor_totals, donor_clinical,
                            top_mask_override=None):
        if top_mask_override is not None:
            _global_top = top_mask_override
        else:
            score_col = f"score_{mkey}"
            sc_vals   = obs_df[score_col].values
            thr       = np.percentile(sc_vals, (1-TOP_PCT)*100)
            _global_top = sc_vals >= thr

        rows = []
        for st in eligible_st:
            st_top = _global_top & (obs_df[stype_col]==st).values
            if int(st_top.sum()) < 20: continue
            per_d = obs_df[st_top].groupby(donor_col).size().rename("n_top_st")
            prop  = pd.concat([all_donor_totals,per_d],axis=1).fillna(0)
            prop["pct"]     = 100*prop["n_top_st"]/prop["n_all"]
            prop            = prop.join(donor_clinical,how="left")

            for label,ck,test_type,expected,green_finding in CLINICAL_TESTS:
                valid = prop[ck].notna() & prop["pct"].notna()
                n_v   = int(valid.sum())
                base  = {"supertype":st,"variable":label,"test":test_type,
                         "rho_or_stat":np.nan,"p_raw":np.nan,"n_donors":n_v,
                         "green_expected":expected,"green_finding":green_finding,
                         "note":""}
                if n_v < 5:
                    base["note"] = f"only {n_v} donors"
                    rows.append(base); continue
                x = prop.loc[valid,ck].values.astype(float)
                y = prop.loc[valid,"pct"].values
                if test_type=="spearman":
                    rho,p = spearmanr(x,y)
                    base.update({"test":"Spearman rho","rho_or_stat":round(rho,4),
                                  "p_raw":round(p,6)})
                else:
                    g0,g1 = y[x==0],y[x==1]
                    if len(g0)<2 or len(g1)<2:
                        base["note"]="insufficient groups"
                    else:
                        stat,p = mannwhitneyu(g1,g0,alternative="greater")
                        base.update({"test":"Mann-Whitney U","rho_or_stat":round(stat,1),
                                     "p_raw":round(p,6),
                                     "note":f"n_dem={len(g1)}, n_nodm={len(g0)}"})
                for a in adnc_order:
                    sub = prop.loc[valid&(prop[adnc_col]==a),"pct"].values
                    base[f"mean_pct_{a.replace(' ','_')}"] = (
                        round(sub.mean(),3) if len(sub)>0 else np.nan)
                rows.append(base)

        df  = pd.DataFrame(rows)
        if df.empty:
            return df   # no eligible supertypes
        out = []
        for st in df.supertype.unique():
            sub = df[df.supertype==st].copy()
            vp  = sub["p_raw"].notna()
            if vp.sum()>0:
                _,p_adj,_,_ = multipletests(sub.loc[vp,"p_raw"].values, method="fdr_bh")
                sub.loc[vp,"p_BH"] = p_adj
            else:
                sub["p_BH"] = np.nan
            out.append(sub)
        df = pd.concat(out,ignore_index=True)
        df["sig_raw"] = df["p_raw"].apply(
            lambda p: ("***" if p<0.001 else "**" if p<0.01 else "*" if p<0.05 else "ns")
            if pd.notna(p) else "")
        df["sig_BH"]  = df["p_BH"].apply(
            lambda p: ("***" if p<0.001 else "**" if p<0.01 else "*" if p<0.05 else "ns")
            if pd.notna(p) else "")
        return df

    eligible_st = [st for st in stype_all.index
                   if int((top_mask & (obs_full[stype_col]==st).values).sum()) >= 20]
    print(f"\n  Eligible supertypes (M1 baseline): {eligible_st}")

    all_clinical = {}
    for mkey in ALL_METHODS:
        print(f"\n  [{mkey}] running clinical suite ...")
        if mkey in INTERSECTION_PAIRS:
            # M4 uses the canonical mask; M5 and M6 are recomputed.
            if mkey == "M4_M1andM2":
                _int_mask = M4_TOP_MASK
            else:
                _pa, _pb  = INTERSECTION_PAIRS[mkey]
                _mask_a   = score_cols[_pa] >= np.percentile(score_cols[_pa], (1-TOP_PCT)*100)
                _mask_b   = score_cols[_pb] >= np.percentile(score_cols[_pb], (1-TOP_PCT)*100)
                _int_mask = _mask_a & _mask_b
            _elig_int = [st for st in stype_all.index
                         if int((_int_mask & (obs_full[stype_col]==st).values).sum()) >= 20]
            print(f"    Intersection mask n={int(_int_mask.sum()):,}  "
                  f"eligible supertypes: {_elig_int}")
            res = _run_clinical_suite(obs_full, mkey, _elig_int,
                                       all_donor_totals, donor_clinical,
                                       top_mask_override=_int_mask)
        else:
            res = _run_clinical_suite(obs_full, mkey, eligible_st,
                                       all_donor_totals, donor_clinical)
        all_clinical[mkey] = res
        fname = os.path.join(OUT_DIR, f"S5_clinical_{mkey}.csv")
        res.to_csv(fname, index=False)
        bh_sig = (res["sig_BH"].isin(["*","**","***"])).sum()
        print(f"  Saved {os.path.basename(fname)}  (BH-sig hits: {bh_sig})")


    # ===== SECTION 9: Head-to-head comparison
    print("\n" + "="*70)
    print("SECTION 9: Head-to-head scoring comparison - M1-M6")
    print("="*70)
    np.random.seed(42); random.seed(42)

    WINNER_VARS_POSITIVE = [
        "ADNC rank", "Braak stage", "CERAD score",
        "Thal phase", "APOE e4 dose", "CAA score",
    ]
    WINNER_VAR_NEG   = "MMSE score"
    WINNER_SUPERTYPE = "Micro-PVM_3-SEAAD"

    compare_rows = []
    for mkey in ALL_METHODS:
        res = all_clinical[mkey]
        valid = res[res["p_raw"].notna()]

        pvm3 = valid[valid["supertype"] == WINNER_SUPERTYPE]
        rho_pool_pvm3 = []
        for var in WINNER_VARS_POSITIVE:
            v_rows = pvm3[pvm3["variable"] == var]["rho_or_stat"].dropna()
            rho_pool_pvm3.extend(v_rows.tolist())
        mmse_rho_pvm3 = pvm3[pvm3["variable"] == WINNER_VAR_NEG]["rho_or_stat"].dropna()
        rho_pool_pvm3.extend(abs(mmse_rho_pvm3).tolist())
        mean_winner_rho = (round(float(np.mean(rho_pool_pvm3)), 4)
                           if rho_pool_pvm3 else np.nan)

        pos     = valid[valid["green_expected"] == "pos"]
        rho_abs = pos["rho_or_stat"].dropna().abs()
        compare_rows.append({
            "method"             : mkey,
            "label"              : SCORER_LABELS[mkey],
            "n_BH_sig"           : int((res["sig_BH"].isin(["*","**","***"])).sum()),
            "n_raw_sig"          : int((res["sig_raw"].isin(["*","**","***"])).sum()),
            "mean_winner_rho"    : mean_winner_rho,
            "mean_abs_rho_pos"   : round(rho_abs.mean(),4)   if len(rho_abs)>0 else np.nan,
            "median_abs_rho_pos" : round(rho_abs.median(),4) if len(rho_abs)>0 else np.nan,
            "n_correct_dir"      : int((pos["rho_or_stat"]>0).sum()),
            "n_pos_tests"        : len(pos),
        })

    compare_df = pd.DataFrame(compare_rows)
    compare_df.to_csv(os.path.join(OUT_DIR,"S6_method_comparison.csv"),index=False)

    print(f"\n  Head-to-head (positive-expected tests only):")
    print(compare_df.to_string(index=False))

    best_method = compare_df.loc[compare_df["mean_winner_rho"].idxmax(),"method"]
    print(f"\n  WINNER by mean_winner_rho : {best_method}  ({SCORER_LABELS[best_method]})")
    for _,_cr in compare_df.iterrows():
        print(f"    {_cr['method']:<15}: mean_winner_rho={_cr['mean_winner_rho']:.4f}  n_BH_sig={_cr['n_BH_sig']}")


    # ===== SECTION 10: Final report
    print("\n" + "="*70)
    print("SECTION 10: Part 2 Final Report")
    print("="*70)
    best_row=compare_df.loc[compare_df["mean_winner_rho"].idxmax()]
    print(f"""
Part 2 complete - neocortical arm
{'='*70}
RELEASE   : 2026-06-22 Multiregion  |  X = ln(UP10K+1), CP10K
FILE      : {FULL_H5AD}
FLAGGED   : {'KEPT' if KEEP_FLAGGED_NUCLEI else 'DROPPED'}
GUARDRAIL : {present_nc}  ({n_after:,} cells retained)
DONORS    : {_n_donors_arm} contributing to this arm
TOP-{TOP_PCT*100:.1f}%    : threshold={top_thresh:.4f}  n={n_top:,}

PRIMARY ENRICHMENT [M1]:
{enrich_df.head(5)[['supertype','n_total','pct_top','enrichment_fold','boot_ci_lo','boot_ci_hi']].to_string(index=False)}

ADNC [M1]:
{adnc_df.head(5)[['supertype','ADNC_rho','ADNC_p']].to_string(index=False)}

DONOR STAGING [M1]:
  ADNC  rho = {rho_adnc:+.4f}  p = {p_adnc:.6f}
  Braak rho = {rho_braak:+.4f}  p = {p_braak:.6f}  (n={n_bv})

METHOD COMPARISON:
{compare_df[['method','n_BH_sig','mean_abs_rho_pos','n_correct_dir']].to_string(index=False)}

WINNER : {best_row['method']}  ({best_row['label']})
  BH-sig  : {best_row['n_BH_sig']}
  mean_winner_rho: {best_row['mean_winner_rho']:.4f}
  mean_abs_rho   : {best_row['mean_abs_rho_pos']:.4f}
  Dir OK  : {best_row['n_correct_dir']}/{best_row['n_pos_tests']}
{'='*70}
""")

    print(f"\nAll outputs in {OUT_DIR}/:")
    for fn in sorted(os.listdir(OUT_DIR)):
        if not fn.startswith("."): print(f"  {fn}")


    # ===== SECTION 11: Export the namespace for later parts (joblib)
    print("\n" + "="*70)
    print(f"SECTION 11: Export canonical namespace → {os.path.join(OUT_DIR, 'Part2_neocortex_env.joblib')}")
    print("="*70)

    import joblib as _joblib

    _env_export = {
        # Core cell metadata
        "obs_full"       : obs_full,          # neocortical-cell DataFrame with all scores
                                              # and obs_full["is_M4_top"] column
        # Score arrays, aligned to obs_full rows
        "score_cols"     : score_cols,        # dict: mkey -> np.ndarray
        # Canonical M4 mask
        "M4_TOP_MASK"    : M4_TOP_MASK,       # bool array, shape (n_cells,)
        "N_M4_TOP"       : N_M4_TOP,          # int
        "M4_THR_M1"      : M4_THR_M1,        # float, M1 score threshold
        "M4_THR_M2"      : M4_THR_M2,        # float, M2 score threshold
        "TOP_PCT"        : TOP_PCT,           # cutoff of this pass
        # UMAP coordinates (if available)
        "umap_coords"    : umap_coords if has_umap else None,
        # Clinical and scoring outputs
        "all_clinical"   : all_clinical,
        "compare_df"     : compare_df,
        "enrich_df"      : enrich_df,
        "adnc_df"        : adnc_df,
        "donor_prop"     : donor_prop,
        "donor_clinical" : donor_clinical,
        "eligible_st"    : eligible_st,
        # Scorer metadata
        "SCORERS"        : SCORERS,
        "SCORER_LABELS"  : SCORER_LABELS,
        "INTERSECTION_PAIRS": INTERSECTION_PAIRS,
        "ALL_METHODS"    : ALL_METHODS,
        # Column names
        "donor_col"      : donor_col,
        "stype_col"      : stype_col,
        "region_col"     : region_col,
        "adnc_col"       : adnc_col,
        "braak_col"      : braak_col,
        # Clinical ordering
        "adnc_order"     : adnc_order,
        "adnc_rank_map"  : adnc_rank_map,
        "adnc_pal"       : adnc_pal,
        "CLINICAL_TESTS" : CLINICAL_TESTS,
        "NEOCORTEX_REGIONS": NEOCORTEX_REGIONS,
        # Provenance
        "RELEASE"        : "2026-06-22 Multiregion",
        "FULL_H5AD"      : FULL_H5AD,
        "NORM_TARGET_SUM": NORM_TARGET_SUM,
        "NORM_CEILING"   : NORM_CEILING,
        "KEEP_FLAGGED_NUCLEI": KEEP_FLAGGED_NUCLEI,
        "n_donors_arm"   : _n_donors_arm,
        # Ordinal maps
        "BRAAK_MAP"      : BRAAK_MAP,
        "CERAD_MAP"      : CERAD_MAP,
        "THAL_MAP"       : THAL_MAP,
        "CAA_MAP"        : CAA_MAP,
        # Column names used by donor_clinical
        "cerad_col"      : cerad_col,
        "thal_col"       : thal_col,
        "caa_col"        : caa_col,
        "apoe_col"       : apoe_col,
        "mmse_col"       : mmse_col,
        "age_col"        : age_col,
        "pmi_col"        : pmi_col,
        "rin_col"        : rin_col,
        "ph_col"         : ph_col,
        "weight_col"     : weight_col,
        "severe_col"     : severe_col,
    }

    # Namespace file in this cutoff's folder.
    _env_path = os.path.join(OUT_DIR, "Part2_neocortex_env.joblib")
    _joblib.dump(_env_export, _env_path, compress=3)
    print(f"  Saved: {os.path.abspath(_env_path)}")
    print(f"  Keys exported: {sorted(_env_export.keys())}")

    # == Sweep bookkeeping ========================================================
    if _IS_PRIMARY:
        # Copy in the part folder, which later parts read.
        _primary_env = os.path.join(ARM_DIR, os.path.basename(_env_path))
        _joblib.dump(_env_export, _primary_env, compress=3)
        print(f"\n  [PRIMARY] also written to:")
        print(f"    {os.path.abspath(_primary_env)}")

    # Per-supertype %top at THIS cutoff, for the sweep summary and the rule.
    _pct_now = {st: 100.0 * int((M4_TOP_MASK & _st_masks[st]).sum())
                     / max(_n_tot_st[st], 1)
                for st in _st_masks}
    _oth_now = {st: v for st, v in _pct_now.items()
                if st not in RULE_ENRICHED_STYPES}
    _oth_gd  = {st: v for st, v in _oth_now.items() if _n_tot_st[st] >= RULE_MIN_N}
    _worst_raw = max(_oth_now, key=_oth_now.get) if _oth_now else None
    _worst_gd  = max(_oth_gd,  key=_oth_gd.get)  if _oth_gd  else None
    _vrow = verdict_df[(verdict_df["method"] == SCAN_METHOD) &
                       (np.isclose(verdict_df["top_pct"], TOP_PCT * 100, atol=1e-6))]

    _SWEEP_RESULTS[round(TOP_PCT * 100, 3)] = {
        "top_pct"          : round(TOP_PCT * 100, 3),
        "is_primary"       : _IS_PRIMARY,
        "n_top_M1"         : int(n_top),
        "N_M4_TOP"         : int(N_M4_TOP),
        "M4_thr_M1"        : round(float(M4_THR_M1), 6),
        "M4_thr_M2"        : round(float(M4_THR_M2), 6),
        "pct_PVM3"         : round(_pct_now.get("Micro-PVM_3-SEAAD", np.nan), 3),
        "pct_PVM4"         : round(_pct_now.get("Micro-PVM_4-SEAAD", np.nan), 3),
        "worst_other"      : _worst_raw,
        "worst_other_pct"  : round(_oth_now[_worst_raw], 3) if _worst_raw else np.nan,
        "worst_other_n"    : _n_tot_st.get(_worst_raw, np.nan),
        "worst_other_gd"   : _worst_gd,
        "worst_other_gd_pct": round(_oth_gd[_worst_gd], 3) if _worst_gd else np.nan,
        "rule_pass_raw"    : (bool(_vrow.iloc[0]["pass_raw"]) if len(_vrow) else np.nan),
        "rule_pass_guarded": (bool(_vrow.iloc[0]["pass_guarded"]) if len(_vrow) else np.nan),
        "margin_cells_raw" : (int(_vrow.iloc[0]["margin_cells_raw"]) if len(_vrow) else np.nan),
        "margin_cells_gd"  : (int(_vrow.iloc[0]["margin_cells_guarded"]) if len(_vrow) else np.nan),
        "donor_rho_ADNC"   : round(float(rho_adnc), 4),
        "donor_rho_Braak"  : round(float(rho_braak), 4),
        "best_method"      : best_method,
        "mean_winner_rho"  : float(best_row["mean_winner_rho"]),
        "out_dir"          : OUT_DIR,
    }
    print(f"\n  [SWEEP {TOP_PCT*100:.1f}%] M4 n={N_M4_TOP:,} | PVM_3 "
          f"{_SWEEP_RESULTS[round(TOP_PCT*100,3)]['pct_PVM3']:.1f}% | PVM_4 "
          f"{_SWEEP_RESULTS[round(TOP_PCT*100,3)]['pct_PVM4']:.1f}% | worst other "
          f"{_worst_raw} {_SWEEP_RESULTS[round(TOP_PCT*100,3)]['worst_other_pct']:.1f}% "
          f"(n={_n_tot_st.get(_worst_raw, 0):,}) | rule raw="
          f"{_SWEEP_RESULTS[round(TOP_PCT*100,3)]['rule_pass_raw']} "
          f"guarded={_SWEEP_RESULTS[round(TOP_PCT*100,3)]['rule_pass_guarded']}")



# %% SECTION 12: Cross-cutoff sweep summary
print("\n" + "=" * 74)
print("SECTION 12: Cutoff sweep summary")
print("=" * 74)

_sw_df = pd.DataFrame(_SWEEP_RESULTS).T.sort_values("top_pct")
_sw_df.to_csv(os.path.join(ARM_DIR, "S7_sweep_summary.csv"), index=False)
print(_sw_df.to_string(index=False))
print(f"\n  Saved {os.path.join(ARM_DIR, 'S7_sweep_summary.csv')}")
print(f"  Total sweep time: {time.time()-_t_sweep:.0f}s")

print(f"\n  Live namespace now holds the PRIMARY cutoff "
      f"({PRIMARY_TOP_PCT*100:.1f}%):")
print(f"    TOP_PCT      = {TOP_PCT}")
print(f"    OUT_DIR      = {OUT_DIR}")
print(f"    N_M4_TOP     = {N_M4_TOP:,}")
print(f"    M4_TOP_MASK  = bool array, {M4_TOP_MASK.sum():,} True of "
      f"{len(M4_TOP_MASK):,}")
