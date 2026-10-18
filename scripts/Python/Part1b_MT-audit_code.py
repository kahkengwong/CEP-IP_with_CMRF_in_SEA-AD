# =============================================================================
# Part1b_MT-audit_code.py
# SEA-AD CEP-IP/CMRF, Part 1b, mitochondrial audit (Figure 1G, Table S1)
#   Pre-check of technical confounds before Part 2: matrix type and layers,
#   traces of ambient correction, the mitochondrial genes present, mito% and
#   depth against pathology, MT-CO3 abundance against ACTB, and the MT-CO3 to
#   cDFG correlation conditioned on mito genes outside Module 1.
#
# Usage (from the SEA-AD folder, or set SEAAD_ROOT), after Part 1a
#   python3 Part1b_MT-audit_code.py
#
# Outputs (Part1_PreChecks/): mito gene inventory (CSV) and report.
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
OUT_DIR = "Part1_PreChecks"   # created under the SEA-AD folder
os.makedirs(OUT_DIR, exist_ok=True)

# == Module 1 definition (as in R_config.R) ===================================
GOI = "MT-CO3"
DFG = ["MT-ATP6", "MT-CO2", "MT-ND3", "MT-ND4",
       "MT-CYB", "MT-ND1", "MT-ND2", "MT-CO1"]
MODULE1 = [GOI] + DFG                      # 9 genes
REF_GENE = "ACTB"

# == obs column names: defaults, resolved against the file below ==============
stype_col  = "Supertype"
region_col = "Brain Region"
donor_col  = "Donor ID"
adnc_col   = "Overall AD neuropathological Change"
braak_col  = "Braak"
TARGET_STYPE = "Micro-PVM_3-SEAAD"

# Normalisation constant (CP10K); X cannot exceed log1p(NORM_TARGET_SUM).
NORM_TARGET_SUM = 1e4          # CP10K, per Part1a and uns['X_normalization']

# Documented SEA-AD supertype renames, reported if TARGET_STYPE is missing.
KNOWN_RENAMES = {
    "VLMC_2":              "Pericyte_1",
    "VLMC_2_1-SEAAD":      "SMC-SEAAD",
    "VLMC_2_2-SEAAD":      "Pericyte_2-SEAAD",
    "Micro-PVM_2_2-SEAAD": "Lymphocyte",
    "Micro-PVM_1_1-SEAAD": "Monocyte",
}


def find_col(cols, keywords):
    """Resolve one obs column name. Exact match wins, then case-insensitive
    substring (shortest hit, to avoid matching a longer derived column).
    Returns None when nothing matches."""
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


# Braak stages are Roman numerals in SEA-AD (same map as Part 2).
BRAAK_MAP = {"Braak 0": 0, "Braak I": 1, "Braak II": 2, "Braak III": 3,
             "Braak IV": 4, "Braak V": 5, "Braak VI": 6}


def parse_braak(s):
    """Map Braak strings to integers. Roman numerals first (the SEA-AD
    convention), then an Arabic-digit fallback for other encodings.
    Returns (values, list of unmapped distinct strings)."""
    txt = pd.Series(s).astype(str).str.strip()
    v = txt.map(BRAAK_MAP)
    unmapped = v.isna()
    if unmapped.any():
        arabic = pd.to_numeric(
            txt[unmapped].str.extract(r"(\d+)")[0], errors="coerce")
        v.loc[unmapped] = arabic
    still = sorted(set(txt[v.isna()])) if v.isna().any() else []
    return pd.to_numeric(v, errors="coerce").values, still

_report = []
def say(*a):
    line = " ".join(str(x) for x in a)
    print(line)
    _report.append(line)

def hdr(t):
    say("")
    say("=" * 74)
    say("  " + t)
    say("=" * 74)


# =============================================================================
# SECTION A - backed pass, metadata only
# =============================================================================
hdr("Q1a/Q2 - Structure and provenance (metadata only, backed mode)")

ad = sc.read_h5ad(FULL_H5AD, backed="r")
say(f"  shape                : {ad.shape}")
say(f"  X dtype (declared)   : {ad.X.dtype if ad.X is not None else '(not exposed in backed mode)'}")
_lk = [k for k in (ad.layers.keys() if ad.layers else []) if k is not None]
say(f"  layers present       : {_lk if _lk else '(none)'}")
say(f"  raw present          : {ad.raw is not None}")
say(f"  obsm keys            : {list(ad.obsm.keys())}")
say(f"  uns keys             : {list(ad.uns.keys())}")


# =============================================================================
# Resolve obs columns against the file; stop if a critical one is missing
# =============================================================================
hdr("Q0 - obs column resolution and TARGET_STYPE assertion")

_oc = list(ad.obs.columns)
_defaults = {
    "stype_col":  stype_col,
    "region_col": region_col,
    "donor_col":  donor_col,
    "adnc_col":   adnc_col,
    "braak_col":  braak_col,
}
_resolved = {
    "stype_col":  find_col(_oc, [stype_col, "Supertype", "supertype"]),
    "region_col": find_col(_oc, [region_col, "Brain Region", "brain_region",
                                 "Region", "region", "tissue"]),
    "donor_col":  find_col(_oc, [donor_col, "Donor ID", "donor_id", "Donor"]),
    "adnc_col":   find_col(_oc, [adnc_col, "AD neuropath", "ADNC"]),
    "braak_col":  find_col(_oc, [braak_col, "Braak", "braak"]),
}

say(f"  obs has {len(_oc)} columns")
say("")
for _k in _defaults:
    _d, _v = _defaults[_k], _resolved[_k]
    _tag = ("ok" if _v == _d else
            "RENAMED" if _v is not None else "*** MISSING ***")
    say(f"    {_k:<11} default={_d:<40} found={str(_v):<40} {_tag}")

stype_col  = _resolved["stype_col"]  or stype_col
region_col = _resolved["region_col"] or region_col
donor_col  = _resolved["donor_col"]  or donor_col
adnc_col   = _resolved["adnc_col"]   or adnc_col
braak_col  = _resolved["braak_col"]  or braak_col


def _bail(msg):
    """Write the report so far, then stop."""
    say("")
    say("  " + "!" * 70)
    say(f"  ABORTING: {msg}")
    say("  " + "!" * 70)
    say("")
    say("  full obs column list:")
    for _i, _c in enumerate(_oc):
        say(f"    [{_i:>3}] {_c}")
    with open(os.path.join(OUT_DIR, "S1b_precheck_report.txt"), "w",
              encoding="utf-8") as _f:
        _f.write("\n".join(_report))
    print(f"\n  report written to {OUT_DIR}/S1b_precheck_report.txt")
    raise SystemExit(msg)


# Critical columns: the run stops if one is missing.
_bad = [k for k in ("stype_col", "region_col", "donor_col")
        if _resolved[k] is None]
if _bad:
    _bail(f"cannot resolve critical obs columns {_bad}")

# ADNC and Braak are optional; without them the Q4 pathology tests are skipped.
for _k in ("adnc_col", "braak_col"):
    if _resolved[_k] is None:
        say("")
        say(f"  [WARN] {_k} unresolved; Q4 pathology correlations skipped.")

# Does the target supertype still exist under this exact name?
_st_vals = ad.obs[stype_col].astype(str)
_st_set = set(_st_vals.unique())
say("")
say(f"  distinct supertypes in '{stype_col}' : {len(_st_set)}")
if TARGET_STYPE in _st_set:
    say(f"  '{TARGET_STYPE}' present : YES  "
        f"(n={int((_st_vals == TARGET_STYPE).sum()):,})")
else:
    say(f"  '{TARGET_STYPE}' present : NO")
    say("")
    say("  microglial / PVM-like supertypes actually present:")
    for _s in sorted(_st_set):
        if "PVM" in _s or "icro" in _s or _s in KNOWN_RENAMES.values():
            say(f"    {_s:<40} n={int((_st_vals == _s).sum()):>9,}")
    say("")
    say("  documented renames to check first (old -> new):")
    for _a, _b in KNOWN_RENAMES.items():
        _mark = "  <- present in file" if _b in _st_set else ""
        say(f"    {_a:<24} -> {_b}{_mark}")
    _bail(f"TARGET_STYPE '{TARGET_STYPE}' not present in {stype_col}")

# Print the uns provenance values (X_normalization documents the CP10K scale).
say("")
say("  --- uns provenance values (printed in full) ---")
_PROV_KEYS = ["X_normalization", "log1p", "title", "batch_condition",
              "default_embedding", "normalization", "schema_version"]
for _k in (ad.uns.keys() if ad.uns else []):
    if any(_p.lower() in str(_k).lower() for _p in _PROV_KEYS):
        try:
            _v = ad.uns[_k]
        except Exception as _e:
            _v = f"(unreadable: {type(_e).__name__})"
        _s = str(_v)
        say(f"    uns['{_k}'] = {_s[:400]}" + (" ..." if len(_s) > 400 else ""))
say("")
say("  --- other uns keys ---")
for _k in (ad.uns.keys() if ad.uns else []):
    if not any(_p.lower() in str(_k).lower() for _p in _PROV_KEYS):
        try:
            _v = ad.uns[_k]
            _d = (f"{type(_v).__name__}, len={len(_v)}"
                  if hasattr(_v, "__len__") else type(_v).__name__)
        except Exception:
            _d = "(unreadable)"
        say(f"    uns['{_k}'] : {_d}")


# SEA-AD inclusion flags: 'Used in analysis' and 'Neurotypical reference'.
hdr("Q0b - SEA-AD inclusion flags")
USED_COL = find_col(_oc, ["Used in analysis", "used_in_analysis",
                          "Used in Analysis"])
NREF_COL = find_col(_oc, ["Neurotypical reference", "neurotypical"])

if USED_COL is None:
    say("  'Used in analysis' not present in this file - nothing to filter on.")
else:
    _u = ad.obs[USED_COL].astype(str).str.strip()
    say(f"  '{USED_COL}' value counts:")
    for _val, _n in _u.value_counts(dropna=False).items():
        say(f"    {str(_val):<20} {int(_n):>9,}  ({100*_n/len(_u):.3f}%)")
    _false = _u.str.lower().isin(["false", "0", "no", "n", "nan", "none"])
    if _false.any():
        say("")
        say(f"  {int(_false.sum()):,} nuclei not flagged for analysis "
            f"({100*_false.mean():.3f}%)")
        say("      Per region:")
        _rr = ad.obs[region_col].astype(str)
        for _r in sorted(set(_rr)):
            _m = (_rr == _r).values
            say(f"        {_r:<6} {int((_false.values & _m).sum()):>8,}"
                f" / {int(_m.sum()):>8,}")
        say(f"      In '{TARGET_STYPE}': "
            f"{int((_false.values & (_st_vals == TARGET_STYPE).values).sum()):,}")
    else:
        say("  All nuclei are flagged for analysis.")

if NREF_COL is None:
    say("")
    say("  'Neurotypical reference' not present.")
else:
    _nr = ad.obs[NREF_COL].astype(str).str.strip().str.lower()
    _isref = _nr.isin(["true", "1", "yes", "y"])
    say("")
    say(f"  '{NREF_COL}': {int(_isref.sum()):,} nuclei flagged as reference")
    if _isref.any():
        say("  [WARN] reference nuclei present.")
    else:
        say("  No reference nuclei.")

# uns['log1p'] is the strongest available provenance marker
if ad.uns and "log1p" in ad.uns:
    say("")
    say(f"  uns['log1p']         : {ad.uns['log1p']}")

say("")
say("  --- ambient-correction keyword hunt (strings only) ---")
KEYS = ["cellbender", "soupx", "decontx", "ambient", "background",
        "remove-background", "scran", "sctransform", "normali", "pooled"]
hits = []
def scan_text(where, text):
    t = str(text).lower()
    for k in KEYS:
        if k in t:
            hits.append((where, k, str(text)[:200]))

for k in (ad.uns.keys() if ad.uns else []):
    scan_text(f"uns key '{k}'", k)
    try:
        scan_text(f"uns['{k}'] value", ad.uns[k])
    except Exception:
        pass
for c in ad.obs.columns:
    scan_text(f"obs column '{c}'", c)
for c in (ad.var.columns if ad.var is not None else []):
    scan_text(f"var column '{c}'", c)

if hits:
    seen = set()
    for where, k, snippet in hits:
        if (where, k) in seen:
            continue
        seen.add((where, k))
        say(f"  [HIT] {k:<18} in {where}")
        say(f"        {snippet}")
else:
    say("  no ambient/normalisation keywords in uns/obs/var names or uns values")

say("")
say("  --- all obs columns ---")
for c in ad.obs.columns:
    say(f"    {c}")

QC_PAT = ["mito", "mt_", "mt-", "pct_counts", "n_counts", "ncount", "n_genes",
          "nfeature", "total_counts", "doublet", "ribo", "fraction",
          "complexity", "qc"]
qc_cols = [c for c in ad.obs.columns if any(p in c.lower() for p in QC_PAT)]
say("")
say(f"  QC-looking obs columns : {qc_cols if qc_cols else '(none)'}")

say("")
say("  --- var columns ---")
for c in (ad.var.columns if ad.var is not None else []):
    say(f"    {c}")

var_names = list(ad.var_names)
obs = ad.obs.copy()
layer_names = _lk
try:
    ad.file.close()
except Exception:
    pass
del ad; gc.collect()


# =============================================================================
# Q3 - mito gene inventory (var_names only)
# =============================================================================
hdr("Q3 - Mitochondrial gene inventory in the feature space")

mt_all = sorted([g for g in var_names if g.startswith("MT-")])
say(f"  genes starting 'MT-'  : {len(mt_all)}")
say(f"    {mt_all}")

missing_mod = [g for g in MODULE1 if g not in var_names]
say("")
say(f"  Module 1 present      : {len(MODULE1) - len(missing_mod)}/{len(MODULE1)}")
if missing_mod:
    say(f"  [WARN] MISSING        : {missing_mod}")

mt_disjoint = [g for g in mt_all if g not in MODULE1]
say("")
say(f"  DISJOINT mito axis    : {len(mt_disjoint)} genes")
say(f"    {mt_disjoint}")

say("")
for key in ["MT-RNR1", "MT-RNR2", "MT-ND5", "MT-ND4L", "MT-ND6", "MT-ATP8"]:
    say(f"    {key:<10} present={key in var_names}")

if not mt_disjoint:
    say("")
    say("  [WARN] no disjoint MT genes; Q6 cannot run.")


# =============================================================================
# SECTION B - full load, then all numeric work
# =============================================================================
hdr("Loading full matrix into memory")

ad = sc.read_h5ad(FULL_H5AD)
say(f"  loaded                : {ad.shape}")

COUNT_LAYER = None
for cand in ["UMIs", "counts", "raw_counts"]:
    if ad.layers and cand in ad.layers:
        COUNT_LAYER = cand
        break
say(f"  count layer           : {COUNT_LAYER if COUNT_LAYER else 'NONE FOUND'}")


# == Q1b: numeric sniff =======================================================
hdr("Q1b - Matrix classification (numeric, in-memory)")

def sniff(M, name):
    if M is None:
        say(f"  {name:<22} [None -- not present]")
        return None
    n_s = min(3000, M.shape[0])
    rows = np.sort(np.random.choice(M.shape[0], n_s, replace=False))
    S = M[rows]
    S = S.toarray() if sp.issparse(S) else np.asarray(S)
    S = np.asarray(S, dtype=np.float64)        # kills object dtype
    nz = S[S != 0]
    if nz.size == 0:
        say(f"  {name:<22} all zero in sample of {n_s}")
        return None
    integer_like = bool(np.allclose(nz, np.rint(nz), atol=1e-8))
    say(f"  {name:<22} min={S.min():.4f}  max={S.max():.4f}  "
        f"mean_nz={nz.mean():.4f}  integer_like={integer_like}   "
        f"[SAMPLED n={n_s:,}]")
    if integer_like:
        say(f"  {'':22}   -> RAW COUNTS")
    elif S.max() < 30:
        say(f"  {'':22}   -> LOG-scale normalised")
    else:
        say(f"  {'':22}   -> LINEAR normalised (not logged)")
    return integer_like

x_is_counts = sniff(ad.X, "X")
for ln in layer_names:
    sniff(ad.layers[ln], f"layers['{ln}']")

# Exact global extremes (the sniff above uses a 3,000-cell sample).
say("")
say("  --- exact global extremes (not sampled) ---")


def exact_minmax(M, name):
    if M is None:
        return None
    try:
        mx = float(M.max())
        mn = float(M.min())
    except Exception as e:
        say(f"  {name:<22} could not compute exact extremes: "
            f"{type(e).__name__}: {e}")
        return None
    say(f"  {name:<22} exact min={mn:.4f}  exact max={mx:.4f}")
    return mx


X_MAX_EXACT = exact_minmax(ad.X, "X")
for ln in layer_names:
    exact_minmax(ad.layers[ln], f"layers['{ln}']")

if X_MAX_EXACT is not None:
    _ceil = float(np.log1p(NORM_TARGET_SUM))
    say("")
    say(f"  NORM_CEILING = log1p({NORM_TARGET_SUM:,.0f}) = {_ceil:.4f}")
    say(f"  headroom     = {_ceil - X_MAX_EXACT:+.4f}")
    if X_MAX_EXACT > _ceil:
        say("  [FAIL] X exceeds the ceiling implied by the target_sum.")
    else:
        say("  [OK] X is under the ceiling.")
    # The same maximum on the CP100K scale, for comparison with CP100K exports.
    _alt = float(np.log1p(np.expm1(X_MAX_EXACT) * 10.0))
    say("")
    say(f"  the same maximum expressed on the CP100K scale : {_alt:.4f}")

if x_is_counts:
    say("")
    say("  [WARN] X looks like raw counts, not normalised expression.")


# == Extract the genes needed, in one block slice =============================
need = [g for g in dict.fromkeys(MODULE1 + mt_all + [REF_GENE]) if g in var_names]
need_idx = [var_names.index(g) for g in need]
say("")
say(f"  extracting {len(need)} genes in one block slice ...")

Craw = ad.layers[COUNT_LAYER] if COUNT_LAYER else ad.X
Xn = ad.X

total_counts = np.asarray(Craw.sum(axis=1)).ravel().astype(float)
n_genes_det = np.asarray((Craw > 0).sum(axis=1)).ravel().astype(float)

mt_idx = [var_names.index(g) for g in mt_all]
mt_counts = np.asarray(Craw[:, mt_idx].sum(axis=1)).ravel().astype(float)
mito_pct = np.where(total_counts > 0, 100.0 * mt_counts / total_counts, np.nan)

Cblk = Craw[:, need_idx]
Cblk = np.asarray(Cblk.toarray() if sp.issparse(Cblk) else Cblk, dtype=float)
Eblk = Xn[:, need_idx]
Eblk = np.asarray(Eblk.toarray() if sp.issparse(Eblk) else Eblk, dtype=float)
Cn = {g: Cblk[:, i] for i, g in enumerate(need)}
E = {g: Eblk[:, i] for i, g in enumerate(need)}
del ad, Cblk, Eblk; gc.collect()
say("  done")


# =============================================================================
# Q4 - mito% and depth, and association with pathology
# =============================================================================
hdr("Q4 - mito% and depth; association with pathology")

def q(v, name, unit=""):
    v = np.asarray(v, dtype=float)
    v = v[~np.isnan(v)]
    if v.size == 0:
        say(f"  {name:<18} (no data)")
        return
    say(f"  {name:<18} min={v.min():.3f}  q25={np.percentile(v,25):.3f}  "
        f"median={np.median(v):.3f}  q75={np.percentile(v,75):.3f}  "
        f"q99={np.percentile(v,99):.3f}  max={v.max():.3f} {unit}")

say("  --- ALL nuclei ---")
q(mito_pct, "mito%", "%")
q(total_counts, "total counts")
q(n_genes_det, "genes detected")

is_target = (obs[stype_col].astype(str).values == TARGET_STYPE)
say("")
say(f"  --- {TARGET_STYPE} (n={int(is_target.sum())}) ---")
q(mito_pct[is_target], "mito%", "%")
q(total_counts[is_target], "total counts")

say("")
say("  --- per region (NeoC vs NonNeoC asymmetry check) ---")
_NEOC = ["DFC", "MTG", "AnG", "STG", "V1C", "ITG"]
_NONNEOC = ["HIP", "MEC", "LEC", "FI"]
_rv = obs[region_col].astype(str).values
say(f"    {'region':<8}{'arm':<9}{'n_nuclei':>10}{'n_PVM3':>9}"
    f"{'mito%_med':>11}{'depth_med':>11}")
for _r in sorted(set(_rv)):
    _m = (_rv == _r)
    _arm = ("NeoC" if _r in _NEOC else
            "NonNeoC" if _r in _NONNEOC else "OTHER")
    _mm = mito_pct[_m]
    _mm = _mm[~np.isnan(_mm)]
    say(f"    {_r:<8}{_arm:<9}{int(_m.sum()):>10,}"
        f"{int((_m & is_target).sum()):>9,}"
        f"{(np.median(_mm) if _mm.size else np.nan):>11.3f}"
        f"{np.median(total_counts[_m]):>11.0f}")
_extra = [r for r in set(_rv) if r not in _NEOC + _NONNEOC]
if _extra:
    say(f"    [WARN] unexpected region values present: {sorted(_extra)}")

say("")
say("  --- donor level: mito% against pathology ---")
if _resolved["adnc_col"] is None and _resolved["braak_col"] is None:
    say("  Skipped: neither ADNC nor Braak resolved in obs (see Q0).")
    gg = pd.DataFrame()
else:
    adnc_rank = {"Not AD": 0, "Low": 1, "Intermediate": 2, "High": 3}
    _cols = {"donor": obs[donor_col].astype(str).values,
             "mito_pct": mito_pct,
             "depth": total_counts}
    if _resolved["adnc_col"] is not None:
        _cols["adnc"] = obs[adnc_col].astype(str).map(adnc_rank).values
    if _resolved["braak_col"] is not None:
        _bv, _bad_braak = parse_braak(obs[braak_col])
        _cols["braak"] = _bv
        _n_ok = int(np.isfinite(_bv).sum())
        say(f"  Braak parsed          : {_n_ok:,}/{len(_bv):,} nuclei")
        if _bad_braak:
            say(f"  [WARN] unmapped Braak values: {_bad_braak[:10]}")
        if _n_ok == 0:
            say("  [WARN] no Braak value parsed; Braak correlation not run.")
    dfd = pd.DataFrame(_cols)
    _agg = dict(mito_pct=("mito_pct", "median"), depth=("depth", "median"))
    for _v in ("adnc", "braak"):
        if _v in dfd.columns:
            _agg[_v] = (_v, "first")
    gg = dfd.groupby("donor").agg(**_agg)
    if "adnc" in gg.columns:
        gg = gg.dropna(subset=["adnc"])
    say(f"  donors in this table  : {len(gg)}")
    if len(gg) not in (43, 84):
        say(f"  [NOTE] n={len(gg)} is neither 84 nor 43.")
    for var in [v for v in ("adnc", "braak") if v in gg.columns]:
        for yv in ["mito_pct", "depth"]:
            ok = gg[[var, yv]].dropna()
            if len(ok) > 5 and ok[var].nunique() > 1:
                r, p = spearmanr(ok[var], ok[yv])
                flag = ("  <-- CONFOUND RISK" if (p < 0.05 and yv == "mito_pct")
                        else "")
                say(f"    {yv:<9} vs {var:<6}: rho={r:+.4f}  p={p:.4g}  "
                    f"n={len(ok)}{flag}")
            else:
                # Report rows that cannot be computed.
                say(f"    {yv:<9} vs {var:<6}: NOT COMPUTED  "
                    f"(n={len(ok)} usable donors, "
                    f"{ok[var].nunique() if len(ok) else 0} distinct levels)")

# Cross-check against any SEA-AD mito column in obs.
for c in qc_cols:
    v = pd.to_numeric(obs[c], errors="coerce").values
    if np.isfinite(v).sum() > 100 and "mito" in c.lower():
        r, p = spearmanr(v[np.isfinite(v)], mito_pct[np.isfinite(v)],
                         nan_policy="omit")
        say(f"  cross-check: computed mito% vs obs['{c}']  rho={r:+.4f} p={p:.4g}")


# =============================================================================
# Q5 - MT-CO3 abundance against ACTB
# =============================================================================
hdr("Q5 - MT-CO3 abundance against ACTB")

# Abundance over all target nuclei (*_all) and over the analysed nuclei
# with GOI > 0 (*_analysed), the set rho is computed on.
sel = is_target & (Cn[GOI] > 0)          # global mask, reused by Q6
analysed_t = Cn[GOI][is_target] > 0      # same mask, indexed into the subset

rows = []
for g_ in [GOI, REF_GENE] + DFG + mt_disjoint:
    if g_ not in Cn:
        continue
    c, e = Cn[g_][is_target], E[g_][is_target]
    if g_ == GOI:
        rho = np.nan
    else:
        rho, _ = spearmanr(E[GOI][sel], E[g_][sel])
    rows.append(dict(gene=g_, in_module1=g_ in MODULE1,
                     is_mito=g_.startswith("MT-"),
                     detect_pct=100.0 * (c > 0).mean(),
                     mean_counts=c.mean(), max_counts=c.max(),
                     mean_norm=e.mean(), max_norm=e.max(),
                     detect_pct_analysed=100.0 * (c[analysed_t] > 0).mean(),
                     mean_counts_analysed=c[analysed_t].mean(),
                     rho_with_GOI=rho))
inv = pd.DataFrame(rows)
inv["ratio_analysed_over_all"] = inv.mean_counts_analysed / inv.mean_counts
say(f"  ({TARGET_STYPE}, n={int(is_target.sum())};  "
    f"analysed, {GOI} > 0, n={int(analysed_t.sum())})")
say(inv[["gene", "in_module1", "detect_pct", "detect_pct_analysed",
         "mean_counts", "mean_counts_analysed", "max_counts",
         "rho_with_GOI"]].to_string(index=False))
say("")
_g = inv[inv.gene == GOI]
if len(_g):
    _mc, _xc = float(_g.mean_counts.iloc[0]), float(_g.max_counts.iloc[0])
    say(f"  {GOI}: mean {_mc:.2f} counts per nucleus, max {_xc:.0f}, detected in "
        f"{float(_g.detect_pct.iloc[0]):.1f}% of {TARGET_STYPE}.")


# =============================================================================
# Q5b - is concordance with the GOI explained by abundance?
# =============================================================================
hdr("Q5b - is concordance with the GOI explained by abundance?")

# Fit rho ~ log10(mean counts) over the non-GOI mito genes and read ACTB off
# the same line, on both abundance axes.


def _fit_trend(x, y):
    slope, intercept = np.polyfit(x, y, 1)
    pred = intercept + slope * x
    r2 = 1.0 - ((y - pred) ** 2).sum() / ((y - y.mean()) ** 2).sum()
    return slope, intercept, r2, pred


_mito = inv[inv.is_mito & inv.rho_with_GOI.notna() & (inv.mean_counts > 0)]
_ref = inv[inv.gene == REF_GENE]
_trend = {}

if len(_mito) >= 5 and len(_ref) == 1:
    _y = _mito.rho_with_GOI.values
    for _label, _col in [("all", "mean_counts"),
                         ("analysed", "mean_counts_analysed")]:
        _x = np.log10(_mito[_col].values)
        _slope, _int, _r2, _pred = _fit_trend(_x, _y)
        _xr = float(np.log10(_ref[_col].iloc[0]))
        _yr = float(_ref.rho_with_GOI.iloc[0])
        _pred_ref = _int + _slope * _xr
        _worst = float(np.max(np.abs(_y - _pred)))
        _worst_gene = _mito.gene.values[int(np.argmax(np.abs(_y - _pred)))]
        _trend[_label] = dict(slope=_slope, intercept=_int, r2=_r2,
                              pred_ref=_pred_ref, resid=_yr - _pred_ref,
                              worst=_worst, worst_gene=_worst_gene)
        # per-gene predictions and residuals on this axis
        inv[f"pred_rho_{_label}"] = _int + _slope * np.log10(
            inv[_col].where(inv[_col] > 0))
        inv[f"resid_{_label}"] = inv.rho_with_GOI - inv[f"pred_rho_{_label}"]

    say(f"  fitted on {len(_mito)} non-GOI mitochondrial genes; "
        f"rho is the same in both fits.")
    say("")
    say(f"  {'':<34}{'all nuclei':>14}{'analysed':>14}")
    say("  " + "-" * 62)
    for _k, _lbl, _f in [("slope", "slope", "{:+.4f}"),
                         ("intercept", "intercept", "{:+.4f}"),
                         ("r2", "R2", "{:.4f}"),
                         ("pred_ref", f"{REF_GENE} predicted rho", "{:+.4f}"),
                         ("resid", f"{REF_GENE} residual", "{:+.4f}"),
                         ("worst", "worst mito residual", "{:+.4f}")]:
        say(f"  {_lbl:<34}{_f.format(_trend['all'][_k]):>14}"
            f"{_f.format(_trend['analysed'][_k]):>14}")
    say(f"  {'ratio to worst mito residual':<34}"
        f"{abs(_trend['all']['resid']) / _trend['all']['worst']:>14.2f}"
        f"{abs(_trend['analysed']['resid']) / _trend['analysed']['worst']:>14.2f}")
    say("")
    say(f"  {REF_GENE} observed rho : "
        f"{float(_ref.rho_with_GOI.iloc[0]):+.4f}  (same in both)")
    say(f"  worst mito residual carried by: {_trend['all']['worst_gene']} "
        f"(all) / {_trend['analysed']['worst_gene']} (analysed)")

    # The same trend against detection rate.
    for _label, _col in [("all", "detect_pct"),
                         ("analysed", "detect_pct_analysed")]:
        _xd = _mito[_col].values
        if np.ptp(_xd) == 0:
            continue
        _, _, _r2d, _ = _fit_trend(_xd, _y)
        say(f"  vs detection% ({_label:<8}): R2 {_r2d:.4f}")

    say("")
    _a = _trend["analysed"]
    if abs(_a["resid"]) >= 3.0 * _a["worst"]:
        say(f"  {REF_GENE} is {abs(_a['resid']) / _a['worst']:.1f}x the worst mito residual off the trend")
    else:
        say(f"  {REF_GENE} is not clearly off the trend "
            f"({abs(_a['resid']) / _a['worst']:.1f}x the worst mito residual).")
else:
    say(f"  [skipped] need >=5 mito genes with rho and exactly one {REF_GENE} "
        f"row; have {len(_mito)} and {len(_ref)}")


# =============================================================================
# Q6 - conditioning on a disjoint mito axis
# =============================================================================
hdr("Q6 - conditioning on a disjoint mito axis")

# sel is defined in Q5.
say(f"  cells: {TARGET_STYPE} with {GOI} > 0   n={int(sel.sum())}")

goi = E[GOI][sel]
dfg_have = [g_ for g_ in DFG if g_ in E]
cdfg = np.mean(np.vstack([E[g_][sel] for g_ in dfg_have]), axis=0)
say(f"  cDFG built from {len(dfg_have)}/{len(DFG)} DFGs")

r0, p0 = spearmanr(goi, cdfg)
say("")
say(f"  UNCONDITIONAL    rho({GOI}, cDFG) = {r0:+.4f}   p={p0:.4g}")

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

avail_dis = [g_ for g_ in mt_disjoint if g_ in E]
if avail_dis:
    camb = np.mean(np.vstack([E[g_][sel] for g_ in avail_dis]), axis=0)
    r1, p1 = partial_spearman(goi, cdfg, camb)
    say(f"  | cAMB ({len(avail_dis)} disjoint MT)          rho = {r1:+.4f}  "
        f"p={p1:.4g}   drop={r0-r1:+.4f}")

rrna = [g_ for g_ in ["MT-RNR1", "MT-RNR2"] if g_ in E]
if rrna:
    crrna = np.mean(np.vstack([E[g_][sel] for g_ in rrna]), axis=0)
    r2, p2 = partial_spearman(goi, cdfg, crrna)
    say(f"  | mt-rRNA {str(rrna):<26} rho = {r2:+.4f}  p={p2:.4g}   "
        f"drop={r0-r2:+.4f}")

r3, p3 = partial_spearman(goi, cdfg, total_counts[sel])
say(f"  | total counts                     rho = {r3:+.4f}  p={p3:.4g}   "
    f"drop={r0-r3:+.4f}")
r4, p4 = partial_spearman(goi, cdfg, np.nan_to_num(mito_pct[sel]))
say(f"  | mito% (circular, reference only) rho = {r4:+.4f}  p={p4:.4g}")

if avail_dis:
    both = np.vstack([camb, total_counts[sel]])
    r5, p5 = partial_spearman(goi, cdfg, both)
    say(f"  | cAMB + total counts             rho = {r5:+.4f}  p={p5:.4g}   "
        f"drop={r0-r5:+.4f}")

say("")
say("  --- ambient directional prediction: strength vs depth ---")
tc = total_counts[sel]
qs = np.quantile(tc, [0, .25, .5, .75, 1.0])
for i in range(4):
    m = (tc >= qs[i]) & (tc <= qs[i + 1])
    if m.sum() > 30:
        rr, _ = spearmanr(goi[m], cdfg[m])
        say(f"    depth Q{i+1} (n={int(m.sum()):5d}, median={np.median(tc[m]):9.0f}): "
            f"rho={rr:+.4f}")


# =============================================================================
# Q7 - implied size factors against mito%
# =============================================================================
hdr("Q7 - implied size factors vs mito% (did normalisation absorb it?)")

try:
    idx_sel = np.where(sel)[0]
    ns = min(3000, idx_sel.size)
    pick = np.random.choice(idx_sel, ns, replace=False)
    Cs = np.vstack([Cn[g_][pick] for g_ in need]).T
    Es = np.vstack([E[g_][pick] for g_ in need]).T
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where((Cs > 0) & (Es > 0), Cs / np.expm1(Es), np.nan)
    valid = np.sum(~np.isnan(ratio), axis=1)
    sf = np.where(valid >= 3, np.nanmedian(ratio, axis=1), np.nan)
    ok = ~np.isnan(sf)
    say(f"  implied size factors for {int(ok.sum())}/{ns} nuclei")
    if ok.sum() > 50:
        q(sf[ok], "implied sf")
        r_m, p_m = spearmanr(sf[ok], mito_pct[pick][ok])
        say(f"  rho(size factor, mito%) = {r_m:+.4f}  p={p_m:.4g}")
        r_d, p_d = spearmanr(sf[ok], total_counts[pick][ok])
        say(f"  rho(size factor, depth) = {r_d:+.4f}  p={p_d:.4g}")
        say("")
        if abs(r_m) >= 0.60:
            say(f"  |rho| = {abs(r_m):.4f} with mito%: strong coupling.")
        elif abs(r_m) >= 0.25:
            say(f"  |rho| = {abs(r_m):.4f} with mito%: moderate coupling.")
        else:
            say(f"  |rho| = {abs(r_m):.4f} with mito%: weak coupling.")
        if abs(r_d) > 0.999:
            say(f"  rho with depth is {r_d:+.4f}: size factor equals depth.")
except Exception as ex:
    say(f"  [skipped] {type(ex).__name__}: {ex}")


# =============================================================================
# Save
# =============================================================================
hdr("Saving")

inv.to_csv(os.path.join(OUT_DIR, "S1b_mito_gene_inventory.csv"), index=False)
say(f"  saved: {OUT_DIR}/S1b_mito_gene_inventory.csv")

# Column names used by this run.
say("")
say("  columns used by this run:")
for _k, _v in [("supertype", stype_col), ("region", region_col),
               ("donor", donor_col), ("ADNC", adnc_col),
               ("Braak", braak_col)]:
    say(f"    {_k:<11} {_v}")

with open(os.path.join(OUT_DIR, "S1b_precheck_report.txt"), "w",
          encoding="utf-8") as f:
    f.write("\n".join(_report))
say(f"  saved: {OUT_DIR}/S1b_precheck_report.txt")

say("")
say("=" * 74)
say("  PART 1b COMPLETE")
say("=" * 74)
