# =============================================================================
# Part1e_TechQC_code.py
# SEA-AD CEP-IP/CMRF, Part 1e, technical-confound checks
#   Tests technical sources that could produce the MT-CO3 module without the
#   biology, before Part 2 selects cells: library batch and the intronic:exonic
#   ambient index, an exogenous soup index, cross-supertype specificity at
#   matched depth, threshold-shape screens, an abundance-matched non-mito
#   panel null, doublets, post-mortem covariates, chemistry and region effects.
#
# Usage (from the SEA-AD folder, or set SEAAD_ROOT), after Part 1d
#   python3 Part1e_TechQC_code.py
#
# Outputs (Part1_PreChecks/): report, per-library, per-supertype and
#   per-region tables and screen results (CSV), figures (S1e_figures/, PNG).
# =============================================================================

import os
import gc
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import spearmanr, rankdata, mannwhitneyu
import scanpy as sc

# == Plotting =================================================================
# Agg writes PNGs without opening windows; each figure is saved before plt.show().
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({
    "figure.dpi": 110,
    "savefig.dpi": 200,
    "savefig.bbox": "tight",
    "font.size": 9,
    "axes.titlesize": 10,
    "axes.labelsize": 9,
    "axes.grid": True,
    "grid.alpha": 0.25,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "legend.frameon": False,
})

SHOW_FIGS = True        # set False for a headless / batch run
FIG_DPI = 200

# Colour scheme
CLR_TARGET = "#14606E"      # navy-teal: the target supertype / the module
CLR_OTHER = "#9E9E9E"       # neutral grey: comparison groups
CLR_NEOC = "#7C6BAF"        # muted violet: NeoC arm
CLR_NONNEOC = "#14606E"     # navy-teal: NonNeoC arm
CLR_LINE = "#333333"        # fitted lines, axes annotation
CLR_ACCENT = "#1A9850"      # green: reference lines and the ACTB marker

np.random.seed(42)
sc.settings.verbosity = 1

os.chdir(os.environ.get("SEAAD_ROOT", os.getcwd()))

FULL_H5AD = ("SEAAD_Immune_HIP_MEC_LEC_ITG_MTG_FI_STG_DFC_AnG_V1C_"
             "RNAseq_final-nuclei.2026-06-22.h5ad")
DONOR_XLSX = os.path.join("MetaData",
                          "sea-ad_cohort_donor_metadata_072524.xlsx")
OUT_DIR = "Part1_PreChecks"
os.makedirs(OUT_DIR, exist_ok=True)

# Files this script no longer writes; F0 warns if a copy is still in OUT_DIR.
SUPERSEDED_OUTPUTS = [
    ("S1e_ip_proxy.csv",
     "superseded R(x) proxy - replaced by F5e / S1e_median_specialness.csv"),
    ("S1e_P13_ip_proxy_displacement.png",
     "figure for the same superseded proxy"),
]

# == Module 1 definition (must match R_config.R) ==============================
GOI = "MT-CO3"
DFG = ["MT-ATP6", "MT-CO2", "MT-ND3", "MT-ND4",
       "MT-CYB", "MT-ND1", "MT-ND2", "MT-CO1"]
MODULE1 = [GOI] + DFG
REF_GENE = "ACTB"
TARGET_STYPE = "Micro-PVM_3-SEAAD"

# == Published CEP-IP inflection point (Part 10b) =============================
# The published IP is the IPRS R(x) = 0 crossing; Si = 0 (the median of the
# GOI) is kept as a comparator.
USE_SI0_IP      = False
IP_RX0          = {"NeoC": 4.613502, "NonNeoC": 4.732708}
IP_RX0_QUANTILE = {"NeoC": 0.5390,   "NonNeoC": 0.4839}
IP_SI0          = {"NeoC": 4.501971, "NonNeoC": 4.785060}
IP_SI0_QUANTILE = {"NeoC": 0.5000,   "NonNeoC": 0.5003}

# == Exogenous (soup) panel: genes expected to be absent from immune nuclei ===
# Their level in an immune nucleus is a per-nucleus ambient index. Lineages are
# kept separate because phagocytosis would load neuronal and oligo genes.
EXO_PANELS = {
    "neuronal":   ["SNAP25", "RBFOX3", "SYT1", "NRGN", "STMN2", "MEG3"],
    "oligo":      ["PLP1", "MBP", "MOBP", "MOG", "ST18", "CNP"],
    "astro":      ["AQP4", "GFAP", "SLC1A2", "ALDH1L1", "GJA1"],
    "opc":        ["PDGFRA", "VCAN", "CSPG4"],
    "endo":       ["CLDN5", "FLT1", "PECAM1"],
}
EXO_ALL = sorted({g for v in EXO_PANELS.values() for g in v})
# immune positive controls, to confirm the nuclei really are immune
IMMUNE_MARKERS = ["PTPRC", "CSF1R", "C1QA", "C1QB", "P2RY12", "CX3CR1"]

NEOC_REGIONS    = ["DFC", "MTG", "AnG", "STG", "V1C", "ITG"]
NONNEOC_REGIONS = ["HIP", "MEC", "LEC", "FI"]

# minimum GOI-positive target cells for a group-level rho to be reported
MIN_CELLS_GROUP = 40
N_PERM = 500            # permutations for F4

ADNC_RANK = {"Not AD": 0, "Low": 1, "Intermediate": 2, "High": 3}
BRAAK_MAP = {"Braak 0": 0, "Braak I": 1, "Braak II": 2, "Braak III": 3,
             "Braak IV": 4, "Braak V": 5, "Braak VI": 6}

FIG_DIR = os.path.join(OUT_DIR, "S1e_figures")
os.makedirs(FIG_DIR, exist_ok=True)

# Result holders, set so the plotting section runs when a section was skipped.
libdf = usable = stdf = regdf = dsum = None
obs_sd = None
null_sd = None
region_gaps = None
lib_gaps = None
gene_stats = None
AMB = None
matched = None
segdf = None
seg_curves = {}
null_rhos = None
POOL_E = {}
pool_stats = None
EXO = None
exo_by_panel = {}
ipdf = None
ip_curves = {}
MIN_N_COMPARE = 500

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


def find_col(cols, keywords):
    """Exact match first, then case-insensitive substring (shortest hit)."""
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


def safe_rho(x, y, min_n=10):
    """Spearman that returns (nan, nan) on degenerate input. Use min_n=5 for
    group-level tests (one point per supertype, region or arm)."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < min_n or len(set(x[m])) < 3 or len(set(y[m])) < 3:
        return np.nan, np.nan
    r, p = spearmanr(x[m], y[m])
    return float(r), float(p)


def partial_rho(x, y, *z, min_n=10):
    """Spearman partial correlation, returning (nan, nan) on degenerate input
    instead of raising. Thin wrapper over partial_spearman."""
    arrs = [np.asarray(a, dtype=float) for a in (x, y) + z]
    m = np.all([np.isfinite(a) for a in arrs], axis=0)
    if m.sum() < min_n:
        return np.nan, np.nan
    try:
        r, p = partial_spearman(arrs[0][m], arrs[1][m], [a[m] for a in arrs[2:]])
        return float(r), float(p)
    except Exception:
        return np.nan, np.nan


def partial_spearman(x, y, z):
    """Spearman partial correlation of x,y given covariate(s) z."""
    z = np.atleast_2d(np.asarray(z, dtype=float))
    if z.shape[0] > z.shape[1]:
        z = z.T
    R = np.column_stack([rankdata(x), rankdata(y)] + [rankdata(c) for c in z])
    sd = R.std(0)
    sd[sd == 0] = 1.0
    R = (R - R.mean(0)) / sd
    Z = np.column_stack([np.ones(len(R)), R[:, 2:]])
    ex = R[:, 0] - Z @ np.linalg.lstsq(Z, R[:, 0], rcond=None)[0]
    ey = R[:, 1] - Z @ np.linalg.lstsq(Z, R[:, 1], rcond=None)[0]
    return spearmanr(ex, ey)


def fmt(v, nd=4):
    return "   n/a  " if v is None or not np.isfinite(v) else f"{v:+.{nd}f}"


_figs_made = []


def finish_fig(fig, name, note=""):
    """Save the figure to PNG, then show it."""
    path = os.path.join(FIG_DIR, f"{name}.png")
    try:
        fig.savefig(path, dpi=FIG_DPI, bbox_inches="tight", facecolor="white")
        _figs_made.append(path)
        say(f"  [figure] {name}.png" + (f"  - {note}" if note else ""))
    except Exception as e:
        say(f"  [figure] FAILED to save {name}: {type(e).__name__}: {e}")
    if SHOW_FIGS:
        plt.show()          # renders into Spyder's Plots pane
    else:
        plt.close(fig)


def boot_rho_ci(x, y, n_boot=250, seed=0):
    """Percentile bootstrap 95% CI for Spearman rho."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    n = len(x)
    if n < 30:
        return np.nan, np.nan
    rr = np.random.default_rng(seed)
    out = []
    for _ in range(n_boot):
        idx = rr.integers(0, n, n)
        r, _p = safe_rho(x[idx], y[idx])
        if np.isfinite(r):
            out.append(r)
    if len(out) < 20:
        return np.nan, np.nan
    return float(np.quantile(out, 0.025)), float(np.quantile(out, 0.975))


def quartile_rho(values, label):
    """rho(GOI, cDFG) within quartiles of `values`, over the analysed set."""
    v = values[ANALYSED]
    g, c = GOI_E[ANALYSED], CDFG[ANALYSED]
    ok = np.isfinite(v)
    if ok.sum() < 200:
        return None
    qs = np.quantile(v[ok], [0, .25, .5, .75, 1.0])
    rows = []
    # Half-open bins [lo, hi), top bin closed, so each cell is counted once.
    for i in range(4):
        m = (ok & (v >= qs[i]) & (v < qs[i + 1])) if i < 3 else \
            (ok & (v >= qs[i]) & (v <= qs[i + 1]))
        if m.sum() > 30:
            r, p = safe_rho(g[m], c[m])
            lo, hi = boot_rho_ci(g[m], c[m], n_boot=200, seed=i)
            rows.append(dict(q=f"Q{i+1}", n=int(m.sum()),
                             median=float(np.median(v[m])),
                             rho=r, lo=lo, hi=hi))
    return pd.DataFrame(rows) if rows else None


# =============================================================================
hdr("F0 - load, resolve columns, technical-covariate inventory")

# The IP definition is printed first because F5e depends on it.
_ipdef = "Si_norm = 0 (median GOI)" if USE_SI0_IP else "IPRS R(x) = 0 crossing"
_ipcmp = "R(x) = 0 crossing" if USE_SI0_IP else "Si_norm = 0 (= median GOI)"
say(f"  PUBLISHED IP definition assumed here : {_ipdef}")
say(f"    (USE_SI0_IP = {USE_SI0_IP}; comparator = {_ipcmp})")
if not USE_SI0_IP:
    say(f"    published IP  NeoC {IP_RX0['NeoC']:.6f} at u="
        f"{IP_RX0_QUANTILE['NeoC']:.4f}   "
        f"NonNeoC {IP_RX0['NonNeoC']:.6f} at u={IP_RX0_QUANTILE['NonNeoC']:.4f}")
    say(f"    comparator    NeoC {IP_SI0['NeoC']:.6f} at u="
        f"{IP_SI0_QUANTILE['NeoC']:.4f}   "
        f"NonNeoC {IP_SI0['NonNeoC']:.6f} at u={IP_SI0_QUANTILE['NonNeoC']:.4f}")

_stale = [(f, why) for f, why in SUPERSEDED_OUTPUTS
          if os.path.exists(os.path.join(OUT_DIR, f))]
if _stale:
    say("")
    say("  Superseded outputs present in " + OUT_DIR + ":")
    for f, why in _stale:
        say(f"      {f:<40} {why}")
say("")

ad = sc.read_h5ad(FULL_H5AD, backed="r")
say(f"  file  : {FULL_H5AD}")
say(f"  shape : {ad.shape}")
_oc = list(ad.obs.columns)

C = {
    "stype":  find_col(_oc, ["Supertype", "supertype"]),
    "region": find_col(_oc, ["Brain Region", "brain_region", "Region"]),
    "donor":  find_col(_oc, ["Donor ID", "donor_id", "Donor"]),
    "lib":    find_col(_oc, ["library_prep", "library"]),
    "adnc":   find_col(_oc, ["Overall AD neuropathological Change", "ADNC"]),
    "braak":  find_col(_oc, ["Braak", "braak"]),
    "intron": find_col(_oc, ["GEX_Reads_mapped_confidently_to_intronic_regions",
                             "intronic"]),
    "exon":   find_col(_oc, ["GEX_Reads_mapped_confidently_to_exonic_regions",
                             "exonic"]),
    "interg": find_col(_oc, ["GEX_Reads_mapped_confidently_to_intergenic_regions",
                             "intergenic"]),
    "sat":    find_col(_oc, ["GEX_sequencing_saturation", "saturation"]),
    "dup":    find_col(_oc, ["GEX_Percent_duplicates", "duplicates"]),
    "doublet": find_col(_oc, ["Doublet score", "doublet"]),
    "mitofrac": find_col(_oc, ["Fraction mitochondrial UMIs"]),
    "nUMI":   find_col(_oc, ["Number of UMIs"]),
    "nGene":  find_col(_oc, ["Genes detected"]),
    "method": find_col(_oc, ["method", "Method"]),
    "pmi":    find_col(_oc, ["PMI"]),
    "rin":    find_col(_oc, ["RIN"]),
    "ph":     find_col(_oc, ["Brain pH"]),
    "weight": find_col(_oc, ["Fresh Brain Weight"]),
    "severe": find_col(_oc, ["Severely Affected Donor"]),
    "sex":    find_col(_oc, ["Sex"]),
}
say("")
say("  resolved columns:")
for k, v in C.items():
    mark = "" if v is not None else "   <- ABSENT"
    say(f"    {k:<9} -> {str(v)}{mark}")

_crit = [k for k in ("stype", "region", "donor", "lib") if C[k] is None]
if _crit:
    say("")
    say(f"  cannot proceed without {_crit}. Full obs column list:")
    for i, c in enumerate(_oc):
        say(f"    [{i:>3}] {c}")
    with open(os.path.join(OUT_DIR, "S1e_techqc_report.txt"), "w",
              encoding="utf-8") as f:
        f.write("\n".join(_report))
    raise SystemExit(f"Part1e: missing critical columns {_crit}")

obs = ad.obs.copy()
var_names = list(ad.var_names)
_lk = [k for k in (ad.layers.keys() if ad.layers else []) if k is not None]
say("")
say(f"  layers: {_lk}")

# Columns expected to be constant
for _k in ("Organism", "Genome", "alignment"):
    _c = find_col(_oc, [_k])
    if _c:
        _v = sorted(set(obs[_c].astype(str)))
        say(f"  {_c:<12} distinct values: {_v[:5]}"
            + ("" if len(_v) == 1 else "   <- NOT constant, note it"))

# Three pass/fail QC flags
say("")
say("  --- pass/fail QC flags ---")
_any_fail = False
for _k in ("rna_amplification_pass_fail", "library_prep_pass_fail",
           "experiment_component_failed"):
    _c = find_col(_oc, [_k])
    if _c is None:
        say(f"  {_k:<34} absent")
        continue
    _vc = obs[_c].astype(str).str.strip().value_counts(dropna=False)
    say(f"  {_c}:")
    for _v, _n in _vc.items():
        say(f"      {str(_v):<24} {int(_n):>9,}  ({100*_n/len(obs):.3f}%)")
    _bad = obs[_c].astype(str).str.strip().str.lower().isin(
        ["fail", "failed", "true", "yes", "y", "1"])
    if _bad.any():
        _any_fail = True
        say(f"      {int(_bad.sum()):,} nuclei on a failed component")
        _rr = obs.loc[_bad.values, C["region"]].astype(str).value_counts()
        say("      by region: " + ", ".join(
            f"{k} {int(v):,}" for k, v in _rr.items()))
if _any_fail:
    say("")
    say("  >>> Nuclei from failed components are present.")
else:
    say("  >>> No failed components among these nuclei.")

try:
    ad.file.close()
except Exception:
    pass
del ad
gc.collect()


# =============================================================================
hdr("F1 - library_prep structure: is batch separable from donor and region?")

lib = obs[C["lib"]].astype(str)
don = obs[C["donor"]].astype(str)
reg = obs[C["region"]].astype(str)
sty = obs[C["stype"]].astype(str)

N_LIBS_TOTAL = int(lib.nunique())        # every library in obs (F3 quotes it)
say(f"  distinct libraries : {N_LIBS_TOTAL:,}")
say(f"  distinct donors    : {don.nunique():,}")
say(f"  distinct regions   : {reg.nunique():,}")
say("")
_n = lib.value_counts()
say(f"  nuclei per library : min {_n.min():,}  median {int(_n.median()):,}  "
    f"max {_n.max():,}")
say(f"  libraries with <50 nuclei : {int((_n < 50).sum())}")

# is a library ever shared across donors, or across regions?
_ld = obs.groupby(lib, observed=True)[C["donor"]].nunique()
_lr = obs.groupby(lib, observed=True)[C["region"]].nunique()
say("")
say(f"  libraries spanning >1 donor  : {int((_ld > 1).sum())} of {len(_ld)}")
say(f"  libraries spanning >1 region : {int((_lr > 1).sum())} of {len(_lr)}")
say("")
if int((_ld > 1).sum()) == 0:
    say("  >>> Every library belongs to one donor: library is nested in donor.")
else:
    say("  >>> Some libraries span several donors: library and donor are partly crossed.")

_dl = obs.groupby(don, observed=True)[C["lib"]].nunique()
say("")
say(f"  libraries per donor : min {_dl.min()}  median {int(_dl.median())}  "
    f"max {_dl.max()}")
if _dl.max() > 1:
    say("  >>> Some donors have more than one library (within-donor comparison in F3).")

# region x library
say("")
say("  libraries per region:")
for r in sorted(set(reg)):
    _m = (reg == r).values
    say(f"    {r:<6} {obs.loc[_m, C['lib']].nunique():>5} libraries, "
        f"{int(_m.sum()):>8,} nuclei")


# =============================================================================
hdr("F2 - intronic:exonic ambient index")

AMB = None
if C["intron"] is None or C["exon"] is None:
    say("  intronic and/or exonic GEX columns absent - no ambient index.")
else:
    _in = pd.to_numeric(obs[C["intron"]], errors="coerce")
    _ex = pd.to_numeric(obs[C["exon"]], errors="coerce")
    say(f"  {C['intron']}")
    say(f"    range {np.nanmin(_in):.6g} .. {np.nanmax(_in):.6g}")
    say(f"  {C['exon']}")
    say(f"    range {np.nanmin(_ex):.6g} .. {np.nanmax(_ex):.6g}")
    _scale = ("fraction 0-1" if np.nanmax([_in.max(), _ex.max()]) <= 1.5
              else "percent 0-100" if np.nanmax([_in.max(), _ex.max()]) <= 100
              else "raw read counts")
    say(f"  interpreted as    : {_scale}")
    with np.errstate(divide="ignore", invalid="ignore"):
        AMB = (_in / _ex).values
    _ok = np.isfinite(AMB)
    say("")
    say(f"  intronic:exonic ratio over nuclei: min {np.nanmin(AMB):.4f}  "
        f"median {np.nanmedian(AMB):.4f}  max {np.nanmax(AMB):.4f}")
    say(f"  finite values: {int(_ok.sum()):,}/{len(AMB):,}")
    say("")
    say("  A lower intronic:exonic ratio means more ambient RNA.")

    # Is the index constant within a library? It should be (per-library metric).
    _g = pd.Series(AMB).groupby(lib.values, observed=True).nunique()
    say("")
    say(f"  distinct ratio values per library: max {int(_g.max())}")
    if int(_g.max()) == 1:
        say("  >>> Constant within library (one value per library).")
    else:
        say("  >>> Not constant within library.")

    # Does the ambient index track pathology?
    if C["adnc"] is not None:
        _lib_amb = pd.Series(AMB).groupby(don.values, observed=True).median()
        _dadnc = pd.to_numeric(
            obs.groupby(don, observed=True)[C["adnc"]].first().map(ADNC_RANK),
            errors="coerce")
        _j = pd.concat([_lib_amb.rename("amb"), _dadnc.rename("adnc")],
                       axis=1).dropna()
        r, p = safe_rho(_j["amb"], _j["adnc"])
        say("")
        say(f"  donor-level ambient index vs ADNC: rho={fmt(r)}  p={p:.4g}  "
            f"n={len(_j)}")
        if np.isfinite(p) and p < 0.05:
            say("  >>> The ambient index tracks pathology.")
        else:
            say("  >>> No pathology association in the ambient index.")


# =============================================================================
hdr("Loading the gene block (one full load, then freed)")

ad = sc.read_h5ad(FULL_H5AD)
COUNT_LAYER = None
for cand in ["UMIs", "counts", "raw_counts"]:
    if ad.layers and cand in ad.layers:
        COUNT_LAYER = cand
        break
say(f"  count layer : {COUNT_LAYER}")

mt_all = sorted([g for g in var_names if g.startswith("MT-")])
need = [g for g in dict.fromkeys(MODULE1 + mt_all + [REF_GENE]
                                 + EXO_ALL + IMMUNE_MARKERS)
        if g in var_names]
need_idx = [var_names.index(g) for g in need]

Craw = ad.layers[COUNT_LAYER] if COUNT_LAYER else ad.X
total_counts = np.asarray(Craw.sum(axis=1)).ravel().astype(float)
mt_idx = [var_names.index(g) for g in mt_all]
mt_counts = np.asarray(Craw[:, mt_idx].sum(axis=1)).ravel().astype(float)
mito_pct = np.where(total_counts > 0, 100.0 * mt_counts / total_counts, np.nan)

Cblk = Craw[:, need_idx]
Cblk = np.asarray(Cblk.toarray() if sp.issparse(Cblk) else Cblk, dtype=float)
Eblk = ad.X[:, need_idx]
Eblk = np.asarray(Eblk.toarray() if sp.issparse(Eblk) else Eblk, dtype=float)
Cn = {g: Cblk[:, i] for i, g in enumerate(need)}
E = {g: Eblk[:, i] for i, g in enumerate(need)}
del Cblk, Eblk
gc.collect()
say(f"  extracted {len(need)} genes")

# == Matched-gene pool for F5d, taken in the same load =========================
# Non-mito genes matched to Module 1 on abundance and detection.
POOL_E = {}
pool_stats = None
try:
    _is_t = (sty.values == TARGET_STYPE)
    _goi_pos_tmp = (Cn[GOI] > 0)
    _an = _is_t & _goi_pos_tmp
    say(f"  building the abundance-matched gene pool on "
        f"{int(_an.sum()):,} analysed cells ...")
    Csub = ad.layers[COUNT_LAYER][_an] if COUNT_LAYER else ad.X[_an]
    if not sp.issparse(Csub):
        Csub = sp.csr_matrix(Csub)
    Csub = Csub.tocsc()
    _n = Csub.shape[0]
    _mean = np.asarray(Csub.mean(axis=0)).ravel()
    _det = np.asarray((Csub > 0).sum(axis=0)).ravel() / max(_n, 1) * 100.0
    pool_stats = pd.DataFrame({"gene": var_names,
                               "mean_counts": _mean,
                               "detect_pct": _det})
    pool_stats["is_mito"] = pool_stats.gene.str.startswith("MT-")
    # Candidates: non-mito, non-module, detected somewhere
    cand = pool_stats[(~pool_stats.is_mito) & (pool_stats.detect_pct > 1)]
    picked = set()
    MATCH_K = 150         # candidates kept per module gene
    for g in MODULE1:
        if g not in pool_stats.gene.values:
            continue
        row = pool_stats[pool_stats.gene == g].iloc[0]
        d = np.abs(cand.detect_pct - row.detect_pct) / 100.0
        m = np.abs(np.log10(cand.mean_counts + 1e-3)
                   - np.log10(row.mean_counts + 1e-3))
        score = d + m
        picked.update(cand.iloc[np.argsort(score.values)[:MATCH_K]].gene.tolist())
    picked = sorted(picked)
    if picked:
        pidx = [var_names.index(g) for g in picked]
        Pb = ad.X[:, pidx]
        Pb = np.asarray(Pb.toarray() if sp.issparse(Pb) else Pb, dtype=np.float32)
        POOL_E = {g: Pb[:, i] for i, g in enumerate(picked)}
        del Pb
    say(f"  matched-gene pool: {len(POOL_E)} candidate genes")
    del Csub
except Exception as _e:
    say(f"  [pool skipped] {type(_e).__name__}: {_e}")

del ad
gc.collect()

dfg_have = [g for g in DFG if g in E]
GOI_E = E[GOI]
CDFG = np.mean(np.vstack([E[g] for g in dfg_have]), axis=0)
is_target = (sty.values == TARGET_STYPE)
goi_pos = (Cn[GOI] > 0)
ANALYSED = is_target & goi_pos
say(f"  {TARGET_STYPE}: {int(is_target.sum()):,} nuclei, "
    f"{int(ANALYSED.sum()):,} with {GOI} > 0")

r_all, p_all = safe_rho(GOI_E[ANALYSED], CDFG[ANALYSED])
say(f"  pooled rho({GOI}, cDFG) over the analysed set = {fmt(r_all)}  "
    f"p={p_all:.4g}")


def group_rho(mask_extra=None):
    """rho(GOI, cDFG) within the analysed set, optionally further masked."""
    m = ANALYSED if mask_extra is None else (ANALYSED & mask_extra)
    if m.sum() < MIN_CELLS_GROUP:
        return np.nan, np.nan, int(m.sum())
    r, p = safe_rho(GOI_E[m], CDFG[m])
    return r, p, int(m.sum())


# =============================================================================
hdr("F2b - exogenous soup index: a within-library ambient measure")

# Foreign-lineage transcripts inside an immune nucleus give a per-nucleus soup
# index; the intronic:exonic index is constant within a library.
EXO = None
exo_by_panel = {}
try:
    say("  immune positive controls:")
    for g in IMMUNE_MARKERS:
        if g in Cn:
            say(f"    {g:<9} detected in {100*(Cn[g][ANALYSED] > 0).mean():>5.1f}% "
                f"of analysed nuclei, mean {Cn[g][ANALYSED].mean():.2f} counts")
    say("")
    say("  exogenous panels:")
    say(f"    {'panel':<10}{'genes found':>12}{'det%':>8}{'mean counts':>13}"
        f"{'% of library':>14}")
    for pname, genes in EXO_PANELS.items():
        have = [g for g in genes if g in Cn]
        if not have:
            say(f"    {pname:<10}{'0':>12}   (none of {genes} in the gene space)")
            continue
        mat = np.vstack([Cn[g] for g in have])
        tot = mat.sum(axis=0)
        exo_by_panel[pname] = tot
        say(f"    {pname:<10}{len(have):>12}"
            f"{100*(tot[ANALYSED] > 0).mean():>8.1f}"
            f"{tot[ANALYSED].mean():>13.3f}"
            f"{100*np.nanmedian(tot[ANALYSED] / np.maximum(total_counts[ANALYSED],1)):>14.4f}")
    if exo_by_panel:
        EXO_COUNTS = np.sum(np.vstack(list(exo_by_panel.values())), axis=0)
        EXO = 100.0 * EXO_COUNTS / np.maximum(total_counts, 1)
        say("")
        say(f"  combined soup index (% of a nucleus's UMIs that are foreign):")
        say(f"    median {np.nanmedian(EXO[ANALYSED]):.4f}%   "
            f"q75 {np.nanquantile(EXO[ANALYSED],.75):.4f}%   "
            f"q99 {np.nanquantile(EXO[ANALYSED],.99):.4f}%   "
            f"max {np.nanmax(EXO[ANALYSED]):.3f}%")
        say(f"    nuclei with zero foreign transcripts: "
            f"{100*(EXO_COUNTS[ANALYSED] == 0).mean():.1f}%")

        # Is it constant within a library? If not, it sees within-library variation.
        _within = pd.Series(EXO[ANALYSED]).groupby(
            lib.values[ANALYSED], observed=True).std()
        _between = pd.Series(EXO[ANALYSED]).groupby(
            lib.values[ANALYSED], observed=True).median().std()
        say("")
        say(f"    SD within libraries (median) {np.nanmedian(_within):.4f}   "
            f"SD between library medians {_between:.4f}")
        if np.nanmedian(_within) > _between:
            say("    >>> Varies more within libraries than between them.")

        say("")
        say("  --- does the module track the soup index? ---")
        r, p = safe_rho(EXO[ANALYSED], GOI_E[ANALYSED])
        say(f"    soup index vs {GOI:<12}      {fmt(r)}  p={p:.3g}")
        r, p = safe_rho(EXO[ANALYSED], CDFG[ANALYSED])
        say(f"    soup index vs cDFG            {fmt(r)}  p={p:.3g}")
        r, p = safe_rho(EXO[ANALYSED], mito_pct[ANALYSED])
        say(f"    soup index vs mito%           {fmt(r)}  p={p:.3g}")
        r, p = safe_rho(EXO[ANALYSED], total_counts[ANALYSED])
        say(f"    soup index vs depth           {fmt(r)}  p={p:.3g}")

        r0, _ = safe_rho(GOI_E[ANALYSED], CDFG[ANALYSED])
        rc, pc = partial_rho(GOI_E[ANALYSED], CDFG[ANALYSED], EXO[ANALYSED])
        say("")
        say(f"    module rho {fmt(r0)} -> {fmt(rc)} conditioned on the soup "
            f"index  (drop {r0 - rc:+.4f})")

        # The directional test, per nucleus
        say("")
        say("  --- module strength by soup quartile (per NUCLEUS) ---")
        qe = quartile_rho(EXO, "soup index")
        if qe is not None:
            say(f"    {'quartile':<9}{'n':>9}{'median soup%':>14}{'rho':>9}"
                f"{'95% CI':>20}")
            for _, rr in qe.iterrows():
                say(f"    {rr['q']:<9}{int(rr['n']):>9,}{rr['median']:>14.4f}"
                    f"{rr['rho']:>9.3f}"
                    f"   [{rr['lo']:+.3f}, {rr['hi']:+.3f}]")
            _lo, _hi = float(qe.iloc[0].rho), float(qe.iloc[-1].rho)
            say("")
            say(f"    Q1 (least soup) {_lo:+.3f}  ->  Q4 (most soup) {_hi:+.3f}"
                f"   change {_hi - _lo:+.3f}")
            say("")
            if _hi - _lo > 0.05:
                say("    >>> The module is stronger in the soupiest nuclei.")
            elif _lo - _hi > 0.05:
                say("    >>> The module is weaker in the soupiest nuclei.")
            else:
                say("    >>> Flat across the soup range.")

        # Lineage consistency: soup should load every lineage equally
        if len(exo_by_panel) >= 3:
            say("")
            say("  --- lineage consistency (soup vs phagocytosis) ---")
            say(f"    {'panel':<10}{'rho with cDFG':>16}{'rho with mito%':>16}")
            for pname, tot in exo_by_panel.items():
                frac = 100.0 * tot / np.maximum(total_counts, 1)
                r1, _ = safe_rho(frac[ANALYSED], CDFG[ANALYSED])
                r2, _ = safe_rho(frac[ANALYSED], mito_pct[ANALYSED])
                say(f"    {pname:<10}{fmt(r1,3):>16}{fmt(r2,3):>16}")
    else:
        say("  no exogenous genes found in the gene space - index unavailable")
except Exception as _e:
    say(f"  [skipped] {type(_e).__name__}: {_e}")


# =============================================================================
hdr("F3 - module strength per library vs the ambient index")

rows = []
for L in sorted(set(lib[ANALYSED])):
    m = (lib.values == L)
    r, p, n = group_rho(m)
    _mm = mito_pct[ANALYSED & m]
    rows.append(dict(
        library=L, n_analysed=n,
        rho=r, p=p,
        amb_index=(np.nanmedian(AMB[m]) if AMB is not None else np.nan),
        donor=don[m].iloc[0] if m.sum() else "",
        region="|".join(sorted(set(reg[m]))),
        median_depth=float(np.nanmedian(total_counts[ANALYSED & m]))
        if (ANALYSED & m).sum() else np.nan,
        median_mito=float(np.nanmedian(_mm)) if _mm.size else np.nan,
        median_goi=float(np.nanmedian(GOI_E[ANALYSED & m]))
        if (ANALYSED & m).sum() else np.nan,
    ))
libdf = pd.DataFrame(rows)
libdf.to_csv(os.path.join(OUT_DIR, "S1e_per_library.csv"), index=False)
usable = libdf[libdf["rho"].notna()]
say(f"  libraries with >={MIN_CELLS_GROUP} analysed cells: {len(usable)} of "
    f"{len(libdf)}")
# The three library counts used in this report, stated together.
say(f"    denominators: {N_LIBS_TOTAL:,} libraries exist in obs; "
    f"{len(libdf):,} contribute >=1 analysed nucleus;")
say(f"                  {len(usable):,} clear the >={MIN_CELLS_GROUP}-cell "
    f"floor and carry a per-library rho.")
say(f"  saved: {OUT_DIR}/S1e_per_library.csv")

if len(usable) >= 8:
    say("")
    say(f"  per-library rho: min {usable.rho.min():+.4f}  "
        f"q25 {usable.rho.quantile(.25):+.4f}  "
        f"median {usable.rho.median():+.4f}  "
        f"q75 {usable.rho.quantile(.75):+.4f}  "
        f"max {usable.rho.max():+.4f}")
    say(f"  SD across libraries : {usable.rho.std():.4f}")
    say(f"  libraries with rho < 0.3 : {int((usable.rho < 0.3).sum())}")
    say(f"  libraries with rho < 0   : {int((usable.rho < 0).sum())}")

    if AMB is not None and usable["amb_index"].notna().sum() >= 8:
        d = usable.dropna(subset=["amb_index", "rho"]).copy()
        r, p = safe_rho(d["amb_index"], d["rho"], min_n=8)
        say("")
        say("  --- the directional ambient test (RAW) ---")
        say(f"  rho( intronic:exonic , module strength ) = {fmt(r)}  "
            f"p={p:.4g}  n={len(d)}")

        # The ratio also tracks mito fraction at library level, so both
        # directions are tested before a verdict.
        say("")
        say("  --- mediation: ambient, or mito fraction, or depth? ---")
        r_am, p_am = safe_rho(d["amb_index"], d["median_mito"], min_n=8)
        r_ad, p_ad = safe_rho(d["amb_index"], d["median_depth"], min_n=8)
        r_md, p_md = safe_rho(d["median_mito"], d["rho"], min_n=8)
        r_dd, p_dd = safe_rho(d["median_depth"], d["rho"], min_n=8)
        say(f"    ambient index vs library mito%     {fmt(r_am)}  p={p_am:.3g}")
        say(f"    ambient index vs library depth     {fmt(r_ad)}  p={p_ad:.3g}")
        say(f"    library mito% vs module strength   {fmt(r_md)}  p={p_md:.3g}")
        say(f"    library depth vs module strength   {fmt(r_dd)}  p={p_dd:.3g}")
        say("")
        r_a_m, p_a_m = partial_rho(d["amb_index"], d["rho"],
                                   d["median_mito"], min_n=8)
        r_m_a, p_m_a = partial_rho(d["median_mito"], d["rho"],
                                   d["amb_index"], min_n=8)
        r_a_d, p_a_d = partial_rho(d["amb_index"], d["rho"],
                                   d["median_depth"], min_n=8)
        r_a_md, p_a_md = partial_rho(d["amb_index"], d["rho"],
                                     d["median_mito"], d["median_depth"],
                                     min_n=8)
        say(f"    ambient -> module | mito%          {fmt(r_a_m)}  p={p_a_m:.3g}")
        say(f"    ambient -> module | depth          {fmt(r_a_d)}  p={p_a_d:.3g}")
        say(f"    ambient -> module | mito% + depth  {fmt(r_a_md)}  p={p_a_md:.3g}")
        say(f"    mito%   -> module | ambient        {fmt(r_m_a)}  p={p_m_a:.3g}")
        AMB_SURVIVES = (np.isfinite(r_a_md) and abs(r_a_md) > 0.15
                        and np.isfinite(p_a_md) and p_a_md < 0.05)
        MITO_SURVIVES = (np.isfinite(r_m_a) and abs(r_m_a) > 0.15
                         and np.isfinite(p_m_a) and p_m_a < 0.05)

        say("")
        say("  --- sensitivity to minimum library size ---")
        say(f"    {'min cells':>10}{'libraries':>11}{'rho':>10}{'p':>12}")
        for cut in (MIN_CELLS_GROUP, 60, 80, 120, 200):
            s_ = d[d.n_analysed >= cut]
            if len(s_) >= 8:
                rc, pc = safe_rho(s_["amb_index"], s_["rho"], min_n=8)
                say(f"    {cut:>10}{len(s_):>11}{fmt(rc,3):>10}{pc:>12.3g}")
            else:
                say(f"    {cut:>10}{len(s_):>11}       too few")

        spread = np.nan
        try:
            qq = pd.qcut(d["amb_index"], 4,
                         labels=["Q1 most ambient", "Q2", "Q3",
                                 "Q4 least ambient"])
            say("")
            say("  --- effect size across the ambient range ---")
            say(f"    {'quartile':<18}{'libs':>6}{'med rho':>10}"
                f"{'med mito%':>11}{'med depth':>11}")
            for k, g in d.groupby(qq, observed=True):
                say(f"    {str(k):<18}{len(g):>6}{g.rho.median():>10.4f}"
                    f"{g.median_mito.median():>11.3f}"
                    f"{g.median_depth.median():>11.0f}")
            spread = float(d.loc[qq == "Q1 most ambient", "rho"].median()
                           - d.loc[qq == "Q4 least ambient", "rho"].median())
            say(f"    spread across the range : {spread:+.4f}")
            say(f"    libraries below rho 0.50: "
                f"{int((d.rho < 0.5).sum())} of {len(d)};  "
                f"lowest {d.rho.min():+.4f}")
        except Exception:
            pass

        say("")
        if not np.isfinite(r):
            say("  >>> not estimable")
        elif AMB_SURVIVES:
            say("  >>> Ambient association survives conditioning on mito% and depth.")
        elif MITO_SURVIVES and abs(r) > 0.15:
            say("  >>> Mediated by mito fraction, not by the ambient index.")
        elif abs(r) < 0.15:
            say("  >>> No association of consequence.")
        else:
            say("  >>> Raw association present; neither candidate survives conditioning.")
        if np.isfinite(spread):
            say("")
            say(f"      EFFECT SIZE: median per-library rho moves {spread:+.3f}")
            say("      across the range of the index.")

        # Cell-level version, conditioning the module on the ambient index
        _m = ANALYSED & np.isfinite(AMB)
        if _m.sum() > 200:
            r0, _ = safe_rho(GOI_E[_m], CDFG[_m])
            rc, pc = partial_spearman(GOI_E[_m], CDFG[_m], AMB[_m])
            say("")
            say(f"  cell-level: rho = {fmt(r0)} unconditional, {fmt(rc)} "
                f"conditioned on the ambient index  (drop {r0 - rc:+.4f})")

    # Within-donor library comparison
    if _dl.max() > 1:
        multi = [d for d in set(don[ANALYSED]) if _dl.get(d, 0) > 1]
        spreads = []
        for d in multi:
            sub = usable[usable["donor"] == d]
            if len(sub) >= 2:
                spreads.append(sub.rho.max() - sub.rho.min())
        if spreads:
            say("")
            say(f"  --- within-donor library spread ({len(spreads)} donors with "
                f">=2 usable libraries) ---")
            say(f"  rho range within a donor: median {np.median(spreads):.4f}  "
                f"max {np.max(spreads):.4f}")
else:
    say("  too few usable libraries for a per-library analysis")


# =============================================================================
hdr("F4 - permutation null: is between-library variance more than chance?")

if len(usable) >= 10:
    obs_sd = float(usable.rho.std())
    idx = np.where(ANALYSED)[0]
    sizes = usable.set_index("library")["n_analysed"].to_dict()
    lib_of = lib.values
    null_sd = []
    for _i in range(N_PERM):
        perm = np.random.permutation(idx)
        pos = 0
        rr = []
        for L, k in sizes.items():
            take = perm[pos:pos + k]
            pos += k
            if len(take) >= MIN_CELLS_GROUP:
                r, _ = safe_rho(GOI_E[take], CDFG[take])
                if np.isfinite(r):
                    rr.append(r)
        if len(rr) >= 5:
            null_sd.append(np.std(rr))
    if null_sd:
        null_sd = np.array(null_sd)
        pval = float((null_sd >= obs_sd).mean())
        say(f"  observed SD of per-library rho : {obs_sd:.4f}")
        say(f"  null SD (labels shuffled, {len(null_sd)} perms): "
            f"median {np.median(null_sd):.4f}  q95 {np.quantile(null_sd,.95):.4f}")
        say(f"  permutation p (observed >= null) : {pval:.4f}")
        say("")
        if pval < 0.05:
            say("  >>> Between-library variation exceeds chance.")
        else:
            say("  >>> Between-library variation is within chance.")
else:
    say("  too few usable libraries for a permutation test")


# =============================================================================
hdr("F5 - cross-supertype specificity (ambient is a droplet property)")

rows = []
for S in sorted(set(sty)):
    m = (sty.values == S) & goi_pos
    n = int(m.sum())
    if n < MIN_CELLS_GROUP:
        rows.append(dict(supertype=S, n_total=int((sty.values == S).sum()),
                         n_goi_pos=n, rho=np.nan, p=np.nan,
                         goi_detect_pct=np.nan, median_depth=np.nan,
                         median_mito=np.nan))
        continue
    r, p = safe_rho(GOI_E[m], CDFG[m])
    rows.append(dict(
        supertype=S,
        n_total=int((sty.values == S).sum()),
        n_goi_pos=n,
        rho=r, p=p,
        goi_detect_pct=100.0 * goi_pos[(sty.values == S)].mean(),
        median_depth=float(np.nanmedian(total_counts[m])),
        median_mito=float(np.nanmedian(mito_pct[m])),
    ))
stdf = pd.DataFrame(rows).sort_values("n_total", ascending=False)
stdf.to_csv(os.path.join(OUT_DIR, "S1e_per_supertype.csv"), index=False)
say(f"  {'supertype':<26}{'n_total':>9}{'n_GOI+':>8}{'rho':>9}"
    f"{'det%':>8}{'depth':>8}{'mito%':>7}")
for _, a in stdf.iterrows():
    say(f"  {a.supertype:<26}{int(a.n_total):>9,}{int(a.n_goi_pos):>8,}"
        f"{fmt(a.rho, 3):>9}{a.goi_detect_pct:>8.1f}"
        f"{(a.median_depth if np.isfinite(a.median_depth) else 0):>8.0f}"
        f"{(a.median_mito if np.isfinite(a.median_mito) else np.nan):>7.2f}")
say(f"  saved: {OUT_DIR}/S1e_per_supertype.csv")

# Comparison supertypes need at least MIN_N_COMPARE GOI-positive cells.
MIN_N_COMPARE = 500
_tgt = stdf[stdf.supertype == TARGET_STYPE]
_oth_all = stdf[(stdf.supertype != TARGET_STYPE) & stdf.rho.notna()]
_oth = _oth_all[_oth_all.n_goi_pos >= MIN_N_COMPARE]
_dropped = _oth_all[_oth_all.n_goi_pos < MIN_N_COMPARE]
if len(_dropped):
    say("")
    say(f"  excluded from the comparison (n < {MIN_N_COMPARE} GOI+ cells): "
        + ", ".join(f"{a.supertype} (n={int(a.n_goi_pos)}, rho {a.rho:+.3f})"
                    for _, a in _dropped.iterrows()))

# Is the supertype ordering explained by detection rate? (min_n=5: few supertypes)
if len(_oth_all) >= 4:
    _big = stdf[stdf.n_goi_pos >= MIN_N_COMPARE]
    _rd, _pd2 = safe_rho(_big.goi_detect_pct, _big.rho, min_n=5)
    _rz, _pz = safe_rho(_big.median_depth, _big.rho, min_n=5)
    say("")
    say(f"  across the {len(_big)} supertypes with n >= {MIN_N_COMPARE}:")
    say(f"    module rho vs GOI detection%  {fmt(_rd)}  p={_pd2:.3g}")
    say(f"    module rho vs median depth    {fmt(_rz)}  p={_pz:.3g}")
    say(f"    (n = {len(_big)} supertypes. With this n a correlation of even\n"
        f"     0.7 is not significant at 0.05, so treat the detectability\n"
        f"     ordering as SUGGESTIVE and let the per-gene R2 0.99 and the\n"
        f"     F5d matched-panel null carry the argument.)")
    if np.isfinite(_rd) and _rd > 0.6:
        say("    >>> The supertype ordering follows the detection-rate ordering.")

if len(_tgt) and len(_oth):
    _t = float(_tgt.rho.iloc[0])
    say("")
    say(f"  target rho {_t:+.4f} vs other supertypes (n >= {MIN_N_COMPARE}): "
        f"min {_oth.rho.min():+.4f}  median {_oth.rho.median():+.4f}  "
        f"max {_oth.rho.max():+.4f}")
    say(f"  best other is {_oth.loc[_oth.rho.idxmax(), 'supertype']} at "
        f"{_oth.rho.max():+.4f} (n={int(_oth.loc[_oth.rho.idxmax(), 'n_goi_pos']):,})")
    # Flag a comparator that beats the target only on a small, deep group.
    _bi = _oth.rho.idxmax()
    _bn = int(_oth.loc[_bi, "n_goi_pos"])
    _bd = float(_oth.loc[_bi, "median_depth"])
    _td = float(_tgt.median_depth.iloc[0])
    if _bn < 2000 or _bd > 1.3 * _td:
        say(f"    CAVEAT: that comparator has only {_bn:,} GOI+ cells and a\n"
            f"    median depth of {_bd:,.0f} against the target's {_td:,.0f}.\n"
            f"    Both inflate rho, so the raw ordering is fragile. F5b\n"
            f"    repeats it at matched depth and is the version to quote.")
    say("")
    if _t - _oth.rho.max() > 0.15:
        say("  >>> The target supertype is clearly above every other.")
    elif _t - _oth.rho.median() > 0.10:
        say("  >>> The target is above the median but overlaps the top of the others.")
    else:
        say("  >>> The module is not specific to the target supertype.")


# =============================================================================
hdr("F5b - cross-supertype comparison at matched sequencing depth")

# rho tracks detectability, so the comparison is repeated by direct
# standardisation onto the target's depth distribution.
matched = None
try:
    tgt_mask = ANALYSED
    edges = np.quantile(total_counts[tgt_mask], [0, .2, .4, .6, .8, 1.0])
    edges[-1] = np.inf
    tgt_w = []
    for i in range(5):
        tgt_w.append(float(((total_counts[tgt_mask] >= edges[i]) &
                            (total_counts[tgt_mask] < edges[i + 1])).mean()))
    tgt_w = np.array(tgt_w)
    say(f"  depth strata from the target's quintiles: "
        + ", ".join(f"{edges[i]:,.0f}" for i in range(5)))
    say("")
    rows = []
    for S in sorted(set(sty)):
        m_all = (sty.values == S) & goi_pos
        if m_all.sum() < MIN_N_COMPARE:
            continue
        strat_rho, strat_n, used_w = [], [], []
        for i in range(5):
            mb = m_all & (total_counts >= edges[i]) & (total_counts < edges[i + 1])
            if mb.sum() >= 40:
                rb, _ = safe_rho(GOI_E[mb], CDFG[mb])
                if np.isfinite(rb):
                    strat_rho.append(rb)
                    strat_n.append(int(mb.sum()))
                    used_w.append(tgt_w[i])
        if len(strat_rho) >= 3:
            w = np.array(used_w) / np.sum(used_w)
            rows.append(dict(
                supertype=S, n=int(m_all.sum()),
                strata_used=len(strat_rho),
                rho_raw=float(stdf.loc[stdf.supertype == S, "rho"].iloc[0]),
                rho_depth_matched=float(np.sum(w * np.array(strat_rho))),
                min_stratum_n=int(min(strat_n))))
    matched = pd.DataFrame(rows)
    if len(matched):
        matched["shift"] = matched.rho_depth_matched - matched.rho_raw
        matched = matched.sort_values("rho_depth_matched", ascending=False)
        say(f"  {'supertype':<24}{'n':>9}{'raw':>9}{'matched':>10}"
            f"{'shift':>9}{'strata':>8}")
        for _, a in matched.iterrows():
            mark = "  <- target" if a.supertype == TARGET_STYPE else ""
            say(f"  {a.supertype:<24}{int(a.n):>9,}{a.rho_raw:>9.3f}"
                f"{a.rho_depth_matched:>10.3f}{a['shift']:>+9.3f}"
                f"{int(a.strata_used):>8}{mark}")
        matched.to_csv(os.path.join(OUT_DIR, "S1e_supertype_depth_matched.csv"),
                       index=False)
        say(f"  saved: {OUT_DIR}/S1e_supertype_depth_matched.csv")
        _tm = matched[matched.supertype == TARGET_STYPE]
        _om = matched[matched.supertype != TARGET_STYPE]
        if len(_tm) and len(_om):
            gap = float(_tm.rho_depth_matched.iloc[0]
                        - _om.rho_depth_matched.max())
            say("")
            say(f"  target minus best other, AT MATCHED DEPTH: {gap:+.4f}")
            say("")
            if gap > 0.15:
                say("  >>> The target leads at matched depth.")
            elif gap > 0.03:
                say("  >>> The target leads slightly at matched depth.")
            else:
                say("  >>> The target does not lead at matched depth.")
except Exception as _e:
    say(f"  [skipped] {type(_e).__name__}: {_e}")


# =============================================================================
hdr("F5c - is a threshold shape specific to the target supertype?")

# Model-free screen for a threshold-like shape and whether it is specific: a
# two-segment fit to the binned cDFG-vs-GOI curve. It is a screen, not the
# CEP-IP inflection point.
say("  Model-free segmented screen (a shape screen, not the CEP-IP IP).")
say("")


def segmented_screen(g, c, n_bins=24):
    """Bin GOI into quantile bins, fit one line vs two, return the break."""
    if len(g) < 400:
        return None
    qs = np.unique(np.quantile(g, np.linspace(0, 1, n_bins + 1)))
    if len(qs) < 8:
        return None
    xb, yb = [], []
    for i in range(len(qs) - 1):
        m = (g >= qs[i]) & (g < qs[i + 1] if i < len(qs) - 2 else g <= qs[i + 1])
        if m.sum() >= 15:
            xb.append(float(np.mean(g[m])))
            yb.append(float(np.mean(c[m])))
    xb, yb = np.array(xb), np.array(yb)
    if len(xb) < 8:
        return None
    b1, a1 = np.polyfit(xb, yb, 1)
    rss1 = float(np.sum((yb - (a1 + b1 * xb)) ** 2))
    # A two-segment fit beats a line for any curve, so it is also compared with
    # a quadratic of similar flexibility.
    q = np.polyfit(xb, yb, 2)
    rss_q = float(np.sum((yb - np.polyval(q, xb)) ** 2))
    best = None
    for k in range(3, len(xb) - 3):
        xk = xb[k]
        left, right = xb <= xk, xb > xk
        if left.sum() < 3 or right.sum() < 3:
            continue
        bl, al = np.polyfit(xb[left], yb[left], 1)
        br, ar = np.polyfit(xb[right], yb[right], 1)
        rss2 = float(np.sum((yb[left] - (al + bl * xb[left])) ** 2)
                     + np.sum((yb[right] - (ar + br * xb[right])) ** 2))
        if best is None or rss2 < best["rss2"]:
            best = dict(brk=float(xk), rss2=rss2, slope_lo=float(bl),
                        slope_hi=float(br))
    if best is None or rss1 <= 0:
        return None
    best["rss1"] = rss1
    best["rss_quad"] = rss_q
    best["rss_drop"] = float(1 - best["rss2"] / rss1)
    best["rss_drop_vs_quad"] = float(1 - best["rss2"] / rss_q) if rss_q > 0 else np.nan
    best["slope_ratio"] = (best["slope_hi"] / best["slope_lo"]
                           if abs(best["slope_lo"]) > 1e-9 else np.nan)
    best["n_bins"] = int(len(xb))
    best["xb"], best["yb"] = xb, yb
    return best


seg_rows, seg_curves = [], {}
for S in sorted(set(sty)):
    m = (sty.values == S) & goi_pos
    if m.sum() < MIN_N_COMPARE:
        continue
    res = segmented_screen(GOI_E[m], CDFG[m])
    if res is None:
        continue
    seg_curves[S] = (res["xb"], res["yb"], res["brk"])
    seg_rows.append(dict(supertype=S, n=int(m.sum()),
                         breakpoint=res["brk"],
                         rss_drop=res["rss_drop"],
                         rss_drop_vs_quad=res["rss_drop_vs_quad"],
                         slope_lo=res["slope_lo"],
                         slope_hi=res["slope_hi"],
                         slope_ratio=res["slope_ratio"]))
segdf = pd.DataFrame(seg_rows)
if len(segdf):
    segdf = segdf.sort_values("rss_drop", ascending=False)
    say(f"  {'supertype':<24}{'n':>9}{'break':>8}{'vs line':>9}"
        f"{'vs quad':>9}{'slope lo':>10}{'slope hi':>10}")
    for _, a in segdf.iterrows():
        mark = "  <- target" if a.supertype == TARGET_STYPE else ""
        say(f"  {a.supertype:<24}{int(a.n):>9,}{a.breakpoint:>8.2f}"
            f"{100*a.rss_drop:>8.1f}%{100*a.rss_drop_vs_quad:>8.1f}%"
            f"{a.slope_lo:>10.3f}{a.slope_hi:>10.3f}{mark}")
    say("")
    say("  'vs line'  = RSS removed relative to a straight line (2 params)")
    say("  'vs quad'  = RSS removed relative to a quadratic (3 params)")
    segdf.to_csv(os.path.join(OUT_DIR, "S1e_segmented_screen.csv"), index=False)
    say(f"  saved: {OUT_DIR}/S1e_segmented_screen.csv")
    _ts = segdf[segdf.supertype == TARGET_STYPE]
    _os = segdf[segdf.supertype != TARGET_STYPE]
    if len(_ts) and len(_os):
        say("")
        say(f"  target RSS drop {100*float(_ts.rss_drop.iloc[0]):.1f}% vs best "
            f"other {100*float(_os.rss_drop.max()):.1f}%")
        say(f"  target break at GOI = {float(_ts.breakpoint.iloc[0]):.3f}; "
            f"other breaks span {_os.breakpoint.min():.3f} to "
            f"{_os.breakpoint.max():.3f}")
        say("")
        _tq = float(_ts.rss_drop_vs_quad.iloc[0])
        _oq = float(_os.rss_drop_vs_quad.max())
        say(f"  against a QUADRATIC: target {100*_tq:.1f}% vs best other "
            f"{100*_oq:.1f}%")
        say("")
        if _tq < 0.10:
            say("  >>> No threshold detected in the target on this screen.")
        elif _tq - _oq > 0.10:
            say("  >>> A threshold shape is stronger in the target than in any other supertype.")
        else:
            say("  >>> A threshold shape exists but is not distinctive to the target.")
else:
    say("  no supertype had enough cells for the screen")


# =============================================================================
hdr("F5e - is the median a distinguished point on the curve?")

# The published IP is the IPRS R(x) = 0 crossing (u = 0.5390 NeoC, 0.4839
# NonNeoC); Si = 0 is a comparator at the median by construction. A two-segment
# fit is placed at every candidate quantile to ask (a) whether the median is a
# distinguished point and (b) whether the published quantile is preferred. A
# central break tends to minimise RSS for any convex curve, and the extreme
# quantiles cannot be evaluated.
say(f"  USE_SI0_IP = {USE_SI0_IP}.")
if USE_SI0_IP:
    say("  Published IP = Si_norm=0 crossing = median GOI.")
else:
    say("  Published IP = R(x)=0 crossing, estimated from the fitted curve:")
    say(f"    NeoC    {IP_RX0['NeoC']:.6f}  at u = {IP_RX0_QUANTILE['NeoC']:.4f}")
    say(f"    NonNeoC {IP_RX0['NonNeoC']:.6f}  at u = "
        f"{IP_RX0_QUANTILE['NonNeoC']:.4f}")
    say("  Comparator = Si_norm=0 = median GOI, u = 0.5000 / 0.5003.")
say("")

# Candidate grid, plus the two published IP quantiles
Q_GRID_NOMINAL = np.round(np.arange(0.10, 0.905, 0.05), 3)
Q_PUBLISHED = sorted({round(float(v), 4) for v in IP_RX0_QUANTILE.values()}) \
    if not USE_SI0_IP else []
Q_GRID = np.round(np.unique(np.concatenate(
    [Q_GRID_NOMINAL, np.array(Q_PUBLISHED, dtype=float)])), 4) \
    if Q_PUBLISHED else Q_GRID_NOMINAL
say(f"  candidate grid: {Q_GRID_NOMINAL.min():.2f} to "
    f"{Q_GRID_NOMINAL.max():.2f} step 0.05 "
    f"({len(Q_GRID_NOMINAL)} points)"
    + (", plus the published IP quantiles "
       + ", ".join(f"{q:.4f}" for q in Q_PUBLISHED)
       + f" -> {len(Q_GRID)} candidates offered" if Q_PUBLISHED else ""))
say("")


def median_specialness(g, c, n_bins=30, extra_q=None):
    """Place a two-segment fit at every candidate quantile; report where the
    data prefer the break, and how the median ranks among the candidates."""
    extra_q = extra_q or {}
    if len(g) < 500:
        return None
    qs = np.unique(np.quantile(g, np.linspace(0, 1, n_bins + 1)))
    if len(qs) < 12:
        return None
    xb, yb = [], []
    for i in range(len(qs) - 1):
        m = (g >= qs[i]) & (g < qs[i + 1] if i < len(qs) - 2 else g <= qs[i + 1])
        if m.sum() >= 15:
            xb.append(float(np.mean(g[m])))
            yb.append(float(np.mean(c[m])))
    xb, yb = np.array(xb), np.array(yb)
    if len(xb) < 12:
        return None
    prof = []
    for q in Q_GRID:
        xq = float(np.quantile(g, q))
        lo, hi = xb <= xq, xb > xq
        if lo.sum() < 4 or hi.sum() < 4:
            prof.append(dict(q=q, xq=xq, rss=np.nan, dslope=np.nan))
            continue
        bl, al = np.polyfit(xb[lo], yb[lo], 1)
        br, ar = np.polyfit(xb[hi], yb[hi], 1)
        rss = float(np.sum((yb[lo] - (al + bl * xb[lo])) ** 2)
                    + np.sum((yb[hi] - (ar + br * xb[hi])) ** 2))
        prof.append(dict(q=q, xq=xq, rss=rss, dslope=float(br - bl)))
    pf = pd.DataFrame(prof).dropna(subset=["rss"])
    if len(pf) < 6 or 0.50 not in set(pf.q):
        return None
    q_best = float(pf.loc[pf.rss.idxmin(), "q"])
    rss_best = float(pf.rss.min())
    row50 = pf[pf.q == 0.50].iloc[0]
    # Rank of the median among candidate break points (1 = best)
    rank50 = int((pf.rss < row50.rss).sum()) + 1
    # How much worse the median is than the best, relative to the spread
    spread = float(pf.rss.max() - pf.rss.min())
    penalty = ((row50.rss - rss_best) / spread) if spread > 0 else 0.0
    # Second penalty against the central candidates only (q 0.25-0.75)
    pc = pf[(pf.q >= 0.25) & (pf.q <= 0.75)]
    spread_c = float(pc.rss.max() - pc.rss.min()) if len(pc) >= 3 else 0.0
    penalty_c = ((row50.rss - float(pc.rss.min())) / spread_c) \
        if spread_c > 0 else np.nan
    # Score the published quantiles as well
    pub = {}
    for _lab, _q in extra_q.items():
        _near = pf.iloc[(pf.q - _q).abs().values.argmin()]
        pub[_lab] = dict(q=float(_near.q),
                         rank=int((pf.rss < _near.rss).sum()) + 1,
                         penalty=float((_near.rss - rss_best) / spread)
                         if spread > 0 else 0.0)
    return dict(median_GOI=float(np.median(g)), q_best=q_best,
                x_best=float(pf.loc[pf.rss.idxmin(), "xq"]),
                rank_of_median=rank50, n_candidates=int(len(pf)),
                n_offered=int(len(Q_GRID)),
                q_eval_lo=float(pf.q.min()), q_eval_hi=float(pf.q.max()),
                rss_penalty=float(penalty),
                rss_penalty_central=float(penalty_c),
                published=pub,
                dslope_at_median=float(row50.dslope),
                dslope_at_best=float(pf.loc[pf.rss.idxmin(), "dslope"]),
                profile=pf)


ip_rows, ip_curves = [], {}
for S in sorted(set(sty)):
    m = (sty.values == S) & goi_pos
    if m.sum() < MIN_N_COMPARE:
        continue
    res = median_specialness(GOI_E[m], CDFG[m],
                             extra_q=({} if USE_SI0_IP else
                                      dict(IP_RX0_QUANTILE)))
    if res is None:
        continue
    ip_curves[S] = res["profile"]
    ip_rows.append(dict(supertype=S, n=int(m.sum()),
                        median_GOI=res["median_GOI"],
                        q_preferred=res["q_best"],
                        x_preferred=res["x_best"],
                        rank_of_median=res["rank_of_median"],
                        n_candidates=res["n_candidates"],
                        n_offered=res["n_offered"],
                        q_eval_lo=res["q_eval_lo"],
                        q_eval_hi=res["q_eval_hi"],
                        rss_penalty=res["rss_penalty"],
                        rss_penalty_central=res["rss_penalty_central"],
                        **{f"rank_at_u_{k}": v["rank"]
                           for k, v in res["published"].items()},
                        **{f"penalty_at_u_{k}": v["penalty"]
                           for k, v in res["published"].items()},
                        dslope_at_median=res["dslope_at_median"]))
ipdf = pd.DataFrame(ip_rows)
if len(ipdf):
    ipdf = ipdf.sort_values("rank_of_median")
    say(f"  {'supertype':<24}{'n':>8}{'median':>9}{'q preferred':>13}"
        f"{'median rank':>13}{'penalty':>10}")
    for _, a in ipdf.iterrows():
        mark = "  <- target" if a.supertype == TARGET_STYPE else ""
        say(f"  {a.supertype:<24}{int(a.n):>8,}{a.median_GOI:>9.3f}"
            f"{a.q_preferred:>13.2f}"
            f"{int(a.rank_of_median):>7} of {int(a.n_candidates):<4}"
            f"{a.rss_penalty:>10.3f}{mark}")
    ipdf.to_csv(os.path.join(OUT_DIR, "S1e_median_specialness.csv"),
                index=False)
    say(f"  saved: {OUT_DIR}/S1e_median_specialness.csv")
    say("")
    say("  'q preferred' = the quantile where a two-segment fit works best.")
    say("  'median rank' = where q=0.50 places among the candidates, 1 = best.")
    say("  'penalty'     = how much worse q=0.50 is than the best break,")
    say("                  scaled 0 (as good as the best) to 1 (the worst).")
    say("  'central'     = the same penalty scaled over q = 0.25-0.75 only.")
    say("")
    say(f"    {'supertype':<24}{'penalty':>9}{'central':>9}"
        f"{'evaluated':>12}{'range':>14}")
    for _, a in ipdf.iterrows():
        _pc = a.get("rss_penalty_central", np.nan)
        say(f"    {a.supertype:<24}{a.rss_penalty:>9.3f}"
            f"{(_pc if np.isfinite(_pc) else np.nan):>9.3f}"
            f"{int(a.n_candidates):>7} of {int(a.n_offered):<4}"
            f"{a.q_eval_lo:>7.2f}-{a.q_eval_hi:<6.2f}")

    # Does the published quantile do well, not just the median?
    _pubcols = [c for c in ipdf.columns if c.startswith("rank_at_u_")]
    if _pubcols:
        say("")
        say("  --- the published IP quantile, scored on the same grid ---")
        say(f"    {'supertype':<24}" + "".join(
            f"{c.replace('rank_at_u_', 'rank ' ):>14}" for c in _pubcols))
        for _, a in ipdf.iterrows():
            say(f"    {a.supertype:<24}" + "".join(
                f"{int(a[c]):>7} of {int(a.n_candidates):<4}"
                for c in _pubcols))
    _t = ipdf[ipdf.supertype == TARGET_STYPE]
    if len(_t):
        a = _t.iloc[0]
        say("")
        _lbl = "published IP" if USE_SI0_IP else "comparator (Si=0)"
        say(f"  TARGET: data prefer a break at q = {a.q_preferred:.2f} "
            f"(GOI = {a.x_preferred:.3f}); the {_lbl} sits at q = 0.50 "
            f"(GOI = {a.median_GOI:.3f})")
        say(f"          q=0.50 ranks {int(a.rank_of_median)} of "
            f"{int(a.n_candidates)}, penalty {a.rss_penalty:.3f}")
        say("")
        if a.rank_of_median <= 2 and a.rss_penalty < 0.15:
            say("  >>> The data prefer a break at the median.")
        elif a.rss_penalty > 0.5:
            say("  >>> The data prefer a break well away from the median.")
        else:
            say("  >>> Intermediate: the median is not the best break but not a bad one.")
    _o = ipdf[ipdf.supertype != TARGET_STYPE]
    if len(_t) and len(_o):
        say("")
        say(f"  comparators prefer q = {_o.q_preferred.min():.2f} to "
            f"{_o.q_preferred.max():.2f}; their median ranks run "
            f"{int(_o.rank_of_median.min())} to "
            f"{int(_o.rank_of_median.max())}")
        # Count supertypes whose preferred q falls in 0.45-0.55
        _inband = int(((ipdf.q_preferred >= 0.45)
                       & (ipdf.q_preferred <= 0.55)).sum())
        say(f"  across ALL {len(ipdf)} supertypes the preferred q spans "
            f"{ipdf.q_preferred.min():.2f} to {ipdf.q_preferred.max():.2f}; "
            f"{_inband} of {len(ipdf)} fall in 0.45-0.55")
        if _inband < len(ipdf):
            say("  supertypes outside 0.45-0.55: " + ", ".join(
                f"{r.supertype} ({r.q_preferred:.2f})"
                for _, r in ipdf.iterrows()
                if not (0.45 <= r.q_preferred <= 0.55)))
else:
    say("  no supertype had enough cells for the test")


# =============================================================================
hdr("F5d - abundance-matched non-mito panel null")

# Random 9-gene panels of non-mitochondrial genes matched to Module 1 on mean
# counts and detection rate: how often does a panel reach the observed rho?
N_NULL_PANELS = 300
null_rhos = None
try:
    if not POOL_E or pool_stats is None:
        say("  matched-gene pool unavailable - skipped")
    else:
        ps = pool_stats.set_index("gene")
        say(f"  pool: {len(POOL_E)} non-mitochondrial candidate genes")
        say(f"  {'module gene':<10}{'mean':>9}{'det%':>8}   "
            f"{'matched pool mean':>18}{'det%':>8}")
        # The matching window widens until each module gene has candidates; the
        # tolerance used is reported.
        TOLS = [(8, 0.30), (15, 0.50), (25, 0.80), (40, 1.20)]
        matched_for, tol_used = {}, {}
        for g in MODULE1:
            if g not in ps.index:
                continue
            row = ps.loc[g]
            cands = []
            for dt, lt in TOLS:
                cands = [c for c in POOL_E
                         if abs(ps.loc[c, "detect_pct"] - row.detect_pct) < dt
                         and abs(np.log10(ps.loc[c, "mean_counts"] + 1e-3)
                                 - np.log10(row.mean_counts + 1e-3)) < lt]
                if len(cands) >= 5:
                    tol_used[g] = (dt, lt)
                    break
            matched_for[g] = cands
            if cands:
                mm = np.mean([ps.loc[c, "mean_counts"] for c in cands])
                dd = np.mean([ps.loc[c, "detect_pct"] for c in cands])
                dt, lt = tol_used.get(g, TOLS[-1])
                say(f"  {g:<10}{row.mean_counts:>9.2f}{row.detect_pct:>8.1f}   "
                    f"{mm:>18.2f}{dd:>8.1f}   n={len(cands)}  "
                    f"tol +-{dt}pp, {lt} log10")
            else:
                say(f"  {g:<10}{row.mean_counts:>9.2f}{row.detect_pct:>8.1f}   "
                    f"{'NO MATCH at any tolerance':>26}")
        _worst = max([lt for _, lt in tol_used.values()], default=np.nan)
        if np.isfinite(_worst) and _worst > 0.5:
            say("")
            say(f"  [WARN] the widest tolerance needed was {_worst} log10 units")
            say("         (about {:.0f}x in abundance)".format(10**_worst))
        usable_dfg = [g for g in DFG if matched_for.get(g)]
        if not matched_for.get(GOI) or len(usable_dfg) < 5:
            say("")
            say("  Too few matched genes to build a null.")
        else:
            rng = np.random.default_rng(42)
            vals = []
            for _ in range(N_NULL_PANELS):
                gsel = rng.choice(matched_for[GOI])
                dsel = [rng.choice(matched_for[g]) for g in usable_dfg]
                if gsel in dsel:
                    continue
                gvec = POOL_E[gsel][ANALYSED]
                cvec = np.mean(np.vstack([POOL_E[d][ANALYSED] for d in dsel]),
                               axis=0)
                rr, _ = safe_rho(gvec, cvec)
                if np.isfinite(rr):
                    vals.append(rr)
            null_rhos = np.array(vals)
            if null_rhos.size >= 30:
                # Not named 'obs', which is the global obs DataFrame
                obs_rho = r_all
                pval = float((null_rhos >= obs_rho).mean())
                say("")
                say(f"  observed Module 1 rho          : {obs_rho:+.4f}")
                say(f"  matched non-mito null ({null_rhos.size} panels):")
                say(f"    min {null_rhos.min():+.4f}   q25 "
                    f"{np.quantile(null_rhos,.25):+.4f}   median "
                    f"{np.median(null_rhos):+.4f}   q75 "
                    f"{np.quantile(null_rhos,.75):+.4f}   max "
                    f"{null_rhos.max():+.4f}")
                say(f"  panels reaching the observed rho: "
                    f"{int((null_rhos >= obs_rho).sum())}/{null_rhos.size}"
                    f"   empirical p = {pval:.4f}")
                say("")
                if pval < 0.01 and (obs_rho - np.median(null_rhos)) > 0.2:
                    say("  >>> Matched non-mito panels do not reach the observed value.")
                elif pval < 0.05:
                    say("  >>> The observed value exceeds most matched panels.")
                else:
                    say("  >>> Matched non-mito panels reach the observed value routinely.")
                pd.DataFrame({"null_rho": null_rhos}).to_csv(
                    os.path.join(OUT_DIR, "S1e_matched_panel_null.csv"),
                    index=False)
                say(f"  saved: {OUT_DIR}/S1e_matched_panel_null.csv")
            else:
                say("  too few valid panels built")
except Exception as _e:
    say(f"  [skipped] {type(_e).__name__}: {_e}")


# =============================================================================
hdr("F6 - doublet score")

if C["doublet"] is None:
    say("  no doublet-score column")
else:
    dbl = pd.to_numeric(obs[C["doublet"]], errors="coerce").values
    say(f"  all nuclei     : median {np.nanmedian(dbl):.4f}  "
        f"q99 {np.nanquantile(dbl, .99):.4f}")
    say(f"  {TARGET_STYPE:<22}: median "
        f"{np.nanmedian(dbl[is_target]):.4f}  "
        f"q99 {np.nanquantile(dbl[is_target], .99):.4f}")
    r, p = safe_rho(dbl[ANALYSED], GOI_E[ANALYSED])
    say(f"  doublet score vs {GOI}   : rho={fmt(r)}  p={p:.4g}")
    r, p = safe_rho(dbl[ANALYSED], CDFG[ANALYSED])
    say(f"  doublet score vs cDFG   : rho={fmt(r)}  p={p:.4g}")
    _m = ANALYSED & np.isfinite(dbl)
    if _m.sum() > 200:
        r0, _ = safe_rho(GOI_E[_m], CDFG[_m])
        rc, _ = partial_spearman(GOI_E[_m], CDFG[_m], dbl[_m])
        say(f"  module conditioned on doublet score: {fmt(r0)} -> {fmt(rc)}  "
            f"(drop {r0 - rc:+.4f})")
    # Module strength in high vs low doublet halves
    _fin = ANALYSED & np.isfinite(dbl)
    if _fin.sum() > 400:
        med = np.nanmedian(dbl[_fin])
        lo = _fin & (dbl <= med)
        hi = _fin & (dbl > med)
        rl, _ = safe_rho(GOI_E[lo], CDFG[lo])
        rh, _ = safe_rho(GOI_E[hi], CDFG[hi])
        say(f"  rho in low-doublet half  {fmt(rl)}  (n={int(lo.sum()):,})")
        say(f"  rho in high-doublet half {fmt(rh)}  (n={int(hi.sum()):,})")


# =============================================================================
hdr("F7 - post-mortem covariates at donor level")

dsum = pd.DataFrame({
    "donor": don.values,
    "mito": mito_pct,
    "depth": total_counts,
})
dsum = dsum.groupby("donor", observed=True).median()
# Per-donor module strength
_pd = {}
for d in sorted(set(don[ANALYSED])):
    m = (don.values == d)
    r, p, n = group_rho(m)
    if np.isfinite(r):
        _pd[d] = r
dsum["module_rho"] = pd.Series(_pd)

for key, lab in [("pmi", "PMI"), ("rin", "RIN"), ("ph", "Brain pH"),
                 ("weight", "Fresh Brain Weight")]:
    if C[key] is None:
        say(f"  {lab:<20} column absent")
        continue
    v = obs.groupby(don, observed=True)[C[key]].first()
    v = pd.to_numeric(v, errors="coerce")
    dsum[key] = v
    for yv in ("mito", "module_rho"):
        j = dsum[[key, yv]].dropna()
        r, p = safe_rho(j[key], j[yv])
        flag = "  <-- ALTERNATIVE EXPLANATION" if (
            np.isfinite(p) and p < 0.05 and yv == "mito") else ""
        say(f"  {lab:<20} vs {yv:<11}: rho={fmt(r)}  p={p:.4g}  "
            f"n={len(j)}{flag}")

if C["adnc"] is not None:
    # pd.to_numeric: the mapped Categorical has no ordering
    dsum["adnc"] = pd.to_numeric(
        obs.groupby(don, observed=True)[C["adnc"]].first().map(ADNC_RANK),
        errors="coerce")
    j = dsum[["adnc", "mito"]].dropna()
    r, p = safe_rho(j["adnc"], j["mito"])
    say("")
    say(f"  reference: ADNC vs mito       : rho={fmt(r)}  p={p:.4g}  n={len(j)}")
    # Condition on PMI, RIN and pH, each and together.
    say("")
    say("  --- does any degradation covariate explain the ADNC-mito link? ---")
    _r0, _p0 = safe_rho(dsum["adnc"], dsum["mito"])
    say(f"    ADNC vs mito, unconditioned        {fmt(_r0)}  p={_p0:.4g}")
    _avail = [k for k in ("pmi", "rin", "ph", "weight") if k in dsum.columns]
    for k in _avail:
        j = dsum[["adnc", "mito", k]].dropna()
        if len(j) > 20:
            rc, pc = partial_rho(j["adnc"], j["mito"], j[k])
            say(f"    ADNC vs mito | {k:<20}{fmt(rc)}  p={pc:.4g}  n={len(j)}")
    if len(_avail) >= 2:
        j = dsum[["adnc", "mito"] + _avail].dropna()
        if len(j) > 25:
            rc, pc = partial_rho(j["adnc"], j["mito"],
                                 *[j[k] for k in _avail])
            say(f"    ADNC vs mito | ALL ({'+'.join(_avail)})".ljust(38)
                + f"{fmt(rc)}  p={pc:.4g}  n={len(j)}")
            _drop = abs(_r0) - abs(rc) if np.isfinite(rc) else np.nan
            say("")
            if np.isfinite(_drop) and _drop > 0.20:
                say("    >>> Degradation largely explains the ADNC-mito link.")
            elif np.isfinite(_drop) and _drop > 0.08:
                say("    >>> Degradation explains part of the ADNC-mito link.")
            else:
                say("    >>> The ADNC-mito link survives adjustment for every degradation covariate.")

    # Do the degradation covariates themselves track pathology?
    say("")
    say("  --- are the covariates themselves pathology-linked? ---")
    for k in _avail:
        j = dsum[["adnc", k]].dropna()
        rc, pc = safe_rho(j["adnc"], j[k])
        flag = "  <- on the causal path" if (np.isfinite(pc) and pc < 0.05) else ""
        say(f"    ADNC vs {k:<22}{fmt(rc)}  p={pc:.4g}  n={len(j)}{flag}")


# =============================================================================
hdr("F8 - chemistry: singleome vs multiome")

if C["method"] is None:
    say("  no method column")
else:
    meth = obs[C["method"]].astype(str)
    say("  nuclei per method:")
    for k, v in meth.value_counts().items():
        say(f"    {k:<28} {int(v):>9,}")
    _lv = sorted(set(meth[ANALYSED]))
    if len(_lv) > 1:
        say("")
        for k in _lv:
            m = (meth.values == k)
            r, p, n = group_rho(m)
            say(f"    {k:<28} rho={fmt(r)}  p={p:.4g}  n={n:,}  "
                f"median depth {np.nanmedian(total_counts[ANALYSED & m]):.0f}")
        # Donors with both chemistries, if any
        _dm = obs.groupby(don, observed=True)[C["method"]].nunique()
        say(f"  donors with more than one method: {int((_dm > 1).sum())}")
    else:
        say(f"  only one method present ({_lv}); nothing to compare")


# =============================================================================
hdr("F9 - Severely Affected Donor sensitivity")

if C["severe"] is None or C["adnc"] is None:
    say("  needs both 'Severely Affected Donor' and ADNC - one is absent")
else:
    sev = obs.groupby(don, observed=True)[C["severe"]].first().astype(str)
    sev_y = set(sev[sev.str.upper().str.startswith("Y")].index)
    say(f"  donors flagged Severely Affected : {len(sev_y)}")
    j = dsum[["adnc", "mito"]].dropna()
    r_all_, p_all_ = safe_rho(j["adnc"], j["mito"])
    j2 = j.loc[[d for d in j.index if d not in sev_y]]
    r_ex, p_ex = safe_rho(j2["adnc"], j2["mito"])
    say(f"  ADNC vs mito, all donors        : rho={fmt(r_all_)}  "
        f"p={p_all_:.4g}  n={len(j)}")
    say(f"  ADNC vs mito, excluding severe  : rho={fmt(r_ex)}  "
        f"p={p_ex:.4g}  n={len(j2)}")
    say("")
    # The change is judged relative to the original value.
    _rel = (abs(r_all_ - r_ex) / abs(r_all_)
            if np.isfinite(r_all_) and abs(r_all_) > 1e-9 else np.nan)
    if np.isfinite(_rel):
        say(f"  change: {r_all_ - r_ex:+.4f} absolute, "
            f"{100*_rel:.1f}% of the original magnitude, from "
            f"{len(sev_y)} of {len(j)} donors "
            f"({100*len(sev_y)/max(len(j),1):.0f}% of the cohort)")
    if np.isfinite(r_all_) and np.isfinite(r_ex):
        if np.isfinite(_rel) and _rel > 0.20:
            say("")
            say("  >>> A large share of the association comes from the severely affected donors.")
        elif abs(r_all_ - r_ex) > 0.15:
            say("  >>> The association depends substantially on the severely affected donors.")
        else:
            say("  >>> The association is not carried by the severely affected donors.")
    # Module strength in severe vs not
    _sm = np.isin(don.values, list(sev_y))
    r1, p1, n1 = group_rho(_sm)
    r2, p2, n2 = group_rho(~_sm)
    say("")
    say(f"  module rho, severe donors      : {fmt(r1)}  n={n1:,}")
    say(f"  module rho, other donors       : {fmt(r2)}  n={n2:,}")


# =============================================================================
hdr("F10 - region-level module strength and leave-one-region-out per arm")

rows = []
for R in sorted(set(reg)):
    m = (reg.values == R)
    r, p, n = group_rho(m)
    rows.append(dict(
        region=R,
        arm=("NeoC" if R in NEOC_REGIONS else
             "NonNeoC" if R in NONNEOC_REGIONS else "OTHER"),
        n_nuclei=int(m.sum()),
        n_analysed=n, rho=r, p=p,
        n_donors=int(obs.loc[m, C["donor"]].nunique()),
        n_libs=int(obs.loc[m, C["lib"]].nunique()),
        median_depth=float(np.nanmedian(total_counts[ANALYSED & m]))
        if (ANALYSED & m).sum() else np.nan,
        median_mito=float(np.nanmedian(mito_pct[ANALYSED & m]))
        if (ANALYSED & m).sum() else np.nan,
    ))
regdf = pd.DataFrame(rows)
regdf.to_csv(os.path.join(OUT_DIR, "S1e_per_region.csv"), index=False)
say(f"  {'region':<7}{'arm':<9}{'analysed':>10}{'rho':>9}{'donors':>8}"
    f"{'libs':>6}{'depth':>8}{'mito%':>7}")
for _, a in regdf.iterrows():
    say(f"  {a.region:<7}{a.arm:<9}{int(a.n_analysed):>10,}{fmt(a.rho,3):>9}"
        f"{int(a.n_donors):>8}{int(a.n_libs):>6}"
        f"{(a.median_depth if np.isfinite(a.median_depth) else 0):>8.0f}"
        f"{(a.median_mito if np.isfinite(a.median_mito) else np.nan):>7.2f}")
say(f"  saved: {OUT_DIR}/S1e_per_region.csv")

# Does region-level rho track region-level depth?
_rg = regdf.dropna(subset=["rho", "median_depth"])
if len(_rg) >= 5:
    _rr, _rp = safe_rho(_rg.median_depth, _rg.rho, min_n=5)
    _rm, _rmp = safe_rho(_rg.median_mito, _rg.rho, min_n=5)
    say("")
    say("  --- does region-level rho track region-level DEPTH? ---")
    say(f"  region rho vs median depth : {fmt(_rr)}  p={_rp:.3g}  "
        f"n={len(_rg)} regions")
    say(f"  region rho vs median mito% : {fmt(_rm)}  p={_rmp:.3g}  "
        f"n={len(_rg)} regions")
    if _rg.median_mito.notna().sum() >= 5:
        _pr, _pp = partial_spearman(_rg.median_depth, _rg.rho, _rg.median_mito)
        say(f"  ...depth vs rho | mito%    : {fmt(_pr)}  p={_pp:.3g}")
    say("")
    say("  --- arm means at region level (the depth-vs-arm question) ---")
    say(f"    {'arm':<10}{'regions':>9}{'mean rho':>10}{'mean depth':>12}")
    for _arm, _g in _rg.groupby("arm"):
        say(f"    {_arm:<10}{len(_g):>9}{_g.rho.mean():>10.4f}"
            f"{_g.median_depth.mean():>12.0f}")

for arm, regs in (("NeoC", NEOC_REGIONS), ("NonNeoC", NONNEOC_REGIONS)):
    present = [r for r in regs if r in set(reg)]
    if not present:
        continue
    inarm = np.isin(reg.values, present)
    r_full, p_full, n_full = group_rho(inarm)
    say("")
    say(f"  {arm}: pooled rho {fmt(r_full)}  n={n_full:,}")
    say(f"    leave-one-region-out:")
    for R in present:
        m = inarm & (reg.values != R)
        r, p, n = group_rho(m)
        _sh = 100.0 * (ANALYSED & inarm & (reg.values == R)).sum() / max(
            (ANALYSED & inarm).sum(), 1)
        say(f"      without {R:<5} rho={fmt(r)}  n={n:>7,}  "
            f"(that region was {_sh:5.1f}% of the arm)")


# =============================================================================
hdr("F10b - is the NeoC / NonNeoC difference just a depth difference?")

# The arms differ in depth and region-level rho tracks depth, so the arm
# contrast is also computed at matched depth.
try:
    inN = np.isin(reg.values, NEOC_REGIONS)
    inX = np.isin(reg.values, NONNEOC_REGIONS)
    rN, _, nN = group_rho(inN)
    rX, _, nX = group_rho(inX)
    say(f"  raw:  NeoC {fmt(rN)} (n={nN:,}, median depth "
        f"{np.nanmedian(total_counts[ANALYSED & inN]):,.0f})")
    say(f"        NonNeoC {fmt(rX)} (n={nX:,}, median depth "
        f"{np.nanmedian(total_counts[ANALYSED & inX]):,.0f})")
    say(f"  raw difference: {rN - rX:+.4f}")

    # Direct standardisation onto a common depth distribution
    both = ANALYSED & (inN | inX)
    edges = np.quantile(total_counts[both], np.linspace(0, 1, 6))
    edges[-1] = np.inf
    say("")
    say(f"  {'depth stratum':<22}{'NeoC rho':>10}{'n':>9}"
        f"{'NonNeoC rho':>13}{'n':>9}")
    accN, accX, wts = [], [], []
    for i in range(5):
        band = (total_counts >= edges[i]) & (total_counts < edges[i + 1])
        mN, mX = ANALYSED & inN & band, ANALYSED & inX & band
        if mN.sum() >= 100 and mX.sum() >= 100:
            a, _ = safe_rho(GOI_E[mN], CDFG[mN])
            b, _ = safe_rho(GOI_E[mX], CDFG[mX])
            w = float((both & band).sum())
            accN.append(a); accX.append(b); wts.append(w)
            say(f"  {edges[i]:>8,.0f}-{edges[i+1]:<12,.0f}"[:22].ljust(22)
                + f"{a:>10.3f}{int(mN.sum()):>9,}{b:>13.3f}{int(mX.sum()):>9,}")
    if len(accN) >= 3:
        w = np.array(wts) / np.sum(wts)
        sN = float(np.sum(w * np.array(accN)))
        sX = float(np.sum(w * np.array(accX)))
        say("")
        say(f"  depth-standardised: NeoC {sN:+.4f}   NonNeoC {sX:+.4f}   "
            f"difference {sN - sX:+.4f}")
        say(f"  the raw difference was {rN - rX:+.4f}")
        _shrink = 1 - abs(sN - sX) / max(abs(rN - rX), 1e-9)
        say("")
        if max(abs(rN - rX), abs(sN - sX)) < 0.02:
            say("  >>> No arm difference, raw or standardised.")
        elif _shrink > 0.5:
            say("  >>> Most of the arm difference is a depth difference.")
        elif _shrink > 0.2:
            say("  >>> Part of the arm difference is depth.")
        else:
            say("  >>> The arm difference survives depth standardisation.")

    # Paired within-donor arm comparison
    say("")
    say("  --- paired within-donor arm comparison ---")
    pairs = []
    for dnr in sorted(set(don[ANALYSED])):
        md = (don.values == dnr)
        a, _, na = group_rho(md & inN)
        b, _, nb = group_rho(md & inX)
        if np.isfinite(a) and np.isfinite(b):
            pairs.append((dnr, a, b))
    if len(pairs) >= 8:
        arr = np.array([[a, b] for _, a, b in pairs])
        diff = arr[:, 0] - arr[:, 1]
        try:
            from scipy.stats import wilcoxon
            _w, _pw = wilcoxon(arr[:, 0], arr[:, 1])
        except Exception:
            _pw = np.nan
        say(f"  {len(pairs)} donors contribute analysable cells to BOTH arms")
        say(f"  median within-donor difference (NeoC - NonNeoC): "
            f"{np.median(diff):+.4f}")
        say(f"  donors with NeoC > NonNeoC: {int((diff > 0).sum())}/{len(diff)}"
            f"   Wilcoxon p = {_pw:.4g}")
    else:
        say("  too few donors contribute to both arms for a paired test")
except Exception as _e:
    say(f"  [skipped] {type(_e).__name__}: {_e}")


# =============================================================================
hdr("F11 - GOI-zero nuclei: is the analysed subset biased?")

z = is_target & (~goi_pos)
say(f"  {TARGET_STYPE}: {int(ANALYSED.sum()):,} analysed, "
    f"{int(z.sum()):,} excluded for {GOI} == 0 "
    f"({100.0*z.sum()/max(is_target.sum(),1):.2f}%)")
if z.sum() > 50:
    for nm, v in (("depth", total_counts), ("mito%", mito_pct),
                  ("cDFG", CDFG)):
        a, b = v[ANALYSED], v[z]
        a, b = a[np.isfinite(a)], b[np.isfinite(b)]
        if a.size > 20 and b.size > 20:
            try:
                _, pu = mannwhitneyu(a, b)
            except Exception:
                pu = np.nan
            say(f"    {nm:<7} analysed median {np.median(a):>10.3f}   "
                f"excluded median {np.median(b):>10.3f}   p={pu:.3g}")


# =============================================================================
hdr("F12 - is the ACTB residual stable across regions and libraries?")

say("  rho(GOI, gene) per region for ACTB and two abundance-matched mitochondrial genes.")
say("")
COMP = [g for g in ("MT-ND5", "MT-ATP8", "MT-ND4L") if g in E]
say(f"  {'region':<7}{'ACTB':>9}" + "".join(f"{g:>10}" for g in COMP)
    + f"{'gap':>9}")
_gap_rows = []
for R in sorted(set(reg)):
    m = ANALYSED & (reg.values == R)
    if m.sum() < MIN_CELLS_GROUP:
        continue
    r_actb, _ = safe_rho(GOI_E[m], E[REF_GENE][m]) if REF_GENE in E else (np.nan, np.nan)
    rr = []
    for g in COMP:
        r_, _ = safe_rho(GOI_E[m], E[g][m])
        rr.append(r_)
    gap = (np.nanmean(rr) - r_actb) if np.isfinite(r_actb) else np.nan
    _gap_rows.append(dict(region=R, arm=("NeoC" if R in NEOC_REGIONS else
                                         "NonNeoC" if R in NONNEOC_REGIONS
                                         else "OTHER"),
                          n=int(m.sum()), actb=r_actb,
                          matched_mito=float(np.nanmean(rr)), gap=gap))
    say(f"  {R:<7}{fmt(r_actb,3):>9}"
        + "".join(f"{fmt(x,3):>10}" for x in rr)
        + f"{fmt(gap,3):>9}")
gapdf = pd.DataFrame(_gap_rows)
gaps = region_gaps = np.array(
    [g for g in gapdf["gap"].values if np.isfinite(g)])
if gaps.size:
    say("")
    say(f"  gap (mean matched-mito rho minus ACTB rho): min {gaps.min():+.3f}  "
        f"median {np.median(gaps):+.3f}  max {gaps.max():+.3f}")
    # Judged relatively: same sign everywhere and no region collapsing.
    _all_pos = bool((gaps > 0).all())
    _ratio = (gaps.min() / np.median(gaps)) if np.median(gaps) > 0 else np.nan
    say(f"  regions with a positive gap : {int((gaps > 0).sum())}/{gaps.size}")
    say(f"  weakest region as a fraction of the median gap : "
        f"{_ratio:.2f}" if np.isfinite(_ratio) else "")

    # Library-level version
    lib_gaps = []
    for L in sorted(set(lib[ANALYSED])):
        m = ANALYSED & (lib.values == L)
        if m.sum() < MIN_CELLS_GROUP:
            continue
        ra, _ = safe_rho(GOI_E[m], E[REF_GENE][m]) if REF_GENE in E \
            else (np.nan, np.nan)
        rr = [safe_rho(GOI_E[m], E[g][m])[0] for g in COMP]
        if np.isfinite(ra) and np.isfinite(np.nanmean(rr)):
            lib_gaps.append(np.nanmean(rr) - ra)
    if lib_gaps:
        lib_gaps = np.array(lib_gaps)
        say(f"  libraries with a positive gap : "
            f"{int((lib_gaps > 0).sum())}/{lib_gaps.size} "
            f"({100*(lib_gaps > 0).mean():.1f}%)")
        say(f"  library-level gap: median {np.median(lib_gaps):+.3f}  "
            f"IQR {np.quantile(lib_gaps,.25):+.3f} to "
            f"{np.quantile(lib_gaps,.75):+.3f}")

    say("")
    if _all_pos and np.isfinite(_ratio) and _ratio > 0.4:
        say("  >>> The gap is positive and stable in every region.")
    elif _all_pos:
        say("  >>> The gap is positive everywhere but varies in size.")
    else:
        say("  >>> The gap reverses in at least one region.")


# =============================================================================
hdr("F13 - figures")

say(f"  PNGs -> {FIG_DIR}")

# Per-gene statistics inside the analysed set, needed by P6
gene_rows = []
for g in [REF_GENE] + [x for x in mt_all if x != GOI]:
    if g not in E:
        continue
    r, p = safe_rho(GOI_E[ANALYSED], E[g][ANALYSED])
    gene_rows.append(dict(
        gene=g,
        in_module1=g in MODULE1,
        is_mito=g.startswith("MT-"),
        mean_counts=float(Cn[g][ANALYSED].mean()),
        detect_pct=100.0 * float((Cn[g][ANALYSED] > 0).mean()),
        rho=r))
gene_stats = pd.DataFrame(gene_rows)


# == P1: the headline ambient test ============================================
try:
    if usable is not None and len(usable) >= 8 and \
            usable["amb_index"].notna().sum() >= 8:
        d = usable.dropna(subset=["amb_index", "rho"])
        fig, (axa, axm, axb) = plt.subplots(
            1, 3, figsize=(15.5, 4.6),
            gridspec_kw={"width_ratios": [1.7, 1.7, 1]})

        sizes = 12 + 120 * (d.n_analysed / d.n_analysed.max())
        # Points coloured by library mito%, which mediates the raw slope.
        sc = axa.scatter(d.amb_index, d.rho, s=sizes, c=d.median_mito,
                         cmap="viridis", alpha=0.85, edgecolor="white",
                         linewidth=0.5, zorder=3)
        cb = fig.colorbar(sc, ax=axa, pad=0.02)
        cb.set_label("library median mito%", fontsize=8)
        if len(d) >= 3:
            b, a = np.polyfit(d.amb_index, d.rho, 1)
            xs = np.linspace(d.amb_index.min(), d.amb_index.max(), 50)
            axa.plot(xs, a + b * xs, color=CLR_LINE, lw=1.8, zorder=4)
        _r, _p = safe_rho(d.amb_index, d.rho, min_n=8)
        axa.axhline(r_all, color=CLR_ACCENT, ls="--", lw=1.2, zorder=2)
        axa.set_xlabel("library intronic : exonic ratio\n"
                       "(LOWER = more ambient contamination)")
        axa.set_ylabel(f"module strength within library\nrho({GOI}, cDFG)")
        axa.set_title(f"P1  RAW: Spearman {_r:+.3f}  (p = {_p:.3g}, "
                      f"n = {len(d)})\n"
                      "ambient predicts a downward slope - but see P1b")
        axa.text(0.03, 0.05,
                 "colour = mito%\nnote the vertical\ncolour gradient",
                 transform=axa.transAxes, fontsize=7.5, va="bottom",
                 color=CLR_LINE,
                 bbox=dict(boxstyle="round,pad=0.35", fc="white",
                           ec=CLR_LINE, lw=0.5, alpha=0.85))

        # Partial-residual panel: the same relationship after removing mito%
        rr = np.column_stack([rankdata(d.amb_index), rankdata(d.rho),
                              rankdata(d.median_mito)])
        rr = (rr - rr.mean(0)) / rr.std(0)
        Zc = np.column_stack([np.ones(len(rr)), rr[:, 2]])
        ex = rr[:, 0] - Zc @ np.linalg.lstsq(Zc, rr[:, 0], rcond=None)[0]
        ey = rr[:, 1] - Zc @ np.linalg.lstsq(Zc, rr[:, 1], rcond=None)[0]
        _rp, _pp = spearmanr(ex, ey)
        axm.scatter(ex, ey, s=sizes, c=CLR_TARGET, alpha=0.55,
                    edgecolor="white", linewidth=0.5, zorder=3)
        bb, aa = np.polyfit(ex, ey, 1)
        xs2 = np.linspace(ex.min(), ex.max(), 30)
        axm.plot(xs2, aa + bb * xs2, color=CLR_LINE, lw=1.8, zorder=4)
        axm.axhline(0, color="black", lw=0.8)
        axm.set_xlabel("ambient index, residual after mito%")
        axm.set_ylabel("module strength, residual after mito%")
        axm.set_title(f"P1b  AFTER removing library mito%: {_rp:+.3f} "
                      f"(p = {_pp:.3g})\n"
                      "if this is flat, the raw slope was mito%, not ambient")

        axb.hist(d.rho, bins=max(8, int(np.sqrt(len(d)))),
                 color=CLR_TARGET, alpha=0.8, edgecolor="white")
        axb.axvline(0, color="black", lw=1.0)
        axb.axvline(r_all, color=CLR_ACCENT, ls="--", lw=1.4)
        axb.set_xlabel(f"rho({GOI}, cDFG) per library")
        axb.set_ylabel("libraries")
        axb.set_title(f"P1c  {int((d.rho < 0.5).sum())} of {len(d)} libraries\n"
                      f"below 0.5; lowest {d.rho.min():+.3f}")
        fig.tight_layout()
        finish_fig(fig, "S1e_P1_ambient_index_vs_module",
                   "the headline directional test")
    else:
        say("  [figure] P1 skipped - no usable ambient index")
except Exception as e:
    say(f"  [figure] P1 FAILED: {type(e).__name__}: {e}")


# == P2: permutation null on between-library variance =========================
try:
    if null_sd is not None and len(np.atleast_1d(null_sd)) > 5:
        ns = np.atleast_1d(null_sd)
        fig, ax = plt.subplots(figsize=(6.4, 4.2))
        ax.hist(ns, bins=25, color=CLR_OTHER, alpha=0.8, edgecolor="white",
                label=f"null: library labels shuffled ({len(ns)} perms)")
        ax.axvline(obs_sd, color=CLR_TARGET, lw=2.2,
                   label=f"observed SD {obs_sd:.4f}")
        ax.axvline(np.quantile(ns, 0.95), color=CLR_LINE, ls="--", lw=1.2,
                   label=f"null 95th pct {np.quantile(ns,0.95):.4f}")
        _pv = float((ns >= obs_sd).mean())
        ax.set_xlabel("SD of module strength across libraries")
        ax.set_ylabel("permutations")
        ax.set_title("P2  Between-library variation vs chance\n"
                     f"permutation p = {_pv:.3f}"
                     f"  ({'library effect present' if _pv < 0.05 else 'no detectable library effect'})")
        ax.legend(loc="upper right", fontsize=8)
        fig.tight_layout()
        finish_fig(fig, "S1e_P2_library_permutation_null")
    else:
        say("  [figure] P2 skipped - permutation test did not run")
except Exception as e:
    say(f"  [figure] P2 FAILED: {type(e).__name__}: {e}")


# == P3: cross-supertype specificity ==========================================
try:
    if stdf is not None and stdf["rho"].notna().sum() >= 2:
        d = stdf.dropna(subset=["rho"]).sort_values("rho")
        # Supertypes below MIN_N_COMPARE cells are drawn faded.
        low_n = d.n_goi_pos < MIN_N_COMPARE
        cols = [CLR_TARGET if s_ == TARGET_STYPE else CLR_OTHER
                for s_ in d.supertype]
        alphas = [0.35 if lo else 1.0 for lo in low_n]
        fig, (axa, axb) = plt.subplots(1, 2, figsize=(12, 4.6))

        ypos = np.arange(len(d))
        for y, (_, a), c, al in zip(ypos, d.iterrows(), cols, alphas):
            axa.barh(y, a.rho, color=c, edgecolor="white", height=0.72,
                     alpha=al, hatch=("///" if a.n_goi_pos < MIN_N_COMPARE
                                      else None))
        axa.set_yticks(ypos)
        axa.set_yticklabels(
            [f"{s_}  (n={int(n):,})" + ("  underpowered"
                                        if n < MIN_N_COMPARE else "")
             for s_, n in zip(d.supertype, d.n_goi_pos)])
        axa.axvline(0, color="black", lw=1.0)
        axa.set_xlabel(f"rho({GOI}, cDFG)")
        axa.set_title("P3  Ambient is a droplet property, not a cell-type one\n"
                      "so an ambient module would appear in every supertype")

        # rho against detectability, the known confounder
        big = d[~low_n.values]
        axb.scatter(d.goi_detect_pct[low_n.values], d.rho[low_n.values],
                    s=55, facecolor="none", edgecolor=CLR_OTHER,
                    linewidth=1.2, zorder=3, label=f"n < {MIN_N_COMPARE}")
        axb.scatter(big.goi_detect_pct, big.rho, s=80,
                    c=[CLR_TARGET if s_ == TARGET_STYPE else CLR_OTHER
                       for s_ in big.supertype],
                    edgecolor="white", linewidth=0.8, zorder=4)
        for _, a in d.iterrows():
            axb.annotate(a.supertype.replace("-SEAAD", "*"),
                         (a.goi_detect_pct, a.rho), fontsize=7,
                         xytext=(4, 4), textcoords="offset points")
        # min_n=5: only a handful of supertypes
        if len(big) >= 4:
            _rd, _pd_ = safe_rho(big.goi_detect_pct, big.rho, min_n=5)
            if len(big) >= 3:
                bb2, aa2 = np.polyfit(big.goi_detect_pct, big.rho, 1)
                xs3 = np.linspace(d.goi_detect_pct.min(),
                                  d.goi_detect_pct.max(), 30)
                axb.plot(xs3, aa2 + bb2 * xs3, color=CLR_LINE, lw=1.4,
                         zorder=2)
            axb.set_title("P3b  ...or is it just detectability?\n"
                          f"rho(detection%, module rho) = {fmt(_rd,3).strip()}"
                          f"  (p = {_pd_:.3g}, n = {len(big)} well-powered)")
            axb.legend(loc="lower right", fontsize=7.5)
        axb.set_xlabel(f"% of nuclei with {GOI} > 0")
        axb.set_ylabel(f"rho({GOI}, cDFG)")
        fig.tight_layout()
        finish_fig(fig, "S1e_P3_cross_supertype_specificity")
    else:
        say("  [figure] P3 skipped - too few supertypes with a usable rho")
except Exception as e:
    say(f"  [figure] P3 FAILED: {type(e).__name__}: {e}")


# == P4: depth and mito% quartile gradients ===================================
try:
    qd = quartile_rho(total_counts, "depth")
    qm = quartile_rho(mito_pct, "mito%")
    if qd is not None or qm is not None:
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.3), sharey=True)
        for ax, q, lab, note in (
                (axes[0], qd, "total UMI depth",
                 "ambient predicts HIGHEST rho in Q1"),
                (axes[1], qm, "mitochondrial fraction",
                 "circular - reference only")):
            if q is None:
                ax.set_visible(False)
                continue
            x = np.arange(len(q))
            yerr = np.vstack([q.rho - q.lo, q.hi - q.rho])
            ax.errorbar(x, q.rho, yerr=yerr, fmt="o-", color=CLR_TARGET,
                        capsize=4, lw=1.6, ms=7, zorder=3)
            ax.set_xticks(x)
            # Bracket access: r.median would return the Series method
            ax.set_xticklabels(
                [f"{r['q']}\nn={r['n']:,}\nmed {r['median']:,.0f}"
                 for _, r in q.iterrows()], fontsize=7)
            ax.set_xlabel(f"quartile of {lab}")
            ax.set_title(f"P4  by {lab}\n{note}")
        axes[0].set_ylabel(f"rho({GOI}, cDFG)  with 95% bootstrap CI")
        fig.tight_layout()
        finish_fig(fig, "S1e_P4_quartile_gradients",
                   "rising with depth argues against ambient")
    else:
        say("  [figure] P4 skipped - not enough cells for quartiles")
except Exception as e:
    say(f"  [figure] P4 FAILED: {type(e).__name__}: {e}")


# == P5: per-region forest plot and leave-one-region-out ======================
try:
    if regdf is not None and regdf["rho"].notna().sum() >= 3:
        d = regdf.dropna(subset=["rho"]).copy()
        ci = [boot_rho_ci(GOI_E[ANALYSED & (reg.values == R)],
                          CDFG[ANALYSED & (reg.values == R)],
                          n_boot=200, seed=i)
              for i, R in enumerate(d.region)]
        d["lo"] = [c[0] for c in ci]
        d["hi"] = [c[1] for c in ci]
        d = d.sort_values(["arm", "rho"])

        fig, (axa, axb) = plt.subplots(
            1, 2, figsize=(12, 4.8),
            gridspec_kw={"width_ratios": [1.35, 1]})

        ypos = np.arange(len(d))
        cols = [CLR_NEOC if a == "NeoC" else CLR_NONNEOC for a in d.arm]
        for y, (_, a), c in zip(ypos, d.iterrows(), cols):
            if np.isfinite(a.lo) and np.isfinite(a.hi):
                axa.plot([a.lo, a.hi], [y, y], color=c, lw=2.0, zorder=2)
            axa.plot(a.rho, y, "o", color=c, ms=8, zorder=3)
        axa.set_yticks(ypos)
        axa.set_yticklabels([f"{a.region}  {a.arm}  n={int(a.n_analysed):,}"
                             f"  {int(a.n_donors)}d/{int(a.n_libs)}L"
                             for _, a in d.iterrows()], fontsize=8)
        axa.axvline(r_all, color=CLR_ACCENT, ls="--", lw=1.2,
                    label=f"pooled {r_all:+.3f}")
        axa.set_xlabel(f"rho({GOI}, cDFG)  with 95% bootstrap CI")
        axa.set_title("P5  Module strength by region\n"
                      "d = donors contributing, L = libraries")
        axa.legend(loc="lower right", fontsize=8)

        # Leave-one-region-out per arm
        bars, labs, bcols = [], [], []
        for arm, regs, c in (("NeoC", NEOC_REGIONS, CLR_NEOC),
                             ("NonNeoC", NONNEOC_REGIONS, CLR_NONNEOC)):
            present = [x for x in regs if x in set(reg)]
            if not present:
                continue
            inarm = np.isin(reg.values, present)
            r_full, _, _ = group_rho(inarm)
            for R in present:
                r_wo, _, _ = group_rho(inarm & (reg.values != R))
                share = 100.0 * (ANALYSED & inarm & (reg.values == R)).sum() \
                    / max((ANALYSED & inarm).sum(), 1)
                if np.isfinite(r_wo) and np.isfinite(r_full):
                    bars.append(r_wo - r_full)
                    labs.append(f"{arm}: −{R}  ({share:.0f}% of arm)")
                    bcols.append(c)
        if bars:
            yp = np.arange(len(bars))
            axb.barh(yp, bars, color=bcols, edgecolor="white", height=0.72)
            axb.axvline(0, color="black", lw=1.0)
            axb.set_yticks(yp)
            axb.set_yticklabels(labs, fontsize=8)
            axb.set_xlabel("change in pooled arm rho when the region is dropped")
            axb.set_title("P5b  Leave-one-region-out\n"
                          "a large bar on a large share = the arm IS that region")
        else:
            axb.set_visible(False)
        fig.tight_layout()
        finish_fig(fig, "S1e_P5_region_forest_and_LORO")
    else:
        say("  [figure] P5 skipped - too few regions with a usable rho")
except Exception as e:
    say(f"  [figure] P5 FAILED: {type(e).__name__}: {e}")


# == P6: the detectability trend and the ACTB residual ========================
try:
    gs = gene_stats.dropna(subset=["rho"]).copy()
    mito = gs[gs.is_mito & (gs.mean_counts > 0)]
    ref = gs[gs.gene == REF_GENE]
    if len(mito) >= 5 and len(ref) == 1:
        x = np.log10(mito.mean_counts.values)
        y = mito.rho.values
        b, a = np.polyfit(x, y, 1)
        pred = a + b * x
        r2 = 1 - ((y - pred) ** 2).sum() / ((y - y.mean()) ** 2).sum()
        xr = float(np.log10(ref.mean_counts.iloc[0]))
        yr = float(ref.rho.iloc[0])
        pred_ref = a + b * xr
        resid = yr - pred_ref
        worst = float(np.max(np.abs(y - pred)))

        fig, ax = plt.subplots(figsize=(7.6, 5.4))
        xs = np.linspace(min(x.min(), xr) - 0.15, x.max() + 0.15, 60)
        ax.plot(xs, a + b * xs, color=CLR_LINE, lw=1.6, zorder=2,
                label=f"mito trend: rho = {a:.3f} + {b:.3f}·log10(counts)"
                      f"   R² = {r2:.4f}")
        inmod = mito.in_module1.values
        ax.scatter(x[inmod], y[inmod], s=85, c=CLR_TARGET, zorder=4,
                   edgecolor="white", linewidth=0.8,
                   label="Module 1 mito genes")
        ax.scatter(x[~inmod], y[~inmod], s=85, c=CLR_OTHER, zorder=4,
                   edgecolor="white", linewidth=0.8,
                   label="non-module mito genes")
        for xi, yi, g in zip(x, y, mito.gene):
            ax.annotate(g.replace("MT-", ""), (xi, yi), fontsize=7,
                        xytext=(5, -9), textcoords="offset points")
        # ACTB, off the line
        ax.scatter([xr], [yr], s=190, marker="D", c=CLR_ACCENT, zorder=5,
                   edgecolor="black", linewidth=0.9, label=REF_GENE)
        ax.plot([xr, xr], [yr, pred_ref], color=CLR_ACCENT, ls="--", lw=1.6,
                zorder=3)
        ax.scatter([xr], [pred_ref], s=60, facecolor="none",
                   edgecolor=CLR_ACCENT, linewidth=1.4, zorder=5)
        ax.annotate(f"{REF_GENE} residual {resid:+.3f}\n"
                    f"({abs(resid)/worst:.1f}× the worst mito residual)",
                    xy=(xr, (yr + pred_ref) / 2),
                    xytext=(18, 0), textcoords="offset points",
                    fontsize=8.5, color=CLR_ACCENT, va="center",
                    arrowprops=dict(arrowstyle="-", color=CLR_ACCENT, lw=0.9))
        ax.set_xlabel(f"log10( mean counts per nucleus )   "
                      f"within {TARGET_STYPE}, {GOI} > 0")
        ax.set_ylabel(f"rho( {GOI} , gene )")
        ax.set_title("P6  Concordance is explained by abundance for mito genes,\n"
                     f"but {REF_GENE} sits far below the same trend - "
                     "so the shared factor is mito-SPECIFIC")
        ax.legend(loc="lower right", fontsize=8)
        fig.tight_layout()
        finish_fig(fig, "S1e_P6_detectability_trend_and_ACTB",
                   "the specificity argument, likely a paper figure")
    else:
        say("  [figure] P6 skipped - not enough genes")
except Exception as e:
    say(f"  [figure] P6 FAILED: {type(e).__name__}: {e}")


# == P7: post-mortem and pathology covariates at donor level ==================
try:
    if dsum is not None and "mito" in dsum.columns:
        panels = [(k, lab) for k, lab in
                  (("pmi", "PMI (hours)"), ("rin", "RIN"),
                   ("ph", "Brain pH"), ("adnc", "ADNC rank"))
                  if k in dsum.columns and dsum[k].notna().sum() > 10]
        if panels:
            ncol = min(4, len(panels))
            fig, axes = plt.subplots(2, ncol, figsize=(3.1 * ncol, 6.4),
                                     squeeze=False)
            for j, (k, lab) in enumerate(panels):
                for i, yv, ylab in ((0, "mito", "donor median mito%"),
                                    (1, "module_rho", "donor module rho")):
                    ax = axes[i][j]
                    if yv not in dsum.columns:
                        ax.set_visible(False)
                        continue
                    jj = dsum[[k, yv]].dropna()
                    if len(jj) < 8:
                        ax.set_visible(False)
                        continue
                    _r, _p = safe_rho(jj[k], jj[yv])
                    sig = (_p < 0.05) if np.isfinite(_p) else False
                    ax.scatter(jj[k], jj[yv], s=26,
                               c=(CLR_TARGET if sig else CLR_OTHER),
                               alpha=0.75, edgecolor="white", linewidth=0.5)
                    if len(jj) >= 4:
                        bb, aa = np.polyfit(jj[k], jj[yv], 1)
                        xs = np.linspace(jj[k].min(), jj[k].max(), 30)
                        ax.plot(xs, aa + bb * xs, color=CLR_LINE, lw=1.2)
                    ax.set_title(f"{lab}\nrho {_r:+.3f}  p {_p:.3g}"
                                 f"{'  *' if sig else ''}", fontsize=8.5)
                    if j == 0:
                        ax.set_ylabel(ylab, fontsize=8.5)
                    if i == 1:
                        ax.set_xlabel(lab, fontsize=8.5)
            fig.suptitle("P7  Donor-level alternatives: does degradation "
                         "(PMI, RIN, pH) explain the mito signal,\n"
                         "or does pathology (ADNC)?   "
                         "navy-teal = p < 0.05,  n = donors", fontsize=9.5)
            fig.tight_layout(rect=[0, 0, 1, 0.94])
            finish_fig(fig, "S1e_P7_donor_covariates")
        else:
            say("  [figure] P7 skipped - no donor covariates available")
    else:
        say("  [figure] P7 skipped - donor summary not built")
except Exception as e:
    say(f"  [figure] P7 FAILED: {type(e).__name__}: {e}")


# == P8: is the analysed subset biased? =======================================
try:
    z = is_target & (~goi_pos)
    if z.sum() > 50:
        fig, axes = plt.subplots(1, 3, figsize=(12.4, 3.9))
        for ax, v, lab, logx in (
                (axes[0], total_counts, "total UMI depth", True),
                (axes[1], mito_pct, "mitochondrial %", True),
                (axes[2], CDFG, "cDFG composite", False)):
            aa = v[ANALYSED]
            bb = v[z]
            aa = aa[np.isfinite(aa) & (aa > (0 if logx else -np.inf))]
            bb = bb[np.isfinite(bb) & (bb > (0 if logx else -np.inf))]
            if aa.size < 20 or bb.size < 20:
                ax.set_visible(False)
                continue
            lo = min(aa.min(), bb.min())
            hi = max(aa.max(), bb.max())
            bins = (np.logspace(np.log10(max(lo, 1e-3)), np.log10(hi), 45)
                    if logx else np.linspace(lo, hi, 45))
            ax.hist(aa, bins=bins, density=True, color=CLR_TARGET, alpha=0.6,
                    label=f"analysed  n={aa.size:,}")
            ax.hist(bb, bins=bins, density=True, color=CLR_OTHER, alpha=0.6,
                    label=f"{GOI}==0, excluded  n={bb.size:,}")
            if logx:
                ax.set_xscale("log")
            ax.set_xlabel(lab)
            ax.set_ylabel("density")
            ax.legend(fontsize=7.5)
        fig.suptitle("P8  The analysed set excludes "
                     f"{GOI}-zero nuclei - how different are they?   "
                     "a depth shift is a selection the methods must state",
                     fontsize=9.5)
        fig.tight_layout(rect=[0, 0, 1, 0.92])
        finish_fig(fig, "S1e_P8_analysed_vs_excluded")
    else:
        say("  [figure] P8 skipped - too few excluded nuclei")
except Exception as e:
    say(f"  [figure] P8 FAILED: {type(e).__name__}: {e}")


# == P9: the segmented threshold screen across supertypes =====================
try:
    if seg_curves and len(seg_curves) >= 2:
        fig, (axa, axb) = plt.subplots(1, 2, figsize=(12.4, 4.8))
        order = [TARGET_STYPE] + [k for k in seg_curves if k != TARGET_STYPE]
        for S in order:
            if S not in seg_curves:
                continue
            xb, yb, brk = seg_curves[S]
            is_t = (S == TARGET_STYPE)
            axa.plot(xb, yb, "o-", ms=3.5,
                     color=(CLR_TARGET if is_t else CLR_OTHER),
                     lw=(2.2 if is_t else 1.1),
                     alpha=(1.0 if is_t else 0.6),
                     zorder=(5 if is_t else 2),
                     label=(S if is_t else None))
            axa.axvline(brk, color=(CLR_TARGET if is_t else CLR_OTHER),
                        ls=":", lw=(1.6 if is_t else 0.8),
                        alpha=(0.9 if is_t else 0.4),
                        zorder=(4 if is_t else 1))
        axa.set_xlabel(f"{GOI} (binned, ln(UP10K+1))")
        axa.set_ylabel("mean cDFG in bin")
        axa.set_title("P9  Binned module curve per supertype\n"
                      "dotted lines = fitted breakpoints; target in colour")
        axa.legend(loc="upper left", fontsize=8)

        if segdf is not None and len(segdf):
            dd = segdf.sort_values("rss_drop")
            cc = [CLR_TARGET if s_ == TARGET_STYPE else CLR_OTHER
                  for s_ in dd.supertype]
            yy = np.arange(len(dd))
            axb.barh(yy, 100 * dd.rss_drop, color=cc, edgecolor="white",
                     height=0.7)
            axb.set_yticks(yy)
            axb.set_yticklabels([f"{s_}  (break {b:.2f})"
                                 for s_, b in zip(dd.supertype,
                                                  dd.breakpoint)],
                                fontsize=8)
            axb.set_xlabel("% of RSS removed by allowing one breakpoint")
            axb.set_title("P9b  How much better is two lines than one?\n"
                          "specificity would mean the target stands out")
        fig.tight_layout()
        finish_fig(fig, "S1e_P9_segmented_threshold_screen",
                   "screen only - NOT the CEP-IP inflection point")
    else:
        say("  [figure] P9 skipped - segmented screen produced no curves")
except Exception as e:
    say(f"  [figure] P9 FAILED: {type(e).__name__}: {e}")


# == P10: depth-matched supertype comparison ==================================
try:
    if matched is not None and len(matched) >= 3:
        dd = matched.sort_values("rho_depth_matched")
        yy = np.arange(len(dd))
        fig, ax = plt.subplots(figsize=(8.2, 4.4))
        for y, (_, a) in zip(yy, dd.iterrows()):
            c = CLR_TARGET if a.supertype == TARGET_STYPE else CLR_OTHER
            ax.plot([a.rho_raw, a.rho_depth_matched], [y, y], color=c,
                    lw=1.2, alpha=0.6, zorder=2)
            ax.plot(a.rho_raw, y, "o", color=c, ms=7, mfc="white",
                    mew=1.6, zorder=3)
            ax.plot(a.rho_depth_matched, y, "o", color=c, ms=9, zorder=4)
        ax.set_yticks(yy)
        ax.set_yticklabels([f"{s_}  (n={int(n):,})"
                            for s_, n in zip(dd.supertype, dd.n)], fontsize=8)
        ax.set_xlabel(f"rho({GOI}, cDFG)")
        ax.set_title("P10  Raw (open) vs depth-standardised (filled)\n"
                     "standardised onto the target's own depth distribution")
        fig.tight_layout()
        finish_fig(fig, "S1e_P10_supertype_depth_matched")
    else:
        say("  [figure] P10 skipped - depth matching produced no table")
except Exception as e:
    say(f"  [figure] P10 FAILED: {type(e).__name__}: {e}")


# == P11: the abundance-matched panel null ====================================
try:
    if null_rhos is not None and null_rhos.size >= 30:
        fig, ax = plt.subplots(figsize=(7.4, 4.4))
        ax.hist(null_rhos, bins=30, color=CLR_OTHER, alpha=0.85,
                edgecolor="white",
                label=f"non-mito panels matched on abundance\nand detection "
                      f"({null_rhos.size} draws)")
        ax.axvline(r_all, color=CLR_TARGET, lw=2.4,
                   label=f"observed Module 1  {r_all:+.3f}")
        ax.axvline(float(np.median(null_rhos)), color=CLR_LINE, ls="--",
                   lw=1.3, label=f"null median  {np.median(null_rhos):+.3f}")
        _pv = float((null_rhos >= r_all).mean())
        ax.set_xlabel(f"rho( GOI , 8-gene composite )")
        ax.set_ylabel("panels")
        ax.set_title("P11  Does an abundance-matched NON-mito panel reach it?\n"
                     f"empirical p = {_pv:.4f}   "
                     "(this is the control detectability cannot explain away)")
        ax.legend(loc="upper left", fontsize=8)
        fig.tight_layout()
        finish_fig(fig, "S1e_P11_matched_panel_null",
                   "the decisive specificity control")
    else:
        say("  [figure] P11 skipped - no null distribution")
except Exception as e:
    say(f"  [figure] P11 FAILED: {type(e).__name__}: {e}")


# == P12: the exogenous soup index ============================================
try:
    if EXO is not None and np.isfinite(EXO[ANALYSED]).sum() > 500:
        qe = quartile_rho(EXO, "soup")
        fig, (axa, axb) = plt.subplots(1, 2, figsize=(11.6, 4.4))
        v = EXO[ANALYSED]
        v = v[np.isfinite(v)]
        axa.hist(v[v > 0], bins=60, color=CLR_TARGET, alpha=0.85,
                 edgecolor="white")
        axa.set_yscale("log")
        axa.set_xlabel("foreign-lineage transcripts, % of the nucleus's UMIs")
        axa.set_ylabel("nuclei (log scale)")
        axa.set_title("P12  Per-nucleus soup index\n"
                      f"{100*(EXO[ANALYSED] == 0).mean():.1f}% of nuclei carry "
                      "none at all")
        if qe is not None:
            x = np.arange(len(qe))
            yerr = np.vstack([qe.rho - qe.lo, qe.hi - qe.rho])
            axb.errorbar(x, qe.rho, yerr=yerr, fmt="o-", color=CLR_TARGET,
                         capsize=4, lw=1.8, ms=8)
            axb.set_xticks(x)
            axb.set_xticklabels([f"{r['q']}\nn={r['n']:,}\n{r['median']:.3f}%"
                                 for _, r in qe.iterrows()], fontsize=7.5)
            axb.set_xlabel("quartile of the per-nucleus soup index")
            axb.set_ylabel(f"rho({GOI}, cDFG)  with 95% CI")
            axb.set_title("P12b  Ambient predicts this to RISE left to right\n"
                          "tested WITHIN library, unlike P1")
        fig.tight_layout()
        finish_fig(fig, "S1e_P12_exogenous_soup_index",
                   "the within-library ambient test")
    else:
        say("  [figure] P12 skipped - no soup index")
except Exception as e:
    say(f"  [figure] P12 FAILED: {type(e).__name__}: {e}")


# == P13: where do the data want the break? ===================================
try:
    if ip_curves and len(ip_curves) >= 2:
        fig, (axa, axb) = plt.subplots(1, 2, figsize=(12.6, 4.8))
        for S, pf in ip_curves.items():
            is_t = (S == TARGET_STYPE)
            # Each RSS profile scaled to 0-1 so supertypes are comparable
            rr = pf.rss.values
            rr = (rr - rr.min()) / max(rr.max() - rr.min(), 1e-12)
            axa.plot(pf.q, rr, "-o", ms=(4 if is_t else 2.5),
                     lw=(2.4 if is_t else 1.0),
                     color=(CLR_TARGET if is_t else CLR_OTHER),
                     alpha=(1.0 if is_t else 0.5),
                     zorder=(5 if is_t else 2),
                     label=(S if is_t else None))
        # Marks the comparator (q = 0.50) and the published IP quantiles.
        axa.axvline(0.50, color=CLR_ACCENT, ls="--", lw=1.8, zorder=3,
                    label="comparator Si=0 (q=0.50)")
        axa.text(0.505, 0.99, " comparator Si=0\n (q=0.50, median)",
                 color=CLR_ACCENT, fontsize=8, va="top")
        if not USE_SI0_IP:
            for _i, (_arm, _u) in enumerate(sorted(IP_RX0_QUANTILE.items())):
                axa.axvline(_u, color=CLR_NEOC, lw=1.5, zorder=3,
                            ls=("-." if _i == 0 else ":"),
                            label=f"published IP, {_arm} (u={_u:.3f})")
            axa.text(0.02, 0.02,
                     "published IP = R(x)=0 crossing (USE_SI0_IP=False);\n"
                     "the median line is the comparator, not the IP",
                     transform=axa.transAxes, fontsize=7.5,
                     color="#444444", va="bottom")
        axa.set_xlabel("candidate break point, as a quantile of GOI")
        axa.set_ylabel("fit error, scaled 0 (best) to 1 (worst)")
        axa.set_title("P13  Where do the data want the break?\n"
                      "a broad flat basin means the location is only "
                      "weakly determined")
        axa.legend(loc="upper right", fontsize=7.5)

        if ipdf is not None and len(ipdf):
            dd = ipdf.sort_values("rss_penalty", ascending=False)
            yy = np.arange(len(dd))
            cc = [CLR_TARGET if s_ == TARGET_STYPE else CLR_OTHER
                  for s_ in dd.supertype]
            # Published-IP penalty beside the comparator's, on the full-grid scale.
            _pubcol = next((c for c in dd.columns
                            if c.startswith("penalty_at_u_")), None)
            if _pubcol:
                axb.barh(yy + 0.19, dd.rss_penalty, color=cc,
                         edgecolor="white", height=0.36,
                         label="comparator, q=0.50")
                axb.barh(yy - 0.19, dd[_pubcol], color=cc, alpha=0.45,
                         edgecolor="white", height=0.36, hatch="//",
                         label=f"published IP, u="
                               f"{_pubcol.replace('penalty_at_u_', '')}")
                for _y, _v in zip(yy + 0.19, dd.rss_penalty):
                    axb.text(_v + 0.012, _y, f"{_v:.3f}", va="center",
                             fontsize=7, color="#444444")
                for _y, _v in zip(yy - 0.19, dd[_pubcol]):
                    axb.text(_v + 0.012, _y, f"{_v:.3f}", va="center",
                             fontsize=7, color="#444444")
                axb.legend(loc="lower right", fontsize=7.5)
            else:
                axb.barh(yy, dd.rss_penalty, color=cc, edgecolor="white",
                         height=0.7)
                for _y, _v in zip(yy, dd.rss_penalty):
                    axb.text(_v + 0.012, _y, f"{_v:.3f}", va="center",
                             fontsize=7, color="#444444")
            axb.set_yticks(yy)
            axb.set_yticklabels(
                [f"{s_}  (prefers q={q:.2f})"
                 for s_, q in zip(dd.supertype, dd.q_preferred)], fontsize=8)
            # The 0-1 axis is kept so the penalties are not exaggerated.
            axb.set_xlim(0, 1)
            axb.set_xlabel("how much worse that break is than the best one\n"
                           "0 = as good as the best, 1 = the worst on the grid")
            axb.set_title("P13b  Penalty at the comparator vs at the "
                          "published IP\n"
                          "all bars near zero = no break point is "
                          "distinguished")
        fig.tight_layout()
        finish_fig(fig, "S1e_P13_is_the_median_special",
                   "the reviewer objection, tested directly")
    else:
        say("  [figure] P13 skipped - no quantile profiles")
except Exception as e:
    say(f"  [figure] P13 FAILED: {type(e).__name__}: {e}")


say("")
say(f"  {len(_figs_made)} figure(s) written to {FIG_DIR}")


# =============================================================================
hdr("F14 - SUMMARY")

say(f"  pooled module rho                    : {fmt(r_all)}")
if len(usable) >= 8:
    say(f"  libraries analysed                   : {len(usable)}")
    say(f"  per-library rho, median (SD)         : "
        f"{usable.rho.median():+.4f} ({usable.rho.std():.4f})")
say(f"  library nested within donor          : "
    f"{'YES - batch and donor inseparable' if int((_ld > 1).sum()) == 0 else 'no'}")
if AMB is not None and len(usable) >= 8 and usable['amb_index'].notna().sum() >= 8:
    _dd = usable.dropna(subset=["amb_index", "rho"])
    _r, _p = safe_rho(_dd["amb_index"], _dd["rho"], min_n=8)
    _rc, _pc = partial_rho(_dd["amb_index"], _dd["rho"],
                           _dd["median_mito"], _dd["median_depth"], min_n=8)
    _rm, _pm = partial_rho(_dd["median_mito"], _dd["rho"],
                           _dd["amb_index"], min_n=8)
    say(f"  ambient index vs module, RAW         : {fmt(_r)} (p={_p:.3g})")
    say(f"  ...conditioned on mito% + depth      : {fmt(_rc)} (p={_pc:.3g})")
    say(f"  mito% vs module | ambient index      : {fmt(_rm)} (p={_pm:.3g})")
if matched is not None and len(matched):
    _tm = matched[matched.supertype == TARGET_STYPE]
    _om = matched[matched.supertype != TARGET_STYPE]
    if len(_tm) and len(_om):
        say(f"  target minus best other, matched depth: "
            f"{float(_tm.rho_depth_matched.iloc[0] - _om.rho_depth_matched.max()):+.4f}")
if segdf is not None and len(segdf):
    _ts = segdf[segdf.supertype == TARGET_STYPE]
    _os = segdf[segdf.supertype != TARGET_STYPE]
    if len(_ts) and len(_os):
        say(f"  threshold screen, target vs best other: "
            f"{100*float(_ts.rss_drop.iloc[0]):.1f}% vs "
            f"{100*float(_os.rss_drop.max()):.1f}% RSS removed")
say(f"  IP definition assumed        : "
    f"{'Si=0 (median)' if USE_SI0_IP else 'R(x)=0 crossing'}"
    f"   [USE_SI0_IP={USE_SI0_IP}]")
if len(_tgt) and len(_oth):
    say(f"  target rho minus best other (n>={MIN_N_COMPARE}): "
        f"{float(_tgt.rho.iloc[0]) - float(_oth.rho.max()):+.4f}")

with open(os.path.join(OUT_DIR, "S1e_techqc_report.txt"), "w",
          encoding="utf-8") as f:
    f.write("\n".join(_report))
say("")
say(f"  saved: {OUT_DIR}/S1e_techqc_report.txt")
say("")
say("=" * 74)
say("  PART 1e COMPLETE")
say("=" * 74)
