# =============================================================================
# Part1f_Ordinals_code.py
# SEA-AD CEP-IP/CMRF, Part 1f, clinical ordinals and donor-level shape test
#   Audits how each clinical ordinal (ADNC, Braak, Thal, CERAD, CAA) parses,
#   builds a donor table with tissue-quality covariates and their complete
#   cases, and compares a linear against a single-step fit of donor mito%
#   (and module strength, if Part 1e ran) on each pathology axis.
#
# Usage (from the SEA-AD folder, or set SEAAD_ROOT), after Part 1e
#   python3 Part1f_Ordinals_code.py
#
# Outputs (Part1_PreChecks/): ordinal levels, donor table and shape tests
#   (CSV), and the report.
# =============================================================================

import os
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import spearmanr, rankdata

np.random.seed(42)

# == Paths ===================================================================
os.chdir(os.environ.get("SEAAD_ROOT", os.getcwd()))

FULL_H5AD = ("SEAAD_Immune_HIP_MEC_LEC_ITG_MTG_FI_STG_DFC_AnG_V1C_"
             "RNAseq_final-nuclei.2026-06-22.h5ad")
OUT_DIR = "Part1_PreChecks"
os.makedirs(OUT_DIR, exist_ok=True)
PER_LIBRARY_CSV = os.path.join(OUT_DIR, "S1e_per_library.csv")

# == Ordinal maps (as in Part 2 and R_config.R) ==============================
BRAAK_MAP = {"Braak 0": 0, "Braak I": 1, "Braak II": 2, "Braak III": 3,
             "Braak IV": 4, "Braak V": 5, "Braak VI": 6}
ADNC_MAP = {"Not AD": 0, "Low": 1, "Intermediate": 2, "High": 3}
# Thal and CERAD: candidate maps; G1 lists the actual level strings.
THAL_MAP = {"Thal 0": 0, "Thal 1": 1, "Thal 2": 2, "Thal 3": 3,
            "Thal 4": 4, "Thal 5": 5,
            "Thal phase 0": 0, "Thal phase 1": 1, "Thal phase 2": 2,
            "Thal phase 3": 3, "Thal phase 4": 4, "Thal phase 5": 5}
CERAD_MAP = {"Absent": 0, "None": 0, "Sparse": 1, "Moderate": 2,
             "Frequent": 3}

ORDINALS = [
    ("ADNC", ["Overall AD neuropathological Change", "ADNC"], ADNC_MAP),
    ("Braak", ["Braak", "Braak stage"], BRAAK_MAP),
    ("Thal", ["Thal", "Thal phase"], THAL_MAP),
    ("CERAD", ["CERAD score", "CERAD"], CERAD_MAP),
    ("CAA", ["Overall CAA Score", "CAA"], None),
]

# Donor-level covariates: Tier 1 (four-covariate model) and Tier 2 additions.
TIER1_COVARS = [("PMI", ["PMI"]), ("Age", ["Age at Death"])]
TIER2_COVARS = [("RIN", ["RIN"]),
                ("pH", ["Brain pH"]),
                ("Weight", ["Fresh Brain Weight"])]

MAX_UNMAPPED_PCT = 5.0     # flag threshold for a parser

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


def find_col(cols, candidates):
    cols = list(cols)
    for c in candidates:
        if c in cols:
            return c
    low = {str(c).lower(): c for c in cols}
    for c in candidates:
        if str(c).lower() in low:
            return low[str(c).lower()]
    for c in candidates:
        for k, v in low.items():
            if str(c).lower() in k:
                return v
    return None


def safe_rho(x, y, min_n=8):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < min_n:
        return np.nan, np.nan, int(ok.sum())
    r, p = spearmanr(x[ok], y[ok])
    return float(r), float(p), int(ok.sum())


def partial_rho(x, y, Z, min_n=8):
    """Spearman partial correlation of x,y given the columns of Z."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    Z = np.asarray(Z, dtype=float)
    if Z.ndim == 1:
        Z = Z[:, None]
    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(Z).all(axis=1)
    if ok.sum() < min_n + Z.shape[1]:
        return np.nan, np.nan, int(ok.sum())
    R = np.column_stack([rankdata(x[ok]), rankdata(y[ok])]
                        + [rankdata(Z[ok, j]) for j in range(Z.shape[1])])
    sd = R.std(0)
    sd[sd == 0] = 1.0
    R = (R - R.mean(0)) / sd
    D = np.column_stack([np.ones(len(R)), R[:, 2:]])
    ex = R[:, 0] - D @ np.linalg.lstsq(D, R[:, 0], rcond=None)[0]
    ey = R[:, 1] - D @ np.linalg.lstsq(D, R[:, 1], rcond=None)[0]
    r, p = spearmanr(ex, ey)
    return float(r), float(p), int(ok.sum())


# =============================================================================
hdr("G0 - load obs (backed: no expression matrix is read)")

ad = sc.read_h5ad(FULL_H5AD, backed="r")
obs = ad.obs.copy()
say(f"  file  : {FULL_H5AD}")
say(f"  shape : {ad.shape}")
say(f"  obs   : {obs.shape[0]:,} rows x {obs.shape[1]} columns")

_oc = list(obs.columns)
donor_col = find_col(_oc, ["Donor ID", "donor_id", "Donor"])
region_col = find_col(_oc, ["Brain Region", "Region"])
mito_col = find_col(_oc, ["Fraction mitochondrial UMIs"])
numi_col = find_col(_oc, ["Number of UMIs"])
ngene_col = find_col(_oc, ["Genes detected"])
severe_col = find_col(_oc, ["Severely Affected Donor"])

if donor_col is None:
    raise SystemExit("FATAL: no donor column resolved. Cannot proceed.")

say("")
say("  resolved columns:")
for nm, cc in [("donor", donor_col), ("region", region_col),
               ("mito frac", mito_col), ("nUMI", numi_col),
               ("nGene", ngene_col), ("severe", severe_col)]:
    say(f"    {nm:<12} -> {cc}")

n_donors = obs[donor_col].astype(str).nunique()
say("")
say(f"  distinct donors : {n_donors}")


# =============================================================================
hdr("G1 - ordinal parsing audit")

say("  Per ordinal: levels, a digit-only parse, and the declared map.")
say("")

level_rows = []
parsed = {}
ordinal_status = {}

for name, cands, cmap in ORDINALS:
    col = find_col(_oc, cands)
    say("-" * 74)
    if col is None:
        say(f"  {name}: NO COLUMN FOUND (tried {cands})")
        ordinal_status[name] = "missing"
        continue

    # Donor attributes, so one value per donor.
    per_donor = obs.groupby(donor_col, observed=True)[col].first()
    is_cat = isinstance(obs[col].dtype, pd.CategoricalDtype)
    txt = per_donor.astype(str).str.strip()
    levels = txt.value_counts().sort_index()

    say(f"  {name}  (column '{col}')")
    say(f"    dtype: {obs[col].dtype}   Categorical: {is_cat}")
    say(f"    distinct levels across {len(per_donor)} donors: {len(levels)}")

    # A digit-only parser, for comparison
    naive = pd.to_numeric(txt.str.extract(r"(\d+)")[0], errors="coerce")
    naive_fail = int(naive.isna().sum())
    naive_pct = 100.0 * naive_fail / max(len(txt), 1)

    # The declared map
    if cmap:
        mapped = pd.to_numeric(txt.map(cmap), errors="coerce")
    else:
        mapped = pd.Series(np.nan, index=txt.index)
    map_fail = int(mapped.isna().sum())
    map_pct = 100.0 * map_fail / max(len(txt), 1)

    say("")
    say(f"    {'level':<28}{'donors':>8}{'digit parse':>13}{'declared map':>14}")
    for lv, n in levels.items():
        _nv = naive[txt == lv].iloc[0] if (txt == lv).any() else np.nan
        _mv = mapped[txt == lv].iloc[0] if (txt == lv).any() else np.nan
        say(f"    {lv[:27]:<28}{n:>8}"
            f"{('NaN' if not np.isfinite(_nv) else f'{_nv:.0f}'):>13}"
            f"{('NaN' if not np.isfinite(_mv) else f'{_mv:.0f}'):>14}")
        level_rows.append(dict(ordinal=name, column=col, level=lv, n_donors=int(n),
                               naive_digit_parse=(np.nan if not np.isfinite(_nv)
                                                  else float(_nv)),
                               declared_map=(np.nan if not np.isfinite(_mv)
                                             else float(_mv)),
                               is_categorical=is_cat))

    say("")
    say(f"    naive digit parser fails on : {naive_fail}/{len(txt)} donors "
        f"({naive_pct:.1f}%)")
    if cmap:
        say(f"    declared map fails on       : {map_fail}/{len(txt)} donors "
            f"({map_pct:.1f}%)")
    else:
        say("    declared map                : none supplied for this field")

    # .map() on a Categorical returns a Categorical
    if is_cat:
        raw_mapped = txt.astype("category").map(cmap) if cmap else None
        if raw_mapped is not None:
            still_cat = isinstance(raw_mapped.dtype, pd.CategoricalDtype)
            say(f"    .map() without pd.to_numeric returns Categorical: "
                f"{still_cat}")

    # Verdict for this field
    if cmap is None:
        ordinal_status[name] = "no map"
        say("")
        say(f"    >>> no map declared for {name}")
    elif map_pct == 0:
        ordinal_status[name] = "ok"
        parsed[name] = mapped
        say("")
        say(f"    >>> clean: all {len(txt)} donors map")
    elif map_pct <= MAX_UNMAPPED_PCT:
        ordinal_status[name] = "partial"
        parsed[name] = mapped
        _un = sorted(set(txt[mapped.isna()]))
        say("")
        say(f"    >>> partial: {map_fail} donor(s) unmapped: {_un[:8]}")
    else:
        ordinal_status[name] = "FAIL"
        parsed[name] = mapped
        _un = sorted(set(txt[mapped.isna()]))
        say("")
        say(f"    >>> fail: {map_pct:.1f}% of donors unmapped {_un[:8]}")

say("-" * 74)
say("")
say("  SUMMARY")
for name, _, _ in ORDINALS:
    say(f"    {name:<8} {ordinal_status.get(name, 'not reached')}")
_bad = [k for k, v in ordinal_status.items() if v in ("FAIL", "missing")]
if _bad:
    say("")
    say(f"  {len(_bad)} ordinal(s) not usable: {_bad}")

if level_rows:
    pd.DataFrame(level_rows).to_csv(
        os.path.join(OUT_DIR, "S1f_ordinal_levels.csv"), index=False)
    say(f"  saved: {OUT_DIR}/S1f_ordinal_levels.csv")


# =============================================================================
hdr("G2 - donor table and covariate completeness")

don = obs[donor_col].astype(str)
dtab = pd.DataFrame(index=sorted(don.unique()))
dtab.index.name = "donor"

for name in parsed:
    dtab[name] = parsed[name].reindex(dtab.index).values

for nm, cands in TIER1_COVARS + TIER2_COVARS:
    c = find_col(_oc, cands)
    if c is None:
        say(f"  [WARN] covariate '{nm}' not found (tried {cands})")
        dtab[nm] = np.nan
        continue
    dtab[nm] = pd.to_numeric(
        obs.groupby(don, observed=True)[c].first().reindex(dtab.index),
        errors="coerce").values

# Donor summaries from the nuclei
if mito_col is not None:
    dtab["mito_pct"] = pd.to_numeric(
        obs.groupby(don, observed=True)[mito_col].median()
        .reindex(dtab.index), errors="coerce").values * 100.0
if numi_col is not None:
    dtab["depth"] = pd.to_numeric(
        obs.groupby(don, observed=True)[numi_col].median()
        .reindex(dtab.index), errors="coerce").values
dtab["n_nuclei"] = obs.groupby(don, observed=True).size().reindex(dtab.index).values
if severe_col is not None:
    _sv = obs.groupby(don, observed=True)[severe_col].first().astype(str)
    dtab["severe"] = _sv.reindex(dtab.index).str.upper().str[0].eq("Y").values

# Module strength per donor, from the Part 1e per-library CSV
MODULE_AVAILABLE = False
if os.path.exists(PER_LIBRARY_CSV):
    lib = pd.read_csv(PER_LIBRARY_CSV)
    if {"donor", "rho", "n_analysed"} <= set(lib.columns):
        _u = lib.dropna(subset=["rho"])
        # Each library weighted by its analysed cells
        _w = (_u.assign(wr=_u.rho * _u.n_analysed)
                .groupby("donor")
                .agg(wr=("wr", "sum"), n=("n_analysed", "sum"),
                     n_libs=("rho", "size")))
        dtab["module_rho"] = (_w.wr / _w.n).reindex(dtab.index).values
        dtab["module_n_libs"] = _w.n_libs.reindex(dtab.index).values
        MODULE_AVAILABLE = True
        say(f"  module strength per donor taken from {PER_LIBRARY_CSV}")
        say(f"    donors with a module rho: "
            f"{int(dtab['module_rho'].notna().sum())} of {len(dtab)}")
else:
    say(f"  [note] {PER_LIBRARY_CSV} not found - module strength unavailable.")

dtab.reset_index().to_csv(os.path.join(OUT_DIR, "S1f_donor_table.csv"),
                          index=False)
say("")
say(f"  saved: {OUT_DIR}/S1f_donor_table.csv   ({dtab.shape[0]} donors x "
    f"{dtab.shape[1]} columns)")

say("")
say("  per-variable completeness:")
say(f"    {'variable':<16}{'present':>9}{'missing':>9}")
for c in dtab.columns:
    if dtab[c].dtype == bool:
        continue
    say(f"    {c:<16}{int(dtab[c].notna().sum()):>9}"
        f"{int(dtab[c].isna().sum()):>9}")

_t1 = [c for c, _ in TIER1_COVARS if c in dtab.columns]
_t2 = [c for c, _ in TIER2_COVARS if c in dtab.columns]
say("")
say("  complete cases for the covariate sets:")
for lbl, cols in [("Tier 1", _t1),
                  ("Tier 1 + Tier 2", _t1 + _t2)]:
    if not cols:
        continue
    n_ok = int(dtab[cols].notna().all(axis=1).sum())
    say(f"    {lbl:<28} {n_ok:>3} of {len(dtab)} donors   ({cols})")

for name in parsed:
    cols = _t1 + _t2 + [name]
    cols = [c for c in cols if c in dtab.columns]
    n_ok = int(dtab[cols].notna().all(axis=1).sum())
    say(f"      ...with {name:<8} as the outcome: {n_ok:>3} donors")



# =============================================================================
hdr("G3 - shape test: linear against single step")

say("  linear: y ~ stage;  step: y ~ 1[stage > cut], best interior cut; compared by R2.")
say("")


def shape_test(y, s):
    """Compare a linear-in-stage fit against the best single-step fit."""
    y = np.asarray(y, dtype=float)
    s = np.asarray(s, dtype=float)
    ok = np.isfinite(y) & np.isfinite(s)
    y, s = y[ok], s[ok]
    if len(y) < 12 or len(np.unique(s)) < 3:
        return None
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    if ss_tot <= 0:
        return None
    # Linear in stage
    b, a = np.polyfit(s, y, 1)
    r2_lin = 1.0 - float(np.sum((y - (a + b * s)) ** 2)) / ss_tot
    # Best single step over the interior cuts
    cuts = np.unique(s)[:-1]
    best = None
    for c in cuts:
        g = s > c
        if g.sum() < 4 or (~g).sum() < 4:
            continue
        pred = np.where(g, y[g].mean(), y[~g].mean())
        r2 = 1.0 - float(np.sum((y - pred) ** 2)) / ss_tot
        if best is None or r2 > best[1]:
            best = (float(c), r2, float(y[g].mean() - y[~g].mean()),
                    int(g.sum()), int((~g).sum()))
    if best is None:
        return None
    cut, r2_step, delta, n_hi, n_lo = best
    rho, p, n = safe_rho(s, y)
    return dict(n=n, rho=rho, p=p, r2_linear=r2_lin, r2_step=r2_step,
                r2_gain=r2_step - r2_lin, best_cut=cut, step_delta=delta,
                n_below=n_lo, n_above=n_hi)


responses = [("mito_pct", "donor median mito%")]
if MODULE_AVAILABLE:
    responses.append(("module_rho", "donor module strength"))

adj_cols = [c for c in _t2 + _t1 if c in dtab.columns]
shape_rows = []

for ycol, ylab in responses:
    if ycol not in dtab.columns:
        continue
    say("-" * 74)
    say(f"  response: {ylab}  ({ycol})")
    say("")
    for name in parsed:
        if name not in dtab.columns:
            continue
        res = shape_test(dtab[ycol], dtab[name])
        if res is None:
            say(f"    {name:<8} not enough donors or stages")
            continue

        # Per-stage means
        _d = dtab[[ycol, name]].dropna()
        _g = _d.groupby(name)[ycol].agg(["size", "mean"])
        say(f"    {name}:  rho={res['rho']:+.4f}  p={res['p']:.3g}  "
            f"n={res['n']}")
        say(f"      {'stage':>7}{'donors':>8}{'mean':>10}")
        for st, row in _g.iterrows():
            say(f"      {st:>7.0f}{int(row['size']):>8}{row['mean']:>10.3f}")
        say(f"      linear R2 {res['r2_linear']:.3f}   "
            f"best-step R2 {res['r2_step']:.3f}   "
            f"gain {res['r2_gain']:+.3f}")
        say(f"      best cut: above stage {res['best_cut']:.0f}  "
            f"({res['n_below']} below vs {res['n_above']} above), "
            f"step size {res['step_delta']:+.3f}")

        # The same comparison on tissue-quality residuals
        rho_adj = p_adj = np.nan
        n_adj = 0
        if adj_cols:
            rho_adj, p_adj, n_adj = partial_rho(
                dtab[name], dtab[ycol], dtab[adj_cols].values)
            say(f"      adjusted for {adj_cols}: rho={rho_adj:+.4f} "
                f"p={p_adj:.3g} n={n_adj}")

        if res["r2_gain"] > 0.05:
            say("      >>> step-like")
        elif res["r2_gain"] > 0.0:
            say("      >>> marginally step-like (treated as linear)")
        else:
            say("      >>> linear")
        say("")

        shape_rows.append(dict(response=ycol, axis=name, **res,
                               rho_adjusted=rho_adj, p_adjusted=p_adj,
                               n_adjusted=n_adj,
                               adjusted_for=";".join(adj_cols)))

if shape_rows:
    pd.DataFrame(shape_rows).to_csv(
        os.path.join(OUT_DIR, "S1f_shape_tests.csv"), index=False)
    say("-" * 74)
    say(f"  saved: {OUT_DIR}/S1f_shape_tests.csv")


# =============================================================================
hdr("G4 - verdict")

_fail = [k for k, v in ordinal_status.items() if v in ("FAIL", "missing")]
_nomap = [k for k, v in ordinal_status.items() if v == "no map"]

say(f"  ordinals clean and usable          : "
    f"{[k for k, v in ordinal_status.items() if v == 'ok']}")
if _nomap:
    say(f"  ordinals with no declared map      : {_nomap}")
if _fail:
    say(f"  ordinals not usable as they stand  : {_fail}")

if _t1 + _t2:
    say(f"  complete cases, 7-covariate model   : "
        f"{int(dtab[[c for c in _t1 + _t2 if c in dtab.columns]].notna().all(axis=1).sum())}"
        f" of {len(dtab)} donors")

if shape_rows:
    _sr = pd.DataFrame(shape_rows)
    _step = _sr[_sr.r2_gain > 0.05]
    say(f"  shape tests run                     : {len(_sr)}")
    say(f"  showing a step-like (flat-then-rising) shape : {len(_step)}")
    if len(_step):
        for _, r in _step.iterrows():
            say(f"      {r['response']} vs {r['axis']}: gain "
                f"{r['r2_gain']:+.3f}, cut above stage {r['best_cut']:.0f}")
    else:
        say("      none - every relationship is as well described by a line.")


with open(os.path.join(OUT_DIR, "S1f_ordinals_report.txt"), "w",
          encoding="utf-8") as f:
    f.write("\n".join(_report))
say("")
say(f"  saved: {OUT_DIR}/S1f_ordinals_report.txt")
