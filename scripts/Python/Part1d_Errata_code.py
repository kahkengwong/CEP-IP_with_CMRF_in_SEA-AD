# =============================================================================
# Part1d_Errata_code.py
# SEA-AD CEP-IP/CMRF, Part 1d, audit of obs against the SEA-AD specimen errata
#   Checks whether the corrections to donor ID, sex and brain region in the
#   SEA-AD Sage change log (06/25/2026) are present in the 2026-06-22 h5ad,
#   checks donor IDs and sex against the 84-donor table, lists the regions and
#   reports donor x region coverage. Reads obs only, in backed mode.
#
# Usage (from the SEA-AD folder, or set SEAAD_ROOT), after Part 1c
#   python3 Part1d_Errata_code.py
#
# Outputs (Part1_PreChecks/): audit report, specimen audit and donor x region
#   coverage tables (CSV).
# =============================================================================

import os
import warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import scanpy as sc

np.random.seed(42)
sc.settings.verbosity = 1

os.chdir(os.environ.get("SEAAD_ROOT", os.getcwd()))

FULL_H5AD = ("SEAAD_Immune_HIP_MEC_LEC_ITG_MTG_FI_STG_DFC_AnG_V1C_"
             "RNAseq_final-nuclei.2026-06-22.h5ad")
DONOR_XLSX = os.path.join("MetaData",
                          "sea-ad_cohort_donor_metadata_072524.xlsx")
OUT_DIR = "Part1_PreChecks"
os.makedirs(OUT_DIR, exist_ok=True)

# Region grouping, must match Part 2
NEOC_REGIONS    = ["DFC", "MTG", "AnG", "STG", "V1C", "ITG"]
NONNEOC_REGIONS = ["HIP", "MEC", "LEC", "FI"]
ALL_EXPECTED    = NEOC_REGIONS + NONNEOC_REGIONS
CAUDATE_CODES   = ["CaH", "CaB", "CaT"]

# Regions where all 84 donors were profiled, per the SEA-AD documentation
FULL_COVERAGE_REGIONS = ["MTG", "DFC", "MEC"]

# =============================================================================
# Errata table, from sea-ad_sage_change_log_063026.pdf
# (specimen_id, needed_fix, donor_correct, donor_org, sex_corr, sex_org,
#  tissue_corr, tissue_org)
# =============================================================================
ERRATA = [
    ("M2AT_220121_206_D01",      "donor",        "H21.33.047", "H18.30.002", "M", "M", "MTG", "MTG"),
    ("M1AT_220902_207_B01",      "donor",        "H20.33.015", "H20.33.013", "M", "M", "V1C", "V1C"),
    ("M1AT_220902_207_C01",      "donor",        "H20.33.013", "H20.33.015", "M", "M", "V1C", "V1C"),
    ("M2AT_221007_221_A01",      "donor",        "H21.33.019", "H21.33.003", "M", "M", "V1C", "V1C"),
    ("M2AT_221007_221_C01",      "donor",        "H21.33.003", "H21.33.019", "M", "M", "V1C", "V1C"),
    ("M2AT_221024_206_A01",      "donor",        "H20.33.031", "H20.33.028", "F", "F", "V1C", "V1C"),
    ("M2AT_221021_206_C01",      "donor",        "H20.33.017", "H20.33.015", "M", "M", "CaH", "CaH"),
    ("M2AT_221021_206_B01",      "donor",        "H20.33.015", "H20.33.017", "M", "M", "CaH", "CaH"),
    ("1406_B09",                 "donor",        "H20.33.015", "H20.33.013", "M", "M", "V1C", "V1C"),
    ("1406_C09",                 "donor",        "H20.33.013", "H20.33.015", "M", "M", "V1C", "V1C"),
    ("1476_C06",                 "donor",        "H20.33.025", "H20.33.015", "M", "M", "CaH", "CaH"),
    ("1476_A06",                 "donor",        "H20.33.015", "H20.33.025", "M", "M", "CaH", "CaH"),
    ("M1HX_230112_401_2_A01",    "donor;sex",    "H20.33.020", "H20.33.031", "M", "F", "STG", "STG"),
    ("M1HX_230112_403_2_B01",    "donor;sex",    "H20.33.031", "H20.33.020", "F", "M", "STG", "STG"),
    ("M1HX_230119_431_2_A01",    "donor;sex",    "H21.33.003", "H21.33.002", "M", "F", "STG", "STG"),
    ("M1HX_230119_433_2_B01",    "donor;sex",    "H21.33.002", "H21.33.003", "F", "M", "STG", "STG"),
    ("M1HX_230130_401_1_A01",    "donor;sex",    "H21.33.004", "H21.33.018", "M", "F", "STG", "STG"),
    ("M1HX_230130_403_1_B01",    "donor;sex",    "H21.33.018", "H21.33.004", "F", "M", "STG", "STG"),
    ("M1TX_191118_107_F01",      "donor;sex",    "H18.30.002", "H200.1023",  "M", "F", "V1C", "V1C"),
    ("M1TX_191118_107_E01",      "donor;sex",    "H18.30.002", "H200.1023",  "M", "F", "V1C", "V1C"),
    ("M1TX_191118_107_D01",      "donor;sex",    "H18.30.002", "H200.1023",  "M", "F", "V1C", "V1C"),
    ("M2AT_221102_221_D01",      "donor;sex",    "H20.33.033", "H20.33.037", "M", "F", "V1C", "V1C"),
    ("M2AT_221102_221_B01",      "donor;sex",    "H20.33.037", "H20.33.033", "F", "M", "V1C", "V1C"),
    ("M2AT_221114_221_D01",      "donor;sex",    "H20.33.046", "H19.33.004", "M", "F", "CaH", "CaH"),
    ("M2AT_221114_221_C01",      "donor;sex",    "H19.33.004", "H20.33.046", "F", "M", "CaH", "CaH"),
    ("M2TX_210727_207_E01",      "donor;tissue", "H20.33.016", "H21.33.037", "F", "F", "DFC", "MTG"),
    ("M1TX_210723_159_G01",      "donor;tissue", "H21.33.037", "H20.33.016", "F", "F", "MTG", "DFC"),
    ("M2AT_221028_221_A01",      "tissue",       "H20.33.027", "H20.33.027", "F", "F", "CaH", "V1C"),
    ("M2AT_221028_221_C01",      "tissue",       "H20.33.001", "H20.33.001", "M", "M", "CaH", "V1C"),
    ("M2AT_221028_221_B01",      "tissue",       "H20.33.019", "H20.33.019", "F", "F", "CaH", "V1C"),
    ("M2TX_201106_201_B01_B01",  "tissue",       "H19.30.002", "H19.30.002", "M", "M", "CaB", "CaH"),
    ("M2TX_201120_202_B01_B01",  "tissue",       "H19.30.001", "H19.30.001", "M", "M", "CaB", "CaH"),
    ("M2TX_201106_201_A01_A01",  "tissue",       "H19.30.002", "H19.30.002", "M", "M", "CaB", "CaH"),
    ("M1TX_181114_001_B01_B01",  "tissue",       "H18.30.002", "H18.30.002", "M", "M", "CaB", "CaH"),
    ("M2TX_201120_202_A01_A01",  "tissue",       "H19.30.001", "H19.30.001", "M", "M", "CaB", "CaH"),
    ("M1TX_181114_001_A01_A01",  "tissue",       "H18.30.002", "H18.30.002", "M", "M", "CaB", "CaH"),
    ("988_D05",                  "tissue",       "H18.30.001", "H18.30.001", "F", "F", "CaT", "CaH"),
    ("897_E03",                  "tissue",       "H20.30.002", "H20.30.002", "M", "M", "CaB", "CaH"),
    ("891_A06",                  "tissue",       "H18.30.001", "H18.30.001", "F", "F", "CaB", "CaH"),
    ("988_C05",                  "tissue",       "H20.30.002", "H20.30.002", "M", "M", "CaT", "CaH"),
    ("930_E05",                  "tissue",       "H20.30.001", "H20.30.001", "M", "M", "CaB", "CaH"),
    ("980_C02",                  "tissue",       "H20.30.001", "H20.30.001", "M", "M", "CaT", "CaH"),
]
ERRATA_COLS = ["specimen_id", "needed_fix", "donor_correct", "donor_org",
               "sex_correct", "sex_org", "tissue_correct", "tissue_org"]
err = pd.DataFrame(ERRATA, columns=ERRATA_COLS)

# Specimens whose CORRECT tissue is one of the ten cortical regions studied
err["in_scope"] = err["tissue_correct"].isin(ALL_EXPECTED)
# Specimens wrongly labelled as an in-scope region (should be removed)
err["wrongly_in_scope"] = (err["tissue_org"].isin(ALL_EXPECTED)
                           & ~err["tissue_correct"].isin(ALL_EXPECTED))

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


def find_col(cols, keywords, exact_first=True):
    """Resolve a column name from candidate keywords. Exact match wins,
    then case-insensitive substring. Returns None if nothing matches."""
    cols = list(cols)
    if exact_first:
        for kw in keywords:
            for c in cols:
                if c == kw:
                    return c
    for kw in keywords:
        hits = [c for c in cols if kw.lower() in str(c).lower()]
        if hits:
            return sorted(hits, key=len)[0]
    return None


def norm_id(s):
    return str(s).strip().upper()


# =============================================================================
hdr("E0 - obs inventory and column resolution")
say(f"  file : {FULL_H5AD}")
ad = sc.read_h5ad(FULL_H5AD, backed="r")
obs = ad.obs.copy()
say(f"  shape: {ad.shape}")
say(f"  obs  : {obs.shape[0]:,} rows x {obs.shape[1]} columns")
say("")
say("  full obs column list:")
for i, c in enumerate(obs.columns):
    say(f"    [{i:>3}] {c}")

say("")
donor_col  = find_col(obs.columns, ["Donor ID", "donor_id", "Donor", "donor"])
region_col = find_col(obs.columns, ["Brain Region", "brain_region",
                                    "Region", "region", "tissue"])
stype_col  = find_col(obs.columns, ["Supertype", "supertype"])
sex_col    = find_col(obs.columns, ["Sex", "sex", "gender"])
lib_col    = find_col(obs.columns, ["library_prep", "library", "specimen",
                                    "sample_id", "sample"])

say("  resolved columns:")
for nm, c in [("donor", donor_col), ("region", region_col),
              ("supertype", stype_col), ("sex", sex_col),
              ("library/specimen", lib_col)]:
    say(f"    {nm:<18} -> {c}")

CRITICAL = {"donor": donor_col, "region": region_col, "supertype": stype_col}
missing_crit = [k for k, v in CRITICAL.items() if v is None]
if missing_crit:
    say("")
    say(f"  critical columns unresolved: {missing_crit}")
    with open(os.path.join(OUT_DIR, "S1d_errata_audit.txt"), "w",
              encoding="utf-8") as f:
        f.write("\n".join(_report))
    raise SystemExit("critical obs columns unresolved")

# =============================================================================
hdr("E1 - can the errata specimen IDs be matched to obs at all?")
err_ids = set(norm_id(x) for x in err["specimen_id"])
say(f"  errata specimens: {len(err_ids)}")
say(f"  of which CORRECT tissue is in the ten studied regions: "
    f"{int(err['in_scope'].sum())}")
say(f"  of which are WRONGLY labelled as a studied region: "
    f"{int(err['wrongly_in_scope'].sum())}")
say("")

# Find the obs column holding the errata specimen IDs: exact match first,
# substring match only on non-numeric columns.
match_report = []
skipped_numeric = []
for c in obs.columns:
    try:
        vals = obs[c].astype(str)
    except Exception:
        continue
    uniq = set(norm_id(v) for v in vals.unique())
    if len(uniq) > 20000:      # not an ID column
        continue
    exact = err_ids & uniq

    # is this a numeric column, or one whose values are all bare numbers?
    _is_num = pd.api.types.is_numeric_dtype(obs[c])
    _all_numeric_like = all(
        v.replace(".", "", 1).replace("-", "", 1).isdigit()
        for v in list(uniq)[:200] if v not in ("NAN", "NONE", ""))
    numericish = _is_num or _all_numeric_like

    sub = set()
    if not exact and not numericish:
        for e in err_ids:
            for u in uniq:
                # require an underscore-bearing ID for substring matching;
                # bare-number specimen IDs are matched by exact only
                if "_" in e and len(u) > 3 and (e in u or u in e):
                    sub.add(e)
                    break
    if not exact and numericish:
        skipped_numeric.append(c)
    if exact or sub:
        match_report.append((c, len(uniq), len(exact), len(sub)))

if skipped_numeric:
    say(f"  numeric columns excluded from substring matching: "
        f"{len(skipped_numeric)} "
        f"(e.g. {skipped_numeric[:4]})")
    say("")

if match_report:
    say("  columns containing errata specimen IDs:")
    say(f"    {'column':<40} {'n_uniq':>8} {'exact':>7} {'substr':>7}")
    for c, nu, ne, ns in sorted(match_report, key=lambda t: -t[2]):
        say(f"    {c:<40} {nu:>8,} {ne:>7} {ns:>7}")
    ID_COL = sorted(match_report, key=lambda t: -t[2])[0][0]
    ID_EXACT = sorted(match_report, key=lambda t: -t[2])[0][2] > 0
    say("")
    say(f"  specimen key: '{ID_COL}' "
        f"({'exact match' if ID_EXACT else 'substring match'})")
else:
    ID_COL = None
    ID_EXACT = False
    say("  no obs column contains an errata specimen ID; E2 skipped")
    if lib_col:
        say(f"  example values of '{lib_col}':")
        for v in list(obs[lib_col].astype(str).unique())[:5]:
            say(f"    {v}")

# =============================================================================
hdr("E2 - does obs match the corrected or the original values?")
CORRECTIONS_APPLIED = None
MATCHED_SPECIMENS = set()
if ID_COL is None:
    say("  skipped - no specimen join available (see E1)")
else:
    ids_norm = obs[ID_COL].astype(str).map(norm_id)
    rows = []
    for _, r in err.iterrows():
        sid = norm_id(r["specimen_id"])
        if ID_EXACT:
            m = (ids_norm == sid).values
        else:
            m = ids_norm.str.contains(sid, regex=False).values
        n = int(m.sum())
        if n == 0:
            continue
        MATCHED_SPECIMENS.add(r["specimen_id"])
        sub = obs.loc[m]
        d_obs = sorted(set(sub[donor_col].astype(str)))
        t_obs = sorted(set(sub[region_col].astype(str)))
        s_obs = (sorted(set(sub[sex_col].astype(str)))
                 if sex_col else ["<no sex col>"])
        d_ok = (len(d_obs) == 1 and d_obs[0] == r["donor_correct"])
        d_bad = (len(d_obs) == 1 and d_obs[0] == r["donor_org"])
        t_ok = (len(t_obs) == 1 and t_obs[0] == r["tissue_correct"])
        t_bad = (len(t_obs) == 1 and t_obs[0] == r["tissue_org"])
        rows.append(dict(specimen=r["specimen_id"], fix=r["needed_fix"],
                         n_nuclei=n,
                         donor_obs="|".join(d_obs),
                         donor_correct=r["donor_correct"],
                         donor_org=r["donor_org"],
                         donor_state=("CORRECT" if d_ok else
                                      "UNCORRECTED" if d_bad else "OTHER"),
                         tissue_obs="|".join(t_obs),
                         tissue_correct=r["tissue_correct"],
                         tissue_org=r["tissue_org"],
                         tissue_state=("CORRECT" if t_ok else
                                       "UNCORRECTED" if t_bad else "OTHER"),
                         sex_obs="|".join(s_obs),
                         sex_correct=r["sex_correct"],
                         sex_org=r["sex_org"]))
    if not rows:
        say("  none of the 42 errata specimens has nuclei in this file")
        CORRECTIONS_APPLIED = True
    else:
        aud = pd.DataFrame(rows)
        say(f"  {len(aud)} errata specimens have nuclei here "
            f"({int(aud['n_nuclei'].sum()):,} nuclei total)")
        say("")
        say(f"    {'specimen':<26}{'fix':<14}{'n':>7}  {'donor':<12}"
            f"{'tissue':<12}")
        for _, a in aud.iterrows():
            say(f"    {a['specimen']:<26}{a['fix']:<14}{a['n_nuclei']:>7,}  "
                f"{a['donor_state']:<12}{a['tissue_state']:<12}")
        n_unc = int(((aud["donor_state"] == "UNCORRECTED") |
                     (aud["tissue_state"] == "UNCORRECTED")).sum())
        n_cor = int(((aud["donor_state"] == "CORRECT") &
                     (aud["tissue_state"] == "CORRECT")).sum())
        n_oth = len(aud) - n_cor - n_unc
        say("")
        say(f"  specimens matching the corrected metadata   : {n_cor}")
        say(f"  specimens matching the original metadata    : {n_unc}")
        say(f"  specimens matching neither (indeterminate)  : {n_oth}")
        # A specimen matching neither value needs inspection, so it is not a pass.
        if n_unc > 0:
            CORRECTIONS_APPLIED = False
        elif n_oth > 0:
            CORRECTIONS_APPLIED = None
        else:
            CORRECTIONS_APPLIED = True
        say("")
        if n_oth and n_unc == 0:
            say(f"  indeterminate: {n_oth} specimen(s) match neither value")
        elif n_unc == 0:
            say("  corrections are in the file")
        else:
            say("  corrections are not in the file")
        aud.to_csv(os.path.join(OUT_DIR, "S1d_errata_specimen_audit.csv"),
                   index=False)
        say(f"  saved: {OUT_DIR}/S1d_errata_specimen_audit.csv")

        # nuclei affected, split by arm
        in_neoc = aud["tissue_obs"].isin(NEOC_REGIONS)
        in_non = aud["tissue_obs"].isin(NONNEOC_REGIONS)
        say("")
        say(f"  affected nuclei in NeoC regions    : "
            f"{int(aud.loc[in_neoc, 'n_nuclei'].sum()):,}")
        say(f"  affected nuclei in NonNeoC regions : "
            f"{int(aud.loc[in_non, 'n_nuclei'].sum()):,}")


# =============================================================================
# E2b - in-scope errata specimens with no nuclei here, classified by assay
# =============================================================================
hdr("E2b - in-scope errata specimens with no nuclei here")


def assay_class(sid):
    """Guess the assay from the specimen-ID prefix. SEA-AD specimen names
    encode it: *AT* = ATAC/multiome, *TX*/*HX* = 10x RNA, a bare number =
    plate-based (SSv4). Only RNA-class specimens could appear in this file."""
    pre = str(sid).split("_")[0]
    if pre.isdigit():
        return "plate"
    if "AT" in pre:
        return "ATAC"
    if "TX" in pre or "HX" in pre:
        return "RNA"
    return "other"


_donors_obs = set(obs[donor_col].astype(str))
_scope = err[err["in_scope"]].copy()
_scope["assay"] = _scope["specimen_id"].map(assay_class)
_scope["tested"] = _scope["specimen_id"].isin(MATCHED_SPECIMENS)
_untested = _scope[~_scope["tested"]]

say(f"  in-scope errata specimens (corrected tissue is a studied region) : "
    f"{len(_scope)}")
say(f"    tested in E2 (joined to obs)                                   : "
    f"{int(_scope['tested'].sum())}")
say(f"    not testable (no nuclei under this specimen ID)                : "
    f"{len(_untested)}")
say("")

# Join rate by assay class.
say("  join rate by assay class (in-scope specimens):")
say(f"    {'assay':<10}{'in errata':>11}{'joined to obs':>15}")
for _a, _grp in _scope.groupby("assay"):
    say(f"    {_a:<10}{len(_grp):>11}{int(_grp['tested'].sum()):>15}")
say("")

if len(_untested) == 0:
    say("  every in-scope errata specimen was tested")
else:
    say(f"    {'specimen':<26}{'assay':<8}{'fix':<12}{'tissue':<7}"
        f"{'donor(correct)':<15}{'in file?':<9}{'nuclei':>9}")
    _resid = []
    for _, r in _untested.iterrows():
        _dc, _tc = str(r["donor_correct"]), str(r["tissue_correct"])
        _present = _dc in _donors_obs
        _n_there = int(((obs[donor_col].astype(str) == _dc)
                        & (obs[region_col].astype(str) == _tc)).sum()) \
            if _present else 0
        _flag = (r["assay"] == "RNA" and _present and _n_there > 0)
        say(f"    {str(r['specimen_id']):<26}{str(r['assay']):<8}"
            f"{str(r['needed_fix']):<12}{_tc:<7}{_dc:<15}"
            f"{('yes' if _present else 'no'):<9}{_n_there:>9,}"
            f"{'  <- residual' if _flag else ''}")
        if _flag:
            _resid.append(str(r["specimen_id"]))

    _n_nonrna = int((_untested["assay"] != "RNA").sum())
    _n_absent = int(((_untested["assay"] == "RNA")
                     & ~_untested["donor_correct"].astype(str)
                     .isin(_donors_obs)).sum())
    say("")
    say(f"  not an RNA specimen, so cannot be in this file      : {_n_nonrna}")
    say(f"  RNA specimen but its corrected donor is absent here : {_n_absent}")
    say(f"  residual (RNA, donor present, unverified)           : {len(_resid)}")
    if not _resid:
        say(f"  verified on {int(_scope['tested'].sum())} of "
            f"{len(_scope)} in-scope specimens; no residual")
    else:
        say("  residual present: flagged specimens were not verified by E2")

    # Specimens wrongly labelled as a studied region; E5 checks they are absent.
    _wrong = err[err["wrongly_in_scope"]]
    if len(_wrong):
        say("")
        say(f"  specimens wrongly labelled as a studied region: {len(_wrong)}")
        for _, r in _wrong.iterrows():
            say(f"    {str(r['specimen_id']):<26}"
                f"{str(r['tissue_org'])} (labelled) -> "
                f"{str(r['tissue_correct'])} (really)")

# =============================================================================
hdr("E3 - donor-ID sanity against the 84-donor cohort table")
if not os.path.exists(DONOR_XLSX):
    say(f"  donor table not found at {DONOR_XLSX}")
    dtab = None
else:
    dtab = pd.read_excel(DONOR_XLSX)
    say(f"  donor table: {dtab.shape[0]} rows x {dtab.shape[1]} columns")
    say(f"  unique Donor ID in table : {dtab['Donor ID'].nunique()}")

obs_donors = set(obs[donor_col].astype(str).unique())
say(f"  unique donors in obs     : {len(obs_donors)}")

if dtab is not None:
    tab_donors = set(dtab["Donor ID"].astype(str))
    orphans = sorted(obs_donors - tab_donors)
    say(f"  donors in obs but not in the 84-donor table : {len(orphans)}")
    for o in orphans:
        n = int((obs[donor_col].astype(str) == o).sum())
        regs = sorted(set(obs.loc[obs[donor_col].astype(str) == o,
                                  region_col].astype(str)))
        note = ""
        if o == "H200.1023":
            note = "  <- MALFORMED ID named in the errata"
        elif o.startswith(("H18.30", "H19.30", "H20.30")):
            note = "  <- young neurotypical REFERENCE donor, not the aged cohort"
        say(f"    {o:<16} {n:>8,} nuclei  regions={regs}{note}")
    if orphans:
        say("  these donors cannot be joined to the donor metadata")
    else:
        say("  all obs donors are in the 84-donor table")

# =============================================================================
hdr("E4 - sex consistency per donor")
if sex_col is None:
    say("  no sex column in obs - skipped")
elif dtab is None:
    say("  no donor table - skipped")
else:
    smap = (dtab.set_index(dtab["Donor ID"].astype(str))["Sex"]
            .astype(str).str.strip().str.upper().str[0].to_dict())
    g = (obs[[donor_col, sex_col]].astype(str)
         .groupby(donor_col)[sex_col]
         .agg(lambda s: sorted(set(s.str.strip().str.upper().str[0]))))
    multi = {d: v for d, v in g.items() if len(v) > 1}
    say(f"  donors with more than one sex value in obs : {len(multi)}")
    for d, v in list(multi.items())[:10]:
        say(f"    {d:<16} {v}")
    mism = []
    for d, v in g.items():
        if len(v) == 1 and d in smap and v[0] != smap[d]:
            mism.append((d, v[0], smap[d]))
    say(f"  donors whose obs sex disagrees with the table : {len(mism)}")
    for d, a, b in mism[:20]:
        n = int((obs[donor_col].astype(str) == d).sum())
        regs = sorted(set(obs.loc[obs[donor_col].astype(str) == d,
                                  region_col].astype(str)))
        say(f"    {d:<16} obs={a}  table={b}  {n:>8,} nuclei  regions={regs}")
    if mism or multi:
        say("  sex mismatch: check these donors against the errata table")
    else:
        say("  sex is consistent within donors and with the table")

# =============================================================================
hdr("E5 - region inventory")
regs = sorted(set(obs[region_col].astype(str)))
say(f"  distinct region values ({len(regs)}): {regs}")
say("")
missing_reg = [r for r in ALL_EXPECTED if r not in regs]
extra_reg = [r for r in regs if r not in ALL_EXPECTED]
say(f"  expected regions absent : {missing_reg if missing_reg else 'none'}")
say(f"  unexpected regions      : {extra_reg if extra_reg else 'none'}")
caud = [r for r in regs if r in CAUDATE_CODES]
if caud:
    n = int(obs[region_col].astype(str).isin(caud).sum())
    say(f"  caudate nuclei present: {caud}, {n:,} nuclei")
else:
    say("  no caudate codes present")
say("")
say("  nuclei per region:")
vc = obs[region_col].astype(str).value_counts()
for r in ALL_EXPECTED:
    arm = "NeoC" if r in NEOC_REGIONS else "NonNeoC"
    say(f"    {r:<6} {arm:<8} {int(vc.get(r, 0)):>9,}")
for r in extra_reg:
    say(f"    {r:<6} {'OTHER':<8} {int(vc.get(r, 0)):>9,}")

# =============================================================================
hdr("E6 - donor x region coverage")
ct = pd.crosstab(obs[donor_col].astype(str), obs[region_col].astype(str))
ct = ct.reindex(columns=[r for r in ALL_EXPECTED if r in ct.columns],
                fill_value=0)
ct["n_regions"] = (ct > 0).sum(axis=1)
ct.to_csv(os.path.join(OUT_DIR, "S1d_donor_region_coverage.csv"))
say(f"  donors x regions matrix saved: "
    f"{OUT_DIR}/S1d_donor_region_coverage.csv")
say("")
say("  donors with >=1 nucleus per region:")
for r in ALL_EXPECTED:
    if r in ct.columns:
        nd = int((ct[r] > 0).sum())
        arm = "NeoC" if r in NEOC_REGIONS else "NonNeoC"
        flag = ""
        if r in FULL_COVERAGE_REGIONS:
            flag = "  <- documented as all 84 donors"
        say(f"    {r:<6} {arm:<8} {nd:>4} donors{flag}")
say("")
full_sets = [set(ct.index[ct[r] > 0]) for r in FULL_COVERAGE_REGIONS
             if r in ct.columns]
if full_sets:
    full_donors = set.union(*full_sets)
    say(f"  donors appearing in MTG/DFC/MEC (union) : {len(full_donors)}")
spec_regions = [r for r in ALL_EXPECTED if r not in FULL_COVERAGE_REGIONS
                and r in ct.columns]
if spec_regions:
    spec_sets = [set(ct.index[ct[r] > 0]) for r in spec_regions]
    spec_union = set.union(*spec_sets)
    spec_inter = set.intersection(*spec_sets)
    say(f"  donors appearing in any AD-Spectrum region (union) : "
        f"{len(spec_union)}")
    say(f"  donors appearing in all AD-Spectrum regions        : "
        f"{len(spec_inter)}")
say("")
say("  distribution of regions-per-donor:")
for k, v in ct["n_regions"].value_counts().sort_index().items():
    say(f"    {int(k):>2} regions : {int(v):>4} donors")

# == Coverage on the nuclei flagged 'Used in analysis' ==========================
_used_col = find_col(obs.columns, ["Used in analysis", "used_in_analysis"])
say("")
if _used_col is None:
    say("  'Used in analysis' column absent")
else:
    _u = obs[_used_col].astype(str).str.strip().str.lower()
    _keep = ~_u.isin(["false", "0", "no", "n", "nan", "none", ""])
    say(f"  '{_used_col}': {int(_keep.sum()):,} of {len(_keep):,} nuclei flagged "
        f"for analysis ({100*_keep.mean():.3f}%)")
    if _keep.all():
        say("  all nuclei flagged")
    else:
        ct2 = pd.crosstab(obs.loc[_keep.values, donor_col].astype(str),
                          obs.loc[_keep.values, region_col].astype(str))
        ct2 = ct2.reindex(columns=[r for r in ALL_EXPECTED
                                   if r in ct2.columns], fill_value=0)
        say("")
        say("  donor counts, all nuclei vs flagged nuclei:")
        say(f"    {'region':<8}{'all':>8}{'flagged':>9}{'delta':>7}")
        for r in ALL_EXPECTED:
            if r in ct.columns:
                _a = int((ct[r] > 0).sum())
                _b = int((ct2[r] > 0).sum()) if r in ct2.columns else 0
                _d = _b - _a
                say(f"    {r:<8}{_a:>8}{_b:>9}{_d:>+7}")
        ct2.to_csv(os.path.join(OUT_DIR,
                                "S1d_donor_region_coverage_flagged.csv"))
        say("")
        say(f"  saved: {OUT_DIR}/S1d_donor_region_coverage_flagged.csv")

# =============================================================================
hdr("E7 - VERDICT")
say(f"  errata corrections applied to this file : {CORRECTIONS_APPLIED}")
say(f"  caudate contamination                   : {bool(caud)}")
say(f"  orphan donors (not in the 84)           : "
    f"{len(orphans) if dtab is not None else 'unknown'}")
say("")
if CORRECTIONS_APPLIED is True and not caud and (
        dtab is not None and not orphans):
    say("  verdict: clean")
elif CORRECTIONS_APPLIED is False:
    say("  verdict: uncorrected donor or region labels present")
else:
    say("  verdict: partial; see E1-E6")

with open(os.path.join(OUT_DIR, "S1d_errata_audit.txt"), "w",
          encoding="utf-8") as f:
    f.write("\n".join(_report))
say("")
say(f"  saved: {OUT_DIR}/S1d_errata_audit.txt")
say("")
say("=" * 74)
say("  PART 1d COMPLETE")
say("=" * 74)
