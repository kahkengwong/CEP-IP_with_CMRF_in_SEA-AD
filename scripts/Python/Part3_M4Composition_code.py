# =============================================================================
# Part3_M4Composition_code.py
# SEA-AD CEP-IP/CMRF, Part 3, supertype composition of the M4 set (Figure 1B)
#   Computes the supertype composition of the M4 selected set with the Part 2
#   Section 3 code, after reproducing the M1 file on disk, and the composition
#   under all six scorers. Reads the h5ad obs only.
#
# Usage (from the SEA-AD folder, or set SEAAD_ROOT), after Part 2
#   python3 Part3_M4Composition_code.py
#
# Outputs (<arm>/Part2_<arm>_M4/top_3.0pct/): <arm>_S1_supertype_enrichment_M4.csv,
#   <arm>_S7_scorer_composition.csv; log Part3_M4Composition_log.txt.
# =============================================================================

import os
import random
from importlib.metadata import version

import numpy as np
import pandas as pd
import anndata as ad

# == Configuration ============================================================
_BASE_CANDIDATES = [os.environ.get("SEAAD_ROOT", os.getcwd())]
BASE = next((p for p in _BASE_CANDIDATES if os.path.isdir(p)), None)
if BASE is None:
    raise SystemExit(f"SEA-AD folder not found in any of: {_BASE_CANDIDATES}")

H5AD = os.path.join(BASE, "SEAAD_Immune_HIP_MEC_LEC_ITG_MTG_FI_STG_DFC_AnG_"
                          "V1C_RNAseq_final-nuclei.2026-06-22.h5ad")
ARM_FOLDER = {"NeoC": "[1]_NeoC", "NonNeoC": "[2]_Non-NeoC"}
TOP_PCT = 0.030
N_BOOT_ENR = 2000                       # Part 2 Section 3
DONOR_COL = "Donor ID"                  # Part 2 donor_col
STYPE_COL = "Supertype"                 # Part 2 stype_col

INTERSECTION_PAIRS = {"M4_M1andM2": ("M1_FC", "M2_Hybrid"),
                      "M5_M1andM3": ("M1_FC", "M3_pRank"),
                      "M6_M2andM3": ("M2_Hybrid", "M3_pRank")}
SCORERS = ["M1_FC", "M2_Hybrid", "M3_pRank",
           "M4_M1andM2", "M5_M1andM3", "M6_M2andM3"]

# Reported values, used as checks only.
EXPECT_N_ARM = {"NeoC": 134088, "NonNeoC": 106563}
EXPECT_M4 = {
    "NeoC": {"Micro-PVM_4-SEAAD": 357, "Micro-PVM_3-SEAAD": 2626,
             "Micro-PVM_2_3-SEAAD": 227, "Micro-PVM_1": 26,
             "Micro-PVM_2_1-SEAAD": 4, "Micro-PVM_2": 170,
             "Lymphocyte": 2, "Monocyte": 0},
    "NonNeoC": {"Micro-PVM_4-SEAAD": 103, "Micro-PVM_3-SEAAD": 2329,
                "Micro-PVM_2_3-SEAAD": 144, "Micro-PVM_2_1-SEAAD": 2,
                "Micro-PVM_1": 8, "Micro-PVM_2": 118,
                "Lymphocyte": 0, "Monocyte": 0}}
# Micro-PVM_4-SEAAD % in the top 3.0% per scorer, from the Part 2 figures
EXPECT_PVM4_PCT = {
    "NeoC":    {"M1_FC": 30.7, "M2_Hybrid": 31.5, "M3_pRank": 40.7,
                "M4_M1andM2": 26.9, "M5_M1andM3": 30.3, "M6_M2andM3": 29.5},
    "NonNeoC": {"M1_FC": 26.2, "M2_Hybrid": 27.4, "M3_pRank": 35.3,
                "M4_M1andM2": 23.3, "M5_M1andM3": 25.6, "M6_M2andM3": 26.2}}
# Monocyte nuclei in the top set, from the same figures
EXPECT_MONO = {"NeoC":    {"M2_Hybrid": 1, "M3_pRank": 2, "M6_M2andM3": 1},
               "NonNeoC": {"M2_Hybrid": 2, "M3_pRank": 3, "M6_M2andM3": 2}}

LOG = []


def say(s=""):
    print(s)
    LOG.append(str(s))


def check(ok, what):
    say(f"  [{'PASS' if ok else 'FAIL'}] {what}")
    if not ok:
        write_log()
        raise SystemExit(f"Check failed: {what}")


def write_log():
    with open(os.path.join(BASE, "Part3_M4Composition_log.txt"), "w",
              encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")


# == Part 2 selection and Section 3 code ======================================
def top_mask(scores, mkey):
    """Part 2 _scorer_top: base scorers by percentile, M4-M6 as intersections."""
    if mkey in INTERSECTION_PAIRS:
        a, b = INTERSECTION_PAIRS[mkey]
        return top_mask(scores, a) & top_mask(scores, b)
    s = scores[f"score_{mkey}"].to_numpy()
    return s >= np.percentile(s, (1 - TOP_PCT) * 100)


def section3(obs, mask, observed):
    """Part 2 Section 3, unchanged except for the explicit groupby 'observed'."""
    np.random.seed(42); random.seed(42)
    stype_all = obs[STYPE_COL].value_counts()
    stype_top = obs[mask][STYPE_COL].value_counts()
    rows = []
    for st in stype_all.index:
        n_st_all = int(stype_all.get(st, 0))
        n_st_top = int(stype_top.get(st, 0))
        enr = (n_st_top / max(n_st_all, 1)) / TOP_PCT
        st_top_msk = mask & (obs[STYPE_COL] == st).values
        n_st_don = int(obs.loc[st_top_msk, DONOR_COL].astype(str).nunique())
        st_mask_all = (obs[STYPE_COL] == st).values
        per_d_all = obs[st_mask_all].groupby(DONOR_COL, observed=observed).size()
        per_d_top = obs[st_mask_all & mask].groupby(DONOR_COL, observed=observed).size()
        per_d = pd.DataFrame({"n_all": per_d_all, "n_top": per_d_top}).fillna(0)
        boot = np.empty(N_BOOT_ENR)
        for bi in range(N_BOOT_ENR):
            bidx = np.random.choice(len(per_d), len(per_d), replace=True)
            b_all = per_d["n_all"].iloc[bidx].sum()
            b_top = per_d["n_top"].iloc[bidx].sum()
            boot[bi] = (b_top / max(b_all, 1)) / TOP_PCT
        lo, hi = np.percentile(boot, [2.5, 97.5])
        rows.append({"supertype": st, "n_total": n_st_all, "n_top": n_st_top,
                     "n_donors_top": n_st_don,
                     "pct_top": round(100 * n_st_top / max(n_st_all, 1), 2),
                     "enrichment_fold": round(enr, 3),
                     "boot_ci_lo": round(lo, 3), "boot_ci_hi": round(hi, 3)})
    return (pd.DataFrame(rows)
            .sort_values("enrichment_fold", ascending=False)
            .reset_index(drop=True))


def same(a, b, cols):
    a = a.set_index("supertype")[cols].sort_index()
    b = b.set_index("supertype")[cols].sort_index()
    return a.index.equals(b.index) and bool((a.values == b.values).all())


# == Run ======================================================================
say("Part3_M4Composition_code.py")
say(f"BASE : {BASE}")
say(f"pandas {pd.__version__}, numpy {np.__version__}, anndata {version('anndata')}")
if not os.path.isfile(H5AD):
    raise SystemExit(f"h5ad not found: {H5AD}")
adata = ad.read_h5ad(H5AD, backed="r")
obs_all = adata.obs[[DONOR_COL, STYPE_COL]].copy()
adata.file.close()
say(f"h5ad obs: {len(obs_all):,} nuclei")

for arm, folder in ARM_FOLDER.items():
    say("")
    say("=" * 70)
    say(f"  {arm}")
    say("=" * 70)
    top_dir = os.path.join(BASE, folder, f"Part2_{arm}_M4", "top_3.0pct")
    s0_path = os.path.join(top_dir, "S0_cell_scores_all_methods.csv")
    s1_candidates = [os.path.join(top_dir, "S1_supertype_enrichment.csv"),
                     os.path.join(top_dir, f"{arm}_S1_supertype_enrichment.csv")]
    s1_path = next((p for p in s1_candidates if os.path.isfile(p)), None)
    if not os.path.isfile(s0_path) or s1_path is None:
        raise SystemExit(f"missing input in {top_dir}")
    say(f"  scores : {s0_path}")
    say(f"  M1 file: {s1_path}")

    scores = pd.read_csv(s0_path, index_col=0, float_precision="round_trip")
    check(len(scores) == EXPECT_N_ARM[arm],
          f"{len(scores):,} nuclei in S0 (expected {EXPECT_N_ARM[arm]:,})")
    check(scores.index.isin(obs_all.index).all(), "every S0 barcode is in the h5ad obs")
    obs = obs_all.loc[scores.index]

    # Step 2: reproduce the M1 file on disk
    m1_disk = pd.read_csv(s1_path)
    m1_mask = top_mask(scores, "M1_FC")
    base_cols = ["n_total", "n_top", "n_donors_top", "pct_top", "enrichment_fold"]
    ci_cols = ["boot_ci_lo", "boot_ci_hi"]
    m1_obs_false = section3(obs, m1_mask, observed=False)
    check(same(m1_obs_false, m1_disk, base_cols),
          "M1 counts, donors, % and fold reproduce the file on disk exactly")
    observed_mode = None
    if same(m1_obs_false, m1_disk, ci_cols):
        observed_mode = False
    else:
        m1_obs_true = section3(obs, m1_mask, observed=True)
        if same(m1_obs_true, m1_disk, ci_cols):
            observed_mode = True
    check(observed_mode is not None,
          "M1 bootstrap CIs reproduce the file on disk exactly "
          f"(groupby observed={observed_mode})")

    # Step 3: M4 with the same code, seed and grouping
    m4_mask = top_mask(scores, "M4_M1andM2")
    m4 = section3(obs, m4_mask, observed=observed_mode)
    got = dict(zip(m4["supertype"], m4["n_top"]))
    check(int(m4_mask.sum()) == sum(EXPECT_M4[arm].values()),
          f"{int(m4_mask.sum()):,} nuclei in the M4 set")
    check(got == EXPECT_M4[arm], "M4 per-supertype counts match Figure 1B")
    out1 = os.path.join(top_dir, f"{arm}_S1_supertype_enrichment_M4.csv")
    m4.to_csv(out1, index=False)
    for _, r in m4.iterrows():
        say(f"    {r.supertype:<22} {r.n_top:>5} / {r.n_total:>6}  donors {r.n_donors_top:>3}  "
            f"{r.pct_top:>6.2f}%  {r.enrichment_fold:>6.3f}x  "
            f"[{r.boot_ci_lo:.3f}, {r.boot_ci_hi:.3f}]")
    say(f"  saved: {out1}")

    # Step 4: six-scorer composition
    rows = []
    for mk in SCORERS:
        msk = top_mask(scores, mk)
        n_all = obs[STYPE_COL].value_counts()
        n_top = obs[msk][STYPE_COL].value_counts()
        for st in n_all.index:
            nt, na = int(n_top.get(st, 0)), int(n_all[st])
            rows.append({"arm": arm, "scorer": mk, "supertype": st,
                         "n_total": na, "n_top": nt,
                         "pct_top": round(100 * nt / na, 2),
                         "enrichment_fold": round(nt / na / TOP_PCT, 3)})
    comp = pd.DataFrame(rows)
    for mk in SCORERS:
        # Compared at the figure's precision (one decimal place).
        r = comp[(comp.scorer == mk) & (comp.supertype == "Micro-PVM_4-SEAAD")].iloc[0]
        v = 100 * r.n_top / r.n_total
        check(f"{v:.1f}" == f"{EXPECT_PVM4_PCT[arm][mk]:.1f}",
              f"{mk}: Micro-PVM_4-SEAAD {v:.3f}% (Part 2 figure {EXPECT_PVM4_PCT[arm][mk]}%)")
    for mk, n in EXPECT_MONO[arm].items():
        v = int(comp[(comp.scorer == mk) & (comp.supertype == "Monocyte")].n_top.iloc[0])
        check(v == n, f"{mk}: {v} Monocyte nuclei in the top set (Part 2 figure {n})")
    out2 = os.path.join(top_dir, f"{arm}_S7_scorer_composition.csv")
    comp.to_csv(out2, index=False)
    say(f"  saved: {out2}")

say("")
say("All checks passed.")
write_log()
