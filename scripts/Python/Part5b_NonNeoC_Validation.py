# =============================================================================
# Part5b_NonNeoC_Validation.py
# SEA-AD CEP-IP/CMRF, Part 5b, Mic.13 equivalence validation, non-neocortical arm (Figure 1)
#   Tests whether the M4-top PVM_3 cells match Mic.13 (Green et al.) using the
#   Part 5 DEGs: hypergeometric overlap with the held-out Mic.13 markers (genes
#   not used by the Part 2 scorer) and with the DAM2 signature (Keren-Shaul et
#   al. 2017). In-scorer Mic.13 genes are reported separately as circular.
#
# Usage (from the SEA-AD folder, or set SEAAD_ROOT), after Part5_NonNeoC_DEGs.py
#   python3 Part5b_NonNeoC_Validation.py
#
# Outputs ([2]_Non-NeoC/Part5b_NonNeoC_Validation/): Part5b_NonNeoC_Mic13_validation.xlsx.
# =============================================================================

import os, warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

from openpyxl import Workbook
from openpyxl.styles import PatternFill, Font, Alignment, Border, Side
from openpyxl.utils import get_column_letter

import random as _random_stdlib
GLOBAL_SEED = 42
np.random.seed(GLOBAL_SEED)
_random_stdlib.seed(GLOBAL_SEED)

# =============================================================================
# SECTION 0: Arm config and folder convention
# =============================================================================
# BASE_DIR: the SEA-AD folder.
_BASE_CANDIDATES = [os.environ.get("SEAAD_ROOT", os.getcwd())]
BASE_DIR = next((p for p in _BASE_CANDIDATES if os.path.isdir(p)), None)
if BASE_DIR is None:
    raise SystemExit(f"SEA-AD folder not found in any of: {_BASE_CANDIDATES}")
ARM_FOLDER = r"[2]_Non-NeoC"
ARM_TAG    = "NonNeoC"
ARM_NAME   = "non-neocortical"
PART_TAG   = "Part5b"
DIR_TAG    = "Validation"
ARM_ROOT   = os.path.join(BASE_DIR, ARM_FOLDER)
ARM_DIR    = os.path.join(ARM_ROOT, f"{PART_TAG}_{ARM_TAG}_{DIR_TAG}")
OUT_DIR    = ARM_DIR
os.makedirs(ARM_DIR, exist_ok=True)
print(f"[OUT] arm folder : {ARM_ROOT}")
print(f"[OUT] this part  : {ARM_DIR}")

print("=" * 70)
print("Part 5b - Mic.13 Equivalence Validation  [NonNeoC arm]")
print("=" * 70)

# == Part 5 DEG input =========================================================
PART5_DIR_TAG   = "DEGs"
PART5_DIR       = os.path.join(ARM_ROOT, f"Part5_{ARM_TAG}_{PART5_DIR_TAG}")
# Accept the original and the arm-prefixed workbook names.
_DEG_NAMES = ["DEG_master_results.xlsx", f"{ARM_TAG}-DEG_master_results.xlsx",
              f"{ARM_TAG}_DEG_master_results.xlsx"]
DEG_MASTER_FILE = next((os.path.join(PART5_DIR, n) for n in _DEG_NAMES
                        if os.path.isfile(os.path.join(PART5_DIR, n))),
                       os.path.join(PART5_DIR, "DEG_master_results.xlsx"))
# Fallback: if running from the same directory as the DEG file
if not os.path.isfile(DEG_MASTER_FILE):
    DEG_MASTER_FILE = os.path.join(BASE_DIR, "DEG_master_results.xlsx")
if not os.path.isfile(DEG_MASTER_FILE):
    raise FileNotFoundError(
        f"\n  FATAL: DEG master results not found at:\n"
        f"    {os.path.join(ARM_ROOT, f'Part5_{ARM_TAG}_{PART5_DIR_TAG}', 'DEG_master_results.xlsx')}\n"
        f"    {os.path.join(BASE_DIR, 'DEG_master_results.xlsx')}\n"
        f"  Run Part5_{ARM_TAG}_DEGs.py first.")
print(f"\n  DEG input: {DEG_MASTER_FILE}")

# == Mic.13 gene list (Green et al. 2024, Extended Data Fig 2c) ===============
MIC13_GENES_FILE = os.path.join(BASE_DIR, "Mic13_genes.txt")
if not os.path.isfile(MIC13_GENES_FILE):
    MIC13_GENES_FILE = "Mic13_genes.txt"
print(f"  Mic.13 gene list: {MIC13_GENES_FILE}")

# == Arm regions (for reading the right DEG sheets) ===========================
ARM_REGIONS  = ["HIP", "MEC", "LEC", "FI"]
REGION_ORDER = ARM_REGIONS

# =============================================================================
# SECTION 1: Define the scorer panel and the held-out set
# =============================================================================
print("\n" + "=" * 70)
print("SECTION 1: Scorer panel vs held-out gene split")
print("=" * 70)

# Scorer panel: the keys of Part 2's MIC13_UP (85) and MIC13_DOWN (48).
# No mitochondrial gene is scored.
PART2_MIC13_UP = {
    "ADAMTS17", "ADARB1", "ADK", "ALCAM", "APOE", "ARID5B", "ASAH1", "ASAP1", "ATG7",
    "BCL2", "BICD1", "CACNA1A", "CADM1", "CD83", "CDK6", "CHSY3", "CLEC5A", "CPM",
    "CTSD", "CYTH3", "DENND4C", "DIRC3", "DPYD", "DSCAM", "DTNA", "ELL2", "FLT1",
    "FMNL2", "FNDC3B", "FNIP2", "FRMD4B", "GAS7", "GLDN", "GPNMB", "HDAC9", "IPCEF1",
    "ITGAX", "LIMS1", "LIPA", "LNCAROD", "LPL", "LPP", "MAFB", "MALAT1", "MAML2",
    "MITF", "MSR1", "MTSS1", "MYO1E", "NCK2", "NHSL1", "NPL", "OLR1", "PADI2", "PDE3B",
    "PDE8A", "PLXDC2", "PPARG", "PRKCE", "PTPRG", "PTPRM", "RALGAPA2", "RASGEF1B",
    "RASGRP3", "RGCC", "RTTN", "SAMD4A", "SAMSN1", "SCIN", "SDK1", "SLC11A1", "SLC1A3",
    "SOCS6", "SPP1", "STARD13", "SYNDIG1", "TANC2", "TGFBR1", "TMEM117", "TPRG1",
    "TRHDE", "WIPF3", "XYLT1", "ZNF618", "ZNF804A"
}
PART2_MIC13_DOWN = {
    "ABCC4", "ADAM28", "AOAH", "CCDC26", "CPED1", "CSGALNACT1", "CTTNBP2", "CX3CR1",
    "ELMO1", "EPB41L3", "F13A1", "FAM13A", "FAM177B", "FCHSD2", "FOXP2", "FRMD4A",
    "FYB1", "GRID2", "HAMP", "HSPA1A", "IGSF21", "IL6ST", "INPP5D", "ITPR2", "MAF",
    "MIR646HG", "MIR99AHG", "MRC1", "MS4A6A", "MVB12B", "NAMPT", "NAV3", "PICALM",
    "PREX1", "PRKCA", "RASGEF1C", "RPS27", "SERPINB9", "SLC9A9", "SRGAP1", "SRGAP2",
    "SRGAP2B", "SRGN", "SSBP2", "ST6GALNAC3", "TBXAS1", "TLN2", "TMEM163"
}
assert len(PART2_MIC13_UP) == 85 and len(PART2_MIC13_DOWN) == 48, "Part 2 scorer lists changed"

# Module 1 mitochondrial genes, for reporting only: Part 2 does not score them
MODULE1_MT_GENES = {
    "MT-CO3",    # GOI
    "MT-CO1", "MT-CO2", "MT-ATP6", "MT-CYB", "MT-ND1", "MT-ND2",
    "MT-ND3", "MT-ND4",
}

# Full scorer gene set = every gene Part 2 scores
SCORER_GENES = PART2_MIC13_UP | PART2_MIC13_DOWN
print(f"  Scorer genes (circular, {len(SCORER_GENES)}): "
      f"{len(PART2_MIC13_UP)} MIC13_UP + {len(PART2_MIC13_DOWN)} MIC13_DOWN (Part 2); no MT genes")

# == Load the 95-gene Mic.13 marker list ======================================
mic13_df = pd.read_csv(MIC13_GENES_FILE, sep="\t")
mic13_df = mic13_df[mic13_df.iloc[:, 0] == "Mic.13"].copy()
mic13_df.columns = ["state", "gene", "mean_exp", "pct_exp"]
MIC13_ALL = set(mic13_df["gene"].values)
print(f"\n  Mic.13 marker genes loaded: {len(MIC13_ALL)}")

# Split into in-scorer (circular) and held-out (the real test)
MIC13_IN_SCORER = MIC13_ALL & SCORER_GENES
MIC13_HELD_OUT  = MIC13_ALL - SCORER_GENES
print(f"  In scorer (circular): {len(MIC13_IN_SCORER)}  "
      f"{sorted(MIC13_IN_SCORER)}")
print(f"  Held out (real test): {len(MIC13_HELD_OUT)}")

# == DAM2 signature (Keren-Shaul et al. 2017, Table S3) =======================
# Core DAM2 genes, upregulated in stage 2 disease-associated microglia.
DAM2_GENES = {
    "APOE", "SPP1", "CST7", "LPL", "TREM2", "AXL", "CTSD",
    "ITGAX", "GPNMB", "CLEC7A", "LILRB4", "TYROBP", "CTSB",
    "CTSL", "LIPA", "CD9", "TIMP2", "CCL6", "IGFBP1",
}
DAM2_HELD_OUT = DAM2_GENES - SCORER_GENES
print(f"\n  DAM2 signature genes: {len(DAM2_GENES)}")
print(f"  DAM2 held out (not in scorer): {len(DAM2_HELD_OUT)}  "
      f"{sorted(DAM2_HELD_OUT)}")


# =============================================================================
# SECTION 2: Load Part 5 DEG results
# =============================================================================
print("\n" + "=" * 70)
print("SECTION 2: Load DEG results from Part 5")
print("=" * 70)

xl = pd.ExcelFile(DEG_MASTER_FILE)
print(f"  Sheets: {xl.sheet_names}")

# Collect all PVM_3 SIG sheets (one per region)
deg_results = {}
for sheet in xl.sheet_names:
    if sheet.startswith("PVM_3_") and sheet.endswith("_SIG"):
        region = sheet.replace("PVM_3_", "").replace("_SIG", "")
        df = xl.parse(sheet)
        deg_results[region] = df
        print(f"    {sheet}: {len(df)} sig DEGs")

# Also load the overlap sheet (genes sig in multiple regions)
if "PVM_3_overlap" in xl.sheet_names:
    overlap_df = xl.parse("PVM_3_overlap")
    print(f"    PVM_3_overlap: {len(overlap_df)} genes")
else:
    overlap_df = pd.DataFrame()
    print("    PVM_3_overlap: not found")

xl.close()

# Build union of all significant DEGs across regions
all_sig_genes = set()
all_sig_up    = set()
all_sig_down  = set()
for region, df in deg_results.items():
    genes = set(df["gene"].values)
    all_sig_genes |= genes
    if "direction" in df.columns:
        _dir = df["direction"].astype(str).str.upper()
        all_sig_up   |= set(df[_dir.str.contains("UP")]["gene"].values)
        all_sig_down |= set(df[_dir.str.contains("DOWN")]["gene"].values)

print(f"\n  Union of sig DEGs across {len(deg_results)} regions:")
print(f"    Total unique: {len(all_sig_genes)}")
print(f"    UP in >=1 region: {len(all_sig_up)}")
print(f"    DOWN in >=1 region: {len(all_sig_down)}")

# Genes sig in >= N regions, counted from the per-region _SIG sheets
# (p_adj < 0.05, |log2FC| >= 0.25 and pct_fg >= 0.15).
sig_region_count = {}
for region, df in deg_results.items():
    for g in set(df["gene"].values):
        sig_region_count[g] = sig_region_count.get(g, 0) + 1
for min_n in [2, 3, 4]:
    n = sum(1 for v in sig_region_count.values() if v >= min_n)
    print(f"    Sig in >= {min_n} regions: {n}")
if len(overlap_df) > 0 and "n_regions" in overlap_df.columns:
    _n_ungated = int((overlap_df["n_regions"] >= 4).sum())
    _n_gated = sum(1 for v in sig_region_count.values() if v >= 4)
    print(f"    sig in >= 4 regions: {_n_gated} gated (used), "
          f"{_n_ungated} in PVM_3_overlap without the pct_fg gate (not used)")

# Get total genes tested (from any ALL sheet) for the hypergeometric universe
for sheet in xl.sheet_names:
    if sheet.endswith("_ALL"):
        xl2 = pd.ExcelFile(DEG_MASTER_FILE)
        all_tested = xl2.parse(sheet)
        UNIVERSE_SIZE = len(all_tested)
        xl2.close()
        print(f"\n  Universe (genes tested): {UNIVERSE_SIZE}  (from {sheet})")
        break


# =============================================================================
# SECTION 3: Hypergeometric overlap tests
# =============================================================================
print("\n" + "=" * 70)
print("SECTION 3: Hypergeometric overlap - Mic.13 held-out genes")
print("=" * 70)

def hypergeom_test(gene_set, deg_set, universe_n, label):
    """One-sided hypergeometric test for enrichment.
    
    H0: the overlap between gene_set and deg_set is no greater than expected
        by chance, given universe_n total genes.
    
    Returns dict with counts and p-value.
    """
    K = len(gene_set)                          # genes in the reference set
    n = len(deg_set)                           # genes in the DEG list
    k = len(gene_set & deg_set)                # overlap
    expected = K * n / max(universe_n, 1)
    fold = k / max(expected, 1e-10)
    # P(X >= k) under hypergeometric
    p = hypergeom.sf(k - 1, universe_n, K, n)
    
    overlap_genes = sorted(gene_set & deg_set)
    
    result = {
        "test": label,
        "reference_set_size": K,
        "deg_list_size": n,
        "universe_size": universe_n,
        "overlap_count": k,
        "expected_overlap": round(expected, 2),
        "fold_enrichment": round(fold, 2),
        "p_value": p,
        "overlap_genes": ", ".join(overlap_genes),
    }
    
    sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
    print(f"  [{label}]")
    print(f"    Reference: {K} genes | DEGs: {n} | Universe: {universe_n}")
    print(f"    Overlap: {k} (expected {expected:.1f}, fold {fold:.1f}x)")
    print(f"    p = {p:.4e}  {sig}")
    if overlap_genes:
        print(f"    Genes: {', '.join(overlap_genes[:20])}"
              + (f" ... (+{len(overlap_genes)-20})" if len(overlap_genes) > 20 else ""))
    
    return result


# == TEST 1a: held-out Mic.13 genes vs all sig DEGs ===========================
print("\n--- TEST 1a: Held-out Mic.13 vs all sig DEGs (union across regions) ---")
test_results = []
r = hypergeom_test(MIC13_HELD_OUT, all_sig_genes, UNIVERSE_SIZE,
                   "Mic.13 held-out vs all sig DEGs")
test_results.append(r)

# == TEST 1b: held-out Mic.13 vs UP DEGs only =================================
print("\n--- TEST 1b: Held-out Mic.13 vs UP DEGs only ---")
r = hypergeom_test(MIC13_HELD_OUT, all_sig_up, UNIVERSE_SIZE,
                   "Mic.13 held-out vs UP DEGs")
test_results.append(r)

# == TEST 1c: in-scorer Mic.13 genes (circular) ===============================
print("\n--- TEST 1c: In-scorer Mic.13 vs all sig DEGs (circular) ---")
r = hypergeom_test(MIC13_IN_SCORER, all_sig_genes, UNIVERSE_SIZE,
                   "Mic.13 IN-SCORER vs all sig DEGs (CIRCULAR)")
test_results.append(r)

# == TEST 1d: per-region overlap ==============================================
print("\n--- TEST 1d: Held-out Mic.13 vs sig DEGs per region ---")
for region in sorted(deg_results.keys()):
    df = deg_results[region]
    reg_genes = set(df["gene"].values)
    r = hypergeom_test(MIC13_HELD_OUT, reg_genes, UNIVERSE_SIZE,
                       f"Mic.13 held-out vs {region} sig DEGs")
    test_results.append(r)

# == TEST 1e: genes sig in >= 4 regions =======================================
print(f"\n--- TEST 1e: Held-out Mic.13 vs genes sig in >= 4 regions ---")
robust_genes = {g for g, v in sig_region_count.items() if v >= 4}
r = hypergeom_test(MIC13_HELD_OUT, robust_genes, UNIVERSE_SIZE,
                   "Mic.13 held-out vs robust DEGs (>=4 regions)")
test_results.append(r)


# =============================================================================
# SECTION 4: DAM2 overlap test
# =============================================================================
print("\n" + "=" * 70)
print("SECTION 4: DAM2 signature overlap (Keren-Shaul et al. 2017)")
print("=" * 70)

print("\n--- TEST 2a: DAM2 held-out vs all sig DEGs ---")
r = hypergeom_test(DAM2_HELD_OUT, all_sig_genes, UNIVERSE_SIZE,
                   "DAM2 held-out vs all sig DEGs")
test_results.append(r)

print("\n--- TEST 2b: DAM2 held-out vs UP DEGs ---")
r = hypergeom_test(DAM2_HELD_OUT, all_sig_up, UNIVERSE_SIZE,
                   "DAM2 held-out vs UP DEGs")
test_results.append(r)


# =============================================================================
# SECTION 5: Detailed gene-level table
# =============================================================================
print("\n" + "=" * 70)
print("SECTION 5: Gene-level annotation table")
print("=" * 70)

# For each Mic.13 gene: is it in-scorer or held-out? Is it a sig DEG?
# In how many regions? Direction? Log2FC range?
gene_rows = []
for _, row in mic13_df.iterrows():
    gene = row["gene"]
    in_scorer = gene in SCORER_GENES
    is_sig = gene in all_sig_genes
    is_up  = gene in all_sig_up
    is_down = gene in all_sig_down
    
    # Determine direction: if in both, use majority of log2FC signs
    if is_up and is_down:
        direction = "MIXED"
    elif is_up:
        direction = "UP"
    elif is_down:
        direction = "DOWN"
    elif is_sig:
        direction = "SIG"  # sig but direction column was empty
    else:
        direction = "-"
    
    # Count regions where this gene is significant
    n_regions_sig = 0
    regions_sig = []
    log2fc_values = {}
    for region, df in deg_results.items():
        if gene in df["gene"].values:
            n_regions_sig += 1
            regions_sig.append(region)
            fc = df[df["gene"] == gene]["log2FC"].values
            if len(fc) > 0:
                log2fc_values[region] = fc[0]
    
    gene_rows.append({
        "gene": gene,
        "pct_exp_mic13": round(row["pct_exp"], 3),
        "mean_exp_mic13": round(row["mean_exp"], 3),
        "in_scorer": "YES" if in_scorer else "no",
        "status": "IN-SCORER (circular)" if in_scorer else "HELD-OUT",
        "is_sig_deg": "YES" if is_sig else "no",
        "direction": direction,
        "n_regions_sig": n_regions_sig,
        "regions_sig": ", ".join(regions_sig) if regions_sig else "-",
        "log2fc_range": (f"{min(log2fc_values.values()):.3f} to "
                        f"{max(log2fc_values.values()):.3f}"
                        if log2fc_values else "-"),
    })

gene_table = pd.DataFrame(gene_rows)
n_held_out_sig = ((gene_table["status"] == "HELD-OUT") & 
                  (gene_table["is_sig_deg"] == "YES")).sum()
n_held_out_total = (gene_table["status"] == "HELD-OUT").sum()
n_scorer_sig = ((gene_table["status"] == "IN-SCORER (circular)") & 
                (gene_table["is_sig_deg"] == "YES")).sum()
n_scorer_total = (gene_table["status"] == "IN-SCORER (circular)").sum()

print(f"  Held-out Mic.13 genes: {n_held_out_sig}/{n_held_out_total} are sig DEGs")
print(f"  In-scorer Mic.13 genes: {n_scorer_sig}/{n_scorer_total} are sig DEGs (circular)")

# Print the held-out genes that ARE sig DEGs
held_out_hits = gene_table[(gene_table["status"] == "HELD-OUT") & 
                           (gene_table["is_sig_deg"] == "YES")]
print(f"\n  Held-out Mic.13 genes recovered as sig DEGs ({len(held_out_hits)}):")
for _, row in held_out_hits.iterrows():
    print(f"    {row['gene']:<15} dir={row['direction']}  "
          f"n_regions={row['n_regions_sig']}  "
          f"log2FC={row['log2fc_range']}")


# =============================================================================
# SECTION 6: Excel export
# =============================================================================
print("\n" + "=" * 70)
print("SECTION 6: Excel export")
print("=" * 70)

# == Style helpers ============================================================
_tn  = Side(style="thin", color="BBBBBB")
_bdr = Border(left=_tn, right=_tn, top=_tn, bottom=_tn)
_bf  = Font(name="Calibri", size=9)
_hf  = Font(name="Calibri", bold=True, size=10, color="FFFFFF")
_hfl = PatternFill("solid", fgColor="1F3864")
_af  = PatternFill("solid", fgColor="EFF3FB")
_ctr = Alignment(horizontal="center", vertical="center")
_lft = Alignment(horizontal="left",   vertical="center")
_sig_fill = PatternFill("solid", fgColor="C6EFCE")
_circ_fill = PatternFill("solid", fgColor="FFC7CE")

def _write_header(ws, headers, widths, row=1):
    for ci, (h, w) in enumerate(zip(headers, widths), 1):
        c = ws.cell(row=row, column=ci, value=h)
        c.font = _hf; c.fill = _hfl; c.border = _bdr; c.alignment = _ctr
        ws.column_dimensions[get_column_letter(ci)].width = w
    ws.row_dimensions[row].height = 24

# Arm prefix in the file name.
fp = os.path.join(OUT_DIR, f"Part5b_{ARM_TAG}_Mic13_validation.xlsx")
wb = Workbook(); wb.remove(wb.active)

# Sheet 1: Hypergeometric test results
ws1 = wb.create_sheet("hypergeometric_tests")
_write_header(ws1,
    ["Test", "Reference set", "DEG list", "Universe", "Overlap",
     "Expected", "Fold enrichment", "p-value", "Sig", "Overlap genes"],
    [40, 14, 12, 12, 10, 10, 16, 14, 8, 80])
for ri, r in enumerate(test_results, 2):
    sig = ("***" if r["p_value"] < 0.001 else "**" if r["p_value"] < 0.01 
           else "*" if r["p_value"] < 0.05 else "ns")
    vals = [r["test"], r["reference_set_size"], r["deg_list_size"],
            r["universe_size"], r["overlap_count"], r["expected_overlap"],
            r["fold_enrichment"], r["p_value"], sig, r["overlap_genes"]]
    for ci, v in enumerate(vals, 1):
        c = ws1.cell(ri, ci, v)
        c.font = _bf; c.border = _bdr
        c.alignment = _lft if ci in (1, 10) else _ctr
        if ci == 9 and sig != "ns":
            c.fill = _sig_fill
        if "CIRCULAR" in str(r["test"]):
            if ci == 1: c.fill = _circ_fill
    if ri % 2 == 0:
        for ci in (1, 2, 3): ws1.cell(ri, ci).fill = _af

# Sheet 2: Gene-level table
ws2 = wb.create_sheet("mic13_gene_detail")
_write_header(ws2,
    ["Gene", "pct_exp (Mic.13)", "mean_exp (Mic.13)", "In scorer?",
     "Status", "Sig DEG?", "Direction", "N regions sig", "Regions", "log2FC range"],
    [14, 16, 16, 12, 22, 10, 10, 14, 30, 18])
for ri, (_, row) in enumerate(gene_table.iterrows(), 2):
    vals = [row["gene"], row["pct_exp_mic13"], row["mean_exp_mic13"],
            row["in_scorer"], row["status"], row["is_sig_deg"],
            row["direction"], row["n_regions_sig"], row["regions_sig"],
            row["log2fc_range"]]
    for ci, v in enumerate(vals, 1):
        c = ws2.cell(ri, ci, v)
        c.font = _bf; c.border = _bdr
        c.alignment = _lft if ci in (1, 5, 9, 10) else _ctr
        if row["status"] == "IN-SCORER (circular)" and ci == 5:
            c.fill = _circ_fill
        elif row["is_sig_deg"] == "YES" and ci == 6:
            c.fill = _sig_fill
    if ri % 2 == 0:
        ws2.cell(ri, 1).fill = _af

wb.save(fp)
print(f"  Saved {fp}")


# =============================================================================
# SECTION 7: Final report
# =============================================================================
print("\n" + "=" * 70)
print("SECTION 7: Part 5b Final Report  [NeoC arm]")
print("=" * 70)

# Extract key numbers for the report
mic13_held_test = [r for r in test_results 
                   if r["test"] == "Mic.13 held-out vs all sig DEGs"][0]
mic13_up_test   = [r for r in test_results 
                   if r["test"] == "Mic.13 held-out vs UP DEGs"][0]
dam2_test       = [r for r in test_results 
                   if r["test"] == "DAM2 held-out vs all sig DEGs"][0]

print(f"""
Part 5b COMPLETE - Mic.13 Equivalence Validation  [{ARM_TAG} arm]
{'='*70}
DEG source      : {DEG_MASTER_FILE}
Mic.13 list     : {len(MIC13_ALL)} genes (Green et al. 2024, Extended Data Fig 2c)
Scorer genes    : {len(SCORER_GENES)} ({len(PART2_MIC13_UP)} MIC13_UP + {len(PART2_MIC13_DOWN)} MIC13_DOWN, Part 2)

GENE SPLIT:
  In-scorer (circular) : {len(MIC13_IN_SCORER)} genes
  Held-out (real test) : {len(MIC13_HELD_OUT)} genes

TEST 1 - Mic.13 held-out genes among Part 5 sig DEGs:
  All DEGs : {mic13_held_test['overlap_count']}/{len(MIC13_HELD_OUT)} recovered, \
fold={mic13_held_test['fold_enrichment']}x, p={mic13_held_test['p_value']:.4e}
  UP DEGs  : {mic13_up_test['overlap_count']}/{len(MIC13_HELD_OUT)} recovered, \
fold={mic13_up_test['fold_enrichment']}x, p={mic13_up_test['p_value']:.4e}

TEST 2 - DAM2 signature (Keren-Shaul et al. 2017):
  All DEGs : {dam2_test['overlap_count']}/{len(DAM2_HELD_OUT)} recovered, \
fold={dam2_test['fold_enrichment']}x, p={dam2_test['p_value']:.4e}
{'='*70}
""")

print(f"\nAll outputs in {OUT_DIR}/:")
for fn in sorted(os.listdir(OUT_DIR)):
    if not fn.startswith("."):
        print(f"  {fn}")
print(f"\nPart 5b ({ARM_NAME}) complete.")
