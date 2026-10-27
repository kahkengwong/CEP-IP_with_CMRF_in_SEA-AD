# R13_Pineda_count_checks.py
# Count checks of external validation 1 (Pineda et al. 2024, BA9 batch 210526), run after
# MapMyCells and before R14. Counts only: no correlation, screen, GAM, IP, or association is computed
# here. This is the Pineda counterpart of R13_Rexach_count_checks.py, with the same checks and
# outputs, plus the doublet shares of section C and the expected r_s of section E.
#
# Checks:
#   A. MapMyCells files: all present (each R12 part, Deep Generative and SEA-AD Hierarchical
#      Mapping), each covering exactly the nuclei of its R12 part, with no duplicates or extras.
#      The algorithm and codebase version are taken from each csv header.
#   B. Taxonomy names: Micro-PVM_3-SEAAD and Micro-PVM_1 written exactly as in the
#      SEA-AD release.
#   C. Agreement: Deep Generative against Hierarchical Mapping at class, subclass,
#      and supertype, within Microglia-PVM, and against the authors' microglia (CellType "Micro" in
#      the authors' per-library metadata); scDblFinder doublets per library.
#   D. Micro-PVM_3 nuclei per donor, by diagnosis, with donors below 5, 10, 20, and 50 nuclei. The
#      arm is Micro-PVM_3 by Deep Generative Mapping.
#   E. Depth of the PVM_3 arm (Micro-PVM_3 with MT-CO3 > 0): median UMIs, genes, MT
#      UMIs, and MT-CO3 UMIs per nucleus, the %MT range, the HKG8 arm count (PVM_3 arm nuclei with
#      ACTB > 0), and the r_s expected at this depth from SEA-AD Stats_15H.
#   F. Export for the R steps: the raw counts of every nucleus that is Micro-PVM_3 in either mapping
#      (all genes, genes x nuclei), with a nucleus table that flags each mapping and whether the
#      nucleus is in the authors' published QC.
#
# Usage (WSL, the same venv as R12):
#   cd /mnt/c/Users/Wong/Desktop/SEA-AD/External_Dataset_1
#   ~/venv_synapse/bin/python R13_Pineda_count_checks.py
# Reads, from the script's folder: R11_nuclei/ (nuclei, genes, and donors), R12_mapmycells_input/
# R12_Pineda_nuclei.csv.gz, the MapMyCells csv files in R13_mapmycells_output/, counts/<LIB>_raw_
# feature_bc_matrix.h5, Stats_15H.csv (SEA-AD, Figure 15H source data), and the authors' metadata
# <LIB>_col_metadata.tsv in AUTHOR_DIR (searched in subfolders too).
# Outputs, in the folder R13_count_checks:
#   R13_Pineda_count_checks.xlsx   sheets Files, Taxonomy_names, Agreement, Microglia_DGM_vs_HM,
#                                  Author_vs_DGM, Doublets, Donors, Arm_depth, Expected_rs,
#                                  Settings, and Session
#   R13_Pineda_MicroPVM3_counts.mtx.gz   raw counts, genes x nuclei (MatrixMarket)
#   R13_Pineda_MicroPVM3_genes.csv.gz    gene_id and symbol, in matrix row order
#   R13_Pineda_MicroPVM3_nuclei.csv.gz   one row per nucleus, in matrix column order
#   R13_Pineda_log.txt                   console output, with progress every 30 s

import glob
import gzip
import importlib.metadata as im
import os
import platform
import re
import sys
import time

import h5py
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

try:
    import openpyxl
    from openpyxl.styles import Font
except ImportError:
    sys.exit("Install openpyxl first: ~/venv_synapse/bin/pip install -q openpyxl")

# The data folder: PINEDA_ROOT when set (the container sets it to /data), otherwise the script folder.
HERE = os.environ.get("PINEDA_ROOT", os.path.dirname(os.path.abspath(__file__)))
# Folder with the authors' per-library metadata (<LIB>_col_metadata.tsv, Synapse). Searched with
# subfolders. If no file is found, the comparisons with the authors' labels and QC are left out.
AUTHOR_DIR = os.path.join(HERE, "Pineda2024_metadata")
MAP_DIR = os.path.join(HERE, "R13_mapmycells_output")
NUC_CSV = os.path.join(HERE, "R12_mapmycells_input", "R12_Pineda_nuclei.csv.gz")
R11_NUC = os.path.join(HERE, "R11_nuclei", "R11_Pineda_nuclei.csv.gz")
R11_GENES = os.path.join(HERE, "R11_nuclei", "R11_Pineda_genes.csv.gz")
STATS_15H = os.path.join(HERE, "Stats_15H.csv")
COUNTS = os.path.join(HERE, "counts")
OUT = os.path.join(HERE, "R13_count_checks")
os.makedirs(OUT, exist_ok=True)
SET = {"arm_supertype": "Micro-PVM_3-SEAAD", "pvm1_supertype": "Micro-PVM_1",
       "microglia_subclass": "Microglia-PVM", "author_microglia": "Micro",
       "goi": "MT-CO3", "hkg8_goi": "ACTB", "donor_floor": 5,
       "mt_report_cuts": [1, 5, 10], "progress_s": 30,
       # Analysis set: only nuclei kept by the authors' published QC (Pineda et al. 2024, STAR Methods:
       # >= 700 genes, <= 20% mtRNA, and per-cluster removal of low-quality, doublet, high-ambient,
       # and damaged nuclei), as for Rexach. Sections D to F use this set; A to C describe all nuclei.
       "author_qc_only": True}
MT13 = ["MT-ND1", "MT-ND2", "MT-CO1", "MT-CO2", "MT-ATP8", "MT-ATP6", "MT-CO3", "MT-ND3",
        "MT-ND4L", "MT-ND4", "MT-ND5", "MT-ND6", "MT-CYB"]
ALGS = {"DeepGenerativeMapping": "DGM", "SEA-ADHierarchicalMapping": "HM"}

log = open(os.path.join(OUT, "R13_Pineda_log.txt"), "w", encoding="utf-8")
pd.set_option("display.width", 220)
pd.set_option("display.max_rows", 300)
pd.set_option("display.max_columns", 30)


def say(*parts):
    s = " ".join(str(p) for p in parts)
    print(s, flush=True)
    log.write(s + "\n")
    log.flush()


def stop(msg):
    say("STOPPED:", msg)
    log.close()
    sys.exit(1)


t0 = time.time()
started = time.strftime("%Y-%m-%d %H:%M:%S")
say("R13_Pineda_count_checks.py started", started)

for p in (NUC_CSV, R11_NUC, R11_GENES):
    if not os.path.exists(p):
        stop("input not found: " + p)
kept = pd.read_csv(NUC_CSV, dtype={"donor_id": str, "barcode": str, "cell_id": str})
say("R12 kept nuclei: %d | donors %d" % (len(kept), kept["donor_id"].nunique()))

# ---- A. MapMyCells files ------------------------------------------------------------------
say("\n== A. MapMyCells files")
pat = re.compile(r"part(\d+)_.*_(DeepGenerativeMapping|SEA-ADHierarchicalMapping)_UTC_\d+\.csv$")
files = {}
for p in sorted(glob.glob(os.path.join(MAP_DIR, "*.csv"))):
    m = pat.search(os.path.basename(p))
    if m:
        key = (int(m.group(1)), ALGS[m.group(2)])
        if key in files:
            stop("two csv files for part %d %s" % key)
        files[key] = p
parts = sorted(kept["part"].unique())
missing = [(p, a) for p in parts for a in ("DGM", "HM") if (p, a) not in files]
if missing:
    stop("missing MapMyCells csv files: " + ", ".join("part%d %s" % k for k in missing))


def read_mapping(path):
    with open(path, encoding="utf-8") as fh:
        header = [ln.strip() for ln in fh if ln.startswith("#")]
    d = pd.read_csv(path, comment="#")
    alg = next((h for h in header if h.startswith("# algorithm")), "")
    return d, alg


maps = {"DGM": [], "HM": []}
file_rows = []
for (p, a), path in sorted(files.items()):
    d, alg = read_mapping(path)
    exp = set(kept.loc[kept["part"] == p, "cell_id"])
    got = d["cell_id"].astype(str)
    row = {"part": p, "algorithm": a, "file": os.path.basename(path), "header": alg,
           "nuclei": len(d), "expected": len(exp), "duplicates": int(got.duplicated().sum()),
           "missing": len(exp - set(got)), "extra": len(set(got) - exp)}
    row["ok"] = row["nuclei"] == row["expected"] and row["duplicates"] == 0 \
        and row["missing"] == 0 and row["extra"] == 0
    file_rows.append(row)
    say("part %d %s: %d nuclei (expected %d) | duplicates %d | missing %d | extra %d | %s" %
        (p, a, row["nuclei"], row["expected"], row["duplicates"], row["missing"], row["extra"],
         "OK" if row["ok"] else "NOT OK"))
    say("   ", alg)
    maps[a].append(d)
fr = pd.DataFrame(file_rows)
if not fr["ok"].all():
    stop("a MapMyCells file does not match its R12 part (see above)")


def levels(d):
    """class, subclass, supertype name and probability columns of a MapMyCells csv."""
    out = {}
    for lv in ("class", "subclass", "supertype"):
        name = lv + "_name"
        probs = [c for c in d.columns if c.startswith(lv + "_") and c.endswith("probability")]
        if name not in d.columns:
            stop("column %s not found; columns are %s" % (name, list(d.columns)))
        out[lv] = (name, probs[0] if probs else None)
    return out


mp = {}
for a in ("DGM", "HM"):
    d = pd.concat(maps[a], ignore_index=True)
    lv = levels(d)
    cols = {"cell_id": d["cell_id"].astype(str)}
    for k, (name, prob) in lv.items():
        cols["%s_%s" % (a, k)] = d[name].astype(str)
        if prob:
            cols["%s_%s_prob" % (a, k)] = pd.to_numeric(d[prob], errors="coerce")
    mp[a] = pd.DataFrame(cols)
    say("%s probability columns used: %s" % (a, ", ".join(p for _, p in lv.values() if p)))
cells = kept.merge(mp["DGM"], on="cell_id", how="left").merge(mp["HM"], on="cell_id", how="left")
if cells["DGM_supertype"].isna().any() or cells["HM_supertype"].isna().any():
    stop("some kept nuclei have no mapping")

# ---- Authors' labels and QC (per-library metadata) -----------------------------------------
say("\n== Authors' metadata")
meta_files = {}
for p in glob.glob(os.path.join(AUTHOR_DIR, "**", "*_col_metadata.tsv"), recursive=True):
    lib = os.path.basename(p)[:-len("_col_metadata.tsv")]
    meta_files.setdefault(lib, p)
libs = sorted(cells["library_id"].unique())
have_meta = [l for l in libs if l in meta_files]
say("Authors' metadata files found for %d of %d libraries (folder %s)" % (len(have_meta), len(libs), AUTHOR_DIR))
# Object columns, so that strings and True/False can be set (NaN stays for libraries without metadata).
cells["author_celltype"] = pd.Series(np.nan, index=cells.index, dtype=object)
cells["in_author_qc"] = pd.Series(np.nan, index=cells.index, dtype=object)
if have_meta:
    cells["bc16"] = cells["barcode"].str.split("-").str[0]
    am_rows = []
    for lib in have_meta:
        m = pd.read_csv(meta_files[lib], sep="\t", dtype=str)
        if "Barcode" not in m.columns or "CellType" not in m.columns:
            stop("%s: columns Barcode and CellType are needed; found %s" % (os.path.basename(meta_files[lib]),
                                                                           list(m.columns)))
        m["bc16"] = m["Barcode"].str.split("-").str[0]
        if m["bc16"].duplicated().any():
            stop("%s: duplicated 16-base barcodes in the authors' metadata" % lib)
        k = cells["library_id"] == lib
        lab = cells.loc[k, "bc16"].map(m.set_index("bc16")["CellType"])
        cells.loc[k, "author_celltype"] = lab.values.astype(object)
        cells.loc[k, "in_author_qc"] = lab.notna().values.astype(object)
        am_rows.append({"library_id": lib, "author_nuclei": len(m), "author_microglia":
                        int((m["CellType"] == SET["author_microglia"]).sum()),
                        "our_kept": int(k.sum()), "our_kept_in_author_set": int(lab.notna().sum())})
    am_tab = pd.DataFrame(am_rows)
    say("Kept nuclei in the authors' set: %d of %d (libraries with metadata)" %
        (int(am_tab["our_kept_in_author_set"].sum()), int(am_tab["our_kept"].sum())))
    say("Authors' nuclei found among our kept nuclei: %d of %d" %
        (int(am_tab["our_kept_in_author_set"].sum()), int(am_tab["author_nuclei"].sum())))
    cells = cells.drop(columns=["bc16"])
else:
    am_tab = pd.DataFrame({"note": ["no authors' metadata found; author comparisons left out"]})

# ---- B. Taxonomy names -----------------------------------------------------------------------
say("\n== B. Taxonomy names")
tn = []
for a in ("DGM", "HM"):
    st = set(cells[a + "_supertype"])
    for nm in (SET["arm_supertype"], SET["pvm1_supertype"]):
        tn.append({"algorithm": a, "name": nm, "present": nm in st})
    say("%s Microglia-PVM supertypes: %s" % (a, ", ".join(sorted(
        cells.loc[cells[a + "_subclass"] == SET["microglia_subclass"], a + "_supertype"].unique()))))
tn = pd.DataFrame(tn)
say(tn.to_string(index=False))

# ---- C. Agreement ----------------------------------------------------------------------------
say("\n== C. Agreement")
agr = []
for lv in ("class", "subclass", "supertype"):
    same = (cells["DGM_" + lv] == cells["HM_" + lv])
    agr.append({"comparison": "DGM vs HM, all kept nuclei", "level": lv, "nuclei": len(cells),
                "agree": int(same.sum()), "pct_agree": round(100 * same.mean(), 2)})
mgd = cells["DGM_subclass"] == SET["microglia_subclass"]
same = cells.loc[mgd, "DGM_supertype"] == cells.loc[mgd, "HM_supertype"]
agr.append({"comparison": "DGM vs HM, DGM Microglia-PVM nuclei", "level": "supertype",
            "nuclei": int(mgd.sum()), "agree": int(same.sum()), "pct_agree": round(100 * same.mean(), 2)})
arm_d = cells["DGM_supertype"] == SET["arm_supertype"]
arm_h = cells["HM_supertype"] == SET["arm_supertype"]
both = arm_d & arm_h
agr.append({"comparison": "Micro-PVM_3: DGM nuclei also Micro-PVM_3 in HM", "level": "supertype",
            "nuclei": int(arm_d.sum()), "agree": int(both.sum()),
            "pct_agree": round(100 * both.sum() / max(int(arm_d.sum()), 1), 2)})
agr.append({"comparison": "Micro-PVM_3: HM nuclei also Micro-PVM_3 in DGM", "level": "supertype",
            "nuclei": int(arm_h.sum()), "agree": int(both.sum()),
            "pct_agree": round(100 * both.sum() / max(int(arm_h.sum()), 1), 2)})
if have_meta:
    inq = cells["in_author_qc"] == True  # noqa: E712 (NaN for libraries without metadata)
    am = inq & (cells["author_celltype"] == SET["author_microglia"])
    for a in ("DGM", "HM"):
        x = (cells[a + "_subclass"] == SET["microglia_subclass"]) & inq
        agr.append({"comparison": "%s Microglia-PVM nuclei (in the authors' set) that the authors call Micro" % a,
                    "level": "subclass", "nuclei": int(x.sum()), "agree": int((x & am).sum()),
                    "pct_agree": round(100 * (x & am).sum() / max(int(x.sum()), 1), 2)})
        agr.append({"comparison": "Authors' Micro nuclei that %s calls Microglia-PVM" % a,
                    "level": "subclass", "nuclei": int(am.sum()), "agree": int((x & am).sum()),
                    "pct_agree": round(100 * (x & am).sum() / max(int(am.sum()), 1), 2)})
        y = (cells[a + "_supertype"] == SET["arm_supertype"]) & inq
        agr.append({"comparison": "%s Micro-PVM_3 nuclei (in the authors' set) that the authors call Micro" % a,
                    "level": "supertype", "nuclei": int(y.sum()), "agree": int((y & am).sum()),
                    "pct_agree": round(100 * (y & am).sum() / max(int(y.sum()), 1), 2)})
    for a in ("DGM", "HM"):
        y = cells[a + "_supertype"] == SET["arm_supertype"]
        agr.append({"comparison": "%s Micro-PVM_3 nuclei in the authors' published QC" % a, "level": "QC",
                    "nuclei": int((y & cells["in_author_qc"].notna()).sum()), "agree": int((y & inq).sum()),
                    "pct_agree": round(100 * (y & inq).sum() / max(int((y & cells["in_author_qc"].notna()).sum()), 1), 2)})
agr = pd.DataFrame(agr)
say(agr.to_string(index=False))
if "DGM_supertype_prob" in cells.columns:
    pr = cells.loc[arm_d, "DGM_supertype_prob"]
    say("DGM Micro-PVM_3 supertype probability: median %.3f | below 0.5: %d of %d" %
        (pr.median(), int((pr < 0.5).sum()), len(pr)))
for a in ("DGM", "HM"):
    sc = cells[a + "_subclass"] == "Sst Chodl"
    say("%s Sst Chodl nuclei: %d | authors' labels for them: %s" % (a, int(sc.sum()), "; ".join(
        "%s %d" % kv for kv in cells.loc[sc, "author_celltype"].value_counts().head(4).items())))
sc_d = cells["DGM_subclass"] == "Sst Chodl"
say("DGM Sst Chodl nuclei, HM subclass: %s" % "; ".join(
    "%s %d" % kv for kv in cells.loc[sc_d, "HM_subclass"].value_counts().head(5).items()))
if sc_d.any():
    say("DGM Sst Chodl nuclei: median genes %.0f against %.0f for all kept nuclei" %
        (cells.loc[sc_d, "n_genes"].median(), cells["n_genes"].median()))

xt_mg = pd.crosstab(cells.loc[mgd | (cells["HM_subclass"] == SET["microglia_subclass"]), "DGM_supertype"],
                    cells.loc[mgd | (cells["HM_subclass"] == SET["microglia_subclass"]), "HM_supertype"])
xt_mg = xt_mg.loc[:, xt_mg.sum() > 0]
xt_auth = pd.crosstab(cells["author_celltype"].fillna("not in the authors' set"), cells["DGM_subclass"])
say("\nMicroglia-PVM supertypes, DGM (rows) against HM (columns):")
say(xt_mg.to_string())

# Doublets per library, from the R02 nuclei table: share of scDblFinder doublets
# among the called nuclei that pass the 500-gene floor.
r11 = pd.read_csv(R11_NUC, dtype={"donor_id": str, "barcode": str, "cell_id": str},
                  usecols=["library_id", "donor_id", "disease", "pass_gene_floor", "scDblFinder_class", "keep"])
fl = r11[r11["pass_gene_floor"]]
dbl = fl.groupby("library_id").agg(donor_id=("donor_id", "first"), disease=("disease", "first"),
                                   nuclei_floor=("keep", "size"),
                                   doublets=("scDblFinder_class", lambda s: int((s == "doublet").sum())),
                                   kept=("keep", "sum")).reset_index()
dbl["pct_doublets"] = (100 * dbl["doublets"] / dbl["nuclei_floor"]).round(2)
dbl = dbl.sort_values(["disease", "donor_id"])
say("\nscDblFinder doublets per library (share of nuclei above the gene floor): median %.2f%% | range %.2f%% to %.2f%%" %
    (dbl["pct_doublets"].median(), dbl["pct_doublets"].min(), dbl["pct_doublets"].max()))
say(dbl.groupby("disease")["pct_doublets"].agg(["count", "median", "min", "max"]).round(2).to_string())

# ---- Analysis set: the authors' published QC ------------------------------------------------------
qc = cells["in_author_qc"] == True  # noqa: E712 (NaN for libraries without metadata)
if SET["author_qc_only"]:
    if len(have_meta) != len(libs):
        stop("the authors' metadata are missing for %d libraries, so the authors' QC cannot be applied" %
             (len(libs) - len(have_meta)))
    n_pvm3_before = int(arm_d.sum())
    mgd, arm_d, arm_h = mgd & qc, arm_d & qc, arm_h & qc
    both = arm_d & arm_h
    say("\n== Analysis set: nuclei in the authors' published QC")
    say("Kept nuclei in the authors' QC: %d of %d | Micro-PVM_3 (DGM): %d of %d" %
        (int(qc.sum()), len(cells), int(arm_d.sum()), n_pvm3_before))

# ---- D. Micro-PVM_3 per donor ------------------------------------------------------------------
say("\n== D. Micro-PVM_3 nuclei per donor (arm: Deep Generative Mapping%s)" %
    (", authors' QC" if SET["author_qc_only"] else ""))
dn = cells.groupby("donor_id").agg(disease=("disease", "first"), sex=("sex", "first"),
                                   library_id=("library_id", "first"),
                                   nuclei_kept=("cell_id", "size")).reset_index()
dn["nuclei_author_qc"] = dn["donor_id"].map(cells[qc].groupby("donor_id").size()).fillna(0).astype(int)
masks = [("microglia_pvm_dgm", mgd), ("micro_pvm3_dgm", arm_d), ("micro_pvm3_hm", arm_h),
         ("micro_pvm3_both", both), ("micro_pvm1_dgm", (cells["DGM_supertype"] == SET["pvm1_supertype"]) & (qc | (not SET["author_qc_only"])))]
if have_meta:
    masks.append(("author_microglia", am))
for col, mask in masks:
    dn[col] = dn["donor_id"].map(cells[mask].groupby("donor_id").size()).fillna(0).astype(int)

# ---- E. Read counts for the union of Micro-PVM_3 nuclei (both mappings) ---------------------------
genes = pd.read_csv(R11_GENES)
gene_ids = np.asarray(genes["gene_id"].astype(str).tolist(), dtype=object)
symbols = np.asarray(genes["symbol"].astype(str).tolist(), dtype=object)
sel = cells[arm_d | arm_h].copy().sort_values(["library_id", "cell_id"]).reset_index(drop=True)
say("\nReading counts of %d nuclei that are Micro-PVM_3 in either mapping, from %d raw .h5 files" %
    (len(sel), sel["library_id"].nunique()))
blocks, last = [], time.time()
for lib, s in sel.groupby("library_id", sort=True):
    f5 = os.path.join(COUNTS, lib + "_raw_feature_bc_matrix.h5")
    if not os.path.exists(f5):
        stop("raw .h5 not found: " + f5)
    with h5py.File(f5, "r") as f:
        g = f["matrix"]
        ids = np.asarray([x.decode() for x in g["features"]["id"][:]], dtype=object)
        if not np.array_equal(ids, gene_ids):
            stop("%s: the gene IDs differ from R11" % lib)
        ng_total = int(g["shape"][0])
        bc = pd.Index([x.decode() for x in g["barcodes"][:]])
        pos = bc.get_indexer(s["barcode"].tolist())
        if (pos < 0).any():
            stop("%s: barcodes not found in the raw .h5" % lib)
        indptr = g["indptr"][:]
        order = np.argsort(pos)
        sp_pos = pos[order]
        breaks = np.where(np.diff(sp_pos) != 1)[0] + 1
        dparts, iparts, lens = [], [], []
        for run in np.split(sp_pos, breaks):
            a0, b0 = indptr[run[0]], indptr[run[-1] + 1]
            dparts.append(g["data"][a0:b0]); iparts.append(g["indices"][a0:b0])
            lens.append(np.diff(indptr[run[0]:run[-1] + 2]))
    m = sp.csr_matrix((np.concatenate(dparts).astype(np.int32), np.concatenate(iparts),
                       np.concatenate([[0], np.cumsum(np.concatenate(lens))])), shape=(len(sp_pos), ng_total))
    inv = np.empty_like(order); inv[order] = np.arange(len(order))
    m = m[inv]
    m.eliminate_zeros()
    m.sort_indices()
    blocks.append(m)
    if time.time() - last >= SET["progress_s"]:
        say("  library %s read | %.0f s" % (lib, time.time() - t0))
        last = time.time()
M = sp.vstack(blocks, format="csr")
gi = {s: i for i, s in enumerate(symbols)}
sel["total_umis_check"] = np.asarray(M.sum(axis=1)).ravel()
if not np.array_equal(sel["total_umis_check"].values, sel["total_umis"].values):
    stop("total UMIs read here differ from R02")
sel["mt_co3_umis"] = M[:, gi[SET["goi"]]].toarray().ravel()
sel["actb_umis"] = M[:, gi[SET["hkg8_goi"]]].toarray().ravel()
sel["is_pvm3_dgm"] = sel["DGM_supertype"] == SET["arm_supertype"]
sel["is_pvm3_hm"] = sel["HM_supertype"] == SET["arm_supertype"]
sel["in_pvm3_arm"] = sel["is_pvm3_dgm"] & (sel["mt_co3_umis"] > 0)
sel["in_hkg8_arm"] = sel["in_pvm3_arm"] & (sel["actb_umis"] > 0)

arm = sel[sel["in_pvm3_arm"]]
dn["pvm3_arm"] = dn["donor_id"].map(arm.groupby("donor_id").size()).fillna(0).astype(int)
dn["hkg8_arm"] = dn["donor_id"].map(sel[sel["in_hkg8_arm"]].groupby("donor_id").size()).fillna(0).astype(int)
dn["pct_pvm3_mtco3_pos"] = (100 * dn["pvm3_arm"] / dn["micro_pvm3_dgm"].replace(0, np.nan)).round(1)
dn["at_donor_floor"] = dn["micro_pvm3_dgm"] >= SET["donor_floor"]
dn = dn.sort_values(["disease", "donor_id"])
say(dn.to_string(index=False))
say("\nMicro-PVM_3 (DGM) per donor by diagnosis:")
say(dn.groupby("disease")["micro_pvm3_dgm"].agg(["count", "sum", "median", "min", "max"]).to_string())
for flr in (5, 10, 20, 50):
    below = dn.loc[dn["micro_pvm3_dgm"] < flr, "donor_id"].tolist()
    say("Donors with fewer than %d Micro-PVM_3 nuclei: %d %s" % (flr, len(below), ", ".join(below)))
el = dn[dn["at_donor_floor"]]
say("Donors at the floor (>= %d Micro-PVM_3 nuclei): %d | by diagnosis: %s" % (
    SET["donor_floor"], len(el), "; ".join("%s %d" % kv for kv in el["disease"].value_counts().items())))

say("\n== E. Depth of the PVM_3 arm (Micro-PVM_3 by DGM with MT-CO3 > 0)")
n_pvm3 = int(sel["is_pvm3_dgm"].sum())
say("Micro-PVM_3 (DGM): %d | with MT-CO3 > 0 (PVM_3 arm): %d (%.1f%%) | HKG8 arm (also ACTB > 0): %d" %
    (n_pvm3, len(arm), 100 * len(arm) / max(n_pvm3, 1), int(sel["in_hkg8_arm"].sum())))
qc = ("total_umis", "n_genes", "mt_umis", "mt_co3_umis", "pct_mt")
dep = pd.DataFrame({
    "measure": ["total UMIs", "genes detected", "MT UMIs", "MT-CO3 UMIs", "%MT"],
    "median": [arm[c].median() for c in qc],
    "p10": [arm[c].quantile(0.1) for c in qc],
    "p90": [arm[c].quantile(0.9) for c in qc],
    "max": [arm[c].max() for c in qc]})
dep[["median", "p10", "p90", "max"]] = dep[["median", "p10", "p90", "max"]].round(3)
say(dep.to_string(index=False))
say("PVM_3 arm nuclei above %s%% MT: %s (%s)" % (
    ", ".join(str(c) for c in SET["mt_report_cuts"]),
    ", ".join("%d" % int((arm["pct_mt"] > c).sum()) for c in SET["mt_report_cuts"]),
    "within the authors' QC, which removed nuclei above 20% MT" if SET["author_qc_only"] else "no %MT cut"))
say("MT-CO3 UMIs in the PVM_3 arm, share of nuclei with 1, 2, 3, and 4 or more: %s" % ", ".join(
    "%.1f%%" % (100 * v) for v in [(arm["mt_co3_umis"] == 1).mean(), (arm["mt_co3_umis"] == 2).mean(),
                                   (arm["mt_co3_umis"] == 3).mean(), (arm["mt_co3_umis"] >= 4).mean()]))

# Expected r_s at this depth: the SEA-AD PVM_3 depth ladder of Figure 15H
# (Stats_15H, level "depth", median r_s of MT-CO3 over its partners). The row used is the largest
# target depth at or below the median total UMIs of the Pineda PVM_3 arm. It sets no threshold.
exp_rs = pd.DataFrame({"note": ["Stats_15H.csv not found in %s; expected r_s not computed" % HERE]})
if os.path.exists(STATS_15H):
    s15 = pd.read_csv(STATS_15H)
    lad = s15[(s15["arm"] == "PVM_3") & (s15["level"] == "depth") &
              (s15["quantity"] == "rs_median_over_partners")].copy()
    lad["target"] = lad["measure"].str.extract(r"target (\d+)")[0].astype(float)
    med = float(arm["total_umis"].median())
    below = lad[lad["target"] <= med]
    row = below.loc[below["target"].idxmax()] if len(below) else lad.loc[lad["target"].idxmin()]
    exp_rs = pd.DataFrame([{"pineda_arm_median_umis": med, "seaad_target_depth": row["target"],
                            "expected_rs": row["estimate"], "expected_rs_lo": row["lo"],
                            "expected_rs_hi": row["hi"], "seaad_observed_median_depth":
                            float(s15.loc[(s15["arm"] == "PVM_3") & (s15["level"] == "observed") &
                                          (s15["quantity"] == "median_depth"), "estimate"].iloc[0]),
                            "rule": "largest SEA-AD PVM_3 target depth at or below the Pineda arm median UMIs"}])
    say("Expected r_s (SEA-AD Stats_15H): arm median %.0f UMIs -> target %d: r_s %.3f (95%% CI %.3f to %.3f)" %
        (med, int(row["target"]), row["estimate"], row["lo"], row["hi"]))
else:
    say("Stats_15H.csv not found; copy it to %s for the expected r_s" % HERE)

# ---- F. Export --------------------------------------------------------------------------------------
say("\n== F. Export of the Micro-PVM_3 nuclei (either mapping)")
with gzip.open(os.path.join(OUT, "tmp_R13_Pineda_MicroPVM3_counts.mtx.gz"), "wb") as fh:
    scipy.io.mmwrite(fh, M.T.tocoo(), field="integer")
os.replace(os.path.join(OUT, "tmp_R13_Pineda_MicroPVM3_counts.mtx.gz"),
           os.path.join(OUT, "R13_Pineda_MicroPVM3_counts.mtx.gz"))
pd.DataFrame({"gene_id": gene_ids, "symbol": symbols}).to_csv(
    os.path.join(OUT, "R13_Pineda_MicroPVM3_genes.csv.gz"), index=False, compression="gzip")
keep_cols = ["cell_id", "donor_id", "library_id", "disease", "sex", "author_celltype", "in_author_qc",
             "total_umis", "n_genes", "mt_umis", "pct_mt", "mt_co3_umis", "actb_umis",
             "DGM_subclass", "DGM_supertype", "DGM_supertype_prob", "HM_subclass", "HM_supertype",
             "HM_supertype_prob", "is_pvm3_dgm", "is_pvm3_hm", "in_pvm3_arm", "in_hkg8_arm"]
sel[[c for c in keep_cols if c in sel.columns]].to_csv(
    os.path.join(OUT, "R13_Pineda_MicroPVM3_nuclei.csv.gz"), index=False, compression="gzip")
say("Matrix: %d genes x %d nuclei | genes with duplicate symbols: %d" %
    (M.shape[1], M.shape[0], int(pd.Series(symbols).duplicated().sum())))

settings = pd.DataFrame({"setting": list(SET.keys()) + ["counts", "authors' metadata"],
                         "value": [str(v) for v in SET.values()] +
                                  ["counts/<LIB>_raw_feature_bc_matrix.h5",
                                   "%d libraries, %s" % (len(have_meta), AUTHOR_DIR)]})
pk = ["h5py", "numpy", "pandas", "scipy", "openpyxl"]
session = pd.DataFrame({"item": ["python", "platform", "started", "finished"] + pk,
                        "value": [sys.version.split()[0], platform.platform(), started,
                                  time.strftime("%Y-%m-%d %H:%M:%S")] + [im.version(x) for x in pk]})


def write_xlsx(path, sheets):
    """Plain workbook: bold header and a frozen top row, no fills or borders."""
    wb = openpyxl.Workbook()
    wb.remove(wb.active)
    for name, df in sheets:
        ws = wb.create_sheet(name)
        ws.append(list(df.columns))
        for c in ws[1]:
            c.font = Font(bold=True)
        for r in df.itertuples(index=False):
            ws.append([None if (isinstance(v, float) and np.isnan(v)) else
                       (v.item() if hasattr(v, "item") else v) for v in r])
        ws.freeze_panes = "A2"
    wb.save(path + ".tmp")
    os.replace(path + ".tmp", path)


write_xlsx(os.path.join(OUT, "R13_Pineda_count_checks.xlsx"), [
    ("Files", fr), ("Taxonomy_names", tn), ("Agreement", agr),
    ("Microglia_DGM_vs_HM", xt_mg.reset_index()), ("Author_vs_DGM", xt_auth.reset_index()),
    ("Author_metadata", am_tab), ("Doublets", dbl), ("Donors", dn), ("Arm_depth", dep),
    ("Expected_rs", exp_rs), ("Settings", settings), ("Session", session)])
say("\nDone in %.0f s" % (time.time() - t0))
log.close()
