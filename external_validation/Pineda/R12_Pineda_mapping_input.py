# R12_Pineda_mapping_input.py
# Prepares the kept nuclei of external validation 1 (Pineda et al. 2024, BA9 batch 210526) for
# MapMyCells. This is the Pineda counterpart of R12_Rexach_mapping_input.py.
# No correlation or other statistic is computed here.
#
# Steps:
#   1. Reads the R11 nuclei table and keeps the nuclei that R02 kept (called by emptyDrops, at least
#      500 genes, and scDblFinder singlet) and that are in the authors' published QC (barcode present
#      in the authors' <LIB>_col_metadata.tsv; Pineda et al. 2024 STAR Methods: >= 700 genes, <= 20%
#      mtRNA, and per-cluster removal of low-quality, doublet, high-ambient, and damaged nuclei), as
#      for Rexach. The QC is applied here, before MapMyCells, so that only these nuclei are mapped.
#   2. Reads their raw counts from counts/<LIB>_raw_feature_bc_matrix.h5 (Cell Ranger 10.1.0), and
#      checks the gene IDs against R11 and the total UMIs and genes detected of every nucleus against
#      R02.
#   3. MapMyCells input: the kept nuclei as gzip h5ad files, split by donor into parts below the 2 GB
#      upload limit of the MapMyCells web app. X holds the raw integer counts, genes are Ensembl IDs
#      (var index), and cells are the R02 cell IDs (<LIB>_<barcode>, obs index), both written as plain
#      strings that the MapMyCells server can read.
#
# Usage (WSL, the same venv as the Rexach scripts):
#   cd /mnt/c/Users/Wong/Desktop/SEA-AD/External_Dataset_1
#   ~/venv_synapse/bin/python R12_Pineda_mapping_input.py
# Inputs, in the script's folder: R11_nuclei/R11_Pineda_nuclei.csv.gz, R11_nuclei/R11_Pineda_genes.csv.gz,
#   counts/<LIB>_raw_feature_bc_matrix.h5, and the authors' metadata Pineda2024_metadata/<LIB>_col_metadata.tsv
#   (searched in subfolders too).
# Outputs, in the folder R12_mapmycells_input:
#   R12_Pineda_MapMyCells_part<k>.h5ad   MapMyCells input, one file per part
#   R12_Pineda_nuclei.csv.gz             one row per kept nucleus, with its QC and part
#   R12_Pineda_summary.xlsx              sheets Donors, Parts, Settings, and Session
#   R12_Pineda_log.txt                   console output, with progress every 30 s

import glob
import importlib.metadata as im
import os
import platform
import sys
import time

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

try:
    import openpyxl
    from openpyxl.styles import Font
except ImportError:
    sys.exit("Install openpyxl first: ~/venv_synapse/bin/pip install -q openpyxl")

SET = {
    "part_cells": 50000,         # nuclei per MapMyCells part (whole donors are kept together)
    "max_bytes": 1_950_000_000,  # below the 2 GB MapMyCells upload limit
    "progress_s": 30,
    "author_qc_only": True,      # map only the nuclei in the authors' published QC
}
MT13 = ["MT-ND1", "MT-ND2", "MT-CO1", "MT-CO2", "MT-ATP8", "MT-ATP6", "MT-CO3", "MT-ND3",
        "MT-ND4L", "MT-ND4", "MT-ND5", "MT-ND6", "MT-CYB"]

# The data folder: PINEDA_ROOT when set (the container sets it to /data), otherwise the script folder.
HERE = os.environ.get("PINEDA_ROOT", os.path.dirname(os.path.abspath(__file__)))
F_NUC = os.path.join(HERE, "R11_nuclei", "R11_Pineda_nuclei.csv.gz")
F_GENES = os.path.join(HERE, "R11_nuclei", "R11_Pineda_genes.csv.gz")
COUNTS = os.path.join(HERE, "counts")
AUTHOR_DIR = os.path.join(HERE, "Pineda2024_metadata")
OUT = os.path.join(HERE, "R12_mapmycells_input")
os.makedirs(OUT, exist_ok=True)
log = open(os.path.join(OUT, "R12_Pineda_log.txt"), "w", encoding="utf-8")
pd.set_option("display.width", 200)
pd.set_option("display.max_rows", 200)
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


def raw_h5(lib):
    return os.path.join(COUNTS, lib + "_raw_feature_bc_matrix.h5")


def read_library(lib, barcodes, gene_ids):
    """Raw counts of the given barcodes of one Cell Ranger raw .h5, as a CSR matrix (nuclei x genes)
    in the order of `barcodes`. Cell Ranger stores the matrix by barcode (CSC with genes as rows), so
    the columns of one barcode are one contiguous slice."""
    with h5py.File(raw_h5(lib), "r") as f:
        g = f["matrix"]
        ids = np.asarray([x.decode() for x in g["features"]["id"][:]], dtype=object)
        if not np.array_equal(ids, gene_ids):
            stop("%s: the gene IDs differ from R11" % lib)
        shape = tuple(int(x) for x in g["shape"][:])
        bc = pd.Index([x.decode() for x in g["barcodes"][:]])
        pos = bc.get_indexer(barcodes)
        if (pos < 0).any():
            stop("%s: %d kept barcodes not found in the raw .h5" % (lib, int((pos < 0).sum())))
        indptr = g["indptr"][:]
        order = np.argsort(pos)
        sp_pos = pos[order]
        breaks = np.where(np.diff(sp_pos) != 1)[0] + 1
        dparts, iparts, lens = [], [], []
        for run in np.split(sp_pos, breaks):
            a, b = indptr[run[0]], indptr[run[-1] + 1]
            dparts.append(g["data"][a:b])
            iparts.append(g["indices"][a:b])
            lens.append(np.diff(indptr[run[0]:run[-1] + 2]))
    data = np.concatenate(dparts)
    if data.size and not np.all(np.equal(np.mod(data, 1), 0)):
        stop("%s: non-integer values in the raw .h5" % lib)
    m = sp.csr_matrix((data.astype(np.int32), np.concatenate(iparts),
                       np.concatenate([[0], np.cumsum(np.concatenate(lens))])),
                      shape=(len(sp_pos), shape[0]))
    inv = np.empty_like(order)
    inv[order] = np.arange(len(order))
    m = m[inv]
    m.eliminate_zeros()
    m.sort_indices()
    return m


t0 = time.time()
started = time.strftime("%Y-%m-%d %H:%M:%S")
say("R12_Pineda_mapping_input.py started", started)
say("Folder:", HERE)

# ---- 1. Kept nuclei -----------------------------------------------------------------------------
for p in (F_NUC, F_GENES):
    if not os.path.exists(p):
        stop("R11 output not found: " + p)
nuc = pd.read_csv(F_NUC, dtype={"donor_id": str, "barcode": str, "cell_id": str})
genes = pd.read_csv(F_GENES)
gene_ids = np.asarray(genes["gene_id"].astype(str).tolist(), dtype=object)
symbols = genes["symbol"].astype(str).values
if len(set(gene_ids)) != len(gene_ids):
    stop("the gene IDs are not unique")
i_mt = np.where(np.isin(symbols, MT13))[0]
if len(i_mt) != 13:
    stop("expected 13 MT genes, found %d" % len(i_mt))
say("\n== 1. Nuclei and genes")
say("R11 called nuclei: %d | kept by R02: %d | genes: %d | gene IDs starting with ENSG: %d" %
    (len(nuc), int(nuc["keep"].sum()), len(gene_ids), sum(s.startswith("ENSG") for s in gene_ids)))
kept = nuc[nuc["keep"]].copy()
libs = sorted(kept["library_id"].unique())
missing = [l for l in libs if not os.path.exists(raw_h5(l))]
if missing:
    stop("%d raw .h5 files missing in counts/, for example %s" % (len(missing), missing[0]))
per_donor_libs = kept.groupby("donor_id")["library_id"].nunique()
if (per_donor_libs != 1).any() or kept.groupby("library_id")["donor_id"].nunique().ne(1).any():
    stop("one library per donor and one donor per library is not true")
say("Libraries: %d | donors: %d | one library per donor: True" % (len(libs), per_donor_libs.shape[0]))

# ---- Authors' published QC -----------------------------------------------------------------------
qc_rows = []
if SET["author_qc_only"]:
    meta = {}
    for p in glob.glob(os.path.join(AUTHOR_DIR, "**", "*_col_metadata.tsv"), recursive=True):
        meta.setdefault(os.path.basename(p)[:-len("_col_metadata.tsv")], p)
    no_meta = [l for l in libs if l not in meta]
    if no_meta:
        stop("the authors' metadata are missing for %d libraries (for example %s) in %s" %
             (len(no_meta), no_meta[0], AUTHOR_DIR))
    kept["bc16"] = kept["barcode"].str.split("-").str[0]
    inq = pd.Series(False, index=kept.index)
    for lib in libs:
        m = pd.read_csv(meta[lib], sep="\t", dtype=str, usecols=["Barcode"])
        b16 = m["Barcode"].str.split("-").str[0]
        if b16.duplicated().any():
            stop("%s: duplicated 16-base barcodes in the authors' metadata" % lib)
        k = kept["library_id"] == lib
        inq[k] = kept.loc[k, "bc16"].isin(set(b16)).values
        qc_rows.append({"library_id": lib, "author_nuclei": len(m), "our_kept": int(k.sum()),
                        "our_kept_in_author_qc": int(inq[k].sum()),
                        "author_nuclei_not_in_our_kept": len(m) - int(inq[k].sum())})
    n_before = len(kept)
    kept = kept[inq].drop(columns=["bc16"]).copy()
    qc_tab = pd.DataFrame(qc_rows)
    say("Authors' published QC: %d of %d kept nuclei are in the authors' QC set (%d libraries) | "
        "authors' nuclei not among our kept nuclei: %d of %d" %
        (len(kept), n_before, len(libs), int(qc_tab["author_nuclei_not_in_our_kept"].sum()),
         int(qc_tab["author_nuclei"].sum())))
    if kept.groupby("donor_id").size().reindex(per_donor_libs.index).isna().any():
        stop("a donor has no nuclei left after the authors' QC")

# ---- Parts: whole donors, in donor order, up to part_cells nuclei each ------------------------------
donor_n = kept.groupby("donor_id").size().sort_index()
part_of, part, n_in = {}, 1, 0
for d, n in donor_n.items():
    if n_in > 0 and n_in + n > SET["part_cells"]:
        part, n_in = part + 1, 0
    part_of[d] = part
    n_in += n
kept["part"] = kept["donor_id"].map(part_of)
n_parts = part
say("MapMyCells parts:", n_parts, "| up to", SET["part_cells"], "nuclei each, whole donors kept together")

# ---- 2 and 3. Read counts per part, check against R02, and write the MapMyCells input ----------------
say("\n== 2. Counts, checks against R02, and MapMyCells files")
part_rows, last = [], time.time()
for p in range(1, n_parts + 1):
    sub = kept[kept["part"] == p].sort_values(["library_id", "cell_id"])
    blocks, ids = [], []
    for lib, s in sub.groupby("library_id", sort=True):
        m = read_library(lib, s["barcode"].tolist(), gene_ids)
        tot = np.asarray(m.sum(axis=1)).ravel()
        ng = np.diff(m.indptr)
        mt = np.asarray(m[:, i_mt].sum(axis=1)).ravel()
        if not (np.array_equal(tot, s["total_umis"].values) and np.array_equal(ng, s["n_genes"].values)
                and np.array_equal(mt, s["mt_umis"].values)):
            stop("%s: total UMIs, genes detected, or MT UMIs differ from R02" % lib)
        blocks.append(m)
        ids.extend(s["cell_id"].tolist())
        if time.time() - last >= SET["progress_s"]:
            say("  part %d: library %s read | %.0f s" % (p, lib, time.time() - t0))
            last = time.time()
    X = sp.vstack(blocks, format="csr")
    # The cell and gene IDs are written as plain (object) strings with no other obs or var columns.
    # With pandas 3, anndata otherwise writes them as "nullable-string-array", which the older
    # anndata on the MapMyCells server cannot read (found with the Rexach upload).
    adata = ad.AnnData(X=X, obs=pd.DataFrame(index=pd.Index(list(ids), dtype=object)),
                       var=pd.DataFrame(index=pd.Index(list(gene_ids), dtype=object)))
    name = "R12_Pineda_MapMyCells_part%d.h5ad" % p
    tmp = os.path.join(OUT, "tmp_" + name)
    say("  part %d: writing %s (%d nuclei, %d donors) | %.0f s" %
        (p, name, X.shape[0], sub["donor_id"].nunique(), time.time() - t0))
    adata.write_h5ad(tmp, compression="gzip")
    with h5py.File(tmp, "r") as chk:
        enc = {k: chk[k + "/" + chk[k].attrs["_index"]].attrs.get("encoding-type") for k in ("obs", "var")}
    if any(v != "string-array" for v in enc.values()):
        stop("part %d has index encodings %s; MapMyCells needs string-array" % (p, enc))
    os.replace(tmp, os.path.join(OUT, name))
    size = os.path.getsize(os.path.join(OUT, name))
    say("  part %d: %.2f GB" % (p, size / 1e9))
    if size > SET["max_bytes"]:
        say("  WARNING: part %d is above the upload limit. Lower part_cells in SET and run again." % p)
    part_rows.append({"part": p, "file": name, "donors": ";".join(sorted(sub["donor_id"].unique())),
                      "n_donors": sub["donor_id"].nunique(), "nuclei": X.shape[0],
                      "bytes": size, "below_limit": size <= SET["max_bytes"]})
    del adata, X, blocks

kept.drop(columns=["keep"]).to_csv(os.path.join(OUT, "R12_Pineda_nuclei.csv.gz"), index=False,
                                   compression="gzip")

# ---- Summaries --------------------------------------------------------------------------------------
say("\n== 3. Summary")
say("Kept nuclei: %d | median UMIs %.0f | median genes %.0f | median %%MT %.2f | max %%MT %.2f" %
    (len(kept), kept["total_umis"].median(), kept["n_genes"].median(), kept["pct_mt"].median(),
     kept["pct_mt"].max()))
dn = kept.groupby("donor_id").agg(disease=("disease", "first"), sex=("sex", "first"),
                                  library_id=("library_id", "first"), part=("part", "first"),
                                  nuclei_kept=("cell_id", "size"), median_umis=("total_umis", "median"),
                                  median_genes=("n_genes", "median"),
                                  median_pct_mt=("pct_mt", "median")).reset_index()
dn = dn.sort_values(["disease", "donor_id"])
say("\nDonors:")
say(dn.to_string(index=False))
pr = pd.DataFrame(part_rows)
settings = pd.DataFrame({"setting": list(SET.keys()) + ["counts", "MT genes", "kept nuclei"],
                         "value": [str(v) for v in SET.values()] +
                                  ["counts/<LIB>_raw_feature_bc_matrix.h5", ", ".join(MT13),
                                   "R02 keep: called (emptyDrops FDR <= 0.01), at least 500 genes, scDblFinder singlet"
                                   + ("; and in the authors' published QC" if SET["author_qc_only"] else "")]})
pk = ["anndata", "h5py", "numpy", "pandas", "scipy", "openpyxl"]
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
    tmp = path + ".tmp"
    wb.save(tmp)
    os.replace(tmp, path)


write_xlsx(os.path.join(OUT, "R12_Pineda_summary.xlsx"),
           [("Donors", dn), ("Parts", pr)] + ([("Author_QC", pd.DataFrame(qc_rows))] if qc_rows else []) +
           [("Settings", settings), ("Session", session)])
say("\nUpload each R12_Pineda_MapMyCells_part<k>.h5ad to MapMyCells (10x Human MTG SEA-AD,"
    " CCN20230505), once with Deep Generative Mapping and once with Hierarchical Mapping.")
say("Put the result csv files in R13_mapmycells_output/.")
say("Done in %.0f s" % (time.time() - t0))
log.close()
