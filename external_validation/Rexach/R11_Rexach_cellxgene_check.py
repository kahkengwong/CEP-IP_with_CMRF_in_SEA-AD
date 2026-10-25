# R11_Rexach_cellxgene_check.py
# Checks the CELLxGENE copy of Rexach et al. 2024 (collection c53573b2-eff4-4c5e-9ad0-b24d422dfd9b)
# before it is used as external dataset 2. Counts and structure only: no correlation or other
# statistic of the readings is computed here, because the Rexach readings are not written yet.
#
# Checks:
#   1. File structure: where the raw counts are (raw/X or X), their shape, and that they are integers.
#   2. obs columns, and the tissue, disease, assay, and cell-type labels with their counts.
#   3. Insula: donors per diagnosis, assay, and sex.
#   4. The 13 MT genes, MT-CO3, and its eight partners in the gene table.
#   5. Insula microglia (CELLxGENE cell_type "microglial cell"): nuclei per donor, and per nucleus
#      total UMIs, genes detected, MT UMIs, and %MT (to confirm the authors' 8% MT cut).
#
# Usage (Spyder, or a terminal with Python 3 and anndata, h5py, numpy, pandas):
#   python R11_Rexach_cellxgene_check.py <path to the .h5ad>
# Without an argument, it uses the only .h5ad file in the script's folder.
# Outputs, in the folder of the .h5ad:
#   R11_Rexach_check_log.txt
#   R11_Rexach_microglia_by_donor.csv   one row per insula donor

import glob
import importlib.metadata as im
import os
import sys
import time

import anndata as ad
import h5py
import numpy as np
import pandas as pd

MT13 = ["MT-ND1", "MT-ND2", "MT-CO1", "MT-CO2", "MT-ATP8", "MT-ATP6", "MT-CO3", "MT-ND3",
        "MT-ND4L", "MT-ND4", "MT-ND5", "MT-ND6", "MT-CYB"]
DFG8 = ["MT-ATP6", "MT-CO1", "MT-CO2", "MT-CYB", "MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4"]
MICROGLIA = "microglial cell"

if len(sys.argv) > 1:
    H5AD = sys.argv[1]
else:
    here = os.path.dirname(os.path.abspath(__file__))
    found = glob.glob(os.path.join(here, "*.h5ad"))
    if len(found) != 1:
        sys.exit("Give the .h5ad path as the argument (found %d .h5ad files here)." % len(found))
    H5AD = found[0]
OUT = os.path.dirname(os.path.abspath(H5AD))
log = open(os.path.join(OUT, "R11_Rexach_check_log.txt"), "w", encoding="utf-8")
pd.set_option("display.width", 200)
pd.set_option("display.max_rows", 200)
pd.set_option("display.max_columns", 30)


def say(*parts):
    s = " ".join(str(p) for p in parts)
    print(s)
    log.write(s + "\n")
    log.flush()


def table(series, title):
    say("\n" + title)
    say(series.value_counts(dropna=False).to_string())


t0 = time.time()
say("R11_Rexach_cellxgene_check.py started", time.strftime("%Y-%m-%d %H:%M:%S"))
say("File:", H5AD, "|", round(os.path.getsize(H5AD) / 1e9, 2), "GB")
say("anndata", im.version("anndata"), "| h5py", im.version("h5py"))

# ---- 1. Structure ----------------------------------------------------------------
f = h5py.File(H5AD, "r")


def matrix_info(path):
    if path not in f:
        return None
    g = f[path]
    enc = g.attrs.get("encoding-type", "dense") if isinstance(g, h5py.Group) else "dense"
    shape = tuple(g.attrs["shape"]) if isinstance(g, h5py.Group) else g.shape
    dtype = g["data"].dtype if isinstance(g, h5py.Group) else g.dtype
    return {"path": path, "encoding": str(enc), "shape": shape, "dtype": str(dtype)}


say("\n== 1. Structure")
say("Top-level keys:", list(f.keys()))
info_x, info_raw = matrix_info("X"), matrix_info("raw/X")
say("X:", info_x)
say("raw/X:", info_raw)

A = ad.read_h5ad(H5AD, backed="r")  # obs and var in memory; matrices stay on disk
obs = A.obs
say("Cells (obs):", obs.shape[0], "| genes in X (var):", A.var.shape[0])


def read_rows(path, rows):
    """Reads the given rows of a CSR matrix in the file, as a scipy CSR matrix."""
    import scipy.sparse as sp
    g = f[path]
    if g.attrs.get("encoding-type") != "csr_matrix":
        raise ValueError("%s is %s; this script reads csr_matrix only" %
                         (path, g.attrs.get("encoding-type")))
    indptr = g["indptr"][:]
    ncol = int(g.attrs["shape"][1])
    rows = np.sort(np.asarray(rows))
    data_parts, ind_parts, lens = [], [], []
    # Contiguous runs of rows are read in one slice each.
    breaks = np.where(np.diff(rows) != 1)[0] + 1
    for run in np.split(rows, breaks):
        a, b = indptr[run[0]], indptr[run[-1] + 1]
        data_parts.append(g["data"][a:b])
        ind_parts.append(g["indices"][a:b])
        lens.append(np.diff(indptr[run[0]:run[-1] + 2]))
    lens = np.concatenate(lens)
    new_indptr = np.concatenate([[0], np.cumsum(lens)])
    return sp.csr_matrix((np.concatenate(data_parts), np.concatenate(ind_parts), new_indptr),
                         shape=(len(rows), ncol)), rows


# Where are the raw counts? raw/X if present, else X. Checked on 2,000 random cells.
rng = np.random.default_rng(42)
probe = rng.choice(obs.shape[0], size=min(2000, obs.shape[0]), replace=False)
counts_path, var_counts = None, None
for path, var in (("raw/X", A.raw.var if A.raw is not None else None), ("X", A.var)):
    if path not in f or var is None:
        continue
    m, _ = read_rows(path, probe)
    vals = m.data
    is_int = bool(np.all(np.equal(np.mod(vals, 1), 0))) if vals.size else False
    say("%s on 2,000 random cells: integer values %s | min %s | max %s | median total per cell %s"
        % (path, is_int, vals.min() if vals.size else "NA", vals.max() if vals.size else "NA",
           np.median(np.asarray(m.sum(axis=1)).ravel())))
    if is_int and counts_path is None:
        counts_path, var_counts = path, var
say("Raw counts used from:", counts_path)
if counts_path is None:
    say("STOPPED: no integer count matrix found.")
    sys.exit(1)

# ---- 2. obs labels ----------------------------------------------------------------
say("\n== 2. obs columns")
say(", ".join(obs.columns))
for col in ("tissue", "disease", "assay", "suspension_type", "cell_type"):
    if col in obs.columns:
        table(obs[col], "Counts by %s (all cells):" % col)

# ---- 3. Insula --------------------------------------------------------------------
say("\n== 3. Insula")
tissue = obs["tissue"].astype(str)
ins = tissue.str.contains("insula", case=False)
say("Insula cells:", int(ins.sum()), "| tissue label(s):", sorted(tissue[ins].unique()))
oi = obs[ins]
donor_col = "donor_id" if "donor_id" in obs.columns else None
say("Donor column:", donor_col)
join = lambda s: ";".join(sorted(s.astype(str).unique()))
dd = oi.groupby(donor_col, observed=True).agg(disease=("disease", join), assay=("assay", join),
                                              nuclei=("disease", "size"))
if "sex" in oi.columns:
    dd["sex"] = oi.groupby(donor_col, observed=True)["sex"].agg(join)
say("Insula donors:", dd.shape[0])
table(dd["disease"], "Insula donors by disease:")
table(dd["assay"], "Insula donors by assay:")
if "sex" in dd.columns:
    table(dd["sex"], "Insula donors by sex:")
table(oi["cell_type"], "Insula cells by cell_type:")
for col in obs.columns:
    if col.lower() in ("celltype", "cell_type_author", "clusters", "cluster", "subcluster",
                       "major_celltype", "broad_celltype", "subclass", "class"):
        table(oi[col], "Insula cells by %s (author label):" % col)

# ---- 4. Genes ---------------------------------------------------------------------
say("\n== 4. Genes in the count matrix")
sym_col = "feature_name" if "feature_name" in var_counts.columns else None
symbols = (var_counts[sym_col] if sym_col else var_counts.index.to_series()).astype(str).values
say("Gene symbol column:", sym_col or "(index)", "| genes:", len(symbols))
for g in MT13:
    say("  %-8s %s" % (g, "present" if g in symbols else "MISSING"))
say("MT-CO3 and its eight partners all present:", all(g in symbols for g in ["MT-CO3"] + DFG8))

# ---- 5. Insula microglia ------------------------------------------------------------
say("\n== 5. Insula microglia (cell_type '%s')" % MICROGLIA)
mg = ins & (obs["cell_type"].astype(str) == MICROGLIA)
rows_mg = np.where(mg.values)[0]
say("Insula microglia:", len(rows_mg))
if len(rows_mg):
    m, rows_sorted = read_rows(counts_path, rows_mg)
    i_mt = np.where(np.isin(symbols, MT13))[0]
    tot = np.asarray(m.sum(axis=1)).ravel()
    mt = np.asarray(m[:, i_mt].sum(axis=1)).ravel()
    ngenes = np.diff(m.indptr)
    pct = 100 * mt / np.maximum(tot, 1)
    cells = pd.DataFrame({"donor": obs[donor_col].astype(str).values[rows_sorted],
                          "disease": obs["disease"].astype(str).values[rows_sorted],
                          "total_umis": tot, "n_genes": ngenes, "mt_umis": mt, "pct_mt": pct})
    q = np.percentile(pct, [50, 90, 99, 100])
    say("Per-nucleus %%MT: median %.2f | 90th %.2f | 99th %.2f | max %.2f" % tuple(q))
    say("Nuclei with %%MT above 8: %d" % int((pct > 8).sum()))
    say("Per-nucleus total UMIs: median %.0f | genes detected: median %.0f | MT UMIs: median %.0f"
        % (np.median(tot), np.median(ngenes), np.median(mt)))
    by = cells.groupby("donor").agg(disease=("disease", "first"), microglia=("total_umis", "size"),
                                    median_umis=("total_umis", "median"),
                                    median_genes=("n_genes", "median"),
                                    median_mt_umis=("mt_umis", "median"),
                                    median_pct_mt=("pct_mt", "median"),
                                    max_pct_mt=("pct_mt", "max")).reset_index()
    by = by.sort_values(["disease", "donor"])
    say("\nInsula microglia by donor:")
    say(by.round(2).to_string(index=False))
    say("\nMicroglia per donor by disease (median, min, max):")
    say(by.groupby("disease")["microglia"].agg(["count", "median", "min", "max"]).to_string())
    by.to_csv(os.path.join(OUT, "R11_Rexach_microglia_by_donor.csv"), index=False)

f.close()
say("\nDone in %.0f s" % (time.time() - t0))
log.close()
