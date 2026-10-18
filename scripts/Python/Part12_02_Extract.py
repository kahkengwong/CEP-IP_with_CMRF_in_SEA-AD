# =============================================================================
# Part12_02_Extract.py
# SEA-AD CEP-IP/CMRF, Part 12 rerun, step 02, extraction from the h5ad
#   Reads the h5ad and Part12_01_config.json, keeps the cells of the two
#   supertypes the arms are drawn from (Micro-PVM_3-SEAAD and Micro-PVM_1) and
#   writes per-cell quantities, raw counts for the genes later R steps need,
#   a per-supertype gene summary and the matched gene sets. No labels, arm or
#   dataset assignment, statistics or figures; those are done in R.
#
# Usage (set SEAAD_ROOT, or run from the SEA-AD folder)
#   python3 Part12_rerun/Part12_02_Extract.py
#
# Outputs (Part12_rerun/02_Extract/): 02_cells.csv.gz, 02_counts.mtx.gz,
#   02_counts_genes.tsv, 02_counts_cells.tsv, 02_gene_summary.csv,
#   02_matched_sets.csv and Part12_02_Extract_log.txt, where step 04 reads them.
# =============================================================================

import gzip
import io
import json
import os
import platform
import re
import sys
import time

import numpy as np
import pandas as pd
import scipy
import scipy.io
import scipy.sparse as sp

STEP = "02"
T0 = time.time()

# == Root and config ==========================================================
ROOT = os.environ.get("SEAAD_ROOT", os.getcwd())
CFG_PATH = os.path.join(ROOT, "Part12_rerun", "Part12_01_config.json")
if not os.path.isfile(CFG_PATH):
    sys.exit("aborted: config not found: %s" % CFG_PATH)
with open(CFG_PATH, encoding="utf-8") as fh:
    CFG = json.load(fh)

RERUN = os.path.join(ROOT, CFG["rerun_dir"])
OUT = os.path.join(RERUN, "02_Extract")
TAB = OUT
os.makedirs(OUT, exist_ok=True)
LOG = os.path.join(OUT, "Part12_02_Extract_log.txt")
LEDGER = os.path.join(ROOT, CFG["ledger_file"])
SCRIPT = os.path.abspath(__file__)

F_CELLS = os.path.join(TAB, "02_cells.csv.gz")
F_MTX = os.path.join(TAB, "02_counts.mtx.gz")
F_GENES = os.path.join(TAB, "02_counts_genes.tsv")
F_BCS = os.path.join(TAB, "02_counts_cells.tsv")
F_GSUM = os.path.join(TAB, "02_gene_summary.csv")
F_SETS = os.path.join(TAB, "02_matched_sets.csv")
OUTPUTS = [F_CELLS, F_MTX, F_GENES, F_BCS, F_GSUM, F_SETS]

# Remove stale outputs and rejected copies first.
for f in OUTPUTS:
    for p in (f, re.sub(r"(\.[a-z]+)(\.gz)?$", r"_REJECTED\1\2", f)):
        if os.path.exists(p):
            os.remove(p)

# == Log and checks ===========================================================
_log = open(LOG, "w", encoding="ascii", errors="replace")
N = {"PASS": 0, "FAIL": 0, "WARN": 0}
_ID = [0]
_DONE = [False]


def say(s=""):
    print(s)
    if not _log.closed:
        _log.write(s + "\n")
        _log.flush()


def chk(ok, desc, detail="", warn_only=False):
    _ID[0] += 1
    tag = "PASS" if ok else ("WARN" if warn_only else "FAIL")
    N[tag] += 1
    say("[%s] %s-%02d %s | %s" % (tag, STEP, _ID[0], desc, detail))
    return bool(ok)


def note(desc, detail=""):
    say("[NOTE] %s | %s" % (desc, detail))


def rel(p):
    return os.path.relpath(p, ROOT).replace("\\", "/")


def finish():
    if _DONE[0]:
        return
    _DONE[0] = True
    verdict = "PASS" if N["FAIL"] == 0 else "FAIL"
    dur = (time.time() - T0) / 60.0
    say("")
    if verdict == "PASS":
        say("Outputs:")
        for f in OUTPUTS:
            say("  %s | %.1f MB"
                % (rel(f), os.path.getsize(f) / 2**20))
    else:
        for f in OUTPUTS:
            if os.path.exists(f):
                os.replace(f, re.sub(r"(\.[a-z]+)(\.gz)?$",
                                     r"_REJECTED\1\2", f))
    say("Checks: %d PASS, %d FAIL, %d WARN | duration %.1f min"
        % (N["PASS"], N["FAIL"], N["WARN"], dur))
    say("VERDICT: %s" % verdict)
    _log.close()
    with open(LEDGER, "a", encoding="ascii") as fh:
        fh.write("%s\t%s\tPart12_02_Extract.py\t%s\t%s\tPASS=%d FAIL=%d WARN=%d"
                 "\t%.1f min\n" % (time.strftime("%Y-%m-%d %H:%M"), STEP,
                                   RUN_ID, verdict, N["PASS"], N["FAIL"],
                                   N["WARN"], dur))
    if verdict == "FAIL":
        sys.exit("aborted: see %s" % LOG)


def abort(msg):
    chk(False, "run stopped", msg)
    finish()


def _on_error(etype, value, tb):
    """On an unexpected error, close the log, rename outputs, write the ledger."""
    import traceback
    if etype is SystemExit:
        return sys.__excepthook__(etype, value, tb)
    say("")
    say("Unexpected error:")
    for line in traceback.format_exception(etype, value, tb):
        say(line.rstrip())
    try:
        abort("%s: %s" % (etype.__name__, value))
    except SystemExit:
        pass


sys.excepthook = _on_error


def num(v):
    """Numeric vector from an obs column of any storage type."""
    return pd.to_numeric(pd.Series(v).astype(str), errors="coerce").values


RUN_ID = "%s-%s" % (time.strftime("%Y%m%d-%H%M"),
                    "windows" if os.name == "nt" else "fedora")

# == Header ===================================================================
import anndata as ad  # noqa: E402
from importlib.metadata import version as _pkgver  # noqa: E402

say("Part 12 rerun -- step 02, extraction from the h5ad")
say("Script: Part12_02_Extract.py")
say("Run ID: %s" % RUN_ID)
say("Start: %s" % time.strftime("%Y-%m-%d %H:%M:%S"))
say("Machine: %s | Python %s" % (platform.node(), platform.python_version()))
say("Packages: numpy %s | pandas %s | scipy %s | anndata %s"
    % (np.__version__, pd.__version__, scipy.__version__, _pkgver("anndata")))
say("Config: %s | version %s"
    % (rel(CFG_PATH), CFG["config_version"]))
# Config version check, before the h5ad is read.
NEED_VERSION = 3
NEED_KEYS = ["rerun_dir", "ledger_file", "inputs", "arms",
             "datasets", "obs_columns", "library_fields", "genes", "panels",
             "nuclear_panels", "inference"]
_v_ok = str(CFG.get("config_version", "")).isdigit() and \
    int(CFG["config_version"]) >= NEED_VERSION
chk(_v_ok, "config is version %d or later" % NEED_VERSION,
    "found version %s" % CFG.get("config_version"))
_miss = [k for k in NEED_KEYS if k not in CFG]
chk(not _miss, "config has every section this step reads",
    "missing: %s" % (", ".join(_miss) or "none"))
if N["FAIL"]:
    abort("replace Part12_01_config.json with the current version")
H5AD = os.path.join(ROOT, CFG["inputs"]["h5ad"]["path"])
if not chk(os.path.isfile(H5AD), "h5ad present", rel(H5AD)):
    abort("h5ad missing")
say("Input: %s | %.1f GB"
    % (rel(H5AD), os.path.getsize(H5AD) / 2**30))
SEED = int(CFG["inference"]["master_seed"]) + \
    int(CFG["inference"]["step_codes"][STEP])
say("Seed: %d (master %d + step code %d)"
    % (SEED, CFG["inference"]["master_seed"],
       CFG["inference"]["step_codes"][STEP]))
say("")

# == Open the h5ad and resolve obs columns by exact name ======================
say("SECTION A -- obs columns and cells")
A = ad.read_h5ad(H5AD, backed="r")
obs = A.obs
low = {str(c).lower(): c for c in obs.columns}
OC = {}
for key, name in CFG["obs_columns"].items():
    if key.startswith("_") or key == "raw_layer":
        continue
    OC[key] = low.get(name.lower())
    chk(OC[key] is not None, "obs column '%s' present" % name,
        str(OC[key]) if OC[key] is not None else "not found")
LAYER = CFG["obs_columns"]["raw_layer"]
chk(LAYER in A.layers.keys(), "raw count layer '%s' present" % LAYER,
    ", ".join(str(k) for k in A.layers.keys() if k is not None))
if N["FAIL"]:
    abort("the h5ad does not carry a required column or layer")

bcs_all = pd.Index(A.obs_names.astype(str))
if not chk(bcs_all.is_unique, "barcodes are unique in the h5ad",
           "%d barcodes" % len(bcs_all)):
    abort("duplicate barcodes in the h5ad")

SUPERTYPES = sorted({a["supertype"] for a in CFG["arms"].values()})
st_all = obs[OC["supertype"]].astype(str).values
for s in SUPERTYPES:
    n = int((st_all == s).sum())
    chk(n > 0, "supertype '%s' present" % s, "%d cells" % n)
if N["FAIL"]:
    abort("a configured supertype is not in the h5ad")
keep = np.where(np.isin(st_all, SUPERTYPES))[0]
note("cells kept (the two supertypes)", "%d of %d" % (keep.size, len(bcs_all)))

# == Raw counts of the kept cells =============================================
say("")
say("SECTION B -- raw counts")
sub = A[keep]
mat = sub.layers[LAYER]
mat = mat.to_memory() if hasattr(mat, "to_memory") else mat
mat = sp.csr_matrix(mat)
genes = pd.Index(A.var_names.astype(str))
chk(genes.is_unique, "gene names are unique", "%d genes" % len(genes))
_int_ok = bool(np.all(mat.data >= 0) and np.all(mat.data == np.rint(mat.data)))
chk(_int_ok, "raw layer holds non-negative integers",
    "%d non-zero entries | stored as %s" % (mat.nnz, mat.data.dtype))
if not _int_ok:
    abort("the raw layer is not integer counts")
mat.data = np.rint(mat.data).astype(np.int32)
tot = np.asarray(mat.sum(axis=1, dtype=np.int64)).ravel().astype(np.float64)
# Column-compressed copy for fast per-gene column selection.
matc = mat.tocsc()
ngen = np.asarray((mat > 0).sum(axis=1)).ravel()
O = sub.obs
if not chk(np.array_equal(np.asarray(O.index.astype(str)),
                          np.asarray(bcs_all[keep])),
           "kept cells are in the same order as the barcode list",
           "%d cells" % len(O)):
    abort("the subset changed the cell order")
obs_umi = num(O[OC["n_umis"]])
d_umi = np.nanmax(np.abs(obs_umi - tot))
chk(d_umi < 0.5, "layer totals equal the obs UMI count",
    "max |difference| %.3g" % d_umi)


def gidx(names):
    return [genes.get_loc(g) for g in names if g in genes]


def rowsum(names):
    """Per-cell sum over the genes present; NaN when none is present."""
    idx = gidx(names)
    if not idx:
        return np.full(mat.shape[0], np.nan)
    return np.asarray(matc[:, idx].sum(axis=1)).ravel().astype(np.float64)


def missing(names):
    return [g for g in names if g not in genes]


# == %MT on the 13 named genes ================================================
say("")
say("SECTION C -- per-cell quantities")
MT13 = CFG["genes"]["mt_coding"]
chk(not missing(MT13), "all 13 MT protein-coding genes present",
    "missing: %s" % (", ".join(missing(MT13)) or "none"))
mt_prefix = [g for g in genes if g.startswith("MT-")]
note("genes with the MT- prefix in the object",
     "%d: %s" % (len(mt_prefix), ", ".join(mt_prefix)))
i_mt = gidx(MT13)
mt_umis = rowsum(MT13)
den = np.where(tot > 0, tot, np.nan)
mito_pct = 100.0 * mt_umis / den
obs_mf = num(O[OC["mito_fraction"]])
d_mf = np.nanmax(np.abs(obs_mf * 100.0 - mito_pct))
chk(d_mf < 1e-3, "13-gene %MT equals the obs mitochondrial fraction x 100",
    "max |difference| %.3g percentage points | %d MT- genes in the object"
    % (d_mf, len(mt_prefix)), warn_only=True)

C = pd.DataFrame({
    "barcode": bcs_all[keep],
    "donor_id": O[OC["donor"]].astype(str).values,
    "region": O[OC["region"]].astype(str).values,
    "supertype": O[OC["supertype"]].astype(str).values,
    "library": O[OC["library"]].astype(str).values,
    "doublet_score": num(O[OC["doublet"]]),
    "total_umis": tot,
    "n_genes": ngen,
    "mt_umis": mt_umis,
    "mito_pct": mito_pct,
})

for k in ("donor_id", "region", "supertype", "library"):
    bad = int(C[k].isin(["", "nan", "None", "NaN"]).sum())
    chk(bad == 0, "no missing %s in the kept cells" % k, "%d missing" % bad)
if N["FAIL"]:
    abort("a kept cell has no donor, region, supertype or library")

# Own-cytoplasm shares, per cent of non-MT UMIs.
P = CFG["panels"]
den_nonmt = np.where(tot - mt_umis > 0, tot - mt_umis, np.nan)
for pan in ("own_er", "ribosome", "nuclear_retained"):
    gs = P[pan]
    chk(not missing(gs), "panel %s genes present" % pan,
        "missing: %s" % (", ".join(missing(gs)) or "none"), warn_only=True)
    C["%s_fracnonmt" % pan] = 100.0 * rowsum(gs) / den_nonmt
C["cyto_share"] = C["own_er_fracnonmt"] + C["ribosome_fracnonmt"]

# Other-lineage shares (fraction of total UMIs), one column set per dataset guard.
OL = P["other_lineage"]
GUARD = P["circularity_guard"]
for ds in CFG["datasets"]:
    drop = set(GUARD[ds])
    union = []
    for lin, gs in OL.items():
        use = [g for g in gs if g not in drop]
        union += use
        C["soup_%s__%s" % (lin, ds)] = rowsum(use) / den
    C["soup_combined__%s" % ds] = rowsum(sorted(set(union))) / den
    note("other-lineage genes used for %s" % ds,
         "%d (guard removed: %s)" % (len(set(union) & set(genes)),
                                     ", ".join(sorted(drop & {g for v in
                                               OL.values() for g in v}))
                                     or "none"))
chk(not missing([g for v in OL.values() for g in v]),
    "other-lineage panel genes present",
    "missing: %s" % (", ".join(missing([g for v in OL.values() for g in v]))
                     or "none"), warn_only=True)

# Library fields per cell; a missing field is a warning.
for key, name in CFG["library_fields"].items():
    if key.startswith("_"):
        continue
    col = low.get(name.lower())
    if chk(col is not None, "library field '%s' present" % name,
           str(col) if col else "not in obs", warn_only=True):
        v = O[col]
        C["lib_%s" % key] = (v.astype(str).values if key == "method"
                             else num(v))

# == Per-supertype gene summary and matched sets ==============================
say("")
say("SECTION D -- gene summary and matched gene sets")
rng = np.random.default_rng(SEED)
st = C["supertype"].values
GS = []
g_mean, g_det = {}, {}
for s in SUPERTYPES:
    r = np.where(st == s)[0]
    m = mat[r]
    g_mean[s] = np.asarray(m.mean(axis=0)).ravel()
    g_det[s] = np.asarray((m > 0).sum(axis=0)).ravel() / float(len(r))
    GS.append(pd.DataFrame({"supertype": s, "gene": genes,
                            "mean_count": g_mean[s], "detection": g_det[s]}))
GS = pd.concat(GS, ignore_index=True)

# Genes that no matched set may use: every configured gene and every MT gene.
NP = CFG["nuclear_panels"]
nuc = {}
for k, rx in NP["prefix"].items():
    nuc[k] = sorted(g for g in genes if re.match(rx, g))
for k, gs in NP["list"].items():
    nuc[k] = [g for g in gs if g in genes]
for k, parts in NP["unions"].items():
    nuc[k] = sorted({g for p in parts for g in nuc[p]})
for k in nuc:
    chk(len(nuc[k]) > 0, "nuclear panel %s resolved" % k,
        "%d genes" % len(nuc[k]), warn_only=True)
chk(not (set(nuc["nuclear_mito_all"]) &
         set(MT13 + CFG["genes"]["cdfg"] + CFG["genes"]["chkg8"] + ["ACTB"])),
    "nuclear panels do not overlap either composite", "checked")

named = (set(MT13) | set(CFG["genes"]["cdfg"]) | set(CFG["genes"]["chkg8"])
         | {"ACTB"} | {g for k in ("own_er", "ribosome", "nuclear_retained")
                       for g in P[k]}
         | {g for v in OL.values() for g in v} | set(P["microglia"])
         | {g for v in GUARD.values() if isinstance(v, list) for g in v}
         | {g for v in nuc.values() for g in v})
excl_prefix = tuple(P["ctrl_pct"]["exclude_prefixes"])
SETS = []


def pool_for(s):
    return np.array([i for i, g in enumerate(genes)
                     if g not in named and not g.startswith(excl_prefix)
                     and g_det[s][i] > 0.01])


def nearest(s, pool, targets, k):
    """For each target gene, the k pool genes closest on log1p mean and detection."""
    out = []
    for i in targets:
        d = (np.abs(np.log1p(g_mean[s][pool]) - np.log1p(g_mean[s][i]))
             + np.abs(g_det[s][pool] - g_det[s][i]))
        out.append(pool[np.argsort(d, kind="stable")][:k])
    return out


ctrl_obs = np.full(len(C), np.nan)
ctrl_null = np.full(len(C), np.nan)
ctrl_null_sd = np.full(len(C), np.nan)
ref_genes = set()
K = int(P["ctrl_pct"]["pool_k"])
NDRAW = int(P["ctrl_pct"]["n_draws"])
for s in SUPERTYPES:
    r = np.where(st == s)[0]
    m = mat[r].tocsc()
    pool = pool_for(s)
    # ctrl_pct: one matched gene per MT gene, then NDRAW random matched sets.
    cand = nearest(s, pool, i_mt, K)
    pick, used = [], set()
    for c in cand:
        free = [int(x) for x in c if int(x) not in used]
        if not free:
            continue
        pick.append(free[0])
        used.add(free[0])
    chk(len(pick) == len(i_mt), "%s ctrl_pct matched one gene per MT gene" % s,
        "%d genes | pool %d" % (len(pick), len(pool)))
    SETS += [{"supertype": s, "set": "ctrl_pct", "draw": 0,
              "gene": genes[j]} for j in pick]
    dn = den[r]
    ctrl_obs[r] = 100.0 * np.asarray(m[:, pick].sum(axis=1)).ravel() / dn
    draws = np.zeros((len(r), NDRAW), dtype=np.float32)
    for b in range(NDRAW):
        used, pk = set(), []
        for c in cand:
            cc = [int(x) for x in c if int(x) not in used] or \
                [int(x) for x in c]
            j = int(rng.choice(cc))
            pk.append(j)
            used.add(j)
        SETS += [{"supertype": s, "set": "ctrl_pct_null", "draw": b + 1,
                  "gene": genes[j]} for j in pk]
        draws[:, b] = 100.0 * np.asarray(m[:, pk].sum(axis=1)).ravel() / dn
    ctrl_null[r] = draws.mean(axis=1)
    ctrl_null_sd[r] = draws.std(axis=1)
    # 18A reference: reference_n pool genes nearest to random nuclear_mito_all genes.
    tgt = [genes.get_loc(g) for g in nuc["nuclear_mito_all"]]
    if not tgt:
        abort("no nuclear_mito_all gene is in the object")
    cand = nearest(s, pool, tgt, K)
    ref, used = [], set()
    order = rng.permutation(len(tgt))
    t = 0
    while len(ref) < int(NP["reference_n"]) and t < 50 * len(tgt):
        c = cand[order[t % len(tgt)]]
        free = [int(x) for x in c if int(x) not in used]
        if free:
            j = free[0]
            ref.append(j)
            used.add(j)
        t += 1
    chk(len(ref) == int(NP["reference_n"]),
        "%s 18A reference has %d genes" % (s, NP["reference_n"]),
        "%d picked" % len(ref))
    SETS += [{"supertype": s, "set": "reference_18A", "draw": 0,
              "gene": genes[j]} for j in ref]
    ref_genes |= {genes[j] for j in ref}
    SETS += [{"supertype": s, "set": "nuclear_%s" % k, "draw": 0, "gene": g}
             for k in sorted(nuc) for g in nuc[k]]
C["ctrl_pct"] = ctrl_obs
C["ctrl_pct_null_mean"] = ctrl_null
C["ctrl_pct_null_sd"] = ctrl_null_sd
SETS = pd.DataFrame(SETS)

# == Raw counts for the genes later R steps need ==============================
say("")
say("SECTION E -- raw counts for later steps")
ctrl_genes = set(SETS.loc[SETS["set"] == "ctrl_pct", "gene"])
want = (set(MT13) | {"ACTB"} | set(CFG["genes"]["chkg8"])
        | {g for k in ("own_er", "ribosome", "nuclear_retained")
           for g in P[k]}
        | {g for v in OL.values() for g in v} | set(P["microglia"])
        | {g for v in nuc.values() for g in v} | ref_genes | ctrl_genes)
cols = sorted(genes.get_loc(g) for g in want if g in genes)
note("genes exported", "%d" % len(cols))
X = matc[:, cols].tocoo()

# == Write and read back ======================================================
say("")
say("SECTION F -- write and verify")
# mtime=0 keeps the gzip header free of a timestamp.
with open(F_CELLS, "wb") as raw, \
        gzip.GzipFile(fileobj=raw, mode="wb", compresslevel=6, mtime=0) as gz, \
        io.TextIOWrapper(gz, encoding="utf-8", newline="") as fh:
    C.to_csv(fh, index=False, float_format="%.10g")
with open(F_MTX, "wb") as raw, \
        gzip.GzipFile(fileobj=raw, mode="wb", compresslevel=6, mtime=0) as fh:
    scipy.io.mmwrite(fh, X, field="integer")
pd.Series(genes[cols]).to_csv(F_GENES, index=False, header=False)
C["barcode"].to_csv(F_BCS, index=False, header=False)
GS.to_csv(F_GSUM, index=False, float_format="%.10g")
SETS.to_csv(F_SETS, index=False)

back = pd.read_csv(F_CELLS, compression="gzip", usecols=["barcode"])
chk(len(back) == len(C) and back["barcode"].is_unique,
    "cell table written and read back", "%d rows" % len(back))
with gzip.open(F_MTX, "rb") as fh:
    Xb = scipy.io.mmread(fh)
chk(Xb.shape == (len(C), len(cols)) and Xb.nnz == X.nnz,
    "count matrix written and read back",
    "%d cells x %d genes | %d non-zero" % (Xb.shape[0], Xb.shape[1], Xb.nnz))
say("")
say("Cells per supertype:")
for s in SUPERTYPES:
    r = st == s
    say("  %-20s %7d cells | %4d donors | %2d regions | median UMIs %.0f | "
        "median %%MT %.3f"
        % (s, int(r.sum()), C.loc[r, "donor_id"].nunique(),
           C.loc[r, "region"].nunique(), float(np.median(tot[r])),
           float(np.nanmedian(mito_pct[r]))))
A.file.close()
finish()
