# =============================================================================
# Part7_DecontX.py
# SEA-AD CEP-IP/CMRF, Part 7, DecontX ambient RNA estimate (Table S6)
#   Runs DecontX on all nuclei of the immune h5ad, per library, with supertype
#   as the cluster label and without empty droplets. Reports the contamination
#   per nucleus and rho(MT-CO3, cDFG) before and after correction.
#
# Usage (from the SEA-AD folder, or set SEAAD_ROOT)
#   python3 Part7_DecontX.py
#   DECONTX_WORKERS=6 python3 Part7_DecontX.py   # fix the number of workers
#   Run with python3 directly (the workers are spawned processes). An
#   interrupted run resumes from Part7_DecontX/_checkpoint.json; delete that
#   file for a full rerun.
#
# Outputs (Part7_DecontX/): contamination per nucleus, module and mito
#   diagnostics (CSV), corrected h5ad files and log.
# =============================================================================

# == One BLAS thread per worker; set before numpy is imported ==============
import os
os.environ["OMP_NUM_THREADS"]       = "1"
os.environ["MKL_NUM_THREADS"]       = "1"
os.environ["OPENBLAS_NUM_THREADS"]  = "1"
os.environ["NUMBA_NUM_THREADS"]     = "1"
os.environ["NUMBA_THREADING_LAYER"] = "workqueue"

import sys
import time
import gc
import json
import warnings
import numpy as np
import pandas as pd
from scipy.sparse import issparse, csr_matrix, vstack as sparse_vstack
from scipy.sparse import save_npz, load_npz
from scipy.stats import spearmanr
from multiprocessing import shared_memory
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing

# == Config ==============================================================
H5AD_PATH = (
    os.environ.get("SEAAD_ROOT", os.getcwd()) + "/"
    "SEAAD_Immune_HIP_MEC_LEC_ITG_MTG_FI_STG_DFC_AnG_V1C_"
    "RNAseq_final-nuclei.2026-06-22.h5ad"
)
OUT_DIR       = os.environ.get("SEAAD_ROOT", os.getcwd()) + "/Part7_DecontX"
TMP_DIR       = os.environ.get("SEAAD_ROOT", os.getcwd()) + "/Part7_DecontX/_tmp_batches"
CHECKPOINT    = os.environ.get("SEAAD_ROOT", os.getcwd()) + "/Part7_DecontX/_checkpoint.json"
TARGET_STYPE  = "Micro-PVM_3-SEAAD"
CLUSTER_KEY   = "Supertype"
BATCH_KEY     = "library_prep"
COUNTS_LAYER  = "UMIs"
GOI           = "MT-CO3"
MODULE_GENES  = [
    "MT-CO3", "MT-ATP6", "MT-CO2", "MT-ND3", "MT-ND4",
    "MT-CYB", "MT-ND1", "MT-ND2", "MT-CO1",
]
DFG_GENES     = [g for g in MODULE_GENES if g != GOI]
ALL_MT_GENES  = [
    "MT-CO3", "MT-ATP6", "MT-CO2", "MT-ND3", "MT-ND4",
    "MT-CYB", "MT-ND1", "MT-ND2", "MT-CO1",
    "MT-ND5", "MT-ND4L", "MT-ATP8", "MT-ND6",
]
NEOC_REGIONS  = {"DFC", "MTG", "AnG", "STG", "V1C", "ITG"}

N_WORKERS           = 18      # upper bound; reduced at runtime to fit RAM
MEM_BUDGET_FRAC     = 0.60    # fraction of MemAvailable the pool may use
DENSE_COPIES        = 5.0     # dense working copies decontx holds per batch
MIN_CELLS_PER_BATCH = 20
MIN_CLUSTERS        = 2
DECONTX_MAX_ITER    = 500
DECONTX_CONVERGENCE = 0.001


# ========================================================================
# Worker functions, at module level so that spawned workers can load them
# ========================================================================

_W = {}  # Worker-local state, populated by _init_worker


def _init_worker(shm_names, arr_shapes, arr_dtypes, matrix_shape,
                 cluster_shm_name, cluster_len, var_names_json):
    """Attach a spawned worker to the shared CSR matrix and cluster labels."""
    global _W
    warnings.filterwarnings("ignore")

    # Attach to the three CSR SharedMemory blocks
    _W["shm_list"] = []
    arrays = []
    for name, shape, dtype_str in zip(shm_names, arr_shapes, arr_dtypes):
        shm = shared_memory.SharedMemory(name=name)
        _W["shm_list"].append(shm)
        arr = np.ndarray(shape, dtype=np.dtype(dtype_str), buffer=shm.buf)
        arrays.append(arr)

    _W["csr"] = csr_matrix(
        (arrays[0], arrays[1], arrays[2]),
        shape=matrix_shape, copy=False,
    )

    # Attach to cluster labels SharedMemory (encoded as fixed-width bytes)
    shm_cl = shared_memory.SharedMemory(name=cluster_shm_name)
    _W["shm_list"].append(shm_cl)
    # Decode fixed-width byte strings back to numpy string array
    raw = np.ndarray((cluster_len,), dtype="S64", buffer=shm_cl.buf)
    _W["clusters"] = np.array([s.decode("utf-8") for s in raw])

    # Gene names for building AnnData var
    _W["var_names"] = json.loads(var_names_json)


def _run_batch(task):
    """Process one batch in a spawn worker. Returns contam + saves npz."""
    import anndata as ad
    import decontx as dcx

    batch_label = task["batch"]
    batch_idx   = np.array(task["idx"])
    tmp_dir     = task["tmp_dir"]
    cluster_key = task["cluster_key"]
    max_iter    = task["max_iter"]
    convergence = task["convergence"]
    min_cells   = task["min_cells"]
    min_clusters = task["min_clusters"]

    batch_id   = task["batch_id"]  # unique sequential index, no collisions
    n_cells    = len(batch_idx)
    clusters   = _W["clusters"][batch_idx]
    n_clusters = len(np.unique(clusters))
    npz_path   = os.path.join(tmp_dir, f"batch_{batch_id:04d}.npz")

    if n_cells < min_cells or n_clusters < min_clusters:
        dcx_block = csr_matrix(_W["csr"][batch_idx, :])
        save_npz(npz_path, dcx_block)
        return {
            "batch": batch_label, "n_cells": n_cells,
            "n_clusters": n_clusters, "skipped": True,
            "contam": np.zeros(n_cells, dtype=np.float64).tolist(),
            "npz_path": npz_path,
        }

    n_genes_full = _W["csr"].shape[1]
    counts_batch = _W["csr"][batch_idx, :].copy()

    # Genes with no count in this batch are dropped before DecontX (they do
    # not change the result) and restored afterwards.
    keep = np.unique(counts_batch.indices)
    counts_sub = csr_matrix(counts_batch[:, keep])
    del counts_batch

    var_df = pd.DataFrame(index=[_W["var_names"][j] for j in keep])
    adata_b = ad.AnnData(
        X=counts_sub,
        obs=pd.DataFrame(
            {cluster_key: clusters},
            index=[str(i) for i in range(n_cells)],
        ),
        var=var_df,
    )

    dcx.decontx(
        adata_b,
        cluster_key=cluster_key,
        max_iter=max_iter,
        convergence=convergence,
    )

    sub = csr_matrix(adata_b.layers["decontX_counts"])
    contam_out = adata_b.obs["decontX_contamination"].values.tolist()
    del adata_b
    gc.collect()

    # Map back onto the full gene axis.
    dcx_block = csr_matrix(
        (sub.data, keep[sub.indices], sub.indptr),
        shape=(n_cells, n_genes_full),
    )
    save_npz(npz_path, dcx_block)

    del sub, dcx_block, counts_sub
    gc.collect()

    return {
        "batch": batch_label, "n_cells": n_cells,
        "n_clusters": n_clusters, "skipped": False,
        "contam": contam_out, "npz_path": npz_path,
    }


# ========================================================================
# Helper functions
# ========================================================================

def rss_gb():
    try:
        with open("/proc/self/status") as f:
            for line in f:
                if line.startswith("VmRSS:"):
                    return int(line.split()[1]) / 1_048_576
    except Exception:
        return -1.0
    return -1.0


def mem_available_gb():
    try:
        with open("/proc/meminfo") as f:
            for line in f:
                if line.startswith("MemAvailable:"):
                    return int(line.split()[1]) / 1_048_576
    except Exception:
        pass
    return -1.0


def find_col(obs, candidates):
    for c in candidates:
        if c in obs.columns:
            return c
    raise KeyError(f"None of {candidates} found in obs")


def compute_module_rho(count_mat, var_names, goi, dfgs):
    gene_idx = {g: i for i, g in enumerate(var_names)}
    missing = [g for g in [goi] + dfgs if g not in gene_idx]
    if missing:
        return np.nan, 0
    # Only the nine module columns are densified.
    cols = [gene_idx[goi]] + [gene_idx[g] for g in dfgs]
    if issparse(count_mat):
        row_sums = np.asarray(count_mat.sum(axis=1)).ravel().astype(np.float64)
        sub = np.asarray(count_mat[:, cols].todense(), dtype=np.float64)
    else:
        arr = np.asarray(count_mat, dtype=np.float64)
        row_sums = arr.sum(axis=1)
        sub = arr[:, cols]
    row_sums[row_sums == 0] = 1.0
    norm = np.log1p(sub / row_sums[:, None] * 10_000)
    goi_vals = norm[:, 0]
    dfg_vals = norm[:, 1:].mean(axis=1)
    pos = goi_vals > 0
    n_pos = int(pos.sum())
    if n_pos < 10:
        return np.nan, n_pos
    rho, _ = spearmanr(goi_vals[pos], dfg_vals[pos])
    return rho, n_pos


def fmt_eta(seconds):
    h = int(seconds // 3600)
    m = int((seconds % 3600) // 60)
    return f"{h}h{m:02d}m"


def save_checkpoint(completed, contam_dict):
    obj = {"completed": list(completed),
           "contam": {b: c for b, c in contam_dict.items()}}
    tmp = CHECKPOINT + ".tmp"
    with open(tmp, "w") as f:
        json.dump(obj, f)
    os.replace(tmp, CHECKPOINT)


def load_checkpoint():
    if not os.path.exists(CHECKPOINT):
        return set(), {}
    with open(CHECKPOINT, "r") as f:
        obj = json.load(f)
    return set(obj["completed"]), obj["contam"]


class _Tee:
    """Duplicate a stream to both the terminal and a log file."""

    def __init__(self, terminal, logfile):
        self.terminal = terminal
        self.logfile = logfile
        self.encoding = getattr(terminal, "encoding", "utf-8")

    def write(self, message):
        self.terminal.write(message)
        self.logfile.write(message)
        self.logfile.flush()

    def flush(self):
        self.terminal.flush()
        self.logfile.flush()

    def fileno(self):
        return self.terminal.fileno()

    def isatty(self):
        return self.terminal.isatty()


# ========================================================================
# Main
# ========================================================================

if __name__ == "__main__":
    import shutil

    os.makedirs(OUT_DIR, exist_ok=True)
    os.makedirs(TMP_DIR, exist_ok=True)

    # -- Built-in log file: captures everything printed to stdout/stderr --
    _log_ts = time.strftime("%Y%m%d_%H%M%S")
    LOG_PATH = os.path.join(OUT_DIR, f"decontx_log_{_log_ts}.txt")
    _logfile = open(LOG_PATH, "w")
    sys.stdout = _Tee(sys.__stdout__, _logfile)
    sys.stderr = _Tee(sys.__stderr__, _logfile)

    t0 = time.time()
    print("=" * 70)
    print("Part7_DecontX.py")
    print("Part 7, DecontX without empty droplets, spawned worker processes")
    print("=" * 70)
    print(f"  RSS at start: {rss_gb():.1f} GB")

    import anndata as ad
    import scanpy as sc

    # == Load ==
    print(f"\nLoading h5ad ...")
    adata_full = sc.read_h5ad(H5AD_PATH)
    print(f"  Shape: {adata_full.shape[0]:,} x {adata_full.shape[1]:,}")
    print(f"  Loaded in {time.time() - t0:.0f}s,  RSS: {rss_gb():.1f} GB")

    # == Extract counts ==
    assert COUNTS_LAYER in adata_full.layers
    counts_csr = csr_matrix(adata_full.layers[COUNTS_LAYER])
    print(f"  CSR: {counts_csr.shape},  nnz: {counts_csr.nnz:,},  "
          f"data: {counts_csr.data.nbytes / 1e9:.2f} GB")

    # == Resolve obs columns ==
    cluster_col = find_col(adata_full.obs, [CLUSTER_KEY, "supertype"])
    batch_col   = find_col(adata_full.obs, [BATCH_KEY, "Library_prep"])
    region_col  = find_col(adata_full.obs, ["Brain Region", "brain_region"])

    clusters_arr = np.asarray(adata_full.obs[cluster_col].values, dtype=str)
    batches_arr  = np.asarray(adata_full.obs[batch_col].values, dtype=str)
    regions_arr  = np.asarray(adata_full.obs[region_col].values, dtype=str)
    obs_index    = adata_full.obs.index.copy()
    var_names    = list(adata_full.var_names)
    var_df       = adata_full.var[[]].copy()

    unique_batches = sorted(np.unique(batches_arr))
    supertypes     = sorted(np.unique(clusters_arr))
    print(f"  Clusters: {len(supertypes)},  Batches: {len(unique_batches)}")
    assert TARGET_STYPE in clusters_arr

    batch_indices = {}
    for b in unique_batches:
        batch_indices[b] = np.where(batches_arr == b)[0]

    n_skip = sum(1 for b in unique_batches
                 if len(batch_indices[b]) < MIN_CELLS_PER_BATCH
                 or len(np.unique(clusters_arr[batch_indices[b]])) < MIN_CLUSTERS)
    print(f"  To process: {len(unique_batches) - n_skip}  (skip {n_skip})")

    # Free the full AnnData
    del adata_full
    gc.collect()
    print(f"  Released adata,  RSS: {rss_gb():.1f} GB")

    # == Choose worker count to fit RAM ==================================
    # Memory per worker is estimated from the five largest batches.
    print(f"\n--- Sizing the worker pool ---")
    sizes = np.array([len(batch_indices[b]) for b in unique_batches])
    probe = [unique_batches[i] for i in np.argsort(sizes)[-5:]]
    peak_bytes = 0.0
    for b in probe:
        blk = counts_csr[batch_indices[b], :]
        n_g = int(np.unique(blk.indices).size)
        peak_bytes = max(peak_bytes, blk.shape[0] * n_g * 8.0 * DENSE_COPIES)
        del blk
    gc.collect()

    avail_gb = mem_available_gb()
    shm_gb = (counts_csr.data.nbytes + counts_csr.indices.nbytes
              + counts_csr.indptr.nbytes) / 1e9
    budget_gb = max(avail_gb * MEM_BUDGET_FRAC - shm_gb, 2.0)
    fit = int(budget_gb * 1e9 / max(peak_bytes, 1.0))
    N_WORKERS = int(os.environ.get("DECONTX_WORKERS", 0)) or \
                max(1, min(N_WORKERS, fit))

    print(f"  Largest batch: {int(sizes.max()):,} cells, "
          f"peak estimate {peak_bytes / 1e9:.2f} GB/worker")
    print(f"  MemAvailable {avail_gb:.1f} GB, shared CSR {shm_gb:.1f} GB, "
          f"budget {budget_gb:.1f} GB")
    print(f"  Workers set to {N_WORKERS}")

    # == Create SharedMemory blocks ======================================
    print(f"\n--- Setting up SharedMemory ---")

    # 1. CSR arrays: data, indices, indptr
    shm_data = shared_memory.SharedMemory(
        create=True, size=counts_csr.data.nbytes)
    shm_indices = shared_memory.SharedMemory(
        create=True, size=counts_csr.indices.nbytes)
    shm_indptr = shared_memory.SharedMemory(
        create=True, size=counts_csr.indptr.nbytes)

    # Copy into shared memory
    np.ndarray(counts_csr.data.shape, dtype=counts_csr.data.dtype,
               buffer=shm_data.buf)[:] = counts_csr.data
    np.ndarray(counts_csr.indices.shape, dtype=counts_csr.indices.dtype,
               buffer=shm_indices.buf)[:] = counts_csr.indices
    np.ndarray(counts_csr.indptr.shape, dtype=counts_csr.indptr.dtype,
               buffer=shm_indptr.buf)[:] = counts_csr.indptr

    shm_names  = [shm_data.name, shm_indices.name, shm_indptr.name]
    arr_shapes = [counts_csr.data.shape, counts_csr.indices.shape,
                  counts_csr.indptr.shape]
    arr_dtypes = [str(counts_csr.data.dtype), str(counts_csr.indices.dtype),
                  str(counts_csr.indptr.dtype)]
    matrix_shape = counts_csr.shape

    # 2. Cluster labels: encode as fixed-width S64 bytes
    clusters_bytes = np.array(clusters_arr, dtype="S64")
    shm_clusters = shared_memory.SharedMemory(
        create=True, size=clusters_bytes.nbytes)
    np.ndarray(clusters_bytes.shape, dtype="S64",
               buffer=shm_clusters.buf)[:] = clusters_bytes

    # 3. Gene names as JSON (sent to initializer)
    var_names_json = json.dumps(var_names)

    print(f"  SharedMemory blocks created:")
    print(f"    data:     {shm_data.name}  ({counts_csr.data.nbytes / 1e9:.2f} GB)")
    print(f"    indices:  {shm_indices.name}  ({counts_csr.indices.nbytes / 1e9:.2f} GB)")
    print(f"    indptr:   {shm_indptr.name}")
    print(f"    clusters: {shm_clusters.name}")

    # Free the parent's copy of the CSR
    del counts_csr, clusters_bytes
    gc.collect()
    print(f"  Released parent CSR copy,  RSS: {rss_gb():.1f} GB")

    # == Checkpoint ======================================================
    # Map batch labels to sequential IDs
    batch_to_id = {b: i for i, b in enumerate(unique_batches)}

    completed_batches, batch_contam = load_checkpoint()
    batch_contam = {b: c for b, c in batch_contam.items()}
    remaining = [b for b in unique_batches if b not in completed_batches]
    if completed_batches:
        # Verify npz files still exist
        for b in list(completed_batches):
            npz = os.path.join(TMP_DIR, f"batch_{batch_to_id[b]:04d}.npz")
            if not os.path.exists(npz):
                completed_batches.discard(b)
                batch_contam.pop(b, None)
        remaining = [b for b in unique_batches if b not in completed_batches]
        print(f"\n  Resuming: {len(completed_batches)} done, "
              f"{len(remaining)} remaining")
    else:
        print(f"\n  Starting fresh: {len(remaining)} batches")

    # == Spawn pool ======================================================
    print(f"\n{'=' * 70}")
    print("Running DecontX")
    print(f"  MODE:       no empty droplets")
    print(f"  WORKERS:    {N_WORKERS}")
    print(f"  BATCHES:    {len(unique_batches)} total, "
          f"{len(remaining)} remaining")
    print(f"{'=' * 70}")

    t1 = time.time()
    ctx = multiprocessing.get_context("spawn")
    done = len(completed_batches)
    skipped = 0
    batch_times = []

    # Build task list
    tasks = []
    for b in remaining:
        tasks.append({
            "batch":        b,
            "batch_id":     batch_to_id[b],
            "idx":          batch_indices[b].tolist(),
            "tmp_dir":      TMP_DIR,
            "cluster_key":  CLUSTER_KEY,
            "max_iter":     DECONTX_MAX_ITER,
            "convergence":  DECONTX_CONVERGENCE,
            "min_cells":    MIN_CELLS_PER_BATCH,
            "min_clusters": MIN_CLUSTERS,
        })

    with ProcessPoolExecutor(
        max_workers=N_WORKERS,
        mp_context=ctx,
        initializer=_init_worker,
        initargs=(shm_names, arr_shapes, arr_dtypes, matrix_shape,
                  shm_clusters.name, len(clusters_arr), var_names_json),
    ) as pool:
        futures = {pool.submit(_run_batch, t): t["batch"] for t in tasks}

        for future in as_completed(futures):
            batch_label = futures[future]
            t_batch = time.time()
            try:
                res = future.result()
                batch_contam[res["batch"]] = res["contam"]
                completed_batches.add(res["batch"])
                done += 1
                if res["skipped"]:
                    skipped += 1
            except Exception as e:
                print(f"  ERROR on {batch_label}: {e}")
                idx = batch_indices[batch_label]
                batch_contam[batch_label] = [0.0] * len(idx)
                bid = batch_to_id[batch_label]
                npz_path = os.path.join(TMP_DIR, f"batch_{bid:04d}.npz")
                # Read from shared memory to save original counts
                d = np.ndarray(arr_shapes[0], dtype=arr_dtypes[0],
                               buffer=shm_data.buf)
                ii = np.ndarray(arr_shapes[1], dtype=arr_dtypes[1],
                                buffer=shm_indices.buf)
                p = np.ndarray(arr_shapes[2], dtype=arr_dtypes[2],
                               buffer=shm_indptr.buf)
                tmp_csr = csr_matrix((d, ii, p), shape=matrix_shape, copy=False)
                save_npz(npz_path, csr_matrix(tmp_csr[idx, :]))
                del tmp_csr, d, ii, p
                completed_batches.add(batch_label)
                done += 1
                skipped += 1

            # Checkpoint every 20 batches
            if done % 20 == 0:
                save_checkpoint(completed_batches, batch_contam)

            # Progress
            batch_times.append(time.time() - t_batch)
            n_rem = len(unique_batches) - done
            avg_t = np.mean(batch_times[-50:])
            eta = avg_t * n_rem / max(N_WORKERS, 1)
            if done % 20 == 0 or done == len(unique_batches):
                elapsed = time.time() - t1
                print(f"  [{done:>4d}/{len(unique_batches)}]  "
                      f"{elapsed / 60:.1f} min  "
                      f"RSS={rss_gb():.1f} GB  "
                      f"free={mem_available_gb():.1f} GB  "
                      f"ETA≈{fmt_eta(eta)}")

    save_checkpoint(completed_batches, batch_contam)
    elapsed_dcx = time.time() - t1
    print(f"\nDecontX finished: {done} batches in {elapsed_dcx / 60:.1f} min "
          f"({skipped} skipped),  RSS: {rss_gb():.1f} GB")

    # == Reassemble from disk ============================================
    # The final CSR is allocated once and each batch is written into its rows.
    print(f"\n--- Reassembling from disk ---")
    t_re = time.time()

    n_total = matrix_shape[0]
    n_genes = matrix_shape[1]

    contam_full = np.zeros(n_total, dtype=np.float64)
    for b in unique_batches:
        contam_full[batch_indices[b]] = np.asarray(batch_contam[b],
                                                   dtype=np.float64)

    # Pass 1: per-row nnz
    row_nnz = np.zeros(n_total, dtype=np.int64)
    for i, b in enumerate(unique_batches):
        blk = load_npz(os.path.join(TMP_DIR, f"batch_{i:04d}.npz"))
        row_nnz[batch_indices[b]] = np.diff(blk.indptr)
        del blk
        if (i + 1) % 200 == 0:
            print(f"    Scanned {i + 1}/{len(unique_batches)}")

    indptr = np.zeros(n_total + 1, dtype=np.int64)
    np.cumsum(row_nnz, out=indptr[1:])
    total_nnz = int(indptr[-1])
    del row_nnz

    idx_dtype = np.int64 if total_nnz > np.iinfo(np.int32).max else np.int32
    starts   = indptr.copy()                       # int64 working copy
    indptr   = indptr.astype(idx_dtype, copy=False)
    data_out = np.empty(total_nnz, dtype=np.float64)
    ind_out  = np.empty(total_nnz, dtype=idx_dtype)

    # Pass 2: scatter
    for i, b in enumerate(unique_batches):
        blk = load_npz(os.path.join(TMP_DIR, f"batch_{i:04d}.npz"))
        idx = batch_indices[b]
        counts = np.diff(blk.indptr)
        dest = (np.repeat(starts[idx], counts)
                + np.arange(len(blk.data), dtype=np.int64)
                - np.repeat(blk.indptr[:-1].astype(np.int64), counts))
        data_out[dest] = blk.data
        ind_out[dest]  = blk.indices
        del blk, dest, counts
        if (i + 1) % 200 == 0:
            print(f"    Loaded {i + 1}/{len(unique_batches)}")

    del starts
    dcx_full = csr_matrix((data_out, ind_out, indptr),
                          shape=(n_total, n_genes))
    del data_out, ind_out, indptr
    gc.collect()

    print(f"  Reassembled: {dcx_full.shape}, nnz={dcx_full.nnz:,}")
    print(f"  Took {time.time() - t_re:.1f}s,  RSS: {rss_gb():.1f} GB")

    # Raw counts as a view on the shared memory, which is released only after
    # the h5ad files are written.
    counts_csr = csr_matrix(
        (np.ndarray(arr_shapes[0], dtype=arr_dtypes[0], buffer=shm_data.buf),
         np.ndarray(arr_shapes[1], dtype=arr_dtypes[1], buffer=shm_indices.buf),
         np.ndarray(arr_shapes[2], dtype=arr_dtypes[2], buffer=shm_indptr.buf)),
        shape=matrix_shape, copy=False,
    )

    # Build output AnnData
    obs_df = pd.DataFrame({
        CLUSTER_KEY:             clusters_arr,
        BATCH_KEY:               batches_arr,
        "Brain Region":          regions_arr,
        "arm":                   ["NeoC" if r in NEOC_REGIONS else "NonNeoC"
                                  for r in regions_arr],
        "decontX_contamination": contam_full,
    }, index=obs_index)

    adata_dcx = ad.AnnData(X=counts_csr, obs=obs_df, var=var_df)
    adata_dcx.layers["decontX_counts"] = dcx_full
    del dcx_full, counts_csr; gc.collect()

    contam = adata_dcx.obs["decontX_contamination"]

    # == Contamination summary ===========================================
    print(f"\n{'=' * 70}")
    print(f"Contamination summary (all {adata_dcx.shape[0]:,} cells)")
    print(f"{'=' * 70}")
    print(f"  mean   {contam.mean():.4f}")
    print(f"  median {contam.median():.4f}")
    print(f"  q25    {contam.quantile(0.25):.4f}")
    print(f"  q75    {contam.quantile(0.75):.4f}")
    print(f"  q95    {contam.quantile(0.95):.4f}")
    print(f"  max    {contam.max():.4f}")

    print(f"\n--- Per-supertype contamination ---")
    for st in supertypes:
        mask = adata_dcx.obs[CLUSTER_KEY] == st
        c = contam[mask]
        print(f"  {st:30s}  n={int(mask.sum()):>7,}  "
              f"mean={c.mean():.4f}  median={c.median():.4f}  "
              f"q95={c.quantile(0.95):.4f}")

    # == Mito-gene diagnostic ============================================
    print(f"\n{'=' * 70}")
    print(f"Mito-gene contamination diagnostic")
    print(f"{'=' * 70}")

    var_names_list = list(adata_dcx.var_names)
    gim = {g: i for i, g in enumerate(var_names_list)}

    REF = ["ACTB", "B2M", "GAPDH", "TMSB4X", "FTL"]
    cg = [g for g in ALL_MT_GENES + REF if g in gim]
    cg_cols = [gim[g] for g in cg]
    sub_pos = {g: k for k, g in enumerate(cg)}

    pvm3_mask = adata_dcx.obs[CLUSTER_KEY] == TARGET_STYPE
    pvm3_idx = np.where(pvm3_mask.values)[0]
    # Only the diagnostic genes are densified.
    rp_d = np.asarray(adata_dcx.X[pvm3_idx, :][:, cg_cols].todense(),
                      dtype=np.float64)
    dp_d = np.asarray(
        adata_dcx.layers["decontX_counts"][pvm3_idx, :][:, cg_cols].todense(),
        dtype=np.float64)

    mdr = []

    print(f"\n  Per-gene in {TARGET_STYPE} (n={len(pvm3_idx):,}):")
    print(f"  {'gene':12s}  {'before':>8s}  {'after':>8s}  "
          f"{'retained':>8s}  {'removed':>8s}  {'module':>8s}")
    print(f"  {'-'*12}  {'-'*8}  {'-'*8}  {'-'*8}  {'-'*8}  {'-'*8}")

    for g in cg:
        gi = sub_pos[g]
        rm = rp_d[:, gi].mean(); dm = dp_d[:, gi].mean()
        ret = dm / rm if rm > 0 else np.nan
        rem = 1.0 - ret if rm > 0 else np.nan
        im = "MODULE" if g in MODULE_GENES else ""
        mdr.append({"gene": g, "in_module": g in MODULE_GENES,
                     "is_mito": g.startswith("MT-"),
                     "mean_raw": rm, "mean_decontX": dm,
                     "fraction_retained": ret, "fraction_removed": rem})
        print(f"  {g:12s}  {rm:8.3f}  {dm:8.3f}  "
              f"{ret:8.4f}  {rem:8.4f}  {im:>8s}")

    for label, filt in [("Module", lambda r: r["in_module"]),
                        ("Non-mod mito", lambda r: r["is_mito"] and not r["in_module"]),
                        ("Non-mito ref", lambda r: not r["is_mito"])]:
        vals = [r["fraction_retained"] for r in mdr
                if filt(r) and not np.isnan(r["fraction_retained"])]
        if vals:
            print(f"    {label:16s} mean {np.mean(vals):.4f}  "
                  f"[{min(vals):.4f}, {max(vals):.4f}]")

    # MT-CO3 quartile
    print(f"\n--- Contamination by MT-CO3 quartile ({TARGET_STYPE}) ---")
    pc = contam[pvm3_mask].values
    if GOI in sub_pos:
        gr = rp_d[:, sub_pos[GOI]]
        for lb, mq in [("= 0", gr == 0), ("> 0", gr > 0)]:
            if mq.sum() > 0:
                print(f"  MT-CO3 {lb}  (n={mq.sum():>6,}):  "
                      f"mean={pc[mq].mean():.4f}  med={np.median(pc[mq]):.4f}")
        gp = gr > 0
        if gp.sum() > 10:
            qe = np.percentile(gr[gp], [0, 25, 50, 75, 100])
            for qi, ql in enumerate(["Q1", "Q2", "Q3", "Q4"]):
                lo, hi = qe[qi], qe[qi + 1]
                qm = ((gr[gp] >= lo) & (gr[gp] < hi)) if qi < 3 \
                    else ((gr[gp] >= lo) & (gr[gp] <= hi))
                if qm.sum() > 0:
                    print(f"  {ql} (n={qm.sum():>6,}, {lo:.0f}-{hi:.0f}):  "
                          f"mean={pc[gp][qm].mean():.4f}  "
                          f"med={np.median(pc[gp][qm]):.4f}")

    del rp_d, dp_d; gc.collect()

    # == Module diagnostic ===============================================
    print(f"\n{'=' * 70}")
    print(f"Module diagnostic: rho(MT-CO3, cDFG) before and after")
    print(f"{'=' * 70}")

    diag_rows = []
    for st in supertypes:
        mask = adata_dcx.obs[CLUSTER_KEY] == st
        idx = np.where(mask.values)[0]
        rb, nb = compute_module_rho(adata_dcx.X[idx, :], var_names_list, GOI, DFG_GENES)
        ra, na = compute_module_rho(adata_dcx.layers["decontX_counts"][idx, :],
                                    var_names_list, GOI, DFG_GENES)
        d = ra - rb if not (np.isnan(ra) or np.isnan(rb)) else np.nan
        diag_rows.append({"supertype": st, "target": st == TARGET_STYPE,
                          "n_GOIpos": nb, "rho_before": rb,
                          "rho_after": ra, "delta": d})
        m = " <<<" if st == TARGET_STYPE else ""
        print(f"  {st:30s}  n={nb:>7,}  "
              f"before={rb:+.4f}  after={ra:+.4f}  delta={d:+.4f}{m}")

    print(f"\n--- {TARGET_STYPE} per arm ---")
    for arm in ["NeoC", "NonNeoC"]:
        mask = ((adata_dcx.obs[CLUSTER_KEY] == TARGET_STYPE) &
                (adata_dcx.obs["arm"] == arm))
        idx = np.where(mask.values)[0]
        if len(idx) == 0: continue
        rb, nb = compute_module_rho(adata_dcx.X[idx, :], var_names_list, GOI, DFG_GENES)
        ra, na = compute_module_rho(adata_dcx.layers["decontX_counts"][idx, :],
                                    var_names_list, GOI, DFG_GENES)
        d = ra - rb if not (np.isnan(ra) or np.isnan(rb)) else np.nan
        print(f"  {arm:8s}  n={nb:>6,}  before={rb:+.4f}  after={ra:+.4f}  "
              f"delta={d:+.4f}")
        diag_rows.append({"supertype": f"{TARGET_STYPE}_{arm}", "target": True,
                          "n_GOIpos": nb, "rho_before": rb,
                          "rho_after": ra, "delta": d})

    pr = [r for r in diag_rows if r["supertype"] == TARGET_STYPE][0]
    dd, rra = pr["delta"], pr["rho_after"]

    print(f"\n{'=' * 70}")
    if np.isnan(dd):
        print("  VERDICT: Could not compute")
    elif abs(dd) < 0.05:
        print(f"  VERDICT: module survives (delta={dd:+.4f})")
    elif rra < 0.30:
        print(f"  VERDICT: module collapses (rho_after={rra:.4f})")
    else:
        print(f"  VERDICT: module weakened but persists "
              f"(delta={dd:+.4f}, rho_after={rra:.4f})")
    print(f"{'=' * 70}")

    # == Save outputs ====================================================
    print(f"\n--- Saving outputs ---")

    adata_dcx.obs[[CLUSTER_KEY, BATCH_KEY, "Brain Region", "arm",
                    "decontX_contamination"]].to_csv(
        os.path.join(OUT_DIR, "decontX_contamination_all.csv"))
    pd.DataFrame(diag_rows).to_csv(
        os.path.join(OUT_DIR, "decontX_module_diagnostic.csv"), index=False)
    pd.DataFrame(mdr).to_csv(
        os.path.join(OUT_DIR, "decontX_mito_contamination.csv"), index=False)

    ap = adata_dcx[adata_dcx.obs[CLUSTER_KEY] == TARGET_STYPE].copy()
    ap.write_h5ad(os.path.join(OUT_DIR, "PVM3_decontX.h5ad"))
    print(f"  PVM3_decontX.h5ad  ({ap.shape[0]:,} cells)")
    del ap; gc.collect()

    adata_dcx.write_h5ad(os.path.join(OUT_DIR, "AllImmune_decontX.h5ad"))
    print(f"  AllImmune_decontX.h5ad  ({adata_dcx.shape[0]:,} cells)")

    # == Cleanup SharedMemory ============================================
    # adata_dcx.X uses these buffers, so this stays after the h5ad writes.
    del adata_dcx; gc.collect()
    shm_data.close();    shm_data.unlink()
    shm_indices.close(); shm_indices.unlink()
    shm_indptr.close();  shm_indptr.unlink()
    shm_clusters.close(); shm_clusters.unlink()

    shutil.rmtree(TMP_DIR, ignore_errors=True)
    if os.path.exists(CHECKPOINT):
        os.remove(CHECKPOINT)

    # == Summary =========================================================
    total_min = (time.time() - t0) / 60
    print(f"\n{'=' * 70}")
    print(f"SUMMARY")
    print(f"{'=' * 70}")
    print(f"  Mode:              no empty droplets")
    print(f"  Workers:           {N_WORKERS}")
    print(f"  Cells:             {n_total:,}")
    print(f"  Batches:           {len(unique_batches)}")
    print(f"  Mean contamination:{contam.mean():.4f}")
    print(f"  Module rho before: {pr['rho_before']:+.4f}")
    print(f"  Module rho after:  {pr['rho_after']:+.4f}")
    print(f"  Delta:             {dd:+.4f}")
    print(f"  DecontX time:      {elapsed_dcx / 60:.1f} min")
    print(f"  Total time:        {total_min:.1f} min")
    print(f"  Final RSS:         {rss_gb():.1f} GB")
    print(f"{'=' * 70}")
    print(f"  Log file: {LOG_PATH}")
    print("Done.")

    # -- Restore streams and close log --
    sys.stdout = sys.__stdout__
    sys.stderr = sys.__stderr__
    _logfile.close()
