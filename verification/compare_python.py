#!/usr/bin/env python3
# =============================================================================
# compare_python.py
# Compares the container rerun of Python Parts 1 to 8 (/data) with the original
# results of the Fedora runs (mounted read-only at /ref).
#   1) Tables (CSV, CSV.GZ, XLSX): compare_rerun.py is run with WIN_ROOT=/ref on
#      the top-level folders found in both, with the same rules as the R check.
#   2) Matrix Market files (.mtx, .mtx.gz) and lists (.tsv, .tsv.gz), which
#      compare_rerun.py counts by presence only: compared byte by byte first;
#      if the bytes differ, entry by entry (row and column indices exactly,
#      values within 1e-6 relative), or line by line for the lists.
#
# Usage (inside the container): python /data/compare_python.py
# Outputs (/data/compare/):
#   compare_python_tables_log.txt, _detail.csv, _cells.csv   (part 1)
#   compare_python_matrices_log.txt, _detail.csv             (part 2)
# =============================================================================
import gzip
import os
import re
import subprocess
import sys
import time

import numpy as np
import pandas as pd

NEW = os.environ.get("SEAAD_ROOT", "/data")
REF = os.environ.get("REF_ROOT", "/ref")
OUT = os.path.join(NEW, "compare")
SKIP = {"compare", "console_logs", "console_logs_py", "pyrerun_scripts", "MetaData"}
REL_TOL, ABS_TOL, LARGE_TOL = 1e-6, 1e-9, 1e-3
CHUNK = 5_000_000
os.makedirs(OUT, exist_ok=True)

parts = sorted(d for d in os.listdir(NEW)
               if os.path.isdir(os.path.join(NEW, d)) and d not in SKIP
               and os.path.isdir(os.path.join(REF, d)))

# ---- 1) Tables, with compare_rerun.py ----------------------------------------
print("Part 1: tables (compare_rerun.py), folders: " + ", ".join(parts), flush=True)
env = dict(os.environ, WIN_ROOT=REF, SEAAD_ROOT=NEW, COMPARE_CELLS="1")
rc = subprocess.run([sys.executable, os.path.join(NEW, "compare_rerun.py"), *parts], env=env).returncode
if rc != 0:
    sys.exit("compare_rerun.py stopped with exit code %d" % rc)
tag = re.sub(r"[^A-Za-z0-9_.-]", "_", "_".join(parts))[:80]
for kind, ext in (("log", "txt"), ("detail", "csv"), ("cells", "csv")):
    src = os.path.join(OUT, "compare_%s_%s.%s" % (tag, kind, ext))
    dst = os.path.join(OUT, "compare_python_tables_%s.%s" % (kind, ext))
    if os.path.exists(src):
        os.replace(src, dst)
p_log = os.path.join(OUT, "compare_python_tables_log.txt")
txt = open(p_log, encoding="utf-8").read()
txt = txt.replace("The reference results were produced on Windows; the rerun used the container on",
                  "The reference results are the original Python runs on Fedora; the rerun used the container on")
txt = txt.replace("compare_%s_detail.csv" % tag, "compare_python_tables_detail.csv")
txt = txt.replace("compare_%s_cells.csv" % tag, "compare_python_tables_cells.csv")
open(p_log, "w", encoding="utf-8").write(txt)

# ---- 2) Matrix Market files and lists ----------------------------------------
LOG = os.path.join(OUT, "compare_python_matrices_log.txt")
DET = os.path.join(OUT, "compare_python_matrices_detail.csv")
_log = open(LOG, "w", encoding="utf-8")


def say(s=""):
    print(s, flush=True)
    _log.write(s + "\n")
    _log.flush()


def opener(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p, "r")


def same_bytes(a, b):
    if os.path.getsize(a) != os.path.getsize(b):
        return False
    with open(a, "rb") as fa, open(b, "rb") as fb:
        while True:
            x, y = fa.read(1 << 24), fb.read(1 << 24)
            if x != y:
                return False
            if not x:
                return True


def mtx_header(p):
    n = 0
    with opener(p) as fh:
        for line in fh:
            n += 1
            if not line.startswith("%"):
                return n, tuple(int(v) for v in line.split())
    return n, None


def compare_mtx(a, b):
    r = dict(n_entries=0, n_index_diff=0, n_val_diff=0, n_outside=0, max_rel=0.0)
    ha, da = mtx_header(a)
    hb, db = mtx_header(b)
    r["dims_rerun"], r["dims_reference"] = str(da), str(db)
    if da != db:
        r["status"] = "DIFFERENT SIZE"
        return r
    kw = dict(sep=" ", header=None, comment=None, chunksize=CHUNK, engine="c")
    ra = pd.read_csv(a, skiprows=ha, **kw)
    rb = pd.read_csv(b, skiprows=hb, **kw)
    for ca, cb in zip(ra, rb):
        r["n_entries"] += len(ca)
        ia, ib = ca.iloc[:, :2].to_numpy(), cb.iloc[:, :2].to_numpy()
        r["n_index_diff"] += int((ia != ib).any(axis=1).sum())
        if ca.shape[1] > 2:
            va = ca.iloc[:, 2].to_numpy(dtype=float)
            vb = cb.iloc[:, 2].to_numpy(dtype=float)
            d = np.abs(va - vb)
            m = np.maximum(np.abs(va), np.abs(vb))
            rel = np.divide(d, m, out=np.zeros_like(d), where=m > 0)
            bad = (d > ABS_TOL) & (rel > REL_TOL)
            r["n_val_diff"] += int(bad.sum())
            r["n_outside"] += int((bad & (rel > LARGE_TOL)).sum())
            if bad.any():
                r["max_rel"] = max(r["max_rel"], float(rel[bad].max()))
    if r["n_index_diff"]:
        r["status"] = "ENTRIES DIFFER"
    elif r["n_outside"]:
        r["status"] = "OUTSIDE 0.1%"
    elif r["n_val_diff"]:
        r["status"] = "WITHIN 0.1%"
    else:
        r["status"] = "SAME"
    return r


def compare_list(a, b):
    with opener(a) as fa, opener(b) as fb:
        la, lb = fa.read().splitlines(), fb.read().splitlines()
    r = dict(n_entries=len(la), dims_rerun=str(len(la)), dims_reference=str(len(lb)))
    if len(la) != len(lb):
        r["status"] = "DIFFERENT SIZE"
    elif la == lb:
        r["status"] = "SAME"
    elif sorted(la) == sorted(lb):
        r["status"] = "SAME VALUES, OTHER ORDER"
    else:
        r["n_text_diff"] = sum(x != y for x, y in zip(la, lb))
        r["status"] = "TEXT DIFFERS"
    return r


t0 = time.time()
say("compare_python.py, part 2 (matrices and lists)   %s" % time.strftime("%Y-%m-%d %H:%M:%S"))
say("Rerun    : %s (container)" % NEW)
say("Reference: %s (original Python runs on Fedora)" % REF)
say("Rules    : bytes first; if they differ, indices exactly and values within "
    "%.0e relative (%.0e absolute)" % (REL_TOL, ABS_TOL))
say("")
EXT_MTX = (".mtx", ".mtx.gz")
EXT_LIST = (".tsv", ".tsv.gz")
rows = []
for part in parts:
    pn, pr = os.path.join(NEW, part), os.path.join(REF, part)
    fn = {os.path.relpath(os.path.join(d, f), pn) for d, _, fs in os.walk(pn) for f in fs
          if f.endswith(EXT_MTX + EXT_LIST)}
    fr = {os.path.relpath(os.path.join(d, f), pr) for d, _, fs in os.walk(pr) for f in fs
          if f.endswith(EXT_MTX + EXT_LIST)}
    if not fn and not fr:
        continue
    say("== %s" % part)
    for f in sorted(fn | fr):
        if f not in fr:
            r = dict(status="ONLY IN RERUN")
        elif f not in fn:
            r = dict(status="ONLY IN REFERENCE")
        elif same_bytes(os.path.join(pn, f), os.path.join(pr, f)):
            r = dict(status="SAME (IDENTICAL BYTES)")
        elif f.endswith(EXT_MTX):
            r = compare_mtx(os.path.join(pn, f), os.path.join(pr, f))
        else:
            r = compare_list(os.path.join(pn, f), os.path.join(pr, f))
        r.update(part=part, file=f)
        rows.append(r)
        extra = ""
        if r.get("n_entries"):
            extra = " | %d entries" % r["n_entries"]
            if r.get("n_val_diff"):
                extra += ", %d values beyond 1e-6, %d beyond 0.1%%, max rel %.2e" % (
                    r["n_val_diff"], r["n_outside"], r["max_rel"])
            if r.get("n_index_diff"):
                extra += ", %d entries at other positions" % r["n_index_diff"]
        if r["status"] == "DIFFERENT SIZE":
            extra = " | rerun %s, reference %s" % (r["dims_rerun"], r["dims_reference"])
        say("   [%s] %s%s" % (r["status"], f, extra))
    say("")

pd.DataFrame(rows).to_csv(DET, index=False)
say("== Overall")
for k, v in pd.Series([r["status"] for r in rows]).value_counts().items():
    say("   %-26s %d" % (k, v))
say("")
say("Tables: %s" % p_log)
say("Detail: %s" % DET)
say("Done in %.1f min" % ((time.time() - t0) / 60))
