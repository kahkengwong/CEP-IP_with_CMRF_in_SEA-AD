#!/usr/bin/env python3
# Compares the Part 8 gam_data tables of the Python rerun with the Windows copy
# used by the R analyses, where the two differ in their number of columns.
# For each table: which columns are only in one of the two, and whether every
# column found in both holds the same values (rows matched in file order;
# numbers within 1e-6 relative, text exactly). Read-only.
# Output: printed on screen and saved to ~/Desktop/SEA-AD_pyrerun/compare/compare_gamdata_log.txt

import glob
import os
import time

import numpy as np
import pandas as pd

NEW = os.path.expanduser("~/Desktop/SEA-AD_pyrerun/R_export_CP10K")
REF = os.path.expanduser("~/Desktop/SEA-AD_Windows/SEA-AD/R_export_CP10K")
OUT = os.path.expanduser("~/Desktop/SEA-AD_pyrerun/compare/compare_gamdata_log.txt")
REL, ABS = 1e-6, 1e-9

lines = [f"compare_gamdata.py   {time.strftime('%Y-%m-%d %H:%M:%S')}",
         "Rerun    : SEA-AD_pyrerun/R_export_CP10K (container)",
         "Reference: SEA-AD_Windows/SEA-AD/R_export_CP10K (copy used by the R analyses)", ""]


def out(t=""):
    print(t)
    lines.append(t)


n_ok = n_bad = 0
extra_sets = {}
for pn in sorted(glob.glob(os.path.join(NEW, "*", "*", "gam_data", "*_selected.csv"))):
    rel = os.path.relpath(pn, NEW)
    pr = os.path.join(REF, rel)
    if not os.path.exists(pr):
        out(f"   [ONLY IN RERUN] {rel}")
        continue
    a = pd.read_csv(pn, dtype=str, keep_default_na=False)
    b = pd.read_csv(pr, dtype=str, keep_default_na=False)
    only_new = [c for c in a.columns if c not in b.columns]
    only_ref = [c for c in b.columns if c not in a.columns]
    extra_sets.setdefault(tuple(only_ref), []).append(rel)
    bad = []
    if len(a) != len(b):
        bad.append(f"rows {len(a)} vs {len(b)}")
    else:
        for c in [c for c in a.columns if c in b.columns]:
            x = pd.to_numeric(a[c], errors="coerce").to_numpy()
            y = pd.to_numeric(b[c], errors="coerce").to_numpy()
            num = ~np.isnan(x) & ~np.isnan(y)
            if num.any():
                d = np.abs(x[num] - y[num])
                m = np.maximum(np.abs(x[num]), np.abs(y[num]))
                diff = (d > ABS) & (d > REL * m)
                txt = (a[c][~num].to_numpy() != b[c][~num].to_numpy())
                k = int(diff.sum()) + int(txt.sum())
            else:
                k = int((a[c].to_numpy() != b[c].to_numpy()).sum())
            if k:
                bad.append(f"{c}: {k} cells differ")
    status = "SAME IN SHARED COLUMNS" if not bad and not only_new else "DIFFERS"
    n_ok += status.startswith("SAME")
    n_bad += not status.startswith("SAME")
    out(f"   [{status}] {rel} | {a.shape[1]} vs {b.shape[1]} columns"
        + (f" | only in rerun: {only_new}" if only_new else "")
        + (f" | {'; '.join(bad)}" if bad else ""))

out("")
out("Columns only in the reference, by group of tables:")
for cols, files in extra_sets.items():
    out(f"   {len(files)} tables: {list(cols)}")
out("")
out(f"Tables with every shared column the same: {n_ok}; with differences: {n_bad}")
os.makedirs(os.path.dirname(OUT), exist_ok=True)
open(OUT, "w", encoding="utf-8").write("\n".join(lines) + "\n")
print(f"\nSaved: {OUT}")
