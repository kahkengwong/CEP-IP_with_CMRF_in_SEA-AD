#!/usr/bin/env python3
"""
build_SuppTable_S6.py

Builds the data sheets of Supplementary Table S6 (DecontX ambient-RNA
removal) from the three tables written by Part7_DecontX.py. Sheet A
summarises the per-nucleus contamination estimate overall, by supertype, by
brain region and, for Micro-PVM_3-SEAAD, by region set (n, mean, median,
quartiles, 95th percentile and maximum; quantiles are linear interpolation,
the same as R quantile type 7). Sheet B is the per-gene table and sheet D the
module diagnostic, both copied as written by Part 7. Sheet C summarises the
retained fraction of sheet B for the module genes, the other mitochondrial
genes and the non-mitochondrial reference genes. Sheet E lists every nucleus
with its supertype, library, brain region, region set and contamination
estimate, in the order of the Part 7 file. Numbers are written with 15
significant digits, as openxlsx wrote them in the published table, so the
values match the published sheets. finalize_supp_tables.py adds the Contents
sheet and saves the published file.

Usage:
    python build_SuppTable_S6.py
    python build_SuppTable_S6.py --root /data/SEA-AD
Root: --root, else $SEAAD_ROOT, else the current folder.

Inputs (all written by scripts/Python/Part7_DecontX.py):
    <root>/Part7_DecontX/decontX_contamination_all.csv
        one row per nucleus: exp_component_name, Supertype, library_prep,
        Brain Region, arm, decontX_contamination
    <root>/Part7_DecontX/decontX_mito_contamination.csv
        per gene in Micro-PVM_3-SEAAD: mean before and after, retained and
        removed fraction
    <root>/Part7_DecontX/decontX_module_diagnostic.csv
        per supertype and per region set of Micro-PVM_3-SEAAD: module rs
        before and after decontamination
Output:
    <root>/Manuscript/Supp_Tables/build/SuppTable_S6_DecontX.xlsx
        sheets A_contamination_summary, B_gene_level, C_gene_set_summary,
        D_module_diagnostic, E_per_nucleus
    <root>/Manuscript/Supp_Tables/build/build_SuppTable_S6_log.txt
        checks and row counts, ending with VERDICT: PASS or VERDICT: FAIL.
        On FAIL no workbook is written.

Formatting: Calibri 11 (openpyxl default), bold header row, frozen top row,
no fills, no borders.
"""

import argparse
import datetime as dt
import math
import os
import sys

import numpy as np
import pandas as pd
from openpyxl import Workbook
from openpyxl.cell import WriteOnlyCell
from openpyxl.styles import Font

TARGET = "Micro-PVM_3-SEAAD"
CT = "decontX_contamination"
STYPE = "Supertype"
PER_NUCLEUS_COLS = ["exp_component_name", "Supertype", "library_prep",
                    "Brain Region", "arm", "decontX_contamination"]
TEXT_COLS = PER_NUCLEUS_COLS[:-1]

# Reference values used only in checks.
N_NUCLEI = 240651                 # nuclei in the 2026-06-22 immune h5ad
PVM3_DELTA, PVM3_DELTA_TOL = -0.0155, 5e-4
# Values printed by the Part 7 run of 2026-07-30; shown in the log, not checked.
LOG_REF = {"n": 240651, "mean": 0.1947, "median": 0.1405, "q25": 0.0345,
           "q75": 0.2945, "q95": 0.5971, "max": 1.0000}


##### Log #####
class Log:
    def __init__(self):
        self.lines, self.fails = [], 0

    def say(self, s=""):
        self.lines.append(str(s))
        print(s, flush=True)

    def chk(self, ok, what):
        self.say(f"  [{'PASS' if ok else 'FAIL'}] {what}")
        if not ok:
            self.fails += 1
        return bool(ok)

    def write(self, path):
        self.say("")
        self.say(f"VERDICT: {'PASS' if self.fails == 0 else 'FAIL'}")
        with open(path, "w", encoding="utf-8") as f:
            f.write("\n".join(self.lines) + "\n")


##### Helpers #####
def cell_value(v):
    """The value as openxlsx wrote it: missing values blank, numbers with 15
    significant digits (whole numbers as integers), booleans as booleans."""
    if v is None:
        return None
    if isinstance(v, (bool, np.bool_)):
        return bool(v)
    if isinstance(v, (int, np.integer)):
        return int(v)
    if isinstance(v, (float, np.floating)):
        v = float(v)
        if not math.isfinite(v):
            return None
        v = float(f"{v:.15g}")
        return int(v) if v.is_integer() and abs(v) < 1e15 else v
    if isinstance(v, str):
        return v
    if pd.isna(v):
        return None
    return v


def summarise(df, label_cols):
    """n, mean and quantiles of contamination, by the given grouping (sorted
    group labels, missing labels kept as their own group)."""
    if not label_cols:
        v = df[CT]
        return pd.DataFrame([{
            "n_cells": int(len(v)), "mean": v.mean(), "median": v.median(),
            "q25": v.quantile(0.25), "q75": v.quantile(0.75),
            "q95": v.quantile(0.95), "max": v.max(),
        }])
    g = df.groupby(label_cols, dropna=False)[CT]
    out = pd.DataFrame({
        "n_cells": g.size(),
        "mean": g.mean(),
        "median": g.median(),
        "q25": g.quantile(0.25),
        "q75": g.quantile(0.75),
        "q95": g.quantile(0.95),
        "max": g.max(),
    })
    return out.reset_index()


def contamination_summary(cells):
    overall = summarise(cells, [])
    overall.insert(0, "group", "all cells")
    overall.insert(1, "grouping", "all cells")
    by_stype = summarise(cells, [STYPE]).rename(columns={STYPE: "group"})
    by_stype.insert(1, "grouping", "supertype")
    by_rg = summarise(cells, ["Brain Region"]).rename(columns={"Brain Region": "group"})
    by_rg.insert(1, "grouping", "brain region")
    tgt = cells[cells[STYPE] == TARGET]
    by_arm = summarise(tgt, ["arm"]).rename(columns={"arm": "group"})
    by_arm.insert(1, "grouping", "%s by region set" % TARGET)
    return pd.concat([overall, by_stype, by_rg, by_arm], ignore_index=True)


def gene_set_summary(mito):
    grp = []
    in_mod = mito["in_module"].astype(bool)
    is_mito = mito["is_mito"].astype(bool)
    for label, mask in [("Module genes (9)", in_mod),
                        ("Mitochondrial, not in module", is_mito & ~in_mod),
                        ("Non-mitochondrial reference", ~is_mito)]:
        v = mito.loc[mask, "fraction_retained"].dropna()
        if len(v):
            grp.append({"gene_set": label, "n_genes": len(v),
                        "mean_fraction_retained": v.mean(),
                        "min_fraction_retained": v.min(),
                        "max_fraction_retained": v.max()})
    return pd.DataFrame(grp, columns=["gene_set", "n_genes", "mean_fraction_retained",
                                      "min_fraction_retained", "max_fraction_retained"])


def rows_of(df):
    """Data rows of a frame as lists of written values."""
    cols = [df[c].tolist() for c in df.columns]
    return [[cell_value(v) for v in rec] for rec in zip(*cols)]


def per_nucleus_rows(cells):
    cols = [cells[c].tolist() for c in TEXT_COLS]
    vals = [cell_value(v) for v in cells[CT].to_numpy(dtype=float)]
    cols = [[None if (v is None or (isinstance(v, float) and math.isnan(v))) else v
             for v in c] for c in cols]
    return [list(rec) for rec in zip(*cols, vals)]


def has_dash(rows, header):
    bad = ("\u2013", "\u2014")
    for v in header:
        if isinstance(v, str) and any(b in v for b in bad):
            return True
    for r in rows:
        for v in r:
            if isinstance(v, str) and any(b in v for b in bad):
                return True
    return False


##### Main #####
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[1])
    ap.add_argument("--root", default=os.environ.get("SEAAD_ROOT", os.getcwd()))
    a = ap.parse_args()
    root = os.path.abspath(a.root)
    p7 = os.path.join(root, "Part7_DecontX")
    out_dir = os.path.join(root, "Manuscript", "Supp_Tables", "build")
    os.makedirs(out_dir, exist_ok=True)
    out = os.path.join(out_dir, "SuppTable_S6_DecontX.xlsx")
    log_path = os.path.join(out_dir, "build_SuppTable_S6_log.txt")
    f_cells = os.path.join(p7, "decontX_contamination_all.csv")
    f_mito = os.path.join(p7, "decontX_mito_contamination.csv")
    f_diag = os.path.join(p7, "decontX_module_diagnostic.csv")

    L = Log()
    L.say("build_SuppTable_S6.py")
    L.say(f"started {dt.datetime.now().isoformat(timespec='seconds')}")
    L.say(f"root    {root}")
    L.say("")
    L.say("Inputs")
    ok = all([L.chk(os.path.isfile(f), f"found {os.path.relpath(f, root)}")
              for f in (f_cells, f_mito, f_diag)])
    if not ok:
        L.write(log_path)
        sys.exit(1)

    # ---- per-nucleus file
    # Sheet E: numbers parsed exactly (round_trip), as data.table::fread reads them.
    cells = pd.read_csv(f_cells, dtype={c: str for c in TEXT_COLS},
                        float_precision="round_trip")
    L.say("")
    L.say("Per-nucleus file")
    L.say(f"  {len(cells):,} rows, columns: {', '.join(map(str, cells.columns))}")
    if not L.chk(all(c in cells.columns for c in PER_NUCLEUS_COLS),
                 "has the six columns " + ", ".join(PER_NUCLEUS_COLS)):
        L.write(log_path)
        sys.exit(1)
    cells = cells[PER_NUCLEUS_COLS]
    cells[CT] = pd.to_numeric(cells[CT], errors="coerce")
    n_bad = int((~np.isfinite(cells[CT].to_numpy(dtype=float))).sum())
    L.chk(n_bad == 0, f"every contamination value is a finite number ({n_bad} are not)")
    L.chk(len(cells) == N_NUCLEI, f"{len(cells):,} nuclei (expected {N_NUCLEI:,})")
    L.chk(not cells["exp_component_name"].duplicated().any(), "each nucleus listed once")
    L.chk(cells["arm"].isin(["NeoC", "NonNeoC"]).all(), "arm is NeoC or NonNeoC")
    valid = cells[np.isfinite(cells[CT].to_numpy(dtype=float))]

    # ---- sheets A to D
    # Sheet A: the published values were computed from the file read with the
    # default pandas parser, which keeps only the first 16 decimals of a long
    # number such as 0.000232978071498691 (relative error below 1e-12). The
    # file is read again that way so that sheet A matches the published one.
    cells_a = pd.read_csv(f_cells, usecols=[STYPE, "Brain Region", "arm", CT],
                          dtype={STYPE: str, "Brain Region": str, "arm": str})
    cells_a[CT] = pd.to_numeric(cells_a[CT], errors="coerce")
    cells_a = cells_a[cells_a[CT].notna()]
    A = contamination_summary(cells_a)
    mito = pd.read_csv(f_mito)
    diag = pd.read_csv(f_diag)
    C = gene_set_summary(mito)

    L.say("")
    L.say("Sheets A to D")
    need_b = ["gene", "in_module", "is_mito", "mean_raw", "mean_decontX",
              "fraction_retained", "fraction_removed"]
    need_d = ["supertype", "target", "n_GOIpos", "rho_before", "rho_after", "delta"]
    L.chk(list(mito.columns) == need_b, "gene table columns: " + ", ".join(need_b))
    L.chk(list(diag.columns) == need_d, "module diagnostic columns: " + ", ".join(need_d))
    o = A[A["grouping"] == "all cells"].iloc[0]
    L.chk(int(o["n_cells"]) == N_NUCLEI, f"sheet A covers {N_NUCLEI:,} nuclei")
    L.say(f"  {'':8s} {'this run':>12s} {'Part 7 log':>12s}")
    L.say(f"  {'n':8s} {int(o['n_cells']):12d} {LOG_REF['n']:12d}")
    for k in ("mean", "median", "q25", "q75", "q95", "max"):
        L.say(f"  {k:8s} {float(o[k]):12.4f} {LOG_REF[k]:12.4f}")
    if "supertype" in diag.columns and "delta" in diag.columns:
        dsel = diag.loc[diag["supertype"] == TARGET, "delta"]
        dval = float(dsel.iloc[0]) if len(dsel) else float("nan")
        L.chk(len(dsel) == 1 and abs(dval - PVM3_DELTA) < PVM3_DELTA_TOL,
              f"module rs change for {TARGET} is {dval:+.4f} (expected {PVM3_DELTA:+.4f} "
              f"within {PVM3_DELTA_TOL})")
    L.chk(len(C) == 3, f"sheet C has the three gene sets ({len(C)} rows)")

    sheets = [("A_contamination_summary", list(A.columns), rows_of(A)),
              ("B_gene_level", list(mito.columns), rows_of(mito)),
              ("C_gene_set_summary", list(C.columns), rows_of(C)),
              ("D_module_diagnostic", list(diag.columns), rows_of(diag))]

    # ---- sheet E, and sheet A recomputed from it
    E_rows = per_nucleus_rows(cells)
    sheets.append(("E_per_nucleus", list(PER_NUCLEUS_COLS), E_rows))
    L.say("")
    L.say("Sheet A against the per-nucleus values")
    A_rows = sheets[0][2]
    stats = ["mean", "median", "q25", "q75", "q95", "max"]
    x = valid[CT].to_numpy(dtype=float)
    n_ok, n_rows = 0, 0
    for r in A_rows:
        rec = dict(zip(A.columns, r))
        if rec["grouping"] == "all cells":
            v = x
        elif rec["grouping"] == "supertype":
            v = x[(valid[STYPE] == rec["group"]).to_numpy()]
        elif rec["grouping"] == "brain region":
            v = x[(valid["Brain Region"] == rec["group"]).to_numpy()]
        else:
            v = x[((valid[STYPE] == TARGET) & (valid["arm"] == rec["group"])).to_numpy()]
        s = [v.mean(), np.median(v), np.quantile(v, 0.25), np.quantile(v, 0.75),
             np.quantile(v, 0.95), v.max()] if len(v) else [np.nan] * 6
        good = rec["n_cells"] == len(v) and all(
            rec[k] is not None and abs(rec[k] - sv) < 1e-9 for k, sv in zip(stats, s))
        n_rows += 1
        n_ok += good
        if not good:
            L.chk(False, f"sheet A row {rec['grouping']} / {rec['group']} reproduced")
    L.chk(n_ok == n_rows, f"all {n_rows} rows of sheet A reproduced from sheet E (difference < 1e-9)")
    L.chk(not any(has_dash(r, h) for _, h, r in sheets), "no en or em dash in any sheet")

    L.say("")
    L.say("Rows per sheet")
    for name, _, rows in sheets:
        L.say(f"  {name:26s} {len(rows):,}")

    if L.fails:
        L.say("")
        L.say(f"{L.fails} check(s) failed; no workbook written.")
        L.write(log_path)
        sys.exit(1)

    # ---- write
    wb = Workbook(write_only=True)
    bold = Font(bold=True)
    for name, hdr, rows in sheets:
        ws = wb.create_sheet(name)
        ws.freeze_panes = "A2"
        hcells = []
        for h in hdr:
            c = WriteOnlyCell(ws, value=str(h))
            c.font = bold
            hcells.append(c)
        ws.append(hcells)
        for r in rows:
            ws.append(r)
    wb.save(out)
    L.say("")
    L.say(f"saved {os.path.relpath(out, root)} ({os.path.getsize(out) / 1e6:.1f} MB)")
    L.write(log_path)


if __name__ == "__main__":
    main()
