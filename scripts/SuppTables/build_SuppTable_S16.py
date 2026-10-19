#!/usr/bin/env python3
"""
build_SuppTable_S16.py

Builds SuppTable_Fig18.xlsx (Figure 18 and Supplementary Figures 14 and 15 of
the Part 12 rerun) from the stats and plot-data files written by
Part12_Fig18_Analysis.r (sections 22 to 24), the figure log written by
Part12_Fig18_Plots.r (section 25), the config (Part12_01_config.json), the unit
table (05_units.csv) with the step 05 log, and the 16A stats file with its log.
finalize_supp_tables.py turns this workbook into Supplementary Table S16.

Usage:
    python build_SuppTable_S16.py
    python build_SuppTable_S16.py --root /data/Part12_rerun
    python build_SuppTable_S16.py --root <folder> --out <folder>/SuppTable_Fig18.xlsx

Without --root the script uses $SEAAD_ROOT/Part12_rerun, or Part12_rerun in the
current folder when SEAAD_ROOT is not set. --root is searched recursively; every
input file must occur exactly once under it, so the script works on the
Part12_rerun folder tree or on a flat folder holding copies of the files.
Without --out the workbook goes to <root>/10_Fig18/SuppTable_Fig18.xlsx (the
folder is created if it does not exist). A build log is written beside the
workbook.

Numbers are parsed and averaged exactly as R's read.csv and mean() do.

The script computes no new statistic. Before it writes anything it checks:
  B-01  the three section logs end with VERDICT: PASS and 0 FAIL
  B-02  every input file is listed in the outputs of the log of the section that wrote it
  B-03  the three analysis logs name one script, one config version and one labels version
  B-04  every stats and plot-data row carries the run ID of its log, and every
        stats row the config and labels version of its log
  B-05  row counts equal the counts the logs report, and the plot-data files have
        the rows the design gives (arms x region sets x panels; arms x effects x
        donor grid; the region rows of 18E)
  B-06  every fdr equals BH recomputed within its family
  B-07  Plot_18A_verdicts agrees with Stats_18A (coefficients, VIF, FDRs, donors),
        and every verdict follows from them
  B-08  the 18A verdict tallies (M4 primary rows per panel, and every arm's rows)
        equal the section 22 log
  B-09  every 18B AUC, best cut, Youden J and kappa recomputed from the Fig18BC_cells
        values equals its stats row, with the same cell and donor counts
  B-10  every 18C donor AUC and cell count recomputed from the Fig18BC_cells values
        equals Plot_18C_donors, and every 18C median and count follows from them
  B-11  every 18D crossing and rate at the largest donor count follows from the
        curves; each arm's donor share recomputed from the unit table equals the 16A
        estimate, and every curve's target share is that share times the effect
  B-12  every 18E retention, unit count and donor share recomputed from the unit
        table equals its plot-data row, and every 18E stats row equals the plot data
  B-13  the figure log is PASS and was started after the last analysis section
  B-14  no text written to the workbook contains an en dash or an em dash
  B-15  the config exists and has the config version of the section logs, the unit
        table and Stats_16A are listed in the outputs of their logs, which are PASS,
        and the unit table carries the labels version of the section logs
Any FAIL stops the build and no workbook is written.

Not recomputed, because they depend on inputs this workbook does not carry or on
R's random streams: the 18A coefficients (the donor panel scores are computed from
the step 02 count matrix), the 18C bootstrap intervals and p values, the 18D
detection rates and the 18E split-half rs.

Formatting: Arial throughout, bold header row, frozen top row. No fills, no
borders, no colours. The verdict column of Fig18A_verdicts, the band column of
Fig18B_readout and the reading column of Fig18C_readout are Excel formulas; the
file is flagged to recalculate on open.
"""

import argparse
import datetime as dt
import json
import os
import pathlib
import re
import sys

import numpy as np
import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Alignment, Font
from openpyxl.utils import get_column_letter

# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------
SECTIONS = [
    ("22", "18A_NuclearPanels", ["Stats_18A.csv", "Plot_18A_verdicts.csv", "Plot_18A_panels.csv"]),
    ("23", "18BC_MTOverlap", ["Stats_18B.csv", "Stats_18C.csv", "Plot_18B_cells.csv", "Plot_18C_donors.csv"]),
    ("24", "18DE_DonorPower", ["Stats_18D.csv", "Stats_18E.csv", "Plot_18D_curves.csv", "Plot_18E_regions.csv"]),
]
FIG_LOG = "Part12_Fig18_log.txt"
CONFIG = "Part12_01_config.json"
UNITS = "05_units.csv"
UNITS_LOG = "*05_Units*_log.txt"
S16A = ("Stats_16A.csv", "Part12_16AB_DonorStructure_log.txt")
ARMS = ["M4", "PVM_3", "PVM_1", "HKG8"]
PROV_COLS = ["run_id", "config_version", "labels_version"]
TOL = 1e-9
DASHES = ("\u2013", "\u2014")
VK18 = ["second only", "Both", "NE panel only", "Neither", "Unreadable", "Untestable"]
CHAR_COLS = ("panel", "arm", "level", "quantity", "measure", "family", "run_id", "labels_version",
             "donor_key", "donor_id", "dataset", "region", "predictor", "outcome", "class", "fitted",
             "verdict", "is_module", "libs")

FONT = Font(name="Arial", size=10)
FONT_B = Font(name="Arial", size=10, bold=True)
FMT_P = "0.00E+00"


# ---------------------------------------------------------------------------
# Check ledger and shared helpers (as build_SuppTable_S15.py)
# ---------------------------------------------------------------------------
class Ledger:
    def __init__(self):
        self.lines, self.n = [], {"PASS": 0, "FAIL": 0, "NOTE": 0}

    def chk(self, cid, ok, what, detail=""):
        tag = "PASS" if ok else "FAIL"
        self.n[tag] += 1
        self._say(f"[{tag}] {cid} {what} | {detail}")

    def note(self, what, detail=""):
        self.n["NOTE"] += 1
        self._say(f"[NOTE] {what} | {detail}")

    def _say(self, s):
        self.lines.append(s)
        print(s)


def find(root, name):
    hits = [p for p in root.rglob(name) if p.is_file()]
    if len(hits) != 1:
        sys.exit(f"ERROR: expected exactly one {name} under {root}, found {len(hits)}")
    return hits[0]


_NUM = re.compile(r"([+-]?)(\d*)(?:\.(\d*))?(?:[eE]([+-]?\d+))?$")


def r_strtod(s):
    """One number as R's read.csv parses it (R_strtod5): the digits are
    accumulated in long double, scaled by a power of ten built by repeated
    squaring in long double, then rounded to double. pandas' own parser and
    Python's float() both differ from read.csv in the last bit for some values,
    which changes the order of values that nearly tie after centring."""
    s = s.strip()
    if s in ("", "NA"):
        return np.nan
    low = s.lower()
    if low in ("inf", "+inf"):
        return np.inf
    if low == "-inf":
        return -np.inf
    if low == "nan":
        return np.nan
    m = _NUM.match(s)
    if not m or not (m.group(2) or m.group(3)):
        raise ValueError(s)
    sign, ip, fp, ex = m.group(1), m.group(2) or "", m.group(3) or "", m.group(4)
    digits = ip + fp
    if len(digits.lstrip("0")) >= 18:
        return float(s)
    expn = -len(fp) + (int(ex) if ex else 0)
    ld = np.longdouble
    ans = ld(0)
    for ch in digits:
        ans = ld(10) * ans + ld(ord(ch) - 48)
    if expn != 0 and ans != 0:
        n, p10, fac = abs(expn), ld(10), ld(1)
        while n:
            if n & 1:
                fac *= p10
            n >>= 1
            p10 *= p10
        ans = ans / fac if expn < 0 else ans * fac
    v = float(ans)
    return -v if sign == "-" else v


def r_mean(v):
    """R's mean() of a double vector: long double sum in order, divided by n,
    plus the long double correction pass, rounded to double."""
    ld = np.longdouble
    n = len(v)
    s = ld(0)
    for a in v:
        s += ld(a)
    s /= n
    if np.isfinite(float(s)):
        t = ld(0)
        for a in v:
            t += ld(a) - s
        s += t / n
    return float(s)


def read_csv_like_r(path, char_cols=()):
    """Read a CSV with R read.csv missing-value rules, so complete-case filters
    match the analysis: "NA" is missing in every column; an empty field is
    missing in a numeric column and an empty string in a character column."""
    d = pd.read_csv(path, dtype=str, keep_default_na=False, na_values=[])
    for c in d.columns:
        v = d[c].where(d[c] != "NA")
        if c in char_cols:
            d[c] = v
            continue
        try:
            d[c] = pd.Series([np.nan if x is None or (isinstance(x, float) and np.isnan(x)) else r_strtod(x)
                              for x in v], index=d.index, dtype=float)
        except ValueError:
            d[c] = v
    return d


def as_bool(s):
    return s.map({"TRUE": True, "FALSE": False})


# ---------------------------------------------------------------------------
# Recomputation (section 19 and section 20 rules)
def ranks(x):
    return pd.Series(np.asarray(x, float)).rank().values


def parse_log(path):
    """Reads a Part 12 section log. Logs written by older copies of the R
    scripts add a field after the file size and after the config version; the
    patterns stop before it, so both forms give the same result."""
    txt = path.read_text(encoding="utf-8", errors="replace").replace("\r\n", "\n")
    g = lambda pat: (re.search(pat, txt, re.M).group(1) if re.search(pat, txt, re.M) else None)
    outs = set(re.findall(r"^\s+\S*/([^/\s|]+) \| [\d.]+ MB", txt, re.M))
    inputs = {lab.strip(): p for lab, p in
              re.findall(r"^\s+input (.+?)\s{2,}(\S.*?) \| [\d.]+ MB", txt, re.M)}
    return {
        "path": path,
        "script": g(r"^Script: (\S+)"),
        "run_id": g(r"^Run ID: (\S+)"),
        "start": g(r"^Start: (.+)$"),
        "config_version": g(r"^Config: .*? \| version (\S+)"),
        "labels_version": g(r"^Labels version: (\S+)"),
        "checks": g(r"^Checks: (.+?) \|"),
        "n_fail": int(g(r"^Checks: \d+ PASS, (\d+) FAIL") or -1),
        "verdict": g(r"^VERDICT: (\S+)"),
        "outputs": outs,
        "inputs": inputs,
        "rows_line": g(r"^\[PASS\] \S+ stats (?:and plot )?files? written(?: and read back)? \| (.+)$"),
        "text": txt,
    }


def bh(p):
    """Benjamini-Hochberg as R p.adjust(method = "BH"), NaN left as NaN."""
    p = np.asarray(p, dtype=float)
    out = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    m = ok.sum()
    if m == 0:
        return out
    pv = p[ok]
    order = np.argsort(pv)[::-1]
    rk = np.arange(m, 0, -1)
    q = np.minimum.accumulate(pv[order] * m / rk)
    res = np.empty(m)
    res[order] = np.minimum(q, 1.0)
    out[ok] = res
    return out


def close(a, b, tol=TOL):
    a, b = np.asarray(a, float), np.asarray(b, float)
    both_nan = np.isnan(a) & np.isnan(b)
    both_inf = np.isinf(a) & np.isinf(b) & (np.sign(a) == np.sign(b))
    with np.errstate(invalid="ignore"):
        return bool(np.all(both_nan | both_inf | (np.abs(a - b) <= tol)))


class Worst:
    """Largest difference over many compared values, and the count that differ."""
    def __init__(self):
        self.worst, self.at, self.n, self.bad = 0.0, "", 0, 0

    def add(self, v, ref, where, tol=TOL):
        self.n += 1
        v, ref = float(v), float(ref)
        if np.isnan(v) or np.isnan(ref) or np.isinf(v) or np.isinf(ref):
            if not close([v], [ref]):
                self.bad += 1
                self.at = self.at or where
            return
        if abs(v - ref) > self.worst:
            self.worst, self.at = abs(v - ref), where
        if abs(v - ref) >= tol:
            self.bad += 1

    def count(self, ok, where):
        self.n += 1
        if not ok:
            self.bad += 1
            self.at = self.at or where

    def detail(self):
        return (f"{self.n} values | {self.bad} differ | largest difference {self.worst:.1e}"
                + (f" ({self.at})" if self.at else ""))


def clean(v):
    """Excel-safe cell value: NaN and Inf to empty, numpy scalars to Python."""
    if v is None:
        return None
    if isinstance(v, (bool, np.bool_)):
        return bool(v)
    if isinstance(v, (float, np.floating)) and not np.isfinite(v):
        return None
    if isinstance(v, np.integer):
        return int(v)
    if isinstance(v, np.floating):
        return float(v)
    return v


def write_table(wb, name, header, rows, widths=None, p_cols=(), fmt_cols=None, wrap=False):
    ws = wb.create_sheet(name)
    ws.append(header)
    for r in rows:
        ws.append([clean(v) for v in r])
    for c in ws[1]:
        c.font = FONT_B
    for row in ws.iter_rows(min_row=2):
        for c in row:
            c.font = FONT
            if wrap:
                c.alignment = Alignment(wrap_text=True, vertical="top")
    fmt_cols = dict(fmt_cols or {})
    for col in p_cols:
        fmt_cols[col] = FMT_P
    for col, fmt in fmt_cols.items():
        j = header.index(col) + 1
        for i in range(2, ws.max_row + 1):
            ws.cell(i, j).number_format = fmt
    for j, h in enumerate(header, start=1):
        if widths and h in widths:
            w = widths[h]
        else:
            vals = [len(str(h))] + [len(str(ws.cell(i, j).value)) for i in range(2, min(ws.max_row, 200) + 1)
                                    if ws.cell(i, j).value is not None
                                    and not str(ws.cell(i, j).value).startswith("=")]
            w = min(max(vals) + 2, 40)
        ws.column_dimensions[get_column_letter(j)].width = w
    ws.freeze_panes = "A2"
    return ws


def stats_sheet(wb, name, df):
    """A stats file as a sheet, provenance columns left out."""
    cols = [c for c in df.columns if c not in PROV_COLS]
    write_table(wb, name, cols, df[cols].itertuples(index=False, name=None), p_cols=("p", "fdr"))



# ---------------------------------------------------------------------------
# Recomputation (section 23 and section 24 rules)
# ---------------------------------------------------------------------------
def auc(x, pos):
    """p12_auc(): the rank-sum AUC of x for membership, ties by midranks."""
    x = np.asarray(x, float); pos = np.asarray(pos, bool)
    n1, n0 = int(pos.sum()), int((~pos).sum())
    if n1 < 1 or n0 < 1:
        return np.nan
    r = ranks(x)
    return float((r[pos].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def best_cut(x, pos):
    """p12_best_cut(): the cut maximising Youden's J, cells above the cut
    called Module, and Cohen's kappa at that cut."""
    x = np.asarray(x, float); pos = np.asarray(pos, bool)
    n1, n0 = int(pos.sum()), int((~pos).sum())
    if n1 < 1 or n0 < 1:
        return {"cut": np.nan, "kappa": np.nan, "youden": np.nan}
    o = np.argsort(x, kind="stable"); xs = x[o]; ps = pos[o]
    tp = n1 - np.cumsum(ps); fp = n0 - np.cumsum(~ps)
    last = np.append(xs[1:] != xs[:-1], True)
    tp, fp, cuts = tp[last], fp[last], xs[last]
    j = tp / n1 + (1 - fp / n0) - 1
    i = int(np.argmax(j))
    a, b = tp[i], fp[i]; cc = n0 - b; d = n1 - a; n = n1 + n0
    acc = (a + cc) / n
    pe = ((a + b) / n) * ((a + d) / n) + ((cc + d) / n) * ((b + cc) / n)
    return {"cut": float(cuts[i]), "kappa": float((acc - pe) / (1 - pe)) if pe < 1 else np.nan,
            "youden": float(j[i])}


def share(y, d):
    """p12_share(): between-donor over total sum of squares, donor entered first."""
    y = np.asarray(y, float)
    m = r_mean(y)
    sst = float(np.sum((y - m) ** 2))
    if not sst > 0:
        return np.nan
    s = pd.Series(y).groupby(np.asarray(d), sort=False).agg(["sum", "size"])
    return float(((s["sum"] ** 2 / s["size"]).sum() - len(y) * m ** 2) / sst)


def verdict18(fdr_s, fdr_n, vif, fitted, alpha, vr):
    """p12_18_verdict() for one row."""
    if not fitted:
        return "Untestable"
    if not (np.isfinite(vif) and vif < vr):
        return "Unreadable"
    ss = bool(np.isfinite(fdr_s) and fdr_s < alpha)
    sn = bool(np.isfinite(fdr_n) and fdr_n < alpha)
    return "second only" if ss and not sn else "NE panel only" if sn and not ss else "Both" if ss and sn else "Neither"


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    default_root = pathlib.Path(os.environ.get("SEAAD_ROOT") or os.getcwd()) / "Part12_rerun"
    ap.add_argument("--root", default=str(default_root),
                    help="the Part12_rerun folder (default: $SEAAD_ROOT/Part12_rerun, else ./Part12_rerun)")
    ap.add_argument("--out", default=None, help="output xlsx (default: <root>/10_Fig18/SuppTable_Fig18.xlsx)")
    args = ap.parse_args()
    root = pathlib.Path(args.root).resolve()
    if not root.is_dir():
        sys.exit(f"ERROR: root folder not found: {root}")
    out = pathlib.Path(args.out).resolve() if args.out else root / "10_Fig18" / "SuppTable_Fig18.xlsx"
    out.parent.mkdir(parents=True, exist_ok=True)
    L = Ledger()
    L._say("build_SuppTable_S16.py")
    L._say(f"Start: {dt.datetime.now():%Y-%m-%d %H:%M:%S} | root {root}")
    L._say("")

    # ---- config ----------------------------------------------------------------------
    cfg_path = find(root, CONFIG)
    cfg = json.loads(cfg_path.read_text(encoding="utf-8"))
    F18 = cfg["figure18"]; PW = F18["power"]
    PK = [p["key"] for p in F18["panels"]]; PL = {p["key"]: p["label"] for p in F18["panels"]}
    OUT = {k: o["column"] for k, o in cfg["clinical"]["outcomes"].items()}
    OCL = {k: o["class"] for k, o in cfg["clinical"]["outcomes"].items()}
    DATASETS = list(cfg["datasets"].keys())
    EFF = [float(e) for e in PW["effects"]]; GRID = [int(n) for n in PW["n_grid"]]
    DETECT = float(PW["detect"])
    min_d = int(cfg["floors"]["donors_per_test"])
    minc = int(cfg["floors"]["cells_per_group_per_donor"])
    alpha = float(cfg["inference"]["alpha"]); vr = float(cfg["floors"]["vif_readable"])
    B_ = int(cfg["inference"]["bootstrap_B"])
    BANDS = [float(b) for b in F18["auc_bands"]]
    share_floor = int(PW["share_floor_donors"])

    # ---- logs and their outputs ----------------------------------------------------
    L._say("SECTION A -- logs and outputs")
    logs, files, D = {}, {}, {}
    for step, folder, names in SECTIONS:
        lg = parse_log(find(root, f"Part12_{folder}_log.txt"))
        logs[step] = lg
        L.chk(f"B-01.{step}", lg["verdict"] == "PASS" and lg["n_fail"] == 0,
              f"section {step} log is PASS with 0 FAIL", f"{lg['checks']} | {lg['verdict']}")
        for nm in names:
            p = find(root, nm)
            files[nm] = (step, p)
            L.chk(f"B-02.{nm}", nm in lg["outputs"], f"{nm} is listed in the outputs of the section {step} log",
                  lg["path"].name)
            D[nm] = read_csv_like_r(p, char_cols=CHAR_COLS)
    one = lambda k: sorted({str(lg[k]) for lg in logs.values()})
    L.chk("B-03", len(one("script")) == 1 and len(one("config_version")) == 1 and len(one("labels_version")) == 1,
          "one analysis script, config version and labels version across sections 22 to 24",
          f"{one('script')} | config {one('config_version')} | labels {one('labels_version')}")
    for nm, (step, _) in files.items():
        d, lg = D[nm], logs[step]
        ok = set(d["run_id"]) == {lg["run_id"]}
        if nm.startswith("Stats_"):
            ok = (ok and set(d["config_version"].astype(int).astype(str)) == {lg["config_version"]}
                  and set(d["labels_version"]) == {lg["labels_version"]})
        L.chk(f"B-04.{nm}", ok, f"{nm} rows carry the log's run ID"
              + (", config and labels version" if nm.startswith("Stats_") else ""),
              f"{lg['run_id']} | {lg['config_version']} | {lg['labels_version']}")

    S18A, V, PAN = D["Stats_18A.csv"], D["Plot_18A_verdicts.csv"], D["Plot_18A_panels.csv"]
    S18B, S18C, CELLS, DON = D["Stats_18B.csv"], D["Stats_18C.csv"], D["Plot_18B_cells.csv"], D["Plot_18C_donors.csv"]
    S18D, S18E, CUR, REG = D["Stats_18D.csv"], D["Stats_18E.csv"], D["Plot_18D_curves.csv"], D["Plot_18E_regions.csv"]
    V = V.rename(columns={"class": "oclass"})
    V["fitted"] = as_bool(V["fitted"]); CELLS["is_module"] = as_bool(CELLS["is_module"])

    # ---- config, unit table and 16A --------------------------------------------------
    L.chk("B-15.config", cfg_path.is_file() and [str(cfg.get("config_version"))] == one("config_version"),
          f"{CONFIG} exists and has the config version of the section logs",
          f"{cfg_path.name} | version {cfg.get('config_version')}")
    hits = [p for p in root.rglob(UNITS_LOG) if p.is_file()]
    if len(hits) != 1:
        sys.exit(f"ERROR: expected exactly one file matching {UNITS_LOG} under {root}, found {len(hits)}")
    up_log = parse_log(hits[0])
    up = find(root, UNITS)
    L.chk("B-15.units", UNITS in up_log["outputs"] and up_log["verdict"] == "PASS",
          f"{UNITS} is listed in the outputs of its log, which is PASS", up_log["path"].name)
    U = read_csv_like_r(up, char_cols=CHAR_COLS)
    L.chk("B-15.labels", set(U["labels_version"]) == set(one("labels_version")),
          "unit table carries the labels version of the section logs", ", ".join(sorted(set(U["labels_version"]))))
    a16_path, a16_log = find(root, S16A[0]), parse_log(find(root, S16A[1]))
    L.chk("B-15.16A", S16A[0] in a16_log["outputs"] and a16_log["verdict"] == "PASS",
          f"{S16A[0]} is listed in the outputs of its log, which is PASS", a16_log["path"].name)
    A16 = read_csv_like_r(a16_path, char_cols=CHAR_COLS)

    # ---- row counts -------------------------------------------------------------------
    rl = {k: logs[k]["rows_line"] or "" for k in logs}
    m22 = re.search(r"18A (\d+) rows \| (\d+) verdicts", rl["22"])
    m23 = re.search(r"18B (\d+) \| 18C (\d+) rows \| (\d+) donor AUCs", rl["23"])
    m24 = re.search(r"18D (\d+) \| 18E (\d+) rows", rl["24"])
    mc = re.search(r"labels carry this labels version \| (\d+) cells", logs["23"]["text"])
    ok = (bool(m22 and m23 and m24 and mc)
          and int(m22.group(1)) == len(S18A) and int(m22.group(2)) == len(V)
          and int(m23.group(1)) == len(S18B) and int(m23.group(2)) == len(S18C) and int(m23.group(3)) == len(DON)
          and int(m24.group(1)) == len(S18D) and int(m24.group(2)) == len(S18E) and int(mc.group(1)) == len(CELLS))
    L.chk("B-05.logs", ok, "stats, verdict, donor-AUC and cell row counts equal the logs",
          f"{rl['22']} || {rl['23']} || {rl['24']} || cells {len(CELLS)}")
    UM = U[U["arm"].isin(ARMS) & np.isfinite(U["module_post_pct"].astype(float))]
    ND = {a: UM[UM.arm == a]["donor_key"].nunique() for a in ARMS}
    n_cur = sum(len(EFF) * sum(1 for n in GRID if n <= ND[a]) for a in ARMS)
    n_reg = 0
    for a in ARMS:
        top = int(UM[UM.arm == a]["donor_key"].value_counts().max())
        n_reg += max(3, top) - 3 + 1
    design = {"Plot_18A_panels.csv": (len(PAN), len(ARMS) * len(DATASETS) * len(PK)),
              "Plot_18A_verdicts.csv": (len(V), len(ARMS) * len(DATASETS) * len(PK) * 2 * len(OUT)),
              "Plot_18D_curves.csv": (len(CUR), n_cur), "Plot_18E_regions.csv": (len(REG), n_reg),
              "Stats_18E.csv": (len(S18E), 3 * n_reg)}
    L.chk("B-05.design", all(a == b for a, b in design.values()), "plot-data row counts equal the design",
          " | ".join(f"{k} {a} of {b}" for k, (a, b) in design.items()))

    for nm in [n for n in D if n.startswith("Stats_")]:
        d = D[nm]
        ok = all(close(bh(g["p"].values), g["fdr"].values, 1e-12) for _, g in d.groupby("family", sort=False))
        L.chk(f"B-06.{nm}", ok, f"{nm} fdr equals BH within family", f"{d['family'].nunique()} families")
    L._say("")

    # ---- 18A ----------------------------------------------------------------------------
    L._say("SECTION B -- panel consistency")
    Si = S18A.set_index(["arm", "level", "measure", "quantity"])
    W = Worst()
    for r in V.itertuples():
        ms = f"{r.predictor}|{r.panel}|{r.outcome}"
        g = lambda q, col="estimate": Si.loc[(r.arm, r.dataset, ms, q), col]
        where = f"18A {r.arm} {r.dataset} {ms}"
        W.add(r.beta_second, g("beta_second"), where + " beta_second")
        W.add(r.beta_panel, g("beta_panel"), where + " beta_panel")
        W.add(r.vif, g("vif"), where + " vif")
        W.add(r.fdr_second, g("beta_second", "fdr"), where + " fdr_second")
        W.add(r.fdr_panel, g("beta_panel", "fdr"), where + " fdr_panel")
        W.add(r.n_donors, g("beta_second", "n_donors"), where + " n")
        W.count(r.fitted == bool(np.isfinite(r.beta_second) and np.isfinite(r.beta_panel)), where + " fitted")
        W.count(r.verdict == verdict18(r.fdr_second, r.fdr_panel, r.vif, r.fitted, alpha, vr), where + " verdict")
        W.count(r.oclass == OCL.get(r.outcome), where + " class")
    L.chk("B-07", W.bad == 0 and W.worst < TOL,
          "Plot_18A_verdicts agrees with Stats_18A and every verdict follows from it", W.detail())

    t22 = logs["22"]["text"]
    ok, n_lines = True, 0
    for sp in ("module", "mito"):
        blk = re.search(r"second predictor: " + sp + r"\n(.*?)(?:\n   second predictor:|\n  M4 )", t22, re.S)
        got = dict(re.findall(r"^\s+(\S+)\s+(second only \d+, Both \d+, NE panel only \d+, Neither \d+, "
                              r"Unreadable \d+, Untestable \d+)$", blk.group(1) if blk else "", re.M))
        for k in PK:
            g = V[(V.arm == "M4") & (V.predictor == sp) & (V.panel == k) & (V.oclass == "cognition")]
            mine = ", ".join(f"{v} {int((g.verdict == v).sum())}" for v in VK18)
            ok &= got.get(k) == mine; n_lines += 1
    for a in ARMS:
        for sp in ("module", "mito"):
            g = V[(V.arm == a) & (V.predictor == sp)]
            mine = ", ".join(f"{v} {int((g.verdict == v).sum())}" for v in VK18)
            mine += f" | max VIF {np.nanmax(g.vif.values):.2f}"
            ok &= bool(re.search(rf"^\s+{re.escape(a)}\s+{sp}\s+all rows: {re.escape(mine)}$", t22, re.M)); n_lines += 1
    L.chk("B-08", ok, "the 18A verdict tallies equal the section 22 log", f"{n_lines} tally lines")

    # ---- 18B and 18C, from the cells -------------------------------------------------
    CS = CELLS.copy()
    CS["_a"] = CS["arm"].map({a: i for i, a in enumerate(ARMS)})
    CS = CS[CS["_a"].notna()].sort_values("_a", kind="stable").reset_index(drop=True)
    CS = CS[["arm", "dataset", "donor_key", "is_module", "mito_pct"]]
    Bi = S18B.set_index(["arm", "level", "quantity"])
    W = Worst()
    for a in ARMS:
        for ds in DATASETS:
            d = CS[(CS.arm == a) & (CS.dataset == ds)]
            x, pos = d["mito_pct"].astype(float).values, d["is_module"].astype(bool).values
            bc = best_cut(x, pos)
            vals = {"auc": auc(x, pos), "kappa": bc["kappa"], "best_cut": bc["cut"], "youden": bc["youden"]}
            for q, v in vals.items():
                r = Bi.loc[(a, ds, q)]
                W.add(v, r["estimate"], f"18B {a} {ds} {q}")
                W.count(int(r["n_units"]) == len(d) and int(r["n_donors"]) == d["donor_key"].nunique(),
                        f"18B {a} {ds} {q} counts")
    L.chk("B-09", W.bad == 0 and W.worst < TOL,
          "every 18B AUC, best cut, Youden J and kappa recomputed from the Fig18BC_cells values equals its row",
          W.detail())

    W = Worst()
    Di = DON.set_index(["arm", "dataset", "donor_key"])
    mine = {}
    for a in ARMS:
        for ds in DATASETS:
            d = CS[(CS.arm == a) & (CS.dataset == ds)]
            for u in pd.unique(d["donor_key"]):
                z = d[d.donor_key == u]
                pos = z["is_module"].astype(bool).values
                n1, n0 = int(pos.sum()), int((~pos).sum())
                if n1 < minc or n0 < minc:
                    continue
                mine[(a, ds, u)] = (n1, n0, auc(z["mito_pct"].astype(float).values, pos))
    W.count(set(mine) == set(Di.index) and len(mine) == len(DON), "18C donor set")
    for k, (n1, n0, v) in mine.items():
        if k in Di.index:
            r = Di.loc[k]
            W.count(int(r["n_module"]) == n1 and int(r["n_nonmodule"]) == n0, f"18C {k} counts")
            W.add(v, r["auc"], f"18C {k} auc")
    for r in S18C.itertuples():
        g = DON[(DON.arm == r.arm) & (DON.dataset == r.level)]
        v = g["auc"].values
        W.add(float(np.median(v)) if len(v) >= min_d else np.nan, r.estimate, f"18C {r.arm} {r.level} median")
        W.count(len(v) == r.n_donors and int((g.n_module + g.n_nonmodule).sum()) == r.n_units,
                f"18C {r.arm} {r.level} counts")
    L.chk("B-10", W.bad == 0 and W.worst < TOL,
          "every 18C donor AUC and count recomputed from the Fig18BC_cells values equals Plot_18C_donors, "
          "and every 18C median and count follows from them", W.detail())

    # ---- 18D and 18E, from the unit table --------------------------------------------
    US = UM.copy()
    US["_a"] = US["arm"].map({a: i for i, a in enumerate(ARMS)})
    US = US.sort_values("_a", kind="stable").reset_index(drop=True)[["arm", "dataset", "donor_key", "region",
                                                                     "module_post_pct"]]
    W = Worst()
    OBS = {}
    for a in ARMS:
        u = US[US.arm == a]
        dons = sorted(pd.unique(u["donor_key"]))
        dcode = u["donor_key"].map({k: i for i, k in enumerate(dons)}).values
        OBS[a] = share(u["module_post_pct"].astype(float).values, dcode)
        a16 = A16[(A16.arm == a) & (A16.quantity == "module_post_pct")]["estimate"]
        W.add(OBS[a], a16.iloc[0] if len(a16) == 1 else np.nan, f"18D {a} observed share against 16A")
        for e in EFF:
            c = CUR[(CUR.arm == a) & (np.abs(CUR.effect - e) < 1e-9)].sort_values("n_donors")
            for t in c["target_share"].values:
                W.add(t, OBS[a] * e, f"18D {a} {e} target share")
            for tag, col in (("R1", "r1_rate"), ("R2", "r2_rate"), ("R3", "r3_rate")):
                okr = np.isfinite(c[col].values) & (c[col].values >= DETECT)
                cross = float(c["n_donors"].values[np.argmax(okr)]) if okr.any() else np.nan
                lv = f"{e:.2f}x"
                s1 = S18D[(S18D.arm == a) & (S18D.level == lv) & (S18D.quantity == "n_at_80") & (S18D.measure == tag)]
                s2 = S18D[(S18D.arm == a) & (S18D.level == lv) & (S18D.quantity == "rate_at_max_n") & (S18D.measure == tag)]
                W.add(cross, s1["estimate"].iloc[0], f"18D {a} {lv} {tag} crossing")
                W.add(c[col].values[-1], s2["estimate"].iloc[0], f"18D {a} {lv} {tag} rate at max n")
                W.count(int(s1["n_donors"].iloc[0]) == int(c["n_donors"].max()), f"18D {a} {lv} {tag} max n")
    L.chk("B-11", W.bad == 0 and W.worst < TOL,
          "every 18D crossing and rate at the largest donor count follows from the curves, each observed "
          "share equals the 16A estimate and every target share is that share times the effect", W.detail())

    W = Worst()
    Ri = REG.set_index(["arm", "min_regions"])
    for a in ARMS:
        u = US[US.arm == a]
        dons = sorted(pd.unique(u["donor_key"]))
        dcode = u["donor_key"].map({k: i for i, k in enumerate(dons)}).values
        nreg = np.bincount(dcode, minlength=len(dons))
        for k in range(3, max(3, int(nreg.max())) + 1):
            keep = np.where(nreg >= k)[0]
            sel = np.isin(dcode, keep)
            sh = share(u["module_post_pct"].astype(float).values[sel], dcode[sel]) if len(keep) >= share_floor else np.nan
            r = Ri.loc[(a, k)]
            where = f"18E {a} {k}+"
            W.count(int(r["donors"]) == len(keep) and int(r["donors_all"]) == len(dons)
                    and int(r["units"]) == int(sel.sum()), where + " counts")
            W.add(len(keep) / len(dons), r["donors_frac"], where + " donors_frac")
            W.add(sh, r["donor_share"], where + " donor_share")
    for r in S18E.itertuples():
        k = int(r.level.split("+")[0])
        g = Ri.loc[(r.arm, k)]
        W.add(r.estimate, g[r.quantity], f"18E {r.arm} {k}+ {r.quantity} stats against plot")
        W.count(int(r.n_units) == int(g["units"]) and int(r.n_donors) == int(g["donors"]), f"18E {r.arm} {k}+ counts")
    L.chk("B-12", W.bad == 0 and W.worst < TOL,
          "every 18E retention, unit count and donor share recomputed from the unit table equals its row, "
          "and every 18E stats row equals the plot data", W.detail())

    fl = parse_log(find(root, FIG_LOG))
    last = max(dt.datetime.strptime(lg["start"], "%Y-%m-%d %H:%M:%S") for lg in logs.values())
    fstart = dt.datetime.strptime(fl["start"], "%Y-%m-%d %H:%M:%S")
    L.chk("B-13", fl["verdict"] == "PASS" and fl["n_fail"] == 0 and fstart >= last,
          "figure log is PASS and started after the last analysis section",
          f"{fl['script']} | {fl['start']}")
    L._say("")

    # ---- workbook ------------------------------------------------------------------
    wb = Workbook()
    wb.remove(wb.active)
    ov = [
        ["Fig 18A, Supp Fig 14", "Can a nuclear-encoded mitochondrial panel replace the Module?",
         "outcome ~ cov4 + second predictor (Module post-IP% or %MT) + nuclear panel score, standardised, OLS, per arm, region set and panel, for the 9 outcomes",
         "22", "Stats_18A.csv; Plot_18A_verdicts.csv; Plot_18A_panels.csv", "Fig18A_verdicts; Fig18A_stats; Fig18A_panels", "Notes 2 to 5 and 11"],
        ["Fig 18B, Supp Fig 15", "Does any %MT cut reproduce the Module cells?",
         "Cells pooled over donors: AUC of %MT for Module membership, best cut by Youden's J, kappa at that cut",
         "23", "Stats_18B.csv; Plot_18B_cells.csv", "Fig18B_readout; Fig18B_stats; Fig18BC_cells", "Notes 2, 6 and 11"],
        ["Fig 18C", "Does %MT distinguish Module cells within a donor?",
         "Per donor with 20 or more Module and non-Module cells, the AUC of %MT; median over donors against 0.5",
         "23", "Stats_18C.csv; Plot_18C_donors.csv", "Fig18C_readout; Fig18C_stats; Fig18C_donors; Fig18BC_cells", "Notes 2, 3 and 7"],
        ["Fig 18D", "How many donors does each effect need?",
         "Detection rate of rules R1 to R3 against donor count, donors subsampled, Module post-IP% rescaled to 1.0x, 0.75x and 0.5x the observed donor share",
         "24", "Stats_18D.csv; Plot_18D_curves.csv", "Fig18D_stats; Fig18D_curves; Fig18DE_units", "Notes 2 and 8"],
        ["Fig 18E", "What does raising the region floor cost?",
         "Donors kept, donor share and split-half rs against the minimum number of regions per donor",
         "24", "Stats_18E.csv; Plot_18E_regions.csv", "Fig18E_stats; Fig18E_regions; Fig18DE_units", "Notes 2 and 9"],
    ]
    write_table(wb, "Overview", ["Panel", "Question", "Quantity", "Section", "Source file", "Sheets", "See"], ov,
                wrap=True, widths={"Panel": 14, "Question": 40, "Quantity": 60, "Section": 9,
                                   "Source file": 34, "Sheets": 38, "See": 18})
    cov = ", ".join(cfg["clinical"]["cov4"].values())
    notes = [
        ["1", "Arms", "M4 = Mic13_PVM3, the masked primary arm. PVM_3 = PVM3_all, the unmasked reference. PVM_1 = PVM1_all, the second supertype. HKG8 = the housekeeping control construct. 18A, 18B and 18C keep the region sets NeoC and NonNeoC apart; 18D and 18E pool them within donor, as 16A and 16B do."],
        ["2", "Inference", f"18A: OLS p, as in the Figure 12 clinical model (no interval is written). 18C: donor bootstrap, B = {B_:,}, 95% percentile interval, p = min(1, 2(k+1)/(B+1)) for the median against 0.5, smallest value 2/(B+1) = {2 / (B_ + 1):.2e}. 18B, 18D and 18E carry no test. fdr is Benjamini-Hochberg within the family column; rows whose family ends in 'descriptive' carry no p."],
        ["3", "Floors", f"A unit is one donor in one region on one arm with at least {int(cfg['floors']['cells_per_unit'])} cells of the arm. 18A: donors with at least {minc} labelled cells of the arm in the region set, and {min_d} donors per model. 18C: donors with at least {minc} Module and {minc} non-Module cells, and {min_d} donors per median. A row below a floor has an empty estimate and is untestable."],
        ["4", "Fig 18A model", f"Per arm, region set, panel and second predictor: outcome ~ {cov} + second predictor + panel score, all standardised; APOE e4 dose leaves the covariates when it is the outcome. The panel score is, per cell, the mean over the panel's genes (detected in at least 10% of the arm's cells in that region set, at least {int(F18['min_panel_genes'])} genes) of log1p(count / total UMIs x 10,000), and per donor the median over the donor's cells. The 15 nuclear-encoded panels and the matched reference (the null) are listed in Fig18A_panels. VIF of the two predictors; {vr:g} or more is unreadable. BH within arm, region set, panel and predictor over the 9 outcomes, each coefficient separately."],
        ["5", "Fig 18A verdicts", f"Per row: Untestable (not fitted), Unreadable (VIF >= {vr:g}), otherwise 'second only' (the Module, or %MT, carries the outcome and the panel does not), NE panel only, Both or Neither, from the two FDRs at {alpha:g}. The primary rows are the M4 cognition rows. A real panel that wins no more rows than the matched reference has shown nothing. In Fig18A_verdicts the verdict column is an Excel formula over the same row. The donor panel scores come from the step 02 count matrix, which this workbook does not carry, so the builder does not refit 18A; it checks the verdicts and tallies against the stats file and the log (checks B-07, B-08)."],
        ["6", "Fig 18B", f"Cells pooled over donors, per arm and region set. AUC = probability that a Module cell has higher %MT than a non-Module cell (ties counted half). Best cut: the %MT value maximising Youden's J, cells above it called Module; kappa at that cut. Bands on the real arms: AUC below {BANDS[0]:.2f} does not recover the Module cells; {BANDS[0]:.2f} to {BANDS[1]:.2f} partial; {BANDS[1]:.2f} or more close to a threshold (Fig18B_readout, Excel formula). Fig18BC_cells holds every cell's %MT and Module membership; the builder recomputes every value from it (check B-09)."],
        ["7", "Fig 18C", f"Per donor with at least {minc} cells on each side, the AUC of %MT for Module membership (Fig18C_donors, recomputed by the builder from Fig18BC_cells, check B-10). Median over donors with a donor-bootstrap interval; BH within region set over the four arms. The reading column of Fig18C_readout is an Excel formula: untestable, does not distinguish (fdr >= {alpha:g}), or the direction of the median."],
        ["8", "Fig 18D", f"Per arm, donors subsampled without replacement at {', '.join(str(n) for n in GRID)} donors (up to the arm's count), {int(PW['draws'])} draws per point, unit Module post-IP% rescaled so the donor share is {', '.join(f'{e:g}' for e in EFF)} times the observed share. R1: donor share above the 95th percentile of its within-region shuffle null. R2: split-half rs at least {float(PW['split_rule']):.2f}. R3: MMSE ~ cov4 + donor Module post-IP%, Module coefficient OLS p < {alpha:g}. n_at_80 is the first donor count with detection rate at least {DETECT:.0%}; empty means not reached. The builder checks the crossings against the curves and the observed shares against 16A from Fig18DE_units (check B-11); the rates themselves come from random draws in R and are not recomputed. Figure 18D draws R3 on M4 and PVM_3 at 1.0x and 0.75x: the donor requirement is read from R3, the clinical association rule. R1, R2, PVM_1, HKG8 and 0.5x are reported in this workbook only. When n_at_80 equals the largest donor count of the arm, the draws at that point hold all or nearly all of the arm's donors, so the rate there is close to a single result on the observed data."],
        ["9", "Fig 18E", f"Per arm and minimum regions per donor: donors kept, the donor share on the kept donors ({share_floor}-donor floor) and the split-half rs ({min_d}-donor floor). Kept donors at a high minimum are those sampled in many regions, not a random subset. The builder recomputes retention, unit counts and donor shares from Fig18DE_units (check B-12); the split-half rs uses random splits and is not recomputed."],
        ["10", "Old values", "The archived run is not compared here (decision D20 dropped, revision 49). Results are reported as they stand."],
        ["11", "Supplementary Figures 14 and 15", "Figures 18A and 18B draw NeoC. Supplementary Figure 14 draws the 18A verdicts for NonNeoC (panel A: Module post-IP% with each panel; panel B: %MT with each panel), and Supplementary Figure 15 the 18B densities for NonNeoC. Their values are the NonNeoC rows of Fig18A_verdicts, Fig18A_stats, Fig18B_readout and Fig18B_stats."],
    ]
    write_table(wb, "Notes", ["Note", "Subject", "Detail"], notes, wrap=True, widths={"Note": 7, "Subject": 30, "Detail": 120})

    gl = [
        ["column", "panel", "Figure panel the row belongs to."],
        ["column", "arm", "M4, PVM_3, PVM_1 or HKG8 (Note 1)."],
        ["column", "level", "18A, 18B, 18C: the region set; 18D: the effect (1.00x, 0.75x, 0.50x); 18E: the minimum regions per donor."],
        ["column", "quantity, measure", "What the row estimates; 18A measure = predictor|panel|outcome; 18D measure = the rule."],
        ["column", "n_units, n_donors", "18B, 18C: cells and donors; 18D: largest donor count on the curve; 18E: units and donors kept."],
        ["column", "estimate, lo, hi, p, fdr, family", "Estimate, 95% interval where one exists, p, Benjamini-Hochberg value within family (Note 2)."],
        ["18A", "beta_second, beta_panel, vif", "Standardised OLS coefficients of the second predictor and the panel score; VIF of the two predictors."],
        ["18A", "predictor", "module (Module post-IP%) or mito (log1p %MT), the second predictor."],
        ["18A", "n_genes, n_in_matrix, n_detected", "In Fig18A_panels: the panel's genes, those in the count matrix, those detected in at least 10% of the arm's cells."],
        ["18B", "auc, best_cut, youden, kappa", "AUC of %MT for Module membership; best %MT cut (per cent); Youden's J and Cohen's kappa at it (Note 6)."],
        ["18C", "median_within_donor_auc", "Median over donors of the within-donor AUC (Note 7)."],
        ["18C", "n_module, n_nonmodule, auc", "In Fig18C_donors: the donor's Module and non-Module cells and its AUC."],
        ["18D", "n_at_80, rate_at_max_n", "First donor count with detection rate at least 80%; detection rate at the largest donor count."],
        ["18D", "r1_rate, r2_rate, r3_rate, r2_evaluable, r3_evaluable, draws, target_share", "In Fig18D_curves: detection rates, the share of draws in which R2 and R3 could be computed, valid draws, and the rescaled donor share."],
        ["18E", "donors_frac, donor_share, split_half", "Share of the arm's donors kept, donor share on them, split-half rs on them."],
        ["cells", "is_module, mito_pct", "In Fig18BC_cells: Module membership and %MT (13 protein-coding MT genes) of each labelled cell."],
        ["units", "post_ip_pct", "In Fig18DE_units: unit Module post-IP%, per cent of the arm's supertype cells."],
    ]
    write_table(wb, "Glossary", ["Scope", "Term", "Meaning"], gl, wrap=True, widths={"Scope": 10, "Term": 46, "Meaning": 100})

    # Fig18A verdicts, verdict as a formula over the row
    hA = ["arm", "region_set", "predictor", "panel", "outcome", "class", "n_donors", "fitted", "vif", "beta_second",
          "beta_panel", "fdr_second", "fdr_panel", "verdict"]
    cA = {h: get_column_letter(j + 1) for j, h in enumerate(hA)}
    rows = []
    for i, r in enumerate(V.itertuples(), start=2):
        c = lambda h: f"{cA[h]}{i}"
        ss = f"AND(ISNUMBER({c('fdr_second')}),{c('fdr_second')}<{alpha})"
        sn = f"AND(ISNUMBER({c('fdr_panel')}),{c('fdr_panel')}<{alpha})"
        f = (f'=IF(NOT({c("fitted")}),"Untestable",IF(NOT(AND(ISNUMBER({c("vif")}),{c("vif")}<{vr:g})),"Unreadable",'
             f'IF(AND({ss},NOT({sn})),"second only",IF(AND({sn},NOT({ss})),"NE panel only",'
             f'IF(AND({ss},{sn}),"Both","Neither")))))')
        rows.append([r.arm, r.dataset, r.predictor, PL.get(r.panel, r.panel), r.outcome, r.oclass, r.n_donors,
                     bool(r.fitted), r.vif, r.beta_second, r.beta_panel, r.fdr_second, r.fdr_panel, f])
    write_table(wb, "Fig18A_verdicts", hA, rows, p_cols=("fdr_second", "fdr_panel"),
                fmt_cols={k: "0.000" for k in ("vif", "beta_second", "beta_panel")}, widths={"verdict": 14})
    stats_sheet(wb, "Fig18A_stats", S18A)
    pcols = [c for c in PAN.columns if c != "run_id"]
    P2 = PAN[pcols].copy(); P2.insert(3, "panel_label", P2["panel"].map(PL))
    write_table(wb, "Fig18A_panels", list(P2.columns), P2.itertuples(index=False, name=None))

    # Fig18B readout
    hB = ["arm", "region_set", "n_cells", "n_donors", "auc", "band", "best_cut_pct", "youden", "kappa"]
    cB = {h: get_column_letter(j + 1) for j, h in enumerate(hB)}
    rows = []
    for i, (a, ds) in enumerate([(a, ds) for a in ARMS for ds in DATASETS], start=2):
        g = lambda q: Bi.loc[(a, ds, q)]
        ac = f"{cB['auc']}{i}"
        rows.append([a, ds, g("auc")["n_units"], g("auc")["n_donors"], g("auc")["estimate"],
                     f'=IF(NOT(ISNUMBER({ac})),"",IF({ac}<{BANDS[0]},"does not recover",IF({ac}<{BANDS[1]},"partial","close to a threshold")))',
                     g("best_cut")["estimate"], g("youden")["estimate"], g("kappa")["estimate"]])
    write_table(wb, "Fig18B_readout", hB, rows, fmt_cols={k: "0.000" for k in ("auc", "youden", "kappa")} |
                {"best_cut_pct": "0.00"}, widths={"band": 22})
    stats_sheet(wb, "Fig18B_stats", S18B)

    # Fig18C readout
    hC = ["arm", "region_set", "n_donors", "n_cells", "median_auc", "lo", "hi", "p", "fdr", "reading"]
    cC = {h: get_column_letter(j + 1) for j, h in enumerate(hC)}
    rows = []
    for i, r in enumerate(S18C.itertuples(), start=2):
        e, fd = f"{cC['median_auc']}{i}", f"{cC['fdr']}{i}"
        rows.append([r.arm, r.level, r.n_donors, r.n_units, r.estimate, r.lo, r.hi, r.p, r.fdr,
                     f'=IF(NOT(ISNUMBER({e})),"untestable (fewer than {min_d} donors)",IF(AND(ISNUMBER({fd}),{fd}<{alpha}),'
                     f'IF({e}>0.5,"Module cells carry more %MT","Module cells carry less %MT"),"does not distinguish"))'])
    write_table(wb, "Fig18C_readout", hC, rows, p_cols=("p", "fdr"),
                fmt_cols={k: "0.000" for k in ("median_auc", "lo", "hi")}, widths={"reading": 32})
    stats_sheet(wb, "Fig18C_stats", S18C)
    dc = [c for c in DON.columns if c != "run_id"]
    write_table(wb, "Fig18C_donors", ["arm", "region_set"] + dc[2:], DON[dc].itertuples(index=False, name=None))

    stats_sheet(wb, "Fig18D_stats", S18D)
    cc = [c for c in CUR.columns if c != "run_id"]
    write_table(wb, "Fig18D_curves", cc, CUR[cc].itertuples(index=False, name=None))
    stats_sheet(wb, "Fig18E_stats", S18E)
    rc = [c for c in REG.columns if c != "run_id"]
    write_table(wb, "Fig18E_regions", rc, REG[rc].itertuples(index=False, name=None))
    write_table(wb, "Fig18BC_cells", ["arm", "region_set", "donor_key", "is_module", "mito_pct"],
                ([a, ds, k, bool(m), x] for a, ds, k, m, x in CS.itertuples(index=False, name=None)))
    write_table(wb, "Fig18DE_units", ["arm", "region_set", "donor_key", "region", "post_ip_pct"],
                US.itertuples(index=False, name=None), fmt_cols={"post_ip_pct": "0.000"})

    order = ["Overview", "Notes", "Glossary",
             "Fig18A_verdicts", "Fig18A_stats", "Fig18A_panels",
             "Fig18B_readout", "Fig18B_stats",
             "Fig18C_readout", "Fig18C_stats", "Fig18C_donors",
             "Fig18D_stats", "Fig18D_curves", "Fig18E_stats", "Fig18E_regions",
             "Fig18BC_cells", "Fig18DE_units"]
    if sorted(order) != sorted(wb.sheetnames):
        sys.exit(f"ERROR: sheet order list does not match the sheets written: {wb.sheetnames}")
    wb._sheets = [wb[n] for n in order]
    wb.active = 0

    texts = []
    for ws in wb.worksheets:
        for row in ws.iter_rows(values_only=True):
            texts += [v for v in row if isinstance(v, str)]
    bad = [t for t in texts if any(d_ in t for d_ in DASHES)]
    L.chk("B-14", not bad, "no en dash or em dash in the workbook text", f"{len(bad)} found")

    L._say("")
    L._say(f"Checks: {L.n['PASS']} PASS, {L.n['FAIL']} FAIL, {L.n['NOTE']} NOTE")
    log_path = out.with_name(out.stem + "_build_log.txt")
    if L.n["FAIL"]:
        L._say("VERDICT: FAIL | workbook not written")
        log_path.write_text("\n".join(L.lines) + "\n", encoding="utf-8")
        sys.exit(1)
    wb.calculation.fullCalcOnLoad = True
    wb.save(out)
    L._say(f"Output: {out}")
    L._say("VERDICT: PASS")
    log_path.write_text("\n".join(L.lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
