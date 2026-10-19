#!/usr/bin/env python3
"""
compare_supp_tables.py

Compares rebuilt Supplementary Tables with the published ones (Zenodo), sheet
by sheet and cell by cell. Sheet names and their order are compared first.
In each sheet, numbers agree when they are within 1e-6 (relative), text must
be identical, and formula cells are compared by their formula text (a
trailing ".0" on a number, which Excel removes when it saves a file, is
ignored). Fonts and widths are not compared.

Each sheet gets one result:
  SAME                     every cell agrees
  SAME VALUES, OTHER ORDER the same rows (matched on the first column) in another order
  TEXT ONLY                numbers agree, some text or formula text differs
  WITHIN 0.1%              some numbers differ by more than 1e-6 but all by at most 0.1%
  OUTSIDE 0.1%             some numbers differ by more than 0.1%, or a cell is
                           empty in one table only
  DIFFERENT SIZE           the sheets have different numbers of rows or columns
  ONLY IN REBUILT / ONLY IN PUBLISHED   the sheet or file exists on one side only

Usage:
    python compare_supp_tables.py --rebuilt <folder> --published <folder> [--tables 13-17]
Defaults: --rebuilt <root>/Manuscript/Supp_Tables, --published /ref/Supp_Tables,
root from $SEAAD_ROOT or the current folder.
Outputs (in the rebuilt folder): compare_supp_tables_log.txt and
compare_supp_tables_cells.csv (each differing cell, up to 200 per sheet).
"""

import argparse
import csv
import glob
import math
import os
import pathlib
import re
import sys

import openpyxl

REL_TOL, ABS_TOL, LARGE_TOL = 1e-6, 1e-9, 1e-3
MAX_CELLS = 200
ORDER = ["SAME", "SAME VALUES, OTHER ORDER", "TEXT ONLY", "WITHIN 0.1%", "OUTSIDE 0.1%", "DIFFERENT SIZE",
         "ONLY IN REBUILT", "ONLY IN PUBLISHED"]
# Known causes, matched on "S<n> <sheet>" and printed under the matching sheet.
EXPLAIN = [
    (r"^S10 A_cells$",
     "UMAP_1 and UMAP_2 come from the Validation2 UMAP, which is not reproducible between computers "
     "(the container rerun of Part 10e Validation2 already differs from the original run)."),
    (r"^S10 (I_Fig8_summary|J_Fig8_tests)$",
     "Neighbour enrichment, kBET and the PrePost projection come from the Validation2 UMAP and subsampling "
     "(the container rerun of Part 10e Validation2 already differs from the original run)."),
    (r"^S3 [CD]_pooled_screen_",
     "best_rs in the rerun of Part 4 (PVM_3_pooled_screen_scatter.csv) differs from the original run "
     "by less than 1e-6 (relative), and rows with near-tied values change places in the sorted sheet."),
    (r"^S13 P_Fig12H_stats$",
     "Depth thinning uses seeded binomial draws. The container rerun of Stats_15H.csv "
     "already differs from the original run in these cells (verification/compare_all_log.txt)."),
    (r"^S19 I_Checks$",
     "The C1 checks record the date each source workbook was last modified."),
]


def find(folder, t):
    hits = glob.glob(os.path.join(folder, f"Supplementary_Table_S{t}_*.xlsx"))
    hits = [h for h in hits if not os.path.basename(h).startswith("~$")]
    return hits[0] if len(hits) == 1 else None


def rows(ws):
    out = [list(r) for r in ws.iter_rows(values_only=True)]
    while out and all(v is None for v in out[-1]):
        out.pop()
    width = max((max((j + 1 for j, v in enumerate(r) if v is not None), default=0) for r in out), default=0)
    return [r[:width] + [None] * (width - len(r)) for r in out]


def norm_formula(f):
    return re.sub(r"(?<![\w.])(\d+)\.0(?![\d])", r"\1", f.replace(" ", ""))


def is_num(v):
    return isinstance(v, (int, float)) and not isinstance(v, bool)


def cell_cmp(a, b):
    """Returns None when equal, else (kind, relative difference)."""
    if a is None and b is None:
        return None
    if is_num(a) and is_num(b):
        if math.isnan(a) and math.isnan(b):
            return None
        d = abs(a - b)
        if d <= ABS_TOL or d <= REL_TOL * max(abs(a), abs(b)):
            return None
        return ("num", d / max(abs(a), abs(b)))
    if isinstance(a, str) and isinstance(b, str) and a.startswith("=") and b.startswith("="):
        return None if norm_formula(a) == norm_formula(b) else ("text", None)
    if a is None or b is None:
        return ("empty", None)
    if is_num(a) or is_num(b):
        try:
            return cell_cmp(float(a), float(b))
        except (TypeError, ValueError):
            return ("empty", None)
    return None if str(a) == str(b) else ("text", None)


def compare_sheet(ra, rb):
    st, diffs, note = compare_cells(ra, rb)
    if st in ("SAME", "DIFFERENT SIZE") or len(ra) < 3:
        return st, diffs, note
    # Same rows in another order: rows matched on the first column when it is a
    # unique key in both tables (near-tied values can sort differently).
    ka, kb = [r[0] for r in ra[1:]], [r[0] for r in rb[1:]]
    if len(set(ka)) == len(ka) and set(ka) == set(kb):
        pos = {k: i for i, k in enumerate(kb)}
        rb2 = [rb[0]] + [rb[1 + pos[k]] for k in ka]
        st2, _, note2 = compare_cells(ra, rb2)
        if st2 in ("SAME", "WITHIN 0.1%"):
            moved = sum(1 for x, y in zip(ka, kb) if x != y)
            return ("SAME VALUES, OTHER ORDER", diffs,
                    f"{moved} rows in another order" + (f"; matched rows: {note2}" if note2 else ""))
    return st, diffs, note


def compare_cells(ra, rb):
    if len(ra) != len(rb) or (ra and len(ra[0]) != len(rb[0])):
        return "DIFFERENT SIZE", [], f"{len(ra)}x{len(ra[0]) if ra else 0} vs {len(rb)}x{len(rb[0]) if rb else 0}"
    diffs, worst, kinds = [], 0.0, set()
    for i, (x, y) in enumerate(zip(ra, rb)):
        for j, (a, b) in enumerate(zip(x, y)):
            r = cell_cmp(a, b)
            if r:
                kinds.add(r[0])
                if r[0] == "num":
                    worst = max(worst, r[1])
                diffs.append((i + 1, j + 1, a, b))
    if not diffs:
        return "SAME", [], ""
    if "empty" in kinds or worst > LARGE_TOL:
        status = "OUTSIDE 0.1%"
    elif "num" in kinds:
        status = "WITHIN 0.1%"
    else:
        status = "TEXT ONLY"
    return status, diffs, f"{len(diffs)} cells differ" + (f", largest relative difference {worst:.2e}" if worst else "")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    root = os.environ.get("SEAAD_ROOT", os.getcwd())
    ap.add_argument("--rebuilt", default=os.path.join(root, "Manuscript", "Supp_Tables"))
    ap.add_argument("--published", default="/ref/Supp_Tables")
    ap.add_argument("--tables", default="1-19")
    args = ap.parse_args()
    tables = []
    for part in args.tables.split(","):
        a, _, b = part.partition("-")
        tables += list(range(int(a), int(b or a) + 1))
    log, cells, tally = [], [], {}

    def say(s=""):
        log.append(s)
        print(s, flush=True)

    say(f"compare_supp_tables.py\n  rebuilt:   {args.rebuilt}\n  published: {args.published}\n")
    for t in tables:
        fa, fb = find(args.rebuilt, t), find(args.published, t)
        if not fa or not fb:
            st = "ONLY IN REBUILT" if fa else "ONLY IN PUBLISHED"
            if not fa and not fb:
                say(f"S{t}: not found in either folder")
                continue
            say(f"S{t}: {st}")
            tally[st] = tally.get(st, 0) + 1
            continue
        wa = openpyxl.load_workbook(fa)
        wb = openpyxl.load_workbook(fb)
        name_note = "" if os.path.basename(fa) == os.path.basename(fb) else f" (file names differ: {os.path.basename(fa)})"
        order_note = "" if wa.sheetnames == wb.sheetnames else "  sheet order or names differ"
        say(f"S{t}: {os.path.basename(fb)}{name_note}{order_note}")
        for s in list(dict.fromkeys(wb.sheetnames + wa.sheetnames)):
            if s not in wa.sheetnames or s not in wb.sheetnames:
                st = "ONLY IN PUBLISHED" if s in wb.sheetnames else "ONLY IN REBUILT"
                note = ""
            else:
                st, diffs, note = compare_sheet(rows(wa[s]), rows(wb[s]))
                for (i, j, a, b) in diffs[:MAX_CELLS]:
                    cells.append([f"S{t}", s, i, openpyxl.utils.get_column_letter(j), a, b])
            tally[st] = tally.get(st, 0) + 1
            say(f"    {st:<24} {s}" + (f"  ({note})" if note else ""))
            if st != "SAME":
                for pat, why in EXPLAIN:
                    if re.search(pat, f"S{t} {s}"):
                        say(f"                             {why}")
    say("")
    say("Sheets per result:")
    for k in ORDER:
        if k in tally:
            say(f"  {k:<24} {tally[k]}")
    outdir = pathlib.Path(args.rebuilt)
    (outdir / "compare_supp_tables_log.txt").write_text("\n".join(log) + "\n", encoding="utf-8")
    with open(outdir / "compare_supp_tables_cells.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["table", "sheet", "row", "column", "rebuilt", "published"])
        w.writerows(cells)
    say(f"Log: {outdir / 'compare_supp_tables_log.txt'}")


if __name__ == "__main__":
    main()
