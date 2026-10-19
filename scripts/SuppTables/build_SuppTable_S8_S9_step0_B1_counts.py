#!/usr/bin/env python3
"""
build_SuppTable_S8_S9_step0_B1_counts.py

Writes B1_group_donor_counts.csv: the cells and donors of the four pre-IP and
post-IP groups of the SKCS-core, per arm and region set. build_SuppTable_S8.py
and build_SuppTable_S9.py read this file and check their inputs against it.
Run this script before them.

Usage:
    python3 build_SuppTable_S8_S9_step0_B1_counts.py
    python3 build_SuppTable_S8_S9_step0_B1_counts.py --root /data
--root defaults to $SEAAD_ROOT, else the current folder.

Inputs (read only; nothing else is modified):
  Part10b_SKCS/Part10b_SKCS_<arm>/Part10b_CP10K_<tag>_SKCS_<region set>_EP_TREP_values.csv
      (Part10b_CP10K_SKCS_TREP_code.r) barcode and group_4 of every SKCS-core cell
  R_export_CP10K/<export arm>/<neocortex|nonneocortex>/gam_data/*_selected.csv
      (Part8_Data-export-CP10K_MT-CO3_8DFG.py) barcode and donor of every exported cell

Output:
  <root>/Manuscript/Supp_Tables/build/B1_group_donor_counts.csv
      arm, dataset, group, cells, donors, unmatched_barcodes
      groups: Pre-IP TREP, Pre-IP non-TREP, Post-IP TREP, Post-IP non-TREP,
              Pre-IP (all), Post-IP (all), SKCS-core (all)
  <root>/Manuscript/Supp_Tables/build/B1_group_donor_counts_log.txt
      ends with VERDICT: PASS or VERDICT: FAIL

Checks before writing (the file is not written if any fails):
  C1  every SKCS-core cell has a donor (unmatched_barcodes = 0)
  C2  cells per group equal the reported values (REF_GROUPS, Figure 3 and
      Supplementary Figure S7)
  C3  SKCS-core donors equal the reported values (REF_CORE_DONORS)
  C4  the four groups add up to the SKCS-core; pre-IP and post-IP (all) are
      the sums of their two groups
"""
import argparse
import glob
import os
import sys

import pandas as pd

_ap = argparse.ArgumentParser(description="Write B1_group_donor_counts.csv for Supplementary Tables S8 and S9")
_ap.add_argument("--root", default=os.environ.get("SEAAD_ROOT") or os.getcwd(),
                 help="SEA-AD folder (default: $SEAAD_ROOT, else the current folder)")
ROOT = os.path.abspath(_ap.parse_args().root)
if not os.path.isdir(ROOT):
    sys.exit("SEA-AD folder not found: %s" % ROOT)
OUT_DIR = os.path.join(ROOT, "Manuscript", "Supp_Tables", "build")
OUT = os.path.join(OUT_DIR, "B1_group_donor_counts.csv")
LOG = os.path.join(OUT_DIR, "B1_group_donor_counts_log.txt")
SKIP = ("Archived", "archive", "R_export_scran", "Code_GitHub", "superseded")

ARMS = ["M4", "PVM_3", "PVM_1", "HKG8"]
DATASETS = ["NeoC", "NonNeoC"]
EXPORT_DIR = {"M4": "Mic13_PVM3", "PVM_3": "PVM3_all", "PVM_1": "PVM1_all", "HKG8": "HKG8_matched"}
EXPORT_DS = {"NeoC": "neocortex", "NonNeoC": "nonneocortex"}
GROUPS = ["Pre-IP TREP", "Pre-IP non-TREP", "Post-IP TREP", "Post-IP non-TREP"]

# Values already reported (build_SuppTable_S7.py REF_GROUPS; build_SuppTable_S8.py
# EXPECTED_CORE_DONORS). Used as checks only.
REF_GROUPS = {
    ("M4", "NeoC"): (640, 31, 668, 80), ("M4", "NonNeoC"): (582, 44, 620, 57),
    ("PVM_3", "NeoC"): (3781, 101, 3836, 689), ("PVM_3", "NonNeoC"): (3103, 116, 3214, 508),
    ("PVM_1", "NeoC"): (787, 1, 921, 153), ("PVM_1", "NonNeoC"): (743, 6, 804, 124),
    ("HKG8", "NeoC"): (71, 0, 147, 32), ("HKG8", "NonNeoC"): (43, 0, 116, 40),
}
REF_CORE_DONORS = {("M4", "NeoC"): 74, ("M4", "NonNeoC"): 79, ("PVM_3", "NeoC"): 84,
                   ("PVM_3", "NonNeoC"): 83, ("PVM_1", "NeoC"): 82, ("PVM_1", "NonNeoC"): 82,
                   ("HKG8", "NeoC"): 55, ("HKG8", "NonNeoC"): 56}

os.makedirs(OUT_DIR, exist_ok=True)
lines, fails = [], []


def say(s=""):
    print(s)
    lines.append(s)


def check(ok, what):
    say("  [%s] %s" % ("PASS" if ok else "FAIL", what))
    if not ok:
        fails.append(what)


def keep(p):
    return not any(s in p.replace("\\", "/") for s in SKIP)


say("build_SuppTable_S8_S9_step0_B1_counts.py")
say("Root: %s" % ROOT)
rows = []
for arm in ARMS:
    for ds in DATASETS:
        say("")
        say("%s %s" % (arm, ds))
        fn = "Part10b_CP10K_%s_SKCS_%s_EP_TREP_values.csv" % (arm, ds)
        hits = [p for p in glob.glob(os.path.join(ROOT, "Part10b_SKCS", "**", fn), recursive=True) if keep(p)]
        if len(hits) != 1:
            check(False, "exactly one %s under Part10b_SKCS (found %d)" % (fn, len(hits)))
            continue
        e = pd.read_csv(hits[0], usecols=["barcode", "group_4"])
        say("  EP_TREP: %s (%d cells)" % (hits[0], len(e)))
        gd_dir = os.path.join(ROOT, "R_export_CP10K", EXPORT_DIR[arm], EXPORT_DS[ds], "gam_data")
        files = sorted(glob.glob(os.path.join(gd_dir, "*_selected.csv")))
        if not files:
            check(False, "gam_data *_selected.csv files in %s" % gd_dir)
            continue
        gd = pd.concat([pd.read_csv(f, usecols=["barcode", "donor"]) for f in files], ignore_index=True)
        gd = gd.drop_duplicates("barcode")
        say("  gam_data: %d files, %d barcodes" % (len(files), len(gd)))
        e["donor"] = e["barcode"].map(gd.set_index("barcode")["donor"])
        check(set(e["group_4"]) <= set(GROUPS), "group labels are the four groups (%s)" % sorted(set(e["group_4"])))
        defs = [(g, e["group_4"] == g) for g in GROUPS]
        defs += [("Pre-IP (all)", e["group_4"].str.startswith("Pre-IP")),
                 ("Post-IP (all)", e["group_4"].str.startswith("Post-IP")),
                 ("SKCS-core (all)", pd.Series(True, index=e.index))]
        got = {}
        for g, m in defs:
            sub = e[m]
            got[g] = dict(cells=int(len(sub)), donors=int(sub["donor"].dropna().nunique()),
                          unmatched=int(sub["donor"].isna().sum()))
            rows.append(dict(arm=arm, dataset=ds, group=g, cells=got[g]["cells"],
                             donors=got[g]["donors"], unmatched_barcodes=got[g]["unmatched"]))
            say("    %-18s cells %6d  donors %3d  unmatched %d" % (g, got[g]["cells"], got[g]["donors"],
                                                                  got[g]["unmatched"]))
        check(got["SKCS-core (all)"]["unmatched"] == 0, "C1 every SKCS-core cell has a donor")
        check(tuple(got[g]["cells"] for g in GROUPS) == REF_GROUPS[(arm, ds)],
              "C2 group cells %s = reported %s" % (tuple(got[g]["cells"] for g in GROUPS), REF_GROUPS[(arm, ds)]))
        check(got["SKCS-core (all)"]["donors"] == REF_CORE_DONORS[(arm, ds)],
              "C3 SKCS-core donors %d = reported %d" % (got["SKCS-core (all)"]["donors"], REF_CORE_DONORS[(arm, ds)]))
        check(sum(got[g]["cells"] for g in GROUPS) == got["SKCS-core (all)"]["cells"]
              and got["Pre-IP (all)"]["cells"] == got["Pre-IP TREP"]["cells"] + got["Pre-IP non-TREP"]["cells"]
              and got["Post-IP (all)"]["cells"] == got["Post-IP TREP"]["cells"] + got["Post-IP non-TREP"]["cells"],
              "C4 group sums")

say("")
if fails or len(rows) != len(ARMS) * len(DATASETS) * 7:
    say("FAIL: %d check(s) failed; %s not written" % (len(fails), OUT))
    if os.path.exists(OUT):
        # a copy from an earlier run would otherwise be read by the S8 and S9 builders
        os.remove(OUT)
        say("Removed the copy from an earlier run: %s" % OUT)
    code = 1
else:
    pd.DataFrame(rows).to_csv(OUT, index=False)
    say("Written: %s (%d rows)" % (OUT, len(rows)))
    say("All checks passed.")
    code = 0
say("VERDICT: %s" % ("PASS" if code == 0 else "FAIL"))
with open(LOG, "w", encoding="utf-8") as f:
    f.write("\n".join(lines) + "\n")
print("Log: %s" % LOG)
sys.exit(code)
