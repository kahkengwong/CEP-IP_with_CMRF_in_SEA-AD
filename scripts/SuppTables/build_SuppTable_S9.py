"""
build_SuppTable_S9.py

Builds SuppTable_S9_Fig4_Copula.xlsx, the workbook behind Supplementary Table
S9: the IPs, the S_i and K_i gates and the worked cells A to P in copula
coordinates (Figure 4, Supplementary Figures S8 to S10). The Part 10b outputs
are also read and checked, because the worked cells are matched to the Part
10b SKCS-core cells and to the Part 10 score CSVs, whose barcodes identify
them in Supplementary Table S7, sheet A. finalize_supp_tables.py turns this
workbook into Supplementary Table S9; it keeps the data sheets and replaces
the Contents sheet written here.

Usage:
    python3 build_SuppTable_S8_S9_step0_B1_counts.py --root /data   (first)
    python3 build_SuppTable_S9.py --root /data [--out <folder>]
--root defaults to $SEAAD_ROOT, else the current folder.

Inputs, under --root (archived and superseded folders are skipped):
    Part10b_IPRS/Part10b_IPRS_<arm>/Part10b_IPRS_CP10K_<tag>_results.xlsx
        (Part10b_CP10K_IPRS_code.r) IP_Results, Residual_Predictions
    Part10b_SKCS/Part10b_SKCS_<arm>/Part10b_CP10K_<tag>_SKCS_<region set>_EP_TREP_values.csv
        (Part10b_CP10K_SKCS_TREP_code.r)
    Part10/Part10_<arm>/Part10_<arm>_<region set>/Part10_CP10K_<tag>_Module1_<GOI>_<composite>_<region set>_Pooled_SKCS_scores.csv
        (Part10_CP10K_SKCS_code.r) barcodes of the worked cells
    Part10c_Copula1/Part10c_Copula1_<arm>/Part10c_CP10K_<tag>_Copula1_console_log_<time>.txt
        (Part10c_CP10K_Copula1_code.r) IP_u and IP_v conversions; newest log
    Part10d_Copula2/Part10d_Copula2_<arm>/Part10d_CP10K_<tag>_Copula2_console_log_<time>.txt
        (Part10d_CP10K_Copula2_code.r) gate conversions; newest log
    Manuscript/Supp_Tables/build/B1_group_donor_counts.csv
        (build_SuppTable_S8_S9_step0_B1_counts.py)
<arm> is Mic13_PVM3, PVM_3, PVM_1 or HKG8; <tag> is M4, PVM_3, PVM_1 or HKG8.

Gates (sheet B_copula_gates). The Part 10d log prints each gate count k with
its copula coordinate, as k/n in the original run and as k/(n + 1) in the
deposited scripts. Both forms are read. The sheet gives u_gate and v_gate as
k/(n + 1). The columns "Form in the log", "u_gate as logged" and "v_gate as
logged" are always written in the form of the published table, cells/n, with
k/n computed from the logged k and n and rounded to four decimals as the
original log printed it. When the log prints k/n, the computed value is
checked against the logged one.

Output:
    <root>/Manuscript/Supp_Tables/build/SuppTable_S9_Fig4_Copula.xlsx
        sheets Contents, A_copula_IP, B_copula_gates, C_worked_cells
    <root>/Manuscript/Supp_Tables/build/SuppTable_S9_Fig4_Copula_build_log.txt
        ends with VERDICT: PASS or VERDICT: FAIL

Every value is copied from these files without modification, except the copula
coordinates, which are computed from the logged counts, and the shares in
sheet C. The script stops if any check fails, and the workbook is then not
written.

Requirements: Python 3.9+, pandas, openpyxl.
"""

import argparse
import atexit
import glob
import io
import os
import re
import sys

import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Font
from openpyxl.utils import get_column_letter

#####
# Settings
#####
SETS = ["M4", "PVM_3", "PVM_1", "HKG8"]
REGION_SETS = ["NeoC", "NonNeoC"]
GOI = {"M4": "MT-CO3", "PVM_3": "MT-CO3", "PVM_1": "MT-CO3", "HKG8": "ACTB"}
COMPOSITE = {"M4": "cDFG", "PVM_3": "cDFG", "PVM_1": "cDFG", "HKG8": "cHKG8"}
GROUPS = ["Pre-IP TREP", "Pre-IP non-TREP", "Post-IP TREP", "Post-IP non-TREP"]
SKIP = ("Archived", "archive", "R_export_scran", "Code_GitHub", "superseded")
P_ZERO_TEXT = "< 2.2 × 10⁻¹⁶"

# Values already reported for Figure 2 and Supplementary Figure 5
# (Supplementary Methods S4, Table 3); this table must agree with them.
EXPECTED_TREP_ALL = {("M4", "NeoC"): 1565, ("M4", "NonNeoC"): 1499,
                     ("PVM_3", "NeoC"): 9072, ("PVM_3", "NonNeoC"): 7904,
                     ("PVM_1", "NeoC"): 2199, ("PVM_1", "NonNeoC"): 1990,
                     ("HKG8", "NeoC"): 160, ("HKG8", "NonNeoC"): 208}
EXPECTED_CORE_DONORS = {("M4", "NeoC"): 74, ("M4", "NonNeoC"): 79,
                        ("PVM_3", "NeoC"): 84, ("PVM_3", "NonNeoC"): 83,
                        ("PVM_1", "NeoC"): 82, ("PVM_1", "NonNeoC"): 82,
                        ("HKG8", "NeoC"): 55, ("HKG8", "NonNeoC"): 56}

# Worked cells as printed in Figure 4 and Supplementary Figures 8 to 10.
# Expression space: (panel, set, region set, x_i, rank x_i, y_i, rank y_i, group)
EXPR_CELLS = {
    "A": ("Figure 4A", "M4", "NeoC", 1.8209, 725, 0.9344, 459, "Pre-IP TREP"),
    "B": ("Figure 4B", "HKG8", "NeoC", 2.7345, 1346, 1.5275, 1032, "Post-IP non-TREP"),
    "E": ("Supplementary Figure S8A", "PVM_3", "NeoC", 3.2135, 11232, 2.4900, 10686, "Post-IP TREP"),
    "F": ("Supplementary Figure S8B", "PVM_1", "NeoC", 2.7748, 2499, 1.5316, 1834, "Post-IP non-TREP"),
    "I": ("Supplementary Figure S9A", "M4", "NonNeoC", 2.2005, 782, 1.1033, 419, "Pre-IP TREP"),
    "J": ("Supplementary Figure S9B", "PVM_3", "NonNeoC", 2.9878, 7922, 2.6856, 8766, "Post-IP TREP"),
    "K": ("Supplementary Figure S9C", "PVM_1", "NonNeoC", 3.2983, 2061, 1.8387, 1460, "Post-IP non-TREP"),
    "L": ("Supplementary Figure S9D", "HKG8", "NonNeoC", 0.6701, 45, 0.7885, 211, "Pre-IP TREP"),
}
# Score space: (panel, set, region set, S_i, rank S_i, K_i, rank K_i, in SKCS-core)
SCORE_CELLS = {
    "C": ("Figure 4C", "M4", "NeoC", 0.5883, 1705, 0.9079, 1939, True),
    "D": ("Figure 4D", "HKG8", "NeoC", 0.3686, 1146, 0.3627, 1022, False),
    "G": ("Supplementary Figure S8C", "PVM_3", "NeoC", 0.4058, 9202, 0.6663, 7585, True),
    "H": ("Supplementary Figure S8D", "PVM_1", "NeoC", 0.0023, 644, 0.3031, 311, False),
    "M": ("Supplementary Figure S10A", "M4", "NonNeoC", 0.4167, 1382, 0.6596, 860.5, True),
    "N": ("Supplementary Figure S10B", "PVM_3", "NonNeoC", 0.4121, 7799, 0.5064, 2323, True),
    "O": ("Supplementary Figure S10C", "PVM_1", "NonNeoC", 0.5601, 2127, 0.8668, 2319.5, True),
    "P": ("Supplementary Figure S10D", "HKG8", "NonNeoC", 0.3674, 864, 0.3933, 781, False),
}

#####
# Helpers
#####
FAILURES = []


def check(ok, message):
    """Record a failed check; the workbook is not written if any check fails."""
    if not ok:
        FAILURES.append(message)


def find_input(root, part_dir, pattern, label, newest=False):
    """One input file under root/part_dir (searched recursively).

    Archived and superseded folders are skipped. Data files must occur
    exactly once. For console logs (newest=True) the log with the latest time
    stamp in its name is used, because each run of a part writes a new log
    and overwrites its other outputs.
    """
    hits = sorted(p for p in glob.glob(os.path.join(root, part_dir, "**", pattern), recursive=True)
                  if not any(s in p.replace("\\", "/") for s in SKIP))
    if not hits:
        sys.exit(f"Input not found: {label} ({part_dir}/**/{pattern})")
    if newest:
        hits.sort(key=os.path.basename)
        path = hits[-1]
        print(f"  {label}: {path}" + (f"  (+{len(hits) - 1} older logs)" if len(hits) > 1 else ""))
        return path
    if len(hits) != 1:
        sys.exit(f"Expected one {label} file, found {len(hits)}: " + "; ".join(hits))
    print(f"  {label}: {hits[0]}")
    return hits[0]


class _Tee(io.TextIOBase):
    """Writes to the console and keeps a copy for the build log."""

    def __init__(self, stream, buf):
        self.stream, self.buf = stream, buf

    def write(self, s):
        self.stream.write(s)
        self.buf.write(s)
        return len(s)

    def flush(self):
        self.stream.flush()


BUILD = {"log": None, "written": False}
_LOG_BUF = io.StringIO()


def _write_build_log():
    """Saves the console output beside the workbook, with the verdict last."""
    sys.stdout.write("VERDICT: %s\n" % ("PASS" if BUILD["written"] else "FAIL"))
    if BUILD["log"]:
        os.makedirs(os.path.dirname(BUILD["log"]), exist_ok=True)
        with open(BUILD["log"], "w", encoding="utf-8") as fh:
            fh.write(_LOG_BUF.getvalue())


def kn4(k, n):
    """k/n rounded to four decimals as the original Part 10d log printed it (sprintf %.4f)."""
    return float("%.4f" % (k / n))


def read_log(path):
    with open(path, encoding="utf-8", errors="replace") as fh:
        return fh.read().replace("\r", "")


#####
# Read inputs
#####
def main():
    ap = argparse.ArgumentParser(description="Build SuppTable_S9_Fig4_Copula.xlsx")
    ap.add_argument("--root", default=os.environ.get("SEAAD_ROOT") or os.getcwd(),
                    help="SEA-AD folder (default: $SEAAD_ROOT, else the current folder)")
    ap.add_argument("--out", default=None,
                    help="folder for the workbook (default: <root>/Manuscript/Supp_Tables/build)")
    args = ap.parse_args()
    root = os.path.abspath(args.root)
    if not os.path.isdir(root):
        sys.exit(f"Root folder not found: {root}")
    build_dir = os.path.join(root, "Manuscript", "Supp_Tables", "build")
    out_dir = args.out or build_dir
    os.makedirs(out_dir, exist_ok=True)
    out = os.path.join(out_dir, "SuppTable_S9_Fig4_Copula.xlsx")
    BUILD["log"] = os.path.join(out_dir, "SuppTable_S9_Fig4_Copula_build_log.txt")
    sys.stdout = _Tee(sys.stdout, _LOG_BUF)
    sys.stderr = _Tee(sys.stderr, _LOG_BUF)
    atexit.register(_write_build_log)
    print("build_SuppTable_S9.py")
    print(f"Root: {root}")
    print(f"Output: {out}")

    print("Inputs:")
    iprs, rx, ep, logs_c, logs_d = {}, {}, {}, {}, {}
    for s in SETS:
        p = find_input(root, "Part10b_IPRS", f"Part10b_IPRS_CP10K_{s}_results.xlsx", f"IPRS {s}")
        iprs[s] = pd.read_excel(p, "IP_Results")
        rx[s] = pd.read_excel(p, "Residual_Predictions")
        logs_c[s] = read_log(find_input(root, "Part10c_Copula1", f"Part10c_CP10K_{s}_Copula1_console_log_*.txt",
                                        f"Part10c log {s}", newest=True))
        logs_d[s] = read_log(find_input(root, "Part10d_Copula2", f"Part10d_CP10K_{s}_Copula2_console_log_*.txt",
                                        f"Part10d log {s}", newest=True))
        for r in REGION_SETS:
            p = find_input(root, "Part10b_SKCS", f"Part10b_CP10K_{s}_SKCS_{r}_EP_TREP_values.csv",
                           f"EP_TREP {s} {r}")
            ep[(s, r)] = pd.read_csv(p)
    b1 = os.path.join(build_dir, "B1_group_donor_counts.csv")
    if not os.path.isfile(b1):
        sys.exit(f"Input not found: {b1} (run build_SuppTable_S8_S9_step0_B1_counts.py first)")
    print(f"  Donor counts: {b1}")
    donors = pd.read_csv(b1)
    # Pooled score CSVs of Part 10, to give every worked cell its barcode in
    # Supplementary Table S7, sheet A (cells outside the SKCS-core are not in Part 10b).
    scores = {}
    for s in SETS:
        for r in REGION_SETS:
            p = find_input(root, "Part10", f"Part10_CP10K_{s}_Module1_*_{r}_Pooled_SKCS_scores.csv",
                           f"Scores {s} {r}")
            t = pd.read_csv(p)
            scores[(s, r)] = t[~t["is_GOI_neg"].astype(bool)].reset_index(drop=True)

    #####
    # Part10c: n, IP_u and IP_v per set and region set
    #####
    rx_conv = re.compile(
        r"IP conversion: IP_x = ([\d.]+) -> IP_y = f\(IP_x\) = ([\d.]+)\s+\|\s+boundary IP_v = "
        r"Fc\(IP_y\) = #\{\S+ <= ([\d.]+)\}/\(n\+1\) = (\d+)/(\d+) = ([\d.]+)\s+\|\s+IP_u = "
        r"Fg\(IP_x\) = (\d+)/(\d+) = ([\d.]+)")
    conv = {}
    for s in SETS:
        for m in rx_conv.finditer(logs_c[s]):
            ipx = float(m.group(1))
            ds = None
            for r in REGION_SETS:
                row = iprs[s].loc[iprs[s]["Dataset"] == r].iloc[0]
                if abs(row["IP_Value"] - ipx) < 6e-5:
                    ds = r
            check(ds is not None, f"Part10c log {s}: IP_x {ipx} matches no IPRS dataset")
            if ds is None:
                continue
            n_plus_1 = int(m.group(5))
            conv[(s, ds)] = dict(ipx=ipx, ipy=float(m.group(2)), n=n_plus_1 - 1,
                                 k_v=int(m.group(4)), ip_v_log=float(m.group(6)),
                                 k_u=int(m.group(7)), ip_u_log=float(m.group(9)))
    for s in SETS:
        for r in REGION_SETS:
            check((s, r) in conv, f"Part10c log {s}: no IP conversion for {r}")

    #####
    # Part10d: gate counts per set and region set (logged as k/n)
    #####
    rx_gate = re.compile(r"Gate conversion: S_i_norm > 0\.05 -> u = (\d+)/(\d+) = ([\d.]+)\s+\|\s+"
                         r"K_i > 0\.5 -> v = (\d+)/(\d+) = ([\d.]+)")
    rx_core = re.compile(r"Non-SKCS cells: (\d+) \| SKCS-core cells: (\d+)")
    gates = {}
    for s in SETS:
        g_hits = list(rx_gate.finditer(logs_d[s]))
        c_hits = list(rx_core.finditer(logs_d[s]))
        check(len(g_hits) == 2 and len(c_hits) == 2, f"Part10d log {s}: expected two gate and two core lines")
        for g, c in zip(g_hits, c_hits):
            # Logs written before the k/(n + 1) change print k/n; later logs print k/(n + 1).
            den = int(g.group(2))
            ds, form = None, None
            for r in REGION_SETS:
                if (s, r) in conv and conv[(s, r)]["n"] == den:
                    ds, form = r, "cells/n"
                elif (s, r) in conv and conv[(s, r)]["n"] + 1 == den:
                    ds, form = r, "cells/(n + 1)"
            check(ds is not None, f"Part10d log {s}: denominator {den} matches no Part10c dataset")
            if ds is None:
                continue
            n = conv[(s, ds)]["n"]
            gates[(s, ds)] = dict(n=n, den=den, form=form, k_s=int(g.group(1)), u_log=float(g.group(3)),
                                  k_k=int(g.group(4)), v_log=float(g.group(6)),
                                  non_core=int(c.group(1)), core=int(c.group(2)))
            print(f"  Part10d log {s} {ds}: gates logged as {form}")

    #####
    # Checks across sources
    #####
    for s in SETS:
        for r in REGION_SETS:
            row = iprs[s].loc[iprs[s]["Dataset"] == r].iloc[0]
            comps = [row[c] for c in ["C1_Precision", "C2_Statistical", "C3_Abruptness",
                                      "C4_GAM_Ribbon", "C5_Smoothness"]]
            check(abs(sum(comps) / 5 - row["IPRS"]) < 1e-6, f"IPRS {s} {r}: not the mean of C1 to C5")
            tier = ("Very Strong" if row["IPRS"] >= 0.8 else "Strong" if row["IPRS"] >= 0.6 else
                    "Moderate" if row["IPRS"] >= 0.4 else "Weak" if row["IPRS"] >= 0.2 else "Failed")
            check(tier == row["IPRS_Tier"], f"IPRS {s} {r}: tier {row['IPRS_Tier']} against {tier}")
            check(int(row["N_Purple_Cells"]) == EXPECTED_TREP_ALL[(s, r)],
                  f"IPRS {s} {r}: TREP cells {row['N_Purple_Cells']} against Figure 2")
            check(int(row["NegPos_Crossings"]) == 1, f"IPRS {s} {r}: crossings != 1")
            cv = conv[(s, r)]
            check(abs(cv["ipy"] - row["IP_Y_Value"]) < 6e-5, f"{s} {r}: IP_y in Part10c log != IPRS")
            check(abs(cv["k_v"] / (cv["n"] + 1) - cv["ip_v_log"]) < 6e-5, f"{s} {r}: IP_v arithmetic")
            check(abs(cv["k_u"] / (cv["n"] + 1) - cv["ip_u_log"]) < 6e-5, f"{s} {r}: IP_u arithmetic")
            e = ep[(s, r)]
            gt = gates[(s, r)]
            check(gt["core"] == len(e), f"{s} {r}: SKCS-core in Part10d log != EP_TREP rows")
            check(gt["core"] + gt["non_core"] == cv["n"], f"{s} {r}: core + non-core != n")
            check(abs(gt["k_s"] / gt["den"] - gt["u_log"]) < 6e-5, f"{s} {r}: logged u_gate arithmetic")
            check(abs(gt["k_k"] / gt["den"] - gt["v_log"]) < 6e-5, f"{s} {r}: logged v_gate arithmetic")
            if gt["form"] == "cells/n":
                check(kn4(gt["k_s"], gt["n"]) == gt["u_log"], f"{s} {r}: u_gate k/n != logged value")
                check(kn4(gt["k_k"], gt["n"]) == gt["v_log"], f"{s} {r}: v_gate k/n != logged value")
            counts = e["group_4"].value_counts()
            for g in GROUPS:
                d = donors[(donors.arm == s) & (donors.dataset == r) & (donors.group == g)]
                check(len(d) == 1 and int(d.cells.iloc[0]) == int(counts.get(g, 0)),
                      f"{s} {r} {g}: donor table cells != EP_TREP")
            d = donors[(donors.arm == s) & (donors.dataset == r) & (donors.group == "SKCS-core (all)")]
            check(int(d.donors.iloc[0]) == EXPECTED_CORE_DONORS[(s, r)], f"{s} {r}: SKCS-core donors != Figure 2")
            check(int(donors[(donors.arm == s) & (donors.dataset == r)].unmatched_barcodes.sum()) == 0,
                  f"{s} {r}: barcodes without donor")

    #####
    # Sheet A: IP detection and IPRS
    #####
    a_rows = []
    for s in SETS:
        for r in REGION_SETS:
            row = iprs[s].loc[iprs[s]["Dataset"] == r].iloc[0]
            p = row["Stat_P_Raw"]
            a_rows.append({
                "Set": s, "Region set": r, "GOI": GOI[s], "Composite": COMPOSITE[s],
                "Cells with GOI >0": conv[(s, r)]["n"], "Donors": int(row["N_Donors"]),
                "DE, donor-aware (fraction)": row["DE_Full_Pop"],
                "DE sizing the TREP set (fraction)": row["DE_Purple_Sizing"],
                "TREP cells": int(row["N_Purple_Cells"]),
                "IP_x": row["IP_Value"], "IP_x 95% CI lower": row["CI_Lower_95"],
                "IP_x 95% CI upper": row["CI_Upper_95"], "CI width": row["CI_Width"],
                "GOI range of the R(x) grid": row["GOI_Range"],
                "CI method (lower; upper)": f"{row['CI_Lower_Method'].capitalize()}; {row['CI_Upper_Method'].capitalize()}",
                "Negative-to-positive crossings": int(row["NegPos_Crossings"]),
                "IP_y = f(IP_x)": row["IP_Y_Value"],
                "IP_y lower (IP_x CI mapped)": row["IP_Y_CI_Lower"],
                "IP_y upper (IP_x CI mapped)": row["IP_Y_CI_Upper"],
                "Test (TREP residuals at or above IP_x against below; one-sided)": {"MannWhitneyU": "Mann-Whitney U", "Welch_t": "Welch t"}.get(row["Stat_Test"], row["Stat_Test"]),
                "p-value": P_ZERO_TEXT if p == 0 else p,
                "C1 precision": row["C1_Precision"], "C2 test": row["C2_Statistical"],
                "C3 slope at IP_x": row["C3_Abruptness"], "C4 band width": row["C4_GAM_Ribbon"],
                "C5 smoothness and sign": row["C5_Smoothness"],
                "|R'(IP_x)|": row["Deriv_at_IP"], "R(x) range on grid": row["Rx_Amplitude"],
                "Mean 95% band width": row["Rx_MeanRibbon"], "R(x) total variation": row["Rx_TotalVar"],
                "IPRS": row["IPRS"], "IPRS tier": row["IPRS_Tier"],
                "Flag": "" if pd.isna(row["IPRS_Flag"]) else row["IPRS_Flag"],
            })
    sheet_a = pd.DataFrame(a_rows)

    #####
    # Sheet B: R(x) on the grid (every tenth of 2,000 points, as exported)
    #####
    b_parts = []
    for s in SETS:
        t = rx[s].copy()
        t.insert(0, "Set", s)
        t = t.rename(columns={"Dataset": "Region set", "Gene_Set": "Composite", "GOI": "GOI value",
                              "r_hat": "R(x)", "se": "SE", "r_upper": "R(x) + 1.96 SE",
                              "r_lower": "R(x) - 1.96 SE"})
        b_parts.append(t)
    sheet_b = pd.concat(b_parts, ignore_index=True)

    #####
    # Sheet C: cells and donors in the four groups of the SKCS-core
    #####
    c_rows = []
    for s in SETS:
        for r in REGION_SETS:
            sub = donors[(donors.arm == s) & (donors.dataset == r)]
            core_n = int(sub[sub.group == "SKCS-core (all)"].cells.iloc[0])
            for g in GROUPS + ["Pre-IP (all)", "Post-IP (all)", "SKCS-core (all)"]:
                x = sub[sub.group == g].iloc[0]
                c_rows.append({"Set": s, "Region set": r, "Group": g, "Cells": int(x.cells),
                               "Donors": int(x.donors),
                               "Share of SKCS-core (%)": round(100 * x.cells / core_n, 1)})
    sheet_c = pd.DataFrame(c_rows)

    #####
    # Sheet D: IP_x and IP_y in copula coordinates
    #####
    d_rows = []
    for s in SETS:
        for r in REGION_SETS:
            cv = conv[(s, r)]
            d_rows.append({"Set": s, "Region set": r, "n": cv["n"], "n + 1": cv["n"] + 1,
                           "IP_x": cv["ipx"], "Cells with GOI ≤ IP_x": cv["k_u"],
                           "IP_u = cells/(n + 1)": round(cv["k_u"] / (cv["n"] + 1), 4),
                           "IP_y": cv["ipy"], "Cells with composite ≤ IP_y": cv["k_v"],
                           "IP_v = cells/(n + 1)": round(cv["k_v"] / (cv["n"] + 1), 4)})
    sheet_d = pd.DataFrame(d_rows)

    #####
    # Sheet E: gates in copula coordinates
    #####
    e_rows = []
    for s in SETS:
        for r in REGION_SETS:
            gt = gates[(s, r)]
            n = gt["n"]
            f = ((n - 1) / (n + 1)) ** 2
            core_below = int((ep[(s, r)]["S_value"] * f <= 0.05).sum())
            e_rows.append({"Set": s, "Region set": r, "n": n, "n + 1": n + 1,
                           "Cells with S_i ≤ 0.05 (copula form)": gt["k_s"],
                           "Cells with K_i ≤ 0.5": gt["k_k"],
                           "u_gate = cells/(n + 1)": round(gt["k_s"] / (n + 1), 4),
                           "v_gate = cells/(n + 1)": round(gt["k_k"] / (n + 1), 4),
                           # published form (cells/n), whatever form this run's log used
                           "Form in the log": "cells/n",
                           "u_gate as logged": kn4(gt["k_s"], n),
                           "v_gate as logged": kn4(gt["k_k"], n),
                           "SKCS-core cells": gt["core"], "Non-SKCS-core cells": gt["non_core"],
                           "Factor ((n - 1)/(n + 1))²": round(f, 5),
                           "SKCS-core cells at or below the S_i gate in copula form": core_below})
    sheet_e = pd.DataFrame(e_rows)

    #####
    # Sheet F: worked cells, arithmetic and agreement with the SKCS-core table
    #####
    f_rows = []
    for c, (panel, s, r, x, rx_, y, ry, group) in EXPR_CELLS.items():
        n = conv[(s, r)]["n"]
        e = ep[(s, r)]
        hit = e[(e["GOI_expr"] - x).abs() < 5e-5]
        check(len(hit) == 1 and hit["group_4"].iloc[0] == group,
              f"Cell {c}: GOI {x} not found once in the {s} {r} SKCS-core with group {group}")
        sc = scores[(s, r)]
        bc = sc.loc[((sc["GOI_expr"] - x).abs() < 5e-5) & ((sc["cDFG_expr"] - y).abs() < 5e-5), "barcode"]
        check(len(bc) == 1, f"Cell {c}: not found once in the {s} {r} score CSV")
        f_rows.append({"Cell": c, "Panel": panel, "Set": s, "Region set": r,
                       "Barcode": bc.iloc[0] if len(bc) == 1 else None, "Space": "Expression",
                       "n": n, "First quantity": GOI[s], "First value": x, "First rank": rx_,
                       "u = rank/(n + 1)": round(rx_ / (n + 1), 4),
                       "Second quantity": COMPOSITE[s], "Second value": y, "Second rank": ry,
                       "v = rank/(n + 1)": round(ry / (n + 1), 4), "Group": group,
                       "Check": "GOI value and group found in Supplementary Table S7, sheet A"})
    for c, (panel, s, r, si, rs, ki, rk, in_core) in SCORE_CELLS.items():
        n = conv[(s, r)]["n"]
        gt = gates[(s, r)]
        u_sk, v_sk = rs / (n + 1), rk / (n + 1)
        u_g, v_g = gt["k_s"] / (n + 1), gt["k_k"] / (n + 1)
        inside = (u_sk > u_g) and (v_sk > v_g)
        check(inside == in_core, f"Cell {c}: copula position disagrees with its SKCS-core label")
        check((si > 0.05 and ki > 0.5) == in_core, f"Cell {c}: scores disagree with its SKCS-core label")
        sc = scores[(s, r)]
        f_sc = ((n - 1) / (n + 1)) ** 2
        bc = sc.loc[((sc["spear_norm"] * f_sc - si).abs() < 5e-5) & ((sc["kend_contrib"] - ki).abs() < 5e-5),
                    "barcode"]
        check(len(bc) == 1, f"Cell {c}: S_i {si} and K_i {ki} not found once in the {s} {r} score CSV")
        note = "Outside the SKCS-core; found in Supplementary Table S7, sheet A"
        group = "Outside the SKCS-core"
        if in_core:
            f = ((n - 1) / (n + 1)) ** 2
            e = ep[(s, r)]
            hit = e[((e["S_value"] * f - si).abs() < 5e-5) & ((e["K_value"] - ki).abs() < 5e-5)]
            check(len(hit) == 1, f"Cell {c}: S_i {si} and K_i {ki} not found once in the {s} {r} SKCS-core")
            if len(hit) == 1:
                group = hit["group_4"].iloc[0]
                note = f"Found in Supplementary Table S7, sheet A (S_i of Supplementary Methods S4 {hit['S_value'].iloc[0]:.4f})"
        f_rows.append({"Cell": c, "Panel": panel, "Set": s, "Region set": r,
                       "Barcode": bc.iloc[0] if len(bc) == 1 else None, "Space": "Score",
                       "n": n, "First quantity": "S_i (copula form)", "First value": si, "First rank": rs,
                       "u = rank/(n + 1)": round(u_sk, 4),
                       "Second quantity": "K_i", "Second value": ki, "Second rank": rk,
                       "v = rank/(n + 1)": round(v_sk, 4), "Group": group, "Check": note})
    sheet_f = pd.DataFrame(f_rows)

    #####
    # Sheet G: every SKCS-core cell
    #####
    g_parts = []
    for s in SETS:
        for r in REGION_SETS:
            t = ep[(s, r)].copy()
            t.insert(0, "Region set", r)
            t.insert(0, "Set", s)
            g_parts.append(t)
    sheet_g = pd.concat(g_parts, ignore_index=True).rename(columns={
        "barcode": "Barcode", "region": "Region", "GOI_expr": "GOI (log1p CP10K)", "S_value": "S_i",
        "K_value": "K_i", "EP_value": "EP", "EP_rank": "EP rank", "is_TREP": "TREP",
        "IP_status": "IP status", "group_4": "Group"})

    #####
    # Contents
    #####
    contents = pd.DataFrame([
        ["A_copula_IP", "IP_x and IP_y converted to copula coordinates (IP_u, IP_v)",
         "Figure 4A, 4B; Supplementary Figures 8A, 8B, 9", "Part 10c logs"],
        ["B_copula_gates", "S_i and K_i gates converted to copula coordinates (u_gate, v_gate)",
         "Figure 4C, 4D; Supplementary Figures 8C, 8D, 10", "Part 10d logs"],
        ["C_worked_cells", "Worked cells A to P: values, ranks and copula coordinates, with checks",
         "Figure 4; Supplementary Figures 8 to 10", "Figures; sheets A and B; Supplementary Table [Fig2], A_cells"],
        ["", "", "", ""],
        ["Notes", "1) GOI and composite: MT-CO3 and cDFG (M4, PVM_3, PVM_1); ACTB and cHKG8 (HKG8). "
                  "Expression values are log1p CP10K.", "", ""],
        ["", "2) Sheet B gives the gate coordinates as cells/(n + 1), the denominator of the cell coordinates "
             "and of the figures; the column 'Form in the log' shows how each Part 10d log printed them. "
             "The counts are those of the logs.", "", ""],
        ["", "3) S_i in sheets B and C (copula form) is 4(u - 0.5)(v - 0.5) from the expression-space copula "
             "coordinates, equal to ((n - 1)/(n + 1))² times the S_i of Supplementary Table [Fig2], "
             "A_cells.", "", ""],
        ["", "4) Barcode identifies each worked cell in Supplementary Table [Fig2], A_cells, which holds "
             "the copula coordinates of every cell.", "", ""],
    ], columns=["Sheet", "Contents", "Panels", "Source"])

    #####
    # Stop or write
    #####
    if FAILURES:
        print("\nChecks failed; workbook not written:")
        for f in FAILURES:
            print("  " + f)
        sys.exit(1)

    sheets = [("Contents", contents), ("A_copula_IP", sheet_d), ("B_copula_gates", sheet_e),
              ("C_worked_cells", sheet_f)]
    wb = Workbook()
    wb.remove(wb.active)
    plain, bold = Font(name="Arial", size=10), Font(name="Arial", size=10, bold=True)
    for name, df in sheets:
        df = df.rename(columns={"Set": "Arm"})
        ws = wb.create_sheet(name)
        ws.append(list(df.columns))
        for rec in df.itertuples(index=False):
            ws.append([None if (isinstance(v, float) and pd.isna(v)) else v for v in rec])
        for row in ws.iter_rows():
            for cell in row:
                cell.font = bold if cell.row == 1 else plain
        ws.freeze_panes = "A2"
        for i, col in enumerate(df.columns, start=1):
            width = max(len(str(col)), *(len(str(v)) for v in df.iloc[:200, i - 1])) if len(df) else len(str(col))
            ws.column_dimensions[get_column_letter(i)].width = min(max(10, width + 2), 60)
    wb.save(out)
    BUILD["written"] = True
    print(f"\nAll checks passed. Written: {out}")
    for name, df in sheets:
        print(f"  {name}: {len(df)} rows")


if __name__ == "__main__":
    main()
