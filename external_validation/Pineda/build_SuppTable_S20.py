#!/usr/bin/env python3
"""
build_SuppTable_S20.py  (2026-10-27)

Builds Supplementary Table S20 (external validation 1, Pineda et al. 2024, BA9) from
R16_Pineda_figure_source_data.xlsx, the table written by R16_Pineda_main_figure.R, the two
mapping comparisons of R13b, and the arm of the earlier run without the authors' QC (R13c).
The A to I sheets hold every value that each panel of the main figure plots, so the figure
can be redrawn from this table. No value is computed. The script:
  - copies the 15 sheets of the R16 source data; sheet letters follow the figure panels;
  - checks them against the R14, R15, and R17 workbooks they came from;
  - adds sheet J (the MapMyCells comparisons: the same file mapped twice, and the nuclei of
    the authors' QC mapped on their own against all kept nuclei mapped together) and sheet K
    (the PVM_3 arm of the earlier run without the authors' QC, by %MT band and by donor);
  - adds the Contents sheet.

Usage (WSL):
    python3 build_SuppTable_S20.py \\
        --root /mnt/c/Users/Wong/Desktop/SEA-AD/External_Dataset_1 \\
        --supp /mnt/c/Users/Wong/Desktop/SEA-AD/Manuscript/Supp_Tables \\
        --figure 19

--figure is the number of the main figure in the manuscript. Output (default
<supp>/Supplementary_Table_S20_Fig<figure>_Pineda.xlsx) and build_SuppTable_S20_log.txt
beside it. An earlier S20 workbook with a different name is renamed with the suffix _vOri.
Any FAIL stops the build before the workbook is written.

Formatting as the other supplementary tables: Arial 10, bold header row, frozen top row;
Contents with a 12 pt bold title. No fills, borders, or colours.
"""

import argparse
import datetime as dt
import math
import pathlib
import sys

import openpyxl
from openpyxl.styles import Alignment, Font

TABLE = "S20"
SHEETS = ["A_DFG_replication", "B_GAM_cells", "B_GAM_curves", "C_CEPIP_cells", "C_CEPIP_curve",
          "C_CEPIP_IP", "D_Expression_space", "E_Pitcher_cells", "E_Pitcher_gates", "F_Donors",
          "G_Contrasts", "H_LODO_donors", "H_LODO_summary", "I_Power_curves", "I_Power_points"]
EXTRA = ["J_Mapping_agreement", "J_Mapping_probability", "K_NoAuthorQC_bands", "K_NoAuthorQC_donors"]
ALL = SHEETS + EXTRA
DESCRIPTION = {
    "A_DFG_replication": ("Panel A. Spearman r_s and Kendall tau_b of MT-CO3 with each of the eight DFGs, "
                          "the tiers passed, and the r_s expected at this depth (SEA-AD, PVM_3 nuclei "
                          "downsampled to 2,000 UMIs). Tier 1 is r_s >= 0.6 and tau_b >= 0.5. Tier 2 is "
                          "r_s >= 0.5 and tau_b >= 0.4."),
    "B_GAM_cells": ("Panel B. MT-CO3, cDFG, and cHKG9 of each PVM_3 arm cell, as ln(CP10K + 1). The plot "
                    "adds a horizontal jitter of 0.01 to the points."),
    "B_GAM_curves": "Panel B. GAM curves of cDFG and cHKG9 on MT-CO3, with 95% bands and the deviance explained (DE).",
    "C_CEPIP_cells": "Panel C. Cells of the CEP-IP plot and their class (TREP or non-TREP).",
    "C_CEPIP_curve": "Panel C. GAM curve that defines the IP, with its 95% band.",
    "C_CEPIP_IP": "Panel C. IP_x with its 95% CI, IP_y, IPRS and its tier, the TREP counts, and the DE of the IP curve.",
    "D_Expression_space": ("Panel D. Class of each PVM_3 arm cell as plotted (Module cells in the four groups, "
                           "Recovered, Original only, and Excluded), with the CMRF class, module membership, "
                           "MT-CO3, and cDFG. IP_x and IP_y are in sheet C_CEPIP_IP."),
    "E_Pitcher_cells": "Panel E. The same cells in copula space (u_SK and v_SK).",
    "E_Pitcher_gates": "Panel E. SKCS-core gates on u_SK and v_SK.",
    "F_Donors": ("Panel F. Module post-IP% per donor (Module post-IP cells / all Micro-PVM_3 nuclei x 100), "
                 "with diagnosis, age at death, sex (F = 0, M = 1), and PMI, for the donors with at least "
                 "5 Micro-PVM_3 nuclei. FTLD_226 has no age at death and is not in panels G to I."),
    "G_Contrasts": ("Panel G. Diagnosis contrasts: standardized beta of Module post-IP% with 95% CI for FTLD and "
                    "ALS against PN, adjusted for age at death and sex (cov2), and with PMI added (cov3), with BH "
                    "FDR over the four contrasts."),
    "H_LODO_donors": ("Panel H (post hoc). Beta of each contrast with one donor left out (leave-one-donor-out, "
                      "LODO), the change from the full fit, DFBETAS, and Cook's distance."),
    "H_LODO_summary": ("Panel H (post hoc). Observed beta, the number of LODO betas above 0, the share of "
                       "10,000 donor bootstrap resamples (stratified by diagnosis) above 0, and the SEA-AD "
                       "PVM_3 cognitive-status betas (cov4)."),
    "I_Power_curves": ("Panel I (post hoc). Power to detect the SEA-AD NeoC beta against the number of donors. "
                       "The Pineda SE is scaled by the residual df, with the same group ratio and covariates."),
    "I_Power_points": "Panel I (post hoc). Current donors, and donors needed for 80% power.",
    "J_Mapping_agreement": ("MapMyCells labels compared nucleus by nucleus (R13b). 'Repeat run': all R02 kept "
                            "nuclei mapped twice with the same files. 'QC before mapping': the nuclei of the "
                            "authors' QC mapped on their own, against the same nuclei mapped together with all "
                            "R02 kept nuclei. Changed: nuclei whose label differs between the two runs."),
    "J_Mapping_probability": "MapMyCells probabilities of the same comparisons: absolute difference between the two runs.",
    "K_NoAuthorQC_bands": ("PVM_3 arm of the earlier run without the authors' QC (all R02 kept nuclei mapped, "
                           "Micro-PVM_3 by Deep Generative Mapping with MT-CO3 > 0), by %MT band, with the "
                           "nuclei in and removed by the authors' QC (R13c)."),
    "K_NoAuthorQC_donors": "The same arm by donor, with the nuclei above 5%, 10%, and 20% MT and in the authors' QC.",
}
TITLE = "Supplementary Table S20. External validation 1 (Pineda et al. 2024, BA9)."

FONT = Font(name="Arial", size=10)
FONT_B = Font(name="Arial", size=10, bold=True)
FONT_T = Font(name="Arial", size=12, bold=True)
DASHES = ("–", "—")
TOL = 1e-9


class Log:
    def __init__(self):
        self.lines, self.fail = [], 0

    def say(self, s=""):
        self.lines.append(s)
        print(s)

    def check(self, ok, what):
        self.say(f"  [{'PASS' if ok else 'FAIL'}] {what}")
        self.fail += 0 if ok else 1
        return ok


def read_rows(path, sheet):
    """All rows of a sheet as lists. Workbooks written by openxlsx carry no reliable
    dimension, so the dimension is reset first."""
    wb = openpyxl.load_workbook(path, read_only=True)
    ws = wb[sheet]
    ws.reset_dimensions()
    rows = [list(r) for r in ws.iter_rows(values_only=True)]
    wb.close()
    while rows and all(v is None for v in rows[-1]):
        rows.pop()
    return rows


def sheet_names(path):
    wb = openpyxl.load_workbook(path, read_only=True)
    out = wb.sheetnames
    wb.close()
    return out


def records(rows):
    """Rows as dicts keyed by the header."""
    h = rows[0]
    return [dict(zip(h, r)) for r in rows[1:]]


def num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return math.nan


def same(a, b):
    a, b = num(a), num(b)
    return (math.isnan(a) and math.isnan(b)) or abs(a - b) <= TOL


def as_bool(v):
    return v in (True, 1, "TRUE", "True", "true")


def guardrail_ok(path):
    """Finished, and the FAIL count is 0."""
    txt = path.read_text(encoding="utf-8", errors="replace").splitlines()
    fails = [x for x in txt if x.strip().startswith("FAIL") and ":" in x]
    n = fails[-1].split(":")[1].strip() if fails else "not found"
    finished = not any(x.strip().startswith("Result: the run did not finish") for x in txt)
    return finished and n == "0", n


def main():
    ap = argparse.ArgumentParser(description="Build Supplementary Table S20 from the R16 source data.")
    ap.add_argument("--root", required=True, help="External_Dataset_1 folder")
    ap.add_argument("--supp", required=True, help="Supp_Tables folder")
    ap.add_argument("--figure", required=True, help="number of the main figure, for example 19")
    ap.add_argument("--out", default=None,
                    help="output xlsx (default: <supp>/Supplementary_Table_S20_Fig<figure>_Pineda.xlsx)")
    args = ap.parse_args()
    root, supp = pathlib.Path(args.root), pathlib.Path(args.supp)
    fig = str(args.figure).strip()
    out = pathlib.Path(args.out) if args.out else supp / f"Supplementary_Table_{TABLE}_Fig{fig}_Pineda.xlsx"
    subtitle = (f"Values behind Figure {fig}. Sheet letters follow the panels of Figure {fig}. "
                "Panels H and I are post hoc analyses. Sheets J and K support the Methods.")
    F = {"r16": root / "R16_Figure_Pineda" / "R16_Pineda_figure_source_data.xlsx",
         "r16_gr": root / "R16_Figure_Pineda" / "R16_Pineda_Guardrail_Log.txt",
         "r14": root / "R14_core" / "R14_Pineda_core.xlsx",
         "r14_gr": root / "R14_core" / "R14_Pineda_Guardrail_Log.txt",
         "r15": root / "R15_contrasts" / "R15_Pineda_contrasts.xlsx",
         "r15_gr": root / "R15_contrasts" / "R15_Pineda_Guardrail_Log.txt",
         "r17": root / "R17_posthoc" / "R17_Pineda_posthoc_beta.xlsx",
         "r17_gr": root / "R17_posthoc" / "R17_Pineda_Guardrail_Log.txt",
         "r13b_rep": root / "R13b_mapping_comparison_repeat" / "R13b_Pineda_mapping_comparison.xlsx",
         "r13b_qc": root / "R13b_mapping_comparison_QC" / "R13b_Pineda_mapping_comparison.xlsx",
         "r13c": root / "R13c_noAuthorQC_arm" / "R13c_Pineda_noAuthorQC_arm.xlsx"}
    L = Log()
    L.say(f"build_SuppTable_S20.py | {dt.datetime.now():%Y-%m-%d %H:%M:%S}")
    L.say(f"Root:   {root}")
    L.say(f"Output: {out}")
    L.say(f"Figure: {fig}")
    L.say("")
    log_path = out.with_name("build_SuppTable_S20_log.txt")

    def stop():
        L.say(f"VERDICT: FAIL | {L.fail} check(s) failed; workbook not written")
        log_path.parent.mkdir(parents=True, exist_ok=True)
        log_path.write_text("\n".join(L.lines) + "\n", encoding="utf-8")
        sys.exit(1)

    # ---- inputs ---------------------------------------------------------------------------
    L.say("SECTION A -- inputs")
    miss = [str(p) for p in F.values() if not p.exists()]
    L.check(not miss, "all input files present" + "".join(f"\n      missing: {m}" for m in miss))
    L.check(supp.is_dir(), f"--supp folder exists ({supp})")
    L.check(fig.isdigit(), f"--figure is a number ({fig})")
    if L.fail:
        stop()
    for k in ("r14_gr", "r15_gr", "r17_gr", "r16_gr"):
        ok, n = guardrail_ok(F[k])
        L.check(ok, f"{F[k].name}: finished, FAIL count {n}")
    names = sheet_names(F["r16"])
    L.check(names == SHEETS, f"R16 source data has the {len(SHEETS)} expected sheets in order"
            + ("" if names == SHEETS else f" | found {names}"))
    if L.fail:
        stop()
    data = {s: read_rows(F["r16"], s) for s in SHEETS}
    # J: the two R13b comparisons, one row block each, with the comparison named in the first column.
    for sheet, src_sheet in (("J_Mapping_agreement", "Agreement"), ("J_Mapping_probability", "Probability")):
        rows = []
        for label, key in (("Repeat run", "r13b_rep"), ("QC before mapping", "r13b_qc")):
            r = read_rows(F[key], src_sheet)
            if not rows:
                rows.append(["comparison"] + r[0])
            rows += [[label] + x for x in r[1:]]
        data[sheet] = rows
    data["K_NoAuthorQC_bands"] = read_rows(F["r13c"], "Bands")
    data["K_NoAuthorQC_donors"] = read_rows(F["r13c"], "Donors")
    for s in ALL:
        L.say(f"  {s:20} {len(data[s]) - 1:6} rows, {len(data[s][0])} columns")
    L.say("")

    # ---- checks against the source workbooks ----------------------------------------------------
    L.say("SECTION B -- values agree with R14, R15, and R17")
    R = {s: records(data[s]) for s in SHEETS}

    dfg14 = {r["DFG"]: r for r in records(read_rows(F["r14"], "02_DFGs"))}
    bad = [r["DFG"] for r in R["A_DFG_replication"]
           if r["DFG"] not in dfg14 or not (same(r["rs"], dfg14[r["DFG"]]["rs"])
                                           and same(r["tau_b"], dfg14[r["DFG"]]["tau_b"])
                                           and as_bool(r["passes_tier1"]) == as_bool(dfg14[r["DFG"]]["passes_tier1"])
                                           and as_bool(r["passes_tier2"]) == as_bool(dfg14[r["DFG"]]["passes_tier2"]))]
    L.check(len(R["A_DFG_replication"]) == 8 and not bad,
            f"A: 8 DFGs, r_s, tau_b, and tiers equal R14 02_DFGs" + (f" | differ: {bad}" if bad else ""))

    ids = {s: [r["cell_id"] for r in R[s]] for s in ("B_GAM_cells", "C_CEPIP_cells", "D_Expression_space",
                                                     "E_Pitcher_cells")}
    n_cells = len(ids["B_GAM_cells"])
    L.check(len(set(ids["B_GAM_cells"])) == n_cells, f"B: {n_cells:,} cells, cell IDs unique")
    L.check(all(set(v) == set(ids["B_GAM_cells"]) for v in ids.values()),
            "B, C, D, and E hold the same cells")
    cls_d = {r["cell_id"]: r["class"] for r in R["D_Expression_space"]}
    L.check(all(cls_d.get(r["cell_id"]) == r["class"] for r in R["E_Pitcher_cells"]),
            "D and E give each cell the same class")
    ipr = R["C_CEPIP_IP"][0] if len(R["C_CEPIP_IP"]) == 1 else {}
    n_trep = sum(r["class"] == "TREP" for r in R["C_CEPIP_cells"])
    L.check(len(R["C_CEPIP_IP"]) == 1 and same(ipr.get("TREP"), n_trep)
            and same(ipr.get("non_TREP"), len(R["C_CEPIP_cells"]) - n_trep),
            f"C: TREP {n_trep:,} and non-TREP {len(R['C_CEPIP_cells']) - n_trep:,} equal the counts in C_CEPIP_IP")
    L.check(sorted({r["set"] for r in R["B_GAM_curves"]}) == ["cDFG", "cHKG9"], "B: curves for cDFG and cHKG9")
    L.check(len(R["E_Pitcher_gates"]) == 2, "E: two SKCS-core gates")

    don15 = {r["donor"]: r for r in records(read_rows(F["r15"], "04_Donors_PVM_3"))}
    bad = [r["donor"] for r in R["F_Donors"]
           if r["donor"] not in don15 or not same(r["pct_PostIP"], don15[r["donor"]]["pct_PostIP"])]
    L.check(len(R["F_Donors"]) == len(don15) and not bad,
            f"F: {len(R['F_Donors'])} donors, Module post-IP% equal to R15 04_Donors_PVM_3"
            + (f" | differ: {bad}" if bad else ""))

    con15 = {r["contrast"]: r for r in records(read_rows(F["r15"], "01_Contrasts")) if r["arm"] == "PVM_3"}
    cols = ("n", "beta_std", "se", "ci_low", "ci_high", "p_value", "FDR")
    bad = [r["contrast"] for r in R["G_Contrasts"]
           if r["contrast"] not in con15 or not all(same(r[c], con15[r["contrast"]][c]) for c in cols)]
    L.check(len(R["G_Contrasts"]) == 4 and not bad, "G: 4 contrasts, beta, SE, CI, p, and FDR equal R15 01_Contrasts"
            + (f" | differ: {bad}" if bad else ""))

    lodo17 = {(r["contrast"], r["donor"]): r for r in records(read_rows(F["r17"], "05_LODO_donors"))}
    bad = [(r["contrast"], r["donor"]) for r in R["H_LODO_donors"]
           if (r["contrast"], r["donor"]) not in lodo17
           or not same(r["beta_without"], lodo17[(r["contrast"], r["donor"])]["beta_without"])]
    L.check(len(R["H_LODO_donors"]) == len(lodo17) and not bad,
            f"H: {len(lodo17)} LODO betas equal R17 05_LODO_donors" + (f" | differ: {bad[:5]}" if bad else ""))
    sum17 = {r["contrast"]: r for r in records(read_rows(F["r17"], "04_LODO_summary"))}
    boot17 = {r["contrast"]: r for r in records(read_rows(F["r17"], "06_Bootstrap"))}
    sea17 = {r["set"]: r for r in records(read_rows(F["r17"], "10_SEAAD_reference"))}
    bad = [r["contrast"] for r in R["H_LODO_summary"]
           if not (r["contrast"] in sum17 and same(r["lodo_positive"], sum17[r["contrast"]]["lodo_positive"])
                   and same(r["beta_std"], sum17[r["contrast"]]["beta_std"])
                   and same(r["bootstrap_share_above_0"], boot17[r["contrast"]]["share_above_0"])
                   and same(r["seaad_beta_NeoC"], sea17["NeoC"]["beta"])
                   and same(r["seaad_beta_NonNeoC"], sea17["NonNeoC"]["beta"]))]
    L.check(len(R["H_LODO_summary"]) == 4 and not bad,
            "H: LODO counts, bootstrap shares, and SEA-AD betas equal R17" + (f" | differ: {bad}" if bad else ""))

    need17 = {r["contrast"]: r for r in records(read_rows(F["r17"], "09_Donors_needed")) if r["seaad_set"] == "NeoC"}
    pts = [r for r in R["I_Power_points"] if r["point"] == "donors for 80% power"]
    bad = [r["contrast"] for r in pts
           if r["contrast"] not in need17 or not same(r["donors"], need17[r["contrast"]]["donors_needed"])]
    L.check(len(pts) == 4 and not bad, "I: donors needed for 80% power equal R17 09_Donors_needed (NeoC beta)"
            + (f" | differ: {bad}" if bad else ""))
    L.check(all(same(r["target_beta"], sea17["NeoC"]["beta"]) for r in R["I_Power_curves"] + R["I_Power_points"]),
            "I: the target beta is the SEA-AD NeoC beta")

    RJ = records(data["J_Mapping_agreement"])
    rep_hm = [r for r in RJ if r["comparison"] == "Repeat run" and r["algorithm"] == "HM" and r["level"] == "supertype"]
    qc_dgm = [r for r in RJ if r["comparison"] == "QC before mapping" and r["algorithm"] == "DGM" and r["level"] == "supertype"]
    L.check(len(rep_hm) == 1 and len(qc_dgm) == 1, "J: both comparisons present, with DGM and HM rows")
    if rep_hm and qc_dgm:
        L.say(f"  J: repeat run HM supertypes the same for {num(rep_hm[0]['pct_same']):.4f}% | QC before mapping "
              f"DGM supertypes changed for {100 - num(qc_dgm[0]['pct_same']):.2f}%")
    kb, kd = records(data["K_NoAuthorQC_bands"]), records(data["K_NoAuthorQC_donors"])
    n_arm = sum(num(r["arm_nuclei"]) for r in kb)
    n20 = sum(num(r["arm_nuclei"]) for r in kb if str(r["pct_mt_band"]).startswith(">20"))
    L.check(same(n_arm, sum(num(r["arm_nuclei"]) for r in kd)), f"K: bands and donors hold the same {n_arm:,.0f} arm nuclei")
    L.check(same(n_arm, 2795) and same(n20, 239),
            f"K: {n_arm:,.0f} arm nuclei and {n20:,.0f} above 20% MT, as cited in the Methods (2,795 and 239)")
    L.check(set(DESCRIPTION) == set(ALL), "every sheet has a description")

    txt = [v for rows in data.values() for r in rows for v in r if isinstance(v, str)]
    txt += [TITLE, subtitle] + list(DESCRIPTION.values())
    L.check(not [v for v in txt if any(d in v for d in DASHES)], "no en dash or em dash in the text")
    L.say("")
    if L.fail:
        stop()

    # ---- earlier S20 workbooks ------------------------------------------------------------------
    L.say("SECTION C -- earlier S20 workbooks")
    old = [p for p in sorted(supp.glob(f"Supplementary_Table_{TABLE}_*.xlsx"))
           if p.resolve() != out.resolve() and "_vOri" not in p.stem and not p.name.startswith("~$")]
    for p in old:
        v = p.with_name(p.stem + "_vOri.xlsx")
        if L.check(not v.exists(), f"{v.name} is not already taken"):
            p.rename(v)
            L.say(f"  {p.name} -> {v.name}")
    if not old:
        L.say("  none")
    if out.exists():
        L.say(f"  {out.name} exists and is replaced")
    L.say("")
    if L.fail:
        stop()

    # ---- write ----------------------------------------------------------------------------------
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "Contents"
    ws["A1"] = TITLE
    ws["A1"].font = FONT_T
    ws["A2"] = subtitle
    ws["A2"].font = FONT
    for j, hdr in enumerate(["Sheet", "Rows", "Description"], start=1):
        ws.cell(4, j, hdr).font = FONT_B
    for i, s in enumerate(ALL, start=5):
        for j, v in enumerate([s, len(data[s]) - 1, DESCRIPTION[s]], start=1):
            c = ws.cell(i, j, v)
            c.font = FONT
            c.alignment = Alignment(vertical="top", wrap_text=(j == 3))
        ws.cell(i, 2).number_format = "#,##0"
    ws.column_dimensions["A"].width = 26
    ws.column_dimensions["B"].width = 9
    ws.column_dimensions["C"].width = 110
    ws.freeze_panes = "A5"
    for s in ALL:
        w = wb.create_sheet(s)
        for r in data[s]:
            w.append(r)
        for row in w.iter_rows():
            for c in row:
                c.font = FONT
        for c in w[1]:
            c.font = FONT_B
        w.freeze_panes = "A2"
    wb.save(out)

    # ---- read back ------------------------------------------------------------------------------
    L.say("SECTION D -- read back")
    back = openpyxl.load_workbook(out, read_only=True)
    L.check(back.sheetnames == ["Contents"] + ALL, "output sheets in order")
    nd = 0
    for s in ALL:
        rows = [list(r) for r in back[s].iter_rows(values_only=True)]
        for a, b in zip(rows, data[s]):
            nd += sum(1 for u, v in zip(a, b) if not (u == v or (u is None and v in (None, ""))))
        nd += abs(len(rows) - len(data[s]))
    back.close()
    L.check(nd == 0, f"every sheet reads back as written ({sum(len(d) - 1 for d in data.values()):,} data rows)")
    L.say("")
    L.say(f"VERDICT: {'FAIL' if L.fail else 'PASS'} | {out.name}")
    log_path.write_text("\n".join(L.lines) + "\n", encoding="utf-8")
    sys.exit(1 if L.fail else 0)


if __name__ == "__main__":
    main()
