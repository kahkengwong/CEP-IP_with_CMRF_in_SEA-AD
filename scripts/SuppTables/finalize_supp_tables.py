#!/usr/bin/env python3
"""
finalize_supp_tables.py

Writes Supplementary Tables S1 to S19 in their published form from the
workbooks written by the builders and the analysis scripts. No value is
computed here. For each table the script:
  - opens the source workbook and keeps the listed sheets, renamed and in
    the published order (other sheets, such as Overview or Notes, are dropped);
  - applies the listed text edits (column names, labels and wording);
  - adds the Contents sheet: title, subtitle and sheet descriptions from
    supp_tables_contents.json, with the number of data rows counted here;
  - saves Manuscript/Supp_Tables/Supplementary_Table_S<n>_<name>.xlsx.

Formula cells (the verdict and class columns of S13 to S17) stay formulas.
Excel computes their values when the file is opened.

Every edit states how many cells it must change. If a sheet or an edit
target is not found as expected, that table is not written and the run
ends with VERDICT: FAIL. The other tables are still written.

Usage:
    python finalize_supp_tables.py                 # all tables with a spec
    python finalize_supp_tables.py --tables 12
    python finalize_supp_tables.py --tables 13-17,19
Root: --root, else $SEAAD_ROOT, else the current folder.
Log:  <root>/Manuscript/Supp_Tables/finalize_supp_tables_log.txt
"""

import argparse
import io
import json
import os
import pathlib
import posixpath
import re
import sys
import zipfile

import openpyxl
from openpyxl.styles import Alignment, Font

HERE = pathlib.Path(__file__).resolve().parent
CONTENTS_JSON = HERE / "supp_tables_contents.json"
OUT_SUB = pathlib.Path("Manuscript") / "Supp_Tables"

F10 = Font(name="Arial", size=10)
F10B = Font(name="Arial", size=10, bold=True)
F12B = Font(name="Arial", size=12, bold=True)


##### Table specifications #####
# src     source workbook, relative to the root
# sheets  (source sheet, published sheet) in the published order
# restyle True: Arial 10 in every data cell, bold header row, top row frozen
# edits   text edits, applied in order after the sheets are renamed:
#   ("header",  sheet, old, new)                one column name
#   ("replace", sheets, old, new, exact, n)     cell text; exact=True matches the
#                                               whole cell, False a substring; n is
#                                               the number of cells that must change
#                                               (None: any number, including 0)
#   ("formula", sheets, old, new, n_min)        text inside formulas, at least n_min cells
#   ("map",     sheets, column, mapping, keep)  every value of the column must be a
#                                               key of mapping or start with a prefix
#                                               in keep
#   ("strip_root", sheet, column)               removes the absolute root from paths
# width_a, width_d  column widths of the Contents sheet (A default 30; D only when set)

S13_SHEETS = [("Fig15A_readout", "A_Fig12A_readout"), ("Fig15A_stats", "B_Fig12A_stats"),
              ("Fig15B_readout", "C_Fig12B_readout"), ("Fig15B_removed", "D_Fig12B_removed"),
              ("Fig15B_stats", "E_Fig12B_stats"), ("Fig15AB_donors", "F_Fig12AB_donors"),
              ("Fig15C_stats", "G_Fig12C_stats"), ("Fig15D_stats", "H_Fig12D_stats"),
              ("Fig15D_tertiles", "I_Fig12D_tertiles"), ("Fig15E_stats", "J_Fig12E_stats"),
              ("Fig15E_samples", "K_Fig12E_samples"), ("Fig15F_stats", "L_Fig12F_stats"),
              ("Fig15F_samples", "M_Fig12F_samples"), ("Fig15G_stats", "N_Fig12G_stats"),
              ("Fig15G_tertiles", "O_Fig12G_tertiles"), ("Fig15H_stats", "P_Fig12H_stats"),
              ("SuppFig11_stats", "Q_SuppFig15_stats")]

S14_SHEETS = [("Fig16A_readout", "A_Fig13A_readout"), ("Fig16A_stats", "B_Fig13A_stats"),
              ("Fig16B_stats", "C_Fig13B_stats"), ("Fig16C_readout", "D_Fig13C_readout"),
              ("Fig16C_stats", "E_Fig13C_stats"), ("Fig16D_stats", "F_Fig13D_stats"),
              ("Fig16D_donors", "G_Fig13D_donors"), ("Fig16E_readout", "H_Fig13E_readout"),
              ("Fig16E_stats", "I_Fig13E_stats"), ("Fig16F_stats", "J_Fig13F_stats"),
              ("Fig16G_conditions", "K_Fig13G_conditions"), ("Fig16G_stats", "L_Fig13G_stats"),
              ("Fig16_units", "M_Fig13_units")]

S15_SHEETS = [("Fig17A_readout", "A_Fig14A_readout"), ("Fig17A_stats", "B_Fig14A_stats"),
              ("Fig17B_stats", "C_Fig14B_stats"), ("Fig17C_verdicts", "D_Fig14C_verdicts"),
              ("Fig17C_stats", "E_Fig14C_stats"), ("Fig17D_verdicts", "F_Fig14D_verdicts"),
              ("Fig17D_stats", "G_Fig14D_stats"), ("Fig17D_donors", "H_Fig14D_donors"),
              ("Fig17E_readout", "I_Fig14E_readout"), ("Fig17E_stats", "J_Fig14E_stats"),
              ("FigS13_stats", "K_SuppFig17_stats"), ("FigS13_measures", "L_SuppFig17_measures"),
              ("FigS13_counts", "M_SuppFig17_counts"), ("Fig17_units", "N_Fig14_units")]

S16_SHEETS = [("Fig18A_verdicts", "A_Fig15A_verdicts"), ("Fig18A_stats", "B_Fig15A_stats"),
              ("Fig18A_panels", "C_Fig15A_panels"), ("Fig18B_readout", "D_Fig15B_readout"),
              ("Fig18B_stats", "E_Fig15B_stats"), ("Fig18C_readout", "F_Fig15C_readout"),
              ("Fig18C_stats", "G_Fig15C_stats"), ("Fig18C_donors", "H_Fig15C_donors"),
              ("Fig18BC_cells", "I_Fig15BC_cells"), ("Fig18D_stats", "J_Fig15D_stats"),
              ("Fig18D_curves", "K_Fig15D_curves"), ("Fig18E_stats", "L_Fig15E_stats"),
              ("Fig18E_regions", "M_Fig15E_regions"), ("Fig18DE_units", "N_Fig15DE_units")]

S17_SHEETS = [("Fig19A_readout", "A_Fig16A_readout"), ("Fig19A_16A_stats", "B_Fig16A_13A_stats"),
              ("Fig19A_16B_stats", "C_Fig16A_13B_stats"), ("Fig19A_16F_stats", "D_Fig16A_13F_stats"),
              ("Fig19B_counts", "E_Fig16B_counts"),
              ("Fig19B_readout_histology", "F_Fig16B_readout_histology"),
              ("Fig19B_readout_clinical", "G_Fig16B_readout_clinical"),
              ("Fig19B_15A_stats", "H_Fig16B_12A_stats"), ("Fig19B_17A_stats", "I_Fig16B_14A_stats"),
              ("Fig19B_17B_stats", "J_Fig16B_14B_stats"), ("Fig19B_17C_stats", "K_Fig16B_14C_stats"),
              ("Fig19B_17D_stats", "L_Fig16B_14D_stats"), ("Fig19B_18A_stats", "M_Fig16B_15A_stats"),
              ("Fig19C_clinical", "N_Fig16C_clinical"), ("Verdicts_17C", "O_Verdicts_14C"),
              ("Verdicts_17D", "P_Verdicts_14D"), ("Verdicts_18A", "Q_Verdicts_15A"),
              ("Influence", "R_Influence"), ("Influence_summary", "S_Influence_summary"),
              ("SA_donors", "T_SA_donors"), ("Frames", "U_Frames"), ("Null_draws", "V_Null_draws")]

S18_SHEETS = [("Fig1A", "A_M4share_cognition"), ("Fig1B", "B_State_vs_MT"),
              ("Fig1C", "C_State_MT_verdicts"), ("Fig1D", "D_Fraction_MT_MMSE"),
              ("Fig2A", "E_Random_cognition"), ("Fig2B", "F_Random_MT"), ("Fig2C", "G_Random_verdicts"),
              ("T1_Models_13", "H_Models"), ("T2_Verdicts_13", "I_Verdicts"),
              ("T3_Bootstrap_13", "J_Bootstrap"), ("T4_Decomposition_13", "K_Decomposition"),
              ("T5_Collinearity_13", "L_Collinearity"), ("T6_PostHoc_13", "M_PostHoc"),
              ("T7_Specs_13", "N_Specs"), ("T8_RandomSubset_13b", "O_RandomSubset"),
              ("T9_LabelAgreement_13b", "P_LabelAgreement"), ("T10_Design_13b", "Q_Design"),
              ("T11_Donors", "R_Donors"), ("T12_Readings", "S_Readings")]

S19_SHEETS = [("Tests", "A_Tests"), ("Cells", "B_Cells"), ("Lead_result", "C_Lead_result"),
              ("Lead_tally", "D_Lead_tally"), ("Both_region_tally", "E_Both_region_tally"),
              ("Boxed", "F_Boxed"), ("Models", "G_Models"), ("Tallies", "H_Tallies"),
              ("Checks", "I_Checks")]

# Fig17_Summary source labels (analysis-run names) -> published table and sheet.
# Labels of S12 are already in published form and pass unchanged.
S19_SOURCE = {
    "Fig15 Fig15A_readout": "S13 A_Fig12A_readout",
    "Fig15 Fig15B_readout": "S13 C_Fig12B_readout",
    "Fig16 Fig16A_readout": "S14 A_Fig13A_readout",
    "Fig16 Fig16B_stats": "S14 C_Fig13B_stats",
    "Fig16 Fig16E_readout": "S14 H_Fig13E_readout",
    "Fig16 Fig16F_stats": "S14 J_Fig13F_stats",
    "Fig17 Fig17D_verdicts": "S15 F_Fig14D_verdicts",
    "Fig18 Fig18A_verdicts": "S16 A_Fig15A_verdicts",
    "Fig18 Fig18B_readout": "S16 D_Fig15B_readout",
    "Fig18 Fig18C_readout": "S16 F_Fig15C_readout",
    "Fig18 Fig18B_readout; Fig18C_readout": "S16 D_Fig15B_readout; F_Fig15C_readout",
    "Fig19 Fig19A_readout": "S17 A_Fig16A_readout",
    "Fig19 Fig19C_clinical": "S17 N_Fig16C_clinical",
    "Part13 C_Models (S03, S07)": "S18 H_Models (S03, S07)",
    "Part13b D_Summary (bA)": "S18 O_RandomSubset (bA)",
}

S12_SHEETS = [(s, s) for s in [
    "A_donor_frame", "B_Fig9A_variance", "C_Fig9B_cognition", "D_Fig9C_pathology",
    "E_Fig9D_mod_vs_nonmod", "F_Fig10A_benchmark", "G_Fig10B_competitors", "H_Fig10CD_decomposition",
    "I_Fig11A_boundary_assoc", "J_Fig11B_boundary_curves", "K_Fig11BC_cells", "L_Fig11C_scorecard",
    "M_Fig11C_tests", "N_Fig11C_overlap", "O_Fig11BC_boundary_summary", "P_SuppFig14_cov4_cov7",
    "Q_SuppFig14_leave_one_in", "R_SuppFig14_without_SA", "S_SuppFig14_tissue_disease"]]

BUILD = "Manuscript/Supp_Tables/build/"

S1_SHEETS = [(s, s) for s in [
    "A_cutoff_rule_verdict", "B_marker_dilution", "C_mixture", "D_bootstrap_stability", "E_arm_matching",
    "F_scorer_overlap", "G_selection_by_supertype", "H_composition_by_scorer", "I_method_comparison",
    "J_sweep_summary"]]
S2_SHEETS = [(s, s) for s in ["A_deg_summary", "B_hypergeometric", "C_mic13_gene_detail"]]

S3_SHEETS = [(s, s) for s in ["A_mito_gene_summary", "B_region_screen", "C_pooled_screen_NeoC",
                               "D_pooled_screen_NonNeoC", "E_GOI_composite"]]
S4_SHEETS = [(s, s) for s in ["A_population_tally", "B_dfg_counts", "C_top_GOIs", "D_pooled_cross_population",
                               "E_pooled_tally"]]

S5_SHEETS = [(s, s) for s in [
    "A_mito_gene_inventory", "B_trend_fits", "C_within_donor_rho", "D_errata_specimen_audit",
    "E_matched_panel_null", "F_per_library", "G_ambient_mediation", "H_ambient_index_quartiles",
    "I_per_region", "J_region_LORO", "K_supertype_specificity", "L_segmented_screen",
    "M_median_specialness", "N_ordinal_levels", "O_donor_module_rho", "P_shape_tests"]]

S6_SHEETS = [(s, s) for s in ["A_contamination_summary", "B_gene_level", "C_gene_set_summary",
                               "D_module_diagnostic", "E_per_nucleus"]]

S7_SHEETS = [(s, s) for s in ["A_cells", "B_HKG_pitcher_cells", "C_GAM", "D_Si_mixture",
                               "E_Si_boundaries", "F_Si_posterior", "G_Si_rank_tests", "H_pitcher_summary"]]
S8_SHEETS = [(s, s) for s in ["A_IP_detection", "B_Rx_grid"]]
S9_SHEETS = [(s, s) for s in ["A_copula_IP", "B_copula_gates", "C_worked_cells"]]

S10_SHEETS = [(s, s) for s in ["A_cells", "B_band", "C_gate", "D_gate_trajectory", "E_fold_RMSE",
                                "F_candidates", "G_tests", "H_resid_by_class", "I_Fig8_summary", "J_Fig8_tests"]]

S11_SHEETS = [(s, s) for s in [
    "A_Contrasts", "B_Control_test", "B_Summary", "C_Axis_check", "D_Control_validity",
    "E_Survivors_regions", "F_Neuronal_summary", "F_Neuronal_genes", "G_UMI_summary", "G_UMI_donors",
    "G_IP_split", "H_Panel_summary", "H_Panel_donors"]]

TABLES = {
    1: dict(src=BUILD + "SuppTable_S1_CellSelection.xlsx", sheets=S1_SHEETS, restyle=True, width_a=28, edits=[]),
    2: dict(src=BUILD + "SuppTable_S2_Mic13Validation.xlsx", sheets=S2_SHEETS, restyle=True, width_a=28, edits=[
        ("replace", ["A_deg_summary"], ">=", "\u2265", False, 1),
        ("replace", ["B_hypergeometric"], "(>=4 regions)", "(\u22654 regions)", False, 2),
    ]),
    3: dict(src=BUILD + "SuppTable_S3_MT-CO3_selection.xlsx", sheets=S3_SHEETS, restyle=True, width_a=28, edits=[
        ("header", "A_mito_gene_summary", "regions_listed_top_in_Part6", "regions_listed_top_S4B"),
    ]),
    4: dict(src=BUILD + "SuppTable_S4_GOI_DFG_screen.xlsx", sheets=S4_SHEETS, restyle=True, width_a=28, edits=[]),
    5: dict(src=BUILD + "SuppTable_S5_PreChecks.xlsx", sheets=S5_SHEETS, restyle=True, width_a=28, edits=[
        ("replace", ["B_trend_fits"], "analysed", "analyzed", True, 1),
        ("header", "G_ambient_mediation", "partialled_on", "partialed_on"),
    ]),
    6: dict(src=BUILD + "SuppTable_S6_DecontX.xlsx", sheets=S6_SHEETS, width_a=28, edits=[]),
    7: dict(src=BUILD + "SuppTable_S7_Fig2_SKCS.xlsx", sheets=S7_SHEETS, edits=[]),
    8: dict(src=BUILD + "SuppTable_S8_Fig3_IP.xlsx", sheets=S8_SHEETS, edits=[]),
    9: dict(src=BUILD + "SuppTable_S9_Fig4_Copula.xlsx", sheets=S9_SHEETS, edits=[]),
    10: dict(src=BUILD + "SuppTable_S10_Fig5_CMRF.xlsx", sheets=S10_SHEETS, edits=[
        ("replace", ["C_gate"], "dRMSE <= BV", "dRMSE \u2264 BV", True, 8),
    ]),
    11: dict(src=BUILD + "SuppTable_S11_DEGs_ModvsNonMod.xlsx", sheets=S11_SHEETS, restyle=True, edits=[
        ("header", "G_IP_split", "part1_module_post", "module_post_de_run"),
        ("header", "G_IP_split", "part1_module_pre", "module_pre_de_run"),
        ("header", "G_IP_split", "part1_nonmodule_post", "nonmodule_post_de_run"),
        ("header", "G_IP_split", "part1_nonmodule_pre", "nonmodule_pre_de_run"),
        ("header", "G_IP_split", "agrees_with_part1", "agrees_with_de_run"),
    ]),
    12: dict(src=BUILD + "SuppTable_S12_Prognosis_&_Ablation.xlsx", sheets=S12_SHEETS, width_d=8, edits=[
        ("replace", ["L_Fig11C_scorecard"], "three real groups", "three main arms", True, 4),
        ("replace", ["G_Fig10B_competitors"], "rival only", "competitor only", True, 57),
        ("replace", ["G_Fig10B_competitors"], "not readable", "unreadable", True, 119),
    ]),
    13: dict(src="Part12_rerun/07_Fig15/SuppTable_Fig15.xlsx", sheets=S13_SHEETS, edits=[]),
    14: dict(src="Part12_rerun/08_Fig16/SuppTable_Fig16.xlsx", sheets=S14_SHEETS, edits=[
        ("header", "B_Fig13A_stats", "lo", "null_median"),
        ("header", "B_Fig13A_stats", "hi", "null_q95"),
    ]),
    15: dict(src="Part12_rerun/09_Fig17/SuppTable_Fig17.xlsx", sheets=S15_SHEETS, edits=[]),
    16: dict(src="Part12_rerun/10_Fig18/SuppTable_Fig18.xlsx", sheets=S16_SHEETS, edits=[]),
    17: dict(src="Part12_rerun/11_Fig19/SuppTable_Fig19.xlsx", sheets=S17_SHEETS, edits=[
        ("formula", ["A_Fig16A_readout", "F_Fig16B_readout_histology", "G_Fig16B_readout_clinical"],
         "without the eleven", "without the 11", 1),
        ("header", "B_Fig16A_13A_stats", "lo", "null_median"),
        ("header", "B_Fig16A_13A_stats", "hi", "null_q95"),
        ("header", "A_Fig16A_readout", "null_est_lo", "null_est_at_q_tail"),
        ("header", "A_Fig16A_readout", "null_est_hi", "null_est_at_q95"),
    ]),
    18: dict(src="Part13b_M4RandomSubset/Part13_13b_CP10K_M4_SuppTable.xlsx", sheets=S18_SHEETS,
             restyle=True, edits=[
        ("replace", ["S_Readings"], "Pre-set reading", "Pre-specified reading", True, 1),
        ("replace", ["S_Readings"], "Pre-set text", "Pre-specified text", True, 2),
        ("replace", ["S_Readings"], "POST-HOC reading (fraction vs %MT, NeoC MMSE, p1)",
         "Post-hoc reading (fraction vs %MT, NeoC MMSE, p1)", True, 1),
    ]),
    19: dict(src="Manuscript/Fig17/Preprint/Fig17_Summary.xlsx", sheets=S19_SHEETS, restyle=True, edits=[
        ("map", ["A_Tests", "B_Cells", "G_Models"], "source", S19_SOURCE, ("S12 ",)),
        ("strip_root", "I_Checks", "detail"),
        ("replace", ["I_Checks"], "C9 9A vs 12A unadjusted (open item c)", "C9 9A vs 12A unadjusted", True, 1),
        ("replace", ["B_Cells", "G_Models"], "both", "Both", True, None),
        ("replace", ["A_Tests", "B_Cells"], ">=", "\u2265", False, None),
        ("replace", ["A_Tests", "B_Cells"], "'not readable' = unreadable", "'unreadable' = unreadable", False, None),
    ]),
}


##### Log #####
class Log:
    def __init__(self):
        self.lines, self.fails = [], 0

    def say(self, s=""):
        self.lines.append(s)
        print(s, flush=True)

    def chk(self, ok, what):
        self.say(f"    [{'PASS' if ok else 'FAIL'}] {what}")
        if not ok:
            self.fails += 1
        return ok


class TableFail(Exception):
    pass


##### Helpers #####
def load_book(L, path):
    """Opens a workbook. Some openxlsx versions write a relationship to a
    drawing part that is not in the file; openpyxl cannot open such a file, so
    those relationships are dropped from an in-memory copy first."""
    try:
        return openpyxl.load_workbook(path)
    except KeyError:
        pass
    zin = zipfile.ZipFile(path)
    names = set(zin.namelist())
    buf, n = io.BytesIO(), 0
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zout:
        for item in zin.infolist():
            data = zin.read(item.filename)
            if item.filename.endswith(".rels"):
                base = posixpath.dirname(posixpath.dirname(item.filename))
                txt = data.decode("utf-8")

                def keep(m):
                    nonlocal n
                    tgt = re.search(r'Target="([^"]+)"', m.group(0)).group(1)
                    full = posixpath.normpath(posixpath.join(base, tgt)).lstrip("/")
                    if 'TargetMode="External"' in m.group(0) or full in names:
                        return m.group(0)
                    n += 1
                    return ""
                data = re.sub(r"<Relationship [^>]*/>", keep, txt).encode("utf-8")
            zout.writestr(item, data)
    L.say(f"    note: {n} links to parts missing from the file were dropped before reading")
    buf.seek(0)
    return openpyxl.load_workbook(buf)


def last_row(ws):
    """Last row holding a value (openxlsx sheets can report a larger dimension)."""
    for r in range(ws.max_row, 0, -1):
        if any(c.value is not None for c in ws[r]):
            return r
    return 0


def header(ws):
    return [c.value for c in ws[1]]


def need(L, ok, what):
    if not L.chk(ok, what):
        raise TableFail(what)


def edit_header(L, wb, sheet, old, new):
    hits = [c for c in wb[sheet][1] if c.value == old]
    need(L, len(hits) == 1, f"{sheet}: column '{old}' found once, renamed '{new}'")
    hits[0].value = new


def edit_replace(L, wb, sheets, old, new, exact, n_expected):
    n = 0
    for s in sheets:
        for row in wb[s].iter_rows():
            for c in row:
                v = c.value
                if not isinstance(v, str) or v.startswith("="):
                    continue
                if exact and v == old:
                    c.value = new
                    n += 1
                elif not exact and old in v:
                    c.value = v.replace(old, new)
                    n += 1
    ok = n_expected is None or n == n_expected
    need(L, ok, f"{', '.join(sheets)}: '{old}' -> '{new}' in {n} cells"
                + ("" if n_expected is None else f" (expected {n_expected})"))


def edit_formula(L, wb, sheets, old, new, n_min):
    for s in sheets:
        n = 0
        for row in wb[s].iter_rows():
            for c in row:
                v = c.value
                if isinstance(v, str) and v.startswith("=") and old in v:
                    c.value = v.replace(old, new)
                    n += 1
        need(L, n >= n_min, f"{s}: '{old}' -> '{new}' in {n} formulas (at least {n_min})")


def edit_map(L, wb, sheets, column, mapping, keep):
    for s in sheets:
        ws = wb[s]
        h = header(ws)
        need(L, h.count(column) == 1, f"{s}: column '{column}' found once")
        j = h.index(column) + 1
        n, bad = 0, set()
        for r in range(2, last_row(ws) + 1):
            c = ws.cell(r, j)
            v = c.value
            if v is None:
                continue
            if v in mapping:
                c.value = mapping[v]
                n += 1
            elif not str(v).startswith(keep):
                bad.add(str(v))
        need(L, not bad, f"{s}: {n} '{column}' labels relabelled"
                         + (f"; labels with no published name: {sorted(bad)}" if bad else ""))


def edit_strip_root(L, wb, sheet, column, root):
    ws = wb[sheet]
    h = header(ws)
    need(L, h.count(column) == 1, f"{sheet}: column '{column}' found once")
    j = h.index(column) + 1
    # The root as this run sees it, then any absolute folder in front of a path
    # that starts with Manuscript/ or a Part folder (other machines and the container).
    rootp = str(root).replace("\\", "/").rstrip("/") + "/"
    abs_prefix = re.compile(r"(?<=\()(?:[A-Za-z]:)?/[^()\s]*?/(?=(?:Manuscript|Part\d+\w*)/)")
    n = 0
    for r in range(2, last_row(ws) + 1):
        c = ws.cell(r, j)
        if isinstance(c.value, str):
            v = abs_prefix.sub("", c.value.replace(rootp, ""))
            if v != c.value:
                c.value = v
                n += 1
    L.chk(True, f"{sheet}: root folder removed from {n} paths")


def restyle(ws):
    for row in ws.iter_rows():
        for c in row:
            c.font = F10
    for c in ws[1]:
        c.font = F10B
    ws.freeze_panes = "A2"


def write_contents(wb, spec_text, width_a, width_d=None):
    ws = wb.create_sheet("Contents", 0)
    ws["A1"] = spec_text["title"]
    ws["A1"].font = F12B
    ws["A2"] = spec_text["subtitle"]
    ws["A2"].font = F10
    for j, h in enumerate(["Sheet", "Rows", "Description"], 1):
        ws.cell(4, j, h).font = F10B
    for i, d in enumerate(spec_text["sheets"], 5):
        n = last_row(wb[d["sheet"]]) - 1
        for j, v in enumerate([d["sheet"], n, d["description"]], 1):
            c = ws.cell(i, j, v)
            c.font = F10
            c.alignment = Alignment(vertical="top", wrap_text=(j == 3))
        ws.cell(i, 2).number_format = "#,##0"
    ws.column_dimensions["A"].width = width_a
    ws.column_dimensions["B"].width = 9
    ws.column_dimensions["C"].width = 110
    if width_d:
        ws.column_dimensions["D"].width = width_d
    ws.freeze_panes = "A5"


##### One table #####
def finalize(L, t, spec, text, root, out_dir):
    src = root / spec["src"]
    L.say(f"S{t}  {spec['src']}")
    need(L, src.is_file(), f"source workbook found")
    wb = load_book(L, src)
    have = wb.sheetnames
    want = [s for s, _ in spec["sheets"]]
    missing = [s for s in want if s not in have]
    need(L, not missing, f"{len(want)} source sheets found" + (f"; missing {missing}" if missing else ""))
    dropped = [s for s in have if s not in want]
    for s in dropped:
        del wb[s]
    if dropped:
        L.say(f"    sheets dropped: {', '.join(dropped)}")
    for s, new in spec["sheets"]:
        wb[s].title = new
    wb._sheets = [wb[new] for _, new in spec["sheets"]]
    for e in spec["edits"]:
        kind = e[0]
        if kind == "header":
            edit_header(L, wb, *e[1:])
        elif kind == "replace":
            edit_replace(L, wb, *e[1:])
        elif kind == "formula":
            edit_formula(L, wb, *e[1:])
        elif kind == "map":
            edit_map(L, wb, *e[1:])
        elif kind == "strip_root":
            edit_strip_root(L, wb, *e[1:], root)
        else:
            raise ValueError(kind)
    if spec.get("restyle"):
        for ws in wb.worksheets:
            restyle(ws)
    published = [new for _, new in spec["sheets"]]
    need(L, [d["sheet"] for d in text["sheets"]] == published,
         "sheet list matches the Contents text in supp_tables_contents.json")
    write_contents(wb, text, spec.get("width_a", 30), spec.get("width_d"))
    wb.active = 0
    for ws in wb.worksheets:
        ws.sheet_view.tabSelected = ws.title == "Contents"
    out = out_dir / text["file"]
    wb.save(out)
    L.say(f"    written: {out.relative_to(root)} ({len(published)} sheets + Contents)")


def parse_tables(arg):
    if not arg:
        return sorted(TABLES)
    out = []
    for part in arg.split(","):
        a, _, b = part.partition("-")
        out += list(range(int(a), int(b or a) + 1))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--root", default=os.environ.get("SEAAD_ROOT", os.getcwd()))
    ap.add_argument("--tables", default="", help="for example 12 or 13-17,19 (default: all)")
    args = ap.parse_args()
    root = pathlib.Path(args.root).resolve()
    out_dir = root / OUT_SUB
    out_dir.mkdir(parents=True, exist_ok=True)
    texts = json.loads(CONTENTS_JSON.read_text(encoding="utf-8"))
    L = Log()
    L.say(f"finalize_supp_tables.py | root {root}")
    tables = parse_tables(args.tables)
    unknown = [t for t in tables if t not in TABLES]
    if unknown:
        L.say(f"No spec yet for: {', '.join(f'S{t}' for t in unknown)}")
    written, failed = [], []
    for t in tables:
        if t not in TABLES:
            continue
        try:
            finalize(L, t, TABLES[t], texts[f"S{t}"], root, out_dir)
            written.append(t)
        except TableFail:
            failed.append(t)
            L.say(f"    S{t} not written")
    L.say("")
    L.say(f"Written: {', '.join(f'S{t}' for t in written) or 'none'}")
    if failed:
        L.say(f"Not written: {', '.join(f'S{t}' for t in failed)}")
    verdict = "PASS" if not failed and not unknown else "FAIL"
    L.say(f"VERDICT: {verdict}")
    (out_dir / "finalize_supp_tables_log.txt").write_text("\n".join(L.lines) + "\n", encoding="utf-8")
    sys.exit(0 if verdict == "PASS" else 1)


if __name__ == "__main__":
    main()
