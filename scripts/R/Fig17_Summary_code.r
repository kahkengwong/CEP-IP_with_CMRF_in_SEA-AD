# =============================================================================
# Fig17_Summary_code.r
# SEA-AD CEP-IP/CMRF, summary of all tests of Module post-IP% (Figure 17)
#   Reads the closed supplementary workbooks (S12, SuppTable_Fig15 to Fig19 and
#   the Part 13 and 13b tables), classifies every model behind each test by the
#   fixed rules below and draws the summary figure. No analysis is rerun.
#   FIG_MODE "preprint" (default) labels it Figure 17; "journal" labels it Figure 12.
#
# Usage (from the SEA-AD folder, or set SEAAD_ROOT), after Parts 10h, 12 and 13
#   source("Fig17_Summary_code.r")
#
# Outputs (Manuscript/Fig17/Preprint/ or Manuscript/Fig17/Journal/): figure
#   (PDF, PNG, SVG), Fig17_Summary.xlsx, cells, models, lead-result and plotted
#   CSVs, log.
#
# Terms: a TEST is one line of the grid (21 tests); a MODEL is one estimate
# behind a test, for one arm, region set, outcome and variant.
#
# Rules (fixed before the grid is filled):
#   G1  Significant = FDR < 0.05 in the BH family of the source, as built.
#   G2  A model holds when significant with the same sign as the cov4
#       (test 1) coefficient, or when the source verdict named in TESTDEF
#       says so.
#   G3  Tests 2-13 are read for an outcome only if test 1 holds for that
#       outcome in that arm and region set; otherwise the outcome is not
#       applicable.
#   G4  Unreadable = flagged by the source or VIF >= 5 in a two-predictor
#       model. Untestable = no estimate, or inconclusive by the source.
#       Both are excluded from n.
#   G5  k = models that hold, n = readable models. Holds k = n > 0; Partial
#       0 < k < n; Lost k = 0 < n; Unreadable n = 0 with an unreadable model;
#       Untestable n = 0 otherwise; Not applicable = no applicable outcome.
#   G6  Tally per column (tests 1-13): tests held over tests read (n > 0).
# =============================================================================

suppressPackageStartupMessages({
    library(ggplot2)
    library(readxl)
    library(grid)
})

# == Configuration ===========================================================
if (!exists("FIG17_ROOT")) FIG17_ROOT <- Sys.getenv("SEAAD_ROOT", getwd())
ROOT     <- FIG17_ROOT
# Figure numbering: "preprint" (Batteries are Figures 12-16, this is Figure 17) or
# "journal" (Batteries are Supplementary Figures S20-S24, this is Figure 12).
if (!exists("FIG_MODE")) FIG_MODE <- "preprint"
if (!FIG_MODE %in% c("preprint", "journal")) stop("FIG_MODE must be \"preprint\" or \"journal\"")
PRE      <- FIG_MODE == "preprint"
FIG_NO   <- if (PRE) "Figure 17" else "Figure 12"
# Battery panel labels by mode: Fig 12-16 (preprint) or Supplementary Figure S20-S24 (journal)
BP <- function(x) {
    if (PRE) return(x)
    n <- as.integer(sub("^(1[2-6]).*$", "\\1", x))
    paste0("S", n + 8, sub("^1[2-6]", "", x))
}
OUT_DIR  <- file.path(ROOT, "Manuscript", "Fig17", if (PRE) "Preprint" else "Journal")
SCRIPT   <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e)
    file.path(ROOT, "Manuscript", "Fig17", "Fig17_Summary_code.r"))
FONT     <- "Arial"
ALPHA    <- 0.05
AUC_BAND <- 0.60            # AUC below 0.60 = not recovering
VIF_MAX  <- 5
OUTCOMES <- c("MMSE score", "CASI score")          # grid and first two figures
PATH_OUT <- c("Thal phase", "Braak stage")         # amyloid, tau
OUTCOMES_ALL <- c(OUTCOMES, PATH_OUT)
# Tests whose MMSE and CASI cells are pooled into one lead cell (none).
MERGE_COG <- character(0)
# Tests left out of every tally so all arms are counted on the same tests (check C22).
TALLY_EXCLUDE <- c("T7", "T11")
TALLY_ARMS    <- c("M4", "PVM_3", "PVM_1")
in_tally <- function(t) TESTDEF$block[match(t, TESTDEF$test)] == 1 & !t %in% TALLY_EXCLUDE
ARMS     <- c("M4", "PVM_3", "PVM_1", "HKG8")
RSETS    <- c("NeoC", "NonNeoC")
LEADS    <- data.frame(arm = c("M4", "PVM_3", "PVM_1"), region_set = "NeoC",
                       stringsAsFactors = FALSE)   # lead-result column pairs
HL_ARM   <- "M4"            # the lead pair shown in the portrait figure
HL_RS    <- "NeoC"
ARM_COL  <- c(M4 = "#3E3480", PVM_3 = "#1A9C78", PVM_1 = "#F2A343", HKG8 = "#9E9E9E")

SRC <- list(
    S12 = list(path = file.path(ROOT, "Manuscript", "Supp_Tables",
                                "Supplementary_Table_S12_Prognosis_&_Ablation.xlsx")),
    F15 = list(path = file.path(ROOT, "Part12_rerun", "07_Fig15", "SuppTable_Fig15.xlsx")),
    F16 = list(path = file.path(ROOT, "Part12_rerun", "08_Fig16", "SuppTable_Fig16.xlsx")),
    F17 = list(path = file.path(ROOT, "Part12_rerun", "09_Fig17", "SuppTable_Fig17.xlsx")),
    F18 = list(path = file.path(ROOT, "Part12_rerun", "10_Fig18", "SuppTable_Fig18.xlsx")),
    F19 = list(path = file.path(ROOT, "Part12_rerun", "11_Fig19", "SuppTable_Fig19.xlsx")),
    # Part 13 (M4 share) and Part 13b (random PVM_3 subsets).
    P13 = list(path = file.path(ROOT, "Part13_M4Share", "Part13_CP10K_M4Share_SuppTable.xlsx")),
    P13B = list(path = file.path(ROOT, "Part13b_M4RandomSubset",
                                 "Part13b_CP10K_M4RandomSubset_SuppTable.xlsx")))
SRC_LABEL <- c(S12 = "SuppTable_S12", F15 = "SuppTable_Fig15", F16 = "SuppTable_Fig16",
               F17 = "SuppTable_Fig17", F18 = "SuppTable_Fig18", F19 = "SuppTable_Fig19",
               P13 = "Part13 SuppTable", P13B = "Part13b SuppTable")

# Test definitions in display order. 'fig' uses manuscript numbers.
TESTDEF <- data.frame(
    # T20-T21 are the M4-only tests, shown as 14-15 (block 3); T14-T19 are shown as 16-21.
    test = paste0("T", 1:21),
    number = c(1:13, 16:21, 14:15),
    block = c(rep(1L, 13), rep(2L, 6), rep(3L, 2)),
    group = c(rep("Clinical association and adjustment", 4),
              rep("Against alternative explanations", 4),
              rep("Batteries 1, 3, 4 and 5: the clinical association retested", 5),
              rep("Batteries 2, 4 and 5: properties of Module post-IP% (per arm, NeoC + NonNeoC)", 6),
              rep("M4 against PVM_3 (M4 only)", 2)),
    battery = c(rep("", 8), "B1", "B1", "B3", "B4", "B5", "B2", "B2", "B2", "B2", "B4", "B5", "", ""),
    label = c("Associated at cov4",
              "Holds with pathology in the model",
              "Holds with tissue quality (cov7)",
              "Holds with each of 18 added covariates",
              "Module exceeds non-Module",
              "Holds against cluster proportion",
              "Holds beyond cell count x gene score",
              "Holds for all four boundaries",
              "Holds with contamination in the model",
              "Holds after doublet removal (3 cuts)",
              "Module carries it with %MT in the model",
              "Holds against 15 nuclear panels + reference",
              "Holds without the 11 SA donors",
              "Donor-structured (raw and TQR)",
              "Repeatable between regions",
              "Offset share within the 10 pp margin",
              "No association with library intronic rate",
              "Not recovered by per-cell %MT",
              "Donor structure holds without SA donors",
              "Holds with the M4 share in the model",
              "M4 cells exceed random PVM_3 cells"),
    # Battery panels follow FIG_MODE (12A in the preprint, S20A in the journal).
    fig = c("9A", "9B", "S14A", "S14E", "9D", "10B", "10C", "11A",
            BP(c("12A", "12B", "14D", "15A", "16C", "13A", "13B", "13E", "13F", "15B-C", "16A")),
            if (PRE) c("S20A", "S21A") else c("Table S18", "Table S18")),
    source = c("S12 B_Fig9A_variance", "S12 C_Fig9B_cognition", "S12 P_SuppFig14_cov4_cov7",
               "S12 Q_SuppFig14_leave_one_in", "S12 E_Fig9D_mod_vs_nonmod",
               "S12 G_Fig10B_competitors", "S12 H_Fig10CD_decomposition",
               "S12 I_Fig11A_boundary_assoc", "Fig15 Fig15A_readout", "Fig15 Fig15B_readout",
               "Fig17 Fig17D_verdicts", "Fig18 Fig18A_verdicts", "Fig19 Fig19C_clinical",
               "Fig16 Fig16A_readout", "Fig16 Fig16B_stats", "Fig16 Fig16E_readout",
               "Fig16 Fig16F_stats", "Fig18 Fig18B_readout; Fig18C_readout",
               "Fig19 Fig19A_readout", "Part13 C_Models (S03, S07)", "Part13b D_Summary (bA)"),
    variants = c("1", "+ADNC, +BTC, +ALL", "cov7", "18 added covariates", "1",
                 "Cluster proportion (supertype)", "score x block (cov4, +BTC)",
                 "4 boundaries x 2 blocks", "1", "doublet score 0.15, 0.20, 0.30", "1",
                 "15 panels + matched reference", "1",
                 "Module post-IP%, TQR", "module_post_pct, module_resid_tq4", "1",
                 "between donors, within donor", "pooled AUC and within-donor AUC per region set",
                 "16A and 16B Module quantities", "cov4 + M4 share; +BTC + M4 share (cognition)",
                 "M4 cells (PVM_3 labels) against 10,000 random PVM_3 sets per region set"),
    rule = c("FDR < 0.05 (no entry condition)",
             "FDR < 0.05, same sign as test 1",
             "FDR < 0.05, same sign as test 1",
             "p_added < 0.05, same sign as test 1 (as Supp Fig 14E)",
             "verdict = 'module > non-module'",
             "verdict_shown in {module only, both}; 'not readable' = unreadable",
             "verdict = 'HOLDS'; 'unreadable' = unreadable",
             "FDR < 0.05, same sign as test 1",
             "fdr_adjusted < 0.05, beta_adjusted same sign as test 1",
             "fdr < 0.05, same sign as test 1",
             "VIF < 5 and fdr_module < 0.05, same sign as test 1 (Module only or Both); VIF >= 5 = unreadable",
             "fdr_second < 0.05, same sign as test 1 (Module only or Both)",
             "class = Holds; Lost (either) = fails; Neither = untestable",
             "fdr < 0.05 against the shuffle null",
             "fdr < 0.05 and split-half rs > 0",
             "verdict = 'within the margin'; underpowered or untestable = untestable",
             "Module rs fdr >= 0.05 (no association) holds; missing = untestable",
             "pooled AUC < 0.60; within-donor median AUC fdr >= 0.05; missing = untestable",
             "read where fdr_all < 0.05; fdr_wo < 0.05 with the same sign",
             "FDR < 0.05, same sign as test 1; largest VIF >= 5 = unreadable",
             "FDR < 0.05 (BH within region set over three cognition outcomes) and the M4 coefficient beyond the random-set median in the direction of test 1"),
    stringsAsFactors = FALSE)
TALLY_SUFFIX <- sprintf("1-13 without %s", paste(TESTDEF$number[match(TALLY_EXCLUDE, TESTDEF$test)],
                                                 collapse = ", "))

# == Logging =================================================================
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
LOG <- file.path(OUT_DIR, "Fig17_Summary_log.txt")
sink(LOG, split = TRUE)
on.exit_sink <- function() while (sink.number() > 0) sink()

CHECKS <- data.frame(check = character(), verdict = character(), detail = character(),
                     stringsAsFactors = FALSE)
add_check <- function(id, ok, detail, warn_only = FALSE) {
    v <- if (isTRUE(ok)) "PASS" else if (warn_only) "WARN" else "FAIL"
    CHECKS[nrow(CHECKS) + 1, ] <<- c(id, v, detail)
    cat(sprintf("  [%s] %s: %s\n", v, id, detail))
}

tryCatch({

cat(sprintf("%s summary (Fig17_Summary_code.r); FIG_MODE = %s\n", FIG_NO, FIG_MODE))
cat(sprintf("Run: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat(sprintf("Script: %s\n", SCRIPT))
cat(sprintf("R %s; ggplot2 %s; readxl %s\n", getRversion(),
            packageVersion("ggplot2"), packageVersion("readxl")))
cat(sprintf("Root: %s\n\n", ROOT))

# == 1. Sources ==============================================================
cat("== 1. Sources\n")
for (k in names(SRC)) {
    p <- SRC[[k]]$path
    if (!file.exists(p)) stop("missing source: ", p)
    add_check(paste0("C1 source ", SRC_LABEL[[k]]), TRUE,
              sprintf("present, modified %s (%s)", format(file.mtime(p), "%Y-%m-%d %H:%M"), p))
}

rd <- function(key, sheet) {
    x <- as.data.frame(read_excel(SRC[[key]]$path, sheet = sheet, guess_max = 1e6))
    cat(sprintf("  read %-16s %-28s %6d rows\n", SRC_LABEL[[key]], sheet, nrow(x)))
    x
}
B   <- rd("S12", "B_Fig9A_variance")
C   <- rd("S12", "C_Fig9B_cognition")
D   <- rd("S12", "D_Fig9C_pathology")
E   <- rd("S12", "E_Fig9D_mod_vs_nonmod")
Fb  <- rd("S12", "F_Fig10A_benchmark")
G   <- rd("S12", "G_Fig10B_competitors")
H   <- rd("S12", "H_Fig10CD_decomposition")
I   <- rd("S12", "I_Fig11A_boundary_assoc")
M   <- rd("S12", "M_Fig11C_tests")
P   <- rd("S12", "P_SuppFig14_cov4_cov7")
Q   <- rd("S12", "Q_SuppFig14_leave_one_in")
F15A <- rd("F15", "Fig15A_readout")
F15B <- rd("F15", "Fig15B_readout")
F16A <- rd("F16", "Fig16A_readout")
F16B <- rd("F16", "Fig16B_stats")
F16C <- rd("F16", "Fig16C_readout")
F16D <- rd("F16", "Fig16D_stats")
F16E <- rd("F16", "Fig16E_readout")
F16F <- rd("F16", "Fig16F_stats")
F17D <- rd("F17", "Fig17D_verdicts")
F18A <- rd("F18", "Fig18A_verdicts")
F18P <- rd("F18", "Fig18A_panels")
F18B <- rd("F18", "Fig18B_readout")
F18C <- rd("F18", "Fig18C_readout")
F19A <- rd("F19", "Fig19A_readout")
F19C <- rd("F19", "Fig19C_clinical")
V17D <- rd("F19", "Verdicts_17D")
V18A <- rd("F19", "Verdicts_18A")
P13M <- rd("P13", "C_Models")
P13S <- rd("P13B", "D_Summary")

# == 2. Models ===============================================================
cat("\n== 2. Models\n")
sig <- function(f) !is.na(f) & f < ALPHA
same <- function(b, ref) !is.na(b) & !is.na(ref) & sign(b) == sign(ref)
# All three inputs are expanded to the same length before ifelse().
st <- function(holds, unreadable = FALSE, untestable = FALSE) {
    n <- length(holds)
    holds <- rep_len(holds, n); unreadable <- rep_len(unreadable, n)
    untestable <- rep_len(untestable, n)
    out <- ifelse(untestable, "untestable", ifelse(unreadable, "unreadable",
                  ifelse(holds, "holds", "fails")))
    stopifnot(length(out) == n)
    out
}
mk <- function(test, d, region_set, outcome, variant, estimate, fdr, status, src_sheet) {
    if (nrow(d) == 0) return(NULL)
    data.frame(test = test, arm = d$arm, region_set = region_set, outcome = outcome,
               variant = as.character(variant), estimate = as.numeric(estimate),
               fdr = as.numeric(fdr), status = status, source = src_sheet,
               stringsAsFactors = FALSE)
}

# T1 and the reference sign (entry condition G4)
b1 <- B[B$population == "Module post-IP" & B$outcome %in% OUTCOMES_ALL, ]
REF <- data.frame(arm = b1$arm, region_set = b1$region_set, outcome = b1$outcome,
                  beta_ref = b1$beta_std, entry = sig(b1$FDR), stringsAsFactors = FALSE)
with_ref <- function(d, rs_col = "region_set") {
    d$region_set_ <- d[[rs_col]]
    merge(d, setNames(REF, c("arm", "region_set_", "outcome", "beta_ref", "entry")),
          by = c("arm", "region_set_", "outcome"), all.x = TRUE)
}
T <- list()
T$T1 <- mk("T1", b1, b1$region_set, b1$outcome, "cov4", b1$beta_std, b1$FDR,
           st(sig(b1$FDR), untestable = is.na(b1$FDR)), "S12 B_Fig9A_variance")

d <- with_ref(C[C$block %in% c("+ADNC", "+BTC", "+ALL") & C$outcome %in% OUTCOMES_ALL, ])
T$T2 <- mk("T2", d, d$region_set_, d$outcome, d$block, d$beta_std, d$FDR,
           st(sig(d$FDR) & same(d$beta_std, d$beta_ref), untestable = is.na(d$FDR)),
           "S12 C_Fig9B_cognition")
d <- with_ref(D[D$block == "+MMSE + CASI" & D$outcome %in% PATH_OUT, ])
T$T2 <- rbind(T$T2, mk("T2", d, d$region_set_, d$outcome, d$block, d$beta_std, d$FDR,
           st(sig(d$FDR) & same(d$beta_std, d$beta_ref), untestable = is.na(d$FDR)),
           "S12 D_Fig9C_pathology"))

d <- with_ref(P[P$block == "cov7" & P$outcome %in% OUTCOMES_ALL, ])
T$T3 <- mk("T3", d, d$region_set_, d$outcome, "cov7", d$beta_std, d$FDR,
           st(sig(d$FDR) & same(d$beta_std, d$beta_ref), untestable = is.na(d$FDR)),
           "S12 P_SuppFig14_cov4_cov7")

d <- with_ref(Q[Q$outcome %in% OUTCOMES_ALL, ])
T$T4 <- mk("T4", d, d$region_set_, d$outcome, d$added_covariate, d$beta_added, d$p_added,
           st(sig(d$p_added) & same(d$beta_added, d$beta_ref), untestable = is.na(d$p_added)),
           "S12 Q_SuppFig14_leave_one_in (p, not FDR)")

d <- with_ref(E[E$outcome %in% OUTCOMES_ALL, ])
T$T5 <- mk("T5", d, d$region_set_, d$outcome, "module vs non-module", d$difference, d$FDR,
           st(d$verdict == "module > non-module", untestable = is.na(d$verdict)),
           "S12 E_Fig9D_mod_vs_nonmod")

d <- with_ref(G[G$competitor == "Cluster proportion (supertype)" & G$outcome %in% OUTCOMES_ALL, ])
T$T6 <- mk("T6", d, d$region_set_, d$outcome, d$competitor_abbrev, d$beta_module, d$FDR_module,
           st(d$verdict_shown %in% c("module only", "both"),
              unreadable = d$verdict_shown == "not readable", untestable = is.na(d$verdict_shown)),
           "S12 G_Fig10B_competitors")

d <- with_ref(H[H$outcome %in% OUTCOMES_ALL, ])
T$T7 <- mk("T7", d, d$region_set_, d$outcome, paste(d$score, d$block), d$beta_module,
           d$FDR_module, st(d$verdict == "HOLDS", unreadable = d$verdict == "unreadable",
                            untestable = is.na(d$verdict)), "S12 H_Fig10CD_decomposition")

d <- with_ref(I[I$outcome %in% OUTCOMES_ALL, ])
T$T8 <- mk("T8", d, d$region_set_, d$outcome, paste(d$boundary, d$block), d$beta_std, d$FDR,
           st(sig(d$FDR) & same(d$beta_std, d$beta_ref), untestable = is.na(d$FDR)),
           "S12 I_Fig11A_boundary_assoc")
d <- with_ref(M[M$model == "module_only" & M$predictor == "Module post-IP %" &
                M$block == "cov4" & M$outcome %in% PATH_OUT, ])
T$T8 <- rbind(T$T8, mk("T8", d, d$region_set_, d$outcome, paste(d$boundary, d$block),
           d$beta_std, d$FDR, st(sig(d$FDR) & same(d$beta_std, d$beta_ref),
                                 untestable = is.na(d$FDR)), "S12 M_Fig11C_tests"))

d <- with_ref(F15A[F15A$outcome %in% OUTCOMES_ALL, ])
T$T9 <- mk("T9", d, d$region_set_, d$outcome, "contamination-adjusted", d$beta_adjusted,
           d$fdr_adjusted, st(sig(d$fdr_adjusted) & same(d$beta_adjusted, d$beta_ref),
                              untestable = is.na(d$fdr_adjusted)), "Fig15 Fig15A_readout")

d <- with_ref(F15B[F15B$outcome %in% OUTCOMES_ALL, ])
T$T10 <- mk("T10", d, d$region_set_, d$outcome, sprintf("doublet %.2f", d$threshold), d$beta,
            d$fdr, st(sig(d$fdr) & same(d$beta, d$beta_ref), untestable = is.na(d$fdr)),
            "Fig15 Fig15B_readout")

d <- with_ref(F17D[F17D$outcome %in% OUTCOMES_ALL, ])
T$T11 <- mk("T11", d, d$region_set_, d$outcome, "Module + %MT", d$beta_module, d$fdr_module,
            st(sig(d$fdr_module) & same(d$beta_module, d$beta_ref),
               unreadable = !is.na(d$vif) & d$vif >= VIF_MAX, untestable = is.na(d$fdr_module)),
            "Fig17 Fig17D_verdicts")

d <- with_ref(F18A[F18A$predictor == "module" & F18A$outcome %in% OUTCOMES_ALL, ])
T$T12 <- mk("T12", d, d$region_set_, d$outcome, d$panel, d$beta_second, d$fdr_second,
            st(sig(d$fdr_second) & same(d$beta_second, d$beta_ref),
               unreadable = !is.na(d$vif) & d$vif >= VIF_MAX, untestable = is.na(d$fdr_second)),
            "Fig18 Fig18A_verdicts")

d <- with_ref(F19C[F19C$outcome %in% OUTCOMES_ALL, ])
T$T13 <- mk("T13", d, d$region_set_, d$outcome, "without SA donors", d$est_wo, d$fdr_wo,
            st(d$class == "Holds", untestable = is.na(d$class) | d$class == "Neither"),
            "Fig19 Fig19C_clinical")

# Tests 16-21, IDs T14-T19 (per arm; region_set "both" unless the source splits it)
d <- F16A[F16A$quantity %in% c("Module post-IP%", "Module post-IP%, TQR"), ]
T$T14 <- mk("T14", d, "both", "", d$quantity, d$donor_share_pct, d$fdr,
           st(sig(d$fdr), untestable = is.na(d$fdr)), "Fig16 Fig16A_readout")

d <- F16B[F16B$quantity %in% c("module_post_pct", "module_resid_tq4"), ]
T$T15 <- mk("T15", d, "both", "", d$quantity, d$estimate, d$fdr,
           st(sig(d$fdr) & d$estimate > 0, untestable = is.na(d$fdr)), "Fig16 Fig16B_stats")

d <- F16E
T$T16 <- mk("T16", d, "both", "", "offset share", d$difference_pp, d$fdr,
           st(d$verdict == "within the margin",
              untestable = is.na(d$verdict) | d$verdict != "within the margin"),
           "Fig16 Fig16E_readout")

d <- F16F[F16F$quantity == "module" & grepl("\\|rs$", F16F$family), ]
T$T17 <- mk("T17", d, "both", "", d$level, d$estimate, d$fdr,
           st(!is.na(d$fdr) & d$fdr >= ALPHA, untestable = is.na(d$fdr)), "Fig16 Fig16F_stats")

d1 <- F18B
d2 <- F18C
T$T18 <- rbind(
    mk("T18", d1, d1$region_set, "", paste("pooled AUC", d1$region_set), d1$auc, NA,
       st(!is.na(d1$auc) & d1$auc < AUC_BAND, untestable = is.na(d1$auc)),
       "Fig18 Fig18B_readout"),
    mk("T18", d2, d2$region_set, "", paste("within-donor AUC", d2$region_set), d2$median_auc,
       d2$fdr, st(!is.na(d2$fdr) & d2$fdr >= ALPHA, untestable = is.na(d2$median_auc)),
       "Fig18 Fig18C_readout"))
T$T18$region_set <- "both"

d <- F19A[F19A$panel %in% c("16A", "16B") &
          F19A$quantity %in% c("module_post_pct", "module_resid_tq4") & sig(F19A$fdr_all), ]
T$T19 <- mk("T19", d, "both", "", paste(d$panel, d$quantity), d$est_wo, d$fdr_wo,
           st(sig(d$fdr_wo) & same(d$est_wo, d$est_all), untestable = is.na(d$fdr_wo)),
           "Fig19 Fig19A_readout")

# Tests 14-15: M4 only, per region set and outcome, entry condition G3
for (k in c("n", "beta", "fdr", "vif_full")) P13M[[k]] <- suppressWarnings(as.numeric(P13M[[k]]))
# Largest VIF of each fitted model; a VIF the workbook left empty (Inf) counts as unreadable
vm13 <- aggregate(list(vif_max = P13M$vif_full),
                  by = list(spec = P13M$spec, dataset = P13M$dataset, outcome = P13M$outcome),
                  FUN = function(v) if (anyNA(v)) Inf else max(v))
d <- P13M[P13M$spec %in% c("S03", "S07") & P13M$predictor == "module" &
          P13M$outcome %in% OUTCOMES_ALL, ]
d$vif_max <- vm13$vif_max[match(paste(d$spec, d$dataset, d$outcome),
                                paste(vm13$spec, vm13$dataset, vm13$outcome))]
d$arm <- "M4"; d$region_set <- d$dataset
d <- with_ref(d)
T$T20 <- mk("T20", d, d$region_set_, d$outcome,
            ifelse(d$spec == "S03", "cov4 + M4 share", "+BTC + M4 share"), d$beta, d$fdr,
            st(sig(d$fdr) & same(d$beta, d$beta_ref),
               unreadable = is.na(d$vif_max) | d$vif_max >= VIF_MAX, untestable = is.na(d$fdr)),
            "Part13 C_Models (S03, S07)")
for (k in c("observed_M4_PVM3labels", "M4_own_labels", "null_median", "fdr"))
    P13S[[k]] <- suppressWarnings(as.numeric(P13S[[k]]))
d <- P13S[P13S$readout %in% paste0("bA_", c("MMSE", "CASI")) & P13S$outcome %in% OUTCOMES_ALL, ]
d$arm <- "M4"; d$region_set <- d$dataset
d <- with_ref(d)
d$diff <- d$observed_M4_PVM3labels - d$null_median
T$T21 <- mk("T21", d, d$region_set_, d$outcome, "M4 cells vs 10,000 random PVM_3 sets", d$diff,
            d$fdr, st(sig(d$fdr) & same(d$diff, d$beta_ref), untestable = is.na(d$fdr)),
            "Part13b D_Summary (bA)")

MODELS <- do.call(rbind, T)
rownames(MODELS) <- NULL
MODELS <- merge(MODELS, REF[, c("arm", "region_set", "outcome", "entry")],
                by = c("arm", "region_set", "outcome"), all.x = TRUE)
blk_of <- function(t) TESTDEF$block[match(t, TESTDEF$test)]
MODELS$applicable <- ifelse(MODELS$test == "T1" | blk_of(MODELS$test) == 2,
                            TRUE, !is.na(MODELS$entry) & MODELS$entry)
MODELS <- MODELS[order(match(MODELS$test, TESTDEF$test), match(MODELS$arm, ARMS),
                       MODELS$region_set, MODELS$outcome, MODELS$variant), ]
for (r in TESTDEF$test) {
    x <- MODELS[MODELS$test == r & MODELS$applicable, ]
    cat(sprintf("  %-4s %4d applicable models: holds %4d, fails %4d, unreadable %3d, untestable %3d\n",
                r, nrow(x), sum(x$status == "holds"), sum(x$status == "fails"),
                sum(x$status == "unreadable"), sum(x$status == "untestable")))
}

# == 3. Cells ================================================================
cat("\n== 3. Cells\n")
app_out <- function(a, rs) intersect(OUTCOMES,
    REF$outcome[REF$arm == a & REF$region_set == rs & REF$entry])
cell_one <- function(r, a, rs) {
    blk <- blk_of(r)
    if (blk == 3 && a != "M4")
        return(data.frame(test = r, block = blk, arm = a, region_set = rs, outcomes_applicable = "",
                          k = 0L, n = 0L, unreadable = 0L, untestable = 0L, missing_outcomes = "",
                          verdict = "M4 only", stringsAsFactors = FALSE))
    if (blk %in% c(1, 3)) {
        app <- if (r == "T1") OUTCOMES else app_out(a, rs)
        x <- MODELS[MODELS$test == r & MODELS$arm == a & MODELS$region_set == rs &
                    MODELS$outcome %in% app, ]
        n_app <- length(app)
        missing_out <- setdiff(app, x$outcome)
    } else {
        x <- MODELS[MODELS$test == r & MODELS$arm == a, ]
        app <- character(0)
        n_app <- nrow(x)
        missing_out <- character(0)
    }
    k <- sum(x$status == "holds"); n <- sum(x$status %in% c("holds", "fails"))
    nu <- sum(x$status == "unreadable"); nt <- sum(x$status == "untestable")
    v <- if (n_app == 0) "not applicable" else if (n > 0) {
        if (k == n) "holds" else if (k == 0) "lost" else "partial"
    } else if (nu > 0) "unreadable" else "untestable"
    data.frame(test = r, block = blk, arm = a, region_set = rs,
               outcomes_applicable = paste(app, collapse = "; "), k = k, n = n,
               unreadable = nu, untestable = nt,
               missing_outcomes = paste(missing_out, collapse = "; "),
               verdict = v, stringsAsFactors = FALSE)
}
CELLS <- do.call(rbind, c(
    lapply(TESTDEF$test[TESTDEF$block %in% c(1, 3)], function(r)
        do.call(rbind, lapply(ARMS, function(a)
            do.call(rbind, lapply(RSETS, function(rs) cell_one(r, a, rs)))))),
    lapply(TESTDEF$test[TESTDEF$block == 2], function(r)
        do.call(rbind, lapply(ARMS, function(a) cell_one(r, a, "both"))))))
CELLS$kn <- ifelse(CELLS$n > 0 & !CELLS$verdict %in% c("not applicable", "M4 only"),
                   paste0(CELLS$k, "/", CELLS$n), "")
CELLS <- merge(CELLS, TESTDEF[, c("test", "number", "label", "battery", "fig", "source", "rule")],
               by = "test")
CELLS <- CELLS[order(CELLS$number, match(CELLS$arm, ARMS), CELLS$region_set), ]
print(table(test = factor(CELLS$test, TESTDEF$test), verdict = CELLS$verdict))

cat("\n  Grid (k/n; verdict initial h/p/l, u = unreadable or untestable, - = not applicable):\n")
code <- c(holds = "h", partial = "p", lost = "l", unreadable = "u", untestable = "u",
          `not applicable` = "-", `M4 only` = ".")
for (r in TESTDEF$test[order(TESTDEF$number)]) {
    x <- CELLS[CELLS$test == r, ]
    cat(sprintf("  %-4s %s\n", r, paste(sprintf("%s:%-6s", paste(x$arm, x$region_set),
        paste0(code[x$verdict], " ", x$kn)), collapse = " ")))
}

TALLY <- do.call(rbind, lapply(ARMS, function(a) do.call(rbind, lapply(RSETS, function(rs) {
    x <- CELLS[in_tally(CELLS$test) & CELLS$arm == a & CELLS$region_set == rs, ]
    data.frame(arm = a, region_set = rs, tests_held = sum(x$verdict == "holds"),
               tests_partial = sum(x$verdict == "partial"), tests_lost = sum(x$verdict == "lost"),
               tests_read = sum(x$n > 0), stringsAsFactors = FALSE)
}))))
TALLY$label <- paste0(TALLY$tests_held, "/", TALLY$tests_read)

# Lead-result cells: one outcome per cell, tests 1-13
lead_cell <- function(r, a, rs, o) {
    if (blk_of(r) == 3 && a != "M4")
        return(data.frame(test = r, arm = a, region_set = rs, outcome = o, k = 0L, n = 0L,
                          unreadable = 0L, verdict = "M4 only", stringsAsFactors = FALSE))
    ent <- REF$entry[REF$arm == a & REF$region_set == rs & REF$outcome == o]
    app <- r == "T1" || isTRUE(ent)
    x <- MODELS[MODELS$test == r & MODELS$arm == a & MODELS$region_set == rs &
                MODELS$outcome == o, ]
    if (!app) x <- x[0, ]
    k <- sum(x$status == "holds"); n <- sum(x$status %in% c("holds", "fails"))
    nu <- sum(x$status == "unreadable")
    never <- !any(MODELS$test == r & MODELS$outcome == o)
    v <- if (never) "not run" else if (!app) "not applicable" else if (n > 0) {
        if (k == n) "holds" else if (k == 0) "lost" else "partial"
    } else if (nu > 0) "unreadable" else "untestable"
    data.frame(test = r, arm = a, region_set = rs, outcome = o, k = k, n = n,
               unreadable = nu, verdict = v, stringsAsFactors = FALSE)
}
HEAD_SETS <- expand.grid(arm = ARMS, region_set = RSETS, stringsAsFactors = FALSE)
HEAD <- do.call(rbind, lapply(seq_len(nrow(HEAD_SETS)), function(i)
    do.call(rbind, lapply(OUTCOMES_ALL, function(o)
        do.call(rbind, lapply(TESTDEF$test[TESTDEF$block %in% c(1, 3)], function(r)
            lead_cell(r, HEAD_SETS$arm[i], HEAD_SETS$region_set[i], o)))))))
# Pool MMSE and CASI for the MERGE_COG tests; per-outcome values kept for the table.
HEAD$k_own <- HEAD$k; HEAD$n_own <- HEAD$n; HEAD$verdict_own <- HEAD$verdict; HEAD$merged <- FALSE
for (r in MERGE_COG) for (a in ARMS) for (rs in RSETS) {
    i <- which(HEAD$test == r & HEAD$arm == a & HEAD$region_set == rs & HEAD$outcome %in% OUTCOMES)
    if (length(i) != length(OUTCOMES)) next
    v <- HEAD$verdict[i]; k <- sum(HEAD$k[i]); n <- sum(HEAD$n[i]); nu <- sum(HEAD$unreadable[i])
    vv <- if (all(v == "not run")) "not run" else if (all(v %in% c("not applicable", "not run")))
        "not applicable" else if (n > 0) { if (k == n) "holds" else if (k == 0) "lost" else "partial" } else
        if (nu > 0) "unreadable" else "untestable"
    HEAD$k[i] <- k; HEAD$n[i] <- n; HEAD$unreadable[i] <- nu; HEAD$verdict[i] <- vv; HEAD$merged[i] <- TRUE
}
# NeoC and NonNeoC combined (region_set 'Both'): k, n and unreadable are summed;
# a region set where test 1 fails contributes nothing.
both_one <- function(x) {
    v <- x$verdict; k <- sum(x$k); n <- sum(x$n); nu <- sum(x$unreadable)
    vv <- if (all(v == "M4 only")) "M4 only" else if (all(v == "not run")) "not run" else
        if (n > 0) { if (k == n) "holds" else if (k == 0) "lost" else "partial" } else
        if (nu > 0) "unreadable" else if (any(v == "untestable")) "untestable" else "not applicable"
    data.frame(test = x$test[1], arm = x$arm[1], region_set = "Both", outcome = x$outcome[1],
               k = k, n = n, unreadable = nu, verdict = vv, k_own = k, n_own = n,
               verdict_own = vv, merged = FALSE,
               regions_read = sum(x$n > 0), stringsAsFactors = FALSE)
}
HEAD$regions_read <- as.integer(HEAD$n > 0)
HB <- HEAD[HEAD$region_set %in% RSETS, ]
HEAD <- rbind(HEAD, do.call(rbind, lapply(split(HB, paste(HB$test, HB$arm, HB$outcome)), both_one)))
rownames(HEAD) <- NULL
HEAD_SETS <- rbind(HEAD_SETS, data.frame(arm = ARMS, region_set = "Both", stringsAsFactors = FALSE))
HEAD$kn <- ifelse(HEAD$n > 0, paste0(HEAD$k, "/", HEAD$n), "")
HEAD$number <- TESTDEF$number[match(HEAD$test, TESTDEF$test)]
HTALLY <- do.call(rbind, lapply(seq_len(nrow(HEAD_SETS)), function(i)
    do.call(rbind, lapply(OUTCOMES_ALL, function(o) {
        # tallies over tests 1-13 without TALLY_EXCLUDE
        x <- HEAD[HEAD$arm == HEAD_SETS$arm[i] & HEAD$region_set == HEAD_SETS$region_set[i] &
                  HEAD$outcome == o & in_tally(HEAD$test), ]
        data.frame(arm = HEAD_SETS$arm[i], region_set = HEAD_SETS$region_set[i], outcome = o,
                   tests_held = sum(x$verdict == "holds"),
                   tests_partial = sum(x$verdict == "partial"),
                   tests_lost = sum(x$verdict == "lost"), tests_read = sum(x$n > 0),
                   # each readable test weighted equally: mean of its k/n
                   mean_pct_per_test = if (any(x$n > 0))
                       100 * mean(x$k[x$n > 0] / x$n[x$n > 0]) else NA,
                   # every model weighted equally: tests with more variants weigh more
                   pooled_k = sum(x$k), pooled_n = sum(x$n),
                   pooled_pct = if (sum(x$n) > 0) 100 * sum(x$k) / sum(x$n) else NA,
                   stringsAsFactors = FALSE)
    }))))
HTALLY$label <- paste0(HTALLY$tests_held, "/", HTALLY$tests_read)
# One decimal; a value below 100 never displays as 100
pct_lab <- function(v) {
    r <- round(v, 1)
    r[!is.na(v) & v < 100 & r >= 100] <- 99.9
    ifelse(is.na(v), "", ifelse(r >= 100, "100", sprintf("%.1f", r)))
}
HTALLY$label_pct <- pct_lab(HTALLY$mean_pct_per_test)
# Tests held in both region sets, per arm and outcome (tests 1-13)
BTALLY <- do.call(rbind, lapply(ARMS, function(a) do.call(rbind, lapply(OUTCOMES_ALL, function(o) {
    x1 <- HEAD[HEAD$arm == a & HEAD$outcome == o & HEAD$region_set == "NeoC" & in_tally(HEAD$test), ]
    x2 <- HEAD[HEAD$arm == a & HEAD$outcome == o & HEAD$region_set == "NonNeoC" & in_tally(HEAD$test), ]
    x2 <- x2[match(x1$test, x2$test), ]
    rd <- x1$n > 0 & x2$n > 0
    data.frame(arm = a, outcome = o, both_held = sum(rd & x1$verdict == "holds" & x2$verdict == "holds"),
               both_read = sum(rd), stringsAsFactors = FALSE)
}))))
BTALLY$label <- paste0(BTALLY$both_held, "/", BTALLY$both_read)
cat(sprintf("\n  Tests held in both region sets (%s):\n", TALLY_SUFFIX))
print(reshape(BTALLY[, c("arm", "outcome", "label")], idvar = "arm", timevar = "outcome", direction = "wide"),
      row.names = FALSE)
cat("\n  Lead results by outcome (tests 1-15; held/read over tests 1-13):\n")
for (i in seq_len(nrow(HTALLY))) {
    h <- HTALLY[i, ]
    x <- HEAD[HEAD$arm == h$arm & HEAD$region_set == h$region_set & HEAD$outcome == h$outcome, ]
    x <- x[order(x$number), ]
    cat(sprintf("  %-6s %-7s %-11s %s | held %d/%d\n", h$arm, h$region_set, h$outcome,
                paste(sprintf("%d:%s", x$number, ifelse(x$kn == "", ifelse(x$verdict == "not run", "nr",
                      ifelse(x$verdict == "M4 only", ".", substr(x$verdict, 1, 1))), x$kn)),
                      collapse = " "), h$tests_held, h$tests_read))
    cat(sprintf("  %-26s mean %% per test %5.1f; pooled %d/%d = %5.1f%%\n", "",
                h$mean_pct_per_test, h$pooled_k, h$pooled_n, h$pooled_pct))
}
cat("\n  Tally, tests 1-13 (held / read):\n")
print(TALLY[, c("arm", "region_set", "tests_held", "tests_partial", "tests_lost", "tests_read")],
      row.names = FALSE)

# == 4. Checks ===============================================================
cat("\n== 4. Checks\n")
add_check("C2 test 1 models", nrow(T$T1) == length(ARMS) * length(RSETS) * length(OUTCOMES_ALL),
          sprintf("%d models (expected %d)", nrow(T$T1),
                  length(ARMS) * length(RSETS) * length(OUTCOMES_ALL)))
add_check("C3 grid complete", nrow(CELLS) == 13 * 8 + 6 * 4 + 2 * 8,
          sprintf("%d cells (expected 144)", nrow(CELLS)))
mo <- CELLS[CELLS$missing_outcomes != "", ]
add_check("C4 applicable outcomes present", nrow(mo) == 0,
          if (nrow(mo) == 0) "every applicable outcome has models" else
              paste(sprintf("%s %s %s: %s", mo$test, mo$arm, mo$region_set, mo$missing_outcomes),
                    collapse = " | "), warn_only = TRUE)

pc <- merge(b1[, c("arm", "region_set", "outcome", "beta_std")],
            P[P$block == "cov4", c("arm", "region_set", "outcome", "beta_std")],
            by = c("arm", "region_set", "outcome"))
add_check("C5 9A vs S14 cov4 beta", nrow(pc) == nrow(b1) &&
          max(abs(pc$beta_std.x - pc$beta_std.y)) < 1e-6,
          sprintf("%d models, max |diff| %.2e", nrow(pc), max(abs(pc$beta_std.x - pc$beta_std.y))))
fc <- merge(b1[, c("arm", "region_set", "outcome", "beta_std")],
            Fb[Fb$predictor_abbrev == "Mod", c("arm", "region_set", "outcome", "beta_std")],
            by = c("arm", "region_set", "outcome"))
add_check("C6 9A vs 10A Module beta", nrow(fc) == nrow(b1) &&
          max(abs(fc$beta_std.x - fc$beta_std.y)) < 1e-3,
          sprintf("%d models, max |diff| %.2e (10A rounded to 4 dp)", nrow(fc),
                  max(abs(fc$beta_std.x - fc$beta_std.y))))

verdict_pair <- function(f1, f2, l1, l2) {
    ifelse(sig(f1) & sig(f2), "Both", ifelse(sig(f1), l1, ifelse(sig(f2), l2, "Neither")))
}
v17 <- F17D
v17$v_re <- ifelse(!is.na(v17$vif) & v17$vif >= VIF_MAX, "Unreadable",
                   verdict_pair(v17$fdr_module, v17$fdr_mito, "Module only", "%MT only"))
m17 <- merge(v17[, c("arm", "region_set", "outcome", "v_re")],
             setNames(V17D[, c("arm", "dataset", "outcome", "verdict_all")],
                      c("arm", "region_set", "outcome", "verdict_all")),
             by = c("arm", "region_set", "outcome"))
add_check(sprintf("C7 %s verdicts recomputed", BP("14D")), nrow(m17) == 72 && all(m17$v_re == m17$verdict_all),
          sprintf("%d of %d agree with Fig19 Verdicts_17D", sum(m17$v_re == m17$verdict_all),
                  nrow(m17)))

pmap <- unique(F18P[, c("panel", "panel_label")])
v18 <- F18A
v18$panel_key <- pmap$panel[match(v18$panel, pmap$panel_label)]
v18$v_re <- verdict_pair(v18$fdr_second, v18$fdr_panel, "second only", "NE panel only")
m18 <- merge(v18[, c("arm", "region_set", "predictor", "panel_key", "outcome", "v_re")],
             setNames(V18A[, c("arm", "dataset", "predictor", "panel", "outcome", "verdict_all")],
                      c("arm", "region_set", "predictor", "panel_key", "outcome", "verdict_all")),
             by = c("arm", "region_set", "predictor", "panel_key", "outcome"))
add_check(sprintf("C8 %s verdicts recomputed", BP("15A")), nrow(m18) == 2304 && all(m18$v_re == m18$verdict_all),
          sprintf("%d of %d agree with Fig19 Verdicts_18A", sum(m18$v_re == m18$verdict_all),
                  nrow(m18)))

# Figure 9A against the Battery 1 unadjusted model, all 72 models
b9 <- B[B$population == "Module post-IP", c("arm", "region_set", "outcome", "n_donors", "FDR")]
o9 <- merge(b9, F15A[, c("arm", "region_set", "outcome", "n_donors", "fdr")],
            by = c("arm", "region_set", "outcome"))
dif <- o9[sig(o9$FDR) != sig(o9$fdr), ]
add_check(sprintf("C9 9A vs %s unadjusted (open item c)", BP("12A")), TRUE,
          sprintf(paste0("9A significant %d of %d; ", BP("12A"), " unadjusted %d of %d; %d differ%s"),
                  sum(sig(o9$FDR)), nrow(o9), sum(sig(o9$fdr)), nrow(o9), nrow(dif),
                  if (nrow(dif) > 0) paste0(": ", paste(sprintf(
                      paste0("%s %s %s (9A FDR %.3f, n %d; ", BP("12A"), " FDR %.3f, n %d)"), dif$arm,
                      dif$region_set, dif$outcome, dif$FDR, as.integer(dif$n_donors.x),
                      dif$fdr, as.integer(dif$n_donors.y)), collapse = "; ")) else ""))

# Not-per-arm values for the figure legend
kap <- F16C[F16C$region_set == "pooled", ]
wd <- F16D[F16D$quantity == "within_median_raw", ]
k_m4 <- kap$kappa[kap$arm == "M4"]; k_hk <- kap$kappa[kap$arm == "HKG8"]
r_mt <- wd$estimate[wd$measure == "MT-CO3 vs cDFG"]
r_ac <- wd$estimate[wd$measure == "ACTB vs cHKG8"]
add_check("C10 legend values", length(k_m4) == 1 && length(k_hk) == 1 &&
          length(r_mt) == 1 && length(r_ac) == 1,
          sprintf("kappa M4 %.3f, HKG8 %.3f; within-donor rs %.3f vs %.3f",
                  k_m4[1], k_hk[1], r_mt[1], r_ac[1]))
FIGREF <- function(x) if (PRE) paste("Figure", x) else paste("Supplementary Figure", BP(x))
FOOT <- sprintf(paste0("Not per arm: Module labels of M4 and HKG8 against PVM_3, Cohen's ",
                       "kappa %.2f and %.2f (", FIGREF("13C"), "); within-donor rs of MT-CO3 ",
                       "with cDFG on PVM_3, median %.2f (ACTB with cHKG8: %.2f; ", FIGREF("13D"),
                       "); external data, Gerrits et al. (", FIGREF("12E"), "-", BP("12H"), ")."),
                k_m4[1], k_hk[1], r_mt[1], r_ac[1])
cat(sprintf("  Legend sentence: %s\n", FOOT))

# Independent recounts: these do not use st() and catch any verdict recycling.
app_m <- function(test) MODELS[MODELS$test == test & MODELS$applicable, ]
n_t1 <- sum(sig(b1$FDR))
add_check("C11 recount test 1", sum(app_m("T1")$status == "holds") == n_t1,
          sprintf("%d holds; direct count %d", sum(app_m("T1")$status == "holds"), n_t1))
d2 <- with_ref(C[C$block %in% c("+ADNC", "+BTC", "+ALL") & C$outcome %in% OUTCOMES, ])
d2 <- d2[!is.na(d2$entry) & d2$entry, ]
n_t2 <- sum(d2$FDR < ALPHA & sign(d2$beta_std) == sign(d2$beta_ref), na.rm = TRUE)
m2 <- app_m("T2"); m2 <- m2[m2$outcome %in% OUTCOMES, ]
add_check("C12 recount test 2", sum(m2$status == "holds") == n_t2,
          sprintf("%d holds; direct count %d", sum(m2$status == "holds"), n_t2))
d9 <- with_ref(F15A[F15A$outcome %in% OUTCOMES_ALL, ])
d9 <- d9[!is.na(d9$entry) & d9$entry, ]
n_t9 <- sum(d9$fdr_adjusted < ALPHA & sign(d9$beta_adjusted) == sign(d9$beta_ref), na.rm = TRUE)
add_check("C13 recount test 9", sum(app_m("T9")$status == "holds") == n_t9,
          sprintf("%d holds; direct count %d", sum(app_m("T9")$status == "holds"), n_t9))
mc <- merge(b1[, c("arm", "region_set", "outcome", "beta_std")],
            M[M$model == "module_only" & M$predictor == "Module post-IP %" &
              M$block == "cov4" & M$boundary == "CMRF",
              c("arm", "region_set", "outcome", "beta_std")],
            by = c("arm", "region_set", "outcome"))
add_check("C15 9A vs 11C CMRF module-only cov4 beta", nrow(mc) == nrow(b1) &&
          max(abs(mc$beta_std.x - mc$beta_std.y)) < 1e-3,
          sprintf("%d of %d models matched, max |diff| %.2e (confirms the sheet M reading)",
                  nrow(mc), nrow(b1), if (nrow(mc)) max(abs(mc$beta_std.x - mc$beta_std.y)) else NA),
          warn_only = TRUE)
# C16: test 14 recounted from the Part 13 table without st().
d16 <- P13M[P13M$spec %in% c("S03", "S07") & P13M$predictor == "module" & P13M$outcome %in% OUTCOMES_ALL, ]
d16$vif_max <- vm13$vif_max[match(paste(d16$spec, d16$dataset, d16$outcome),
                                  paste(vm13$spec, vm13$dataset, vm13$outcome))]
d16 <- merge(d16, REF[REF$arm == "M4", ], by.x = c("dataset", "outcome"), by.y = c("region_set", "outcome"))
d16 <- d16[d16$entry, ]
n16 <- sum(d16$fdr < ALPHA & sign(d16$beta) == sign(d16$beta_ref) & d16$vif_max < VIF_MAX, na.rm = TRUE)
add_check("C16 recount test 14", sum(app_m("T20")$status == "holds") == n16,
          sprintf("%d holds; direct count %d", sum(app_m("T20")$status == "holds"), n16))
# C17: Part 13 S01 (Module post-IP% alone, cov4) is the Figure 9A model, so the
# Part 13 workbook must hold the same M4 coefficients as S12 sheet B.
c17 <- merge(b1[b1$arm == "M4", c("region_set", "outcome", "beta_std")],
             setNames(P13M[P13M$spec == "S01" & P13M$predictor == "module", c("dataset", "outcome", "beta")],
                      c("region_set", "outcome", "beta")), by = c("region_set", "outcome"))
add_check("C17 Part 13 S01 vs 9A (M4)", nrow(c17) == sum(b1$arm == "M4") &&
          max(abs(c17$beta_std - c17$beta)) < 1e-3,
          sprintf("%d of %d models matched, max |diff| %.2e", nrow(c17), sum(b1$arm == "M4"),
                  if (nrow(c17)) max(abs(c17$beta_std - c17$beta)) else NA))
# C18: Part 13b's own-label M4 coefficient (MA) is again the Figure 9A model.
c18 <- merge(b1[b1$arm == "M4" & b1$outcome %in% OUTCOMES, c("region_set", "outcome", "beta_std")],
             setNames(P13S[P13S$readout %in% c("bA_MMSE", "bA_CASI"), c("dataset", "outcome", "M4_own_labels")],
                      c("region_set", "outcome", "own")), by = c("region_set", "outcome"))
add_check("C18 Part 13b own-label MA vs 9A (M4)", nrow(c18) == 4 && max(abs(c18$beta_std - c18$own)) < 1e-3,
          sprintf("%d of 4 models matched, max |diff| %.2e", nrow(c18),
                  if (nrow(c18)) max(abs(c18$beta_std - c18$own)) else NA))
# C19: the both-region tally recounted from the models (G3 applicability, G5 rules,
# MMSE and CASI pooled for MERGE_COG tests), without HEAD.
c19 <- do.call(rbind, lapply(ARMS, function(a) do.call(rbind, lapply(OUTCOMES_ALL, function(o) {
    held <- 0L; rdn <- 0L
    for (r in TESTDEF$test[in_tally(TESTDEF$test)]) {
        outs <- if (r %in% MERGE_COG && o %in% OUTCOMES) OUTCOMES else o
        st2 <- sapply(RSETS, function(rs) {
            x <- MODELS[MODELS$test == r & MODELS$arm == a & MODELS$region_set == rs &
                        MODELS$outcome %in% outs & MODELS$applicable, ]
            k <- sum(x$status == "holds"); n <- sum(x$status %in% c("holds", "fails"))
            c(read = n > 0, held = n > 0 && k == n)
        })
        if (all(st2["read", ] == 1)) { rdn <- rdn + 1L; if (all(st2["held", ] == 1)) held <- held + 1L }
    }
    data.frame(arm = a, outcome = o, held = held, read = rdn)
}))))
m19 <- merge(c19, BTALLY, by = c("arm", "outcome"))
add_check("C19 both-region tally recount", nrow(m19) == length(ARMS) * length(OUTCOMES_ALL) &&
          all(m19$held == m19$both_held & m19$read == m19$both_read),
          sprintf("%d of %d arm x outcome cells agree", sum(m19$held == m19$both_held & m19$read == m19$both_read),
                  nrow(m19)))
# C20: the combined cells recounted from the models (G3 applicability, G5 rules),
# without HEAD.
hb20 <- HEAD[HEAD$region_set == "Both" & HEAD$verdict != "M4 only", ]
c20 <- mapply(function(r, a, o) {
    x <- MODELS[MODELS$test == r & MODELS$arm == a & MODELS$outcome == o &
                MODELS$region_set %in% RSETS & MODELS$applicable, ]
    c(sum(x$status == "holds"), sum(x$status %in% c("holds", "fails")))
}, hb20$test, hb20$arm, hb20$outcome)
ok20 <- c20[1, ] == hb20$k & c20[2, ] == hb20$n
add_check("C20 combined cells recount", all(ok20),
          sprintf("%d of %d combined cells agree (k and n from the models)", sum(ok20), length(ok20)))
# C21: in the combined figure a full circle read in both region sets is a test held
# in both; the both-region tally must equal that count of drawn cells.
c21 <- do.call(rbind, lapply(ARMS, function(a) do.call(rbind, lapply(OUTCOMES_ALL, function(o) {
    x <- HEAD[HEAD$region_set == "Both" & HEAD$arm == a & HEAD$outcome == o & in_tally(HEAD$test), ]
    data.frame(arm = a, outcome = o, held = sum(x$regions_read == 2 & x$verdict == "holds"),
               read = sum(x$regions_read == 2))
}))))
m21 <- merge(c21, BTALLY, by = c("arm", "outcome"))
add_check("C21 combined glyphs vs both-region tally", nrow(m21) == length(ARMS) * length(OUTCOMES_ALL) &&
          all(m21$held == m21$both_held & m21$read == m21$both_read),
          sprintf("%d of %d arm x outcome cells agree", sum(m21$held == m21$both_held & m21$read == m21$both_read),
                  nrow(m21)))
# C22: TALLY_EXCLUDE rederived from the cells: block-1 tests (other than test 1)
# with no reading (unreadable, untestable or not run) in at least one M4, PVM_3 or
# PVM_1 cell of NeoC or NonNeoC.
h22 <- HEAD[HEAD$region_set %in% RSETS & HEAD$arm %in% TALLY_ARMS & blk_of(HEAD$test) == 1 &
            HEAD$test != "T1" & HEAD$verdict %in% c("unreadable", "untestable", "not run"), ]
ex22 <- TESTDEF$test[TESTDEF$test %in% h22$test]
add_check("C22 tests left out of the tallies", setequal(ex22, TALLY_EXCLUDE),
          sprintf("derived {%s}; TALLY_EXCLUDE {%s}", paste(ex22, collapse = ", "),
                  paste(TALLY_EXCLUDE, collapse = ", ")))
cat("  Cells without a reading, tests left out of the tallies (for the legend):\n")
for (r in TALLY_EXCLUDE) {
    x <- HEAD[HEAD$test == r & HEAD$region_set %in% RSETS & HEAD$arm %in% TALLY_ARMS, ]
    tb <- table(paste(x$arm, ifelse(x$outcome %in% OUTCOMES, "cognition", "pathology")), x$verdict)
    cat(sprintf("  %s (test %d):\n", r, TESTDEF$number[TESTDEF$test == r])); print(tb)
}
var_n <- tapply(seq_len(nrow(MODELS)), MODELS$test, function(i)
    length(unique(MODELS$status[i])))
add_check("C14 verdict spread", TRUE, paste(sprintf("%s %d", names(var_n), var_n),
                                            collapse = ", "))

# == 5. Figures ==============================================================
cat("\n== 5. Figures\n")
PT <- function(x) x / .pt
LEAD_LIGHT <- c(M4 = "#ECE9F6", PVM_3 = "#E1F2EC", PVM_1 = "#FDF0DD", HKG8 = "#EEEEEE")
# Darker text shades so the lead headers and tallies stay legible on the light bands
LEAD_TEXT  <- c(M4 = "#3E3480", PVM_3 = "#137A5E", PVM_1 = "#B8690A", HKG8 = "#555555")
# Solid title-box fills (white text); PVM_1 and HKG8 deepened for contrast with white
HEAD_FILL  <- c(M4 = "#3E3480", PVM_3 = "#1A9C78", PVM_1 = "#E08A1E", HKG8 = "#7A7A7A")

circle <- function(cx, cy, r, from = 0, to = 2 * pi, n = 60) {
    t <- seq(from, to, length.out = n)
    data.frame(x = cx + r * cos(t), y = cy + r * sin(t))
}
# Fill fraction for quartile glyphs: all, >=75%, 25% to <75%, <25%, none
quart_frac <- function(k, n) {
    f <- k / n
    ifelse(is.na(f), NA, ifelse(f >= 1, 1, ifelse(f >= 0.75, 0.75, ifelse(f >= 0.25, 0.5,
           ifelse(f > 0, 0.25, 0)))))
}
glyph_polys <- function(id, verdict, cx, cy, col, rad, qf = NA) {
    out <- list()
    if (verdict %in% c("holds", "partial", "lost", "not applicable")) {
        if (verdict == "holds") {
            p <- circle(cx, cy, rad); p$id <- paste0(id, "f"); p$fill <- col; out$f <- p
        }
        if (verdict == "partial" && !is.na(qf)) {
            # wedge filled anticlockwise from 3 o'clock (angle 0, unit-circle convention)
            arc <- if (qf == 0.5) circle(cx, cy, rad, pi / 2, 3 * pi / 2) else
                circle(cx, cy, rad, 0, 2 * pi * qf)
            p <- rbind(data.frame(x = cx, y = cy), arc,
                       data.frame(x = cx, y = cy))
            p$id <- paste0(id, "q"); p$fill <- col; out$h <- p
        }
        if (verdict == "partial" && is.na(qf)) {
            p <- rbind(circle(cx, cy, rad, pi / 2, 3 * pi / 2), data.frame(x = cx, y = cy + rad))
            p$id <- paste0(id, "h"); p$fill <- col; out$h <- p
        }
        o <- circle(cx, cy, rad); o$id <- paste0(id, "o"); o$fill <- NA; out$o <- o
    }
    res <- do.call(rbind, out)
    if (!is.null(res)) res$col <- col
    res
}

# Layout specifications. Grid positions are in data units; the page width
# sets the scale. Font sizes are in pt.
SPEC_PORTRAIT <- list(
    name = "portrait", base = "Fig17_Summary", PAGE_W = 7.5, block2 = TRUE,
    leads = LEADS[LEADS$arm == HL_ARM & LEADS$region_set == HL_RS, ],
    X_NUM = 0.9, X_LAB = 1.25, X_BAT = 13.2, X_FIG = 14.8, LEAD_X0 = 16.6, LEAD_SUB = 1.6,
    LEAD_GAP = 0, X0 = 20.4, ARM_W = 3.55, SUB_W = 1.7, RH = 1.22, RAD = 0.30,
    TXT_LAB = 8.5, TXT_GRP = 9, TXT_SUB = 7.5, TXT_RS = 6.8, TXT_KN = 6.5, TXT_HEAD = 8.5,
    KEY_X = c(0.6, 7.6, 14.6),
    FB_H = 0.74, FB_GAP = 0.22, FB_ROWGAP = 0.30, FB_BAND = 0.17, FB_TXT = 7.5,
    FB_STAGE = 6.5, FB_TAG = 6, FB_TOP = 0.30, FB_HEAD = 9,
    flow_text = c("SEA-AD immune\nsnRNA-seq: 240,651\nnuclei, 84 donors",
                  "M4: Mic.13-like\nPVM_3 cells;\nMT-CO3 + 8 DFGs",
                  "GAM fit and\nSKCS-core\nper cell",
                  "Inflection point\n(IP) and IPRS",
                  "Copula space\n(S_i, K_i)",
                  "CMRF band and\nRSk-fCV gate:\nModule pre/post-IP",
                  "Validation:\nGAM refit, UMAP",
                  "Module post-IP%\nper donor, 4 arms\nx 2 region sets",
                  "Clinical association,\nbenchmark,\nboundaries",
                  "Batteries 1-5:\ncontamination,\nartifact, %MT, SA"),
    feeds = c(rep("", 8), "Tests 1-8", "Tests 9-21"))
SPEC_LANDSCAPE <- list(
    name = "landscape", base = "Fig17_Summary_landscape", PAGE_W = 11, block2 = FALSE,
    leads = LEADS,
    X_NUM = 1.0, X_LAB = 1.4, X_BAT = 14.0, X_FIG = 15.8, LEAD_X0 = 17.95, LEAD_SUB = 1.8,
    LEAD_GAP = 0.6, X0 = 31.0, ARM_W = 3.6, SUB_W = 1.75, RH = 1.2, RAD = 0.31,
    TXT_LAB = 9.5, TXT_GRP = 10, TXT_SUB = 8, TXT_RS = 7.5, TXT_KN = 7, TXT_HEAD = 9.5,
    KEY_X = c(0.6, 9.0, 17.4),
    FB_H = 0.72, FB_GAP = 0.32, FB_ROWGAP = 0.28, FB_BAND = 0.20, FB_TXT = 9,
    FB_STAGE = 7.5, FB_TAG = 7, FB_TOP = 0.34, FB_HEAD = 10,
    flow_text = c("SEA-AD immune snRNA-seq\n240,651 nuclei, 84 donors",
                  "M4: Mic.13-like PVM_3 cells\nMT-CO3 + 8 DFGs",
                  "GAM fit and SKCS-core\nper cell",
                  "Inflection point (IP)\nand IPRS",
                  "Copula space\n(S_i, K_i)",
                  "CMRF band and RSk-fCV gate\nModule pre-IP and post-IP",
                  "Validation\nGAM refit, UMAP",
                  "Module post-IP% per donor\n4 arms x 2 region sets",
                  "Clinical association,\nbenchmark, boundaries",
                  "Batteries 1-5: contamination,\nartifact, %MT, SA donors"),
    feeds = c(rep("", 8), "Tests 1-8", "Tests 9-13"))

# Fields are replaced one by one: modifyList() would try to merge the 'leads'
# data frames column by column.
SPEC_PATHOLOGY <- SPEC_LANDSCAPE
PATH_UPDATE <- list(
    name = "landscape_pathology", base = "Fig17_Summary_landscape_pathology",
    leads = data.frame(arm = ARMS, region_set = "NeoC", stringsAsFactors = FALSE),
    lead_outcomes = OUTCOMES_ALL, arm_grid = FALSE, battery_col = FALSE,
    X_FIG = 14.5, LEAD_X0 = 17.2, LEAD_SUB = 1.7, LEAD_GAP = 0.5, RH = 1.2, GRID_TOP = 0.45,
    KEY_X = c(0.6, 11.5, 22.4),
    group_override = c("Batteries 1, 3, 4 and 5: the clinical association retested" =
                           "Batteries 1, 3, 4 & 5: Clinical association rechallenged",
                       "Against alternative explanations" =
                           "Module vs non-Module, competitors, ablations"),
    label_override = c(T1 = "Associated at cov4 (age, sex, PMI, APOE \u03b54)",
                       T2 = "Holds with the other domain in the model",
                       T5 = "Module outperforms non-Module",
                       T8 = "Holds for all four boundaries (ablations)"),
    fig_override = c(T2 = "9B-C", T3 = "S14A-B", T8 = "11A, 11C"),
    TXT_LAB = 10.5, TXT_GRP = 10.5, TXT_HEAD = 10.5, TXT_LEG = 10, TXT_FIGCOL = 10.5,
    TXT_PANEL = 12, TXT_NUM = 10, TXT_OUT = 8, TXT_TALLY = 9.5,
    MERGE_NOT_RUN = TRUE, MARK_CHAR = "-", MARK_SIZE = 16, MARK_COL = "grey20",
    MARK_LABEL = "Not assessable",
    FIGCOL_COL = "black", LEAD_BOX = TRUE, flow_style = "metro", TXT_KN = 8,
    HEAD_BOX_FILL = "#002060", TXT_STRIP = 12, HEAD_DROP = 0.45, STRIP_FILL = "white",
    STRIP_LW = 0.3, GLYPH_COL = "grey20", STRIP_COL = "grey20", HEAD_BOX_EXTRA = 0.15,
    HEAD_LABEL_UP = 0.08,
    show_title = FALSE, QUARTILES = TRUE, TALLY_PCT = TRUE, TALLY_BOTH = FALSE, QKEY_X = c(0.6, 8.6, 17.4, 25.6, 32.4),
    QLINE2_X = c(9.0, 24.6))
for (k in names(PATH_UPDATE)) SPEC_PATHOLOGY[[k]] <- PATH_UPDATE[[k]]
OUT_SHORT <- c("MMSE score" = "MMSE", "CASI score" = "CASI", "Thal phase" = "Thal",
               "Braak stage" = "Braak")

build_grid <- function(S) {
    if (is.null(S$lead_outcomes)) S$lead_outcomes <- OUTCOMES
    if (is.null(S$arm_grid)) S$arm_grid <- TRUE
    if (is.null(S$battery_col)) S$battery_col <- TRUE
    # Text-size fields; unset ones fall back to the shared sizes of the spec
    for (k in c("TXT_LEG", "TXT_FIGCOL", "TXT_NUM", "TXT_OUT"))
        if (is.null(S[[k]])) S[[k]] <- S$TXT_SUB
    if (is.null(S$TXT_TALLY)) S$TXT_TALLY <- S$TXT_LAB
    if (is.null(S$TXT_PANEL)) S$TXT_PANEL <- 11
    if (is.null(S$MERGE_NOT_RUN)) S$MERGE_NOT_RUN <- FALSE
    if (is.null(S$MARK_CHAR)) S$MARK_CHAR <- "\u00D7"
    if (is.null(S$MARK_SIZE)) S$MARK_SIZE <- 10
    if (is.null(S$MARK_COL)) S$MARK_COL <- "grey45"
    if (is.null(S$MARK_LABEL)) S$MARK_LABEL <- "Unreadable, untestable or inconclusive"
    if (is.null(S$FIGCOL_COL)) S$FIGCOL_COL <- "grey30"
    if (is.null(S$LEAD_BOX)) S$LEAD_BOX <- FALSE
    LO <- S$lead_outcomes; nk <- length(LO)
    tdef <- TESTDEF
    if (!is.null(S$group_override))
        for (k in names(S$group_override)) tdef$group[tdef$group == k] <- S$group_override[[k]]
    if (!is.null(S$label_override))
        tdef$label[match(names(S$label_override), tdef$test)] <- unname(S$label_override)
    if (!is.null(S$fig_override))
        tdef$fig[match(names(S$fig_override), tdef$test)] <- unname(S$fig_override)
    tests_shown <- tdef[if (S$block2) TRUE else tdef$block %in% c(1, 3), ]
    # Block 3 (M4 only) sits after block 1 and above the tally
    tests_shown <- tests_shown[order(match(tests_shown$block, c(1, 3, 2)), tests_shown$number), ]
    GRID_TOP <- if (is.null(S$GRID_TOP)) 0.6 else S$GRID_TOP
    nl <- nrow(S$leads)
    # pairs are spaced by one pair width plus LEAD_GAP
    lead_x <- function(i, o) S$LEAD_X0 + (i - 1) * ((nk - 1) * S$LEAD_SUB + 1.9 + S$LEAD_GAP) +
        (match(o, LO) - 1) * S$LEAD_SUB
    col_x <- function(a, rs) S$X0 + (match(a, ARMS) - 1) * S$ARM_W +
        ifelse(rs == "both", S$SUB_W / 2, (match(rs, RSETS) - 1) * S$SUB_W)
    arm_mid <- function(a) S$X0 + (match(a, ARMS) - 1) * S$ARM_W + S$SUB_W / 2
    X_MAX <- if (S$arm_grid) S$X0 + 3 * S$ARM_W + S$SUB_W + 1.0 else
        lead_x(nl, LO[nk]) + 0.95 + 0.4
    RH <- S$RH; RAD <- S$RAD

    lines <- list(); y <- -0.2
    add_line <- function(type, ...) { lines[[length(lines) + 1]] <<- list(type = type, y = y, ...) }
    y <- y - 1.4; add_line("colhead")
    groups <- unique(tests_shown$group)
    for (gname in groups) {
        if (tests_shown$block[tests_shown$group == gname][1] == 2) {
            y <- y - RH - 0.25; add_line("tally"); y <- y - 0.3
        }
        y <- y - 1.1; add_line("group", text = gname)
        for (r in tests_shown$test[tests_shown$group == gname]) {
            y <- y - RH; add_line("row", test = r)
        }
    }
    # Tally rows in display order; the first carries the separator line
    TALS <- if (!is.null(S$TALLIES)) S$TALLIES else
        c("tally", if (isTRUE(S$TALLY_PCT)) "tally_pct", if (isTRUE(S$TALLY_BOTH)) "tally_both")
    if (!S$block2) for (j in seq_along(TALS)) {
        y <- y - if (j == 1) RH + 0.25 else RH * 0.85; add_line(TALS[j], first = j == 1)
    }
    y <- y - 1.35; add_line("legend")
    Y_MIN <- y - 1.5 - if (!is.null(S$LEGEND_LINE3)) 1.05 else 0
    row_y <- unlist(lapply(Filter(function(l) l$type == "row", lines),
                           function(l) setNames(l$y, l$test)))

    P <- CELLS[CELLS$test %in% tests_shown$test, ]
    if (!S$arm_grid) P <- P[0, ]
    P$x <- col_x(P$arm, P$region_set); P$y <- row_y[P$test]; P$gy <- P$y + 0.20
    P$col <- if (is.null(S$GLYPH_COL)) unname(ARM_COL[P$arm]) else rep(S$GLYPH_COL, nrow(P))
    H <- merge(HEAD[HEAD$outcome %in% LO, ], data.frame(arm = S$leads$arm, region_set = S$leads$region_set,
                                lead_i = seq_len(nl), stringsAsFactors = FALSE),
               by = c("arm", "region_set"))
    H$x <- lead_x(H$lead_i, H$outcome); H$y <- row_y[H$test]; H$gy <- H$y + 0.20
    H$col <- if (is.null(S$GLYPH_COL)) unname(ARM_COL[H$arm]) else rep(S$GLYPH_COL, nrow(H))
    # A pooled MMSE + CASI cell is drawn once, centred over the two columns
    if (all(OUTCOMES %in% LO)) {
        mm <- H$merged & H$outcome == OUTCOMES[1]
        H$x[mm] <- (lead_x(H$lead_i[mm], OUTCOMES[1]) + lead_x(H$lead_i[mm], OUTCOMES[2])) / 2
        H <- H[!(H$merged & H$outcome == OUTCOMES[2]), ]
    }
    polys <- do.call(rbind, c(
        lapply(seq_len(nrow(P)), function(i)
            glyph_polys(paste0("c", i), P$verdict[i], P$x[i], P$gy[i], P$col[i], RAD)),
        lapply(seq_len(nrow(H)), function(i)
            glyph_polys(paste0("h", i), H$verdict[i], H$x[i], H$gy[i], H$col[i], RAD,
                        if (isTRUE(S$QUARTILES)) quart_frac(H$k[i], H$n[i]) else NA))))
    marks <- rbind(P[P$verdict %in% c("unreadable", "untestable"), c("x", "gy")],
                   H[H$verdict %in% c("unreadable", "untestable",
                                      if (S$MERGE_NOT_RUN) "not run"), c("x", "gy")])
    marks$sym <- rep(S$MARK_CHAR, nrow(marks))
    KN <- rbind(P[P$kn != "", c("x", "y", "kn")], H[H$kn != "", c("x", "y", "kn")])
    # One 'M4 only' note per non-M4 column or strip, centred on the block 3 rows
    b3 <- tests_shown$test[tests_shown$block == 3]
    M4ONLY <- data.frame(x = numeric(0), y = numeric(0))
    if (length(b3)) {
        yb3 <- mean(row_y[b3]) + 0.20
        xo <- if (S$arm_grid) sapply(setdiff(ARMS, "M4"), arm_mid) else {
            i_o <- which(S$leads$arm != "M4")
            (lead_x(i_o, LO[1]) + lead_x(i_o, LO[nk])) / 2 }
        if (length(xo)) M4ONLY <- data.frame(x = unname(xo), y = yb3)
    }

    ycol <- Filter(function(l) l$type == "colhead", lines)[[1]]$y
    HDROP <- if (is.null(S$HEAD_DROP)) 0 else S$HEAD_DROP     # headers moved down
    HBOX  <- if (is.null(S$HEAD_BOX_EXTRA)) 0 else S$HEAD_BOX_EXTRA  # taller title box
    HUP   <- if (is.null(S$HEAD_LABEL_UP)) 0 else S$HEAD_LABEL_UP
    tl <- Filter(function(l) grepl("^tally", l$type), lines)
    ytal_count <- tl[[1]]$y                 # first tally row
    ytal <- tl[[length(tl)]]$y              # last tally row
    BAND <- data.frame(arm = S$leads$arm, region_set = S$leads$region_set,
                       xmin = lead_x(seq_len(nl), LO[1]) - 0.95,
                       xmax = lead_x(seq_len(nl), LO[nk]) + 0.95,
                       ymin = ytal - 0.42, ymax = ycol + 0.98 - HDROP + HBOX, stringsAsFactors = FALSE)
    BAND$fill <- if (is.null(S$STRIP_FILL)) unname(LEAD_LIGHT[BAND$arm]) else S$STRIP_FILL
    BAND$col <- if (is.null(S$STRIP_COL)) unname(ARM_COL[BAND$arm]) else S$STRIP_COL
    BAND$tcol <- unname(LEAD_TEXT[BAND$arm])

    g <- ggplot() + coord_fixed(xlim = c(0, X_MAX), ylim = c(Y_MIN, GRID_TOP), expand = FALSE) +
        theme_void(base_family = FONT) + theme(plot.margin = margin(0, 0, 0, 0)) +
        geom_rect(data = BAND, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                                   fill = fill, colour = col),
                  linewidth = if (is.null(S$STRIP_LW)) 0.6 else S$STRIP_LW)
    if (any(LO %in% PATH_OUT) && any(LO %in% OUTCOMES)) {
        j <- max(which(LO %in% OUTCOMES))
        DIV <- data.frame(x = lead_x(seq_len(nl), LO[j]) + S$LEAD_SUB / 2,
                          ymin = ytal_count + 0.68, ymax = ycol - 0.15 - HDROP, col = BAND$col)
        g <- g + geom_segment(data = DIV, aes(x = x, xend = x, y = ymin, yend = ymax,
                                              colour = col), linetype = "dotted",
                              linewidth = 0.45)
    }
    for (l in lines) {
        if (isTRUE(l$first))
            g <- g + annotate("segment", x = 0.1, xend = X_MAX - 0.1, y = l$y + 0.68,
                              yend = l$y + 0.68, colour = "grey30", linewidth = 0.4)
        if (l$type == "colhead") {
            hy <- l$y + 0.5 - HDROP
            hd <- data.frame(arm = ARMS, xmin = S$X0 + (seq_along(ARMS) - 1) * S$ARM_W - 0.6,
                             xmax = S$X0 + (seq_along(ARMS) - 1) * S$ARM_W + S$SUB_W + 0.6,
                             ymin = hy - 0.45, ymax = hy + 0.45, y = hy,
                             fill = unname(ARM_COL[ARMS]),
                             tcol = ifelse(ARMS == "HKG8", "black", "white"))
            sh <- expand.grid(arm = ARMS, rs = RSETS, stringsAsFactors = FALSE)
            sh$x <- col_x(sh$arm, sh$rs); sh$y <- hy - 0.95
            lh <- data.frame(x = (BAND$xmin + BAND$xmax) / 2, y = hy + HUP,
                             label = ifelse(BAND$region_set == "Both", BAND$arm,
                                            paste(BAND$arm, BAND$region_set)), col = BAND$tcol)
            if (S$LEAD_BOX) {
                # Solid title box spanning the strip and joined to its top edge
                TB <- data.frame(xmin = BAND$xmin, xmax = BAND$xmax, ymin = hy - 0.42,
                                 ymax = BAND$ymax,
                                 fill = if (is.null(S$HEAD_BOX_FILL)) unname(HEAD_FILL[BAND$arm])
                                        else S$HEAD_BOX_FILL)
                g <- g + geom_rect(data = TB, aes(xmin = xmin, xmax = xmax, ymin = ymin,
                                                  ymax = ymax, fill = fill, colour = fill),
                                   linewidth = 0.6)
                lh$col <- "white"
            }
            ls <- expand.grid(i = seq_len(nl), o = LO, stringsAsFactors = FALSE)
            ls$x <- lead_x(ls$i, ls$o); ls$y <- hy - 0.95
            ls$label <- unname(OUT_SHORT[ls$o])
            if (S$arm_grid)
                g <- g + geom_rect(data = hd, aes(xmin = xmin, xmax = xmax, ymin = ymin,
                                                  ymax = ymax, fill = fill), colour = NA) +
                    geom_text(data = hd, aes(x = (xmin + xmax) / 2, y = y, label = arm,
                                             colour = tcol),
                              fontface = "bold", size = PT(S$TXT_HEAD), family = FONT) +
                    geom_text(data = sh, aes(x = x, y = y, label = rs), size = PT(S$TXT_RS),
                              family = FONT)
            hx <- if (S$battery_col) c(S$X_BAT, S$X_FIG) else S$X_FIG
            hl <- if (S$battery_col) c("Battery", "Figure") else "Figure"
            g <- g + annotate("text", x = hx, y = hy - 0.95, label = hl, size = PT(S$TXT_FIGCOL),
                              family = FONT, colour = S$FIGCOL_COL) +
                geom_text(data = lh, aes(x = x, y = y, label = label, colour = col),
                          fontface = "bold",
                          size = PT(if (is.null(S$TXT_STRIP)) S$TXT_HEAD else S$TXT_STRIP),
                          family = FONT) +
                geom_text(data = ls, aes(x = x, y = y, label = label), fontface = "bold",
                          size = PT(S$TXT_OUT), family = FONT)
        }
        if (l$type == "group") {
            g <- g + annotate("text", x = 0.1, y = l$y, label = l$text, hjust = 0,
                              fontface = "bold", size = PT(S$TXT_GRP), family = FONT) +
                annotate("segment", x = 0.1, xend = X_MAX - 0.1, y = l$y - 0.45,
                         yend = l$y - 0.45, colour = "grey70", linewidth = 0.3)
        }
        if (l$type == "row") {
            rr <- tdef[tdef$test == l$test, ]
            ty <- l$y + 0.20
            g <- g + annotate("text", x = S$X_NUM, y = ty, label = rr$number, hjust = 1,
                              size = PT(S$TXT_NUM), family = FONT, colour = "grey40",
                              fontface = "bold") +
                annotate("text", x = S$X_LAB, y = ty, label = rr$label, hjust = 0,
                         size = PT(S$TXT_LAB), family = FONT) +
                annotate("text", x = S$X_BAT, y = ty,
                         label = if (S$battery_col) rr$battery else "", size = PT(S$TXT_FIGCOL),
                         family = FONT, colour = S$FIGCOL_COL) +
                annotate("text", x = S$X_FIG, y = ty, label = rr$fig, size = PT(S$TXT_FIGCOL),
                         family = FONT, colour = S$FIGCOL_COL)
        }
        if (l$type == "tally") {
            td <- transform(TALLY, x = col_x(arm, region_set), y = l$y + 0.17)
            if (!S$arm_grid) td <- td[0, ]
            ht <- merge(HTALLY[HTALLY$outcome %in% LO, ], data.frame(arm = S$leads$arm, region_set = S$leads$region_set,
                                           lead_i = seq_len(nl), stringsAsFactors = FALSE),
                        by = c("arm", "region_set"))
            ht$x <- lead_x(ht$lead_i, ht$outcome); ht$y <- l$y + 0.17
            ht$col <- if (is.null(S$GLYPH_COL)) unname(LEAD_TEXT[ht$arm]) else rep(S$GLYPH_COL, nrow(ht))
            g <- g + annotate("text", x = S$X_LAB, y = l$y + 0.17,
                         label = sprintf("Tests held / tests read (%s)", TALLY_SUFFIX), hjust = 0,
                         fontface = "bold", size = PT(S$TXT_LAB), family = FONT) +
                geom_text(data = td, aes(x = x, y = y, label = label), fontface = "bold",
                          size = PT(S$TXT_TALLY), family = FONT) +
                geom_text(data = ht, aes(x = x, y = y, label = label, colour = col),
                          fontface = "bold", size = PT(S$TXT_TALLY), family = FONT)
        }
        if (l$type == "tally_both") {
            hb <- merge(BTALLY[BTALLY$outcome %in% LO, ],
                        data.frame(arm = S$leads$arm, lead_i = seq_len(nl), stringsAsFactors = FALSE), by = "arm")
            hb$x <- lead_x(hb$lead_i, hb$outcome); hb$y <- l$y + 0.17
            hb$col <- if (is.null(S$GLYPH_COL)) unname(LEAD_TEXT[hb$arm]) else rep(S$GLYPH_COL, nrow(hb))
            g <- g + annotate("text", x = S$X_LAB, y = l$y + 0.17,
                              label = sprintf("Tests held in both region sets (%s)", TALLY_SUFFIX), hjust = 0,
                              fontface = "bold", size = PT(S$TXT_LAB), family = FONT) +
                geom_text(data = hb, aes(x = x, y = y, label = label, colour = col),
                          fontface = "bold", size = PT(S$TXT_TALLY), family = FONT)
        }
        if (l$type == "tally_pct") {
            hp <- merge(HTALLY[HTALLY$outcome %in% LO, ],
                        data.frame(arm = S$leads$arm, region_set = S$leads$region_set,
                                   lead_i = seq_len(nl), stringsAsFactors = FALSE),
                        by = c("arm", "region_set"))
            hp$x <- lead_x(hp$lead_i, hp$outcome); hp$y <- l$y + 0.17
            hp$col <- if (is.null(S$GLYPH_COL)) unname(LEAD_TEXT[hp$arm]) else rep(S$GLYPH_COL, nrow(hp))
            if (!is.null(S$BOX_MIN)) {
                # every mean % at or above BOX_MIN (the '>=75% of models hold' cut),
                # judged on the displayed value
                hp$shown <- suppressWarnings(as.numeric(hp$label_pct))
                bx <- hp[!is.na(hp$shown) & hp$shown >= S$BOX_MIN, ]
                if (nrow(bx))
                    g <- g + geom_rect(data = bx, aes(xmin = x - S$LEAD_SUB * 0.42, xmax = x + S$LEAD_SUB * 0.42,
                                                      ymin = y - 0.36, ymax = y + 0.36),
                                       fill = NA, colour = S$BOX_COL, linewidth = 0.5)
                BOXED[[S$name]] <<- bx[, c("arm", "region_set", "outcome", "label_pct")]
            }
            g <- g + annotate("text", x = S$X_LAB, y = l$y + 0.17,
                              label = "Mean % of models held per test (same tests)", hjust = 0,
                              fontface = "bold", size = PT(S$TXT_LAB), family = FONT) +
                geom_text(data = hp, aes(x = x, y = y, label = label_pct, colour = col),
                          fontface = "bold", size = PT(S$TXT_TALLY), family = FONT)
        }
        if (l$type == "legend" && isTRUE(S$QUARTILES)) {
            ly <- l$y
            kx <- S$QKEY_X
            KCOL <- if (is.null(S$GLYPH_COL)) "grey30" else S$GLYPH_COL
            kv <- c("holds", "partial", "partial", "partial", "lost")
            kq <- c(NA, 0.75, 0.5, 0.25, NA)
            kp <- do.call(rbind, lapply(1:5, function(i)
                glyph_polys(paste0("k", i), kv[i], kx[i], ly,
                            if (is.null(S$GLYPH_COL)) "grey30" else S$GLYPH_COL, RAD, kq[i])))
            g <- g + annotate("segment", x = 0.1, xend = X_MAX - 0.1, y = ly + 0.8,
                              yend = ly + 0.8, colour = "grey30", linewidth = 0.4) +
                geom_polygon(data = kp[!is.na(kp$fill), ], aes(x, y, group = id),
                             fill = KCOL, colour = NA) +
                geom_path(data = kp[is.na(kp$fill), ], aes(x, y, group = id),
                          colour = KCOL, linewidth = 0.4) +
                annotate("text", x = kx + 0.6, y = ly,
                         label = c("All models hold", "\u226575% of models hold",
                                   "25% to <75% hold", "<25% hold", "No model holds"),
                         hjust = 0, size = PT(S$TXT_LEG), family = FONT) +
                annotate("text", x = kx[1], y = ly - 1.05, label = S$MARK_CHAR,
                         size = PT(S$MARK_SIZE), family = FONT, colour = S$MARK_COL) +
                annotate("text", x = kx[1] + 0.6, y = ly - 1.05, label = S$MARK_LABEL,
                         hjust = 0, size = PT(S$TXT_LEG), family = FONT) +
                annotate("text", x = S$QLINE2_X[1], y = ly - 1.05,
                         label = "Open circle without k/n: nothing to retest (HKG8)",
                         hjust = 0, size = PT(S$TXT_LEG), family = FONT) +
                annotate("text", x = S$QLINE2_X[2], y = ly - 1.05,
                         label = "k/n: models that hold / readable models (variants)",
                         hjust = 0, size = PT(S$TXT_LEG), family = FONT)
            if (!is.null(S$LEGEND_LINE3))
                g <- g + annotate("text", x = kx[1] - 0.3, y = ly - 2.1, label = S$LEGEND_LINE3,
                                  hjust = 0, size = PT(S$TXT_LEG), family = FONT)
        }
        if (l$type == "legend" && !isTRUE(S$QUARTILES)) {
            ly <- l$y
            key <- data.frame(v = c("holds", "partial", "lost"), x = S$KEY_X)
            kp <- do.call(rbind, lapply(1:3, function(i)
                glyph_polys(paste0("k", i), key$v[i], key$x[i], ly, "grey30", RAD)))
            g <- g + annotate("segment", x = 0.1, xend = X_MAX - 0.1, y = ly + 0.8,
                              yend = ly + 0.8, colour = "grey30", linewidth = 0.4) +
                geom_polygon(data = kp[!is.na(kp$fill), ], aes(x, y, group = id),
                             fill = "grey30", colour = NA) +
                geom_path(data = kp[is.na(kp$fill), ], aes(x, y, group = id),
                          colour = "grey30", linewidth = 0.4) +
                annotate("text", x = key$x + 0.6, y = ly,
                         label = c("All models hold", "Some models hold",
                                   "No model holds; without k/n: nothing to retest (HKG8)"),
                         hjust = 0, size = PT(S$TXT_LEG), family = FONT) +
                annotate("text", x = S$KEY_X[1], y = ly - 1.05, label = S$MARK_CHAR,
                         size = PT(S$MARK_SIZE), family = FONT, colour = S$MARK_COL) +
                annotate("text", x = S$KEY_X[1] + 0.6, y = ly - 1.05,
                         label = S$MARK_LABEL, hjust = 0,
                         size = PT(S$TXT_LEG), family = FONT) +
                annotate("text", x = S$KEY_X[3], y = ly - 1.05,
                         label = if (nk > 2) "k/n: models that hold / readable models (variants)"
                                 else "k/n: models that hold / readable models (MMSE, CASI x variants)",
                         hjust = 0, size = PT(S$TXT_LEG), family = FONT)
            nr <- sort(unique(H$number[H$verdict == "not run"]))
            if (length(nr) > 0 && !S$MERGE_NOT_RUN)
                g <- g + annotate("text", x = S$KEY_X[2], y = ly - 1.05,
                                  label = sprintf("Blank: not run for pathology (test %s)",
                                                  paste(nr, collapse = ", ")),
                                  hjust = 0, size = PT(S$TXT_LEG), family = FONT)
        }
    }
    g <- g +
        geom_polygon(data = polys[!is.na(polys$fill), ], aes(x, y, group = id, fill = fill),
                     colour = NA) +
        geom_path(data = polys[is.na(polys$fill), ], aes(x, y, group = id, colour = col),
                  linewidth = 0.5) +
        geom_text(data = KN, aes(x = x, y = y - (0.38 + (S$TXT_KN - 7) * 0.03), label = kn),
                  size = PT(S$TXT_KN),
                  family = FONT, colour = "grey25") +
        geom_text(data = marks, aes(x = x, y = gy, label = sym), size = PT(S$MARK_SIZE),
                  family = FONT, colour = S$MARK_COL) +
        # white label so the dotted outcome divider does not run through the text
        do.call(geom_label, c(list(data = M4ONLY, mapping = aes(x = x, y = y), label = "M4 only",
                                   size = PT(S$TXT_LAB), family = FONT, colour = "grey40",
                                   fontface = "italic", fill = "white",
                                   label.padding = unit(0.08, "lines")),
                              if (utils::packageVersion("ggplot2") >= "3.5.0") list(linewidth = 0)
                              else list(label.size = 0))) +
        scale_fill_identity() + scale_colour_identity() +
        annotate("text", x = 0.1, y = GRID_TOP - 0.35,
                 label = if (is.null(S$panel_label)) "(B)" else S$panel_label,
                 hjust = 0, fontface = "bold",
                 size = PT(S$TXT_PANEL), family = FONT)
    # Panel title after '(B)', aligned with 'Machinery' in panel A
    if (!is.null(S$panel_note))
        g <- g + annotate("text", x = 0.1 + (METRO$head / 72 * 2.3) * X_MAX / S$PAGE_W,
                          y = GRID_TOP - 0.35, label = S$panel_note, hjust = 0, fontface = "bold",
                          size = PT(METRO$head - 1), family = FONT)
    plotted <- rbind(
        data.frame(figure = rep(S$name, nrow(P)), part = rep("grid", nrow(P)), number = P$number, test = P$test, arm = P$arm,
                   region_set = P$region_set, outcome = rep("", nrow(P)), verdict = P$verdict, k = P$k,
                   n = P$n, x = P$x, y = P$y, stringsAsFactors = FALSE),
        data.frame(figure = rep(S$name, nrow(H)), part = rep("lead", nrow(H)), number = H$number, test = H$test, arm = H$arm,
                   region_set = H$region_set, outcome = H$outcome, verdict = H$verdict, k = H$k,
                   n = H$n, x = H$x, y = H$y, stringsAsFactors = FALSE))
    list(g = g, X_MAX = X_MAX, Y_MIN = Y_MIN, GRID_TOP = GRID_TOP, plotted = plotted)
}

# == Flowchart (drawn with grid, units in inches) ============================
STAGE <- data.frame(
    stage = c("Data and selection", "Per-cell model", "Module cells", "Donor quantity", "Tests"),
    dark  = c("#2E5E9E", "#A23B5D", "#5E7322", "#33414E", "#9B2C2C"),
    light = c("#E3ECF7", "#F7E1E8", "#EDF1DA", "#E3E7EB", "#F8E2E2"),
    stringsAsFactors = FALSE)
FLOW_STAGE <- c("Data and selection", "Data and selection", "Per-cell model", "Per-cell model",
                "Per-cell model", "Module cells", "Module cells", "Donor quantity", "Tests", "Tests")
FLOW_FIG <- c("Fig 1", "Fig 1", "Fig 2", "Fig 3", "Fig 4", "Figs 5-7", "Fig 8", "Supp Fig 13",
              "Figs 9-11, S14", if (PRE) "Figs 12-16" else "S15-S24")
flow_height <- function(S) S$FB_TOP + 2 * S$FB_H + S$FB_ROWGAP + 0.08 -
    if (identical(S$flow_style, "metro")) METRO$trim else 0
draw_flow <- function(S) {
    W <- S$PAGE_W; FH <- flow_height(S)
    FB_W <- (W - 0.2 - 4 * S$FB_GAP) / 5
    cx <- 0.1 + (rep(1:5, 2) - 1) * (FB_W + S$FB_GAP) + FB_W / 2
    top <- FH - S$FB_TOP - (rep(1:2, each = 5) - 1) * (S$FB_H + S$FB_ROWGAP)
    pushViewport(viewport(xscale = c(0, W), yscale = c(0, FH)))
    nat <- function(v) unit(v, "native")
    grid.text("(A)", x = nat(0.1), y = nat(FH - 0.13), just = "left",
              gp = gpar(fontsize = 11, fontface = "bold", fontfamily = FONT))
    grid.text("Machinery", x = nat(0.45), y = nat(FH - 0.13), just = "left",
              gp = gpar(fontsize = S$FB_HEAD, fontface = "bold", fontfamily = FONT))
    arr <- arrow(length = unit(0.07, "in"), type = "closed")
    agp <- gpar(col = "grey35", fill = "grey35", lwd = 1.1)
    for (i in 1:10) {
        st <- STAGE[STAGE$stage == FLOW_STAGE[i], ]
        bot <- top[i] - S$FB_H
        grid.roundrect(x = nat(cx[i]), y = nat(top[i] - S$FB_H / 2), width = nat(FB_W),
                       height = nat(S$FB_H), r = unit(0.05, "in"),
                       gp = gpar(fill = st$light, col = st$dark, lwd = 1.2))
        grid.roundrect(x = nat(cx[i]), y = nat(top[i] - S$FB_BAND / 2), width = nat(FB_W),
                       height = nat(S$FB_BAND), r = unit(0.05, "in"),
                       gp = gpar(fill = st$dark, col = st$dark))
        grid.rect(x = nat(cx[i]), y = nat(top[i] - S$FB_BAND * 0.75), width = nat(FB_W),
                  height = nat(S$FB_BAND / 2), gp = gpar(fill = st$dark, col = st$dark))
        grid.text(FLOW_STAGE[i], x = nat(cx[i]), y = nat(top[i] - S$FB_BAND / 2),
                  gp = gpar(fontsize = S$FB_STAGE, fontface = "bold", col = "white",
                            fontfamily = FONT))
        body_mid <- (top[i] - S$FB_BAND + bot + 0.14) / 2
        dq <- FLOW_STAGE[i] == "Donor quantity"
        grid.text(S$flow_text[i], x = nat(cx[i]), y = nat(body_mid),
                  gp = gpar(fontsize = S$FB_TXT, fontfamily = FONT, lineheight = 0.95,
                            col = if (dq) "black" else "grey10",
                            fontface = if (dq) "bold" else "plain"))
        grid.text(FLOW_FIG[i], x = nat(cx[i] + FB_W / 2 - 0.05), y = nat(bot + 0.08),
                  just = "right", gp = gpar(fontsize = S$FB_TAG, col = "grey35", fontfamily = FONT))
        if (S$feeds[i] != "")
            grid.text(S$feeds[i], x = nat(cx[i] - FB_W / 2 + 0.05), y = nat(bot + 0.08),
                      just = "left", gp = gpar(fontsize = S$FB_TAG, fontface = "bold",
                                               col = st$dark, fontfamily = FONT))
        if (i %% 5 != 0) {
            ym <- top[i] - S$FB_H / 2
            grid.lines(x = nat(c(cx[i] + FB_W / 2 + 0.02, cx[i] + FB_W / 2 + S$FB_GAP - 0.02)),
                       y = nat(c(ym, ym)), arrow = arr, gp = agp)
        }
    }
    mid <- top[5] - S$FB_H - S$FB_ROWGAP / 2
    grid.lines(x = nat(c(cx[5], cx[5], cx[6], cx[6])),
               y = nat(c(top[5] - S$FB_H - 0.02, mid, mid, top[6] + 0.02)), arrow = arr, gp = agp)
    popViewport()
}

# == Metro panel A ===========================================================
METRO <- list(
    dark  = c("#2E5E9E", "#A23B5D", "#5E7322", "#33414E", "#6B5B95"),
    stage = c("Data and GOI/DFGs\nidentification",
              "SKCS-core, pre-IP & post-IP,\nand transform into copula space",
              "Module cells\nrecovery & validation",
              "Donor\nquantity",
              "Clinical tests, competitors,\nablation, Batteries"),
    step_stage = c(1, 1, 2, 2, 2, 3, 3, 4, 5, 5),
    line1 = c("SEA-AD immune", "M4 cells;", "GAM fit,", "Inflection", "Copula space,",
              "CMRF band,", "Validation:", "Module post-IP%", "Clinical", "Batteries"),
    line2 = c("snRNA-seq", "MT-CO3/DFGs", "SKCS-core", "point, IPRS", "pitcher plot",
              "RSk-fCV gate", "Pitcher, UMAP", "per donor", "tests", "1-5"),
    fig = c("Fig 1", "Fig 1", "Fig 2", "Fig 3", "Fig 4", "Figs 5-7", "Fig 8", "S13",
            "Figs 9-11, S14", if (PRE) "Figs 12-16" else "S15-S24"),
    # stage text 10 pt so the fifth stage text fits its box
    head = 12, lab = 10, figsz = 9.5, ring = "grey20", ring_lwd = 1.8, num = 9, tab = 10, r = 0.15, lwd = 5, yl = 1.00,
    # station and stage-box offsets; the panel is 0.12 in shorter (trim)
    tab_y = 0.33, tab_h = 0.46, pad = 0.42, pad1 = 0.55, trim = 0.12,
    # the fifth stage box runs to the arrow tip
    last_to_tip = TRUE)
draw_flow_metro <- function(S) {
    M <- METRO; W <- S$PAGE_W; FH <- flow_height(S)
    pushViewport(viewport(xscale = c(0, W), yscale = c(0, FH)))
    nat <- function(v) unit(v, "native")
    gp_t <- function(size, face = "plain", col = "grey10", lh = 1)
        gpar(fontsize = size, fontface = face, col = col, fontfamily = FONT, lineheight = lh)
    # '(A)' and 'Machinery' are added in the figure layout, not drawn here
    xs <- seq(0.65, W - 0.65, length.out = 10); yl <- M$yl; r <- M$r
    st <- M$step_stage; dk <- M$dark
    # Lines are drawn as dark-grey bands with the coloured line on top, so each
    # carries the same outline width as the station circles.
    xt <- W - 0.1; L <- 0.2; hh <- 0.1; xn <- xt - L * 0.62
    seg_x0 <- c(xs[1:9], xs[10]); seg_x1 <- c(xs[2:10], xn + 0.01)
    seg_col <- c(dk[st[2:10]], dk[5])
    for (j in seq_along(seg_x0))
        grid.lines(x = nat(c(seg_x0[j], seg_x1[j])), y = nat(c(yl, yl)),
                   gp = gpar(col = M$ring, lwd = M$lwd + 2 * M$ring_lwd, lineend = "butt"))
    for (j in seq_along(seg_x0))
        grid.lines(x = nat(c(seg_x0[j], seg_x1[j])), y = nat(c(yl, yl)),
                   gp = gpar(col = seg_col[j], lwd = M$lwd, lineend = "butt"))
    # notched (barbed) arrowhead: tip, two back corners, notch
    grid.polygon(x = nat(c(xt, xt - L, xn, xt - L)), y = nat(c(yl, yl + hh, yl, yl - hh)),
                 gp = gpar(fill = dk[5], col = M$ring, lwd = M$ring_lwd, linejoin = "mitre"))
    lh <- M$lab / 72 * 1.18
    for (i in 1:10) {
        cx <- xs[i]
        grid.circle(x = nat(cx), y = nat(yl), r = unit(r, "in"),
                    gp = gpar(fill = dk[st[i]], col = M$ring, lwd = M$ring_lwd))
        grid.text(as.character(as.roman(i)), x = nat(cx), y = nat(yl),
                  gp = gp_t(M$num, "bold", "white"))
        y2 <- yl + r + 0.07 + lh * 0.5
        grid.text(M$line1[i], x = nat(cx), y = nat(y2 + lh), gp = gp_t(M$lab, "bold"))
        grid.text(M$line2[i], x = nat(cx), y = nat(y2), gp = gp_t(M$lab, "bold"))
        grid.text(M$fig[i], x = nat(cx), y = nat(yl - r - 0.11),
                  gp = gp_t(M$figsz, "bold", dk[st[i]]))
    }
    for (k in 1:5) {
        idx <- which(st == k)
        pad <- if (length(idx) == 1) M$pad1 else M$pad
        xa <- xs[min(idx)] - pad; xb <- xs[max(idx)] + pad
        if (k == 5 && isTRUE(M$last_to_tip)) xb <- xt
        grid.roundrect(x = nat((xa + xb) / 2), y = nat(M$tab_y), width = nat(xb - xa),
                       height = nat(M$tab_h), r = unit(0.05, "in"),
                       gp = gpar(fill = dk[k], col = M$ring, lwd = M$ring_lwd))
        grid.text(M$stage[k], x = nat((xa + xb) / 2), y = nat(M$tab_y),
                  gp = gp_t(M$tab, "bold", "white", 1.05))
    }
    popViewport()
}

# == Pages ===================================================================
out_path <- function(base, ext) {
    p <- paste0(base, ext)
    locked <- file.exists(p) && !isTRUE(tryCatch({
        con <- file(p, open = "ab"); close(con); TRUE }, error = function(e) FALSE,
        warning = function(w) FALSE))
    if (locked) {
        p <- paste0(base, "_", format(Sys.time(), "%Y%m%d-%H%M%S"), ext)
        cat(sprintf("  target open in another program; writing %s\n", basename(p)))
    }
    p
}
# Main figure: NeoC and NonNeoC combined.
SPEC_BOTH <- SPEC_PATHOLOGY
SPEC_BOTH$name <- "landscape_pathology_BothRegions"
SPEC_BOTH$base <- "Fig17_Summary_landscape_pathology_BothRegions"
SPEC_BOTH$leads <- data.frame(arm = ARMS, region_set = "Both", stringsAsFactors = FALSE)
SPEC_BOTH$TALLIES <- c("tally_both", "tally_pct")
# No panel label or panel note; both are added in the figure layout
SPEC_BOTH$panel_label <- ""
SPEC_BOTH$fig_override <- c(SPEC_PATHOLOGY$fig_override, T12 = paste0(BP("15A"), ", S18"))
SPEC_BOTH$BOX_MIN <- 75
SPEC_BOTH$BOX_COL <- "#002060"              # header navy
SPEC_BOTH$LEGEND_LINE3 <- sprintf(paste0("Tallies count the same tests for every arm: 1-13 without %s ",
    "(not assessable for PVM_3 and PVM_1). Box: mean %% of 75 or more."),
    paste(TESTDEF$number[match(TALLY_EXCLUDE, TESTDEF$test)], collapse = " and "))
# Strips start further right so "S23A, S18" fits before the M4 strip
SPEC_BOTH$LEAD_X0 <- 18.0
SPEC_BOTH$X_FIG <- 15.3
SPEC_PATHOLOGY$no_flow <- TRUE
SPEC_PATHOLOGY$LEGEND_LINE3 <- sprintf(paste0("Tallies count the same tests for every arm: 1-13 without %s ",
    "(not assessable for PVM_3 and PVM_1)."),
    paste(TESTDEF$number[match(TALLY_EXCLUDE, TESTDEF$test)], collapse = " and "))
SPEC_PATHOLOGY$panel_label <- ""            # numbered when placed as a supplementary figure
SPEC_PATH_NONNEOC <- SPEC_PATHOLOGY
SPEC_PATH_NONNEOC$name <- "landscape_pathology_NonNeoC"
SPEC_PATH_NONNEOC$base <- "Fig17_Summary_landscape_pathology_NonNeoC"
SPEC_PATH_NONNEOC$leads <- data.frame(arm = ARMS, region_set = "NonNeoC", stringsAsFactors = FALSE)
SPEC_PATH_NONNEOC$no_flow <- TRUE           # panel B only
SPEC_PATH_NONNEOC$panel_label <- ""         # numbered when placed as a supplementary figure
# Supplementary Figure S23A shows NeoC; its NonNeoC panels are Supplementary Figure S18
SPEC_PATH_NONNEOC$fig_override <- c(SPEC_PATHOLOGY$fig_override, T12 = "S18")
# Figures written by this run; the other specs remain defined and can be added here.
RENDER_SPECS <- list(SPEC_BOTH)
PLOTTED <- list()
BOXED <- list()
for (S in RENDER_SPECS) {
    G <- build_grid(S)
    PLOTTED[[S$name]] <- G$plotted
    TITLE_H <- if (isFALSE(S$show_title)) 0.06 else 0.34
    FH <- if (isTRUE(S$no_flow)) 0 else flow_height(S)
    grid_h <- S$PAGE_W * (G$GRID_TOP - G$Y_MIN) / G$X_MAX
    PAGE_H <- TITLE_H + FH + grid_h + 0.05
    draw_page <- function() {
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(3, 1,
            heights = unit(c(TITLE_H, FH, grid_h), "in"), widths = unit(S$PAGE_W, "in"))))
        if (!isFALSE(S$show_title))
            grid.text(FIG_NO, x = unit(0.1, "in"), just = "left",
                      gp = gpar(fontsize = 12, fontface = "bold", fontfamily = FONT),
                      vp = viewport(layout.pos.row = 1))
        if (!isTRUE(S$no_flow)) {
            pushViewport(viewport(layout.pos.row = 2))
            if (identical(S$flow_style, "metro")) draw_flow_metro(S) else draw_flow(S)
            popViewport()
        }
        print(G$g, vp = viewport(layout.pos.row = 3))
        popViewport()
    }
    base <- file.path(OUT_DIR, S$base)
    render <- function(open_fn, ext) {
        p <- out_path(base, ext)
        ok <- tryCatch({ open_fn(p); draw_page(); dev.off(); TRUE },
                       error = function(e) { if (dev.cur() > 1) dev.off()
                           cat(sprintf("  %s not written: %s\n", ext, conditionMessage(e)))
                           FALSE })
        if (ok) cat(sprintf("  wrote %s (%.2f x %.2f in)\n", basename(p), S$PAGE_W, PAGE_H))
    }
    render(function(p) grDevices::cairo_pdf(p, width = S$PAGE_W, height = PAGE_H), ".pdf")
    render(function(p) grDevices::png(p, width = S$PAGE_W, height = PAGE_H, units = "in",
                                      res = 300, type = "cairo"), ".png")
    render(function(p) {
        if (requireNamespace("svglite", quietly = TRUE))
            svglite::svglite(p, width = S$PAGE_W, height = PAGE_H)
        else grDevices::svg(p, width = S$PAGE_W, height = PAGE_H)
    }, ".svg")
}

# Tests 16-21 (IDs T14-T19) for the landscape legend (not drawn in that figure)
cat("\n  Tests 16-21 (per arm, NeoC + NonNeoC; for the landscape legend):\n")
for (r in TESTDEF$test[TESTDEF$block == 2]) {
    x <- CELLS[CELLS$test == r, ]
    cat(sprintf("  %2d %-42s %s\n", x$number[1], x$label[1],
                paste(sprintf("%s %s%s", x$arm, x$verdict,
                              ifelse(x$n > 0, paste0(" (", x$k, "/", x$n, ")"), "")),
                      collapse = "; ")))
}

# == 6. Tables ===============================================================
cat("\n== 6. Tables\n")
CELLS_OUT <- CELLS[, c("number", "test", "label", "battery", "fig", "block", "arm",
                       "region_set", "outcomes_applicable", "k", "n", "unreadable",
                       "untestable", "verdict", "source", "rule")]
MODELS_OUT <- MODELS
MODELS_OUT$number <- TESTDEF$number[match(MODELS_OUT$test, TESTDEF$test)]
MODELS_OUT <- MODELS_OUT[, c("number", "test", "arm", "region_set", "outcome", "variant",
                             "estimate", "fdr", "status", "applicable", "source")]
PLOT_OUT <- do.call(rbind, PLOTTED)
write.csv(CELLS_OUT,  file.path(OUT_DIR, "Fig17_cells.csv"), row.names = FALSE)
write.csv(MODELS_OUT, file.path(OUT_DIR, "Fig17_models.csv"), row.names = FALSE)
write.csv(PLOT_OUT,   file.path(OUT_DIR, "Fig17_plotted.csv"), row.names = FALSE)
HEAD_OUT <- HEAD[, c("number", "test", "arm", "region_set", "outcome", "k", "n", "unreadable",
                     "verdict", "regions_read")]
HEAD_OUT <- HEAD_OUT[order(match(HEAD_OUT$arm, ARMS), HEAD_OUT$outcome, HEAD_OUT$number), ]
write.csv(HEAD_OUT, file.path(OUT_DIR, "Fig17_lead_result.csv"), row.names = FALSE)
SHEETS <- list(Tests = TESTDEF, Cells = CELLS_OUT, Lead_result = HEAD_OUT, Models = MODELS_OUT,
               Tallies = TALLY[, c("arm", "region_set", "tests_held", "tests_partial",
                                   "tests_lost", "tests_read")],
               Lead_tally = HTALLY[, c("arm", "region_set", "outcome", "tests_held",
                                       "tests_partial", "tests_lost", "tests_read",
                                       "mean_pct_per_test", "pooled_k", "pooled_n",
                                       "pooled_pct")],
               Both_region_tally = BTALLY[, c("arm", "outcome", "both_held", "both_read")],
               Boxed = do.call(rbind, c(list(data.frame(figure = character(0), arm = character(0),
                                                        region_set = character(0), outcome = character(0),
                                                        label_pct = character(0))),
                                         Map(function(n, b) if (nrow(b)) cbind(figure = n, b) else NULL,
                                             names(BOXED), BOXED))),
               Checks = CHECKS)
xp <- out_path(file.path(OUT_DIR, "Fig17_Summary"), ".xlsx")
if (requireNamespace("openxlsx", quietly = TRUE)) {
    wb <- openxlsx::createWorkbook()
    for (s in names(SHEETS)) {
        openxlsx::addWorksheet(wb, s)
        openxlsx::writeData(wb, s, SHEETS[[s]],
                            headerStyle = openxlsx::createStyle(textDecoration = "bold"))
        openxlsx::freezePane(wb, s, firstRow = TRUE)
    }
    openxlsx::saveWorkbook(wb, xp, overwrite = TRUE)
    cat(sprintf("  wrote %s\n", basename(xp)))
} else if (requireNamespace("writexl", quietly = TRUE)) {
    writexl::write_xlsx(SHEETS, xp)
    cat(sprintf("  wrote %s (writexl; no frozen pane)\n", basename(xp)))
} else cat("  xlsx not written: install openxlsx or writexl; CSV files written\n")
if (file.exists(file.path(OUT_DIR, "Fig17_tests.csv")))
    cat("  note: Fig17_tests.csv is superseded by Fig17_models.csv\n")

cat(sprintf("\n== Summary: %d PASS, %d WARN, %d FAIL\n", sum(CHECKS$verdict == "PASS"),
            sum(CHECKS$verdict == "WARN"), sum(CHECKS$verdict == "FAIL")))
cat(sprintf("VERDICT: %s\n", if (any(CHECKS$verdict == "FAIL")) "FAIL" else "PASS"))
cat(sprintf("Log: %s\n", LOG))

}, error = function(e) {
    cat(sprintf("\nERROR: %s\n", conditionMessage(e)))
}, finally = on.exit_sink())
