##############################################################################
# Part13c_CP10K_M4_SuppConsolidate_code.r
#
# CEP-IP SEA-AD - Part 13c: two supplementary figures and one supplementary
# table from the finished Part 13 and Part 13b runs. Nothing is refitted: every
# value is read from the two SuppTable workbooks and the Part 13b draws file.
#
#   Supplementary Figure 1 (Part 13, M4 share and Module state against %MT)
#     (A) Module post-IP% beta without and with the M4 share (cov4, +BTC, cov4 + %MT)
#     (B) Module state and %MT of the same M4 cells in one model (Fig 14D, fraction)
#     (C) Verdict counts: Fig 14D, Fig 14D + M4 share, fraction + %MT
#     (D) Within-M4 fraction against %MT per donor, shaded by MMSE
#   Supplementary Figure 2 (Part 13b, M4 cells against random PVM_3 cells)
#     (A) Module post-IP% beta on cognition against 10,000 random subsets
#     (B) r with %MT and the mean fraction against the random subsets
#     (C) Figure 14D verdicts of the random subsets and of the M4 cells
#   Supplementary Table: one sheet per figure panel (the plotted values), then
#     the complete Part 13 and Part 13b results. No guardrail or log sheets.
#   In the preprint, Supplementary Figures 1 and 2 are Supplementary Figures S20
#   and S21, and finalize_supp_tables.py turns the table into Supplementary
#   Table S18.
#
# Usage (from the SEA-AD folder, or set SEAAD_ROOT), after Part 13 and Part 13b
#   source("Part13c_CP10K_M4_SuppConsolidate_code.r")
#
# INPUTS (found under CEP_HOME)
#   Part13_M4Share/Part13_CP10K_M4Share_SuppTable.xlsx
#   Part13b_M4RandomSubset/Part13b_CP10K_M4RandomSubset_SuppTable.xlsx
#   Part13b_M4RandomSubset/Part13b_CP10K_M4RandomSubset_null_draws.csv
#
# OUTPUTS (Part13b_M4RandomSubset)
#   Part13_13b_CP10K_M4_SuppFig1_M4Share.pdf / .png        180 x 270 mm
#   Part13_13b_CP10K_M4_SuppFig2_RandomSubset.pdf / .png   180 x 200 mm
#   Part13_13b_CP10K_M4_SuppTable.xlsx                     plain sheets
#   Part13_13b_CP10K_M4_legend_draft.txt
#   Part13_13b_CP10K_M4_Guardrail_Log.txt, _console_log_<stamp>.txt
##############################################################################

cfg <- function(name, default) {
    o <- if (exists("P13C_CFG", envir = globalenv(), inherits = FALSE)) get("P13C_CFG", envir = globalenv()) else list()
    if (name %in% names(o)) o[[name]] else default
}
CEP_HOME   <- cfg("CEP_HOME", Sys.getenv("SEAAD_ROOT", getwd()))
OUT_DIR    <- file.path(CEP_HOME, "Part13b_M4RandomSubset")
P13_XLSX   <- cfg("P13_XLSX", file.path(CEP_HOME, "Part13_M4Share", "Part13_CP10K_M4Share_SuppTable.xlsx"))
P13B_XLSX  <- cfg("P13B_XLSX", file.path(OUT_DIR, "Part13b_CP10K_M4RandomSubset_SuppTable.xlsx"))
P13B_DRAWS <- cfg("P13B_DRAWS", file.path(OUT_DIR, "Part13b_CP10K_M4RandomSubset_null_draws.csv"))
PREFIX     <- "Part13_13b_CP10K_M4"
GUARD_STRICT <- cfg("GUARD_STRICT", TRUE)
VIF_MAX <- 5; ALPHA <- 0.05
DATASETS <- c("NeoC", "NonNeoC")
COG_LAB  <- c("Cognitive status", "MMSE score", "CASI score")
SHORT    <- c("Cognitive status" = "CogS", "MMSE score" = "MMSE", "CASI score" = "CASI")

suppressPackageStartupMessages({ library(openxlsx); library(ggplot2); library(patchwork) })
if (!dir.exists(CEP_HOME)) stop("CEP_HOME does not exist: ", CEP_HOME, call. = FALSE)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
f_out <- function(x) file.path(OUT_DIR, paste0(PREFIX, x))
OUT_F1 <- f_out("_SuppFig1_M4Share"); OUT_F2 <- f_out("_SuppFig2_RandomSubset")
OUT_XLSX <- f_out("_SuppTable.xlsx"); OUT_LEG <- f_out("_legend_draft.txt"); LEDGER <- f_out("_Guardrail_Log.txt")
CONSOLE <- f_out(sprintf("_console_log_%s.txt", format(Sys.time(), "%Y%m%d_%H%M%S")))
.stale <- c(paste0(OUT_F1, c(".pdf", ".png")), paste0(OUT_F2, c(".pdf", ".png")), OUT_XLSX, OUT_LEG, LEDGER)
invisible(file.remove(.stale[file.exists(.stale)]))
.con <- file(CONSOLE, open = "wt"); sink(.con, split = TRUE)
.old_width <- options(width = 200)
cat(strrep("=", 78), "\n  Part 13c - consolidated supplementary figures and table\n  ", format(Sys.time()),
    "\n  Root : ", CEP_HOME, "\n  Write: ", OUT_DIR, "\n", strrep("=", 78), "\n", sep = "")

.gr <- list()
gr <- function(id, status, detail) {
    .gr[[length(.gr) + 1L]] <<- data.frame(id = id, status = status, detail = detail, stringsAsFactors = FALSE)
    cat(sprintf("  [%-4s]  %-10s %s\n", status, id, detail))
    if (status == "FAIL" && GUARD_STRICT) stop(sprintf("%s FAILED: %s", id, detail), call. = FALSE)
}
num <- function(x) suppressWarnings(as.numeric(x))
fmt <- function(x, d = 3) ifelse(is.na(x), "NA", formatC(x, format = "f", digits = d))
# FDR labels: 3 decimals from 0.001; below 0.001 as m x 10 with a superscript
# exponent (manuscript style), never e notation
# (a session that is not UTF-8 falls back to "^-4")
sup_digits <- function(z) if (isTRUE(l10n_info()[["UTF-8"]]))
    chartr("-0123456789", "\u207b\u2070\u00b9\u00b2\u00b3\u2074\u2075\u2076\u2077\u2078\u2079", z) else paste0("^", z)
fmt_p <- function(x) {
    e <- ifelse(is.na(x) | x <= 0, 0, floor(log10(x)))
    small <- sprintf("%.2fx10%s", x / 10^e, sup_digits(as.character(e)))
    ifelse(is.na(x), "NA", ifelse(x < 0.001, small, formatC(x, format = "f", digits = 3)))
}
short_cog <- function(lab) ifelse(lab == "Cognitive status", "CogS\u2020", SHORT[lab])

# ---- palettes and theme (as Part 13 and Part 13b) ----------------------------
COL <- list(original = "#3B3B3B", variant = "#6A3D9A", strip = "#262626", grid = "#EBEBEB",
            m4 = "#1B9E77", pvm3 = "#F2A950", own = "#3B3B3B",
            hist_fill = "#F9D8A8", hist_border = "#B9762A")   # histogram style V2, chosen 2026-10-10
HIST_BINS <- 45L
VERDICT_FILL <- c("Module only" = "#2D2D2D", "Second only" = "#1275D8", "Both" = "#7B7B7B",
                  "Neither" = "#D9D9D9", "Unreadable" = "#FBE7A2")
# Font sizes (pt). Axis tick numbers stay at 6 pt; every other text was enlarged
# on 2026-10-10 for legibility (TXT), including the category labels on the y-axis
# of panels 1A-1C and 2C (TXT$ycat). GEOM_PT converts pt to the mm size of geom_text.
TXT <- list(base = 9, strip = 8.5, title = 9.5, axis_title = 8.5, ticks = 6, ycat = 8,
            legend_title = 8.5, legend_text = 8, tag = 11, in_panel = 8)
GEOM_PT <- function(pt) pt / ggplot2::.pt
th <- theme_bw(base_size = TXT$base) +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = COL$grid, linewidth = 0.25),
          strip.background = element_rect(fill = COL$strip, colour = COL$strip),
          strip.text = element_text(colour = "white", face = "bold", size = TXT$strip),
          plot.title = element_text(face = "bold", size = TXT$title, hjust = 0),
          axis.title = element_text(size = TXT$axis_title), axis.text = element_text(size = TXT$ticks, colour = "black"),
          legend.position = "bottom", legend.title = element_text(size = TXT$legend_title, face = "bold"),
          legend.text = element_text(size = TXT$legend_text), legend.key.size = unit(3.4, "mm"),
          legend.margin = margin(0, 0, 0, 0), plot.tag = element_text(face = "bold", size = TXT$tag))
# category labels on a y-axis (text, not numbers)
th_ycat <- theme(axis.text.y = element_text(size = TXT$ycat, colour = "black"))
# geom_label border: `linewidth` from ggplot2 3.5.0, `label.size` before.
label_noborder <- if (utils::packageVersion("ggplot2") >= "3.5.0") list(linewidth = 0) else list(label.size = 0)

run_13c <- function() {
    # ---- STEP 1: inputs ---------------------------------------------------------
    cat("\nSTEP 1 -- inputs\n")
    for (f in c(P13_XLSX, P13B_XLSX, P13B_DRAWS)) if (!file.exists(f)) stop("Input not found: ", f, call. = FALSE)
    for (f in c(P13_XLSX, P13B_XLSX, P13B_DRAWS)) cat(sprintf("  %s | %.1f MB\n", f, file.size(f) / 2^20))
    rd <- function(f, s) { x <- openxlsx::read.xlsx(f, sheet = s); if (is.null(x)) stop(basename(f), ": sheet ", s, " is empty.", call. = FALSE); x }
    need13 <- c("A_Specs", "B_Frame", "C_Models", "D_Verdicts", "E_Decomposition", "F_Bootstrap", "G_Collinearity",
                "H_Readings", "L_PostHoc")
    need13b <- c("A_Design", "B_Frame", "D_Summary", "E_Verdicts_MB", "F_Label_Agreement", "H_Readings")
    miss <- c(setdiff(need13, getSheetNames(P13_XLSX)), setdiff(need13b, getSheetNames(P13B_XLSX)))
    if (length(miss)) stop("Sheets missing: ", paste(miss, collapse = ", "),
                           ". The Part 13 table must come from the run with the post-hoc models (sheet L_PostHoc).", call. = FALSE)
    T13 <- lapply(stats::setNames(need13, need13), function(s) rd(P13_XLSX, s))
    T13b <- lapply(stats::setNames(need13b, need13b), function(s) rd(P13B_XLSX, s))
    DR <- utils::read.csv(P13B_DRAWS, check.names = FALSE, stringsAsFactors = FALSE)
    gr("GR-13c-01", "PASS", sprintf("inputs read: Part 13 table %d sheets, Part 13b table %d sheets, %d draws",
                                    length(T13), length(T13b), nrow(DR)))
    M <- T13$C_Models; V <- T13$D_Verdicts
    for (k in c("n", "beta", "ci_lo", "ci_hi", "fdr", "vif_full", "fdr_second", "vif_max", "p_value"))
        for (nm in c("M", "V")) { z <- get(nm); if (k %in% names(z)) { z[[k]] <- num(z[[k]]); assign(nm, z) } }

    # ---- STEP 2: the same donors and values in both runs --------------------------
    cat("\nSTEP 2 -- Part 13 and Part 13b describe the same donors\n")
    F13 <- T13$B_Frame; F13b <- T13b$B_Frame
    k13 <- paste(F13$dataset, F13$donor_id); k13b <- paste(F13b$dataset, F13b$donor_id)
    j <- match(k13, k13b)
    same_set <- setequal(k13, k13b) && !anyNA(j)
    dv <- function(a, b) { a <- num(a); b <- num(b); if (!identical(is.na(a), is.na(b))) return(Inf)
        if (all(is.na(a))) 0 else max(abs(a - b), na.rm = TRUE) }
    d_share <- if (same_set) dv(F13$share, F13b$m4_share[j]) else Inf
    d_mt    <- if (same_set) dv(F13$mt, F13b$mt_m4[j]) else Inf
    d_frac  <- if (same_set) dv(F13$frac, F13b$frac_m4_own_labels[j]) else Inf
    d_post  <- if (same_set) dv(F13$pct_post, F13b$x_m4_own_labels[j]) else Inf
    gr("GR-13c-02", if (max(d_share, d_mt, d_frac, d_post) <= 1e-9) "PASS" else "FAIL",
       sprintf(paste("same %d donor rows in both runs (%s); M4 share, %%MT, within-M4 fraction and Module post-IP%% agree",
                     "(max |diff| %s, %s, %s, %s)"), nrow(F13), same_set, format(d_share), format(d_mt), format(d_frac),
               format(d_post)))

    # ---- STEP 3: Supplementary Figure 1 data ----------------------------------------
    cat("\nSTEP 3 -- Supplementary Figure 1 (Part 13)\n")
    vmax <- stats::aggregate(list(vif_max = M$vif_full), by = list(spec = M$spec, dataset = M$dataset, outcome = M$outcome),
                             FUN = function(v) if (anyNA(v)) Inf else max(v))
    # (A) Module post-IP% without and with the M4 share
    setA <- data.frame(spec = c("S04", "S06", "S05", "S07", "S09", "S10"),
                       set = rep(c("cov4", "+BTC", "cov4 + %MT"), each = 2),
                       series = rep(c("Without M4 share", "With M4 share"), 3), stringsAsFactors = FALSE)
    A <- merge(M[M$predictor == "module" & M$outcome %in% COG_LAB & M$spec %in% setA$spec, ], setA, by = "spec")
    A <- merge(A, vmax, by = c("spec", "dataset", "outcome"), all.x = TRUE)
    # (B) Module state and %MT in one model
    setB <- data.frame(spec = c("S09", "S15", "S16"), set = c("Module post-IP% + %MT (Fig 14D)", "fraction + %MT, cov4",
                       "fraction + %MT, +BTC"), state = c("module", "frac", "frac"), stringsAsFactors = FALSE)
    B <- merge(M[M$spec %in% setB$spec & M$outcome %in% COG_LAB, ], setB, by = "spec")
    B <- B[B$predictor == B$state | B$predictor == "mt", ]
    B$series <- ifelse(B$predictor == "mt", "%MT of the M4 cells", "Module state (Module post-IP% or within-M4 fraction)")
    B <- merge(B, vmax, by = c("spec", "dataset", "outcome"), all.x = TRUE)
    # (C) verdicts
    C <- V[V$spec %in% c("S09", "S10", "S15") & V$group %in% c("Cognition", "Pathology", "APOE"), ]
    C$row <- paste0(C$dataset, c(S09 = ", Module (Fig 14D)", S10 = ", Module + M4 share", S15 = ", fraction")[C$spec])
    # (D) per donor
    mm <- if ("MMSE_num" %in% names(F13)) "MMSE_num" else grep("^MMSE", names(F13), value = TRUE)[1]
    D <- data.frame(dataset = F13$dataset, donor_id = F13$donor_id, frac = num(F13$frac), mt = num(F13$mt),
                    mmse = num(F13[[mm]]), stringsAsFactors = FALSE)
    D <- D[is.finite(D$frac) & is.finite(D$mt), ]
    nA <- nrow(A); nB <- nrow(B); nC <- nrow(C)
    gr("GR-13c-03", if (nA == 36L && nB == 36L && nC == 54L && nrow(D) > 0 && !anyNA(A$beta) && !anyNA(B$beta)) "PASS" else "FAIL",
       sprintf(paste("Figure 1 rows from the Part 13 table: panel A %d (expected 36), panel B %d (36), panel C %d (54 =",
                     "3 models x 2 region sets x 9 outcomes), panel D %d donors; MMSE column '%s'"), nA, nB, nC, nrow(D), mm))

    # ---- STEP 4: Supplementary Figure 2 data ----------------------------------------
    cat("\nSTEP 4 -- Supplementary Figure 2 (Part 13b)\n")
    S <- T13b$D_Summary
    for (k in setdiff(names(S), c("dataset", "readout", "outcome", "bh_family"))) S[[k]] <- num(S[[k]])
    roA <- paste0("bA_", SHORT); roB <- c("r_x_mt", "r_frac_mt", "mean_frac")
    DR$valid <- as.logical(DR$valid)
    chk <- do.call(rbind, lapply(DATASETS, function(ds) do.call(rbind, lapply(c(roA, roB), function(ro) {
        v <- DR[DR$dataset == ds & DR$valid, ro]; v <- v[is.finite(v)]
        s <- S[S$dataset == ds & S$readout == ro, ]
        data.frame(ds = ds, ro = ro, d = max(abs(c(stats::median(v) - s$null_median,
                   unname(stats::quantile(v, .025)) - s$null_q025, unname(stats::quantile(v, .975)) - s$null_q975))),
                   n = length(v), n_tab = s$draws_used)
    }))))
    gr("GR-13c-04", if (all(chk$d <= 1e-9) && all(chk$n == chk$n_tab)) "PASS" else "FAIL",
       sprintf(paste("the draws file reproduces the Part 13b summary for the 12 plotted readouts (median, 2.5 and 97.5",
                     "percentiles, max |diff| %s; draws used %s)"), format(max(chk$d)), paste(unique(chk$n), collapse = "/")))
    VB <- T13b$E_Verdicts_MB
    for (k in grep("^null_", names(VB), value = TRUE)) VB[[k]] <- num(VB[[k]])
    vsum <- rowSums(VB[, grep("^null_", names(VB))])
    gr("GR-13c-05", if (nrow(VB) == 6L && all(abs(vsum - 1) < 1e-9)) "PASS" else "FAIL",
       sprintf("Figure 2C: %d rows (expected 6), verdict shares of the random subsets sum to 1 (max |diff| %s)",
               nrow(VB), format(max(abs(vsum - 1)))))

    # ---- STEP 5: figures --------------------------------------------------------------
    cat("\nSTEP 5 -- figures\n")
    ok1 <- tryCatch({ fig1(A, B, C, D); TRUE }, error = function(e) { cat("  FIGURE 1 FAILED: ", conditionMessage(e), "\n"); FALSE })
    ok2 <- tryCatch({ fig2(S, DR, VB); TRUE }, error = function(e) { cat("  FIGURE 2 FAILED: ", conditionMessage(e), "\n"); FALSE })
    .lc <- file(OUT_LEG, open = "w", encoding = "UTF-8"); writeLines(c(LEG1, "", LEG2), .lc); close(.lc)

    # ---- STEP 6: consolidated table -----------------------------------------------------
    cat("\nSTEP 6 -- supplementary table\n")
    keepA <- c("set", "series", "spec", "dataset", "outcome", "n", "beta", "ci_lo", "ci_hi", "fdr", "vif_max")
    keepB <- c("set", "series", "spec", "dataset", "outcome", "predictor", "n", "beta", "ci_lo", "ci_hi", "fdr", "vif_max")
    keepC <- c("row", "spec", "dataset", "outcome", "group", "predictor", "n", "beta", "fdr", "beta_second", "fdr_second",
               "vif_max", "verdict")
    C_out <- C[, intersect(keepC, names(C))]
    C_out$verdict <- sub("^Second only$", "%MT only", sub("^Module only$", "State only", C_out$verdict))
    r_tab <- do.call(rbind, lapply(DATASETS, function(ds) { z <- D[D$dataset == ds, ]
        data.frame(dataset = ds, n = nrow(z), r_fraction_mt = stats::cor(z$frac, z$mt)) }))
    SA <- S[S$readout %in% roA, ]; SB <- S[S$readout %in% roB, ]
    cols_s <- c("dataset", "readout", "outcome", "observed_M4_PVM3labels", "M4_own_labels", "null_median", "null_q025",
                "null_q975", "p_two", "bh_family", "fdr")
    frame <- merge(F13, F13b[, setdiff(names(F13b), c("n_total"))], by = c("donor_id", "dataset"), all = TRUE)
    readings <- rbind(data.frame(part = "Part 13", T13$H_Readings, stringsAsFactors = FALSE),
                      data.frame(part = "Part 13b", T13b$H_Readings, stringsAsFactors = FALSE))
    README <- data.frame(sheet = c("Fig1A", "Fig1B", "Fig1C", "Fig1D", "Fig2A", "Fig2B", "Fig2C",
                                   "T1_Models_13", "T2_Verdicts_13", "T3_Bootstrap_13", "T4_Decomposition_13",
                                   "T5_Collinearity_13", "T6_PostHoc_13", "T7_Specs_13", "T8_RandomSubset_13b",
                                   "T9_LabelAgreement_13b", "T10_Design_13b", "T11_Donors", "T12_Readings"),
        content = c("Supplementary Figure 1A: standardized beta of Module post-IP% on cognition without and with the M4 share (cov4, +BTC, cov4 + %MT), 95% CI, FDR in the family of the figure the model extends, largest full-model VIF.",
                    "Supplementary Figure 1B: standardized beta of the Module state and of %MT in one model (Fig 14D model; within-M4 fraction at cov4 and +BTC).",
                    "Supplementary Figure 1C: verdicts at FDR < 0.05 for the Module state and %MT in one model, all nine clinical outcomes (Fig 14D; Fig 14D + M4 share; fraction + %MT). Unreadable = largest VIF >= 5.",
                    "Supplementary Figure 1D: within-M4 fraction, %MT of the M4 cells and MMSE per donor (donors with at least 5 M4 cells), with Pearson r per region set.",
                    "Supplementary Figure 2A: Module post-IP% beta on cognition, M4 cells under PVM_3 labels, M4 cells under their own labels (reference), and the random-subset median and 95% range; two-sided empirical p and BH FDR.",
                    "Supplementary Figure 2B: r of Module post-IP% with %MT, r of the fraction with %MT and the mean fraction, as in 2A.",
                    "Supplementary Figure 2C: Figure 14D verdicts of the M4 cells and the share of random subsets with each verdict.",
                    "Part 13, every fitted model, one row per predictor (n, beta, SE, p, 95% CI, FDR, dR2, VIF).",
                    "Part 13, every verdict (Module against M4 share, CP or %MT; fraction against %MT).",
                    "Part 13, donor bootstrap of the change in the Module coefficient when the M4 share is added.",
                    "Part 13, Figure 10C/10D decomposition models with and without the M4 share.",
                    "Part 13, Pearson and Spearman r between the M4 share, the Module state and %MT.",
                    "Part 13, post-hoc models of the within-M4 fraction and %MT, cognition.",
                    "Part 13, every model specification and its BH family.",
                    "Part 13b, every readout against the random subsets (empirical p, BH FDR and family).",
                    "Part 13b, Module post-IP labels of the M4 cells under the M4 arm's own labels and PVM_3 labels.",
                    "Part 13b, pool, labels, strata, draws, seed and models.",
                    "Per donor and region set, Part 13 and Part 13b donor values in one row.",
                    "Pre-set and post-hoc readings of both parts and their results."),
        stringsAsFactors = FALSE)
    wb <- list(README = README,
               Fig1A = A[order(A$dataset, A$outcome, A$set, A$series), intersect(keepA, names(A))],
               Fig1B = B[order(B$dataset, B$outcome, B$spec, B$predictor), intersect(keepB, names(B))],
               Fig1C = C_out[order(C_out$dataset, C_out$spec, C_out$outcome), ],
               Fig1D = rbind(transform(D, r_fraction_mt = NA_real_),
                             data.frame(dataset = r_tab$dataset, donor_id = "Pearson r (all donors in the panel)", frac = NA,
                                        mt = NA, mmse = NA, r_fraction_mt = r_tab$r_fraction_mt)),
               Fig2A = SA[, cols_s], Fig2B = SB[, cols_s], Fig2C = VB,
               T1_Models_13 = T13$C_Models, T2_Verdicts_13 = T13$D_Verdicts, T3_Bootstrap_13 = T13$F_Bootstrap,
               T4_Decomposition_13 = T13$E_Decomposition, T5_Collinearity_13 = T13$G_Collinearity,
               T6_PostHoc_13 = T13$L_PostHoc, T7_Specs_13 = T13$A_Specs, T8_RandomSubset_13b = S,
               T9_LabelAgreement_13b = T13b$F_Label_Agreement, T10_Design_13b = T13b$A_Design,
               T11_Donors = frame, T12_Readings = readings)
    wb <- lapply(wb, function(x) { for (k in names(x)) if (is.numeric(x[[k]])) x[[k]][!is.finite(x[[k]]) & !is.na(x[[k]])] <- NA; x })
    openxlsx::write.xlsx(wb, OUT_XLSX, overwrite = TRUE, firstRow = TRUE)
    n_back <- length(getSheetNames(OUT_XLSX))
    gr("GR-13c-06", if (ok1 && ok2 && n_back == length(wb) && all(file.exists(c(paste0(OUT_F1, c(".pdf", ".png")),
        paste0(OUT_F2, c(".pdf", ".png")), OUT_LEG)))) "PASS" else "FAIL",
       sprintf("outputs written: two figures (pdf, png), table with %d sheets, legend draft", n_back))
    for (f in c(paste0(OUT_F1, c(".pdf", ".png")), paste0(OUT_F2, c(".pdf", ".png")), OUT_XLSX, OUT_LEG))
        cat(sprintf("  %s | %.1f MB\n", basename(f), file.size(f) / 2^20))
    "completed"
}

# ---- FIGURE 1 -----------------------------------------------------------------
coef_panel <- function(fa, sets, series_cols, title, xlab) {
    rowdef <- expand.grid(set = sets, outcome = COG_LAB, stringsAsFactors = FALSE)
    rowdef$y <- seq_len(nrow(rowdef)) + (match(rowdef$outcome, COG_LAB) - 1) * 0.6
    rowdef$label <- paste0(short_cog(rowdef$outcome), ", ", rowdef$set)
    fa <- merge(fa, rowdef[, c("set", "outcome", "y")], by = c("set", "outcome"))
    fl <- ifelse(fa$outcome == "Cognitive status", -1, 1)
    fa$b <- fa$beta * fl
    fa$lo <- ifelse(fl < 0, -fa$ci_hi, fa$ci_lo); fa$hi <- ifelse(fl < 0, -fa$ci_lo, fa$ci_hi)
    fa$series <- factor(fa$series, names(series_cols))
    fa$yy <- fa$y + ifelse(fa$series == names(series_cols)[1], -0.2, 0.2)
    # ASCII keys for the shape scale; the label carries the symbol, so the match
    # does not depend on the session's text encoding.
    fa$readable <- factor(ifelse(is.finite(fa$vif_max) & fa$vif_max < VIF_MAX, "ok", "unreadable"), c("ok", "unreadable"))
    fa$fillc <- ifelse(!is.na(fa$fdr) & fa$fdr < ALPHA & fa$readable == "ok",
                       unname(series_cols[as.character(fa$series)]), "white")
    ggplot(fa) +
        geom_vline(xintercept = 0, linetype = "dashed", colour = "grey45", linewidth = 0.3) +
        geom_segment(aes(x = lo, xend = hi, y = yy, yend = yy, colour = series), linewidth = 0.45) +
        geom_point(aes(x = b, y = yy, colour = series, fill = fillc, shape = readable), size = 1.7, stroke = 0.5) +
        scale_colour_manual(values = series_cols, name = NULL) + scale_fill_identity() +
        scale_shape_manual(values = c(ok = 21, unreadable = 23), labels = c(ok = "VIF < 5", unreadable = "VIF \u22655 (unreadable)"),
                           name = NULL, drop = FALSE) +
        scale_y_reverse(breaks = rowdef$y, labels = rowdef$label, expand = expansion(add = 0.6)) +
        facet_grid(. ~ dataset) + labs(x = xlab, y = NULL, title = title) + th + th_ycat +
        guides(colour = guide_legend(order = 1, override.aes = list(fill = unname(series_cols), shape = 21)),
               shape = guide_legend(order = 2, override.aes = list(fill = "white", colour = "black"))) +
        theme(legend.box = "vertical", legend.spacing.y = unit(0.5, "mm"))
}
fig1 <- function(A, B, C, D) {
    pA <- coef_panel(A, c("cov4", "+BTC", "cov4 + %MT"),
                     c("Without M4 share" = COL$original, "With M4 share" = COL$variant),
                     "Module post-IP% \u03b2, without and with the M4 share", "Standardized \u03b2 of Module post-IP% (95% CI)")
    pB <- coef_panel(B, c("Module post-IP% + %MT (Fig 14D)", "fraction + %MT, cov4", "fraction + %MT, +BTC"),
                     c("Module state (Module post-IP% or within-M4 fraction)" = COL$variant, "%MT of the M4 cells" = COL$original),
                     "Module state and %MT of the same M4 cells in one model", "Standardized \u03b2 in one model (95% CI)")
    # (C) verdict counts
    C$grp <- factor(C$group, c("Cognition", "Pathology", "APOE"))
    cnt <- stats::aggregate(list(k = rep(1, nrow(C))), by = list(grp = C$grp, row = C$row, verdict = C$verdict), sum)
    labs_v <- c("Module only" = "State only", "Second only" = "%MT only", "Both" = "Both", "Neither" = "Neither",
                "Unreadable" = "Unreadable")
    cnt$verdict <- factor(labs_v[cnt$verdict], labs_v)
    fills <- stats::setNames(unname(VERDICT_FILL), labs_v)
    rows <- c("NeoC, Module (Fig 14D)", "NeoC, Module + M4 share", "NeoC, fraction",
              "NonNeoC, Module (Fig 14D)", "NonNeoC, Module + M4 share", "NonNeoC, fraction")
    cnt$row <- factor(cnt$row, rev(rows))
    cnt$txt <- ifelse(cnt$verdict %in% c("State only", "%MT only", "Both"), "white", "black")
    pC <- ggplot(cnt, aes(x = k, y = row, fill = verdict, group = verdict)) +
        geom_col(position = position_stack(reverse = TRUE), colour = "black", linewidth = 0.25, width = 0.7, show.legend = TRUE) +
        geom_text(aes(label = k, colour = txt), position = position_stack(vjust = 0.5, reverse = TRUE), size = GEOM_PT(TXT$in_panel)) +
        scale_colour_identity() +
        scale_fill_manual(values = fills, drop = FALSE, name = "Verdict (FDR<0.05; Unreadable, VIF \u22655)") +
        scale_x_continuous(breaks = function(l) seq(0, floor(l[2]), 1), expand = expansion(mult = c(0, 0.02))) +
        facet_grid(. ~ grp, scales = "free_x", space = "free_x") +
        labs(x = "Number of outcomes", y = NULL, title = "Module state vs %MT, one model") + th + th_ycat +
        theme(panel.grid.major.y = element_blank(), panel.spacing.x = unit(2, "mm")) +
        guides(fill = guide_legend(nrow = 2, title.position = "top"))
    # (D) per donor
    lab <- do.call(rbind, lapply(split(D, D$dataset), function(d)
        data.frame(dataset = d$dataset[1], txt = sprintf("r = %.2f, n = %d", stats::cor(d$frac, d$mt), nrow(d)))))
    pD <- ggplot(D, aes(x = frac, y = mt)) +
        geom_point(aes(fill = mmse), shape = 21, colour = "grey30", size = 1.5, stroke = 0.25) +
        geom_text(data = lab, aes(label = txt), x = -Inf, y = Inf, hjust = -0.08, vjust = 1.4, size = GEOM_PT(TXT$in_panel), inherit.aes = FALSE) +
        scale_fill_gradient(low = COL$variant, high = "#F2F2F2", na.value = "white", name = "MMSE") +
        facet_grid(. ~ dataset) +
        labs(x = "Within-M4 fraction (%)", y = "%MT of the M4 cells (log1p)", title = "Per donor (\u22655 M4 cells)") + th +
        theme(legend.key.width = unit(5, "mm"))
    # C and D each take the full width (larger fonts did not fit side by side)
    fig <- pA / pB / pC / (pD + theme(legend.position = "right")) +
        plot_layout(heights = c(1, 1, 0.66, 0.55)) + plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")")
    dev_pdf <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
    ggsave(paste0(OUT_F1, ".pdf"), fig, width = 180, height = 270, units = "mm", device = dev_pdf)
    ggsave(paste0(OUT_F1, ".png"), fig, width = 180, height = 270, units = "mm", dpi = 300, bg = "white")
}

# ---- FIGURE 2 -----------------------------------------------------------------
hist_panel <- function(S, DR, ros, labeller, title, xlab, flip_cogs) {
    LINES <- c("M4 cells, PVM_3 labels" = COL$m4, "M4 cells, own labels" = COL$own)
    h <- do.call(rbind, lapply(DATASETS, function(ds) do.call(rbind, lapply(ros, function(ro) {
        v <- DR[DR$dataset == ds & DR$valid, ro]; fl <- if (flip_cogs && grepl("_CogS$", ro)) -1 else 1
        data.frame(dataset = ds, readout = ro, value = fl * v[is.finite(v)], stringsAsFactors = FALSE) }))))
    s <- S[S$readout %in% ros, ]
    fl <- ifelse(flip_cogs & grepl("_CogS$", s$readout), -1, 1)
    ln <- rbind(data.frame(dataset = s$dataset, readout = s$readout, value = fl * s$observed_M4_PVM3labels,
                           series = "M4 cells, PVM_3 labels"),
                data.frame(dataset = s$dataset, readout = s$readout, value = fl * s$M4_own_labels, series = "M4 cells, own labels"))
    ln$series <- factor(ln$series, names(LINES))
    tx <- data.frame(dataset = s$dataset, readout = s$readout, txt = sprintf("FDR = %s", fmt_p(s$fdr)))
    lv <- unname(labeller[ros]); flev <- paste0(rep(lv, each = 2), ", ", rep(DATASETS, length(lv)))
    addf <- function(z) { z$facet <- factor(paste0(labeller[z$readout], ", ", z$dataset), flev); z }
    h <- addf(h); ln <- addf(ln); tx <- addf(tx)
    ggplot(h, aes(x = value)) +
        geom_histogram(bins = HIST_BINS, fill = COL$hist_fill, colour = COL$hist_border, linewidth = 0.2) +
        geom_vline(data = ln, aes(xintercept = value, colour = series, linetype = series), linewidth = 0.55) +
        do.call(geom_label, c(list(data = tx, mapping = aes(label = txt), x = -Inf, y = Inf, hjust = -0.05, vjust = 1.15,
                                   size = GEOM_PT(TXT$in_panel), fill = "white", label.padding = unit(0.6, "mm"), inherit.aes = FALSE),
                              label_noborder)) +
        scale_colour_manual(values = LINES, name = NULL, drop = FALSE) +
        scale_linetype_manual(values = c("M4 cells, PVM_3 labels" = "solid", "M4 cells, own labels" = "22"), name = NULL, drop = FALSE) +
        facet_wrap(~ facet, ncol = 2, scales = "free") +
        labs(x = xlab, y = "Random subsets (draws)", title = title) + th
}
fig2 <- function(S, DR, VB) {
    labA <- stats::setNames(short_cog(COG_LAB), paste0("bA_", SHORT))
    pA <- hist_panel(S, DR, names(labA), labA, "Module post-IP% on cognition against random PVM_3 subsets (as Figure 9A)",
                     "Standardized \u03b2 of Module post-IP%", TRUE)
    labB <- c(r_x_mt = "r, Module post-IP% with %MT", r_frac_mt = "r, fraction with %MT", mean_frac = "Mean fraction (%)")
    pB <- hist_panel(S, DR, names(labB), labB, "Collinearity with %MT and the fraction level", "Value", FALSE)
    # (C) verdict shares of the random subsets; the M4 verdict printed beside each bar
    vv <- c(null_module_only = "Module only", null_mt_only = "%MT only", null_both = "Both", null_neither = "Neither",
            null_unreadable = "Unreadable")
    L <- do.call(rbind, lapply(names(vv), function(k) data.frame(dataset = VB$dataset, outcome = VB$outcome,
                                                                 verdict = vv[[k]], share = VB[[k]], stringsAsFactors = FALSE)))
    L$row <- factor(paste0(VB$dataset, ", ", short_cog(VB$outcome))[match(paste(L$dataset, L$outcome), paste(VB$dataset, VB$outcome))],
                    rev(paste0(rep(DATASETS, each = 3), ", ", short_cog(COG_LAB))))
    L$verdict <- factor(L$verdict, vv)
    fills <- stats::setNames(unname(VERDICT_FILL), vv)
    m4 <- data.frame(row = factor(paste0(VB$dataset, ", ", short_cog(VB$outcome)), levels(L$row)),
                     txt = paste0("M4 cells: ", VB$observed_M4_PVM3labels))
    pC <- ggplot(L, aes(x = share, y = row, fill = verdict)) +
        geom_col(position = position_stack(reverse = TRUE), colour = "black", linewidth = 0.25, width = 0.7) +
        geom_text(data = m4, aes(x = 1.02, y = row, label = txt), hjust = 0, size = GEOM_PT(TXT$in_panel), inherit.aes = FALSE) +
        scale_fill_manual(values = fills, drop = FALSE, name = "Verdict of the random subsets (FDR<0.05)") +
        scale_x_continuous(limits = c(0, 1.45), breaks = seq(0, 1, 0.25), labels = function(x) paste0(100 * x, "%"),
                           expand = expansion(mult = c(0, 0))) +
        labs(x = "Share of random subsets", y = NULL,
             title = "Module post-IP% and %MT in one model (Figure 14D): verdicts of the random subsets") + th + th_ycat +
        theme(panel.grid.major.y = element_blank()) + guides(fill = guide_legend(nrow = 1, title.position = "top"))
    fig <- (pA + theme(legend.position = "none")) / pB / pC + plot_layout(heights = c(1, 1, 0.42)) +
        plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")")
    dev_pdf <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
    ggsave(paste0(OUT_F2, ".pdf"), fig, width = 180, height = 200, units = "mm", device = dev_pdf)
    ggsave(paste0(OUT_F2, ".png"), fig, width = 180, height = 200, units = "mm", dpi = 300, bg = "white")
}

# ---- LEGENDS ------------------------------------------------------------------
LEG1 <- c(
    "Supplementary Figure [n]. The M4 share and the Module state against %MT on M4. M4 share: the donor's M4 cells (M4-selected Micro-PVM_3 cells with MT-CO3 >0) as a percentage of the donor's Micro-PVM_3 cells, per region set. Module post-IP% has the same denominator and equals the M4 share x the within-M4 fraction (the donor's Module post-IP cells as a percentage of its M4 cells; donors with at least 5 M4 cells). %MT: log1p of the donor mean of per-cell %MT over the M4 cells, as in Figure 14D. Outcome and predictors were standardized after the complete-case filter; cov4, age at death, sex, PMI and APOE \u03b54 dose; +BTC, cov4 plus Braak, Thal and CERAD. FDR, BH within region set and predictor over the outcomes of the figure each model extends (Figure 9B: the cognition outcomes per block; Figure 14D: the nine clinical outcomes). Unreadable, a full-model VIF \u22655. \u2020Cognitive status is sign-reversed so that lower \u03b2 means worse cognition throughout.",
    "(A) Standardized \u03b2 of Module post-IP% on cognition without (dark) and with (violet) the M4 share, at cov4 and +BTC (as Figure 9B) and at cov4 with %MT (as Figure 14D). With the M4 share in the model, the coefficient is the association of a larger Module post-IP fraction of the M4 cells at a fixed M4 share. Bars, 95% CIs from the t distribution; filled symbols, FDR <0.05 in a readable model; diamonds, unreadable (never filled).",
    "(B) The Module state and %MT of the same M4 cells in one model: Module post-IP% with %MT (Figure 14D), and the within-M4 fraction with %MT at cov4 (BH over the nine clinical outcomes) and at +BTC (BH over the cognition outcomes). The fraction models were added after the M4 share models had been read (post-hoc). A Module post-IP cell has cDFG, the mean of the eight mitochondrially encoded partner genes, at or above IP_y, so the fraction and %MT partly measure the same reads.",
    "(C) Outcomes counted by verdict at FDR <0.05 for the Module state and %MT in one model, over the nine clinical outcomes: Module post-IP% with %MT (Figure 14D), the same with the M4 share added, and the within-M4 fraction with %MT. State only, Module post-IP% or the fraction only.",
    "(D) Within-M4 fraction against %MT of the M4 cells per donor, shaded by MMSE (white, missing MMSE); Pearson r and the number of donors. Values for every panel, the donor bootstrap of the change in \u03b2, the Figure 10C/10D decomposition models and the remaining collinearity values are in the supplementary table.")
LEG2 <- c(
    "Supplementary Figure [n+1]. Are the M4 results a property of the M4 cells? Each donor's M4 cells were compared with random PVM_3 cells of the same number, drawn without replacement within each donor and region from the PVM_3 arm cells (Micro-PVM_3 cells with MT-CO3 >0, the M4 cells included), so the M4 share is the same in every draw (10,000 draws per region set). Every cell, M4 or random, carries its PVM_3 arm labels (PVM_3 Module and PVM_3 IP_y). Module post-IP%, Module post-IP cells in the set as a percentage of the donor's Micro-PVM_3 cells; fraction, Module post-IP cells as a percentage of the cells in the set (donors with at least 5 M4 cells); %MT, log1p of the mean per-cell %MT over the set, as in Figure 14D. Models at cov4, outcome and predictors standardized after the complete-case filter. Histograms (light amber), the random subsets; solid green line, the M4 cells under PVM_3 labels; dashed dark line, the M4 cells under the M4 arm's own labels (its own Module and IP_y), shown for reference only. FDR, BH-corrected two-sided empirical p of the M4 cells under PVM_3 labels against the random subsets; the empirical p is min(1, 2 x the smaller of (1 + draws at or above) / (B + 1) and (1 + draws at or below) / (B + 1)); BH within region set over the three cognition outcomes (panel A) and over the three panel B readouts. The random subsets are matched on cell number per donor and region only, so a difference can reflect any property that separates Mic.13-like cells from other PVM_3 cells, including sequencing depth.",
    "(A) Standardized \u03b2 of Module post-IP% on cognition, as Figure 9A. \u2020Cognitive status is sign-reversed so that lower \u03b2 means worse cognition throughout.",
    "(B) Donor-level Pearson r of Module post-IP% with %MT (donors with M4 cells) and of the fraction with %MT (donors with at least 5 M4 cells), and the mean fraction over donors with at least 5 M4 cells.",
    "(C) Module post-IP% and %MT in one model, as Figure 14D (BH over the nine clinical outcomes; unreadable, a full-model VIF \u22655): share of the random subsets with each verdict on the cognition outcomes, and the verdict of the M4 cells under PVM_3 labels. Values for every panel and the remaining readouts are in the supplementary table.")

# ---- RUN ----------------------------------------------------------------------
.res <- tryCatch(withCallingHandlers(run_13c(), warning = function(w) {
    cat("  WARNING: ", conditionMessage(w), "\n", sep = ""); invokeRestart("muffleWarning") }),
    error = function(e) { cat("\nRUN STOPPED: ", conditionMessage(e), "\n", sep = ""); paste("STOPPED -", conditionMessage(e)) })
G <- if (length(.gr)) do.call(rbind, .gr) else data.frame(id = character(0), status = character(0), detail = character(0))
writeLines(c(strrep("=", 72), "Part 13c - consolidated supplementary figures and table - Guardrail Log",
             paste("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), paste("Write:", OUT_DIR), strrep("=", 72), "",
             sprintf("  [%-4s]  %-10s %s", G$status, G$id, G$detail), "", strrep("=", 72),
             sprintf("  PASS %d | WARN %d | FAIL %d", sum(G$status == "PASS"), sum(G$status == "WARN"), sum(G$status == "FAIL")),
             paste("  Run :", .res), strrep("=", 72)), LEDGER)
cat("\nGuardrail ledger : ", LEDGER, "\nConsole log      : ", CONSOLE, "\nRun: ", .res, "\n", sep = "")
options(.old_width)
sink(); close(.con)
