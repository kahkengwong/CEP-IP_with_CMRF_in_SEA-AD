# R16_Rexach_main_figure.R
# Panels of the main figure for external validation 2 (Rexach et al. 2024, insula), in the plot
# types, colours, and styles of the main manuscript figures (Figures 2, 4, 8, 9, and 15).
#
#   (A) DFG replication (Section 3): r_s against tau_b of MT-CO3 with each of the eight DFGs, with
#       the tier 1 and tier 2 regions and the r_s expected at this depth (SEA-AD Stats_15H).
#   (B) GAM plot (as Figure 2): MT-CO3 against cDFG and against cHKG9 on the same PVM_3 arm cells.
#   (C) CEP-IP plot (as Figure 2): TREP and non-TREP cells of the IP curve, with IP_x, IP_y, and IPRS.
#   (D) Expression space (as Figures 4 and 8): Module cells in the four groups, Recovered, Original
#       only, and Excluded cells.
#   (E) Pitcher plot in copula space (as Figure 8A): the same classes on u_SK and v_SK, with the
#       SKCS-core gates.
#   (F) Module post-IP% per donor by diagnosis (box and points, as Figure 15C).
#   (G) Diagnosis contrasts (forest plot, as Figure 9B): standardized beta and 95% CI, adjusted for
#       age at death, sex, and PMI (cov3). Filled symbols mark BH < 0.05; open symbols mark BH >= 0.05.
#   (H) Post hoc (R17): leave-one-donor-out betas of each contrast, with the observed beta, the share
#       of bootstrap resamples above 0, and the SEA-AD PVM_3 cognitive-status betas.
#   (I) Post hoc (R17): power to detect the SEA-AD NeoC beta (0.465) against the number of donors,
#       with the current donors and the donors needed for 80% power marked.
# The HKG8 control (CEP-IP plot) and the %MT densities are in the supplementary material.
# The panels are written one by one (PDF and PNG) for assembly, and also as one draft figure.
#
# Guardrails (ledger R16_Figure_Rexach/R16_Rexach_Guardrail_Log.txt, PASS / WARN / FAIL):
#   GR-01  input files present             GR-04  panel data consistent with R14 and R15
#   GR-02  R14, R15, and R17 finished,     GR-05  R warnings during the run
#          no FAIL                         GR-06  output files written
#   GR-03  counts equal the R14 results
#
# Inputs, in External_Dataset_2: R14_readings/ (results .rds, readings .xlsx, guardrail log),
#   R15_readings/ (contrasts .xlsx, guardrail log), R17_posthoc/ (post hoc .xlsx, guardrail log), and
#   R13_count_checks/R13_Rexach_MicroPVM3_nuclei.csv.gz.
# Outputs, in R16_Figure_Rexach/: R16_A_... to R16_I_... (PDF and PNG), R16_Rexach_Figure_draft
#   (PDF and PNG), R16_Rexach_figure_source_data.xlsx, R16_Rexach_log.txt, and the guardrail log.
#   The source data holds every value each panel plots, in sheets named by panel letter;
#   build_SuppTable_S21.py turns it into Supplementary Table S21.
#
# Run in R on Windows:
#   source("C:/Users/Wong/Desktop/SEA-AD/External_Dataset_2/R16_Rexach_main_figure.R")

# The data folder: REXACH_ROOT when set (the container sets it to /data), otherwise the Windows folder.
ROOT <- Sys.getenv("REXACH_ROOT", "C:/Users/Wong/Desktop/SEA-AD/External_Dataset_2")

# r_s expected at the depth of the Rexach arm: SEA-AD Stats_15H, PVM_3 downsampled to 3,000 UMIs
# (median r_s over the eight partners 0.630, 95% CI 0.626 to 0.632); readings count check 2.
EXPECTED_RS <- c(estimate = 0.630, lo = 0.626, hi = 0.632)
# Counts of the R14 run (R14_Rexach_log.txt, 2026-10-25), checked in GR-03.
EXPECT <- list(arm = 3417L, module = 2578L, post_ip = 1341L, recovered = 430L, floor_donors = 33L)

# Colours of the manuscript figures and of the original scripts (Part10_CP10K_SKCS_code.r for the GAM
# and CEP-IP plots; Part10e_CP10K_Validation1_code.r for the CMRF classes).
COL <- list(
    trep = "#4B0082", non_trep = "#C0C0C0", gam_dfg = "#4B0082", gam_hkg = "#4F75DE",
    curve_ip = "#FFCC99", ip_line = "#1A237E",
    four = c("Pre-IP TREP" = "#E69F00", "Pre-IP non-TREP" = "#1B3A8C",
             "Post-IP TREP" = "#009E73", "Post-IP non-TREP" = "#D94F4F"),
    recovered = "#7C3AED", original_only = "#8B6914", excluded = "#C0C0C0", gate = "#009E73",
    plum = "#A3408A", expected = "#7F7F7F",
    module_density = "#7F7F7F", other_density = "#C9A27E",
    box_pts = "#F0A850", box_fill = "#EBE3C0", forest = "#782888", header = "#262626",
    lodo = "#B48CC0", seaad = "#7F7F7F",
    power = c("All three vs Control" = "#782888", "AD vs Control" = "#0072B2",
              "bvFTD vs Control" = "#D55E00", "PSP vs Control" = "#7F7F7F"))

IN14 <- file.path(ROOT, "R14_readings")
IN15 <- file.path(ROOT, "R15_readings")
OUT <- file.path(ROOT, "R16_Figure_Rexach")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
LOG <- file.path(OUT, "R16_Rexach_log.txt")
GR_LOG <- file.path(OUT, "R16_Rexach_Guardrail_Log.txt")
FILES_IN <- c(r14_rds = file.path(IN14, "R14_Rexach_results.rds"),
              r14_xlsx = file.path(IN14, "R14_Rexach_readings.xlsx"),
              r14_gr = file.path(IN14, "R14_Rexach_Guardrail_Log.txt"),
              r15_xlsx = file.path(IN15, "R15_Rexach_contrasts.xlsx"),
              r15_gr = file.path(IN15, "R15_Rexach_Guardrail_Log.txt"),
              r17_xlsx = file.path(ROOT, "R17_posthoc", "R17_Rexach_posthoc_beta.xlsx"),
              r17_gr = file.path(ROOT, "R17_posthoc", "R17_Rexach_Guardrail_Log.txt"),
              nuclei = file.path(ROOT, "R13_count_checks", "R13_Rexach_MicroPVM3_nuclei.csv.gz"))

need <- c("cepip", "ggplot2", "gridExtra", "openxlsx")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Install these packages first: ", paste(miss, collapse = ", "))
suppressPackageStartupMessages({ library(cepip); library(ggplot2) })

say <- function(...) cat(sprintf(...), "\n", sep = "")
GR <- new.env()
GR$log <- list(); GR$warnings <- character(0)
gr <- function(id, status, detail) {
    if (length(detail) != 1L || is.na(detail)) detail <- "[detail unavailable]"
    say("  [%s]  %-16s  %s", status, id, detail)
    GR$log[[length(GR$log) + 1L]] <- list(id = id, status = status, detail = detail)
}
gr_count <- function(s) sum(vapply(GR$log, function(x) x$status == s, logical(1)))
gr_write <- function(finished) {
    n_p <- gr_count("PASS"); n_w <- gr_count("WARN"); n_f <- gr_count("FAIL")
    writeLines(c(strrep("=", 72), "R16 Rexach main figure - Guardrail Log",
                 paste("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), strrep("=", 72), "",
                 vapply(GR$log, function(x) sprintf("  [%s]  %-16s  %s", x$status, x$id, x$detail),
                        character(1)),
                 "", strrep("=", 72), sprintf("  Total : %d", n_p + n_w + n_f),
                 sprintf("  PASS  : %d", n_p), sprintf("  WARN  : %d", n_w), sprintf("  FAIL  : %d", n_f),
                 if (!finished) "\n  Result: the run did not finish; see R16_Rexach_log.txt"
                 else if (n_f == 0 && n_w == 0) "\n  Result: no failures, no warnings"
                 else if (n_f == 0) "\n  Result: no failures; review the warnings above"
                 else "\n  Result: failures present", strrep("=", 72)), GR_LOG)
}
stop_if_fail <- function() if (gr_count("FAIL") > 0)
    stop("Guardrail failure; no figure was drawn. See the guardrail log.", call. = FALSE)

# ---- Style ---------------------------------------------------------------------------------------
# Colour helpers of Fig1_vector_code.r (main Figure 1).
darker <- function(hex, f = 0.72) {
    m <- grDevices::col2rgb(hex) * f
    grDevices::rgb(m[1, ], m[2, ], m[3, ], maxColorValue = 255)
}
lighter <- function(hex, f) {
    m <- grDevices::col2rgb(hex) * (1 - f) + 255 * f
    grDevices::rgb(m[1, ], m[2, ], m[3, ], maxColorValue = 255)
}
# Base R lwd to ggplot linewidth (mm) and point stroke.
lw <- function(lwd) lwd * 0.75 / ggplot2::.pt
st <- function(lwd) lwd * 0.75 / ggplot2::.stroke * 2
# Panel A follows Figure 1D: plum #A3408A, region fill lighter(plum, 0.9), points with a fill
# lighter(plum, 0.35) and a border darker(plum, 0.6) of lwd 1.0, dashed gates (lwd 0.8, grey25),
# grey60 leader lines (lwd 0.4), and genes below both tiers as open circles with grey40 labels.
# Tier 1 is a slightly stronger plum and tier 2 a slightly lighter plum.
A_STY <- list(
    reg1 = lighter(COL$plum, 0.84), reg2 = lighter(COL$plum, 0.93),
    fill1 = lighter(COL$plum, 0.15), edge1 = darker(COL$plum, 0.6),
    fill2 = lighter(COL$plum, 0.6), edge2 = darker(COL$plum, 0.8),
    gate1 = "grey25", gate2 = "grey55", lab1 = darker(COL$plum, 0.6), lab2 = darker(COL$plum, 0.9))
theme_fig <- function(bs = 9) {
    theme_minimal(base_size = bs, base_family = "sans") +
        theme(panel.grid.major = element_line(colour = "gray92", linewidth = 0.3),
              panel.grid.minor = element_blank(),
              panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4),
              axis.title = element_text(size = bs, face = "bold"),
              axis.text = element_text(size = bs - 1.5, colour = "black"),
              strip.background = element_rect(fill = COL$header, colour = COL$header),
              strip.text = element_text(colour = "white", face = "bold", size = bs + 0.5,
                                        margin = margin(3, 3, 3, 3)),
              plot.tag = element_text(face = "bold", size = bs + 2),
              plot.tag.position = c(0.005, 0.995),
              legend.text = element_text(size = bs - 1.5), legend.title = element_blank(),
              legend.key.size = unit(0.32, "cm"), legend.background = element_blank(),
              plot.background = element_rect(fill = "white", colour = NA),
              plot.margin = margin(14, 6, 4, 12))
}
# A dark header strip (as in Figures 9 and 15) from a one-level facet.
header <- function(text) facet_wrap(vars(!!text))
legend_in <- function(x, y, just = c(x, y)) {
    if (utils::packageVersion("ggplot2") >= "3.5.0")
        theme(legend.position = "inside", legend.position.inside = c(x, y), legend.justification = just)
    else theme(legend.position = c(x, y), legend.justification = just)
}
fmt_n <- function(n) formatC(n, format = "d", big.mark = ",")
curve_df <- function(fit, n = 300L) {
    xs <- seq(min(fit$data$GOI), max(fit$data$GOI), length.out = n)
    pr <- cepip:::predict_population(fit$model, xs, fit$data, se = TRUE)
    data.frame(x = xs, y = as.numeric(pr$fit), lo = as.numeric(pr$fit) - 1.96 * as.numeric(pr$se.fit),
               hi = as.numeric(pr$fit) + 1.96 * as.numeric(pr$se.fit))
}
lab_gene <- function(g) bquote(bolditalic(.(g)))
save_panel <- function(p, name, w, h) {
    for (ext in c("pdf", "png")) {
        f <- file.path(OUT, paste0(name, ".", ext))
        if (ext == "pdf") grDevices::pdf(f, width = w, height = h, useDingbats = FALSE)
        else grDevices::png(f, width = w, height = h, units = "in", res = 600)
        print(p); grDevices::dev.off()
    }
    file.path(OUT, paste0(name, c(".pdf", ".png")))
}

main <- function() {
t0 <- Sys.time()
say("R16_Rexach_main_figure.R started %s | ggplot2 %s", format(t0), packageVersion("ggplot2"))

# ---- Guardrails on the inputs ----------------------------------------------------------------------
say("\n== Guardrails on the inputs")
ok <- file.exists(FILES_IN)
gr("GR-01-inputs", if (all(ok)) "PASS" else "FAIL",
   if (all(ok)) "R14, R15, and R17 outputs and the R13 nuclei table are present"
   else paste("missing:", paste(basename(FILES_IN[!ok]), collapse = ", ")))
stop_if_fail()
fails <- vapply(FILES_IN[c("r14_gr", "r15_gr", "r17_gr")], function(f) {
    g <- readLines(f, warn = FALSE)
    n <- suppressWarnings(as.integer(sub(".*FAIL\\s*:\\s*", "", grep("^\\s*FAIL\\s*:", g, value = TRUE))))
    if (any(grepl("did not finish", g)) || length(n) != 1L) NA_integer_ else n
}, integer(1))
gr("GR-02-upstream", if (all(fails %in% 0L)) "PASS" else "FAIL",
   sprintf("R14 FAIL count %s | R15 FAIL count %s | R17 FAIL count %s (NA = did not finish)",
           fails[1], fails[2], fails[3]))
stop_if_fail()

r14 <- readRDS(FILES_IN[["r14_rds"]])
res <- r14$res$PVM_3; hkg8 <- r14$res$HKG8; hkc <- r14$res$HKG_control
nuc <- utils::read.csv(gzfile(FILES_IN[["nuclei"]]), stringsAsFactors = FALSE)
cl <- res$cmrf$cells
grp <- cepip_groups(res)
got <- c(arm = nrow(cl), module = sum(cl$in_module), post_ip = sum(grp$group == "PostIP"),
         recovered = sum(cl$cell_class == "Recovered"), floor_donors = length(r14$floor_donors))
exp_v <- unlist(EXPECT)
gr("GR-03-counts", if (all(got == exp_v)) "PASS" else "WARN",
   paste(sprintf("%s %d (%d)", names(got), got, exp_v), collapse = " | "))

dfg <- openxlsx::read.xlsx(FILES_IN[["r14_xlsx"]], "02_S3_DFGs")
con <- openxlsx::read.xlsx(FILES_IN[["r15_xlsx"]], "01_Contrasts")
con <- con[con$arm == "PVM_3", ]
don <- openxlsx::read.xlsx(FILES_IN[["r15_xlsx"]], "04_Donors_PVM_3")
r17_obs <- openxlsx::read.xlsx(FILES_IN[["r17_xlsx"]], "02_Observed_cov3")
r17_lodo <- openxlsx::read.xlsx(FILES_IN[["r17_xlsx"]], "05_LODO_donors")
r17_boot <- openxlsx::read.xlsx(FILES_IN[["r17_xlsx"]], "06_Bootstrap")
r17_need <- openxlsx::read.xlsx(FILES_IN[["r17_xlsx"]], "09_Donors_needed")
r17_sea <- openxlsx::read.xlsx(FILES_IN[["r17_xlsx"]], "10_SEAAD_reference")
cons_ok <- nrow(dfg) == 8L && nrow(con) == 4L && nrow(don) == length(r14$floor_donors) &&
    all(abs(don$pct_PostIP - r14$donors$PVM_3$pct_PostIP[match(don$donor, r14$donors$PVM_3$donor)]) < 1e-9) &&
    all(res$ep$id %in% nuc$cell_id)
r17_dev <- max(abs(r17_obs$beta_std - con$beta_std[match(r17_obs$contrast, con$contrast)]))
r17_ok <- nrow(r17_obs) == 4L && is.finite(r17_dev) && r17_dev < 1e-8 && nrow(r17_sea) == 2L &&
    nrow(r17_lodo) == sum(r17_obs$n)
gr("GR-04-consistency", if (cons_ok && r17_ok) "PASS" else "FAIL",
   sprintf("DFG rows %d | PVM_3 contrasts %d | donors %d, Module post-IP%% equal to R14 | arm cells found in the R13 table %s | R17 betas equal R15 (max diff %.1e), LODO rows %d",
           nrow(dfg), nrow(con), nrow(don), all(res$ep$id %in% nuc$cell_id), r17_dev, nrow(r17_lodo)))
stop_if_fail()

src <- list(); files <- character(0)

# ---- (A) DFG replication ----------------------------------------------------------------------------
dfg$tier <- ifelse(dfg$passes_tier1, "Tier 1", ifelse(dfg$passes_tier2, "Tier 2 only", "Neither"))
# Gene labels sit in a column to the right of the points, joined by thin lines, so that close
# points do not overlap.
dfg <- dfg[order(-dfg$tau_b, -dfg$rs), ]
lab_x <- max(dfg$rs) + 0.06
dfg$lab_y <- seq(max(dfg$tau_b) + 0.02, min(dfg$tau_b) - 0.02, length.out = nrow(dfg))
xr <- c(min(0.45, min(dfg$rs) - 0.03), lab_x + 0.11)
yr <- c(min(0.33, min(dfg$tau_b) - 0.05), max(0.66, max(dfg$tau_b) + 0.06))
dfg$tier <- factor(dfg$tier, levels = c("Tier 1", "Tier 2 only", "Neither"))
pA <- ggplot(dfg, aes(x = rs, y = tau_b)) +
    annotate("rect", xmin = 0.5, xmax = Inf, ymin = 0.4, ymax = Inf, fill = A_STY$reg2) +
    annotate("rect", xmin = 0.6, xmax = Inf, ymin = 0.5, ymax = Inf, fill = A_STY$reg1) +
    annotate("rect", xmin = EXPECTED_RS[["lo"]], xmax = EXPECTED_RS[["hi"]], ymin = -Inf, ymax = Inf,
             fill = COL$expected, alpha = 0.25) +
    geom_vline(xintercept = EXPECTED_RS[["estimate"]], colour = COL$expected, linetype = "dotted",
               linewidth = lw(0.8)) +
    annotate("segment", x = c(0.5, 0.5), xend = c(0.5, Inf), y = c(0.4, 0.4), yend = c(Inf, 0.4),
             colour = A_STY$gate2, linetype = "dashed", linewidth = lw(0.8)) +
    annotate("segment", x = c(0.6, 0.6), xend = c(0.6, Inf), y = c(0.5, 0.5), yend = c(Inf, 0.5),
             colour = A_STY$gate1, linetype = "dashed", linewidth = lw(0.8)) +
    geom_segment(aes(xend = lab_x - 0.005, yend = lab_y), colour = "grey60", linewidth = lw(0.4)) +
    geom_point(aes(fill = tier, colour = tier), shape = 21, size = 2.4, stroke = st(1.0)) +
    geom_text(data = dfg[dfg$tier != "Neither", ], aes(x = lab_x, y = lab_y, label = DFG),
              fontface = "italic", size = 2.7, hjust = 0, colour = "black") +
    geom_text(data = dfg[dfg$tier == "Neither", ], aes(x = lab_x, y = lab_y, label = DFG),
              fontface = "italic", size = 2.7, hjust = 0, colour = "grey40") +
    annotate("text", x = 0.6, y = yr[2], label = "Tier 1", hjust = -0.1, vjust = 1.2, size = 2.9,
             fontface = "bold", colour = A_STY$lab1) +
    annotate("text", x = 0.5, y = yr[2], label = "Tier 2", hjust = -0.1, vjust = 1.2, size = 2.9,
             fontface = "bold", colour = A_STY$lab2) +
    annotate("text", x = EXPECTED_RS[["estimate"]], y = yr[2], label = "Expected~r[s]", parse = TRUE,
             hjust = -0.08, vjust = 2.6, size = 2.6, colour = "gray30") +
    scale_fill_manual(values = c("Tier 1" = A_STY$fill1, "Tier 2 only" = A_STY$fill2, "Neither" = "white"),
                      labels = function(l) sprintf("%s (n=%d)", l, as.integer(table(dfg$tier)[l])),
                      name = NULL) +
    scale_colour_manual(values = c("Tier 1" = A_STY$edge1, "Tier 2 only" = A_STY$edge2, "Neither" = "grey40"),
                        labels = function(l) sprintf("%s (n=%d)", l, as.integer(table(dfg$tier)[l])),
                        name = NULL) +
    coord_cartesian(xlim = xr, ylim = yr) +
    labs(x = expression(bold("Spearman") ~ bold(r[s])), y = expression(bold("Kendall") ~ bold(tau[b])),
         tag = "(A)") +
    header(sprintf("DFG replication: %d of 8 at tier 1", sum(dfg$passes_tier1))) +
    theme_fig() + legend_in(0.02, 0.02, c(0, 0))
files <- c(files, save_panel(pA, "R16_A_DFG_replication", 3.3, 3.2))
src[["A_DFG_replication"]] <- data.frame(dfg[, setdiff(names(dfg), "lab_y")],
                                         expected_rs = EXPECTED_RS[["estimate"]],
                                         expected_rs_lo = EXPECTED_RS[["lo"]], expected_rs_hi = EXPECTED_RS[["hi"]])

# ---- (B) GAM plot (styles of Part10_CP10K_SKCS_code.r, GAM curve plots) -----------------------------------
fd <- res$fit$data; fh <- hkc$fit$data
m <- match(fd$id, fh$id)
dB <- rbind(data.frame(x = fd$GOI, y = fd$Expression, set = "cDFG"),
            data.frame(x = fd$GOI, y = fh$Expression[m], set = "cHKG9"))
cB <- rbind(cbind(curve_df(res$fit), set = "cDFG"), cbind(curve_df(hkc$fit), set = "cHKG9"))
lab_B <- c(cDFG = sprintf("cDFG: DE=%.1f%%", 100 * res$fit$de),
           cHKG9 = sprintf("cHKG9 (control): DE=%.1f%%", 100 * max(0, hkc$fit$de)))
cols_B <- c(cDFG = COL$gam_dfg, cHKG9 = COL$gam_hkg)
pB <- ggplot() +
    geom_point(data = dB, aes(x = x, y = y, colour = set),
               position = position_jitter(width = 0.01, height = 0, seed = 123), alpha = 0.25, size = 0.45) +
    geom_line(data = cB, aes(x = x, y = y, colour = set), linewidth = 1) +
    geom_ribbon(data = cB, aes(x = x, ymin = lo, ymax = hi, fill = set), alpha = 0.20, colour = NA) +
    scale_colour_manual(values = cols_B, labels = lab_B) +
    scale_fill_manual(values = cols_B, guide = "none") +
    guides(colour = guide_legend(override.aes = list(size = 1.8, alpha = 1, linewidth = 1))) +
    labs(x = lab_gene("MT-CO3"), y = expression(bold("cDFG / cHKG9")), tag = "(B)") +
    header(sprintf("GAM plot: PVM_3 arm (n=%s)", fmt_n(nrow(fd)))) +
    theme_fig() + legend_in(0.02, 0.98, c(0, 1))
files <- c(files, save_panel(pB, "R16_B_GAM_plot", 3.3, 3.2))
# Cell values without the plot's horizontal jitter (width 0.01, seed 123).
src[["B_GAM_cells"]] <- data.frame(cell_id = fd$id, donor = fd$donor, MT_CO3 = fd$GOI, cDFG = fd$Expression,
                                   cHKG9 = fh$Expression[m])
src[["B_GAM_curves"]] <- data.frame(set = cB$set, MT_CO3 = cB$x, fit = cB$y, ci95_lo = cB$lo, ci95_hi = cB$hi,
                                    DE = ifelse(cB$set == "cDFG", res$fit$de, hkc$fit$de))

# ---- (C) CEP-IP plot (styles of Part10_CP10K_SKCS_code.r, CEP-IP plot) -----------------------------------------
# TREP #4B0082 and non-TREP #C0C0C0 points (alpha 0.45), and the GAM curve and its 95% band in #FFCC99
# (line, then the band at alpha 0.2), on the curve that defines the IP.
# Boxed inset of plotmath lines, right-aligned, in the bottom right corner of the panel.
ip_box_grob <- function(lines, size_mm) {
    th <- gridExtra::ttheme_minimal(core = list(
        fg_params = list(parse = TRUE, hjust = 1, x = 0.97, fontsize = size_mm * ggplot2::.pt),
        bg_params = list(fill = NA, col = NA), padding = grid::unit(c(4, 1.6), "pt")))
    tg <- gridExtra::tableGrob(matrix(lines, ncol = 1), theme = th, rows = NULL, cols = NULL)
    tg <- gtable::gtable_add_grob(tg, grid::roundrectGrob(r = grid::unit(2, "pt"),
                                  gp = grid::gpar(fill = "white", col = NA)),
                                  t = 1, b = nrow(tg), l = 1, r = ncol(tg), z = -Inf, name = "background")
    tg <- gtable::gtable_add_grob(tg, grid::roundrectGrob(r = grid::unit(2, "pt"),
                                  gp = grid::gpar(fill = NA, col = "black", lwd = 0.6)),
                                  t = 1, b = nrow(tg), l = 1, r = ncol(tg), z = Inf, name = "border")
    # annotation_custom() replaces the viewport of the grob it is given, so the corner viewport sits
    # on an inner gTree and the outer gTree takes the panel viewport.
    inner <- grid::gTree(children = grid::gList(tg),
                         vp = grid::viewport(x = grid::unit(1, "npc") - grid::unit(3, "pt"),
                                             y = grid::unit(3, "pt"), just = c(1, 0),
                                             width = sum(tg$widths), height = sum(tg$heights)))
    grid::gTree(children = grid::gList(inner))
}
cepip_panel <- function(r, xg, yl, title, tag, ip_lines) {
    ip <- r$ip; e <- ip$ep
    d <- data.frame(x = e$x, y = e$y, cls = ifelse(ip$is_top, "TREP", "non-TREP"))
    d <- d[order(d$cls == "TREP"), ]
    cv <- curve_df(ip$fit)
    labs_c <- c(TREP = sprintf("TREP (n=%s)", fmt_n(sum(ip$is_top))),
                "non-TREP" = sprintf("non-TREP (n=%s)", fmt_n(sum(!ip$is_top))))
    # Plotmath lines, so that x in IPx and y in IPy are italic.
    box <- if (isTRUE(ip$has_ip))
        c(sprintf('"DE=%.1f%%"', 100 * ip$fit$de),
          sprintf('IP*italic(x)*"=%.2f (95%% CI: %.2f-%.2f)"', ip$ip_x, ip$ip_ci[1], ip$ip_ci[2]),
          sprintf('IP*italic(y)*"=%.2f"', ip$ip_y),
          sprintf('"IPRS=%.3f (%s tier)"', ip$iprs$IPRS, ip$iprs$tier))
    else c(sprintf('"DE=%.1f%%"', 100 * max(0, ip$fit$de)), '"No IP"',
           sprintf('"IPRS=%.3f (%s tier)"', ip$iprs$IPRS, ip$iprs$tier))
    p <- ggplot()
    if (ip_lines && isTRUE(ip$has_ip) && all(is.finite(ip$ip_ci)))
        p <- p + annotate("rect", xmin = ip$ip_ci[1], xmax = ip$ip_ci[2], ymin = -Inf, ymax = Inf,
                          fill = COL$ip_line, alpha = 0.10)
    p <- p +
        geom_point(data = d, aes(x = x, y = y, colour = cls), size = 0.9, alpha = 0.45) +
        geom_line(data = cv, aes(x = x, y = y), colour = COL$curve_ip, linewidth = 1) +
        geom_ribbon(data = cv, aes(x = x, ymin = lo, ymax = hi), fill = COL$curve_ip, alpha = 0.2)
    if (ip_lines && isTRUE(ip$has_ip))
        p <- p + geom_vline(xintercept = ip$ip_x, colour = COL$ip_line, linetype = "dashed", linewidth = 0.5) +
            geom_hline(yintercept = ip$ip_y, colour = "black", linetype = "dashed", linewidth = 0.5)
    p + annotation_custom(ip_box_grob(box, 2.3), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
        scale_colour_manual(values = c(TREP = COL$trep, "non-TREP" = COL$non_trep),
                            limits = c("TREP", "non-TREP"), labels = labs_c) +
        guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
        labs(x = lab_gene(xg), y = yl, tag = tag) +
        header(title) + theme_fig() + legend_in(0.02, 0.98, c(0, 1))
}
ip <- res$ip
pC <- cepip_panel(res, "MT-CO3", expression(bold("cDFG")), "CEP-IP plot: PVM_3 arm", "(C)", TRUE)
files <- c(files, save_panel(pC, "R16_C_CEPIP_plot", 3.3, 3.2))
src[["C_CEPIP_cells"]] <- data.frame(cell_id = ip$ep$id, donor = ip$ep$donor, MT_CO3 = ip$ep$x, cDFG = ip$ep$y,
                                     class = ifelse(ip$is_top, "TREP", "non-TREP"))
cvC <- curve_df(ip$fit)
src[["C_CEPIP_curve"]] <- data.frame(MT_CO3 = cvC$x, fit = cvC$y, ci95_lo = cvC$lo, ci95_hi = cvC$hi)
src[["C_CEPIP_IP"]] <- data.frame(IP_x = ip$ip_x, IP_x_CI_lower = ip$ip_ci[1], IP_x_CI_upper = ip$ip_ci[2],
                                  IP_y = ip$ip_y, IPRS = ip$iprs$IPRS, tier = ip$iprs$tier,
                                  TREP = sum(ip$is_top), non_TREP = sum(!ip$is_top), DE_IP_curve = ip$fit$de)

# ---- (D, E) Module classes in expression and copula space (styles of Part10e_CP10K_Validation1_code.r) -------
# OrigRef cells: four-group circles (alpha 0.58). Recovered: violet #7C3AED diamonds with a white border
# (shape 23, stroke 0.15, alpha 0.93). Original only: tan #8B6914 inverted triangles with a white border
# (shape 25, alpha 0.88). Excluded: gray #C0C0C0 (alpha 0.45), as the non-TREP cells of panel C. Drawn in that order from the back:
# Excluded, Original only, OrigRef (non-TREP before TREP), Recovered.
lev <- c(names(COL$four), "Recovered", "Original only", "Excluded")
cl$cls <- ifelse(cl$cell_class == "Original + Refined", cl$group_orig, cl$cell_class)
cl$cls <- factor(cl$cls, levels = lev)
draw <- c("Excluded" = 1, "Original only" = 2, "Pre-IP non-TREP" = 3, "Post-IP non-TREP" = 4,
          "Pre-IP TREP" = 5, "Post-IP TREP" = 6, "Recovered" = 7)
cl <- cl[order(draw[as.character(cl$cls)]), ]
cnt <- table(cl$cls)
lab_E <- sprintf("%s (n=%s)", lev, fmt_n(as.integer(cnt[lev])))
v_col <- c(COL$four, "Recovered" = "white", "Original only" = "white", "Excluded" = COL$excluded)
v_fill <- c(COL$four, "Recovered" = COL$recovered, "Original only" = COL$original_only, "Excluded" = COL$excluded)
v_shape <- c(16, 16, 16, 16, 23, 25, 19)
# Excluded matches the panel C points (shape 19, size 0.9, stroke 0.5). This layer has stroke 0.15,
# so its size is raised by the stroke difference: 0.9 + 0.35 * .stroke / 2 / .pt = 1.13.
v_size <- c(0.8, 0.8, 0.8, 0.8, 1.5, 1.4, 0.9 + 0.35 * ggplot2::.stroke / 2 / ggplot2::.pt)
v_alpha <- c(0.58, 0.58, 0.58, 0.58, 0.93, 0.88, 0.45)
class_layer <- function(xv, yv) list(
    geom_point(data = cl, aes(x = .data[[xv]], y = .data[[yv]], colour = cls, fill = cls, shape = cls,
                              size = cls, alpha = cls), stroke = 0.15),
    scale_colour_manual(values = v_col, limits = lev, labels = lab_E, drop = FALSE, name = NULL),
    scale_fill_manual(values = v_fill, limits = lev, labels = lab_E, drop = FALSE, name = NULL),
    scale_shape_manual(values = v_shape, limits = lev, labels = lab_E, drop = FALSE, name = NULL),
    scale_size_manual(values = v_size, limits = lev, labels = lab_E, drop = FALSE, name = NULL),
    scale_alpha_manual(values = v_alpha, limits = lev, labels = lab_E, drop = FALSE, name = NULL),
    guides(colour = guide_legend(override.aes = list(size = c(2, 2, 2, 2, 2.4, 2.2, 2), alpha = 1,
                                                     colour = c(COL$four, "white", "white", COL$excluded)))))
pD <- ggplot() + class_layer("x", "y") +
    geom_vline(xintercept = ip$ip_x, colour = COL$ip_line, linetype = "dashed", linewidth = 0.4) +
    geom_hline(yintercept = ip$ip_y, colour = "black", linetype = "dashed", linewidth = 0.5) +
    labs(x = lab_gene("MT-CO3"), y = expression(bold("cDFG")), tag = "(D)") +
    header(sprintf("Expression space: Module n=%s", fmt_n(sum(cl$in_module)))) +
    theme_fig() + legend_in(0.995, 0.005, c(1, 0))
files <- c(files, save_panel(pD, "R16_D_expression_space", 3.3, 3.2))
gs <- res$cmrf$gates_SK
pE <- ggplot() + class_layer("u_SK", "v_SK") +
    geom_vline(xintercept = gs[["u"]], colour = COL$gate, linewidth = 0.6) +
    geom_hline(yintercept = gs[["v"]], colour = COL$gate, linewidth = 0.6) +
    annotate("text", x = gs[["u"]], y = 0, label = 'bold(S[i]*" > 0.05")', parse = TRUE, hjust = -0.08, vjust = 0,
             size = 2.8, fontface = "bold") +
    annotate("text", x = 0, y = gs[["v"]], label = 'bold(K[i]*" > 0.5")', parse = TRUE, hjust = 0, vjust = -0.4,
             size = 2.8, fontface = "bold") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(x = expression(bold(u[SK])), y = expression(bold(v[SK])), tag = "(E)") +
    header("Pitcher plot (copula space)") + theme_fig() + theme(legend.position = "none")
files <- c(files, save_panel(pE, "R16_E_pitcher_copula", 3.3, 3.2))
src[["D_Expression_space"]] <- data.frame(cell_id = cl$id, donor = cl$donor, class = as.character(cl$cls),
                                          cell_class = cl$cell_class, in_module = cl$in_module,
                                          MT_CO3 = cl$x, cDFG = cl$y)
src[["E_Pitcher_cells"]] <- data.frame(cell_id = cl$id, class = as.character(cl$cls), u_SK = cl$u_SK, v_SK = cl$v_SK)
src[["E_Pitcher_gates"]] <- data.frame(gate = c("S_i > 0.05 (u_SK)", "K_i > 0.5 (v_SK)"),
                                       value = c(gs[["u"]], gs[["v"]]))

# ---- (F) Module post-IP% per donor by diagnosis ------------------------------------------------------------
lv_dx <- c("Control", "AD", "bvFTD", "PSP")
don$dx <- factor(don$dx, levels = lv_dx)
nlab <- sprintf("n=%d", as.integer(table(don$dx)))
set.seed(42)
pF <- ggplot(don, aes(x = dx, y = pct_PostIP)) +
    geom_boxplot(width = 0.44, outlier.shape = NA, fill = COL$box_fill, colour = "black", linewidth = 0.4) +
    geom_jitter(width = 0.096, height = 0, shape = 21, fill = COL$box_pts, colour = "black", stroke = 0.25,
                size = 1.8, alpha = 0.9) +
    annotate("text", x = seq_along(lv_dx), y = -Inf, label = nlab, vjust = -0.6, size = 3.0, colour = "gray35") +
    scale_y_continuous(limits = c(-6, max(don$pct_PostIP) * 1.08)) +
    labs(x = NULL, y = expression(bold("Module post-IP% per donor")), tag = "(F)") +
    header("Module post-IP% by diagnosis") + theme_fig() +
    theme(axis.text.x = element_text(face = "bold", size = 8.5, colour = "black"))
files <- c(files, save_panel(pF, "R16_F_donor_by_diagnosis", 3.3, 3.2))
src[["F_Donors"]] <- don[, intersect(c("donor", "dx", "pct_PostIP", "n_PostIP", "n_total", "age_at_death", "sex", "pmi"),
                                      names(don))]

# ---- (G) Diagnosis contrasts --------------------------------------------------------------------------------
lab_ct <- c(All_vs_Control = "All three vs Control", AD_vs_Control = "AD vs Control",
            bvFTD_vs_Control = "bvFTD vs Control", PSP_vs_Control = "PSP vs Control")
con$label <- factor(lab_ct[con$contrast], levels = rev(lab_ct))
con$sig <- !is.na(con$BH_over_4_contrasts) & con$BH_over_4_contrasts < 0.05
con$ptxt <- sprintf("n=%d, p=%s", con$n, formatC(con$p_value, format = "f", digits = 3))
lim <- max(abs(c(con$ci_low, con$ci_high)), na.rm = TRUE) * 1.15
pG <- ggplot(con, aes(y = label, x = beta_std)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "gray45", linewidth = 0.4) +
    geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0, colour = COL$forest, linewidth = 0.6) +
    geom_point(aes(fill = sig), shape = 21, size = 2.4, colour = COL$forest, stroke = 0.8) +
    geom_text(aes(x = lim, label = ptxt), hjust = 1, vjust = -0.9, size = 2.6, colour = "gray25") +
    scale_fill_manual(values = c("TRUE" = COL$forest, "FALSE" = "white"), guide = "none") +
    coord_cartesian(xlim = c(-lim, lim)) +
    labs(x = expression(bold("Standardized") ~ bolditalic(beta) ~ bold("(95% CI)")), y = NULL, tag = "(G)") +
    header("Diagnosis contrasts (cov3)") + theme_fig()
files <- c(files, save_panel(pG, "R16_G_diagnosis_forest", 3.3, 3.2))
src[["G_Contrasts"]] <- data.frame(con[, c("contrast", "n", "n_disease", "n_control", "beta_std", "se", "ci_low",
                                           "ci_high", "p_value", "BH_over_4_contrasts")], model = "cov3")

# ---- (H) Post hoc: leave-one-donor-out (LODO) stability (R17) -------------------------------------------------------------
# Each small point is the beta with one donor left out (leave-one-donor-out, LODO). The diamond is the observed beta. The grey band
# spans the two SEA-AD PVM_3 cognitive-status betas (NeoC and NonNeoC, cov4).
r17_lodo$label <- factor(lab_ct[r17_lodo$contrast], levels = rev(lab_ct))
r17_obs$label <- factor(lab_ct[r17_obs$contrast], levels = rev(lab_ct))
r17_boot$label <- factor(lab_ct[r17_boot$contrast], levels = rev(lab_ct))
pos <- tapply(r17_lodo$beta_without > 0, r17_lodo$contrast, sum)
r17_obs$htxt <- sprintf("LODO>0: %d of %d; bootstrap>0: %.0f%%", as.integer(pos[r17_obs$contrast]), r17_obs$n,
                        100 * r17_boot$share_above_0[match(r17_obs$contrast, r17_boot$contrast)])
sea_rng <- range(r17_sea$beta)
limH <- max(abs(c(r17_lodo$beta_without, r17_obs$beta_std, sea_rng)), na.rm = TRUE) * 1.25
set.seed(7)
pH <- ggplot() +
    annotate("rect", xmin = sea_rng[1], xmax = sea_rng[2], ymin = -Inf, ymax = Inf, fill = COL$seaad, alpha = 0.25) +
    annotate("text", x = mean(sea_rng), y = Inf, label = "SEA-AD", vjust = 1.4, size = 2.6, colour = "gray30",
             fontface = "bold") +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "gray45", linewidth = 0.4) +
    geom_jitter(data = r17_lodo, aes(x = beta_without, y = label), height = 0.12, width = 0, shape = 21,
                fill = COL$lodo, colour = "white", stroke = 0.2, size = 1.6, alpha = 0.85) +
    geom_point(data = r17_obs, aes(x = beta_std, y = label), shape = 23, size = 2.6, fill = COL$forest,
               colour = "black", stroke = 0.4) +
    geom_label(data = r17_obs, aes(x = limH, y = label, label = htxt), hjust = 1, nudge_y = -0.36, size = 2.4,
               colour = "gray25", fill = "white", label.size = 0, label.padding = unit(0.08, "lines")) +
    coord_cartesian(xlim = c(-limH * 0.6, limH)) +
    labs(x = expression(bold("Standardized") ~ bolditalic(beta) ~ bold("(one donor left out)")), y = NULL,
         tag = "(H)") +
    header("Leave-one-donor-out (post hoc)") + theme_fig()
files <- c(files, save_panel(pH, "R16_H_LODO_stability", 3.3, 3.2))
src[["H_LODO_donors"]] <- r17_lodo[, c("contrast", "donor", "dx", "pct_PostIP", "beta_without", "change",
                                       "dfbetas", "cooks_d")]
src[["H_LODO_summary"]] <- data.frame(contrast = r17_obs$contrast, n = r17_obs$n, beta_std = r17_obs$beta_std,
                                      lodo_positive = as.integer(pos[r17_obs$contrast]),
                                      bootstrap_share_above_0 = r17_boot$share_above_0[match(r17_obs$contrast, r17_boot$contrast)],
                                      seaad_beta_NeoC = r17_sea$beta[r17_sea$set == "NeoC"],
                                      seaad_beta_NonNeoC = r17_sea$beta[r17_sea$set == "NonNeoC"])

# ---- (I) Post hoc: power at the SEA-AD beta (R17) ---------------------------------------------------------------------
# Power of the cov3 model to detect the SEA-AD NeoC beta, from the Rexach SE scaled by the residual df,
# with the same group ratio and covariates (as R17 2c). Filled circle: current donors. Open circle:
# donors needed for 80% power.
k_par <- 5L
pow_at <- function(n, se0, df0, beta) {
    dfn <- n - k_par; if (dfn < 1) return(0)
    se <- se0 * sqrt(df0 / dfn); tc <- stats::qt(0.975, dfn); ncp <- beta / se
    1 - stats::pt(tc, dfn, ncp) + stats::pt(-tc, dfn, ncp)
}
b_tgt <- r17_sea$beta[r17_sea$set == "NeoC"]
ns <- 8:90
cvI <- do.call(rbind, lapply(seq_len(nrow(r17_obs)), function(k) data.frame(
    label = unname(lab_ct[r17_obs$contrast[k]]), n = ns,
    power = vapply(ns, pow_at, numeric(1), se0 = r17_obs$se[k], df0 = r17_obs$df_resid[k], beta = b_tgt))))
ptI <- do.call(rbind, lapply(seq_len(nrow(r17_obs)), function(k) {
    nd <- r17_need$donors_needed[r17_need$contrast == r17_obs$contrast[k] & r17_need$seaad_set == "NeoC"]
    data.frame(label = unname(lab_ct[r17_obs$contrast[k]]), kind = c("now", "needed"), n = c(r17_obs$n[k], nd),
               power = c(pow_at(r17_obs$n[k], r17_obs$se[k], r17_obs$df_resid[k], b_tgt),
                         pow_at(nd, r17_obs$se[k], r17_obs$df_resid[k], b_tgt)))
}))
cvI$label <- factor(cvI$label, levels = lab_ct); ptI$label <- factor(ptI$label, levels = lab_ct)
need_chk <- vapply(seq_len(nrow(r17_obs)), function(k) {
    pw <- vapply(ns, pow_at, numeric(1), se0 = r17_obs$se[k], df0 = r17_obs$df_resid[k], beta = b_tgt)
    ns[which(pw >= 0.80)[1]]
}, numeric(1))
need_r17 <- r17_need$donors_needed[match(paste(r17_obs$contrast, "NeoC"), paste(r17_need$contrast, r17_need$seaad_set))]
gr("GR-04b-power", if (isTRUE(all(need_chk == need_r17))) "PASS" else "WARN",
   sprintf("donors needed recomputed here equal R17 (NeoC beta %.3f): %s", b_tgt,
           paste(sprintf("%s %s/%s", r17_obs$contrast, need_chk, need_r17), collapse = " | ")))
pI <- ggplot(cvI, aes(x = n, y = power, colour = label)) +
    geom_hline(yintercept = 0.8, linetype = "dashed", colour = "gray45", linewidth = 0.4) +
    geom_line(linewidth = 0.7) +
    geom_point(data = ptI[ptI$kind == "now", ], size = 2.2, shape = 21, aes(fill = label), colour = "black",
               stroke = 0.3) +
    geom_point(data = ptI[ptI$kind == "needed", ], size = 2.2, shape = 21, fill = "white", stroke = 0.7) +
    annotate("text", x = max(ns), y = 0.8, label = "80%", hjust = 1, vjust = -0.5, size = 2.6, colour = "gray30") +
    annotate("text", x = min(ns), y = 1, label = "Filled: current donors\nOpen: donors for 80% power",
             hjust = 0, vjust = 1, size = 2.5, colour = "gray25", lineheight = 0.95) +
    scale_colour_manual(values = COL$power, name = NULL) +
    guides(colour = guide_legend(override.aes = list(shape = NA, linewidth = 0.9))) +
    scale_fill_manual(values = COL$power, guide = "none") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(x = expression(bold("Donors in the contrast")),
         y = expression(bold("Power at SEA-AD") ~ bolditalic(beta) ~ bold("= 0.465")), tag = "(I)") +
    header("Power analysis (post hoc)") + theme_fig() + legend_in(0.98, 0.02, c(1, 0))
files <- c(files, save_panel(pI, "R16_I_power_SEAAD", 3.3, 3.2))
ct_of <- stats::setNames(names(lab_ct), lab_ct)
src[["I_Power_curves"]] <- data.frame(contrast = unname(ct_of[as.character(cvI$label)]), donors = cvI$n,
                                      power = cvI$power, target_beta = b_tgt)
src[["I_Power_points"]] <- data.frame(contrast = unname(ct_of[as.character(ptI$label)]),
                                      point = ifelse(ptI$kind == "now", "current donors", "donors for 80% power"),
                                      donors = ptI$n, power = ptI$power, target_beta = b_tgt)

# ---- Draft figure and source data -------------------------------------------------------------------------------
say("\n== Draft figure")
g <- gridExtra::arrangeGrob(pA, pB, pC, pD, pE, pF, pG, pH, pI, ncol = 3,
                            top = grid::textGrob("External validation 2 (Rexach et al. 2024, insula)",
                                                 gp = grid::gpar(fontface = "bold", fontsize = 12)))
for (ext in c("pdf", "png")) {
    f <- file.path(OUT, paste0("R16_Rexach_Figure_draft.", ext))
    if (ext == "pdf") grDevices::pdf(f, width = 10, height = 10.2, useDingbats = FALSE)
    else grDevices::png(f, width = 10, height = 10.2, units = "in", res = 300)
    grid::grid.newpage(); grid::grid.draw(g); grDevices::dev.off()
    files <- c(files, f)
}
wb <- openxlsx::createWorkbook()
bold <- openxlsx::createStyle(textDecoration = "bold")
for (nm in names(src)) {
    sh <- substr(nm, 1, 31)
    openxlsx::addWorksheet(wb, sh)
    openxlsx::writeData(wb, sh, src[[nm]], headerStyle = bold)
    openxlsx::freezePane(wb, sh, firstRow = TRUE)
}
xl <- file.path(OUT, "R16_Rexach_figure_source_data.xlsx")
openxlsx::saveWorkbook(wb, xl, overwrite = TRUE)
files <- c(files, xl)

gr("GR-05-R_warnings", if (length(GR$warnings)) "WARN" else "PASS",
   if (length(GR$warnings)) sprintf("%d R warning(s): %s", length(GR$warnings), paste(unique(GR$warnings), collapse = " || "))
   else "no R warnings during the run")
gr("GR-06-outputs", if (all(file.exists(files))) "PASS" else "FAIL",
   sprintf("%d of %d files written", sum(file.exists(files)), length(files)))
say("\nGuardrail tally: %d PASS / %d WARN / %d FAIL", gr_count("PASS"), gr_count("WARN"), gr_count("FAIL"))
say("Done in %.0f s", as.numeric(difftime(Sys.time(), t0, units = "secs")))
invisible(TRUE)
}

# Panels of an earlier run are removed first, so renamed or failed panels cannot be left behind.
old_out <- list.files(OUT, pattern = "^R16_([A-I]_.*|Rexach_Figure_draft)\\.(pdf|png)$|^R16_Rexach_figure_source_data\\.xlsx$",
                      full.names = TRUE)
invisible(file.remove(old_out))
zz <- file(LOG, open = "wt")
sink(zz, split = TRUE)
finished <- FALSE
tryCatch({
    withCallingHandlers(main(),
        warning = function(w) {
            GR$warnings <- c(GR$warnings, conditionMessage(w))
            say("R warning: %s", conditionMessage(w))
            invokeRestart("muffleWarning")
        },
        error = function(e) say("ERROR: %s", conditionMessage(e)))
    finished <- TRUE
}, finally = {
    gr_write(finished)
    say("Guardrail log written: %s", GR_LOG)
    sink(); close(zz)
})
