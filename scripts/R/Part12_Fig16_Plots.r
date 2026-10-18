# =============================================================================
# Part12_Fig16_Plots.r
# SEA-AD CEP-IP/CMRF, Part 12, Figure 16 plots (Figure 16)
#   Draws the seven Figure 16 panels from the stats and plot-data files written
#   by sections 13 to 17 of Part12_Fig16_Analysis.r, after checking them
#   against their logs. Computes no statistic.
#
# Usage (from the SEA-AD folder, or set SEAAD_ROOT), after Part12_Fig16_Analysis.r
#   source("Part12_rerun/Part12_Fig16_Plots.r")
#
# Outputs (Part12_rerun/08_Fig16/Fig16/): figures (PDF, PNG), log.
# =============================================================================

source(file.path(Sys.getenv("SEAAD_ROOT", getwd()), "Part12_rerun", "Part12_03_Config.r"))

SCRIPT <- "Part12_Fig16_Plots.r"
FIG_DIR <- "08_Fig16"   # every Figure 16 section folder sits inside this one

# =============================================================================
# SECTION 18: Figure 16
# Panels 16A to 16G; arms in the order M4, PVM_3, PVM_1, HKG8, NeoC before
# NonNeoC. Filled marker and * at FDR < 0.05, ** at FDR < 0.01, "ns" otherwise.
# =============================================================================
p12_section_18 <- function() {
    STEP <- "18"
    D <- p12_step_dirs(STEP, "Fig16", parent = FIG_DIR)
    p12_clean(D$out)
    p12_start(STEP, SCRIPT, D)

    PAL <- P12_CFG$palette
    ARM_COL <- unlist(PAL$arms); GRP_COL <- unlist(PAL$groups)
    INK <- PAL$neutral$dark; MID <- PAL$neutral$mid; LIGHT <- PAL$neutral$light
    ACCENT <- PAL$neutral$accent
    ARMS <- P12_ARMS
    B <- as.integer(P12_CFG$inference$bootstrap_B)
    stars <- function(q) ifelse(!is.finite(q), "", ifelse(q < 0.01, "**", ifelse(q < 0.05, "*", "")))
    band <- function(q) ifelse(!is.finite(q), "", ifelse(q < 0.05, stars(q), "ns"))

    ## Sizes (absolute; every figure sets par(cex = 1)).
    T_TITLE <- 0.92; T_AXIS <- 0.75; T_AXTITLE <- 0.83; T_LAB <- 0.70; T_KEY <- 0.70
    T_NOTE <- 0.56; T_VAL <- 0.72; PT <- 1.3; PT_LWD <- 1.2; WH_LWD <- 1.2
    GREY_T <- "grey30"; BAND <- "#F2F2F2"; ZERO <- "#808080"; NULL_COL <- "#C0392B"
    QLAB <- c(module_post_pct = "Module post-IP%", module_resid_tq4 = "Module post-IP%, TQR",
              mito_log1p = "%MT (plain)", contamination = "Contamination",
              ctrl_pct_null = "Matched share, null")
    QCOL <- c(module_post_pct = ACCENT, module_resid_tq4 = ACCENT, mito_log1p = INK,
              contamination = MID, ctrl_pct_null = LIGHT)
    QPCH <- c(module_post_pct = 21, module_resid_tq4 = 23, mito_log1p = 21, contamination = 24,
              ctrl_pct_null = 22)
    QEDGE <- c(module_post_pct = ACCENT, module_resid_tq4 = ACCENT, mito_log1p = "black",
               contamination = "black", ctrl_pct_null = "black")
    ## 16B bar fills; TQR is the accent purple tinted 60% towards white.
    QBAR <- c(module_post_pct = ACCENT, module_resid_tq4 = "#C3B1D7", mito_log1p = "#A6A6A6",
              contamination = "#5E5E5E", ctrl_pct_null = "#A6A6A6")
    QHATCH <- c(module_post_pct = FALSE, module_resid_tq4 = FALSE, mito_log1p = FALSE,
                contamination = FALSE, ctrl_pct_null = TRUE)

    fig_dev <- function(stem, w, h, draw) {
        ok <- TRUE
        for (kind in c("pdf", "png")) {
            f <- file.path(D$out, paste0(stem, ".", kind))
            p12_output(f)
            if (kind == "pdf") grDevices::pdf(f, width = w, height = h, useDingbats = FALSE)
            else grDevices::png(f, width = w, height = h, units = "in", res = 300,
                                type = if (.Platform$OS.type == "windows") "windows" else "cairo")
            r <- tryCatch({ par(cex = 1, family = "sans"); draw(); TRUE }, error = function(e) {
                p12_say(sprintf("    %s %s: %s", stem, kind, conditionMessage(e))); FALSE })
            grDevices::dev.off()
            ok <- ok && r
        }
        p12_chk(ok, sprintf("%s drawn", stem), "pdf and png")
    }
    rd <- function(folder, file) utils::read.csv(p12_upstream_file(file.path(FIG_DIR, folder), file),
                                                 check.names = FALSE, stringsAsFactors = FALSE)
    ## Greedy word wrap on measured widths.
    fit_lines <- function(txt, w_in, cex, font = 1) {
        wid <- function(z) strwidth(z, units = "inches", cex = cex, font = font) / par("cex")
        unlist(lapply(txt, function(z) {
            w <- strsplit(z, " ", fixed = TRUE)[[1]]; out <- character(0); cur <- ""
            for (x in w) {
                cand <- if (nzchar(cur)) paste(cur, x) else x
                if (nzchar(cur) && wid(cand) > w_in) { out <- c(out, cur); cur <- x } else cur <- cand
            }
            c(out, cur)
        }))
    }
    ## Title at the top of the device and a small note at the bottom, both on a
    ## device-wide overlay drawn last.
    frame_text <- function(title, note = NULL, key = NULL) {
        din <- par("din")
        tl <- fit_lines(title, din[1] - 0.25, T_TITLE, 2)
        nl <- if (is.null(note)) character(0) else fit_lines(note, din[1] - 0.25, T_NOTE)
        par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mai = c(0, 0, 0, 0), new = TRUE, xpd = NA,
            cex = 1)
        plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "",
             xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
        iy <- 1 / din[2]
        for (k in seq_along(tl)) text(0.5, 1 - (0.08 + 0.19 * (k - 1)) * iy, tl[k], adj = c(0.5, 1),
                                      cex = T_TITLE, font = 2)
        for (k in seq_along(nl)) text(0.5, (0.07 + 0.12 * (length(nl) - k)) * iy, nl[k],
                                      adj = c(0.5, 0), cex = T_NOTE, col = GREY_T)
        if (!is.null(key)) {
            ## `overlay`: a second legend at the same place (hatching over a key box).
            ov <- key$overlay; key$overlay <- NULL
            args <- c(list(x = 0.5, y = 1 - (0.10 + 0.19 * length(tl)) * iy, xjust = 0.5, yjust = 1,
                           horiz = is.null(key$ncol), bty = "n", cex = T_KEY), key)
            do.call(legend, args)
            if (!is.null(ov)) do.call(legend, utils::modifyList(args, ov))
        }
    }
    top_in <- function(title) { din <- par("din"); 0.14 + 0.19 * length(fit_lines(title, din[1] - 0.25, T_TITLE, 2)) }
    note_in <- function(note) { din <- par("din"); 0.10 + 0.12 * length(fit_lines(note, din[1] - 0.25, T_NOTE)) }
    ## Hatched rectangle (published style for TQR and the matched share).
    hrect <- function(x0, y0, x1, y1, fill, hatch, lwd = 0.8) {
        rect(x0, y0, x1, y1, col = fill, border = "black", lwd = lwd)
        if (hatch) rect(x0, y0, x1, y1, col = "#2D2D2D", border = NA, density = 22, angle = 45,
                        lwd = 0.7)
    }
    whisker <- function(x0, x1, y, h, horiz = TRUE) {
        if (!is.finite(x0) || !is.finite(x1)) return(invisible(NULL))
        if (horiz) { segments(x0, y, x1, y, lwd = WH_LWD); segments(c(x0, x1), y - h, c(x0, x1), y + h, lwd = WH_LWD) }
        else { segments(y, x0, y, x1, lwd = WH_LWD); segments(y - h, c(x0, x1), y + h, c(x0, x1), lwd = WH_LWD) }
    }

    main <- function() {
        p12_say("SECTION A -- inputs")
        A <- rd("16AB_DonorStructure", "Stats_16A.csv")
        Bt <- rd("16AB_DonorStructure", "Stats_16B.csv")
        C <- rd("16C_LabelAgreement", "Stats_16C.csv")
        Dd <- rd("16D_WithinDonorRho", "Stats_16D.csv")
        Dp <- rd("16D_WithinDonorRho", "Plot_16D_donors.csv")
        E <- rd("16E_Offset", "Stats_16E.csv")
        Ff <- rd("16FG_LibraryFields", "Stats_16F.csv")
        G <- rd("16FG_LibraryFields", "Stats_16G.csv")
        wm <- Dd$estimate[Dd$quantity == "within_median_raw" & Dd$measure == "MT-CO3 vs cDFG"]
        pm <- stats::median(Dp$rs_raw[Dp$pair == "MT-CO3 vs cDFG"], na.rm = TRUE)
        p12_chk(isTRUE(abs(wm - pm) < 1e-9), "16D plot data agree with the stats file",
                sprintf("within-donor median %.6f against %.6f", pm, wm))
        labs <- c(QLAB, ARMS)
        p12_chk(p12_ascii_ok(labs), "figure labels are ASCII", sprintf("%d labels", length(labs)))
        ## Figure text is typed in this file, so the whole file is checked too.
        src <- readLines(file.path(P12_RERUN, SCRIPT), warn = FALSE)
        p12_chk(p12_ascii_ok(src), "figure script text is ASCII", sprintf("%d lines", length(src)))
        if (P12$n[["FAIL"]] > 0L) p12_abort("checks failed; figures not drawn")
        p12_say("")
        p12_say("SECTION B -- figures")

        # == 16A =============================================================
        fig_dev("Fig16A_donor_share", 6.93, 5.4, function() {
            ttl <- "Is Module post-IP% donor-structured?"
            note <- paste("Donor share of variance across donor-region units (donor first) against",
                          "10,000 shuffles of donor labels within region; FDR BH within arm over the",
                          "five quantities; filled marker and * FDR < 0.05, ** FDR < 0.01, ns FDR >= 0.05.")
            ## Rows sorted by the gap over the null 95th percentile; untestable last.
            gap <- A$estimate - A$hi
            v <- A[order(-ifelse(is.finite(gap), gap, -Inf)), ]
            n <- nrow(v)
            lab <- sprintf("%s %s", v$arm, QLAB[v$quantity])
            lab_in <- max(strwidth(lab, units = "inches", cex = T_LAB)) + 0.12
            par(mai = c(0.55 + note_in(note), lab_in, top_in(ttl) + 0.42, 0.42), cex = 1)
            plot(NA, xlim = c(0, 100), ylim = c(0.5, n + 0.5), yaxt = "n", xaxt = "n", xlab = "",
                 ylab = "", xaxs = "i", yaxs = "i")
            abline(v = seq(0, 100, 20), col = "grey92")
            axis(1, at = seq(0, 100, 20), cex.axis = T_AXIS, mgp = c(3, 0.4, 0), tcl = -0.25)
            mtext("Donor share of variance (%)", side = 1, line = 1.5, font = 2, cex = T_AXTITLE)
            y <- rev(seq_len(n))
            ## A faint rule along each row, under the data.
            abline(h = y, col = "grey90", lwd = 0.7)
            text(par("usr")[1] - diff(grconvertX(c(0, 0.06), "inches", "user")), y, lab,
                 adj = c(1, 0.5), cex = T_LAB, xpd = NA)
            for (i in seq_len(n)) {
                r <- v[i, ]
                if (is.finite(r$hi)) segments(100 * r$hi, y[i] - 0.32, 100 * r$hi, y[i] + 0.32,
                                              col = NULL_COL, lwd = 1.6)
                if (!is.finite(r$estimate)) {
                    text(50, y[i], "untestable", cex = T_LAB, col = MID, font = 3); next
                }
                ## Arrow stops at the marker edge; none when the gap is shorter than the head.
                if (is.finite(r$hi)) {
                    rad <- diff(grconvertX(c(0, 0.075), "inches", "user"))
                    x0 <- 100 * r$hi; x1 <- 100 * r$estimate
                    if (abs(x1 - x0) > rad + diff(grconvertX(c(0, 0.08), "inches", "user")))
                        arrows(x0, y[i], x1 - sign(x1 - x0) * rad, y[i], length = 0.07, angle = 25,
                               col = "#5A5A5A", lwd = 1.1)
                }
                sig <- is.finite(r$fdr) && r$fdr < 0.05
                points(100 * r$estimate, y[i], pch = QPCH[[r$quantity]], cex = PT, lwd = PT_LWD,
                       col = QEDGE[[r$quantity]], bg = if (sig) QCOL[[r$quantity]] else "white")
                if (is.finite(r$hi)) {
                    ## Right of tick and point, or left when it would run past the axis end.
                    lb <- sprintf("%+.0f pp%s", 100 * (r$estimate - r$hi), if (i == 1L) " over null" else "")
                    w <- strwidth(lb, cex = T_NOTE + 0.04)
                    xr <- 100 * max(r$estimate, r$hi) + 2
                    if (xr + w <= 99) text(xr, y[i], lb, adj = c(0, 0.5), cex = T_NOTE + 0.04, col = GREY_T)
                    else text(100 * min(r$estimate, r$hi) - 2, y[i], lb, adj = c(1, 0.5), cex = T_NOTE + 0.04,
                              col = GREY_T)
                }
            }
            u <- par("usr")
            text(u[2] + diff(grconvertX(c(0, 0.22), "inches", "user")), y, band(v$fdr), cex = T_LAB,
                 xpd = NA, font = 2)
            text(u[2] + diff(grconvertX(c(0, 0.22), "inches", "user")), n + 0.9, "FDR", cex = T_LAB,
                 xpd = NA, font = 2)
            box()
            frame_text(ttl, note, key = list(legend = c(QLAB, "Null 95th percentile"), ncol = 3,
                                             pch = c(QPCH, NA), pt.bg = c(QCOL, NA), col = c(QEDGE, NULL_COL),
                                             lty = c(rep(NA, 5), 1), lwd = c(rep(PT_LWD, 5), 1.6),
                                             pt.cex = 1.1))
        })

        # == 16B =============================================================
        fig_dev("Fig16B_split_half", 6.16, 4.8, function() {
            ttl <- "Does a donor keep its Module post-IP% across regions?"
            note <- paste("Split-half rs, median over 200 random splits of each donor's units;",
                          "bars 95% donor-bootstrap intervals; * FDR < 0.05, ** FDR < 0.01, BH within arm.")
            v <- Bt[order(-ifelse(is.finite(Bt$estimate), Bt$estimate, -Inf)), ]
            n <- nrow(v)
            lab <- sprintf("%s %s", v$arm, QLAB[v$quantity])
            lab_in <- max(strwidth(lab, units = "inches", cex = T_LAB)) + 0.12
            par(mai = c(0.55 + note_in(note), lab_in, top_in(ttl) + 0.30, 0.25), cex = 1)
            xl <- c(floor(min(-0.2, v$lo, na.rm = TRUE) / 0.2 - 1e-9) * 0.2, max(1, v$hi, na.rm = TRUE))
            plot(NA, xlim = xl, ylim = c(0.4, n + 0.6), yaxt = "n", xaxt = "n", xlab = "", ylab = "",
                 xaxs = "i", yaxs = "i")
            ## No vertical grid; the zero line is the only vertical reference.
            tk <- round(seq(xl[1], 1, 0.2), 1)
            axis(1, at = tk, labels = format(tk), cex.axis = T_AXIS, mgp = c(3, 0.4, 0), tcl = -0.25)
            abline(v = 0, col = ZERO)
            mtext("Split-half rs", side = 1, line = 1.5, font = 2, cex = T_AXTITLE)
            y <- rev(seq_len(n))
            text(par("usr")[1] - diff(grconvertX(c(0, 0.06), "inches", "user")), y, lab,
                 adj = c(1, 0.5), cex = T_LAB, xpd = NA)
            for (i in seq_len(n)) {
                r <- v[i, ]
                if (!is.finite(r$estimate)) { text(0.4, y[i], "untestable", cex = T_LAB, col = MID, font = 3); next }
                hrect(0, y[i] - 0.42, r$estimate, y[i] + 0.42, QBAR[[r$quantity]], QHATCH[[r$quantity]],
                      lwd = 0.55)
                whisker(r$lo, r$hi, y[i], 0.12)
                s <- stars(r$fdr)
                if (nzchar(s)) text(max(r$hi, r$estimate, na.rm = TRUE), y[i], s, pos = 4, offset = 0.2,
                                    cex = T_VAL, font = 2)
            }
            box()
            qs <- c("module_post_pct", "mito_log1p", "module_resid_tq4")
            ## Hatched entries: the bar's fill with dark hatching over it.
            frame_text(ttl, note, key = list(legend = QLAB[qs], fill = QBAR[qs], border = "black",
                                             overlay = list(fill = "#2D2D2D",
                                                            density = ifelse(QHATCH[qs], 55, 0),
                                                            border = NA, text.col = NA)))
        })

        # == 16C =============================================================
        fig_dev("Fig16C_label_agreement", 4.58, 4.68, function() {
            ttl <- "Do M4 and HKG8 label the same cells as PVM_3?"
            kp <- C[C$quantity == "kappa" & C$level == "pooled", ]
            note <- sprintf(paste("Shared cells in regions with 20 or more; pooled kappa M4 %.3f [%.3f, %.3f]",
                                  "(threshold > 0.60), HKG8 %.3f [%.3f, %.3f] (threshold |kappa| < 0.20).",
                                  "Segments too narrow to label are given above the bar, left to right."),
                            kp$estimate[kp$arm == "M4"], kp$lo[kp$arm == "M4"], kp$hi[kp$arm == "M4"],
                            kp$estimate[kp$arm == "HKG8"], kp$lo[kp$arm == "HKG8"], kp$hi[kp$arm == "HKG8"])
            ## Rows from the top: HKG8 NeoC, HKG8 NonNeoC, M4 NeoC, M4 NonNeoC.
            rows <- expand.grid(ds = P12_DATASETS, arm = c("HKG8", "M4"), stringsAsFactors = FALSE)
            sh <- c("share_both", "share_test_only", "share_ref_only", "share_neither")
            M <- t(vapply(seq_len(nrow(rows)), function(i) vapply(sh, function(q) {
                v <- C$estimate[C$arm == rows$arm[i] & C$level == rows$ds[i] & C$quantity == q]
                if (length(v) && is.finite(v)) 100 * v else 0 }, 0), numeric(4)))
            kv <- vapply(seq_len(nrow(rows)), function(i) {
                v <- C$estimate[C$arm == rows$arm[i] & C$level == rows$ds[i] & C$quantity == "kappa"]
                if (length(v)) v[1] else NA_real_ }, 0)
            COL4 <- c("#1B9E77", "#F2A950", "#B4451F", "#B8B8B8")
            lab <- sprintf("%s vs\nPVM_3\n%s", rows$arm, rows$ds)
            par(mai = c(0.55 + note_in(note), 0.75, top_in(ttl) + 0.52, 0.62), cex = 1)
            n <- nrow(M)
            plot(NA, xlim = c(0, 100), ylim = c(0.4, n + 0.6), axes = FALSE, xlab = "", ylab = "",
                 xaxs = "i", yaxs = "i")
            for (i in seq_len(n)) {
                y <- n - i + 1; x0 <- 0
                small <- character(0); small_x <- numeric(0)
                for (j in 1:4) {
                    ## Bars 0.61 of the row pitch.
                    rect(x0, y - 0.305, x0 + M[i, j], y + 0.305, col = COL4[j], border = "black", lwd = 0.8)
                    ## Percentage at cex 0.95, else 0.78; if neither fits, printed above the bar.
                    pl <- sprintf("%.0f%%", M[i, j])
                    fits <- function(cx) strwidth(pl, cex = cx, font = 2) + strwidth("m", cex = 0.5) <= M[i, j]
                    cx <- if (fits(0.95)) 0.95 else if (fits(0.78)) 0.78 else NA_real_
                    if (is.finite(cx)) text(x0 + M[i, j] / 2, y, pl, col = "white", font = 2, cex = cx)
                    else if (M[i, j] > 0) { small <- c(small, pl); small_x <- c(small_x, x0 + M[i, j] / 2) }
                    x0 <- x0 + M[i, j]
                }
                if (length(small)) {
                    txt <- paste(small, collapse = " / ")
                    hw <- strwidth(txt, cex = T_NOTE + 0.06) / 2
                    xc <- min(max(mean(small_x), hw), 100 - hw)
                    text(xc, y + 0.335, txt, adj = c(0.5, 0), cex = T_NOTE + 0.06, col = INK, xpd = NA)
                }
            }
            axis(1, at = c(0, 25, 50, 75, 100), labels = c("0", "25%", "50%", "75%", "100%"),
                 cex.axis = T_AXIS, mgp = c(3, 0.4, 0), tcl = -0.25)
            mtext("Cells present in both arms (%)", side = 1, line = 1.5, font = 2, cex = T_AXTITLE)
            text(-2, rev(seq_len(n)), lab, adj = c(1, 0.5), cex = T_LAB, font = 2, xpd = NA)
            text(103, n + 0.62, "kappa:", adj = c(0, 0), font = 2, cex = T_LAB, xpd = NA)
            text(103, rev(seq_len(n)), ifelse(is.finite(kv), sprintf("%.3f", kv), "NA"), adj = c(0, 0.5),
                 cex = 0.85, font = ifelse(rows$arm == "M4", 2, 1), xpd = NA)
            frame_text(ttl, note, key = list(legend = c("Both arms", "Test arm only", "Reference only",
                                                        "Neither"), fill = COL4, border = "black",
                                             ncol = 2))
        })

        # == 16D =============================================================
        fig_dev("Fig16D_within_donor_rs", 4.4, 4.4, function() {
            ttl <- "Does the MT-CO3 to cDFG relationship hold inside single donors?"
            note <- paste("PVM_3 arm. Pre-set marks on the within-donor median: 0.75 and 0.40.",
                          "Dashed outlines: control bars where they overlap the within-donor bars.")
            r <- Dp$rs_raw[Dp$pair == "MT-CO3 vs cDFG" & is.finite(Dp$rs_raw)]
            rc <- Dp$rs_raw[Dp$pair == "ACTB vs cHKG8" & is.finite(Dp$rs_raw)]
            g <- function(q, m) Dd$estimate[Dd$quantity == q & Dd$measure == m]
            pooled <- g("pooled_raw", "MT-CO3 vs cDFG"); between <- g("between_raw", "MT-CO3 vs cDFG")
            wmed <- g("within_median_raw", "MT-CO3 vs cDFG")
            br <- seq(min(c(r, rc, pooled, between), na.rm = TRUE),
                      max(c(r, rc, pooled, between), na.rm = TRUE), length.out = 25)
            h1 <- graphics::hist(r, breaks = br, plot = FALSE)
            h2 <- graphics::hist(rc, breaks = br, plot = FALSE)
            par(mai = c(0.55 + note_in(note), 0.62, top_in(ttl) + 0.05, 0.15), cex = 1)
            ## y axis to 21 donors (legend clear of the bars), higher only if needed.
            plot(NA, xlim = range(br), ylim = c(0, max(21, max(h1$counts, h2$counts) * 1.05)), las = 1,
                 xlab = "", ylab = "", cex.axis = T_AXIS, mgp = c(3, 0.45, 0), tcl = -0.25, yaxs = "i")
            plot(h2, col = "#B8B8B8", border = "black", add = TRUE)
            plot(h1, col = grDevices::adjustcolor("#1B9E77", 0.9), border = "black", add = TRUE)
            ## Control bars outlined again on top, dashed, so shared bins stay visible.
            plot(h2, col = NA, border = "#2D2D2D", lty = 2, lwd = 1.1, add = TRUE)
            abline(v = pooled, col = "#B4451F", lwd = 2)
            abline(v = between, col = "#F2A950", lwd = 2)
            abline(v = wmed, col = "black", lty = 2, lwd = 1.8)
            mtext("rs, MT-CO3 against cDFG", side = 1, line = 1.5, font = 2, cex = T_AXTITLE)
            mtext("Donors", side = 2, line = 1.9, font = 2, cex = T_AXTITLE)
            legend("topleft", bty = "o", bg = "white", box.col = "grey80", cex = T_NOTE + 0.06,
                   fill = c("#1B9E77", "#B8B8B8", NA, NA, NA), border = c("black", "black", NA, NA, NA),
                   col = c(NA, NA, "#B4451F", "#F2A950", "black"), lty = c(NA, NA, 1, 1, 2),
                   lwd = c(NA, NA, 2, 2, 1.8),
                   legend = c(sprintf("Within-donor rs (%d donors)", length(r)),
                              "Control panel (same donors)",
                              sprintf("Pooled (%+.4f)", pooled), sprintf("Between-donor (%+.4f)", between),
                              sprintf("Within-donor (median %+.4f)", wmed)))
            frame_text(ttl, note)
        })

        # == 16E =============================================================
        ## HKG8 is untestable in 16E (too few donors per branch) and named in the note.
        EA <- setdiff(ARMS, "HKG8")
        rh <- E[E$arm == "HKG8" & E$quantity == "d_offset_share_pp", ]
        np <- E$n_donors[E$arm == "HKG8" & E$level == "post-IP" & E$quantity == "offset_share_pct"]
        nq <- E$n_donors[E$arm == "HKG8" & E$level == "pre-IP" & E$quantity == "offset_share_pct"]
        hk_note <- if (nrow(rh) && !is.finite(rh$estimate))
            sprintf(" HKG8 not shown: untestable (pre-IP %s and post-IP %s donors; 10 needed).",
                    if (length(nq)) nq else "NA", if (length(np)) np else "NA") else " HKG8 not shown."
        fig_dev("Fig16E_offset_share", 3.69, 4.255, function() {
            ttl <- "Is the offset share equal between the branches?"
            note <- paste0("Post-IP minus pre-IP with its 90% donor-bootstrap interval, which the 10 pp ",
                           "margin reads; * FDR < 0.05, ** FDR < 0.01 for non-zero, BH over the arms. ",
                           "a = median intercept over donors (MT UMIs at zero depth).", hk_note)
            ## Pre-IP drawn first (left), then post-IP.
            GC <- c(pre = GRP_COL[["Module pre-IP"]], post = GRP_COL[["Module post-IP"]])
            par(mai = c(0.40 + note_in(note), 0.62, top_in(ttl) + 0.30, 0.10), cex = 1)
            xl <- c(0.55, length(EA) + 0.45)
            ## Bar width and offset fixed in inches.
            upi <- par("pin")[1] / diff(xl)
            half <- 0.169 / upi; off <- 0.189 / upi
            ## y range leaves 0.62 in for the text and 0.20 in for the bar label.
            ph <- par("pin")[2]
            sh <- E$estimate[E$quantity == "offset_share_pct"]
            smax <- if (any(is.finite(sh))) max(sh[is.finite(sh)], 0) else 100
            ytop <- max(20, smax * ph / (ph - 0.82))
            xpos <- seq_along(EA)
            plot(NA, xlim = xl, ylim = c(0, ytop), xaxt = "n", yaxt = "n", xlab = "", ylab = "",
                 xaxs = "i", yaxs = "i")
            yin <- function(inch) inch * ytop / ph
            yt <- seq(0, ytop, 20)
            ## Faint grid under the bars, horizontal lines below the text band only.
            abline(h = yt[yt > 0 & yt < ytop - yin(0.60)], col = "grey90", lwd = 0.6)
            abline(v = seq_len(length(EA) - 1) + 0.5, col = "grey90", lwd = 0.6)
            axis(2, at = yt, las = 1, cex.axis = T_AXIS, mgp = c(3, 0.45, 0), tcl = -0.25)
            mtext("Offset share of MT UMIs (%)", side = 2, line = 1.9, font = 2, cex = T_AXTITLE)
            for (i in seq_along(EA)) {
                a <- EA[i]
                for (j in 1:2) {
                    lv <- c("pre-IP", "post-IP")[j]
                    s <- E$estimate[E$arm == a & E$level == lv & E$quantity == "offset_share_pct"]
                    ic <- E$estimate[E$arm == a & E$level == lv & E$quantity == "intercept"]
                    xx <- xpos[i] + (j - 1.5) * 2 * off
                    if (length(s) && is.finite(s)) {
                        rect(xx - half, 0, xx + half, max(s, 0), col = GC[[j]], border = "black", lwd = 0.8)
                        text(xx, max(s, 0) + yin(0.09), sprintf("%.0f%%", s), cex = T_NOTE + 0.12, col = INK,
                             font = 2)
                        if (s > yin(0.50) && is.finite(ic)) text(xx, s / 2, sprintf("a=%.0f", ic),
                                                                 cex = T_NOTE + 0.16, srt = 90, col = "white")
                    }
                }
                mtext(a, side = 1, at = xpos[i], line = 0.4, font = 2, cex = T_AXTITLE)
                r <- E[E$arm == a & E$quantity == "d_offset_share_pp", ]
                r9 <- E[E$arm == a & E$quantity == "d_offset_share_pp_90", ]
                if (nrow(r) && is.finite(r$estimate)) {
                    VL <- c("within the margin" = "within 10 pp",
                            "non-zero but within the margin" = "non-zero, within 10 pp",
                            "non-zero and not shown to be within the margin" = "non-zero, not within 10 pp",
                            "underpowered: neither non-zero nor shown to be within the margin" = "underpowered")
                    text(xpos[i], ytop - yin(0.14), sprintf("%+.1f pp%s", r$estimate, stars(r$fdr)),
                         cex = T_VAL + 0.08, font = 2)
                    text(xpos[i], ytop - yin(0.33), sprintf("90%% [%+.1f, %+.1f]", r9$lo, r9$hi),
                         cex = T_NOTE + 0.04)
                    text(xpos[i], ytop - yin(0.50), if (r9$measure %in% names(VL)) VL[[r9$measure]] else r9$measure,
                         cex = T_NOTE + 0.02)
                } else text(xpos[i], ytop - yin(0.33), "untestable", cex = T_NOTE + 0.04, col = MID, font = 3)
            }
            box()
            frame_text(ttl, note, key = list(legend = c("Module pre-IP", "Module post-IP"), fill = GC,
                                             border = "black"))
        })

        # == 16F =============================================================
        FQ <- c(module = "Module post-IP%", ncells = "Log cell count", mito = "%MT")
        FC <- c(module = ACCENT, ncells = "#7F7F7F", mito = "#BFBFBF")
        fig_dev("Fig16F_intronic_rate", 4.8, 6.67, function() {
            ttl <- "Does Module post-IP% track the library intronic rate?"
            note <- paste("Spearman rs with the library intronic mapping rate; bars 95% donor-bootstrap",
                          "intervals; * FDR < 0.05, ** FDR < 0.01, BH within arm and level.")
            R <- Ff[Ff$quantity %in% names(FQ), ]
            LVS <- c("between donors", "within donor")
            ## Panel heights in proportion to each y range plus room for the text, so
            ## one rs unit has the same length in both panels.
            dr <- lapply(LVS, function(lv) {
                d <- R[R$level == lv, ]
                c(min(c(d$lo, d$estimate, 0), na.rm = TRUE), max(c(d$hi, d$estimate, 0), na.rm = TRUE))
            })
            star_end <- lapply(LVS, function(lv) {
                d <- R[R$level == lv & is.finite(R$fdr) & R$fdr < 0.05, ]
                c(lo = if (any(d$estimate < 0)) 0.17 else 0.05, hi = if (any(d$estimate >= 0)) 0.17 else 0.05)
            })
            TEXT_EST <- 0.34
            layout(matrix(1:2, 2), heights = vapply(1:2, function(i) diff(dr[[i]]) + TEXT_EST +
                                                        sum(star_end[[i]]), 0))
            par(cex = 1)
            par(omi = c(0.30 + note_in(note), 0.55, top_in(ttl) + 0.25, 0.30), cex = 1)
            for (li in 1:2) {
                lv <- LVS[li]
                d <- R[R$level == lv, ]
                par(mai = c(if (li == 2L) 0.10 else 0.05, 0.05, 0.05, 0.05))
                w <- 0.24
                xl <- c(0.45, length(ARMS) + 0.55)
                plot.new()
                cc <- Ff[Ff$level == lv & Ff$quantity == "arm_minus_hkg8", ]
                cc <- cc[order(match(cc$arm, ARMS)), ]
                txt <- if (nrow(cc)) fit_lines(paste("Module post-IP% - control:",
                                                     paste(sprintf("%s %+.2f [%+.2f, %+.2f]", cc$arm, cc$estimate,
                                                                   cc$lo, cc$hi), collapse = "; ")),
                                               par("pin")[1] - 0.2, T_NOTE, 1) else character(0)
                text_in <- if (length(txt)) 0.08 + 0.13 * length(txt) else 0.05
                se <- star_end[[li]]
                k <- diff(dr[[li]]) / (par("pin")[2] - text_in - se[["lo"]] - se[["hi"]])  # units per inch
                ylim <- c(dr[[li]][1] - se[["lo"]] * k, dr[[li]][2] + (text_in + se[["hi"]]) * k)
                plot.window(xlim = xl, ylim = ylim, yaxs = "i")
                yt <- pretty(dr[[li]]); yt <- yt[yt >= ylim[1] + 1e-9 & yt <= ylim[2] - (text_in * k)]
                axis(2, at = yt, las = 1, cex.axis = T_AXIS, mgp = c(3, 0.45, 0), tcl = -0.25)
                abline(h = 0, col = "black", lwd = 0.8)
                for (i in seq_along(ARMS)) for (j in seq_along(FQ)) {
                    r <- d[d$arm == ARMS[i] & d$quantity == names(FQ)[j], ]
                    if (!nrow(r) || !is.finite(r$estimate)) next
                    xx <- i + (j - 2) * w
                    rect(xx - w * 0.45, 0, xx + w * 0.45, r$estimate, col = FC[[j]], border = "black",
                         lwd = 0.7)
                    whisker(r$lo, r$hi, xx, 0.05, horiz = FALSE)
                    s <- stars(r$fdr)
                    if (nzchar(s)) text(xx, if (r$estimate < 0) r$lo else r$hi, s,
                                        pos = if (r$estimate < 0) 1 else 3, offset = 0.15, cex = T_VAL,
                                        font = 2)
                }
                for (i in seq_along(ARMS))
                    if (!any(is.finite(d$estimate[d$arm == ARMS[i]])))
                        text(i, 0.05, "untestable", srt = 90, adj = c(0, 0.5), cex = T_NOTE, col = MID)
                u <- par("usr")
                if (length(txt))
                    text(mean(u[1:2]), u[4] - 0.05 * k, paste(txt, collapse = "\n"),
                         adj = c(0.5, 1), cex = T_NOTE, font = 2)
                mtext(c("Between donors", "Within donor")[li], side = 4, line = 0.35, font = 2,
                      cex = T_LAB + 0.05)
                box()
                if (li == 2L) {
                    axis(1, at = seq_along(ARMS), labels = FALSE, tcl = -0.25)
                    mtext(ARMS, side = 1, at = seq_along(ARMS), line = 0.5, font = 2, cex = T_AXTITLE)
                }
            }
            mtext("rs with the library intronic mapping rate", side = 2, outer = TRUE, line = 1.7,
                  font = 2, cex = T_AXTITLE)
            frame_text(ttl, note, key = list(legend = FQ, fill = FC, border = "black"))
        })

        # == 16G, one file per arm ===========================================
        GQ <- c(module = "Module", mito = "%MT", ncells = "Log cell count", depth = "Median UMIs")
        ## 16G draws the accent purple as a 25% tint towards white.
        tint <- function(col, f) {
            m <- grDevices::col2rgb(col)[, 1]
            grDevices::rgb(t(round(m + (255 - m) * f)), maxColorValue = 255)
        }
        ACCENT_G <- tint(ACCENT, 0.25)
        GCOL <- c(module = ACCENT_G, mito = "#7F7F7F", ncells = "#BFBFBF", depth = INK)
        GPCH <- c(module = 21, mito = 22, ncells = 24, depth = 23)
        CQ <- c(module_minus_depth = "Module - depth", module_minus_mito = "Module - %MT",
                arm_minus_hkg8 = "Module - control")
        CCOL <- c(module_minus_depth = INK, module_minus_mito = "#7F7F7F", arm_minus_hkg8 = ACCENT_G)
        fields <- unique(G$measure)
        FL <- c(intronic = "Intronic\nrate", antisense = "Anti-\nsense", neun_fans = "NeuN\nfraction",
                cdna_400 = "cDNA\n>400 bp", pcr_cycles = "PCR\ncycles", saturation = "Satu-\nration")
        LAB_Q <- LAB_F <- NA_real_   # label widths in inches, measured on each 16G device
        forest_block <- function(g, qs, labs, cols, pchs, xl, show_lab) {
            per <- length(qs); gap <- 0.6
            ypos <- function(pi, qi) (length(fields) - pi) * (per + gap) + (per - qi) + 1
            ytop <- ypos(1, 1) + 0.7
            plot(NA, xlim = xl, ylim = c(0.3, ytop), xaxt = "n", yaxt = "n", xlab = "", ylab = "",
                 xaxs = "i", yaxs = "i")
            for (pi in seq_along(fields)) if (pi %% 2L == 0L)
                rect(xl[1], ypos(pi, per) - 0.5 - gap / 2, xl[2], ypos(pi, 1) + 0.5 + gap / 2, col = BAND,
                     border = NA)
            abline(v = 0, col = ZERO)
            ## Ticks at least 6% inside the panel edges.
            tk <- pretty(xl, 5); tk <- tk[tk > xl[1] + 0.06 * diff(xl) & tk < xl[2] - 0.06 * diff(xl)]
            axis(1, at = tk, cex.axis = T_AXIS, mgp = c(3, 0.4, 0), tcl = -0.25)
            for (pi in seq_along(fields)) {
                yy <- vapply(seq_len(per), function(qi) ypos(pi, qi), 0)
                for (qi in seq_len(per)) {
                    q <- qs[qi]
                    r <- g[g$measure == fields[pi] & g$quantity == q, ]
                    if (!nrow(r) || !is.finite(r$estimate)) next
                    if (is.finite(r$lo)) segments(r$lo, yy[qi], r$hi, yy[qi], lwd = WH_LWD, col = "black")
                    sig <- is.finite(r$fdr) && r$fdr < 0.05
                    ## Markers that do not pass FDR < 0.05 are drawn open and 15% smaller.
                    points(r$estimate, yy[qi], pch = pchs[[q]], cex = if (sig) PT else 0.85 * PT,
                           lwd = PT_LWD, col = "black", bg = if (sig) cols[[q]] else "white")
                    s <- stars(r$fdr)
                    if (nzchar(s)) text(min(if (is.finite(r$hi)) r$hi else r$estimate, xl[2]), yy[qi], s,
                                        pos = 4, offset = 0.15, cex = T_VAL, font = 2)
                }
                if (show_lab) {
                    u <- par("usr"); gx <- diff(grconvertX(c(0, 1), "inches", "user"))
                    text(u[1] - 0.05 * gx, yy, labs, adj = c(1, 0.5), cex = T_LAB - 0.04, xpd = NA)
                    segments(u[1] - (LAB_Q + 0.08) * gx, min(yy) - 0.3, u[1] - (LAB_Q + 0.08) * gx,
                             max(yy) + 0.3, xpd = NA)
                    text(u[1] - (LAB_Q + 0.14) * gx, mean(range(yy)), FL[[fields[pi]]], adj = c(1, 0.5),
                         cex = T_LAB, font = 2, xpd = NA)
                }
            }
            box()
        }
        for (a in ARMS) {
            ## HKG8 has no arm - HKG8 contrast row.
            cq <- if (a == "HKG8") CQ[names(CQ) != "arm_minus_hkg8"] else CQ
            ## Device size: T0 and B0 are the rs and contrast panel-row heights; the rs
            ## row is drawn 20% taller.
            nf0 <- length(fields)
            h1 <- nf0 * (length(GQ) + 0.6); h2 <- nf0 * (length(cq) + 0.6)
            L0 <- (2.4 + 1.05 * nf0) - 1.22
            T0 <- L0 * h1 / (h1 + h2); B0 <- L0 * h2 / (h1 + h2)
            fig_dev(sprintf("Fig16G_library_fields_%s", a), 5.78, 1.59 + 1.2 * T0 + B0, function() {
                ttl <- sprintf("Do sort purity or library construction order Module post-IP%%? (%s)", a)
                note <- paste("Spearman rs with the field (top) and differences between two rs values",
                              "on the same donor resamples (bottom); bars 95% donor-bootstrap intervals;",
                              "filled marker and * FDR < 0.05, ** FDR < 0.01, BH within arm and level.",
                              if (a == "HKG8") "" else "Module - control is the arm minus HKG8.")
                g0 <- G[G$arm == a, ]
                pad <- function(r) r + c(-0.05, 0.12) * diff(r)
                xl1 <- pad(range(c(g0$lo[g0$quantity %in% names(GQ)], g0$hi[g0$quantity %in% names(GQ)],
                                   -0.6, 0.6), na.rm = TRUE))
                xl2 <- pad(range(c(g0$lo[g0$quantity %in% names(cq)], g0$hi[g0$quantity %in% names(cq)],
                                   -0.5, 0.5), na.rm = TRUE))
                layout(matrix(1:4, 2, byrow = TRUE), heights = c(1.2 * T0, B0))
                ## layout() sets par(cex = 0.83); reset before measuring label widths.
                par(cex = 1)
                LAB_Q <<- max(strwidth(c(GQ, CQ), units = "inches", cex = T_LAB - 0.04)) + 0.08
                LAB_F <<- max(unlist(lapply(strsplit(FL, "\n"), strwidth, units = "inches", cex = T_LAB,
                                            font = 2)))
                par(omi = c(0.25 + note_in(note), LAB_Q + LAB_F + 0.30, top_in(ttl) + 0.30, 0.12), cex = 1)
                for (row in 1:2) for (li in 1:2) {
                    lv <- c("between donors", "within donor")[li]
                    g <- g0[g0$level == lv, ]
                    par(mai = c(0.42, 0.05, if (row == 1L) 0.22 else 0.38, 0.12))
                    if (row == 1L) forest_block(g, names(GQ), GQ, GCOL, GPCH, xl1, li == 1L)
                    else forest_block(g, names(cq), cq, CCOL, c(module_minus_depth = 21,
                                                               module_minus_mito = 22, arm_minus_hkg8 = 21),
                                      xl2, li == 1L)
                    mtext(if (row == 1L) c("Between donors", "Within donor")[li] else "Contrasts",
                          side = 3, line = 0.2, font = 2, cex = T_AXTITLE, col = if (row == 2L) MID else "black")
                    mtext(if (row == 1L) "rs against the field" else "Difference in rs",
                          side = 1, line = 1.4, font = 2, cex = T_AXTITLE - 0.05)
                    if (!any(is.finite(g$estimate)))
                        text(mean(par("usr")[1:2]), mean(par("usr")[3:4]), "untestable: fewer than 10 donors",
                             col = MID, cex = T_LAB, font = 3)
                }
                frame_text(ttl, note, key = list(legend = GQ, pch = GPCH, pt.bg = GCOL, col = "black",
                                                 pt.cex = 1.1))
            })
        }
        p12_say("")
    }

    p12_run(main)
}

# =============================================================================
# RUN
# =============================================================================
local({
    v <- p12_section_18()
    if (!identical(v, "PASS"))
        stop(sprintf("section 18 ended with verdict %s. See Part12_Fig16_log.txt.", v), call. = FALSE)
    message("Figure 16: section 18 PASS. Figures in Part12_rerun/08_Fig16/Fig16")
})
