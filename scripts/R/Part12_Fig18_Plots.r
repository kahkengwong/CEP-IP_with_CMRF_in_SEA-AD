# =============================================================================
# Part12_Fig18_Plots.r
# SEA-AD CEP-IP/CMRF, Part 12, Figure 18 plots (Figure 18, Supplementary Figure 17)
#   Draws the Figure 18 and Supplementary Figure 17 panels from the stats and
#   plot-data files written by sections 22 to 24 of Part12_Fig18_Analysis.r,
#   after checking them against their logs. Computes no statistic.
#
# Usage (from the SEA-AD folder, or set SEAAD_ROOT), after Part12_Fig18_Analysis.r
#   source("Part12_rerun/Part12_Fig18_Plots.r")
#
# Outputs (Part12_rerun/10_Fig18/Fig18/): figures (PDF, PNG), log.
# =============================================================================

source(file.path(Sys.getenv("SEAAD_ROOT", getwd()), "Part12_rerun", "Part12_03_Config.r"))

SCRIPT <- "Part12_Fig18_Plots.r"
ANALYSIS <- "Part12_Fig18_Analysis.r"
FIG_DIR <- "10_Fig18"   # every Figure 18 section folder sits inside this one
if (!isTRUE(as.integer(P12_CFG$config_version) >= 16L))
    stop(sprintf("config is version %s; Figure 18 needs version 16 or later",
                 P12_CFG$config_version), call. = FALSE)
F18 <- P12_CFG$figure18
PW18 <- F18$power

# == Run from Part12_rerun =====================================================
local({
    src <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
    want <- file.path(P12_RERUN, SCRIPT)
    if (!is.null(src) && !identical(tolower(normalizePath(src, winslash = "/", mustWork = FALSE)),
                                    tolower(normalizePath(want, winslash = "/", mustWork = FALSE))))
        stop(sprintf("run this file from Part12_rerun: save it there, then source(\"%s\")", want),
             call. = FALSE)
})

# == Helpers copied from Part12_Fig18_Analysis.r ==============================
VK18 <- c("second only", "Both", "NE panel only", "Neither", "Unreadable", "Untestable")
p12_18_verdict <- function(fdr_s, fdr_n, vif, fitted) {
    alpha <- as.numeric(P12_CFG$inference$alpha); vr <- as.numeric(P12_CFG$floors$vif_readable)
    ss <- is.finite(fdr_s) & fdr_s < alpha; sn <- is.finite(fdr_n) & fdr_n < alpha
    ifelse(!fitted, "Untestable",
           ifelse(!(is.finite(vif) & vif < vr), "Unreadable",
                  ifelse(ss & !sn, "second only",
                         ifelse(!ss & sn, "NE panel only",
                                ifelse(ss & sn, "Both", "Neither")))))
}

## Named assignments of an R file as parsed language, for exact comparison.
p12_18_defs <- function(path, what) {
    out <- list()
    if (!file.exists(path)) return(out)
    for (e in parse(path, keep.source = FALSE))
        if (is.call(e) && identical(e[[1L]], as.name("<-")) && is.name(e[[2L]]) &&
            as.character(e[[2L]]) %in% what) out[[as.character(e[[2L]])]] <- e[[3L]]
    out
}


# =============================================================================
# SECTION 25: Figure 18 and Supplementary Figure 17
# Panels 18A to 18E and S17, from the stats and plot files of sections 22 to
# 24. Checks first that the 18A verdicts, 18C medians and 18D crossings follow
# from those files.
# =============================================================================
p12_section_25 <- function() {
    STEP <- "25"
    D <- p12_step_dirs(STEP, "Fig18", parent = FIG_DIR)
    p12_clean(D$out)
    p12_start(STEP, SCRIPT, D, pkgs = c("ggplot2", "gtable"))

    PAL <- P12_CFG$palette; C18 <- F18$colours
    ARMS <- P12_ARMS; ARM_COL <- unlist(PAL$arms)
    INK <- "#3B3B3B"; MID <- "#7F7F7F"; PALE <- "#BFBFBF"; GRID <- C18$grid
    HEAD <- PAL$header
    VFILL <- c("second only" = "#2D2D2D", "Both" = "#7B7B7B", "NE panel only" = C18$ne_only,
               "Neither" = "#D9D9D9", "Unreadable" = "#FBE7A2", "Untestable" = "#F2F2F2")
    DETECT <- as.numeric(PW18$detect); PLAN <- as.numeric(PW18$plan_effect)
    EPS <- as.numeric(F18$eps_mt)
    ARM_LABEL <- c(M4 = "M4 (primary)", PVM_3 = "PVM_3 (reference)", PVM_1 = "PVM_1", HKG8 = "HKG8 (control)")
    effx <- function(e) if (abs(e * 10 - round(e * 10)) > 1e-9) sprintf("%.2fx", e) else sprintf("%.1fx", e)

    fig_dev <- function(stem, w, h, draw) {
        ok <- TRUE
        for (kind in c("pdf", "png")) {
            f <- file.path(D$out, paste0(stem, ".", kind))
            p12_output(f)
            if (kind == "pdf") grDevices::pdf(f, width = w, height = h, useDingbats = FALSE)
            else grDevices::png(f, width = w, height = h, units = "in", res = 300,
                                type = if (.Platform$OS.type == "windows") "windows" else "cairo")
            r <- tryCatch({ draw(); TRUE }, error = function(e) {
                p12_say(sprintf("    %s %s: %s", stem, kind, conditionMessage(e))); FALSE })
            grDevices::dev.off()
            ok <- ok && r
        }
        p12_chk(ok, sprintf("%s drawn", stem), "pdf and png")
    }
    rd <- function(folder, file) utils::read.csv(p12_upstream_file(file.path(FIG_DIR, folder), file),
                                                 check.names = FALSE, stringsAsFactors = FALSE)
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
    ## Label colour: white or dark grey by WCAG contrast.
    txt_col <- function(h) {
        f <- grDevices::col2rgb(h)[, 1] / 255
        lin <- ifelse(f <= 0.04045, f / 12.92, ((f + 0.055) / 1.055)^2.4)
        lum <- 0.2126 * lin[1] + 0.7152 * lin[2] + 0.0722 * lin[3]
        if (1.05 / (lum + 0.05) >= (lum + 0.05) / 0.0685) "white" else "grey15"
    }
    overlay <- function() {
        par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE, xpd = NA)
        plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "",
             xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
    }

    main <- function() {
        p12_say("SECTION A -- inputs")
        ## Stats must come from the current analysis file, with the same verdict rule.
        for (fd in c("18A_NuclearPanels", "18BC_MTOverlap", "18DE_DonorPower")) {
            lg <- file.path(P12_RERUN, FIG_DIR, fd, sprintf("Part12_%s_log.txt", fd))
            sl <- if (file.exists(lg)) sub("\r$", "", grep("^Script: ", readLines(lg, warn = FALSE),
                                                             value = TRUE)[1]) else NA_character_
            p12_chk(identical(sl, sprintf("Script: %s", ANALYSIS)),
                    sprintf("%s log was written by %s", fd, ANALYSIS),
                    if (is.na(sl)) "log or Script line not found" else sl)
        }
        sh <- c("VK18", "p12_18_verdict")
        da <- p12_18_defs(file.path(P12_RERUN, ANALYSIS), sh)
        dp <- p12_18_defs(file.path(P12_RERUN, SCRIPT), sh)
        same <- vapply(sh, function(n) !is.null(da[[n]]) && !is.null(dp[[n]]) &&
                           identical(da[[n]], dp[[n]]), TRUE)
        p12_chk(all(same), sprintf("shared helpers identical to %s", ANALYSIS),
                paste(sprintf("%s %s", sh, ifelse(same, "same", "DIFFERS or not found")),
                      collapse = " | "))
        S18A <- rd("18A_NuclearPanels", "Stats_18A.csv")
        V18A <- rd("18A_NuclearPanels", "Plot_18A_verdicts.csv")
        S18B <- rd("18BC_MTOverlap", "Stats_18B.csv")
        S18C <- rd("18BC_MTOverlap", "Stats_18C.csv")
        CELLS <- rd("18BC_MTOverlap", "Plot_18B_cells.csv")
        DON <- rd("18BC_MTOverlap", "Plot_18C_donors.csv")
        S18D <- rd("18DE_DonorPower", "Stats_18D.csv")
        CUR <- rd("18DE_DonorPower", "Plot_18D_curves.csv")
        REG <- rd("18DE_DonorPower", "Plot_18E_regions.csv")
        CELLS$is_module <- as.logical(CELLS$is_module)

        ## Closure checks.
        g <- function(q) S18A[S18A$quantity == q, ][match(
            paste(V18A$arm, V18A$dataset, paste(V18A$predictor, V18A$panel, V18A$outcome, sep = "|")),
            with(S18A[S18A$quantity == q, ], paste(arm, level, measure))), ]
        v2 <- p12_18_verdict(g("beta_second")$fdr, g("beta_panel")$fdr, g("vif")$estimate,
                             is.finite(g("beta_second")$estimate) & is.finite(g("beta_panel")$estimate))
        p12_chk(identical(v2, V18A$verdict), "18A verdicts follow from Stats_18A", sprintf("%d rows", nrow(V18A)))
        med <- vapply(seq_len(nrow(S18C)), function(i) {
            v <- DON$auc[DON$arm == S18C$arm[i] & DON$dataset == S18C$level[i]]
            if (length(v) >= as.integer(P12_CFG$floors$donors_per_test)) stats::median(v) else NA_real_ }, 0)
        p12_chk(all((is.na(med) & is.na(S18C$estimate)) | abs(med - S18C$estimate) < 1e-12, na.rm = TRUE) &&
                    all(is.na(med) == is.na(S18C$estimate)),
                "18C medians follow from the donor AUCs", sprintf("%d rows", nrow(S18C)))
        cr <- S18D[S18D$quantity == "n_at_80", ]
        cr2 <- vapply(seq_len(nrow(cr)), function(i) {
            d <- CUR[CUR$arm == cr$arm[i] & sprintf("%.2fx", CUR$effect) == cr$level[i], ]
            d <- d[order(d$n_donors), ]
            col <- c(R1 = "r1_rate", R2 = "r2_rate", R3 = "r3_rate")[[cr$measure[i]]]
            ok <- is.finite(d[[col]]) & d[[col]] >= DETECT
            if (any(ok)) d$n_donors[which(ok)[1]] else NA_real_ }, 0)
        p12_chk(identical(is.na(cr2), is.na(cr$estimate)) && all(cr2 == cr$estimate, na.rm = TRUE),
                "18D crossings follow from the curves", sprintf("%d crossings", nrow(cr)))
        labs <- c(vapply(F18$panels, function(p) p$label, ""), names(VFILL), ARM_LABEL, ARMS)
        p12_chk(p12_ascii_ok(labs), "figure labels are ASCII", sprintf("%d labels", length(labs)))
        ## Figure text is typed in this file, so the whole file is checked.
        src <- readLines(file.path(P12_RERUN, SCRIPT), warn = FALSE)
        p12_chk(p12_ascii_ok(src), "figure script text is ASCII", sprintf("%d lines", length(src)))
        if (P12$n[["FAIL"]] > 0L) p12_abort("checks failed; figures not drawn")
        p12_say("")
        p12_say("SECTION B -- figures")

        # == 18A: the two duel figures ==============================================
        PAN <- F18$panels
        PK <- vapply(PAN, function(p) p$key, ""); PL <- vapply(PAN, function(p) p$label, "")
        CLS <- c(cognition = "Cognition", pathology = "Pathology", genotype = "APOE")
        OCL <- vapply(P12_CFG$clinical$outcomes, function(o) o$class, "")
        ncl <- vapply(names(CLS), function(k) sum(OCL == k), 0)
        GAP <- 0.7
        ## 18A: one bar per panel row, one region set per file (NeoC main,
        ## NonNeoC supplementary); panel height 3.352 in, then 10% taller.
        gg18a <- function(sp, ds) {
            d <- V18A[V18A$predictor == sp & V18A$dataset == ds, ]
            own <- if (sp == "mito") "%MT only" else "Module only"
            lab <- c("second only" = own, "Both" = "Both", "NE panel only" = "NE panel only",
                     "Neither" = "Neither", "Unreadable" = "Unreadable", "Untestable" = "Untestable")
            pred <- if (sp == "mito") "%MT" else "Module post-IP%"
            starts <- cumsum(c(0, ncl[-length(ncl)])) + (seq_along(ncl) - 1) * GAP
            rows <- list()
            for (a in ARMS) for (r in seq_along(PK)) for (ci in seq_along(CLS)) {
                z <- d[d$arm == a & d$panel == PK[r] & d$class == names(CLS)[ci], ]
                tb <- table(factor(z$verdict, VK18)); x0 <- starts[ci]
                for (k in VK18) {
                    n <- tb[[k]]; if (!n) next
                    rows[[length(rows) + 1L]] <- data.frame(arm = a, ypos = length(PK) - r + 1, xmin = x0,
                                                            xmax = x0 + n, n = n, verdict = k,
                                                            stringsAsFactors = FALSE)
                    x0 <- x0 + n
                }
            }
            tal <- do.call(rbind, rows)
            tal$arm <- factor(tal$arm, levels = ARMS)
            tal$verdict <- factor(tal$verdict, levels = names(VFILL))
            tal$xmid <- (tal$xmin + tal$xmax) / 2
            tal$txt_col <- vapply(as.character(tal$verdict), function(k) txt_col(VFILL[[k]]), "")
            ytop <- length(PK) + 0.75
            blk <- data.frame(xstart = starts, xend = starts + ncl, class = unname(CLS))
            shade <- blk[seq_len(nrow(blk)) %% 2L == 0L, , drop = FALSE]
            clab <- data.frame(xmid = (blk$xstart + blk$xend) / 2, class = blk$class,
                               arm = factor(ARMS[1], levels = ARMS))
            rn <- d[d$arm != "HKG8", ]; hk <- d[d$arm == "HKG8", ]
            other <- setdiff(P12_DATASETS, ds)
            wrap <- function(z, w) paste(strwrap(z, width = w), collapse = "\n")
            sub <- sprintf(paste("outcome ~ cov4 + nuclear panel + %s, both predictors standardised, FDR within",
                                 "arm, region set and panel over the outcomes. %s only: each bar sums to the %d",
                                 "outcomes, ordered by outcome class. %s is in the %s figure."),
                           pred, ds, length(OCL), other,
                           if (ds == P12_DATASETS[1]) "supplementary" else "main")
            cap <- c(sprintf("Segments run left to right: %s. The classes are set apart by a gap, so there is no x axis; the numerals are the counts.",
                             paste(sprintf("%s (%d)", CLS, ncl), collapse = ", ")),
                     "The matched reference panel is the null: a nuclear panel that wins no more rows than it has shown nothing.",
                     sprintf("VIF of the two predictors: at most %.2f (ceiling %g). Real arms: NE panel only on %d of %d rows; control: %d of %d.",
                             suppressWarnings(max(d$vif, na.rm = TRUE)), as.numeric(P12_CFG$floors$vif_readable),
                             sum(rn$verdict == "NE panel only"), nrow(rn), sum(hk$verdict == "NE panel only"),
                             nrow(hk)))
            p <- ggplot2::ggplot(tal) +
                ggplot2::geom_rect(data = shade, ggplot2::aes(xmin = xstart - GAP / 2, xmax = xend + GAP / 2,
                                                              ymin = 0.4, ymax = ytop + 0.7),
                                   inherit.aes = FALSE, fill = C18$class_shade, colour = NA) +
                ggplot2::geom_rect(ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ypos - 0.34, ymax = ypos + 0.34,
                                                fill = verdict), colour = "black", linewidth = 0.3) +
                ## Counts at size 3.5.
                ggplot2::geom_text(ggplot2::aes(x = xmid, y = ypos, label = n, colour = txt_col), size = 3.5) +
                ggplot2::geom_text(data = clab, ggplot2::aes(x = xmid, y = ytop, label = class),
                                   inherit.aes = FALSE, size = 2.5, colour = C18$class_label, fontface = "bold",
                                   vjust = 0) +
                ggplot2::scale_colour_identity() +
                ggplot2::scale_fill_manual(values = VFILL, breaks = names(VFILL), labels = unname(lab[names(VFILL)]),
                                           name = NULL, drop = FALSE) +
                ggplot2::scale_x_continuous(breaks = NULL, limits = c(NA, max(starts + ncl)),
                                            expand = ggplot2::expansion(mult = c(0.01, 0.07))) +
                ggplot2::scale_y_continuous(breaks = seq_along(PK), labels = rev(PL), limits = c(0.4, ytop + 0.7),
                                            expand = ggplot2::expansion(add = 0)) +
                ggplot2::facet_grid(. ~ arm) +
                ggplot2::labs(title = sprintf("%s vs nuclear-encoded panels in one model (cov4; clinical; %s)", pred, ds),
                              subtitle = wrap(sub, round(13 * 9.5)),
                              caption = paste(vapply(cap, wrap, "", w = round(11 * 9.5)), collapse = "\n"),
                              x = NULL, y = NULL) +
                ggplot2::theme_bw(base_size = 10) +
                ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 12),
                               plot.subtitle = ggplot2::element_text(size = 8, colour = "grey30"),
                               plot.caption = ggplot2::element_text(size = 7.5, colour = "grey30", hjust = 0),
                               strip.background = ggplot2::element_rect(fill = HEAD),
                               strip.text = ggplot2::element_text(colour = "white", face = "bold", size = 9),
                               panel.grid.minor = ggplot2::element_blank(),
                               legend.position = "bottom",
                               panel.spacing = grid::unit(14, "pt"),
                               plot.background = ggplot2::element_rect(fill = "white", colour = NA))
            ## Panel height, and the device height that holds the rest.
            grDevices::pdf(NULL, width = 9.5, height = 7.7)
            g <- ggplot2::ggplotGrob(p)
            pr <- unique(g$layout$t[grepl("^panel", g$layout$name)])
            g$heights[pr] <- grid::unit(3.352, "in")
            h <- sum(grid::convertHeight(g$heights, "in", valueOnly = TRUE))
            ## 10% taller; the extra height goes to the panels.
            g$heights[pr] <- grid::unit(3.352 + 0.10 * h, "in")
            h <- 1.10 * h
            grDevices::dev.off()
            list(g = g, h = h)
        }
        for (ds in P12_DATASETS) for (sp in c("module", "mito")) {
            stem <- if (ds == P12_DATASETS[1]) sprintf("Fig18A_duel_%s", sp)
                    else sprintf("FigS18%s_duel_%s_%s", if (sp == "module") "A" else "B", sp, ds)
            gg <- gg18a(sp, ds)
            fig_dev(stem, 9.5, gg$h, function() grid::grid.draw(gg$g))
        }

        # == 18B: densities ===========================================================
        TCL <- -0.25; MGP <- c(2.2, 0.55, 0); CEX_AXIS <- 0.95; CEX_MAIN <- 1.15
        COL_MOD <- C18$module_cells; COL_NON <- C18$nonmodule_cells
        dens_pair <- function(a, ds) {
            d <- CELLS[CELLS$arm == a & CELLS$dataset == ds, ]
            x1 <- d$mito_pct[d$is_module]; x0 <- d$mito_pct[!d$is_module]
            if (length(x1) < 5L || length(x0) < 5L) {
                plot.new(); title(main = sprintf("%s  %s", a, ds), col.main = INK); return(invisible(NULL))
            }
            l1 <- log10(x1 + EPS); l0 <- log10(x0 + EPS)
            lo <- stats::quantile(c(l1, l0), 0.002); hi <- stats::quantile(c(l1, l0), 0.999)
            d1 <- stats::density(l1, from = lo, to = hi); d0 <- stats::density(l0, from = lo, to = hi)
            ymax <- max(c(d1$y, d0$y)) * 1.10
            plot(NA, xlim = c(lo, hi), ylim = c(0, ymax), xaxt = "n", yaxt = "n",
                 xlab = "%MT per cell (log scale)", ylab = "density")
            tk <- c(0, 0.1, 0.3, 1, 3, 10, 30)
            tk <- tk[log10(tk + EPS) >= lo & log10(tk + EPS) <= hi]
            abline(v = log10(tk + EPS), col = GRID, lwd = 1.1)
            abline(h = pretty(c(0, ymax), 5), col = GRID, lwd = 1.1)
            polygon(c(d0$x, rev(d0$x)), c(d0$y, rep(0, length(d0$y))), col = paste0(COL_NON, "59"),
                    border = COL_NON, lwd = 1.8)
            polygon(c(d1$x, rev(d1$x)), c(d1$y, rep(0, length(d1$y))), col = paste0(COL_MOD, "45"),
                    border = COL_MOD, lwd = 1.8)
            s <- function(q) S18B$estimate[S18B$arm == a & S18B$level == ds & S18B$quantity == q]
            if (is.finite(s("best_cut"))) abline(v = log10(s("best_cut") + EPS), lty = 2, col = INK, lwd = 1.5)
            mtext(sprintf("AUC %.3f   kappa %.3f   best cut %.2f%%", s("auc"), s("kappa"), s("best_cut")),
                  side = 3, line = 0.2, cex = 0.7, col = MID)
            axis(1, at = log10(tk + EPS), labels = tk, tcl = TCL, mgp = MGP)
            axis(2, at = pretty(c(0, ymax), 5), las = 1, tcl = TCL, mgp = MGP)
            box(col = INK)
            title(main = sprintf("%s  %s", a, ds), cex.main = CEX_MAIN, col.main = INK, line = 1.9)
        }
        ## Main figure NeoC, supplementary NonNeoC; one row of panels each.
        dens_fig <- function(ds) function() {
            ## Bottom outer margin 3.4 lines, so the key sits below the x-axis titles.
            par(mfrow = c(1, length(ARMS)), mar = c(4.2, 4.2, 3.6, 0.8), oma = c(3.4, 0.4, 5.4, 0.4),
                col.axis = INK, col.lab = INK, fg = INK, tcl = TCL, mgp = MGP, cex.axis = CEX_AXIS)
            for (a in ARMS) dens_pair(a, ds)
            sel <- S18B$arm != "HKG8" & S18B$level == ds
            kr <- range(S18B$estimate[S18B$quantity == "kappa" & sel], na.rm = TRUE)
            cr <- range(S18B$estimate[S18B$quantity == "best_cut" & sel], na.rm = TRUE)
            ar <- range(S18B$estimate[S18B$quantity == "auc" & sel], na.rm = TRUE)
            mtext(sprintf("Does any %%MT cut reproduce the Module cells? (%s)", ds), outer = TRUE, side = 3,
                  line = 3.6, font = 2, cex = 0.78)
            mtext(sprintf(paste("Real arms: AUC %.2f to %.2f, best-cut kappa %.2f to %.2f; the best cut falls at",
                                "%.2f to %.2f per cent %%MT."), ar[1], ar[2], kr[1], kr[2], cr[1], cr[2]),
                  outer = TRUE, side = 3, line = 2.3, cex = 0.72, col = MID)
            mtext(paste("Cells pooled over donors (descriptive); the within-donor AUC in 18C is the test.",
                        "Dashed line: the single best %MT cut (Youden's J).",
                        if (ds == P12_DATASETS[1]) sprintf("%s is in the supplementary figure.", P12_DATASETS[2])
                        else sprintf("%s is in the main figure.", P12_DATASETS[1])),
                  outer = TRUE, side = 3, line = 1.2, cex = 0.72, col = MID)
            par(fig = c(0, 1, 0, 1), mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0), new = TRUE)
            plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
            legend("bottom", horiz = TRUE, bty = "n", cex = 0.8,
                   fill = c(paste0(COL_MOD, "45"), paste0(COL_NON, "59")), border = c(COL_MOD, COL_NON),
                   legend = c("Module cells", "non-Module cells"))
        }
        ## Panel plot region 2.564 x 2.521 in; one margin line is 0.132 in.
        H18B <- (8 - (1.4 + 5.4) * 0.132) / 2 + (3.4 + 5.4) * 0.132
        fig_dev("Fig18B_mito_densities", 13, H18B, dens_fig(P12_DATASETS[1]))
        fig_dev("FigS18C_mito_densities_NonNeoC", 13, H18B, dens_fig(P12_DATASETS[2]))

        # == 18C: within-donor AUC ====================================================
        fig_dev("Fig18C_within_donor_auc", 8, 5.6, function() {
            par(mfrow = c(1, length(P12_DATASETS)), mar = c(6.2, 4.4, 3.6, 0.8), oma = c(0.4, 0.4, 4.2, 0.4),
                col.axis = INK, col.lab = INK, fg = INK, tcl = TCL, mgp = MGP, cex.axis = CEX_AXIS, las = 1)
            set.seed(p12_seed("23"))
            for (ds in P12_DATASETS) {
                plot(NA, xlim = c(0.70, length(ARMS) + 0.30), ylim = c(0, 1.18), xaxt = "n", yaxt = "n",
                     xlab = "", ylab = "AUC of %MT, within donor")
                abline(h = seq(0, 1, 0.2), col = GRID, lwd = 1.1)
                abline(v = seq_along(ARMS), col = GRID, lwd = 1.1)
                abline(h = 0.5, lty = 2, col = PALE, lwd = 1.8)
                axis(2, at = seq(0, 1, 0.2), las = 1, tcl = TCL, mgp = MGP)
                box(col = INK)
                for (j in seq_along(ARMS)) {
                    v <- DON$auc[DON$arm == ARMS[j] & DON$dataset == ds]; v <- v[is.finite(v)]
                    r <- S18C[S18C$arm == ARMS[j] & S18C$level == ds, ]
                    if (length(v)) {
                        points(jitter(rep(j, length(v)), amount = 0.16), v, pch = 21,
                               bg = paste0(unname(ARM_COL[ARMS[j]]), "80"), col = unname(ARM_COL[ARMS[j]]),
                               cex = 1.15, lwd = 0.8)
                        if (length(v) >= 2L)
                            boxplot(v, at = j, add = TRUE, boxwex = 0.52, axes = FALSE,
                                    col = paste0(unname(ARM_COL[ARMS[j]]), "1A"), border = "#000000",
                                    outline = FALSE, lwd = 2.1, medlwd = 3.0, whisklty = 1, staplewex = 0.42)
                    }
                    text(j, 0.042, sprintf("n=%d", length(v)), cex = 0.85, col = MID)
                    text(j, 1.11, if (nrow(r) && is.finite(r$fdr))
                        sprintf("%.3f%s", r$fdr, p12_stars(r$fdr)) else "untestable",
                         cex = if (nrow(r) && is.finite(r$fdr)) 0.85 else 0.7,
                         col = if (nrow(r) && is.finite(r$fdr)) "#000000" else MID)
                }
                mtext("FDR=", side = 2, at = 1.11, las = 1, adj = 1, line = 0.3, cex = 0.85, col = "#000000")
                for (j in seq_along(ARMS)) mtext(ARMS[j], side = 1, at = j, line = 0.7, cex = 0.86, col = INK)
                title(main = ds, cex.main = CEX_MAIN, col.main = INK)
            }
            re <- S18C$estimate[S18C$arm != "HKG8"]; ce <- S18C$estimate[S18C$arm == "HKG8"]
            mtext("Does %MT distinguish Module cells within a donor?", outer = TRUE, side = 3, line = 2.4,
                  font = 2, cex = 0.8)
            mtext(sprintf("Median within-donor AUC %s on the real arms, %s on the control.",
                          if (any(is.finite(re))) sprintf("%.2f to %.2f", min(re, na.rm = TRUE),
                                                          max(re, na.rm = TRUE)) else "untestable",
                          if (any(is.finite(ce))) sprintf("%.2f to %.2f", min(ce, na.rm = TRUE),
                                                          max(ce, na.rm = TRUE)) else "untestable"),
                  outer = TRUE, side = 3, line = 1.3, cex = 0.7, col = MID)
            mtext(sprintf(paste("One point per donor with %d or more cells on each side; dashed line is chance.",
                                "FDR: donor bootstrap of the median against 0.5, BH within region set."),
                          as.integer(P12_CFG$floors$cells_per_group_per_donor)),
                  outer = TRUE, side = 3, line = 0.3, cex = 0.7, col = MID)
        })

        # == 18D and S17: power curves ================================================
        RULES <- c("r1_rate", "r2_rate", "r3_rate")
        RULE_TAG <- c(r1_rate = "R1", r2_rate = "R2", r3_rate = "R3")
        RULE_LAB <- c(r1_rate = "R1  donor share clears its null", r2_rate = "R2  split-half >= 0.50",
                      r3_rate = "R3  clinical association (MMSE, cov4)")
        RULE_COL <- c(r1_rate = C18$rule_r1, r2_rate = C18$rule_r2, r3_rate = C18$rule_r3)
        MGP2 <- c(2.4, 0.55, 0); CEX_MAIN2 <- 1.05; XTK <- c(10, 20, 40, 60, 80)
        XMAX <- max(CUR$n_donors, na.rm = TRUE) * 1.04
        cross_of <- function(a, e, rule) {
            r <- S18D$estimate[S18D$arm == a & S18D$level == sprintf("%.2fx", e) & S18D$quantity == "n_at_80" &
                                   S18D$measure == rule]
            if (length(r)) r[1] else NA_real_
        }
        ## `rules`: the rules drawn. 18D draws R3 only; Supplementary Figure 17 all three.
        curve_panel <- function(a, e, show_ylab, rules = RULES) {
            d <- CUR[CUR$arm == a & abs(CUR$effect - e) < 1e-9, ]; d <- d[order(d$n_donors), ]
            ## Crossings go in the note line under the title.
            xing <- vapply(rules, function(rr) {
                nx <- cross_of(a, e, RULE_TAG[[rr]])
                nn <- if (is.na(nx)) "n.r." else as.integer(nx)
                if (length(rules) == 1L) sprintf("%s donors", nn) else sprintf("%s %s", RULE_TAG[[rr]], nn) }, "")
            note <- if (nrow(d)) sprintf("target share %.2f; curves run to %d donors; %.0f%% at %s",
                                         d$target_share[1], as.integer(max(d$n_donors)), 100 * DETECT,
                                         paste(xing, collapse = ", ")) else "no rows"
            plot(NA, xlim = c(0, XMAX), ylim = c(0, 1.06), xaxt = "n", yaxt = "n", xlab = "donors",
                 ylab = if (show_ylab) "detection rate" else "")
            abline(h = seq(0, 1, 0.2), col = GRID, lwd = 1.1); abline(v = XTK, col = GRID, lwd = 1.1)
            abline(h = DETECT, lty = 2, col = PALE, lwd = 1.8)
            drew <- FALSE
            if (nrow(d)) {
                if ("r2_rate" %in% rules) {
                    gap <- suppressWarnings(max(abs(d$r2_evaluable - d$r2_rate), na.rm = TRUE))
                    drew <- is.finite(gap) && gap > 0.15
                    if (drew) lines(d$n_donors, d$r2_evaluable, lty = 3, lwd = 2.2, col = RULE_COL[["r2_rate"]])
                }
                for (rr in rules) {
                    ok <- is.finite(d[[rr]])
                    lines(d$n_donors[ok], d[[rr]][ok], col = RULE_COL[[rr]], lwd = 3.0)
                    points(d$n_donors[ok], d[[rr]][ok], pch = 16, col = RULE_COL[[rr]], cex = 1.35)
                    nx <- cross_of(a, e, RULE_TAG[[rr]])
                    if (is.finite(nx)) points(nx, d[[rr]][d$n_donors == nx][1], pch = 23, bg = C18$crossing_fill,
                                              col = "#000000", cex = 2.1, lwd = 0.9)
                }
            }
            axis(1, at = XTK, tcl = TCL, mgp = MGP2)
            axis(2, at = seq(0, 1, 0.2), las = 1, tcl = TCL, mgp = MGP2)
            box(col = INK)
            title(main = sprintf("%s  %s", ARM_LABEL[[a]], effx(e)), cex.main = CEX_MAIN2, col.main = INK,
                  line = 1.7)
            mtext(note, side = 3, line = 0.2, cex = 0.60, col = MID)
            drew
        }
        curve_figure <- function(arms, effs, headers, rules = RULES) function() {
            par(mfrow = c(1, 4), mar = c(4.2, 4.0, 4.0, 0.9), oma = c(3.0, 0.6, 5.0, 0.4), col.axis = INK,
                col.lab = INK, fg = INK, tcl = TCL, mgp = MGP2, cex.axis = CEX_AXIS)
            any_eval <- FALSE
            for (i in 1:4) {
                if (i > length(arms) || is.na(arms[i])) { plot.new(); next }
                any_eval <- curve_panel(arms[i], effs[i], i == 1L, rules) || any_eval
            }
            for (k in seq_along(headers))
                mtext(headers[[k]], outer = TRUE, side = 3, line = c(3.3, 2.1, 1.0)[k], font = if (k == 1L) 2 else 1,
                      cex = if (k == 1L) 0.82 else 0.72, col = if (k == 1L) "#000000" else MID)
            par(fig = c(0, 1, 0, 1), mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0), new = TRUE)
            plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
            lab <- if (length(rules) == 1L) sub("^R[0-9]  ", "", RULE_LAB[rules]) else RULE_LAB[rules]
            legend("bottom", horiz = TRUE, bty = "n", cex = 0.72, seg.len = 2.2,
                   lwd = c(rep(3, length(rules)), if (any_eval) 2.2), lty = c(rep(1, length(rules)), if (any_eval) 3),
                   col = c(unname(RULE_COL[rules]), if (any_eval) RULE_COL[["r2_rate"]]), text.col = INK,
                   legend = c(unname(lab), if (any_eval) "R2 evaluable (the test could be run)"))
        }
        fmt_n <- function(x) if (is.na(x)) "not reached" else sprintf("%d donors", as.integer(x))
        ## 18D: the clinical rule only; its crossing at the planning effect is the
        ## donor requirement.
        n_m4 <- cross_of("M4", PLAN, "R3"); n_p3 <- cross_of("PVM_3", PLAN, "R3")
        hdr_main <- list(
            sprintf("The clinical association (MMSE) sets the requirement: about %s on M4 and %s on PVM_3",
                    fmt_n(n_m4), if (is.na(n_p3)) "not reached" else as.integer(n_p3)),
            sprintf(paste("Detection of MMSE ~ cov4 + donor Module post-IP%% (OLS p < 0.05) against donor count,",
                          "at 1.0x and %s the observed effect. Donor share and split-half rules:",
                          "SuppTable_Fig18."), effx(PLAN)),
            sprintf(paste("Donors subsampled from this cohort, %d draws per point; the donor share rescaled to each",
                          "effect. Yellow diamond: first donor count at %.0f%% detection. Simulation, no test."),
                    as.integer(PW18$draws), 100 * DETECT))
        fig_dev("Fig18D_power_curves", 14, 4.5,
                curve_figure(c("M4", "M4", "PVM_3", "PVM_3"), c(1.0, PLAN, 1.0, PLAN), hdr_main,
                             rules = "r3_rate"))
        for (a in ARMS) {
            hdr <- list(sprintf("%s: donor-count curves at 1.0x, %s and 0.5x the observed effect", ARM_LABEL[[a]],
                                effx(PLAN)),
                        sprintf("At %s: R1 %s, R2 %s, R3 %s.", effx(PLAN), fmt_n(cross_of(a, PLAN, "R1")),
                                fmt_n(cross_of(a, PLAN, "R2")), fmt_n(cross_of(a, PLAN, "R3"))),
                        if (a == "HKG8") "This is the control: its R1 curve should not reach 80% at any donor count." else
                            sprintf("Curves run to %d of this arm's %d donors (the last grid point at or below it).",
                                    as.integer(max(CUR$n_donors[CUR$arm == a])),
                                    as.integer(REG$donors_all[REG$arm == a][1])))
            fig_dev(sprintf("FigS17_power_curves_%s", a), 14, 4.5,
                    curve_figure(c(a, a, a, NA), c(1.0, PLAN, 0.5, NA), hdr))
        }

        # == 18E: region check ========================================================
        REG_COL <- c(frac = C18$region_frac, share = C18$region_share, split = C18$region_split)
        REG_LAB <- c(frac = "donors retained", share = "donor share", split = "split-half rs")
        ## 4.785 x 7.5 in.
        fig_dev("Fig18E_region_check", 4.35 * 1.10, 7.5, function() {
            xr <- range(REG$min_regions, na.rm = TRUE); xt <- seq(xr[1], xr[2], by = 1)
            fl <- as.integer(P12_CFG$floors$donors_per_test)
            w1 <- strwrap("What does raising the region floor cost?", width = 54)
            rf <- REG[REG$arm == "M4", ]; rf <- rf[order(rf$min_regions), ]
            w2 <- unlist(lapply(c(
                if (nrow(rf) >= 2L) sprintf("On M4, 3+ keeps %d of %d donors at share %.2f; 4+ keeps %d.",
                                            as.integer(rf$donors[1]), as.integer(rf$donors_all[1]), rf$donor_share[1],
                                            as.integer(rf$donors[2])) else "",
                sprintf("Open points use fewer than %d donors. Dashed grey line is the 0.50 split-half rule.", fl),
                "Kept donors at a high minimum are those sampled in many regions, not a random subset."),
                strwrap, width = 66))
            w2 <- w2[nzchar(w2)]
            nH <- length(w1) + length(w2)
            par(mfrow = c(length(ARMS), 1), oma = c(6.4, 4.6, 1.4 + 1.15 * nH, 0.6), col.axis = INK,
                col.lab = INK, fg = INK, tcl = TCL, mgp = MGP2, cex.axis = CEX_AXIS)
            for (k in seq_along(ARMS)) {
                a <- ARMS[k]
                d <- REG[REG$arm == a, ]; d <- d[order(d$min_regions), ]
                par(mar = c(if (k == length(ARMS)) 0.9 else 0.5, 0.6, 1.9, 0.8))
                plot(NA, xlim = c(xr[1] - 0.3, xr[2] + 0.3), ylim = c(0, 1.06), xaxt = "n", yaxt = "n",
                     xlab = "", ylab = "")
                abline(h = seq(0, 1, 0.2), col = GRID, lwd = 1.1); abline(v = xt, col = GRID, lwd = 1.1)
                abline(h = 0.5, lty = 2, col = PALE, lwd = 1.8)
                thin <- d$donors < fl
                for (v in c("donors_frac", "donor_share", "split_half")) {
                    nm <- c(donors_frac = "frac", donor_share = "share", split_half = "split")[[v]]
                    ok <- is.finite(d[[v]]); if (!any(ok)) next
                    lines(d$min_regions[ok], d[[v]][ok], col = REG_COL[[nm]], lwd = 3.0)
                    sel <- ok & !thin
                    if (any(sel)) points(d$min_regions[sel], d[[v]][sel], pch = 16, col = REG_COL[[nm]], cex = 1.35)
                    sel <- ok & thin
                    if (any(sel)) points(d$min_regions[sel], d[[v]][sel], pch = 21, bg = "#FFFFFF",
                                         col = REG_COL[[nm]], cex = 1.25, lwd = 2.0)
                }
                if (k == length(ARMS)) axis(1, at = xt, tcl = TCL, mgp = MGP2)
                axis(2, at = seq(0, 1, 0.5), las = 1, tcl = TCL, mgp = MGP2)
                box(col = INK)
                ns <- sum(is.finite(d$split_half))
                sub <- if (!nrow(d)) "" else if (ns == 0L) "split-half never runs" else
                    if (ns < nrow(d)) sprintf("split-half stops at %d+",
                                              as.integer(max(d$min_regions[is.finite(d$split_half)])) + 1L) else ""
                mtext(sprintf("%s  (%d donors)", ARM_LABEL[[a]], as.integer(if (nrow(d)) d$donors_all[1] else 0L)),
                      side = 3, line = 0.7, adj = 0, font = 2, cex = 0.72, col = INK)
                if (nzchar(sub)) mtext(sub, side = 3, line = 0.7, adj = 1, cex = 0.62, col = MID)
            }
            mtext("minimum regions per donor", side = 1, outer = TRUE, line = 2.4, cex = 0.78, col = INK)
            mtext("fraction, or correlation", side = 2, outer = TRUE, line = 2.4, cex = 0.78, col = INK)
            ln <- nH + 0.1
            for (t in w1) { mtext(t, outer = TRUE, side = 3, line = ln, font = 2, cex = 0.78); ln <- ln - 1.15 }
            for (t in w2) { mtext(t, outer = TRUE, side = 3, line = ln, cex = 0.64, col = MID); ln <- ln - 1.15 }
            par(fig = c(0, 1, 0, 1), mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0), new = TRUE)
            plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
            legend("bottom", horiz = TRUE, bty = "n", cex = 0.72, seg.len = 2.2, lwd = 3.0, lty = 1,
                   col = unname(REG_COL), text.col = INK, legend = unname(REG_LAB))
        })
        p12_say("")
    }

    p12_run(main)
}
# =============================================================================
# RUN
# =============================================================================
local({
    v <- p12_section_25()
    if (!identical(v, "PASS"))
        stop(sprintf("section 25 ended with verdict %s. See Part12_Fig18_log.txt.", v), call. = FALSE)
    message("Figure 18: section 25 PASS. Figures in Part12_rerun/10_Fig18/Fig18")
})
