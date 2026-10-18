# =============================================================================
# Part12_Fig17_Plots.r
# SEA-AD CEP-IP/CMRF, Part 12, Figure 17 plots (Figure 17, Supplementary Figure 12)
#   Draws the Figure 17 and Supplementary Figure 12 panels from the stats and
#   plot-data files written by sections 19 and 20 of Part12_Fig17_Analysis.r,
#   after checking them against their logs. Computes no statistic.
#
# Usage (from the SEA-AD folder, or set SEAAD_ROOT), after Part12_Fig17_Analysis.r
#   source("Part12_rerun/Part12_Fig17_Plots.r")
#
# Outputs (Part12_rerun/09_Fig17/Fig17/): figures (PDF, PNG), log.
# =============================================================================

source(file.path(Sys.getenv("SEAAD_ROOT", getwd()), "Part12_rerun", "Part12_03_Config.r"))

SCRIPT <- "Part12_Fig17_Plots.r"
ANALYSIS <- "Part12_Fig17_Analysis.r"
FIG_DIR <- "09_Fig17"   # every Figure 17 section folder sits inside this one
## Config version 15 or later holds the colours drawn here.
if (!isTRUE(as.integer(P12_CFG$config_version) >= 15L))
    stop(sprintf("config is version %s; the Figure 17 plots need version 15 or later",
                 P12_CFG$config_version), call. = FALSE)

# == Helpers shared with Part12_Fig17_Analysis.r ==============================

## Verdict of a two-predictor row (17C, 17D).
p12_17_verdict <- function(fdr_m, fdr_t, vif, fitted) {
    alpha <- as.numeric(P12_CFG$inference$alpha)
    vr <- as.numeric(P12_CFG$floors$vif_readable)
    sm <- is.finite(fdr_m) & fdr_m < alpha
    st <- is.finite(fdr_t) & fdr_t < alpha
    ifelse(!fitted, "Untestable",
           ifelse(!(is.finite(vif) & vif < vr), "Unreadable",
                  ifelse(sm & !st, "Module only",
                         ifelse(!sm & st, "%MT only",
                                ifelse(sm & st, "Both", "Neither")))))
}

## Measure to stain family key, and family key to label and block.
p12_17_families <- function() {
    F <- P12_CFG$histology$families
    list(of = unlist(lapply(F, function(f) stats::setNames(rep(f$key, length(f$measures)),
                                                           f$measures))),
         keys = vapply(F, function(f) f$key, ""),
         label = stats::setNames(vapply(F, function(f) f$label, ""), vapply(F, function(f) f$key, "")),
         block = stats::setNames(vapply(F, function(f) f$block, ""), vapply(F, function(f) f$key, "")),
         n = stats::setNames(vapply(F, function(f) length(f$measures), 0L),
                             vapply(F, function(f) f$key, "")))
}

LEVELS <- c("within donor", "between donors")

## Parsed definitions of the objects in `what`, for an exact comparison.
p12_17_defs <- function(path, what) {
    out <- list()
    if (!file.exists(path)) return(out)
    for (e in parse(path, keep.source = FALSE))
        if (is.call(e) && identical(e[[1L]], as.name("<-")) && is.name(e[[2L]]) &&
            as.character(e[[2L]]) %in% what) out[[as.character(e[[2L]])]] <- e[[3L]]
    out
}


# =============================================================================
# SECTION 21: Figure 17 and Supplementary Figure 12
# Panels 17A to 17E and S12A, S12B; arms in the order M4, PVM_3, PVM_1, HKG8,
# NeoC before NonNeoC. A point is filled when its FDR is below 0.05.
# =============================================================================
p12_section_21 <- function() {
    STEP <- "21"
    D <- p12_step_dirs(STEP, "Fig17", parent = FIG_DIR)
    p12_clean(D$out)
    p12_start(STEP, SCRIPT, D, pkgs = c("ggplot2", "gtable"))

    PAL <- P12_CFG$palette
    INK <- PAL$neutral$dark; MID <- PAL$neutral$mid; LIGHT <- PAL$neutral$light
    ACCENT <- PAL$neutral$accent
    ## Strip colour for every Figure 17 panel: RGB 0, 32, 96 (#002060).
    HEAD <- grDevices::rgb(0, 32, 96, maxColorValue = 255)
    VCOL <- unlist(PAL$verdict); GCOL <- unlist(PAL$grain); VLINE <- PAL$vif_line
    HI <- P12_CFG$histology
    QLAB <- unlist(HI$quantities)
    QCOL <- c(module_post_pct = ACCENT, mito_log1p = INK, cyto_share = MID, contamination = MID,
              log_n_cells = LIGHT, abs_module_minus_abs_mito = INK)
    ARMS <- P12_ARMS
    VK <- c("Module only", "Both", "%MT only", "Neither", "Unreadable", "Untestable")
    VR <- as.numeric(P12_CFG$floors$vif_readable)
    R_LINE <- sqrt(1 - 1 / VR)
    B <- as.integer(P12_CFG$inference$bootstrap_B)

    ## Forest sizes as printed (base size reset to 1 after layout()).
    F_PT <- 1.45 * 0.66; F_LWD <- 1.4; F_PT_LWD <- 1.4; F_STAR <- 1.45 * 0.66
    F_AXIS <- 0.67; F_ROWLAB <- 0.67; F_XLAB <- 0.83; F_STRIP <- 0.75; F_TITLE <- 1.0
    F_CAP <- 0.62; F_KEY <- 0.67; F_KEY_PT <- 1.0; F_BAND <- "#F4F4F4"; F_ZERO <- "#808080"
    ## Inset notes and tick-label minus signs; S12 changes them via with_forest_style().
    F_NOTE_CEX <- 0.55; F_NOTE_COL <- MID; F_HYPHEN <- FALSE
    ## ggplot sizes as base R cex at pointsize 12: 12 pt = 1.0, 10 pt = 0.83, 9 pt = 0.75.
    G_TITLE <- 1.0; G_SUB <- 0.67; G_AXIS <- 0.67; G_AXTITLE <- 0.83; G_STRIP <- 0.75
    G_CAP <- 0.62; G_COUNT <- 0.71; G_KEY <- 0.67; G_BAR_H <- 0.68; G_BORDER <- 0.85
    G_GRID <- "grey92"; G_FRAME <- "grey20"; G_TEXT <- "grey30"
    STRIP_IN <- 0.22

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
    wrap <- function(s, w) vapply(s, function(z) paste(strwrap(z, w), collapse = "\n"), "")
    ## Lines of text that fit across `w_in` inches at size `cex` (greedy wrap).
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
    ## Footer height in inches: caption lines, then the key above them.
    foot_in <- function(n_lines, key, key_rows = 1) 0.12 + 0.15 * n_lines +
        (if (key) 0.28 + 0.14 * key_rows else 0.04)
    ## Footer: caption lines and star note at fixed offsets, the key above them.
    footer <- function(lines, key = NULL, cap_cex = F_CAP, key_cex = F_KEY) {
        ## Drawn last on a device, so the layout parameters are not restored.
        par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mai = c(0, 0, 0, 0), new = TRUE, xpd = NA,
            cex = 1)
        plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "",
             xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
        iy <- 1 / par("din")[2]
        for (i in seq_along(lines))
            text(0.5, (0.12 + 0.15 * (length(lines) - i)) * iy, lines[i], cex = cap_cex, col = G_TEXT)
        if (!is.null(key))
            do.call(legend, c(list(x = 0.5, y = (0.12 + 0.15 * length(lines) + 0.02) * iy,
                                   xjust = 0.5, yjust = 0, horiz = is.null(key$ncol), bty = "n",
                                   cex = key_cex), key))
    }
    ## Strips attached to the plot region.
    strip_top <- function(label, cex, col = HEAD) {
        u <- par("usr"); dy <- diff(grconvertY(c(0, STRIP_IN), "inches", "user"))
        rect(u[1], u[4], u[2], u[4] + dy, col = col, border = col, xpd = NA)
        text(mean(u[1:2]), u[4] + dy / 2, label, col = "white", font = 2, cex = cex, xpd = NA)
    }
    strip_right <- function(label, cex) {
        u <- par("usr"); dx <- diff(grconvertX(c(0, STRIP_IN), "inches", "user"))
        rect(u[2], u[3], u[2] + dx, u[4], col = HEAD, border = HEAD, xpd = NA)
        text(u[2] + dx / 2, mean(u[3:4]), label, col = "white", font = 2, srt = 270, cex = cex,
             xpd = NA)
    }
    ## Forest x range: symmetric, max(`floor`, widest interval) x 1.12.
    ## Title above subtitle lines, centred on the device, drawn on an overlay.
    heading <- function(tl, sl = character(0), tcex = G_TITLE, scex = G_SUB) {
        par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mai = c(0, 0, 0, 0), new = TRUE, xpd = NA,
            cex = 1)
        plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "",
             xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
        iy <- 1 / par("din")[2]
        y <- 1 - 0.10 * iy
        for (k in seq_along(tl)) {
            text(0.5, y, tl[k], adj = c(0.5, 1), cex = tcex, font = 2); y <- y - 0.22 * iy
        }
        for (k in seq_along(sl)) {
            text(0.5, y, sl[k], adj = c(0.5, 1), cex = scex, col = G_TEXT); y <- y - 0.15 * iy
        }
    }
    fx <- function(v, floor = 0.4) { m <- max(floor, abs(v[is.finite(v)])) * 1.12; c(-m, m) }
    ## As fx(), with room for a star beyond the widest interval end.
    fx_star <- function(v, floor = 0.4) { m <- max(floor, abs(v[is.finite(v)])) * 1.12 + 0.15; c(-m, m) }

    ## One forest panel; rows top to bottom, a row's series over 0.74 of the row.
    ## `ticks` fixes x ticks; `star_xpd` lets stars enter the right margin;
    ## `star_nudge = FALSE` keeps each star at the end of its own interval.
    forest <- function(rows, qs, get, xl, labels = NULL, note = NULL, ticks = NULL, star_xpd = FALSE,
                       star_nudge = TRUE) {
        n <- length(rows); ns <- length(qs); step <- 0.74 / max(1, ns)
        plot(NA, xlim = xl, ylim = c(0.4, n + 0.6), xaxt = "n", yaxt = "n", xlab = "", ylab = "",
             xaxs = "i", yaxs = "i")
        for (i in seq_len(n)) if (i %% 2L == 1L)
            rect(xl[1], n - i + 0.5, xl[2], n - i + 1.5, col = F_BAND, border = NA)
        abline(v = 0, col = F_ZERO)
        tk <- if (is.null(ticks)) pretty(xl, n = 4) else ticks
        tk <- tk[tk > xl[1] + 1e-9 & tk < xl[2] - 1e-9]
        ## With F_HYPHEN, negative tick labels use a soft hyphen (the PDF device
        ## draws "-" as a long minus sign).
        tl <- format(tk, trim = TRUE)
        if (F_HYPHEN) tl <- gsub("-", if (identical(names(grDevices::dev.cur()), "pdf")) "\u00ad" else "-", tl)
        axis(1, at = tk, labels = tl, cex.axis = F_AXIS, mgp = c(3, 0.35, 0), tcl = -0.25)
        for (i in seq_len(n)) {
            y0 <- n - i + 1
            any_r <- FALSE
            placed <- list()   # stars already drawn in this row, as c(x0, x1, y)
            for (j in seq_len(ns)) {
                r <- get(rows[i], qs[j])
                if (is.null(r) || !nrow(r) || !is.finite(r$estimate[1])) next
                any_r <- TRUE
                y <- y0 + ((ns + 1) / 2 - j) * step
                cq <- QCOL[[qs[j]]]
                if (is.finite(r$lo[1]))
                    segments(max(r$lo[1], xl[1]), y, min(r$hi[1], xl[2]), y, col = cq, lwd = F_LWD)
                sig <- is.finite(r$fdr[1]) && r$fdr[1] < 0.05
                points(r$estimate[1], y, pch = if (sig) 19L else 1L, col = cq, bg = "white",
                       cex = F_PT, lwd = F_PT_LWD)
                st <- p12_stars(r$fdr[1])
                if (nzchar(st)) {
                    at <- if (is.finite(r$hi[1])) r$hi[1] else r$estimate[1]
                    xs <- min(at, xl[2] - 0.10 * abs(xl[2]))
                    ## A star that would overlap one already drawn (< 0.095 in apart) moves right.
                    off <- 0.18 * strwidth("m", cex = F_STAR)
                    w <- strwidth(st, cex = F_STAR, font = 2)
                    if (star_nudge) for (b in placed)
                        if (abs(diff(grconvertY(c(y, b[3]), "user", "inches"))) < 0.095 &&
                            xs + off < b[2] && xs + off + w > b[1])
                            xs <- b[2] + 0.05 * w - off
                    text(xs, y, st, pos = 4, offset = 0.18, col = cq, cex = F_STAR, font = 2,
                         xpd = if (star_xpd) NA else FALSE)
                    placed[[length(placed) + 1L]] <- c(xs + off, xs + off + w, y)
                }
            }
            if (!any_r) text(0, y0, "untestable", col = F_NOTE_COL, cex = F_AXIS, font = 3)
            if (!is.null(note) && nzchar(note[i]))
                text(xl[1] + 0.02 * diff(xl), y0 + 0.43, note[i], adj = c(0, 0.5), cex = F_NOTE_CEX,
                     col = F_NOTE_COL)
        }
        if (!is.null(labels))
            text(xl[1] - 0.03 * diff(xl), rev(seq_len(n)), labels, adj = c(1, 0.5), cex = F_ROWLAB,
                 font = 2, xpd = NA)
        box()
    }
    ## Forests: one row of arm panels per block, a block being
    ## list(rows, labels, get = function(arm, row, q), level, note = function(arm)).
    forest_fig <- function(blocks, qs, xl, xlab, title, caption, key = TRUE) {
        din <- par("din")
        lab_in <- max(unlist(lapply(blocks, function(b)
            strwidth(unlist(strsplit(b$labels, "\n")), units = "inches", cex = F_ROWLAB,
                     font = 2)))) + 0.18
        star <- sprintf("* FDR < 0.05, ** FDR < 0.01 (BH). Filled point: FDR < 0.05. Bars: 95%% donor-bootstrap intervals, B = %s.",
                        format(B, big.mark = ","))
        lines <- c(fit_lines(caption, din[1] - 0.6, F_CAP), star)
        bot <- foot_in(length(lines), key)
        tl <- fit_lines(title, din[1] - 0.4, F_TITLE, 2)
        top <- 0.12 + 0.22 * length(tl)
        lev_in <- if (any(vapply(blocks, function(b) !is.null(b$level), TRUE))) 0.24 else 0
        mt <- STRIP_IN + 0.04 + lev_in; mb <- 0.42; ml1 <- lab_in; ml <- 0.08; mr <- 0.18
        nr <- vapply(blocks, function(b) length(b$rows), 0)
        Hp <- din[2] - top - bot - length(blocks) * (mt + mb)
        heights <- nr / sum(nr) * Hp + mt + mb
        Wp <- din[1] - ml1 - 3 * ml - 4 * mr
        widths <- c(Wp / 4 + ml1 + mr, rep(Wp / 4 + ml + mr, 3))
        layout(matrix(seq_len(4 * length(blocks)), length(blocks), byrow = TRUE),
               widths = widths, heights = heights)
        par(omi = c(bot, 0, top, 0), family = "sans", cex = 1)
        for (bi in seq_along(blocks)) {
            b <- blocks[[bi]]
            for (i in seq_along(ARMS)) {
                a <- ARMS[i]
                par(mai = c(mb, if (i == 1L) ml1 else ml, mt, mr))
                forest(b$rows, qs, function(row, q) b$get(a, row, q), xl,
                       labels = if (i == 1L) b$labels else NULL,
                       note = if (is.null(b$note)) NULL else b$note(a))
                strip_top(a, F_STRIP)
                mtext(xlab, side = 1, line = 1.35, cex = F_XLAB)
                if (i == 1L && !is.null(b$level)) {
                    u <- par("usr")
                    text(grconvertX(0.02, "nfc", "user"),
                         u[4] + diff(grconvertY(c(0, STRIP_IN + 0.06), "inches", "user")),
                         b$level, adj = c(0, 0), font = 2, cex = F_XLAB, xpd = NA)
                }
            }
        }
        heading(tl)
        footer(lines, if (key) list(legend = QLAB[qs], pch = 19, pt.cex = F_KEY_PT,
                                    col = QCOL[qs]) else NULL)
    }

    ## Draws f() with forest fonts x k, notes in note_col and hyphen minus signs.
    with_forest_style <- function(k, note_col, f) {
        env <- environment(forest)
        nm <- c("F_AXIS", "F_ROWLAB", "F_XLAB", "F_STRIP", "F_TITLE", "F_CAP", "F_KEY", "F_STAR", "F_PT",
                "F_KEY_PT", "F_NOTE_CEX")
        old <- mget(c(nm, "F_NOTE_COL", "F_HYPHEN"), envir = env)
        on.exit(list2env(old, envir = env))
        for (v in nm) assign(v, old[[v]] * k, envir = env)
        assign("F_NOTE_COL", note_col, envir = env)
        assign("F_HYPHEN", TRUE, envir = env)
        f()
    }

    ## 17B: 2 x 4 panels with margins fixed in inches, so each plot region is
    ## 1.561 x 2.229 in and each panel carries its own row labels.
    forest_mfrow <- function(blocks, qs, xl, xlab, title, caption) {
        LN <- 0.132
        par(mfrow = c(length(blocks), length(ARMS)))
        par(mai = c(4.2, 8.5, 3, 2.4) * LN, omi = c(6.8, 0, 2.6, 0) * LN, family = "sans", cex = 1)
        for (bi in seq_along(blocks)) {
            b <- blocks[[bi]]
            for (i in seq_along(ARMS)) {
                a <- ARMS[i]
                forest(b$rows, qs, function(row, q) b$get(a, row, q), xl, labels = b$labels,
                       ticks = c(-0.5, 0, 0.5), star_xpd = TRUE, star_nudge = FALSE)
                strip_top(a, F_STRIP)
                mtext(xlab, side = 1, line = 1.35, cex = F_XLAB)
                if (i == 1L) {
                    u <- par("usr")
                    text(grconvertX(0.02, "nfc", "user"),
                         u[4] + diff(grconvertY(c(0, STRIP_IN + 0.04), "inches", "user")),
                         b$level, adj = c(0, 0), font = 2, cex = F_XLAB, xpd = NA)
                }
            }
        }
        star <- sprintf("* FDR < 0.05, ** FDR < 0.01 (BH). Filled point: FDR < 0.05. Bars: 95%% donor-bootstrap intervals, B = %s.",
                        format(B, big.mark = ","))
        lines <- c(fit_lines(caption, par("din")[1] - 0.6, F_CAP), star)
        heading(fit_lines(title, par("din")[1] - 0.4, F_TITLE, 2))
        footer(lines, list(legend = QLAB[qs], pch = 19, pt.cex = F_KEY_PT, col = QCOL[qs]))
    }

    ## 17C and 17D in ggplot2: theme_bw base size 10, bars of height 0.68 with
    ## black borders, counts in white or grey15 by contrast, facet_grid(arm ~ block),
    ## strip fill RGB 0, 32, 96, and a fixed panel height (`panel_h`).
    ## `x_min_range` widens a named column's x range so its strip label fits;
    ## `gap_x` is the column spacing in points.
    gg_duel <- function(tab, cols, collab, bars, barlab, xlab, title, sub, caption, panel_h,
                        title_w, sub_w, x_min_range = NULL, gap_x = 14) {
        tal <- as.data.frame(table(arm = factor(tab$arm, levels = ARMS),
                                   col = factor(tab$col, levels = cols),
                                   bar = factor(tab$bar, levels = rev(bars)),
                                   verdict = factor(tab$verdict, levels = VK)))
        tal <- tal[tal$Freq > 0, , drop = FALSE]
        tal <- tal[order(tal$arm, tal$col, tal$bar, tal$verdict), , drop = FALSE]
        key <- interaction(tal$arm, tal$col, tal$bar, drop = TRUE)
        tal$xmax <- stats::ave(tal$Freq, key, FUN = cumsum)
        tal$xmin <- tal$xmax - tal$Freq
        tal$xmid <- tal$xmax - tal$Freq / 2
        tal$ypos <- as.integer(tal$bar)
        ## White or grey15, whichever contrasts more (WCAG).
        lum <- function(h) {
            f <- grDevices::col2rgb(h)[, 1] / 255
            l <- ifelse(f <= 0.04045, f / 12.92, ((f + 0.055) / 1.055)^2.4)
            sum(c(0.2126, 0.7152, 0.0722) * l)
        }
        tal$txt_col <- vapply(as.character(tal$verdict), function(k) {
            L <- lum(VCOL[[k]]); if (1.05 / (L + 0.05) >= (L + 0.05) / 0.0685) "white" else "grey15" }, "")
        wrap <- function(z, w) paste(strwrap(z, width = w), collapse = "\n")
        p <- ggplot2::ggplot(tal) +
            ggplot2::geom_rect(ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ypos - 0.34,
                                            ymax = ypos + 0.34, fill = verdict),
                               colour = "black", linewidth = 0.3) +
            ggplot2::geom_text(ggplot2::aes(x = xmid, y = ypos, label = Freq, colour = txt_col),
                               size = 3.0) +
            (if (length(x_min_range))
                ggplot2::geom_blank(data = data.frame(col = factor(names(x_min_range), levels = cols),
                                                      x = unname(x_min_range), y = 1),
                                    ggplot2::aes(x = x, y = y))) +
            ggplot2::scale_colour_identity() +
            ggplot2::scale_fill_manual(values = VCOL[VK], breaks = VK, name = NULL, drop = FALSE) +
            ggplot2::scale_x_continuous(breaks = function(l) {
                b <- unique(round(pretty(l))); b[is.finite(b) & b >= 0] }) +
            ggplot2::scale_y_continuous(breaks = seq_along(bars), labels = rev(barlab),
                                        expand = ggplot2::expansion(add = 0.6)) +
            ggplot2::facet_grid(arm ~ col, scales = "free_x", space = "free_x",
                                labeller = ggplot2::labeller(col = stats::setNames(
                                    vapply(collab, wrap, "", w = 13), cols))) +
            ggplot2::labs(title = wrap(title, title_w), subtitle = wrap(sub, sub_w),
                          caption = paste(vapply(caption, wrap, "", w = 45), collapse = "\n"),
                          x = xlab, y = NULL) +
            ggplot2::theme_bw(base_size = 10) +
            ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 12),
                           plot.subtitle = ggplot2::element_text(size = 8, colour = "grey30"),
                           plot.caption = ggplot2::element_text(size = 7.5, colour = "grey30", hjust = 0),
                           strip.background = ggplot2::element_rect(fill = HEAD),
                           strip.text = ggplot2::element_text(colour = "white", face = "bold", size = 9),
                           panel.grid.minor = ggplot2::element_blank(),
                           legend.position = "bottom",
                           panel.spacing = grid::unit(14, "pt"),
                           panel.spacing.x = grid::unit(gap_x, "pt"),
                           plot.background = ggplot2::element_rect(fill = "white", colour = NA))
        g <- ggplot2::ggplotGrob(p)
        rows <- unique(g$layout$t[grepl("^panel", g$layout$name)])
        g$heights[rows] <- grid::unit(panel_h, "in")
        g <- gtable::gtable_add_rows(g, grid::unit(1, "null"), -1)
        ## No grid.newpage(): fig_dev() already opens a fresh device.
        grid::grid.draw(g)
    }


    main <- function() {
        p12_say("SECTION A -- inputs")
        ## Stats logs written by the current analysis file; shared helpers identical.
        for (fd in c("17AB_Histology", "17CE_ModuleVsMT")) {
            lg <- file.path(P12_RERUN, FIG_DIR, fd, sprintf("Part12_%s_log.txt", fd))
            sl <- if (file.exists(lg)) sub("\r$", "", grep("^Script: ", readLines(lg, warn = FALSE),
                                                             value = TRUE)[1]) else NA_character_
            p12_chk(identical(sl, sprintf("Script: %s", ANALYSIS)),
                    sprintf("%s log was written by %s", fd, ANALYSIS),
                    if (is.na(sl)) "log or Script line not found" else sl)
        }
        sh <- c("p12_17_families", "p12_17_verdict")
        da <- p12_17_defs(file.path(P12_RERUN, ANALYSIS), sh)
        dp <- p12_17_defs(file.path(P12_RERUN, SCRIPT), sh)
        same <- vapply(sh, function(n) !is.null(da[[n]]) && !is.null(dp[[n]]) &&
                           identical(da[[n]], dp[[n]]), TRUE)
        p12_chk(all(same), sprintf("shared helpers identical to %s", ANALYSIS),
                paste(sprintf("%s %s", sh, ifelse(same, "same", "DIFFERS")), collapse = " | "))
        A <- rd("17AB_Histology", "Stats_17A.csv")
        Bt <- rd("17AB_Histology", "Stats_17B.csv")
        S <- rd("17AB_Histology", "Stats_S12.csv")
        SM <- rd("17AB_Histology", "Stats_S12_measures.csv")
        SC <- rd("17AB_Histology", "Plot_S12_counts.csv")
        C <- rd("17CE_ModuleVsMT", "Stats_17C.csv")
        Dd <- rd("17CE_ModuleVsMT", "Stats_17D.csv")
        E <- rd("17CE_ModuleVsMT", "Stats_17E.csv")
        CV <- rd("17CE_ModuleVsMT", "Plot_17C_verdicts.csv")
        DV <- rd("17CE_ModuleVsMT", "Plot_17D_verdicts.csv")
        FAM <- p12_17_families()
        IBA <- vapply(HI$iba1_17A, function(z) z$measure, "")
        IBL <- vapply(HI$iba1_17A, function(z) z$label, "")
        CPSL <- stats::setNames(vapply(HI$cps, function(z) z$label, ""),
                                vapply(HI$cps, function(z) z$column, ""))
        CPS <- names(CPSL); LOCAL <- CPS[vapply(HI$cps, function(z) z$scope, "") == "local"]

        ## Drawn verdicts and family medians follow from the stats files.
        gc <- function(q) C[C$quantity == q, ][match(paste(CV$arm, CV$level, CV$measure),
                                                     with(C[C$quantity == q, ], paste(arm, level, measure))), ]
        v2 <- p12_17_verdict(gc("beta_module")$fdr, gc("beta_mito")$fdr, gc("vif")$estimate,
                             is.finite(gc("beta_module")$estimate) & is.finite(gc("beta_mito")$estimate))
        p12_chk(identical(v2, CV$verdict), "17C verdicts follow from Stats_17C",
                sprintf("%d rows", nrow(CV)))
        gd <- function(q) Dd[Dd$quantity == q, ][match(paste(DV$arm, DV$dataset, DV$outcome),
                                                       with(Dd[Dd$quantity == q, ], paste(arm, level, measure))), ]
        v3 <- p12_17_verdict(gd("beta_module")$fdr, gd("beta_mito")$fdr, gd("vif")$estimate,
                             is.finite(gd("beta_module")$estimate) & is.finite(gd("beta_mito")$estimate))
        p12_chk(identical(v3, DV$verdict), "17D verdicts follow from Stats_17D",
                sprintf("%d rows", nrow(DV)))
        SA <- S[S$quantity != "abs_module_minus_abs_mito", ]
        med <- vapply(seq_len(nrow(SA)), function(i) {
            ms <- names(FAM$of)[FAM$of == SA$measure[i]]
            v <- SM$estimate[SM$arm == SA$arm[i] & SM$level == SA$level[i] &
                                 SM$quantity == SA$quantity[i] & SM$measure %in% ms]
            v <- v[is.finite(v)]; if (length(v)) stats::median(v) else NA_real_ }, 0)
        same <- (is.na(med) & is.na(SA$estimate)) | abs(med - SA$estimate) < 1e-12
        p12_chk(all(same, na.rm = TRUE) && !anyNA(same), "S12 family medians follow from the per-measure rs",
                sprintf("%d family rows", nrow(SA)))
        nt <- vapply(seq_len(nrow(SC)), function(i) {
            ms <- names(FAM$of)[FAM$of == SC$family[i]]
            sum(is.finite(SM$estimate[SM$arm == SC$arm[i] & SM$level == SC$level[i] &
                                          SM$quantity == SC$quantity[i] & SM$measure %in% ms])) }, 0)
        p12_chk(all(nt == SC$n_tested), "S12 tested counts follow from the per-measure rs",
                sprintf("%d rows", nrow(SC)))
        labs <- c(IBL, CPSL, FAM$label, FAM$block, QLAB, VK, names(GCOL))
        p12_chk(p12_ascii_ok(labs), "every figure label is ASCII", sprintf("%d labels", length(labs)))
        ## Figure text is typed in this file, so the whole file is checked too.
        src <- readLines(file.path(P12_RERUN, SCRIPT), warn = FALSE)
        p12_chk(p12_ascii_ok(src), "figure script text is ASCII", sprintf("%d lines", length(src)))
        if (P12$n[["FAIL"]] > 0L) p12_abort("checks failed; figures not drawn")
        p12_say("")
        p12_say("SECTION B -- figures")

        # == 17A ================================================================
        fig_dev("Fig17A_microglial_shape", 9, 5.6 * 1.2 * 0.95, function() {
            QD <- HI$drawn[["17A"]]
            g <- A[A$level == "within donor", ]
            forest_fig(list(list(rows = IBA, labels = wrap(IBL, 22),
                                 get = function(a, m, q) g[g$arm == a & g$measure == m & g$quantity == q, ])),
                       QD, fx(unlist(g[g$quantity %in% QD, c("lo", "hi")])), "Spearman rs (within donor)",
                       "Does microglial shape track Module post-IP%, or just %MT?",
                       paste("Within donor: every variable centred on the donor; NeoC and NonNeoC pooled.",
                             "BH within arm and quantity over the 6 IBA1 measures."))
        })

        # == 17B ================================================================
        fig_dev("Fig17B_cps", 13.0 * 0.95 * 0.98 * 1.03, 7.6, function() {
            QD <- HI$drawn[["17B"]]
            BET <- c(CPS[!(CPS %in% LOCAL)], LOCAL)
            getb <- function(lv) function(a, m, q)
                Bt[Bt$level == lv & Bt$arm == a & Bt$measure == m & Bt$quantity == q, ]
            forest_mfrow(list(list(rows = LOCAL, labels = CPSL[LOCAL], get = getb("within donor"),
                                   level = "Within donor"),
                              list(rows = BET, labels = CPSL[BET], get = getb("between donors"),
                                   level = "Between donors")),
                         ## x range capped at -0.95 to +0.95.
                         QD, c(-0.95, 0.95),
                         "Spearman rs", "Does nearby pathology track Module post-IP%, or just %MT?",
                         paste("CPS: continuous pseudo-progression score. Top: within donor, regions centred",
                               "on the donor (the global scores are constant within donor). Bottom: donor means.",
                               "BH within arm, level and quantity over the scores of that level."))
        })

        # == 17C ================================================================
        real <- CV[CV$arm != "HKG8", ]
        n_mt <- sum(real$verdict == "%MT only"); n_fl <- sum(real$verdict == "%MT only" & real$module_sign_flip)
        fig_dev("Fig17C_duel_histology", 14 * 0.4125 * 0.95, (3.2 + 0.86 * 4) * 1.05 * 1.03 * 1.01,
                function() {
            BL <- HI$blocks
            tab <- data.frame(arm = CV$arm, col = CV$block, bar = CV$level, verdict = CV$verdict)
            gg_duel(tab, BL, BL, LEVELS,
                    c("Within donor, histology units", "Between donors, histology units"), "measures",
                    "Module post-IP% vs %MT in one model (histology)",
                    paste("rank(measure) ~ rank(Module post-IP%) + rank(%MT), no covariates; NeoC and NonNeoC",
                          "pooled. Verdict at FDR < 0.05, BH within arm, level and predictor over the 28",
                          "measures. A bar sums to the measures in that block."),
                    c(sprintf("Unreadable: VIF >= %g, the two cannot be separated in the joint model.", VR),
                      "Untestable: fewer than 10 donors with 3 or more units.",
                      if (n_mt) sprintf(paste("On %d of the %d %%MT-only rows on the real arms the Module",
                                              "coefficient reversed sign against its marginal rs:",
                                              "shared variance the model could not allocate."), n_fl, n_mt)
                      else "No %MT-only row on the real arms."),
                    panel_h = 0.481 * 1.05 * 1.03 * 1.01, title_w = 36, sub_w = 72)
        })

        # == 17D ================================================================
        fig_dev("Fig17D_duel_clinical", 11 * 0.375, (3.2 + 0.86 * 4) * 1.06, function() {
            CLS <- c(cognition = "Cognition", pathology = "Pathology", genotype = "APOE")
            tab <- data.frame(arm = DV$arm, col = DV$class, bar = DV$dataset, verdict = DV$verdict)
            gg_duel(tab, names(CLS), unname(CLS), P12_DATASETS, P12_DATASETS, "outcomes",
                    "Module post-IP% vs %MT in one model (clinical)",
                    paste("outcome ~ cov4 + Module post-IP% + %MT, standardised; region sets fitted",
                          "separately. Verdict at FDR < 0.05, BH within arm, region set and predictor over",
                          "the 9 outcomes."),
                    c(sprintf("Unreadable: VIF >= %g. cov4: age, APOE e4 dose, sex, PMI.", VR),
                      "Each arm carries the %MT of its own cells."),
                    panel_h = 0.593 * 1.06, title_w = 34, sub_w = 72, x_min_range = c(genotype = 1.55),
                    gap_x = 7)
        })

        # == 17E ================================================================
        fig_dev("Fig17E_collinearity", 11 * 0.4125 * 1.10, 2.6 + 1.05 * 4, function() {
            RW <- data.frame(level = c("within donor HU", "between donors HU", rep("between donors CG", 2)),
                             measure = c("pooled", "pooled", P12_DATASETS),
                             label = c("Within donor (HU)", "Between donors (HU)",
                                       sprintf("Between donors (CG), %s", P12_DATASETS)),
                             stringsAsFactors = FALSE)
            g <- E[E$quantity == "abs_r", ]; gv <- E[E$quantity == "vif", ]
            din <- par("din")
            tl <- fit_lines("Collinearity between Module post-IP% and %MT", din[1] - 0.3, G_TITLE, 2)
            sl <- fit_lines(sprintf(paste("Right of the dashed line the joint model cannot separate the",
                                          "two (VIF >= %g, |r| %.3f)."), VR, R_LINE), din[1] - 0.3, G_SUB)
            lines <- fit_lines(c(paste("HU: histology units joined to a neuropathology row (Spearman,",
                                       "NeoC and NonNeoC pooled, as 17C). CG: clinical grain, donor",
                                       "table with cov4 complete (Pearson, as 17D)."),
                                 "95% donor-bootstrap intervals are in SuppTable_Fig17 (Fig17E_readout)."),
                               din[1] - 0.3, G_CAP)
            top <- 0.10 + 0.22 * length(tl) + 0.15 * length(sl) + 0.06
            bot <- 0.24 + foot_in(length(lines), TRUE)
            lab_in <- max(strwidth(RW$label, units = "inches", cex = G_AXIS)) + 0.10
            ph <- (din[2] - top - bot - 3 * 0.08 - 0.26) / 4
            layout(matrix(1:4, 4), heights = c(rep(ph + 0.08, 3), ph + 0.26))
            par(omi = c(bot, lab_in, top, 0.05), family = "sans", cex = 1)
            for (ai in seq_along(ARMS)) {
                a <- ARMS[ai]
                par(mai = c(if (ai == 4L) 0.22 else 0.04, 0.03, 0.04, 0.03 + STRIP_IN))
                nr <- nrow(RW)
                plot(NA, xlim = c(0, 1.45), ylim = c(0.4, nr + 0.6), axes = FALSE, xlab = "", ylab = "",
                     xaxs = "i", yaxs = "i")
                abline(v = seq(0, 1, 0.5), col = G_GRID, lwd = 0.8)
                abline(v = R_LINE, col = VLINE, lwd = 0.6 * 2.845, lty = 2)
                for (i in seq_len(nr)) {
                    y <- nr - i + 1
                    r <- g[g$arm == a & g$level == RW$level[i] & g$measure == RW$measure[i], ]
                    v <- gv[gv$arm == a & gv$level == RW$level[i] & gv$measure == RW$measure[i], ]
                    cg <- GCOL[[RW$level[i]]]
                    if (!nrow(r) || !is.finite(r$estimate)) {
                        text(0.03, y, "untestable: fewer than 10 donors", adj = 0, cex = G_AXIS,
                             col = MID, font = 3)
                        next
                    }
                    ## Lollipop only; the intervals are in the SuppTable.
                    segments(0, y, r$estimate, y, col = cg, lwd = 1.1 * 2.845)
                    points(r$estimate, y, pch = 19, col = cg, cex = 1.3)
                    ## Label on a white box, so the dashed VIF line passes behind it.
                    vl <- if (is.finite(v$estimate)) sprintf("VIF=%.1f", v$estimate) else "VIF=Inf"
                    vcex <- 0.82
                    vx <- r$estimate + diff(grconvertX(c(0, 0.09), "inches", "user"))
                    pad <- diff(grconvertX(c(0, 0.015), "inches", "user"))
                    vw <- strwidth(vl, cex = vcex); vh <- strheight(vl, cex = vcex)
                    rect(vx - pad, y - 0.8 * vh, vx + vw + pad, y + 0.8 * vh, col = "white", border = NA)
                    text(vx, y, vl, adj = c(0, 0.5), cex = vcex, col = "black")
                }
                box(col = G_FRAME, lwd = 1)
                if (ai == 4L) axis(1, at = seq(0, 1, 0.5), cex.axis = G_AXIS, col = G_FRAME,
                                   col.axis = G_TEXT, mgp = c(3, 0.2, 0), tcl = -0.2, lwd = 0.8)
                text(par("usr")[1] - diff(grconvertX(c(0, 0.06), "inches", "user")), rev(seq_len(nr)),
                     RW$label, adj = c(1, 0.5), cex = G_AXIS, col = G_TEXT, xpd = NA)
                strip_right(a, G_STRIP)
            }
            mtext("|r|, Module post-IP% against %MT", outer = TRUE, side = 1, line = 0.1,
                  cex = G_AXTITLE)
            heading(tl, sl)
            footer(lines, list(legend = c("Within donor (HU)", "Between donors (HU)", "Between donors (CG)"),
                               pch = 19, pt.cex = 1.2, col = unname(GCOL)), cap_cex = G_CAP, key_cex = G_KEY)
        })

        # == S12A and S12B =======================================================
        FL <- sprintf("%s (%d)", FAM$label, FAM$n)
        fig_dev("FigS12A_stain_families", 9, 6.0, function() {
            QD <- HI$drawn[["S12"]]
            g <- S[S$level == "within donor" & S$quantity %in% QD, ]
            note <- function(a) {
                cn <- SC[SC$arm == a & SC$level == "within donor" & SC$quantity == "module_post_pct", ]
                cn <- cn[match(FAM$keys, cn$family), ]
                ifelse(cn$n_tested < cn$n_total, sprintf("tested %d of %d", cn$n_tested, cn$n_total), "")
            }
            ## Fonts x 1.25, inset notes in grey10, hyphen minus signs.
            with_forest_style(1.25, "grey10", function()
            forest_fig(list(list(rows = FAM$keys, labels = FL, note = note,
                                 get = function(a, f, q) g[g$arm == a & g$measure == f & g$quantity == q, ])),
                       QD, fx(c(g$lo, g$hi)), expression(bold("Family median r"[s])),
                       "Stain families against Module post-IP%, %MT, contamination and log cell count",
                       paste("Within donor. Points are the median Spearman rs over the family's tested measures;",
                             "(n) is the measures in the family. BH within arm and quantity over the 7 families.")))
        })
        fig_dev("FigS12B_module_minus_mt", 9, 5.4, function() {
            g <- S[S$level == "within donor" & S$quantity == "abs_module_minus_abs_mito", ]
            with_forest_style(1.25, "grey10", function()
            forest_fig(list(list(rows = FAM$keys, labels = FL,
                                 get = function(a, f, q) g[g$arm == a & g$measure == f, ])),
                       "abs_module_minus_abs_mito", fx(c(g$lo, g$hi), 0.25),
                       ## Short axis title; the caption states the difference.
                       expression(bold("Difference in |r"[s] * "|")),
                       "Does Module post-IP% or %MT track each stain family more strongly?",
                       paste("Within donor. Difference = |rs| Module minus |rs| %MT; family median of the paired",
                             "difference, on the same donor resamples; negative means %MT tracks the family more",
                             "strongly. BH within arm",
                             "over the 7 families."), key = FALSE))
        })
        p12_say("")
    }

    p12_run(main)
}
# =============================================================================
# Run
# =============================================================================
local({
    v <- p12_section_21()
    if (!identical(v, "PASS"))
        stop(sprintf("section 21 ended with verdict %s. See Part12_Fig17_log.txt.", v), call. = FALSE)
    message("Figure 17: section 21 PASS. Figures in Part12_rerun/09_Fig17/Fig17")
})
