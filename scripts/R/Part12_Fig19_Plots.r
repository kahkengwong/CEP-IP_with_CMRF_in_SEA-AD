# =============================================================================
# Part12_Fig19_Plots.r
# SEA-AD CEP-IP/CMRF, Part 12, Figure 19 plots (Figure 19)
#   Draws Figures 19A to 19C and the supplementary verdict table from the
#   readout, null and verdict files of sections 26 to 29 of
#   Part12_Fig19_Analysis.r, after checking them against their logs and
#   recomputing the classes from the null files. Computes no statistic.
#
# Usage (from the SEA-AD folder, or set SEAAD_ROOT), after Part12_Fig19_Analysis.r
#   source("Part12_rerun/Part12_Fig19_Plots.r")
#
# Outputs (Part12_rerun/11_Fig19/Fig19/): figures (PDF, PNG), log.
# =============================================================================


source(file.path(Sys.getenv("SEAAD_ROOT", getwd()), "Part12_rerun", "Part12_03_Config.r"))

SCRIPT <- "Part12_Fig19_Plots.r"
ANALYSIS <- "Part12_Fig19_Analysis.r"
FIG_DIR <- "11_Fig19"   # every Figure 19 section folder sits inside this one
if (!isTRUE(as.integer(P12_CFG$config_version) >= 16L))
    stop(sprintf("config is version %s; Figure 19 needs version 16 or later",
                 P12_CFG$config_version), call. = FALSE)
B5_CFG_PATH <- file.path(P12_RERUN, "Part12_Fig19_config.json")
if (!file.exists(B5_CFG_PATH)) stop("Battery 5 config not found: ", B5_CFG_PATH, call. = FALSE)
B5 <- jsonlite::fromJSON(B5_CFG_PATH, simplifyVector = TRUE, simplifyDataFrame = FALSE,
                         simplifyMatrix = FALSE)

## This file must be run from its copy in Part12_rerun.
local({
    src <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
    want <- file.path(P12_RERUN, SCRIPT)
    if (!is.null(src) && !identical(tolower(normalizePath(src, winslash = "/", mustWork = FALSE)),
                                    tolower(normalizePath(want, winslash = "/", mustWork = FALSE))))
        stop(sprintf("run this file from Part12_rerun: save it there, then source(\"%s\")", want),
             call. = FALSE)
})

# == Helpers, identical copies of those in Part12_Fig19_Analysis.r ==========
b5_change <- function(e_removed, e) {
    s <- sign(e); s[!is.finite(s) | s == 0] <- NA_real_
    s * (e_removed - e)
}

b5_readout <- function(R, N) {
    tail_q <- as.numeric(B5$null$tail)
    alpha <- as.numeric(P12_CFG$inference$alpha)
    R$delta_sa <- b5_change(R$est_wo, R$est_all)
    DN <- if (nrow(R)) b5_change(N, matrix(R$est_all, nrow(N), ncol(N))) else N
    q <- function(v, pr) { v <- v[is.finite(v)]; if (length(v) < 20L) NA_real_ else
        stats::quantile(v, pr, names = FALSE, type = 7) }
    R$null_n <- if (nrow(R)) rowSums(is.finite(DN)) else integer(0)
    R$null_q_tail <- if (nrow(R)) apply(DN, 1L, q, tail_q) else numeric(0)
    R$null_q50 <- if (nrow(R)) apply(DN, 1L, q, 0.5) else numeric(0)
    R$null_q95 <- if (nrow(R)) apply(DN, 1L, q, 0.95) else numeric(0)
    R$null_est_lo <- R$est_all + sign(R$est_all) * R$null_q_tail
    R$null_est_hi <- R$est_all + sign(R$est_all) * R$null_q95
    R$p_null <- if (nrow(R)) vapply(seq_len(nrow(R)), function(i) {
        v <- DN[i, ]; v <- v[is.finite(v)]
        if (length(v) < 20L || !is.finite(R$delta_sa[i])) NA_real_ else
            (1 + sum(v <= R$delta_sa[i] + 1e-12)) / (1 + length(v))
    }, 0) else numeric(0)
    sig_a <- R$tested_all & is.finite(R$fdr_all) & R$fdr_all < alpha
    sig_w <- R$tested_wo & is.finite(R$fdr_wo) & R$fdr_wo < alpha
    same <- sign(R$est_all) == sign(R$est_wo)
    beyond <- is.finite(R$delta_sa) & is.finite(R$null_q_tail) & R$delta_sa < R$null_q_tail
    R$class <- ifelse(!R$tested_all, "Untestable",
               ifelse(sig_a & !R$tested_wo, "Untestable without the eleven",
               ifelse(sig_a & sig_w & same, "Holds",
               ifelse(sig_a & sig_w, "Reversed",
               ifelse(sig_a & !R$null_ok, "Lost, null untestable",
               ifelse(sig_a & beyond, "Lost, beyond random removal",
               ifelse(sig_a, "Lost, within random removal",
               ifelse(sig_w, "Gained", "Neither"))))))))
    R
}

## Parsed definitions of the objects named in `what`, for an exact comparison.
b5_defs <- function(path, what) {
    out <- list()
    if (!file.exists(path)) return(out)
    for (e in parse(path, keep.source = FALSE))
        if (is.call(e) && identical(e[[1L]], as.name("<-")) && is.name(e[[2L]]) &&
            as.character(e[[2L]]) %in% what) out[[as.character(e[[2L]])]] <- e[[3L]]
    out
}

# =============================================================================
# SECTION 31: Figure 19
# Panels 19A to 19C and the supplementary verdict cross-tabulation.
# =============================================================================
p12_section_31 <- function() {
    STEP <- "31"
    D <- p12_step_dirs(STEP, "Fig19", parent = FIG_DIR)
    p12_clean(D$out)
    p12_start(STEP, SCRIPT, D, pkgs = "ggplot2")

    PAL <- P12_CFG$palette
    INK <- PAL$neutral$dark; MID <- PAL$neutral$mid; PALE <- PAL$neutral$light
    ## Class colours for 19A and 19B: rose tones for the lost classes, blue-purple for Gained.
    KCOL <- unlist(B5$class_colours)
    KCOL[c("Reversed", "Lost, beyond random removal", "Lost, within random removal",
           "Lost, null untestable", "Gained")] <- c("#7E3B4A", "#B5566A", "#DDA3AF", "#F0D4DA", "#7474C1")
    KLEV <- B5$classes
    ARM_LABEL <- c(M4 = "M4 (primary)", PVM_3 = "PVM_3 (reference)", PVM_1 = "PVM_1", HKG8 = "HKG8 (control)")
    QLAB <- c(module_post_pct = "Module post-IP%", module_resid_tq4 = "Module post-IP%, TQR",
              mito_log1p = "%MT (log1p)", contamination = "Contamination",
              ctrl_pct_null = "Matched control, null", module = "Module post-IP%",
              mito = "%MT (log1p)", ncells = "Log cell count", depth = "Log median depth",
              module_minus_mito = "|Module rs| - |%MT rs|", module_minus_depth = "Module - depth",
              arm_minus_hkg8 = "Arm - HKG8")
    PTITLE <- c("16A" = "16A donor share of variance", "16B" = "16B split-half rs",
                "16F" = "16F rs with intronic rate, between donors")
    ## A label as a plotmath string with each standalone "rs" as r[s]; bold = TRUE for strip titles.
    b5_rs_pm <- function(x, bold = FALSE) vapply(x, function(z) {
        parts <- strsplit(z, "\\brs\\b", perl = TRUE)[[1]]
        if (grepl("\\brs$", z, perl = TRUE)) parts <- c(parts, "")
        q <- sprintf("\"%s\"", gsub("\"", "\\\\\"", parts))
        if (bold) q <- sprintf("bold(%s)", q)
        paste(q, collapse = if (bold) " * bold(r[s]) * " else " * r[s] * ")
    }, "", USE.NAMES = FALSE)
    PTITLE_X <- stats::setNames(b5_rs_pm(PTITLE, bold = TRUE), names(PTITLE))
    ## axis numbers: R's default format, with 0.00 written as 0
    b5_num <- function(b) { l <- format(b, trim = TRUE); l[is.finite(b) & abs(b) < 1e-12] <- "0"; l }
    PANELS <- c("15A", "16A", "16B", "16F", "17A", "17B", "17C", "17D", "18A")   # figure order, top to bottom in 19B
    PLAB <- c("16A" = "16A donor share", "16B" = "16B split-half", "16F" = "16F intronic, between",
              "17A" = "17A IBA1, within", "17B" = "17B CPS, between", "17C" = "17C two-predictor, within",
              "15A" = "15A clinical", "17D" = "17D two-predictor, clinical", "18A" = "18A nuclear panels")
    theme19 <- ggplot2::theme_bw(base_size = 9) +
        ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                       panel.grid.major.y = ggplot2::element_blank(),
                       strip.background = ggplot2::element_rect(fill = "white", colour = INK),
                       strip.text = ggplot2::element_text(face = "bold", colour = INK),
                       plot.title = ggplot2::element_text(face = "bold", size = 10, colour = INK),
                       plot.subtitle = ggplot2::element_text(size = 8, colour = MID),
                       plot.caption = ggplot2::element_text(size = 7, colour = MID, hjust = 0),
                       legend.position = "bottom", legend.title = ggplot2::element_blank())

    fig_dev <- function(stem, w, h, plot) {
        ok <- TRUE
        for (kind in c("pdf", "png")) {
            f <- file.path(D$out, paste0(stem, ".", kind))
            p12_output(f)
            if (kind == "pdf") grDevices::cairo_pdf(f, width = w, height = h)
            else grDevices::png(f, width = w, height = h, units = "in", res = 300,
                                type = if (.Platform$OS.type == "windows") "windows" else "cairo")
            r <- tryCatch({ print(plot); TRUE }, error = function(e) {
                p12_say(sprintf("    %s %s: %s", stem, kind, conditionMessage(e))); FALSE })
            grDevices::dev.off()
            ok <- ok && r
        }
        p12_chk(ok, sprintf("%s drawn", stem), "pdf and png")
    }
    rd <- function(folder, file) utils::read.csv(p12_upstream_file(file.path(FIG_DIR, folder), file),
                                                 check.names = FALSE, stringsAsFactors = FALSE)

    main <- function() {
        p12_say("SECTION A -- inputs")
        for (fd in c("19_SADonors", "19_UnitPanels", "19_Histology", "19_Clinical")) {
            lg <- file.path(P12_RERUN, FIG_DIR, fd, sprintf("Part12_%s_log.txt", fd))
            sl <- if (file.exists(lg)) sub("\r$", "", grep("^Script: ", readLines(lg, warn = FALSE),
                                                             value = TRUE)[1]) else NA_character_
            p12_chk(identical(sl, sprintf("Script: %s", ANALYSIS)),
                    sprintf("%s log was written by %s", fd, ANALYSIS),
                    if (is.na(sl)) "log or Script line not found" else sl)
        }
        sh <- c("b5_change", "b5_readout")
        da <- b5_defs(file.path(P12_RERUN, ANALYSIS), sh)
        dp <- b5_defs(file.path(P12_RERUN, SCRIPT), sh)
        same <- vapply(sh, function(n) !is.null(da[[n]]) && !is.null(dp[[n]]) &&
                           identical(da[[n]], dp[[n]]), TRUE)
        p12_chk(all(same), sprintf("shared helpers identical to %s", ANALYSIS),
                paste(sprintf("%s %s", sh, ifelse(same, "same", "DIFFERS or not found")), collapse = " | "))
        if (P12$n[["FAIL"]] > 0L) p12_abort("input checks failed; nothing drawn")

        RD <- list(Units = rd("19_UnitPanels", "Readout_19_Units.csv"),
                   Histology = rd("19_Histology", "Readout_19_Histology.csv"),
                   Clinical = rd("19_Clinical", "Readout_19_Clinical.csv"))
        NL <- list(Units = rd("19_UnitPanels", "Null_19_Units.csv"),
                   Histology = rd("19_Histology", "Null_19_Histology.csv"),
                   Clinical = rd("19_Clinical", "Null_19_Clinical.csv"))
        V17C <- rd("19_Histology", "Verdicts_19_17C.csv")
        V17D <- rd("19_Clinical", "Verdicts_19_17D.csv")
        V18A <- rd("19_Clinical", "Verdicts_19_18A.csv")
        FR <- rd("19_SADonors", "Frames_19.csv")
        p12_say("")

        p12_say("SECTION B -- the readouts follow from the null files")
        keep_in <- c("row_id", "panel", "arm", "level", "quantity", "measure", "frame", "k_removed",
                     "est_all", "fdr_all", "n_donors_all", "est_wo", "lo_wo", "hi_wo", "fdr_wo",
                     "n_donors_wo", "tested_all", "tested_wo", "null_ok")
        for (nm in names(RD)) {
            R <- RD[[nm]]; L <- NL[[nm]]
            nd <- max(c(0L, L$draw))
            N <- matrix(NA_real_, nrow(R), nd)
            i <- match(L$row_id, R$row_id)
            ok_map <- !anyNA(i) && !anyDuplicated(paste(L$row_id, L$draw))
            if (ok_map) N[cbind(i, L$draw)] <- L$estimate
            X <- b5_readout(R[, keep_in], N)
            num <- c("delta_sa", "null_q_tail", "null_q50", "null_q95", "p_null")
            dq <- vapply(num, function(v) {
                d <- abs(X[[v]] - R[[v]]); d[is.na(X[[v]]) & is.na(R[[v]])] <- 0
                if (any(!is.finite(d))) Inf else max(c(0, d)) }, 0)
            p12_chk(ok_map && identical(X$class, R$class) && all(dq < 1e-10) &&
                        nd == as.integer(B5$null$draws),
                    sprintf("%s: classes and null quantiles recomputed from the null file", nm),
                    sprintf("%d rows | %d draws | classes %s | largest quantile difference %.1e",
                            nrow(R), nd, if (identical(X$class, R$class)) "same" else "DIFFER", max(dq)))
            p12_chk(all(R$class %in% KLEV), sprintf("%s: every class is a listed class", nm),
                    paste(sort(unique(R$class)), collapse = ", "))
        }
        AL <- do.call(rbind, lapply(RD, function(R) R[, c("panel", "arm", "level", "quantity", "measure",
                                                          "frame", "k_removed", "est_all", "fdr_all",
                                                          "est_wo", "fdr_wo", "null_est_lo", "null_est_hi",
                                                          "class")]))
        rownames(AL) <- NULL
        p12_chk(setequal(unique(AL$panel), PANELS), "the readouts hold the nine retested panels",
                paste(sort(unique(AL$panel)), collapse = ", "))
        if (P12$n[["FAIL"]] > 0L) p12_abort("readout checks failed; nothing drawn")
        p12_say("")

        p12_say("SECTION C -- figures")
        alpha <- as.numeric(P12_CFG$inference$alpha)
        kline <- paste(sprintf("%s %s", sub("^U\\|", "", FR$frame[grepl("^U\\|", FR$frame)]),
                               FR$k[grepl("^U\\|", FR$frame)]), collapse = ", ")

        # == 19A ==============================================================
        A <- RD$Units
        A$arm_f <- factor(A$arm, levels = P12_ARMS)
        A$q_lab <- QLAB[A$quantity]
        ## One y level per panel and row; arms in arm order, quantities in QLAB order.
        A$row_key <- paste(A$panel, A$arm, A$q_lab, sep = "|")
        ord <- unlist(lapply(c("16A", "16B", "16F"), function(p) {
            g <- A[A$panel == p, ]
            g$row_key[order(match(g$arm, P12_ARMS), match(g$quantity, names(QLAB)))]
        }))
        A$row_f <- factor(A$row_key, levels = rev(ord))
        ## panel titles as plotmath so "rs" prints with a subscript s
        A$panel_f <- factor(PTITLE_X[A$panel], levels = PTITLE_X)
        sa_n <- sum(A$class %in% c("Reversed", "Lost, beyond random removal"))
        s_all <- sum(is.finite(A$fdr_all) & A$fdr_all < alpha)

        ## Dumbbells: all donors deep teal, without the eleven mid teal. Rows that
        ## change class carry a label; rows that hold, or are significant in neither run, do not.
        T13_DARK <- "#2E6E6B"; T13_MID <- "#6FB0AB"; T13_BAR <- "#C3DEDB"
        A$lab13 <- ifelse(A$class == "Lost, beyond random removal", "lost",
                   ifelse(A$class %in% c("Lost, within random removal", "Lost, null untestable"),
                          "lost, within null",
                   ifelse(A$class == "Gained", "gained",
                   ifelse(A$class == "Reversed", "reversed",
                   ifelse(A$class == "Untestable without the eleven", "untestable", "")))))
        A$bold13 <- A$lab13 %in% c("lost", "reversed")
        ## Label above the dumbbell, centred on it.
        A$xlab13 <- (A$est_all + A$est_wo) / 2
        A$xlab13[!is.finite(A$xlab13)] <- A$est_all[!is.finite(A$xlab13)]
        ## A label that would cross the panel edge is aligned to the edge.
        A$hj13 <- 0.5
        for (p in unique(A$panel)) {
            i <- which(A$panel == p)
            r <- range(c(A$est_all[i], A$est_wo[i]), na.rm = TRUE); w <- diff(r)
            half <- ifelse(A$lab13[i] == "lost, within null", 0.2, 0.08) * w
            hi <- r[2] + 0.08 * w; lo <- r[1] - 0.08 * w
            right <- A$xlab13[i] + half > hi; left <- A$xlab13[i] - half < lo
            A$hj13[i[right]] <- 1; A$xlab13[i[right]] <- hi
            A$hj13[i[left & !right]] <- 0; A$xlab13[i[left & !right]] <- lo
        }
        D13 <- rbind(data.frame(A[, c("row_f", "panel_f")], x = A$est_all, run = "All donors"),
                     data.frame(A[, c("row_f", "panel_f")], x = A$est_wo, run = "Without the eleven"))
        D13$run <- factor(D13$run, levels = c("All donors", "Without the eleven"))
        gA2 <- ggplot2::ggplot(A, ggplot2::aes(y = row_f)) +
            ggplot2::geom_segment(ggplot2::aes(x = est_all, xend = est_wo, yend = row_f), colour = T13_BAR,
                                  linewidth = 1.6, na.rm = TRUE) +
            ggplot2::geom_point(data = D13, ggplot2::aes(x = x, y = row_f, fill = run), shape = 21,
                                colour = "white", size = 3.3, stroke = 0.5, na.rm = TRUE) +
            ggplot2::geom_text(data = A[nzchar(A$lab13) & A$bold13, ],
                               ggplot2::aes(x = xlab13, y = row_f, label = lab13, hjust = hj13), vjust = 0,
                               position = ggplot2::position_nudge(y = 0.24),
                               fontface = "bold", size = 3.6, colour = T13_DARK) +
            ggplot2::geom_text(data = A[nzchar(A$lab13) & !A$bold13, ],
                               ggplot2::aes(x = xlab13, y = row_f, label = lab13, hjust = hj13), vjust = 0,
                               position = ggplot2::position_nudge(y = 0.24),
                               size = 3.1, colour = MID) +
            ggplot2::scale_fill_manual(values = c("All donors" = T13_DARK, "Without the eleven" = T13_MID)) +
            ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.1, 0.1)), labels = b5_num) +
            ggplot2::scale_y_discrete(labels = function(x) parse(text = b5_rs_pm(sub("^[^|]*\\|([^|]*)\\|", "\\1  ", x))),
                                      expand = ggplot2::expansion(add = c(0.6, 1.0))) +
            ggplot2::facet_wrap(~ panel_f, scales = "free", ncol = 3, labeller = ggplot2::label_parsed) +
            ggplot2::labs(x = "estimate", y = NULL,
                          title = "19A. Is the donor structure carried by the eleven severely affected donors?",
                          subtitle = sprintf(paste0("Rows significant with all donors: %d. Lost beyond random removal ",
                                                    "or reversed without the eleven: %d."), s_all, sa_n),
                          caption = paste0("lost: significant with all donors, not without, and the change is below ",
                                           "the 5th percentile of ", B5$null$draws, " random removals of as many ",
                                           "High-ADNC donors.\nlost, within null: lost but within that range. ",
                                           "No label: significant in both runs, or in neither. SA donors per arm: ",
                                           kline, ".")) +
            ggplot2::theme_bw(base_size = 9) +
            ggplot2::theme(panel.grid = ggplot2::element_blank(),
                           panel.grid.major.y = ggplot2::element_line(colour = "#ECECEC", linewidth = 0.35),
                           axis.text.y = ggplot2::element_text(colour = "black", size = 9),
                           axis.text.x = ggplot2::element_text(size = 9),
                           axis.ticks.y = ggplot2::element_blank(),
                           strip.background = ggplot2::element_rect(fill = "white", colour = INK),
                           strip.text = ggplot2::element_text(face = "bold", colour = INK),
                           plot.title = ggplot2::element_text(face = "bold", size = 10, colour = INK),
                           plot.subtitle = ggplot2::element_text(size = 8, colour = MID),
                           plot.caption = ggplot2::element_text(size = 7, colour = MID, hjust = 0),
                           legend.position = "bottom", legend.title = ggplot2::element_blank()) +
            ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(size = 3.6)))
        fig_dev("Fig19A", 13, 8.6, gA2)

        # == 19B ==============================================================
        CT <- stats::aggregate(list(n = rep(1L, nrow(AL))),
                               list(panel = AL$panel, arm = AL$arm, class = AL$class), sum)
        CT$panel_f <- factor(PLAB[CT$panel], levels = rev(PLAB[PANELS]))
        CT$arm_f <- factor(ARM_LABEL[CT$arm], levels = ARM_LABEL[P12_ARMS])
        CT$class_f <- factor(CT$class, levels = rev(KLEV))
        tot <- table(factor(AL$class, KLEV))
        NT <- stats::aggregate(list(n = CT$n), list(panel_f = CT$panel_f, arm_f = CT$arm_f), sum)
        ## Count inside a segment of at least LAB_MIN of its bar, otherwise just above it.
        LAB_MIN <- 0.07
        DARKK <- c("Holds", "Reversed", "Lost, beyond random removal", "Gained")
        LB <- do.call(rbind, lapply(split(CT, list(CT$panel, CT$arm), drop = TRUE), function(g) {
            g <- g[order(match(g$class, KLEV)), ]
            sh <- g$n / sum(g$n); x1 <- cumsum(sh); x0 <- x1 - sh
            g$mid <- (x0 + x1) / 2; g$share <- sh
            g$inside <- sh >= LAB_MIN
            xs <- g$mid; last <- -Inf
            for (i in which(!g$inside)) { xs[i] <- max(g$mid[i], last + 0.045); last <- xs[i] }
            g$xlab <- xs
            g
        }))
        LB$tcol <- ifelse(LB$class %in% DARKK, "white", INK)
        ## numeric row positions for the leaders (the bar half-width is 0.36)
        LB$ypos <- as.numeric(LB$panel_f)
        gB <- ggplot2::ggplot(CT, ggplot2::aes(y = panel_f, x = n, fill = class_f)) +
            ggplot2::geom_col(position = "fill", width = 0.72, colour = "black", linewidth = 0.2) +
            ggplot2::geom_text(data = LB[LB$inside, ], ggplot2::aes(y = panel_f, x = mid, label = n, colour = tcol),
                               inherit.aes = FALSE, size = 3.5) +
            ggplot2::geom_segment(data = LB[!LB$inside, ], ggplot2::aes(x = mid, xend = xlab, y = ypos + 0.36,
                                                                         yend = ypos + 0.45),
                                  inherit.aes = FALSE, colour = MID, linewidth = 0.25) +
            ggplot2::geom_text(data = LB[!LB$inside, ], ggplot2::aes(y = ypos + 0.46, x = xlab, label = n),
                               inherit.aes = FALSE, size = 1.9, colour = INK, vjust = 0) +
            ggplot2::scale_colour_identity() +
            ggplot2::geom_text(data = NT, ggplot2::aes(y = panel_f, x = 1.02, label = n), inherit.aes = FALSE,
                               hjust = 0, size = 2.3, colour = MID) +
            ggplot2::scale_fill_manual(values = KCOL, breaks = KLEV[KLEV %in% CT$class]) +
            ggplot2::scale_x_continuous(breaks = c(0, 0.5, 1), expand = ggplot2::expansion(mult = c(0, 0))) +
            ggplot2::coord_cartesian(xlim = c(0, 1.14), clip = "off") +
            ggplot2::facet_wrap(~ arm_f, nrow = 1) +
            ggplot2::labs(x = "share of retested rows (numbers: rows in each class; total at the right of each bar)", y = NULL,
                          title = "19B. Retested rows by class, with and without the eleven",
                          subtitle = paste(sprintf("%s %d", names(tot)[tot > 0], tot[tot > 0]), collapse = " | "),
                          caption = paste0("Class rules as pre-set in Part12_Preset_Readings_Battery5.md. ",
                                           "Significant: FDR < 0.05. Lost beyond random removal: attenuation ",
                                           "below the 5th percentile of random removals of as many High-ADNC donors.")) +
            theme19 + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 10.5, colour = "black"),
                                     axis.text.x = ggplot2::element_text(size = 8.5, colour = "black")) +
            ggplot2::guides(fill = ggplot2::guide_legend(nrow = 2, reverse = FALSE))
        fig_dev("Fig19B", 13, 5.21, gB)

        # == Supplementary: the verdict cross-tabulation ======================
        VK17 <- c("Module only", "%MT only", "Both", "Neither", "Unreadable", "Untestable")
        VK18 <- c("second only", "Both", "NE panel only", "Neither", "Unreadable", "Untestable")
        vt <- function(V, panel, lev) {
            tb <- as.data.frame(table(factor(V$arm, P12_ARMS), factor(V$verdict_all, lev),
                                      factor(V$verdict_wo, lev)))
            names(tb) <- c("arm", "all", "wo", "n")
            tb$panel <- panel
            tb
        }
        VT <- rbind(vt(V17C, "17C within donor", VK17), vt(V17D, "17D clinical", VK17),
                    vt(V18A, "18A nuclear panels", VK18))
        p12_chk(sum(VT$n) == nrow(V17C) + nrow(V17D) + nrow(V18A), "verdict counts cover every model",
                sprintf("%d models", sum(VT$n)))
        VT$arm_f <- factor(ARM_LABEL[as.character(VT$arm)], levels = ARM_LABEL[P12_ARMS])
        VT$panel_f <- factor(VT$panel, levels = c("17C within donor", "17D clinical", "18A nuclear panels"))
        VT$diag <- as.character(VT$all) == as.character(VT$wo)
        ## Every cell carries its count; an empty combination shows 0 in grey.
        VT$lab <- as.character(VT$n)
        ## plasma runs dark (few models) to yellow (many): white text on the darker low end
        VT$dark <- sqrt(VT$n) < 0.55 * sqrt(max(VT$n))
        VT$tcol <- ifelse(VT$n == 0, "#A6A6A6", ifelse(VT$dark, "white", INK))
        VT$tsize <- ifelse(VT$n == 0, 2.0, 3.4)   # zeros in small grey type
        VT$all <- factor(VT$all, levels = rev(unique(c(VK17, VK18))))
        VT$wo <- factor(VT$wo, levels = unique(c(VK17, VK18)))
        moved <- sum(VT$n[!VT$diag])
        gC <- ggplot2::ggplot(VT, ggplot2::aes(x = wo, y = all)) +
            ggplot2::geom_tile(ggplot2::aes(fill = ifelse(n > 0, n, NA)), colour = "black", linewidth = 0.2) +
            ggplot2::geom_tile(data = VT[VT$diag & VT$n > 0, ], fill = NA, colour = "black", linewidth = 0.9) +
            ggplot2::geom_text(ggplot2::aes(label = lab, colour = tcol, size = tsize), show.legend = FALSE) +
            ggplot2::scale_size_identity() +
            ggplot2::scale_fill_viridis_c(option = "plasma", trans = "sqrt", na.value = "white",
                                          name = "models") +
            ggplot2::scale_colour_identity() +
            ## arm and panel on two strip lines, so the titles fit the narrower figure
            ggplot2::facet_wrap(~ arm_f + panel_f, ncol = 3, scales = "free") +
            ggplot2::labs(x = "verdict without the eleven", y = "verdict with all donors",
                          title = "Do the Battery 3 and 4 verdicts change without the eleven?",
                          subtitle = sprintf("Models whose verdict changes: %d of %d.", moved, sum(VT$n)),
                          caption = paste0("Numbers are models; 0 marks a combination with no model. Thick outline: verdict unchanged. Verdict rules of 17C, 17D and 18A; ",
                                           "FDR recomputed without the eleven within the same families.")) +
            theme19 + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 9.5,
                                                                          colour = "black"),
                                     axis.text.y = ggplot2::element_text(size = 9.5, colour = "black"),
                                     axis.title = ggplot2::element_text(colour = "black"),
                                     strip.text = ggplot2::element_text(face = "bold", colour = "black", size = 9.5),
                                     panel.grid = ggplot2::element_blank(),
                                     legend.position = "right", legend.title = ggplot2::element_text())
        fig_dev("FigS19_verdicts", 9.75, 10, gC)

        # == 19C: FDR of the clinical association =============================
        ## The 15A unadjusted Module coefficient per arm, region set and outcome.
        OC <- P12_CFG$clinical$outcomes
        OUTS <- names(OC)
        OCL <- vapply(OC, function(o) o$class, "")
        CQ <- RD$Clinical[RD$Clinical$panel == "15A" & RD$Clinical$quantity == "module_beta_unadjusted", ]
        p12_chk(nrow(CQ) == length(P12_ARMS) * length(P12_DATASETS) * length(OUTS) && all(CQ$measure %in% OUTS),
                "19C holds one unadjusted 15A row per arm, region set and outcome",
                sprintf("%d rows", nrow(CQ)))
        RUNS <- c("All donors", "Without the eleven")
        FL <- rbind(data.frame(CQ[, c("arm", "level", "measure", "class")], run = RUNS[1], fdr = CQ$fdr_all),
                    data.frame(CQ[, c("arm", "level", "measure", "class")], run = RUNS[2], fdr = CQ$fdr_wo))
        FL$y <- -log10(pmax(FL$fdr, 1e-300))
        FL$x <- factor(FL$measure, levels = OUTS)
        FL$arm_f <- factor(ARM_LABEL[FL$arm], levels = ARM_LABEL[P12_ARMS])
        FL$run <- factor(FL$run, levels = RUNS)
        sig_a <- sum(CQ$fdr_all < alpha, na.rm = TRUE)
        sig_b <- sum(CQ$fdr_all < alpha & CQ$fdr_wo < alpha, na.rm = TRUE)
        brk <- cumsum(table(factor(OCL, unique(OCL))))
        brk <- brk[-length(brk)] + 0.5
        gC <- ggplot2::ggplot(FL, ggplot2::aes(x = x, y = y, group = run, colour = run, linetype = run)) +
            ## declared first so the separators sit on the outcome axis; short outcome names
            ggplot2::scale_x_discrete(labels = function(x) {
                sh <- c("Cognitive status" = "CS", "MMSE score" = "MMSE", "CASI score" = "CASI", "CERAD score" = "CERAD",
                        "Braak stage" = "Braak", "Thal phase" = "Thal", "ADNC rank" = "ADNC",
                        "CAA score" = "CAA", "APOE e4 dose" = "APOE")
                ifelse(x %in% names(sh), sh[x], x) }) +
            ggplot2::geom_vline(xintercept = brk, colour = PALE, linewidth = 0.35, linetype = "dotted") +
            ggplot2::geom_hline(yintercept = -log10(alpha), colour = MID, linewidth = 0.45, linetype = "dashed") +
            ggplot2::geom_line(linewidth = 0.7, na.rm = TRUE) +
            ggplot2::geom_point(ggplot2::aes(shape = run), size = 2.4, stroke = 0.8, na.rm = TRUE) +
            ggplot2::scale_colour_manual(values = stats::setNames(c(INK, PAL$neutral$accent), RUNS)) +
            ggplot2::scale_linetype_manual(values = stats::setNames(c("solid", "22"), RUNS)) +
            ## both runs as filled circles (an open circle would read as not significant)
            ggplot2::scale_shape_manual(values = stats::setNames(c(16, 16), RUNS)) +
            ggplot2::facet_grid(level ~ arm_f) +
            ggplot2::labs(x = NULL, y = "-log10 FDR",
                          title = "19C. Does the clinical association hold without the eleven?",
                          subtitle = sprintf(paste0("Module post-IP%% coefficient, cov4 (Figure 12A model). ",
                                                    "Rows significant with all donors: %d; still significant ",
                                                    "without the eleven: %d."), sig_a, sig_b),
                          caption = paste0("FDR is BH within arm and region set over the nine outcomes, as in 15A. ",
                                           "Dashed line: FDR 0.05; points above it are significant. Dotted lines ",
                                           "separate cognition, pathology and genotype. CS: cognitive status.")) +
            theme19 + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 10,
                                                                          colour = "black"),
                                     axis.text.y = ggplot2::element_text(size = 9, colour = "black"),
                                     panel.grid.major.x = ggplot2::element_blank(),
                                     panel.grid.major.y = ggplot2::element_line(colour = "#ECECEC", linewidth = 0.35))
        fig_dev("Fig19C", 11.7, 5.85, gC)
        p12_say("")

        p12_say("SECTION D -- result lines")
        p12_say(sprintf("  19A: significant with all donors %d | reversed or lost beyond random removal %d", s_all, sa_n))
        for (p in PANELS) {
            g <- AL[AL$panel == p, ]
            tb <- table(factor(g$class, KLEV))
            p12_say(sprintf("  %-4s %s", p, paste(sprintf("%s %d", names(tb)[tb > 0], tb[tb > 0]), collapse = " | ")))
        }
        p12_say(sprintf("  19C: 15A unadjusted rows significant with all donors %d, still significant without the eleven %d",
                        sig_a, sig_b))
        for (a in P12_ARMS) for (ds in P12_DATASETS) {
            g <- CQ[CQ$arm == a & CQ$level == ds, ]
            p12_say(sprintf("    %-6s %-8s %s", a, ds, paste(sprintf("%s %s -> %s", g$measure, p12_fmt(g$fdr_all),
                                                                   p12_fmt(g$fdr_wo)), collapse = " | ")))
        }
        p12_say(sprintf("  Supplementary verdicts: changed %d of %d", moved, sum(VT$n)))
    }

    p12_run(main)
}

# =============================================================================
# RUN
# =============================================================================
local({
    v <- p12_section_31()
    if (!identical(v, "PASS"))
        stop(sprintf("section 31 ended with verdict %s. See Part12_Fig19_log.txt.", v), call. = FALSE)
    message("Figure 19: section 31 PASS. Figures in Part12_rerun/11_Fig19/Fig19")
})
