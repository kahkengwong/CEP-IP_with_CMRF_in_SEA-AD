# =============================================================================
# Part12_Fig15_Plots.r
# SEA-AD CEP-IP/CMRF, Part 12, step 12: Figure 15 plots (Figure 15, Supplementary Figure 10)
#   Draws Figure 15 (panels A to H) and Supplementary Figure 10 from the stats
#   and plot files written by steps 07 to 11 of Part12_Fig15_Analysis.r. Each
#   input is checked against its step log. No statistic is computed here.
#
# Usage (from the SEA-AD folder, or set SEAAD_ROOT), after Part12_Fig15_Analysis.r
#   source("Part12_rerun/Part12_Fig15_Plots.r")
#
# Outputs (Part12_rerun/07_Fig15/Fig15/): figures (PDF, PNG), log.
# =============================================================================

source(file.path(Sys.getenv("SEAAD_ROOT", getwd()), "Part12_rerun", "Part12_03_Config.r"))

SCRIPT <- "Part12_Fig15_Plots.r"
FIG_DIR <- "07_Fig15"   # every Figure 15 section folder sits inside this one

# =============================================================================
# SECTION 12: Figure 15 (panels A to H) and Supplementary Figure 10
#   15A  dR2 before against after the contamination covariate
#   15B  the same, one panel per doublet threshold
#   15C  within-donor contamination contrasts by cell group
#   15D  change in rs from the cleanest contamination tertile
#   15E  ambient against nuclear MT fraction per sample
#   15F  raw and partial rs per DFG gene
#   15G  rs by soup-index tertile
#   15H  depth-thinning curve with the Gerrits band
#   S10  post-IP minus pre-IP per metric, per arm
# Arms are drawn in the order M4, PVM_3, PVM_1, HKG8. Significance: filled
# symbol and * at FDR < 0.05, ** at FDR < 0.01.
# =============================================================================
p12_section_12 <- function() {
    STEP <- "12"
    D <- p12_step_dirs(STEP, "Fig15", parent = FIG_DIR)
    p12_clean(c(D$tables, D$figures))
    p12_start(STEP, SCRIPT, D)

    PAL <- P12_CFG$palette
    ARM_COL <- unlist(PAL$arms)
    G_DARK <- PAL$neutral$dark; G_MID <- PAL$neutral$mid; G_LIGHT <- PAL$neutral$light
    G_RULE <- PAL$neutral$rule; ACCENT <- PAL$neutral$accent; BAND <- PAL$neutral$band
    ARMS <- P12_ARMS

    # Significance mark from an FDR.
    stars <- function(q) ifelse(!is.finite(q), "", ifelse(q < 0.01, "**", ifelse(q < 0.05, "*", "")))
    # An FDR from a p value at the resampling floor is a bound, written with "<=".
    fdr_txt <- function(q, p = NULL) {
        out <- p12_fmt(q, p12_boot_floor())
        if (is.null(p)) return(out)
        at <- is.finite(p) & is.finite(q) & p <= p12_boot_floor() * (1 + 1e-9)
        ifelse(at & !startsWith(out, "<="), paste0("<=", sprintf("%.2e", q)), out)
    }
    n_donors_txt <- function(n) sprintf("%d donor%s", n, ifelse(n == 1, "", "s"))
    # Greedy word wrap to a width in inches on the open device; one element of
    # `txt` per paragraph.
    wrap_in <- function(txt, width_in, cex) {
        out <- character(0)
        for (para in txt) {
            cur <- ""
            for (w in strsplit(para, " ", fixed = TRUE)[[1]]) {
                nxt <- if (nzchar(cur)) paste(cur, w) else w
                if (nzchar(cur) && strwidth(nxt, units = "inches", cex = cex) > width_in) {
                    out <- c(out, cur); cur <- w
                } else cur <- nxt
            }
            out <- c(out, cur)
        }
        out
    }

    # Four significance states for the before/after plots.
    STATE_PCH <- c(both = 21, lost = 4, gained = 24, never = 1)
    STATE_LAB <- c(both = "significant before and after", lost = "lost significance",
                   gained = "gained significance", never = "never significant")
    state_of <- function(fb, fa) {
        sb <- is.finite(fb) & fb < 0.05; sa <- is.finite(fa) & fa < 0.05
        ifelse(sb & sa, "both", ifelse(sb, "lost", ifelse(sa, "gained", "never")))
    }
    state_points <- function(x, y, st, arm, cex) {
        pchv <- unname(STATE_PCH[st])
        points(x, y, pch = pchv, cex = cex, lwd = 1.4, bg = ARM_COL[arm],
               col = ifelse(pchv == 4, ARM_COL[arm], "black"))
    }
    state_legend <- function(st, where = "topleft", cex = 0.72) {
        k <- names(STATE_PCH)
        n <- vapply(k, function(z) sum(st == z), integer(1))
        legend(where, bty = "n", cex = cex, pt.cex = 1.2, pch = unname(STATE_PCH[k]),
               pt.bg = "grey70", col = "black", legend = sprintf("%s (n = %d)", STATE_LAB[k], n))
    }
    arm_legend <- function(where = "bottomright", cex = 0.78) {
        legend(where, bty = "n", cex = cex, pt.cex = 1.3, pch = 21, pt.bg = ARM_COL[ARMS],
               col = "black", legend = ARMS)
    }

    # Writes one figure as PDF and PNG. An error while drawing is a FAIL.
    fig_dev <- function(stem, w, h, draw) {
        ok <- TRUE
        for (kind in c("pdf", "png")) {
            f <- file.path(D$figures, paste0(stem, ".", kind))
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

    main <- function() {
        p12_say("SECTION A -- inputs")
        A <- rd("15AB_ContamDoublets", "Stats_15A.csv")
        B <- rd("15AB_ContamDoublets", "Stats_15B.csv")
        C <- rd("15CD_DecontXStrata", "Stats_15C.csv")
        Dd <- rd("15CD_DecontXStrata", "Stats_15D.csv")
        Dp <- rd("15CD_DecontXStrata", "Plot_15D_tertiles.csv")
        E <- rd("15EG_GerritsAmbient", "Stats_15E.csv")
        Fs <- rd("15EG_GerritsAmbient", "Stats_15F.csv")
        G <- rd("15EG_GerritsAmbient", "Stats_15G.csv")
        Pe <- rd("15EG_GerritsAmbient", "Plot_15E_samples.csv")
        Pf <- rd("15EG_GerritsAmbient", "Plot_15F_samples.csv")
        Pg <- rd("15EG_GerritsAmbient", "Plot_15G_tertiles.csv")
        H <- rd("15H_DepthThinning", "Stats_15H.csv")
        S10 <- rd("S10_Cargo", "Stats_S10.csv")
        lv <- unique(c(A$labels_version, B$labels_version, C$labels_version, Dd$labels_version,
                       H$labels_version, S10$labels_version))
        p12_chk(identical(lv, P12_CFG$labels$labels_version),
                "every stats file carries this labels version", paste(lv, collapse = ", "))
        PART <- P12_CFG$genes$cdfg
        p12_say("")

        p12_say("SECTION B -- figures")

        # == 15A =============================================================
        key <- function(X) paste(X$arm, X$level, X$measure)
        ub <- A[A$quantity == "module_beta_unadjusted", ]
        ab <- A[A$quantity == "module_beta_contam_adjusted", ]
        ux <- A[A$quantity == "dR2_unadjusted_same_donors", ]
        ay <- A[A$quantity == "dR2_contam_adjusted", ]
        d15a <- data.frame(arm = ux$arm, x = ux$estimate, y = ay$estimate[match(key(ux), key(ay))],
                           fb = ub$fdr[match(key(ux), key(ub))], fa = ab$fdr[match(key(ux), key(ab))])
        d15a <- d15a[is.finite(d15a$x) & is.finite(d15a$y), ]
        d15a <- d15a[order(match(d15a$arm, ARMS)), ]
        fig_dev("Fig15A_contamination_dR2", 7.2, 6.6, function() {
            op <- par(mar = c(4.8, 5.0, 4.4, 1.4), family = "sans"); on.exit(par(op))
            st <- state_of(d15a$fb, d15a$fa)
            lim <- c(0, max(c(d15a$x, d15a$y)) * 1.08)
            plot(NA, xlim = lim, ylim = lim, las = 1, font.lab = 2,
                 xlab = "Module post-IP% dR2 before contamination is added",
                 ylab = "Module post-IP% dR2 after contamination is added")
            abline(0, 1, col = "grey35", lwd = 1.6, lty = 2)
            grid(col = "grey92", lty = 1)
            state_points(d15a$x, d15a$y, st, d15a$arm, 1.3)
            mtext("Module post-IP% dR2 before and after adjusting for donor-mean contamination",
                  side = 3, line = 1.8, font = 2, cex = 1.0)
            mtext(sprintf(paste0("supertype non-Module cells | Figure 12A patients | %d rows | ",
                                 "%d of %d significant rows lost"),
                          nrow(d15a), sum(st == "lost"), sum(st %in% c("both", "lost"))),
                  side = 3, line = 0.4, cex = 0.7, col = "grey30")
            state_legend(st); arm_legend()
        })

        # == 15B =============================================================
        ud <- A[A$quantity == "dR2_unadjusted", ]
        cuts <- sort(unique(as.numeric(sub("^module_beta_doublet_", "",
                                           B$quantity[grepl("^module_beta_doublet_", B$quantity)]))))
        fig_dev("Fig15B_doublet_dR2", 4.4 * length(cuts) + 0.8, 5.2, function() {
            op <- par(mfrow = c(1, length(cuts)), mar = c(4.6, 4.6, 4.0, 1.0), oma = c(0, 0, 3.0, 0),
                      family = "sans"); on.exit(par(op))
            allv <- c(ud$estimate, B$estimate[grepl("^dR2_doublet_", B$quantity)])
            lim <- c(0, max(allv, na.rm = TRUE) * 1.08)
            for (cu in cuts) {
                tag <- sprintf("%.2f", cu)
                bb <- B[B$quantity == paste0("module_beta_doublet_", tag), ]
                by <- B[B$quantity == paste0("dR2_doublet_", tag), ]
                d <- data.frame(arm = ud$arm, x = ud$estimate, y = by$estimate[match(key(ud), key(by))],
                                fb = ub$fdr[match(key(ud), key(ub))],
                                fa = bb$fdr[match(key(ud), key(bb))])
                d <- d[is.finite(d$x) & is.finite(d$y), ]
                d <- d[order(match(d$arm, ARMS)), ]
                st <- state_of(d$fb, d$fa)
                plot(NA, xlim = lim, ylim = lim, las = 1, font.lab = 2,
                     xlab = "Module post-IP% dR2, all nuclei",
                     ylab = if (cu == cuts[1]) "Module post-IP% dR2, doublet-like nuclei removed" else "")
                abline(0, 1, col = "grey35", lwd = 1.6, lty = 2)
                grid(col = "grey92", lty = 1)
                state_points(d$x, d$y, st, d$arm, 1.5)
                mtext(sprintf("threshold %s", tag), side = 3, line = 1.5, font = 2, cex = 0.95)
                mtext(sprintf("%d rows | %d of %d significant rows lost", nrow(d), sum(st == "lost"),
                              sum(st %in% c("both", "lost"))), side = 3, line = 0.3, cex = 0.68,
                      col = "grey30")
                state_legend(st, cex = 0.68)
                if (cu == cuts[1]) arm_legend()
            }
            mtext("Module post-IP% dR2 before and after removing doublet-like nuclei, Figure 12A patients",
                  side = 3, outer = TRUE, line = 1.0, font = 2, cex = 1.1)
        })

        # == 15C =============================================================
        CON <- c(post_minus_pre = "Module post-IP vs Module pre-IP",
                 post_minus_nonmodule = "Module post-IP vs non-Module",
                 pre_minus_nonmodule = "Module pre-IP vs non-Module")
        cc <- C[C$quantity %in% names(CON), ]
        cc <- cc[order(match(cc$arm, ARMS), match(cc$quantity, names(CON))), ]
        at_floor <- startsWith(fdr_txt(cc$fdr, cc$p), "<=")
        p12_chk(all(cc$fdr[at_floor] < 0.001), "15C FDR values printed as < 0.001 are below 0.001",
                sprintf("%d floor-derived values, largest %s", sum(at_floor),
                        if (any(at_floor)) sprintf("%.2e", max(cc$fdr[at_floor])) else "none"))
        fig_dev("Fig15C_contamination_by_group", 12.0, 6.14, function() {
            layout(matrix(c(1, 2), nrow = 2), heights = c(6.2, 1.0))
            op <- par(mar = c(5.4, 21.0, 4.0, 13.0), xpd = FALSE); on.exit({ par(op); layout(1) })
            ypos <- rev(seq_len(nrow(cc)))
            xr <- range(c(0, cc$estimate, cc$lo, cc$hi), na.rm = TRUE)
            xpad <- max(0.01, diff(xr) * 0.10)
            xtk <- pretty(c(xr[1] - xpad, xr[2] + xpad))
            plot(NA, xlim = range(xtk), ylim = c(0.55, nrow(cc) + 0.45),
                 axes = FALSE, xlab = "", ylab = "", yaxs = "i")
            brk <- which(cc$arm[-1] != cc$arm[-nrow(cc)])
            if (length(brk)) abline(h = ypos[brk] - 0.5, col = G_RULE, lty = 3)
            abline(v = 0, col = G_DARK, lwd = 1.2)
            ok <- is.finite(cc$estimate)
            col <- ARM_COL[cc$arm]
            if (any(ok)) {
                # One line per row: the 95% interval, drawn at a single width.
                segments(cc$lo[ok], ypos[ok], cc$hi[ok], ypos[ok], col = col[ok], lwd = 2)
                sig <- is.finite(cc$fdr) & cc$fdr < 0.05
                points(cc$estimate[ok], ypos[ok], pch = 21, cex = 1.8, lwd = 1.4,
                       bg = ifelse(sig[ok], col[ok], "white"), col = G_DARK)
            }
            if (any(!ok)) text(0, ypos[!ok], "untestable", cex = 0.8, col = G_MID, adj = c(-0.15, 0.5))
            axis(1, at = xtk, col = G_DARK, col.axis = G_DARK, cex.axis = 0.95)
            axis(2, at = ypos, labels = CON[cc$quantity], las = 1, tick = FALSE, line = -0.5,
                 cex.axis = 0.92)
            x_arm <- grconvertX(0.015, from = "ndc", to = "user")
            for (a in unique(cc$arm))
                text(x_arm, mean(ypos[cc$arm == a]), labels = a, adj = 0, font = 2, cex = 0.98,
                     col = G_DARK, xpd = NA)
            axis(4, at = ypos, las = 1, tick = FALSE, line = -0.5, cex.axis = 0.88, col.axis = G_DARK,
                 # A floor-derived FDR is a bound below 0.001.
                 labels = ifelse(ok, sprintf("%s   FDR %s %s", n_donors_txt(cc$n_donors),
                                             ifelse(startsWith(fdr_txt(cc$fdr, cc$p), "<="),
                                                    "< 0.001",
                                                    paste("=", fdr_txt(cc$fdr, cc$p))),
                                             stars(cc$fdr)),
                                 sprintf("%s   untestable", n_donors_txt(cc$n_donors))))
            mtext("median over donors of the within-donor difference in DecontX contamination",
                  side = 1, line = 2.6, cex = 0.98)
            mtext("<- first group cleaner        first group dirtier ->", side = 1, line = 3.9,
                  cex = 0.86, col = G_MID)
            mtext("Contamination by cell group, within donor", side = 3, line = 2.2, adj = 0,
                  font = 2, cex = 1.12)
            mtext(paste("NeoC and NonNeoC pooled within donor; bars are 95% donor-bootstrap",
                        "intervals; BH within arm over the three contrasts."),
                  side = 3, line = 1.1, adj = 0, cex = 0.80, col = G_MID)
            mtext(sprintf("A donor enters a contrast with at least %d cells in both groups; an arm with fewer than %d donors is untestable.",
                          as.integer(P12_CFG$floors$cells_per_group_per_donor),
                          as.integer(P12_CFG$floors$donors_per_test)),
                  side = 3, line = 0.2, adj = 0, cex = 0.80, col = G_MID)
            par(mar = c(0, 0, 0, 0)); plot.new()
            legend("left", legend = ARMS, pch = 15, col = ARM_COL[ARMS], pt.cex = 1.8,
                   horiz = TRUE, bty = "n", cex = 0.92, text.col = G_DARK, inset = c(0.02, 0))
            legend("right", legend = c("FDR < 0.05", "ns"), pch = 21, pt.bg = c(G_DARK, "white"),
                   col = G_DARK, pt.cex = 1.6, horiz = TRUE, bty = "n", cex = 0.92,
                   text.col = G_DARK, inset = c(0.02, 0))
        })

        # == 15D =============================================================
        dt <- Dd[Dd$quantity == "rs_T3_minus_T1", ]
        Dp <- Dp[order(match(Dp$arm, ARMS), Dp$tertile), ]
        Dp$d_rs <- Dp$rs - Dp$rs[match(paste(Dp$arm, 1L), paste(Dp$arm, Dp$tertile))]
        # The drawn T3 - T1 must equal the tested estimate.
        dchk <- vapply(dt$arm[is.finite(dt$estimate)], function(a)
            abs(Dp$d_rs[Dp$arm == a & Dp$tertile == 3] - dt$estimate[dt$arm == a]), 0)
        p12_chk(all(dchk < 1e-9), "15D plot data agree with the stats file",
                sprintf("max |difference| %.2g over %d arms", max(c(0, dchk)), length(dchk)))
        fig_dev("Fig15D_rs_by_contamination", 6.3, 6.4, function() {
            layout(matrix(c(1, 2), nrow = 2), heights = c(6.2, 1.1))
            x_hi <- max(Dp$contamination_median, na.rm = TRUE)
            xlim <- c(-0.035 * x_hi, x_hi * 1.28)
            yr <- range(c(Dp$d_rs, dt$lo, dt$hi, 0), na.rm = TRUE)
            ypad <- diff(yr) * 0.10
            # Display offset only: arms within 2.4% of the x range are spread sideways.
            t3x <- vapply(ARMS, function(a) {
                v <- Dp$contamination_median[Dp$arm == a & Dp$tertile == 3]
                if (length(v) == 1L) v else NA_real_ }, 0)
            dw <- 0.012 * diff(xlim)
            dodge <- setNames(rep(0, length(ARMS)), ARMS)
            o <- ARMS[order(t3x)]; o <- o[is.finite(t3x[o])]
            cl <- cumsum(c(1, diff(t3x[o]) >= 2 * dw))
            for (g in unique(cl)) {
                m <- o[cl == g]
                if (length(m) > 1L) dodge[m] <- (seq_along(m) - (length(m) + 1) / 2) * 2 * dw
            }
            # Title and subtitle start at the left edge and wrap to the figure width.
            SUB_CEX <- 0.78
            subs <- wrap_in(c(
                "Each series starts at zero at its own cleanest tertile; an ambient-driven correlation would rise.",
                "Tertiles cut within donor, NeoC and NonNeoC pooled. GOI MT-CO3 (ACTB on HKG8).",
                "Bars are 95% donor-bootstrap intervals on T3 minus T1; * FDR < 0.05, ** FDR < 0.01.",
                paste0("BH over the four arms.",
                       if (any(dodge != 0)) " Coincident T3 points are offset sideways." else "")),
                width_in = 0.94 * par("din")[1], cex = SUB_CEX)
            ns <- length(subs)
            op <- par(mar = c(4.6, 5.4, 1.5 + 0.95 * ns + 0.3, 1.2)); on.exit({ par(op); layout(1) })
            plot(NA, xlim = xlim, ylim = c(yr[1] - ypad, yr[2] + ypad), axes = FALSE, xlab = "",
                 ylab = "", xaxs = "i", yaxs = "i")
            abline(h = 0, col = G_DARK, lwd = 1.2)
            axis(1, col = G_DARK, col.axis = G_DARK); axis(2, las = 1, col = G_DARK, col.axis = G_DARK)
            box(col = G_DARK)
            for (a in ARMS) {
                s <- Dp[Dp$arm == a, ]
                r <- dt[dt$arm == a, ]
                if (!nrow(s) || !nrow(r)) next
                testable <- is.finite(r$estimate)
                sx <- s$contamination_median + ifelse(s$tertile == 3, dodge[[a]], 0)
                lines(sx, s$d_rs, col = ARM_COL[a], lwd = 2.8, lty = if (testable) 1 else 2)
                xe <- sx[s$tertile == 3]
                if (testable) arrows(xe, r$lo, xe, r$hi, code = 3, angle = 90, length = 0.035,
                                     col = ARM_COL[a], lwd = 1.6)
                points(sx[s$tertile < 3], s$d_rs[s$tertile < 3], pch = 19, cex = 1.0,
                       col = ARM_COL[a])
                points(xe, s$d_rs[s$tertile == 3], pch = 21, cex = 1.9, lwd = 1.4,
                       bg = if (testable && r$fdr < 0.05) ARM_COL[a] else "white", col = G_DARK)
                # The left member of an offset pair is labelled below and to the left.
                lab_t3 <- if (testable) sprintf("%+.3f %s", r$estimate, stars(r$fdr)) else "untestable"
                if (dodge[[a]] < 0)
                    text(xe, s$d_rs[s$tertile == 3], labels = lab_t3, adj = c(1.12, 1.6), cex = 0.82,
                         font = 2, col = ARM_COL[a])
                else
                    text(xe, s$d_rs[s$tertile == 3], labels = lab_t3, pos = 4, offset = 0.6, cex = 0.82,
                         font = 2, col = ARM_COL[a])
            }
            mtext("median DecontX contamination of the tertile", side = 1, line = 2.8, cex = 0.95)
            # Two lines, so the label stays within the plot height at this width.
            mtext("change in rs(GOI, composite)", side = 2, line = 4.4, cex = 0.95)
            mtext("from the cleanest tertile", side = 2, line = 3.4, cex = 0.95)
            x_left <- grconvertX(0.03, from = "ndc", to = "user")
            mtext("Co-expression across contamination tertiles of the Module cells", side = 3,
                  line = 0.3 + 0.95 * ns + 0.25, at = x_left, adj = 0, font = 2, cex = 1.0)
            for (i in seq_len(ns))
                mtext(subs[i], side = 3, line = 0.3 + 0.95 * (ns - i), at = x_left, adj = 0,
                      cex = SUB_CEX, col = G_MID)
            # Key in two rows: arms, then the T3 significance symbols.
            par(mar = c(0, 0, 0, 0)); plot.new()
            legend(0.5, 0.98, xjust = 0.5, yjust = 1, legend = ARMS, col = ARM_COL[ARMS], lty = 1,
                   lwd = 2.8, pch = 21, pt.bg = ARM_COL[ARMS], pt.cex = 1.3, horiz = TRUE, bty = "n",
                   cex = 0.88, text.col = G_DARK, seg.len = 2.0,
                   text.width = max(strwidth(ARMS, cex = 0.88)) * 1.10)
            legend(0.5, 0.50, xjust = 0.5, yjust = 1, legend = c("T3, FDR < 0.05", "T3, ns"),
                   pch = 21, pt.bg = c(G_DARK, "white"), pt.cex = 1.3, col = G_DARK, horiz = TRUE,
                   bty = "n", cex = 0.88, text.col = G_DARK)
        })

        # == 15E =============================================================
        e1 <- E[E$quantity == "log_ambient_over_nuclear_mito", ]
        fig_dev("Fig15E_ambient_vs_nuclear", 6.4, 6.2, function() {
            op <- par(mar = c(4.8, 5.0, 4.6, 1.2), family = "sans"); on.exit(par(op))
            x <- 100 * Pe$ambient_mito_frac; y <- 100 * Pe$nuclear_mito_frac
            lim <- c(0, max(c(x, y), na.rm = TRUE) * 1.06)
            plot(NA, xlim = lim, ylim = lim, las = 1, bty = "l", font.lab = 2,
                 xlab = "ambient MT fraction, % (barcodes of 1 to 100 UMIs)",
                 ylab = "nuclear MT fraction, % (called nuclei)")
            abline(0, 1, lty = 2, lwd = 1.6, col = G_DARK)
            points(x, y, pch = 21, cex = 1.6, bg = G_MID, col = "black", lwd = 0.6)
            legend("topleft", bty = "n", cex = 0.8, lty = 2, lwd = 1.6, col = G_DARK,
                   legend = "equal ambient and nuclear fraction")
            mtext("Nuclear against ambient mitochondrial fraction (Gerrits)", side = 3, line = 2.0,
                  font = 2, cex = 1.05)
            mtext(sprintf("%d samples | %d below the line | median ambient:nuclear ratio %.2fx | FDR %s %s",
                          nrow(Pe), sum(y < x, na.rm = TRUE), exp(e1$estimate), fdr_txt(e1$fdr),
                          stars(e1$fdr)),
                  side = 3, line = 0.6, cex = 0.72, col = "grey40")
        })

        # == 15F =============================================================
        f_raw <- Fs[Fs$quantity == "rs_raw" & Fs$measure == "median over 8 partners", ]
        e_med <- stats::median(Pe$log_ratio, na.rm = TRUE)
        p12_chk(abs(e_med - e1$estimate) < 1e-9, "15E plot data agree with the stats file",
                sprintf("median log ratio %.6f against %.6f", e_med, e1$estimate))
        f_par <- Fs[Fs$quantity == "rs_partial" & Fs$measure == "median over 8 partners", ]
        fig_dev("Fig15F_raw_vs_partial", 8.0, 5.6, function() {
            op <- par(mar = c(6.0, 5.0, 5.0, 1.2), family = "sans"); on.exit(par(op))
            at <- seq_along(PART); w <- 0.18
            raw <- lapply(PART, function(g) Pf$rs_raw[Pf$partner == g & is.finite(Pf$rs_raw)])
            par_ <- lapply(PART, function(g) Pf$rs_partial[Pf$partner == g & is.finite(Pf$rs_partial)])
            yr <- range(c(0, unlist(raw), unlist(par_)), na.rm = TRUE)
            plot(NA, xlim = c(0.5, length(PART) + 0.5), ylim = yr, axes = FALSE, xlab = "", ylab = "")
            abline(h = 0, col = "grey", lwd = 0.8)
            boxplot(raw, at = at - w, boxwex = 0.30, add = TRUE, axes = FALSE, outline = FALSE,
                    col = G_LIGHT, border = G_DARK, whisklty = 1, medcol = G_DARK, medlwd = 1.6)
            boxplot(par_, at = at + w, boxwex = 0.30, add = TRUE, axes = FALSE, outline = FALSE,
                    col = ACCENT, border = G_DARK, whisklty = 1, medcol = "white", medlwd = 1.6)
            axis(1, at = at, labels = PART, las = 2, cex.axis = 0.85)
            axis(2, las = 1)
            box(bty = "l")
            mtext("rs of DFG genes with MT-CO3", side = 2, line = 3.4, font = 2)
            mtext("Raw and partial co-expression per DFG gene (Gerrits)", side = 3, line = 3.0,
                  font = 2, cex = 1.05)
            mtext(paste("Raw: MT-CO3-positive called nuclei. Partial: soup index and log depth",
                        "removed. One value per sample."), side = 3, line = 1.8, cex = 0.72,
                  col = "grey40")
            mtext(sprintf("Median over samples: raw %.3f, FDR %s %s | partial %.3f, FDR %s %s | partial / raw %.2f",
                          f_raw$estimate, fdr_txt(f_raw$fdr), stars(f_raw$fdr), f_par$estimate,
                          fdr_txt(f_par$fdr), stars(f_par$fdr),
                          Fs$estimate[Fs$quantity == "partial_over_raw"]),
                  side = 3, line = 0.6, cex = 0.72, col = "grey40")
            legend("bottomright", bty = "n", cex = 0.8, fill = c(G_LIGHT, ACCENT), border = G_DARK,
                   legend = c("raw", "partial"))
        })

        # == 15G =============================================================
        g1 <- G[G$quantity == "rs_top_minus_bottom" & G$measure == "median over 8 partners", ]
        fig_dev("Fig15G_rs_by_soup_tertile", 5.28, 5.6, function() {
            op <- par(mar = c(5.0, 5.0, 4.6, 1.2), family = "sans"); on.exit(par(op))
            vals <- lapply(1:3, function(t) Pg$rs[Pg$tertile == t & is.finite(Pg$rs)])
            yr <- range(c(0, unlist(vals)), na.rm = TRUE)
            plot(NA, xlim = c(0.5, 3.5), ylim = yr, axes = FALSE, xlab = "", ylab = "")
            abline(h = 0, col = "grey", lwd = 0.8)
            boxplot(vals, at = 1:3, boxwex = 0.5, add = TRUE, axes = FALSE, outline = FALSE,
                    col = c(G_LIGHT, G_MID, ACCENT), border = G_DARK, whisklty = 1, medcol = G_DARK, medlwd = 1.6)
            med <- vapply(vals, function(v) if (length(v)) stats::median(v) else NA_real_, 0)
            lines(1:3, med, col = G_DARK, lwd = 1.8)
            points(1:3, med, pch = 18, cex = 1.6, col = G_DARK)
            axis(1, at = 1:3, labels = c("T1\nleast soup", "T2", "T3\nmost soup"), padj = 0.5)
            axis(2, las = 1)
            box(bty = "l")
            mtext("rs of DFG genes with MT-CO3", side = 2, line = 3.4, font = 2)
            mtext("Co-expression by soup-index tertile (Gerrits)", side = 3, line = 2.2, font = 2,
                  cex = 1.05)
            mtext(sprintf("T3 minus T1, median over samples: %+.3f [%+.3f, %+.3f], FDR %s %s",
                          g1$estimate, g1$lo, g1$hi, fdr_txt(g1$fdr), stars(g1$fdr)),
                  side = 3, line = 1.0, cex = 0.72, col = "grey40")
            mtext(sprintf("%d sample-by-gene values per tertile; tertiles cut within sample",
                          length(vals[[1]])), side = 3, line = 0.1, cex = 0.72, col = "grey40")
        })

        # == 15H =============================================================
        hr <- H[H$quantity == "rs_median_over_partners" & H$level %in% c("depth", "depth+mito"), ]
        hk <- H[H$quantity == "rs_housekeeping" & H$level %in% c("depth", "depth+mito"), ]
        hr$target <- as.numeric(sub("^target ", "", hr$measure))
        hk$target <- as.numeric(sub("^target ", "", hk$measure))
        obs_rs <- H$estimate[H$quantity == "rs_median_over_partners" & H$level == "observed"]
        obs_hk <- H$estimate[H$quantity == "rs_housekeeping" & H$level == "observed"]
        obs_x <- H$estimate[H$quantity == "median_depth" & H$level == "observed"]
        g_dep <- H$estimate[H$quantity == "reference_depth"]
        band <- H[H$quantity == "reference_rs_raw_iqr", ]
        # A target within 2% of the Gerrits depth stays in Stats_15H but is not drawn.
        near <- setdiff(unique(hr$target[abs(hr$target - g_dep) / g_dep < 0.02]), g_dep)
        if (length(near)) {
            hr <- hr[!hr$target %in% near, ]
            hk <- hk[!hk$target %in% near, ]
            p12_note("15H targets left out of the drawing",
                     sprintf("%s, within 2%% of the Gerrits depth %g; kept in Stats_15H",
                             paste(near, collapse = ", "), g_dep))
        }
        # The MT match does not change the housekeeping pair; drawn once when identical.
        hk_d <- hk[hk$level == "depth", ]; hk_m <- hk[hk$level == "depth+mito", ]
        hk_diff <- abs(hk_d$estimate[order(hk_d$target)] - hk_m$estimate[order(hk_m$target)])
        hk_same <- nrow(hk_d) == nrow(hk_m) && nrow(hk_d) > 0L && all(hk_diff < 1e-12)
        if (hk_same) p12_note("15H housekeeping series identical with and without the MT match",
                              "drawn as one line")
        fig_dev("Fig15H_depth_curve", 6.24, 6.0, function() {
            # Title and subtitle centred; the subtitle wraps to the figure width.
            SUB_CEX <- 0.72
            subs <- wrap_in(c(
                "Median over 20 thinning repeats; bars span the repeats.",
                paste("Solid: depth only; dashed: depth and Gerrits' MT fraction.",
                      "Cells below a target are not thinned."),
                "Grey: housekeeping control (ACTB with cHKG8) on the same cells and draws."),
                width_in = 0.94 * par("din")[1], cex = SUB_CEX)
            ns <- length(subs)
            op <- par(mar = c(6.0, 5.0, 1.6 + 0.95 * ns + 0.3, 1.2), family = "sans")
            on.exit(par(op))
            xs <- c(hr$target, obs_x)
            allv <- c(hr$lo, hr$hi, hk$estimate, obs_rs, obs_hk, band$lo, band$hi)
            yl <- c(min(0, allv, na.rm = TRUE), max(allv, na.rm = TRUE) * 1.14)
            plot(NA, xlim = c(max(xs) * 1.25, min(xs) * 0.8), ylim = yl, log = "x", axes = FALSE,
                 xlab = "", ylab = "")
            rect(10^par("usr")[1], band$lo, 10^par("usr")[2], band$hi, col = BAND, border = NA)
            text(10^par("usr")[1], band$hi, adj = c(-0.02, -0.4), cex = 0.7, col = ACCENT,
                 labels = "Gerrits per-sample raw rs, interquartile range")
            # Dashed series use LTY_DASH (dash 6, gap 2).
            LTY_DASH <- "62"
            sty <- list(list(d = hr, lev = "depth", col = ACCENT, lty = 1),
                        list(d = hr, lev = "depth+mito", col = ACCENT, lty = LTY_DASH),
                        list(d = hk, lev = "depth", col = G_MID, lty = 1))
            if (!hk_same) sty[[4]] <- list(d = hk, lev = "depth+mito", col = G_MID, lty = LTY_DASH)
            for (s in sty) {
                k <- s$d[s$d$level == s$lev, ]
                k <- k[order(k$target), ]
                # Every series gets its range over the repeats, as the subtitle states.
                segments(k$target, k$lo, k$target, k$hi, col = s$col, lwd = 1)
                lines(k$target, k$estimate, col = s$col, lty = s$lty, lwd = 2.2)
                points(k$target, k$estimate, pch = 19, cex = 0.7, col = s$col)
            }
            one <- function(v) length(v) == 1L && is.finite(v)
            if (one(obs_x)) {
                if (one(obs_rs)) points(obs_x, obs_rs, pch = 1, cex = 1.8, lwd = 2.2, col = ACCENT)
                if (one(obs_hk)) points(obs_x, obs_hk, pch = 1, cex = 1.8, lwd = 2.2, col = G_MID)
                text(obs_x, yl[2] * 0.97, "observed", cex = 0.7, col = "grey40")
            }
            tk <- sort(unique(hr$target))
            near <- abs(tk - g_dep) / g_dep < 0.02 & tk != g_dep
            tk <- tk[!near]
            axis(1, at = tk, labels = ifelse(tk == g_dep, sprintf("%d\n(Gerrits)", as.integer(tk)),
                                             as.character(as.integer(tk))),
                 las = 2, cex.axis = 0.75, padj = 0.3)
            axis(2, las = 1)
            box(bty = "l")
            mtext("target depth, UMIs per nucleus (depth falls to the right)", side = 1, line = 4.2,
                  font = 2)
            mtext("Median rs of DFG genes with MT-CO3", side = 2, line = 3.4, font = 2)
            x_mid <- grconvertX(0.5, from = "ndc", to = "user")
            mtext("Depth and mitochondrial content against the Gerrits range (PVM_3 arm)", side = 3,
                  line = 0.3 + 0.95 * ns + 0.3, at = x_mid, font = 2, cex = 0.95)
            for (i in seq_len(ns))
                mtext(subs[i], side = 3, line = 0.3 + 0.95 * (ns - i), at = x_mid, cex = SUB_CEX,
                      col = "grey40")
            if (hk_same)
                legend("topright", bty = "n", cex = 0.72, col = c(ACCENT, ACCENT, G_MID),
                       lty = c("solid", LTY_DASH, "solid"), lwd = 2.2,
                       legend = c("Module genes, depth only", "Module genes, depth + MT",
                                  "housekeeping, depth only = depth + MT"))
            else
                legend("topright", bty = "n", cex = 0.72, col = c(ACCENT, ACCENT, G_MID, G_MID),
                       lty = c("solid", LTY_DASH, "solid", LTY_DASH), lwd = 2.2,
                       legend = c("Module genes, depth only", "Module genes, depth + MT",
                                  "housekeeping, depth only", "housekeeping, depth + MT"))
        })

        # == S10 =============================================================
        MET <- c(own_er = "Own-ER share", ribosome = "Ribosomal share", neat1 = "NEAT1 share",
                 other_lineage = "Other-lineage index")
        sx <- S10[S10$quantity == "post_minus_pre", ]
        fig_dev("FigS10_own_vs_other_lineage", 12.0, 5.2, function() {
            op <- par(mfrow = c(1, length(MET)), mar = c(4.0, 4.6, 3.6, 0.8), oma = c(0, 0, 3.4, 0),
                      family = "sans"); on.exit(par(op))
            for (k in names(MET)) {
                d <- sx[sx$measure == k, ]
                d <- d[match(ARMS, d$arm), ]
                v <- d$estimate
                yr <- range(c(0, v, d$lo, d$hi), na.rm = TRUE)
                yr <- yr + c(-1, 1) * diff(yr) * 0.15
                bp <- barplot(ifelse(is.finite(v), v, 0), names.arg = ARMS, col = NA, border = NA,
                              ylim = yr, las = 1, cex.names = 0.85, axes = FALSE,
                              ylab = if (k == names(MET)[1]) "post-IP minus pre-IP, percentage points" else "",
                              font.lab = 2)
                grid(nx = NA, ny = NULL, col = "grey92", lty = 1)
                barplot(ifelse(is.finite(v), v, 0), col = ARM_COL[ARMS], border = "black",
                        add = TRUE, axes = TRUE, names.arg = rep("", length(ARMS)), las = 1)
                abline(h = 0, lwd = 1.2)
                ok <- is.finite(v)
                ci <- ok & is.finite(d$lo) & is.finite(d$hi) & d$hi > d$lo
                if (any(ci)) arrows(bp[ci], d$lo[ci], bp[ci], d$hi[ci], code = 3, angle = 90,
                                    length = 0.03, lwd = 1.2)
                if (any(!ok)) text(bp[!ok], 0, "untestable", srt = 90, adj = c(-0.1, 0.5),
                                   cex = 0.75, col = G_MID)
                lab <- ifelse(ok, stars(d$fdr), "")
                ytxt <- ifelse(is.finite(d$hi) & v >= 0, d$hi, d$lo)
                if (any(ok)) text(bp[ok], ytxt[ok], lab[ok], pos = ifelse(v[ok] >= 0, 3, 1),
                             cex = 1.1, font = 2)
                mtext(MET[[k]], side = 3, line = 1.2, font = 2, cex = 0.9)
                mtext(sprintf("donors: %s", paste(d$n_donors, collapse = " / ")), side = 3,
                      line = 0.2, cex = 0.65, col = "grey40")
            }
            mtext("Own cytoplasm and other lineages, Module post-IP minus Module pre-IP", side = 3,
                  outer = TRUE, line = 1.8, font = 2, cex = 1.1)
            mtext(paste("Median over donors of the donor difference; bars are 95% donor-bootstrap",
                        "intervals; * FDR < 0.05, ** FDR < 0.01, BH within arm over the four metrics;",
                        "NeoC and NonNeoC pooled within donor."),
                  side = 3, outer = TRUE, line = 0.4, cex = 0.72, col = "grey40")
        })
        p12_say("")
    }

    p12_run(main)
}

# =============================================================================
# RUN
# =============================================================================
local({
    v <- p12_section_12()
    if (!identical(v, "PASS"))
        stop(sprintf("section 12 ended with verdict %s. See Part12_Fig15_log.txt.", v), call. = FALSE)
    message("Figure 15: section 12 PASS. Figures in Part12_rerun/07_Fig15/Fig15")
})
