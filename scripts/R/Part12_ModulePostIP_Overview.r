# =============================================================================
# Part12_ModulePostIP_Overview.r
# SEA-AD CEP-IP/CMRF, Part 12, Module post-IP% overview (Supplementary Figure S13)
#   Panel A: cell groups as a share of the supertype, pooled over donors, per
#   arm and region set. Panel B: Module post-IP% per donor (histogram and
#   density) with skewness, and G1 without the severely affected (SA) donors
#   against the Battery 5 random removals. Descriptive; no test.
#
# Usage (from the SEA-AD folder, or set SEAAD_ROOT), after steps 04 and 05 and
# section 26 of Part12_Fig19_Analysis.r
#   source("Part12_rerun/Part12_ModulePostIP_Overview.r")
#
# Outputs (Part12_rerun/SF_ModulePostIP/): figures (PDF, PNG), tables (CSV), log.
# =============================================================================

source(file.path(Sys.getenv("SEAAD_ROOT", getwd()), "Part12_rerun", "Part12_03_Config.r"))

SCRIPT <- "Part12_ModulePostIP_Overview.r"
STEP <- "SF"   # own step label; no random numbers are drawn

p12_section_modulepostip <- function() {
    D <- p12_step_dirs(STEP, "ModulePostIP")
    p12_clean(D$tables)
    p12_start(STEP, SCRIPT, D)

    PAL <- P12_CFG$palette
    ARM_COL <- unlist(PAL$arms)
    ## Panel B lines for HKG8 use the neutral mid grey (the arm grey is too light).
    LINE_COL <- ARM_COL; LINE_COL[["HKG8"]] <- PAL$neutral$mid
    ARMS <- P12_ARMS; SETS <- P12_DATASETS
    GROUPS <- c("Module post-IP", "Module pre-IP", "non-Module post-IP", "non-Module pre-IP",
                "Outside the arm's cells")

    ## pdf() draws "-" as a long minus; hy() swaps it for the soft hyphen on
    ## the PDF device so labels show a true hyphen.
    ON_PDF <- FALSE
    hy <- function(x) if (ON_PDF) gsub("-", "\u00ad", x, fixed = TRUE) else x

    fig_dev <- function(stem, w, h, draw) {
        ok <- TRUE
        for (kind in c("pdf", "png")) {
            f <- file.path(D$figures, paste0(stem, ".", kind))
            p12_output(f)
            ON_PDF <<- kind == "pdf"
            if (kind == "pdf") grDevices::pdf(f, width = w, height = h, useDingbats = FALSE)
            else grDevices::png(f, width = w, height = h, units = "in", res = 300,
                                type = if (isTRUE(capabilities("cairo"))) "cairo" else "windows")   # cairo keeps semi-transparent fills
            r <- tryCatch({ draw(); TRUE }, error = function(e) {
                p12_say(sprintf("    %s %s: %s", stem, kind, conditionMessage(e))); FALSE })
            grDevices::dev.off()
            ok <- ok && r
        }
        p12_chk(ok, sprintf("%s drawn", stem), "pdf and png")
    }

    main <- function() {
        # == A. Inputs =======================================================
        p12_say("SECTION A -- inputs")
        CU <- utils::read.csv(p12_upstream_file("04_Labels", "04_module_cuts.csv"),
                              check.names = FALSE, stringsAsFactors = FALSE)
        DN <- utils::read.csv(p12_upstream_file("05_Units", "05_donors.csv"),
                              check.names = FALSE, stringsAsFactors = FALSE,
                              colClasses = c(donor_key = "character", donor_id = "character"))
        need_cu <- c("arm", "dataset", "labels_version", "reproduced_post", "n_v1", "n_module",
                     "n_nonmodule_post", "n_nonmodule_pre")
        need_dn <- c("arm", "dataset", "donor_key", "n_cells", "n_post", "n_super",
                     "module_post_pct", "labels_version", "ADNC_rank")
        miss <- c(setdiff(need_cu, names(CU)), setdiff(need_dn, names(DN)))
        if (length(miss)) p12_abort(paste("input lacks:", paste(miss, collapse = ", ")))
        p12_chk(all(CU$labels_version == P12_CFG$labels$labels_version) &&
                    all(DN$labels_version == P12_CFG$labels$labels_version),
                "cut table and donor table carry this labels version",
                P12_CFG$labels$labels_version)
        p12_chk(identical(P12_CFG$quantities$module_post_pct_denominator, "supertype"),
                "Module post-IP% denominator is the supertype",
                P12_CFG$quantities$module_post_pct_denominator)
        p12_say("")

        # == B. Cell groups per arm and region set ===========================
        p12_say("SECTION B -- cell groups")
        GR <- list()
        for (a in ARMS) for (ds in SETS) {
            cu <- CU[CU$arm == a & CU$dataset == ds, ]
            dn <- DN[DN$arm == a & DN$dataset == ds, ]
            if (nrow(cu) != 1L || !nrow(dn)) p12_abort(sprintf("no rows for %s %s", a, ds))
            post <- cu$reproduced_post; pre <- cu$n_module - cu$reproduced_post
            nmp <- cu$n_nonmodule_post; nmq <- cu$n_nonmodule_pre
            sup <- sum(dn$n_super); out <- sup - cu$n_v1
            k <- sprintf("%s %s", a, ds)
            p12_chk(cu$n_module + nmp + nmq == cu$n_v1,
                    sprintf("%s labelled cells = Module + non-Module post + non-Module pre", k),
                    sprintf("%d = %d + %d + %d", cu$n_v1, cu$n_module, nmp, nmq))
            p12_chk(sum(dn$n_post) == post,
                    sprintf("%s Module post-IP cells in the cut table = sum over donors", k),
                    sprintf("%d against %d", post, sum(dn$n_post)))
            p12_chk(sum(dn$n_cells) == cu$n_v1,
                    sprintf("%s labelled cells in the cut table = sum over donors", k),
                    sprintf("%d against %d", cu$n_v1, sum(dn$n_cells)))
            p12_chk(out >= 0, sprintf("%s supertype cells cover the labelled cells", k),
                    sprintf("%d supertype, %d labelled", sup, cu$n_v1))
            n <- c(post, pre, nmp, nmq, out)
            GR[[k]] <- data.frame(arm = a, region_set = ds, group = GROUPS, n_cells = n,
                                  n_supertype = sup, pct_supertype = 100 * n / sup,
                                  n_donors = nrow(dn))
        }
        GR <- do.call(rbind, GR); rownames(GR) <- NULL
        f_gr <- file.path(D$tables, "SF_cell_groups.csv")
        p12_chk(p12_ascii_ok(GR), "cell-group table is ASCII", "checked")
        p12_output(f_gr); utils::write.csv(GR, f_gr, row.names = FALSE)
        p12_say("")

        # == C. Module post-IP% per donor ====================================
        p12_say("SECTION C -- Module post-IP% per donor")
        dv <- 100 * DN$n_post / DN$n_super
        p12_chk(all(abs(dv - DN$module_post_pct) < 1e-9, na.rm = TRUE) && !anyNA(DN$module_post_pct),
                "donor Module post-IP% = 100 x n_post / n_super",
                sprintf("%d donor rows | max |difference| %.2g", nrow(DN),
                        max(abs(dv - DN$module_post_pct), na.rm = TRUE)))
        SM <- list()
        for (a in ARMS) for (ds in SETS) {
            v <- sort(DN$module_post_pct[DN$arm == a & DN$dataset == ds])
            q <- stats::quantile(v, c(0.25, 0.5, 0.75), names = FALSE, type = 7)
            ## G1: bias-adjusted sample skewness. Quartile (Bowley) skewness:
            ## (Q3 + Q1 - 2 x median) / (Q3 - Q1). Descriptive only.
            n <- length(v); dv0 <- v - mean(v)
            g1 <- mean(dv0^3) / mean(dv0^2)^1.5
            G1 <- if (n > 2) sqrt(n * (n - 1)) / (n - 2) * g1 else NA_real_
            bow <- if (q[3] > q[1]) (q[3] + q[1] - 2 * q[2]) / (q[3] - q[1]) else NA_real_
            SM[[paste(a, ds)]] <- data.frame(arm = a, region_set = ds, n_donors = length(v),
                                             median = q[2], mean = mean(v), q1 = q[1], q3 = q[3],
                                             min = min(v), max = max(v), donors_at_0 = sum(v == 0),
                                             skewness_g1 = G1, skewness_quartile = bow)
        }
        SM <- do.call(rbind, SM); rownames(SM) <- NULL
        for (i in seq_len(nrow(SM)))
            p12_note(sprintf("%s %s", SM$arm[i], SM$region_set[i]),
                     sprintf(paste("n %d | median %.2f | mean %.2f | IQR %.2f-%.2f | range %.2f-%.2f |",
                                   "at 0%%: %d | skewness G1 %.2f | quartile skewness %.2f"),
                             SM$n_donors[i], SM$median[i], SM$mean[i], SM$q1[i], SM$q3[i],
                             SM$min[i], SM$max[i], SM$donors_at_0[i], SM$skewness_g1[i],
                             SM$skewness_quartile[i]))
        f_sm <- file.path(D$tables, "SF_donor_summary.csv")
        p12_output(f_sm); utils::write.csv(SM, f_sm, row.names = FALSE)
        p12_say("")

        # == C2. SA donors and the tail ======================================
        ## Descriptive. Reuses the Battery 5 SA donors, frames and random-removal
        ## draws (section 26 of Part12_Fig19_Analysis.r). G1 is recomputed without
        ## the frame's k SA donors and after each draw that removes k other
        ## High-ADNC donors; p_low = (1 + draws with change <= SA change) / (1 + draws).
        p12_say("SECTION C2 -- SA donors and the tail (descriptive)")
        B5_DIR <- "11_Fig19/19_SADonors"
        SAt <- utils::read.csv(p12_upstream_file(B5_DIR, "SA_donors_19.csv"), check.names = FALSE,
                               colClasses = c(donor_key = "character", donor_id = "character"))
        FRt <- utils::read.csv(p12_upstream_file(B5_DIR, "Frames_19.csv"), check.names = FALSE,
                               stringsAsFactors = FALSE)
        DRt <- utils::read.csv(p12_upstream_file(B5_DIR, "NullDraws_19.csv"), check.names = FALSE,
                               colClasses = c(frame = "character", removed = "character"))
        SA <- SAt$donor_key[as.logical(SAt$severely_affected)]
        p12_chk(length(SA) == 11L, "Battery 5 lists eleven SA donors", sprintf("%d listed", length(SA)))
        p12_note("SA donors", paste(sort(SAt$donor_id[SAt$donor_key %in% SA]), collapse = ", "))
        skew_g1 <- function(v) {
            n <- length(v); d0 <- v - mean(v)
            if (n < 3) return(NA_real_)
            sqrt(n * (n - 1)) / (n - 2) * mean(d0^3) / mean(d0^2)^1.5
        }
        TL <- list()
        for (a in ARMS) for (ds in SETS) {
            fr <- paste("DN", a, ds, sep = "|")
            dn <- DN[DN$arm == a & DN$dataset == ds, ]
            v <- dn$module_post_pct; names(v) <- dn$donor_key
            is_sa <- dn$donor_key %in% SA; k <- sum(is_sa)
            fi <- FRt[FRt$frame == fr, ]
            p12_chk(nrow(fi) == 1L && fi$n_donors == nrow(dn) && fi$k == k,
                    sprintf("%s: Battery 5 frame has the same donors and SA count", fr),
                    sprintf("donors %d against %d | k %d against %d", if (nrow(fi)) fi$n_donors else NA,
                            nrow(dn), if (nrow(fi)) fi$k else NA, k))
            g_all <- skew_g1(v); g_wo <- skew_g1(v[!is_sa])
            rm <- DRt$removed[DRt$frame == fr]
            g_null <- vapply(strsplit(rm, ";", fixed = TRUE), function(x) skew_g1(v[!(names(v) %in% x)]), 0)
            ok_null <- length(g_null) > 0L
            ch_null <- g_null - g_all
            top <- v >= stats::quantile(v, 0.9, type = 7)
            TL[[fr]] <- data.frame(
                arm = a, region_set = ds, n_donors = length(v), k_sa = k, n_draws = length(g_null),
                g1_all = g_all, g1_without_sa = g_wo, change = g_wo - g_all,
                null_g1_p05 = if (ok_null) stats::quantile(g_null, 0.05, names = FALSE) else NA_real_,
                null_g1_p95 = if (ok_null) stats::quantile(g_null, 0.95, names = FALSE) else NA_real_,
                p_low = if (ok_null) (1 + sum(ch_null <= g_wo - g_all)) / (1 + sum(is.finite(ch_null))) else NA_real_,
                sa_in_top10 = sum(top & is_sa), n_top10 = sum(top))
            r <- TL[[fr]]
            p12_note(sprintf("%s %s", a, ds),
                     sprintf("k %d | G1 %.2f all, %.2f without SA | random removal 5-95%% %.2f-%.2f (%d draws) | p_low %.3f | SA in top 10%%: %d of %d",
                             k, g_all, g_wo, r$null_g1_p05, r$null_g1_p95, r$n_draws, r$p_low,
                             r$sa_in_top10, r$n_top10))
        }
        TL <- do.call(rbind, TL); rownames(TL) <- NULL
        f_tl <- file.path(D$tables, "SF_sa_tail.csv")
        p12_output(f_tl); utils::write.csv(TL, f_tl, row.names = FALSE)
        p12_say("")

        # == D. Figures ======================================================
        p12_say("SECTION D -- figures")
        ## Panel A: heatmap on one single-hue scale (cream to dark slate), chosen to
        ## avoid hues used elsewhere in the paper. Scale 0-40%; values above 40% take
        ## the top colour (key marks 40+). Labels are white above 25%.
        ramp <- function(cols, lo = 0, hi = 1) {
            full <- grDevices::colorRampPalette(cols)(1001)
            full[round(seq(1 + 1000 * lo, 1 + 1000 * hi, length.out = 101))]
        }
        HEAT <- list(khaki = list(col = ramp(c("#FBF8EE", "#E4D9A9", "#9A9A6A", "#2F4A4A"))))
        LAYOUT <- list(list(stem = "SF_A_cell_groups", w = 7.4, h = 4.6, asp = NA))
        VMAX <- 40
        TXT_WHITE_ABOVE <- 25
        heat_col <- HEAT$khaki$col
        col_of <- function(v) heat_col[1 + round(100 * min(v, VMAX) / VMAX)]
        for (lay in LAYOUT) {
        fig_dev(lay$stem, lay$w, lay$h, function() {
            rows <- paste(rep(ARMS, each = length(SETS)), rep(SETS, times = length(ARMS)))
            nr <- length(rows); nc <- length(GROUPS); GAP <- 0.18
            y0 <- (seq_len(nr) - 1) + GAP * ((seq_len(nr) - 1) %/% length(SETS))
            ytot <- max(y0) + 1
            graphics::par(mar = c(0.5, 9, 4.4, 9), family = "sans", xpd = NA)
            graphics::plot.new()
            graphics::plot.window(xlim = c(0, nc), ylim = c(ytot, 0), xaxs = "i", yaxs = "i", asp = lay$asp)
            for (i in seq_len(nr)) for (j in seq_len(nc)) {
                r <- GR[paste(GR$arm, GR$region_set) == rows[i] & GR$group == GROUPS[j], ]
                v <- r$pct_supertype; fill <- col_of(v)
                graphics::rect(j - 1, y0[i], j, y0[i] + 1, col = fill, border = "black", lwd = 1.8)
                graphics::text(j - 0.5, y0[i] + 0.5, sprintf("%.1f", v), cex = 0.95,
                               col = if (v > TXT_WHITE_ABOVE) "white" else PAL$neutral$dark,
                               font = if (j == 1) 2 else 1)
            }
            ## Column labels in bold; the Module post-IP column slightly larger.
            graphics::text(seq_len(nc) - 0.5, -0.15, hy(c("Module\npost-IP", "Module\npre-IP",
                           "non-Module\npost-IP", "non-Module\npre-IP", "Outside the\narm's cells")),
                           adj = c(0.5, 0), cex = c(0.85, 0.72, 0.72, 0.72, 0.72), font = 2)
            ## Row labels are the region sets; each arm's pair of rows is joined by
            ## a square bracket with the arm name to its left.
            graphics::text(-0.1, y0 + 0.5, rep(SETS, times = length(ARMS)), adj = c(1, 0.5), cex = 0.8)
            bx <- -0.1 - max(graphics::strwidth(SETS, cex = 0.8)) - 0.1
            tick <- 0.06
            for (k in seq_along(ARMS)) {
                top <- y0[(k - 1) * length(SETS) + 1] + 0.18
                bot <- y0[k * length(SETS)] + 1 - 0.18
                graphics::segments(bx, top, bx, bot, lwd = 1.2)
                graphics::segments(bx, c(top, bot), bx + tick, c(top, bot), lwd = 1.2)
                mid <- (top + bot) / 2
                graphics::segments(bx - tick, mid, bx, mid, lwd = 1.2)   # short tick from the midpoint
                graphics::text(bx - tick - 0.04, mid, hy(ARMS[k]), adj = c(1, 0.5), cex = 0.85, font = 2)
            }
            sup <- GR$n_supertype[GR$group == GROUPS[1]]
            graphics::text(nc + 0.12, y0 + 0.5, format(sup, big.mark = ","), adj = c(0, 0.5), cex = 0.75)
            graphics::text(nc + 0.12, -0.15, "Supertype\ncells", adj = c(0, 0), cex = 0.72)
            ## colour key
            kx <- nc + 1.05; kw <- 0.18
            for (s in 0:99) graphics::rect(kx, ytot - ytot * s / 100, kx + kw, ytot - ytot * (s + 1) / 100,
                                           col = heat_col[s + 1], border = NA)
            graphics::rect(kx, 0, kx + kw, ytot, border = "black", lwd = 0.6)
            tk <- seq(0, VMAX, by = 10)
            graphics::text(kx + kw + 0.06, ytot - ytot * tk / VMAX,
                           c(tk[-length(tk)], paste0(VMAX, "+")), adj = c(0, 0.5), cex = 0.65)
            graphics::text(kx + kw + 0.45, ytot / 2, hy("% of supertype cells"), srt = 90, cex = 0.7)
            graphics::title(main = hy("Cell groups as a share of the supertype, pooled over donors"),
                            adj = 0, line = 3.1, cex.main = 0.95)
        })
        }

        ## Panel B: share of donors per bin (%), NeoC filled and NonNeoC outline bars
        ## in the arm colour, with kernel density curves reflected at 0% and 100%
        ## (bw.nrd0), scaled to the bar axis. About 12-15 bins per arm.
        tint <- function(col, a) grDevices::adjustcolor(col, alpha.f = a)
        shade <- function(col, f) {   # f < 1 darker, f > 1 towards white
            m <- grDevices::col2rgb(col) / 255
            m <- if (f <= 1) m * f else m + (1 - m) * (f - 1)
            grDevices::rgb(m[1], m[2], m[3])
        }
        refl_density <- function(v, xs) {
            bw <- stats::bw.nrd0(v)
            d <- stats::density(c(v, -v, 200 - v), bw = bw, from = 0, to = 100, n = 1024)
            stats::approx(d$x, 3 * d$y, xs)$y
        }
        ## Per region set: bar fill, border, border line type, curve colour, curve line type.
        SCHEMES <- list(
            v5 = list(name = "filled against outline bars, arm colour", f = function(a) list(
                NeoC = list(fill = tint(LINE_COL[[a]], 0.45), border = shade(LINE_COL[[a]], 0.7), lty = 1,
                            curve = shade(LINE_COL[[a]], 0.7), clty = 1),
                NonNeoC = list(fill = NA, border = shade(LINE_COL[[a]], 0.55), lty = 2,
                               curve = shade(LINE_COL[[a]], 0.55), clty = 2))))
        for (vk in names(SCHEMES)) {
            sc <- SCHEMES[[vk]]
            fig_dev("SF_B_donor_histogram", 7.4, 6.6, function() {
                ## Larger text except tick labels; x tick labels sit closer to the ticks.
                graphics::par(mfrow = c(2, 2), mar = c(3.6, 3.9, 2.4, 0.8), mgp = c(2.2, 0.35, 0),
                              family = "sans", las = 1, tcl = -0.2, cex.lab = 1.1, font.lab = 2,
                              cex.axis = 0.8)
                for (a in ARMS) {
                    st <- sc$f(a)
                    va <- DN$module_post_pct[DN$arm == a]
                    w <- c(0.5, 1, 2, 2.5, 5, 10)[which.min(abs(max(va) / c(0.5, 1, 2, 2.5, 5, 10) - 13))]
                    br <- seq(0, w * max(1, ceiling(max(va) / w - 1e-9)), by = w)
                    xs <- seq(0, max(br), length.out = 400)
                    H <- lapply(SETS, function(ds) {
                        v <- DN$module_post_pct[DN$arm == a & DN$dataset == ds]
                        list(v = v, pct = 100 * graphics::hist(v, breaks = br, right = FALSE,
                                                               include.lowest = TRUE, plot = FALSE)$counts / length(v),
                             dens = 100 * w * refl_density(v, xs))
                    }); names(H) <- SETS
                    ymax <- max(unlist(lapply(H, function(h) c(h$pct, h$dens))))
                    graphics::plot(NA, xlim = range(br), ylim = c(0, ymax * 1.6), xaxs = "i", yaxs = "i",
                                   xlab = "", ylab = "Donors (%)", bty = "l", xaxt = "n")
                    graphics::axis(1, mgp = c(2.2, 0.15, 0))
                    graphics::title(xlab = hy("Module post-IP% (%)"), line = 2.0)
                    graphics::abline(h = pretty(c(0, ymax)), col = PAL$neutral$rule, lwd = 0.5)
                    for (ds in SETS) {
                        s <- st[[ds]]
                        graphics::rect(br[-length(br)], 0, br[-1], H[[ds]]$pct, col = s$fill,
                                       border = s$border, lty = s$lty, lwd = 0.9)
                    }
                    for (ds in SETS)
                        graphics::lines(xs, H[[ds]]$dens, col = st[[ds]]$curve, lwd = 2, lty = st[[ds]]$clty)
                    lab <- vapply(SETS, function(ds) {
                        s <- SM[SM$arm == a & SM$region_set == ds, ]
                        sprintf("%s (n=%d; median %.2f)", ds, s$n_donors, s$median)
                    }, "")
                    graphics::legend("topright", legend = hy(lab), bty = "n", cex = 1.0, y.intersp = 1.5,
                                     fill = vapply(SETS, function(ds) {
                                         f <- st[[ds]]$fill; if (is.na(f)) "white" else f }, ""),
                                     border = vapply(SETS, function(ds) st[[ds]]$border, ""),
                                     lty = vapply(SETS, function(ds) st[[ds]]$clty, 1),
                                     col = vapply(SETS, function(ds) st[[ds]]$curve, ""),
                                     lwd = 2, seg.len = 2.2, merge = TRUE)
                    ## Shape, descriptive: skewness G1 and quartile skewness per region set.
                    sk <- SM[SM$arm == a, ]
                    sk <- sk[match(SETS, sk$region_set), ]
                    u <- graphics::par("usr")
                    xr <- u[2] - 0.03 * (u[2] - u[1]); y1 <- u[3] + 0.70 * (u[4] - u[3])
                    lh <- 1.5 * graphics::strheight("Ag", cex = 0.92)
                    graphics::text(xr, y1, hy("Skewness (G1; quartile):"), adj = c(1, 1),
                                   cex = 0.92, font = 2, col = "black")
                    graphics::text(xr, y1 - lh * seq_along(SETS),
                                   hy(sprintf("%s (%.2f; %.2f)", SETS, sk$skewness_g1, sk$skewness_quartile)),
                                   adj = c(1, 1), cex = 0.92, col = "black")
                    tl <- TL[TL$arm == a, ]; tl <- tl[match(SETS, tl$region_set), ]
                    y2 <- y1 - lh * (length(SETS) + 1.4)
                    graphics::text(xr, y2, hy("G1 without SA (random 5-95%):"), adj = c(1, 1),
                                   cex = 0.92, font = 2, col = "black")
                    graphics::text(xr, y2 - lh * seq_along(SETS),
                                   hy(sprintf("%s k=%d (%.2f; %.2f-%.2f)", SETS, tl$k_sa, tl$g1_without_sa,
                                              tl$null_g1_p05, tl$null_g1_p95)),
                                   adj = c(1, 1), cex = 0.92, col = "black")
                    graphics::title(main = hy(sprintf("%s (bin %s pp)", a, format(w))), adj = 0,
                                    cex.main = 1.25, line = 0.7)
                }
            })
        }
        p12_note("HKG8 colour in panel B", sprintf("%s (neutral mid grey) in place of %s",
                                                        PAL$neutral$mid, ARM_COL[["HKG8"]]))
    }
    p12_run(main)
}

p12_section_modulepostip()
