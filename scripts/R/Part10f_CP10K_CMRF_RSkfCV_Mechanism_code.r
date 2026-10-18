# =============================================================================
# Part10f_CP10K_CMRF_RSkfCV_Mechanism_code.r
# SEA-AD CEP-IP/CMRF, Part 10f, RSk-fCV Stage 2 gate mechanism (Figure 10)
#   Reads the saved Part 10e CMRF results (runs no CMRF) and draws how the
#   RSk-fCV Stage 2 gate works. A candidate is removed only when both hold:
#   paired Wilcoxon p < alpha, and dRMSE = median(RMSE_B - RMSE_A) /
#   median(RMSE_A) > BV, where BV = IQR(RMSE_A) / median(RMSE_A).
#   Arm Mic13_PVM3 (M4), NeoC.
#
# Usage (from the SEA-AD folder, or set SEAAD_ROOT), after Part 10e
#   source("Part10f_CP10K_CMRF_RSkfCV_Mechanism_code.r")
#
# Outputs (Part10f_RSkfCV_Mechanism/Part10f_RSkfCV_Mechanism_<arm>/): figure
#   pages (multi-page PDF, one PNG per page) and workbook.
# =============================================================================

suppressPackageStartupMessages({
    library(ggplot2); library(grid); library(gridExtra); library(writexl)
})

# == Paths and options ========================================================
if (!exists("CEP_HOME")) CEP_HOME <- Sys.getenv("SEAAD_ROOT", getwd())
SHOW <- TRUE
SEED <- 42L         # must match Part 10e's gate seed
# Per-fold bar panels: "repeat1" (repeat 1's five folds) or "mean" (over 40 repeats)
BAR_ROWS <- "repeat1"

# One entry per figure
TARGETS <- list(
    list(tag = "M4 NeoC",   dir = "Part10e_CMRF_Mic13_PVM3",
         rds = "Part10e_CP10K_M4_CMRF_all_results.rds",   ds = "NeoC",
         sub = "Mic13_PVM3", pre = "Part10e_CP10K_M4_RSkfCV"),
    NULL
)
TARGETS <- Filter(Negate(is.null), TARGETS)

OUT_ROOT <- file.path(CEP_HOME, "Part10f_RSkfCV_Mechanism")
PNG_W <- 3300; PNG_H <- 2600; PNG_RES <- 200   # per-page PNG geometry
PDF_W <- 15;   PDF_H <- 12

# == Template cosmetics =======================================================
theme_copula <- function(base_size = 10) {
    theme_minimal(base_size = base_size) +
        theme(
            plot.title       = element_text(face = "bold", size = base_size + 1),
            plot.subtitle    = element_text(size = base_size - 1, color = "gray30"),
            axis.title       = element_text(face = "bold", size = base_size - 1),
            axis.text        = element_text(size = base_size - 2),
            panel.grid.major = element_line(color = "gray90"),
            panel.grid.minor = element_line(color = "gray95"),
            panel.border     = element_rect(color = "black", fill = NA,
                                            linewidth = 0.35),
            plot.background  = element_rect(fill = "white", color = NA),
            panel.background = element_rect(fill = "white", color = NA)
        )
}
TH <- function() theme_copula(base_size = 9.5)

COL_A     <- "#009E73"    # Curve A / OriRef / alpha and BV lines
COL_B     <- "#1A237E"    # Curve B / navy trajectory / PAVA step
COL_CAND  <- "#7C3AED"    # candidates
COL_TEST  <- "#D94F4F"    # test cells, worst cell, removal region
COL_TRAIN <- rgb(  0, 130, 176, maxColorValue = 255)
COL_TESTB <- rgb(238, 207, 206, maxColorValue = 255)
COL_ZERO  <- "#8B4513"    # zero reference line on every dRMSE plot
DECILE_COLS <- c(
    D1  = rgb(  7,   5,  27, maxColorValue = 255),
    D2  = rgb( 40,  11,  83, maxColorValue = 255),
    D3  = rgb( 82,  14, 109, maxColorValue = 255),
    D4  = rgb(120,  28, 109, maxColorValue = 255),
    D5  = rgb(159,  42,  99, maxColorValue = 255),
    D6  = rgb(195,  59,  79, maxColorValue = 255),
    D7  = rgb(225,  86,  53, maxColorValue = 255),
    D8  = rgb(245, 125,  21, maxColorValue = 255),
    D9  = rgb(252, 168,  13, maxColorValue = 255),
    D10 = rgb(246, 215,  70, maxColorValue = 255))

# ocean.curl for the heatmaps (pals when installed, else its anchor colours)
OCEAN_CURL <- if (requireNamespace("pals", quietly = TRUE)) pals::ocean.curl(256) else
    c("#151D44", "#1D4F6B", "#2F8C7C", "#6FBE96", "#B7DDBA", "#F7F5EF",
      "#F2CDB2", "#E29A74", "#C96A5B", "#A33E55", "#6B1A3B", "#341E43")
TILE_BORDER <- "black"   # thin outline on every heatmap cell
TILE_LW     <- 0.12
COL_CAND_D  <- "#4C1D95"   # darker violet, candidate dot borders
COL_REMOVED <- "#000000"   # removed bars, ranking panels
COL_BAR_ACC <- "#5B21B6"   # accepted bars, saturated violet
COL_STICK   <- "#8A9BAD"   # residual sticks on the test-fold panel
COL_DIFF    <- "#D96A6A"   # Q and R boxplot+jitter fill (salmon)

# Boxplot median line width 0.88 (median.linewidth from ggplot2 4.0.0, else fatten)
BOX_LW      <- 0.55
BOX_MED_LW  <- 0.88
BOX_MED_ARG <- if (utils::packageVersion("ggplot2") >= "4.0.0")
    list(median.linewidth = BOX_MED_LW) else list(fatten = BOX_MED_LW / BOX_LW)
geom_box_std <- function(...) do.call(geom_boxplot, c(list(...), BOX_MED_ARG))

# == Part 10e helpers, copied verbatim so folds reproduce exactly =============
iso_predict <- function(u_train, v_train, u_new) {
    ord <- order(u_train); u_s <- u_train[ord]; n_tr <- length(u_s)
    iso_v <- isoreg(u_s, v_train[ord])$yf
    idx <- findInterval(u_new, u_s, rightmost.closed = TRUE)
    idx[idx < 1L] <- 1L; idx[idx > n_tr] <- n_tr
    iso_v[idx]
}
strat_folds <- function(u_set, K, rep_seed) {
    n_set <- length(u_set)
    n_strat <- max(1L, min(10L, floor(n_set / K)))
    u_rank <- rank(u_set, ties.method = "first")
    strat  <- pmin(ceiling(u_rank / n_set * n_strat), n_strat)
    fa <- integer(n_set)
    for (s in seq_len(n_strat)) {
        idx_s <- which(strat == s); set.seed(rep_seed + s)
        # sample.int keeps single-cell strata correct
        sh <- idx_s[sample.int(length(idx_s))]
        fa[sh] <- ((seq_along(sh) - 1L) %% K) + 1L
    }
    lapply(seq_len(K), function(k) fa == k)
}
run_cv <- function(u, v, K, R, seed = SEED, offset = 0L) {
    unlist(lapply(seq_len(R), function(r) {
        tm <- strat_folds(u, K, seed + offset + r * 1000L)
        vapply(seq_len(K), function(k) {
            te <- tm[[k]]; tr <- !te
            if (sum(tr) < 2L || sum(te) < 1L) return(NA_real_)
            sqrt(mean((v[te] - iso_predict(u[tr], v[tr], u[te]))^2))
        }, numeric(1))
    }), use.names = FALSE)
}
# Step data for a geom_step curve
step_df <- function(u_tr, v_tr) {
    o <- order(u_tr)
    data.frame(u = u_tr[o], y = isoreg(u_tr[o], v_tr[o])$yf)
}

# =============================================================================
# Per-target pages
# =============================================================================
make_pages <- function(RDS_PATH, DS, TAG) {

# == Load =====================================================================
stopifnot(file.exists(RDS_PATH))
res <- readRDS(RDS_PATH)[[DS]]
df  <- res$df_refined
s2  <- res$stage2_res
stopifnot(!is.null(df), !is.null(s2))

K <- s2$k_folds; R <- s2$n_repeats; ALPHA <- s2$alpha; BV <- s2$noise_floor
# Row indices come from stage2_res, which is in the gate's cell order
if (is.null(s2$anchor_row) || is.null(s2$cand_row))
    stop("stage2_res has no anchor_row/cand_row; rerun Part 10e.", call. = FALSE)
a_idx <- s2$anchor_row
c_idx <- s2$cand_row
stopifnot(length(c_idx) == s2$n_candidates)
if (!identical(sort(c_idx), sort(which(!df$is_orig_skcs & df$pass_stage1))))
    stop("stage2_res$cand_row does not match df's Stage 1 candidates.",
         call. = FALSE)
if (!identical(sort(a_idx),
               sort(which(df$cell_class == "Original + Refined"))))
    stop("stage2_res$anchor_row does not match df's Orig+Ref cells.",
         call. = FALSE)

u_a <- df$u_all[a_idx]; v_a <- df$v_all[a_idx]
u_c <- df$u_all[c_idx]; v_c <- df$v_all[c_idx]
n_a <- length(u_a); n_c <- length(u_c)

rem_ord  <- s2$removal_order; kept <- s2$pass_stage2
resid_pv <- s2$resid_pava_cand
p_iter   <- s2$p_per_iter; d_iter <- s2$effect_size_per_iter
n_it     <- length(p_iter)

# The gate closes when either condition fails; GATE_BY names which one did
GATE_P_OK <- s2$wilcox_p_final >= ALPHA          # p no longer significant
GATE_D_OK <- s2$relative_delta_final <= BV       # effect within background
# Braces needed for the top-level if/else chain
GATE_BY <- {
    if (GATE_P_OK && GATE_D_OK) "both conditions"
    else if (GATE_P_OK)         "the Wilcoxon p"
    else if (GATE_D_OK)         "the dRMSE effect size"
    else                        "neither (iteration cap reached)"
}

cat(sprintf(paste0("Figure 10 inputs -- %s\n",
                   "  OriRef %d | candidates %d | removed %d | kept %d\n",
                   "  iterations %d | final p=%.4f | dRMSE=%.4f | BV=%.4f\n",
                   "  CV: %d folds x %d repeats = %d fold-RMSEs, alpha=%.2f\n"),
            TAG, n_a, n_c, s2$n_removed, sum(kept), n_it,
            s2$wilcox_p_final, s2$relative_delta_final, BV, K, R, K * R, ALPHA))

# == Recompute the two iterations the figure shows ============================
rmse_A <- s2$rmse_grp1a
u_p1 <- c(u_a, u_c);       v_p1 <- c(v_a, v_c)
u_pN <- c(u_a, u_c[kept]); v_pN <- c(v_a, v_c[kept])
# Curve B at closure: the stored vector the gate tested
rmse_BN <- s2$rmse_final_union
# Curve B at iteration 1 is rebuilt; the check below tests the rebuild at closure
cat("  reconstructing Curve B, iteration 1 ...\n")
rmse_B1 <- run_cv(u_p1, v_p1, K, R)
cat("  verifying the reconstruction against the stored closure union ...\n")
.chk <- run_cv(u_pN, v_pN, K, R)
.dev <- max(abs(.chk - rmse_BN), na.rm = TRUE)
RECON_OK <- is.finite(.dev) && .dev < 1e-10
cat(sprintf("    max |recomputed - stored| at closure = %.3e -> %s\n", .dev,
            if (RECON_OK) "reproduces exactly" else
            "does not reproduce"))

wp <- function(b) suppressWarnings(wilcox.test(b, rmse_A, paired = TRUE,
                                   alternative = "greater", exact = FALSE))$p.value
p1 <- wp(rmse_B1); pN <- wp(rmse_BN)
# relative_delta as in Part 10e: median(B - A) / median(A), NAs dropped first
a_valid <- rmse_A[!is.na(rmse_A)]
med_A   <- median(a_valid)
IQR_A   <- as.numeric(IQR(a_valid))          # BV * med_A, the absolute threshold
rel_delta <- function(rB) {
    b <- rB[!is.na(rB)]
    median(b - a_valid[seq_along(b)]) / med_A
}
d1 <- rel_delta(rmse_B1); dN <- rel_delta(rmse_BN)
.d_chk <- rel_delta(.chk)
cat(sprintf(paste0("  dRMSE check: recomputed closure %.4f vs Part 10e %.4f",
                   " (diff %.2e)\n"),
            .d_chk, s2$relative_delta_final,
            .d_chk - s2$relative_delta_final))
cat(sprintf("  iter 1  : p=%.3g  dRMSE=%.4f\n  iter %d : p=%.4f  dRMSE=%.4f\n",
            p1, d1, n_it, pN, dN))

# == Paired tests for the figure ==============================================
# Paired on (repeat, fold); Shapiro-Wilk picks paired t or Wilcoxon; BH over 3
paired_test <- function(x, y) {
    d <- x - y
    d <- d[is.finite(d)]
    sw <- suppressWarnings(shapiro.test(d))$p.value
    if (is.finite(sw) && sw > 0.05) {
        tt <- t.test(x, y, paired = TRUE)
        list(name = "paired t-test", p = tt$p.value, shapiro = sw)
    } else {
        wt <- suppressWarnings(wilcox.test(x, y, paired = TRUE, exact = FALSE))
        list(name = "Wilcoxon signed-rank", p = wt$p.value, shapiro = sw)
    }
}
d_first <- (rmse_B1 - rmse_A) / med_A
d_last  <- (rmse_BN - rmse_A) / med_A
t_F  <- paired_test(rmse_B1, rmse_A)     # RMSE_B vs RMSE_A, iteration 1
t_M  <- paired_test(rmse_BN, rmse_A)     # RMSE_B vs RMSE_A, final iteration
t_R3 <- paired_test(d_first, d_last)     # dRMSE, first vs final iteration
FDR  <- p.adjust(c(t_F$p, t_M$p, t_R3$p), method = "BH")
names(FDR) <- c("F", "M", "R3")
fdr_txt <- function(k) {
    if (FDR[[k]] < 1e-3) sprintf("%.2e", FDR[[k]]) else sprintf("%.4f", FDR[[k]])
}
cat(sprintf(paste0("  Paired tests (BH over 3):\n",
                   "    F  %-22s p=%.3g  FDR=%.3g  (Shapiro p=%.3g)\n",
                   "    M  %-22s p=%.3g  FDR=%.3g  (Shapiro p=%.3g)\n",
                   "    R3 %-22s p=%.3g  FDR=%.3g  (Shapiro p=%.3g)\n"),
            t_F$name,  t_F$p,  FDR[["F"]],  t_F$shapiro,
            t_M$name,  t_M$p,  FDR[["M"]],  t_M$shapiro,
            t_R3$name, t_R3$p, FDR[["R3"]], t_R3$shapiro))

fold_id <- factor(rep(seq_len(K), times = R))
rep_id  <- factor(rep(seq_len(R), each = K))

# =============================================================================
# A: K-fold design    B: stratification by decile
# =============================================================================
n_pool <- n_a + n_c; n_te <- round(n_pool / K); n_tr <- n_pool - n_te
gap <- 0.06; rect_l <- list(); lab_l <- list()
for (r in 1:K) {
    rect_l[[length(rect_l) + 1L]] <- data.frame(
        xmin = r - 0.47, xmax = r + 0.47, ymin = r - 0.42, ymax = r + 0.42,
        Role = "TEST")
    lab_l[[length(lab_l) + 1L]] <- data.frame(x = r, y = r, label = n_te,
                                              Role = "TEST")
    if (r > 1) rect_l[[length(rect_l) + 1L]] <- data.frame(
        xmin = 0.53, xmax = r - 0.47 - gap, ymin = r - 0.42, ymax = r + 0.42,
        Role = "TRAIN")
    if (r < K) rect_l[[length(rect_l) + 1L]] <- data.frame(
        xmin = r + 0.47 + gap, xmax = K + 0.47, ymin = r - 0.42, ymax = r + 0.42,
        Role = "TRAIN")
    lx <- if (r == 1) (r + 0.47 + gap + K + 0.47) / 2
          else if (r == K) (0.53 + r - 0.47 - gap) / 2
          else if ((r - 1) >= (K - r)) (0.53 + r - 0.47 - gap) / 2
          else (r + 0.47 + gap + K + 0.47) / 2
    lab_l[[length(lab_l) + 1L]] <- data.frame(x = lx, y = r, label = n_tr,
                                              Role = "TRAIN")
}
d_rect <- do.call(rbind, rect_l); d_lab <- do.call(rbind, lab_l)
pA <- ggplot() +
    geom_rect(data = d_rect, aes(xmin = xmin, xmax = xmax, ymin = ymin,
                                 ymax = ymax, fill = Role),
              colour = "white", linewidth = 0.8) +
    geom_text(data = d_lab, aes(x = x, y = y, label = label), size = 3.5,
              fontface = "bold",
              colour = ifelse(d_lab$Role == "TRAIN", "white", "grey30")) +
    scale_fill_manual(values = c(TEST = COL_TESTB, TRAIN = COL_TRAIN),
                      name = NULL) +
    scale_x_continuous(breaks = NULL, expand = expansion(mult = 0.04)) +
    scale_y_continuous(breaks = 1:K, labels = paste("Fold", 1:K),
                       expand = expansion(mult = 0.08)) +
    labs(title = "A: K-fold CV (K=5) -- pool split into folds",
         subtitle = sprintf(paste0("Pool = OriRef (%d) + candidates (%d) = %d cells.\n",
                                   "Test (pink, n~%d), Train (blue, n~%d)."),
                            n_a, n_c, n_pool, n_te, n_tr),
         x = NULL, y = NULL) +
    TH() + theme(panel.grid = element_blank(), panel.border = element_blank(),
                 axis.ticks = element_blank(), legend.position = "bottom",
                 legend.key.size = unit(0.4, "cm"))

brks <- quantile(u_p1, probs = seq(0, 1, by = 0.1))
brks[1] <- 0; brks[length(brks)] <- 1
pB <- ggplot(data.frame(u = u_p1, v = v_p1,
                        decile = cut(u_p1, breaks = brks,
                                     labels = paste0("D", 1:10),
                                     include.lowest = TRUE)),
             aes(u, v, colour = decile)) +
    geom_point(size = 0.7, alpha = 0.65) +
    geom_vline(xintercept = brks[-c(1, length(brks))], colour = "grey60",
               linetype = "dotted", linewidth = 0.3) +
    scale_colour_manual(values = DECILE_COLS, name = "Decile") +
    coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(title = "B: Stratified -- balanced folds across copula",
         subtitle = paste0("Dotted lines = decile boundaries. Each decile contributes\n",
                           "cells proportionally to every fold."),
         x = "u = Fx(GOI)", y = "v = Fc(cDFG)") +
    TH() + theme(panel.grid.minor = element_blank(), legend.position = "bottom",
                 legend.key.size = unit(0.3, "cm"),
                 legend.text = element_text(size = 7))

# =============================================================================
# Curve panels: A alone, B alone, both on the same test fold
# =============================================================================
curve_set <- function(u_pool, v_pool, tag) {
    te <- strat_folds(u_pool, K, SEED + 1000L)[[1]]
    tr <- !te; is_or <- seq_along(u_pool) <= n_a
    trA <- tr & is_or
    dfA <- step_df(u_pool[trA], v_pool[trA])
    dfB <- step_df(u_pool[tr],  v_pool[tr])

    o_i <- which(is_or); c_i <- which(!is_or)
    dA_or <- data.frame(u = u_pool[o_i], v = v_pool[o_i],
                        f = iso_predict(u_pool[trA], v_pool[trA], u_pool[o_i]))
    dB_or <- data.frame(u = u_pool[o_i], v = v_pool[o_i],
                        f = iso_predict(u_pool[tr], v_pool[tr], u_pool[o_i]))
    dB_ca <- if (length(c_i))
        data.frame(u = u_pool[c_i], v = v_pool[c_i],
                   f = iso_predict(u_pool[tr], v_pool[tr], u_pool[c_i]))
        else data.frame(u = numeric(0), v = numeric(0), f = numeric(0))
    fA <- iso_predict(u_pool[trA], v_pool[trA], u_pool[te])
    fB <- iso_predict(u_pool[tr],  v_pool[tr],  u_pool[te])
    d_te <- data.frame(u = u_pool[te], v = v_pool[te], fA = fA, fB = fB)
    rA <- sqrt(mean((v_pool[te] - fA)^2)); rB <- sqrt(mean((v_pool[te] - fB)^2))

    list(
        # == C1 / J1: OriRef only =============================================
        A = ggplot() +
            geom_segment(data = dA_or, aes(u, v, xend = u, yend = f),
                         colour = COL_A, alpha = 0.70, linewidth = 0.35) +
            geom_point(data = dA_or, aes(u, v), colour = "grey65",
                       size = 1.0, alpha = 0.35, shape = 16) +
            geom_step(data = dfA, aes(u, y), colour = COL_A, linewidth = 1.0,
                      direction = "hv") +
            coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
            labs(title = sprintf("%s1: Curve A -- train on OriRef only", tag),
                 subtitle = sprintf(paste0("Green step = Curve A; grey dots = the %d OriRef cells.\n",
                                           "Sticks = v - PAVA_A(u). Held-out fold RMSE_A = %.4f.\nSCHEMATIC: both curves are scored on one shared pool fold.\nThe GATE does not do this -- its Curve A is the OriRef-only CV over OriRef folds."),
                                    nrow(dA_or), rA),
                 x = "u = Fx(GOI)", y = "v = Fc(cDFG)") +
            TH() + theme(panel.grid.minor = element_blank()),

        # == C2 / J2: OriRef + candidates =====================================
        B = ggplot() +
            geom_segment(data = dB_or, aes(u, v, xend = u, yend = f),
                         colour = "grey55", alpha = 0.65, linewidth = 0.35) +
            geom_point(data = dB_or, aes(u, v), colour = "grey65",
                       size = 1.0, alpha = 0.32, shape = 16) +
            geom_segment(data = dB_ca, aes(u, v, xend = u, yend = f),
                         colour = COL_CAND, alpha = 0.65, linewidth = 0.35) +
            geom_point(data = dB_ca, aes(u, v), colour = COL_CAND,
                       size = 1.2, alpha = 0.35, shape = 16) +
            geom_step(data = dfB, aes(u, y), colour = COL_B, linewidth = 1.0,
                      direction = "hv") +
            coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
            labs(title = sprintf("%s2: Curve B -- train on OriRef + candidates", tag),
                 subtitle = sprintf(paste0("Navy step = Curve B; grey = OriRef, violet = the %d candidates.\n",
                                           "Sticks = v - PAVA_B(u). Held-out fold RMSE_B = %.4f.\nSCHEMATIC: both curves are scored on one shared pool fold.\nThe GATE does not do this -- its Curve A is the OriRef-only CV over OriRef folds."),
                                    nrow(dB_ca), rB),
                 x = "u = Fx(GOI)", y = "v = Fc(cDFG)") +
            TH() + theme(panel.grid.minor = element_blank()),

        # == C3 / J3: the shared test fold ====================================
        C = ggplot() +
            geom_segment(data = d_te, aes(u, v, xend = u, yend = fB),
                         colour = COL_STICK, alpha = 0.60, linewidth = 0.50) +
            geom_point(data = d_te, aes(u, v), colour = COL_TEST, size = 1.4,
                       alpha = 0.75, shape = 17) +
            geom_step(data = dfA, aes(u, y), colour = COL_A, linewidth = 1.0,
                      direction = "hv") +
            geom_step(data = dfB, aes(u, y), colour = COL_B, linewidth = 1.0,
                      direction = "hv") +
            coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
            labs(title = sprintf("%s3: Both curves on same test set", tag),
                 subtitle = sprintf(paste0("Red triangles = the %d held-out cells, sticks to Curve B.\n",
                                           "Curve A -> RMSE_A = %.4f, Curve B -> RMSE_B = %.4f.\nSCHEMATIC: both curves are scored on one shared pool fold.\nThe GATE does not do this -- its Curve A is the OriRef-only CV over OriRef folds."),
                                    nrow(d_te), rA, rB),
                 x = "u = Fx(GOI)", y = "v = Fc(cDFG)") +
            TH() + theme(panel.grid.minor = element_blank()))
}
cv1 <- curve_set(u_p1, v_p1, "C")
cvN <- curve_set(u_pN, v_pN, "J")

# =============================================================================
# Per-fold RMSE bars, RMSE heatmaps, RMSE boxplots
# =============================================================================
bar_panel <- function(rB, ttl) {
    if (identical(BAR_ROWS, "repeat1")) {
        vA   <- rmse_A[as.integer(as.character(rep_id)) == 1L]
        vB   <- rB[as.integer(as.character(rep_id)) == 1L]
        note <- "Repeat 1 only -- the heatmap row the callout marks."
    } else {
        vA   <- as.numeric(tapply(rmse_A, fold_id, mean, na.rm = TRUE))
        vB   <- as.numeric(tapply(rB,     fold_id, mean, na.rm = TRUE))
        note <- "Mean over the 40 repeats."
    }
    d <- data.frame(
        fold  = factor(rep(1:K, 2)),
        curve = rep(c("Curve A (OriRef)", "Curve B (+ cand)"), each = K),
        rmse  = c(as.numeric(vA), as.numeric(vB)))
    ggplot(d, aes(fold, rmse, fill = curve)) +
        geom_col(position = position_dodge(width = 0.7), width = 0.6) +
        scale_fill_manual(values = c("Curve A (OriRef)" = COL_A,
                                     "Curve B (+ cand)" = COL_B), name = NULL) +
        labs(title = ttl,
             subtitle = paste0(note,
                               "\nEach fold produces one RMSE_A and one RMSE_B.\n",
                               "Navy near green means the candidates cost little fit."),
             x = "Fold", y = "RMSE") +
        TH() + theme(legend.position = "bottom",
                     legend.key.size = unit(0.3, "cm"),
                     panel.grid.minor = element_blank())
}
pD <- bar_panel(rmse_B1, "D: Per-fold RMSE, first iteration")
pK <- bar_panel(rmse_BN, sprintf("K: Per-fold RMSE, iteration %d", n_it))

heat_panel <- function(rB, ttl, sub) {
    d <- data.frame(rep_id = rep(rep_id, 2), fold_id = rep(fold_id, 2),
                    rmse = c(rmse_A, rB),
                    curve = factor(rep(c("RMSE_A (OriRef)", "RMSE_B (Union)"),
                                       each = K * R),
                                   levels = c("RMSE_A (OriRef)", "RMSE_B (Union)")))
    ggplot(d, aes(fold_id, rep_id, fill = rmse)) +
        geom_tile(colour = TILE_BORDER, linewidth = TILE_LW) +
        facet_wrap(~ curve) +
        scale_fill_gradientn(colours = OCEAN_CURL, name = "RMSE") +
        scale_x_discrete(expand = expansion(0)) +
        scale_y_discrete(expand = expansion(0)) +
        labs(title = ttl, subtitle = sub, x = "Fold (1-5)", y = "Repeat (1-40)") +
        TH() + theme(panel.grid = element_blank(),
                     panel.border = element_blank(),
                     axis.ticks = element_blank(),
                     axis.text.y = element_text(size = 5),
                     strip.text = element_text(face = "bold", size = 8))
}
pE <- heat_panel(rmse_B1, "E: 200 fold-RMSEs, first iteration",
                 paste0("Each tile = one fold's RMSE. A systematically warmer right\n",
                        "panel means the union fits worse than OriRef alone."))
pL <- heat_panel(rmse_BN, sprintf("L: 200 fold-RMSEs, iteration %d", n_it),
                 sprintf(paste0("The same comparison once the gate has closed, on %s.\n",
                                "The right panel stays warmer here; removal narrows the\n",
                                "gap without having to erase it."), GATE_BY))


# == FB / MB: the two RMSE distributions, R3 cosmetics ========================
box_panel <- function(rB, pv, ttl, sub, tst, fdr_key) {
    d <- data.frame(curve = factor(rep(c("Curve A (OriRef)", "Curve B (+ cand)"),
                                       each = length(rmse_A)),
                                   levels = c("Curve A (OriRef)", "Curve B (+ cand)")),
                    rmse = c(rmse_A, rB))
    ggplot(d, aes(curve, rmse)) +
        geom_jitter(aes(fill = curve), width = 0.20, size = 1.7, alpha = 0.38,
                    shape = 21, stroke = 0.30, colour = "grey35") +
        geom_box_std(aes(fill = curve), width = 0.55, colour = "black",
                     linewidth = BOX_LW, outlier.shape = NA, alpha = 0.35) +
        scale_fill_manual(values = c(COL_A, COL_CAND), guide = "none") +
        labs(title = ttl,
             subtitle = sprintf(paste0("%s\nGate statistic: one-sided paired Wilcoxon (H1: RMSE_B > RMSE_A) p = %s.\n",
                                       "Two-sided %s on the same pairs: FDR = %s (BH over the 3 paired tests here).\n",
                                       "The two curves use different fold partitions (A over OriRef, B over the pool);\n",
                                       "they are paired by (repeat, fold) index, not by shared test cells."),
                                sub, format(pv, scientific = TRUE, digits = 3),
                                tst$name, fdr_txt(fdr_key)),
             x = NULL, y = "Per-fold isotonic RMSE") +
        TH() + theme(panel.grid.minor = element_blank())
}
pFB <- box_panel(rmse_B1, p1, "FB: RMSE distributions, first iteration",
                 "The full candidate set raises the per-fold RMSE.", t_F, "F")
pMB <- box_panel(rmse_BN, pN, sprintf("MB: RMSE distributions, iteration %d", n_it),
                 sprintf("At gate closure. The gate closed on %s.", GATE_BY),
                 t_M, "M")

# == F / M: the effect-size condition, drawn ==================================
# dRMSE > BV reduces to median(RMSE_B - RMSE_A) > IQR(RMSE_A)
diff_panel <- function(rB, pv, ttl, sub, tst, fdr_key, col_fill, note = "") {
    d   <- data.frame(g = "", d = rB - rmse_A)
    d   <- d[is.finite(d$d), ]
    med <- median(d$d)
    n_hi <- sum(d$d > IQR_A)
    ok  <- med > IQR_A
    ggplot(d, aes(g, d)) +
        geom_hline(yintercept = 0, colour = COL_ZERO, linetype = "solid",
                   linewidth = 0.4) +
        geom_hline(yintercept = IQR_A, colour = COL_A, linetype = "dashed",
                   linewidth = 0.7) +
        geom_jitter(width = 0.20, size = 1.7, alpha = 0.38, shape = 21,
                    stroke = 0.30, colour = "grey35", fill = col_fill) +
        geom_box_std(width = 0.55, colour = "black", linewidth = BOX_LW,
                     outlier.shape = NA, alpha = 0.35, fill = col_fill) +
        annotate("segment", x = 0.55, xend = 1.45, y = med, yend = med,
                 colour = "black", linewidth = 1.3) +
        annotate("text", x = 1.48, y = med, hjust = 0, size = 2.9,
                 fontface = "bold",
                 label = sprintf("median = %.5f", med)) +
        annotate("text", x = 1.48, y = IQR_A, hjust = 0, size = 2.9,
                 colour = COL_A, fontface = "italic",
                 label = sprintf("IQR(RMSE_A) = %.5f", IQR_A)) +
        annotate("text", x = 1.48, y = 0, hjust = 0, size = 2.7,
                 colour = COL_ZERO, label = "no difference") +
        scale_x_discrete(expand = expansion(add = c(0.6, 1.9))) +
        labs(title = ttl,
             subtitle = sprintf(paste0("%s%s\n",
                                       "EFFECT: median(B - A) = %.5f vs IQR(RMSE_A) = %.5f -> %s. %d/%d folds above the line.\n",
                                       "CONSISTENCY: one-sided paired Wilcoxon p = %s; two-sided %s FDR = %s."),
                                sub, note, med, IQR_A,
                                if (ok) "MATERIAL, removal continues"
                                else "WITHIN BACKGROUND, gate closes",
                                n_hi, nrow(d),
                                format(pv, scientific = TRUE, digits = 3),
                                tst$name, fdr_txt(fdr_key)),
             x = NULL, y = "RMSE_B - RMSE_A, per fold") +
        TH() + theme(panel.grid.minor = element_blank(),
                     axis.text.x = element_blank(),
                     axis.ticks.x = element_blank())
}
pF <- diff_panel(rmse_B1, p1, "F: First iteration",
                 "All 352 candidates in the pool.",
                 t_F, "F", COL_DIFF,
                 if (RECON_OK) "" else
                     " Curve B here is RECONSTRUCTED and did not reproduce the stored closure vector.")
pM <- diff_panel(rmse_BN, pN, sprintf("M: Iteration %d (gate closure)", n_it),
                 sprintf("The %d accepted candidates. Curve B is Part 10e's stored final union.",
                         sum(kept)),
                 t_M, "M", COL_DIFF)

# == CI: is the margin bigger than resampling error? ==========================
# Fold-level bootstrap of the margin, resampling the paired folds together
boot_margin <- function(rB, B = 2000L) {
    ok <- is.finite(rmse_A) & is.finite(rB)
    a <- rmse_A[ok]; b <- rB[ok]; n <- length(a)
    set.seed(SEED)
    m <- vapply(seq_len(B), function(i) {
        k <- sample.int(n, n, replace = TRUE)
        median(b[k] - a[k]) - IQR(a[k])
    }, numeric(1))
    list(draws = m, est = median(b - a) - IQR(a),
         lo = unname(quantile(m, 0.025)), hi = unname(quantile(m, 0.975)))
}
bm1 <- boot_margin(rmse_B1); bmN <- boot_margin(rmse_BN)
ci1 <- c(est = bm1$est, lo = bm1$lo, hi = bm1$hi)
ciN <- c(est = bmN$est, lo = bmN$lo, hi = bmN$hi)
STAGE_LV <- c("First iteration", sprintf("Iteration %d", n_it))
d_draw <- data.frame(
    stage = factor(rep(STAGE_LV, c(length(bm1$draws), length(bmN$draws))),
                   levels = STAGE_LV),
    m = c(bm1$draws, bmN$draws))
d_ci <- data.frame(stage = factor(STAGE_LV, levels = STAGE_LV),
                   est = c(ci1[["est"]], ciN[["est"]]),
                   lo  = c(ci1[["lo"]],  ciN[["lo"]]),
                   hi  = c(ci1[["hi"]],  ciN[["hi"]]))
cat(sprintf(paste0("  Margin median(B-A) - IQR(RMSE_A), fold bootstrap 2000x:\n",
                   "    iteration 1  %+.5f [%+.5f, %+.5f]\n",
                   "    iteration %-2d %+.5f [%+.5f, %+.5f]\n"),
            ci1[["est"]], ci1[["lo"]], ci1[["hi"]],
            n_it, ciN[["est"]], ciN[["lo"]], ciN[["hi"]]))

pCI <- ggplot(d_draw, aes(stage, m)) +
    geom_hline(yintercept = 0, colour = COL_ZERO, linetype = "solid",
               linewidth = 0.4) +
    geom_jitter(aes(fill = stage), width = 0.22, size = 1.7, alpha = 0.10,
                shape = 21, stroke = 0.20, colour = "grey35") +
    geom_errorbar(data = d_ci, aes(x = stage, ymin = lo, ymax = hi),
                  width = 0.14, colour = "black", linewidth = 0.55,
                  inherit.aes = FALSE) +
    geom_point(data = d_ci, aes(stage, est, fill = stage), size = 3.4,
               shape = 21, stroke = 0.55, colour = "grey35",
               inherit.aes = FALSE) +
    scale_fill_manual(values = c(COL_A, COL_CAND), guide = "none") +
    labs(title = "CI: margin over background variability",
         subtitle = sprintf(paste0("Margin = median(RMSE_B - RMSE_A) - IQR(RMSE_A). Above the brown line the\n",
                                   "effect-size condition holds; at or below it the gate closes. Faint dots =\n",
                                   "the 2000 bootstrap draws; solid dot = estimate, bar = 95%% percentile CI.\n",
                                   "Fold-level bootstrap, 2000 resamples: %+.5f [%+.5f, %+.5f] -> %+.5f [%+.5f, %+.5f]."),
                            ci1[["est"]], ci1[["lo"]], ci1[["hi"]],
                            ciN[["est"]], ciN[["lo"]], ciN[["hi"]]),
         x = NULL, y = "median(B - A) - IQR(RMSE_A)") +
    TH() + theme(panel.grid.minor = element_blank())

# =============================================================================
# G / N: PAVA residual ranking     H: PAVA on the initial union
# =============================================================================
pava_step <- step_df(u_p1, v_p1)
rank_panel <- function(hl, ttl, sub) {
    o <- order(resid_pv, decreasing = TRUE)
    ggplot(data.frame(rank = seq_along(o), r = resid_pv[o], hl = hl[o]),
           aes(rank, r, fill = hl)) +
        geom_col(width = 1, colour = "grey35", linewidth = 0.005) +
        scale_fill_manual(values = c(`FALSE` = COL_BAR_ACC, `TRUE` = COL_REMOVED),
                          guide = "none") +
        labs(title = ttl, subtitle = sub,
             x = "Rank, descending |v - PAVA(u)|", y = "|v - PAVA(u)|") +
        TH() + theme(panel.grid.minor = element_blank())
}
hl1 <- rep(FALSE, n_c); if (length(rem_ord)) hl1[rem_ord[1]] <- TRUE
pG <- rank_panel(hl1, "G: Candidates ranked by PAVA residual",
                 paste0("PAVA computed once on the initial union, so removal is\n",
                        "decoupled from spine admission. Dark bar = removed this iteration."))
pN_ <- rank_panel(!kept, "N: Final ranking",
                  sprintf("Removed %d/%d (%.1f%%, dark bars); accepted in purple.",
                          s2$n_removed, n_c, 100 * s2$n_removed / n_c))

# X marks the first and last removed cells
pH <- ggplot() +
    geom_point(data = data.frame(u = u_a, v = v_a), aes(u, v),
               colour = "grey75", size = 0.3, alpha = 0.3) +
    geom_point(data = data.frame(u = u_c, v = v_c), aes(u, v),
               shape = 21, size = 1.6, stroke = 0.35,
               fill = COL_CAND, colour = COL_CAND_D, alpha = 0.28) +
    geom_step(data = pava_step, aes(u, y), colour = COL_B, linewidth = 0.6,
              alpha = 0.7, direction = "hv") +
    coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(title = "H: PAVA precomputed on the initial union",
         subtitle = sprintf(paste0("Navy step = PAVA on OriRef + all %d candidates, computed once.\n",
                                   "Violet = candidates. X = first and last cell removed."), n_c),
         x = "u = Fx(GOI)", y = "v = Fc(cDFG)") +
    TH() + theme(panel.grid.minor = element_blank())
if (length(rem_ord)) {
    mk_i <- unique(c(rem_ord[1], rem_ord[length(rem_ord)]))
    d_mk <- data.frame(u = u_c[mk_i], v = v_c[mk_i],
                       f = iso_predict(u_p1, v_p1, u_c[mk_i]),
                       lab = c("first removed", "last removed")[seq_along(mk_i)])
    pH <- pH +
        geom_segment(data = d_mk, aes(u, v, xend = u, yend = f),
                     colour = "black", linewidth = 0.5) +
        geom_point(data = d_mk, aes(u, v), colour = COL_TEST, size = 3.4,
                   shape = 4, stroke = 1.6) +
        geom_text(data = d_mk, aes(u, v, label = lab), colour = COL_TEST,
                  fontface = "bold", size = 2.4,
                  vjust = c(2.2, -1.6)[seq_len(nrow(d_mk))],
                  hjust = c(1.05, -0.05)[seq_len(nrow(d_mk))])
}

# =============================================================================
# I: gate convergence on p
# =============================================================================
d_cv <- data.frame(iter = seq_len(n_it), p = p_iter)
pI <- ggplot(d_cv, aes(iter, p)) +
    geom_hline(yintercept = ALPHA, colour = COL_A, linetype = "dashed",
               linewidth = 0.7) +
    annotate("text", x = 1, y = ALPHA + 0.008,
             label = sprintf("alpha = %.2f", ALPHA), colour = COL_A,
             hjust = 0, size = 2.8, fontface = "italic") +
    geom_line(colour = COL_B, linewidth = 0.6) +
    geom_point(colour = COL_B, size = 0.8, alpha = 0.6) +
    geom_point(data = d_cv[n_it, ], colour = COL_A, size = 3.0, shape = 16) +
    annotate("text", x = n_it, y = p_iter[n_it] + 0.012,
             label = sprintf("p=%.4f\n%d accepted", p_iter[n_it], sum(kept)),
             colour = COL_A, hjust = 0.5, size = 2.5, fontface = "bold") +
    scale_x_continuous(breaks = seq(0, n_it, by = max(1, n_it %/% 10))) +
    labs(title = "I: Gate convergence -- p-value vs iteration",
         subtitle = sprintf(paste0("Each iteration removes 1 candidate (highest PAVA residual) + re-runs the 200-fold CV.\n",
                                   "%d removals. Green dashed = alpha. Green dot = gate closed, on %s --\n",
                                   "p need not have reached alpha -- see the dRMSE page."),
                            s2$n_removed, GATE_BY),
         x = "Iteration (each = 1 removal + 200-fold RSk-fCV)",
         y = "Wilcoxon p-value") +
    TH() + theme(panel.grid.minor = element_blank())

# =============================================================================
# R1-R4: the dRMSE / background-variability arm of the gate
# =============================================================================
q1 <- quantile(rmse_A, 0.25, na.rm = TRUE); q3 <- quantile(rmse_A, 0.75, na.rm = TRUE)
pR1 <- ggplot(data.frame(r = rmse_A), aes(r)) +
    geom_histogram(bins = 40, fill = COL_A, colour = "white", linewidth = 0.2,
                   alpha = 0.75) +
    geom_vline(xintercept = med_A, colour = COL_B, linewidth = 0.7) +
    geom_vline(xintercept = c(q1, q3), colour = COL_B, linetype = "dashed",
               linewidth = 0.5) +
    labs(title = "R1: Background variability of the reference",
         subtitle = sprintf(paste0("RMSE_A alone, one value per fold. Navy solid = median, dashed = quartiles.\n",
                                   "BV = IQR/median = %.5f/%.5f = %.4f (%.2f%%)."),
                            q3 - q1, med_A, BV, 100 * BV),
         x = "RMSE_A", y = "Count") +
    TH() + theme(panel.grid.minor = element_blank())

pR2 <- ggplot(data.frame(iter = seq_len(n_it), d = d_iter), aes(iter, d)) +
    geom_ribbon(aes(ymin = -Inf, ymax = BV), fill = COL_A, alpha = 0.10) +
    geom_hline(yintercept = 0, colour = COL_ZERO, linetype = "solid",
               linewidth = 0.4) +
    geom_hline(yintercept = BV, colour = COL_A, linetype = "dashed",
               linewidth = 0.7) +
    annotate("text", x = 1, y = BV, label = sprintf("BV = %.4f", BV),
             colour = COL_A, hjust = 0, vjust = -0.6, size = 2.8,
             fontface = "italic") +
    geom_line(colour = COL_B, linewidth = 0.6) +
    geom_point(colour = COL_B, size = 0.8, alpha = 0.6) +
    geom_point(data = data.frame(iter = n_it, d = d_iter[n_it]),
               colour = COL_A, size = 3.0, shape = 16) +
    annotate("text", x = n_it, y = d_iter[n_it],
             label = sprintf("dRMSE = %.4f", d_iter[n_it]), colour = COL_A,
             hjust = 1.05, vjust = 1.8, size = 2.5, fontface = "bold") +
    scale_x_continuous(breaks = seq(0, n_it, by = max(1, n_it %/% 10))) +
    labs(title = "R2: Effect size vs iteration",
         subtitle = sprintf(paste0("dRMSE = median(RMSE_B - RMSE_A) / median(RMSE_A).\n",
                                   "Shaded = inside the reference's own fold-to-fold spread.\n",
                                   "Gate closed on %s."), GATE_BY),
         x = "Iteration (each = 1 removal + 200-fold RSk-fCV)", y = "dRMSE") +
    TH() + theme(panel.grid.minor = element_blank())

d_r3 <- data.frame(
    stage = factor(rep(c("First iteration", sprintf("Iteration %d", n_it)),
                       each = K * R),
                   levels = c("First iteration", sprintf("Iteration %d", n_it))),
    d = c(d_first, d_last))
pR3 <- ggplot(d_r3, aes(stage, d)) +
    geom_hline(yintercept = 0, colour = COL_ZERO, linetype = "solid",
               linewidth = 0.4) +
    geom_hline(yintercept = BV, colour = COL_A, linetype = "dashed",
               linewidth = 0.7) +
    geom_jitter(aes(fill = stage), width = 0.20, size = 1.7, alpha = 0.38,
                shape = 21, stroke = 0.30, colour = "grey35") +
    geom_box_std(aes(fill = stage), width = 0.55, colour = "black",
                 linewidth = BOX_LW, outlier.shape = NA, alpha = 0.35) +
    scale_fill_manual(values = c(COL_A, COL_CAND), guide = "none") +
    labs(title = "R3: Paired per-fold differences",
         subtitle = sprintf(paste0("The 200 values the median is taken over. Medians %.4f -> %.4f\n",
                                   "against BV %.4f (green dashed). Paired on fold identity:\n",
                                   "%s, FDR = %s (BH over the 3 paired tests in this figure)."),
                            d1, dN, BV, t_R3$name, fdr_txt("R3")),
         x = NULL, y = "(RMSE_B - RMSE_A) / median(RMSE_A)") +
    TH() + theme(panel.grid.minor = element_blank())

d_r4 <- data.frame(p = pmax(p_iter, 1e-300), d = d_iter, iter = seq_len(n_it))
# Zoom from the 2nd percentile of p to just past alpha
r4_lo   <- as.numeric(quantile(d_r4$p, 0.02, na.rm = TRUE))
r4_hi   <- max(ALPHA * 5, max(d_r4$p, na.rm = TRUE))
r4_noff <- sum(d_r4$p < r4_lo, na.rm = TRUE)
pR4 <- ggplot(d_r4, aes(p, d)) +
    annotate("rect", xmin = 1e-300, xmax = ALPHA, ymin = BV, ymax = Inf,
             fill = COL_TEST, alpha = 0.10) +
    geom_hline(yintercept = 0, colour = COL_ZERO, linetype = "solid",
               linewidth = 0.4) +
    geom_vline(xintercept = ALPHA, colour = COL_A, linetype = "dashed",
               linewidth = 0.7) +
    geom_hline(yintercept = BV, colour = COL_A, linetype = "dashed",
               linewidth = 0.7) +
    geom_path(colour = COL_B, linewidth = 0.4, alpha = 0.7) +
    geom_point(aes(colour = iter), size = 1.2) +
    geom_point(data = d_r4[n_it, ], colour = COL_A, size = 3.0, shape = 16) +
    scale_colour_viridis_c(option = "plasma", name = "Iteration") +
    scale_x_log10() +
    # coord_cartesian zooms without dropping points
    coord_cartesian(xlim = c(r4_lo, r4_hi)) +
    labs(title = "R4: Decision plane -- both conditions",
         subtitle = sprintf(paste0("Removal needs p < alpha AND dRMSE > BV (shaded).\n",
                                   "Path runs iteration 1 -> %d; green dot = gate closed, on %s.%s"),
                            n_it, GATE_BY,
                            if (r4_noff > 0)
                                sprintf(" %d early iteration(s) off-scale left.", r4_noff)
                            else ""),
         x = "Wilcoxon p (log scale)", y = "dRMSE") +
    TH() + theme(panel.grid.minor = element_blank(), legend.position = "bottom",
                 legend.key.width = unit(1, "cm"),
                 legend.key.height = unit(0.3, "cm"))

# =============================================================================
# Output: 2x2 pages, console then PDF
# =============================================================================
pages <- list(
    "CV design and stratification"  = list(pA, pB, cv1$A, cv1$B),
    "First iteration"               = list(cv1$C, pD, pG, pE),
    "PAVA removal and convergence"  = list(pH, pI, pN_, cvN$A),
    "Final iteration"               = list(cvN$B, cvN$C, pK, pL),
    "Effect size vs background variability" = list(pF, pM, pCI),
    "dRMSE arm of the gate"   = list(pR1, pR2, pR4),
    "RMSE and dRMSE distributions"  = list(pFB, pMB, pR3)
)

# == Force each panel while its data are in scope =============================
pages <- lapply(pages, function(ps) lapply(ps, function(g) { force(g); g }))

# == Workbook: the numbers behind every panel =================================
xl <- list()
xl[["Summary"]] <- data.frame(
    arm = TAG, dataset = DS, n_cells = nrow(df),
    n_origref = n_a, n_candidates = n_c,
    n_removed = s2$n_removed, n_recovered = sum(kept),
    iterations = n_it, k_folds = K, n_repeats = R, n_fold_rmse = K * R,
    alpha = ALPHA, wilcox_p_final = s2$wilcox_p_final,
    dRMSE_final = s2$relative_delta_final, BV = BV,
    median_RMSE_A = med_A, IQR_RMSE_A = IQR_A,
    margin_final = ciN[["est"]], margin_lo = ciN[["lo"]], margin_hi = ciN[["hi"]],
    gate_closed_on = GATE_BY,
    curveB_reconstruction_ok = RECON_OK,
    stringsAsFactors = FALSE)
xl[["Fold_RMSE"]] <- data.frame(
    repeat_id = as.integer(as.character(rep_id)),
    fold_id   = as.integer(as.character(fold_id)),
    RMSE_A = rmse_A, RMSE_B_iter1 = rmse_B1, RMSE_B_final = rmse_BN,
    diff_iter1 = rmse_B1 - rmse_A, diff_final = rmse_BN - rmse_A,
    stringsAsFactors = FALSE)
xl[["Gate_Trajectory"]] <- data.frame(
    iteration = seq_len(n_it), wilcox_p = p_iter, dRMSE = d_iter,
    n_active = s2$n_per_iter, BV = BV, alpha = ALPHA,
    stringsAsFactors = FALSE)
.ord <- order(resid_pv, decreasing = TRUE)
xl[["Candidates"]] <- data.frame(
    rank_by_pava = seq_along(.ord),
    candidate_index = .ord,
    resid_pava = resid_pv[.ord],
    accepted = kept[.ord],
    removal_step = match(.ord, rem_ord),
    u = u_c[.ord], v = v_c[.ord],
    stringsAsFactors = FALSE)
xl[["Paired_Tests"]] <- data.frame(
    panel = c("F", "M", "R3"),
    comparison = c("RMSE_B vs RMSE_A, iteration 1",
                   sprintf("RMSE_B vs RMSE_A, iteration %d", n_it),
                   sprintf("dRMSE, iteration 1 vs %d", n_it)),
    test = c(t_F$name, t_M$name, t_R3$name),
    shapiro_p = c(t_F$shapiro, t_M$shapiro, t_R3$shapiro),
    p_value = c(t_F$p, t_M$p, t_R3$p),
    FDR_BH = as.numeric(FDR[c("F", "M", "R3")]),
    stringsAsFactors = FALSE)
xl[["Bootstrap_CI"]] <- data.frame(
    stage = STAGE_LV, estimate = d_ci$est, ci_lo = d_ci$lo, ci_hi = d_ci$hi,
    n_resamples = 2000L, quantity = "median(RMSE_B - RMSE_A) - IQR(RMSE_A)",
    stringsAsFactors = FALSE)

list(pages = pages, xl = xl, s2 = s2, GATE_BY = GATE_BY, n_it = n_it,
     ALPHA = ALPHA, BV = BV)
}

# == Plot letters: A, B, C ... in order of appearance =========================
letter_pages <- function(pages) {
    k <- 0L
    idx <- list()
    out <- lapply(names(pages), function(pg) {
        lapply(pages[[pg]], function(g) {
            k <<- k + 1L
            if (k > length(LETTERS))
                stop("more than 26 panels; the letter scheme needs extending.",
                     call. = FALSE)
            old <- if (!is.null(g$labels$title)) g$labels$title else ""
            body <- sub("^[A-Za-z0-9_]+:[[:space:]]*", "", old)
            idx[[length(idx) + 1L]] <<- data.frame(
                plot = LETTERS[k], page = pg, title = body,
                stringsAsFactors = FALSE)
            g + labs(title = paste0("Plot ", LETTERS[k], ": ", body))
        })
    })
    names(out) <- names(pages)
    list(pages = out, index = do.call(rbind, idx))
}

draw_page <- function(ps, ttl, tag)
    grid.arrange(grobs = ps, ncol = 2, nrow = ceiling(length(ps) / 2),
                 top = textGrob(sprintf("Figure 10 -- %s (%s)", ttl, tag),
                                gp = gpar(fontface = "bold", fontsize = 10)))

for (tg in TARGETS) {
    rds <- file.path(CEP_HOME, "Part10e_CMRF", tg$dir, tg$rds)
    if (!file.exists(rds)) { cat("SKIP (not found):", rds, "\n"); next }
    cat(sprintf("\n================ %s ================\n", tg$tag))
    r <- make_pages(rds, tg$ds, tg$tag)

    lp    <- letter_pages(r$pages)
    pgs   <- lp$pages
    r$xl[["Panel_Index"]] <- lp$index

    out_dir <- file.path(OUT_ROOT, paste0("Part10f_RSkfCV_Mechanism_", tg$sub))
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    base <- file.path(out_dir, paste0(tg$pre, "_", tg$ds))

    # Draw to the console device
    if (isTRUE(SHOW)) {
        for (nm in names(pgs)) {
            cat("  page:", nm, "\n")
            grid.arrange(grobs = pgs[[nm]], ncol = 2,
                         nrow = ceiling(length(pgs[[nm]]) / 2),
                         top = textGrob(sprintf("Figure 10 -- %s (%s)", nm, tg$tag),
                                        gp = gpar(fontface = "bold", fontsize = 10)))
        }
    }
    # Pages kept for redrawing
    FIG10_PAGES <- pgs

    # One multi-page PDF
    pdf_path <- paste0(base, "_Figure10_mechanism.pdf")
    tryCatch({
        pdf(pdf_path, width = PDF_W, height = PDF_H, onefile = TRUE)
        on.exit(try(dev.off(), silent = TRUE), add = TRUE)
        for (nm in names(pgs)) draw_page(pgs[[nm]], nm, tg$tag)
        dev.off(); on.exit()
        cat("  PDF : ", basename(pdf_path), "\n", sep = "")
    }, error = function(e) cat("  PDF FAILED: ", conditionMessage(e), "\n", sep = ""))

    # One PNG per page
    for (i in seq_along(pgs)) {
        png_path <- sprintf("%s_Figure10_page%d.png", base, i)
        tryCatch({
            png(png_path, width = PNG_W, height = PNG_H, res = PNG_RES)
            on.exit(try(dev.off(), silent = TRUE), add = TRUE)
            draw_page(pgs[[i]], names(pgs)[i], tg$tag)
            dev.off(); on.exit()
        }, error = function(e)
            cat("  PNG page ", i, " FAILED: ", conditionMessage(e), "\n", sep = ""))
    }
    cat("  PNG : ", length(pgs), " page file(s)\n", sep = "")

    xlsx_path <- paste0(base, "_RSkfCV_Mechanism.xlsx")
    tryCatch({
        write_xlsx(r$xl, path = xlsx_path)
        cat("  Excel: ", basename(xlsx_path), " (", length(r$xl), " sheets)\n",
            sep = "")
    }, error = function(e)
        cat("  Excel FAILED: ", conditionMessage(e), "\n", sep = ""))

    cat(sprintf("  Folder: %s\n", normalizePath(out_dir, mustWork = FALSE)))
    cat(sprintf(paste0("Gate closed at iteration %d.\n",
                       "  Wilcoxon p = %.3g (alpha %.2f) | dRMSE = %.4f (BV %.4f)\n",
                       "  Closing criterion: %s\n"),
                r$n_it, r$s2$wilcox_p_final, r$ALPHA,
                r$s2$relative_delta_final, r$BV, r$GATE_BY))
}
