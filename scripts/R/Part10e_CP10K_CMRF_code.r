# =============================================================================
# Part10e_CP10K_CMRF_code.r
# SEA-AD CEP-IP/CMRF, Part 10e, Copula Monotone Recovery Filter (CMRF)
#   (Table S10; the CMRF classes behind Figures 6 to 8)
#   For each arm and dataset: fits the spine g(u) on TREP SKCS-core cells,
#   calibrates delta_base, admits cells within the adaptive band (Stage 1),
#   tests the candidates with the RSk-fCV fidelity gate and removes the worst
#   by PAVA residual until the gate closes (Stage 2), then classifies cells as
#   Original + Refined, Original only, Recovered or Excluded.
#
# Usage (from the SEA-AD folder, or set SEAAD_ROOT), after Part 10,
#   Part10b_CP10K_IPRS_code.r, Part10b_CP10K_SKCS_TREP_code.r and Part 10d
#   source("Part10e_CP10K_CMRF_code.r")                         # all arms
#   CEP_ARM <- "PVM1_all"; source("Part10e_CP10K_CMRF_code.r")  # one arm
#   Needs CEP_arms.R, R_config.R and R_config_HKG8.R in the same folder.
#
# Outputs (Part10e_CMRF/Part10e_CMRF_<arm>/): figures (PDF, PNG), SK
#   interaction and RSk-fCV audit workbooks per dataset, verification
#   workbook, results RDS, guardrail ledger and logs.
# =============================================================================


# =============================================================================
# SECTION 0: Libraries and shared arm definitions
# =============================================================================

library(mgcv)
library(Matrix)        # readMM for the normalisation guardrails
library(ggplot2)
library(dplyr)
library(tidyr)
library(writexl)
library(scales)
library(gridExtra)
library(grid)
library(gtable)
library(cowplot)
library(MASS)          # kde2d
library(clinfun)       # jonckheere.test
library(Hmisc)         # hoeffd()
library(energy)        # dcor.test()

set.seed(42)

# CEP_arms.R holds the arm table, the folder layout and the driver.
if (!exists("CEP_HOME")) CEP_HOME <- Sys.getenv("SEAAD_ROOT", getwd())
.p10e_harness <- file.path(CEP_HOME, "CEP_arms.R")
if (!file.exists(.p10e_harness))
    stop("CEP_arms.R not found at: ", .p10e_harness,
         "\n  It holds the arm table shared by Parts 9 onwards.",
         call. = FALSE)
source(.p10e_harness)

PART_TAG     <- "Part10e"        # file-name tag
PART_DIR_TAG <- "Part10e_CMRF"   # results-folder tag

# == Estimator settings for the spine ==========================================
SPINE_COL_10E           <- "#12369E"  # navy: the fitted spine g(u)
IP_U_COL_10E             <- "#008080"  # teal: IP_u marker lines

GAM_DONOR_RE_10E        <- TRUE   # donor random intercept on the spine
DONOR_RE_MIN_LEVELS_10E <- 3L     # below this the fit falls back to cell level
W_DONOR_WEIGHTS_10E     <- TRUE   # inverse-donor weights
GAM_METHOD_10E          <- "REML"
GAM_SELECT_10E          <- TRUE
GAM_GAMMA_10E           <- 1.5
PRSS_PHASE2_MAX_P10E    <- 5L     # phase 1 = 8 iterations, phase 2 = at most 5
# Phase 1 of the PRSS loop (8 independent fits) may run on PSOCK workers;
# phase 2 is sequential. Results are identical either way.
P10E_PARALLEL           <- TRUE
P10E_N_CORES            <- max(1L, parallel::detectCores() - 1L)
# Workers are used only if one sequential run took at least this many seconds.
P10E_PAR_MIN_SECONDS    <- 2.0
PRSS_EARLY_STOP_P10E    <- 20L    # unreachable under the phase-2 cap
MIN_SPINE_CELLS         <- 50L    # below this the spine is not fitted

# == CMRF parameters =========================================================
CMRF_COVERAGE  <- 0.95   # delta_base quantile level: p<0.05 two-tailed equivalent
CMRF_KAPPA     <- 0.50   # adaptive widening at the waist
CMRF_BW        <- 0.10   # Gaussian kernel bandwidth
CMRF_N_REPEATS <- 40L    # RSk-fCV repeats for the Stage 2 fidelity gate
CMRF_K_FOLDS   <- 5L     # RSk-fCV folds (5 x 40 = 200 fold RMSEs)
CMRF_ALPHA     <- 0.05   # Wilcoxon alpha for the gate
CMRF_JT_NPERM  <- 10000L # Jonckheere-Terpstra permutations
CMRF_HD_NPERM  <- 9999L  # Hoeffding D permutations


# =============================================================================
# SECTION 0d: Pure helpers
# =============================================================================
# Helpers that depend only on their arguments.


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

# == theme_vp(): theme for the verification plots VP1 to VP3 ==================
VP_PANEL_BG    <- "white"     # white panel background (default ggplot2 style)
VP_PLOT_BG     <- "white"
VP_BORDER_COL  <- "grey20"
VP_COL_ORIGREF   <- "#009E73"   # Orig+Ref -- Post-IP TREP green (from Plot 7)
VP_COL_UNION   <- "#006B4E"   # Final union (darker green to distinguish from Orig+Ref)
VP_COL_ORIGREF_T1<- "#009E73"   # Orig+Ref alone in T1 (same as Orig+Ref)
VP_COL_RECOV   <- "#7C3AED"   # Recovered -- vivid purple (from Plot 7 RECOV_COL)
VP_COL_ORIONLY <- "#A0522D"   # Original only -- tan/sienna (from Plot 7 Orig only)
VP_COL_EXCL    <- "#8A9BAD"   # Excluded (grey, not in Plot 7)
VP_VIO_OUTLINE <- "grey30"
VP_VIO_LW      <- 0.50
VP_VIO_ALPHA   <- 0.55
VP_BOX_FILL    <- "white"
VP_BOX_COL     <- "black"
VP_BOX_LW      <- 0.55
VP_BOX_W       <- 0.22
VP_OL_SZ       <- 0.6
VP_OL_ALPHA    <- 0.35
VP_DLINE_COL   <- "#1A237E"   # navy for delta_base line and 0.5 hlines
VP_DLINE_LW    <- 0.75
VP_TAG_BG      <- "#D0D8E8"
VP_TAG_COL     <- "black"
VP_TAG_SZ      <- 2.8

theme_vp <- function(base_size = 9.5) {
    theme_grey(base_size = base_size) %+replace%
        theme(
            panel.background  = element_rect(fill = "white", colour = NA),
            panel.border      = element_rect(fill = NA, colour = "grey20",
                                             linewidth = 0.45),
            panel.grid.major  = element_line(colour = "grey85", linewidth = 0.45),
            panel.grid.minor  = element_line(colour = "grey92", linewidth = 0.25),
            axis.title        = element_text(face = "bold", size = base_size - 0.5),
            axis.text         = element_text(size = base_size - 1.5, colour = "grey20"),
            axis.ticks        = element_line(colour = "grey50", linewidth = 0.35),
            plot.title        = element_text(face = "bold", size = base_size + 0.5,
                                             margin = margin(b = 2)),
            plot.subtitle     = element_text(size = base_size - 1.5, colour = "grey30",
                                             margin = margin(b = 3)),
            strip.background  = element_rect(fill = "grey85", colour = NA),
            strip.text        = element_text(face = "bold", size = base_size - 0.5),
            legend.background = element_rect(fill = "white", colour = NA),
            legend.key        = element_rect(fill = "white", colour = NA),
            legend.text       = element_text(size = base_size - 1.5),
            plot.background   = element_rect(fill = VP_PLOT_BG, colour = NA),
            plot.margin       = margin(4, 6, 4, 6)
        )
}

# vp_tag() returns the plot unchanged.
vp_tag <- function(p, label) { p }

save_pp <- function(p, base, pdf_w = 10, pdf_h = 9,
                    png_w = 2200, png_h = 2000, res = 200) {
    # Each render is wrapped, so a failure costs only that figure.
    render_to_device <- function(open_device_fn, draw_fn, what = "figure") {
        tryCatch({
            open_device_fn()
            on.exit(dev.off(), add = TRUE)
            draw_fn()
        }, error = function(e) {
            cat("  [save_pp] ", what, " failed for ", basename(base), ": ",
                conditionMessage(e), "\n", sep = "")
        })
    }
    is_grid_obj <- inherits(p, "grob") || inherits(p, "gtable") ||
        inherits(p, "arrangelist")
    if (is_grid_obj) {
        render_to_device(
            function() pdf(paste0(base, ".pdf"), width = pdf_w, height = pdf_h),
            function() { grid.newpage(); grid.draw(p) }, "PDF")
        render_to_device(
            function() png(paste0(base, ".png"), width = png_w, height = png_h, res = res),
            function() { grid.newpage(); grid.draw(p) }, "PNG")
    } else {
        render_to_device(
            function() pdf(paste0(base, ".pdf"), width = pdf_w, height = pdf_h),
            function() print(p), "PDF")
        render_to_device(
            function() png(paste0(base, ".png"), width = png_w, height = png_h, res = res),
            function() print(p), "PNG")
    }
}

compute_empirical_copula_at <- function(p, q, p_eval, q_eval) {
    n <- length(p)
    vapply(seq_along(p_eval), function(k) {
        sum(p <= p_eval[k] & q <= q_eval[k]) / n
    }, numeric(1L))
}

FOUR_COLS <- c(
    "Pre-IP TREP"      = "#E69F00",
    "Pre-IP non-TREP"  = "#1B3A8C",
    "Post-IP TREP"     = "#009E73",
    "Post-IP non-TREP" = "#D94F4F"
)


##############################################################################
# SPINE FITTING -- CEP-IP spine g(u), mgcv with donor random intercept
##############################################################################
# Unconstrained penalised spline with donor random intercept, inverse-donor
# weights and REML; monotonicity is checked per fit.

DONOR_RE_LABEL_10E <- "s(donor)"
SPINE_LABEL_10E    <- "s(u_trep)"

# Donor as a factor with unused levels dropped.
as_donor_factor_10e <- function(x) droplevels(factor(as.character(x)))

n_donor_levels_10e <- function(x) {
    if (is.null(x)) return(0L)
    nlevels(as_donor_factor_10e(x))
}

donor_re_ok_10e <- function(x)
    isTRUE(GAM_DONOR_RE_10E) && n_donor_levels_10e(x) >= DONOR_RE_MIN_LEVELS_10E

# w = 1 / (cells this donor contributes to this frame), as a plain numeric
# vector; pass weights = w by bare column name.
donor_weights_10e <- function(donor_fac) {
    if (!isTRUE(W_DONOR_WEIGHTS_10E)) return(rep(1, length(donor_fac)))
    f <- as_donor_factor_10e(donor_fac)
    as.numeric(1 / as.numeric(table(f)[f]))
}

p10e_spine_formula <- function(k, use_re) {
    rhs <- sprintf("s(u_trep, bs = \"tp\", k = %d)", as.integer(k))
    if (use_re) rhs <- paste0(rhs, " + s(donor, bs = \"re\")")
    stats::as.formula(paste("v_resp ~", rhs))
}

# The s(u_trep) term located by label (the donor term is a second smooth).
p10e_smooth_idx <- function(model) {
    if (is.null(model$smooth) || length(model$smooth) == 0) return(NA_integer_)
    labs <- vapply(model$smooth, function(s) s$label, character(1))
    i <- which(labs == SPINE_LABEL_10E)
    if (length(i) == 1L) return(i)
    i <- which(labs != DONOR_RE_LABEL_10E)
    if (length(i) >= 1L) return(i[1])
    NA_integer_
}

# first.sp indexes model$sp; the name match and sp[1] are fallbacks.
p10e_smooth_lambda <- function(model, i) {
    if (is.na(i) || is.null(model$sp) || length(model$sp) == 0) return(NA_real_)
    sm <- model$smooth[[i]]
    if (!is.null(sm$first.sp) && sm$first.sp >= 1L &&
        sm$first.sp <= length(model$sp))
        return(unname(model$sp[sm$first.sp]))
    hit <- which(startsWith(names(model$sp), sm$label))
    if (length(hit) > 0) return(unname(model$sp[hit[1]]))
    unname(model$sp[1])
}

# PRSS = penalised residual sum of squares; the penalty comes from the spine
# term.
calc_prss_10e <- function(model) {
    RSS <- sum(stats::residuals(model)^2)
    i   <- p10e_smooth_idx(model)
    if (is.na(i))
        return(list(PRSS = RSS, RSS = RSS,
                    f_double_prime_integral = NA_real_, lambda = NA_real_))
    sm     <- model$smooth[[i]]
    lambda <- p10e_smooth_lambda(model, i)
    if (is.null(sm$S) || length(sm$S) == 0)
        return(list(PRSS = RSS, RSS = RSS,
                    f_double_prime_integral = NA_real_, lambda = lambda))
    S     <- sm$S[[1]]
    beta  <- stats::coef(model)[sm$first.para:sm$last.para]
    f_dpi <- as.numeric(t(beta) %*% S %*% beta)
    list(PRSS = RSS + lambda * f_dpi, RSS = RSS,
         f_double_prime_integral = f_dpi, lambda = lambda)
}

# Population-level curve: a reference donor level is supplied and the random
# intercept excluded, so the prediction is the same for every donor.
spine_predict_pop <- function(model, u_vals) {
    nd <- data.frame(u_trep = as.numeric(u_vals))
    has_re <- !is.null(model$smooth) &&
        any(vapply(model$smooth, function(s) s$label, character(1)) ==
                DONOR_RE_LABEL_10E)
    if (has_re) {
        ref <- model$var.summary[["donor"]]
        if (is.null(ref)) ref <- model$xlevels[["donor"]]
        ref <- if (is.factor(ref)) levels(ref)[1] else as.character(ref)[1]
        nd$donor <- factor(ref, levels = levels(model$model$donor))
        return(as.numeric(stats::predict(model, newdata = nd,
                                         exclude = DONOR_RE_LABEL_10E,
                                         type = "response")))
    }
    as.numeric(stats::predict(model, newdata = nd, type = "response"))
}

# Monotonicity check on the population-level curve.
spine_monotone_check <- function(model, n_grid = 501L, tol = 1e-8) {
    ug <- seq(0, 1, length.out = n_grid)
    gg <- spine_predict_pop(model, ug)
    dd <- diff(gg)
    i  <- which.min(dd)
    list(monotone = all(dd >= -tol),
         min_diff = min(dd),
         worst_u  = ug[i],
         n_grid   = n_grid,
         g_range  = range(gg))
}

# == The PRSS fitting loop ===================================================
# Phase 1: k = 3..max_k in 8 iterations; phase 2 refines for at most phase2_max.
fit_gam_spine_prss <- function(u_trep, v_trep, donor = NULL,
                               label      = "g(u)",
                               num_iters  = 100L,
                               seed_base  = 42L,
                               phase2_max = PRSS_PHASE2_MAX_P10E) {
    n_unique <- length(unique(u_trep))
    max_k    <- min(10L, n_unique - 1L)
    if (max_k < 3L) max_k <- 3L

    donor_f <- if (is.null(donor)) NULL else as_donor_factor_10e(donor)
    use_re  <- donor_re_ok_10e(donor_f)
    n_don   <- n_donor_levels_10e(donor_f)
    df_fit  <- data.frame(u_trep = as.numeric(u_trep),
                          v_resp = as.numeric(v_trep))
    if (!is.null(donor_f)) {
        df_fit$donor <- donor_f
        df_fit$w     <- donor_weights_10e(donor_f)
    } else {
        df_fit$w <- rep(1, nrow(df_fit))
    }
    if (!use_re && !is.null(donor_f))
        cat(sprintf(paste0("    %s: %d donor level(s) is below the minimum of ",
                           "%d, so the fit is cell level.\\n"),
                    label, n_don, DONOR_RE_MIN_LEVELS_10E))

    best_prss   <- Inf
    best_model  <- NULL
    best_k      <- 3L
    best_iter   <- 1L
    n_iter_max  <- min(as.integer(num_iters), 8L + as.integer(phase2_max))
    prss_values <- rep(NA_real_, n_iter_max)
    early_ctr   <- 0L
    early_stop  <- FALSE

    run_one <- function(i, k) {
        set.seed(seed_base + i)
        fit_at <- function(kk)
            mgcv::gam(p10e_spine_formula(kk, use_re),
                      data = df_fit, weights = w,
                      method = GAM_METHOD_10E, select = GAM_SELECT_10E,
                      gamma = GAM_GAMMA_10E)
        gm <- tryCatch(fit_at(k), error = function(e) {
            set.seed(seed_base + i)
            fit_at(max(3L, k - 1L))
        })
        list(model = gm, prss = calc_prss_10e(gm)$PRSS, k = k)
    }

    n_phase1 <- min(8L, n_iter_max)
    k_phase1 <- vapply(seq_len(n_phase1), function(i)
        min(3L + ((i - 1L) %% max(1L, max_k - 2L)), max_k), numeric(1))

    # Remaining phase-1 fits use workers only if the estimated time (from the
    # first fit, timed sequentially) is at least P10E_PAR_MIN_SECONDS.
    run_one_timed <- function(i, k) {
        t0  <- proc.time()[["elapsed"]]
        res <- tryCatch(run_one(i, k), error = function(e) NULL)
        if (!is.null(res)) res$secs <- proc.time()[["elapsed"]] - t0
        res
    }

    phase1_res <- vector("list", n_phase1)
    phase1_res[[1L]] <- run_one_timed(1L, k_phase1[1L])
    t_fit1 <- if (is.null(phase1_res[[1L]])) NA_real_ else phase1_res[[1L]]$secs
    if (!is.null(phase1_res[[1L]]))
        cat(sprintf("    %s phase1 iter 1 | k=%d | PRSS=%.5f | %.1f s\n",
                    label, phase1_res[[1L]]$k, phase1_res[[1L]]$prss, t_fit1))

    rest <- if (n_phase1 >= 2L) seq.int(2L, n_phase1) else integer(0)
    n_workers <- if (isTRUE(P10E_PARALLEL))
        max(1L, min(as.integer(P10E_N_CORES), length(rest))) else 1L
    est_rest  <- if (is.finite(t_fit1)) t_fit1 * length(rest) else 0

    used_par <- FALSE
    if (length(rest) > 0L && n_workers > 1L &&
        est_rest >= P10E_PAR_MIN_SECONDS) {
        cl <- tryCatch(parallel::makePSOCKcluster(n_workers),
                       error = function(e) NULL)
        if (!is.null(cl)) {
            on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
            par_res <- tryCatch({
                parallel::clusterEvalQ(cl, library(mgcv))
                parallel::clusterExport(
                    cl,
                    varlist = c("df_fit", "use_re", "seed_base", "k_phase1",
                                "p10e_spine_formula", "calc_prss_10e",
                                "p10e_smooth_idx", "p10e_smooth_lambda",
                                "SPINE_LABEL_10E", "DONOR_RE_LABEL_10E",
                                "GAM_METHOD_10E", "GAM_SELECT_10E",
                                "GAM_GAMMA_10E"),
                    envir = environment())
                parallel::parLapply(cl, rest, function(i) {
                    kk <- k_phase1[i]
                    t0 <- proc.time()[["elapsed"]]
                    set.seed(seed_base + i)
                    fit_at <- function(kv)
                        mgcv::gam(p10e_spine_formula(kv, use_re),
                                  data = df_fit, weights = w,
                                  method = GAM_METHOD_10E,
                                  select = GAM_SELECT_10E,
                                  gamma  = GAM_GAMMA_10E)
                    gm <- tryCatch(fit_at(kk), error = function(e) {
                        set.seed(seed_base + i)
                        tryCatch(fit_at(max(3L, kk - 1L)), error = function(e2) NULL)
                    })
                    if (is.null(gm)) return(NULL)
                    list(model = gm, prss = calc_prss_10e(gm)$PRSS, k = kk,
                         secs = proc.time()[["elapsed"]] - t0)
                })
            }, error = function(e) {
                cat(sprintf("    %s: parallel phase 1 failed (%s); sequential.\n",
                            label, conditionMessage(e)))
                NULL
            })
            if (!is.null(par_res)) {
                phase1_res[rest] <- par_res
                used_par <- TRUE
                cat(sprintf("    %s: fits 2-%d ran on %d worker(s)\n",
                            label, n_phase1, n_workers))
            }
        }
    } else if (length(rest) > 0L && n_workers > 1L) {
        cat(sprintf("    %s: sequential (one fit %.1f s, remaining %d about %.1f s, below %.1f s)\n",
                    label, t_fit1, length(rest), est_rest,
                    P10E_PAR_MIN_SECONDS))
    }

    if (!used_par && length(rest) > 0L)
        for (i in rest) phase1_res[[i]] <- run_one_timed(i, k_phase1[i])

    for (i in seq_len(n_phase1)) {
        res <- phase1_res[[i]]
        if (!is.null(res)) {
            prss_values[i] <- res$prss
            if (i > 1L)   # iteration 1 was printed as soon as it finished
                cat(sprintf("    %s phase1 iter %d | k=%d | PRSS=%.5f | %.1f s\n",
                            label, i, res$k, res$prss,
                            if (is.null(res$secs)) NA_real_ else res$secs))
            if (res$prss < best_prss) {
                best_prss <- res$prss; best_model <- res$model
                best_k <- res$k; best_iter <- i; early_ctr <- 0L
            } else early_ctr <- early_ctr + 1L
        } else early_ctr <- early_ctr + 1L
    }
    if (is.null(best_model)) {
        cat(sprintf("    %s: every phase-1 iteration failed.\n", label))
        return(list(best_model = NULL, best_k = NA_integer_,
                    best_iter = NA_integer_, best_prss = NA_real_,
                    prss_values = prss_values, actual_iters = 0L,
                    early_stopped = FALSE, max_k = max_k,
                    phase2_max = as.integer(phase2_max),
                    n_iter_max = n_iter_max,
                    donor_re = use_re, n_donors = n_don))
    }

    k_range <- unique(pmax(3L, pmin(max_k,
        if (best_k >= max_k) c(max(3L, max_k - 1L), max_k)
        else if (best_k <= 3L) c(3L, min(4L, max_k))
        else c(best_k - 1L, best_k, min(best_k + 1L, max_k)))))
    cat(sprintf("    %s phase1 best: k=%d | PRSS=%.5f\n",
                label, best_k, best_prss))

    # seq_len() avoids 9L:8L counting backwards when phase2_max = 0.
    for (i in 8L + seq_len(max(0L, n_iter_max - 8L))) {
        if (early_ctr >= PRSS_EARLY_STOP_P10E) {
            cat(sprintf("    %s early stop at iter %d\n", label, i))
            early_stop <- TRUE; break
        }
        # The phase-2 k draw is seeded from (seed_base, i), so the loop is
        # reproducible in any execution order; sample.int() handles length 1.
        set.seed(seed_base + 5000L + i)
        k   <- if (i %% 3L == 0L) k_range[sample.int(length(k_range), 1L)]
               else best_k
        res <- tryCatch(run_one(i, k), error = function(e) NULL)
        if (!is.null(res)) {
            prss_values[i] <- res$prss
            cat(sprintf("    %s refine iter %d | k=%d | PRSS=%.5f\n",
                        label, i, k, res$prss))
            if (res$prss < best_prss) {
                best_prss <- res$prss; best_model <- res$model
                best_k <- k; best_iter <- i; early_ctr <- 0L
                k_range <- unique(pmax(3L, pmin(max_k,
                    if (best_k >= max_k) c(max(3L, max_k - 1L), max_k)
                    else if (best_k <= 3L) c(3L, min(4L, max_k))
                    else c(best_k - 1L, best_k, min(best_k + 1L, max_k)))))
            } else early_ctr <- early_ctr + 1L
        } else early_ctr <- early_ctr + 1L
    }

    actual_iters <- suppressWarnings(max(which(!is.na(prss_values))))
    if (!is.finite(actual_iters)) actual_iters <- n_iter_max

    list(best_model    = best_model,
         best_k        = best_k,
         best_iter     = best_iter,
         best_prss     = best_prss,
         prss_values   = prss_values[seq_len(actual_iters)],
         actual_iters  = actual_iters,
         early_stopped = early_stop,
         max_k         = max_k,
         phase2_max    = as.integer(phase2_max),
         n_iter_max    = n_iter_max,
         donor_re      = use_re,
         n_donors      = n_don)
}

# == Weighted quantile =======================================================
# Donor-weighted quantile; with unit weights it equals the type-7 quantile.
weighted_quantile_10e <- function(x, w, probs) {
    ok <- is.finite(x) & is.finite(w) & w > 0
    x <- x[ok]; w <- w[ok]
    if (length(x) == 0L) return(NA_real_)
    if (length(unique(w)) == 1L)
        return(as.numeric(stats::quantile(x, probs = probs, type = 7,
                                          na.rm = TRUE)))
    o  <- order(x); x <- x[o]; w <- w[o]
    cw <- cumsum(w) / sum(w)
    as.numeric(stats::approx(cw, x, xout = probs, method = "linear",
                             rule = 2, ties = "ordered")$y)
}
##############################################################################
# CMRF HELPERS
##############################################################################

# == RSk-fCV folds and repeats ===============================================
# Top level so PSOCK workers can be given them by name.
p10e_stratified_folds <- function(u_set, K, rep_seed) {
    n_set   <- length(u_set)
    n_strat <- max(1L, min(10L, floor(n_set / K)))
    u_rank  <- rank(u_set, ties.method = "first")
    strat   <- pmin(ceiling(u_rank / n_set * n_strat), n_strat)
    fold_assign <- integer(n_set)
    for (s in seq_len(n_strat)) {
        idx_s <- which(strat == s)
        set.seed(rep_seed + s)
        # sample.int() avoids sample(x) treating a length-1 x as 1:x.
        idx_s_shuf <- idx_s[sample.int(length(idx_s))]
        fold_assign[idx_s_shuf] <- ((seq_along(idx_s_shuf) - 1L) %% K) + 1L
    }
    lapply(seq_len(K), function(k) fold_assign == k)
}

# One RSk-fCV repeat: k_folds RMSEs in fold order. The seed depends only on
# (seed, base_seed_offset, r), so fold k stays paired across curves.
p10e_rskfcv_repeat <- function(r, u_set, v_set, base_seed_offset, k_folds, seed) {
    rep_seed   <- seed + base_seed_offset + r * 1000L
    test_masks <- p10e_stratified_folds(u_set, k_folds, rep_seed)
    vapply(seq_len(k_folds), function(k) {
        test_mask  <- test_masks[[k]]
        train_mask <- !test_mask
        if (sum(train_mask) < 2L || sum(test_mask) < 1L) return(NA_real_)
        yhat <- iso_predict_cmrf(u_set[train_mask], v_set[train_mask],
                                 u_set[test_mask])
        sqrt(mean((v_set[test_mask] - yhat)^2))
    }, numeric(1))
}

# == iso_predict_cmrf(): isotonic (PAVA) prediction at u_new ===================
iso_predict_cmrf <- function(u_train, v_train, u_new) {
    ord   <- order(u_train)
    u_s   <- u_train[ord]
    n_tr  <- length(u_s)
    fit   <- isoreg(u_s, v_train[ord])
    iso_v <- fit$yf
    # Vectorised lookup; the clamps match the end-point branches.
    idx <- findInterval(u_new, u_s, rightmost.closed = TRUE)
    idx[idx < 1L]   <- 1L
    idx[idx > n_tr] <- n_tr
    iso_v[idx]
}

# == calibrate_cmrf_delta() ===================================================
# delta_base = donor-weighted quantile(|v - g(u)|, coverage) over all original
# SKCS-core cells.
calibrate_cmrf_delta <- function(u_cal, v_cal, g_u_cal, w_cal = NULL,
                                 coverage = 0.95) {
    cat(sprintf(paste0("  Calibrating CMRF delta_base: weighted quantile(|v ",
                       "- g(u)|, %.2f) over all original SKCS-core cells ...\n"),
                coverage))

    resid_all <- abs(v_cal - g_u_cal)
    if (is.null(w_cal)) w_cal <- rep(1, length(resid_all))

    if (length(resid_all) < 10L)
        cat("  WARNING: fewer than 10 cells for the delta calibration.\n")

    delta_base <- weighted_quantile_10e(resid_all, w_cal, coverage)

    cat(sprintf(paste0("  Calibration cells: %d | sum(w) = %.4f | %.0fth ",
                       "weighted pctile |v - g(u)| = %.6f (delta_base)\n"),
                length(resid_all), sum(w_cal), 100 * coverage, delta_base))
    cat(sprintf("  Calibration resid range: [%.6f, %.6f] | median=%.6f\n",
                min(resid_all), max(resid_all), median(resid_all)))

    list(delta_base    = delta_base,
         coverage      = coverage,
         n_grp1a       = length(u_cal),
         n_cal         = length(u_cal),
         sum_w_cal     = sum(w_cal),
         weighted      = length(unique(w_cal)) > 1L,
         resid_all     = resid_all,
         w_cal         = w_cal,
         resid_median  = median(resid_all),
         resid_max     = max(resid_all))
}

# == cmrf_stage2_fidelity_gate() ==============================================
# Stage 2: RSk-fCV on Orig+Ref alone (Curve A, once) and on the pool of
# Orig+Ref plus active candidates (Curve B, each iteration), paired by
# (repeat, fold). Removal needs Wilcoxon p < alpha and dRMSE > BV; the active
# candidate with the highest PAVA residual is removed until the gate closes.
cmrf_stage2_fidelity_gate <- function(u_grp1a, v_grp1a,
                                      u_cand, v_cand,
                                      resid_cand,
                                      n_repeats = 40L,
                                      k_folds   = 5L,
                                      alpha     = 0.05,
                                      seed      = 42L) {

    n1a        <- length(u_grp1a)
    n_cnd      <- length(u_cand)
    n_total_folds <- n_repeats * k_folds   # = 200 with defaults (40 x 5)
    cat(sprintf("  CMRF Stage 2 (B2): %d Orig+Ref anchor cells, %d Stage1 candidates\n",
                n1a, n_cnd))
    cat(sprintf("  RSk-fCV: n_repeats=%d, k_folds=%d (%d total folds) | Wilcoxon alpha=%.3f\n",
                n_repeats, k_folds, n_total_folds, alpha))
    cat("  Effect-size gate: calibrated on background variability (Formulation B)\n")

    # Top-level fold builder under its in-function name.
    make_stratified_folds <- p10e_stratified_folds

    # == Run RSk-fCV on a set of (u,v) points ================================
    # Returns n_repeats * k_folds RMSEs; fold (r, k) at (r-1)*k_folds + k.
    # Seeds per repeat keep folds paired; repeats may run on workers.
    .gate_cl <- NULL

    run_rskfcv <- function(u_set, v_set, base_seed_offset) {
        if (length(u_set) < k_folds + 1L) return(rep(NA_real_, n_total_folds))
        res <- NULL
        if (!is.null(.gate_cl))
            res <- tryCatch(
                parallel::parLapply(.gate_cl, seq_len(n_repeats),
                                    p10e_rskfcv_repeat,
                                    u_set = u_set, v_set = v_set,
                                    base_seed_offset = base_seed_offset,
                                    k_folds = k_folds, seed = seed),
                error = function(e) {
                    cat(sprintf("  RSk-fCV: parallel run failed (%s); sequential.\n",
                                conditionMessage(e)))
                    NULL
                })
        if (is.null(res))
            res <- lapply(seq_len(n_repeats), p10e_rskfcv_repeat,
                          u_set = u_set, v_set = v_set,
                          base_seed_offset = base_seed_offset,
                          k_folds = k_folds, seed = seed)
        as.numeric(unlist(res, use.names = FALSE))
    }

    # == Step 1: Orig+Ref-alone RSk-fCV baseline (timed, sequential) ========
    cat("  Step 1: Orig+Ref-alone RSk-fCV baseline ...\n")
    .t0        <- proc.time()[["elapsed"]]
    rmse_grp1a <- run_rskfcv(u_grp1a, v_grp1a, base_seed_offset = 0L)
    .t_cv_seq  <- proc.time()[["elapsed"]] - .t0
    cat(sprintf("  RSk-fCV: one %d-fold run took %.2f s sequentially\n",
                n_total_folds, .t_cv_seq))

    # == Decide on workers, then keep the cluster for the loop ===============
    .gate_workers <- if (isTRUE(P10E_PARALLEL))
        max(1L, min(as.integer(P10E_N_CORES), n_repeats)) else 1L
    if (.gate_workers > 1L && .t_cv_seq >= P10E_PAR_MIN_SECONDS) {
        .gate_cl <- tryCatch(parallel::makePSOCKcluster(.gate_workers),
                             error = function(e) NULL)
        if (!is.null(.gate_cl)) {
            .ok <- tryCatch({
                parallel::clusterExport(
                    cl = .gate_cl,
                    varlist = c("p10e_rskfcv_repeat", "p10e_stratified_folds",
                                "iso_predict_cmrf"),
                    envir = globalenv())
                TRUE
            }, error = function(e) FALSE)
            if (!isTRUE(.ok)) {
                try(parallel::stopCluster(.gate_cl), silent = TRUE)
                .gate_cl <- NULL
            }
        }
        if (is.null(.gate_cl))
            cat("  RSk-fCV: cluster unavailable; the repeats stay sequential.\n")
        else
            cat(sprintf("  RSk-fCV: repeats run on %d worker(s)\n", .gate_workers))
    } else if (.gate_workers > 1L) {
        cat(sprintf("  RSk-fCV: sequential (one run %.2f s, below %.1f s)\n",
                    .t_cv_seq, P10E_PAR_MIN_SECONDS))
    }
    on.exit(if (!is.null(.gate_cl))
                try(parallel::stopCluster(.gate_cl), silent = TRUE), add = TRUE)

    # == Background variability (BV) =========================================
    # noise_floor = IQR(RMSE_A) / median(RMSE_A)
    rmse_a_valid   <- rmse_grp1a[!is.na(rmse_grp1a)]
    median_rmse_a  <- median(rmse_a_valid)
    iqr_rmse_a     <- IQR(rmse_a_valid)
    noise_floor    <- if (median_rmse_a > 0) iqr_rmse_a / median_rmse_a else 0

    cat(sprintf("  Orig+Ref baseline RMSE: mean=%.6f  sd=%.6f  median=%.6f  IQR=%.6f\n",
                mean(rmse_a_valid), sd(rmse_a_valid), median_rmse_a, iqr_rmse_a))
    cat(sprintf("  Background variability BV (IQR/median): %.6f (%.2f%% relative threshold)\n",
                noise_floor, noise_floor * 100))

    # == PAVA removal residuals (once, on the initial union) =================
    # Removal is ranked by |v - PAVA(u)|, independent of the Stage 1 spine.
    u_init_union <- c(u_grp1a, u_cand)
    v_init_union <- c(v_grp1a, v_cand)
    pava_pred_cand <- iso_predict_cmrf(u_init_union, v_init_union, u_cand)
    resid_pava_cand <- abs(v_cand - pava_pred_cand)
    cat(sprintf("  PAVA removal residuals: mean=%.6f, max=%.6f (initial union, n=%d)\n",
                mean(resid_pava_cand), max(resid_pava_cand), length(u_init_union)))
    cat(sprintf("  Spine removal residuals: mean=%.6f, max=%.6f (for comparison)\n",
                mean(resid_cand), max(resid_cand)))

    # Edge case: no candidates
    if (n_cnd == 0L) {
        return(list(pass_stage2      = logical(0L),
                    rmse_grp1a       = rmse_grp1a,
                    rmse_final_union = rmse_grp1a,
                    wilcox_p_final   = NA_real_,
                    n_removed        = 0L,
                    removal_order    = integer(0L),
                    resid_pava_cand  = numeric(0L),
                    p_per_iter       = numeric(0L),
                    n_per_iter       = integer(0L),
                    noise_floor      = noise_floor,
                    relative_delta_final = 0,
                    effect_size_per_iter = numeric(0L),
                    n_candidates     = 0L,
                    n_repeats        = n_repeats,
                    k_folds          = k_folds,
                    n_total_folds    = n_total_folds,
                    alpha            = alpha))
    }

    # == Steps 2 to 7: iterative batch test and removal ======================
    # Remove only if Wilcoxon p < alpha and dRMSE > BV.
    active_mask    <- rep(TRUE, n_cnd)
    removal_order  <- integer(0L)
    p_per_iter     <- numeric(0L)   # p-value at each iteration for convergence plot
    n_per_iter     <- integer(0L)   # n_active at each iteration
    effect_size_per_iter <- numeric(0L)  # relative delta at each iteration
    rmse_union_cur <- NULL
    wilcox_p_cur   <- NA_real_
    relative_delta <- 0
    iter           <- 0L

    repeat {
        iter <- iter + 1L
        active_idx <- which(active_mask)
        n_active   <- length(active_idx)

        # Pool Orig+Ref + remaining active candidates
        u_pool <- c(u_grp1a, u_cand[active_idx])
        v_pool <- c(v_grp1a, v_cand[active_idx])

        # Curve B: CV on the pool, paired with Curve A by (repeat, fold) index.
        rmse_union_cur <- run_rskfcv(u_pool, v_pool,
                                     base_seed_offset = 0L)

        # Paired Wilcoxon, one-tailed: H1 = union RMSE > Orig+Ref RMSE.
        wt <- tryCatch(
            suppressWarnings(
                wilcox.test(rmse_union_cur, rmse_grp1a,
                            paired      = TRUE,
                            alternative = "greater",
                            exact       = FALSE)
            ),
            error = function(e) list(p.value = 1.0)
        )
        wilcox_p_cur <- wt$p.value

        # Relative RMSE increase, paired on folds where both curves have an RMSE.
        .pair_ok       <- !is.na(rmse_union_cur) & !is.na(rmse_grp1a)
        relative_delta <- if (any(.pair_ok)) {
            median(rmse_union_cur[.pair_ok] - rmse_grp1a[.pair_ok]) / median_rmse_a
        } else NA_real_
        if (sum(.pair_ok) < length(rmse_grp1a))
            cat(sprintf(paste0("  NOTE: %d of %d fold(s) dropped from dRMSE",
                               " (NA in Curve A or Curve B)\n"),
                        length(rmse_grp1a) - sum(.pair_ok), length(rmse_grp1a)))

        p_per_iter           <- c(p_per_iter, wilcox_p_cur)
        n_per_iter           <- c(n_per_iter, n_active)
        effect_size_per_iter <- c(effect_size_per_iter, relative_delta)

        # Removal needs both conditions; NA counts as not met.
        sig_worse     <- isTRUE(wilcox_p_cur < alpha)
        exceeds_noise <- isTRUE(relative_delta > noise_floor)
        should_remove <- sig_worse && exceeds_noise

        status_tag <- if (!sig_worse) {
            "[p >= alpha -- accepting]"
        } else if (!exceeds_noise) {
            sprintf("[p<alpha BUT dRMSE=%.4f <= BV=%.4f -- accepting (within background variability)]",
                    relative_delta, noise_floor)
        } else {
            "[WORSE & exceeds noise -- removing worst]"
        }

        cat(sprintf("  Iter %d: n_active=%d | union RMSE mean=%.6f | Wilcoxon p=%.4f | dRMSE=%.4f %s\n",
                    iter, n_active,
                    mean(rmse_union_cur, na.rm = TRUE),
                    wilcox_p_cur, relative_delta, status_tag))

        # Accept if either condition fails (union is OK)
        if (!should_remove) break

        # All candidates already removed or only 1 left -- stop
        if (n_active <= 1L) {
            last_active <- active_idx[1L]
            active_mask[last_active] <- FALSE
            removal_order <- c(removal_order, last_active)
            cat("  All candidates removed by fidelity gate.\n")
            break
        }

        # Remove the active candidate with the highest PAVA residual.
        worst_in_active <- active_idx[which.max(resid_pava_cand[active_idx])]
        active_mask[worst_in_active] <- FALSE
        removal_order <- c(removal_order, worst_in_active)
    }

    pass_stage2 <- active_mask
    n_removed   <- sum(!active_mask)

    cat(sprintf("  Stage 2 result: %d accepted, %d removed | final Wilcoxon p=%.4f | dRMSE=%.4f | BV=%.4f\n",
                sum(pass_stage2), n_removed, wilcox_p_cur, relative_delta, noise_floor))

    list(pass_stage2      = pass_stage2,
         rmse_grp1a       = rmse_grp1a,
         rmse_final_union = rmse_union_cur,
         wilcox_p_final   = wilcox_p_cur,
         n_removed        = n_removed,
         removal_order    = removal_order,
         resid_pava_cand  = resid_pava_cand,
         p_per_iter       = p_per_iter,
         n_per_iter       = n_per_iter,
         noise_floor      = noise_floor,
         relative_delta_final = relative_delta,
         effect_size_per_iter = effect_size_per_iter,
         n_candidates     = n_cnd,
         n_repeats        = n_repeats,
         k_folds          = k_folds,
         n_total_folds    = n_total_folds,
         alpha            = alpha)
}

# == apply_cmrf() =============================================================
# Stage 1: |v - g(u)| <= delta_adapt(u) = delta_base * (1 + kappa * K(u, IP_u, bw)).
# Stage 2: fidelity gate on the non-SKCS Stage 1 candidates; Orig+Ref cells are
# the reference and are not re-filtered.
apply_cmrf <- function(u_all, v_all, g_u_all, is_orig_skcs, IP_u,
                       delta_base, kappa = 0.50, bw = 0.10,
                       n_repeats = 40L, k_folds = 5L, alpha = 0.05,
                       seed = 42L) {

    # == Adaptive delta_adapt(u_i) ===========================================
    K_weight    <- as.numeric(exp(-(u_all - IP_u)^2 / (2 * bw^2)))
    delta_adapt <- as.numeric(delta_base * (1 + kappa * K_weight))

    # == Stage 1: spine band proximity =======================================
    resid_abs <- as.numeric(abs(v_all - g_u_all))
    cond2     <- resid_abs <= delta_adapt

    # == Stage 1 pass (band proximity for all cells) ==========================
    pass_stage1 <- cond2

    # == Stage 1 SKCS anchor: spine-proximal SKCS (band pass) ================
    grp1a_mask <- is_orig_skcs & pass_stage1   # = "Original + Refined"
    u_grp1a    <- u_all[grp1a_mask]
    v_grp1a    <- v_all[grp1a_mask]
    cat(sprintf("  Stage 1 SKCS anchor: %d (band pass) of %d total SKCS\n",
                sum(grp1a_mask), sum(is_orig_skcs)))
    yhat_iso <- iso_predict_cmrf(u_grp1a, v_grp1a, u_all)

    # == Stage 2: fidelity gate on non-SKCS Stage 1 candidates only ===========
    # Orig+Ref cells passing Stage 1 are accepted directly (reference set).
    non_skcs_stage1_idx <- which(!is_orig_skcs & pass_stage1)
    stage2_res <- NULL

    if (length(non_skcs_stage1_idx) > 0L) {
        stage2_res <- cmrf_stage2_fidelity_gate(
            u_grp1a    = u_grp1a,
            v_grp1a    = v_grp1a,
            u_cand     = u_all[non_skcs_stage1_idx],
            v_cand     = v_all[non_skcs_stage1_idx],
            resid_cand = resid_abs[non_skcs_stage1_idx],
            n_repeats  = n_repeats,
            k_folds    = k_folds,
            alpha      = alpha,
            seed       = seed
        )
    }

    # stage2_res vectors are positional; record the rows (of u_all) they refer to.
    if (!is.null(stage2_res)) {
        stage2_res$anchor_row <- which(grp1a_mask)
        stage2_res$cand_row   <- non_skcs_stage1_idx
    }

    # == Assemble per-cell Stage 2 result ====================================
    pass_stage2_all <- logical(length(u_all))
    # Orig+Ref cells that passed Stage 1 always pass Stage 2 (reference set)
    pass_stage2_all[is_orig_skcs & pass_stage1] <- TRUE
    # Non-SKCS candidates: use RSk-fCV gate result
    if (!is.null(stage2_res) && length(non_skcs_stage1_idx) > 0L) {
        pass_stage2_all[non_skcs_stage1_idx] <- stage2_res$pass_stage2
    }

    # == Final CMRF pass = Stage 1 AND Stage 2 ===============================
    pass_cmrf <- pass_stage1 & pass_stage2_all

    # == Classification ======================================================
    cell_class <- dplyr::case_when(
        is_orig_skcs &  pass_cmrf ~ "Original + Refined",
        is_orig_skcs & !pass_cmrf ~ "Original only",
        !is_orig_skcs &  pass_cmrf ~ "Recovered",
        TRUE                       ~ "Excluded"
    )

    is_recovered <- !is_orig_skcs & pass_cmrf

    # == Stage 2 per-cell diagnostics columns =================================
    # Candidates removed by the iterative gate.
    removed_by_stage2_all <- rep(FALSE, length(u_all))
    if (!is.null(stage2_res) && length(non_skcs_stage1_idx) > 0L) {
        removed_idx <- non_skcs_stage1_idx[!stage2_res$pass_stage2]
        removed_by_stage2_all[removed_idx] <- TRUE
    }

    list(g_u_all              = g_u_all,
         resid_abs             = resid_abs,
         delta_adapt           = delta_adapt,
         K_weight              = K_weight,
         yhat_iso              = yhat_iso,
         cond2                 = cond2,
         pass_stage1           = pass_stage1,
         pass_stage2           = pass_stage2_all,
         pass_cmrf             = pass_cmrf,
         cell_class            = cell_class,
         is_recovered          = is_recovered,
         stage2_res            = stage2_res,
         removed_by_stage2     = removed_by_stage2_all)
}


##############################################################################
# SW-GATED STATISTICAL TEST HELPERS
# Used by run_part10e() (Plot 2) and the verification section.
##############################################################################

# == SW helper (subsample to 5000 if n > 5000) ===============================
sw_pval <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) < 3)   return(NA_real_)
    if (length(x) > 5000) x <- sample(x, 5000, replace = FALSE)
    tryCatch(shapiro.test(x)$p.value, error = function(e) NA_real_)
}

# == 2-group SW-gated test ===================================================
two_grp <- function(x1, x2, alternative = "two.sided", paired = FALSE) {
    x1 <- x1[!is.na(x1)]; x2 <- x2[!is.na(x2)]
    if (length(x1) < 2 || length(x2) < 2)
        return(list(test_used="n<2", stat=NA_real_, p=NA_real_,
                    sw_p1=NA_real_, sw_p2=NA_real_, n1=NA_real_,
                    norm1=NA, norm2=NA))
    sw1 <- sw_pval(x1); sw2 <- sw_pval(x2)
    n1  <- !is.na(sw1) && sw1 >= 0.05
    n2  <- !is.na(sw2) && sw2 >= 0.05
    if (isTRUE(n1) && isTRUE(n2)) {
        res <- tryCatch(
            if (paired) t.test(x1, x2, paired=TRUE,  alternative=alternative) else        t.test(x1, x2, paired=FALSE, alternative=alternative),
            error = function(e) NULL)
        tnm <- if (paired) sprintf("Paired t-test (%s)", alternative) else        sprintf("Student t-test (%s)", alternative)
    } else {
        res <- tryCatch(suppressWarnings(
            wilcox.test(x1, x2, paired=paired, alternative=alternative,
                        exact=FALSE)),
            error = function(e) NULL)
        tnm <- if (paired) sprintf("Wilcoxon signed-rank paired (%s)", alternative) else        sprintf("Mann-Whitney U (%s)", alternative)
    }
    list(test_used = tnm,
         stat      = if (is.null(res)) NA_real_ else round(as.numeric(res$statistic), 4),
         p         = if (is.null(res)) NA_real_ else res$p.value,
         sw_p1     = round(sw1, 6), sw_p2 = round(sw2, 6),
         norm1     = n1,            norm2  = n2)
}

# == Dunn's post-hoc test =====================================================
# Pooled mid-ranks, tie-corrected z, two-sided p, BH; `levs` fixes the order.
dunn_bh <- function(v, g, levs = NULL) {
    keep <- !is.na(v) & !is.na(g)
    v <- as.numeric(v[keep]); g <- as.character(g[keep])
    gs <- if (!is.null(levs)) levs[levs %in% unique(g)] else sort(unique(g))
    if (length(gs) < 2) return(NULL)
    N    <- length(v)
    R    <- rank(v)                       # average ranks for ties
    tb   <- as.numeric(table(v))
    Tcor <- sum(tb^3 - tb)                # tie-correction term
    sigma <- (N * (N + 1) / 12) - (Tcor / (12 * (N - 1)))
    Rbar <- vapply(gs, function(k) mean(R[g == k]), numeric(1))
    nn   <- vapply(gs, function(k) sum(g == k),     numeric(1))
    cmp <- character(0); Z <- numeric(0); P <- numeric(0)
    for (i in seq_len(length(gs) - 1L)) for (j in (i + 1L):length(gs)) {
        a <- gs[i]; b <- gs[j]
        z <- (Rbar[a] - Rbar[b]) / sqrt(sigma * (1 / nn[a] + 1 / nn[b]))
        cmp <- c(cmp, paste(a, "-", b)); Z <- c(Z, z); P <- c(P, 2 * pnorm(-abs(z)))
    }
    p_bh <- p.adjust(P, method = "BH")
    data.frame(comparison = cmp, Z = round(Z, 4),
               p_raw = signif(P, 4), p_BH = signif(p_bh, 4),
               stringsAsFactors = FALSE)
}

# == Multi-group SW-gated test (4 classes) ====================================
multi_grp <- function(vals, grps, levs = NULL) {
    df <- data.frame(v = vals, g = as.character(grps), stringsAsFactors=FALSE)
    df <- df[!is.na(df$v) & !is.na(df$g), ]
    if (!is.null(levs)) df <- df[df$g %in% levs, ]
    gsz <- table(df$g)
    if (length(gsz) < 2 || any(gsz < 2))
        return(list(test_used="insufficient groups", stat=NA_real_, p=NA_real_,
                    sw_note="n/a", all_norm=NA, posthoc=NULL, ph_method="n/a"))
    sw_ps   <- sapply(names(gsz), function(g) sw_pval(df$v[df$g == g]))
    norms   <- !is.na(sw_ps) & sw_ps >= 0.05
    all_n   <- all(norms)
    sw_note <- paste(names(gsz), sprintf("SW_p=%s(%s)",
                                         ifelse(is.na(sw_ps), "NA", sprintf("%.4f", sw_ps)),
                                         ifelse(norms, "N", "non-N")), collapse="; ")
    if (all_n) {
        aov_r <- tryCatch(summary(aov(v ~ g, data=df)), error=function(e) NULL)
        overall_p    <- if (!is.null(aov_r)) aov_r[[1]][["Pr(>F)"]][1] else NA_real_
        overall_stat <- if (!is.null(aov_r)) aov_r[[1]][["F value"]][1] else NA_real_
        ph_r <- tryCatch(
            pairwise.t.test(df$v, df$g, p.adjust.method="holm"),
            error=function(e) NULL)
        ph_df <- if (!is.null(ph_r)) {
            pm <- ph_r$p.value
            idx <- which(!is.na(pm), arr.ind=TRUE)
            data.frame(group1=rownames(pm)[idx[,1]],
                       group2=colnames(pm)[idx[,2]],
                       p_holm_sidak=round(as.numeric(pm[idx]),6),
                       stringsAsFactors=FALSE)
        } else NULL
        list(test_used="One-way ANOVA", stat=round(overall_stat,4),
             p=round(overall_p,6), sw_note=sw_note, all_norm=TRUE,
             posthoc=ph_df, ph_method="Holm-Sidak (pairwise t-test)")
    } else {
        kw_r <- tryCatch(kruskal.test(v ~ g, data=df), error=function(e) NULL)
        overall_p    <- if (!is.null(kw_r)) kw_r$p.value   else NA_real_
        overall_stat <- if (!is.null(kw_r)) kw_r$statistic else NA_real_
        ph_df <- tryCatch(dunn_bh(df$v, df$g, levs = levs),
                          error = function(e) NULL)
        list(test_used="Kruskal-Wallis", stat=round(overall_stat,4),
             p=round(overall_p,6), sw_note=sw_note, all_norm=FALSE,
             posthoc=ph_df, ph_method="Dunn's post-hoc (BH-adjusted)")
    }
}


# == fmt_p(): compact p formatter, underflow-safe ==============================
#   fmt_p(4.51e-8) -> "4.51e-08";  fmt_p(0) -> "< 1e-300"
fmt_p <- function(p) {
    if (length(p) == 0 || is.na(p)) return("NA")
    if (p <= 0 || p < 1e-300)       return("< 1e-300")
    formatC(p, format = "e", digits = 2)
}

# Pull a BH-adjusted Dunn p-value for an unordered class pair from a posthoc
# table whose `comparison` column is "<A> - <B>" (dunn.test convention).
get_dunn_bh <- function(ph, a, b) {
    if (is.null(ph)) return(NA_real_)
    if ("comparison" %in% names(ph)) {              # Dunn BH path
        parts <- strsplit(as.character(ph$comparison), " - ", fixed = TRUE)
        idx <- which(vapply(parts, function(pp) setequal(trimws(pp), c(a, b)), logical(1)))
        if (length(idx) >= 1) return(as.numeric(ph$p_BH[idx[1]]))
    } else if (all(c("group1", "group2") %in% names(ph))) {   # ANOVA Holm-Sidak path
        idx <- which((ph$group1 == a & ph$group2 == b) |
                     (ph$group1 == b & ph$group2 == a))
        if (length(idx) >= 1) return(as.numeric(ph$p_holm_sidak[idx[1]]))
    }
    NA_real_
}


# =============================================================================
# SECTION 0e: RDS size control
# =============================================================================
# Fitted-model environments are re-pointed at baseenv() (lossless) and plot
# objects, already written to disk, are dropped. RDS_KEEP_PLOTS = TRUE keeps them.
RDS_KEEP_PLOTS <- FALSE

# Handles gam, glm and lm.
rds_strip_gam_env <- function(m) {
    if (!inherits(m, c("gam", "glm", "lm"))) return(m)
    if (inherits(m$formula, "formula"))      environment(m$formula)      <- baseenv()
    if (inherits(m$pred.formula, "formula")) environment(m$pred.formula) <- baseenv()
    if (!is.null(m$terms))  attr(m$terms,  ".Environment") <- baseenv()
    if (!is.null(m$pterms)) attr(m$pterms, ".Environment") <- baseenv()
    if (!is.null(m$model) && !is.null(attr(m$model, "terms")))
        attr(attr(m$model, "terms"), ".Environment") <- baseenv()
    for (i in seq_along(m$smooth)) {
        if (inherits(m$smooth[[i]]$form, "formula"))
            environment(m$smooth[[i]]$form) <- baseenv()
        if (!is.null(m$smooth[[i]]$terms))
            attr(m$smooth[[i]]$terms, ".Environment") <- baseenv()
    }
    m
}

.rds_is_plot <- function(x)
    inherits(x, c("ggplot", "gg", "ggproto", "plotly", "htmlwidget",
                  "gtable", "grob", "gTree", "gList", "trellis",
                  "recordedplot"))

# Strips every fitted model and, unless RDS_KEEP_PLOTS, replaces plots with
# NULL; names, classes and shape are unchanged.
rds_slim <- function(x, keep_plots = RDS_KEEP_PLOTS) {
    if (inherits(x, c("gam", "glm", "lm"))) return(rds_strip_gam_env(x))
    if (!isTRUE(keep_plots) && .rds_is_plot(x)) return(NULL)
    if (is.list(x) && !is.data.frame(x)) {
        att <- attributes(x)
        out <- lapply(x, rds_slim, keep_plots = keep_plots)
        attributes(out) <- att
        return(out)
    }
    x
}

# Saves through rds_slim() and reports the file size.
rds_save <- function(obj, path, label = "RDS") {
    saveRDS(rds_slim(obj), path)
    cat(sprintf("  %s saved: %s  (%.1f MB%s)\n", label, basename(path),
                file.size(path) / 1e6,
                if (isTRUE(RDS_KEEP_PLOTS)) ", plots kept" else ", plots dropped"))
    invisible(path)
}


# =============================================================================
# SECTION 0b: The pipeline for one arm
# =============================================================================

run_one_arm <- function(arm_tag, config_file, out_dir, file_tag) {

    config_path <- file.path(CEP_HOME, config_file)
    if (!file.exists(config_path))
        stop("Configuration not found: ", config_path, call. = FALSE)

    .p10e_preamble <- character(0)

    assign("CEP_ARM", arm_tag, envir = globalenv())
    on.exit(suppressWarnings(rm("CEP_ARM", envir = globalenv())), add = TRUE)

    .p10e_preamble <- tryCatch(
        utils::capture.output({
            cat("\n", strrep("=", 76), "\n", sep = "")
            cat("  ARM: ", arm_tag, "   config: ", config_file, "\n", sep = "")
            cat(strrep("=", 76), "\n", sep = "")
            source(config_path)
        }),
        error = function(e) {
            cat("\n", strrep("=", 76), "\n", sep = "")
            cat("  ARM: ", arm_tag, "   config: ", config_file, "\n", sep = "")
            cat(strrep("=", 76), "\n", sep = "")
            stop(e)
        })

    .p10e_needed <- c("ARM_TAG", "ARM_DIR", "R_EXPORT_ROOT", "GOI", "DFG",
                      "HKG_GENES", "DFG_COMPOSITE_LABEL", "RUN_SUFFIX",
                      "PRIMARY_STYPE", "NORM_TAG", "NORM_TARGET_SUM",
                      "NORM_CEILING")
    .p10e_absent <- .p10e_needed[!vapply(.p10e_needed, exists, logical(1))]
    if (length(.p10e_absent) > 0)
        stop("The configuration loaded from\n    ", config_path,
             "\n  does not define: ", paste(.p10e_absent, collapse = ", "),
             "\n  This is the signature of a pre-migration config.",
             call. = FALSE)

    if (!identical(ARM_TAG, arm_tag))
        stop("The configuration reports ARM_TAG '", ARM_TAG,
             "' but this arm is '", arm_tag, "'.", call. = FALSE)

    .p10e_preamble <- c(.p10e_preamble, utils::capture.output(
        cat("  Config OK: ", ARM_TAG, " | ", GOI, " + ", length(DFG), " DFG | ",
            NORM_TAG, "\n", sep = "")))

    set.seed(42)

    DATASETS <- c("NonNeoC", "NeoC")

    # == Input paths =========================================================
    PART10_DIR      <- cep_arm_dir_for("Part10", arm_tag)
    P10_FILE_PREFIX <- paste0(cep_out_prefix("Part10", NORM_TAG,
                                             cep_file_tag_for(arm_tag)),
                              "_Module1_", gsub("-", "", GOI), "_",
                              DFG_COMPOSITE_LABEL)
    p10e_scores_candidates <- function(ds) {
        fn <- paste0(P10_FILE_PREFIX, "_", ds, "_Pooled_SKCS_scores.csv")
        c(file.path(PART10_DIR, paste0(basename(PART10_DIR), "_", ds), fn),
          file.path(PART10_DIR, fn))
    }
    p10e_scores_csv <- function(ds) {
        cand <- p10e_scores_candidates(ds)
        hit  <- cand[file.exists(cand)]
        if (length(hit) > 0L) hit[1] else cand[1]
    }

    P10B_SKCS_PREFIX <- paste0(cep_out_prefix("Part10b", NORM_TAG,
                                              cep_file_tag_for(arm_tag)),
                               "_SKCS")
    p10e_ep_candidates <- function(ds) {
        fn <- paste0(P10B_SKCS_PREFIX, "_", ds, "_EP_TREP_values.csv")
        c(file.path(cep_arm_dir_for("Part10b_SKCS", arm_tag), fn),
          file.path(cep_arm_dir_for("Part10b", arm_tag), fn),
          file.path(R_EXPORT_ROOT, "Part10b_SKCS",
                    paste0("Part10b_scran_", ds, "_EP_TREP_values.csv")))
    }
    p10e_ep_csv <- function(ds) {
        cand <- p10e_ep_candidates(ds)
        hit  <- cand[file.exists(cand)]
        if (length(hit) > 0L) hit[1] else cand[1]
    }

    p10e_p10d_candidates <- function() {
        c(file.path(cep_arm_dir_for("Part10d_Copula2", arm_tag),
                    paste0(cep_out_prefix("Part10d", NORM_TAG,
                                          cep_file_tag_for(arm_tag)),
                           "_Copula2_all_results.rds")),
          file.path(R_EXPORT_ROOT, "Part10d_Copula2",
                    "Part10d_scran_Copula2_all_results.rds"))
    }
    p10e_p10d_rds <- function() {
        cand <- p10e_p10d_candidates()
        hit  <- cand[file.exists(cand)]
        if (length(hit) > 0L) hit[1] else cand[1]
    }

    # == Donor lookup from the per-region gam_data CSVs ========================
    # Built once per dataset and cached.
    .p10e_donor_cache <- list()
    p10e_donor_map <- function(ds) {
        if (!is.null(.p10e_donor_cache[[ds]])) return(.p10e_donor_cache[[ds]])
        gam_dir <- if (identical(ds, "NeoC")) NEOC_DIR else NONNC_DIR
        regs    <- if (identical(ds, "NeoC")) NEOC_REGIONS else NONNC_REGIONS
        stype   <- gsub("-", "_", gsub("/", "_", PRIMARY_STYPE))
        bc <- character(0); don <- character(0)
        for (reg in regs) {
            gp <- file.path(gam_dir, "gam_data",
                            paste0(reg, "_", stype, "_selected.csv"))
            if (!file.exists(gp)) next
            df <- read.csv(gp, stringsAsFactors = FALSE, check.names = FALSE)
            if (!("donor" %in% colnames(df)))
                stop(sprintf(paste0("Region %s: gam_data has no donor column, ",
                                    "so the spine cannot carry the donor ",
                                    "random intercept."), reg), call. = FALSE)
            if (!("barcode" %in% colnames(df)))
                stop(sprintf("Region %s: gam_data has no barcode column.", reg),
                     call. = FALSE)
            bc  <- c(bc,  df$barcode)
            don <- c(don, as.character(df$donor))
        }
        if (length(bc) == 0L)
            stop("No gam_data CSVs found for ", ds, " under ",
                 file.path(gam_dir, "gam_data"), call. = FALSE)
        dup <- duplicated(bc)
        if (any(dup)) {
            chk <- tapply(don, bc, function(v) length(unique(v)))
            if (any(chk > 1L))
                stop(sprintf(paste0("%s: %d barcode(s) carry two different ",
                                    "donor labels across region files."),
                             ds, sum(chk > 1L)), call. = FALSE)
            bc  <- bc[!dup]; don <- don[!dup]
        }
        m <- setNames(don, bc)
        .p10e_donor_cache[[ds]] <<- m
        m
    }

    # == Output directory and file prefix ====================================
    OUT_DIR <- out_dir
    dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
    OUT_PREFIX  <- cep_out_prefix(PART_TAG, NORM_TAG, cep_file_tag_for(arm_tag))
    FILE_PREFIX <- paste0(OUT_PREFIX, "_CMRF")

    # == Log file ============================================================
    LOG_MESSAGES_P10E <- TRUE
    LOG_PATH_P10E <- file.path(OUT_DIR,
                               paste0(FILE_PREFIX, "_console_log_",
                                      format(Sys.time(), "%Y%m%d_%H%M%S"),
                                      ".txt"))
    MSG_PATH_P10E <- file.path(OUT_DIR,
                               paste0(FILE_PREFIX, "_console_log_messages.tmp"))

    while (sink.number() > 0) sink()
    if (sink.number(type = "message") != 2L) sink(type = "message")

    .p10e_log_con <- file(LOG_PATH_P10E, open = "wt")
    .p10e_msg_con <- if (isTRUE(LOG_MESSAGES_P10E))
        file(MSG_PATH_P10E, open = "wt") else NULL
    sink(.p10e_log_con, append = FALSE, split = TRUE)
    if (!is.null(.p10e_msg_con)) sink(.p10e_msg_con, type = "message")
    on.exit({
        tryCatch({
            if (!isTRUE(.p10e_completed))
                cat("\n", strrep("!", 72),
                    "\n  Arm did not complete. See \"*** ARM FAILED\" in the console.\n",
                    strrep("!", 72), "\n", sep = "")
            if (sink.number(type = "message") != 2L) sink(type = "message")
            while (sink.number() > 0) sink()
            if (!is.null(.p10e_msg_con)) close(.p10e_msg_con)
            close(.p10e_log_con)
            if (!is.null(.p10e_msg_con) && file.exists(MSG_PATH_P10E)) {
                .msg_txt <- readLines(MSG_PATH_P10E, warn = FALSE)
                if (length(.msg_txt) > 0)
                    cat(c("", strrep("=", 72),
                          "  MESSAGE STREAM (warnings, messages, errors)",
                          strrep("=", 72), .msg_txt, ""),
                        sep = "\n", file = LOG_PATH_P10E, append = TRUE)
                unlink(MSG_PATH_P10E)
            }
        }, error = function(e) invisible(NULL))
    }, add = TRUE)

    .p10e_completed <- FALSE
    .p10e_warn_opt  <- getOption("warn")
    options(warn = 1)
    on.exit(options(warn = .p10e_warn_opt), add = TRUE)
    cat(sprintf("Part10e CMRF - log started: %s\n\n", Sys.time()))
    if (length(.p10e_preamble) > 0)
        cat(paste(.p10e_preamble, collapse = "\n"), "\n\n", sep = "")

    # == IPRS-derived IP values ==============================================
    .iprs_candidates <- c(
        file.path(cep_arm_dir_for("Part10b_SKCS", arm_tag),
                  paste0(cep_out_prefix("Part10b", NORM_TAG,
                                        cep_file_tag_for(arm_tag)),
                         "_IPRS_results.rds")),
        file.path(cep_arm_dir_for("Part10b_IPRS", arm_tag),
                  paste0(cep_out_prefix("Part10b_IPRS", NORM_TAG,
                                        cep_file_tag_for(arm_tag)),
                         "_results.rds")),
        file.path(cep_arm_dir_for("Part10b", arm_tag),
                  paste0(cep_out_prefix("Part10b", NORM_TAG,
                                        cep_file_tag_for(arm_tag)),
                         "_IPRS_results.rds")),
        file.path(R_EXPORT_ROOT, "Part10b_IPRS",
                  "Part10b_scran_IPRS_results.rds"))
    .iprs_hit     <- .iprs_candidates[file.exists(.iprs_candidates)]
    IPRS_RDS_PATH <- if (length(.iprs_hit) > 0L) .iprs_hit[1]
                     else .iprs_candidates[1]

    IPRS_IP_DFG   <- NULL
    IPRS_IP_DFG_Y <- NULL
    ip_from       <- "no IPRS RDS read"
    if (file.exists(IPRS_RDS_PATH)) {
        iprs_rds    <- readRDS(IPRS_RDS_PATH)
        IPRS_IP_DFG <- iprs_rds$inflection_points_DFG
        # IP_y = f(IP_x), read from the IPRS RDS.
        IPRS_IP_DFG_Y <- iprs_rds$inflection_points_DFG_y
        ip_from     <- if (!is.null(iprs_rds$ip_definition))
            iprs_rds$ip_definition else "unrecorded definition"
        cat("  [IPRS] ", basename(IPRS_RDS_PATH), " | definition: ", ip_from,
            "\n", sep = "")
        if (!is.null(IPRS_IP_DFG))
            for (ds in names(IPRS_IP_DFG)) {
                .y_ds <- if (!is.null(IPRS_IP_DFG_Y) &&
                             ds %in% names(IPRS_IP_DFG_Y))
                    suppressWarnings(as.numeric(IPRS_IP_DFG_Y[ds])) else NA_real_
                cat(sprintf("  [IPRS] %-10s IP_x = %.6f | IP_y = %s\n",
                            ds, IPRS_IP_DFG[ds],
                            if (is.finite(.y_ds)) sprintf("%.6f", .y_ds)
                            else "ABSENT (Step 4/5 will be skipped)"))
            }
    } else {
        cat("  [IPRS] RDS not found. Looked in:\n")
        for (.p in .iprs_candidates) cat("      ", .p, "\n", sep = "")
        cat("  [IPRS] No IP and no IP_y: Steps 1-3 run, Step 4/5 skipped.\n")
    }

    cat("\n", strrep("=", 72), "\n")
    cat("  Part10e: Si x Ki copula and CMRF\n")
    cat("  Arm          :", ARM_TAG, "\n")
    cat("  GOI          :", GOI, "\n")
    cat("  DFG          :", DFG_COMPOSITE_LABEL, " - ",
        paste(DFG, collapse = ", "), "\n")
    cat("  Normalisation:", NORM_TAG, " log1p(UMIs/total x ",
        format(NORM_TARGET_SUM, scientific = FALSE), ")\n")
    cat("  Spine        : mgcv s(u) +",
        if (GAM_DONOR_RE_10E) "s(donor, bs = \"re\")" else "(no donor RE)",
        "|", if (W_DONOR_WEIGHTS_10E) "inverse-donor weights" else "uniform weights",
        "| monotonicity checked\n")
    cat("  PRSS         : phase 1 = 8 iterations, phase 2 = at most",
        PRSS_PHASE2_MAX_P10E, "\n")
    cat("  Output dir   :", OUT_DIR, "\n")
    cat(strrep("=", 72), "\n\n")

    # Marker removed when the arm completes.
    INCOMPLETE_MARKER_P10E <- file.path(OUT_DIR, "_RUN_INCOMPLETE.txt")
    writeLines(c(
        "This arm did not finish.",
        "",
        paste("Arm     :", arm_tag),
        paste("Part    :", PART_DIR_TAG),
        paste("Started :", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
        "",
        "Any files in this folder are from a partial run and must not be used.",
        "This file is removed automatically when the arm completes."
    ), INCOMPLETE_MARKER_P10E)


    ##########################################################################
    # SECTION A0: PRE-ANALYSIS GUARDRAILS
    ##########################################################################
    #   GR-01  norm.mtx max within the CP10K ceiling
    #   GR-02  expm1(norm) row sums pin NORM_TARGET_SUM
    #   GR-04  raw count matrix integer-like
    #   GR-05  Part 10b EP/TREP CSV exists
    #   GR-06  Part 10d RDS exists
    #   GR-07  IPRS IP values present

    GR_LOG_10E <- file.path(OUT_DIR, paste0(FILE_PREFIX, "_Guardrail_Log.txt"))
    .gr10e_log  <- list()
    .gr10e_pass <- 0L
    .gr10e_fail <- 0L
    .gr10e_warn <- 0L

    gr10e <- function(check_id, status, detail) {
        # A zero-length detail would vanish from the log.
        if (length(detail) != 1L || is.na(detail))
            detail <- paste0("[detail unavailable: the message was built from a",
                             " missing value and came out zero-length]")
        # An unrecognised status keeps its row.
        tag <- switch(status, PASS = "[PASS]", FAIL = "[FAIL]", WARN = "[WARN]",
                      sprintf("[%s]", as.character(status)[1]))
        cat(sprintf("  %s  %-16s  %s\n", tag, check_id, detail))
        .gr10e_log[[length(.gr10e_log) + 1L]] <<- list(
            id = check_id, status = status, detail = detail)
        if (status == "PASS") .gr10e_pass <<- .gr10e_pass + 1L
        if (status == "FAIL") .gr10e_fail <<- .gr10e_fail + 1L
        if (status == "WARN") .gr10e_warn <<- .gr10e_warn + 1L
    }

    gr10e_write <- function() {
        con <- file(GR_LOG_10E, "wt")
        # Close the connection on error.
        on.exit(try(close(con), silent = TRUE), add = TRUE)
        writeLines(c(
            strrep("=", 72),
            paste0("Part 10e CMRF - Guardrail Log [", NORM_TAG, "]"),
            paste("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
            paste("Arm :", ARM_TAG),
            paste0("Normalisation: adata.X = log1p(UMIs/total x ",
                   format(NORM_TARGET_SUM, scientific = FALSE),
                   "), exported unchanged"),
            strrep("=", 72), ""
        ), con)
        for (item in .gr10e_log)
            writeLines(sprintf("  [%s]  %-16s  %s", item$status, item$id,
                               item$detail), con)
        writeLines(c("", strrep("=", 72),
                     sprintf("  Total : %d",
                             .gr10e_pass + .gr10e_fail + .gr10e_warn),
                     sprintf("  PASS  : %d", .gr10e_pass),
                     sprintf("  WARN  : %d", .gr10e_warn),
                     sprintf("  FAIL  : %d", .gr10e_fail),
                     if (.gr10e_fail == 0 && .gr10e_warn == 0)
                         "\n  Result: no failures, no warnings"
                     else if (.gr10e_fail == 0)
                         "\n  Result: no failures; review the warnings above"
                     else "\n  Result: failures present",
                     strrep("=", 72)), con)
        close(con)
        cat("  Guardrail log written:", GR_LOG_10E, "\n")
    }

    cat("\n", strrep("=", 72), "\n")
    cat("  SECTION A0: PRE-ANALYSIS GUARDRAILS\n")
    cat(strrep("=", 72), "\n\n")

    for (.ds in list(
        list(tag = "NeoC",    dir = NEOC_DIR,  dfg = DFG_NEOC),
        list(tag = "NonNeoC", dir = NONNC_DIR, dfg = DFG_NONNC)
    )) {
        cat("  -- Dataset:", .ds$tag, "--------------------------------------\n")
        norm_mtx_path <- file.path(.ds$dir, "deg_support", "norm.mtx")
        raw_mtx_path  <- file.path(.ds$dir, "counts", "raw_selected",
                                   "matrix.mtx")

        if (file.exists(norm_mtx_path)) {
            X_n   <- as(readMM(norm_mtx_path), "CsparseMatrix")
            n_max <- max(X_n)
            gr10e(paste0("GR-01-", .ds$tag),
                  if (n_max > 0 && n_max <= NORM_CEILING + 1e-6) "PASS" else "FAIL",
                  sprintf("norm.mtx max = %.4f  (ceiling log1p(%s) = %.4f)",
                          n_max, format(NORM_TARGET_SUM, scientific = FALSE),
                          NORM_CEILING))
            n_rows <- min(50L, nrow(X_n))
            rs     <- rowSums(expm1(X_n[seq_len(n_rows), , drop = FALSE]))
            rs_dev <- max(abs(rs - NORM_TARGET_SUM))
            gr10e(paste0("GR-02-", .ds$tag),
                  if (rs_dev < NORM_TARGET_SUM * 1e-6) "PASS" else "FAIL",
                  sprintf("expm1(norm) row sums: mean = %.6f, max |dev| = %.2e, n = %d",
                          mean(rs), rs_dev, n_rows))
            rm(X_n); gc(verbose = FALSE)
        } else {
            gr10e(paste0("GR-01-", .ds$tag), "FAIL",
                  paste("norm.mtx not found:", norm_mtx_path))
            gr10e(paste0("GR-02-", .ds$tag), "FAIL",
                  paste("norm.mtx not found, row-sum identity not checked:",
                        norm_mtx_path))
        }

        if (file.exists(raw_mtx_path)) {
            X_r    <- as(readMM(raw_mtx_path), "CsparseMatrix")
            r_samp <- as.numeric(X_r@x[seq_len(min(10000, length(X_r@x)))])
            is_int <- all(abs(r_samp - round(r_samp)) < 0.01)
            r_max  <- max(X_r)
            gr10e(paste0("GR-04-", .ds$tag),
                  if (is_int && r_max > 20) "PASS" else "FAIL",
                  sprintf("raw max = %.0f, integer-like = %s",
                          r_max, if (is_int) "TRUE" else "FALSE"))
            rm(X_r); gc(verbose = FALSE)
        } else {
            gr10e(paste0("GR-04-", .ds$tag), "FAIL",
                  paste("raw_selected matrix not found:", raw_mtx_path))
        }
    }

    for (ds in DATASETS) {
        .p <- p10e_ep_csv(ds)
        gr10e(paste0("GR-05-", ds),
              if (file.exists(.p)) "PASS" else "FAIL",
              sprintf("Part 10b EP/TREP CSV: %s",
                      if (file.exists(.p))
                          file.path(basename(dirname(.p)), basename(.p))
                      else paste("not found in any location:",
                                 paste(p10e_ep_candidates(ds), collapse = "  |  "))))
    }

    .p10d <- p10e_p10d_rds()
    gr10e("GR-06-Part10d",
          if (file.exists(.p10d)) "PASS" else "FAIL",
          sprintf("Part 10d RDS: %s",
                  if (file.exists(.p10d))
                      file.path(basename(dirname(.p10d)), basename(.p10d))
                  else paste("not found in any location:",
                             paste(p10e_p10d_candidates(), collapse = "  |  "))))

    .have_x <- !is.null(IPRS_IP_DFG)   && length(IPRS_IP_DFG)   > 0L
    .have_y <- !is.null(IPRS_IP_DFG_Y) && length(IPRS_IP_DFG_Y) > 0L
    .y_str  <- if (!.have_y) "none" else
        paste(sprintf("%s=%s", names(IPRS_IP_DFG_Y),
                      vapply(as.numeric(IPRS_IP_DFG_Y),
                             function(z) if (is.finite(z)) sprintf("%.4f", z)
                                         else "NA",
                             character(1))), collapse = ", ")
    .y_gap <- if (!.have_x) character(0) else {
        .xn <- names(IPRS_IP_DFG)
        .xn[vapply(.xn, function(nm) {
            .xv <- suppressWarnings(as.numeric(IPRS_IP_DFG[nm]))
            .yv <- if (.have_y && nm %in% names(IPRS_IP_DFG_Y))
                suppressWarnings(as.numeric(IPRS_IP_DFG_Y[nm])) else NA_real_
            isTRUE(is.finite(.xv)) && !isTRUE(is.finite(.yv))
        }, logical(1))]
    }
    gr10e("GR-07-IPRS",
          if (.have_x && .have_y && length(.y_gap) == 0L) "PASS" else "WARN",
          if (.have_x)
              sprintf(paste0("IPRS IP values (%s): IP_x %s | IP_y %s%s"),
                      ip_from,
                      paste(sprintf("%s=%.4f", names(IPRS_IP_DFG),
                                    as.numeric(IPRS_IP_DFG)), collapse = ", "),
                      .y_str,
                      if (!.have_y)
                          "  -- IP_y absent, Step 4/5 skipped"
                      else if (length(.y_gap) > 0L)
                          sprintf(paste0("  -- no finite IP_y for %s, ",
                                         "Step 4/5 skipped there"),
                                  paste(.y_gap, collapse = ", "))
                      else "")
          else paste0("IPRS IP values not available; Steps 1-3 run and ",
                      "Step 4/5 is skipped"))

    gr10e_write()
    cat(sprintf("\n  Pre-analysis guardrails: %d PASS / %d WARN / %d FAIL\n",
                .gr10e_pass, .gr10e_warn, .gr10e_fail))
    if (.gr10e_fail > 0) {
        cat("  Guardrail failures: review before interpreting results.\n\n")
        warning("Part 10e pre-analysis guardrails: ", .gr10e_fail, " FAIL(s)",
                call. = FALSE)
    } else if (.gr10e_warn > 0) {
        cat(sprintf("  %d passed, %d warning(s), no failures.\n\n",
                    .gr10e_pass, .gr10e_warn))
    } else {
        cat("  All pre-analysis guardrails passed.\n\n")
    }


    ##########################################################################
    # SECTION B: MAIN FUNCTION (runs per dataset)
    ##########################################################################
    run_part10e <- function(dataset_tag) {

        set.seed(42)

        cat("\n", strrep("=", 72), "\n")
        cat("  Part10e CMRF:", dataset_tag, " [", ARM_TAG, "]\n")
        cat(strrep("=", 72), "\n\n")

        all_plots  <- list()
        excel_tabs <- list()

        # == Load Part10d RDS ================================================
        rds_path <- p10e_p10d_rds()
        if (!file.exists(rds_path)) {
            cat("  ERROR: Part 10d RDS not found. Looked in:\n")
            for (.p in p10e_p10d_candidates()) cat("      ", .p, "\n", sep = "")
            return(invisible(NULL))
        }
        all_results_10d <- readRDS(rds_path)
        res10d <- all_results_10d[[dataset_tag]]
        if (is.null(res10d)) {
            cat("  ERROR: No Part10d results for dataset:", dataset_tag, "\n")
            return(invisible(NULL))
        }

        Si_obs  <- res10d$Si_obs
        Ki_obs  <- res10d$Ki_obs
        u       <- res10d$u
        v       <- res10d$v
        n       <- length(Si_obs)

        # == IP_x from Part10b_IPRS (SKCS-core derived) =========================
        # IP_u is computed later from IP_x on the full-population copula axis.
        IP_x          <- NA_real_
        IP_y          <- NA_real_      # IP_y, the fitted cDFG at IP_x
        IP_v          <- NA_real_      # its projection onto the cDFG marginal
        IP_u_10d      <- res10d$IP_u   # keep Part10d value for comparison log
        IP_u          <- NA_real_      # will be set after scores_goi is loaded
        IP_u_source   <- "pending"

        if (!is.null(IPRS_IP_DFG) && dataset_tag %in% names(IPRS_IP_DFG) &&
            !is.na(IPRS_IP_DFG[dataset_tag])) {
            IP_x <- as.numeric(IPRS_IP_DFG[dataset_tag])
            cat(sprintf("  IP_x = %.6f (SKCS-core, from Part10b_IPRS)\n", IP_x))
            if (!is.null(IPRS_IP_DFG_Y) && dataset_tag %in% names(IPRS_IP_DFG_Y) &&
                !is.na(IPRS_IP_DFG_Y[dataset_tag])) {
                IP_y <- as.numeric(IPRS_IP_DFG_Y[dataset_tag])
                cat(sprintf("  IP_y = %.6f = f(IP_x), from Part10b_IPRS\n", IP_y))
            } else {
                stop("Part10b_IPRS supplied IP_x for ", dataset_tag,
                     " but no IP_y. Rerun the IPRS part for this arm.",
                     call. = FALSE)
            }
        } else {
            # Fallback: Part10d IP_u for the band kernel only; Step 4/5 is skipped.
            IP_u        <- IP_u_10d
            IP_u_source <- "Part10d_fallback"
            cat(sprintf(paste0("  WARNING: IPRS IP_x not available for %s; ",
                               "Part10d IP_u = %s kept for the band kernel; Step 4/5 skipped.\n"),
                        dataset_tag,
                        if (is.finite(IP_u)) sprintf("%.4f", IP_u) else "NA"))
        }

        cat(sprintf("  Loaded %d SKCS-core cells from Part10d RDS\n", n))
        if (!is.na(IP_u_10d))
            cat(sprintf("  Part10d IP_u = %.4f (full-pop; for comparison only)\n",
                        IP_u_10d))
        cat("  Spine: fitted in this part\n")

        # == Load Part10b CSV for group labels and raw EPi ===================
        ep_csv <- p10e_ep_csv(dataset_tag)
        if (!file.exists(ep_csv)) {
            cat("  ERROR: Part 10b EP/TREP CSV not found. Looked in:\n")
            for (.p in p10e_ep_candidates(dataset_tag)) cat("      ", .p, "\n", sep = "")
            return(invisible(NULL))
        }
        # check.names = FALSE, as in Part 10b.
        skcs_df <- read.csv(ep_csv, stringsAsFactors = FALSE, check.names = FALSE)

        required_cols <- c("barcode", "IP_status", "is_TREP", "EP_value",
                           "S_value", "K_value")
        missing_cols  <- setdiff(required_cols, names(skcs_df))
        if (length(missing_cols) > 0) {
            cat("  ERROR: Missing columns in Part10b CSV:",
                paste(missing_cols, collapse = ", "), "\n")
            return(invisible(NULL))
        }

        # Align the EP CSV to the Part 10d cells by barcode.
        bc_10d <- res10d$barcode
        if (is.null(bc_10d))
            stop("Part 10d RDS carries no barcode vector.",
                 call. = FALSE)
        .m <- match(bc_10d, skcs_df$barcode)
        if (anyNA(.m))
            stop(sprintf(paste0("%s: %d of %d Part 10d barcode(s) are absent from ",
                                "the Part 10b EP/TREP CSV."),
                         dataset_tag, sum(is.na(.m)), length(bc_10d)), call. = FALSE)
        skcs_df <- skcs_df[.m, , drop = FALSE]
        rownames(skcs_df) <- NULL
        cat(sprintf("  EP/TREP CSV aligned to the Part 10d cells by barcode: %d cells\n",
                    nrow(skcs_df)))

        # Donor matched from the gam_data CSVs by barcode.
        donor_map   <- p10e_donor_map(dataset_tag)
        donor_cells <- unname(donor_map[bc_10d])
        if (anyNA(donor_cells))
            stop(sprintf(paste0("%s: %d of %d barcode(s) have no donor in the ",
                                "gam_data CSVs."),
                         dataset_tag, sum(is.na(donor_cells)), length(bc_10d)),
                 call. = FALSE)
        donor_cells <- as_donor_factor_10e(donor_cells)
        cat(sprintf("  Donors on the SKCS-core cells: %d\n",
                    n_donor_levels_10e(donor_cells)))

        ip_stat  <- skcs_df$IP_status
        # is_TREP must coerce to TRUE/FALSE.
        is_trep  <- as.logical(skcs_df$is_TREP)
        if (any(is.na(is_trep))) {
            .bad <- unique(as.character(skcs_df$is_TREP)[is.na(is_trep)])
            stop(sprintf(paste0("%s: %d of %d is_TREP value(s) did not coerce ",
                                "to TRUE/FALSE (saw: %s). Check ",
                                "<prefix>_%s_EP_TREP_values.csv."),
                        dataset_tag, sum(is.na(is_trep)), length(is_trep),
                        paste(utils::head(.bad, 4), collapse = ", "), dataset_tag),
                 call. = FALSE)
        }
        EP_raw   <- skcs_df$EP_value

        EP_clip_lo <- quantile(EP_raw, 0.02, na.rm = TRUE)
        EP_plot    <- pmax(EP_raw, EP_clip_lo)

        group <- paste(ip_stat, ifelse(is_trep, "TREP", "non-TREP"))
        # An IP_status that forms none of the four groups is a stop.
        .bad_grp <- !group %in% names(FOUR_COLS)
        if (any(.bad_grp)) {
            .bad_lab <- unique(group[.bad_grp])
            stop(sprintf(paste0("%s: %d of %d cell(s) carry an IP_status that ",
                                "does not form one of the four Pre/Post x TREP ",
                                "groups (saw: %s). Check ",
                                "<prefix>_%s_EP_TREP_values.csv."),
                        dataset_tag, sum(.bad_grp), length(group),
                        paste(utils::head(.bad_lab, 4), collapse = ", "),
                        dataset_tag),
                 call. = FALSE)
        }
        group <- factor(group, levels = names(FOUR_COLS))

        cat(sprintf("  Group counts: %s\n",
                    paste(names(table(group)), table(group), sep = "=",
                          collapse = "  ")))

        # =========================================================================
        # STEP 1 -- Rank-transform Si and Ki to uniform marginals p_i, q_i
        # =========================================================================
        cat("\n  STEP 1 -- Rank-transform to copula marginals\n")

        p_i <- rank(Si_obs, ties.method = "average") / (n + 1)
        q_i <- rank(Ki_obs, ties.method = "average") / (n + 1)

        cat(sprintf("  p_i range [%.4f, %.4f] | q_i range [%.4f, %.4f]\n",
                    min(p_i), max(p_i), min(q_i), max(q_i)))

        ks_p <- suppressWarnings(ks.test(p_i, "punif")$p.value)
        ks_q <- suppressWarnings(ks.test(q_i, "punif")$p.value)
        cat(sprintf("  KS test p_i ~ Uniform: p=%.4f  q_i ~ Uniform: p=%.4f\n",
                    ks_p, ks_q))



        # =========================================================================
        # STEP 2a -- Si-Ki copula statistics
        # =========================================================================
        cat("\n  STEP 2a -- Si-Ki copula statistics\n")

        df_sk <- data.frame(
            p_i     = p_i,
            q_i     = q_i,
            Si      = Si_obs,
            Ki      = Ki_obs,
            EP      = EP_plot,
            EP_raw  = EP_raw,
            group   = group,
            is_trep = is_trep,
            ip_stat = ip_stat
        )

        rho_sk <- cor(Si_obs, Ki_obs, method = "spearman")
        tau_sk <- cor(Si_obs, Ki_obs, method = "kendall")
        cat(sprintf("  Spearman rho(Si,Ki)=%.4f | Kendall tau(Si,Ki)=%.4f\n",
                    rho_sk, tau_sk))


        # =========================================================================
        # STEP 2d -- Full-population Si-Ki pitcher plant in copula space
        #
        cat("\n  STEP 2d -- Full-population Si-Ki pitcher plant copula\n")

        scores_csv <- p10e_scores_csv(dataset_tag)

        scores_goi    <- NULL
        p_full        <- NULL
        q_full        <- NULL
        full_group    <- NULL
        df_full_pop   <- NULL
        p_thresh      <- NULL
        q_thresh      <- NULL
        n_full_pop    <- NULL
        n_skcs_found  <- NULL
        skcs_barcodes <- skcs_df$barcode

        # Part 10b Pre/Post labels recomputed on both axes and reported in the
        # guardrail ledger; no label is overwritten.
        ip_axis_check <- list(n = 0L, agree_cdfg = NA_real_,
                              agree_goi = NA_real_, axis = "cDFG at IP_y",
                              note = "not computed")

        if (!file.exists(scores_csv)) {
            cat("  WARNING: Part10 scores CSV not found:", scores_csv, "\n")
            cat("  Plots 5, 6, 7 skipped (require full-pop scores CSV).\n\n")
        } else {

            scores_all <- read.csv(scores_csv, stringsAsFactors = FALSE)
            scores_goi <- scores_all[!is.na(scores_all$spear_contrib), ]
            rownames(scores_goi) <- NULL   # reset after subset to avoid row-name warnings
            n_full_pop <- nrow(scores_goi)
            cat(sprintf("  Full GOI-positive population: %d cells\n", n_full_pop))

            p_full <- rank(scores_goi$spear_norm,
                           ties.method = "average") / (n_full_pop + 1)
            q_full <- rank(scores_goi$kend_contrib,
                           ties.method = "average") / (n_full_pop + 1)

            skcs_label_map <- setNames(
                paste(skcs_df$IP_status,
                      ifelse(as.logical(skcs_df$is_TREP), "TREP", "non-TREP")),
                skcs_df$barcode
            )

            full_group_raw <- skcs_label_map[scores_goi$barcode]
            full_group_raw[is.na(full_group_raw)] <- "Non-selected"
            full_group_raw[!full_group_raw %in% c(names(FOUR_COLS), "Non-selected")] <-
                "Non-selected"

            group_levels_full <- c("Non-selected", names(FOUR_COLS))
            full_group <- factor(full_group_raw, levels = group_levels_full)

            # Recompute Part 10b's Pre/Post labels on both axes (cDFG: cDFG_expr < IP_y).
            .ax_i <- match(skcs_df$barcode, scores_goi$barcode)
            .ax_ok <- !is.na(.ax_i) & !is.na(ip_stat)
            if (!("cDFG_expr" %in% names(scores_goi))) {
                ip_axis_check$note <-
                    "scores CSV has no cDFG_expr column, so the cDFG rule cannot be reproduced"
            } else if (!any(.ax_ok)) {
                ip_axis_check$note <-
                    "no SKCS cell has both an IP_status and a row in the scores CSV"
            } else {
                .ax_cdfg <- as.numeric(scores_goi$cDFG_expr)[.ax_i[.ax_ok]]
                .ax_goi  <- if ("GOI_expr" %in% names(scores_goi))
                    as.numeric(scores_goi$GOI_expr)[.ax_i[.ax_ok]] else NULL
                .ax_obs  <- ip_stat[.ax_ok]
                .lab_cdfg <- if (is.finite(IP_y))
                    ifelse(.ax_cdfg <  IP_y, "Pre-IP", "Post-IP") else NULL
                .lab_goi  <- if (is.finite(IP_x) && !is.null(.ax_goi))
                    ifelse(.ax_goi  <= IP_x, "Pre-IP", "Post-IP") else NULL
                ip_axis_check <- list(
                    n          = sum(.ax_ok),
                    agree_cdfg = if (is.null(.lab_cdfg)) NA_real_
                                 else mean(.lab_cdfg == .ax_obs),
                    agree_goi  = if (is.null(.lab_goi))  NA_real_
                                 else mean(.lab_goi  == .ax_obs),
                    axis       = "cDFG at IP_y",
                    note       = if (sum(.ax_ok) < nrow(skcs_df))
                        sprintf("%d of %d SKCS cell(s) had no IP_status or no scores row",
                                nrow(skcs_df) - sum(.ax_ok), nrow(skcs_df))
                    else "")
                cat(sprintf(paste0("  Part 10b Pre/Post labels reproduce: %s on ",
                                   "the cDFG axis at IP_y, %s on the GOI axis ",
                                   "at IP_x  (n=%d)\n"),
                            if (is.finite(ip_axis_check$agree_cdfg))
                                sprintf("%.4f", ip_axis_check$agree_cdfg) else "NA",
                            if (is.finite(ip_axis_check$agree_goi))
                                sprintf("%.4f", ip_axis_check$agree_goi) else "NA",
                            ip_axis_check$n))
            }

            df_full_pop <- data.frame(
                p_full    = p_full,
                q_full    = q_full,
                group     = full_group,
                is_skcs   = scores_goi$barcode %in% skcs_barcodes,
                barcode   = scores_goi$barcode,
                Si_raw    = scores_goi$spear_norm,
                Ki_raw    = scores_goi$kend_contrib,
                stringsAsFactors = FALSE
            )
            if (all(c("u_full", "v_full") %in% names(scores_goi))) {
                df_full_pop$u_full <- scores_goi$u_full
                df_full_pop$v_full <- scores_goi$v_full
            }

            n_skcs_found  <- sum(df_full_pop$is_skcs)
            rho_full      <- cor(scores_goi$spear_norm, scores_goi$kend_contrib,
                                 method = "spearman")
            p_thresh <- mean(scores_goi$spear_norm <= 0.05)
            q_thresh <- mean(scores_goi$kend_contrib  <= 0.5)

            cat(sprintf("  SKCS matched: %d | Non-SKCS: %d\n",
                        n_skcs_found, sum(!df_full_pop$is_skcs)))
            cat(sprintf("  p_thresh=%.3f | q_thresh=%.3f\n", p_thresh, q_thresh))

            excel_tabs[["Step2d_FullPop_Summary"]] <- data.frame(
                dataset       = dataset_tag,
                n_full_pop    = n_full_pop,
                n_SKCS        = n_skcs_found,
                n_nonSKCS     = sum(!df_full_pop$is_skcs),
                rho_SK_full   = round(rho_full, 6),
                p_thresh_Sinorm05 = round(p_thresh, 6),
                q_thresh_Ki05 = round(q_thresh, 6),
                stringsAsFactors = FALSE
            )
        }


        # =========================================================================
        # STEP 3b -- Per-group Kendall tau and Spearman rho
        #
        cat("\n  STEP 3b -- Per-group dependence statistics\n")

        groups_uniq <- levels(group)
        group_stats <- lapply(groups_uniq, function(g) {
            idx <- group == g
            if (sum(idx) < 10) return(NULL)
            Si_g <- Si_obs[idx]; Ki_g <- Ki_obs[idx]
            data.frame(group = g, n = sum(idx),
                       tau_SK  = round(cor(Si_g, Ki_g, method = "kendall"),  4),
                       rho_SK  = round(cor(Si_g, Ki_g, method = "spearman"), 4),
                       mean_Si = round(mean(Si_g), 4), mean_Ki = round(mean(Ki_g), 4),
                       sd_Si   = round(sd(Si_g),   4), sd_Ki   = round(sd(Ki_g),   4),
                       stringsAsFactors = FALSE)
        })
        group_stats <- do.call(rbind, Filter(Negate(is.null), group_stats))

        overall_row <- data.frame(
            group   = "Overall (SKCS)", n = n,
            tau_SK  = round(tau_sk, 4),  rho_SK  = round(rho_sk, 4),
            mean_Si = round(mean(Si_obs), 4), mean_Ki = round(mean(Ki_obs), 4),
            sd_Si   = round(sd(Si_obs),   4), sd_Ki   = round(sd(Ki_obs),   4),
            stringsAsFactors = FALSE
        )
        group_stats_full <- rbind(overall_row, group_stats)

        # =========================================================================
        # STEP 4 -- CMRF: TWO-STAGE MEMBERSHIP FILTER
        #   Stage 1  band proximity, |v_i - g(u_i)| <= delta_adapt(u_i), with the
        #            SKCS core (Si_norm > 0.05, Ki > 0.5) as the anchor set.
        #   Stage 2  RSk-fCV dual gate on the non-SKCS Stage 1 candidates.
        #   Classification:
        #     SKCS core + pass Stage 1            -> "Original + Refined"
        #     SKCS core + fail Stage 1            -> "Original only"
        #     non-SKCS  + pass Stage 1 and 2      -> "Recovered"
        #     everything else                     -> "Excluded"
        # =========================================================================

        cat("\n  STEP 4 -- CMRF two-stage membership filter\n")

        if (is.null(scores_goi)) {
            cat("  Step 4 skipped: full-population scores CSV not available.\n")
        } else if (!is.finite(IP_y) ||
                   (!is.finite(IP_x) && !is.finite(IP_u_10d))) {
            # No IP_y or no inflection point: Steps 1-3 stand, IP-dependent work is skipped.
            cat("  Step 4/5 skipped: ",
                if (!is.finite(IP_y))
                    "no IP_y, the cDFG-axis classification threshold.\n"
                else "no inflection point available.\n", sep = "")
            cat(sprintf("    IP_x from the IPRS RDS : %s\n",
                        if (is.finite(IP_x)) sprintf("%.6f", IP_x) else "NA"))
            cat(sprintf("    IP_y from the IPRS RDS : %s\n",
                        if (is.finite(IP_y)) sprintf("%.6f", IP_y) else "NA"))
            cat(sprintf("    IP_u from Part 10d     : %s\n",
                        if (is.finite(IP_u_10d)) sprintf("%.6f", IP_u_10d)
                        else "NA"))
            if (!is.finite(IP_y) && !is.finite(IP_x) && is.finite(IP_u_10d))
                cat("  NOTE: Part 10d IP_u present but no IP_x or IP_y; the Part 10d RDS is stale.\n", sep = "")
            cat("  Not written for ", dataset_tag, ":\n", sep = "")
            for (.miss in c("_05_waist_zoom_with_orig_only.pdf/.png",
                            "_06_CMRF_band_uv_space.pdf/.png",
                            "_07_recovery_bar_comparison.pdf/.png",
                            "_01_CMRF_diagnostics.pdf/.png",
                            "_RSk_fCV_Audit.xlsx"))
                cat("      ", FILE_PREFIX, "_", dataset_tag, .miss, "\n", sep = "")
            cat("  Needs ", if (!is.finite(IP_y)) "IP_y (inflection_points_DFG_y)"
                else "an inflection point", " from Part10b_CP10K_IPRS_code.r.\n", sep = "")
        } else {

            # == Build or retrieve (u, v) for full population ================
            has_uv_in_csv <- all(c("u_full", "v_full") %in% names(scores_goi))

            if (has_uv_in_csv) {
                u_all <- scores_goi$u_full
                v_all <- scores_goi$v_full
                cat("  Using u_full/v_full from scores CSV\n")
            } else if (all(c("GOI_expr", "cDFG_expr") %in% names(scores_goi))) {
                u_all <- rank(scores_goi$GOI_expr,
                              ties.method = "average") / (n_full_pop + 1)
                v_all <- rank(scores_goi$cDFG_expr,
                              ties.method = "average") / (n_full_pop + 1)
                cat("  Reconstructed u_all/v_all from GOI_expr/cDFG_expr ranks\n")
            } else {
                u_all <- p_full
                v_all <- q_full
                cat("  WARNING: using Si/Ki ranks as u/v proxies (less accurate)\n")
            }

            # == Compute IP_u from SKCS-core IP_x on full-pop copula axis ==========
            # IP_u = #{GOI_j <= IP_x, j in all GOI+} / (n_full + 1)
            if (!is.na(IP_x) && "GOI_expr" %in% names(scores_goi)) {
                IP_u <- sum(scores_goi$GOI_expr <= IP_x) / (n_full_pop + 1)
                IP_u_source <- "SKCS_IPRS"
                cat(sprintf("  IP_u = %.6f  (IP_x=%.6f projected onto full-pop copula; n=%d)\n",
                            IP_u, IP_x, n_full_pop))
                if (!is.na(IP_u_10d))
                    cat(sprintf("  cf. Part10d IP_u = %.6f  (delta = %.6f)\n",
                                IP_u_10d, IP_u - IP_u_10d))
                # IP_v = Fc(IP_y) on the full-population axis.
                if (is.finite(IP_y) && "cDFG_expr" %in% names(scores_goi)) {
                    IP_v <- sum(scores_goi$cDFG_expr <= IP_y) / (n_full_pop + 1)
                    cat(sprintf(paste0("  IP_v = %.6f  (IP_y=%.6f projected onto",
                                       " full-pop copula; n=%d)\n"),
                                IP_v, IP_y, n_full_pop))
                }
            } else if (is.na(IP_u) && !is.na(IP_u_10d)) {
                IP_u <- IP_u_10d
                IP_u_source <- "Part10d_fallback"
                cat(sprintf("  IP_u = %.6f (Part10d fallback)\n", IP_u))
            }

            # == Fit the spine g(u) on TREP SKCS-core cells ========================
            cat("  Fitting the spine g(u) on TREP SKCS-core cells ...\n")
            trep_idx_e <- which(is_trep)
            if (length(trep_idx_e) < MIN_SPINE_CELLS) {
                cat(sprintf(paste0("  ERROR: only %d TREP SKCS-core cell(s); ",
                                   "minimum is %d. Step 4/5 skipped.\n"),
                            length(trep_idx_e), MIN_SPINE_CELLS))
                return(invisible(NULL))
            }
            gam_result <- fit_gam_spine_prss(
                u_trep     = u[trep_idx_e],
                v_trep     = v[trep_idx_e],
                donor      = donor_cells[trep_idx_e],
                label      = "g(u) CEP-IP spine",
                phase2_max = PRSS_PHASE2_MAX_P10E,
                seed_base  = 42L)
            gam_obj <- gam_result$best_model
            if (is.null(gam_obj)) {
                cat("  ERROR: spine fit failed -- Step 4/5 skipped.\n")
                return(invisible(NULL))
            }
            gam_pred_var  <- "u_trep"
            spine_n_cells <- length(trep_idx_e)
            spine_u_range <- range(u[trep_idx_e])
            cat(sprintf(paste0("    g(u): n=%d | best_k=%d | donor RE=%s (%d donors) ",
                               "| weights=%s | PRSS=%.5f | iters=%d\n"),
                        spine_n_cells, gam_result$best_k,
                        if (gam_result$donor_re) "YES" else "NO",
                        gam_result$n_donors,
                        if (isTRUE(W_DONOR_WEIGHTS_10E)) "inverse donor" else "uniform",
                        gam_result$best_prss, gam_result$actual_iters))
            excel_tabs[["Step4_Spine_PRSS"]] <- data.frame(
                dataset       = dataset_tag,
                iter          = seq_along(gam_result$prss_values),
                prss          = gam_result$prss_values,
                best_iter     = seq_along(gam_result$prss_values) == gam_result$best_iter,
                best_k        = gam_result$best_k,
                max_k         = gam_result$max_k,
                phase2_max    = gam_result$phase2_max,
                n_iter_max    = gam_result$n_iter_max,
                early_stopped = gam_result$early_stopped,
                n_cells       = spine_n_cells,
                donor_re      = gam_result$donor_re,
                n_donors      = gam_result$n_donors,
                stringsAsFactors = FALSE
            )

            # Population-level curve (donor term excluded).
            gam_predict <- function(u_vals) spine_predict_pop(gam_obj, u_vals)

            # Monotonicity check on the population-level curve.
            spine_mono <- spine_monotone_check(gam_obj)
            cat(sprintf(paste0("    monotonicity: min consecutive difference = ",
                               "%.3e over %d grid points -- %s\n"),
                        spine_mono$min_diff, spine_mono$n_grid,
                        if (spine_mono$monotone) "non-decreasing"
                        else sprintf("DECREASING near u = %.4f", spine_mono$worst_u)))

            # Where the IP pair sits relative to the spine (reported only).
            g_at_IP_u  <- if (is.finite(IP_u)) gam_predict(IP_u) else NA_real_
            g_inv_IP_v <- NA_real_
            if (is.finite(IP_v)) {
                .ug <- seq(0, 1, length.out = 2001L)
                .gg <- gam_predict(.ug)
                if (IP_v >= min(.gg) && IP_v <= max(.gg))
                    g_inv_IP_v <- .ug[which.min(abs(.gg - IP_v))]
            }
            cat(sprintf(paste0("    spine at the IP: g(IP_u) = %s (IP_v = %s)",
                               " | g_inv(IP_v) = %s (IP_u = %s)\n"),
                        if (is.finite(g_at_IP_u))  sprintf("%.4f", g_at_IP_u)  else "NA",
                        if (is.finite(IP_v))      sprintf("%.4f", IP_v)      else "NA",
                        if (is.finite(g_inv_IP_v)) sprintf("%.4f", g_inv_IP_v)
                        else "NA (IP_v outside the fitted range)",
                        if (is.finite(IP_u))      sprintf("%.4f", IP_u)      else "NA"))

            # Third subtitle row of every figure that draws the IP.
            ip_coords_line <- paste0(
                sprintf("IP_u = %s | IP_v = %s | IP_x = %s | IP_y = %s | ",
                        if (is.finite(IP_u)) sprintf("%.4f", IP_u) else "NA",
                        if (is.finite(IP_v)) sprintf("%.4f", IP_v) else "NA",
                        if (is.finite(IP_x)) sprintf("%.4f", IP_x) else "NA",
                        if (is.finite(IP_y)) sprintf("%.4f", IP_y) else "NA"),
                sprintf("g(IP_u) = %s | g_inv(IP_v) = %s",
                        if (is.finite(g_at_IP_u))  sprintf("%.4f", g_at_IP_u)  else "NA",
                        if (is.finite(g_inv_IP_v)) sprintf("%.4f", g_inv_IP_v) else "NA"))
            cat("  IP coordinates: ", ip_coords_line, "\n", sep = "")

            # Horizontal boundary at IP_v, which the Pre/Post colours obey; none
            # when IP_v is not finite.
            IP_V_COL_10E <- "#8B1A1A"
            ip_boundary_v <- function(linewidth = 0.55) {
                if (!is.finite(IP_v)) return(list())
                list(geom_hline(yintercept = IP_v, colour = IP_V_COL_10E,
                                linetype = "dashed", linewidth = linewidth))
            }

            # goi_at_ip is reported in Step 3e.
            goi_at_ip <- if ("GOI_expr" %in% names(scores_goi)) {
                quantile(scores_goi$GOI_expr, probs = IP_u, type = 7, na.rm = TRUE)
            } else NA_real_

            # == Step 3e: threshold translation ===================================
            cat("\n  STEP 3e -- Threshold translation summary\n")

            goi_gene      <- "MT-CO3"   # GOI is MT-CO3 for both NeoC and NonNeoC (R_config.R)
            p_ip_full_val <- mean(u_all <= IP_u)
            p_thresh_val  <- mean(scores_goi$spear_norm <= 0.05)
            si_at_pthresh <- quantile(scores_goi$spear_norm,
                                      probs = p_thresh_val, type = 7, na.rm = TRUE)
            goi_at_sinorm05 <- if ("GOI_expr" %in% names(scores_goi)) {
                quantile(scores_goi$GOI_expr, probs = p_thresh_val, type = 7, na.rm = TRUE)
            } else NA_real_
            q_thresh_val  <- mean(scores_goi$kend_contrib <= 0.5)
            ki_at_qthresh <- quantile(scores_goi$kend_contrib,
                                      probs = q_thresh_val, type = 7, na.rm = TRUE)
            q_v05_full    <- mean(v_all <= 0.5)

            df_thresh_translation <- data.frame(
                section = c(
                    rep("(a) IP_u -> GOI value", 3),
                    rep("(b) IP_u -> p_full",    5),
                    rep("(c) Si_norm=0.05 boundary -> p_full", 7),
                    rep("(d) Ki=0.5 boundary -> q_full", 7)
                ),
                step = c(
                    "Input: IP_u (copula u where g(IP_u)=0.5)",
                    "Formula: GOI_IP = quantile(GOI_expr, probs=IP_u, type=7)",
                    sprintf("Result: GOI_IP (%s expression at IP)", goi_gene),
                    "Input: IP_u",
                    "Formula: p_IP_full = mean(u_all <= IP_u)",
                    "Expansion: counts GOI-pos cells with u rank <= IP_u",
                    "Result: p_IP_full",
                    "Interpretation: X% of cells have GOI rank at or below IP rank",
                    "Definition: Si_norm = spear_norm",
                    "Threshold: Si_norm > 0.05",
                    "Formula: p_thresh = mean(spear_norm <= 0.05)",
                    "Result: p_thresh",
                    "Verification: quantile(Si_norm, probs=p_thresh) should be ~0.05",
                    sprintf("GOI value at Si_norm=0.05 [%s]", goi_gene),
                    "Delta: p_thresh - p_IP_full",
                    "Definition: Ki = kend_contrib",
                    "Threshold: Ki > 0.5",
                    "Formula: q_thresh = mean(kend_contrib <= 0.5)",
                    "Result: q_thresh",
                    "Verification: quantile(Ki, probs=q_thresh) should be ~0.5",
                    "Reference: q_v05_full = mean(v_all <= 0.5)",
                    "Delta: q_thresh - q_v05_full"
                ),
                value = c(
                    sprintf("%.6f", IP_u),
                    sprintf("quantile(GOI_expr, probs=%.6f, type=7)", IP_u),
                    sprintf("%.4f", ifelse(is.na(goi_at_ip), NA_real_, goi_at_ip)),
                    sprintf("%.6f", IP_u),
                    "mean(u_all <= IP_u)",
                    sprintf("mean(u_all <= %.6f)", IP_u),
                    sprintf("%.6f", p_ip_full_val),
                    sprintf("%.1f%%", 100 * p_ip_full_val),
                    "spear_norm column",
                    "Si_norm > 0.05; boundary Si_norm = 0.05",
                    "mean(spear_norm <= 0.05)",
                    sprintf("%.6f", p_thresh_val),
                    sprintf("quantile(Si_norm, %.6f) = %.6f", p_thresh_val, si_at_pthresh),
                    sprintf("%.4f", ifelse(is.na(goi_at_sinorm05), NA_real_, goi_at_sinorm05)),
                    sprintf("%.6f - %.6f = %.6f", p_thresh_val, p_ip_full_val,
                            p_thresh_val - p_ip_full_val),
                    "kend_contrib column",
                    "Ki > 0.5; boundary Ki = 0.5",
                    "mean(kend_contrib <= 0.5)",
                    sprintf("%.6f", q_thresh_val),
                    sprintf("quantile(Ki, %.6f) = %.6f", q_thresh_val, ki_at_qthresh),
                    sprintf("mean(v_all <= 0.5) = %.6f", q_v05_full),
                    sprintf("%.6f - %.6f = %.6f", q_thresh_val, q_v05_full,
                            q_thresh_val - q_v05_full)
                ),
                dataset = dataset_tag,
                stringsAsFactors = FALSE
            )
            excel_tabs[["Step3e_Threshold_Translation"]] <- df_thresh_translation
            cat("  Step 3e threshold translation table built.\n")

            # =====================================================================
            # STEP 4 (continued) -- CMRF two-stage filter
            # =====================================================================

            # == Original SKCS membership flag ===============================
            is_orig_skcs <- scores_goi$barcode %in% skcs_barcodes
            cat(sprintf("  Original SKCS-core: %d cells (of %d full-pop)\n",
                        sum(is_orig_skcs), n_full_pop))

            # Calibration weights recomputed on the original SKCS-core cells; other
            # cells keep weight 1.
            donor_full <- unname(p10e_donor_map(dataset_tag)[scores_goi$barcode])
            if (anyNA(donor_full))
                stop(sprintf(paste0("%s: %d full-population barcode(s) have no ",
                                    "donor in the gam_data CSVs."),
                             dataset_tag, sum(is.na(donor_full))), call. = FALSE)
            w_all <- rep(1, length(u_all))
            w_all[is_orig_skcs] <- donor_weights_10e(donor_full[is_orig_skcs])
            cat(sprintf(paste0("  Calibration weights: %d donor(s) over the ",
                               "original SKCS-core; sum(w) = %.4f\n"),
                        n_donor_levels_10e(donor_full[is_orig_skcs]),
                        sum(w_all[is_orig_skcs])))

            # == Pre-compute the spine at all full-population u coordinates ====
            # g(u_i) for every cell (donor term excluded).
            g_u_all_cmrf <- gam_predict(u_all)

            # == Stage 1: delta calibration from original SKCS-core residuals ====
            # delta_base = weighted quantile(|v - g(u)|, 0.95) over all original
            # SKCS-core cells.
            cmrf_delta_cal <- calibrate_cmrf_delta(
                u_cal     = u_all[is_orig_skcs],
                v_cal     = v_all[is_orig_skcs],
                g_u_cal   = g_u_all_cmrf[is_orig_skcs],
                w_cal     = w_all[is_orig_skcs],
                coverage  = CMRF_COVERAGE
            )
            delta <- cmrf_delta_cal$delta_base

            # == Apply CMRF (Stage 1 + Stage 2) ==============================
            # Stage 1 band: |v - g(u)| <= delta_adapt(u)  [all cells]
            #   delta_adapt(u) = delta_base*(1+kappa*K(u,IP_u,bw))
            #   delta_base from the 95th weighted percentile of |v-g(u)|
            #   kappa=0.50, bw=0.10: widens at waist, narrows in arms
            # Stage 2: RSk-fCV fidelity gate on non-SKCS Stage 1 candidates
            cat("\n  Applying CMRF Stage 1 (band proximity for all cells) ...\n")
            cmrf_res <- apply_cmrf(
                u_all        = u_all,
                v_all        = v_all,
                g_u_all      = g_u_all_cmrf,
                is_orig_skcs = is_orig_skcs,
                IP_u         = IP_u,
                delta_base   = delta,
                kappa        = CMRF_KAPPA,
                bw           = CMRF_BW,
                n_repeats    = CMRF_N_REPEATS,
                k_folds      = CMRF_K_FOLDS,
                alpha        = CMRF_ALPHA,
                seed         = 42L
            )

            g_u_all_res      <- cmrf_res$g_u_all        # spine (= g_u_all_cmrf)
            resid_abs        <- cmrf_res$resid_abs       # |v - g(u)|
            delta_adapt      <- cmrf_res$delta_adapt     # per-cell adaptive tolerance
            K_weight_res     <- cmrf_res$K_weight        # Gaussian kernel weights
            yhat_iso         <- cmrf_res$yhat_iso        # PAVA step (diagnostic only)
            cond2_cmrf       <- cmrf_res$cond2           # spine proximity within delta_adapt
            cell_class       <- cmrf_res$cell_class
            is_recovered     <- cmrf_res$is_recovered
            is_refined_skcs  <- cmrf_res$pass_cmrf       # Stage1 AND Stage2 pass
            removed_s2       <- cmrf_res$removed_by_stage2

            cat(sprintf("  delta_base = %.6f (95th pctile, weighted, all original SKCS-core; n=%d)\n",
                        delta, cmrf_delta_cal$n_grp1a))
            cat(sprintf("  Delta cal resid: median=%.6f  max=%.6f\n",
                        cmrf_delta_cal$resid_median, cmrf_delta_cal$resid_max))
            cat(sprintf("  delta_adapt range: [%.6f, %.6f] (kappa=%.2f, bw=%.2f)\n",
                        min(delta_adapt), max(delta_adapt), CMRF_KAPPA, CMRF_BW))
            if (!is.null(cmrf_res$stage2_res)) {
                cat(sprintf("  Stage 2 RSk-fCV Wilcoxon p (final): %.4f | Removed by gate: %d\n",
                            cmrf_res$stage2_res$wilcox_p_final,
                            cmrf_res$stage2_res$n_removed))
            }
            cat(sprintf("  Original + Refined (Orig+Ref anchor):  %d\n",
                        sum(cell_class == "Original + Refined")))
            cat(sprintf("  Original only (failed Stage 1 -- outside band):  %d\n",
                        sum(cell_class == "Original only")))
            cat(sprintf("  Recovered (Grp 1b, passed both stages):%d\n",
                        sum(is_recovered)))
            cat(sprintf("  Excluded:                              %d\n",
                        sum(cell_class == "Excluded")))
            cat(sprintf("  CMRF-core (union):                     %d\n",
                        sum(is_refined_skcs)))

            # == Build df_refined ============================================
            # Recovered cells labelled on the cDFG axis at IP_y, like every other class.
            if (!is.finite(IP_y))
                stop("Cannot label the Recovered cells: IP_y is not available ",
                     "for ", dataset_tag, ".", call. = FALSE)
            if (!"cDFG_expr" %in% names(scores_goi))
                stop("Cannot label the Recovered cells: scores_goi has no ",
                     "cDFG_expr column.", call. = FALSE)
            .cdfg_all <- as.numeric(scores_goi$cDFG_expr)
            if (anyNA(.cdfg_all[is_recovered]))
                stop(sum(is.na(.cdfg_all[is_recovered])), " Recovered cell(s) ",
                     "have no cDFG_expr, so they cannot be classified on the ",
                     "cDFG axis.", call. = FALSE)
            group_refined_raw <- as.character(full_group)
            group_refined_raw[is_recovered] <- ifelse(
                .cdfg_all[is_recovered] < IP_y,
                "Recovered Pre-IP",
                "Recovered Post-IP"
            )
            cat(sprintf(paste0("  Recovered labelled on the cDFG axis at",
                               " IP_y = %.6f: %d Pre-IP / %d Post-IP\n"),
                        IP_y,
                        sum(group_refined_raw[is_recovered] == "Recovered Pre-IP"),
                        sum(group_refined_raw[is_recovered] == "Recovered Post-IP")))

            df_refined <- data.frame(
                barcode          = unname(scores_goi$barcode),
                u_all            = unname(u_all),
                v_all            = unname(v_all),
                p_full           = unname(p_full),
                q_full           = unname(q_full),
                # Upstream scores from Part10 CSV
                region           = unname(scores_goi$region),
                GOI_expr         = unname(scores_goi$GOI_expr),
                cDFG_expr        = unname(scores_goi$cDFG_expr),
                spear_contrib    = unname(scores_goi$spear_contrib),
                spear_norm       = unname(scores_goi$spear_norm),
                kend_contrib     = unname(scores_goi$kend_contrib),
                is_GOI_neg       = unname(scores_goi$is_GOI_neg),
                is_trep          = unname(scores_goi$is_trep),
                # CMRF core outputs
                yhat_iso         = unname(yhat_iso),
                resid_abs        = unname(resid_abs),
                delta_base       = delta,
                delta_adapt      = unname(delta_adapt),
                pass_cmrf        = unname(is_refined_skcs),
                pass_stage1      = unname(cmrf_res$pass_stage1),
                pass_stage2      = unname(cmrf_res$pass_stage2),
                removed_by_stage2 = unname(removed_s2),
                is_orig_skcs     = unname(is_orig_skcs),
                is_refined       = unname(is_refined_skcs),
                is_recovered     = unname(is_recovered),
                cell_class       = unname(cell_class),
                group_orig       = unname(as.character(full_group)),
                group_refined    = unname(group_refined_raw),
                g_u              = unname(g_u_all_cmrf),
                K_density        = unname(K_weight_res),
                cond2_refined    = unname(cond2_cmrf),
                stringsAsFactors = FALSE
            )

            # =====================================================================
            # PLOT 1 -- CMRF diagnostic panel (spine + adaptive delta band)
            # Panel A: spine g(u) + adaptive delta_adapt(u) band in (u,v) space.
            #          PAVA step shown as dotted reference.
            # Panel B: |v - g(u)| GAM residuals vs u, coloured by cell_class,
            #          adaptive delta_adapt(u) line shown.
            # Panel C: residual boxplot per class.
            # Panel D: Orig+Ref GAM residual distribution (delta_base calibration).
            # =====================================================================
            cat("\n  STEP 3f -- CMRF diagnostic plots (spine + adaptive delta band)\n")

            # Compute GAM spine and adaptive tau on a dense grid for smooth ribbon
            u_diag_grid   <- seq(0, 1, by = 0.005)
            g_diag        <- gam_predict(u_diag_grid)
            Kw_diag       <- exp(-(u_diag_grid - IP_u)^2 / (2 * 0.10^2))
            delta_diag      <- delta * (1 + 0.50 * Kw_diag)
            df_band_diag  <- data.frame(
                u     = u_diag_grid,
                g     = g_diag,
                upper = pmin(g_diag + delta_diag, 1),
                lower = pmax(g_diag - delta_diag, 0)
            )

            # PAVA step for reference overlay in Panel A
            df_1a_uv  <- df_refined[df_refined$is_orig_skcs, ]
            ord_1a_uv <- order(df_1a_uv$u_all)
            df_iso_uv <- data.frame(
                u = df_1a_uv$u_all[ord_1a_uv],
                y = df_1a_uv$yhat_iso[ord_1a_uv]
            )

            # Prepare multi-group colouring for Panel A (same scheme as Plot 6)
            RECOV_COL_03F      <- "#7C3AED"
            ORIG_ONLY_COL_03F  <- "#8B6914"
            ZOOM_COLS_03F <- c("Non-selected" = "gray60", FOUR_COLS)

            df_03f_all <- df_refined
            df_03f_all$plot_colour <- dplyr::case_when(
                df_03f_all$is_recovered ~ "Recovered",
                df_03f_all$cell_class == "Original + Refined" &
                    df_03f_all$group_orig == "Pre-IP TREP"      ~ "Pre-IP TREP",
                df_03f_all$cell_class == "Original + Refined" &
                    df_03f_all$group_orig == "Pre-IP non-TREP"  ~ "Pre-IP non-TREP",
                df_03f_all$cell_class == "Original + Refined" &
                    df_03f_all$group_orig == "Post-IP TREP"     ~ "Post-IP TREP",
                df_03f_all$cell_class == "Original + Refined" &
                    df_03f_all$group_orig == "Post-IP non-TREP" ~ "Post-IP non-TREP",
                TRUE                                             ~ "Excluded"
            )

            df_03f_nonsel   <- df_03f_all[!df_03f_all$is_orig_skcs &
                                              !df_03f_all$is_recovered &
                                              df_03f_all$cell_class != "Original only", ]
            df_03f_nonsel$plot_colour <- "Non-selected"
            df_03f_skcs     <- df_03f_all[df_03f_all$cell_class == "Original + Refined", ]
            df_03f_recov    <- df_03f_all[df_03f_all$is_recovered, ]
            df_03f_origonly <- df_03f_all[df_03f_all$cell_class == "Original only", ]

            # Factors with the full level set keep five legend keys when a group is empty.
            df_03f_nonsel$plot_colour <- factor(df_03f_nonsel$plot_colour,
                                                levels = names(ZOOM_COLS_03F))
            df_03f_skcs$plot_colour   <- factor(df_03f_skcs$plot_colour,
                                                levels = names(ZOOM_COLS_03F))

            p_f1 <- ggplot() +
                # Adaptive delta_adapt(u) band around spine
                geom_ribbon(data = df_band_diag,
                            aes(x = u, ymin = lower, ymax = upper),
                            fill = "#607D8B", alpha = 0.20, inherit.aes = FALSE) +
                # Non-selected background
                geom_point(data = df_03f_nonsel,
                           aes(x = u_all, y = v_all, colour = plot_colour),
                           size = 0.35, alpha = 0.30) +
                # Original-only (tan triangles, like Plot 6)
                geom_point(data = df_03f_origonly,
                           aes(x = u_all, y = v_all),
                           colour = "white", fill = ORIG_ONLY_COL_03F,
                           shape = 25, size = 2.0, alpha = 0.85, stroke = 0.5) +
                # Orig+Ref retained (four-group colours, like Plot 6)
                geom_point(data = df_03f_skcs,
                           aes(x = u_all, y = v_all, colour = plot_colour),
                           size = 0.9, alpha = 0.55) +
                # Recovered cells (violet diamonds, top layer, like Plot 6)
                geom_point(data = df_03f_recov,
                           aes(x = u_all, y = v_all),
                           colour = "white", fill = RECOV_COL_03F,
                           shape = 23, size = 2.0, alpha = 0.91, stroke = 0.6) +
                scale_colour_manual(values = ZOOM_COLS_03F, name = "Cell group",
                                    drop = FALSE,
                                    guide = guide_legend(
                                        override.aes = list(
                                            size  = c(1.8, rep(2.2, 4)),
                                            shape = rep(16L, 5),
                                            alpha = c(0.65, rep(0.75, 4))
                                        )
                                    )) +
                # spine g(u), solid navy
                geom_line(data = df_band_diag, aes(x = u, y = g),
                          colour = SPINE_COL_10E, linewidth = 0.9) +
                # PAVA step -- diagnostic reference only (dotted dark teal)
                geom_step(data = df_iso_uv, aes(x = u, y = y),
                          colour = "#004D40", linewidth = 0.5,
                          linetype = "dotted", direction = "hv") +
                # Reference lines: IP_u vertical, IP_v horizontal, v=0.5 horizontal
                geom_vline(xintercept = IP_u, colour = IP_U_COL_10E,
                           linetype = "solid", linewidth = 0.7) +
                ip_boundary_v(linewidth = 0.5) +
                geom_hline(yintercept = 0.5, colour = "#1A237E",
                           linetype = "dashed", linewidth = 0.5) +
                annotate("text", x = IP_u + 0.01, y = 0.04,
                         label = sprintf("IP_u=%.3f", IP_u),
                         colour = IP_U_COL_10E, hjust = 0, size = 2.5, fontface = "bold") +
                # IP_v label
                (if (is.finite(IP_v))
                     annotate("text", x = 0.015, y = IP_v + 0.025,
                              label = sprintf("IP_v=%.3f", IP_v),
                              colour = IP_V_COL_10E, hjust = 0, size = 2.5,
                              fontface = "bold")
                 else NULL) +
                coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
                labs(title = "A: CMRF: spine + adaptive delta band in (u,v)",
                     subtitle = paste0(
                         "Steel band = adaptive delta_adapt(u)=delta_base*(1+0.5*K(u)). ",
                         "Solid navy = spine g(u). Dotted = PAVA step (diag).\n",
                         "Teal solid = IP_u (band-kernel centre). Dark red ",
                         "dashed = IP_v, the boundary the Pre/Post colours ",
                         "are assigned on (cDFG at IP_y).\n", ip_coords_line),
                     x = "u = Fx(GOI)", y = "v = Fc(cDFG)") +
                theme_copula(base_size = 9) +
                theme(aspect.ratio = 1, panel.grid.minor = element_blank(),
                      legend.position = "bottom",
                      legend.key.size = unit(0.35, "cm"),
                      legend.text = element_text(size = 7))

            # Panel B: |v - g(u)| residuals vs u, adaptive delta_adapt(u) line
            # CLASS_COLS_F2 retained here for Panel C (boxplot by cell_class).
            CLASS_COLS_F2 <- c(
                "Original + Refined" = "#009E73",
                "Original only"      = "#8B6914",
                "Recovered"          = "#7C3AED",
                "Excluded"           = "gray65"
            )
            # Colour by the four positional groups; other cells are Non-selected.
            RESID_COLS_F2 <- c(
                "Non-selected"     = "gray65",
                "Pre-IP TREP"      = "#E69F00",
                "Pre-IP non-TREP"  = "#1B3A8C",
                "Post-IP TREP"     = "#009E73",
                "Post-IP non-TREP" = "#D94F4F"
            )
            # Orig+Ref cells for residual diagnostic (reference population)
            df_resid_plot <- df_refined[df_refined$cell_class == "Original + Refined", ]
            df_resid_plot <- df_resid_plot[order(df_resid_plot$u_all), ]

            # Assign plot_colour_b: 4 positional groups for "Original + Refined" cells.
            df_resid_plot$plot_colour_b <- dplyr::case_when(
                df_resid_plot$group_orig == "Pre-IP TREP"      ~ "Pre-IP TREP",
                df_resid_plot$group_orig == "Pre-IP non-TREP"  ~ "Pre-IP non-TREP",
                df_resid_plot$group_orig == "Post-IP TREP"     ~ "Post-IP TREP",
                df_resid_plot$group_orig == "Post-IP non-TREP" ~ "Post-IP non-TREP",
                TRUE                                            ~ "Non-selected"
            )
            df_resid_plot$plot_colour_b <- factor(
                df_resid_plot$plot_colour_b,
                levels = names(RESID_COLS_F2)
            )

            # Background gray layer: all non-Orig+Ref cells (Excluded, Recovered,
            # Original only) plotted as gray points behind the coloured Orig+Ref layer.
            df_resid_bg <- df_refined[df_refined$cell_class != "Original + Refined", ]
            df_resid_bg <- df_resid_bg[order(df_resid_bg$u_all), ]
            df_resid_bg$plot_colour_b <- factor("Non-selected",
                                                levels = names(RESID_COLS_F2))

            # Per-cell adaptive tolerance (recomputed on Orig+Ref subset for counts)
            delta_adapt_ref <- delta * (1 + 0.50 *
                exp(-(df_resid_plot$u_all - IP_u)^2 / (2 * 0.10^2)))

            # Three-zone counts (among Orig+Ref cells):
            #   Zone 1: within delta_base (flat floor)
            #   Zone 2: within adaptive band but outside delta_base
            #   Zone 3: outside both bands
            n_within_base    <- sum(df_resid_plot$resid_abs <= delta)
            n_adapt_only     <- sum(df_resid_plot$resid_abs >  delta &
                                    df_resid_plot$resid_abs <= delta_adapt_ref)
            n_outside_both   <- sum(df_resid_plot$resid_abs >  delta_adapt_ref)
            n_origref_total  <- nrow(df_resid_plot)

            # cDFG is not an axis here, so the Pre/Post boundary is described in words.
            subtitle_b <- paste0(
                "Solid navy = adaptive delta_adapt(u). Dashed = delta_base. Below curve = pass Stage 1. ",
                "Orig+Ref cells coloured (n=", n_origref_total,
                "); gray background = non-Orig+Ref (n=", nrow(df_resid_bg), ").\n",
                "Teal vertical = IP_u (band-kernel centre), NOT the Pre/Post ",
                "boundary: those colours are assigned on cDFG at IP_y, which ",
                "is not an axis of this panel.\n",
                sprintf("Within delta_base: n=%d | Adaptive band only: n=%d | Outside both: n=%d",
                        n_within_base, n_adapt_only, n_outside_both)
            )

            # Adaptive delta_adapt(u) curve for Panel B
            df_tau_line <- data.frame(
                u           = u_diag_grid,
                delta_adapt = delta * (1 + 0.50 * exp(-(u_diag_grid - IP_u)^2 / (2 * 0.10^2)))
            )

            # Split into gray background layer and colored foreground layer so the
            # 4-group cells are rendered on top (same layering approach as Panel A).
            df_resid_gray    <- df_resid_plot[df_resid_plot$plot_colour_b == "Non-selected", ]
            df_resid_colored <- df_resid_plot[df_resid_plot$plot_colour_b != "Non-selected", ]

            p_f2 <- ggplot() +
                # Gray background layer: non-Orig+Ref cells (Excluded, Recovered,
                # Original only) rendered first so coloured cells sit on top.
                geom_point(data = df_resid_bg,
                           aes(x = u_all, y = resid_abs),
                           colour = "gray80", size = 0.3, alpha = 0.25) +
                # Gray foreground layer (if any Orig+Ref cells have no positional group)
                geom_point(data = df_resid_gray,
                           aes(x = u_all, y = resid_abs, colour = plot_colour_b),
                           size = 0.4, alpha = 0.30) +
                # 4-group colored foreground layer (Original + Refined cells only)
                geom_point(data = df_resid_colored,
                           aes(x = u_all, y = resid_abs, colour = plot_colour_b),
                           size = 0.4, alpha = 0.45) +
                # Adaptive delta_adapt(u) curve (navy solid)
                geom_line(data = df_tau_line, aes(x = u, y = delta_adapt),
                          colour = "#1A237E", linewidth = 0.9,
                          linetype = "solid", inherit.aes = FALSE) +
                # delta_base flat reference (navy dashed)
                geom_hline(yintercept = delta, colour = "#1A237E",
                           linewidth = 0.5, linetype = "dashed") +
                geom_vline(xintercept = IP_u, colour = IP_U_COL_10E,
                           linetype = "solid", linewidth = 0.5) +
                annotate("text", x = 0.02, y = delta * 1.08,
                         label = sprintf("delta_base=%.5f", delta),
                         colour = "#1A237E", hjust = 0, vjust = 0, size = 2.4) +
                annotate("text", x = IP_u + 0.01,
                         y = delta * (1 + 0.50) * 1.05,
                         label = sprintf("delta_max=%.5f\n(at IP_u)", delta * 1.50),
                         colour = "#1A237E", hjust = 0, vjust = 0, size = 2.4) +
                scale_colour_manual(values = RESID_COLS_F2, name = "Cell group",
                                    drop = FALSE,
                                    guide = guide_legend(
                                        override.aes = list(
                                            size  = c(1.8, rep(2.2, 4)),
                                            shape = rep(16L, 5),
                                            alpha = c(0.40, rep(0.75, 4))
                                        )
                                    )) +
                labs(title = "B: |v - g(u)| vs u  [Stage 1: GAM residuals]",
                     subtitle = subtitle_b,
                     x = "u = Fx(GOI)", y = "|v - g(u)|") +
                theme_copula(base_size = 9) +
                theme(panel.grid.minor = element_blank(),
                      legend.position = "bottom",
                      legend.key.size = unit(0.35, "cm"),
                      legend.text = element_text(size = 7),
                      plot.subtitle = element_text(size = 6.5, lineheight = 1.3))

            # Panel C: per-class residual boxplot
            # Only include classes that are actually present in df_refined
            present_classes <- intersect(
                c("Original + Refined", "Original only", "Recovered", "Excluded"),
                unique(df_refined$cell_class)
            )
            df_box <- df_refined[df_refined$cell_class %in% present_classes, ]
            df_box$cell_class_f <- factor(df_box$cell_class, levels = present_classes)

            # delta_max = delta_adapt at IP_u (the actual admission criterion peak)
            delta_max_val <- delta * (1 + 0.50)  # kappa=0.50 at IP_u
            p_f3 <- ggplot(df_box, aes(x = cell_class_f, y = resid_abs,
                                       fill = cell_class_f)) +
                geom_boxplot(outlier.size = 0.4, outlier.alpha = 0.3,
                             linewidth = 0.4) +
                # delta_base (dashed navy) -- floor tolerance
                geom_hline(yintercept = delta, colour = "#1A237E",
                           linewidth = 0.7, linetype = "dashed") +
                # delta_max at IP_u (dotted navy) -- peak adaptive tolerance
                geom_hline(yintercept = delta_max_val, colour = "#1A237E",
                           linewidth = 0.5, linetype = "dotted") +
                annotate("text", x = 0.55, y = delta * 1.04,
                         label = sprintf("delta_base=%.4f", delta),
                         colour = "#1A237E", size = 2.3, hjust = 0, vjust = 0) +
                annotate("text", x = 0.55, y = delta_max_val * 1.04,
                         label = sprintf("delta_max=%.4f (at IP_u)", delta_max_val),
                         colour = "#1A237E", size = 2.3, hjust = 0, vjust = 0) +
                scale_fill_manual(values = CLASS_COLS_F2, guide = "none", drop = TRUE) +
                scale_x_discrete(labels = function(x) gsub(" \\(", "\n(", x)) +
                labs(title = "C: Residual distribution by cell class",
                     subtitle = "Navy dashed = delta_base; navy dotted = delta_max at IP_u (actual criterion). Boxes below delta_max pass Stage 1.",
                     x = NULL, y = "|residual|") +
                theme_copula(base_size = 9) +
                theme(panel.grid.minor = element_blank(),
                      axis.text.x = element_text(size = 7, lineheight = 0.85))

            # Panel D: Orig+Ref GAM residual distribution (delta_base calibration)
            # Shows the empirical basis for delta_base (95th percentile marked)
            df_grp1a_resid <- data.frame(
                resid = cmrf_delta_cal$resid_all,
                u     = u_all[is_orig_skcs]
            )
            p_f4 <- ggplot(df_grp1a_resid, aes(x = resid)) +
                geom_histogram(bins = 30, fill = "#64B5F6", alpha = 1.00,
                               colour = "#0D47A1", linewidth = 0.38) +
                geom_vline(xintercept = delta, colour = "#004D40",
                           linewidth = 0.9, linetype = "solid") +
                geom_vline(xintercept = cmrf_delta_cal$resid_median, colour = "#00695C",
                           linewidth = 0.9, linetype = "dashed") +
                annotate("text", x = delta * 1.02, y = Inf,
                         label = sprintf("delta_base=%.5f\n(95th pctile)", delta),
                         colour = "#004D40", vjust = 1.5, hjust = 0, size = 2.6) +
                annotate("text", x = cmrf_delta_cal$resid_median * 1.02, y = Inf,
                         label = sprintf("median=%.5f", cmrf_delta_cal$resid_median),
                         colour = "#00695C", vjust = 3.0, hjust = 0, size = 2.6) +
                labs(title = "D: Orig+Ref GAM residuals [delta_base calibration]",
                     subtitle = sprintf("n=%d Orig+Ref cells | delta_base=%.5f (95th pctile) | delta_max=%.5f at waist",
                                        cmrf_delta_cal$n_grp1a, delta, delta * 1.50),
                     x = "|v - g(u)|  [GAM residual]", y = "Count") +
                theme_copula(base_size = 9) +
                theme(panel.grid.minor = element_blank())

            p3f_grid <- cowplot::plot_grid(p_f1, p_f2, p_f3, p_f4,
                                           nrow = 2, ncol = 2)
            p01 <- cowplot::plot_grid(
                cowplot::ggdraw() + cowplot::draw_label(
                    sprintf("Part10e Plot 1 -- CMRF B+C diagnostics: %s  |  delta_base=%.5f (Orig+Ref 95th pctile)",
                            dataset_tag, delta),
                    fontface = "bold", size = 10, x = 0.01, hjust = 0),
                p3f_grid,
                ncol = 1, rel_heights = c(0.05, 1)
            )

            save_pp(p01,
                    file.path(OUT_DIR, paste0(FILE_PREFIX, "_", dataset_tag,
                                              "_01_CMRF_diagnostics")),
                    pdf_w = 13, pdf_h = 11, png_w = 2800, png_h = 2400)
            cat("  Plot 1 saved (CMRF diagnostics)\n")
            all_plots[["01_CMRF_diagnostics"]] <- p01

            # =====================================================================
            # PLOT 2 -- CMRF Summary Panel (Stage 2 gate, post-hoc, delta cal,
            #             RSk-fCV fold RMSE comparison)
            # =====================================================================
            cat("\n  STEP 3g -- CMRF summary panel (gate, post-hoc, calibration, RMSE)\n")

            s2 <- cmrf_res$stage2_res  # Stage 2 gate results

            # == Panel A: Stage 2 gate waterfall (narrower) ==================
            n_s1_cand <- if (!is.null(s2)) s2$n_candidates else 0L
            n_removed <- if (!is.null(s2)) s2$n_removed    else 0L
            n_accept  <- if (!is.null(s2)) sum(s2$pass_stage2) else 0L
            p_final   <- if (!is.null(s2)) s2$wilcox_p_final  else NA_real_
            gate_pass <- !is.na(p_final) && p_final >= CMRF_ALPHA

            df_gate <- data.frame(
                stage = factor(c("Stage 1\ncandidates", "Removed\nby RSk-fCV", "Recovered"),
                               levels = c("Stage 1\ncandidates", "Removed\nby RSk-fCV", "Recovered")),
                n     = c(n_s1_cand, n_removed, n_accept),
                fill  = c("candidate", "removed", "recovered")
            )
            p_g1 <- ggplot(df_gate, aes(x = stage, y = n, fill = fill)) +
                geom_col(width = 0.55, colour = "black", linewidth = 0.4) +
                geom_text(aes(label = n), vjust = -0.4, size = 3.5, fontface = "bold") +
                scale_fill_manual(values = c("candidate" = "grey80",
                                             "removed"   = "grey50",
                                             "recovered" = "#7C3AED"),
                                  guide = "none") +
                scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
                labs(title = "A: Stage 2 RSk-fCV gate",
                     subtitle = sprintf("Wilcoxon p = %.4f; alpha = %.2f [%s]",
                                        p_final, CMRF_ALPHA,
                                        ifelse(gate_pass, "PASS", "FAIL")),
                     x = NULL, y = "Cell count") +
                theme_copula(base_size = 9.5) +
                theme(panel.grid.minor = element_blank())

            # == Panel B: T2 post-hoc pairwise (SW-gated: KW+Dunn or ANOVA+Holm) ===
            # multi_grp() is shared with the verification section.
            t2_mg_result <- multi_grp(
                df_refined$resid_abs, df_refined$cell_class,
                levs = c("Original + Refined", "Recovered", "Original only", "Excluded")
            )
            t2_posthoc <- t2_mg_result$posthoc
            if (!is.null(t2_posthoc) && nrow(t2_posthoc) > 0) {
                if ("comparison" %in% names(t2_posthoc)) {
                    # Dunn's BH posthoc (non-normal path)
                    df_dunn_plot <- data.frame(
                        comparison = gsub(" - ", "\nvs\n", t2_posthoc$comparison),
                        p_bh       = t2_posthoc$p_BH,
                        sig        = ifelse(t2_posthoc$p_BH < 0.05, "p < 0.05", "ns"),
                        stringsAsFactors = FALSE
                    )
                } else {
                    # Holm-Sidak posthoc (normal path)
                    df_dunn_plot <- data.frame(
                        comparison = paste0(t2_posthoc$group1, "\nvs\n", t2_posthoc$group2),
                        p_bh       = t2_posthoc$p_holm_sidak,
                        sig        = ifelse(t2_posthoc$p_holm_sidak < 0.05, "p < 0.05", "ns"),
                        stringsAsFactors = FALSE
                    )
                }
            } else {
                df_dunn_plot <- data.frame(
                    comparison = "n/a", p_bh = NA_real_, sig = "ns",
                    stringsAsFactors = FALSE
                )
            }
            df_dunn_plot$comparison <- factor(df_dunn_plot$comparison,
                                              levels = df_dunn_plot$comparison)

            p_g2 <- ggplot(df_dunn_plot, aes(x = comparison, y = -log10(p_bh), fill = sig)) +
                geom_col(width = 0.55, colour = "black", linewidth = 0.4) +
                geom_hline(yintercept = -log10(0.05), linetype = "dashed",
                           colour = "#1A237E", linewidth = 0.6) +
                annotate("text", x = Inf, y = -log10(0.05) + 0.08,
                         label = "FDR = 0.05", colour = "#1A237E",
                         hjust = 1.1, size = 2.8, fontface = "italic") +
                geom_text(aes(label = sprintf("%.4f", p_bh)),
                          vjust = -0.3, size = 2.6, fontface = "bold") +
                scale_fill_manual(values = c("p < 0.05" = "#009E73", "ns" = "#8A9BAD"),
                                  guide = "none") +
                labs(title = "B: T2 post-hoc pairwise (Dunn's BH)",
                     subtitle = "Bars above dashed line (FDR<0.05): class pair has significantly different spine residuals,\ni.e. CMRF assigns them distinct fidelity levels. Below: residuals statistically similar.",
                     x = NULL, y = expression(-log[10](p[BH]))) +
                theme_copula(base_size = 9.5) +
                theme(panel.grid.minor = element_blank(),
                      axis.text.x = element_text(size = 6.5, lineheight = 0.85))

            # == Panel C: Delta calibration summary (table-grob) =============
            cal_df <- data.frame(
                Parameter = c("delta_base", "delta_max (at waist)",
                              "Coverage (quantile)", "n_Orig+Ref (total)",
                              "Orig+Ref resid median",
                              "Orig+Ref resid max", "kappa", "bw (Gaussian)"),
                Value = c(sprintf("%.5f", delta),
                          sprintf("%.5f", delta * (1 + CMRF_KAPPA)),
                          sprintf("%.0f%%", CMRF_COVERAGE * 100),
                          as.character(cmrf_delta_cal$n_grp1a),
                          sprintf("%.5f", cmrf_delta_cal$resid_median),
                          sprintf("%.5f", cmrf_delta_cal$resid_max),
                          as.character(CMRF_KAPPA),
                          as.character(CMRF_BW)),
                stringsAsFactors = FALSE
            )

            p_g3 <- ggplot(cal_df, aes(x = 1, y = rev(seq_len(nrow(cal_df))))) +
                geom_text(aes(x = 0.7, label = Parameter), hjust = 1, size = 3.0,
                          fontface = "bold", colour = "grey20") +
                geom_text(aes(x = 0.75, label = Value), hjust = 0, size = 3.0,
                          colour = "#1A237E") +
                xlim(0, 1.5) +
                labs(title = "C: Delta calibration parameters",
                     subtitle = sprintf("CMRF tolerance calibration: %s", dataset_tag)) +
                theme_void(base_size = 9.5) +
                theme(plot.title = element_text(face = "bold", size = 10,
                                                margin = margin(b = 2)),
                      plot.subtitle = element_text(size = 8, colour = "grey30",
                                                   margin = margin(b = 5)),
                      plot.background = element_rect(fill = "white", colour = NA),
                      plot.margin = margin(8, 8, 8, 8))

            # == Panel D: Fold RMSE boxplot (Orig+Ref vs Orig+Ref + Recovered) ====
            if (!is.null(s2) && length(s2$rmse_grp1a) > 0) {
                df_rmse_03g <- data.frame(
                    RMSE  = c(s2$rmse_grp1a, s2$rmse_final_union),
                    Group = factor(rep(c("Orig+Ref\nalone",
                                         "Orig+Ref +\nRecovered"),
                                       each = length(s2$rmse_grp1a)),
                                   levels = c("Orig+Ref\nalone",
                                              "Orig+Ref +\nRecovered"))
                )
                p_g4 <- ggplot(df_rmse_03g, aes(x = Group, y = RMSE, fill = Group)) +
                    geom_boxplot(width = 0.45, colour = "black", linewidth = 0.4,
                                 outlier.shape = NA, alpha = 0.35) +
                    geom_jitter(aes(colour = Group), width = 0.15, size = 0.8,
                                alpha = 0.25, shape = 21, stroke = 0.3,
                                show.legend = FALSE) +
                    scale_fill_manual(values = c("Orig+Ref\nalone"     = "#009E73",
                                                 "Orig+Ref +\nRecovered" = "#7C3AED"),
                                      guide = "none") +
                    scale_colour_manual(values = c("Orig+Ref\nalone"     = "#009E73",
                                                   "Orig+Ref +\nRecovered" = "#7C3AED"),
                                        guide = "none") +
                    labs(title = "D: RSk-fCV fold RMSE",
                         subtitle = sprintf("p = %.4f (n.s.): adding Recovered does not degrade monotone fit.\nOverlapping boxes = union RMSE indistinguishable from Orig+Ref alone (%d paired folds).",
                                            p_final, length(s2$rmse_grp1a)),
                         x = NULL, y = "Isotonic RMSE per fold") +
                    theme_copula(base_size = 9.5) +
                    theme(panel.grid.minor = element_blank())
            } else {
                p_g4 <- ggplot() +
                    annotate("text", x = 0.5, y = 0.5,
                             label = "No Stage 2 candidates", size = 5) +
                    theme_void()
            }

            # Combine: A narrower (~0.7x), BCD equal
            p3g_left  <- cowplot::plot_grid(p_g1, p_g3, nrow = 2, rel_heights = c(1, 1))
            p3g_right <- cowplot::plot_grid(p_g2, p_g4, nrow = 2, rel_heights = c(1, 1))
            p3g_grid  <- cowplot::plot_grid(p3g_left, p3g_right,
                                            ncol = 2, rel_widths = c(0.7, 1))
            p02 <- cowplot::plot_grid(
                cowplot::ggdraw() + cowplot::draw_label(
                    sprintf("Part10e Plot 2 -- CMRF Summary Panel: %s", dataset_tag),
                    fontface = "bold", size = 11, x = 0.01, hjust = 0),
                p3g_grid,
                ncol = 1, rel_heights = c(0.05, 1)
            )

            save_pp(p02,
                    file.path(OUT_DIR, paste0(FILE_PREFIX, "_", dataset_tag,
                                              "_02_CMRF_summary_panel")),
                    pdf_w = 14, pdf_h = 11, png_w = 3000, png_h = 2400)
            cat("  Plot 2 saved (CMRF summary panel)\n")
            all_plots[["02_CMRF_summary_panel"]] <- p02


            # =====================================================================
            # PLOT 3 -- Before vs After CMRF in (u,v) copula space
            # Panel A: BEFORE -- 4-group SKCS colours (like Plot 6)
            # Panel B: AFTER  -- Recovered / Removed / Orig-only
            # =====================================================================
            cat("\n  STEP 3h -- Before vs After CMRF in (u,v) copula space\n")

            idx_recov   <- which(df_refined$is_recovered)
            idx_removed <- which(df_refined$cell_class == "Excluded" & removed_s2)

            # == Panel A: BEFORE CMRF (4-group SKCS colours) =================
            ZOOM_COLS_03H <- c("Non-selected" = "gray60", FOUR_COLS)
            df_03h_all <- df_refined
            df_03h_all$plot_colour <- dplyr::case_when(
                df_03h_all$is_orig_skcs &
                    df_03h_all$group_orig == "Pre-IP TREP"      ~ "Pre-IP TREP",
                df_03h_all$is_orig_skcs &
                    df_03h_all$group_orig == "Pre-IP non-TREP"  ~ "Pre-IP non-TREP",
                df_03h_all$is_orig_skcs &
                    df_03h_all$group_orig == "Post-IP TREP"     ~ "Post-IP TREP",
                df_03h_all$is_orig_skcs &
                    df_03h_all$group_orig == "Post-IP non-TREP" ~ "Post-IP non-TREP",
                TRUE ~ "Non-selected"
            )
            df_03h_nonsel <- df_03h_all[df_03h_all$plot_colour == "Non-selected", ]
            df_03h_skcs   <- df_03h_all[df_03h_all$plot_colour != "Non-selected", ]

            p_h1 <- ggplot() +
                geom_ribbon(data = df_band_diag,
                            aes(x = u, ymin = lower, ymax = upper),
                            fill = "#607D8B", alpha = 0.15, inherit.aes = FALSE) +
                geom_line(data = df_band_diag, aes(x = u, y = g),
                          colour = SPINE_COL_10E, linewidth = 0.7) +
                geom_point(data = df_03h_nonsel,
                           aes(x = u_all, y = v_all),
                           colour = "gray75", size = 0.35, alpha = 0.40) +
                geom_point(data = df_03h_skcs,
                           aes(x = u_all, y = v_all, colour = plot_colour),
                           size = 0.9, alpha = 0.55) +
                scale_colour_manual(values = ZOOM_COLS_03H, name = "Cell group",
                                    drop = FALSE) +
                geom_vline(xintercept = IP_u, colour = IP_U_COL_10E,
                           linetype = "solid", linewidth = 0.5) +
                ip_boundary_v(linewidth = 0.5) +
                geom_hline(yintercept = 0.5, colour = "#1A237E",
                           linetype = "dashed", linewidth = 0.4) +
                coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
                labs(title = "A: BEFORE CMRF -- All SKCS-selected cells (4-group)",
                     subtitle = paste0(
                         "All original SKCS core cells (Si_norm>0.05, Ki>0.5) by ",
                         "Pre/Post-IP and TREP/non-TREP.\nGrey = non-selected ",
                         "cells. Navy = spine g(u); band = adaptive delta.\n",
                         "Teal solid = IP_u (band-kernel centre). Dark red ",
                         "dashed = IP_v, the boundary the Pre/Post colours ",
                         "are assigned on (cDFG at IP_y).\n",
                         ip_coords_line),
                     x = "u = Fx(GOI)", y = "v = Fc(cDFG)") +
                theme_copula(base_size = 9.5) +
                theme(panel.grid.minor = element_blank(), aspect.ratio = 1,
                      legend.position = "bottom", legend.key.size = unit(0.3, "cm"),
                      legend.text = element_text(size = 7))

            # == Panel B: AFTER CMRF (Recovered / Removed / Orig-only) =======
            df_h_origonly <- df_refined[df_refined$cell_class == "Original only", ]
            df_h_recov   <- df_refined[idx_recov, ]
            df_h_removed <- if (length(idx_removed) > 0) df_refined[idx_removed, ] else NULL
            df_h_bg <- df_refined[!df_refined$is_recovered & !removed_s2 &
                                      df_refined$cell_class != "Original only", ]

            p_h2 <- ggplot() +
                geom_ribbon(data = df_band_diag,
                            aes(x = u, ymin = lower, ymax = upper),
                            fill = "#607D8B", alpha = 0.15, inherit.aes = FALSE) +
                geom_line(data = df_band_diag, aes(x = u, y = g),
                          colour = SPINE_COL_10E, linewidth = 0.7) +
                geom_point(data = df_h_bg,
                           aes(x = u_all, y = v_all),
                           colour = "grey85", size = 0.3, alpha = 0.25) +
                {if (nrow(df_h_origonly) > 0)
                    geom_point(data = df_h_origonly,
                               aes(x = u_all, y = v_all),
                               colour = "white", fill = "#8B6914",
                               shape = 25, size = 1.8, alpha = 0.85, stroke = 0.4)
                } +
                geom_point(data = df_h_recov,
                           aes(x = u_all, y = v_all),
                           colour = "white", fill = "#7C3AED",
                           shape = 23, size = 1.8, alpha = 0.8, stroke = 0.5) +
                {if (!is.null(df_h_removed))
                    geom_point(data = df_h_removed,
                               aes(x = u_all, y = v_all),
                               colour = "grey40", fill = "grey50",
                               size = 2.0, alpha = 0.8, shape = 4, stroke = 0.7)
                } +
                geom_vline(xintercept = IP_u, colour = IP_U_COL_10E,
                           linetype = "solid", linewidth = 0.5) +
                geom_hline(yintercept = 0.5, colour = "#1A237E",
                           linetype = "dashed", linewidth = 0.4) +
                coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
                labs(title = "B: AFTER CMRF -- Recovered, Removed, Orig-only",
                     subtitle = sprintf("Purple <> = Recovered (n=%d); grey x = Removed by RSk-fCV (n=%d);\ntan v = Orig-only (n=%d). Spine + adaptive delta band shown.",
                                        length(idx_recov), length(idx_removed),
                                        nrow(df_h_origonly)),
                     x = "u = Fx(GOI)", y = "v = Fc(cDFG)") +
                theme_copula(base_size = 9.5) +
                theme(panel.grid.minor = element_blank(), aspect.ratio = 1)

            p3h_grid <- cowplot::plot_grid(p_h1, p_h2, nrow = 1, ncol = 2)
            p03 <- cowplot::plot_grid(
                cowplot::ggdraw() + cowplot::draw_label(
                    sprintf("Part10e Plot 3 -- Before vs After CMRF: %s", dataset_tag),
                    fontface = "bold", size = 11, x = 0.01, hjust = 0),
                p3h_grid,
                ncol = 1, rel_heights = c(0.05, 1)
            )

            save_pp(p03,
                    file.path(OUT_DIR, paste0(FILE_PREFIX, "_", dataset_tag,
                                              "_03_before_after_CMRF")),
                    pdf_w = 16, pdf_h = 8, png_w = 3400, png_h = 1800)
            cat("  Plot 3 saved (Before vs After CMRF)\n")
            all_plots[["03_before_after_CMRF"]] <- p03

            # =====================================================================
            # PLOT 4 -- RSk-fCV Stage 2 gate: mechanism and evidence
            # 4 panels explaining how the repeated stratified k-fold CV works.
            # All statistical tests follow SW-gated procedure with BH FDR correction.
            # =====================================================================
            cat("\n  STEP 3i -- RSk-fCV mechanism and evidence\n")

            if (!is.null(s2) && length(s2$rmse_grp1a) > 0) {

                n_folds_total <- length(s2$rmse_grp1a)
                rmse_diff_all <- s2$rmse_final_union - s2$rmse_grp1a

                # == Panel A: Recovered vs Removed (SW-gated test + OVL) =====
                if (length(idx_removed) > 0 && length(idx_recov) > 0) {
                    rv_recov <- df_refined$resid_abs[idx_recov]
                    rv_remov <- df_refined$resid_abs[idx_removed]
                    # OVL via kernel density
                    dens_r <- density(rv_recov, from = 0,
                                      to = max(c(rv_recov, rv_remov)) * 1.1, n = 512)
                    dens_m <- density(rv_remov, from = 0,
                                      to = max(c(rv_recov, rv_remov)) * 1.1, n = 512)
                    ovl_min <- pmin(approx(dens_r$x, dens_r$y, dens_m$x)$y, dens_m$y)
                    ovl_min[is.na(ovl_min)] <- 0
                    ovl_val <- min(max(sum(ovl_min * diff(dens_m$x[1:2])), 0), 1)
                    # SW-gated 2-group test
                    sw_rv <- shapiro.test(rv_recov)$p.value
                    sw_rm <- shapiro.test(rv_remov)$p.value
                    both_normal <- sw_rv >= 0.05 & sw_rm >= 0.05
                    if (both_normal) {
                        tt <- t.test(rv_recov, rv_remov, alternative = "less")
                        test_name <- "t-test (less)"
                        test_p <- tt$p.value
                    } else {
                        wt <- wilcox.test(rv_recov, rv_remov, alternative = "less")
                        test_name <- "Mann-Whitney U (less)"
                        test_p <- wt$p.value
                    }
                    df_rv <- data.frame(
                        resid_abs = c(rv_recov, rv_remov),
                        group     = factor(
                            c(rep("Recovered", length(rv_recov)),
                              rep("Removed by\nRSk-fCV", length(rv_remov))),
                            levels = c("Recovered", "Removed by\nRSk-fCV"))
                    )
                    p_i1 <- ggplot(df_rv, aes(x = group, y = resid_abs, fill = group)) +
                        geom_boxplot(width = 0.45, colour = "black", linewidth = 0.4,
                                     outlier.shape = NA, alpha = 0.35) +
                        geom_jitter(aes(colour = group), width = 0.15, size = 1.0,
                                    alpha = 0.35, shape = 21, stroke = 0.4,
                                    show.legend = FALSE) +
                        geom_hline(yintercept = delta, colour = "#1A237E",
                                   linetype = "dashed", linewidth = 0.5) +
                        scale_fill_manual(values = c("Recovered" = "#7C3AED",
                                                     "Removed by\nRSk-fCV" = "grey50"),
                                          guide = "none") +
                        scale_colour_manual(values = c("Recovered" = "#7C3AED",
                                                        "Removed by\nRSk-fCV" = "grey50"),
                                            guide = "none") +
                        labs(title = "A: Spine residual: Recovered vs Removed",
                             subtitle = sprintf("SW-gated %s p = %.2e.\nOVL = %.3f; low OVL = near-complete separation. Boundary overlap is expected\n(RSk-fCV decides the cutoff, not OVL).",
                                                test_name, test_p, ovl_val),
                             x = NULL, y = "|v - g(u)|") +
                        theme_copula(base_size = 9.5) +
                        theme(panel.grid.minor = element_blank())
                } else {
                    p_i1 <- ggplot() +
                        annotate("text", x = 0.5, y = 0.5,
                                 label = "No removed cells", size = 4) +
                        theme_void()
                }

                # == Panel B: Paired fold scatter (Orig+Ref RMSE vs Union RMSE) ====
                df_pair <- data.frame(origref = s2$rmse_grp1a, union = s2$rmse_final_union)
                # Paired Wilcoxon (matches gate's own test)
                wsr <- suppressWarnings(
                    wilcox.test(s2$rmse_final_union, s2$rmse_grp1a,
                                paired = TRUE, alternative = "greater", exact = FALSE))
                pair_p <- wsr$p.value
                p_i2 <- ggplot(df_pair, aes(x = origref, y = union)) +
                    geom_abline(slope = 1, intercept = 0, colour = "#CC7000",
                                linetype = "dashed", linewidth = 0.6) +
                    geom_point(size = 1.2, alpha = 0.45, colour = "#1A237E",
                               shape = 16) +
                    coord_fixed() +
                    labs(title = "B: Paired fold RMSE (test-set evaluation)",
                         subtitle = sprintf("Each dot = one test-fold's isotonic RMSE (train on 4/5, evaluate on 1/5).\nWilcoxon signed-rank (greater) p = %.4f. Dots near orange identity line = no degradation.",
                                            pair_p),
                         x = "Orig+Ref alone RMSE", y = "Orig+Ref + Recovered RMSE") +
                    theme_copula(base_size = 9.5) +
                    theme(panel.grid.minor = element_blank())

                # == Panel C: All Stage 1 candidates ranked by PAVA residual ====
                non_skcs_s1_idx <- which(!df_refined$is_orig_skcs &
                                             cmrf_res$pass_stage1)
                if (length(s2$removal_order) > 0 &&
                    length(non_skcs_s1_idx) >= max(c(s2$removal_order, 0))) {
                    removed_global <- non_skcs_s1_idx[s2$removal_order]
                    accepted_idx   <- non_skcs_s1_idx[s2$pass_stage2]
                    # PAVA residuals: removed in removal order, accepted sorted desc
                    resid_pava_removed  <- s2$resid_pava_cand[s2$removal_order]
                    resid_pava_accepted <- s2$resid_pava_cand[s2$pass_stage2]
                    df_bars <- data.frame(
                        step      = seq_len(length(removed_global) + length(accepted_idx)),
                        resid_pava = c(resid_pava_removed,
                                       sort(resid_pava_accepted, decreasing = TRUE)),
                        type      = factor(c(rep("Removed", length(removed_global)),
                                             rep("Recovered", length(accepted_idx))),
                                           levels = c("Removed", "Recovered"))
                    )
                    p_i3 <- ggplot(df_bars, aes(x = step, y = resid_pava, fill = type)) +
                        geom_col(colour = NA, width = 0.9) +
                        scale_fill_manual(values = c("Removed" = "grey50",
                                                     "Recovered" = "#7C3AED"),
                                          name = NULL) +
                        labs(title = "C: All Stage 1 candidates ranked by PAVA residual",
                             subtitle = "Left grey = removed worst-first by PAVA isotonic residual.\nRight purple = accepted as Recovered. Gate stops when Wilcoxon p >= alpha.",
                             x = "Candidate rank (descending PAVA residual)",
                             y = "|v - PAVA(u)|") +
                        theme_copula(base_size = 9.5) +
                        theme(panel.grid.minor = element_blank(),
                              legend.position = "bottom",
                              legend.key.size = unit(0.3, "cm"))
                } else {
                    p_i3 <- ggplot() + annotate("text", x=.5, y=.5,
                        label="No removal data", size=4) + theme_void()
                }

                # == Panel D: Cumulative mean RMSE difference ================
                df_cum <- data.frame(
                    fold     = seq_along(rmse_diff_all),
                    diff     = rmse_diff_all,
                    cum_mean = cumsum(rmse_diff_all) / seq_along(rmse_diff_all)
                )
                p_i4 <- ggplot(df_cum) +
                    geom_hline(yintercept = 0, colour = "#8B4513",
                               linetype = "solid", linewidth = 0.5) +
                    geom_point(aes(x = fold, y = diff), colour = "grey70",
                               size = 0.6, alpha = 0.35) +
                    geom_line(aes(x = fold, y = cum_mean), colour = "#1A237E",
                              linewidth = 0.8) +
                    labs(title = "D: Cumulative mean dRMSE across folds",
                         subtitle = sprintf("Grey dots = individual fold RMSE differences (noisy per-fold values).\nNavy line = running mean, converging to dmu = %.6f ~ 0 (stable over %d folds).",
                                            mean(rmse_diff_all), n_folds_total),
                         x = "Fold (cumulative)", y = "Running mean dRMSE") +
                    theme_copula(base_size = 9.5) +
                    theme(panel.grid.minor = element_blank())

                p3i_grid <- cowplot::plot_grid(p_i1, p_i2, p_i3, p_i4,
                                               nrow = 2, ncol = 2)
                p04 <- cowplot::plot_grid(
                    cowplot::ggdraw() + cowplot::draw_label(
                        sprintf("Part10e Plot 4 -- RSk-fCV mechanism and evidence: %s",
                                dataset_tag),
                        fontface = "bold", size = 11, x = 0.01, hjust = 0),
                    p3i_grid,
                    ncol = 1, rel_heights = c(0.05, 1)
                )

                save_pp(p04,
                        file.path(OUT_DIR, paste0(FILE_PREFIX, "_", dataset_tag,
                                                  "_04_RSk_fCV_mechanism")),
                        pdf_w = 14, pdf_h = 12, png_w = 3000, png_h = 2600)
                cat("  Plot 4 saved (RSk-fCV mechanism and evidence)\n")
                all_plots[["04_RSk_fCV_mechanism"]] <- p04


                # =============================================================
                # RSk-fCV AUDIT EXCEL (8 sheets per dataset)
                # Full intermediate data for reproducing Wilcoxon from scratch.
                # =============================================================
                cat("\n  RSk-fCV audit Excel export...\n")

                audit_sheets <- list()
                n_reps_a  <- CMRF_N_REPEATS
                n_k_a     <- CMRF_K_FOLDS

                # == Sheet 1: Cell_Metadata ==================================
                safe_col <- function(df, col) {
                    if (col %in% names(df)) df[[col]] else rep(NA, nrow(df))
                }
                n_rows <- nrow(df_refined)
                sh1 <- data.frame(
                    barcode      = df_refined$barcode,
                    dataset      = rep(dataset_tag, n_rows),
                    u            = safe_col(df_refined, "u_all"),
                    v            = safe_col(df_refined, "v_all"),
                    is_skcs      = safe_col(df_refined, "is_orig_skcs"),
                    group_orig   = safe_col(df_refined, "group_orig"),
                    cell_class   = safe_col(df_refined, "cell_class"),
                    is_recovered = safe_col(df_refined, "is_recovered"),
                    g_u          = safe_col(df_refined, "g_u"),
                    resid_abs    = safe_col(df_refined, "resid_abs"),
                    region       = safe_col(df_refined, "region"),
                    GOI_expr     = safe_col(df_refined, "GOI_expr"),
                    cDFG_expr    = safe_col(df_refined, "cDFG_expr"),
                    spear_contrib = safe_col(df_refined, "spear_contrib"),
                    kend_contrib  = safe_col(df_refined, "kend_contrib"),
                    is_GOI_neg   = safe_col(df_refined, "is_GOI_neg"),
                    is_trep      = safe_col(df_refined, "is_trep"),
                    stringsAsFactors = FALSE
                )
                # Add delta_adapt and decile
                Kw_all <- exp(-(df_refined$u_all - IP_u)^2 / (2 * CMRF_BW^2))
                sh1$delta_adapt_u <- delta * (1 + CMRF_KAPPA * Kw_all)
                sh1$stage1_band_pass  <- if (length(cmrf_res$cond2) == n_rows) cmrf_res$cond2 else rep(NA, n_rows)
                u_or_audit <- df_refined$u_all[df_refined$is_orig_skcs]
                brks_audit <- quantile(u_or_audit, probs = seq(0, 1, by = 0.1))
                brks_audit[1] <- 0; brks_audit[length(brks_audit)] <- 1
                # Tied u quantiles: fall back to equal-width deciles.
                if (anyDuplicated(brks_audit) > 0L) {
                    cat("  NOTE: tied u quantiles in the audit deciles;",
                        " using equal-width breaks instead\n", sep = "")
                    brks_audit <- seq(0, 1, length.out = 11L)
                }
                sh1$u_quantile_decile <- as.character(cut(df_refined$u_all,
                    breaks = brks_audit, labels = paste0("D", 1:10),
                    include.lowest = TRUE))
                audit_sheets[["Cell_Metadata"]] <- sh1

                # == Sheet 1c: Stage 1: spine band proximity audit ===========
                u_vec      <- unname(df_refined$u_all)
                v_vec      <- unname(df_refined$v_all)
                g_u_vec    <- unname(safe_col(df_refined, "g_u"))
                resid_vec  <- abs(v_vec - g_u_vec)
                dadapt_vec <- unname(sh1$delta_adapt_u)
                is_skcs_vec <- unname(df_refined$is_orig_skcs)
                sh_c2 <- data.frame(
                    barcode                    = unname(df_refined$barcode),
                    u                          = u_vec,
                    v                          = v_vec,
                    delta_base                 = delta,
                    delta_adapt_u              = dadapt_vec,
                    abs_v_minus_g_u            = resid_vec,
                    resid_within_delta         = resid_vec <= dadapt_vec,
                    pass_stage1                = resid_vec <= dadapt_vec,
                    is_orig_skcs               = is_skcs_vec,
                    stringsAsFactors           = FALSE
                )
                audit_sheets[["Stage1_Band"]] <- sh_c2

                # == Sheet 2: Fold_RMSE (200 rows) ===========================
                rep_id_a  <- rep(seq_len(n_reps_a), each = n_k_a)
                fold_id_a <- rep(seq_len(n_k_a), times = n_reps_a)
                sh2 <- data.frame(
                    repeat_id     = rep_id_a,
                    fold_id       = fold_id_a,
                    fold_global   = seq_along(s2$rmse_grp1a),
                    RMSE_A        = s2$rmse_grp1a,
                    RMSE_B        = s2$rmse_final_union,
                    delta_RMSE    = s2$rmse_final_union - s2$rmse_grp1a
                )
                audit_sheets[["Fold_RMSE"]] <- sh2

                # == Sheet 3: Wilcoxon_Summary ===============================
                diffs_a <- s2$rmse_final_union - s2$rmse_grp1a
                wsr_a <- suppressWarnings(
                    wilcox.test(s2$rmse_final_union, s2$rmse_grp1a,
                                paired = TRUE, alternative = "greater", exact = FALSE))
                sh3 <- data.frame(
                    test_name      = "Wilcoxon signed-rank (greater)",
                    statistic      = as.numeric(wsr_a$statistic),
                    p_value        = wsr_a$p.value,
                    n_pairs        = length(diffs_a),
                    alpha          = CMRF_ALPHA,
                    decision       = ifelse(wsr_a$p.value >= CMRF_ALPHA, "PASS (accept)", "FAIL (remove)"),
                    mean_delta     = mean(diffs_a),
                    median_delta   = median(diffs_a)
                )
                audit_sheets[["Wilcoxon_Summary"]] <- sh3

                # == Sheet 4: Gate_Iteration =================================
                non_skcs_s1_audit <- which(!df_refined$is_orig_skcs &
                                               cmrf_res$pass_stage1)
                if (length(s2$removal_order) > 0 &&
                    length(non_skcs_s1_audit) >= max(s2$removal_order)) {
                    removed_glob <- non_skcs_s1_audit[s2$removal_order]
                    sh4 <- data.frame(
                        step           = seq_along(removed_glob),
                        barcode        = df_refined$barcode[removed_glob],
                        u              = df_refined$u_all[removed_glob],
                        v              = df_refined$v_all[removed_glob],
                        resid_abs_spine = df_refined$resid_abs[removed_glob],
                        resid_pava     = s2$resid_pava_cand[s2$removal_order],
                        n_remaining    = s2$n_candidates - seq_along(removed_glob)
                    )
                } else {
                    sh4 <- data.frame(step    = NA_integer_,
                                      barcode = NA_character_,
                                      note    = "No cells removed by the gate",
                                      stringsAsFactors = FALSE)
                }
                audit_sheets[["Gate_Iteration"]] <- sh4

                # == Sheet: Gate_PValue_Trajectory ==================================
                # Per-iteration p-values; iter = number of removals done.
                if (length(s2$p_per_iter) > 0) {
                    sh_pv <- data.frame(
                        iter     = seq_along(s2$p_per_iter) - 1L,
                        p_value  = s2$p_per_iter,
                        n_active = s2$n_per_iter
                    )
                    audit_sheets[["Gate_PValue_Trajectory"]] <- sh_pv
                }

                # == Sheet 5: One_Fold_Audit (fold 1 of repeat 1) ============
                # Recreate one fold to capture per-cell predictions
                set.seed(42)
                anchor_mask <- df_refined$is_orig_skcs & df_refined$pass_stage1
                union_mask <- anchor_mask | df_refined$is_recovered
                u_union <- df_refined$u_all[union_mask]
                v_union <- df_refined$v_all[union_mask]
                bc_union <- df_refined$barcode[union_mask]
                is_origref_union <- anchor_mask[union_mask]
                # Stratified fold assignment
                dec_audit <- cut(u_union, breaks = brks_audit,
                                 labels = 1:10, include.lowest = TRUE)
                fold_assign_audit <- integer(length(u_union))
                for (d in 1:10) {
                    idx_d <- which(dec_audit == d)
                    .fa <- rep(seq_len(n_k_a), length.out = length(idx_d))
                    fold_assign_audit[idx_d] <- .fa[sample.int(length(.fa))]
                }
                # Fold 1 as test
                test_mask_a <- fold_assign_audit == 1
                train_mask_a <- !test_mask_a
                # Curve A: Orig+Ref train only
                origref_train <- train_mask_a & is_origref_union
                iso_A <- isoreg(u_union[origref_train], v_union[origref_train])
                pred_A <- as.stepfun(iso_A)(u_union[test_mask_a])
                # Curve B: union train
                iso_B <- isoreg(u_union[train_mask_a], v_union[train_mask_a])
                pred_B <- as.stepfun(iso_B)(u_union[test_mask_a])
                # Build audit sheet
                v_test <- v_union[test_mask_a]
                err_A <- v_test - pred_A
                err_B <- v_test - pred_B
                sh5 <- data.frame(
                    barcode    = bc_union[test_mask_a],
                    u          = u_union[test_mask_a],
                    v_actual   = v_test,
                    g_hat_A    = pred_A,
                    g_hat_B    = pred_B,
                    error_A    = err_A,
                    error_B    = err_B,
                    error_sq_A = err_A^2,
                    error_sq_B = err_B^2
                )
                rmse_A_audit <- sqrt(mean(err_A^2))
                rmse_B_audit <- sqrt(mean(err_B^2))
                # Add summary row
                sh5_summary <- data.frame(
                    barcode = "=== SUMMARY ===",
                    u = NA, v_actual = NA, g_hat_A = NA, g_hat_B = NA,
                    error_A = NA, error_B = NA,
                    error_sq_A = mean(err_A^2),
                    error_sq_B = mean(err_B^2)
                )
                sh5_summary2 <- data.frame(
                    barcode = sprintf("RMSE_A=%.6f | RMSE_B=%.6f | dRMSE=%.6f",
                                      rmse_A_audit, rmse_B_audit,
                                      rmse_B_audit - rmse_A_audit),
                    u = NA, v_actual = NA, g_hat_A = NA, g_hat_B = NA,
                    error_A = NA, error_B = NA, error_sq_A = NA, error_sq_B = NA
                )
                # Both curves here are scored on one shared pool fold (schematic).
                sh5_summary3 <- data.frame(
                    barcode = paste0("SCHEMATIC: both curves scored on one shared ",
                                     "pool fold; the gate's RMSE_A is the OriRef-only ",
                                     "CV (see Stage2_Folds)."),
                    u = NA, v_actual = NA, g_hat_A = NA, g_hat_B = NA,
                    error_A = NA, error_B = NA, error_sq_A = NA, error_sq_B = NA
                )
                sh5 <- rbind(sh5, sh5_summary, sh5_summary2, sh5_summary3)
                audit_sheets[["One_Fold_Audit"]] <- sh5

                # == Sheet 6: Fold_Assignments (all 40 repeats) ==============
                set.seed(42)
                fold_matrix <- matrix(NA_integer_, nrow = length(u_union),
                                      ncol = n_reps_a)
                colnames(fold_matrix) <- paste0("repeat_", sprintf("%02d", 1:n_reps_a))
                for (rr in seq_len(n_reps_a)) {
                    fa <- integer(length(u_union))
                    for (d in 1:10) {
                        idx_d <- which(dec_audit == d)
                        .fa <- rep(seq_len(n_k_a), length.out = length(idx_d))
                        fa[idx_d] <- .fa[sample.int(length(.fa))]
                    }
                    fold_matrix[, rr] <- fa
                }
                sh6 <- data.frame(
                    barcode     = bc_union,
                    u           = u_union,
                    v           = v_union,
                    is_origref  = is_origref_union,
                    is_recovered = !is_origref_union,
                    fold_matrix,
                    stringsAsFactors = FALSE
                )
                audit_sheets[["Fold_Assignments"]] <- sh6

                # Write Excel; a failure is reported and the run continues.
                audit_path <- file.path(OUT_DIR,
                    paste0(FILE_PREFIX, "_", dataset_tag, "_RSk_fCV_Audit.xlsx"))
                tryCatch({
                    write_xlsx(audit_sheets, path = audit_path)
                    cat(sprintf("  RSk-fCV audit Excel saved: %s (%d sheets)\n",
                                basename(audit_path), length(audit_sheets)))
                }, error = function(e)
                    cat("  RSk-fCV audit Excel FAILED: ", conditionMessage(e),
                        "\n  (diagnostic output only; the run continues)\n",
                        sep = ""))


            } else {
                cat("  Plot 4 and the RSk-fCV audit Excel skipped",
                    " (no Stage 2 result on this dataset)\n", sep = "")
            }


            # Store delta calibration table for Excel
            excel_tabs[["Step4_CMRF_Delta_Params"]] <- data.frame(
                dataset              = dataset_tag,
                cal_method           = "weighted quantile(|v - g(u)|, 0.95) over ALL original SKCS-core cells",
                coverage             = cmrf_delta_cal$coverage,
                n_grp1a_total        = cmrf_delta_cal$n_grp1a,
                grp1a_resid_median   = round(cmrf_delta_cal$resid_median, 8),
                grp1a_resid_max      = round(cmrf_delta_cal$resid_max,    8),
                delta_base           = round(delta, 8),
                delta_max_at_waist   = round(delta * 1.50, 8),
                kappa                = CMRF_KAPPA,
                bw                   = CMRF_BW,
                stage2_n_repeats     = CMRF_N_REPEATS,
                stage2_k_folds       = CMRF_K_FOLDS,
                stage2_total_folds   = CMRF_N_REPEATS * CMRF_K_FOLDS,
                stage2_alpha         = CMRF_ALPHA,
                spine_method         = "mgcv g(u) fitted in Part 10e on TREP SKCS-core cells; donor RE + inverse-donor weights; PRSS-optimised; monotonicity checked, not imposed",
                adapt_formula        = "delta_adapt(u) = delta_base * (1 + kappa * exp(-(u-IP_u)^2/(2*bw^2)))",
                cv_method            = "RSk-fCV: K=5 folds x 40 repeats=200 total folds; stratification by u-quantile decile (10 strata)",
                justification        = "Stage 1: 95th weighted percentile of the original SKCS-core spine residuals. Stage 2: RSk-fCV (K=5 x 40 repeats); Curve A = OriRef alone (once), Curve B = OriRef + candidates (each iteration), matched by (repeat, fold); paired Wilcoxon; remove the highest PAVA residual until p>=0.05 OR dRMSE<=BV.",
                cv_method_note = "removal requires p<0.05 AND dRMSE > BV = IQR(RMSE_A)/median(RMSE_A)",
                stringsAsFactors     = FALSE
            )

            # == p_thresh / q_thresh for plot reference lines ================
            p_IP_full  <- mean(u_all <= IP_u)

            # =====================================================================
            # PLOT 5 -- All Orig-only & Recovered cells with CMRF tau band
            # In (p_full, q_full) space with PAVA isotonic step curve + tau ribbon.
            # =====================================================================
            cat("\n  STEP 4 -- Plot 5: Orig-only & Recovered with delta band\n")

            ORIG_ONLY_COL_04C <- "#8B6914"
            RECOV_COL         <- "#7C3AED"
            ZOOM_COLS_04C     <- c("Non-selected" = "gray60", FOUR_COLS)

            df_04c_all <- df_refined
            df_04c_all$plot_colour <- dplyr::case_when(
                df_04c_all$is_recovered ~ "Recovered",
                df_04c_all$cell_class == "Original + Refined" &
                    df_04c_all$group_orig == "Pre-IP TREP"      ~ "Pre-IP TREP",
                df_04c_all$cell_class == "Original + Refined" &
                    df_04c_all$group_orig == "Pre-IP non-TREP"  ~ "Pre-IP non-TREP",
                df_04c_all$cell_class == "Original + Refined" &
                    df_04c_all$group_orig == "Post-IP TREP"     ~ "Post-IP TREP",
                df_04c_all$cell_class == "Original + Refined" &
                    df_04c_all$group_orig == "Post-IP non-TREP" ~ "Post-IP non-TREP",
                TRUE                                             ~ "Excluded"
            )

            df_04c_excl     <- df_04c_all[df_04c_all$plot_colour == "Excluded", ]
            df_04c_excl$plot_colour <- "Non-selected"
            df_04c_skcs     <- df_04c_all[df_04c_all$cell_class == "Original + Refined", ]
            df_04c_origonly <- df_refined[df_refined$cell_class == "Original only", ]
            df_04c_recov    <- df_refined[df_refined$is_recovered, ]

            # Plot 5 is in (p_full, q_full) space; see Plot 6 for the band in (u,v).

            p05 <- ggplot() +
                geom_point(data = df_04c_excl,
                           aes(x = p_full, y = q_full, colour = plot_colour),
                           size = 0.55, alpha = 0.40) +
                geom_point(data = df_04c_skcs,
                           aes(x = p_full, y = q_full, colour = plot_colour),
                           size = 1.55, alpha = 0.65) +
                geom_point(data = df_04c_origonly,
                           aes(x = p_full, y = q_full),
                           colour = "white", fill = ORIG_ONLY_COL_04C,
                           shape = 25, size = 2.4, alpha = 0.90, stroke = 0.6) +
                geom_point(data = df_04c_recov,
                           aes(x = p_full, y = q_full),
                           colour = "white", fill = RECOV_COL,
                           shape = 23, size = 2.6, alpha = 0.93, stroke = 0.6) +
                scale_colour_manual(values = ZOOM_COLS_04C, name = "Cell group",
                                    drop = FALSE) +
                geom_vline(xintercept = p_thresh, linetype = "solid",
                           colour = "gray25", linewidth = 0.55) +
                geom_hline(yintercept = q_thresh, linetype = "solid",
                           colour = "gray25", linewidth = 0.55) +
                geom_hline(yintercept = 0.5, linetype = "dotted",
                           colour = "gray45", linewidth = 0.45) +
                annotate("text", x = p_thresh + 0.01, y = 0.03,
                         label = sprintf("Si_norm>0.05\n(p=%.2f)", p_thresh),
                         hjust = 0, vjust = 0, size = 2.4, colour = "gray25") +
                annotate("text", x = 0.02, y = q_thresh + 0.01,
                         label = sprintf("Ki>0.5\n(q=%.2f)", q_thresh),
                         hjust = 0, vjust = 0, size = 2.4, colour = "gray25") +
                annotate("text", x = 0.005, y = 0.97,
                         label = sprintf(
                             "Orig-only (dark-tan v): %d\nRecovered CMRF (violet +): %d\ndelta_adapt_base=%.5f (in u,v space)",
                             nrow(df_04c_origonly), nrow(df_04c_recov), delta),
                         hjust = 0, vjust = 1, size = 2.7, colour = "gray20") +
                coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
                labs(
                    title    = sprintf(
                        "Part10e Plot 5 -- Orig-only & Recovered (CMRF) with delta band: %s",
                        dataset_tag),
                    subtitle = paste0(
                        "CMRF classification in (p_full,q_full) space. ",
                        "Dark-gray = original Si_norm>0.05/Ki>0.5 thresholds. ",
                        "Dark-tan (v) = Orig-only (failed CMRF). ",
                        "Violet (+) = Recovered (CMRF). ",
                        "See Plot 6 for CMRF delta band in (u,v) space.\n",
                        ip_coords_line
                    ),
                    x = "p_full = F_S(Si)  [full-pop Si marginal CDF]",
                    y = "q_full = F_K(Ki)  [full-pop Ki marginal CDF]"
                ) +
                theme_copula() +
                theme(legend.position = "right", aspect.ratio = 1,
                      panel.grid.minor = element_blank())

            save_pp(p05,
                    file.path(OUT_DIR, paste0(FILE_PREFIX, "_", dataset_tag,
                                              "_05_waist_zoom_with_orig_only")),
                    pdf_w = 10, pdf_h = 8, png_w = 2200, png_h = 1800)
            cat("  Plot 5 saved\n")
            all_plots[["05_waist_zoom_with_orig_only"]] <- p05

            # Constants and frames used by Plots 6 and 7.
            RECOV_COL        <- "#7C3AED"
            ORIG_ONLY_COL_4A <- "#8B6914"
            SCALE_COLS_04a   <- c("Non-selected" = "gray60", FOUR_COLS)

            df_04a          <- df_refined
            df_04a_nonsel   <- df_04a[!df_04a$is_orig_skcs & !df_04a$is_recovered, ]
            df_04a_origonly <- df_04a[df_04a$cell_class == "Original only", ]
            df_04a_orig     <- df_04a[df_04a$cell_class == "Original + Refined", ]
            df_04a_recov    <- df_04a[df_04a$is_recovered, ]

            # =====================================================================
            # PLOT 6 -- CMRF in (u, v) GOI-cDFG copula space
            # Shows the tau tolerance band around the PAVA isotonic step curve
            # mapped to (u, v) space, with recovered cells highlighted.
            # =====================================================================
            cat("\n  STEP 4 -- Plot 6: CMRF in (u,v) space\n")

            ZOOM_COLS <- c("Non-selected" = "gray60", FOUR_COLS)

            df_04b_all <- df_refined
            df_04b_all$plot_colour <- dplyr::case_when(
                df_04b_all$is_recovered ~ "Recovered",
                df_04b_all$cell_class == "Original + Refined" &
                    df_04b_all$group_orig == "Pre-IP TREP"      ~ "Pre-IP TREP",
                df_04b_all$cell_class == "Original + Refined" &
                    df_04b_all$group_orig == "Pre-IP non-TREP"  ~ "Pre-IP non-TREP",
                df_04b_all$cell_class == "Original + Refined" &
                    df_04b_all$group_orig == "Post-IP TREP"     ~ "Post-IP TREP",
                df_04b_all$cell_class == "Original + Refined" &
                    df_04b_all$group_orig == "Post-IP non-TREP" ~ "Post-IP non-TREP",
                TRUE                                             ~ "Excluded"
            )

            df_04b_nonsel <- df_04b_all[!df_04b_all$is_orig_skcs &
                                            !df_04b_all$is_recovered, ]
            df_04b_nonsel$plot_colour <- "Non-selected"
            df_04b_skcs   <- df_04b_all[df_04b_all$cell_class == "Original + Refined", ]
            df_04b_recov  <- df_04b_all[df_04b_all$is_recovered, ]
            df_04b_origonly <- df_04b_all[df_04b_all$cell_class == "Original only", ]

            # CMRF: adaptive delta_adapt(u) band around spine
            # Dense grid for smooth ribbon display
            u_band_grid_04b <- seq(0, 1, by = 0.002)
            g_band_04b      <- gam_predict(u_band_grid_04b)
            Kw_04b          <- exp(-(u_band_grid_04b - IP_u)^2 / (2 * 0.10^2))
            tau_04b         <- delta * (1 + 0.50 * Kw_04b)
            df_04b_band <- data.frame(
                u     = u_band_grid_04b,
                g     = g_band_04b,
                upper = pmin(g_band_04b + tau_04b, 1),
                lower = pmax(g_band_04b - tau_04b, 0)
            )
            # PAVA step for diagnostic reference
            df_04b_iso_sorted <- df_refined[order(df_refined$u_all), ]

            # Factors with the full level set keep five legend keys.
            df_04b_nonsel$plot_colour <- factor(df_04b_nonsel$plot_colour,
                                                levels = names(ZOOM_COLS))
            df_04b_skcs$plot_colour   <- factor(df_04b_skcs$plot_colour,
                                                levels = names(ZOOM_COLS))

            p06 <- ggplot() +
                # Non-selected background
                geom_point(data = df_04b_nonsel,
                           aes(x = u_all, y = v_all, colour = plot_colour),
                           size = 0.55, alpha = 0.45) +
                # Adaptive delta_adapt(u) band around spine
                geom_ribbon(data = df_04b_band,
                            aes(x = u, ymin = lower, ymax = upper),
                            fill = "#607D8B", alpha = 0.18, inherit.aes = FALSE) +
                # spine g(u), solid navy
                geom_line(data = df_04b_band, aes(x = u, y = g),
                          colour = SPINE_COL_10E, linewidth = 0.9,
                          inherit.aes = FALSE) +
                # PAVA step -- diagnostic reference only (dotted)
                geom_step(data = df_04b_iso_sorted,
                          aes(x = u_all, y = yhat_iso),
                          colour = "#004D40", linewidth = 0.45,
                          linetype = "dotted", direction = "hv",
                          inherit.aes = FALSE) +
                # Original-only excluded (dark-tan triangles)
                geom_point(data = df_04b_origonly,
                           aes(x = u_all, y = v_all),
                           colour = "white", fill = ORIG_ONLY_COL_4A,
                           shape = 25, size = 2.0, alpha = 0.85, stroke = 0.5) +
                # Orig+Ref retained (four-group colours)
                geom_point(data = df_04b_skcs,
                           aes(x = u_all, y = v_all, colour = plot_colour),
                           size = 1.2, alpha = 0.55) +
                # Recovered (violet diamonds, top layer)
                geom_point(data = df_04b_recov,
                           aes(x = u_all, y = v_all),
                           colour = "white", fill = RECOV_COL,
                           shape = 23, size = 2.0, alpha = 0.91, stroke = 0.6) +
                scale_colour_manual(values = ZOOM_COLS, name = "Cell group",
                                    drop = FALSE,
                                    guide = guide_legend(
                                        override.aes = list(
                                            size  = c(1.8, rep(2.2, 4)),
                                            shape = rep(16L, 5),
                                            alpha = c(0.65, rep(0.75, 4))
                                        )
                                    )) +
                # Reference lines: IP_u vertical, IP_v horizontal, v=0.5 horizontal
                geom_vline(xintercept = IP_u, linetype = "solid",
                           colour = IP_U_COL_10E, linewidth = 0.7) +
                ip_boundary_v(linewidth = 0.7) +
                geom_hline(yintercept = 0.5, linetype = "dashed",
                           colour = "#1A237E", linewidth = 0.55) +
                annotate("text", x = IP_u + 0.01, y = 0.03,
                         label = sprintf("IP_u=%.3f", IP_u),
                         colour = IP_U_COL_10E, hjust = 0, vjust = 0, size = 2.8,
                         fontface = "bold") +
                annotate("text", x = 0.02, y = 0.52,
                         label = "v=0.5", colour = "#1A237E",
                         hjust = 0, size = 2.5, fontface = "italic") +
                # IP_v label (none when IP_v is not finite)
                (if (is.finite(IP_v))
                     annotate("text", x = 0.02, y = IP_v + 0.02,
                              label = sprintf("IP_v=%.3f (Pre/Post boundary)", IP_v),
                              colour = IP_V_COL_10E, hjust = 0, vjust = 0,
                              size = 2.8, fontface = "bold")) +
                labs(
                    title    = sprintf("Part10e Plot 6 -- CMRF filter in (u,v) space: %s",
                                       dataset_tag),
                    subtitle = paste0(
                        sprintf("IP_u=%.3f. ", IP_u),
                        "Solid navy = g(u) spine. ",
                        "Steel band = adaptive delta_adapt(u)=delta_base*(1+0.5*K(u,IP_u,bw)); delta_base=",
                        sprintf("%.5f (Orig+Ref 95th pctile)", delta), "]. ",
                        "Navy dashed = v=0.5 reference. ",
                        "Dotted dark teal = PAVA step (diag). ",
                        "Violet (+) = Recovered (CMRF). Dark-tan (v) = Orig-only excl.\n",
                        ip_coords_line, ""
                    ),
                    x = "u = Fx(GOI)  [GOI marginal CDF]",
                    y = "v = Fc(cDFG)  [cDFG marginal CDF]"
                ) +
                coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
                theme_copula() +
                theme(legend.position = "right", aspect.ratio = 1,
                      panel.grid.minor = element_blank())

            save_pp(p06,
                    file.path(OUT_DIR, paste0(FILE_PREFIX, "_", dataset_tag,
                                              "_06_CMRF_band_uv_space")),
                    pdf_w = 10, pdf_h = 8, png_w = 2200, png_h = 1800)
            cat("  Plot 6 saved\n")
            all_plots[["06_CMRF_band_uv_space"]] <- p06

            # =====================================================================
            # STEP 5 -- QUANTITATIVE COMPARISON: ORIGINAL vs CMRF
            # =====================================================================
            cat("\n  STEP 5 -- Quantitative comparison original vs CMRF\n")

            # == Plot 7: Bar chart of cell counts per group ==================
            orig_counts <- as.data.frame(table(
                group = df_04a_orig$group_orig[df_04a_orig$group_orig %in%
                                                   names(FOUR_COLS)]),
                stringsAsFactors = FALSE)
            names(orig_counts)[names(orig_counts) == "Freq"] <- "count"
            orig_counts$type <- "Original SKCS"

            # Recovered counted on the cDFG axis at IP_y.
            if (!is.finite(IP_y) || !("cDFG_expr" %in% names(df_refined)))
                stop("Recovery counts need IP_y and df_refined$cDFG_expr.",
                     call. = FALSE)
            recov_pre  <- sum(df_refined$is_recovered &
                                  df_refined$cDFG_expr <  IP_y, na.rm = TRUE)
            recov_post <- sum(df_refined$is_recovered &
                                  df_refined$cDFG_expr >= IP_y, na.rm = TRUE)

            recov_counts <- data.frame(
                group = c("Pre-IP\nRecovered", "Post-IP\nRecovered"),
                count = c(recov_pre, recov_post),
                type  = "Recovered",
                stringsAsFactors = FALSE
            )

            # Excluded = "Original only" cells
            df_excl <- df_refined[df_refined$cell_class == "Original only", ]

            # Guard: ip_side and fail_label only assignable on non-empty df_excl
            if (nrow(df_excl) > 0) {
                # cDFG axis at IP_y.
                df_excl$ip_side    <- ifelse(df_excl$cDFG_expr < IP_y,
                                             "Pre-IP", "Post-IP")
                df_excl$fail_label <- paste0(df_excl$ip_side, "\nOrig\nonly")
                excl_tbl <- as.data.frame(table(group = df_excl$fail_label),
                                          stringsAsFactors = FALSE)
                names(excl_tbl)[names(excl_tbl) == "Freq"] <- "count"
                excl_tbl$type <- "Orig only"
            } else {
                excl_tbl <- data.frame(group = character(0), count = integer(0),
                                       type  = character(0), stringsAsFactors = FALSE)
            }

            excl_levels <- c("Pre-IP\nOrig\nonly", "Post-IP\nOrig\nonly")
            excl_all <- data.frame(group = excl_levels, stringsAsFactors = FALSE)
            if (nrow(excl_tbl) > 0) {
                excl_all <- merge(excl_all, excl_tbl[, c("group", "count")],
                                  by = "group", all.x = TRUE)
            } else {
                excl_all$count <- 0L
            }
            excl_all$count[is.na(excl_all$count)] <- 0L
            excl_all$type <- "Orig only"

            df_bar <- dplyr::bind_rows(orig_counts, recov_counts, excl_all)

            df_bar$group <- factor(df_bar$group,
                                   levels = c(names(FOUR_COLS),
                                              "Pre-IP\nRecovered",
                                              "Post-IP\nRecovered",
                                              excl_levels))
            df_bar$type  <- factor(df_bar$type,
                                   levels = c("Original SKCS", "Recovered", "Orig only"))

            df_bar <- df_bar %>%
                dplyr::group_by(type) %>%
                dplyr::mutate(pct = ifelse(sum(count, na.rm = TRUE) > 0,
                                           100 * count / sum(count, na.rm = TRUE), 0)) %>%
                dplyr::ungroup()

            # Pad to equal-width facets (4 bars per facet)
            pad_recov <- data.frame(
                group = factor(excl_levels, levels = levels(df_bar$group)),
                count = NA_real_, pct = NA_real_,
                type  = factor("Recovered", levels = levels(df_bar$type)),
                stringsAsFactors = FALSE
            )
            pad_excl <- data.frame(
                group = factor(c("Pre-IP\nRecovered", "Post-IP\nRecovered"),
                               levels = levels(df_bar$group)),
                count = NA_real_, pct = NA_real_,
                type  = factor("Orig only", levels = levels(df_bar$type)),
                stringsAsFactors = FALSE
            )
            df_bar_padded <- dplyr::bind_rows(df_bar, pad_recov, pad_excl)

            # Recode group labels to two-line
            recode_group_label <- function(x) {
                dplyr::recode(as.character(x),
                              "Pre-IP TREP"         = "Pre-IP\nTREP",
                              "Pre-IP non-TREP"     = "Pre-IP\nnon-TREP",
                              "Post-IP TREP"        = "Post-IP\nTREP",
                              "Post-IP non-TREP"    = "Post-IP\nnon-TREP",
                              "Pre-IP\nRecovered"   = "Pre-IP\nRecovered",
                              "Post-IP\nRecovered"  = "Post-IP\nRecovered",
                              "Pre-IP\nOrig\nonly"  = "Pre-IP\nOrig\nonly",
                              "Post-IP\nOrig\nonly" = "Post-IP\nOrig\nonly"
                )
            }

            new_levels <- c(
                "Pre-IP\nTREP", "Pre-IP\nnon-TREP",
                "Post-IP\nTREP", "Post-IP\nnon-TREP",
                "Pre-IP\nRecovered", "Post-IP\nRecovered",
                "Pre-IP\nOrig\nonly", "Post-IP\nOrig\nonly"
            )

            BAR_COLS_05 <- c(
                setNames(unname(FOUR_COLS["Pre-IP TREP"]),    "Pre-IP\nTREP"),
                setNames(unname(FOUR_COLS["Pre-IP non-TREP"]), "Pre-IP\nnon-TREP"),
                setNames(unname(FOUR_COLS["Post-IP TREP"]),   "Post-IP\nTREP"),
                setNames(unname(FOUR_COLS["Post-IP non-TREP"]),"Post-IP\nnon-TREP"),
                "Pre-IP\nRecovered"      = "#9D5CF0",
                "Post-IP\nRecovered"     = RECOV_COL,
                "Pre-IP\nOrig\nonly"  = "#A0522D",
                "Post-IP\nOrig\nonly" = "#8B6914"
            )

            df_bar_padded$group_lbl <- factor(
                recode_group_label(as.character(df_bar_padded$group)),
                levels = new_levels)
            df_bar$group_lbl <- factor(
                recode_group_label(as.character(df_bar$group)),
                levels = new_levels)

            y_max_data <- max(df_bar$count, na.rm = TRUE)
            y_upper    <- y_max_data * 1.18

            p07 <- ggplot(df_bar_padded,
                          aes(x = group_lbl, y = count, fill = group_lbl)) +
                geom_col(data = df_bar_padded[!is.na(df_bar_padded$count), ],
                         width = 0.65, colour = "black", linewidth = 0.3) +
                geom_text(
                    data = df_bar[!is.na(df_bar$count) & df_bar$count > 0, ],
                    aes(label = sprintf("%d", count)),
                    vjust = -0.2, size = 4.2, fontface = "bold"
                ) +
                scale_fill_manual(values = BAR_COLS_05, name = NULL, drop = FALSE) +
                scale_y_continuous(
                    expand = expansion(mult = c(0.02, 0)),
                    limits = c(0, y_upper)
                ) +
                facet_wrap(~ type, scales = "free_x", nrow = 1) +
                labs(
                    title    = sprintf(
                        "Part10e Plot 7 -- Cell count: Original SKCS vs Recovered (CMRF) vs Orig only: %s",
                        dataset_tag),
                    subtitle = sprintf(
                        "Original+Refined n=%d | Recovered (CMRF) n=%d | Orig only n=%d | CMRF-core n=%d | delta_base=%.5f",
                        nrow(df_04a_orig), sum(is_recovered),
                        nrow(df_excl), sum(is_refined_skcs), delta),
                    x = NULL,
                    y = "Cell count"
                ) +
                theme_copula(base_size = 12) +
                theme(legend.position  = "none",
                      strip.text       = element_text(face = "bold", size = 13),
                      axis.text.x      = element_text(size = 11, face = "bold",
                                                      lineheight = 0.9,
                                                      hjust = 0.5, vjust = 1,
                                                      margin = margin(t = 4)),
                      axis.text.y      = element_text(size = 11),
                      axis.title.y     = element_text(size = 12, face = "bold"),
                      plot.title       = element_text(size = 14, face = "bold",
                                                      margin = margin(b = 4)),
                      plot.subtitle    = element_text(size = 9.5, margin = margin(b = 8)),
                      panel.spacing    = unit(1.2, "lines"),
                      plot.margin      = margin(t = 12, r = 15, b = 10, l = 15))

            save_pp(p07,
                    file.path(OUT_DIR, paste0(FILE_PREFIX, "_", dataset_tag,
                                              "_07_recovery_bar_comparison")),
                    pdf_w = 14, pdf_h = 7, png_w = 3000, png_h = 1600)
            cat("  Plot 7 saved\n")
            all_plots[["07_recovery_bar_comparison"]] <- p07

            # == Excel tabs for Steps 4 and 5 ================================
            skcs_meta <- data.frame(
                barcode   = skcs_df$barcode,
                IP_status = skcs_df$IP_status,
                is_TREP   = as.logical(skcs_df$is_TREP),
                EP_value  = skcs_df$EP_value,
                Si_value  = skcs_df$S_value,
                Ki_value  = skcs_df$K_value,
                stringsAsFactors = FALSE
            )
            # Left join by match() keeps df_refined's row order (stage2_res is positional).
            if (anyDuplicated(skcs_meta$barcode) > 0L)
                stop("skcs_meta has duplicate barcodes; a left join would ",
                     "change the row count of df_refined_full.", call. = FALSE)
            .meta_i   <- match(df_refined$barcode, skcs_meta$barcode)
            .meta_col <- setdiff(names(skcs_meta), "barcode")
            df_refined_full <- cbind(
                df_refined,
                skcs_meta[.meta_i, .meta_col, drop = FALSE])
            rownames(df_refined_full) <- NULL
            stopifnot(identical(df_refined_full$barcode, df_refined$barcode))
            df_refined_full$in_union_skcs <-
                df_refined_full$cell_class == "Original + Refined" |
                df_refined_full$is_recovered
            df_refined_full$dataset_tag_col <- dataset_tag

            # Fixed leading columns, then the CMRF-specific columns, then the rest.
            first_cols <- c(
                "barcode",         "dataset_tag_col",
                "in_union_skcs",   "is_orig_skcs",     "is_recovered",
                "cell_class",      "IP_status",         "is_TREP",
                "group_orig",      "group_refined",
                "u_all",           "v_all",             "K_density",
                "EP_value",        "Si_value",          "Ki_value",
                "g_u",
                "cond2_refined",     "is_refined",
                "p_full",          "q_full",
                # CMRF specific
                "yhat_iso",        "resid_abs",         "delta_base",
                "delta_adapt",     "pass_cmrf",
                "pass_stage1",     "pass_stage2",       "removed_by_stage2"
            )
            first_cols <- first_cols[first_cols %in% names(df_refined_full)]
            rest_cols  <- setdiff(names(df_refined_full), first_cols)
            df_refined_full <- df_refined_full[, c(first_cols, rest_cols)]

            # GOI-negative cells stub
            scores_goi_neg <- scores_all[is.na(scores_all$spear_contrib), ]
            rownames(scores_goi_neg) <- NULL
            n_goi_neg      <- nrow(scores_goi_neg)
            cat(sprintf("  GOI-negative cells: %d (will be appended to Step4_Refined_Cells)\n",
                        n_goi_neg))

            if (n_goi_neg > 0) {
                # Build a stub with every column in df_refined_full, typed correctly.
                # Start from a template of NAs matching df_refined_full's schema.
                goi_neg_n   <- n_goi_neg
                stub_base   <- as.data.frame(
                    lapply(df_refined_full[0, , drop = FALSE], function(col) {
                        # Replicate NA of the correct type for each column
                        rep(
                            switch(class(col)[1],
                                   "logical"   = NA,
                                   "integer"   = NA_integer_,
                                   "numeric"   = NA_real_,
                                   "character" = NA_character_,
                                   NA),
                            goi_neg_n
                        )
                    }),
                    stringsAsFactors = FALSE
                )
                rownames(stub_base) <- NULL

                # Fill in the columns we know
                stub_base$barcode           <- unname(scores_goi_neg$barcode)
                stub_base$dataset_tag_col   <- dataset_tag
                stub_base$cell_class        <- "GOI-negative"
                stub_base$in_union_skcs     <- FALSE
                stub_base$is_orig_skcs      <- FALSE
                stub_base$is_recovered      <- FALSE
                stub_base$pass_cmrf         <- FALSE
                stub_base$is_refined        <- FALSE
                stub_base$cond2_refined     <- FALSE
                stub_base$pass_stage1       <- FALSE
                stub_base$pass_stage2       <- FALSE
                stub_base$removed_by_stage2 <- FALSE
                stub_base$group_orig        <- "GOI-negative"
                stub_base$group_refined     <- "GOI-negative"
                stub_base$delta_base        <- delta

                df_step4_all <- dplyr::bind_rows(df_refined_full, stub_base)
            } else {
                df_step4_all <- df_refined_full
            }

            df_step4_all <- df_step4_all[order(
                is.na(df_step4_all$u_all), df_step4_all$u_all, na.last = TRUE), ]
            rownames(df_step4_all) <- NULL

            excel_tabs[["Step4_Refined_Cells"]] <- df_step4_all
            cat(sprintf("  Step4_Refined_Cells: %d total rows (%d GOI-pos + %d GOI-neg)\n",
                        nrow(df_step4_all),
                        sum(!is.na(df_step4_all$u_all)),
                        sum(df_step4_all$cell_class == "GOI-negative", na.rm = TRUE)))

            df_union <- df_refined_full[df_refined_full$in_union_skcs, ]
            cat(sprintf("  CMRF-core: %d cells (%d Orig+Refined + %d Recovered)\n",
                        nrow(df_union),
                        sum(df_union$cell_class == "Original + Refined"),
                        sum(df_union$is_recovered)))
            excel_tabs[["Step4_Union_SKCS"]] <- df_union

            excel_tabs[["Step4_Boundary_Params"]] <- data.frame(
                dataset           = dataset_tag,
                IP_u              = round(IP_u,              6),
                IP_x              = round(IP_x,              6),
                IP_y              = round(IP_y,              6),
                IP_v              = round(if (exists("IP_v")) IP_v else NA_real_, 6),
                class_axis        = "cDFG at IP_y",
                p_IP_full         = round(p_ip_full_val,     6),
                q_v05_full        = round(q_v05_full,        6),
                delta_base          = round(delta,               8),
                delta_max_waist     = round(delta * 1.50,        8),
                kappa             = 0.50,
                bw_gaussian       = 0.10,
                coverage          = cmrf_delta_cal$coverage,
                n_grp1a           = cmrf_delta_cal$n_grp1a,
                grp1a_resid_median = round(cmrf_delta_cal$resid_median, 8),
                filter_method     = paste0(
                    "CMRF: ",
                    "cond2=|v-g(u)|<=delta_adapt(u) [spine proximity+adaptive; band proximity for all cells]; ",
                    "g(u)=Part10e spine, donor-aware, fitted on TREP SKCS-core; ",
                    "delta_adapt(u)=delta_base*(1+kappa*K(u,IP_u,bw)); ",
                    "delta_base=weighted quantile(all original SKCS-core |v-g(u)|, 0.95) [95th pctile = p<0.05 equiv; data-driven]; ",
                    "Stage2: RSk-fCV (K=5 folds x 40 repeats=200 total folds; stratified by u-quantile decile) paired Wilcoxon"
                ),
                p_thresh_orig     = round(p_thresh,          6),
                q_thresh_orig     = round(q_thresh,          6),
                n_orig_skcs       = sum(is_orig_skcs),
                n_refined_skcs    = sum(is_refined_skcs),
                n_recovered       = sum(is_recovered),
                n_orig_only       = sum(is_orig_skcs & !is_refined_skcs),
                gam_source        = "fitted in Part 10e on TREP SKCS-core cells",
                spine_donor_re    = if (exists("gam_result")) gam_result$donor_re else NA,
                spine_n_donors    = if (exists("gam_result")) gam_result$n_donors else NA_integer_,
                spine_n_cells     = if (exists("spine_n_cells")) spine_n_cells else NA_integer_,
                spine_monotone    = if (exists("spine_mono")) spine_mono$monotone else NA,
                spine_min_diff    = if (exists("spine_mono")) spine_mono$min_diff else NA_real_,
                gam_pred_var      = gam_pred_var,
                stringsAsFactors = FALSE
            )

            excel_tabs[["Step5_Recovery_Summary"]] <- data.frame(
                dataset      = dataset_tag,
                category     = c("Original SKCS total",
                                 "Refined SKCS total (CMRF B+C pass)",
                                 "Recovered (CMRF B+C, new)",
                                 "Dropped (orig only, CMRF B+C fail)",
                                 "Recovered Pre-IP",
                                 "Recovered Post-IP"),
                n_cells      = c(sum(is_orig_skcs),
                                 sum(is_refined_skcs),
                                 sum(is_recovered),
                                 sum(is_orig_skcs & !is_refined_skcs),
                                 recov_pre,
                                 recov_post),
                pct_of_orig  = round(100 * c(
                    1,
                    sum(is_refined_skcs)  / sum(is_orig_skcs),
                    sum(is_recovered)     / sum(is_orig_skcs),
                    sum(is_orig_skcs & !is_refined_skcs) / sum(is_orig_skcs),
                    recov_pre             / sum(is_orig_skcs),
                    recov_post            / sum(is_orig_skcs)
                ), 2),
                stringsAsFactors = FALSE
            )

            # CMRF Orig+Ref residual table (audit trail for delta_base)
            excel_tabs[["Step4_CMRF_Delta_Cal"]] <- data.frame(
                dataset              = dataset_tag,
                cell_idx             = seq_along(cmrf_delta_cal$resid_all),
                grp1a_spine_resid     = round(cmrf_delta_cal$resid_all, 8),
                stringsAsFactors     = FALSE
            )

            # Stage 2 RSk-fCV fidelity gate summary (batch Wilcoxon iterative removal)
            if (!is.null(cmrf_res$stage2_res)) {
                s2 <- cmrf_res$stage2_res
                excel_tabs[["Step4_CMRF_Stage2_Gate"]] <- data.frame(
                    dataset              = dataset_tag,
                    n_grp1a              = sum(is_orig_skcs & cmrf_res$pass_stage1),
                    n_stage1_candidates  = s2$n_candidates,
                    n_repeats            = s2$n_repeats,
                    k_folds              = s2$k_folds,
                    n_total_folds        = s2$n_total_folds,
                    alpha                = s2$alpha,
                    grp1a_rmse_mean      = round(mean(s2$rmse_grp1a, na.rm = TRUE), 8),
                    grp1a_rmse_sd        = round(sd(s2$rmse_grp1a,   na.rm = TRUE), 8),
                    final_union_rmse_mean = round(mean(s2$rmse_final_union, na.rm = TRUE), 8),
                    wilcox_p_final       = signif(s2$wilcox_p_final, 6),
                    # The two terms behind the gate decision.
                    dRMSE_final          = signif(s2$relative_delta_final, 6),
                    BV                   = signif(s2$noise_floor, 6),
                    n_removed_by_gate    = s2$n_removed,
                    n_accepted           = sum(s2$pass_stage2),
                    # Names the criterion that closed the gate.
                    gate_decision        = if (isTRUE(s2$n_removed == 0L))
                        "PASS (no candidate removed)"
                    else if (isTRUE(s2$wilcox_p_final >= s2$alpha))
                        "PARTIAL (closed on the Wilcoxon criterion)"
                    else if (isTRUE(s2$relative_delta_final <= s2$noise_floor))
                        "PARTIAL (closed on the effect-size criterion)"
                    else "PARTIAL (all candidates removed)",
                    stringsAsFactors     = FALSE
                )
                # Per-fold RMSE comparison for audit
                excel_tabs[["Step4_CMRF_Stage2_Folds"]] <- data.frame(
                    dataset          = dataset_tag,
                    rep_fold         = seq_along(s2$rmse_grp1a),
                    rmse_grp1a       = round(s2$rmse_grp1a,       8),
                    rmse_final_union = round(s2$rmse_final_union, 8),
                    cv_method        = "RSk-fCV: each row is one fold of one repeat (K=5 folds x 40 repeats=200 rows); structurally paired between Orig+Ref and Union",
                    stringsAsFactors = FALSE
                )
            }

        }  # end if scores_goi available

        # =========================================================================
        # Assemble Excel output
        # =========================================================================
        cat("\n  Assembling Excel output ...\n")

        excel_tabs[["Step1_Cell_Data"]] <- data.frame(
            dataset  = dataset_tag,
            barcode  = skcs_df$barcode,
            ip_stat  = ip_stat,
            is_trep  = is_trep,
            group    = as.character(group),
            Si_obs   = round(Si_obs, 6),
            Ki_obs   = round(Ki_obs, 6),
            EP_raw   = round(EP_raw, 6),
            EP_plot  = round(EP_plot, 6),
            p_i      = round(p_i, 6),
            q_i      = round(q_i, 6),
            u        = round(u, 6),
            v        = round(v, 6),
            stringsAsFactors = FALSE
        )

        excel_tabs[["Step3_Group_Stats"]] <- group_stats_full
        if (exists("grid_df")) excel_tabs[["Step3_Copula_Grid"]] <- grid_df

        # == EP regression models (additive, interaction, full) ==============
        # These summaries feed the Excel regression sheets below.
        frac_preip <- mean(ip_stat == "Pre-IP", na.rm = TRUE)

        lm_add  <- lm(EP_raw ~ Si_obs + Ki_obs)
        lm_int  <- lm(EP_raw ~ Si_obs * Ki_obs)
        lm_full <- lm(EP_raw ~ Si_obs * Ki_obs + ip_stat + is_trep)
        sum_add  <- summary(lm_add)
        sum_int  <- summary(lm_int)
        sum_full <- summary(lm_full)

        rename_coef_df <- function(df, model_label) {
            df <- as.data.frame(df)
            df$term  <- rownames(df)
            df$model <- model_label
            rownames(df) <- NULL
            on <- names(df); nn <- on
            nn[grepl("Estimate",   on, fixed = TRUE)] <- "estimate"
            nn[grepl("Std. Error", on, fixed = TRUE)] <- "std_error"
            nn[grepl("t value",    on, fixed = TRUE)] <- "t_value"
            nn[grepl("Pr(",        on, fixed = TRUE)] <- "p_value"
            names(df) <- nn
            df[, c("model", "term", "estimate", "std_error", "t_value", "p_value")]
        }
        coef_int  <- rename_coef_df(coef(sum_int),  "EP ~ Si * Ki")
        coef_full <- rename_coef_df(coef(sum_full), "EP ~ Si*Ki + IP + TREP")
        reg_coefs <- rbind(coef_int, coef_full)

        reg_summary <- data.frame(
            model   = c("EP ~ Si + Ki", "EP ~ Si * Ki", "EP ~ Si*Ki + IP + TREP"),
            R2adj   = round(c(sum_add$adj.r.squared,
                              sum_int$adj.r.squared,
                              sum_full$adj.r.squared), 5),
            F_stat  = round(c(sum_add$fstatistic[1],
                              sum_int$fstatistic[1],
                              sum_full$fstatistic[1]), 3),
            n_cells = n,
            stringsAsFactors = FALSE
        )
        excel_tabs[["Step3_Regression_Summary"]] <- reg_summary
        excel_tabs[["Step3_Regression_Coefs"]]   <- reg_coefs

        overall_dep <- data.frame(
            dataset    = dataset_tag,
            n_SKCS     = n,
            rho_SK     = round(rho_sk, 6),
            tau_SK     = round(tau_sk, 6),
            KS_p_pi    = round(ks_p, 6),
            KS_p_qi    = round(ks_q, 6),
            IP_u       = round(IP_u, 6),
            frac_preip = round(frac_preip, 6),
            stringsAsFactors = FALSE
        )
        excel_tabs[["Summary"]] <- overall_dep

        desired_order <- c(
            "Summary",
            "Step1_Cell_Data",
            "Step2d_FullPop_Summary",
            "Step3_Group_Stats",
            "Step3_Copula_Grid",
            "Step3_Regression_Summary",
            "Step3_Regression_Coefs",
            "Step3e_Threshold_Translation",
            "Step4_Union_SKCS",
            "Step4_Refined_Cells",
            "Step4_Boundary_Params",
            "Step4_CMRF_Delta_Params",
            "Step4_CMRF_Delta_Cal",
            "Step4_CMRF_Stage2_Gate",
            "Step4_CMRF_Stage2_Folds",
            "Step5_Recovery_Summary"
        )
        present_tabs   <- desired_order[desired_order %in% names(excel_tabs)]
        remaining_tabs <- setdiff(names(excel_tabs), present_tabs)
        excel_tabs     <- excel_tabs[c(present_tabs, remaining_tabs)]

        xlsx_path <- file.path(OUT_DIR,
                               paste0(FILE_PREFIX, "_", dataset_tag,
                                      "_SK_interaction_analysis.xlsx"))
        tryCatch({
            write_xlsx(excel_tabs, path = xlsx_path)
            cat(sprintf("  Excel saved: %s\n", basename(xlsx_path)))
        }, error = function(e)
            cat("  Excel write FAILED: ", conditionMessage(e),
                "\n  (close the workbook if it is open, then rerun)\n", sep = ""))

        # == Display all plots ===============================================
        cat("  Displaying plots in RStudio ...\n")
        # Close any file device leaked by save_pp before displaying.
        base_dev <- dev.cur()
        leaked_devs <- setdiff(dev.list(), c(1L, base_dev))
        if (length(leaked_devs) > 0L) {
            for (ld in leaked_devs) tryCatch(dev.off(ld), error=function(e) NULL)
            cat(sprintf("  [closed %d leaked device(s) before display loop]\n",
                        length(leaked_devs)))
        }
        for (nm in names(all_plots)) {
            cat("    ->", nm, "\n")
            p_obj <- all_plots[[nm]]
            if (is.null(p_obj)) {
                cat("    [skipped: NULL plot object]\n")
                next
            }
            tryCatch({
                if (inherits(p_obj, "ggplot")) {
                    print(p_obj)
                } else {
                    grid.newpage(); grid.draw(p_obj)
                }
            }, error = function(e) {
                cat(sprintf("    [display error for %s: %s]\n", nm, conditionMessage(e)))
            })
        }

        invisible(list(
            dataset_tag  = dataset_tag,
            n            = n,
            IP_x         = IP_x,
            IP_u         = IP_u,
            IP_u_10d     = IP_u_10d,
            IP_u_source  = IP_u_source,
            ip_axis_check = ip_axis_check,
            p_i = p_i, q_i = q_i,
            Si_obs = Si_obs, Ki_obs = Ki_obs,
            EP_raw = EP_raw, EP_plot = EP_plot,
            group  = group,
            rho_SK = rho_sk, tau_SK = tau_sk,
            group_stats  = group_stats_full,
            grid_df      = if (exists("grid_df")) grid_df else NULL,
            lm_add = lm_add, lm_int = lm_int, lm_full = lm_full,
            # CMRF results
            cmrf_delta_cal      = if (exists("cmrf_delta_cal"))  cmrf_delta_cal  else NULL,
            delta        = if (exists("delta"))    delta    else NULL,
            stage2_res   = if (exists("cmrf_res") && !is.null(cmrf_res$stage2_res)) cmrf_res$stage2_res else NULL,
            df_refined   = if (exists("df_refined_full")) df_refined_full else
                if (exists("df_refined"))      df_refined      else NULL,
            df_union     = if (exists("df_union")) df_union else NULL,
            # Scalar summaries
            p_thresh     = if (exists("p_thresh"))    p_thresh    else NULL,
            q_thresh     = if (exists("q_thresh"))    q_thresh    else NULL,
            n_orig_skcs  = if (exists("is_orig_skcs")) sum(is_orig_skcs) else NULL,
            n_recovered  = if (exists("is_recovered")) sum(is_recovered) else NULL,
            n_refined    = if (exists("is_refined_skcs")) sum(is_refined_skcs) else NULL,
            # Spine for downstream
            gam_obj      = if (exists("gam_obj"))      gam_obj      else NULL,
            gam_pred_var = if (exists("gam_pred_var")) gam_pred_var else NULL,
            spine_prss   = if (exists("gam_result"))
                gam_result[c("best_k", "best_iter", "best_prss", "actual_iters",
                             "early_stopped", "max_k", "phase2_max", "n_iter_max",
                             "donor_re", "n_donors")] else NULL,
            spine_monotone = if (exists("spine_mono")) spine_mono else NULL,
            spine_n_cells  = if (exists("spine_n_cells")) spine_n_cells else NULL,
        g_at_IP_u       = if (exists("g_at_IP_u"))  g_at_IP_u  else NULL,
        g_inv_IP_v      = if (exists("g_inv_IP_v")) g_inv_IP_v else NULL,
        IP_v           = if (exists("IP_v")) IP_v else NULL,
        IP_y           = if (exists("IP_y")) IP_y else NULL,
        # Old names kept as aliases for downstream scripts.
        u_IP           = if (exists("IP_u")) IP_u else NULL,   # deprecated
        v_IP           = if (exists("IP_v")) IP_v else NULL,   # deprecated
        x_IP           = if (exists("IP_x")) IP_x else NULL,   # deprecated
        y_IP           = if (exists("IP_y")) IP_y else NULL,   # deprecated
            spine_u_range  = if (exists("spine_u_range")) spine_u_range else NULL,
            n_donors_core  = if (exists("donor_cells"))
                n_donor_levels_10e(donor_cells) else NULL,
        # Counted over the full level set, so an empty class reports 0.
        cmrf_class_counts = if (exists("cell_class"))
            table(factor(cell_class,
                         levels = c("Original + Refined", "Original only",
                                    "Recovered", "Excluded"))) else NULL,
            excel_tabs   = excel_tabs,
            plot_names   = names(all_plots),
            t2_mg_result = if (exists("t2_mg_result")) t2_mg_result else NULL
        ))
    }



    ##########################################################################
    # SECTION C: EXECUTE
    ##########################################################################

    all_results_10e <- vector("list", length(DATASETS))
    names(all_results_10e) <- DATASETS

    for (ds in DATASETS) {
        # Single-bracket assignment keeps a failed dataset as NULL in the list.
        all_results_10e[ds] <- list(tryCatch(
            run_part10e(ds),
            error = function(e) {
                cat("  ERROR in dataset", ds, ":", conditionMessage(e), "\n")
                NULL
            }
        ))
    }

    # rds_save() strips model environments and drops plot objects.
    rds_path <- file.path(OUT_DIR, paste0(FILE_PREFIX, "_all_results.rds"))
    rds_save(all_results_10e, rds_path, label = "Results object")
    # =============================================================================
    # IP_u COMPARISON LOG -- SKCS-core IP vs Part10d full-population IP
    # =============================================================================
    uip_log_path <- file.path(OUT_DIR, paste0(FILE_PREFIX, "_IP_u_comparison.txt"))
    uip_log <- file(uip_log_path, open = "wt")
    # Close the connection on error.
    on.exit(try(close(uip_log), silent = TRUE), add = TRUE)

    writeLines(c(
        strrep("=", 72),
        sprintf("  Part10e IP_u Comparison: SKCS-core vs Full-population"),
        sprintf("  Generated: %s", Sys.time()),
        strrep("=", 72), ""
    ), uip_log)

    for (ds in DATASETS) {
        res <- all_results_10e[[ds]]
        if (is.null(res)) {
            writeLines(sprintf("  %s: SKIPPED (run failed)", ds), uip_log)
            next
        }

        writeLines(c(
            strrep("-", 72),
            sprintf("  Dataset: %s", ds),
            strrep("-", 72),
            sprintf("  IP_x (SKCS-core, Part10b_IPRS):  %.6f",
                    ifelse(is.null(res$IP_x) || is.na(res$IP_x), NA, res$IP_x)),
            sprintf("  IP_u used (SKCS-core -> full-pop): %.6f  [source: %s]",
                    res$IP_u, res$IP_u_source),
            sprintf("  IP_u Part10d (full-pop computed):   %.6f",
                    ifelse(is.null(res$IP_u_10d) || is.na(res$IP_u_10d), NA, res$IP_u_10d)),
            sprintf("  Delta (used - Part10d):             %.6f",
                    ifelse(is.null(res$IP_u_10d) || is.na(res$IP_u_10d), NA,
                           res$IP_u - res$IP_u_10d)),
            ""
        ), uip_log)

        # == Recompute Stage 1 under Part10d IP_u to compare classifications ====
        df_r <- res$df_refined
        if (!is.null(df_r) && !is.na(res$IP_u_10d) && !is.na(res$IP_u) &&
            abs(res$IP_u - res$IP_u_10d) > 1e-8) {

            u_vec <- df_r$u_all
            v_vec <- df_r$v_all

            # Classify under Part10d IP_u for comparison (IP-arm assignment only)
            cell_class_cur <- df_r$cell_class

            # Category breakdown for current run
            is_rec     <- !is.na(df_r$is_recovered) & df_r$is_recovered
            is_orionly <- !is.na(cell_class_cur) & cell_class_cur == "Original only"

            # IP status for each cell based on current IP_u
            ip_cur <- ifelse(u_vec >= res$IP_u, "Post-IP", "Pre-IP")

            n_rec_pre  <- sum(is_rec & ip_cur == "Pre-IP", na.rm = TRUE)
            n_rec_post <- sum(is_rec & ip_cur == "Post-IP", na.rm = TRUE)
            n_oo_pre   <- sum(is_orionly & ip_cur == "Pre-IP", na.rm = TRUE)
            n_oo_post  <- sum(is_orionly & ip_cur == "Post-IP", na.rm = TRUE)

            # IP status under Part10d IP_u
            ip_10d <- ifelse(u_vec >= res$IP_u_10d, "Post-IP", "Pre-IP")

            n_rec_pre_10d  <- sum(is_rec & ip_10d == "Pre-IP", na.rm = TRUE)
            n_rec_post_10d <- sum(is_rec & ip_10d == "Post-IP", na.rm = TRUE)
            n_oo_pre_10d   <- sum(is_orionly & ip_10d == "Pre-IP", na.rm = TRUE)
            n_oo_post_10d  <- sum(is_orionly & ip_10d == "Post-IP", na.rm = TRUE)

            # Split on the u axis (two IP_u values compared); rows carry "[u]".
            writeLines(c(
                "  Category counts (SKCS-core IP_u vs Part10d IP_u):",
                "",
                "  [u] = split on the GOI-axis copula coordinate u at IP_u,",
                "  not the workbook's Pre/Post split (cDFG at IP_y).",
                "",
                sprintf("                          SKCS IP_u   Part10d IP_u   Diff"),
                sprintf("                          (%.4f)     (%.4f)",
                        res$IP_u, res$IP_u_10d),
                sprintf("    Pre-IP [u] Recovered:  %4d          %4d         %+d",
                        n_rec_pre, n_rec_pre_10d,
                        n_rec_pre - n_rec_pre_10d),
                sprintf("    Post-IP [u] Recovered: %4d          %4d         %+d",
                        n_rec_post, n_rec_post_10d,
                        n_rec_post - n_rec_post_10d),
                sprintf("    Pre-IP [u] Orig-only:  %4d          %4d         %+d",
                        n_oo_pre, n_oo_pre_10d,
                        n_oo_pre - n_oo_pre_10d),
                sprintf("    Post-IP [u] Orig-only: %4d          %4d         %+d",
                        n_oo_post, n_oo_post_10d,
                        n_oo_post - n_oo_post_10d),
                "",
                "  Note: Recovered vs Excluded does not depend on IP_u; only the Pre/Post labelling differs.",
                ""
            ), uip_log)

        } else if (!is.null(df_r)) {
            writeLines("  IP_u values are identical -- no CMRF classification differences.\n",
                       uip_log)
        }
    }

    writeLines(c(
        strrep("=", 72),
        "  COMPARISON COMPLETE",
        strrep("=", 72)
    ), uip_log)
    close(uip_log)
    cat(sprintf("  IP_u comparison log saved: %s\n", basename(uip_log_path)))

    # =============================================================================
    # CMRF VERIFICATION EXCEL
    # =============================================================================
    cat("\n  Building CMRF verification Excel ...\n")

    ref_verify_sheets  <- list()
    all_plots_verify   <- list()   # stores per-dataset 4-panel verification plots

    set.seed(42)   # ensures sw_pval sample() subsampling is reproducible across runs

    for (ds in DATASETS) {
        res <- all_results_10e[[ds]]
        if (is.null(res)) next

        rc   <- res$excel_tabs[["Step4_Refined_Cells"]]
        bp   <- res$excel_tabs[["Step4_Boundary_Params"]]
        dp   <- res$excel_tabs[["Step4_CMRF_Delta_Params"]]
        s2   <- res$excel_tabs[["Step4_CMRF_Stage2_Gate"]]
        s2f  <- res$excel_tabs[["Step4_CMRF_Stage2_Folds"]]
        dcal <- res$excel_tabs[["Step4_CMRF_Delta_Cal"]]

        if (is.null(rc) || is.null(bp)) {
            cat(sprintf("  Skipping %s: refined cells not available\n", ds))
            next
        }

        # == 1. Boundary parameters (from provided data) =====================
        ref_verify_sheets[[paste0(ds, "_Boundary_Params")]] <- bp

        # == 2. Delta calibration parameters =================================
        if (!is.null(dp)) ref_verify_sheets[[paste0(ds, "_Delta_Params")]] <- dp

        # == 3. Stage 2 gate summary =========================================
        if (!is.null(s2)) ref_verify_sheets[[paste0(ds, "_Stage2_Gate")]] <- s2

        # == 4. Stage 2 per-fold RMSE (Orig+Ref vs final union) ==============
        if (!is.null(s2f)) ref_verify_sheets[[paste0(ds, "_Stage2_Folds")]] <- s2f

        # == 5. Delta calibration residuals ==================================
        if (!is.null(dcal)) ref_verify_sheets[[paste0(ds, "_Delta_Cal_Resid")]] <- dcal

        # =========================================================================
        # COMPREHENSIVE STATISTICAL TESTS
        #   SW-gated: t-test or Mann-Whitney U (two groups); ANOVA + Holm-Sidak or
        #   Kruskal-Wallis + Dunn (more groups). BH FDR across tests. JT, Hoeffding D
        #   and dCor are distribution-free.
        # =========================================================================
        rc_goi   <- rc[!is.na(rc$u_all), ]
        IP_u_val <- bp$IP_u                       # band-kernel centre only
        # Pre/Post selection below is on the cDFG axis at IP_y.
        IP_y_val <- if ("IP_y" %in% names(bp))
            suppressWarnings(as.numeric(bp$IP_y)[1]) else NA_real_

        orig_ref    <- rc_goi[rc_goi$cell_class == "Original + Refined", ]
        recovered   <- rc_goi[rc_goi$cell_class == "Recovered", ]
        orig_only   <- rc_goi[rc_goi$cell_class == "Original only", ]
        excl_s1fail <- rc_goi[!rc_goi$is_orig_skcs & !rc_goi$pass_stage1, ]

        stat_rows <- list()

        ph_str <- function(ph) {
            if (is.null(ph) || nrow(ph)==0) return("none")
            if ("comparison" %in% names(ph))
                paste(sprintf("%s p_BH=%.4f",  ph$comparison,   ph$p_BH),        collapse="; ") else
                    paste(sprintf("%s vs %s p_holm=%.4f", ph$group1, ph$group2, ph$p_holm_sidak),
                          collapse="; ")
        }

        # =========================================================================
        # T2: resid_abs -- multi-group across 4 classes (SW-gated: ANOVA or KW)
        # =========================================================================
        # Recomputed from rc_goi, the cell set drawn in VP1.
        t1_mg <- multi_grp(rc_goi$resid_abs, rc_goi$cell_class,
                           levs=c("Original + Refined","Recovered","Original only","Excluded"))
        stat_rows[["T2_resid_multigroup"]] <- data.frame(
            dataset        = ds,
            test_type      = "Multi-group (SW-gated)",
            test           = t1_mg$test_used,
            comparison     = "resid_abs across Original+Refined / Recovered / Original-only / Excluded",
            H0             = "No difference in resid_abs among cell classes",
            normality_SW   = t1_mg$sw_note,
            all_normal     = isTRUE(t1_mg$all_norm),
            overall_stat   = t1_mg$stat,
            p_value        = t1_mg$p,
            posthoc_method = t1_mg$ph_method,
            posthoc_summary = ph_str(t1_mg$posthoc),
            interpretation = ifelse(!is.na(t1_mg$p) && t1_mg$p < 0.05,
                                    "PASS: Significant difference in resid_abs across cell classes",
                                    "NOTE: No significant multi-group residual difference"),
            stringsAsFactors = FALSE)
        if (!is.null(t1_mg$posthoc) && nrow(t1_mg$posthoc) > 0) {
            t1_mg$posthoc$dataset <- ds
            ref_verify_sheets[[paste0(ds, "_T2_Posthoc")]] <- t1_mg$posthoc
        }
        t2_posthoc_canonical <- t1_mg$posthoc   # reused by VP1-B FDR subtitle

        # =========================================================================
        # T3: Jonckheere-Terpstra ordered trend test on resid_abs (distribution-free)
        # Tests H1: resid_abs decreases monotonically across 4 classes ordered by
        # expected monotone fidelity: Excluded > Orig-only > Recovered > Orig+Refined
        # =========================================================================
        {
            jt_vals  <- rc_goi$resid_abs[rc_goi$cell_class %in%
                                             c("Excluded","Original only","Recovered","Original + Refined")]
            jt_grps  <- rc_goi$cell_class[rc_goi$cell_class %in%
                                              c("Excluded","Original only","Recovered","Original + Refined")]
            # Ordered factor: increasing monotone fidelity (decreasing resid_abs expected)
            jt_order <- factor(jt_grps,
                               levels = c("Excluded","Original only","Recovered","Original + Refined"),
                               ordered = TRUE)
            jt_ok    <- !is.na(jt_vals) & !is.na(jt_grps)
            jt_res   <- tryCatch(
                clinfun::jonckheere.test(jt_vals[jt_ok], as.integer(jt_order[jt_ok]),
                                         alternative = "decreasing", nperm = CMRF_JT_NPERM),
                error = function(e) list(statistic = NA_real_, p.value = NA_real_)
            )
            stat_rows[["T3_JT_ordered_trend"]] <- data.frame(
                dataset        = ds,
                test_type      = "Ordered trend (distribution-free)",
                test           = sprintf("Jonckheere-Terpstra (H1: resid_abs decreasing, nperm=%d)", CMRF_JT_NPERM),
                comparison     = "resid_abs ordered: Excluded > Orig-only > Recovered > Orig+Refined",
                H0             = "No monotone decrease in resid_abs across ordered classes",
                n_excluded     = sum(jt_grps == "Excluded",        na.rm = TRUE),
                n_orig_only    = sum(jt_grps == "Original only",   na.rm = TRUE),
                n_recovered    = sum(jt_grps == "Recovered",       na.rm = TRUE),
                n_orig_ref     = sum(jt_grps == "Original + Refined", na.rm = TRUE),
                JT_statistic   = if (is.null(jt_res$statistic)) NA_real_ else
                    round(as.numeric(jt_res$statistic), 4),
                p_value        = if (is.null(jt_res$p.value)) NA_real_ else
                    signif(jt_res$p.value, 6),
                interpretation = ifelse(!is.na(jt_res$p.value) && jt_res$p.value < 0.05,
                                        "PASS: resid_abs decreases monotonically across ordered fidelity classes",
                                        "NOTE: No significant ordered decrease"),
                stringsAsFactors = FALSE)
            cat(sprintf("  %s T3 JT: stat=%.2f p=%.4f\n", ds,
                        if (is.null(jt_res$statistic)) NA_real_ else as.numeric(jt_res$statistic),
                        if (is.null(jt_res$p.value)) NA_real_ else jt_res$p.value))
        }

        # =========================================================================
        # T7: resid_abs Recovered vs Excluded (2-group, SW-gated)
        # =========================================================================
        if (nrow(recovered) > 1 && nrow(excl_s1fail) > 1) {
            tg4 <- two_grp(recovered$resid_abs, excl_s1fail$resid_abs,
                           alternative="two.sided", paired=FALSE)
            stat_rows[["T7_resid_Recov_vs_Excluded"]] <- data.frame(
                dataset      = ds, test_type = "2-group (SW-gated)", test = tg4$test_used,
                comparison   = "resid_abs: Recovered vs Excluded (non-SKCS Stage1-fail)",
                H0           = "No difference in resid_abs (two-sided); direction read from medians",
                n_group1     = nrow(recovered),    n_group2 = nrow(excl_s1fail),
                median_g1    = round(median(recovered$resid_abs, na.rm=TRUE), 6),
                median_g2    = round(median(excl_s1fail$resid_abs, na.rm=TRUE), 6),
                SW_p_g1      = tg4$sw_p1,          SW_p_g2 = tg4$sw_p2,
                normal_g1    = tg4$norm1,           normal_g2 = tg4$norm2,
                statistic    = tg4$stat,            p_value = signif(tg4$p, 4),
                interpretation = ifelse(
                    !is.na(tg4$p) && tg4$p < 0.05 &&
                        median(recovered$resid_abs, na.rm=TRUE) <
                        median(excl_s1fail$resid_abs, na.rm=TRUE),
                    "PASS: Recovered residuals significantly lower than Excluded (two-sided; direction confirmed)",
                    ifelse(!is.na(tg4$p) && tg4$p < 0.05,
                           "NOTE: significant but Recovered NOT lower than Excluded",
                           "NOTE: no significant difference (two-sided)")),
                stringsAsFactors = FALSE)
        }

        # =========================================================================
        # T7: Hoeffding D test for independence in (u,v) within Recovered
        # Distribution-free; permutation p-value (CMRF_HD_NPERM shuffles).
        # =========================================================================
        if (nrow(recovered) > 4) {
            hd_res <- tryCatch({
                u_hd <- recovered$u_all
                v_hd <- recovered$v_all
                D_obs <- Hmisc::hoeffd(u_hd, v_hd)$D[1,2]
                # Permutation: shuffle u, recompute D each time
                set.seed(42L)
                D_perm <- vapply(seq_len(CMRF_HD_NPERM), function(i) {
                    Hmisc::hoeffd(sample(u_hd), v_hd)$D[1,2]
                }, numeric(1))
                # Phipson & Smyth (2010) correction: (b+1)/(nperm+1)
                p_perm <- (sum(D_perm >= D_obs) + 1L) / (CMRF_HD_NPERM + 1L)
                list(D = D_obs, P = p_perm, nperm = CMRF_HD_NPERM)
            }, error = function(e) list(D = NA_real_, P = NA_real_,
                                         nperm = CMRF_HD_NPERM))
            stat_rows[["T7_hoeffding_D_Recovered"]] <- data.frame(
                dataset     = ds,
                test_type   = "Independence test (distribution-free -- no SW; permutation)",
                test        = sprintf("Hoeffding D (H1: u,v not independent; R=%d permutations)", CMRF_HD_NPERM),
                comparison  = "Hoeffding D for (u_all, v_all) within Recovered cells",
                H0          = "u_all and v_all are independent in Recovered cells",
                n_group1    = nrow(recovered),
                D_statistic = round(hd_res$D, 6),
                p_value     = signif(hd_res$P, 4),
                interpretation = ifelse(!is.na(hd_res$P) && hd_res$P < 0.05,
                                        "PASS: Significant dependence in Recovered (u,v) -- monotone association supported",
                                        "NOTE: No significant dependence detected"),
                stringsAsFactors = FALSE)
            cat(sprintf("  %s T7 Hoeffding D=%.4f p=%.4f (permutation, R=%d)\n", ds,
                        if(is.na(hd_res$D)) NA_real_ else hd_res$D,
                        if(is.na(hd_res$P)) NA_real_ else hd_res$P,
                        CMRF_HD_NPERM))
        }

        # =========================================================================
        # T8: Distance correlation test for (u,v) within Recovered
        # energy::dcor.test(); distribution-free; H1: dCor(u,v) > 0.
        # =========================================================================
        if (nrow(recovered) > 4) {
            dc_res <- tryCatch({
                u_dc <- recovered$u_all[!is.na(recovered$u_all) & !is.na(recovered$v_all)]
                v_dc <- recovered$v_all[!is.na(recovered$u_all) & !is.na(recovered$v_all)]
                # R=999 permutations for stable p-value; set seed for reproducibility
                set.seed(42L)
                dc_test <- energy::dcor.test(u_dc, v_dc, R = 999L)
                # unname()[[1]] extracts the values across package versions.
                dcor_val <- unname(dc_test$statistic)[[1]]
                p_val    <- unname(dc_test$p.value)[[1]]
                list(dcor = dcor_val,
                     p    = p_val,
                     n    = length(u_dc))
            }, error = function(e) list(dcor = NA_real_, p = NA_real_,
                                        n    = NA_integer_))
            stat_rows[["T8_dcor_Recovered"]] <- data.frame(
                dataset     = ds,
                test_type   = "Distance correlation (distribution-free -- no SW)",
                test        = "Distance correlation dCor (H1: dCor(u,v) > 0; R=999 permutations)",
                comparison  = "dCor(u_all, v_all) within Recovered cells",
                H0          = "u_all and v_all are independent (dCor = 0) in Recovered cells",
                n_group1    = dc_res$n,
                dCor        = round(dc_res$dcor, 6),
                p_value     = round(dc_res$p,    6),
                interpretation = ifelse(!is.na(dc_res$p) && dc_res$p < 0.05,
                                        "PASS: Significant distance correlation in Recovered (u,v) -- dependence confirmed",
                                        "NOTE: No significant distance correlation detected"),
                stringsAsFactors = FALSE)
            cat(sprintf("  %s T8 dCor=%.4f p=%.4f\n", ds,
                        if(is.na(dc_res$dcor)) NA_real_ else dc_res$dcor,
                        if(is.na(dc_res$p))    NA_real_ else dc_res$p))
        }

        # =========================================================================
        # T1 (CONFIRMATORY GATE): RSk-fCV RMSE paired Wilcoxon (matches gate's own test)
        # =========================================================================
        if (!is.null(s2f) && all(c("rmse_grp1a","rmse_final_union") %in% names(s2f))) {
            cc7 <- complete.cases(s2f$rmse_grp1a, s2f$rmse_final_union)
            xg  <- s2f$rmse_grp1a[cc7]
            xu  <- s2f$rmse_final_union[cc7]
            wsr_t1 <- suppressWarnings(
                wilcox.test(xu, xg, paired = TRUE, alternative = "greater", exact = FALSE))
            stat_rows[["T1_confirmatory_gate"]] <- data.frame(
                dataset    = ds, test_type = "2-group paired (Wilcoxon)",
                test       = "Wilcoxon signed-rank (greater)",
                comparison = "RSk-fCV RMSE: final union vs Orig+Ref alone (K=5 folds x 40 repeats=200 paired fold RMSEs)",
                cv_method  = "RSk-fCV: stratified by u-quantile decile; every cell tested exactly once per repeat",
                H0         = "Final union RMSE <= Orig+Ref RMSE (dual gate: p>=0.05 OR delta<=noise_floor)",
                alpha_used = 0.05,
                n_total_folds = length(xg),
                mean_grp1a = round(mean(xg), 8), mean_union = round(mean(xu), 8),
                statistic  = round(as.numeric(wsr_t1$statistic), 4),
                p_value    = round(wsr_t1$p.value, 6),
                interpretation = ifelse(!is.na(wsr_t1$p.value) && wsr_t1$p.value >= 0.05,
                                        "PASS: Final union RMSE not significantly worse than Orig+Ref (p>=0.05; further relaxed gate)",
                                        "FAIL: Final union RMSE significantly worse than Orig+Ref (p<0.001) -- review Stage 2"),
                stringsAsFactors = FALSE)
        }

        # =========================================================================
        # T12: TREP enrichment Fisher's exact (categorical -- no SW)
        # =========================================================================
        # TREP lookup over all cells with a known is_TREP value.
        trep_src    <- rc[!is.na(rc$barcode) & !is.na(rc$is_TREP),
                          c("barcode","is_TREP")]
        trep_lookup <- setNames(as.logical(trep_src$is_TREP), trep_src$barcode)
        # Post-IP selected on cDFG >= IP_y.
        if (!is.finite(IP_y_val) || !("cDFG_expr" %in% names(rc_goi)))
            stop("The Post-IP contingency test needs IP_y and cDFG_expr.", call. = FALSE)
        recov_postip <- recovered[!is.na(recovered$cDFG_expr) &
                                      recovered$cDFG_expr >= IP_y_val, ]
        excl_postip  <- excl_s1fail[!is.na(excl_s1fail$cDFG_expr) &
                                        excl_s1fail$cDFG_expr >= IP_y_val, ]
        resolve_trep <- function(df) {
            tv  <- df$is_TREP; nm <- is.na(tv)
            if (any(nm) && "barcode" %in% names(df)) tv[nm] <- trep_lookup[df$barcode[nm]]
            tv
        }
        rv <- resolve_trep(recov_postip); ev <- resolve_trep(excl_postip)
        nrt <- sum(rv==TRUE, na.rm=TRUE); nrn <- sum(rv!=TRUE, na.rm=TRUE)
        net <- sum(ev==TRUE, na.rm=TRUE); nen <- sum(ev!=TRUE, na.rm=TRUE)
        if (nrow(recov_postip)>0 && nrow(excl_postip)>0 &&
            (nrt+nrn)>0 && (net+nen)>0) {
            ct11 <- matrix(c(nrt,nrn,net,nen), nrow=2, byrow=TRUE,
                           dimnames=list(c("Recovered_PostIP","Excl_PostIP"),
                                         c("TREP","nonTREP")))
            ft11 <- tryCatch(fisher.test(ct11, alternative="greater"),
                             error=function(e) list(p.value=NA_real_, estimate=NA_real_))
        }

        # =========================================================================
        # Combine all tests; BH FDR across all numeric p-values
        # =========================================================================
        stat_df <- do.call(dplyr::bind_rows, lapply(stat_rows, function(x) x))
        if (nrow(stat_df) > 0 && "p_value" %in% names(stat_df)) {
            pv  <- suppressWarnings(as.numeric(stat_df$p_value))
            ok  <- !is.na(pv)
            pbh <- rep(NA_real_, length(pv))
            if (sum(ok) > 0) pbh[ok] <- p.adjust(pv[ok], method="BH")
            stat_df$p_value_BH_FDR <- round(pbh, 6)
        }
        ref_verify_sheets[[paste0(ds, "_StatTests")]] <- stat_df
        cat(sprintf("  %s: %d statistical tests computed\n", ds, nrow(stat_df)))

        # =========================================================================
        # VERIFICATION PLOTS (VP1, VP3)
        # Plot VP1: Confirmatory gate + ordered fidelity gradient
        #   Panel A: T1 (RSk-fCV RMSE box: Orig+Ref vs Union) -- confirmatory gate
        #   Panel B: T2 (resid_abs 4-class box) -- Dunn BH FDR vs Orig+Ref in subtitle
        #   Panel C: T3 (resid_abs box, classes ordered as in T2) -- JT trend
        # Plot VP3: Secondary monotone signal characterisation
        #   Panel A: T6 (resid_abs Recovered vs Excluded box, two-sided MWU)
        #   Panel B: T7 (Hoeffding D independence (u,v) within Recovered)
        #   Panel C: T8 (distance correlation dCor (u,v) within Recovered)
        # =========================================================================
        cat(sprintf("  %s: building verification plots (VP1, VP3) ...\n", ds))

        VER_COLS <- c(
            "Original + Refined" = VP_COL_ORIGREF,
            "Original only"      = VP_COL_ORIONLY,
            "Recovered"          = VP_COL_RECOV,
            "Excluded"           = VP_COL_EXCL
        )
        CLASS_LEVS_V <- c("Original + Refined", "Recovered",
                          "Original only", "Excluded")
        plot_cells <- rc_goi[rc_goi$cell_class %in% CLASS_LEVS_V, ]
        plot_cells$cell_class_f <- factor(plot_cells$cell_class, levels = CLASS_LEVS_V)

        # =========================================================================
        # PLOT VP1: Confirmatory gate (T1) + ordered fidelity gradient (T2, T3)
        # =========================================================================

        # VP1-A: T1 -- RSk-fCV fold RMSE violin (confirmatory gate)
        t1_row <- stat_rows[["T1_confirmatory_gate"]]
        t1_p_str <- if (!is.null(t1_row)) {
            sprintf("Wilcoxon signed-rank p=%s", format(t1_row$p_value, scientific=TRUE, digits=3))
        } else "n/a"

        if (!is.null(s2f) && all(c("rmse_grp1a","rmse_final_union") %in% names(s2f))) {
            df_rskfcv_vp <- data.frame(
                rep_fold = rep(s2f$rep_fold, 2),
                RMSE     = c(s2f$rmse_grp1a, s2f$rmse_final_union),
                group = rep(c("Orig+Ref alone","Final union\n(Orig+Ref + Recovered)"),
                            each = nrow(s2f))
            )
            df_rskfcv_vp$group <- factor(df_rskfcv_vp$group,
                                         levels = c("Orig+Ref alone","Final union\n(Orig+Ref + Recovered)"))
            vp1A <- ggplot(df_rskfcv_vp, aes(x=group, y=RMSE, fill=group)) +
                geom_boxplot(width=0.45, colour="black", linewidth=VP_BOX_LW,
                             outlier.shape=NA, alpha=0.35) +
                geom_jitter(aes(colour=group), width=0.18, size=1.0, alpha=0.40,
                            shape=21, stroke=0.45,
                            show.legend=FALSE) +
                scale_fill_manual(
                    values=c("Orig+Ref alone"=VP_COL_ORIGREF_T1,
                             "Final union\n(Orig+Ref + Recovered)"=VP_COL_UNION),
                    guide="none") +
                scale_colour_manual(
                    values=c("Orig+Ref alone"=VP_COL_ORIGREF_T1,
                             "Final union\n(Orig+Ref + Recovered)"=VP_COL_UNION),
                    guide="none") +
                labs(title="A: T1 -- RSk-fCV isotonic RMSE (confirmatory gate)",
                     subtitle=sprintf("T1 (%s) | %d repeats x %d folds = %d total fold-RMSEs", t1_p_str,
                                      CMRF_N_REPEATS, CMRF_K_FOLDS,
                                      CMRF_N_REPEATS * CMRF_K_FOLDS),
                     x=NULL, y="Per-fold isotonic RMSE") +
                theme_vp() +
                theme(axis.text.x=element_text(size=8.5, lineheight=0.85))
        } else {
            vp1A <- ggplot() +
                annotate("text", x=0.5, y=0.5, label="RSk-fCV data not available", size=3.5) +
                labs(title="A: T1 RSk-fCV RMSE") + theme_void() + theme_vp()
        }
        vp1A <- vp_tag(vp1A, ds)

        # VP1-B: T2, resid_abs by class; the subtitle carries the omnibus result.
        t2_row  <- stat_rows[["T2_resid_multigroup"]]
        t2_sub  <- if (!is.null(t2_row)) {
            om_p <- fmt_p(t2_row$p_value)
            sprintf("%s omnibus across 4 classes: p %s%s", t2_row$test,
                    ifelse(grepl("^<", om_p), "", "= "), om_p)
        } else "n/a"

        vp1B <- ggplot(plot_cells, aes(x=cell_class_f, y=resid_abs, fill=cell_class_f)) +
            geom_boxplot(width=0.45, colour="black", linewidth=VP_BOX_LW,
                         outlier.shape=NA, alpha=0.35) +
            geom_jitter(aes(colour=cell_class_f), width=0.18, size=1.0, alpha=0.18,
                        shape=21, stroke=0.45,
                        show.legend=FALSE) +
            geom_hline(yintercept=bp$delta_base, colour=VP_DLINE_COL,
                       linetype="dashed", linewidth=VP_DLINE_LW) +
            annotate("text", x=Inf, y=bp$delta_base*1.06,
                     label=sprintf("delta_base\n%.4f", bp$delta_base),
                     colour=VP_DLINE_COL, size=2.3, hjust=1.05, vjust=0,
                     fontface="italic") +
            scale_fill_manual(values=VER_COLS, guide="none") +
            scale_colour_manual(values=VER_COLS, guide="none") +
            scale_x_discrete(labels=c("Original + Refined"="Orig+Ref",
                                      "Recovered"="Recov",
                                      "Original only"="Orig-only",
                                      "Excluded"="Excluded")) +
            labs(title="B: T2 -- Spine residuals by class",
                 subtitle=t2_sub,
                 x=NULL, y="|v - g(u)| (spine residual)") +
            theme_vp() +
            theme(axis.text.x=element_text(size=8))
        vp1B <- vp_tag(vp1B, ds)

        # VP1-C: T3 -- resid_abs boxplot, classes ordered matching T2 (JT trend)
        t3_row <- stat_rows[["T3_JT_ordered_trend"]]
        jt_floor <- 1 / CMRF_JT_NPERM
        jt_p_disp <- if (!is.null(t3_row) && !is.na(t3_row$p_value)) {
            if (t3_row$p_value <= jt_floor)
                sprintf("p < %s (perm floor, nperm=%d)",
                        formatC(jt_floor, format="e", digits=0), CMRF_JT_NPERM)
            else fmt_p(t3_row$p_value)
        } else "n/a"
        t3_str <- if (!is.null(t3_row) && !is.na(t3_row$JT_statistic)) {
            sprintf("Jonckheere-Terpstra JT=%.0f, %s", t3_row$JT_statistic, jt_p_disp)
        } else "JT n/a"
        # Median |resid| per class.
        med_fn <- function(cl) {
            v <- plot_cells$resid_abs[plot_cells$cell_class == cl]
            v <- v[!is.na(v)]
            if (length(v) == 0) return(NA_real_)
            median(v)
        }
        t3_sub <- sprintf("%s\nMedian |resid|:  Orig+Ref %.3f | Recov %.3f | Orig-only %.3f | Excl %.3f",
                          t3_str,
                          med_fn("Original + Refined"), med_fn("Recovered"),
                          med_fn("Original only"),      med_fn("Excluded"))

        jt_class_levs <- c("Original + Refined","Recovered","Original only","Excluded")
        df_jt_bar <- do.call(rbind, lapply(jt_class_levs, function(cl) {
            v <- plot_cells$resid_abs[plot_cells$cell_class == cl]
            v <- v[!is.na(v)]
            if (length(v) == 0) return(NULL)
            data.frame(class=cl, med_resid=median(v), n=length(v),
                       q25=quantile(v,0.25), q75=quantile(v,0.75),
                       stringsAsFactors=FALSE)
        }))
        if (!is.null(df_jt_bar)) {
            df_jt_bar$class_f <- factor(df_jt_bar$class, levels=jt_class_levs)
            JT_BAR_COLS <- c("Excluded"=VP_COL_EXCL,"Original only"=VP_COL_ORIONLY,
                             "Recovered"=VP_COL_RECOV,"Original + Refined"=VP_COL_ORIGREF)
            # Boxplot ordered as in T2 (Orig+Ref first, Excluded last).
            df_jt_box <- plot_cells[plot_cells$cell_class %in% jt_class_levs, ]
            df_jt_box$class_f <- factor(df_jt_box$cell_class, levels=jt_class_levs)
            n_lab <- setNames(df_jt_bar$n, df_jt_bar$class)
            vp1C <- ggplot(df_jt_box, aes(x=class_f, y=resid_abs, fill=class_f)) +
                geom_boxplot(width=0.45, colour="black", linewidth=VP_BOX_LW,
                             outlier.shape=NA, alpha=0.35) +
                geom_jitter(aes(colour=class_f), width=0.18, size=1.0, alpha=0.18,
                            shape=21, stroke=0.45, show.legend=FALSE) +
                geom_hline(yintercept=bp$delta_base, colour=VP_DLINE_COL,
                           linetype="dashed", linewidth=VP_DLINE_LW) +
                annotate("text", x=Inf, y=bp$delta_base*1.06,
                         label=sprintf("delta_base\n%.4f", bp$delta_base),
                         colour=VP_DLINE_COL, size=2.3, hjust=1.05, vjust=0,
                         fontface="italic") +
                scale_fill_manual(values=JT_BAR_COLS, guide="none") +
                scale_colour_manual(values=JT_BAR_COLS, guide="none") +
                scale_x_discrete(labels=function(z) sprintf("%s\nn=%d",
                                     c("Excluded"="Excl","Original only"="Orig-only",
                                       "Recovered"="Recov","Original + Refined"="Orig+Ref")[z],
                                     n_lab[z])) +
                labs(title="C: T3 -- Ordered fidelity gradient (JT)",
                     subtitle=t3_sub,
                     x="Cell class", y="|v - g(u)| (spine residual)") +
                theme_vp() +
                theme(axis.text.x=element_text(size=8, lineheight=0.85))
        } else {
            vp1C <- ggplot() +
                annotate("text", x=0.5, y=0.5, label="JT data not available", size=3.5) +
                labs(title="C: T3 JT ordered trend") + theme_void() + theme_vp()
        }
        vp1C <- vp_tag(vp1C, ds)

        p_vp1 <- cowplot::plot_grid(
            cowplot::ggdraw() + cowplot::draw_label(
                sprintf("Part10e VP1 -- Confirmatory Gate [RSk-fCV] & Ordered Fidelity Gradient: %s  (T1/T2/T3)",
                        ds),
                fontface="bold", size=11, x=0.01, hjust=0),
            cowplot::plot_grid(vp1A, vp1B, vp1C, nrow=1, ncol=3, align="hv"),
            ncol=1, rel_heights=c(0.05, 1)
        )
        vp1_base <- file.path(OUT_DIR, paste0(FILE_PREFIX, "_", ds, "_VP1_ConfGate_FidelityGrad"))
        save_pp(p_vp1, vp1_base,
                pdf_w=16, pdf_h=6, png_w=3400, png_h=1400)
        cat(sprintf("  %s: VP1 (T1/T2/T3) saved\n", ds))
        tryCatch(print(p_vp1),
                 error = function(e) cat(sprintf("  [VP1 display error: %s]\n",
                                                 conditionMessage(e))))

        # =========================================================================
        # PLOT VP3: Secondary monotone signal characterisation (T6/T7/T8)
        #   A: T6 resid Recov vs Excluded (two-sided MWU) | B: T7 Hoeffding D | C: T8 dCor
        # =========================================================================

        # VP3-A: T7 -- resid_abs Recovered vs Excluded boxplot
        t7_row <- stat_rows[["T7_resid_Recov_vs_Excluded"]]
        t7_str <- if (!is.null(t7_row) && !is.na(t7_row$p_value))
            sprintf("T6: %s p=%s", t7_row$test, fmt_p(t7_row$p_value)) else "T6: n/a"

        df_t7 <- rbind(
            if (nrow(recovered)>0)    data.frame(grp="Recovered",resid=recovered$resid_abs, stringsAsFactors=FALSE) else NULL,
            if (nrow(excl_s1fail)>0)  data.frame(grp="Excluded", resid=excl_s1fail$resid_abs, stringsAsFactors=FALSE) else NULL
        )
        if (!is.null(df_t7) && nrow(df_t7)>0) {
            df_t7$grp_f <- factor(df_t7$grp, levels=c("Recovered","Excluded"))
            vp3A <- ggplot(df_t7, aes(x=grp_f, y=resid, fill=grp_f)) +
                geom_boxplot(width=0.45, colour="black", linewidth=VP_BOX_LW,
                             outlier.shape=NA, alpha=0.35) +
                geom_jitter(aes(colour=grp_f), width=0.18, size=1.0, alpha=0.28,
                            shape=21, stroke=0.45,
                            show.legend=FALSE) +
                scale_fill_manual(values=c("Recovered"=VP_COL_RECOV,"Excluded"=VP_COL_EXCL),
                                  guide="none") +
                scale_colour_manual(values=c("Recovered"=VP_COL_RECOV,"Excluded"=VP_COL_EXCL),
                                    guide="none") +
                labs(title="A: T6 -- Spine resid: Recovered vs Excluded",
                     subtitle=t7_str,
                     x=NULL, y="|v - g(u)| (spine residual)") +
                theme_vp()
            vp3A <- vp_tag(vp3A, ds)
        } else {
            vp3A <- ggplot() +
                annotate("text", x=0.5, y=0.5, label="T7 data not available", size=3.5) +
                labs(title="A: T6 resid Recov vs Excl") + theme_void() + theme_vp()
        }

        # VP3-B: T7, Hoeffding D independence test (u,v) in Recovered
        t7_hd_row <- stat_rows[["T7_hoeffding_D_Recovered"]]
        hd_floor  <- 1 / (CMRF_HD_NPERM + 1)
        t7_hd_str <- if (!is.null(t7_hd_row) && !is.na(t7_hd_row$D_statistic)) {
            hd_p_disp <- if (!is.na(t7_hd_row$p_value) && t7_hd_row$p_value <= hd_floor)
                sprintf("p < %s (perm floor, R=%d)",
                        formatC(hd_floor, format="e", digits=0), CMRF_HD_NPERM)
            else sprintf("p=%s", format(t7_hd_row$p_value, scientific=TRUE, digits=3))
            sprintf("T7: Hoeffding D=%.4f %s (n=%d)",
                    t7_hd_row$D_statistic, hd_p_disp, nrow(recovered))
        } else "T7: n/a"

        if (nrow(recovered) > 4) {
            vp3B <- ggplot(recovered, aes(x=u_all, y=v_all)) +
                geom_point(colour=VP_COL_RECOV, size=1.2, alpha=0.55, shape=16) +
                geom_vline(xintercept=IP_u_val, colour=VP_DLINE_COL,
                           linetype="dashed", linewidth=VP_DLINE_LW) +
                geom_hline(yintercept=0.5, colour=VP_DLINE_COL,
                           linetype="dashed", linewidth=VP_DLINE_LW) +
                coord_fixed(xlim=c(0,1), ylim=c(0,1)) +
                labs(title="B: T7 -- Hoeffding D independence (Recovered)",
                     subtitle=t7_hd_str,
                     x="u = Fx(GOI)", y="v = Fc(cDFG)") +
                theme_vp() +
                theme(aspect.ratio=1,
                      panel.grid.minor = element_blank())
            vp3B <- vp_tag(vp3B, ds)
        } else {
            vp3B <- ggplot() +
                annotate("text", x=0.5, y=0.5,
                         label="T7: insufficient Recovered cells", size=3.5) +
                labs(title="B: T7 Hoeffding D") + theme_void() + theme_vp()
        }

        # VP3-C: T8, dCor(u,v) in Recovered; points coloured by spine residual.
        t8_dc_row <- stat_rows[["T8_dcor_Recovered"]]
        dc_floor  <- 1 / (999 + 1)   # dcor.test uses R=999 permutations
        t8_dc_str <- if (!is.null(t8_dc_row) && !is.na(t8_dc_row$dCor)) {
            dc_p_disp <- if (!is.na(t8_dc_row$p_value) && t8_dc_row$p_value <= dc_floor)
                sprintf("p < %s (perm floor, R=999)",
                        formatC(dc_floor, format="e", digits=0))
            else sprintf("p=%s", format(t8_dc_row$p_value, scientific=TRUE, digits=3))
            sprintf("T8: dCor=%.4f %s (n=%d)",
                    t8_dc_row$dCor, dc_p_disp, nrow(recovered))
        } else "T8: n/a"

        if (nrow(recovered) > 4) {
            # Viridis magma palette (5%-95% spectrum)
            vp3C <- ggplot(recovered, aes(x=u_all, y=v_all, colour=resid_abs)) +
                geom_point(size=1.4, alpha=0.75, shape=16) +
                scale_colour_viridis_c(option="magma", begin=0.05, end=0.95,
                                       name="|resid|", limits=c(0, NA),
                                       guide=guide_colorbar(barwidth=0.5,
                                                            barheight=3.5,
                                                            title.theme=element_text(size=7),
                                                            label.theme=element_text(size=6))) +
                geom_vline(xintercept=IP_u_val, colour=VP_DLINE_COL,
                           linetype="dashed", linewidth=VP_DLINE_LW) +
                geom_hline(yintercept=0.5, colour=VP_DLINE_COL,
                           linetype="dashed", linewidth=VP_DLINE_LW) +
                coord_fixed(xlim=c(0,1), ylim=c(0,1)) +
                labs(title="C: T8 -- Distance correlation dCor (Recovered)",
                     subtitle=t8_dc_str,
                     x="u = Fx(GOI)", y="v = Fc(cDFG)") +
                theme_vp() +
                theme(aspect.ratio=1,
                      panel.grid.minor = element_blank(),
                      legend.position="right",
                      legend.key.size=unit(0.35, "cm"))
            vp3C <- vp_tag(vp3C, ds)
        } else {
            vp3C <- ggplot() +
                annotate("text", x=0.5, y=0.5,
                         label="T8: insufficient Recovered cells", size=3.5) +
                labs(title="C: T8 Distance correlation") + theme_void() + theme_vp()
        }

        p_vp3 <- cowplot::plot_grid(
            cowplot::ggdraw() + cowplot::draw_label(
                sprintf("Part10e VP3 -- Secondary Monotone Signal: %s  (T6/T7/T8)", ds),
                fontface="bold", size=11, x=0.01, hjust=0),
            cowplot::plot_grid(vp3A, vp3B, vp3C, nrow=1, ncol=3, align="hv"),
            ncol=1, rel_heights=c(0.05, 1)
        )
        vp3_base <- file.path(OUT_DIR, paste0(FILE_PREFIX, "_", ds, "_VP3_MonotoneSignal"))
        save_pp(p_vp3, vp3_base,
                pdf_w=16, pdf_h=6, png_w=3400, png_h=1400)
        cat(sprintf("  %s: VP3 (T6/T7/T8) saved\n", ds))
        tryCatch(print(p_vp3),
                 error = function(e) cat(sprintf("  [VP3 display error: %s]\n",
                                                 conditionMessage(e))))

        # =========================================================================
        # EXCEL EXPORT: Plot data for T1--T7A (one file per dataset)
        # Allows re-plotting in external software (Prism, Python, etc.)
        # =========================================================================
        vp_excel_sheets <- list()

        # T1: RSk-fCV fold RMSE (Orig+Ref vs Final union)
        if (exists("df_rskfcv_vp") && !is.null(df_rskfcv_vp) && nrow(df_rskfcv_vp) > 0) {
            vp_excel_sheets[["T1_RSk_fCV_RMSE"]] <- df_rskfcv_vp
        }

        # T2: Spine residuals by 4-class (cell-level)
        if (exists("plot_cells") && !is.null(plot_cells) && nrow(plot_cells) > 0) {
            df_t2_export <- data.frame(
                cell_class = as.character(plot_cells$cell_class_f),
                resid_abs  = plot_cells$resid_abs,
                stringsAsFactors = FALSE
            )
            if ("barcode" %in% names(plot_cells)) df_t2_export$barcode <- plot_cells$barcode
            if ("u_all"   %in% names(plot_cells)) df_t2_export$u       <- plot_cells$u_all
            if ("v_all"   %in% names(plot_cells)) df_t2_export$v       <- plot_cells$v_all
            vp_excel_sheets[["T2_Resid_by_class"]] <- df_t2_export
        }

        # T3: JT ordered fidelity gradient (summary per class)
        if (exists("df_jt_bar") && !is.null(df_jt_bar) && nrow(df_jt_bar) > 0) {
            vp_excel_sheets[["T3_JT_ordered_trend"]] <- df_jt_bar
        }

        # T7A (T6 in stat_rows): Spine resid Recovered vs Excluded (cell-level)
        if (exists("df_t7") && !is.null(df_t7) && nrow(df_t7) > 0) {
            vp_excel_sheets[["T7A_Resid_Recov_vs_Excl"]] <- df_t7
        }

        if (length(vp_excel_sheets) > 0) {
            vp_xlsx_path <- file.path(OUT_DIR,
                                      paste0(FILE_PREFIX, "_", ds, "_VP_PlotData_T1_to_T7A.xlsx"))
            tryCatch({
                write_xlsx(vp_excel_sheets, path = vp_xlsx_path)
                cat(sprintf("  %s: VP plot data Excel saved: %s\n", ds,
                            basename(vp_xlsx_path)))
            }, error = function(e)
                cat("  ", ds, ": VP plot data Excel FAILED: ",
                    conditionMessage(e), "\n", sep = ""))
        }

        # Store the plots
        all_plots_verify[[ds]] <- list(vp1=p_vp1, vp3=p_vp3)


    }

    # == Write verification Excel ============================================
    ref_verify_path <- file.path(OUT_DIR, paste0(FILE_PREFIX, "_Verification.xlsx"))
    tryCatch({
        write_xlsx(ref_verify_sheets, path = ref_verify_path)
        cat(sprintf("  CMRF verification Excel saved: %s\n",
                    basename(ref_verify_path)))
    }, error = function(e)
        cat("  CMRF verification Excel FAILED: ", conditionMessage(e),
            "\n  (close the workbook if it is open, then rerun)\n", sep = ""))


    ##########################################################################
    # SECTION C1: POST-ANALYSIS GUARDRAILS
    ##########################################################################
    cat("\n", strrep("#", 72), "\n")
    cat("  SECTION C1: POST-ANALYSIS GUARDRAILS\n")
    cat(strrep("#", 72), "\n\n")

    for (ds in names(all_results_10e)) {
        r <- all_results_10e[[ds]]
        if (is.null(r)) {
            gr10e(paste0("GR-08-", ds), "FAIL", "Dataset returned NULL")
            next
        }

        cc <- r$cmrf_class_counts
        gr10e(paste0("GR-08-", ds),
              if (!is.null(cc)) "PASS" else "WARN",
              if (!is.null(cc))
                  paste0("CMRF classification: ",
                         paste(sprintf("%s=%d", names(cc), as.integer(cc)),
                               collapse = " | "))
              else "CMRF classification: not available")

        gr10e(paste0("GR-09-", ds),
              if (!is.null(r$n_recovered) && r$n_recovered >= 0) "PASS" else "WARN",
              sprintf("Recovered cells = %s | original SKCS-core = %s | refined = %s",
                      if (is.null(r$n_recovered)) "NA" else as.character(r$n_recovered),
                      if (is.null(r$n_orig_skcs)) "NA" else as.character(r$n_orig_skcs),
                      if (is.null(r$n_refined))   "NA" else as.character(r$n_refined)))

        # == GR-10: spine monotonicity ============================================
        sm <- r$spine_monotone
        gr10e(paste0("GR-10-", ds),
              if (is.null(sm)) "WARN" else if (isTRUE(sm$monotone)) "PASS" else "WARN",
              if (is.null(sm)) "spine not fitted, so monotonicity is unknown"
              else sprintf(paste0("spine %s: min consecutive difference = %.3e ",
                                  "over %d grid points%s"),
                           if (isTRUE(sm$monotone)) "non-decreasing" else "DECREASING",
                           sm$min_diff, sm$n_grid,
                           if (isTRUE(sm$monotone)) ""
                           else sprintf(" | steepest decrease near u = %.4f",
                                        sm$worst_u)))

        # == GR-11: donor random intercept and fitting set ========================
        sp <- r$spine_prss
        gr10e(paste0("GR-11-", ds),
              if (is.null(sp)) "WARN"
              else if (isTRUE(sp$donor_re) || !isTRUE(GAM_DONOR_RE_10E)) "PASS"
              else "WARN",
              if (is.null(sp)) "spine not fitted"
              else sprintf(paste0("spine on %s TREP SKCS-core cell(s), u in ",
                                  "[%s, %s] | donor RE = %s (%d donors) | ",
                                  "k = %d | %d iteration(s)"),
                           if (is.null(r$spine_n_cells)) "NA"
                           else as.character(r$spine_n_cells),
                           if (is.null(r$spine_u_range)) "NA"
                           else sprintf("%.4f", r$spine_u_range[1]),
                           if (is.null(r$spine_u_range)) "NA"
                           else sprintf("%.4f", r$spine_u_range[2]),
                           if (isTRUE(sp$donor_re)) "YES" else "NO",
                           sp$n_donors, sp$best_k, sp$actual_iters))

        # == GR-12: axis of the inherited Pre/Post labels =========================
        # The cDFG rule must reproduce Part 10b's labels; fields are length-checked.
        .num1 <- function(z, digits = 6L) if (length(z) != 1L || !is.finite(z))
            "NA" else sprintf(paste0("%.", digits, "f"), z)
        .ok1  <- function(z) length(z) == 1L && is.finite(z)
        .ax <- r$ip_axis_check
        if (!is.list(.ax) || !.ok1(.ax$n) || .ax$n == 0L ||
            !.ok1(.ax$agree_cdfg)) {
            gr10e(paste0("GR-12-", ds), "WARN",
                  sprintf(paste0("Pre/Post labels not comparable: %s | ",
                                 "IP_y = %s | IP_v = %s"),
                          if (!is.list(.ax) || length(.ax$note) != 1L ||
                              is.na(.ax$note) || !nzchar(.ax$note))
                              "no IP_status or no cDFG_expr to compare against"
                          else .ax$note,
                          .num1(r$IP_y), .num1(r$IP_v)))
        } else {
            gr10e(paste0("GR-12-", ds),
                  if (.ax$agree_cdfg >= 0.999) "PASS" else "FAIL",
                  sprintf(paste0("classification axis = %s | IP_y = %s | ",
                                 "IP_v = %s | Part 10b Pre/Post labels ",
                                 "reproduce %.4f on the cDFG axis at IP_y and ",
                                 "%s on the GOI axis at IP_x, over %d cell(s)%s"),
                          if (length(.ax$axis) == 1L) .ax$axis else "cDFG at IP_y",
                          .num1(r$IP_y), .num1(r$IP_v),
                          .ax$agree_cdfg,
                          .num1(.ax$agree_goi, 4L),
                          as.integer(.ax$n),
                          if (.ax$agree_cdfg < 0.999)
                              " -- the EP CSV was not written on the cDFG axis"
                          else ""))
        }
    }

    gr10e_write()
    cat(sprintf("\n  Final guardrail tally: %d PASS / %d WARN / %d FAIL\n",
                .gr10e_pass, .gr10e_warn, .gr10e_fail))
    if (.gr10e_fail > 0) {
        cat("  Failures present -- review the ledger before interpreting.\n")
    } else if (.gr10e_warn > 0) {
        cat(sprintf("  %d passed, %d warning(s), no failures.\n",
                    .gr10e_pass, .gr10e_warn))
    } else {
        cat(sprintf("  All %d checks passed.\n", .gr10e_pass))
    }


    ##########################################################################
    # SECTION C2: ARM SUMMARY
    ##########################################################################
    cat("\n", strrep("=", 72), "\n")
    cat("  Part10e CMRF COMPLETE [", ARM_TAG, "]\n")
    cat("  Output directory:", OUT_DIR, "\n")
    cat("  Files generated:\n")
    .transient <- c("_RUN_INCOMPLETE.txt", basename(MSG_PATH_P10E))
    for (f in setdiff(sort(list.files(OUT_DIR)), .transient)) {
        sz <- round(file.size(file.path(OUT_DIR, f)) / 1e3, 0)
        cat(sprintf("    %-62s  %5d KB\n", f, sz))
    }
    cat("  Load the results with: readRDS(\"", rds_path, "\")\n", sep = "")
    cat(strrep("=", 72), "\n")

    cat("\n  Log written to:", LOG_PATH_P10E,
        if (isTRUE(LOG_MESSAGES_P10E)) " (console output and warnings)"
        else " (console output only)", "\n")
    cat(sprintf("\nPart10e CMRF - log ended: %s\n", Sys.time()))

    .p10e_completed <- TRUE

    if (sink.number(type = "message") != 2L) sink(type = "message")
    while (sink.number() > 0) sink()

    suppressWarnings(file.remove(INCOMPLETE_MARKER_P10E))

    invisible(list(arm  = arm_tag, dir = OUT_DIR,
                   pass = .gr10e_pass, warn = .gr10e_warn,
                   fail = .gr10e_fail))
}


# =============================================================================
# SECTION 0c: Driver
# =============================================================================

cep_run_arms(PART_DIR_TAG, run_one_arm)
