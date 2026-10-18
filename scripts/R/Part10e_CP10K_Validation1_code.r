# =============================================================================
# Part10e_CP10K_Validation1_code.r
# SEA-AD CEP-IP/CMRF, Part 10e, GAM refit on the final Module (Figure 8A, Table S10)
#   For each arm: refits the cDFG GAM (donor random intercept, inverse-donor
#   weights) on the CMRF final Module, fits the cHKG curve on the same cells as
#   a within-arm reference, recomputes EP and TREP, and reports concordance
#   with the Part 10b TREP labels. Draws the CEP-IP, Pitcher and copula-space
#   figures with refitted and original labels.
#
# Usage (from the SEA-AD folder, or set SEAAD_ROOT), after Part10e_CP10K_CMRF_code.r
#   source("Part10e_CP10K_Validation1_code.r")                         # all arms
#   CEP_ARM <- "PVM1_all"; source("Part10e_CP10K_Validation1_code.r")  # one arm
#   Needs CEP_arms.R, R_config.R and R_config_HKG8.R in the same folder.
#
# Outputs (Part10e_Validation1/Part10e_Validation1_<arm>/): figures (PDF, PNG),
#   workbook and results RDS per dataset, cell-count report, guardrail ledger
#   and logs.
# =============================================================================


# =============================================================================
# SECTION 0: Libraries and shared arm definitions
# =============================================================================

library(mgcv)
library(Matrix)        # readMM for norm.mtx and the guardrails
library(ggplot2)
library(dplyr)
library(readxl)
library(writexl)
library(gridExtra)
library(grid)
library(scales)

set.seed(42)

# CEP_arms.R holds the arm table, the folder layout and the driver.
if (!exists("CEP_HOME")) CEP_HOME <- Sys.getenv("SEAAD_ROOT", getwd())
.p10v_harness <- file.path(CEP_HOME, "CEP_arms.R")
if (!file.exists(.p10v_harness))
    stop("CEP_arms.R not found at: ", .p10v_harness,
         "\n  It holds the arm table shared by Parts 9 onwards.", call. = FALSE)
source(.p10v_harness)

PART_TAG     <- "Part10e"              # file-name tag; "_V1" is appended
PART_DIR_TAG <- "Part10e_Validation1"  # results-folder tag

# == Estimator settings, as in Parts 10b and 10e ==============================
GAM_DONOR_RE_V1        <- TRUE
DONOR_RE_MIN_LEVELS_V1 <- 3L
W_DONOR_WEIGHTS_V1     <- TRUE
GAM_GAMMA_V1           <- 1.5
PRSS_PHASE2_MAX_V1     <- 5L    # phase 1 = 8 iterations, phase 2 = at most 5
# Phase 1 of both PRSS loops uses workers only if the timed first fit predicts
# at least P10V_PAR_MIN_SECONDS for the remaining seven fits.
P10V_PARALLEL          <- TRUE
P10V_N_CORES           <- max(1L, parallel::detectCores() - 1L)
P10V_PAR_MIN_SECONDS   <- 2.0
PRSS_EARLY_STOP_V1     <- 20L   # unreachable under the phase-2 cap
MIN_CELLS              <- 50L
IP_CLASS_AXIS_V1       <- "cDFG"   # "GOI" splits Pre/Post on the GOI axis

DONOR_RE_LABEL_V1 <- "s(donor)"
GOI_SMOOTH_LABEL  <- "s(GOI)"


# =============================================================================
# SECTION 0d: Pure helpers
# =============================================================================
# Depend only on their arguments; the config is sourced inside run_one_arm().

# Donor as a factor, and inverse-donor weights w = 1 / (cells per donor).
as_donor_factor_v1 <- function(x) droplevels(factor(as.character(x)))

n_donor_levels_v1 <- function(x) {
    if (is.null(x)) return(0L)
    nlevels(as_donor_factor_v1(x))
}

donor_re_ok_v1 <- function(x)
    isTRUE(GAM_DONOR_RE_V1) && n_donor_levels_v1(x) >= DONOR_RE_MIN_LEVELS_V1

donor_weights_v1 <- function(donor_fac) {
    if (!isTRUE(W_DONOR_WEIGHTS_V1)) return(rep(1, length(donor_fac)))
    f <- as_donor_factor_v1(donor_fac)
    as.numeric(1 / as.numeric(table(f)[f]))
}

p10v_gam_formula <- function(k, use_re, bs) {
    rhs <- sprintf("GOI + s(GOI, bs = \"%s\", k = %d)", bs, as.integer(k))
    if (use_re) rhs <- paste0(rhs, " + s(donor, bs = \"re\")")
    stats::as.formula(paste("Expression ~", rhs))
}

# The s(GOI) term located by label, and its lambda by first.sp.
p10v_smooth_idx <- function(model) {
    if (is.null(model$smooth) || length(model$smooth) == 0) return(NA_integer_)
    labs <- vapply(model$smooth, function(s) s$label, character(1))
    i <- which(labs == GOI_SMOOTH_LABEL)
    if (length(i) == 1L) return(i)
    i <- which(labs != DONOR_RE_LABEL_V1)
    if (length(i) >= 1L) return(i[1])
    NA_integer_
}

p10v_smooth_lambda <- function(model, i) {
    if (is.na(i) || is.null(model$sp) || length(model$sp) == 0) return(NA_real_)
    sm <- model$smooth[[i]]
    if (!is.null(sm$first.sp) && sm$first.sp >= 1L &&
        sm$first.sp <= length(model$sp))
        return(unname(model$sp[sm$first.sp]))
    hit <- which(startsWith(names(model$sp), sm$label))
    if (length(hit) > 0) return(unname(model$sp[hit[1]]))
    unname(model$sp[1])
}

# PRSS = RSS + lambda * f'' penalty of s(GOI); the data argument is unused.
.calculate_prss <- function(model, data = NULL) {
    RSS <- sum(stats::residuals(model)^2)
    i   <- p10v_smooth_idx(model)
    if (is.na(i))
        return(list(PRSS = RSS, RSS = RSS,
                    f_double_prime_integral = NA_real_, lambda = NA_real_))
    sm     <- model$smooth[[i]]
    lambda <- p10v_smooth_lambda(model, i)
    if (is.null(sm$S) || length(sm$S) == 0)
        return(list(PRSS = RSS, RSS = RSS,
                    f_double_prime_integral = NA_real_, lambda = lambda))
    beta  <- stats::coef(model)[sm$first.para:sm$last.para]
    f_dpi <- as.numeric(t(beta) %*% sm$S[[1]] %*% beta)
    list(PRSS = RSS + lambda * f_dpi, RSS = RSS,
         f_double_prime_integral = f_dpi, lambda = lambda)
}

# Population curve: a reference donor level is supplied and s(donor) excluded.
p10v_predict_pop <- function(model, goi_vals, se.fit = FALSE) {
    nd <- data.frame(GOI = as.numeric(goi_vals))
    has_re <- !is.null(model$smooth) &&
        any(vapply(model$smooth, function(s) s$label, character(1)) ==
                DONOR_RE_LABEL_V1)
    if (has_re) {
        ref <- model$var.summary[["donor"]]
        if (is.null(ref)) ref <- model$xlevels[["donor"]]
        ref <- if (is.factor(ref)) levels(ref)[1] else as.character(ref)[1]
        nd$donor <- factor(ref, levels = levels(model$model$donor))
        return(stats::predict(model, newdata = nd,
                              exclude = DONOR_RE_LABEL_V1,
                              type = "response", se.fit = se.fit))
    }
    stats::predict(model, newdata = nd, type = "response", se.fit = se.fit)
}

# Weights for both deviance measures: the frame's w, else the model's prior
# weights, else inverse-donor weights, else 1.
ep_weight_vec_v1 <- function(model, data) {
    if ("w" %in% names(data)) return(as.numeric(data$w))
    if (!is.null(model) && !is.null(model$prior.weights) &&
        length(model$prior.weights) == nrow(data))
        return(as.numeric(model$prior.weights))
    if ("donor" %in% names(data)) return(donor_weights_v1(data$donor))
    rep(1, nrow(data))
}

# Does the model carry the donor random intercept.
p10v_model_has_re <- function(model)
    !is.null(model$smooth) &&
    any(vapply(model$smooth, function(s) s$label, character(1)) ==
            DONOR_RE_LABEL_V1)

# TREP sizing DE: weighted 1 - RSS/TSS against a flat mean, from the
# population prediction (s(donor) excluded). Used only to size the TREP set.
pop_dev_explained_v1 <- function(model, data) {
    if (is.null(model) || is.null(data)) return(NA_real_)
    if (!all(c("GOI", "Expression") %in% names(data))) return(NA_real_)
    y   <- as.numeric(data$Expression)
    fit <- as.numeric(p10v_predict_pop(model, data$GOI))
    if (length(fit) != length(y)) return(NA_real_)
    w   <- ep_weight_vec_v1(model, data)
    if (length(w) != length(y)) w <- rep(1, length(y))
    ok  <- is.finite(y) & is.finite(fit) & is.finite(w) & w > 0
    if (sum(ok) < 2L) return(NA_real_)
    ybar   <- sum(w[ok] * y[ok]) / sum(w[ok])
    ss_tot <- sum(w[ok] * (y[ok] - ybar)^2)
    if (!is.finite(ss_tot) || ss_tot <= 0) return(NA_real_)
    as.numeric(1 - sum(w[ok] * (y[ok] - fit[ok])^2) / ss_tot)
}

# Scalar or NA, so a log line built from a NULL field does not vanish.
.num1_v1 <- function(x) {
    if (is.null(x) || length(x) != 1L) return(NA_real_)
    suppressWarnings(as.numeric(x))
}

.p10v_null_dev_cache <- new.env(hash = TRUE, parent = emptyenv())

get_null_deviance_v1 <- function(data, w, use_re, method = "REML") {
    d <- data.frame(Expression = as.numeric(data$Expression),
                    w          = as.numeric(w))
    if (use_re) d$donor <- as_donor_factor_v1(data$donor)
    # Cache key separates frames, including how y is paired with donor.
    .yy <- as.numeric(d$Expression)
    key <- paste0(nrow(d), "_", round(sum(.yy), 6), "_",
                  round(sum(.yy^2), 6), "_",
                  round(sum(seq_along(.yy) * .yy), 6), "_",
                  if (use_re)
                      paste0(nlevels(d$donor), "_",
                             round(sum(as.integer(d$donor) * .yy), 6))
                  else "0_0", "_",
                  as.integer(use_re), "_", round(sum(d$w^2), 9), "_", method)
    if (exists(key, envir = .p10v_null_dev_cache))
        return(get(key, envir = .p10v_null_dev_cache))
    set.seed(123)
    # weights = w by bare column name, so w must be a column of d.
    nd <- if (use_re)
        stats::deviance(mgcv::gam(Expression ~ 1 + s(donor, bs = "re"),
                                  data = d, weights = w, method = method))
    else
        stats::deviance(mgcv::gam(Expression ~ 1, data = d, weights = w))
    assign(key, nd, envir = .p10v_null_dev_cache)
    nd
}

# Reported DE (donor-aware): the null carries the same donor random intercept
# and weights as the model, fitted with the model's own method.
de_explained_v1 <- function(model, data) {
    if (is.null(model) || is.null(data)) return(NA_real_)
    if (!("Expression" %in% names(data))) return(NA_real_)
    w <- ep_weight_vec_v1(model, data)
    if (length(w) != nrow(data)) return(NA_real_)
    use_re <- p10v_model_has_re(model) && ("donor" %in% names(data)) &&
        donor_re_ok_v1(data$donor)
    method <- if (!is.null(model$method)) model$method else "REML"
    nd <- tryCatch(get_null_deviance_v1(data, w, use_re, method),
                   error = function(e) NA_real_)
    if (!is.finite(nd) || nd <= 0) return(NA_real_)
    dm <- stats::deviance(model)
    if (!is.finite(dm)) return(NA_real_)
    as.numeric((nd - dm) / nd)
}

# Pre/Post boundary layers: IP_y in expression space, IP_v in copula space;
# nothing is drawn when the value is not available.
ip_hline_v1 <- function(y_val, x_lim, size = 2.8, colour = "black") {
    if (!is.finite(y_val)) return(list())
    list(geom_hline(yintercept = y_val, linetype = "dashed",
                    colour = colour, linewidth = 0.7),
         annotate("text", x = x_lim[1] + diff(x_lim) * 0.02,
                  y = y_val, vjust = -0.6,
                  label = sprintf("y[IP]==%.3f", y_val),
                  parse = TRUE, hjust = 0, size = size, colour = colour))
}

ip_hline_cop_v1 <- function(v_val, colour = "gray10") {
    if (!is.finite(v_val)) return(list())
    list(geom_hline(yintercept = v_val, linetype = "dashed",
                    colour = colour, linewidth = 0.4),
         annotate("text", x = 0.02, y = v_val, vjust = -0.6,
                  label = sprintf("v[IP]==%.3f", v_val),
                  parse = TRUE, hjust = 0, size = 2.8, colour = colour))
}

# TREP sizing fraction, clamped to [0, 1].
trep_count_fraction_v1 <- function(de) {
    if (!is.finite(de)) return(NA_real_)
    max(0, min(1, de))
}

# == PRSS fitting loop ========================================================
# Phase 1: 8 fits over k; phase 2: at most PRSS_PHASE2_MAX_V1 refinements
# around the best k, with the k draw seeded from (seed, i).
.fit_full_gam <- function(data, label,
                          donor      = NULL,
                          phase2_max = PRSS_PHASE2_MAX_V1,
                          gamma      = GAM_GAMMA_V1,
                          verbose    = TRUE) {
    set.seed(123)
    unique_x <- length(unique(data$GOI))
    max_k    <- min(10L, unique_x - 1L)
    if (unique_x < 4L) {
        if (verbose) cat("    Skipping: < 4 unique GOI values\n")
        return(NULL)
    }

    donor_f <- if (is.null(donor)) NULL else as_donor_factor_v1(donor)
    use_re  <- donor_re_ok_v1(donor_f)
    n_don   <- n_donor_levels_v1(donor_f)
    df_fit  <- data.frame(GOI = as.numeric(data$GOI),
                          Expression = as.numeric(data$Expression))
    if (!is.null(donor_f)) {
        df_fit$donor <- donor_f
        df_fit$w     <- donor_weights_v1(donor_f)
    } else {
        df_fit$w <- rep(1, nrow(df_fit))
    }
    if (verbose)
        cat(sprintf(paste0("    [%s] n=%d | unique_x=%d | max_k=%d | donor RE=%s",
                           " (%d donors) | weights=%s\n"),
                    label, nrow(df_fit), unique_x, max_k,
                    if (use_re) "YES" else "NO", n_don,
                    if (isTRUE(W_DONOR_WEIGHTS_V1)) "inverse donor" else "uniform"))

    n_iter_max  <- 8L + as.integer(phase2_max)
    prss_values <- rep(NA_real_, n_iter_max)
    best_prss   <- Inf; best_model <- NULL; best_k <- 3L; best_iter <- 1L
    early_ctr   <- 0L; early_stop <- FALSE

    run_one <- function(i, k) {
        set.seed(123L + i)
        # weights = w by bare column name.
        m <- mgcv::gam(p10v_gam_formula(k, use_re, GAM_BS),
                       data = df_fit, weights = w, method = GAM_METHOD,
                       select = GAM_SELECT, gamma = gamma)
        list(model = m, prss = .calculate_prss(m)$PRSS, k = k)
    }

    k_phase1 <- vapply(seq_len(8L), function(i)
        min(3L + ((i - 1L) %% max(1L, max_k - 2L)), max_k), numeric(1))

    # Fit 1 runs here and is timed; workers are used for fits 2-8 if worthwhile.
    run_one_timed <- function(i, k) {
        t0  <- proc.time()[["elapsed"]]
        res <- tryCatch(run_one(i, k), error = function(e) {
            if (verbose) cat("      error iter", i, ":", conditionMessage(e), "\n")
            NULL
        })
        if (!is.null(res)) res$secs <- proc.time()[["elapsed"]] - t0
        res
    }

    phase1_res <- vector("list", 8L)
    phase1_res[[1L]] <- run_one_timed(1L, k_phase1[1L])
    t_fit1 <- if (is.null(phase1_res[[1L]])) NA_real_ else phase1_res[[1L]]$secs
    if (verbose && !is.null(phase1_res[[1L]]))
        cat(sprintf("      init  iter  1 | k=%d | PRSS=%.4f | %.1f s\n",
                    phase1_res[[1L]]$k, phase1_res[[1L]]$prss, t_fit1))

    rest      <- 2:8
    n_workers <- if (isTRUE(P10V_PARALLEL))
        max(1L, min(as.integer(P10V_N_CORES), length(rest))) else 1L
    est_rest  <- if (is.finite(t_fit1)) t_fit1 * length(rest) else 0
    used_par  <- FALSE

    if (n_workers > 1L && est_rest >= P10V_PAR_MIN_SECONDS) {
        cl <- tryCatch(parallel::makePSOCKcluster(n_workers),
                       error = function(e) NULL)
        if (!is.null(cl)) {
            on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
            par_res <- tryCatch({
                parallel::clusterEvalQ(cl, library(mgcv))
                parallel::clusterExport(
                    cl,
                    varlist = c("df_fit", "use_re", "k_phase1", "gamma",
                                "p10v_gam_formula", ".calculate_prss",
                                "p10v_smooth_idx", "p10v_smooth_lambda",
                                "GOI_SMOOTH_LABEL", "DONOR_RE_LABEL_V1",
                                "GAM_BS", "GAM_METHOD", "GAM_SELECT"),
                    envir = environment())
                parallel::parLapply(cl, rest, function(i) {
                    kk <- k_phase1[i]
                    t0 <- proc.time()[["elapsed"]]
                    set.seed(123L + i)
                    m <- tryCatch(
                        mgcv::gam(p10v_gam_formula(kk, use_re, GAM_BS),
                                  data = df_fit, weights = w,
                                  method = GAM_METHOD, select = GAM_SELECT,
                                  gamma = gamma),
                        error = function(e) NULL)
                    if (is.null(m)) return(NULL)
                    list(model = m, prss = .calculate_prss(m)$PRSS, k = kk,
                         secs = proc.time()[["elapsed"]] - t0)
                })
            }, error = function(e) {
                if (verbose)
                    cat(sprintf("      parallel phase 1 failed (%s); sequential.\n",
                                conditionMessage(e)))
                NULL
            })
            if (!is.null(par_res)) {
                phase1_res[rest] <- par_res
                used_par <- TRUE
                if (verbose)
                    cat(sprintf("      fits 2-8 ran on %d worker(s)\n", n_workers))
            }
        }
    } else if (n_workers > 1L && verbose) {
        cat(sprintf(paste0("      staying sequential: one fit %.1f s, remaining 7",
                           " about %.1f s, below %.1f s\n"),
                    t_fit1, est_rest, P10V_PAR_MIN_SECONDS))
    }

    if (!used_par)
        for (i in rest) phase1_res[[i]] <- run_one_timed(i, k_phase1[i])

    for (i in seq_len(8L)) {
        res <- phase1_res[[i]]
        if (!is.null(res)) {
            prss_values[i] <- res$prss
            if (verbose && i > 1L)
                cat(sprintf("      init  iter %2d | k=%d | PRSS=%.4f | %.1f s\n",
                            i, res$k, res$prss,
                            if (is.null(res$secs)) NA_real_ else res$secs))
            if (res$prss < best_prss) {
                best_prss <- res$prss; best_model <- res$model
                best_k <- res$k; best_iter <- i; early_ctr <- 0L
            } else early_ctr <- early_ctr + 1L
        } else early_ctr <- early_ctr + 1L
    }
    if (is.null(best_model)) {
        if (verbose) cat(sprintf("    [%s]: every phase-1 iteration failed.\n", label))
        return(NULL)
    }

    k_rng <- function(bk) unique(pmax(3L, pmin(max_k,
        if (bk >= max_k) c(max(3L, max_k - 1L), max_k)
        else if (bk <= 3L) c(3L, min(4L, max_k))
        else c(bk - 1L, bk, min(bk + 1L, max_k)))))
    k_range <- k_rng(best_k)

    # Phase 2: refinement around the best k.
    for (i in 8L + seq_len(max(0L, n_iter_max - 8L))) {
        if (early_ctr >= PRSS_EARLY_STOP_V1) { early_stop <- TRUE; break }
        set.seed(200L + i)
        k   <- if (i %% 3L == 0L) sample(k_range, 1L) else best_k
        res <- tryCatch(run_one(i, k), error = function(e) NULL)
        if (!is.null(res)) {
            prss_values[i] <- res$prss
            if (res$prss < best_prss) {
                best_prss <- res$prss; best_model <- res$model
                best_k <- k; best_iter <- i; early_ctr <- 0L
                k_range <- k_rng(best_k)
            } else early_ctr <- early_ctr + 1L
        } else early_ctr <- early_ctr + 1L
    }
    if (verbose)
        cat(sprintf("    Final best k=%d PRSS=%.4f  iters=%d  early_stop=%s\n",
                    best_k, best_prss, n_iter_max, early_stop))

    list(best_model = best_model, best_k = best_k, best_prss = best_prss,
         best_iter = best_iter, prss_values = prss_values,
         n_iter_max = n_iter_max, phase2_max = as.integer(phase2_max),
         early_stopped = early_stop, donor_re = use_re, n_donors = n_don,
         data = data, label = label)
}
# == Plot theme and figure savers =============================================
theme_p10b <- function(base_size = 10) {
    theme_minimal(base_size = base_size) +
        theme(
            plot.title       = element_text(face = "bold", size = base_size + 1),
            plot.subtitle    = element_text(size = base_size - 1),
            axis.title       = element_text(face = "bold", size = base_size - 1),
            axis.text        = element_text(size = base_size - 2),
            panel.grid.major = element_line(color = "gray90"),
            panel.grid.minor = element_blank(),
            panel.border     = element_rect(color = "black", fill = NA,
                                            linewidth = 0.35),
            plot.background  = element_rect(fill = "white", color = NA),
            panel.background = element_rect(fill = "white", color = NA)
        )
}

show_plot <- function(p) {
    tryCatch(print(p), error = function(e) invisible(NULL))
}

save_pp <- function(p, base, pdf_w = 10, pdf_h = 8,
                    png_w = 2000, png_h = 1600, res = 200) {
    tryCatch({
        pdf(paste0(base, ".pdf"), width = pdf_w, height = pdf_h)
        print(p); dev.off()
    }, error = function(e) {
        tryCatch(dev.off(), error = function(e2) invisible(NULL))
        cat("  [save_pp] PDF error:", conditionMessage(e), "\n")
    })
    tryCatch({
        png(paste0(base, ".png"), width = png_w, height = png_h, res = res)
        print(p); dev.off()
    }, error = function(e) {
        tryCatch(dev.off(), error = function(e2) invisible(NULL))
        cat("  [save_pp] PNG error:", conditionMessage(e), "\n")
    })
}

save_grob <- function(g, base, pdf_w = 14, pdf_h = 7.5,
                      png_w = 2800, png_h = 1500, res = 200) {
    tryCatch({
        pdf(paste0(base, ".pdf"), width = pdf_w, height = pdf_h)
        grid.draw(g); dev.off()
    }, error = function(e) {
        tryCatch(dev.off(), error = function(e2) invisible(NULL))
        cat("  [save] PDF error:", conditionMessage(e), "\n")
    })
    tryCatch({
        png(paste0(base, ".png"), width = png_w, height = png_h, res = res)
        grid.draw(g); dev.off()
    }, error = function(e) {
        tryCatch(dev.off(), error = function(e2) invisible(NULL))
        cat("  [save] PNG error:", conditionMessage(e), "\n")
    })
}


# =============================================================================
# SECTION 0e: RDS size control
# =============================================================================
# Fitted models are saved with their environments pointed at baseenv(),
# and plot objects are dropped unless RDS_KEEP_PLOTS = TRUE.
RDS_KEEP_PLOTS <- FALSE

# Points the model environments of a gam, glm or lm at baseenv().
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

# Strips fitted models and, unless RDS_KEEP_PLOTS, replaces plots with NULL.
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

    .p10v_preamble <- character(0)
    assign("CEP_ARM", arm_tag, envir = globalenv())
    on.exit(suppressWarnings(rm("CEP_ARM", envir = globalenv())), add = TRUE)

    .p10v_preamble <- tryCatch(
        utils::capture.output({
            cat("\n", strrep("=", 76), "\n", sep = "")
            cat("  ARM: ", arm_tag, "   config: ", config_file, "\n", sep = "")
            cat(strrep("=", 76), "\n", sep = "")
            source(config_path)
        }),
        error = function(e) {
            cat("\n  ARM: ", arm_tag, "   config: ", config_file, "\n", sep = "")
            stop(e)
        })

    .p10v_needed <- c("ARM_TAG", "R_EXPORT_ROOT", "GOI", "DFG", "HKG_GENES",
                      "DFG_COMPOSITE_LABEL", "PRIMARY_STYPE", "NORM_TAG",
                      "NORM_TARGET_SUM", "NORM_CEILING", "GAM_BS", "GAM_METHOD",
                      "GAM_SELECT")
    .p10v_absent <- .p10v_needed[!vapply(.p10v_needed, exists, logical(1))]
    if (length(.p10v_absent) > 0)
        stop("The configuration loaded from\n    ", config_path,
             "\n  does not define: ", paste(.p10v_absent, collapse = ", "),
             call. = FALSE)
    if (!identical(ARM_TAG, arm_tag))
        stop("The configuration reports ARM_TAG '", ARM_TAG,
             "' but this arm is '", arm_tag, "'.", call. = FALSE)

    .p10v_preamble <- c(.p10v_preamble, utils::capture.output(
        cat("  Config OK: ", ARM_TAG, " | ", GOI, " + ", length(DFG), " DFG | ",
            NORM_TAG, "\n", sep = "")))

    set.seed(42)

    DATASETS        <- c("NeoC", "NonNeoC")
    GOI_LABEL_PLAIN <- GOI
    # Axis cap: the CP10K ceiling.
    X_MAX_CAP       <- NORM_CEILING
    # cHKG reference: housekeeping genes minus the GOI. If it equals the DFG
    # set (HKG8 arm), the cHKG curve is skipped.
    HKG_V1          <- setdiff(HKG_GENES[HKG_GENES != "18S rRNA"], GOI)
    HKG_OVERLAP_DFG <- intersect(HKG_V1, DFG)
    HKG_DEGENERATE  <- length(setdiff(HKG_V1, DFG)) == 0L
    if (length(HKG_OVERLAP_DFG) > 0)
        cat(sprintf(paste0("  cHKG reference: %d of %d gene(s) also in the DFG",
                           " set (%s)\n"),
                    length(HKG_OVERLAP_DFG), length(HKG_V1),
                    paste(HKG_OVERLAP_DFG, collapse = ", ")))
    if (HKG_DEGENERATE)
        cat("  cHKG reference equals the DFG set on this arm; cHKG curve skipped.\n")

    # == Input paths ==========================================================
    CMRF_DIR    <- cep_arm_dir_for("Part10e_CMRF", arm_tag)
    CMRF_PREFIX <- paste0(cep_out_prefix("Part10e", NORM_TAG,
                                         cep_file_tag_for(arm_tag)), "_CMRF")
    p10v_xlsx <- function(ds)
        file.path(CMRF_DIR, paste0(CMRF_PREFIX, "_", ds,
                                   "_SK_interaction_analysis.xlsx"))

    # == Output ===============================================================
    OUT_DIR <- out_dir
    dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
    FILE_PREFIX <- paste0(cep_out_prefix(PART_TAG, NORM_TAG,
                                         cep_file_tag_for(arm_tag)), "_V1")

    # == Logging ==============================================================
    LOG_MESSAGES_V1 <- TRUE
    LOG_PATH_V1 <- file.path(OUT_DIR, paste0(FILE_PREFIX, "_console_log_",
                                             format(Sys.time(), "%Y%m%d_%H%M%S"),
                                             ".txt"))
    MSG_PATH_V1 <- file.path(OUT_DIR, paste0(FILE_PREFIX,
                                             "_console_log_messages.tmp"))
    while (sink.number() > 0) sink()
    if (sink.number(type = "message") != 2L) sink(type = "message")
    # Set before on.exit() registers, so the handler can always read it.
    .p10v_completed <- FALSE
    .p10v_log_con <- file(LOG_PATH_V1, open = "wt")
    .p10v_msg_con <- if (isTRUE(LOG_MESSAGES_V1))
        file(MSG_PATH_V1, open = "wt") else NULL
    sink(.p10v_log_con, append = FALSE, split = TRUE)
    if (!is.null(.p10v_msg_con)) sink(.p10v_msg_con, type = "message")
    # The exit handler runs four independent stages.
    on.exit({
        tryCatch({
            if (!isTRUE(.p10v_completed))
                cat("\n", strrep("!", 72),
                    "\n  Arm did not complete; the error is printed to the console.\n",
                    strrep("!", 72), "\n", sep = "")
        }, error = function(e) invisible(NULL))
        tryCatch({
            if (sink.number(type = "message") != 2L) sink(type = "message")
            while (sink.number() > 0) sink()
        }, error = function(e) invisible(NULL))
        tryCatch(if (!is.null(.p10v_msg_con)) close(.p10v_msg_con),
                 error = function(e) invisible(NULL))
        tryCatch(close(.p10v_log_con), error = function(e) invisible(NULL))
        tryCatch({
            if (!is.null(.p10v_msg_con) && file.exists(MSG_PATH_V1)) {
                .msg_txt <- readLines(MSG_PATH_V1, warn = FALSE)
                if (length(.msg_txt) > 0)
                    cat(c("", strrep("=", 72),
                          "  MESSAGE STREAM (warnings, messages, errors)",
                          strrep("=", 72), .msg_txt, ""),
                        sep = "\n", file = LOG_PATH_V1, append = TRUE)
                unlink(MSG_PATH_V1)
            }
        }, error = function(e) invisible(NULL))
    }, add = TRUE)

    .p10v_warn_opt  <- getOption("warn")
    options(warn = 1)
    on.exit(options(warn = .p10v_warn_opt), add = TRUE)
    cat(sprintf("Part10e Validation1 - log started: %s\n\n", Sys.time()))
    if (length(.p10v_preamble) > 0)
        cat(paste(.p10v_preamble, collapse = "\n"), "\n\n", sep = "")

    # == IP read from the IPRS RDS ============================================
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
    IPRS_RDS_PATH <- if (length(.iprs_hit) > 0L) .iprs_hit[1] else .iprs_candidates[1]

    IPRS_IP_X <- NULL; IPRS_IP_Y <- NULL; ip_from <- "no IPRS RDS read"
    if (file.exists(IPRS_RDS_PATH)) {
        .iprs     <- readRDS(IPRS_RDS_PATH)
        IPRS_IP_X <- .iprs$inflection_points_DFG
        IPRS_IP_Y <- .iprs$inflection_points_DFG_y
        ip_from   <- if (!is.null(.iprs$ip_definition)) .iprs$ip_definition
                     else "unrecorded definition"
        cat("  [IPRS] ", basename(IPRS_RDS_PATH), " | definition: ", ip_from,
            "\n", sep = "")
        if (!is.null(IPRS_IP_X))
            for (.d in names(IPRS_IP_X))
                cat(sprintf("  [IPRS] %-10s IP_x = %.6f | IP_y = %s\n", .d,
                            IPRS_IP_X[.d],
                            if (!is.null(IPRS_IP_Y) && .d %in% names(IPRS_IP_Y))
                                sprintf("%.6f", IPRS_IP_Y[.d]) else "NA"))
    } else {
        cat("  [IPRS] RDS not found. Looked in:\n")
        for (.p in .iprs_candidates) cat("      ", .p, "\n", sep = "")
    }

    # == Donor lookup from the per-region gam_data CSVs =======================
    .p10v_donor_cache <- list()
    p10v_donor_map <- function(ds) {
        if (!is.null(.p10v_donor_cache[[ds]])) return(.p10v_donor_cache[[ds]])
        gam_dir <- if (identical(ds, "NeoC")) NEOC_DIR else NONNC_DIR
        regs    <- if (identical(ds, "NeoC")) NEOC_REGIONS else NONNC_REGIONS
        stype   <- gsub("-", "_", gsub("/", "_", PRIMARY_STYPE))
        bc <- character(0); don <- character(0)
        for (reg in regs) {
            gp <- file.path(gam_dir, "gam_data",
                            paste0(reg, "_", stype, "_selected.csv"))
            if (!file.exists(gp)) next
            df <- read.csv(gp, stringsAsFactors = FALSE, check.names = FALSE)
            if (!all(c("barcode", "donor") %in% colnames(df)))
                stop(sprintf("Region %s: gam_data lacks barcode or donor.", reg),
                     call. = FALSE)
            bc  <- c(bc,  df$barcode)
            don <- c(don, as.character(df$donor))
        }
        if (length(bc) == 0L)
            stop("No gam_data CSVs found for ", ds, call. = FALSE)
        # Duplicated barcodes must agree on the donor.
        dup <- duplicated(bc)
        if (any(dup)) {
            .first <- setNames(don[!dup], bc[!dup])
            .clash <- unique(bc[dup][don[dup] != .first[bc[dup]]])
            .clash <- .clash[!is.na(.clash)]
            if (length(.clash) > 0)
                stop(sprintf(paste0("%s: %d barcode(s) appear in more than one",
                                    " region with different donors, e.g. %s."),
                             ds, length(.clash),
                             paste(utils::head(.clash, 3), collapse = ", ")),
                     call. = FALSE)
            cat(sprintf(paste0("    %d duplicated barcode(s) across regions,",
                               " all agreeing on the donor; keeping one each\n"),
                        sum(dup)))
            bc <- bc[!dup]; don <- don[!dup]
        }
        m <- setNames(don, bc)
        .p10v_donor_cache[[ds]] <<- m
        m
    }

    DATASET_META <- list(
        NeoC    = list(xlsx = p10v_xlsx("NeoC"),    gam_dir = NEOC_DIR,
                       regions = NEOC_REGIONS),
        NonNeoC = list(xlsx = p10v_xlsx("NonNeoC"), gam_dir = NONNC_DIR,
                       regions = NONNC_REGIONS)
    )

    cat("\n", strrep("=", 72), "\n")
    cat("  Part10e Validation1 - CMRF final-module validation\n")
    cat("  Arm          :", ARM_TAG, "\n")
    cat("  GOI          :", GOI, "\n")
    cat("  Normalisation:", NORM_TAG, " log1p(UMIs/total x ",
        format(NORM_TARGET_SUM, scientific = FALSE), ")\n")
    cat("  GAM          : GOI + s(GOI) +",
        if (GAM_DONOR_RE_V1) "s(donor, bs = \"re\")" else "(no donor RE)",
        "|", if (W_DONOR_WEIGHTS_V1) "inverse-donor weights" else "uniform weights",
        "\n")
    cat("  Pre/Post-IP  :", IP_CLASS_AXIS_V1, "axis\n")
    cat("  CMRF input   :", CMRF_DIR, "\n")
    cat("  Output dir   :", OUT_DIR, "\n")
    cat(strrep("=", 72), "\n\n")

    # Marker file, removed when the arm completes.
    INCOMPLETE_MARKER_V1 <- file.path(OUT_DIR, "_RUN_INCOMPLETE.txt")
    writeLines(c("This arm did not finish.", "",
                 paste("Arm     :", arm_tag),
                 paste("Part    :", PART_DIR_TAG),
                 paste("Started :", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), "",
                 "Any files in this folder are from a partial run.",
                 "This file is removed automatically when the arm completes."),
               INCOMPLETE_MARKER_V1)


    # =============================================================================
    # SECTION A0: PRE-ANALYSIS GUARDRAILS
    # =============================================================================
    GR_LOG_V1  <- file.path(OUT_DIR, paste0(FILE_PREFIX, "_Guardrail_Log.txt"))
    .grv1_log  <- list()
    .grv1_pass <- 0L; .grv1_fail <- 0L; .grv1_warn <- 0L

    grv1 <- function(check_id, status, detail) {
        # A zero-length detail would vanish from console and ledger.
        if (length(detail) != 1L || is.na(detail))
            detail <- paste0("[detail unavailable: the message was built from a",
                             " missing value and came out zero-length]")
        # An unrecognised status is recorded as FAIL.
        if (!status %in% c("PASS", "FAIL", "WARN")) {
            detail <- paste0("[status '", as.character(status)[1],
                             "' is not PASS/FAIL/WARN] ", detail)
            status <- "FAIL"
        }
        tag <- switch(status, PASS = "[PASS]", FAIL = "[FAIL]", WARN = "[WARN]")
        cat(sprintf("  %s  %-16s  %s\n", tag, check_id, detail))
        .grv1_log[[length(.grv1_log) + 1L]] <<- list(
            id = check_id, status = status, detail = detail)
        if (status == "PASS") .grv1_pass <<- .grv1_pass + 1L
        if (status == "FAIL") .grv1_fail <<- .grv1_fail + 1L
        if (status == "WARN") .grv1_warn <<- .grv1_warn + 1L
    }

    grv1_write <- function() {
        con <- file(GR_LOG_V1, "wt")
        # on.exit() closes the connection on any failure.
        on.exit(try(close(con), silent = TRUE), add = TRUE)
        writeLines(c(strrep("=", 72),
                     paste0("Part 10e Validation1 - Guardrail Log [", NORM_TAG, "]"),
                     paste("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
                     paste("Arm :", ARM_TAG),
                     paste0("Normalisation: adata.X = log1p(UMIs/total x ",
                            format(NORM_TARGET_SUM, scientific = FALSE),
                            "), exported unchanged"),
                     strrep("=", 72), ""), con)
        for (item in .grv1_log)
            writeLines(sprintf("  [%s]  %-16s  %s", item$status, item$id,
                               item$detail), con)
        writeLines(c("", strrep("=", 72),
                     sprintf("  Total : %d", .grv1_pass + .grv1_fail + .grv1_warn),
                     sprintf("  PASS  : %d", .grv1_pass),
                     sprintf("  WARN  : %d", .grv1_warn),
                     sprintf("  FAIL  : %d", .grv1_fail),
                     if (.grv1_fail == 0 && .grv1_warn == 0)
                         "\n  Result: no failures, no warnings"
                     else if (.grv1_fail == 0)
                         "\n  Result: no failures; warnings present"
                     else "\n  Result: failures present",
                     strrep("=", 72)), con)
        cat("  Guardrail log written:", GR_LOG_V1, "\n")
    }

    cat("  SECTION A0 - PRE-ANALYSIS GUARDRAILS\n\n")

    # Positive tests on the CP10K normalisation.
    for (.ds in list(list(tag = "NeoC",    dir = NEOC_DIR),
                     list(tag = "NonNeoC", dir = NONNC_DIR))) {
        norm_mtx_path <- file.path(.ds$dir, "deg_support", "norm.mtx")
        raw_mtx_path  <- file.path(.ds$dir, "counts", "raw_selected", "matrix.mtx")
        if (file.exists(norm_mtx_path)) {
            X_n   <- as(readMM(norm_mtx_path), "CsparseMatrix")
            n_max <- max(X_n)
            grv1(paste0("GR-01-", .ds$tag),
                 if (n_max > 0 && n_max <= NORM_CEILING + 1e-4) "PASS" else "FAIL",
                 sprintf("norm.mtx max = %.4f  (ceiling log1p(%s) = %.4f)",
                         n_max, format(NORM_TARGET_SUM, scientific = FALSE),
                         NORM_CEILING))
            n_rows <- min(50L, nrow(X_n))
            rs     <- rowSums(expm1(X_n[seq_len(n_rows), , drop = FALSE]))
            rs_dev <- max(abs(rs - NORM_TARGET_SUM))
            grv1(paste0("GR-02-", .ds$tag),
                 if (rs_dev < NORM_TARGET_SUM * 1e-6) "PASS" else "FAIL",
                 sprintf("expm1(norm) row sums: mean = %.6f, max |dev| = %.2e, n = %d",
                         mean(rs), rs_dev, n_rows))
            rm(X_n); gc(verbose = FALSE)
        } else {
            grv1(paste0("GR-01-", .ds$tag), "FAIL",
                 paste("norm.mtx not found:", norm_mtx_path))
        }
        if (file.exists(raw_mtx_path)) {
            X_r    <- as(readMM(raw_mtx_path), "CsparseMatrix")
            r_samp <- as.numeric(X_r@x[seq_len(min(10000, length(X_r@x)))])
            is_int <- all(abs(r_samp - round(r_samp)) < 0.01)
            r_max  <- max(X_r)
            grv1(paste0("GR-04-", .ds$tag),
                 if (is_int && r_max > 20) "PASS" else "FAIL",
                 sprintf("raw max = %.0f, integer-like = %s", r_max,
                         if (is_int) "TRUE" else "FALSE"))
            rm(X_r); gc(verbose = FALSE)
        } else {
            grv1(paste0("GR-04-", .ds$tag), "FAIL",
                 paste("raw_selected matrix not found:", raw_mtx_path))
        }
    }

    for (ds in DATASETS) {
        .p <- p10v_xlsx(ds)
        grv1(paste0("GR-03-", ds),
             if (file.exists(.p)) "PASS" else "FAIL",
             sprintf("CMRF workbook: %s",
                     if (file.exists(.p))
                         file.path(basename(dirname(.p)), basename(.p))
                     else paste("not found:", .p)))
    }

    grv1("GR-05-IPRS",
         if (!is.null(IPRS_IP_X) && length(IPRS_IP_X) > 0L) "PASS" else "WARN",
         if (!is.null(IPRS_IP_X) && length(IPRS_IP_X) > 0L)
             sprintf("IP values (%s): %s", ip_from,
                     paste(sprintf("%s=%.4f", names(IPRS_IP_X),
                                   as.numeric(IPRS_IP_X)), collapse = ", "))
         else "IP values not available from the IPRS RDS")

    cat(sprintf("\n  Pre-analysis guardrails: %d PASS / %d WARN / %d FAIL\n\n",
                .grv1_pass, .grv1_warn, .grv1_fail))


    # =============================================================================
    # SECTION B: COLOUR PALETTES
    # =============================================================================
    four_cols <- c(
        "Pre-IP TREP"      = "#E69F00",
        "Pre-IP non-TREP"  = "#1B3A8C",
        "Post-IP TREP"     = "#009E73",
        "Post-IP non-TREP" = "#D94F4F"
    )
    trep_groups    <- c("Pre-IP TREP",     "Post-IP TREP")
    nontrep_groups <- c("Pre-IP non-TREP", "Post-IP non-TREP")

    gs_colors <- c("cDFG" = "#4B0082", "cHKG" = "#4F75DE")
    gs_fills  <- gs_colors

    cep_cols <- c("TREP" = "#4B0082", "non-TREP" = "#C0C0C0")


    # =============================================================================
    # SECTION C: MAIN FUNCTION (runs per dataset)
    # =============================================================================

    run_V1_plots <- function(dataset_tag) {

        cat("\n", strrep("=", 72), "\n")
        cat("  Part10e Validation1:", dataset_tag, " [", ARM_TAG, "]\n")
        cat(strrep("=", 72), "\n\n")

        meta     <- DATASET_META[[dataset_tag]]
        xlsx_path <- meta$xlsx
        gam_dir   <- meta$gam_dir

        # == C1. Load data from the CMRF workbook =============================
        cat("  Loading Excel data ...\n")
        all_cells   <- as.data.frame(read_excel(xlsx_path,
                                                sheet = "Step4_Refined_Cells"))
        # Columns this script needs, checked at the read.
        .v1_need <- c("barcode", "cell_class", "in_union_skcs",
                      "GOI_expr", "cDFG_expr", "spear_norm", "kend_contrib")
        .v1_miss <- setdiff(.v1_need, names(all_cells))
        if (length(.v1_miss) > 0)
            stop("Step4_Refined_Cells in\n  ", xlsx_path,
                 "\n  lacks: ", paste(.v1_miss, collapse = ", "), call. = FALSE)
        if (nrow(all_cells) == 0L)
            stop("Step4_Refined_Cells is empty for ", dataset_tag, ".",
                 call. = FALSE)
        .v1_bc_dup <- sum(duplicated(as.character(all_cells$barcode)))
        if (.v1_bc_dup > 0)
            stop(.v1_bc_dup, " duplicated barcode(s) in Step4_Refined_Cells for ",
                 dataset_tag, ".", call. = FALSE)
        if (anyNA(all_cells$cell_class))
            stop(sum(is.na(all_cells$cell_class)),
                 " cell(s) carry no cell_class in Step4_Refined_Cells for ",
                 dataset_tag, ".",
                 call. = FALSE)
        if (anyNA(all_cells$in_union_skcs))
            stop(sum(is.na(all_cells$in_union_skcs)),
                 " cell(s) carry no in_union_skcs for ", dataset_tag, ".",
                 call. = FALSE)

        # == C2. IP thresholds ================================================
        # IP_x and IP_y from the IPRS RDS.
        IP_x <- if (!is.null(IPRS_IP_X) && dataset_tag %in% names(IPRS_IP_X))
            as.numeric(IPRS_IP_X[dataset_tag]) else NA_real_
        IP_y      <- if (!is.null(IPRS_IP_Y) && dataset_tag %in% names(IPRS_IP_Y))
            as.numeric(IPRS_IP_Y[dataset_tag]) else NA_real_
        cat(sprintf("  IP_x = %s  |  IP_y = %s   (from the IPRS RDS)\n",
                    if (is.finite(IP_x)) sprintf("%.6f", IP_x) else "NA",
                    if (is.finite(IP_y))      sprintf("%.6f", IP_y)      else "NA"))

        # Pre/Post split axis; a missing IP_y on the cDFG axis is a hard stop.
        IP_AXIS_USED <- IP_CLASS_AXIS_V1
        if (identical(IP_CLASS_AXIS_V1, "cDFG") && !is.finite(IP_y))
            stop("IP_CLASS_AXIS_V1 is \"cDFG\" but the IPRS RDS carries no IP_y ",
                 "for ", dataset_tag, ".", call. = FALSE)
        IP_ACTIVE <- if (identical(IP_AXIS_USED, "cDFG")) IP_y else IP_x
        cat(sprintf("  Pre/Post-IP split on %s at %s\n", IP_AXIS_USED,
                    if (is.finite(IP_ACTIVE)) sprintf("%.6f", IP_ACTIVE) else "NA"))

        # == C3. Cell counts ==================================================
        n_total_M4top <- nrow(all_cells)
        n_goi_neg     <- sum(all_cells$cell_class == "GOI-negative")
        n_goi_pos     <- n_total_M4top - n_goi_neg
        n_union       <- sum(all_cells$in_union_skcs == TRUE)
        n_origref     <- sum(all_cells$cell_class == "Original + Refined")
        n_recovered   <- sum(all_cells$cell_class == "Recovered")
        n_orig_only   <- sum(all_cells$cell_class == "Original only")
        n_excluded    <- sum(all_cells$cell_class == "Excluded")
        n_non_module  <- n_goi_pos - n_union

        cat(sprintf("  Total=%d  GOI+=%d  Module=%d  (OrigRef=%d  Recovered=%d)  Gray=%d\n",
                    n_total_M4top, n_goi_pos, n_union, n_origref, n_recovered, n_non_module))

        # == C4. Filter =======================================================
        goi_pos_cells <- all_cells[all_cells$cell_class != "GOI-negative", ]
        module_cells  <- goi_pos_cells[goi_pos_cells$in_union_skcs == TRUE, ]

        # == C5. Fit cDFG GAM =================================================
        cat("\n  Fitting cDFG GAM on final module ...\n")
        # Donor matched by barcode; an unmatched barcode is a hard stop.
        donor_map_ds <- p10v_donor_map(dataset_tag)
        donor_module <- unname(donor_map_ds[module_cells$barcode])
        if (anyNA(donor_module))
            stop(sprintf(paste0("%s: %d of %d module barcode(s) have no donor in ",
                                "the gam_data CSVs."), dataset_tag,
                         sum(is.na(donor_module)), nrow(module_cells)),
                 call. = FALSE)
        cat(sprintf("  Donors on the module cells: %d\n",
                    n_donor_levels_v1(donor_module)))

        df_gam <- data.frame(GOI = module_cells$GOI_expr,
                             Expression = module_cells$cDFG_expr,
                             barcode = module_cells$barcode,
                             donor = as_donor_factor_v1(donor_module),
                             stringsAsFactors = FALSE)
        res_dfg <- .fit_full_gam(df_gam, label = paste0("cDFG_", dataset_tag),
                                 donor = df_gam$donor)
        if (is.null(res_dfg) || is.null(res_dfg$best_model)) {
            cat("  ERROR: cDFG GAM fitting failed\n"); return(NULL)
        }
        gam_model <- res_dfg$best_model
        # model_de (donor-aware) is reported; model_de_pop sizes the TREP set.
        model_de     <- de_explained_v1(gam_model, df_gam)
        model_de_pop <- pop_dev_explained_v1(gam_model, df_gam)
        cat(sprintf("  cDFG DE (donor-aware): %.2f%%\n", model_de * 100))
        cat(sprintf("  TREP sizing DE (flat null, population fit): %.2f%%\n",
                    model_de_pop * 100))

        # == C6. Fit cHKG GAM (from norm.mtx) =================================
        cat("\n  Loading cHKG data from norm.mtx ...\n")
        res_hkg      <- NULL
        df_hkg       <- NULL
        hkg_data_ok  <- FALSE

        mtx_path  <- file.path(gam_dir, "deg_support", "norm.mtx")
        bar_path  <- file.path(gam_dir, "deg_support", "norm_barcodes.tsv")
        feat_path <- file.path(gam_dir, "deg_support", "norm_features.tsv")

        if (all(file.exists(c(mtx_path, bar_path, feat_path)))) {
            tryCatch({
                X_hkg   <- readMM(mtx_path)
                hkg_bc  <- readLines(bar_path, warn = FALSE)
                hkg_ft  <- read.table(feat_path, sep = "\t", header = FALSE,
                                      stringsAsFactors = FALSE)[, 2]
                rownames(X_hkg) <- hkg_bc
                colnames(X_hkg) <- hkg_ft
                hkg_avail <- HKG_V1[HKG_V1 %in% colnames(X_hkg)]
                if (length(hkg_avail) > 0 && !HKG_DEGENERATE) {
                    bc_module <- module_cells$barcode
                    bc_in_mtx <- bc_module[bc_module %in% hkg_bc]
                    if (length(bc_in_mtx) >= MIN_CELLS) {
                        sub_mat  <- as.matrix(
                            X_hkg[bc_in_mtx, hkg_avail, drop = FALSE])
                        hkg_expr <- if (length(hkg_avail) == 1)
                            as.numeric(sub_mat[, 1]) else rowMeans(sub_mat)
                        goi_vals <- module_cells$GOI_expr[
                            match(bc_in_mtx, module_cells$barcode)]
                        keep <- !is.na(goi_vals) & goi_vals > 0
                        df_hkg <- data.frame(
                            GOI = goi_vals[keep],
                            Expression = hkg_expr[keep],
                            barcode = bc_in_mtx[keep],
                            donor = as_donor_factor_v1(
                                unname(donor_map_ds[bc_in_mtx[keep]])),
                            stringsAsFactors = FALSE)
                        hkg_data_ok <- nrow(df_hkg) >= MIN_CELLS
                        cat(sprintf("    cHKG: %d cells, %d genes (%s)\n",
                                    nrow(df_hkg), length(hkg_avail),
                                    paste(hkg_avail, collapse = ", ")))
                    }
                }
                rm(X_hkg); gc(verbose = FALSE)
            }, error = function(e)
                cat("    norm.mtx error:", conditionMessage(e), "\n"))
        } else {
            cat("    norm.mtx not found, cHKG will be skipped\n")
        }

        if (hkg_data_ok) {
            cat("  Fitting cHKG GAM ...\n")
            res_hkg <- .fit_full_gam(df_hkg, label = paste0("cHKG_", dataset_tag),
                                     donor = df_hkg$donor, verbose = FALSE)
            if (!is.null(res_hkg) && !is.null(res_hkg$best_model))
                cat(sprintf("  cHKG DE (donor-aware): %.2f%%\n",
                            de_explained_v1(res_hkg$best_model, df_hkg) * 100))
        }

        # == C7. EP / TREP ====================================================
        cat("\n  Computing EP and TREP ...\n")
        # EP: grand-mean null against the population prediction (s(donor)
        # excluded).
        null_model   <- mgcv::gam(Expression ~ 1, data = df_gam)
        null_fitted  <- as.numeric(fitted(null_model))
        model_fitted <- as.numeric(p10v_predict_pop(gam_model, df_gam$GOI))
        # Both vectors must have one value per row.
        if (length(null_fitted) != nrow(df_gam) ||
            length(model_fitted) != nrow(df_gam))
            stop(sprintf(paste0("%s: the EP vectors are the wrong length ",
                                "(null %d, model %d, rows %d)."),
                         dataset_tag, length(null_fitted),
                         length(model_fitted), nrow(df_gam)), call. = FALSE)
        null_sq      <- (df_gam$Expression - null_fitted)^2
        model_sq     <- (df_gam$Expression - model_fitted)^2
        # EP = 0 when the null residual is zero.
        ep_values    <- ifelse(null_sq > 0, 1 - (model_sq / null_sq), 0)
        sorted_idx   <- order(ep_values, decreasing = TRUE, na.last = TRUE)

        # TREP count = n x sizing fraction (clamped to [0, 1]); a non-finite DE
        # gives zero TREP cells.
        n_ep         <- nrow(df_gam)
        # Sized by model_de_pop, on the same prediction as the EP ranking.
        .trep_frac   <- trep_count_fraction_v1(model_de_pop)
        if (!is.finite(.trep_frac)) {
            cat("  WARN: TREP sizing DE not finite; target_cells set to 0.\n")
            .trep_frac <- 0
        }
        target_cells <- max(0L, min(
            as.integer(round(n_ep * .trep_frac)), n_ep))
        trep_barcodes    <- if (target_cells > 0L)
            df_gam$barcode[sorted_idx[seq_len(target_cells)]] else character(0)
        nontrep_barcodes <- if (target_cells < n_ep)
            df_gam$barcode[sorted_idx[seq.int(target_cells + 1L, n_ep)]]
        else character(0)
        n_trep    <- length(trep_barcodes)
        n_nontrep <- length(nontrep_barcodes)
        cat(sprintf("  TREP=%d  non-TREP=%d\n", n_trep, n_nontrep))

        # == C8. 4-group labels ===============================================
        # Pre/Post on cDFG (< IP_y) or on GOI (<= IP_x), per IP_AXIS_USED.
        ip_status_bc <- setNames(
            if (identical(IP_AXIS_USED, "cDFG"))
                ifelse(module_cells$cDFG_expr <  IP_ACTIVE, "Pre-IP", "Post-IP")
            else
                ifelse(module_cells$GOI_expr  <= IP_ACTIVE, "Pre-IP", "Post-IP"),
            module_cells$barcode)
        trep_flag_bc <- setNames(
            ifelse(module_cells$barcode %in% trep_barcodes, "TREP", "non-TREP"),
            module_cells$barcode)
        four_group_bc <- setNames(paste(ip_status_bc, trep_flag_bc),
                                  module_cells$barcode)
        group_tab <- table(four_group_bc)
        cat("  4-group counts:\n")
        for (g in names(four_cols)) cat(sprintf("    %-22s %d\n", g,
                                                as.integer(ifelse(is.na(group_tab[g]), 0, group_tab[g]))))


        # =============================================================================
        # FIGURE 1: CEP-IP plot and Pitcher plot
        # =============================================================================
        cat("\n  [Fig 1] Building CEP-IP + Pitcher ...\n")

        # == Points data ======================================================
        pts <- data.frame(GOI = df_gam$GOI, Expression = df_gam$Expression,
                          type = "data", se.fit = NA_real_, barcode = df_gam$barcode,
                          four_group = four_group_bc[df_gam$barcode],
                          stringsAsFactors = FALSE)
        draw_ord_map <- c("Pre-IP non-TREP"=1L, "Post-IP non-TREP"=2L,
                          "Pre-IP TREP"=3L, "Post-IP TREP"=4L)
        pts$draw_ord <- draw_ord_map[match(pts$four_group, names(draw_ord_map))]
        pts$draw_ord[is.na(pts$draw_ord)] <- 0L
        pts <- pts[order(pts$draw_ord), ]

        # GAM curve
        x_seq <- seq(min(df_gam$GOI), max(df_gam$GOI), length.out = 1000)
        set.seed(123)
        # Population curve, s(donor) excluded.
        pred <- p10v_predict_pop(gam_model, x_seq, se.fit = TRUE)
        fit_df <- data.frame(GOI = x_seq, Expression = as.numeric(pred$fit),
                             type = "fit", se.fit = as.numeric(pred$se.fit),
                             barcode = NA_character_, four_group = NA_character_,
                             draw_ord = NA_integer_, stringsAsFactors = FALSE)

        # Gray cells
        non_module <- goi_pos_cells[goi_pos_cells$in_union_skcs == FALSE, ]
        gray_pts <- data.frame(GOI = non_module$GOI_expr,
                               Expression = non_module$cDFG_expr,
                               stringsAsFactors = FALSE)

        ABS_PAD <- 0.08
        all_x <- c(pts$GOI, gray_pts$GOI)
        all_y <- c(pts$Expression, gray_pts$Expression)
        x_lim <- c(min(all_x) - ABS_PAD, min(max(all_x) + ABS_PAD, X_MAX_CAP))
        y_lim <- c(max(0, min(all_y) - ABS_PAD), max(all_y) + ABS_PAD)

        all_cols_cep <- c(four_cols, "Non-module" = "#D0D0D0")
        cep_labs <- sapply(names(all_cols_cep), function(g) {
            if (g == "Non-module") sprintf("Non-module (n=%d)", nrow(non_module))
            else { n_g <- as.integer(group_tab[g]); if (is.na(n_g)) n_g <- 0L
            sprintf("%s (n=%d)", g, n_g) }
        })
        de_label <- sprintf("cDFG DE=%.1f%%", model_de * 100)

        # Build subtitle with all group counts
        n_preip_trep    <- as.integer(ifelse(is.na(group_tab["Pre-IP TREP"]), 0, group_tab["Pre-IP TREP"]))
        n_preip_nontrep <- as.integer(ifelse(is.na(group_tab["Pre-IP non-TREP"]), 0, group_tab["Pre-IP non-TREP"]))
        n_postip_trep   <- as.integer(ifelse(is.na(group_tab["Post-IP TREP"]), 0, group_tab["Post-IP TREP"]))
        n_postip_nontrep<- as.integer(ifelse(is.na(group_tab["Post-IP non-TREP"]), 0, group_tab["Post-IP non-TREP"]))
        n_gray          <- nrow(non_module)

        cep_subtitle <- sprintf(
            paste0("%s  |  IP=%.3f\n",
                   "Pre-IP: TREP=%d, non-TREP=%d  |  ",
                   "Post-IP: TREP=%d, non-TREP=%d  |  Gray=%d"),
            de_label, IP_x,
            n_preip_trep, n_preip_nontrep,
            n_postip_trep, n_postip_nontrep, n_gray)

        set.seed(123)
        p_cep <- ggplot() +
            geom_point(data = gray_pts, aes(x = GOI, y = Expression),
                       colour = "#D0D0D0", alpha = 0.25, size = 0.6) +
            geom_point(data = pts[pts$type == "data", ],
                       aes(x = GOI, y = Expression, colour = four_group),
                       position = position_jitter(width=0.01, height=0, seed=123),
                       alpha = 0.50, size = 0.9) +
            geom_line(data = fit_df, aes(x = GOI, y = Expression),
                      colour = "#4B0082", linewidth = 1.3) +
            geom_ribbon(data = fit_df,
                        aes(x = GOI, ymin = Expression - 1.96*se.fit,
                            ymax = Expression + 1.96*se.fit),
                        fill = "#4B0082", alpha = 0.12, colour = NA) +
            # Pre/Post boundary at IP_y (labels assigned on the cDFG axis).
            ip_hline_v1(IP_y, x_lim, size = 2.8) +
            scale_colour_manual(values = all_cols_cep, labels = cep_labs,
                                name = "Cell group",
                                breaks = c("Pre-IP TREP", "Pre-IP non-TREP",
                                           "Post-IP TREP", "Post-IP non-TREP",
                                           "Non-module")) +
            scale_x_continuous(labels = number_format(accuracy=0.1), limits = x_lim) +
            scale_y_continuous(labels = number_format(accuracy=0.1), limits = y_lim) +
            labs(title = "CEP-IP plot (cDFG vs GOI)",
                 subtitle = cep_subtitle,
                 x = bquote(italic(.(GOI))),
                 y = DFG_COMPOSITE_LABEL) +
            theme_p10b(base_size = 10) +
            theme(legend.position = "bottom",
                  legend.text = element_text(size = 7),
                  legend.title = element_text(size = 8, face = "bold"),
                  legend.key.size = unit(0.4, "cm"),
                  plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
                  plot.subtitle = element_text(size = 7, colour = "gray30", hjust = 0.5),
                  aspect.ratio = 1)

        # == Pitcher plot =====================================================
        pitcher_df <- data.frame(S = goi_pos_cells$spear_norm,
                                 K = goi_pos_cells$kend_contrib,
                                 barcode = goi_pos_cells$barcode,
                                 stringsAsFactors = FALSE)
        pitcher_df$group_4c <- "Non-module"
        in_mod <- pitcher_df$barcode %in% module_cells$barcode
        pitcher_df$group_4c[in_mod] <- four_group_bc[pitcher_df$barcode[in_mod]]

        draw_ord_p <- c("Non-module"=0L, "Pre-IP non-TREP"=1L,
                        "Post-IP non-TREP"=2L, "Pre-IP TREP"=3L, "Post-IP TREP"=4L)
        pitcher_df$draw_ord <- draw_ord_p[match(pitcher_df$group_4c, names(draw_ord_p))]
        pitcher_df$draw_ord[is.na(pitcher_df$draw_ord)] <- 0L
        pitcher_df <- pitcher_df[order(pitcher_df$draw_ord), ]

        # Pitcher subtitle with all group counts
        pitch_subtitle <- sprintf(
            paste0("Pre-IP: TREP=%d, non-TREP=%d  |  ",
                   "Post-IP: TREP=%d, non-TREP=%d  |  Gray=%d"),
            n_preip_trep, n_preip_nontrep,
            n_postip_trep, n_postip_nontrep, n_gray)

        p_pitch <- ggplot(pitcher_df, aes(x = S, y = K)) +
            geom_point(data = subset(pitcher_df, group_4c == "Non-module"),
                       colour = "#D0D0D0", alpha = 0.25, size = 0.6) +
            geom_point(data = subset(pitcher_df,
                                     group_4c %in% c(trep_groups, nontrep_groups)),
                       aes(colour = group_4c), alpha = 0.50, size = 0.85) +
            geom_vline(xintercept = 0.05, colour = "#009E73",
                       linetype = "solid", linewidth = 0.6) +
            geom_hline(yintercept = 0.5, colour = "#009E73",
                       linetype = "solid", linewidth = 0.6) +
            scale_colour_manual(values = c(four_cols, "Non-module" = "#D0D0D0"),
                                name = "Cell group") +
            labs(title = expression("Pitcher plot (" * K[i] * " vs " * S[i]^norm * ")"),
                 subtitle = pitch_subtitle,
                 x = expression(S[i]^norm), y = expression(K[i])) +
            theme_p10b(base_size = 10) +
            theme(legend.position = "none",
                  plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
                  plot.subtitle = element_text(size = 7, colour = "gray30", hjust = 0.5),
                  aspect.ratio = 1)

        # == Assemble Figure 1 ================================================
        title_text <- sprintf(
            "%s (%s\u207A; n=%d) \u2014 CMRF Final Module (OriRef+Recovered=%d)",
            dataset_tag, GOI_LABEL_PLAIN, n_goi_pos, n_union)
        header <- arrangeGrob(
            textGrob(title_text, gp = gpar(fontface = "bold", fontsize = 14)),
            linesGrob(x = unit(c(0.05, 0.95), "npc"), y = unit(c(0.3, 0.3), "npc"),
                      gp = gpar(col = "black", lwd = 1.5)),
            ncol = 1, heights = unit(c(0.7, 0.3), "cm"))

        fig1 <- arrangeGrob(header,
                            arrangeGrob(p_cep, p_pitch, ncol = 2,
                                        widths = c(1, 1)),
                            ncol = 1, heights = unit(c(1.2, 10), "cm"))
        grid.newpage(); grid.draw(fig1)
        save_grob(fig1, file.path(OUT_DIR, paste0(
            FILE_PREFIX, "_", dataset_tag, "_CEP-IP_Pitcher_FinalModule")),
            pdf_w = 14, pdf_h = 8.5, png_w = 2800, png_h = 1700)
        cat("  Fig 1 saved\n")


        # =============================================================================
        # FIGURE 2: GAM curves (cDFG + cHKG) + CEP-IP scatter (TREP/non-TREP)
        # =============================================================================
        cat("\n  [Fig 2] Building GAM curves + CEP-IP scatter ...\n")

        # == Panel A: dual GAM curve (cDFG and cHKG) ==========================
        plot_data_gs <- data.frame()

        for (item in list(
            list(res = res_dfg, label = "cDFG"),
            list(res = res_hkg, label = "cHKG"))) {
            if (is.null(item$res) || is.null(item$res$best_model)) next
            m   <- item$res$best_model
            dat <- item$res$data
            pts_gs <- data.frame(GOI = dat$GOI, Expression = dat$Expression,
                                 Gene_Set = item$label, se.fit = NA_real_,
                                 type = "data", stringsAsFactors = FALSE)
            x_seq_gs <- seq(min(dat$GOI), max(dat$GOI), length.out = 1000)
            fit_gs <- tryCatch({
                set.seed(123)
                # Population curve, s(donor) excluded.
                pred_gs <- p10v_predict_pop(m, x_seq_gs, se.fit = TRUE)
                data.frame(GOI = x_seq_gs, Expression = as.numeric(pred_gs$fit),
                           Gene_Set = item$label,
                           se.fit = as.numeric(pred_gs$se.fit),
                           type = "fit", stringsAsFactors = FALSE)
            }, error = function(e) NULL)
            plot_data_gs <- rbind(plot_data_gs, pts_gs)
            if (!is.null(fit_gs)) plot_data_gs <- rbind(plot_data_gs, fit_gs)
        }

        de_dfg_lab <- sprintf("cDFG DE=%.1f%%", model_de * 100)
        de_hkg_lab <- if (!is.null(res_hkg) && !is.null(res_hkg$best_model))
            sprintf("cHKG DE=%.1f%%",
                    de_explained_v1(res_hkg$best_model, df_hkg) * 100)
        else if (HKG_DEGENERATE)
            sprintf("cHKG reference not shown: identical to the %s set on this arm",
                    DFG_COMPOSITE_LABEL)
        else NA_character_
        de_labs_all <- c(de_dfg_lab, de_hkg_lab)[!is.na(c(de_dfg_lab, de_hkg_lab))]

        n_pts_gs <- max(
            sum(plot_data_gs$type == "data" & plot_data_gs$Gene_Set == "cDFG"),
            sum(plot_data_gs$type == "data" & plot_data_gs$Gene_Set == "cHKG"))

        x_data_gs <- plot_data_gs$GOI[plot_data_gs$type == "data"]
        y_data_gs <- plot_data_gs$Expression[plot_data_gs$type == "data"]
        x_lim_gs <- c(min(x_data_gs) - ABS_PAD, min(max(x_data_gs) + ABS_PAD, X_MAX_CAP))
        y_lim_gs <- c(max(0, min(y_data_gs) - ABS_PAD), max(y_data_gs) + ABS_PAD)

        p_gam_dual <- ggplot(plot_data_gs,
                             aes(x = GOI, y = Expression, colour = Gene_Set)) +
            geom_point(data = subset(plot_data_gs, type == "data"),
                       position = position_jitter(width = 0.01, height = 0,
                                                  seed = 123),
                       alpha = 0.25, size = 0.8) +
            geom_line(data = subset(plot_data_gs, type == "fit"),
                      linewidth = 1.2) +
            geom_ribbon(data = subset(plot_data_gs, type == "fit"),
                        aes(ymin = Expression - 1.96 * se.fit,
                            ymax = Expression + 1.96 * se.fit,
                            fill = Gene_Set),
                        alpha = 0.20, colour = NA) +
            scale_colour_manual(values = gs_colors, name = "Gene Set") +
            scale_fill_manual(values = gs_fills, name = "Gene Set") +
            scale_x_continuous(labels = number_format(accuracy = 0.1),
                               limits = x_lim_gs) +
            scale_y_continuous(labels = number_format(accuracy = 0.1),
                               limits = y_lim_gs) +
            labs(title = paste0("CMRF Final Module (n=", n_pts_gs,
                                ", \u03b3=", GAM_GAMMA_V1, ")"),
                 subtitle = paste(de_labs_all, collapse = "  |  "),
                 x = paste0(GOI_LABEL_PLAIN, " Expression (", NORM_TAG, " log1p)"),
                 y = paste0("Gene Set Expression (", NORM_TAG, " log1p)")) +
            theme_p10b(base_size = 9) +
            theme(legend.position = "bottom",
                  legend.text = element_text(size = 8),
                  legend.title = element_text(size = 9, face = "bold"),
                  legend.key.size = unit(0.4, "cm"),
                  plot.title = element_text(face = "bold", size = 10),
                  plot.subtitle = element_text(size = 7),
                  aspect.ratio = 1)
        cat("  GAM dual curve built\n")

        # == Panel B: CEP-IP scatter (TREP and non-TREP) ======================
        pd_cep <- data.frame(
            GOI = df_gam$GOI, Expression = df_gam$Expression,
            Group = ifelse(df_gam$barcode %in% trep_barcodes, "TREP", "non-TREP"),
            stringsAsFactors = FALSE)
        pd_cep$draw_order <- ifelse(pd_cep$Group == "TREP", 2L, 1L)
        pd_cep <- pd_cep[order(pd_cep$draw_order), ]

        set.seed(123)
        # Population curve, s(donor) excluded.
        pred_cep  <- p10v_predict_pop(gam_model, x_seq, se.fit = TRUE)
        pred_cep_df <- data.frame(GOI = x_seq,
                                  fit = as.numeric(pred_cep$fit),
                                  se.fit = as.numeric(pred_cep$se.fit))

        x_range_cep <- range(df_gam$GOI, na.rm = TRUE)
        y_range_cep <- range(df_gam$Expression, na.rm = TRUE)
        # The subtitle carries the TREP sizing DE.
        dev_pct <- round(model_de_pop * 100, 1)

        p_cep_scatter <- ggplot() +
            geom_point(data = pd_cep,
                       aes(x = GOI, y = Expression, color = Group),
                       size = 1.8, alpha = 0.45) +
            geom_line(data = pred_cep_df, aes(x = GOI, y = fit),
                      color = "#FFCC99", linewidth = 1.2) +
            geom_ribbon(data = pred_cep_df,
                        aes(x = GOI,
                            ymin = fit - 1.96 * se.fit,
                            ymax = fit + 1.96 * se.fit),
                        fill = "#FFCC99", alpha = 0.2) +
            # Pre/Post boundary at IP_y.
            ip_hline_v1(IP_y, x_range_cep, size = 2.5) +
            scale_color_manual(values = cep_cols, name = "CEP-IP") +
            theme_p10b(base_size = 9) +
            theme(legend.position = "bottom",
                  legend.text = element_text(size = 8),
                  legend.title = element_text(size = 9, face = "bold"),
                  legend.key.size = unit(0.4, "cm"),
                  plot.title = element_text(face = "bold", size = 10),
                  plot.subtitle = element_text(size = 7),
                  aspect.ratio = 1) +
            labs(title = paste0(GOI_LABEL_PLAIN, " vs cDFG \u2014 CMRF Final Module"),
                 subtitle = sprintf(
                     paste0("TREP sizing DE (flat null, population fit)=%.1f%%",
                            "  |  n=%d  (TREP=%d, non-TREP=%d)"),
                     dev_pct, nrow(df_gam), n_trep, n_nontrep),
                 x = paste0(GOI_LABEL_PLAIN, " Expression (", NORM_TAG, " log1p)"),
                 y = paste0(DFG_COMPOSITE_LABEL, " Expression (", NORM_TAG, " log1p)")) +
            scale_x_continuous(labels = number_format(accuracy = 0.1),
                               limits = c(x_range_cep[1] - ABS_PAD,
                                          min(x_range_cep[2] + ABS_PAD, X_MAX_CAP))) +
            scale_y_continuous(labels = number_format(accuracy = 0.1),
                               limits = c(max(0, y_range_cep[1] - ABS_PAD),
                                          y_range_cep[2] + ABS_PAD))
        cat("  CEP-IP scatter built\n")

        # == Assemble Figure 2 ================================================
        title2_text <- sprintf(
            "%s \u2014 CMRF Final Module: GAM Curves & CEP-IP Scatter (n=%d)",
            dataset_tag, n_union)
        header2 <- arrangeGrob(
            textGrob(title2_text, gp = gpar(fontface = "bold", fontsize = 13)),
            linesGrob(x = unit(c(0.05, 0.95), "npc"), y = unit(c(0.3, 0.3), "npc"),
                      gp = gpar(col = "black", lwd = 1.5)),
            ncol = 1, heights = unit(c(0.7, 0.3), "cm"))

        fig2 <- arrangeGrob(header2,
                            arrangeGrob(p_gam_dual, p_cep_scatter,
                                        ncol = 2, widths = c(1, 1)),
                            ncol = 1, heights = unit(c(1.2, 10), "cm"))
        grid.newpage(); grid.draw(fig2)
        save_grob(fig2, file.path(OUT_DIR, paste0(
            FILE_PREFIX, "_", dataset_tag, "_GAM_cHKG_CEP-IP_scatter")),
            pdf_w = 14, pdf_h = 8.5, png_w = 2800, png_h = 1700)
        cat("  Fig 2 saved\n\n")


        # =============================================================================
        # FIGURE 3: Copula space, cDFG vs GOI (u, v) and Ki vs Si (u_SK, v_SK)
        # =============================================================================
        cat("  [Fig 3] Building copula space plots ...\n")

        set.seed(42)

        # == Copula transform: rank/(n+1) over all GOI-positive cells =========
        n_cop <- nrow(goi_pos_cells)

        # Panel A: cDFG vs GOI in copula space (u = F_GOI(x), v = F_cDFG(y))
        set.seed(42)
        x_all <- goi_pos_cells$GOI_expr
        y_all <- goi_pos_cells$cDFG_expr
        u_cop <- rank(x_all, ties.method = "average") / (n_cop + 1)
        v_cop <- rank(y_all, ties.method = "average") / (n_cop + 1)

        # IP_u and IP_v: IP_x and IP_y projected onto the marginals, via /(n+1).
        IP_u <- if (is.finite(IP_x))
            sum(x_all <= IP_x) / (n_cop + 1) else NA_real_
        IP_v <- if (is.finite(IP_y))
            sum(y_all <= IP_y) / (n_cop + 1) else NA_real_
        cat(sprintf("  Copula (u,v): n=%d  IP_u=%s  IP_v=%s\n", n_cop,
                    if (is.finite(IP_u)) sprintf("%.4f", IP_u) else "NA",
                    if (is.finite(IP_v)) sprintf("%.4f", IP_v) else "NA"))

        # Panel B: Ki vs Si in copula space (u_SK = F(Si), v_SK = F(Ki))
        set.seed(42)
        Si_all <- goi_pos_cells$spear_norm
        Ki_all <- goi_pos_cells$kend_contrib
        u_sk <- rank(Si_all, ties.method = "average") / (n_cop + 1)
        v_sk <- rank(Ki_all, ties.method = "average") / (n_cop + 1)

        # SKCS cuts (spear_norm > 0.05, kend_contrib > 0.5) in copula coordinates.
        SI_NORM_CUT_V1 <- 0.05
        KI_CUT_V1      <- 0.5
        u_Sinorm05_cop <- mean(Si_all <= SI_NORM_CUT_V1)
        v_Ki05_cop     <- mean(Ki_all <= KI_CUT_V1)
        u_Si0_cop      <- u_Sinorm05_cop   # alias kept for downstream readers
        cat(sprintf(paste0("  SKCS cuts: Si_norm > %s, Ki > %s  ->  copula",
                           " (u_SK, v_SK) = (%.4f, %.4f)\n"),
                    format(SI_NORM_CUT_V1), format(KI_CUT_V1),
                    u_Sinorm05_cop, v_Ki05_cop))

        # == Build copula data frame for both panels ==========================
        cop_df <- data.frame(
            u = u_cop, v = v_cop,
            u_sk = u_sk, v_sk = v_sk,
            barcode = goi_pos_cells$barcode,
            stringsAsFactors = FALSE
        )
        cop_df$group <- "Non-module"
        in_mod_cop <- cop_df$barcode %in% module_cells$barcode
        cop_df$group[in_mod_cop] <- four_group_bc[cop_df$barcode[in_mod_cop]]

        draw_ord_cop <- c("Non-module"=0L, "Pre-IP non-TREP"=1L,
                          "Post-IP non-TREP"=2L, "Pre-IP TREP"=3L,
                          "Post-IP TREP"=4L)
        cop_df$draw_ord <- draw_ord_cop[match(cop_df$group, names(draw_ord_cop))]
        cop_df$draw_ord[is.na(cop_df$draw_ord)] <- 0L
        cop_df <- cop_df[order(cop_df$draw_ord), ]

        cop_cols <- c(four_cols, "Non-module" = "gray80")

        # == Panel A: cDFG vs MT-CO3 copula space =============================
        set.seed(123)
        p_cop_uv <- ggplot(cop_df, aes(x = u, y = v)) +
            geom_point(data = subset(cop_df, group == "Non-module"),
                       colour = "gray80", alpha = 0.30, size = 0.6) +
            geom_point(data = subset(cop_df, group != "Non-module"),
                       aes(colour = group), alpha = 0.50, size = 0.9) +
            # Pre/Post boundary at IP_v = Fc(IP_y).
            ip_hline_cop_v1(IP_v) +
            scale_colour_manual(values = cop_cols, name = "Cell group") +
            scale_x_continuous(breaks = seq(0, 1, 0.25)) +
            scale_y_continuous(breaks = seq(0, 1, 0.25)) +
            coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
            labs(title = "Copula space",
                 subtitle = sprintf("u=F(GOI)  v=F(cDFG)  n=%d", n_cop),
                 x = expression(""*"" %->% "" * u),
                 y = expression(""*"" %->% "" * v)) +
            theme_p10b(base_size = 10) +
            theme(legend.position = "bottom",
                  legend.text = element_text(size = 7),
                  legend.title = element_text(size = 8, face = "bold"),
                  legend.key.size = unit(0.35, "cm"),
                  plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
                  plot.subtitle = element_text(size = 7, colour = "gray30",
                                               hjust = 0.5),
                  aspect.ratio = 1)

        # == Panel B: Ki vs Si copula space ===================================
        set.seed(123)
        p_cop_sk <- ggplot(cop_df, aes(x = u_sk, y = v_sk)) +
            geom_point(data = subset(cop_df, group == "Non-module"),
                       colour = "gray80", alpha = 0.30, size = 0.6) +
            geom_point(data = subset(cop_df, group != "Non-module"),
                       aes(colour = group), alpha = 0.50, size = 0.9) +
            geom_vline(xintercept = u_Sinorm05_cop, colour = "#009E73",
                       linetype = "solid", linewidth = 0.6) +
            geom_hline(yintercept = v_Ki05_cop, colour = "#009E73",
                       linetype = "solid", linewidth = 0.6) +
            scale_colour_manual(values = cop_cols, name = "Cell group") +
            scale_x_continuous(breaks = seq(0, 1, 0.25)) +
            scale_y_continuous(breaks = seq(0, 1, 0.25)) +
            coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
            labs(title = "Copula space",
                 subtitle = sprintf(
                     "u[SK]=F(Si_norm)  v[SK]=F(Ki)  n=%d", n_cop),
                 x = expression(""*"" %->% "" * u[SK]),
                 y = expression(""*"" %->% "" * v[SK])) +
            theme_p10b(base_size = 10) +
            theme(legend.position = "none",
                  plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
                  plot.subtitle = element_text(size = 7, colour = "gray30",
                                               hjust = 0.5),
                  aspect.ratio = 1)

        cat("  Copula panels built\n")

        # == Assemble Figure 3 ================================================
        title3_text <- sprintf(
            "%s (%s\u207A, n=%d) \u2014 CMRF Final Module: Copula Space",
            dataset_tag, GOI_LABEL_PLAIN, n_goi_pos)
        header3 <- arrangeGrob(
            textGrob(title3_text, gp = gpar(fontface = "bold", fontsize = 13)),
            linesGrob(x = unit(c(0.05, 0.95), "npc"), y = unit(c(0.3, 0.3), "npc"),
                      gp = gpar(col = "black", lwd = 1.5)),
            ncol = 1, heights = unit(c(0.7, 0.3), "cm"))

        fig3 <- arrangeGrob(header3,
                            arrangeGrob(p_cop_uv, p_cop_sk,
                                        ncol = 2, widths = c(1, 1)),
                            ncol = 1, heights = unit(c(1.2, 10), "cm"))
        grid.newpage(); grid.draw(fig3)
        save_grob(fig3, file.path(OUT_DIR, paste0(
            FILE_PREFIX, "_", dataset_tag, "_Copula_space")),
            pdf_w = 14, pdf_h = 8.5, png_w = 2800, png_h = 1700)
        cat("  Fig 3 saved\n\n")


        # =============================================================================
        # FIGURES 4 & 5: ORIGINAL LABELS (no GAM refit)
        # OrigRef: Part 10b 4-group; Recovered: diamond; Original only: inverted
        # triangle.
        # =============================================================================
        cat("  [Fig 4-5] Building original-label plots (no refit) ...\n")

        set.seed(42)

        # == Colour and shape constants =======================================
        RECOV_COL   <- "#7C3AED"   # violet
        ORIGONLY_COL <- "#8B6914"  # dark tan

        # == Original group labels from the workbook ==========================
        # OrigRef cells: use group_refined (4-group from Part10b)
        # Recovered: "Recovered Pre-IP" / "Recovered Post-IP" (shown as diamonds)
        # Original-only: have Part10b group but shown as tan triangles
        # Excluded: gray background

        origref_cells <- goi_pos_cells[goi_pos_cells$cell_class == "Original + Refined", ]
        recov_cells   <- goi_pos_cells[goi_pos_cells$cell_class == "Recovered", ]
        oo_cells      <- goi_pos_cells[goi_pos_cells$cell_class == "Original only", ]
        excl_cells    <- goi_pos_cells[goi_pos_cells$cell_class == "Excluded", ]

        cat(sprintf("    OrigRef=%d  Recovered=%d  Orig-only=%d  Excluded=%d\n",
                    nrow(origref_cells), nrow(recov_cells),
                    nrow(oo_cells), nrow(excl_cells)))

        # OrigRef original 4-group
        # Concordance uses group_orig (Part 10b's assignment); any difference
        # from group_refined is reported.
        orig_grp_bc <- setNames(as.character(origref_cells$group_orig),
                                origref_cells$barcode)
        .gr_diff <- sum(as.character(origref_cells$group_orig) !=
                            as.character(origref_cells$group_refined),
                        na.rm = TRUE)
        if (.gr_diff > 0)
            cat(sprintf(paste0("    NOTE: group_orig and group_refined differ",
                               " on %d of %d Original+Refined cell(s); the",
                               " concordance uses group_orig.\n"),
                        .gr_diff, nrow(origref_cells)))

        # Draw-order for OrigRef
        origref_plot <- data.frame(
            GOI = origref_cells$GOI_expr, Expression = origref_cells$cDFG_expr,
            S = origref_cells$spear_norm, K = origref_cells$kend_contrib,
            barcode = origref_cells$barcode,
            group_orig = orig_grp_bc[origref_cells$barcode],
            stringsAsFactors = FALSE)
        ord_map <- c("Pre-IP non-TREP"=1L, "Post-IP non-TREP"=2L,
                     "Pre-IP TREP"=3L, "Post-IP TREP"=4L)
        origref_plot$draw_ord <- ord_map[match(origref_plot$group_orig,
                                               names(ord_map))]
        origref_plot$draw_ord[is.na(origref_plot$draw_ord)] <- 0L
        origref_plot <- origref_plot[order(origref_plot$draw_ord), ]
        # A complete factor, so the legend always has four keys.
        origref_plot$group_orig <- factor(origref_plot$group_orig,
                                          levels = names(four_cols))

        # Subtitle with original group counts
        or_tab <- table(origref_plot$group_orig)
        or_preip_t  <- as.integer(ifelse(is.na(or_tab["Pre-IP TREP"]), 0, or_tab["Pre-IP TREP"]))
        or_preip_nt <- as.integer(ifelse(is.na(or_tab["Pre-IP non-TREP"]), 0, or_tab["Pre-IP non-TREP"]))
        or_postip_t <- as.integer(ifelse(is.na(or_tab["Post-IP TREP"]), 0, or_tab["Post-IP TREP"]))
        or_postip_nt<- as.integer(ifelse(is.na(or_tab["Post-IP non-TREP"]), 0, or_tab["Post-IP non-TREP"]))

        orig_sub <- sprintf(
            paste0("Original Part10b labels (no refit)  |  IP=%.3f\n",
                   "Pre-IP: TREP=%d, non-TREP=%d  |  Post-IP: TREP=%d, non-TREP=%d\n",
                   "Recovered=%d (\u25C6)  |  Orig-only=%d (\u25BD)  |  Excluded=%d (gray)"),
            IP_x, or_preip_t, or_preip_nt, or_postip_t, or_postip_nt,
            nrow(recov_cells), nrow(oo_cells), nrow(excl_cells))

        # == Figure 4A: CEP-IP (cDFG vs GOI), original labels =================
        set.seed(123)
        p4_cep <- ggplot() +
            # Excluded (gray background)
            geom_point(data = excl_cells,
                       aes(x = GOI_expr, y = cDFG_expr),
                       colour = "#D0D0D0", alpha = 0.25, size = 0.6) +
            # Original-only (tan inverted triangle)
            geom_point(data = oo_cells,
                       aes(x = GOI_expr, y = cDFG_expr),
                       colour = "white", fill = ORIGONLY_COL,
                       shape = 25, size = 2.2, alpha = 0.88, stroke = 0.15) +
            # OrigRef cells (4-group circles, original labels)
            geom_point(data = origref_plot,
                       aes(x = GOI, y = Expression, colour = group_orig),
                       size = 1.35, alpha = 0.58, shape = 16) +
            # Recovered (purple diamond)
            geom_point(data = recov_cells,
                       aes(x = GOI_expr, y = cDFG_expr),
                       colour = "white", fill = RECOV_COL,
                       shape = 23, size = 2.2, alpha = 0.93, stroke = 0.15) +
            # IP threshold at IP_y (cDFG axis)
            ip_hline_v1(IP_y, x_lim, size = 2.8) +
            scale_colour_manual(
                values = four_cols, name = "OrigRef group", drop = FALSE,
                guide = guide_legend(override.aes = list(
                    size = rep(2.5, 4), shape = rep(16L, 4), alpha = rep(0.75, 4)))) +
            scale_x_continuous(labels = number_format(accuracy = 0.1), limits = x_lim) +
            scale_y_continuous(labels = number_format(accuracy = 0.1), limits = y_lim) +
            labs(title = "CEP-IP plot (cDFG vs GOI) \u2014 Original labels",
                 subtitle = orig_sub,
                 x = bquote(italic(.(GOI))),
                 y = DFG_COMPOSITE_LABEL) +
            theme_p10b(base_size = 10) +
            theme(legend.position = "bottom",
                  legend.text = element_text(size = 7),
                  legend.title = element_text(size = 8, face = "bold"),
                  legend.key.size = unit(0.4, "cm"),
                  plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
                  plot.subtitle = element_text(size = 6.5, colour = "gray30",
                                               hjust = 0.5, lineheight = 1.1),
                  aspect.ratio = 1)
        cat("  Fig 4A: CEP-IP (original labels) built\n")

        # == Figure 4B: Pitcher (Ki vs Si), original labels ===================
        set.seed(123)
        p4_pitch <- ggplot() +
            geom_point(data = excl_cells,
                       aes(x = spear_norm, y = kend_contrib),
                       colour = "#D0D0D0", alpha = 0.25, size = 0.6) +
            geom_point(data = oo_cells,
                       aes(x = spear_norm, y = kend_contrib),
                       colour = "white", fill = ORIGONLY_COL,
                       shape = 25, size = 2.2, alpha = 0.88, stroke = 0.15) +
            geom_point(data = origref_plot,
                       aes(x = S, y = K, colour = group_orig),
                       size = 1.35, alpha = 0.58, shape = 16) +
            geom_point(data = recov_cells,
                       aes(x = spear_norm, y = kend_contrib),
                       colour = "white", fill = RECOV_COL,
                       shape = 23, size = 2.2, alpha = 0.93, stroke = 0.15) +
            geom_vline(xintercept = 0.05, colour = "#009E73",
                       linetype = "solid", linewidth = 0.6) +
            geom_hline(yintercept = 0.5, colour = "#009E73",
                       linetype = "solid", linewidth = 0.6) +
            scale_colour_manual(values = four_cols) +
            labs(title = expression("Pitcher plot (" * K[i] * " vs " * S[i]^norm * ") \u2014 Original labels"),
                 subtitle = orig_sub,
                 x = expression(S[i]^norm), y = expression(K[i])) +
            theme_p10b(base_size = 10) +
            theme(legend.position = "none",
                  plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
                  plot.subtitle = element_text(size = 6.5, colour = "gray30",
                                               hjust = 0.5, lineheight = 1.1),
                  aspect.ratio = 1)
        cat("  Fig 4B: Pitcher (original labels) built\n")

        # == Assemble Figure 4 ================================================
        title4_text <- sprintf(
            "%s (%s\u207A, n=%d) \u2014 Original Part10b Labels (no refit)",
            dataset_tag, GOI_LABEL_PLAIN, n_goi_pos)
        header4 <- arrangeGrob(
            textGrob(title4_text, gp = gpar(fontface = "bold", fontsize = 13)),
            linesGrob(x = unit(c(0.05, 0.95), "npc"), y = unit(c(0.3, 0.3), "npc"),
                      gp = gpar(col = "black", lwd = 1.5)),
            ncol = 1, heights = unit(c(0.7, 0.3), "cm"))
        fig4 <- arrangeGrob(header4,
                            arrangeGrob(p4_cep, p4_pitch, ncol = 2, widths = c(1, 1)),
                            ncol = 1, heights = unit(c(1.2, 10), "cm"))
        grid.newpage(); grid.draw(fig4)
        save_grob(fig4, file.path(OUT_DIR, paste0(
            FILE_PREFIX, "_", dataset_tag, "_OrigLabel_CEP-IP_Pitcher")),
            pdf_w = 14, pdf_h = 8.5, png_w = 2800, png_h = 1700)
        cat("  Fig 4 saved\n")


        # == Figure 5A: Copula (u,v), original labels =========================
        # Reuse cop_df from Fig 3 but re-colour with original labels

        # Build copula data with original group assignment
        cop_orig <- data.frame(
            u = u_cop, v = v_cop, u_sk = u_sk, v_sk = v_sk,
            barcode = goi_pos_cells$barcode,
            cell_class = goi_pos_cells$cell_class,
            stringsAsFactors = FALSE)

        # OrigRef copula subset
        cop_origref <- cop_orig[cop_orig$cell_class == "Original + Refined", ]
        cop_origref$group_orig <- factor(orig_grp_bc[cop_origref$barcode],
                                         levels = names(four_cols))
        cop_origref$draw_ord <- ord_map[match(cop_origref$group_orig, names(ord_map))]
        cop_origref$draw_ord[is.na(cop_origref$draw_ord)] <- 0L
        cop_origref <- cop_origref[order(cop_origref$draw_ord), ]

        cop_recov <- cop_orig[cop_orig$cell_class == "Recovered", ]
        cop_oo    <- cop_orig[cop_orig$cell_class == "Original only", ]
        cop_excl  <- cop_orig[cop_orig$cell_class == "Excluded", ]

        set.seed(123)
        p5_uv <- ggplot() +
            geom_point(data = cop_excl, aes(x = u, y = v),
                       colour = "gray80", alpha = 0.30, size = 0.6) +
            geom_point(data = cop_oo, aes(x = u, y = v),
                       colour = "white", fill = ORIGONLY_COL,
                       shape = 25, size = 2.0, alpha = 0.88, stroke = 0.15) +
            geom_point(data = cop_origref, aes(x = u, y = v, colour = group_orig),
                       size = 1.2, alpha = 0.55, shape = 16) +
            geom_point(data = cop_recov, aes(x = u, y = v),
                       colour = "white", fill = RECOV_COL,
                       shape = 23, size = 2.0, alpha = 0.93, stroke = 0.15) +
            # Pre/Post boundary at IP_v = Fc(IP_y).
            ip_hline_cop_v1(IP_v) +
            scale_colour_manual(values = four_cols, name = "OrigRef group",
                                drop = FALSE) +
            coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
            scale_x_continuous(breaks = seq(0, 1, 0.25)) +
            scale_y_continuous(breaks = seq(0, 1, 0.25)) +
            labs(title = "Copula space \u2014 Original labels",
                 subtitle = sprintf("u=F(GOI)  v=F(cDFG)  n=%d", n_cop),
                 x = expression(""*"" %->% "" * u),
                 y = expression(""*"" %->% "" * v)) +
            theme_p10b(base_size = 10) +
            theme(legend.position = "bottom",
                  legend.text = element_text(size = 7),
                  legend.title = element_text(size = 8, face = "bold"),
                  legend.key.size = unit(0.35, "cm"),
                  plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
                  plot.subtitle = element_text(size = 7, colour = "gray30", hjust = 0.5),
                  aspect.ratio = 1)

        # == Figure 5B: Copula (u_SK, v_SK), original labels ==================
        set.seed(123)
        p5_sk <- ggplot() +
            geom_point(data = cop_excl, aes(x = u_sk, y = v_sk),
                       colour = "gray80", alpha = 0.30, size = 0.6) +
            geom_point(data = cop_oo, aes(x = u_sk, y = v_sk),
                       colour = "white", fill = ORIGONLY_COL,
                       shape = 25, size = 2.0, alpha = 0.88, stroke = 0.15) +
            geom_point(data = cop_origref, aes(x = u_sk, y = v_sk, colour = group_orig),
                       size = 1.2, alpha = 0.55, shape = 16) +
            geom_point(data = cop_recov, aes(x = u_sk, y = v_sk),
                       colour = "white", fill = RECOV_COL,
                       shape = 23, size = 2.0, alpha = 0.93, stroke = 0.15) +
            geom_vline(xintercept = u_Sinorm05_cop, colour = "#009E73",
                       linetype = "solid", linewidth = 0.6) +
            geom_hline(yintercept = v_Ki05_cop, colour = "#009E73",
                       linetype = "solid", linewidth = 0.6) +
            scale_colour_manual(values = four_cols) +
            coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
            scale_x_continuous(breaks = seq(0, 1, 0.25)) +
            scale_y_continuous(breaks = seq(0, 1, 0.25)) +
            labs(title = "Copula space \u2014 Original labels",
                 subtitle = sprintf("u[SK]=F(Si_norm)  v[SK]=F(Ki)  n=%d", n_cop),
                 x = expression(""*"" %->% "" * u[SK]),
                 y = expression(""*"" %->% "" * v[SK])) +
            theme_p10b(base_size = 10) +
            theme(legend.position = "none",
                  plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
                  plot.subtitle = element_text(size = 7, colour = "gray30", hjust = 0.5),
                  aspect.ratio = 1)
        cat("  Fig 5: Copula (original labels) built\n")

        # == Assemble Figure 5 ================================================
        title5_text <- sprintf(
            "%s (%s\u207A, n=%d) \u2014 Original Labels: Copula Space",
            dataset_tag, GOI_LABEL_PLAIN, n_goi_pos)
        header5 <- arrangeGrob(
            textGrob(title5_text, gp = gpar(fontface = "bold", fontsize = 13)),
            linesGrob(x = unit(c(0.05, 0.95), "npc"), y = unit(c(0.3, 0.3), "npc"),
                      gp = gpar(col = "black", lwd = 1.5)),
            ncol = 1, heights = unit(c(0.7, 0.3), "cm"))
        fig5 <- arrangeGrob(header5,
                            arrangeGrob(p5_uv, p5_sk, ncol = 2, widths = c(1, 1)),
                            ncol = 1, heights = unit(c(1.2, 10), "cm"))
        grid.newpage(); grid.draw(fig5)
        save_grob(fig5, file.path(OUT_DIR, paste0(
            FILE_PREFIX, "_", dataset_tag, "_OrigLabel_Copula_space")),
            pdf_w = 14, pdf_h = 8.5, png_w = 2800, png_h = 1700)
        cat("  Fig 5 saved\n\n")

        # == Concordance: refitted TREP vs original Part 10b TREP =============
        # Only for OrigRef cells (which exist in both assignments)
        cat("  Computing TREP concordance (refitted vs original) ...\n")
        origref_bc <- origref_cells$barcode
        v3_label   <- four_group_bc[origref_bc]
        orig_label <- orig_grp_bc[origref_bc]

        concord_df <- data.frame(
            barcode    = origref_bc,
            v3_group   = v3_label,
            orig_group = orig_label,
            stringsAsFactors = FALSE)
        concord_df$agree <- concord_df$v3_group == concord_df$orig_group
        n_agree <- sum(concord_df$agree, na.rm = TRUE)
        n_total_or <- nrow(concord_df)
        # Unlabelled cells are reported; they count as disagreements.
        n_unlabelled <- sum(is.na(concord_df$agree))
        if (n_unlabelled > 0)
            cat(sprintf(paste0("    WARN: %d of %d Original+Refined cell(s)",
                               " carry no label on one side; they count as",
                               " disagreements below.\n"),
                        n_unlabelled, n_total_or))
        # NA when there are no Original+Refined cells.
        pct_agree <- if (n_total_or > 0L)
            round(100 * n_agree / n_total_or, 1) else NA_real_
        cat(sprintf("  Concordance: %d/%d (%s)\n\n", n_agree, n_total_or,
                    if (is.finite(pct_agree)) sprintf("%.1f%%", pct_agree)
                    else "no Original+Refined cells"))

        # Cross-tabulation
        xtab <- table(Original = concord_df$orig_group, Refitted = concord_df$v3_group)


        # =============================================================================
        # SECTION G: AUGMENTED DATA FRAME AND RDS BUNDLE
        # Part 10e sheet with the refitted EP, TREP and group for module cells.
        # =============================================================================
        cat("  Building V1-augmented Step4_Refined_Cells ...\n")

        # == Original Part 10e sheet (full schema) ============================
        rc_orig <- all_cells

        # == Map refitted values to barcodes ==================================
        # EP values
        ep_bc <- setNames(ep_values, df_gam$barcode)

        # Lookup tables for union cells
        v3_ep_lkp     <- ep_bc                       # EP_value
        v3_ip_lkp     <- ip_status_bc                # "Pre-IP" / "Post-IP"
        v3_group_lkp  <- four_group_bc               # 4-group
        v3_istrep_lkp <- setNames(
            as.integer(module_cells$barcode %in% trep_barcodes),
            module_cells$barcode)

        # == Update columns for union cells ===================================
        rc_v3 <- rc_orig
        rc_v3$dataset_tag_col <- dataset_tag

        # Identify union barcodes in the data frame
        bc_in_df <- rc_v3$barcode

        # EP_value: update for union cells
        rc_v3$EP_value <- ifelse(
            bc_in_df %in% names(v3_ep_lkp),
            v3_ep_lkp[bc_in_df],
            rc_v3$EP_value)

        # is_TREP: update for union cells, kept logical.
        rc_v3$is_TREP <- as.logical(ifelse(
            bc_in_df %in% names(v3_istrep_lkp),
            v3_istrep_lkp[bc_in_df],
            rc_v3$is_TREP))
        # is_trep (full-population column) is left as Part 10e wrote it.

        # IP_status: recompute for union cells (same IP threshold)
        rc_v3$IP_status <- ifelse(
            bc_in_df %in% names(v3_ip_lkp),
            v3_ip_lkp[bc_in_df],
            rc_v3$IP_status)

        # group_refined: refitted 4-group for union cells
        rc_v3$group_refined <- ifelse(
            bc_in_df %in% names(v3_group_lkp),
            v3_group_lkp[bc_in_df],
            rc_v3$group_refined)

        # == Column ordering (Part 10e CMRF schema) ===========================
        target_cols <- c(
            "barcode",         "dataset_tag_col",
            "in_union_skcs",   "is_orig_skcs",     "is_recovered",
            "cell_class",      "IP_status",         "is_TREP",
            "group_orig",      "group_refined",
            "u_all",           "v_all",             "K_density",
            "EP_value",        "Si_value",          "Ki_value",
            "g_u",
            "cond2_refined",   "is_refined",
            "p_full",          "q_full",
            "yhat_iso",        "resid_abs",         "delta_base",
            "delta_adapt",     "pass_cmrf",
            "pass_stage1",     "pass_stage2",       "removed_by_stage2",
            "region",          "GOI_expr",          "cDFG_expr",
            "spear_contrib",   "spear_norm",      "kend_contrib",
            "is_GOI_neg",      "is_trep"
        )

        # Keep target cols that exist, then append any extras
        have_cols <- target_cols[target_cols %in% names(rc_v3)]
        rest_cols <- setdiff(names(rc_v3), have_cols)
        rc_v3 <- rc_v3[, c(have_cols, rest_cols)]

        # Sort: GOI-positive first (by u_all), GOI-negative last
        rc_v3 <- rc_v3[order(is.na(rc_v3$u_all), rc_v3$u_all, na.last = TRUE), ]
        rownames(rc_v3) <- NULL

        cat(sprintf("  V1 data frame: %d rows (%d GOI-pos + %d GOI-neg)\n",
                    nrow(rc_v3),
                    sum(rc_v3$cell_class != "GOI-negative", na.rm = TRUE),
                    sum(rc_v3$cell_class == "GOI-negative", na.rm = TRUE)))

        # == Return for guardrail =============================================
        invisible(list(
            dataset_tag = dataset_tag, n_total_M4top = n_total_M4top,
            n_goi_neg = n_goi_neg, n_goi_pos = n_goi_pos, n_union = n_union,
            n_origref = n_origref, n_recovered = n_recovered,
            n_orig_only = n_orig_only, n_excluded = n_excluded,
            n_non_module = n_non_module, n_trep = n_trep, n_nontrep = n_nontrep,
            model_de = model_de, model_de_pop = model_de_pop,
            IP_x = IP_x, group_tab = group_tab,
            # Reported in the post-analysis guardrails
            IP_y = IP_y, IP_axis_used = IP_AXIS_USED, IP_active = IP_ACTIVE,
            donor_re = res_dfg$donor_re, n_donors = res_dfg$n_donors,
            best_k = res_dfg$best_k, prss_values = res_dfg$prss_values,
            n_iter_max = res_dfg$n_iter_max,
            concordance_df = concord_df, concordance_xtab = xtab,
            concordance_pct = pct_agree,
            # Objects for the RDS
            gam_model = gam_model,
            gam_best_k = res_dfg$best_k,
            gam_best_prss = res_dfg$best_prss,
            ep_values = ep_bc,
            trep_barcodes = trep_barcodes,
            nontrep_barcodes = nontrep_barcodes,
            four_group_bc = four_group_bc,
            ip_status_bc = ip_status_bc,
            trep_flag_bc = trep_flag_bc,
            copula_u = setNames(u_cop, goi_pos_cells$barcode),
            copula_v = setNames(v_cop, goi_pos_cells$barcode),
            copula_u_sk = setNames(u_sk, goi_pos_cells$barcode),
            copula_v_sk = setNames(v_sk, goi_pos_cells$barcode),
            IP_u = IP_u, IP_v = IP_v,
            # Old names kept as aliases for downstream readers.
            u_IP = IP_u, v_IP = IP_v,                  # deprecated
            IP_thresh = IP_x,                          # deprecated
            # u_Si0_cop kept as an alias of u_Sinorm05_cop.
            u_Sinorm05_cop = u_Sinorm05_cop, v_Ki05_cop = v_Ki05_cop,
            u_Si0_cop = u_Sinorm05_cop,
            Si_norm_cut = SI_NORM_CUT_V1, Ki_cut = KI_CUT_V1,
            rc_v3 = rc_v3
        ))
    }


    # =============================================================================
    # SECTION D: EXECUTE
    # =============================================================================

    all_results <- vector("list", length(DATASETS))
    names(all_results) <- DATASETS
    for (ds in DATASETS) {
        # Single-bracket assignment keeps a NULL result in the list.
        all_results[ds] <- list(tryCatch(
            run_V1_plots(ds),
            error = function(e) {
                cat("  ERROR in dataset", ds, ":", conditionMessage(e), "\n")
                NULL
            }))
    }


    # =============================================================================
    # SECTION D2: SAVE RDS AND EXCEL
    # =============================================================================
    cat("\n", strrep("=", 72), "\n")
    cat("  SAVING RDS + EXCEL\n")
    cat(strrep("=", 72), "\n\n")

    for (ds in DATASETS) {
        r <- all_results[[ds]]
        if (is.null(r)) { cat("  [", ds, "] skipped (no results)\n"); next }

        # rds_save() strips model environments and drops plot objects.
        rds_path <- file.path(OUT_DIR,
                              paste0(FILE_PREFIX, "_", ds, "_results.rds"))
        rds_save(r, rds_path, label = sprintf("[%s] RDS", ds))

        xlsx_out <- file.path(OUT_DIR, paste0(FILE_PREFIX, "_", ds,
                                              "_SK_interaction_analysis.xlsx"))
        excel_sheets <- list("Step4_Refined_Cells" = r$rc_v3,
                             "Concordance"         = r$concordance_df)
        tryCatch({
            write_xlsx(excel_sheets, path = xlsx_out)
            cat(sprintf("  [%s] Excel saved: %s (%d rows)\n",
                        ds, basename(xlsx_out), nrow(r$rc_v3)))
        }, error = function(e)
            cat("  [", ds, "] Excel error:", conditionMessage(e), "\n"))
    }
    cat("\n")
    # =============================================================================
    # SECTION E: CELL COUNT GUARDRAIL
    # =============================================================================

    cat("\n", strrep("=", 72), "\n")
    cat("  CELL COUNT GUARDRAIL\n")
    cat(strrep("=", 72), "\n\n")

    log_path <- file.path(OUT_DIR, paste0(FILE_PREFIX, "_cell_count_guardrail.txt"))
    log_con <- file(log_path, open = "wt")
    # on.exit() closes the report connection on any failure.
    on.exit(try(close(log_con), silent = TRUE), add = TRUE)
    write_log <- function(...) {
        msg <- paste0(...)
        cat(msg)
        writeLines(msg, con = log_con, sep = "")
    }

    write_log(strrep("=", 72), "\n")
    write_log(paste0("Part10e Validation1 [", ARM_TAG, "] - Cell Count Guardrail Report\n"))
    write_log(sprintf("Generated: %s\n", Sys.time()))
    write_log(strrep("=", 72), "\n\n")

    # Reconciliation identities on the counts from Step4_Refined_Cells; they
    # hold on any arm.
    all_pass <- TRUE
    for (ds in DATASETS) {
        r <- all_results[[ds]]
        if (is.null(r)) {
            write_log(sprintf("[%s]  SKIP - no results\n\n", ds))
            all_pass <- FALSE; next
        }
        write_log(sprintf("[%s] Cell Count Reconciliation  (arm %s, GOI %s)\n",
                          ds, ARM_TAG, GOI))
        write_log(strrep("-", 78), "\n")
        write_log(sprintf("%-52s %10s %10s\n", "Identity", "Observed", "Status"))
        write_log(strrep("-", 78), "\n")

        checks <- list(
            list(sprintf("OrigRef (%d) + Recovered (%d) = module (%d)",
                         r$n_origref, r$n_recovered, r$n_union),
                 r$n_origref + r$n_recovered, r$n_union),
            list(sprintf("the four CMRF classes (%d+%d+%d+%d) = %s-positive (%d)",
                         r$n_origref, r$n_recovered, r$n_orig_only,
                         r$n_excluded, GOI, r$n_goi_pos),
                 r$n_origref + r$n_recovered + r$n_orig_only + r$n_excluded,
                 r$n_goi_pos),
            list(sprintf("%s-positive (%d) + %s-negative (%d) = sheet rows (%d)",
                         GOI, r$n_goi_pos, GOI, r$n_goi_neg, r$n_total_M4top),
                 r$n_goi_pos + r$n_goi_neg, r$n_total_M4top),
            list(sprintf("module (%d) + non-module (%d) = %s-positive (%d)",
                         r$n_union, r$n_non_module, GOI, r$n_goi_pos),
                 r$n_union + r$n_non_module, r$n_goi_pos),
            list(sprintf("TREP (%d) + non-TREP (%d) = module (%d)",
                         r$n_trep, r$n_nontrep, r$n_union),
                 r$n_trep + r$n_nontrep, r$n_union))

        ds_pass <- TRUE
        for (chk in checks) {
            status <- if (identical(as.integer(chk[[2]]),
                                    as.integer(chk[[3]]))) "PASS" else "FAIL"
            if (status == "FAIL") { ds_pass <- FALSE; all_pass <- FALSE }
            write_log(sprintf("%-52s %10d %10s\n",
                              chk[[1]], as.integer(chk[[2]]), status))
        }

        write_log(strrep("-", 78), "\n")
        write_log(sprintf("Dataset %s: %s\n\n", ds,
                          if (ds_pass) "ALL RECONCILE" else "*** DOES NOT RECONCILE ***"))
        write_log(sprintf("  %s-positive (plotted):%s%d\n", GOI,
                          strrep(" ", max(1L, 24L - nchar(GOI))), r$n_goi_pos))
        write_log(sprintf("  Non-module (gray):             %d\n", r$n_non_module))
        write_log(sprintf("  GAM DE (donor-aware):          %.2f%%\n",
                          r$model_de * 100))
        write_log(sprintf("  TREP sizing DE (flat null):    %.2f%%\n",
                          .num1_v1(r$model_de_pop) * 100))
        write_log(sprintf("  IP threshold:                  %.4f\n", r$IP_x))
        write_log(sprintf("  TREP / non-TREP:               %d / %d\n", r$n_trep, r$n_nontrep))
        write_log("\n  4-group breakdown:\n")
        for (g in names(four_cols)) {
            n_g <- as.integer(r$group_tab[g])
            if (is.na(n_g)) n_g <- 0L
            write_log(sprintf("    %-25s %d\n", g, n_g))
        }

        # Concordance: refitted vs original Part 10b TREP assignment
        write_log("\n  TREP Concordance (OrigRef cells: refitted vs original Part10b):\n")
        write_log(sprintf("    Agreement: %d / %d (%s)\n",
                          sum(r$concordance_df$agree, na.rm = TRUE),
                          nrow(r$concordance_df),
                          if (is.null(r$concordance_pct) ||
                              !is.finite(r$concordance_pct)) "not measurable"
                          else sprintf("%.1f%%", r$concordance_pct)))

        xt <- r$concordance_xtab
        rn <- rownames(xt); cn <- colnames(xt)
        write_log(sprintf("    %-22s", "Original \\ Refitted"))
        for (cc in cn) write_log(sprintf(" %18s", cc))
        write_log("\n")
        for (rr in rn) {
            write_log(sprintf("    %-22s", rr))
            for (cc in cn) {
                val <- if (rr %in% rownames(xt) && cc %in% colnames(xt)) xt[rr, cc] else 0L
                write_log(sprintf(" %18d", val))
            }
            write_log("\n")
        }
        write_log("\n\n")
    }

    write_log(strrep("=", 72), "\n")
    write_log(sprintf("OVERALL GUARDRAIL: %s\n",
                      if (all_pass) "ALL COUNTS RECONCILE - PASS"
                      else "*** COUNTS DO NOT RECONCILE ***"))
    write_log(strrep("=", 72), "\n")
    # Closed here so the file is flushed; on.exit() covers the failure path.
    close(log_con)
    cat(sprintf("\n  Guardrail log saved: %s\n", log_path))


    # =============================================================================
    # SECTION F: POST-ANALYSIS GUARDRAILS AND SUMMARY
    # =============================================================================
    cat("\n  POST-ANALYSIS GUARDRAILS\n\n")
    for (ds in DATASETS) {
        r <- all_results[[ds]]
        if (is.null(r)) {
            grv1(paste0("GR-06-", ds), "FAIL", "Dataset returned NULL")
            next
        }
        grv1(paste0("GR-06-", ds), "PASS",
             sprintf(paste0("module = %d cell(s) (OrigRef %d + Recovered %d) | ",
                            "cDFG DE = %.2f%% (donor-aware) | ",
                            "TREP %d / non-TREP %d"),
                     r$n_union, r$n_origref, r$n_recovered,
                     100 * .num1_v1(r$model_de), r$n_trep, r$n_nontrep))

        # Donor random intercept on the refitted curve.
        grv1(paste0("GR-07-", ds),
             if (isTRUE(r$donor_re) || !isTRUE(GAM_DONOR_RE_V1)) "PASS" else "WARN",
             sprintf("cDFG refit: donor RE = %s (%d donors) | k = %s | %s weights",
                     if (isTRUE(r$donor_re)) "YES" else "NO",
                     if (is.null(r$n_donors)) 0L else r$n_donors,
                     if (is.null(r$best_k)) "NA" else as.character(r$best_k),
                     if (W_DONOR_WEIGHTS_V1) "inverse-donor" else "uniform"))

        # Axis of the Pre/Post split; anything other than cDFG is a FAIL.
        grv1(paste0("GR-08-", ds),
             if (identical(r$IP_axis_used, "cDFG")) "PASS" else "FAIL",
             sprintf("Pre/Post split on %s at %s | IP_x = %s | IP_y = %s",
                     if (is.null(r$IP_axis_used)) "NA" else r$IP_axis_used,
                     if (is.null(r$IP_active) || !is.finite(r$IP_active)) "NA"
                     else sprintf("%.6f", r$IP_active),
                     if (is.null(r$IP_x) || !is.finite(r$IP_x)) "NA"
                     else sprintf("%.6f", r$IP_x),
                     if (is.null(r$IP_y) || !is.finite(r$IP_y)) "NA"
                     else sprintf("%.6f", r$IP_y)))

        # The concordance between the refitted TREP call and Part 10b's.
        grv1(paste0("GR-09-", ds),
             if (is.null(r$concordance_pct) ||
                 !is.finite(r$concordance_pct)) "WARN" else "PASS",
             if (is.null(r$concordance_pct) || !is.finite(r$concordance_pct))
                 "TREP concordance not measurable: no Original+Refined cells"
             else sprintf("Refitted vs Part10b TREP: %.1f%% agreement on %d cell(s)",
                          r$concordance_pct, nrow(r$concordance_df)))
    }

    grv1_write()
    cat(sprintf("\n  Final guardrail tally: %d PASS / %d WARN / %d FAIL\n",
                .grv1_pass, .grv1_warn, .grv1_fail))
    if (.grv1_fail > 0) {
        cat("  Failures present.\n")
    } else if (.grv1_warn > 0) {
        cat(sprintf("  %d passed, %d warning(s), no failures.\n",
                    .grv1_pass, .grv1_warn))
    } else {
        cat(sprintf("  All %d checks passed.\n", .grv1_pass))
    }

    cat("\n", strrep("=", 72), "\n")
    cat("  Part10e Validation1 COMPLETE [", ARM_TAG, "]\n")
    cat("  Output directory:", OUT_DIR, "\n")
    cat("  Files generated:\n")
    # Lists the output files, excluding transient ones.
    .transient <- c("_RUN_INCOMPLETE.txt", basename(MSG_PATH_V1))
    for (f in setdiff(sort(list.files(OUT_DIR)), .transient)) {
        sz <- round(file.size(file.path(OUT_DIR, f)) / 1e3, 0)
        cat(sprintf("    %-62s  %5d KB\n", f, sz))
    }
    cat(strrep("=", 72), "\n")

    cat("\n  Log written to:", LOG_PATH_V1, "\n")
    cat(sprintf("\nPart10e Validation1 - log ended: %s\n", Sys.time()))

    .p10v_completed <- TRUE
    if (sink.number(type = "message") != 2L) sink(type = "message")
    while (sink.number() > 0) sink()
    suppressWarnings(file.remove(INCOMPLETE_MARKER_V1))

    invisible(list(arm = arm_tag, dir = OUT_DIR, pass = .grv1_pass,
                   warn = .grv1_warn, fail = .grv1_fail))
}


# =============================================================================
# SECTION 0c: Driver
# =============================================================================

cep_run_arms(PART_DIR_TAG, run_one_arm)
