# =============================================================================
# Part9b_CP10K_SiKi_code.r
# SEA-AD CEP-IP/CMRF, Part 9b, Si_norm floor (Supplementary Figure S4, Table S7)
#   Justifies the Si_norm > 0.05 floor on the primary arm (Mic13_PVM3). SKCS
#   scores are computed from the Part 9 GAM data; the waist near Si_norm = 0 is
#   located by KDE, rank-based bounds and a GMM, and GOI rank dispersion at the
#   threshold is compared with the waist (Kruskal-Wallis, Dunn's, Kendall tau).
#
# Usage (from the SEA-AD folder, or set SEAAD_ROOT), after Part 9
#   source("Part9b_CP10K_SiKi_code.r")
#
# Outputs (Part9b/Part9b_Mic13_PVM3/): figures (PNG), threshold workbook,
#   results RDS and log.
# =============================================================================


# =============================================================================
# SECTION 0: Libraries
# =============================================================================

suppressPackageStartupMessages({
    library(ggplot2)
    library(mclust)
    library(openxlsx)
    library(grid)
    library(dunn.test)
})


# =============================================================================
# SECTION 0a: Shared arm definitions
# =============================================================================

if (!exists("CEP_HOME")) CEP_HOME <- Sys.getenv("SEAAD_ROOT", getwd())
.p9b_harness <- file.path(CEP_HOME, "CEP_arms.R")
if (!file.exists(.p9b_harness))
    stop("CEP_arms.R not found at: ", .p9b_harness,
         call. = FALSE)
source(.p9b_harness)

PART_TAG <- "Part9b"


# =============================================================================
# SECTION 0b: Helper functions (pure, outside run_one_arm)
# =============================================================================

## Shapiro-Wilk with subsampling for large n
safe_shapiro <- function(x, max_n = 5000) {
    x <- x[!is.na(x)]
    n <- length(x)
    if (n < 3) return(list(W = NA, p.value = NA, n_tested = n))
    test_x <- if (n > max_n) sample(x, max_n) else x
    sw <- shapiro.test(test_x)
    list(W = sw$statistic, p.value = sw$p.value, n_tested = length(test_x))
}

## Two-group comparison: Shapiro-Wilk -> t-test or Mann-Whitney U
two_group_test <- function(g1, g2, alternative = "two.sided") {
    sw1 <- safe_shapiro(g1)
    sw2 <- safe_shapiro(g2)
    both_normal <- (!is.na(sw1$p.value) && sw1$p.value >= 0.05) &&
                   (!is.na(sw2$p.value) && sw2$p.value >= 0.05)
    if (both_normal) {
        tt <- t.test(g1, g2, alternative = alternative)
        return(list(method = "t-test", statistic = tt$statistic,
                    p.value = tt$p.value, sw1 = sw1, sw2 = sw2))
    } else {
        mw <- wilcox.test(g1, g2, alternative = alternative)
        return(list(method = "Mann-Whitney U", statistic = mw$statistic,
                    p.value = mw$p.value, sw1 = sw1, sw2 = sw2))
    }
}

## GMM posterior: P(component k | x) for all k
gmm_posterior <- function(x, means, sds, props) {
    wd <- props * dnorm(x, means, sds)
    total <- sum(wd)
    if (total < .Machine$double.eps) return(rep(0, length(means)))
    wd / total
}

## Intersection of two weighted Gaussians - all real roots
gmm_pairwise_intersection <- function(mu1, sd1, pi1, mu2, sd2, pi2) {
    A <- 1/(2*sd2^2) - 1/(2*sd1^2)
    B <- mu1/sd1^2 - mu2/sd2^2
    C <- mu2^2/(2*sd2^2) - mu1^2/(2*sd1^2) - log(pi2/pi1) - log(sd1/sd2)
    if (abs(A) < 1e-12) {
        if (abs(B) < 1e-12) return(numeric(0))
        return(-C / B)
    }
    disc <- B^2 - 4*A*C
    if (disc < 0) return(numeric(0))
    sort(c((-B + sqrt(disc)) / (2*A), (-B - sqrt(disc)) / (2*A)))
}

## Between-means adjacent intersections for a GMM
gmm_adjacent_intersections <- function(params) {
    ord     <- order(params$mean)
    means_s <- params$mean[ord]; sds_s <- params$sd[ord]; props_s <- params$prop[ord]
    result <- data.frame(pair = character(), intersection = numeric(),
                         stringsAsFactors = FALSE)
    for (i in 1:(length(means_s) - 1)) {
        xs <- gmm_pairwise_intersection(means_s[i], sds_s[i], props_s[i],
                                         means_s[i+1], sds_s[i+1], props_s[i+1])
        between <- xs[xs > means_s[i] & xs < means_s[i+1]]
        for (b in between)
            result <- rbind(result, data.frame(
                pair = sprintf("%.4f-%.4f", means_s[i], means_s[i+1]),
                intersection = b, stringsAsFactors = FALSE))
    }
    result
}

## Find x where P(waist) = target via uniroot
find_posterior_boundary <- function(target, waist_idx, means, sds, props,
                                    lo = 0.0, hi = 0.20) {
    f <- function(x) gmm_posterior(x, means, sds, props)[waist_idx] - target
    if (f(lo) < 0) return(lo)
    if (f(hi) > 0) return(NA_real_)
    uniroot(f, c(lo, hi), tol = 1e-8)$root
}

## Extract GMM parameters as data.frame
extract_gmm_params <- function(fit) {
    k <- fit$G
    means <- fit$parameters$mean
    sds   <- sqrt(fit$parameters$variance$sigmasq)
    props <- fit$parameters$pro
    if (length(sds) == 1) sds <- rep(sds, k)
    data.frame(comp = 1:k, mean = means, sd = sds, prop = props)
}


# =============================================================================
# SECTION 0c: Inline SKCS computation (lightweight)
# =============================================================================
# Per-cell Si_norm (rank-product deviation / Si_max) and Ki (tau_b
# contribution) from the GOI and cDFG vectors, as in Part 10.

compute_skcs_p9b <- function(goi_vec, cdfg_vec, label = "",
                              si_thresh = 0.05, ki_thresh = 0.50) {
    n <- length(goi_vec)
    stopifnot(length(cdfg_vec) == n, n >= 3)

    # -- Spearman contribution --
    rk_goi  <- rank(goi_vec,  ties.method = "average")
    rk_cdfg <- rank(cdfg_vec, ties.method = "average")
    mr_goi  <- mean(rk_goi)
    mr_cdfg <- mean(rk_cdfg)
    spear_contrib <- (rk_goi - mr_goi) * (rk_cdfg - mr_cdfg) / (n - 1)
    si_max      <- (n - 1) / 4
    spear_norm  <- spear_contrib / si_max

    global_rho <- cor(goi_vec, cdfg_vec, method = "spearman")
    cat(sprintf("  [%s] Global Spearman rho = %.4f  |  Si_max = %.2f\n",
                label, global_rho, si_max))

    # -- Kendall contribution (tau_b, tie-adjusted) --
    KENDALL_MAX_N <- 2000L
    KENDALL_CHUNK <- 500L

    if (n <= KENDALL_MAX_N) {
        sign_goi  <- outer(goi_vec,  goi_vec,  "-")
        sign_cdfg <- outer(cdfg_vec, cdfg_vec, "-")
        sg <- sign(sign_goi);  sc <- sign(sign_cdfg)
        prod_sc   <- sg * sc
        concordant <- prod_sc > 0;  discordant <- prod_sc < 0
        diag(concordant) <- FALSE;  diag(discordant) <- FALSE
        tied_goi  <- sg == 0; diag(tied_goi)  <- FALSE
        tied_cdfg <- sc == 0; diag(tied_cdfg) <- FALSE
        n0_goi    <- rowSums(tied_goi)
        n0_cdfg   <- rowSums(tied_cdfg)
        Ci_minus_Di <- rowSums(concordant) - rowSums(discordant)
        denom_b <- sqrt((n - 1 - n0_goi) * (n - 1 - n0_cdfg))
        denom_b[denom_b == 0] <- 1
        kend_contrib <- Ci_minus_Di / denom_b
        rm(sign_goi, sign_cdfg, sg, sc, prod_sc,
           concordant, discordant, tied_goi, tied_cdfg)
        gc(verbose = FALSE)
    } else {
        cat(sprintf("  [%s] n=%d > %d -- chunked Kendall (chunk=%d)\n",
                    label, n, KENDALL_MAX_N, KENDALL_CHUNK))
        kend_contrib <- numeric(n)
        n0_goi_vec   <- numeric(n)
        n0_cdfg_vec  <- numeric(n)
        chunk_starts <- seq(1L, n, by = KENDALL_CHUNK)
        for (cs in chunk_starts) {
            ce  <- min(cs + KENDALL_CHUNK - 1L, n)
            idx <- cs:ce
            sg_blk <- sign(outer(goi_vec[idx],  goi_vec,  "-"))
            sc_blk <- sign(outer(cdfg_vec[idx], cdfg_vec, "-"))
            prod   <- sg_blk * sc_blk
            for (ii in seq_along(idx)) prod[ii, idx[ii]] <- 0L
            kend_contrib[idx] <- rowSums(prod > 0L) - rowSums(prod < 0L)
            tied_g <- sg_blk == 0L
            tied_c <- sc_blk == 0L
            for (ii in seq_along(idx)) { tied_g[ii, idx[ii]] <- FALSE
                                          tied_c[ii, idx[ii]] <- FALSE }
            n0_goi_vec[idx]  <- rowSums(tied_g)
            n0_cdfg_vec[idx] <- rowSums(tied_c)
            rm(sg_blk, sc_blk, prod, tied_g, tied_c); gc(verbose = FALSE)
        }
        Ci_minus_Di <- kend_contrib
        denom_b <- sqrt((n - 1 - n0_goi_vec) * (n - 1 - n0_cdfg_vec))
        denom_b[denom_b == 0] <- 1
        kend_contrib <- Ci_minus_Di / denom_b
    }

    n_tied <- sum(denom_b != (n - 1))
    cat(sprintf("  [%s] Using tau_b denominator | cells with ties: %d / %d\n",
                label, n_tied, n))

    global_tau <- cor(goi_vec, cdfg_vec, method = "kendall")
    cat(sprintf("  [%s] Global Kendall tau  = %.4f\n", label, global_tau))

    is_skcs <- (spear_norm > si_thresh) & (kend_contrib > ki_thresh)
    n_skcs  <- sum(is_skcs)
    cat(sprintf("  [%s] SKCS-core: %d (%.1f%%)  |  Non-SKCS: %d (%.1f%%)\n",
                label, n_skcs, 100 * n_skcs / n,
                n - n_skcs, 100 * (n - n_skcs) / n))

    list(spear_norm   = spear_norm,
         kend_contrib = kend_contrib,
         is_skcs      = is_skcs,
         global_rho   = global_rho,
         global_tau   = global_tau)
}


# =============================================================================
# SECTION 0d: The pipeline for one arm
# =============================================================================

run_one_arm <- function(arm_tag, config_file, out_dir, file_tag) {

    config_path <- file.path(CEP_HOME, config_file)
    if (!file.exists(config_path))
        stop("Configuration not found: ", config_path, call. = FALSE)

    cat("\n", strrep("=", 76), "\n", sep = "")
    cat("  ARM: ", arm_tag, "\n", sep = "")
    cat("  Config: ", config_path, "\n", sep = "")
    cat(strrep("=", 76), "\n", sep = "")

    assign("CEP_ARM", arm_tag, envir = globalenv())
    on.exit(suppressWarnings(rm("CEP_ARM", envir = globalenv())), add = TRUE)

    source(config_path)

    # == Confirm the config loaded correctly ===================================
    .p9b_needed <- c("ARM_TAG", "NORM_TAG", "GOI", "DFG")
    .p9b_absent <- .p9b_needed[!vapply(.p9b_needed, exists, logical(1))]
    if (length(.p9b_absent) > 0)
        stop("Configuration loaded from\n    ", config_path,
             "\n  does not define: ", paste(.p9b_absent, collapse = ", "),
             call. = FALSE)

    if (!identical(ARM_TAG, arm_tag))
        stop("ARM_TAG '", ARM_TAG, "' does not match arm '", arm_tag,
             "'. Check ", config_file, ".", call. = FALSE)

    cat("  Config OK: ", ARM_TAG, " | ", GOI, " + ", length(DFG), " DFG | ",
        NORM_TAG, "\n\n", sep = "")

    # == Output paths ==========================================================
    PART9B_DIR <- out_dir
    if (!dir.exists(PART9B_DIR)) dir.create(PART9B_DIR, recursive = TRUE)

    OUT_PREFIX <- cep_out_prefix(PART_TAG, NORM_TAG, file_tag)
    LOG_PATH   <- file.path(PART9B_DIR, paste0(OUT_PREFIX, "_SiKi_Log.txt"))

    # == Log file ==============================================================
    while (sink.number() > 0) sink()
    log_con <- file(LOG_PATH, open = "wt")
    sink(log_con, append = FALSE, split = TRUE)
    on.exit({
        tryCatch({
            while (sink.number() > 0) sink()
            if (isOpen(log_con)) close(log_con)
        }, error = function(e) invisible(NULL))
    }, add = TRUE)


    # == Constants =============================================================
    COL_SI <- "spear_norm"
    COL_KI <- "kend_contrib"
    SKCS_SI_THRESH <- 0.05
    SKCS_KI_THRESH <- 0.50
    DATASETS <- c("NonNeoC", "NeoC")
    GMM_G_RANGE <- 1:6
    PLOT_W <- 10; PLOT_H <- 6; PLOT_DPI <- 150

    # == Plot helper ===========================================================
    save_plot <- function(filename, plot_obj, width = PLOT_W, height = PLOT_H) {
        filepath <- file.path(PART9B_DIR, filename)
        png(filepath, width = width, height = height,
            units = "in", res = PLOT_DPI)
        if (inherits(plot_obj, "gtable") || inherits(plot_obj, "grob")) {
            grid::grid.draw(plot_obj)
        } else {
            print(plot_obj)
        }
        dev.off()
        cat(sprintf("    Saved: %s\n", filename))
    }


    # ==========================================================================
    # Banner
    # ==========================================================================

    cat("\n", strrep("=", 72), "\n")
    cat("  Part9b Si/Ki threshold justification\n")
    cat("  Arm     : ", ARM_TAG, "\n")
    cat("  Norm    : ", NORM_TAG, "\n")
    cat("  Output  : ", PART9B_DIR, "\n")
    cat(strrep("=", 72), "\n")


    # ==========================================================================
    # Load Part9 RDS + compute SKCS inline
    # ==========================================================================
    # Part 9 RDS from the CEP_arms layout.

    arm_subdir   <- CEP_ARM_TABLE$subdir[CEP_ARM_TABLE$arm == arm_tag]
    part9_arm_dir <- cep_arm_dir("Part9", arm_subdir)
    part9_prefix  <- cep_out_prefix("Part9", NORM_TAG, file_tag)

    part9_rds_paths <- list(
        NeoC    = file.path(part9_arm_dir,
                            paste0(part9_prefix, "_NeoC_GAM_models.rds")),
        NonNeoC = file.path(part9_arm_dir,
                            paste0(part9_prefix, "_NonNeoC_GAM_models.rds"))
    )

    cat("\n  Part9 RDS:\n")
    for (nm in names(part9_rds_paths))
        cat("    ", nm, ": ", part9_rds_paths[[nm]], "\n")

    for (nm in names(part9_rds_paths)) {
        if (!file.exists(part9_rds_paths[[nm]]))
            stop("Part9 RDS not found: ", part9_rds_paths[[nm]], call. = FALSE)
    }

    cat("\n  Loading Part9 RDS and computing SKCS inline...\n")

    p9_rds_raw <- list(
        NonNeoC = readRDS(part9_rds_paths[["NonNeoC"]]),
        NeoC    = readRDS(part9_rds_paths[["NeoC"]])
    )

    p9b_data <- list()

    for (ds in DATASETS) {
        p9 <- p9_rds_raw[[ds]]
        rd <- p9$region_data
        regions_ds <- p9$regions

        cat(sprintf("\n  Pooling %s: %d regions (%s)\n",
                    ds, length(regions_ds), paste(regions_ds, collapse = ", ")))

        goi_list  <- list()
        expr_list <- list()
        for (reg in regions_ds) {
            if (is.null(rd[[reg]])) {
                cat(sprintf("    %s: NULL (skipped)\n", reg))
                next
            }
            goi_list[[reg]]  <- rd[[reg]]$goi_vec
            expr_list[[reg]] <- rd[[reg]]$dfg$Expression
            cat(sprintf("    %s: %d GOI-positive cells\n",
                        reg, length(rd[[reg]]$goi_vec)))
        }
        goi_pooled  <- unlist(goi_list,  use.names = FALSE)
        expr_pooled <- unlist(expr_list, use.names = FALSE)
        cat(sprintf("  %s pooled: %d cells\n", ds, length(goi_pooled)))

        skcs <- compute_skcs_p9b(goi_pooled, expr_pooled,
                                  paste0(ds, "_Pooled"),
                                  si_thresh = SKCS_SI_THRESH,
                                  ki_thresh = SKCS_KI_THRESH)

        score_pooled_df <- data.frame(
            GOI_expr     = goi_pooled,
            cDFG_expr    = expr_pooled,
            spear_norm   = skcs$spear_norm,
            kend_contrib = skcs$kend_contrib,
            is_GOI_neg   = FALSE,
            is_skcs      = skcs$is_skcs,
            stringsAsFactors = FALSE
        )

        p9b_data[[ds]] <- list(score_pooled_df = score_pooled_df)
    }

    rm(p9_rds_raw); gc(verbose = FALSE)
    cat("  SKCS computation complete.\n")


    # ==========================================================================
    # Storage
    # ==========================================================================

    all_synth_rows   <- list()
    all_results      <- list()
    all_posterior_df  <- data.frame()
    all_pvals        <- data.frame()


    # ==========================================================================
    # Main loop
    # ==========================================================================

    for (ds in DATASETS) {

        cat("\n\n", strrep("=", 72), "\n")
        cat(sprintf("  Dataset: %s\n", ds))
        cat(strrep("=", 72), "\n")

        # == Load and filter ==
        score_df <- p9b_data[[ds]]$score_pooled_df
        score_df <- score_df[!score_df$is_GOI_neg, ]

        si_norm <- score_df[[COL_SI]]
        ki      <- score_df[[COL_KI]]
        n_cells <- length(si_norm)

        cat(sprintf("\n  %d GOI-positive cells\n", n_cells))
        cat(sprintf("  Si_norm: [%.4f, %.4f]  median=%.4f\n",
                    min(si_norm), max(si_norm), median(si_norm)))
        cat(sprintf("  Ki:      [%.4f, %.4f]  median=%.4f\n",
                    min(ki), max(ki), median(ki)))

        pval_collector <- data.frame(Test = character(), Method = character(),
                                      p_value = numeric(),
                                      stringsAsFactors = FALSE)


        # ======================================================================
        # SECTION 1: Empirical waist boundary (density gap)
        # ======================================================================

        cat("\n  -- Section 1: Empirical waist boundary --\n")

        abs_si <- abs(si_norm)

        kde_data <- abs_si[abs_si <= 0.15]
        kde <- density(kde_data, from = 0, to = 0.15, n = 1024, bw = "SJ")

        dy <- diff(kde$y)
        sign_dy <- sign(dy)
        sign_changes <- diff(sign_dy)

        peak_idx   <- which(sign_changes <= -1) + 1
        trough_idx <- which(sign_changes >= 1) + 1

        waist_bnd <- NA_real_
        waist_bnd_method <- "default"

        if (length(peak_idx) > 0) {
            waist_peak_x <- kde$x[peak_idx[1]]
            waist_peak_y <- kde$y[peak_idx[1]]

            cat(sprintf("    Waist peak: |Si_norm| = %.4f (density = %.2f)\n",
                        waist_peak_x, waist_peak_y))

            troughs_after <- trough_idx[trough_idx > peak_idx[1]]

            if (length(troughs_after) > 0) {
                waist_bnd <- kde$x[troughs_after[1]]
                waist_bnd_density <- kde$y[troughs_after[1]]
                waist_bnd_method <- "density trough"
                cat(sprintf(
                    "    Density trough at |Si_norm| = %.4f (density = %.2f)\n",
                    waist_bnd, waist_bnd_density))
                cat(sprintf("    Density drop: %.1f%% from peak to trough\n",
                            (1 - waist_bnd_density / waist_peak_y) * 100))
            } else {
                d2y <- diff(dy)
                inflect_idx <- which(diff(sign(d2y)) >= 1) + 2
                inflects_after <- inflect_idx[inflect_idx > peak_idx[1] &
                                               inflect_idx <= length(kde$x)]
                if (length(inflects_after) > 0) {
                    waist_bnd <- kde$x[inflects_after[1]]
                    waist_bnd_method <- "density inflection"
                    cat(sprintf(
                        "    No trough; using inflection at |Si_norm| = %.4f\n",
                        waist_bnd))
                }
            }
        }

        if (is.na(waist_bnd)) {
            waist_bnd <- 0.01
            waist_bnd_method <- "default (no gap detected)"
            cat("    Warning: no density gap detected; using default 0.01\n")
        }

        n_waist_emp <- sum(abs_si < waist_bnd)
        cat(sprintf(
            "    >>> Empirical waist boundary: |Si_norm| < %.4f (%s) <<<\n",
            waist_bnd, waist_bnd_method))
        cat(sprintf("    Waist cells: %d (%.1f%% of total)\n",
                    n_waist_emp, n_waist_emp / n_cells * 100))

        # Waist boundary plot
        df_kde <- data.frame(x = kde$x, y = kde$y)

        p_waist <- ggplot(df_kde, aes(x = x, y = y)) +
            geom_line(linewidth = 1, colour = "#2166AC") +
            geom_area(alpha = 0.2, fill = "#2166AC") +
            geom_vline(xintercept = waist_bnd, linetype = "dashed",
                       colour = "#2C6E49", linewidth = 0.7) +
            annotate("text", x = waist_bnd + 0.002, y = max(kde$y) * 0.85,
                     label = sprintf("Empirical = %.4f\n(%s)",
                                     waist_bnd, waist_bnd_method),
                     colour = "#2C6E49", size = 3.2, hjust = 0,
                     fontface = "bold") +
            geom_vline(xintercept = 0.01, linetype = "dotted",
                       colour = "grey50", linewidth = 0.5) +
            annotate("text", x = 0.01 + 0.002, y = max(kde$y) * 0.65,
                     label = "Old default = 0.01",
                     colour = "grey50", size = 2.8, hjust = 0) +
            labs(
                title = sprintf(
                    "%s - Empirical waist boundary via density gap", ds),
                subtitle = sprintf(
                    paste0("|Si_norm| < %.4f defines the waist (%s) | ",
                           "%d cells (%.1f%%) | KDE bandwidth = %.4f (SJ)\n",
                           "The density trough between the waist spike and ",
                           "the broader distribution provides a data-driven ",
                           "boundary, replacing the arbitrary 0.01 convention"),
                    waist_bnd, waist_bnd_method,
                    n_waist_emp, n_waist_emp / n_cells * 100,
                    kde$bw),
                x = "|Si_norm|",
                y = "Density"
            ) +
            coord_cartesian(xlim = c(0, 0.10)) +
            theme_bw(base_size = 12) +
            theme(plot.title = element_text(face = "bold"),
                  plot.subtitle = element_text(size = 9, colour = "grey30"))

        print(p_waist)
        save_plot(sprintf("%s_%s_waist_boundary.png", OUT_PREFIX, ds), p_waist)


        # ======================================================================
        # SECTION 1b: Rank-based waist boundary (nonparametric)
        # ======================================================================

        cat("\n  -- Section 1b: Rank-based waist boundary --\n")

        goi_expr   <- score_df$GOI_expr
        goi_ranks  <- rank(goi_expr, ties.method = "average")
        median_rank <- (n_cells + 1) / 2
        rank_dev   <- abs(goi_ranks - median_rank)

        cat(sprintf("    n_cells = %d | median_rank = %.1f\n",
                    n_cells, median_rank))

        K_A <- 1; K_B <- sqrt(n_cells); K_C <- sqrt(n_cells / 12)

        cat(sprintf("    Option A (resolution limit): K = %.1f ranks\n", K_A))
        cat(sprintf(
            "    Option B (rank noise floor): K = %.1f ranks (sqrt(n))\n", K_B))
        cat(sprintf(
            "    Option C (rank SE):          K = %.1f ranks (sqrt(n/12))\n",
            K_C))

        wb_A_theory <- 2 * K_A / (n_cells - 1)
        wb_B_theory <- 2 * K_B / (n_cells - 1)
        wb_C_theory <- 2 * K_C / (n_cells - 1)

        cat(sprintf("    Theoretical bounds (2K/(n-1)):\n"))
        cat(sprintf("      A: %.6f | B: %.6f | C: %.6f\n",
                    wb_A_theory, wb_B_theory, wb_C_theory))

        cells_in_A <- rank_dev <= K_A
        cells_in_B <- rank_dev <= K_B
        cells_in_C <- rank_dev <= K_C

        wb_A_emp <- if (any(cells_in_A & abs_si > 0))
                        max(abs_si[cells_in_A]) else wb_A_theory
        wb_B_emp <- if (any(cells_in_B & abs_si > 0))
                        max(abs_si[cells_in_B]) else wb_B_theory
        wb_C_emp <- if (any(cells_in_C & abs_si > 0))
                        max(abs_si[cells_in_C]) else wb_C_theory

        n_in_A <- sum(cells_in_A)
        n_in_B <- sum(cells_in_B)
        n_in_C <- sum(cells_in_C)

        cat(sprintf(
            "    Empirical bounds (max|Si| for cells within K ranks):\n"))
        cat(sprintf("      A: %.6f (%d cells, %.1f%%)\n",
                    wb_A_emp, n_in_A, 100 * n_in_A / n_cells))
        cat(sprintf("      B: %.6f (%d cells, %.1f%%)\n",
                    wb_B_emp, n_in_B, 100 * n_in_B / n_cells))
        cat(sprintf("      C: %.6f (%d cells, %.1f%%)\n",
                    wb_C_emp, n_in_C, 100 * n_in_C / n_cells))

        n_recov_A   <- sum(abs_si > wb_A_emp & abs_si < SKCS_SI_THRESH)
        n_recov_B   <- sum(abs_si > wb_B_emp & abs_si < SKCS_SI_THRESH)
        n_recov_C   <- sum(abs_si > wb_C_emp & abs_si < SKCS_SI_THRESH)
        n_recov_KDE <- sum(abs_si > waist_bnd & abs_si < SKCS_SI_THRESH)

        cat(sprintf("\n    Recovery zone size (waist_bnd to 0.05):\n"))
        cat(sprintf("      KDE (%s):       waist=%.4f | recovery = %d cells\n",
                    waist_bnd_method, waist_bnd, n_recov_KDE))
        cat(sprintf(
            "      Option A (K=1):        waist=%.4f | recovery = %d cells\n",
            wb_A_emp, n_recov_A))
        cat(sprintf(
            "      Option B (sqrt(n)):    waist=%.4f | recovery = %d cells\n",
            wb_B_emp, n_recov_B))
        cat(sprintf(
            "      Option C (sqrt(n/12)): waist=%.4f | recovery = %d cells\n",
            wb_C_emp, n_recov_C))

        B_to_thresh_gap <- SKCS_SI_THRESH - wb_B_emp
        cat(sprintf(
            "\n    Option B gap to threshold: %.4f (%.1f%% of threshold)\n",
            B_to_thresh_gap, 100 * B_to_thresh_gap / SKCS_SI_THRESH))
        if (B_to_thresh_gap < 0.01)
            cat(
              "    Warning: Option B waist nearly reaches the Si_norm threshold\n"
            )
        if (wb_B_emp > SKCS_SI_THRESH)
            cat(
              "    Warning: Option B waist exceeds the Si_norm threshold\n"
            )

        wb_comparison <- data.frame(
            Option = c("KDE density trough", "A: Resolution (K=1)",
                       "B: Noise floor (sqrt(n))",
                       "C: Rank SE (sqrt(n/12))"),
            K_ranks = c(NA, K_A, K_B, K_C),
            Theory_bound = round(
                c(NA, wb_A_theory, wb_B_theory, wb_C_theory), 6),
            Empirical_bound = round(
                c(waist_bnd, wb_A_emp, wb_B_emp, wb_C_emp), 6),
            N_waist_cells = c(n_waist_emp, n_in_A, n_in_B, n_in_C),
            Pct_of_total = round(c(100 * n_waist_emp / n_cells,
                                   100 * n_in_A / n_cells,
                                   100 * n_in_B / n_cells,
                                   100 * n_in_C / n_cells), 1),
            N_recovery_zone = c(n_recov_KDE, n_recov_A,
                                n_recov_B, n_recov_C),
            stringsAsFactors = FALSE
        )
        cat("\n    Waist boundary comparison:\n")
        print(wb_comparison, row.names = FALSE)

        # == Rank-based waist boundary plot ==
        p_waist_rank <- ggplot(df_kde, aes(x = x, y = y)) +
            geom_line(linewidth = 1, colour = "#2166AC") +
            geom_area(alpha = 0.15, fill = "#2166AC") +
            geom_vline(xintercept = waist_bnd, linetype = "dashed",
                       colour = "#2C6E49", linewidth = 0.7) +
            annotate("text", x = waist_bnd, y = max(kde$y) * 0.95,
                     label = sprintf("KDE = %.4f", waist_bnd),
                     colour = "#2C6E49", size = 3, hjust = -0.1,
                     fontface = "bold") +
            geom_vline(xintercept = wb_A_emp, linetype = "dotted",
                       colour = "#D6604D", linewidth = 0.6) +
            annotate("text", x = wb_A_emp, y = max(kde$y) * 0.85,
                     label = sprintf("A (K=1) = %.4f", wb_A_emp),
                     colour = "#D6604D", size = 2.8, hjust = -0.1) +
            geom_vline(xintercept = min(wb_B_emp, 0.098), linetype = "longdash",
                       colour = "#4393C3", linewidth = 0.7) +
            annotate("text", x = min(wb_B_emp, 0.096), y = max(kde$y) * 0.75,
                     label = sprintf("B (sqrt(n)) = %.4f", wb_B_emp),
                     colour = "#4393C3", size = 2.8, hjust = -0.1) +
            geom_vline(xintercept = wb_C_emp, linetype = "twodash",
                       colour = "#A0522D", linewidth = 0.7) +
            annotate("text", x = wb_C_emp, y = max(kde$y) * 0.65,
                     label = sprintf("C (sqrt(n/12)) = %.4f", wb_C_emp),
                     colour = "#A0522D", size = 2.8, hjust = -0.1) +
            geom_vline(xintercept = SKCS_SI_THRESH, linetype = "solid",
                       colour = "#6B5B8A", linewidth = 0.5) +
            annotate("text", x = SKCS_SI_THRESH, y = max(kde$y) * 0.55,
                     label = "SKCS = 0.05",
                     colour = "#6B5B8A", size = 2.8, hjust = -0.1) +
            labs(
                title = sprintf(
                    "%s - Waist boundary: KDE vs rank-based options", ds),
                subtitle = sprintf(
                    paste0("KDE=%.4f | A(K=1)=%.4f | B(sqrt(n)=%.0f)=%.4f",
                           " | C(sqrt(n/12)=%.0f)=%.4f\n",
                           "Recovery zone: KDE=%d | A=%d | B=%d | C=%d | ",
                           "B gap to threshold = %.4f"),
                    waist_bnd, wb_A_emp, K_B, wb_B_emp, K_C, wb_C_emp,
                    n_recov_KDE, n_recov_A, n_recov_B, n_recov_C,
                    B_to_thresh_gap),
                x = "|Si_norm|", y = "Density"
            ) +
            coord_cartesian(xlim = c(0, 0.10)) +
            theme_bw(base_size = 12) +
            theme(plot.title = element_text(face = "bold"),
                  plot.subtitle = element_text(size = 9, colour = "grey30"))

        print(p_waist_rank)
        save_plot(sprintf("%s_%s_waist_rankbased.png", OUT_PREFIX, ds),
                  p_waist_rank, width = 10, height = 6)

        # == Scatter: |Si_norm| vs rank deviation ==
        df_rank_si <- data.frame(rank_dev = rank_dev, abs_si = abs_si)

        p_rank_scatter <- ggplot(df_rank_si,
                                 aes(x = rank_dev, y = abs_si)) +
            geom_point(alpha = 0.25, size = 0.6, colour = "grey40") +
            geom_abline(intercept = 0, slope = 2 / (n_cells - 1),
                        colour = "#D6604D", linewidth = 0.8,
                        linetype = "dashed") +
            annotate("text", x = max(rank_dev) * 0.7,
                     y = 2 * max(rank_dev) * 0.7 / (n_cells - 1) + 0.02,
                     label = "|Si|_max = 2d/(n-1)",
                     colour = "#D6604D", size = 3.2, fontface = "italic") +
            geom_vline(xintercept = K_C, linetype = "twodash",
                       colour = "#A0522D", linewidth = 0.5) +
            annotate("text", x = K_C, y = max(abs_si) * 0.95,
                     label = sprintf("C: sqrt(n/12)=%.0f", K_C),
                     colour = "#A0522D", size = 2.8, hjust = -0.1) +
            geom_vline(xintercept = K_B, linetype = "longdash",
                       colour = "#4393C3", linewidth = 0.5) +
            annotate("text", x = K_B, y = max(abs_si) * 0.85,
                     label = sprintf("B: sqrt(n)=%.0f", K_B),
                     colour = "#4393C3", size = 2.8, hjust = -0.1) +
            geom_hline(yintercept = SKCS_SI_THRESH, linetype = "solid",
                       colour = "#6B5B8A", linewidth = 0.5) +
            annotate("text", x = max(rank_dev) * 0.95,
                     y = SKCS_SI_THRESH,
                     label = "Si = 0.05",
                     colour = "#6B5B8A", size = 2.8, vjust = -0.5,
                     hjust = 1) +
            labs(
                title = sprintf(
                    "%s - |Si_norm| vs GOI rank deviation from median", ds),
                subtitle = sprintf(
                    paste0("Theoretical bound (dashed red): ",
                           "|Si|_max = 2d/(n-1) | n = %d\n",
                           "Cells within K ranks of median are ",
                           "mechanically compressed"),
                    n_cells),
                x = "GOI rank deviation from median |d_i|",
                y = "|Si_norm|"
            ) +
            theme_bw(base_size = 12) +
            theme(plot.title = element_text(face = "bold"),
                  plot.subtitle = element_text(size = 9, colour = "grey30"))

        print(p_rank_scatter)
        save_plot(sprintf("%s_%s_rank_vs_Si.png", OUT_PREFIX, ds),
                  p_rank_scatter, width = 10, height = 6)

        rank_waist_results <- list(
            K_A = K_A, K_B = K_B, K_C = K_C,
            wb_A_theory = wb_A_theory, wb_B_theory = wb_B_theory,
            wb_C_theory = wb_C_theory,
            wb_A_emp = wb_A_emp, wb_B_emp = wb_B_emp, wb_C_emp = wb_C_emp,
            n_in_A = n_in_A, n_in_B = n_in_B, n_in_C = n_in_C,
            n_recov_A = n_recov_A, n_recov_B = n_recov_B,
            n_recov_C = n_recov_C,
            B_to_thresh_gap = B_to_thresh_gap,
            wb_comparison = wb_comparison
        )


        # ======================================================================
        # SECTION 2: Si_norm GMM decomposition
        # ======================================================================

        cat("\n  -- Section 2: Si_norm GMM decomposition --\n")

        fit <- Mclust(si_norm, G = GMM_G_RANGE, verbose = FALSE)
        params <- extract_gmm_params(fit)

        cat(sprintf("    BIC-optimal: %d components (BIC=%.1f)\n",
                    fit$G, max(fit$BIC, na.rm = TRUE)))
        for (i in 1:nrow(params))
            cat(sprintf("      C%d: mean=%.4f  sd=%.4f  prop=%.1f%%\n",
                        i, params$mean[i], params$sd[i],
                        params$prop[i] * 100))

        adj_ints <- gmm_adjacent_intersections(params)
        if (nrow(adj_ints) > 0)
            cat(sprintf("    Adjacent intersections: %s\n",
                        paste(round(adj_ints$intersection, 4),
                              collapse = ", ")))

        waist_idx  <- which.min(params$mean)
        waist_mean <- params$mean[waist_idx]
        waist_sd   <- params$sd[waist_idx]
        waist_prop <- params$prop[waist_idx]

        cat(sprintf(
            "    Waist spike: C%d (mean=%.4f, sd=%.4f, prop=%.1f%%)\n",
            params$comp[waist_idx], waist_mean, waist_sd,
            waist_prop * 100))

        # GMM plot
        df_gmm <- data.frame(value = si_norm)
        x_grid <- seq(min(si_norm) - 0.1, max(si_norm) + 0.1,
                      length.out = 500)

        comp_df <- data.frame()
        mix_dens <- rep(0, length(x_grid))
        for (i in 1:nrow(params)) {
            d <- params$prop[i] * dnorm(x_grid, params$mean[i],
                                        params$sd[i])
            mix_dens <- mix_dens + d
            comp_df <- rbind(comp_df, data.frame(
                x = x_grid, density = d,
                Component = sprintf("C%d (mean=%.3f, sd=%.3f, %.0f%%)",
                                    i, params$mean[i], params$sd[i],
                                    params$prop[i] * 100)
            ))
        }
        mix_df <- data.frame(x = x_grid, density = mix_dens)

        p_gmm <- ggplot() +
            geom_histogram(data = df_gmm,
                           aes(x = value, y = after_stat(density)),
                           bins = 50, fill = "#9B8EC4", colour = "white",
                           alpha = 0.5) +
            geom_line(data = mix_df, aes(x = x, y = density),
                      linewidth = 1, colour = "black") +
            geom_line(data = comp_df,
                      aes(x = x, y = density, colour = Component),
                      linewidth = 0.7, linetype = "dashed") +
            geom_vline(xintercept = SKCS_SI_THRESH, linetype = "dashed",
                       colour = "#6B5B8A", linewidth = 0.5) +
            annotate("text", x = SKCS_SI_THRESH + 0.01, y = Inf, vjust = 2,
                     label = "SKCS = 0.05", colour = "#6B5B8A", size = 3,
                     hjust = 0) +
            labs(
                title = sprintf(
                    "%s - Si_norm GMM decomposition (%d components)",
                    ds, fit$G),
                subtitle = sprintf(
                    paste0("Waist spike C%d: mean=%.4f, sd=%.4f, ",
                           "%.1f%% of cells | BIC=%.1f\n",
                           "The waist is a distinct narrow component ",
                           "near Si_norm = 0, capturing cells ",
                           "structurally pinned by GOI-median proximity"),
                    params$comp[waist_idx], waist_mean, waist_sd,
                    waist_prop * 100,
                    max(fit$BIC, na.rm = TRUE)),
                x = "Si_norm", y = "Density"
            ) +
            theme_bw(base_size = 12) +
            theme(plot.title = element_text(face = "bold"),
                  plot.subtitle = element_text(size = 9,
                                               colour = "grey30"),
                  legend.position = "right")

        if (nrow(adj_ints) > 0)
            for (j in 1:nrow(adj_ints))
                p_gmm <- p_gmm +
                    geom_vline(xintercept = adj_ints$intersection[j],
                               linetype = "dotted", colour = "#7A6840",
                               linewidth = 0.5)

        print(p_gmm)
        save_plot(sprintf("%s_%s_gmm_Si_norm.png", OUT_PREFIX, ds), p_gmm)


        # ======================================================================
        # SECTION 3: GMM posterior P(waist | Si_norm)
        # ======================================================================

        cat("\n  -- Section 3: GMM posterior P(waist | Si_norm) --\n")

        cat("    Posterior decay:\n")
        for (x in c(0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.10)) {
            pw <- gmm_posterior(x, params$mean, params$sd,
                                params$prop)[waist_idx]
            marker <- ifelse(x == 0.05, " <- SKCS", "")
            cat(sprintf("      Si_norm = %.2f: P(waist) = %.4f%%%s\n",
                        x, pw * 100, marker))
        }

        p_waist_at_005 <- gmm_posterior(
            0.05, params$mean, params$sd, params$prop)[waist_idx]

        w_pctls <- qnorm(c(0.95, 0.99, 0.999), waist_mean, waist_sd)
        cat(sprintf(
            "    Waist percentiles: 95th=%.4f  99th=%.4f  99.9th=%.4f\n",
            w_pctls[1], w_pctls[2], w_pctls[3]))

        p05_bnd <- find_posterior_boundary(
            0.05, waist_idx, params$mean, params$sd, params$prop)
        p01_bnd <- find_posterior_boundary(
            0.01, waist_idx, params$mean, params$sd, params$prop)

        cat(sprintf("    P(waist) = 5%%  at Si_norm = %.4f\n", p05_bnd))
        cat(sprintf("    P(waist) = 1%%  at Si_norm = %.4f\n", p01_bnd))
        cat(sprintf("    >>> P(waist | 0.05) = %.4f%% <<<\n",
                    p_waist_at_005 * 100))

        waist_grey_int <- if (nrow(adj_ints) > 0) {
            adj_ints$intersection[1]
        } else {
            NA_real_
        }

        # Reconcile waist boundary with GMM
        waist_bnd_raw <- waist_bnd
        if (!is.na(waist_grey_int) && waist_bnd > waist_grey_int) {
            waist_bnd <- waist_grey_int
            waist_bnd_method <- sprintf(
                "capped at GMM waist-grey (density trough was %.4f)",
                waist_bnd_raw)
            cat(sprintf("\n    Waist boundary reconciliation:\n"))
            cat(sprintf("      Density trough: %.4f\n", waist_bnd_raw))
            cat(sprintf("      GMM waist-grey: %.4f (lower, adopted)\n",
                        waist_grey_int))
            cat(sprintf(
                "      >>> Reconciled waist boundary: |Si_norm| < %.4f <<<\n",
                waist_bnd))
            n_waist_emp <- sum(abs(si_norm) < waist_bnd)
            cat(sprintf("      Waist cells: %d (%.1f%% of total)\n",
                        n_waist_emp, n_waist_emp / n_cells * 100))
        } else {
            cat(sprintf(
                "\n    Waist boundary: %.4f (%s) - no GMM cap needed\n",
                waist_bnd, waist_bnd_method))
        }

        # Store posterior curve for cross-dataset plot
        x_seq <- seq(-0.02, 0.12, length.out = 500)
        pw_curve <- sapply(x_seq, function(x)
            gmm_posterior(x, params$mean, params$sd,
                          params$prop)[waist_idx])
        all_posterior_df <- rbind(all_posterior_df,
            data.frame(Si_norm = x_seq, P_waist = pw_curve, Dataset = ds))


        # ======================================================================
        # SECTION 4: Rank-based verification (tests)
        # ======================================================================

        cat("\n  -- Section 4: Rank-based verification --\n")

        goi_rank_pctile <- rank(score_df$GOI_expr,
                                ties.method = "average") / n_cells
        goi_dev <- abs(goi_rank_pctile - 0.5)

        # Bin definitions
        waist_label <- sprintf("Waist: |Si|<%.3f", waist_bnd)
        trans_label <- sprintf("[%.3f, 0.04)", waist_bnd)
        bin_labels <- c(waist_label, trans_label,
                        "[0.04, 0.06)",
                        "[0.08, 0.12)", "[0.15, 0.25)", "[0.30, 0.50)")
        bin_defs <- list()
        bin_defs[[bin_labels[1]]] <- abs(si_norm) < waist_bnd
        bin_defs[[bin_labels[2]]] <- si_norm >= waist_bnd & si_norm < 0.04
        bin_defs[[bin_labels[3]]] <- si_norm >= 0.04 & si_norm < 0.06
        bin_defs[[bin_labels[4]]] <- si_norm >= 0.08 & si_norm < 0.12
        bin_defs[[bin_labels[5]]] <- si_norm >= 0.15 & si_norm < 0.25
        bin_defs[[bin_labels[6]]] <- si_norm >= 0.30 & si_norm < 0.50

        thresh_label <- bin_labels[3]   # "[0.04, 0.06)"

        cat("    GOI rank deviation by Si_norm bin:\n")
        cat(sprintf("    %-22s %5s %8s %8s\n", "Bin", "n", "median", "IQR"))
        for (bname in bin_labels) {
            idx <- bin_defs[[bname]]
            n_bin <- sum(idx)
            if (n_bin < 3) next
            cat(sprintf("    %-22s %5d %8.4f %8.4f\n",
                        bname, n_bin,
                        median(goi_dev[idx]), IQR(goi_dev[idx])))
        }

        waist_dev  <- goi_dev[bin_defs[[waist_label]]]
        thresh_dev <- goi_dev[bin_defs[[thresh_label]]]
        med_waist  <- median(waist_dev)
        med_thresh <- median(thresh_dev)
        med_ratio  <- med_thresh / max(med_waist, 1e-6)

        # Violin data
        df_rv <- data.frame()
        for (bname in bin_labels) {
            idx <- bin_defs[[bname]]
            if (sum(idx) < 3) next
            df_rv <- rbind(df_rv, data.frame(
                GOI_dev = goi_dev[idx], Si_bin = bname,
                stringsAsFactors = FALSE))
        }
        df_rv$Si_bin <- factor(df_rv$Si_bin, levels = bin_labels)

        # Shapiro-Wilk per bin
        cat("\n    Normality (Shapiro-Wilk per bin):\n")
        all_normal <- TRUE
        for (bname in bin_labels) {
            bin_vals <- goi_dev[bin_defs[[bname]]]
            if (length(bin_vals) < 3) next
            sw <- safe_shapiro(bin_vals)
            is_normal <- !is.na(sw$p.value) && sw$p.value >= 0.05
            if (!is_normal) all_normal <- FALSE
            cat(sprintf("      %-22s W=%.4f, p=%.2e -> %s\n",
                        bname, sw$W, sw$p.value,
                        ifelse(is_normal, "NORMAL", "NON-NORMAL")))
        }

        if (all_normal) {
            cat("    All bins normal; Kruskal-Wallis + Dunn's used\n")
        } else {
            cat("    At least one bin non-normal -> ",
                "Kruskal-Wallis + Dunn's post hoc\n")
        }

        # Kruskal-Wallis omnibus
        kw <- kruskal.test(GOI_dev ~ Si_bin, data = df_rv)
        cat(sprintf("\n    Kruskal-Wallis: chi2 = %.2f, df = %d, p = %.2e\n",
                    kw$statistic, kw$parameter, kw$p.value))

        pval_collector <- rbind(pval_collector, data.frame(
            Test = sprintf("%s_rankverif_KW_omnibus", ds),
            Method = "Kruskal-Wallis",
            p_value = kw$p.value, stringsAsFactors = FALSE))

        # Dunn's post hoc with BH correction
        cat("    Dunn's post hoc (BH-corrected):\n")
        dunn_res <- dunn.test(df_rv$GOI_dev, df_rv$Si_bin,
                               method = "bh", alpha = 0.05, altp = TRUE)

        dunn_comp_labels <- dunn_res$comparisons
        waist_thresh_pattern <- paste0(waist_label, " - ", thresh_label)
        thresh_waist_pattern <- paste0(thresh_label, " - ", waist_label)

        dunn_match <- which(dunn_comp_labels == waist_thresh_pattern |
                             dunn_comp_labels == thresh_waist_pattern)

        if (length(dunn_match) > 0) {
            dunn_z     <- dunn_res$Z[dunn_match[1]]
            dunn_p_adj <- dunn_res$altP.adjusted[dunn_match[1]]
            cat(sprintf("\n    Key comparison (Waist vs %s):\n",
                        thresh_label))
            cat(sprintf("      z = %.4f, BH-adjusted p = %.2e\n",
                        dunn_z, dunn_p_adj))
            cat(sprintf("      Waist:  median = %.4f (n=%d)\n",
                        med_waist, length(waist_dev)))
            cat(sprintf("      %s: median = %.4f (n=%d)\n", thresh_label,
                        med_thresh, length(thresh_dev)))
            cat(sprintf("      Ratio: %.1fx more dispersed at threshold\n",
                        med_ratio))
        } else {
            dunn_z <- NA_real_; dunn_p_adj <- NA_real_
            cat("    Warning: waist vs threshold comparison ",
                "not found in Dunn's output\n")
        }

        # Kendall tau
        set.seed(42)
        tau_idx  <- if (n_cells > 5000) sample(n_cells, 5000) else 1:n_cells
        tau_test <- cor.test(si_norm[tau_idx], goi_dev[tau_idx],
                             method = "kendall")
        cat(sprintf("    Kendall tau: %.4f (p = %.2e)\n",
                    tau_test$estimate, tau_test$p.value))

        pval_collector <- rbind(pval_collector, data.frame(
            Test = sprintf("%s_rankverif_kendall_tau", ds),
            Method = "Kendall",
            p_value = tau_test$p.value, stringsAsFactors = FALSE))

        # BH FDR correction
        pval_collector$FDR <- p.adjust(pval_collector$p_value, method = "BH")
        pval_collector$Sig <- ifelse(pval_collector$FDR < 0.05, "*", "")

        fdr_kw  <- pval_collector$FDR[pval_collector$Test ==
                                        sprintf("%s_rankverif_KW_omnibus", ds)]
        fdr_tau <- pval_collector$FDR[pval_collector$Test ==
                                        sprintf("%s_rankverif_kendall_tau", ds)]

        cat(sprintf(
            "    FDR-corrected: KW FDR = %.2e | Kendall FDR = %.2e\n",
            fdr_kw, fdr_tau))


        # ======================================================================
        # SECTION 5: Synthesis
        # ======================================================================

        cat("\n  -- Section 5: Synthesis --\n")

        synth_ds <- data.frame(
            Dataset  = ds,
            Approach = c("SKCS (adopted)",
                         "GMM waist-grey intersection",
                         "GMM P(waist) < 5%",
                         "GMM P(waist) < 1%"),
            Si_thresh = c(
                SKCS_SI_THRESH,
                waist_grey_int,
                p05_bnd,
                p01_bnd
            ),
            Note = c(
                "Adopted threshold",
                "GMM component boundary",
                sprintf("P(waist|0.05) = %.2f%%", p_waist_at_005 * 100),
                sprintf("P(waist|0.05) = %.2f%%", p_waist_at_005 * 100)
            ),
            stringsAsFactors = FALSE
        )

        cat("\n    Threshold comparison:\n")
        print(synth_ds, row.names = FALSE)
        cat(sprintf("\n    Key result: P(waist | Si_norm = 0.05) = %.4f%%\n",
                    p_waist_at_005 * 100))

        all_synth_rows[[ds]] <- synth_ds

        # Synthesis plot
        df_synth <- data.frame(Si_norm = si_norm)

        p_synth <- ggplot(df_synth, aes(x = Si_norm)) +
            geom_histogram(aes(y = after_stat(density)), bins = 50,
                           fill = "#9B8EC4", colour = "white", alpha = 0.7) +
            geom_density(linewidth = 0.8)

        if (!is.na(waist_grey_int))
            p_synth <- p_synth +
                geom_vline(xintercept = waist_grey_int, linetype = "dashed",
                           colour = "darkorange3", linewidth = 0.7) +
                annotate("text", x = waist_grey_int, y = Inf, vjust = -0.5,
                         hjust = -0.05,
                         label = sprintf("GMM waist-grey = %.4f",
                                         waist_grey_int),
                         colour = "darkorange3", size = 3.2,
                         fontface = "bold")

        if (!is.na(p05_bnd))
            p_synth <- p_synth +
                geom_vline(xintercept = p05_bnd, linetype = "dashed",
                           colour = "darkcyan", linewidth = 0.7) +
                annotate("text", x = p05_bnd, y = Inf, vjust = -2.0,
                         hjust = -0.05,
                         label = sprintf("P(waist) < 5%% = %.4f", p05_bnd),
                         colour = "darkcyan", size = 3.2, fontface = "bold")

        p_synth <- p_synth +
            geom_vline(xintercept = SKCS_SI_THRESH, linetype = "dashed",
                       colour = "#6B5B8A", linewidth = 0.7) +
            annotate("text", x = SKCS_SI_THRESH, y = Inf, vjust = -3.5,
                     hjust = -0.05,
                     label = sprintf("SKCS = %.2f", SKCS_SI_THRESH),
                     colour = "#6B5B8A", size = 3.2, fontface = "bold")

        p_synth <- p_synth +
            labs(
                title = sprintf("%s - Si_norm threshold synthesis", ds),
                subtitle = sprintf(
                    paste0("GMM: P(waist|0.05) = %.2f%% | Rank: %.1fx ",
                           "dispersion, KW FDR = %.1e, Dunn's p.adj = %.1e\n",
                           "Two independent lines of evidence support ",
                           "Si_norm > 0.05 as a conservative calibration ",
                           "floor"),
                    p_waist_at_005 * 100, med_ratio,
                    fdr_kw, dunn_p_adj),
                x = "Si_norm", y = "Density"
            ) +
            theme_bw(base_size = 12) +
            theme(plot.title = element_text(face = "bold"),
                  plot.subtitle = element_text(size = 9, colour = "grey30"))

        print(p_synth)
        save_plot(sprintf("%s_%s_synthesis.png", OUT_PREFIX, ds), p_synth)


        # Accumulate p-values
        all_pvals <- rbind(all_pvals, pval_collector)

        # Store results
        all_results[[ds]] <- list(
            n_cells          = n_cells,
            waist_bnd        = waist_bnd,
            waist_bnd_method = waist_bnd_method,
            gmm_fit          = fit,
            gmm_params       = params,
            gmm_ints         = adj_ints,
            waist_idx        = waist_idx,
            waist_grey_int   = waist_grey_int,
            p_waist_at_005   = p_waist_at_005,
            p05_boundary     = p05_bnd,
            p01_boundary     = p01_bnd,
            waist_pctls      = w_pctls,
            rank_verification = list(
                med_waist    = med_waist,
                med_thresh   = med_thresh,
                med_ratio    = med_ratio,
                kw_chi2      = kw$statistic,
                kw_df        = kw$parameter,
                kw_p         = kw$p.value,
                kw_fdr       = fdr_kw,
                dunn_z       = dunn_z,
                dunn_p_adj   = dunn_p_adj,
                tau          = tau_test$estimate,
                tau_p        = tau_test$p.value,
                tau_fdr      = fdr_tau
            ),
            pval_table       = pval_collector,
            synth_table      = synth_ds,
            rank_waist       = rank_waist_results
        )

    }  # end dataset loop


    # ==========================================================================
    # After the loop
    # ==========================================================================

    # ======================================================================
    # SECTION 6: Posterior profile plot (both datasets)
    # ======================================================================

    cat("\n\n  -- Section 6: Posterior profile plot --\n")

    p_post <- ggplot(all_posterior_df,
                     aes(x = Si_norm, y = P_waist, colour = Dataset)) +
        geom_line(linewidth = 1) +

        geom_hline(yintercept = 0.05, linetype = "dotted",
                   colour = "grey50") +
        annotate("text", x = 0.115, y = 0.055, label = "5%",
                 colour = "grey50", size = 3, hjust = 1) +
        geom_hline(yintercept = 0.01, linetype = "dotted",
                   colour = "grey50") +
        annotate("text", x = 0.115, y = 0.015, label = "1%",
                 colour = "grey50", size = 3, hjust = 1) +

        geom_vline(xintercept = SKCS_SI_THRESH, linetype = "dashed",
                   colour = "#6B5B8A", linewidth = 0.5, alpha = 0.7) +
        annotate("text", x = SKCS_SI_THRESH + 0.002, y = 0.50,
                 label = sprintf("SKCS = %.2f", SKCS_SI_THRESH),
                 colour = "#6B5B8A", size = 3, hjust = 0) +

        scale_colour_manual(
            values = c("NonNeoC" = "#2166AC", "NeoC" = "#A0522D")) +
        labs(
            title = "GMM posterior probability of waist-spike membership",
            subtitle = sprintf(
                paste0("P(waist|0.05): NonNeoC = %.2f%%, ",
                       "NeoC = %.4f%% | ",
                       "P < 5%% at Si = %.4f / %.4f | ",
                       "P < 1%% at Si = %.4f / %.4f\n",
                       "At the SKCS threshold, waist contamination ",
                       "probability is negligible in both datasets"),
                all_results$NonNeoC$p_waist_at_005 * 100,
                all_results$NeoC$p_waist_at_005 * 100,
                all_results$NonNeoC$p05_boundary,
                all_results$NeoC$p05_boundary,
                all_results$NonNeoC$p01_boundary,
                all_results$NeoC$p01_boundary),
            x = "Si_norm",
            y = "P(waist component)"
        ) +
        coord_cartesian(xlim = c(-0.02, 0.12), ylim = c(0, 0.70)) +
        theme_bw(base_size = 12) +
        theme(
            plot.title = element_text(face = "bold"),
            plot.subtitle = element_text(size = 9, colour = "grey30"),
            legend.position = c(0.85, 0.85),
            legend.background = element_rect(colour = "grey80")
        )

    print(p_post)
    save_plot(sprintf("%s_posterior_profile.png", OUT_PREFIX), p_post)


    # ======================================================================
    # SECTION 8: Summary + Export
    # ======================================================================

    cat("\n\n", strrep("=", 72), "\n")
    cat("  SUMMARY\n")
    cat(strrep("=", 72), "\n\n")

    cat("  Waist boundary (density gap):\n")
    for (d in DATASETS) {
        r <- all_results[[d]]
        cat(sprintf("     %s: |Si_norm| < %.4f (%s)\n",
                    d, r$waist_bnd, r$waist_bnd_method))
    }

    cat("\n  GMM waist spike:\n")
    for (d in DATASETS) {
        r <- all_results[[d]]
        cat(sprintf(
            "     %s: waist spike C%d (mean=%.4f, sd=%.4f, %.1f%% of cells)\n",
            d, r$gmm_params$comp[r$waist_idx],
            r$gmm_params$mean[r$waist_idx],
            r$gmm_params$sd[r$waist_idx],
            r$gmm_params$prop[r$waist_idx] * 100))
    }

    cat("\n  Rank verification at 0.05:\n")
    for (d in DATASETS) {
        r <- all_results[[d]]$rank_verification
        cat(sprintf(
            "     %s: %.1fx dispersion ratio | KW FDR = %.1e, Dunn's p.adj = %.1e | tau = %.3f\n",
            d, r$med_ratio, r$kw_fdr, r$dunn_p_adj, r$tau))
    }

    cat("\n  GMM posterior at threshold:\n")
    for (d in DATASETS)
        cat(sprintf("     %s: P(waist | Si_norm = 0.05) = %.4f%%\n",
                    d, all_results[[d]]$p_waist_at_005 * 100))


    # == Export ==
    cat("\n", strrep("=", 72), "\n")
    cat("  Export\n")
    cat(strrep("=", 72), "\n")

    synth_all <- do.call(rbind, all_synth_rows)

    wb <- createWorkbook()
    hs <- createStyle(textDecoration = "bold", fgFill = "#DCE6F1",
                      border = "Bottom")

    # Sheet 1: Waist boundary
    waist_summary <- data.frame(
        Dataset = DATASETS,
        Waist_boundary = sapply(DATASETS, function(d)
            all_results[[d]]$waist_bnd),
        Method = sapply(DATASETS, function(d)
            all_results[[d]]$waist_bnd_method),
        N_waist_cells = sapply(DATASETS, function(d) {
            sum(abs(p9b_data[[d]]$score_pooled_df[[COL_SI]][
                !p9b_data[[d]]$score_pooled_df$is_GOI_neg]) <
                all_results[[d]]$waist_bnd)
        }),
        Pct_of_total = sapply(DATASETS, function(d) {
            si <- p9b_data[[d]]$score_pooled_df[[COL_SI]][
                !p9b_data[[d]]$score_pooled_df$is_GOI_neg]
            sum(abs(si) < all_results[[d]]$waist_bnd) / length(si) * 100
        })
    )
    addWorksheet(wb, "Waist_Boundary")
    writeData(wb, "Waist_Boundary", waist_summary)
    addStyle(wb, "Waist_Boundary", hs, rows = 1,
             cols = 1:ncol(waist_summary))
    setColWidths(wb, "Waist_Boundary", cols = 1:ncol(waist_summary),
                 widths = "auto")

    # Sheet 1b: Rank-based waist boundary comparison
    rank_wb_all <- do.call(rbind, lapply(DATASETS, function(d) {
        rw <- all_results[[d]]$rank_waist$wb_comparison
        rw$Dataset <- d
        rw[, c("Dataset", names(rw)[names(rw) != "Dataset"])]
    }))
    addWorksheet(wb, "Waist_RankBased")
    writeData(wb, "Waist_RankBased", rank_wb_all)
    addStyle(wb, "Waist_RankBased", hs, rows = 1,
             cols = 1:ncol(rank_wb_all))
    setColWidths(wb, "Waist_RankBased", cols = 1:ncol(rank_wb_all),
                 widths = "auto")

    # Sheet 2: Threshold comparison
    addWorksheet(wb, "Thresholds")
    writeData(wb, "Thresholds", synth_all)
    addStyle(wb, "Thresholds", hs, rows = 1, cols = 1:ncol(synth_all))
    setColWidths(wb, "Thresholds", cols = 1:ncol(synth_all),
                 widths = "auto")

    # Sheet 3: GMM posterior summary
    post_summary <- data.frame(
        Dataset = DATASETS,
        Waist_mean = sapply(DATASETS, function(d)
            all_results[[d]]$gmm_params$mean[all_results[[d]]$waist_idx]),
        Waist_sd = sapply(DATASETS, function(d)
            all_results[[d]]$gmm_params$sd[all_results[[d]]$waist_idx]),
        Waist_prop_pct = sapply(DATASETS, function(d)
            all_results[[d]]$gmm_params$prop[
                all_results[[d]]$waist_idx] * 100),
        Waist_grey_GMM = sapply(DATASETS, function(d)
            all_results[[d]]$waist_grey_int),
        P_waist_at_005_pct = sapply(DATASETS, function(d)
            all_results[[d]]$p_waist_at_005 * 100),
        P05_boundary = sapply(DATASETS, function(d)
            all_results[[d]]$p05_boundary),
        P01_boundary = sapply(DATASETS, function(d)
            all_results[[d]]$p01_boundary)
    )
    addWorksheet(wb, "GMM_Posterior")
    writeData(wb, "GMM_Posterior", post_summary)
    addStyle(wb, "GMM_Posterior", hs, rows = 1,
             cols = 1:ncol(post_summary))
    setColWidths(wb, "GMM_Posterior", cols = 1:ncol(post_summary),
                 widths = "auto")

    # Sheet 4: Rank verification summary
    rank_summary <- data.frame(
        Dataset = DATASETS,
        Waist_median_GOI_dev = sapply(DATASETS, function(d)
            all_results[[d]]$rank_verification$med_waist),
        Thresh_median_GOI_dev = sapply(DATASETS, function(d)
            all_results[[d]]$rank_verification$med_thresh),
        Dispersion_ratio = sapply(DATASETS, function(d)
            all_results[[d]]$rank_verification$med_ratio),
        KW_chi2 = sapply(DATASETS, function(d)
            all_results[[d]]$rank_verification$kw_chi2),
        KW_p = sapply(DATASETS, function(d)
            all_results[[d]]$rank_verification$kw_p),
        KW_FDR = sapply(DATASETS, function(d)
            all_results[[d]]$rank_verification$kw_fdr),
        Dunn_z = sapply(DATASETS, function(d)
            all_results[[d]]$rank_verification$dunn_z),
        Dunn_p_adj = sapply(DATASETS, function(d)
            all_results[[d]]$rank_verification$dunn_p_adj),
        Kendall_tau = sapply(DATASETS, function(d)
            all_results[[d]]$rank_verification$tau),
        Kendall_p = sapply(DATASETS, function(d)
            all_results[[d]]$rank_verification$tau_p),
        Kendall_FDR = sapply(DATASETS, function(d)
            all_results[[d]]$rank_verification$tau_fdr)
    )
    addWorksheet(wb, "Rank_Verification")
    writeData(wb, "Rank_Verification", rank_summary)
    addStyle(wb, "Rank_Verification", hs, rows = 1,
             cols = 1:ncol(rank_summary))
    setColWidths(wb, "Rank_Verification", cols = 1:ncol(rank_summary),
                 widths = "auto")

    # Sheet 5: All p-values with FDR
    addWorksheet(wb, "P_values_FDR")
    writeData(wb, "P_values_FDR", all_pvals)
    addStyle(wb, "P_values_FDR", hs, rows = 1, cols = 1:ncol(all_pvals))
    setColWidths(wb, "P_values_FDR", cols = 1:ncol(all_pvals),
                 widths = "auto")

    xlsx_path <- file.path(PART9B_DIR,
                           paste0(OUT_PREFIX,
                                  "_SiKi_threshold_comparison.xlsx"))
    saveWorkbook(wb, xlsx_path, overwrite = TRUE)
    cat(sprintf("  Excel: %s\n", basename(xlsx_path)))

    rds_path <- file.path(PART9B_DIR,
                          paste0(OUT_PREFIX, "_SiKi_results.rds"))
    saveRDS(list(all_results  = all_results,
                 arm_tag      = ARM_TAG,
                 norm_tag     = NORM_TAG,
                 goi          = GOI,
                 dfg_genes    = DFG,
                 si_thresh    = SKCS_SI_THRESH,
                 ki_thresh    = SKCS_KI_THRESH),
            rds_path)
    cat(sprintf("  RDS: %s\n", basename(rds_path)))

    cat("\n  Files:\n")
    flist <- list.files(PART9B_DIR, full.names = TRUE)
    for (f in flist) {
        sz <- file.info(f)$size
        cat(sprintf("    %-55s %5.0f KB\n", basename(f), sz / 1024))
    }

    cat("\n", strrep("=", 72), "\n")
    cat("  Part9b complete (", ARM_TAG, ")\n")
    cat(strrep("=", 72), "\n")

    while (sink.number() > 0) sink()

    invisible(list(pass = 2L, warn = 0L, fail = 0L))
}


# =============================================================================
# SECTION 0e: Driver
# =============================================================================
# Runs only the primary arm (Mic13_PVM3).

CEP_ARM <- "Mic13_PVM3"
cep_run_arms(PART_TAG, run_one_arm)
