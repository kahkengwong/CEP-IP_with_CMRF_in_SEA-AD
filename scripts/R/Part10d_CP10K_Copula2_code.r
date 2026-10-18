# =============================================================================
# Part10d_CP10K_Copula2_code.r
# SEA-AD CEP-IP/CMRF, Part 10d, the IP and gates in copula coordinates
#   (Figure 4, Supplementary Figures S8 to S10, Table S9)
#   For each arm: places the inflection point at (IP_u, IP_v) in copula space
#   and tests the copula density at IP_u and IP_v against the flanks, fits the
#   GAM image curve g(u) on TREP cells, draws the synthesis figure, and maps
#   the Ki-Si pair and the SKCS gates into copula coordinates. All statistics
#   are cell level (no donor random intercept or donor weights).
#
# Usage (from the SEA-AD folder, or set SEAAD_ROOT), after Part 10,
#   Part10b_CP10K_IPRS_code.r and Part10b_CP10K_SKCS_TREP_code.r
#   source("Part10d_CP10K_Copula2_code.r")                         # all arms
#   CEP_ARM <- "PVM1_all"; source("Part10d_CP10K_Copula2_code.r")  # one arm
#   Needs CEP_arms.R, R_config.R and R_config_HKG8.R in the same folder.
#
# Outputs (Part10d_Copula2/Part10d_Copula2_<arm>/): figures (PDF, PNG),
#   copula workbook per dataset, results RDS, guardrail ledger and logs.
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
library(viridis)       # turbo palette for the K_i gradient

set.seed(42)

# Install ggnewscale (second fill scale in Step 5) if missing.
if (!requireNamespace("ggnewscale", quietly = TRUE)) {
    cat("  Installing 'ggnewscale' for the dual fill scales in Step 5 ...\n")
    tryCatch(
        install.packages("ggnewscale", quiet = TRUE,
                         repos = "https://cloud.r-project.org"),
        error   = function(e)
            cat("  ggnewscale install failed:", conditionMessage(e), "\n"),
        warning = function(w)
            cat("  ggnewscale install warning:", conditionMessage(w), "\n"))
}

# CEP_arms.R holds the arm table, the folder layout and the driver.
if (!exists("CEP_HOME")) CEP_HOME <- Sys.getenv("SEAAD_ROOT", getwd())
.p10d_harness <- file.path(CEP_HOME, "CEP_arms.R")
if (!file.exists(.p10d_harness))
    stop("CEP_arms.R not found at: ", .p10d_harness,
         "\n  It holds the arm table shared by Parts 9 onwards.",
         call. = FALSE)
source(.p10d_harness)

# PART_DIR_TAG builds the folder name, PART_TAG the file names.
PART_TAG     <- "Part10d"           # file-name tag
PART_DIR_TAG <- "Part10d_Copula2"   # results-folder tag

# Phase-2 cap for the PRSS loop (at most 8 + 5 iterations).
PRSS_PHASE2_MAX_P10D <- 5L


# =============================================================================
# SECTION 0d: Pure helpers
# =============================================================================


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

# Writes one figure as PDF and PNG; a device error skips that figure only.
save_pp <- function(p, base, pdf_w = 10, pdf_h = 9,
                    png_w = 2200, png_h = 2000, res = 200) {
    .draw <- if (inherits(p, "grob") || inherits(p, "gtable") ||
                 inherits(p, "arrangelist"))
        function() { grid.newpage(); grid.draw(p) } else function() print(p)
    .dev_entry <- grDevices::dev.cur()
    .unwind <- function() {
        while (grDevices::dev.cur() != .dev_entry &&
               grDevices::dev.cur() != 1L)
            tryCatch(grDevices::dev.off(),
                     error = function(e2) return(invisible(NULL)))
    }
    tryCatch({
        pdf(paste0(base, ".pdf"), width = pdf_w, height = pdf_h)
        .draw(); dev.off()
    }, error = function(e) {
        .unwind()
        cat("  [save_pp] PDF error for", basename(base), ":",
            conditionMessage(e), "\n")
    })
    tryCatch({
        png(paste0(base, ".png"), width = png_w, height = png_h, res = res)
        .draw(); dev.off()
    }, error = function(e) {
        .unwind()
        cat("  [save_pp] PNG error for", basename(base), ":",
            conditionMessage(e), "\n")
    })
    invisible(NULL)
}

# C_hat_i = #{j: u_j <= u_i and v_j <= v_i} / n, computed in row blocks.
compute_ecdf_copula <- function(u, v, block_budget = 5e6) {
    n <- length(u)
    # For each i, count #{j: u_j <= u_i and v_j <= v_i} including j=i
    blk   <- max(1L, min(n, as.integer(block_budget %/% max(1L, n))))
    C_hat <- numeric(n)
    for (.s in seq(1L, n, by = blk)) {
        idx <- .s:min(.s + blk - 1L, n)
        C_hat[idx] <- rowSums(outer(u[idx], u, "<=") & outer(v[idx], v, "<="))
    }
    C_hat / n
}

# Kendall tau_b per cell from the strict-inequality copula counts:
# Ki = (Ci - Di) / sqrt((n-1-T_u)(n-1-T_v)), computed in row blocks.
compute_ordinal_copula_ki <- function(u, v, block_budget = 5e6) {
    n  <- length(u)
    # Rows per block, so one block matrix holds about `block_budget` logicals.
    blk <- max(1L, min(n, as.integer(block_budget %/% max(1L, n))))

    C_lo_ct <- numeric(n)   # counts, divided by n at the end
    C_hi_ct <- numeric(n)
    D_a_ct  <- numeric(n)
    D_b_ct  <- numeric(n)
    for (.s in seq(1L, n, by = blk)) {
        idx <- .s:min(.s + blk - 1L, n)
        u_lt <- outer(u[idx], u, "<"); u_gt <- outer(u[idx], u, ">")
        v_lt <- outer(v[idx], v, "<"); v_gt <- outer(v[idx], v, ">")
        C_lo_ct[idx] <- rowSums(u_lt & v_lt)
        C_hi_ct[idx] <- rowSums(u_gt & v_gt)
        D_a_ct[idx]  <- rowSums(u_lt & v_gt)
        D_b_ct[idx]  <- rowSums(u_gt & v_lt)
        rm(u_lt, u_gt, v_lt, v_gt)
    }
    Ci <- C_lo_ct + C_hi_ct
    Di <- D_a_ct  + D_b_ct

    # tau_b per-cell tie counts, excluding the cell itself.
    .iu <- match(u, unique(u)); .iv <- match(v, unique(v))
    T_u <- tabulate(.iu, nbins = max(.iu))[.iu] - 1
    T_v <- tabulate(.iv, nbins = max(.iv))[.iv] - 1

    denom_b <- sqrt((n - 1 - T_u) * (n - 1 - T_v))
    denom_b[denom_b == 0] <- 1  # guard
    gc(verbose = FALSE)
    list(C_strict_lower = C_lo_ct / n, C_strict_upper = C_hi_ct / n,
         Ki_formula = (Ci - Di) / denom_b)
}

# Ki from the empirical copula (tau_a identity, finite-n form):
# K_i = (n+1)/(n-1) + [4n/(n-1)]C_hat_i - [2n/(n-1)](u_i+v_i)
compute_ki_formula <- function(u, v, C_hat, exact_n_correction = TRUE) {
    n <- length(u)
    if (exact_n_correction) {
        # Exact finite-n formula
        (n + 1) / (n - 1) + (4 * n / (n - 1)) * C_hat - (2 * n / (n - 1)) * (u + v)
    } else {
        # Large-n approximation
        1 + 4 * C_hat - 2 * u - 2 * v
    }
}

# == PRSS helper: RSS + lambda * integral(f'') for one GAM ====================
calc_prss_copula <- function(gam_model) {
    RSS <- sum(residuals(gam_model)^2)
    # Guard: if no smooth or penalty matrix, return RSS only
    if (is.null(gam_model$smooth) || length(gam_model$smooth) == 0 ||
        is.null(gam_model$smooth[[1]]$S) || length(gam_model$smooth[[1]]$S) == 0 ||
        is.null(gam_model$sp) || length(gam_model$sp) == 0) {
        return(list(PRSS = RSS, RSS = RSS,
                    f_double_prime_integral = NA_real_, lambda = NA_real_))
    }
    lambda <- gam_model$sp[1]
    S      <- gam_model$smooth[[1]]$S[[1]]
    fp     <- gam_model$smooth[[1]]$first.para
    lp     <- gam_model$smooth[[1]]$last.para
    beta   <- coef(gam_model)[fp:lp]
    f_dpi  <- as.numeric(t(beta) %*% S %*% beta)
    list(PRSS = RSS + lambda * f_dpi,
         RSS  = RSS,
         f_double_prime_integral = f_dpi,
         lambda = lambda)
}

# == GAM image curve g(u): PRSS-based fitting =================================
# Fits v_trep ~ s(u_trep, bs = "tp") with REML inside a PRSS search over
# k: phase 1 explores k = 3..max_k, phase 2 refines around the best k.
fit_gam_image_prss <- function(u_trep, v_trep,
                                label      = "g(u)",
                                num_iters  = 100L,
                                seed_base  = 42L,
                                phase2_max = PRSS_PHASE2_MAX_P10D) {
    n_unique <- length(unique(u_trep))
    max_k    <- min(10L, n_unique - 1L)
    if (max_k < 3L) max_k <- 3L          # absolute floor

    best_prss    <- Inf
    best_model   <- NULL
    best_k       <- 3L
    best_iter    <- 1L
    n_iter_max   <- min(as.integer(num_iters), 8L + as.integer(phase2_max))
    prss_values  <- rep(NA_real_, n_iter_max)
    early_ctr    <- 0L
    early_stop   <- FALSE

    run_one <- function(i, k) {
        set.seed(seed_base + i)
        df_fit <- data.frame(u_trep = u_trep, v_resp = v_trep)
        gm <- tryCatch(
            mgcv::gam(v_resp ~ s(u_trep, bs = "tp", k = k),
                      data = df_fit, method = "REML",
                      select = TRUE, gamma = 1.5),
            error = function(e) {
                set.seed(seed_base + i)
                mgcv::gam(v_resp ~ s(u_trep, bs = "tp",
                                     k = max(3L, k - 1L)),
                          data = df_fit, method = "REML",
                          select = TRUE, gamma = 1.5)
            }
        )
        pc <- calc_prss_copula(gm)
        list(model = gm, prss = pc$PRSS, k = k)
    }

    # == Phase 1: exploration (iterations 1-8) ================================
    for (i in seq_len(min(8L, num_iters))) {
        k <- 3L + ((i - 1L) %% max(1L, max_k - 2L))
        k <- min(k, max_k)
        res <- tryCatch(run_one(i, k), error = function(e) NULL)
        if (!is.null(res)) {
            prss_values[i] <- res$prss
            cat(sprintf("    %s phase1 iter %d | k=%d | PRSS=%.5f\n",
                        label, i, k, res$prss))
            if (res$prss < best_prss) {
                best_prss <- res$prss; best_model <- res$model
                best_k <- k; best_iter <- i; early_ctr <- 0L
            } else { early_ctr <- early_ctr + 1L }
        } else {
            prss_values[i] <- NA_real_; early_ctr <- early_ctr + 1L
        }
    }

    # k_range seeded from phase 1 best
    k_range <- unique(pmax(3L, pmin(max_k,
                   if (best_k >= max_k)   c(max(3L, max_k - 1L), max_k) else
                   if (best_k <= 3L)      c(3L, min(4L, max_k)) else
                   c(best_k - 1L, best_k, min(best_k + 1L, max_k)))))
    cat(sprintf("    %s phase1 best: k=%d | PRSS=%.5f\n",
                label, best_k, best_prss))

    # == Phase 2: adaptive refinement, capped at phase2_max ===================
    for (i in 8L + seq_len(max(0L, n_iter_max - 8L))) {
        if (early_ctr >= 20L) {
            cat(sprintf("    %s early stop at iter %d\n", label, i))
            early_stop <- TRUE; break
        }
        k <- if (i %% 3L == 0L) sample(k_range, 1L) else best_k
        res <- tryCatch(run_one(i, k), error = function(e) NULL)
        if (!is.null(res)) {
            prss_values[i] <- res$prss
            cat(sprintf("    %s refine iter %d | k=%d | PRSS=%.5f\n",
                        label, i, k, res$prss))
            if (res$prss < best_prss) {
                best_prss <- res$prss; best_model <- res$model
                best_k <- k; best_iter <- i; early_ctr <- 0L
                k_range <- unique(pmax(3L, pmin(max_k,
                               if (best_k >= max_k) c(max(3L, max_k-1L), max_k) else
                               if (best_k <= 3L)    c(3L, min(4L, max_k)) else
                               c(best_k-1L, best_k, min(best_k+1L, max_k)))))
            } else { early_ctr <- early_ctr + 1L }
        } else {
            prss_values[i] <- NA_real_; early_ctr <- early_ctr + 1L
        }
    }

    actual_iters <- max(which(!is.na(prss_values)), na.rm = TRUE)
    if (!is.finite(actual_iters)) actual_iters <- n_iter_max

    list(best_model   = best_model,
         best_k       = best_k,
         best_iter    = best_iter,
         best_prss    = best_prss,
         phase2_max   = as.integer(phase2_max),
         n_iter_max   = n_iter_max,
         prss_values  = prss_values[seq_len(actual_iters)],
         actual_iters = actual_iters,
         early_stopped = early_stop,
         max_k        = max_k)
}

# Non-parametric tail dependence coefficients via empirical copula
tail_dep_np <- function(u, v, q = 0.10) {
    Cqq   <- mean(u <= q & v <= q)
    lam_L <- if (q > 0) Cqq / q else NA_real_
    Cuq   <- mean(u > (1 - q) & v > (1 - q))
    lam_U <- if (q > 0) Cuq / q else NA_real_
    list(lambda_L = lam_L, lambda_U = lam_U)
}

# Cramer-von Mises distance between two empirical copulas
cvm_copula_dist <- function(u1, v1, u2, v2, grid_size = 30) {
    pts <- seq(0, 1, length.out = grid_size)
    cvm <- 0
    for (uu in pts) {
        for (vv in pts) {
            C1 <- mean(u1 <= uu & v1 <= vv)
            C2 <- mean(u2 <= uu & v2 <= vv)
            cvm <- cvm + (C1 - C2)^2
        }
    }
    cvm / grid_size^2
}

# Four-group colour palette (as in Part 10c)
FOUR_COLS <- c(
    "Pre-IP TREP"      = "#E69F00",
    "Pre-IP non-TREP"  = "#1B3A8C",
    "Post-IP TREP"     = "#009E73",
    "Post-IP non-TREP" = "#D94F4F",
    "SKCS-core"        = "#4B0082",
    "Unknown"          = "gray60"
)

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

    # Output before the log exists is buffered and replayed into it.
    .p10d_preamble <- character(0)

    assign("CEP_ARM", arm_tag, envir = globalenv())
    on.exit(suppressWarnings(rm("CEP_ARM", envir = globalenv())), add = TRUE)

    .p10d_preamble <- tryCatch(
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

    .p10d_needed <- c("ARM_TAG", "ARM_DIR", "R_EXPORT_ROOT", "GOI", "DFG",
                      "HKG_GENES", "DFG_COMPOSITE_LABEL", "RUN_SUFFIX",
                      "PRIMARY_STYPE", "NORM_TAG", "NORM_TARGET_SUM",
                      "NORM_CEILING")
    .p10d_absent <- .p10d_needed[!vapply(.p10d_needed, exists, logical(1))]
    if (length(.p10d_absent) > 0)
        stop("The configuration loaded from\n    ", config_path,
             "\n  does not define: ", paste(.p10d_absent, collapse = ", "),
             call. = FALSE)

    if (!identical(ARM_TAG, arm_tag))
        stop("The configuration reports ARM_TAG '", ARM_TAG,
             "' but this arm is '", arm_tag, "'.", call. = FALSE)

    .p10d_preamble <- c(.p10d_preamble, utils::capture.output(
        cat("  Config OK: ", ARM_TAG, " | ", GOI, " + ", length(DFG), " DFG | ",
            NORM_TAG, "\n", sep = "")))

    set.seed(42)

    # == User parameters ======================================================
    DATASETS  <- c("NonNeoC", "NeoC")
    MIN_CELLS <- 50

    # == Input paths ==========================================================
    P10B_SKCS_PREFIX <- paste0(cep_out_prefix("Part10b", NORM_TAG,
                                              cep_file_tag_for(arm_tag)),
                               "_SKCS")
    p10d_ep_candidates <- function(ds) {
        fn <- paste0(P10B_SKCS_PREFIX, "_", ds, "_EP_TREP_values.csv")
        c(file.path(cep_arm_dir_for("Part10b_SKCS", arm_tag), fn),
          file.path(cep_arm_dir_for("Part10b", arm_tag), fn),
          file.path(R_EXPORT_ROOT, "Part10b_SKCS",
                    paste0("Part10b_scran_", ds, "_EP_TREP_values.csv")))
    }
    p10d_ep_csv <- function(ds) {
        cand <- p10d_ep_candidates(ds)
        hit  <- cand[file.exists(cand)]
        if (length(hit) > 0L) hit[1] else cand[1]
    }

    # == Output directory and file prefix =====================================
    OUT_DIR <- out_dir
    dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

    OUT_PREFIX  <- cep_out_prefix(PART_TAG, NORM_TAG,
                                  cep_file_tag_for(arm_tag))
    FILE_PREFIX <- paste0(OUT_PREFIX, "_Copula2")

    # == Log file =============================================================
    LOG_MESSAGES_P10D <- TRUE
    LOG_PATH_P10D <- file.path(OUT_DIR,
                               paste0(FILE_PREFIX, "_console_log_",
                                      format(Sys.time(), "%Y%m%d_%H%M%S"),
                                      ".txt"))
    MSG_PATH_P10D <- file.path(OUT_DIR,
                               paste0(FILE_PREFIX, "_console_log_messages.tmp"))

    while (sink.number() > 0) sink()
    # sink.number(type = "message") is 2 (stderr) when nothing is diverted.
    if (sink.number(type = "message") != 2L) sink(type = "message")

    .p10d_log_con <- file(LOG_PATH_P10D, open = "wt")
    .p10d_msg_con <- if (isTRUE(LOG_MESSAGES_P10D))
        file(MSG_PATH_P10D, open = "wt") else NULL
    sink(.p10d_log_con, append = FALSE, split = TRUE)
    if (!is.null(.p10d_msg_con)) sink(.p10d_msg_con, type = "message")
    on.exit({
        tryCatch({
            if (!isTRUE(.p10d_completed))
                cat("\n", strrep("!", 72),
                    "\n  Arm did not complete; the error is printed to the console.\n",
                    strrep("!", 72), "\n", sep = "")
            if (sink.number(type = "message") != 2L) sink(type = "message")
            while (sink.number() > 0) sink()
            if (!is.null(.p10d_msg_con)) close(.p10d_msg_con)
            close(.p10d_log_con)
            if (!is.null(.p10d_msg_con) && file.exists(MSG_PATH_P10D)) {
                .msg_txt <- readLines(MSG_PATH_P10D, warn = FALSE)
                if (length(.msg_txt) > 0)
                    cat(c("", strrep("=", 72),
                          "  MESSAGE STREAM (warnings, messages, errors)",
                          strrep("=", 72), .msg_txt, ""),
                        sep = "\n", file = LOG_PATH_P10D, append = TRUE)
                unlink(MSG_PATH_P10D)
            }
        }, error = function(e) invisible(NULL))
    }, add = TRUE)

    .p10d_completed <- FALSE

    .p10d_warn_opt <- getOption("warn")
    options(warn = 1)
    on.exit(options(warn = .p10d_warn_opt), add = TRUE)
    cat(sprintf("Part10d Copula2 - log started: %s\n\n", Sys.time()))

    if (length(.p10d_preamble) > 0)
        cat(paste(.p10d_preamble, collapse = "\n"), "\n\n", sep = "")

    # == DATASET_META =========================================================
    # IP fallbacks are NA; every IP-dependent step is gated on is.finite().
    DATASET_META <- list(
        NonNeoC = list(regions     = NONNC_REGIONS,
                       gam_dir     = NONNC_DIR,
                       IP_thresh   = NA_real_,
                       IP_y_thresh = NA_real_),
        NeoC    = list(regions     = NEOC_REGIONS,
                       gam_dir     = NEOC_DIR,
                       IP_thresh   = NA_real_,
                       IP_y_thresh = NA_real_)
    )

    # == Load IPRS-derived IP values ==========================================
    # Same candidate locations, in the same order, as Parts 10b and 10c.
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
    iprs_rds_path <- if (length(.iprs_hit) > 0L) .iprs_hit[1]
                     else .iprs_candidates[1]
    if (length(.iprs_hit) > 0L && !identical(iprs_rds_path,
                                             .iprs_candidates[1]))
        cat("  [IPRS] Not in this arm's Part10b_SKCS folder; using:",
            iprs_rds_path, "\n")

    ip_from <- "no IPRS RDS read"
    if (file.exists(iprs_rds_path)) {
        iprs_rds <- readRDS(iprs_rds_path)
        # Field presence checked before indexing.
        ip_vec   <- iprs_rds$inflection_points_DFG
        # IP_y = f(IP_x), the fitted cDFG at the inflection point.
        ip_y_vec <- iprs_rds$inflection_points_DFG_y
        ip_from <- if (!is.null(iprs_rds$ip_definition))
            switch(iprs_rds$ip_definition,
                   Si0_full  = "Si=0 median GOI, full population",
                   Si0_core  = "Si=0 median GOI, SKCS-core",
                   IPRS_full = "IPRS R(x)=0, full population",
                   IPRS_core = "IPRS R(x)=0, SKCS-core",
                   iprs_rds$ip_definition)
        else "unrecorded definition (RDS predates ip_definition)"

        if (!is.null(ip_vec) && length(ip_vec) > 0L)
            cat("  [IPRS] IP definition:", ip_from, "\n")
        if (is.null(iprs_rds$ip_definition))
            warning("This IPRS RDS predates ip_definition.", call. = FALSE)

        if (is.null(ip_vec) || length(ip_vec) == 0L) {
            cat("  [IPRS] RDS has no usable inflection_points_DFG field:",
                iprs_rds_path, "\n")
            cat("  [IPRS] Step 4 will be centred on the copula median.\n")
        } else {
            for (ds in DATASETS) {
                ip_val <- suppressWarnings(as.numeric(ip_vec[ds]))
                if (length(ip_val) == 1L && !is.na(ip_val)) {
                    DATASET_META[[ds]]$IP_thresh <- ip_val
                    cat(sprintf("  [IPRS] %-10s IP_thresh set to %.6f (%s)\n",
                                ds, ip_val, ip_from))
                } else {
                    cat(sprintf("  [IPRS] %-10s IP_thresh unavailable in IPRS RDS\n",
                                ds))
                }
                ip_y_val <- if (is.null(ip_y_vec)) NA_real_
                            else suppressWarnings(as.numeric(ip_y_vec[ds]))
                if (length(ip_y_val) == 1L && !is.na(ip_y_val)) {
                    DATASET_META[[ds]]$IP_y_thresh <- ip_y_val
                    cat(sprintf("  [IPRS] %-10s IP_y_thresh set to %.6f  = f(IP_thresh)\n",
                                ds, ip_y_val))
                } else {
                    cat(sprintf("  [IPRS] %-10s IP_y_thresh unavailable in IPRS RDS\n",
                                ds))
                }
            }
        }
    } else {
        cat("  [IPRS] RDS not found. Looked in:\n")
        for (.p in .iprs_candidates) cat("      ", .p, "\n", sep = "")
        cat("  [IPRS] Step 4 will be centred on the copula median.\n")
    }

    cat("\n", strrep("=", 72), "\n")
    cat("  Part10d - Empirical copula: S_i, K_i, EP_i, R(IP_x)=0\n")
    cat("  Arm          :", ARM_TAG, "\n")
    cat("  GOI          :", GOI, "\n")
    cat("  DFG          :", DFG_COMPOSITE_LABEL, " - ",
        paste(DFG, collapse = ", "), "\n")
    cat("  Normalisation:", NORM_TAG, " log1p(UMIs/total x ",
        format(NORM_TARGET_SUM, scientific = FALSE), ")\n")
    cat("  Datasets     :", paste(DATASETS, collapse = ", "), "\n")
    cat("  Part 10b EP  :", dirname(p10d_ep_csv(DATASETS[1])), "\n")
    cat("  Output dir   :", OUT_DIR, "\n")
    cat(strrep("=", 72), "\n\n")

    # Marker file, removed when the arm completes.
    INCOMPLETE_MARKER_P10D <- file.path(OUT_DIR, "_RUN_INCOMPLETE.txt")
    writeLines(c(
        "This arm did not finish.",
        "",
        paste("Arm     :", arm_tag),
        paste("Part    :", PART_DIR_TAG),
        paste("Started :", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
    ), INCOMPLETE_MARKER_P10D)


    # =============================================================================
    # SECTION A0: PRE-ANALYSIS GUARDRAILS
    # =============================================================================

    GR_LOG_10D <- file.path(OUT_DIR, paste0(FILE_PREFIX, "_Guardrail_Log.txt"))

    .gr10d_log  <- list()
    .gr10d_pass <- 0L
    .gr10d_fail <- 0L
    .gr10d_warn <- 0L

    gr10d <- function(check_id, status, detail) {
        # A zero-length detail would print nothing, so it is replaced.
        if (length(detail) != 1L || is.na(detail))
            detail <- paste0("[detail unavailable: the message was built from a",
                             " missing value and came out zero-length]")
        # Unknown status values still produce a tagged row.
        tag <- switch(status, PASS = "[PASS]", FAIL = "[FAIL]", WARN = "[WARN]",
                      sprintf("[%s]", as.character(status)[1]))
        cat(sprintf("  %s  %-16s  %s\n", tag, check_id, detail))
        .gr10d_log[[length(.gr10d_log) + 1L]] <<- list(
            id = check_id, status = status, detail = detail)
        if (status == "PASS") .gr10d_pass <<- .gr10d_pass + 1L
        if (status == "FAIL") .gr10d_fail <<- .gr10d_fail + 1L
        if (status == "WARN") .gr10d_warn <<- .gr10d_warn + 1L
    }

    gr10d_write <- function() {
        con <- file(GR_LOG_10D, "wt")
        on.exit(try(close(con), silent = TRUE), add = TRUE)
        writeLines(c(
            strrep("=", 72),
            paste0("Part 10d Copula2 - Guardrail Log [", NORM_TAG, "]"),
            paste("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
            paste("Arm :", ARM_TAG),
            paste0("Normalisation: adata.X = log1p(UMIs/total x ",
                   format(NORM_TARGET_SUM, scientific = FALSE),
                   "), exported unchanged"),
            strrep("=", 72), ""
        ), con)
        for (item in .gr10d_log)
            writeLines(sprintf("  [%s]  %-16s  %s", item$status, item$id,
                               item$detail), con)
        writeLines(c("", strrep("=", 72),
                     sprintf("  Total : %d",
                             .gr10d_pass + .gr10d_fail + .gr10d_warn),
                     sprintf("  PASS  : %d", .gr10d_pass),
                     sprintf("  WARN  : %d", .gr10d_warn),
                     sprintf("  FAIL  : %d", .gr10d_fail),
                     if (.gr10d_fail == 0 && .gr10d_warn == 0)
                         "\n  Result: no failures, no warnings"
                     else if (.gr10d_fail == 0)
                         "\n  Result: no failures; warnings present"
                     else "\n  Result: failures present",
                     strrep("=", 72)), con)
        close(con)
        cat("  Guardrail log written:", GR_LOG_10D, "\n")
    }

    cat("\n", strrep("=", 72), "\n")
    cat("  SECTION A0 - PRE-ANALYSIS GUARDRAILS\n")
    cat(strrep("=", 72), "\n\n")

    # Positive tests on the CP10K normalisation.
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

            gr10d(paste0("GR-01-", .ds$tag),
                  if (n_max > 0 && n_max <= NORM_CEILING + 1e-6) "PASS"
                  else "FAIL",
                  sprintf("norm.mtx max = %.4f  (ceiling log1p(%s) = %.4f)",
                          n_max,
                          format(NORM_TARGET_SUM, scientific = FALSE),
                          NORM_CEILING))

            n_rows <- min(50L, nrow(X_n))
            rs     <- rowSums(expm1(X_n[seq_len(n_rows), , drop = FALSE]))
            rs_dev <- max(abs(rs - NORM_TARGET_SUM))
            gr10d(paste0("GR-02-", .ds$tag),
                  if (rs_dev < NORM_TARGET_SUM * 1e-6) "PASS" else "FAIL",
                  sprintf("expm1(norm) row sums: mean = %.6f, max |dev| = %.2e, n = %d",
                          mean(rs), rs_dev, n_rows))

            rm(X_n); gc(verbose = FALSE)
        } else {
            gr10d(paste0("GR-01-", .ds$tag), "FAIL",
                  paste("norm.mtx not found:", norm_mtx_path))
            gr10d(paste0("GR-02-", .ds$tag), "FAIL",
                  paste("norm.mtx not found, so the CP10K row-sum identity",
                        "could not be checked:", norm_mtx_path))
        }

        if (file.exists(raw_mtx_path)) {
            X_r    <- as(readMM(raw_mtx_path), "CsparseMatrix")
            r_samp <- as.numeric(X_r@x[seq_len(min(10000, length(X_r@x)))])
            is_int <- all(abs(r_samp - round(r_samp)) < 0.01)
            r_max  <- max(X_r)
            gr10d(paste0("GR-04-", .ds$tag),
                  if (is_int && r_max > 20) "PASS" else "FAIL",
                  sprintf("raw max = %.0f, integer-like = %s",
                          r_max, if (is_int) "TRUE" else "FALSE"))
            rm(X_r); gc(verbose = FALSE)
        } else {
            gr10d(paste0("GR-04-", .ds$tag), "FAIL",
                  paste("raw_selected matrix not found:", raw_mtx_path))
        }
    }

    # == Part 10b EP/TREP CSV exists (FAIL if absent) =========================
    for (ds in DATASETS) {
        .p     <- p10d_ep_csv(ds)
        .found <- file.exists(.p)
        gr10d(paste0("GR-05-", ds),
              if (.found) "PASS" else "FAIL",
              sprintf("Part 10b EP/TREP CSV: %s",
                      if (.found)
                          file.path(basename(dirname(.p)), basename(.p))
                      else paste0(
                          "not found in any location: ",
                          paste(p10d_ep_candidates(ds), collapse = "  |  "))))
    }

    # Written here and at the end, so an aborted arm still leaves a ledger.
    gr10d_write()

    cat(sprintf("\n  Pre-analysis guardrails: %d PASS / %d WARN / %d FAIL\n",
                .gr10d_pass, .gr10d_warn, .gr10d_fail))
    if (.gr10d_fail > 0) {
        cat("  Guardrail failures present\n\n")
        warning("Part 10d pre-analysis guardrails: ", .gr10d_fail, " FAIL(s)",
                call. = FALSE)
    } else if (.gr10d_warn > 0) {
        cat(sprintf("  %d passed, %d warning(s), no failures.\n\n",
                    .gr10d_pass, .gr10d_warn))
    } else {
        cat("  All pre-analysis guardrails passed\n\n")
    }


    # =============================================================================
    # SECTION B: MAIN FUNCTION (runs per dataset)
    # =============================================================================


    run_part10d <- function(dataset_tag) {

        set.seed(42)

        meta      <- DATASET_META[[dataset_tag]]
        gam_dir   <- meta$gam_dir
        regions   <- meta$regions
        IP_thresh   <- meta$IP_thresh
    IP_y_thresh <- meta$IP_y_thresh
    if (is.null(IP_y_thresh)) IP_y_thresh <- NA_real_

        cat("\n", strrep("=", 72), "\n")
        cat("  Part10d Copula2:", dataset_tag, " [", ARM_TAG, "]\n")
        cat(strrep("=", 72), "\n\n")

        all_plots  <- list()
        excel_tabs <- list()

        # == Load GOI-positive GAM data =======================================
        region_dfs <- list()
        for (reg in regions) {
            gp <- file.path(gam_dir, "gam_data",
                            paste0(reg, "_",
                                   gsub("-", "_", gsub("/", "_", PRIMARY_STYPE)),
                                   "_selected.csv"))
            if (!file.exists(gp)) next
            df        <- read.csv(gp, stringsAsFactors = FALSE, check.names = FALSE)
            df$region <- reg
            df        <- df[df$GOI > 0, ]
            if (nrow(df) >= MIN_CELLS) region_dfs[[reg]] <- df
        }
        if (length(region_dfs) == 0) {
            cat("  ERROR: no GAM data found\n"); return(invisible(NULL))
        }
        df_all <- do.call(rbind, region_dfs)
        rownames(df_all) <- NULL
        cat(sprintf("  GOI-positive cells loaded: %d\n", nrow(df_all)))

        # == Load SKCS-core cells from Part10b EP CSV =========================
        ep_csv <- p10d_ep_csv(dataset_tag)
        if (!file.exists(ep_csv)) {
            cat("  ERROR: Part 10b EP/TREP CSV not found. Looked in:\n")
            for (.p in p10d_ep_candidates(dataset_tag)) cat("      ", .p, "\n", sep = "")
            return(invisible(NULL))
        }
        # check.names = FALSE so hyphenated column names are kept.
        skcs_df <- read.csv(ep_csv, stringsAsFactors = FALSE, check.names = FALSE)
        cat(sprintf("  SKCS-core cells loaded: %d\n", nrow(skcs_df)))

        # == Verify required columns ==========================================
        # Column names as they appear in Part10b_*_EP_TREP_values.csv
        required_cols <- c("barcode", "IP_status", "is_TREP", "EP_value", "K_value", "S_value")
        missing_cols  <- setdiff(required_cols, names(skcs_df))
        if (length(missing_cols) > 0) {
            cat("  ERROR: Missing columns in Part10b CSV:",
                paste(missing_cols, collapse = ", "), "\n")
            cat("  Available columns:", paste(names(skcs_df), collapse = ", "), "\n")
            return(invisible(NULL))
        }

        # == Match SKCS cells to expression data ==============================
        in_skcs <- df_all$barcode %in% skcs_df$barcode
        df_skcs <- df_all[in_skcs, ]
        cat(sprintf("  Matched %d / %d SKCS cells to GAM data\n\n",
                    sum(in_skcs), nrow(skcs_df)))
        # EP CSV barcodes with no row in the GAM data are dropped and counted.
        n_ep_unmatched <- sum(!(skcs_df$barcode %in% df_all$barcode))

        # Merge EP/TREP metadata
        df_skcs <- merge(df_skcs,
                         skcs_df[, required_cols],
                         by = "barcode", all.x = FALSE)
        n_SKCS <- nrow(df_skcs)     # cells in SKCS-core (for internal use)
        n_full <- nrow(df_all)      # ALL GOI-positive cells (matches Part10b rank population)
        n      <- n_full            # n used in scale_factor and K_i formula = n_full
        cat(sprintf("  Working dataset: %d SKCS-core cells | %d GOI-positive (full pop)\n\n",
                    n_SKCS, n_full))

        # == S_value (already normalised in Part 10b) =========================
        Si_obs_raw <- df_skcs$S_value   # kept for the workbook
        Si_obs     <- df_skcs$S_value   # use directly -- already normalised in Part10b

        # == Core vectors =====================================================
        x        <- df_skcs$GOI           # GOI expression
        y        <- df_skcs$Expression    # cDFG expression
        Ki_obs   <- df_skcs$K_value       # K_i (already normalised, [-1,+1]) from Part10b
        EP_obs_raw <- df_skcs$EP_value     # raw EP_i (can be very negative for non-TREP cells)
        # Clip EP for plotting/modelling: cells where model is far worse than null
        # (large negative EP) create extreme outliers that compress the y-axis.
        EP_clip_lo <- quantile(df_skcs$EP_value, 0.02, na.rm = TRUE)
        EP_obs     <- pmax(df_skcs$EP_value, EP_clip_lo)  # clipped at 2nd percentile
        cat(sprintf("  EP_value range: [%.1f, %.4f] | clipped at 2nd pct (%.4f): [%.4f, %.4f]\n\n",
                    min(df_skcs$EP_value, na.rm=TRUE), max(df_skcs$EP_value, na.rm=TRUE),
                    EP_clip_lo, min(EP_obs, na.rm=TRUE), max(EP_obs, na.rm=TRUE)))
        ip_stat  <- df_skcs$IP_status     # "Pre-IP" / "Post-IP"
        # Stop if is_TREP does not coerce to logical.
        is_trep  <- as.logical(df_skcs$is_TREP)  # coerce: CSV may read as character
        if (anyNA(is_trep)) {
            .bad <- unique(as.character(df_skcs$is_TREP)[is.na(is_trep)])
            stop(sprintf(paste0("%s: is_TREP in the Part 10b EP/TREP CSV did ",
                                "not coerce to logical for %d of %d row(s). ",
                                "Unrecognised value(s): %s."),
                        dataset_tag, sum(is.na(is_trep)), length(is_trep),
                        paste(utils::head(.bad, 5), collapse = ", ")),
                 call. = FALSE)
        }
        if (any(is.na(is_trep))) {
            .bad <- unique(as.character(df_skcs$is_TREP)[is.na(is_trep)])
            stop(sprintf(paste0("%s: %d of %d is_TREP value(s) did not coerce to ",
                                "TRUE/FALSE (saw: %s). Check ",
                                "<prefix>_%s_EP_TREP_values.csv."),
                        dataset_tag, sum(is.na(is_trep)), length(is_trep),
                        paste(utils::head(.bad, 4), collapse = ", "), dataset_tag),
                 call. = FALSE)
        }

        # Pre/Post labels recomputed on both axes to check the axis Part 10b used.
        ip_axis_check <- list(n = 0L, agree_cdfg = NA_real_, agree_goi = NA_real_)
        .have_lab <- !is.na(ip_stat)
        if (any(.have_lab)) {
            .lab_goi  <- if (is.finite(IP_thresh))
                ifelse(df_skcs$GOI[.have_lab]        <= IP_thresh,
                       "Pre-IP", "Post-IP") else NULL
            .lab_cdfg <- if (is.finite(IP_y_thresh))
                ifelse(df_skcs$Expression[.have_lab] <  IP_y_thresh,
                       "Pre-IP", "Post-IP") else NULL
            ip_axis_check <- list(
                n          = sum(.have_lab),
                agree_cdfg = if (is.null(.lab_cdfg)) NA_real_
                             else mean(.lab_cdfg == ip_stat[.have_lab]),
                agree_goi  = if (is.null(.lab_goi))  NA_real_
                             else mean(.lab_goi  == ip_stat[.have_lab]))
            cat(sprintf(paste0("  Part 10b Pre/Post labels reproduce: %s on the",
                               " cDFG axis at IP_y, %s on the GOI axis at IP_x",
                               "  (n=%d)\n"),
                        if (is.finite(ip_axis_check$agree_cdfg))
                            sprintf("%.4f", ip_axis_check$agree_cdfg) else "NA",
                        if (is.finite(ip_axis_check$agree_goi))
                            sprintf("%.4f", ip_axis_check$agree_goi) else "NA",
                        ip_axis_check$n))
        }

        # Four-group label
        group <- paste(ip_stat, ifelse(is_trep, "TREP", "non-TREP"))
        group[!group %in% names(FOUR_COLS)] <- "Unknown"
        group <- factor(group, levels = names(FOUR_COLS))

        # == Full-population copula coordinates ===============================
        # Ranks over all n_full GOI-positive cells, as in Part 10b, so that
        # S_i = ((n_full+1)^2/(n_full-1)) * (u_i-0.5)(v_i-0.5).
        x_full <- df_all$GOI
        y_full <- df_all$Expression
        u_full_all <- rank(x_full, ties.method = "average") / (n_full + 1)
        v_full_all <- rank(y_full, ties.method = "average") / (n_full + 1)
        # Extract SKCS-core subset by barcode position in df_all
        skcs_idx_in_full <- match(df_skcs$barcode, df_all$barcode)
        # Every SKCS barcode must be in df_all.
        if (any(is.na(skcs_idx_in_full)))
            stop(sprintf(paste0("%s: %d SKCS barcode(s) have no match in the ",
                                "GOI-positive frame."),
                        dataset_tag, sum(is.na(skcs_idx_in_full))),
                 call. = FALSE)
        u <- u_full_all[skcs_idx_in_full]   # full-pop copula u for SKCS cells
        v <- v_full_all[skcs_idx_in_full]   # full-pop copula v for SKCS cells
        cat(sprintf("  Copula coords: u range [%.4f, %.4f] | v range [%.4f, %.4f]\n",
                    min(u), max(u), min(v), max(v)))
        cat(sprintf("  (full-pop ranks: n_full=%d -> u,v not uniform over [0,1]\n\n",
                    n_full))

        # == Non-selected (gray) cells ========================================
        # GOI-positive cells outside the SKCS-core, ranked on the full population.
        non_sel_bc  <- setdiff(df_all$barcode, df_skcs$barcode)
        n_non_sel   <- length(non_sel_bc)
        cat(sprintf("  Non-selected (gray) cells: %d  (n_full=%d - n_SKCS=%d)\n\n",
                    n_non_sel, n_full, n_SKCS))

        # Non-selected cells are counted but not plotted (no Si_obs/Ki_obs available)

        # == Full-population empirical copula C_hat ===========================
        # Over all n_full GOI-positive cells, the population of K_value in Part 10b.
        cat(sprintf("  Computing C_hat over all %d GOI-positive cells (O(n^2)) ...\n",
                    n_full))
        C_hat_full_all <- compute_ecdf_copula(u_full_all, v_full_all)
        C_hat <- C_hat_full_all[skcs_idx_in_full]   # SKCS subset
        cat(sprintf("  C_hat computed. SKCS subset: n=%d | range [%.4f, %.4f]\n\n",
                    length(C_hat), min(C_hat), max(C_hat)))

        # == GAM image curve in copula space (TREP cells only) ================
        trep_idx   <- which(is_trep)
        u_trep     <- u[trep_idx]
        v_trep     <- v[trep_idx]
        cat(sprintf("  Fitting GAM image curve g(u) on %d TREP cells (PRSS phase1+phase2) ...\n",
                    length(trep_idx)))
        # g(u) needs at least 4 distinct u values for a k = 3 basis.
        .n_u_trep <- length(unique(u_trep))
        if (.n_u_trep < 4L)
            stop(sprintf(paste0("%s: the TREP set has only %d cell(s) with %d ",
                                "distinct u value(s), which cannot support the ",
                                "g(u) image curve (a k = 3 thin-plate basis ",
                                "needs at least 4); see <prefix>_%s_",
                                "EP_TREP_values.csv."),
                        dataset_tag, length(trep_idx), .n_u_trep, dataset_tag),
                 call. = FALSE)

        gam_image_result <- fit_gam_image_prss(u_trep, v_trep,
                                                label     = "g(u)",
                                                num_iters = 100L,
                                                seed_base = 42L)
        gam_image  <- gam_image_result$best_model
        if (is.null(gam_image))
            stop(sprintf(paste0("%s: every g(u) GAM fit failed (best PRSS = ",
                                "%s over %d TREP cell(s), %d distinct u)."),
                        dataset_tag,
                        format(gam_image_result$best_prss),
                        length(trep_idx), .n_u_trep),
                 call. = FALSE)
        gam_img_sm <- summary(gam_image)
        cat(sprintf("  g(u) PRSS result: best_k=%d | best_iter=%d | best_PRSS=%.5f | early_stopped=%s\n",
                    gam_image_result$best_k,
                    gam_image_result$best_iter,
                    gam_image_result$best_prss,
                    gam_image_result$early_stopped))
        cat(sprintf("  g(u) GAM summary: edf=%.2f | Dev=%.1f%% | p=%.2e\n\n",
                    gam_img_sm$edf,
                    gam_img_sm$dev.expl * 100,
                    gam_img_sm$s.table[4]))

        # Predict g(u_i) for all cells
        g_u_all <- as.numeric(predict(gam_image,
                                      newdata = data.frame(u_trep = u)))

        # Fine grid for g(u) curve in plots
        u_grid    <- seq(0.001, 0.999, length.out = 500)
        g_u_grid_pred <- predict(gam_image,
                                 newdata   = data.frame(u_trep = u_grid),
                                 se.fit    = TRUE)
        g_u_grid  <- as.numeric(g_u_grid_pred$fit)
        g_u_se    <- as.numeric(g_u_grid_pred$se.fit)
        g_curve_df <- data.frame(
            u       = u_grid,
            g_u     = pmin(1, pmax(0, g_u_grid)),
            g_upper = pmin(1, pmax(0, g_u_grid + 1.96 * g_u_se)),
            g_lower = pmin(1, pmax(0, g_u_grid - 1.96 * g_u_se))
        )

        # IP_u = #{j in all GOI+ cells: x_j <= IP_thresh} / (n_full + 1)
        IP_u <- if (is.finite(IP_thresh))
            sum(x_full <= IP_thresh) / (n_full + 1) else NA_real_

        # Without an IP, Step 4 centres on u_REF = 0.5 and draws no IP line.
    # IP_v = Fc(IP_y), the IP projected onto the cDFG marginal.
    IP_v <- if (is.finite(IP_y_thresh))
        sum(y_full <= IP_y_thresh) / (n_full + 1) else NA_real_
    if (is.finite(IP_v))
        cat(sprintf(paste0("  IP in copula space: (IP_u, IP_v) = (%.4f, %.4f)",
                           "   [IP_y = %.4f]\n"),
                    IP_u, IP_v, IP_y_thresh))
    else
        cat("  IP_y not available from the IPRS RDS, so IP_v is not computed.\n")

    IP_AVAILABLE <- is.finite(IP_u)
        u_REF        <- if (IP_AVAILABLE) IP_u else 0.5
        ip_vline <- function(colour = "#0D1B4B", linewidth = 0.55,
                             linetype = "dashed") {
            if (!IP_AVAILABLE) return(list())
            list(geom_vline(xintercept = IP_u, linetype = linetype,
                            colour = colour, linewidth = linewidth))
        }
        if (IP_AVAILABLE)
            cat(sprintf(paste0("  IP threshold: IP_x = %.3f  ->  IP_u = Fx(IP_x)",
                               " = %.4f  (full-pop n=%d)\n\n"),
                        IP_thresh, IP_u, n_full))
        else
            cat(sprintf(paste0("  IP threshold: not available; Step 4 centred on",
                               " u = 0.5  (full-pop n=%d)\n\n"), n_full))


        # =============================================================================
        # STEP 1: S_i as copula quadrant displacement
        # S_i = ((n+1)^2/(n-1)) * (u_i - 0.5)(v_i - 0.5)
        # =============================================================================
        cat("  ── Step 1: S_i copula computation ──\n")

        # scale_factor = (n_full+1)^2/(n_full-1): uses full-population n
        scale_factor <- (n + 1)^2 / (n - 1)   # n = n_full here
        Si_formula   <- scale_factor * (u - 0.5) * (v - 0.5)
        Si_formula_all <- scale_factor * (u_full_all - 0.5) * (v_full_all - 0.5)

        # =============================================================================
        # STEP 2: K_i as empirical copula formula
        # K_i = (n+1)/(n-1) + [4n/(n-1)]*C_hat_i - [2n/(n-1)]*(u_i + v_i)
        # =============================================================================
        cat("  ── Step 2: K_i copula computation ──\n")

        # == Ki over all n_full GOI-positive cells (tau_b, as in Part 10b) ====
        cat(sprintf("  Computing Ki ordinal copula over ALL %d GOI+ cells ...\n", n_full))
        ki_full_result <- compute_ordinal_copula_ki(u_full_all, v_full_all)
        Ki_formula     <- ki_full_result$Ki_formula[skcs_idx_in_full]
        Ki_formula_all <- ki_full_result$Ki_formula   # ALL n_full cells (for all-cell plots)
        C_strict_lower <- ki_full_result$C_strict_lower[skcs_idx_in_full]
        C_strict_upper <- ki_full_result$C_strict_upper[skcs_idx_in_full]
        Ki_approx      <- Ki_formula   # kept for the workbook
        cat("  Done.\n\n")

        # =============================================================================
        # STEP 4: The inflection point in copula space
        # (i) the IP as (IP_u, IP_v); (ii) the copula density at IP_u against
        # u = 0.25 and u = 0.75, one-sided t-tests on the strips around each.
        # =============================================================================
        cat("  ── Step 4: the inflection point in copula space ──\n")

        # == (i) the IP as a copula coordinate pair ===========================
        cat(sprintf("  (i)  IP_u = Fx(IP_x=%.4f) = %.4f | IP_v = Fc(IP_y=%s) = %s\n",
                    IP_thresh, IP_u,
                    if (is.finite(IP_y_thresh)) sprintf("%.4f", IP_y_thresh) else "NA",
                    if (is.finite(IP_v)) sprintf("%.4f", IP_v) else "NA"))


        # == (ii) Copula density along the diagonal ===========================
        # kde2d on (u, v), limits extended by COP_PAD_4 so contours are not clipped
        COP_PAD_4 <- 0.20
        bw_u_d <- bw.nrd0(u)
        bw_v_d <- bw.nrd0(v)
        dens2d_uv <- kde2d(u, v, n = 150,
                           h    = c(bw_u_d * 2.5, bw_v_d * 2.5),   # wider bw -> halos extend beyond dots
                           lims = c(0 - COP_PAD_4, 1 + COP_PAD_4,
                                    0 - COP_PAD_4, 1 + COP_PAD_4))
        dens_uv_df <- expand.grid(u_g = dens2d_uv$x, v_g = dens2d_uv$y)
        dens_uv_df$density       <- as.vector(dens2d_uv$z)
        dens_uv_df$log1p_density <- log1p(dens_uv_df$density)   # log-normalised for p5

        # Density on the diagonal: sample at v = u (interpolated via nearest grid)
        diag_u    <- seq(0.05, 0.95, by = 0.01)
        diag_dens <- sapply(diag_u, function(uu) {
            iu <- which.min(abs(dens2d_uv$x - uu))
            iv <- which.min(abs(dens2d_uv$y - uu))
            dens2d_uv$z[iu, iv]
        })
        diag_df <- data.frame(u = diag_u, density_on_diagonal = diag_dens)

        # Statistical comparison: density at IP_u vs u=0.25 and u=0.75
        # Extract density values in a strip around each u-slice
        strip_w <- 0.04
        get_strip_density <- function(u_center) {
            dens_uv_df$density[abs(dens_uv_df$u_g - u_center) <= strip_w &
                                   abs(dens_uv_df$u_g - dens_uv_df$v_g) <= 0.15]
        }

        d_at_025   <- get_strip_density(0.25)
        d_at_IP_u   <- get_strip_density(u_REF)
        d_at_075   <- get_strip_density(0.75)

        # One-sided t-test: is density at IP_u lower than at u=0.25 and u=0.75?
        ttest_025_IP_u <- if (length(d_at_IP_u) >= 3 && length(d_at_025) >= 3)
            t.test(d_at_025, d_at_IP_u, alternative = "greater") else NULL
        ttest_075_IP_u <- if (length(d_at_IP_u) >= 3 && length(d_at_075) >= 3)
            t.test(d_at_075, d_at_IP_u, alternative = "greater") else NULL

        # == Cells in each strip ==============================================
        WAIST_MIN_CELLS_10D <- 10L
        # Flank references: a minimum cell count and the uniform-copula density.
        WAIST_FLANK_MIN_10D <- max(20L, as.integer(ceiling(0.03 * length(u))))
        WAIST_UNIF_REF_10D  <- 1 / (1 + 2 * COP_PAD_4)^2
        strip_cells_u <- function(u_center)
            sum(abs(u - u_center) <= strip_w & abs(u - v) <= 0.15)
        n_cells_025  <- strip_cells_u(0.25)
        n_cells_IP_u  <- strip_cells_u(u_REF)
        n_cells_075  <- strip_cells_u(0.75)
        waist_occupied <- n_cells_IP_u >= WAIST_MIN_CELLS_10D

        # == The same test in the v direction, at IP_v ========================
        v_REF <- if (is.finite(IP_v)) IP_v else 0.5
        get_strip_density_v <- function(v_center) {
            dens_uv_df$density[abs(dens_uv_df$v_g - v_center) <= strip_w &
                                   abs(dens_uv_df$u_g - dens_uv_df$v_g) <= 0.15]
        }
        strip_cells_v <- function(v_center)
            sum(abs(v - v_center) <= strip_w & abs(u - v) <= 0.15)

        d_at_v025  <- get_strip_density_v(0.25)
        d_at_IP_v   <- get_strip_density_v(v_REF)
        d_at_v075  <- get_strip_density_v(0.75)
        n_cells_v025 <- strip_cells_v(0.25)
        n_cells_IP_v  <- strip_cells_v(v_REF)
        n_cells_v075 <- strip_cells_v(0.75)
        waist_v_occupied <- n_cells_IP_v >= WAIST_MIN_CELLS_10D

        ttest_025_IP_v <- if (length(d_at_IP_v) >= 3 && length(d_at_v025) >= 3)
            t.test(d_at_v025, d_at_IP_v, alternative = "greater") else NULL
        ttest_075_IP_v <- if (length(d_at_IP_v) >= 3 && length(d_at_v075) >= 3)
            t.test(d_at_v075, d_at_IP_v, alternative = "greater") else NULL

        cat(sprintf("  (ii) Copula density waist test (strip ±%.2f around diagonal):\n", strip_w))
        cat(sprintf("    Mean density at u=0.25: %.5f (n=%d)\n",
                    mean(d_at_025), length(d_at_025)))
        cat(sprintf("    Mean density at IP_u:  %.5f (n=%d)\n",
                    mean(d_at_IP_u), length(d_at_IP_u)))
        cat(sprintf("    Mean density at u=0.75: %.5f (n=%d)\n",
                    mean(d_at_075), length(d_at_075)))
        if (!is.null(ttest_025_IP_u))
            cat(sprintf("    t-test u=0.25 > IP_u:  t=%.3f  p=%.4f\n",
                        ttest_025_IP_u$statistic, ttest_025_IP_u$p.value))
        if (!is.null(ttest_075_IP_u))
            cat(sprintf("    t-test u=0.75 > IP_u:  t=%.3f  p=%.4f\n",
                        ttest_075_IP_u$statistic, ttest_075_IP_u$p.value))
        cat(sprintf(paste0("    cells in strip: %d (u=0.25), %d (IP_u), ",
                           "%d (u=0.75)%s\n"),
                    n_cells_025, n_cells_IP_u, n_cells_075,
                    if (waist_occupied) ""
                    else sprintf(paste0("   fewer than %d cells at",
                                        " IP_u: near-empty region"),
                                 WAIST_MIN_CELLS_10D)))

        cat(sprintf(paste0("\n  (iii) Copula density waist test in the v ",
                           "direction, at v_%s:\n"),
                    if (is.finite(IP_v)) "IP" else "median"))
        cat(sprintf("    Mean density at v=0.25: %.5f (n=%d)\n",
                    mean(d_at_v025), length(d_at_v025)))
        cat(sprintf("    Mean density at IP_v:  %.5f (n=%d)\n",
                    mean(d_at_IP_v), length(d_at_IP_v)))
        cat(sprintf("    Mean density at v=0.75: %.5f (n=%d)\n",
                    mean(d_at_v075), length(d_at_v075)))
        if (!is.null(ttest_025_IP_v))
            cat(sprintf("    t-test v=0.25 > IP_v:  t=%.3f  p=%.4f\n",
                        ttest_025_IP_v$statistic, ttest_025_IP_v$p.value))
        if (!is.null(ttest_075_IP_v))
            cat(sprintf("    t-test v=0.75 > IP_v:  t=%.3f  p=%.4f\n",
                        ttest_075_IP_v$statistic, ttest_075_IP_v$p.value))
        cat(sprintf(paste0("    cells in strip: %d (v=0.25), %d (IP_v), ",
                           "%d (v=0.75)%s\n\n"),
                    n_cells_v025, n_cells_IP_v, n_cells_v075,
                    if (waist_v_occupied) ""
                    else sprintf(paste0("   fewer than %d cells at",
                                        " IP_v: near-empty region"),
                                 WAIST_MIN_CELLS_10D)))

        # Plot 5 scales its log-normalised fill against this maximum.
        cop_log1p_max4 <- max(dens_uv_df$log1p_density)

        # Workbook sheets for Step 4
        excel_tabs[["Step4_DiagDensity"]] <- diag_df   # density along copula diagonal

        # == Estimation basis =================================================
        excel_tabs[["Estimation_Basis"]] <- data.frame(
            dataset  = dataset_tag,
            quantity = c("g(u) image curve",
                         "Ki-Si Spearman rho / Kendall tau",
                         "Ki-Si tail dependence lambda_L, lambda_U",
                         "CvM copula distance and its permutation test",
                         "copula density waist t-tests (u and v directions)"),
            unit     = "cell",
            donor_random_intercept = FALSE,
            inverse_donor_weights  = FALSE,
            note = paste0(
                "Cell level: no donor random intercept or donor weights."),
            stringsAsFactors = FALSE
        )

        excel_tabs[["Step4_IP_CopulaSpace"]] <- data.frame(
            dataset               = dataset_tag,
            IP_x                  = IP_thresh,
            IP_u_Fx_IP_x           = round(IP_u, 5),
            IP_u_vs_05            = round(IP_u - 0.5, 5),
            IP_y                  = round(IP_y_thresh, 6),
            IP_v_Fc_IP_y           = round(IP_v, 5),
            mean_density_u025     = round(mean(d_at_025), 5),
            mean_density_IP_u      = round(mean(d_at_IP_u), 5),
            mean_density_u075     = round(mean(d_at_075), 5),
            ttest_025_IP_u_t       = round(ifelse(is.null(ttest_025_IP_u), NA,
                                                 ttest_025_IP_u$statistic), 4),
            ttest_025_IP_u_p       = round(ifelse(is.null(ttest_025_IP_u), NA,
                                                 ttest_025_IP_u$p.value), 5),
            ttest_075_IP_u_t       = round(ifelse(is.null(ttest_075_IP_u), NA,
                                                 ttest_075_IP_u$statistic), 4),
            ttest_075_IP_u_p       = round(ifelse(is.null(ttest_075_IP_u), NA,
                                                 ttest_075_IP_u$p.value), 5),
            waist_confirmed       = ifelse(
                !is.null(ttest_025_IP_u) && !is.null(ttest_075_IP_u),
                ttest_025_IP_u$p.value < 0.05 & ttest_075_IP_u$p.value < 0.05,
                NA),
            # Cells, not grid points, in each strip
            n_cells_strip_u025    = n_cells_025,
            n_cells_strip_IP_u     = n_cells_IP_u,
            n_cells_strip_u075    = n_cells_075,
            waist_region_occupied = waist_occupied,
            # The same test in the v direction
            mean_density_v025     = round(mean(d_at_v025), 5),
            mean_density_IP_v      = round(mean(d_at_IP_v), 5),
            mean_density_v075     = round(mean(d_at_v075), 5),
            ttest_025_IP_v_t       = round(ifelse(is.null(ttest_025_IP_v), NA,
                                                 ttest_025_IP_v$statistic), 4),
            ttest_025_IP_v_p       = round(ifelse(is.null(ttest_025_IP_v), NA,
                                                 ttest_025_IP_v$p.value), 5),
            ttest_075_IP_v_t       = round(ifelse(is.null(ttest_075_IP_v), NA,
                                                 ttest_075_IP_v$statistic), 4),
            ttest_075_IP_v_p       = round(ifelse(is.null(ttest_075_IP_v), NA,
                                                 ttest_075_IP_v$p.value), 5),
            waist_v_confirmed     = ifelse(
                !is.null(ttest_025_IP_v) && !is.null(ttest_075_IP_v),
                ttest_025_IP_v$p.value < 0.05 & ttest_075_IP_v$p.value < 0.05,
                NA),
            n_cells_strip_v025    = n_cells_v025,
            n_cells_strip_IP_v     = n_cells_IP_v,
            n_cells_strip_v075    = n_cells_v075,
            waist_v_region_occupied = waist_v_occupied,
            note = paste0(
                "IP_u = Fx(IP_x) = #{x<=IP_x}/(n+1); IP_v = Fc(IP_y). ",
                "Waist: density at the IP coordinate lower than flanks ",
                "(one-sided t-test).")
        )


        # =============================================================================
        # STEP 5: Synthesis, pitcher plant as copula geometry
        # Density (background), K_i (colour), EP_i (size), g(u) (dashed), IP_u.
        # =============================================================================
        cat("  ── Step 5: Synthesis -- pitcher plant as copula geometry ──\n")

        # Prepare synthesis data
        synth_df <- data.frame(
            u           = u,
            v           = v,
            Ki_obs      = Ki_obs,
            EP_obs      = EP_obs,
            Si_obs      = Si_obs,
            group       = group,
            is_trep     = is_trep,
            trep_label  = ifelse(is_trep, "TREP cells", "non-TREP cells"),
            stringsAsFactors = FALSE
        )

        # EP size range for plot: [0.5, 3.5] mm
        ep_min  <- min(EP_obs); ep_max <- max(EP_obs)

        # 5 representative EP breaks for size legend, ensuring at least one positive
        # Pad inward slightly so ggplot does not drop boundary breaks
        ep_lo <- ep_min + 0.02 * (ep_max - ep_min)
        ep_hi <- ep_max - 0.02 * (ep_max - ep_min)
        ep_breaks <- round(seq(ep_lo, ep_hi, length.out = 5), 1)
        # Guarantee the highest break is positive (use ep_max if still negative)
        if (all(ep_breaks <= 0)) ep_breaks[5] <- round(ep_hi, 1)

        # Build synthesis plot -- one panel, all cells, Ki as colour, EP as size
        p5 <- ggplot() +
            # Copula density background, log1p-normalised, transparent at low density
            geom_raster(data = dens_uv_df,
                        aes(x = u_g, y = v_g, fill = log1p_density),
                        interpolate = TRUE, na.rm = TRUE) +
            scale_fill_gradientn(
                colours  = c("#FFFFFF00", "#FFFFFF00", "#E8E0FFA0", "#C8B8F8C8",
                             "#9E7AE8E0", "#7B2FBEEE", "#4B0082F8", "#2A0050FF"),
                values   = scales::rescale(c(0, 0.20, 0.35, 0.50,
                                             0.65, 0.78, 0.90, 1.00)),
                limits   = c(0, cop_log1p_max4),
                oob      = scales::squish,
                name     = "Copula\ndensity\n(log1p)",
                na.value = "#FFFFFF00"
            ) +
            # Cell points: K_i as fill colour (turbo 5-95%), EP as size
            ggnewscale::new_scale_fill() +
            geom_point(data = synth_df,
                       aes(x = u, y = v, fill = Ki_obs, size = EP_obs),
                       shape = 21, alpha = 0.72, colour = "gray15", stroke = 0.15) +
            scale_fill_gradientn(
                colours = viridis::turbo(256)[round(seq(0.05 * 255 + 1,
                                                         0.95 * 255 + 1,
                                                         length.out = 9))],
                limits  = c(0.5, 1),
                oob     = scales::squish,
                name    = "K_i"
            ) +
            scale_size_continuous(
                range  = c(0.5, 3.5),
                breaks = ep_breaks,
                name   = expression(EP[i])
            ) +
            guides(size = guide_legend(
                override.aes = list(fill = "gray45", colour = "gray20",
                                    stroke = 0.5, alpha = 1)
            )) +
            # GAM image curve g(u), dark teal dashed
            geom_ribbon(data = g_curve_df,
                        aes(x = u, ymin = g_lower, ymax = g_upper),
                        inherit.aes = FALSE,
                        fill = "#006D6D", alpha = 0.15) +
            geom_line(data = g_curve_df,
                      aes(x = u, y = g_u),
                      inherit.aes = FALSE,
                      colour = "#006D6D", linewidth = 1.0, linetype = "dashed") +
            # IP_u line, dark navy (omitted when no IP)
            ip_vline() +
            geom_vline(xintercept = 0.5, linetype = "dashed",
                       colour = "gray50", linewidth = 0.4) +
            geom_hline(yintercept = 0.5, linetype = "dashed",
                       colour = "gray50", linewidth = 0.4) +
            scale_x_continuous(limits = c(0, 1), expand = c(0.05, 0),
                               breaks = seq(0, 1, 0.25)) +
            scale_y_continuous(limits = c(0, 1), expand = c(0.05, 0),
                               breaks = seq(0, 1, 0.25)) +
            labs(title    = paste0(dataset_tag,
                                   " -- Step 5: Pitcher plant as copula geometry"),
                 subtitle = paste0(
                     "Background: empirical copula density (log1p-normalised).  ",
                     "Point fill: K_i (turbo palette, 0.5 to 1.0).  ",
                     "Point size: EP_i (larger=higher EP).\n",
                     "Dashed dark teal: GAM image g(u).  ",
                     "Dark navy vertical: IP_u."),
                 x = "u = Fx(x_i)  [GOI marginal CDF]",
                 y = "v = Fc(y_i)  [cDFG marginal CDF]") +
            theme_copula() +
            theme(legend.position = "right", aspect.ratio = 1)

        save_pp(p5,
                file.path(OUT_DIR,
                          paste0(FILE_PREFIX, "_", dataset_tag, "_05_synthesis")),
                pdf_w = 11, pdf_h = 10, png_w = 2400, png_h = 2200)
        cat("  Plot 5 saved (synthesis: pitcher plant as copula geometry)\n\n")
        all_plots[["synthesis"]] <- p5

        # =============================================================================
        # STEP 6: Ki-Si copula analysis
        # The copula framework of Part 10c applied to (Si_norm, Ki) over all
        # n_full GOI-positive cells: SKCS-core by IP x TREP group, others gray80.
        # 06a marginal and joint densities, 06b demo cells, 06c unit square,
        # 06d expression-space scatter, 06e copula density, 06f shuffled.
        # =============================================================================
        cat("\n  -- Step 6: Ki-Si Copula Analysis (ALL GOI-positive cells) --\n\n")

        # == 6-0: Ki and Si_norm for all GOI-positive cells ===================
        # Si_norm = 4(u-0.5)(v-0.5), the scale of Part 10b spear_norm; the
        # rank-product form is kept as Si_rankprod.
        SI_NORM_SCALE_10D <- 4          # 4(u-0.5)(v-0.5), range [-1,+1]
        Si_norm     <- SI_NORM_SCALE_10D * (u_full_all - 0.5) * (v_full_all - 0.5)
        Si_rankprod <- scale_factor      * (u_full_all - 0.5) * (v_full_all - 0.5)

        # Ki: ordinal copula formula from Step 2
        # ki_full_result was computed over all n_full cells in Step 2
        Ki_all <- ki_full_result$Ki_formula

        n_ks   <- n_full

        cat(sprintf("  Ki-Si copula: n = %d GOI-positive cells (all)\n", n_ks))
        cat(sprintf("  Si_norm range: [%.4f, %.4f] | Ki range: [%.4f, %.4f]\n",
                    min(Si_norm), max(Si_norm), min(Ki_all), max(Ki_all)))
        cat(sprintf(paste0("  Plots use Si_norm = 4(u-0.5)(v-0.5); rank-product",
                           " form (range [%.1f, %.1f], scale_factor = %.1f)",
                           " kept as Si_rankprod\n"),
                    min(Si_rankprod), max(Si_rankprod), scale_factor))

        # Core vectors
        ks_x   <- Si_norm           # x-axis = Si_norm (normalized)
        ks_y   <- Ki_all           # y-axis = Ki ([-1,+1])

        # Copula coordinates: rank/(n_full+1) over ALL cells
        u_sk <- rank(ks_x, ties.method = "average") / (n_ks + 1)   # F(Si_norm)
        v_sk <- rank(ks_y, ties.method = "average") / (n_ks + 1)   # F(Ki)

        cat(sprintf("  u_sk range: [%.4f, %.4f] | v_sk range: [%.4f, %.4f]\n\n",
                    min(u_sk), max(u_sk), min(v_sk), max(v_sk)))

        # == Group assignment for all cells ===================================
        NON_SKCS_COL <- "gray80"

        ks_group <- rep("Non-SKCS", n_ks)
        skcs_grp_lkp <- setNames(as.character(group), df_skcs$barcode)
        bc_all <- df_all$barcode
        for (i in seq_along(bc_all)) {
            g <- skcs_grp_lkp[bc_all[i]]
            if (!is.na(g)) ks_group[i] <- g
        }

        ks_group_levels <- c("Non-SKCS",
                             "Pre-IP non-TREP", "Post-IP non-TREP",
                             "Pre-IP TREP", "Post-IP TREP")
        ks_group <- factor(ks_group, levels = ks_group_levels)

        COLS_KS <- c(
            "Non-SKCS"         = NON_SKCS_COL,
            "Pre-IP TREP"      = "#E69F00",
            "Pre-IP non-TREP"  = "#1B3A8C",
            "Post-IP TREP"     = "#009E73",
            "Post-IP non-TREP" = "#D94F4F"
        )
        cols_ks <- COLS_KS[names(COLS_KS) %in% levels(ks_group)]

        n_by_ks <- table(ks_group)
        grp_labs_ks <- sapply(names(cols_ks), function(g) {
            cnt <- as.integer(n_by_ks[g])
            if (is.na(cnt)) cnt <- 0L
            sprintf("%s  n=%d", g, cnt)
        })

        cat(sprintf("  Non-SKCS cells: %d | SKCS-core cells: %d\n",
                    sum(ks_group == "Non-SKCS"),
                    sum(ks_group != "Non-SKCS")))

        # == Pre-compute dependence stats =====================================
        rho_ks_orig <- cor(ks_x, ks_y, method = "spearman")
        td_ks_orig  <- tail_dep_np(u_sk, v_sk, q = 0.10)

        set.seed(999)
        ks_y_shuffled <- sample(ks_y)
        u_sk_shuf <- rank(ks_x,          ties.method = "average") / (n_ks + 1)
        v_sk_shuf <- rank(ks_y_shuffled, ties.method = "average") / (n_ks + 1)

        rho_ks_shuf <- cor(ks_x, ks_y_shuffled, method = "spearman")
        td_ks_shuf  <- tail_dep_np(u_sk_shuf, v_sk_shuf, q = 0.10)
        ks_test_Ki  <- suppressWarnings(ks.test(ks_y, ks_y_shuffled))

        cvm_ks_dist <- cvm_copula_dist(u_sk, v_sk, u_sk_shuf, v_sk_shuf,
                                        grid_size = 25)
        set.seed(42)
        cvm_ks_null <- replicate(200, {
            idx1 <- sample(n_ks, n_ks, replace = TRUE)
            idx2 <- sample(n_ks, n_ks, replace = TRUE)
            cvm_copula_dist(u_sk[idx1], v_sk[idx1],
                            u_sk[idx2], v_sk[idx2], grid_size = 20)
        })
        cvm_ks_p <- mean(cvm_ks_null >= cvm_ks_dist)
        cvm_ks_p_label <- if (cvm_ks_p < 0.001) "p < 0.001" else
                          sprintf("p = %.3f", cvm_ks_p)

        cat(sprintf("  Ki-Si original: rho = %.4f | lambda_L = %.4f | lambda_U = %.4f\n",
                    rho_ks_orig, td_ks_orig$lambda_L, td_ks_orig$lambda_U))
        cat(sprintf("  Ki-Si shuffled: rho = %.4f | lambda_L = %.4f | lambda_U = %.4f\n",
                    rho_ks_shuf, td_ks_shuf$lambda_L, td_ks_shuf$lambda_U))
        cat(sprintf("  KS marginals: D = %.4f  p = %.4f\n",
                    ks_test_Ki$statistic, ks_test_Ki$p.value))
        cat(sprintf("  CvM copula distance = %.5f  (%s)\n\n",
                    cvm_ks_dist, cvm_ks_p_label))

        ann_ks_expr <- sprintf(
            "r_s = %.3f\nlambda_L = %.3f\nlambda_U = %.3f",
            rho_ks_orig, td_ks_orig$lambda_L, td_ks_orig$lambda_U)
        ann_ks_shuf <- sprintf(
            "r_s = %.3f\nlambda_L = %.3f\nlambda_U = %.3f",
            rho_ks_shuf, td_ks_shuf$lambda_L, td_ks_shuf$lambda_U)

        # == Copula data frame (ALL cells) ====================================
        ks_copula_df <- data.frame(
            u_sk    = u_sk,
            v_sk    = v_sk,
            Si_norm = ks_x,
            Ki_raw  = ks_y,
            barcode = bc_all,
            group   = ks_group,
            stringsAsFactors = FALSE
        )

        draw_ord_ks <- c(
            "Non-SKCS"         = 0L,
            "Pre-IP non-TREP"  = 1L,
            "Post-IP non-TREP" = 2L,
            "Pre-IP TREP"      = 3L,
            "Post-IP TREP"     = 4L
        )
        ks_copula_df$draw_ord <- draw_ord_ks[
            match(as.character(ks_copula_df$group), names(draw_ord_ks))]
        ks_copula_df$draw_ord[is.na(ks_copula_df$draw_ord)] <- 0L
        ks_copula_df <- ks_copula_df[order(ks_copula_df$draw_ord), ]


        # =============================================================================
        # PLOT 06a: Marginal PDFs + joint density [EXPRESSION SPACE]
        # =============================================================================
        cat("  Building Plot 06a: Ki-Si marginal joint (expression space) ...\n")

        # Axis limits: full data range with 4% padding
        ks_x_pad <- diff(range(ks_x)) * 0.04
        ks_y_pad <- diff(range(ks_y)) * 0.04
        ks_x_lim <- range(ks_x) + c(-ks_x_pad, ks_x_pad)
        ks_y_lim <- range(ks_y) + c(-ks_y_pad, ks_y_pad)

        ks_bw_x <- bw.nrd0(ks_x) * 1.2
        ks_bw_y <- bw.nrd0(ks_y) * 1.2
        ks_dens2d  <- kde2d(ks_x, ks_y, n = 150,
                            h = c(ks_bw_x, ks_bw_y),
                            lims = c(ks_x_lim, ks_y_lim))
        ks_dens_df <- expand.grid(Si = ks_dens2d$x, Ki = ks_dens2d$y)
        ks_dens_df$density <- as.vector(ks_dens2d$z)

        ks_bw_si  <- density(ks_x, bw = "nrd0")
        ks_bw_ki  <- density(ks_y, bw = "nrd0")
        si_dens   <- data.frame(val = ks_bw_si$x, pdf = ks_bw_si$y)
        ki_dens   <- data.frame(val = ks_bw_ki$x, pdf = ks_bw_ki$y)
        si_dens   <- si_dens[si_dens$val >= ks_x_lim[1] & si_dens$val <= ks_x_lim[2], ]
        ki_dens   <- ki_dens[ki_dens$val >= ks_y_lim[1] & ki_dens$val <= ks_y_lim[2], ]

        ks_dens_max <- max(ks_dens_df$density)
        jet_cols_ks <- colorRampPalette(
            c("white", "white", "#E8E8FF", "#C0B0FF",
              "#7B2FBE", "#4B0082", "purple",
              "cyan", "#00CC88", "yellow"))(512)

        p_ks_main <- ggplot(ks_dens_df, aes(x = Si, y = Ki)) +
            geom_raster(aes(fill = density), interpolate = TRUE, na.rm = TRUE) +
            geom_contour(aes(z = density), colour = "black",
                         alpha = 0.45, linewidth = 0.45, bins = 10, na.rm = TRUE) +
            annotate("label", x = ks_x_lim[2] - diff(ks_x_lim) * 0.03,
                     y = ks_y_lim[1] + diff(ks_y_lim) * 0.12,
                     label = ann_ks_expr, hjust = 1, vjust = 0,
                     size = 3.0, fontface = "bold", fill = "white", alpha = 0.80,
                     label.padding = unit(0.3, "lines"), colour = "gray15") +
            scale_fill_gradientn(colours = jet_cols_ks, limits = c(0, ks_dens_max),
                oob = scales::squish, name = "Density", na.value = "white") +
            scale_x_continuous(limits = ks_x_lim, expand = c(0, 0)) +
            scale_y_continuous(limits = ks_y_lim, expand = c(0, 0)) +
            labs(x = "S_i_norm = 4(u-0.5)(v-0.5), range [-1,+1]  (same scale as Part 10b spear_norm; SKCS gate at 0.05)",
                 y = "K_i  (Kendall concordance score)") +
            theme_copula() +
            theme(legend.position = "none", panel.grid.major = element_line(color = "gray85"),
                  panel.grid.minor = element_blank())

        p_ks_top <- ggplot(si_dens, aes(x = val, y = pdf)) +
            geom_area(fill = "#5555BB", alpha = 0.30) +
            geom_line(colour = "#6666CC", linewidth = 1.0) +
            annotate("text", x = ks_x_lim[2] - diff(ks_x_lim) * 0.02,
                     y = max(si_dens$pdf) * 0.85,
                     label = "p(Si_norm) - marginal PDF of S_i_norm",
                     hjust = 1, size = 2.8, colour = "gray30") +
            scale_x_continuous(limits = ks_x_lim, expand = c(0, 0)) +
            scale_y_continuous(expand = c(0.02, 0)) +
            labs(x = NULL, y = "Density") +
            theme_copula() +
            theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                  panel.grid.minor = element_blank())

        p_ks_right <- ggplot(ki_dens, aes(x = pdf, y = val)) +
            geom_ribbon(aes(xmin = 0, xmax = pdf), fill = "#009988", alpha = 0.25) +
            geom_path(colour = "#00BBAA", linewidth = 1.0) +
            annotate("text", x = max(ki_dens$pdf) * 0.50,
                     y = ks_y_lim[1] + diff(ks_y_lim) * 0.92,
                     label = "p(Ki)\nK_i PDF", hjust = 0, size = 2.6, colour = "gray30") +
            scale_y_continuous(limits = ks_y_lim, expand = c(0, 0)) +
            scale_x_continuous(expand = c(0.02, 0)) +
            labs(x = "Density", y = NULL) +
            theme_copula() +
            theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                  panel.grid.minor = element_blank())

        left_col_6a <- plot_grid(p_ks_top, p_ks_main, ncol = 1, align = "v",
                                 axis = "lr", rel_heights = c(0.5, 3))
        right_col_6a <- plot_grid(ggplot() + theme_void(), p_ks_right, ncol = 1,
                                  align = "v", axis = "l", rel_heights = c(0.5, 3))
        g6a_inner <- plot_grid(left_col_6a, right_col_6a, ncol = 2,
                               rel_widths = c(3, 0.5))
        title_6a <- paste0(dataset_tag,
            "  -  Ki-Si joint density with marginal PDFs  [EXPRESSION SPACE]\n",
            "Main: f(S_i_norm, K_i) | Top: p(S_i_norm) | Right: p(K_i)  |  ",
            "All GOI-positive cells (n=", n_ks, ")")
        g6a <- plot_grid(
            ggdraw() + draw_label(title_6a, fontface = "bold", size = 10,
                                  x = 0.5, hjust = 0.5),
            g6a_inner, ncol = 1, rel_heights = c(0.08, 1))

        base6a <- file.path(OUT_DIR, paste0(FILE_PREFIX, "_", dataset_tag,
                                            "_06a_KiSi_marginal_joint"))
        save_pp(g6a, base6a, pdf_w = 10, pdf_h = 10,
                png_w = 2000, png_h = 2000, res = 200)
        cat("  Plot 06a saved (Ki-Si marginal joint)\n")
        all_plots[["KiSi_marginal_joint"]] <- g6a


        # =============================================================================
        # PLOT 06c: Copula unit square [0,1]x[0,1]  (Ki-Si)
        # Non-SKCS in gray80 (bottom), SKCS by IP x TREP group (top)
        # =============================================================================
        cat("  Building Plot 06c: Ki-Si copula unit square ...\n")

        # Si_norm = 0.05 and Ki = 0.5 gates in copula coordinates, k/(n+1)
        n_le_Si005  <- sum(ks_x <= 0.05)
        n_le_Ki05   <- sum(ks_y <= 0.5)
        u_Si005_cop <- n_le_Si005 / (n_ks + 1)  # F(S_i_norm) at the 0.05 gate
        v_Ki05_cop  <- n_le_Ki05  / (n_ks + 1)  # F(K_i) at the 0.5 gate

        gate_conv_line <- sprintf(
            paste0("Gate conversion: S_i_norm > 0.05 -> u = %d/%d = %.4f",
                   "   |   K_i > 0.5 -> v = %d/%d = %.4f   (k/(n+1))"),
            n_le_Si005, n_ks + 1L, u_Si005_cop,
            n_le_Ki05,  n_ks + 1L, v_Ki05_cop)
        cat("  ", gate_conv_line, "\n", sep = "")

        p_ks_copula <- ggplot(ks_copula_df, aes(x = u_sk, y = v_sk, colour = group)) +
            geom_point(data = subset(ks_copula_df, group == "Non-SKCS"),
                       alpha = 0.30, size = 0.4) +
            geom_point(data = subset(ks_copula_df, group != "Non-SKCS"),
                       alpha = 0.30, size = 0.5) +
            # Fixed threshold lines (Si_norm=0.05, Ki=0.5) in copula space
            geom_vline(xintercept = u_Si005_cop, colour = "#009E73",
                       linetype = "solid", linewidth = 0.6) +
            geom_hline(yintercept = v_Ki05_cop, colour = "#009E73",
                       linetype = "solid", linewidth = 0.6) +
            scale_colour_manual(values = cols_ks, labels = grp_labs_ks,
                                name = "Cell group") +
            scale_x_continuous(breaks = seq(0, 1, 0.25),
                               labels = c("0\n(lowest S_i)", "0.25",
                                          "0.5\n(median)", "0.75",
                                          "1\n(highest S_i)")) +
            scale_y_continuous(breaks = seq(0, 1, 0.25),
                               labels = c("0\n(lowest K_i)", "0.25",
                                          "0.5\n(median)", "0.75",
                                          "1\n(highest K_i)")) +
            labs(title = paste0(dataset_tag,
                                "  -  Ki-Si Copula unit square  [0,1] x [0,1]"),
                 subtitle = paste0(
                     "u_sk = F(Si_norm) = rank(S_i_norm) / (n+1)   |   ",
                     "v_sk = F(Ki) = rank(K_i) / (n+1)\n",
                     "Gray = Non-SKCS (K<0.5)  |  Coloured = SKCS-core (K>0.5)  |  ",
                     "All GOI-positive cells, n=", n_ks),
                 x = "u_sk = F(S_i_norm)  [S_i_norm marginal CDF]",
                 y = "v_sk = F(K_i)  [K_i marginal CDF]") +
            theme_copula() +
            theme(legend.position = "right", aspect.ratio = 1)

        save_pp(p_ks_copula, file.path(OUT_DIR, paste0(FILE_PREFIX, "_", dataset_tag,
                "_06c_KiSi_copula_unit_square")),
                pdf_w = 10, pdf_h = 9, png_w = 2000, png_h = 1800)
        cat("  Plot 06c saved (Ki-Si copula unit square)\n")
        all_plots[["KiSi_copula_unit_square"]] <- p_ks_copula


        # =============================================================================
        # PLOT 06b: Copula unit square with 3 demo cells + formula panel
        # =============================================================================
        cat("  Building Plot 06b: Ki-Si copula demo cells ...\n")

        pick_cell_ks <- function(target_u, df) {
            idx <- which.min(abs(df$u_sk - target_u))
            df[idx, , drop = FALSE]
        }
        cell_ks_A <- pick_cell_ks(0.20, ks_copula_df)
        # Cell B: middle of the Pre-IP TREP arm
        preip_trep_df <- subset(ks_copula_df, group == "Pre-IP TREP")
        preip_trep_median_u <- median(preip_trep_df$u_sk)
        cell_ks_B <- pick_cell_ks(preip_trep_median_u, preip_trep_df)
        cell_ks_C <- pick_cell_ks(0.80, ks_copula_df)
        demo_ks_cells <- rbind(cell_ks_A, cell_ks_B, cell_ks_C)
        demo_ks_cells$label <- c("A", "B", "C")

        rank_ks_x <- rank(ks_x, ties.method = "average")
        rank_ks_y <- rank(ks_y, ties.method = "average")
        np1_ks    <- n_ks + 1L

        demo_ks_info <- lapply(seq_len(3), function(i) {
            bc <- demo_ks_cells$barcode[i]
            j  <- which(bc_all == bc)[1]
            list(label = demo_ks_cells$label[i], si = ks_x[j], ki = ks_y[j],
                 rk_si = rank_ks_x[j], rk_ki = rank_ks_y[j],
                 ui = demo_ks_cells$u_sk[i], vi = demo_ks_cells$v_sk[i])
        })

        p_demo_ks_scatter <- ggplot(ks_copula_df,
                                    aes(x = u_sk, y = v_sk, colour = group)) +
            geom_point(alpha = 0.50, size = 0.9) +
            geom_point(data = demo_ks_cells, aes(x = u_sk, y = v_sk),
                       inherit.aes = FALSE, shape = 21, size = 5,
                       fill = "#8B5CF6", colour = "#4C1D95",
                       stroke = 1.2, alpha = 0.95) +
            geom_text(data = demo_ks_cells, aes(x = u_sk, y = v_sk, label = label),
                      inherit.aes = FALSE, fontface = "bold", size = 3.0,
                      colour = "white") +
            scale_colour_manual(values = cols_ks, labels = grp_labs_ks,
                                name = "Cell group") +
            scale_x_continuous(breaks = seq(0, 1, 0.25),
                               labels = c("0\n(lowest S_i)", "0.25",
                                          "0.5\n(median)", "0.75",
                                          "1\n(highest S_i)")) +
            scale_y_continuous(breaks = seq(0, 1, 0.25),
                               labels = c("0", "0.25", "0.5", "0.75", "1")) +
            labs(title = paste0(dataset_tag,
                                "  -  Ki-Si Copula unit square  [0,1] x [0,1]"),
                 subtitle = paste0(
                     "Three representative cells (A, B, C) highlighted\n",
                     "Formulas show how each cell's (u_sk, v_sk) is computed from rank / (n+1)\n",
                     gate_conv_line),
                 x = "u_sk = F(S_i_norm)  [S_i_norm marginal CDF]",
                 y = "v_sk = F(K_i)  [K_i marginal CDF]") +
            theme_copula() +
            theme(legend.position = "right", aspect.ratio = 1)

        formula_ks_lines <- c(
            "   How (u_sk, v_sk) is computed",
            "   ______________________________",
            sprintf("   n = %d cells   (n+1 = %d)", n_ks, np1_ks), "")
        for (ci in demo_ks_info) {
            formula_ks_lines <- c(formula_ks_lines,
                sprintf("   Cell %s", ci$label),
                sprintf("   S_i_norm = %.4f   rank = %.1f", ci$si, ci$rk_si),
                sprintf("   u_sk = %.1f / %d = %.4f", ci$rk_si, np1_ks, ci$ui), "",
                sprintf("   K_i  = %.4f   rank = %.1f", ci$ki, ci$rk_ki),
                sprintf("   v_sk = %.1f / %d = %.4f", ci$rk_ki, np1_ks, ci$vi),
                "   ______________________________", "")
        }
        p_formula_ks <- ggplot() +
            annotate("text", x = 0, y = 0.5,
                     label = paste(formula_ks_lines, collapse = "\n"),
                     hjust = 0, vjust = 0.5, family = "mono", size = 3.2,
                     colour = "gray15") +
            theme_void() +
            coord_cartesian(xlim = c(-0.05, 1), ylim = c(0, 1), clip = "off")

        p_demo_ks <- plot_grid(p_demo_ks_scatter, p_formula_ks,
                               ncol = 2, rel_widths = c(1.1, 0.6), align = "none")
        save_pp(p_demo_ks, file.path(OUT_DIR, paste0(FILE_PREFIX, "_", dataset_tag,
                "_06b_KiSi_copula_demo_cells")),
                pdf_w = 16, pdf_h = 9, png_w = 3200, png_h = 1800)
        cat("  Plot 06b saved (Ki-Si copula demo cells)\n")
        all_plots[["KiSi_copula_demo_cells"]] <- p_demo_ks


        # =============================================================================
        # PLOT 06d: Expression-space scatter (Si vs Ki) coloured by group
        # =============================================================================
        cat("  Building Plot 06d: Ki-Si expression-space scatter ...\n")

        ks_x_pad <- diff(range(ks_x)) * 0.04
        ks_y_pad <- diff(range(ks_y)) * 0.04
        ks_x_lim2 <- c(min(ks_x) - ks_x_pad, max(ks_x) + ks_x_pad)
        ks_y_lim2 <- c(max(-1, min(ks_y) - ks_y_pad), min(1, max(ks_y) + ks_y_pad))

        p_ks_expr <- ggplot(ks_copula_df, aes(x = Si_norm, y = Ki_raw, colour = group)) +
            geom_point(data = subset(ks_copula_df, group == "Non-SKCS"),
                       alpha = 0.30, size = 0.4) +
            geom_point(data = subset(ks_copula_df, group != "Non-SKCS"),
                       alpha = 0.30, size = 0.5) +
            # Fixed threshold lines (Si_norm=0.05, Ki=0.5) in expression space
            geom_vline(xintercept = 0.05, colour = "#009E73",
                       linetype = "solid", linewidth = 0.6) +
            geom_hline(yintercept = 0.5, colour = "#009E73",
                       linetype = "solid", linewidth = 0.6) +
            scale_colour_manual(values = cols_ks, labels = grp_labs_ks,
                                name = "Cell group") +
            scale_x_continuous(limits = ks_x_lim2, expand = c(0.05, 0)) +
            scale_y_continuous(limits = ks_y_lim2, expand = c(0.05, 0)) +
            labs(title = paste0(dataset_tag,
                                "  -  Expression space scatter  [S_i_norm vs K_i]"),
                 subtitle = paste0(
                     "Gray = Non-SKCS (K<0.5)  |  Coloured = SKCS-core  |  ",
                     "Same colouring as Plot 06c\n",
                     "K_i in [-1,+1]  |  S_i_norm in [-1,+1] (the scale the SKCS gate 0.05 is on)  |  n=", n_ks),
                 x = "S_i_norm = 4(u-0.5)(v-0.5), range [-1,+1]  (same scale as Part 10b spear_norm; SKCS gate at 0.05)",
                 y = "K_i  (Kendall concordance score)") +
            theme_copula() +
            theme(legend.position = "right", aspect.ratio = 1)

        save_pp(p_ks_expr, file.path(OUT_DIR, paste0(FILE_PREFIX, "_", dataset_tag,
                "_06d_KiSi_expression_space_scatter")),
                pdf_w = 10, pdf_h = 9, png_w = 2000, png_h = 1800)
        cat("  Plot 06d saved (Ki-Si expression-space scatter)\n")
        all_plots[["KiSi_expression_space_scatter"]] <- p_ks_expr


        # =============================================================================
        # PLOT 06e: Copula joint density with uniform marginal PDFs
        # =============================================================================
        cat("  Building Plot 06e: Ki-Si copula joint density ...\n")

        # Install goftest if missing; KS only when it is unavailable.
        if (!requireNamespace("goftest", quietly = TRUE)) {
            cat("  Installing 'goftest' for Anderson-Darling test ...\n")
            tryCatch(
                install.packages("goftest", quiet = TRUE,
                                 repos = "https://cloud.r-project.org"),
                error   = function(e)
                    cat("  goftest install failed:", conditionMessage(e), "\n"),
                warning = function(w)
                    cat("  goftest install warning:", conditionMessage(w), "\n"))
        }
        has_goftest_ks <- requireNamespace("goftest", quietly = TRUE)

        ks_test_u <- suppressWarnings(ks.test(u_sk, "punif", min = 0, max = 1))
        ks_test_v <- suppressWarnings(ks.test(v_sk, "punif", min = 0, max = 1))

        if (has_goftest_ks) {
            ad_u_sk <- goftest::ad.test(u_sk, "punif", min = 0, max = 1)
            ad_v_sk <- goftest::ad.test(v_sk, "punif", min = 0, max = 1)
            raw_p_ks <- c(ks_test_u$p.value, ks_test_v$p.value,
                          ad_u_sk$p.value, ad_v_sk$p.value)
            bh_p_ks  <- p.adjust(raw_p_ks, method = "BH")
        } else {
            ad_u_sk <- list(statistic = NA_real_, p.value = NA_real_)
            ad_v_sk <- list(statistic = NA_real_, p.value = NA_real_)
            raw_p_ks <- c(ks_test_u$p.value, ks_test_v$p.value, NA_real_, NA_real_)
            bh_p_ks  <- p.adjust(raw_p_ks[!is.na(raw_p_ks)], method = "BH")
            bh_p_ks  <- c(bh_p_ks, NA_real_, NA_real_)
        }
        cat(sprintf("  KS u_sk vs Uniform: D=%.5f p_BH=%.5f\n", ks_test_u$statistic, bh_p_ks[1]))
        cat(sprintf("  KS v_sk vs Uniform: D=%.5f p_BH=%.5f\n\n", ks_test_v$statistic, bh_p_ks[2]))

        COP_PAD_KS <- 0.20
        u_lim_6e <- c(0 - COP_PAD_KS, 1 + COP_PAD_KS)
        v_lim_6e <- c(0 - COP_PAD_KS, 1 + COP_PAD_KS)
        bw_u_sk <- bw.nrd0(u_sk); bw_v_sk <- bw.nrd0(v_sk)

        dens2d_ks_cop <- kde2d(u_sk, v_sk, n = 150, h = c(bw_u_sk, bw_v_sk),
                               lims = c(u_lim_6e, v_lim_6e))
        dens_ks_cop_df <- expand.grid(u_grid = dens2d_ks_cop$x, v_grid = dens2d_ks_cop$y)
        dens_ks_cop_df$density <- as.vector(dens2d_ks_cop$z)

        u_sk_dens <- density(u_sk, bw = "nrd0", from = u_lim_6e[1], to = u_lim_6e[2], n = 512)
        v_sk_dens <- density(v_sk, bw = "nrd0", from = v_lim_6e[1], to = v_lim_6e[2], n = 512)
        u_sk_dens_df <- data.frame(val = u_sk_dens$x, pdf = u_sk_dens$y)
        v_sk_dens_df <- data.frame(val = v_sk_dens$x, pdf = v_sk_dens$y)
        u_sk_dens_df <- u_sk_dens_df[u_sk_dens_df$val >= u_lim_6e[1] &
                                         u_sk_dens_df$val <= u_lim_6e[2], ]
        v_sk_dens_df <- v_sk_dens_df[v_sk_dens_df$val >= v_lim_6e[1] &
                                         v_sk_dens_df$val <= v_lim_6e[2], ]

        cop_ks_dens_max <- max(dens_ks_cop_df$density)
        jet_cols_ks_cop <- colorRampPalette(
            c("white", "white", "#E8E8FF", "#C0B0FF",
              "#7B2FBE", "#4B0082", "purple", "cyan", "#00CC88", "yellow"))(512)
        cop_ks_breaks <- c(0, 0.25, 0.50, 0.75, 1.00)

        p_ks_cop_main <- ggplot(dens_ks_cop_df, aes(x = u_grid, y = v_grid)) +
            geom_raster(aes(fill = density), interpolate = TRUE, na.rm = TRUE) +
            geom_contour(aes(z = density), colour = "black",
                         alpha = 0.50, linewidth = 0.45, bins = 10, na.rm = TRUE) +
            geom_vline(xintercept = c(0, 1), colour = "gray60",
                       linewidth = 0.35, linetype = "dotted") +
            geom_hline(yintercept = c(0, 1), colour = "gray60",
                       linewidth = 0.35, linetype = "dotted") +
            geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                        colour = "white", linewidth = 0.6, alpha = 0.7) +
            scale_fill_gradientn(colours = jet_cols_ks_cop,
                limits = c(0, cop_ks_dens_max), oob = scales::squish,
                name = "Density", na.value = "white") +
            scale_x_continuous(limits = u_lim_6e, expand = c(0, 0), breaks = cop_ks_breaks) +
            scale_y_continuous(limits = v_lim_6e, expand = c(0, 0), breaks = cop_ks_breaks) +
            annotate("label", x = 0.98, y = 0.03,
                     label = paste0(
                         sprintf("r_s = %.3f\nlambda_L = %.3f\nlambda_U = %.3f\n",
                                 rho_ks_orig, td_ks_orig$lambda_L, td_ks_orig$lambda_U),
                         sprintf("CvM: %s\n", cvm_ks_p_label),
                         sprintf("KS & AD: FDR = %.1f (uniform)", min(bh_p_ks, na.rm = TRUE))),
                     hjust = 1, vjust = 0, size = 2.8, fontface = "bold",
                     fill = alpha("black", 0.55), label.padding = unit(0.35, "lines"),
                     colour = "white") +
            labs(x = "u_sk = F(S_i_norm)  [S_i_norm marginal CDF]",
                 y = "v_sk = F(K_i)  [K_i marginal CDF]") +
            theme_copula() +
            theme(legend.position = "none",
                  panel.background = element_rect(fill = "gray98", colour = NA),
                  panel.grid.major = element_line(color = "gray80"),
                  panel.grid.minor = element_blank())

        u_sk_pdf_max <- max(u_sk_dens_df$pdf, 1.0) * 1.30
        ks_u_label_6e <- sprintf("KS p = %.4f (BH: %.4f)", raw_p_ks[1], bh_p_ks[1])
        p_ks_cop_top <- ggplot(u_sk_dens_df, aes(x = val, y = pdf)) +
            geom_hline(yintercept = 1.0, linetype = "dashed", colour = "gray50", linewidth = 0.7) +
            geom_area(fill = "#5555BB", alpha = 0.30) +
            geom_line(colour = "#6666CC", linewidth = 1.0) +
            annotate("text", x = u_lim_6e[2] - 0.01, y = u_sk_pdf_max * 0.90,
                     label = ks_u_label_6e, hjust = 1, size = 2.5, colour = "gray30") +
            scale_x_continuous(limits = u_lim_6e, expand = c(0, 0), breaks = cop_ks_breaks) +
            scale_y_continuous(limits = c(0, u_sk_pdf_max), expand = c(0.02, 0)) +
            labs(x = NULL, y = "Density") +
            theme_copula() +
            theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                  panel.grid.minor = element_blank())

        v_sk_pdf_max <- max(v_sk_dens_df$pdf, 1.0) * 1.30
        ks_v_label_6e <- sprintf("KS p = %.4f\n(BH: %.4f)", raw_p_ks[2], bh_p_ks[2])
        p_ks_cop_right <- ggplot(v_sk_dens_df, aes(x = pdf, y = val)) +
            geom_vline(xintercept = 1.0, linetype = "dashed", colour = "gray50", linewidth = 0.7) +
            geom_ribbon(aes(xmin = 0, xmax = pdf), fill = "#009988", alpha = 0.25) +
            geom_path(colour = "#00BBAA", linewidth = 1.0) +
            annotate("text", x = v_sk_pdf_max * 0.40, y = v_lim_6e[2] - 0.03,
                     label = ks_v_label_6e, hjust = 0, vjust = 1, size = 2.4, colour = "gray30") +
            scale_y_continuous(limits = v_lim_6e, expand = c(0, 0), breaks = cop_ks_breaks) +
            scale_x_continuous(limits = c(0, v_sk_pdf_max), expand = c(0.02, 0)) +
            labs(x = "Density", y = NULL) +
            theme_copula() +
            theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                  panel.grid.minor = element_blank())

        left_col_6e <- plot_grid(p_ks_cop_top, p_ks_cop_main, ncol = 1,
                                 align = "v", axis = "lr", rel_heights = c(0.5, 3))
        right_col_6e <- plot_grid(ggplot() + theme_void(), p_ks_cop_right, ncol = 1,
                                  align = "v", axis = "l", rel_heights = c(0.5, 3))
        g6e_inner <- plot_grid(left_col_6e, right_col_6e, ncol = 2, rel_widths = c(3, 0.5))

        ks_summary_6e <- sprintf(
            "KS u_sk: D=%.4f p=%s  |  KS v_sk: D=%.4f p=%s  (BH; p>0.05 = uniform)",
            ks_test_u$statistic,
            ifelse(bh_p_ks[1] < 0.001, "<0.001", sprintf("%.4f", bh_p_ks[1])),
            ks_test_v$statistic,
            ifelse(bh_p_ks[2] < 0.001, "<0.001", sprintf("%.4f", bh_p_ks[2])))
        title_6e <- paste0(dataset_tag,
            "  -  Ki-Si Copula joint density with marginal PDFs  [COPULA SPACE]\n",
            "Main: C(u_sk,v_sk) kernel density  |  Top: p(u_sk)  |  Right: p(v_sk)  |  ",
            "All GOI-positive cells (n=", n_ks, ")\n",
            "Dashed grey = Uniform(0,1) reference  |  Dashed white = concordance diagonal  |  ",
            "Axis +/- 0.20 padding\n", ks_summary_6e)
        g6e <- plot_grid(
            ggdraw() + draw_label(title_6e, fontface = "bold", size = 8.5,
                                  x = 0.5, hjust = 0.5),
            g6e_inner, ncol = 1, rel_heights = c(0.12, 1))

        base6e <- file.path(OUT_DIR, paste0(FILE_PREFIX, "_", dataset_tag,
                                            "_06e_KiSi_copula_joint_density"))
        save_pp(g6e, base6e, pdf_w = 10, pdf_h = 10,
                png_w = 2000, png_h = 2000, res = 200)
        cat("  Plot 06e saved (Ki-Si copula joint density)\n")
        all_plots[["KiSi_copula_joint_density"]] <- g6e


        # =============================================================================
        # PLOT 06f: Shuffled joint density with marginal PDFs
        # =============================================================================
        cat("  Building Plot 06f: Ki-Si shuffled marginal joint ...\n")

        # Same axis limits as Plot 06a (shuffled Ki has same range)
        ks_x_lim_s <- ks_x_lim
        ks_y_lim_s <- ks_y_lim

        dens2d_6f  <- kde2d(ks_x, ks_y_shuffled, n = 150,
                            h = c(bw.nrd0(ks_x) * 1.2, bw.nrd0(ks_y_shuffled) * 1.2),
                            lims = c(ks_x_lim_s, ks_y_lim_s))
        dens_df_6f <- expand.grid(Si = dens2d_6f$x, Ki = dens2d_6f$y)
        dens_df_6f$density <- as.vector(dens2d_6f$z)

        bw_si_6f <- density(ks_x, bw = "nrd0")
        bw_ki_6f <- density(ks_y_shuffled, bw = "nrd0")
        si_dens_6f <- data.frame(val = bw_si_6f$x, pdf = bw_si_6f$y)
        ki_dens_6f <- data.frame(val = bw_ki_6f$x, pdf = bw_ki_6f$y)
        si_dens_6f <- si_dens_6f[si_dens_6f$val >= ks_x_lim_s[1] &
                                     si_dens_6f$val <= ks_x_lim_s[2], ]
        ki_dens_6f <- ki_dens_6f[ki_dens_6f$val >= ks_y_lim_s[1] &
                                     ki_dens_6f$val <= ks_y_lim_s[2], ]

        dens_max_6f <- max(dens_df_6f$density, na.rm = TRUE)
        jet_cols_6f <- colorRampPalette(
            c("white", "white", "#E8E8FF", "#C0B0FF",
              "#7B2FBE", "#4B0082", "purple", "cyan", "#00CC88", "yellow"))(512)

        p_main_6f <- ggplot(dens_df_6f, aes(x = Si, y = Ki)) +
            geom_raster(aes(fill = density), interpolate = TRUE, na.rm = TRUE) +
            geom_contour(aes(z = density), colour = "black",
                         alpha = 0.45, linewidth = 0.45, bins = 10, na.rm = TRUE) +
            annotate("label", x = ks_x_lim_s[2] - diff(ks_x_lim_s) * 0.03,
                     y = ks_y_lim_s[1] + diff(ks_y_lim_s) * 0.12,
                     label = ann_ks_shuf, hjust = 1, vjust = 0,
                     size = 3.0, fontface = "bold", fill = "white", alpha = 0.80,
                     label.padding = unit(0.3, "lines"), colour = "gray15") +
            scale_fill_gradientn(colours = jet_cols_6f, limits = c(0, dens_max_6f),
                oob = scales::squish, name = "Density", na.value = "white") +
            scale_x_continuous(limits = ks_x_lim_s, expand = c(0, 0)) +
            scale_y_continuous(limits = ks_y_lim_s, expand = c(0, 0)) +
            labs(x = "S_i_norm = 4(u-0.5)(v-0.5), range [-1,+1]  (same scale as Part 10b spear_norm; SKCS gate at 0.05)",
                 y = "K_i  (Kendall concordance score)") +
            theme_copula() +
            theme(legend.position = "none", panel.grid.major = element_line(color = "gray85"),
                  panel.grid.minor = element_blank())

        p_top_6f <- ggplot(si_dens_6f, aes(x = val, y = pdf)) +
            geom_area(fill = "#5555BB", alpha = 0.30) +
            geom_line(colour = "#6666CC", linewidth = 1.0) +
            annotate("text", x = ks_x_lim_s[2] - diff(ks_x_lim_s) * 0.02,
                     y = max(si_dens_6f$pdf) * 0.85,
                     label = "p(Si_norm) - marginal PDF of S_i_norm",
                     hjust = 1, size = 2.8, colour = "gray30") +
            scale_x_continuous(limits = ks_x_lim_s, expand = c(0, 0)) +
            scale_y_continuous(expand = c(0.02, 0)) +
            labs(x = NULL, y = "Density") +
            theme_copula() +
            theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                  panel.grid.minor = element_blank())

        p_right_6f <- ggplot(ki_dens_6f, aes(x = pdf, y = val)) +
            geom_ribbon(aes(xmin = 0, xmax = pdf), fill = "#009988", alpha = 0.25) +
            geom_path(colour = "#00BBAA", linewidth = 1.0) +
            annotate("text", x = max(ki_dens_6f$pdf) * 0.50,
                     y = ks_y_lim_s[1] + diff(ks_y_lim_s) * 0.92,
                     label = "p(Ki)\nK_i PDF", hjust = 0, size = 2.6, colour = "gray30") +
            scale_y_continuous(limits = ks_y_lim_s, expand = c(0, 0)) +
            scale_x_continuous(expand = c(0.02, 0)) +
            labs(x = "Density", y = NULL) +
            theme_copula() +
            theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                  panel.grid.minor = element_blank())

        left_col_6f <- plot_grid(p_top_6f, p_main_6f, ncol = 1, align = "v",
                                 axis = "lr", rel_heights = c(0.5, 3))
        right_col_6f <- plot_grid(ggplot() + theme_void(), p_right_6f, ncol = 1,
                                  align = "v", axis = "l", rel_heights = c(0.5, 3))
        g6f_inner <- plot_grid(left_col_6f, right_col_6f, ncol = 2, rel_widths = c(3, 0.5))
        title_str_6f <- paste0(dataset_tag,
            "  -  SHUFFLED Ki-Si joint density with marginal PDFs  [EXPRESSION SPACE]\n",
            "Same marginal distributions as Plot 06a (identical p(Si_norm) and p(Ki))  |  ",
            "Dependence broken by random permutation of K_i\n",
            sprintf("rho=%.3f | lambda_L=%.3f | lambda_U=%.3f  |  n=%d",
                    rho_ks_shuf, td_ks_shuf$lambda_L, td_ks_shuf$lambda_U, n_ks))
        g6f <- plot_grid(
            ggdraw() + draw_label(title_str_6f, fontface = "bold", size = 10,
                                  x = 0.01, hjust = 0, vjust = 1),
            g6f_inner, ncol = 1, rel_heights = c(0.08, 1))

        save_pp(g6f, file.path(OUT_DIR, paste0(FILE_PREFIX, "_", dataset_tag,
                "_06f_KiSi_shuffled_marginal_joint")),
                pdf_w = 10, pdf_h = 9, png_w = 2000, png_h = 1800)
        cat("  Plot 06f saved (Ki-Si shuffled marginal joint)\n\n")
        all_plots[["KiSi_shuffled_marginal_joint"]] <- g6f

        # == Ki-Si copula summary to the workbook =============================
        excel_tabs[["Step6_KiSi_CopulaSummary"]] <- data.frame(
            dataset              = dataset_tag,
            n_cells_all          = n_ks,
            n_SKCS               = n_SKCS,
            n_nonSKCS            = n_ks - n_SKCS,
            Spearman_rho_KiSi    = round(rho_ks_orig, 4),
            Spearman_rho_shuffled = round(rho_ks_shuf, 4),
            lambda_L_q10         = round(td_ks_orig$lambda_L, 4),
            lambda_U_q10         = round(td_ks_orig$lambda_U, 4),
            lambda_L_shuf_q10    = round(td_ks_shuf$lambda_L, 4),
            lambda_U_shuf_q10    = round(td_ks_shuf$lambda_U, 4),
            CvM_copula_distance  = round(cvm_ks_dist, 6),
            CvM_perm_p           = round(cvm_ks_p, 4),
            KS_marginals_D       = round(ks_test_Ki$statistic, 4),
            KS_marginals_p       = round(ks_test_Ki$p.value, 4),
            KS_u_sk_D            = round(ks_test_u$statistic, 5),
            KS_u_sk_p_BH         = round(bh_p_ks[1], 5),
            KS_v_sk_D            = round(ks_test_v$statistic, 5),
            KS_v_sk_p_BH         = round(bh_p_ks[2], 5),
            stringsAsFactors = FALSE)

        excel_tabs[["Step6_KiSi_CellData"]] <- data.frame(
            barcode  = bc_all,
            group    = as.character(ks_group),
            is_SKCS  = bc_all %in% df_skcs$barcode,
            Si_norm     = round(ks_x, 6),          # 4(u-.5)(v-.5)
            Si_rankprod = round(Si_rankprod, 6),   # scale_factor form
            Ki       = round(ks_y, 6),
            u_sk     = round(u_sk, 5),
            v_sk     = round(v_sk, 5),
            stringsAsFactors = FALSE)

        cat("\n  All Part 10d quantities are cell level (no donor random intercept or weights).\n")

        cat("  Step 6 (Ki-Si copula) complete.\n\n")


        # == PRSS iteration traces to the workbook ============================
        excel_tabs[["GAM_gU_PRSS_iters"]] <- data.frame(
            dataset  = dataset_tag,
            iter     = seq_along(gam_image_result$prss_values),
            prss     = gam_image_result$prss_values,
            best_iter = seq_along(gam_image_result$prss_values) == gam_image_result$best_iter,
            stringsAsFactors = FALSE
        )
        excel_tabs[["GAM_Summary"]] <- data.frame(
            dataset           = dataset_tag,
            gam               = "g(u) image curve",
            best_k            = gam_image_result$best_k,
            best_iter         = gam_image_result$best_iter,
            best_prss         = round(gam_image_result$best_prss, 6),
            actual_iters      = gam_image_result$actual_iters,
            early_stopped     = gam_image_result$early_stopped,
            max_k_allowed     = gam_image_result$max_k,
            phase2_max        = gam_image_result$phase2_max,
            n_iter_max        = gam_image_result$n_iter_max,
            edf               = sum(summary(gam_image)$edf),
            dev_expl          = round(summary(gam_image)$dev.expl, 4),
            gamma             = 1.5,
            select            = TRUE,
            method            = "REML",
            bs                = "tp",
            stringsAsFactors  = FALSE
        )

        # == Export workbook ==================================================
        xlsx_path <- file.path(OUT_DIR,
                               paste0(FILE_PREFIX, "_", dataset_tag,
                                      "_copula_SKCS_analysis.xlsx"))
        write_xlsx(excel_tabs, path = xlsx_path)
        cat(sprintf("  Excel saved: %s\n", basename(xlsx_path)))

        # == Display all plots in RStudio =====================================
        cat("  Displaying plots in RStudio ...\n")
        for (nm in names(all_plots)) {
            cat("    ->", nm, "\n")
            p_obj <- all_plots[[nm]]
            if (inherits(p_obj, "ggplot")) {
                print(p_obj)
            } else {
                grid.newpage(); grid.draw(p_obj)
            }
        }

        invisible(list(
            dataset_tag       = dataset_tag,
            n                 = n,
            u = u, v = v, C_hat = C_hat,
            x = x, y = y,
            Si_obs = Si_obs, Ki_obs = Ki_obs, EP_obs = EP_obs,
            EP_obs_raw = EP_obs_raw,
            ip_stat = ip_stat, is_trep = is_trep, group = group,
            ip_axis_check = ip_axis_check,
            barcode = df_skcs$barcode,
            Si_formula = Si_formula, Ki_formula = Ki_formula,
            g_u_all = g_u_all, IP_u = IP_u, IP_v = IP_v,
            # u_IP and v_IP: alternative field names read by later parts.
            u_IP = IP_u, v_IP = IP_v,
            IP_thresh = IP_thresh, IP_y_thresh = IP_y_thresh,
            IP_available = IP_AVAILABLE, n_full = n_full,
            n_ep_unmatched = n_ep_unmatched,
            # Density waist test
            waist = list(
                mean_density_u025 = mean(d_at_025),
                mean_density_IP_u  = mean(d_at_IP_u),
                mean_density_u075 = mean(d_at_075),
                t_025 = if (is.null(ttest_025_IP_u)) NA_real_
                        else unname(ttest_025_IP_u$statistic),
                p_025 = if (is.null(ttest_025_IP_u)) NA_real_
                        else ttest_025_IP_u$p.value,
                t_075 = if (is.null(ttest_075_IP_u)) NA_real_
                        else unname(ttest_075_IP_u$statistic),
                p_075 = if (is.null(ttest_075_IP_u)) NA_real_
                        else ttest_075_IP_u$p.value,
                n_cells_025 = n_cells_025,
                n_cells_IP_u = n_cells_IP_u,
                n_cells_075 = n_cells_075,
                occupied    = waist_occupied,
                flank_min   = WAIST_FLANK_MIN_10D,
                unif_ref    = WAIST_UNIF_REF_10D),
            # v-direction waist test
            waist_v = list(
                mean_density_v025 = mean(d_at_v025),
                mean_density_IP_v  = mean(d_at_IP_v),
                mean_density_v075 = mean(d_at_v075),
                t_025 = if (is.null(ttest_025_IP_v)) NA_real_
                        else unname(ttest_025_IP_v$statistic),
                p_025 = if (is.null(ttest_025_IP_v)) NA_real_
                        else ttest_025_IP_v$p.value,
                t_075 = if (is.null(ttest_075_IP_v)) NA_real_
                        else unname(ttest_075_IP_v$statistic),
                p_075 = if (is.null(ttest_075_IP_v)) NA_real_
                        else ttest_075_IP_v$p.value,
                n_cells_025 = n_cells_v025,
                n_cells_IP_v = n_cells_IP_v,
                n_cells_075 = n_cells_v075,
                occupied    = waist_v_occupied,
                flank_min   = WAIST_FLANK_MIN_10D,
                unif_ref    = WAIST_UNIF_REF_10D),
            # Ki-Si copula
            u_sk = u_sk, v_sk = v_sk,
            rho_KiSi = rho_ks_orig,
            td_KiSi  = td_ks_orig,
            gam_image         = gam_image,
            gam_image_prss    = gam_image_result[c("best_k", "best_iter",
                                                    "best_prss", "actual_iters",
                                                    "early_stopped", "max_k",
                                                    "prss_values")],
            excel_tabs        = excel_tabs,
            plot_names        = names(all_plots)   # plot objects excluded from RDS (saved to disk)
        ))
    }


    # =============================================================================
    # SECTION C: EXECUTE
    # =============================================================================

    all_results_10d <- vector("list", length(DATASETS))
    names(all_results_10d) <- DATASETS

    for (ds in DATASETS) {
        # Single-bracket assignment keeps a NULL entry for a failed dataset.
        all_results_10d[ds] <- list(tryCatch(
            run_part10d(ds),
            error = function(e) {
                cat("  ERROR in dataset", ds, ":", conditionMessage(e), "\n")
                NULL
            }
        ))
    }

    # Saved through rds_save(), which strips model environments and plots.
    rds_path <- file.path(OUT_DIR, paste0(FILE_PREFIX, "_all_results.rds"))
    rds_save(all_results_10d, rds_path, label = "Results object")


    # =============================================================================
    # SECTION C1: POST-ANALYSIS GUARDRAILS
    # =============================================================================
    cat("\n", strrep("#", 72), "\n")
    cat("  SECTION C1 - POST-ANALYSIS GUARDRAILS\n")
    cat(strrep("#", 72), "\n\n")

    for (ds in names(all_results_10d)) {
        r <- all_results_10d[[ds]]
        if (is.null(r)) {
            gr10d(paste0("GR-06-", ds), "FAIL", "Dataset returned NULL")
            next
        }
        gr10d(paste0("GR-06-", ds),
              if (!is.null(r$n_full) && r$n_full >= MIN_CELLS) "PASS" else "FAIL",
              sprintf("GOI-positive cells = %s  (min = %d)",
                      if (is.null(r$n_full)) "NULL" else as.character(r$n_full),
                      MIN_CELLS))

        # == The IP this run used, and where it came from =====================
        gr10d(paste0("GR-07-", ds),
              if (isTRUE(r$IP_available)) "PASS" else "WARN",
              sprintf("IP_thresh = %s  (%s)",
                      if (is.null(r$IP_thresh) || !is.finite(r$IP_thresh))
                          "NA (IPRS not available; Step 4 centred on u = 0.5)"
                      else sprintf("%.6f", r$IP_thresh),
                      ip_from))

        # == The IP as a copula coordinate pair ===============================
        .u <- r$IP_u; .v <- r$IP_v
        gr10d(paste0("GR-08-", ds),
              if (is.null(.u) || !is.finite(.u)) "WARN"
              else if (is.null(.v) || !is.finite(.v)) "WARN" else "PASS",
              sprintf("IP in copula space: (IP_u, IP_v) = (%s, %s)",
                      if (is.null(.u) || !is.finite(.u)) "NA"
                      else sprintf("%.4f", .u),
                      if (is.null(.v) || !is.finite(.v))
                          "NA (IP_y absent from the IPRS RDS)"
                      else sprintf("%.4f", .v)))

        # == Which axis Part 10b's Pre/Post labels were made on ===============
        .ax <- r$ip_axis_check
        if (is.null(.ax) || is.null(.ax$n) || .ax$n == 0L ||
            !is.finite(.ax$agree_cdfg)) {
            gr10d(paste0("GR-10-", ds), "WARN",
                  paste0("Pre/Post labels not comparable (no IP_status in the ",
                         "EP CSV, or no IP_y from the IPRS RDS)"))
        } else {
            gr10d(paste0("GR-10-", ds),
                  if (.ax$agree_cdfg >= 0.999) "PASS" else "FAIL",
                  sprintf(paste0("Part 10b Pre/Post labels reproduce %.4f on ",
                                 "the cDFG axis at IP_y and %s on the GOI axis ",
                                 "at IP_x, over %d cell(s)%s"),
                          .ax$agree_cdfg,
                          if (is.finite(.ax$agree_goi))
                              sprintf("%.4f", .ax$agree_goi) else "NA",
                          .ax$n,
                          if (.ax$agree_cdfg < 0.999)
                              " -- the EP CSV was not written on the cDFG axis"
                          else ""))
        }

        # == Density waist test at the IP =====================================
        # Verdict: waist (strip occupied, density below both flanks, both tests
        # significant), gap (near-empty strip between two dense flanks), or none.
        # waist and gap PASS; none WARNs.
        .waist_verdict <- function(w, mid, f1, f2, n_mid, n_f1, n_f2) {
            .n1 <- function(z) is.numeric(z) && length(z) == 1L && is.finite(z)
            if (is.null(w)) return(list(kind = "none", status = "WARN",
                                        why = "test not available"))
            if (!.n1(w$p_025) || !.n1(w$p_075))
                return(list(kind = "none", status = "WARN",
                            why = "one or both tests could not be computed"))
            lower <- .n1(mid) && .n1(f1) && .n1(f2) && mid < f1 && mid < f2
            sig   <- w$p_025 < 0.05 && w$p_075 < 0.05
            if (!lower || !sig)
                return(list(kind = "none", status = "WARN",
                            why = if (!lower)
                                      "density at the IP is not below both flanks"
                                  else "one or both tests are not significant"))
            fmin <- if (.n1(w$flank_min)) w$flank_min else 20L
            uref <- if (.n1(w$unif_ref))  w$unif_ref  else 0
            if (isTRUE(w$occupied))
                return(list(kind = "waist", status = "PASS",
                            why = "IP strip occupied; density below both flanks"))
            flanks_dense <- .n1(n_f1) && .n1(n_f2) &&
                            n_f1 >= fmin && n_f2 >= fmin &&
                            f1 > uref && f2 > uref
            if (flanks_dense)
                return(list(kind = "gap", status = "PASS",
                            why = sprintf(paste0("IP strip holds %d cell(s) ",
                                                 "while both flanks exceed %d ",
                                                 "cells and the uniform-copula ",
                                                 "density %.3f: a GAP between ",
                                                 "two dense lobes, stronger ",
                                                 "than a waist"),
                                          as.integer(n_mid), as.integer(fmin),
                                          uref)))
            list(kind = "none", status = "WARN",
                 why = sprintf(paste0("IP strip holds %d cell(s) and the ",
                                      "flanks are not dense enough to call ",
                                      "this a gap (%d and %d cells, densities ",
                                      "%.3f and %.3f, against %d cells and ",
                                      "%.3f): the region is EMPTY rather than ",
                                      "a minimum"),
                              as.integer(n_mid), as.integer(n_f1),
                              as.integer(n_f2), f1, f2,
                              as.integer(fmin), uref))
        }

        .ratio_fmt <- function(flank, waist) {
            if (!is.finite(flank) || !is.finite(waist)) return("NA")
            if (waist <= 0) return("Inf (waist density is 0)")
            r <- flank / waist
            if (r >= 1000) sprintf("%.3g", r) else sprintf("%.2f", r)
        }
        .w <- r$waist
        # EP-CSV to GAM-data barcode coverage
        .nu <- r$n_ep_unmatched
        gr10d(paste0("GR-11-", ds),
              if (is.null(.nu)) "WARN" else if (.nu == 0L) "PASS" else "FAIL",
              if (is.null(.nu))
                  "EP-CSV barcode coverage not recorded for this dataset"
              else if (.nu == 0L)
                  "every Part 10b EP/TREP barcode has a row in the GAM data"
              else
                  sprintf(paste0("%d Part 10b EP/TREP barcode(s) have no row in ",
                                 "the GAM data and were dropped."), .nu))

        # v-direction waist at IP_v (WARN when not confirmed)
        .wv <- r$waist_v
        .vv <- if (is.null(.wv)) list(kind = "none", status = "WARN",
                                      why = "test not available")
               else .waist_verdict(.wv, .wv$mean_density_IP_v,
                                   .wv$mean_density_v025, .wv$mean_density_v075,
                                   .wv$n_cells_IP_v, .wv$n_cells_025,
                                   .wv$n_cells_075)
        gr10d(paste0("GR-12-", ds), .vv$status,
              if (is.null(.wv))
                  "v-direction density waist test not available"
              else sprintf(paste0("v-direction (cDFG marginal, supplementary) ",
                                  "-- %s: density at IP_v = %.4f vs %.4f ",
                                  "(v=0.25) and %.4f (v=0.75) | ratios %s, %s | ",
                                  "t = %.2f, %.2f | p = %.3g, %.3g | ",
                                  "cells in strip %d/%d/%d | %s"),
                           toupper(.vv$kind),
                           .wv$mean_density_IP_v, .wv$mean_density_v025,
                           .wv$mean_density_v075,
                           .ratio_fmt(.wv$mean_density_v025, .wv$mean_density_IP_v),
                           .ratio_fmt(.wv$mean_density_v075, .wv$mean_density_IP_v),
                           .wv$t_025, .wv$t_075, .wv$p_025, .wv$p_075,
                           .wv$n_cells_025, .wv$n_cells_IP_v, .wv$n_cells_075,
                           .vv$why))

        # u-direction waist at IP_u, same verdict rule
        .vu <- if (is.null(.w)) list(kind = "none", status = "WARN",
                                     why = "test not available")
               else .waist_verdict(.w, .w$mean_density_IP_u,
                                   .w$mean_density_u025, .w$mean_density_u075,
                                   .w$n_cells_IP_u, .w$n_cells_025,
                                   .w$n_cells_075)
        gr10d(paste0("GR-09-", ds), .vu$status,
              if (is.null(.w)) "density waist test not available"
              else sprintf(paste0("u-direction (GOI marginal, primary) -- %s: ",
                                  "density at IP_u = %.4f vs %.4f (u=0.25) and ",
                                  "%.4f (u=0.75) | ratios %s, %s | ",
                                  "t = %.2f, %.2f | p = %.3g, %.3g",
                                  " | cells in strip %d/%d/%d | %s"),
                           toupper(.vu$kind),
                           .w$mean_density_IP_u, .w$mean_density_u025,
                           .w$mean_density_u075,
                           .ratio_fmt(.w$mean_density_u025, .w$mean_density_IP_u),
                           .ratio_fmt(.w$mean_density_u075, .w$mean_density_IP_u),
                           .w$t_025, .w$t_075,
                           .w$p_025, .w$p_075,
                           .w$n_cells_025, .w$n_cells_IP_u, .w$n_cells_075,
                           .vu$why))
    }

    gr10d_write()
    cat(sprintf("\n  Final guardrail tally: %d PASS / %d WARN / %d FAIL\n",
                .gr10d_pass, .gr10d_warn, .gr10d_fail))
    if (.gr10d_fail > 0) {
        cat("  Failures present -- review the ledger before interpreting.\n")
    } else if (.gr10d_warn > 0) {
        cat(sprintf("  %d passed, %d warning(s), no failures.\n",
                    .gr10d_pass, .gr10d_warn))
    } else {
        cat(sprintf("  All %d checks passed.\n", .gr10d_pass))
    }


    # =============================================================================
    # SECTION C2: ARM SUMMARY
    # =============================================================================
    cat("\n", strrep("=", 72), "\n")
    cat("  Part10d Copula2 complete [", ARM_TAG, "]\n")
    cat("  Output directory:", OUT_DIR, "\n")
    cat("  Files generated:\n")
    # The marker and the message temporary file are removed below.
    .transient <- c("_RUN_INCOMPLETE.txt", basename(MSG_PATH_P10D))
    for (f in setdiff(sort(list.files(OUT_DIR)), .transient)) {
        sz <- round(file.size(file.path(OUT_DIR, f)) / 1e3, 0)
        cat(sprintf("    %-62s  %5d KB\n", f, sz))
    }
    cat("  Load the results with: readRDS(\"", rds_path, "\")\n", sep = "")
    cat(strrep("=", 72), "\n")

    cat("\n  Log written to:", LOG_PATH_P10D,
        if (isTRUE(LOG_MESSAGES_P10D)) " (console output and warnings)"
        else " (console output only)", "\n")
    cat(sprintf("\nPart10d Copula2 - log ended: %s\n", Sys.time()))

    .p10d_completed <- TRUE

    if (sink.number(type = "message") != 2L) sink(type = "message")
    while (sink.number() > 0) sink()

    # The arm completed, so the marker is removed.
    suppressWarnings(file.remove(INCOMPLETE_MARKER_P10D))

    invisible(list(arm  = arm_tag, dir = OUT_DIR,
                   pass = .gr10d_pass, warn = .gr10d_warn,
                   fail = .gr10d_fail))
}


# =============================================================================
# SECTION 0c: Driver
# =============================================================================

cep_run_arms(PART_DIR_TAG, run_one_arm)
