# =============================================================================
# Part10b_CP10K_SKCS_TREP_code.r
# SEA-AD CEP-IP/CMRF, Part 10b, SKCS fixed thresholds and IP x TREP classes
#   (Table S7)
#   For each arm and dataset (NonNeoC, NeoC): stratifies the non-SKCS cells,
#   fits the donor-aware two-phase PRSS-REML GAM at the fixed SKCS thresholds
#   (S_norm > 0.05, K > 0.5) and at the 50th percentiles, ranks SKCS-core cells
#   by explanatory power (EP) to define TREP cells, and classifies each cell as
#   Pre- or Post-IP on the cDFG axis at IP_y (IP from Part 10b IPRS).
#
# Usage (from the SEA-AD folder, or set SEAAD_ROOT), after Part 10 and
#   Part10b_CP10K_IPRS_code.r
#   source("Part10b_CP10K_SKCS_TREP_code.r")                         # all arms
#   CEP_ARM <- "PVM1_all"; source("Part10b_CP10K_SKCS_TREP_code.r")  # one arm
#   Needs CEP_arms.R and R_config.R (or R_config_HKG8.R) in the same folder.
#
# Outputs (Part10b_SKCS/Part10b_SKCS_<arm>/): EP/TREP values CSV, non-SKCS
#   stratification CSV, diagnostics workbook, cross-dataset summary CSV,
#   figures (PDF, PNG), 3D scatter (HTML), guardrail log, results RDS and log.
# =============================================================================


# =============================================================================
# SECTION 0: Libraries and shared arm definitions
# =============================================================================

library(mgcv)
library(Matrix)       # readMM for norm.mtx
library(ggplot2)
library(dplyr)
library(writexl)
library(scales)
library(gridExtra)
library(MASS)         # kde2d

set.seed(42)

# CEP_arms.R holds the arm table, the folder layout and the driver.
if (!exists("CEP_HOME")) CEP_HOME <- Sys.getenv("SEAAD_ROOT", getwd())
.p10b_harness <- file.path(CEP_HOME, "CEP_arms.R")
if (!file.exists(.p10b_harness))
    stop("CEP_arms.R not found at: ", .p10b_harness,
         "\n  It holds the arm table shared by Parts 9 onwards.",
         call. = FALSE)
source(.p10b_harness)

# PART_TAG names the output files; PART_DIR_TAG names the results folder.
PART_TAG     <- "Part10b"        # file-name tag
PART_DIR_TAG <- "Part10b_SKCS"   # results-folder tag

# Self-tests (bodies in Section Z); neither needs external files.
RUN_SMOKE_TEST  <- FALSE
RUN_STRESS_TEST <- FALSE


# =============================================================================
# SECTION 0d: Pure model helpers and PSOCK workers
# =============================================================================
# Defined at the top level so that PSOCK workers receive them without the arm's
# data; the configuration constants are passed in as arguments.

DONOR_RE_LABEL <- "s(donor)"

# GAM formula, with the donor random intercept when use_re is TRUE.
p10b_gam_formula <- function(k, use_re, bs) {
    rhs <- sprintf("GOI + s(GOI, bs = \"%s\", k = %d)", bs, as.integer(k))
    if (use_re) rhs <- paste0(rhs, " + s(donor, bs = \"re\")")
    stats::as.formula(paste("Expression ~", rhs))
}

# Index of the s(GOI) smooth, located by label.
p10b_goi_smooth_idx <- function(model) {
    if (is.null(model$smooth) || length(model$smooth) == 0)
        return(NA_integer_)
    labs <- vapply(model$smooth, function(s) s$label, character(1))
    i    <- which(labs == "s(GOI)")
    if (length(i) == 1L) return(i)
    i <- which(labs != DONOR_RE_LABEL)
    if (length(i) >= 1L) return(i[1])
    NA_integer_
}

# Smoothing parameter of one smooth term.
p10b_smooth_lambda <- function(model, i) {
    if (is.na(i) || is.null(model$sp) || length(model$sp) == 0)
        return(NA_real_)
    sm <- model$smooth[[i]]
    if (!is.null(sm$first.sp) && sm$first.sp >= 1L &&
        sm$first.sp <= length(model$sp))
        return(unname(model$sp[sm$first.sp]))
    hit <- which(startsWith(names(model$sp), sm$label))
    if (length(hit) > 0) return(unname(model$sp[hit[1]]))
    unname(model$sp[1])
}

# One column of summary()$s.table for one term, by row name.
p10b_smooth_stat <- function(sm_summary, label, column) {
    st <- sm_summary$s.table
    if (is.null(st) || nrow(st) == 0) return(NA_real_)
    i <- match(label, rownames(st))
    if (is.na(i)) i <- 1L
    if (!(column %in% colnames(st))) return(NA_real_)
    unname(st[i, column])
}

p10b_model_has_re <- function(model)
    !is.null(model$smooth) &&
        any(vapply(model$smooth, function(s) s$label, character(1)) ==
                DONOR_RE_LABEL)

# PRSS of the fitted model, from the s(GOI) penalty and coefficients.
p10b_calculate_prss <- function(model) {
    RSS <- sum(stats::residuals(model)^2)
    i   <- p10b_goi_smooth_idx(model)
    if (is.na(i))
        return(list(PRSS = RSS, RSS = RSS,
                    f_double_prime_integral = NA_real_, lambda = NA_real_))
    sm     <- model$smooth[[i]]
    lambda <- p10b_smooth_lambda(model, i)
    S      <- sm$S[[1]]
    beta   <- stats::coef(model)[sm$first.para:sm$last.para]
    f_dpi  <- as.numeric(t(beta) %*% S %*% beta)
    list(PRSS = RSS + lambda * f_dpi, RSS = RSS,
         f_double_prime_integral = f_dpi, lambda = lambda)
}

# == One PRSS iteration, run on a worker ======================================
# Reads the frame and constants from p10b_fit_ctx and returns the PRSS only;
# the master refits the winning iteration.
p10b_prss_worker <- function(job) {
    ctx <- p10b_fit_ctx
    res <- tryCatch({
        set.seed(123L + job$i)
        m <- mgcv::gam(p10b_gam_formula(job$k, ctx$use_re, ctx$bs),
                       data = ctx$data, weights = w,
                       method = ctx$method, select = ctx$select,
                       gamma = ctx$gamma)
        list(ok = TRUE, prss = p10b_calculate_prss(m)$PRSS)
    }, error = function(e) list(ok = FALSE, message = conditionMessage(e)))
    res$i <- job$i
    res$k <- job$k
    res
}


# =============================================================================
# SECTION 0e: RDS size control
# =============================================================================
# Fitted models carry the environments they were created in, which inflates
# the saved RDS; those are re-pointed at baseenv() (lossless). Plots are
# dropped unless RDS_KEEP_PLOTS is TRUE (they are already saved to disk).
RDS_KEEP_PLOTS <- FALSE

rds_strip_gam_env <- function(m) {
    if (!inherits(m, "gam")) return(m)
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

# Strips every fitted model and, unless keep_plots, replaces plots with NULL;
# names and classes are kept.
rds_slim <- function(x, keep_plots = RDS_KEEP_PLOTS) {
    if (inherits(x, "gam")) return(rds_strip_gam_env(x))
    if (!isTRUE(keep_plots) && .rds_is_plot(x)) return(NULL)
    if (is.list(x) && !is.data.frame(x)) {
        att <- attributes(x)
        out <- lapply(x, rds_slim, keep_plots = keep_plots)
        attributes(out) <- att
        return(out)
    }
    x
}

# Saves through rds_slim() and prints the file size.
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

    # Output printed before the log exists is buffered and replayed into it.
    .p10b_preamble <- character(0)

    assign("CEP_ARM", arm_tag, envir = globalenv())
    on.exit(suppressWarnings(rm("CEP_ARM", envir = globalenv())), add = TRUE)

    .p10b_preamble <- tryCatch(
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

    .p10b_needed <- c("ARM_TAG", "ARM_DIR", "R_EXPORT_ROOT", "GOI", "DFG",
                      "HKG_GENES", "DFG_COMPOSITE_LABEL", "RUN_SUFFIX",
                      "PRIMARY_STYPE", "NORM_TAG", "NORM_TARGET_SUM",
                      "NORM_CEILING")
    .p10b_absent <- .p10b_needed[!vapply(.p10b_needed, exists, logical(1))]
    if (length(.p10b_absent) > 0)
        stop("The configuration loaded from\n    ", config_path,
             "\n  does not define: ", paste(.p10b_absent, collapse = ", "),
             call. = FALSE)

    if (!identical(ARM_TAG, arm_tag))
        stop("The configuration reports ARM_TAG '", ARM_TAG,
             "' but this arm is '", arm_tag, "'.", call. = FALSE)

    .p10b_preamble <- c(.p10b_preamble, utils::capture.output(
        cat("  Config OK: ", ARM_TAG, " | ", GOI, " + ", length(DFG), " DFG | ",
            NORM_TAG, "\n", sep = "")))

    set.seed(42)

    # == User parameters ======================================================
    DATASETS <- c("NonNeoC", "NeoC")

    # GAM parameters, as in Parts 9 and 10
    MIN_CELLS            <- 50
    GAM_GAMMA_P10        <- 1.5
    PRSS_N_ITER_P10      <- 100L
    PRSS_EARLY_STOP_P10  <- 20L

    # == Estimator ============================================================
    # Donor random intercept; FALSE gives cell-level fits.
    GAM_DONOR_RE        <- TRUE
    DONOR_RE_MIN_LEVELS <- 3L
    # DONOR_RE_LABEL is defined at the top level (Section 0d) for the workers.

    # Inverse-donor-count weights; FALSE gives every cell w = 1.
    W_DONOR_WEIGHTS     <- TRUE

    # == Parallelism ==========================================================
    # Number of PSOCK workers; 1L runs everything in this process.
    .p10b_detected <- suppressWarnings(parallel::detectCores())
    PAR_CORES_P10B <- if (is.na(.p10b_detected)) 1L
                      else max(1L, as.integer(.p10b_detected) - 1L)
    .p10b_cl <- NULL

    # Phase-2 refinement cap: a fit runs at most 8 + PRSS_PHASE2_MAX_P10B
    # iterations, so PRSS_EARLY_STOP_P10 is not reached.
    PRSS_PHASE2_MAX_P10B <- 5L

    # == Part 10 output directory (the source for this part) ==================
    PART10_DIR <- cep_arm_dir_for("Part10", arm_tag)
    P10_PREFIX <- cep_out_prefix("Part10", NORM_TAG,
                                 cep_file_tag_for(arm_tag))
    P10_FILE_PREFIX <- paste0(P10_PREFIX, "_Module1_", gsub("-", "", GOI), "_",
                              DFG_COMPOSITE_LABEL)

    # Part 10 writes each dataset to <arm folder>/<arm folder name>_<dataset>;
    # the arm folder root is tried second.
    p10b_scores_candidates <- function(ds) {
        fn <- paste0(P10_FILE_PREFIX, "_", ds, "_Pooled_SKCS_scores.csv")
        c(file.path(PART10_DIR, paste0(basename(PART10_DIR), "_", ds), fn),
          file.path(PART10_DIR, fn))
    }
    p10b_scores_csv <- function(ds) {
        cand <- p10b_scores_candidates(ds)
        hit  <- cand[file.exists(cand)]
        if (length(hit) > 0L) hit[1] else cand[1]
    }

    # == Output directories ===================================================
    OUT_DIR <- out_dir
    dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

    OUT_PREFIX  <- cep_out_prefix(PART_TAG, NORM_TAG,
                                  cep_file_tag_for(arm_tag))
    FILE_PREFIX <- paste0(OUT_PREFIX, "_SKCS")

    # == Gene definitions =====================================================
    GOI_10 <- GOI
    DFG_10 <- DFG
    DFG_DISPLAY_LABEL <- DFG_COMPOSITE_LABEL

    # HKG comparison gene set (none when the GOI is itself an HKG gene)
    if (GOI %in% HKG_GENES) {
        HKG_10 <- NULL
        cat("  HKG comparison: SKIPPED (GOI is an HKG gene; DFG = ",
            DFG_COMPOSITE_LABEL, ")\n", sep = "")
    } else {
        HKG_10 <- HKG_GENES
    }

    # == Log file =============================================================
    # Printed output and the message stream both go to the log; LOG_MESSAGES_P10B
    # = FALSE keeps warnings on the console. Messages are appended at the end.
    LOG_MESSAGES_P10B <- TRUE
    LOG_PATH_P10B <- file.path(OUT_DIR, paste0(OUT_PREFIX,
                                               "_SKCS_console_log_",
                                               format(Sys.time(),
                                                      "%Y%m%d_%H%M%S"),
                                               ".txt"))
    MSG_PATH_P10B <- file.path(OUT_DIR,
                               paste0(OUT_PREFIX, "_SKCS_console_log_messages.tmp"))

    while (sink.number() > 0) sink()
    # sink.number(type = "message") is 2 when nothing is diverted.
    if (sink.number(type = "message") != 2L) sink(type = "message")

    .p10b_log_con <- file(LOG_PATH_P10B, open = "wt")
    .p10b_msg_con <- if (isTRUE(LOG_MESSAGES_P10B))
        file(MSG_PATH_P10B, open = "wt") else NULL
    sink(.p10b_log_con, append = FALSE, split = TRUE)
    if (!is.null(.p10b_msg_con)) sink(.p10b_msg_con, type = "message")
    on.exit({
        tryCatch({
            if (!isTRUE(.p10b_completed))
                cat("\n", strrep("!", 72),
                    "\n  ARM DID NOT COMPLETE. The error is printed to the console.\n",
                    strrep("!", 72), "\n", sep = "")
            if (sink.number(type = "message") != 2L) sink(type = "message")
            while (sink.number() > 0) sink()
            if (!is.null(.p10b_msg_con)) close(.p10b_msg_con)
            close(.p10b_log_con)
            # Append the message stream to the log and remove the temporary file.
            if (!is.null(.p10b_msg_con) && file.exists(MSG_PATH_P10B)) {
                .msg_txt <- readLines(MSG_PATH_P10B, warn = FALSE)
                if (length(.msg_txt) > 0)
                    cat(c("", strrep("=", 72),
                          "  MESSAGE STREAM (warnings, messages, errors)",
                          strrep("=", 72), .msg_txt, ""),
                        sep = "\n", file = LOG_PATH_P10B, append = TRUE)
                unlink(MSG_PATH_P10B)
            }
        }, error = function(e) invisible(NULL))
    }, add = TRUE)

    # Set at the end of the arm; the on.exit handler reports an incomplete run.
    .p10b_completed <- FALSE

    # Restore the caller's warn level on exit.
    .p10b_warn_opt <- getOption("warn")
    options(warn = 1)   # surface warnings immediately rather than deferring
    on.exit(options(warn = .p10b_warn_opt), add = TRUE)
    cat(sprintf("Part10b SKCS TREP - log started: %s\n\n",
                Sys.time()))

    # Replay what was printed before the log existed.
    if (length(.p10b_preamble) > 0)
        cat(paste(.p10b_preamble, collapse = "\n"), "\n\n", sep = "")

    # == DATASET_META: regions, GAM data folder and IP thresholds =============
    # IP values come from the IPRS RDS; IP-dependent analysis is skipped when NA.
    DATASET_META <- list(
        NonNeoC = list(regions     = NONNC_REGIONS,
                       gam_dir     = NONNC_DIR,
                       dfg_genes   = DFG_NONNC,
                       IP_thresh   = NA_real_,
                       IP_y_thresh = NA_real_),
        NeoC    = list(regions     = NEOC_REGIONS,
                       gam_dir     = NEOC_DIR,
                       dfg_genes   = DFG_NEOC,
                       IP_thresh   = NA_real_,
                       IP_y_thresh = NA_real_)
    )

    # Axis of the Pre/Post-IP split: "GOI" (GOI <= IP_thresh) or "cDFG"
    # (cDFG < IP_y_thresh = f(IP_thresh)). IP_thresh is the reported IP either way.
    IP_CLASS_AXIS <- "cDFG"
    if (!IP_CLASS_AXIS %in% c("GOI", "cDFG"))
        stop("IP_CLASS_AXIS must be \"GOI\" or \"cDFG\", got: ", IP_CLASS_AXIS,
             call. = FALSE)

    # == Load the IPRS IP values ==============================================
    # The IPRS RDS is searched in these locations, in order; the one used is printed.
    .iprs_candidates <- c(
        # 1. The copy the IPRS part writes into this part's folder.
        file.path(OUT_DIR, paste0(OUT_PREFIX, "_IPRS_results.rds")),
        # 2. The IPRS part's own output.
        file.path(cep_arm_dir_for("Part10b_IPRS", arm_tag),
                  paste0(cep_out_prefix("Part10b_IPRS", NORM_TAG,
                                        cep_file_tag_for(arm_tag)),
                         "_results.rds")),
        # 3. The older Part10b folder.
        file.path(cep_arm_dir_for("Part10b", arm_tag),
                  paste0(OUT_PREFIX, "_IPRS_results.rds")),
        # 4. The scran location (primary arm only).
        file.path(R_EXPORT_ROOT, "Part10b_IPRS",
                  "Part10b_scran_IPRS_results.rds"))

    # Label used when the RDS does not record the IP definition.
    ip_from <- "population unrecorded"

    .iprs_hit    <- .iprs_candidates[file.exists(.iprs_candidates)]
    iprs_rds_path <- if (length(.iprs_hit) > 0L) .iprs_hit[1]
                     else .iprs_candidates[1]
    if (length(.iprs_hit) > 0L && !identical(iprs_rds_path,
                                             .iprs_candidates[1]))
        cat("  [IPRS] Not in this arm's folder; using:", iprs_rds_path, "\n")

    if (file.exists(iprs_rds_path)) {
        iprs_rds <- readRDS(iprs_rds_path)
        # inflection_points_DFG holds the IP definition the IPRS part selected;
        # ip_definition names it.
        ip_vec  <- iprs_rds$inflection_points_DFG
        # Response-axis (cDFG) threshold, IP_y = f(IP_thresh).
        ip_y_vec <- iprs_rds$inflection_points_DFG_y
        ip_from <- if (!is.null(iprs_rds$ip_definition))
            switch(iprs_rds$ip_definition,
                   Si0_full  = "Si=0 median GOI, full population",
                   Si0_core  = "Si=0 median GOI, SKCS-core",
                   IPRS_full = "IPRS R(x)=0, full population",
                   IPRS_core = "IPRS R(x)=0, SKCS-core",
                   iprs_rds$ip_definition)
        else if (identical(iprs_rds$ip_population_for_part10b, "SKCS_core"))
            "SKCS-core (RDS predates ip_definition)"
        else "unrecorded definition (RDS predates ip_definition)"

        if (!is.null(ip_vec) && length(ip_vec) > 0L)
            cat("  [IPRS] IP definition:", ip_from, "\n")
        if (is.null(iprs_rds$ip_definition))
            warning("This IPRS RDS predates ip_definition, so which IP ",
                    "definition it carries cannot be confirmed from the file. ",
                    "Rerun Part10b_CP10K_IPRS_code.r.", call. = FALSE)
        if (is.null(ip_vec) || length(ip_vec) == 0L) {
            cat("  [IPRS] RDS has no usable inflection_points_DFG field:",
                iprs_rds_path, "\n")
            cat("  [IPRS] IP-dependent outputs will be skipped.\n")
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
                    cat(sprintf(paste0("  [IPRS] %-10s IP_y_thresh set to %.6f",
                                       "  = f(IP_thresh), classification axis\n"),
                                ds, ip_y_val))
                } else if (identical(IP_CLASS_AXIS, "cDFG")) {
                    # Without IP_y_thresh the dataset stops at the classification step.
                    cat(sprintf(paste0("  [IPRS] %-10s IP_y_thresh absent from ",
                                       "the IPRS RDS; this dataset will stop ",
                                       "at the classification step.\n"), ds))
                }
            }
        }
    } else {
        cat("  [IPRS] RDS not found. Looked in:\n")
        for (.p in .iprs_candidates) cat("      ", .p, "\n", sep = "")
        cat("  [IPRS] IP-dependent outputs will be skipped.\n")
    }

    cat("\n", strrep("=", 72), "\n")
    cat("  Part10b - SKCS Fixed-Threshold Diagnostic\n")
    cat("  Arm          :", ARM_TAG, "\n")
    cat("  GOI          :", GOI_10, "\n")
    cat("  DFG          :", DFG_COMPOSITE_LABEL, " - ",
        paste(DFG_10, collapse = ", "), "\n")
    cat("  Normalisation:", NORM_TAG, " log1p(UMIs/total x ",
        format(NORM_TARGET_SUM, scientific = FALSE), ")\n")
    cat("  Datasets     :", paste(DATASETS, collapse = ", "), "\n")
    cat("  S threshold  : > 0   (cell pulls rho upward)\n")
    cat("  K threshold  : > 0.5 (concordant with >50% of cell pairs)\n")
    cat("  Estimator    :",
        if (GAM_DONOR_RE)
            sprintf("donor random intercept, s(donor, bs = \"re\"), minimum %d donors",
                    DONOR_RE_MIN_LEVELS)
        else "cell level (donor random intercept disabled)", "\n")
    cat("  Weights      :",
        if (W_DONOR_WEIGHTS)
            "inverse donor cell count, w = 1/n_i (each donor sums to 1)"
        else "uniform, w = 1 (inverse-donor weighting disabled)", "\n")
    cat("  PRSS         : phase 1 = 8 iterations, phase 2 = up to",
        PRSS_PHASE2_MAX_P10B, "(early stop", PRSS_EARLY_STOP_P10, ")\n")
    cat("  Part 10 dir  :", PART10_DIR, "\n")
    cat("  Scores read  :",
        paste(paste0(basename(PART10_DIR), "_", DATASETS), collapse = ", "),
        "(per-dataset folders; arm root as fallback)\n")
    cat("  Output dir   :", OUT_DIR, "\n")

    # One PSOCK cluster per arm, stopped at the end of the arm or on abort.
    if (PAR_CORES_P10B > 1L) {
        .p10b_cl <- parallel::makeCluster(PAR_CORES_P10B)   # PSOCK: no fork
        on.exit({
            if (!is.null(.p10b_cl))
                try(parallel::stopCluster(.p10b_cl), silent = TRUE)
        }, add = TRUE)
        parallel::clusterEvalQ(.p10b_cl, { library(mgcv); NULL })
        parallel::clusterExport(.p10b_cl,
            c("DONOR_RE_LABEL", "p10b_gam_formula", "p10b_goi_smooth_idx",
              "p10b_smooth_lambda", "p10b_smooth_stat", "p10b_model_has_re",
              "p10b_calculate_prss", "p10b_prss_worker"))
        cat("  Parallel     :", PAR_CORES_P10B,
            "PSOCK worker(s) for the PRSS iterations\n")
    } else {
        cat("  Parallel     : disabled (PAR_CORES_P10B = 1); sequential\n")
    }
    cat(strrep("=", 72), "\n\n")

    # Marker file written on entry and removed only on a clean finish.
    INCOMPLETE_MARKER_P10B <- file.path(OUT_DIR, "_RUN_INCOMPLETE.txt")
    writeLines(c(
        "This arm did not finish.",
        "",
        paste("Arm     :", arm_tag),
        paste("Part    :", PART_DIR_TAG),
        paste("Started :", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
        "",
        "Any files in this folder are from a partial run and must not be used."
    ), INCOMPLETE_MARKER_P10B)


    # =============================================================================
    # SECTION A0: PRE-ANALYSIS GUARDRAILS
    # =============================================================================
    # Input checks before any fitting: the CP10K ceiling and row sums of norm.mtx,
    # integer raw counts, and the Part 10 scores CSV for each dataset.

    GR_LOG_P10B <- file.path(OUT_DIR, paste0(OUT_PREFIX,
                                              "_SKCS_Guardrail_Log.txt"))

    .gr10b_log   <- list()
    .gr10b_pass  <- 0L
    .gr10b_fail  <- 0L
    .gr10b_warn  <- 0L

    gr10b <- function(check_id, status, detail) {
        # A zero-length or NA detail is replaced so the row is never lost.
        if (length(detail) != 1L || is.na(detail))
            detail <- paste0("[detail unavailable: the message was built from a",
                             " missing value and came out zero-length]")
        # Any other status is printed as given.
        tag <- switch(status, PASS = "[PASS]", FAIL = "[FAIL]", WARN = "[WARN]",
                      sprintf("[%s]", as.character(status)[1]))
        cat(sprintf("  %s  %-16s  %s\n", tag, check_id, detail))
        .gr10b_log[[length(.gr10b_log) + 1L]] <<- list(
            id = check_id, status = status, detail = detail)
        if (status == "PASS") .gr10b_pass <<- .gr10b_pass + 1L
        if (status == "FAIL") .gr10b_fail <<- .gr10b_fail + 1L
        if (status == "WARN") .gr10b_warn <<- .gr10b_warn + 1L
    }

    gr10b_write <- function() {
        con <- file(GR_LOG_P10B, "wt")
        # Close the connection even if a write fails.
        on.exit(try(close(con), silent = TRUE), add = TRUE)
        writeLines(c(
            strrep("=", 72),
            paste0("Part 10b SKCS - Guardrail Log [", NORM_TAG, "]"),
            paste("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
            paste("Arm :", ARM_TAG),
            paste0("Normalisation: adata.X = log1p(UMIs/total x ",
                   format(NORM_TARGET_SUM, scientific = FALSE),
                   "), exported unchanged"),
            strrep("=", 72), ""
        ), con)
        for (item in .gr10b_log)
            writeLines(sprintf("  [%s]  %-16s  %s", item$status, item$id,
                               item$detail), con)
        writeLines(c("", strrep("=", 72),
                     sprintf("  Total : %d", .gr10b_pass + .gr10b_fail + .gr10b_warn),
                     sprintf("  PASS  : %d", .gr10b_pass),
                     sprintf("  WARN  : %d", .gr10b_warn),
                     sprintf("  FAIL  : %d", .gr10b_fail),
                     if (.gr10b_fail == 0 && .gr10b_warn == 0)
                         "\n  Result: no failures, no warnings"
                     else if (.gr10b_fail == 0)
                         "\n  Result: no failures, warnings present"
                     else "\n  Result: failures present",
                     strrep("=", 72)), con)
        close(con)
        cat("  Guardrail log written:", GR_LOG_P10B, "\n")
    }

    cat("\n", strrep("=", 72), "\n")
    cat("  SECTION A0 - PRE-ANALYSIS GUARDRAILS\n")
    cat(strrep("=", 72), "\n\n")

    # == Normalisation checks, per dataset ====================================
    for (.ds in list(
        list(tag = "NeoC",    dir = NEOC_DIR,  dfg = DFG_NEOC),
        list(tag = "NonNeoC", dir = NONNC_DIR, dfg = DFG_NONNC)
    )) {
        cat("  -- Dataset:", .ds$tag, "--------------------------------------\n")

        norm_mtx_path  <- file.path(.ds$dir, "deg_support", "norm.mtx")
        raw_mtx_path   <- file.path(.ds$dir, "counts", "raw_selected",
                                    "matrix.mtx")

        if (file.exists(norm_mtx_path)) {
            X_n <- as(readMM(norm_mtx_path), "CsparseMatrix")
            n_max <- max(X_n)

            # == Maximum within the CP10K ceiling =============================
            gr10b(paste0("GR-01-", .ds$tag),
                  if (n_max > 0 && n_max <= NORM_CEILING + 1e-6) "PASS"
                  else "FAIL",
                  sprintf("norm.mtx max = %.4f  (ceiling log1p(%s) = %.4f)",
                          n_max,
                          format(NORM_TARGET_SUM, scientific = FALSE),
                          NORM_CEILING))

            # == expm1 row sums equal the target ==============================
            n_rows <- min(50L, nrow(X_n))
            rs     <- rowSums(expm1(X_n[seq_len(n_rows), , drop = FALSE]))
            rs_dev <- max(abs(rs - NORM_TARGET_SUM))
            gr10b(paste0("GR-02-", .ds$tag),
                  if (rs_dev < NORM_TARGET_SUM * 1e-6) "PASS" else "FAIL",
                  sprintf("expm1(norm) row sums: mean = %.6f, max |dev| = %.2e, n = %d",
                          mean(rs), rs_dev, n_rows))

            rm(X_n); gc(verbose = FALSE)
        } else {
            # Both checks fail when norm.mtx is absent.
            gr10b(paste0("GR-01-", .ds$tag), "FAIL",
                  paste("norm.mtx NOT FOUND:", norm_mtx_path))
            gr10b(paste0("GR-02-", .ds$tag), "FAIL",
                  paste("norm.mtx NOT FOUND, so the CP10K row-sum identity",
                        "could not be checked:", norm_mtx_path))
        }

        # == Raw counts are integer-valued ====================================
        if (file.exists(raw_mtx_path)) {
            X_r    <- as(readMM(raw_mtx_path), "CsparseMatrix")
            r_samp <- as.numeric(X_r@x[seq_len(min(10000, length(X_r@x)))])
            is_int <- all(abs(r_samp - round(r_samp)) < 0.01)
            r_max  <- max(X_r)
            gr10b(paste0("GR-04-", .ds$tag),
                  if (is_int && r_max > 20) "PASS" else "FAIL",
                  sprintf("raw max = %.0f, integer-like = %s",
                          r_max, if (is_int) "TRUE" else "FALSE"))
            rm(X_r); gc(verbose = FALSE)
        } else {
            gr10b(paste0("GR-04-", .ds$tag), "FAIL",
                  paste("raw_selected matrix NOT FOUND:", raw_mtx_path))
        }
    }

    # == Part 10 scores CSV exists ============================================
    for (ds in DATASETS) {
        scores_csv_gr <- p10b_scores_csv(ds)
        .found <- file.exists(scores_csv_gr)
        gr10b(paste0("GR-05-", ds),
              if (.found) "PASS" else "FAIL",
              sprintf("Scores CSV: %s",
                      if (.found)
                          file.path(basename(dirname(scores_csv_gr)),
                                    basename(scores_csv_gr))
                      else paste("NOT FOUND in either location:",
                                 paste(p10b_scores_candidates(ds),
                                       collapse = "  |  "))))
    }

    # Written now and again at the end, so an aborted arm still leaves a ledger.
    gr10b_write()

    cat(sprintf("\n  Pre-analysis guardrails: %d PASS / %d WARN / %d FAIL\n",
                .gr10b_pass, .gr10b_warn, .gr10b_fail))
    if (.gr10b_fail > 0) {
        cat("  *** GUARDRAIL FAILURES ***\n\n")
        warning("Part 10b pre-analysis guardrails: ", .gr10b_fail, " FAIL(s)",
                call. = FALSE)
    } else if (.gr10b_warn > 0) {
        cat(sprintf("  %d passed, %d warning(s), no failures.\n\n",
                    .gr10b_pass, .gr10b_warn))
    } else {
        cat("  *** ALL PRE-ANALYSIS GUARDRAILS PASSED ***\n\n")
    }


    # =============================================================================
    # SECTION A0b: DONOR RANDOM INTERCEPT, WEIGHTS, PREDICTION, DE
    # =============================================================================
    # Local copies of the Part 9 and Part 10 helpers (those files are not sourced).

    # Donor as a factor with unused levels dropped.
    as_donor_factor <- function(x) droplevels(factor(as.character(x)))

    # Accepts either a fitting frame or a bare donor vector.
    .donor_vec_of <- function(x) {
        if (is.null(x)) return(NULL)
        if (is.data.frame(x)) {
            if (!("donor" %in% names(x))) return(NULL)
            return(x$donor)
        }
        x
    }

    n_donor_levels <- function(x) {
        v <- .donor_vec_of(x)
        if (is.null(v)) return(0L)
        nlevels(as_donor_factor(v))
    }

    # The random intercept needs at least DONOR_RE_MIN_LEVELS donors; otherwise
    # the cell-level formula is used.
    donor_re_ok <- function(x)
        isTRUE(GAM_DONOR_RE) && n_donor_levels(x) >= DONOR_RE_MIN_LEVELS

    # w = 1 / (cells this donor contributes to this frame), so each donor sums to 1.
    # Recomputed for every fitted subset.
    donor_weights <- function(donor_fac) {
        if (!isTRUE(W_DONOR_WEIGHTS)) return(rep(1, length(donor_fac)))
        f <- as_donor_factor(donor_fac)
        as.numeric(1 / as.numeric(table(f)[f]))
    }

    # Level the donor factor and recompute the weights from this frame.
    prep_fit_frame <- function(data) {
        if ("donor" %in% names(data)) data$donor <- as_donor_factor(data$donor)
        data$w <- if ("donor" %in% names(data)) donor_weights(data$donor)
                  else rep(1, nrow(data))
        data
    }

    # The one place that decides what happens to a frame with no w column.
    ensure_weights <- function(data) {
        if (!("w" %in% names(data))) {
            data$w <- if ("donor" %in% names(data))
                donor_weights(data$donor) else rep(1, nrow(data))
        }
        data
    }

    # Wrappers around the top-level helpers in Section 0d, with the arm's GAM_BS.
    make_gam_formula <- function(k, use_re) p10b_gam_formula(k, use_re, GAM_BS)
    model_has_re     <- p10b_model_has_re

    # Population-level prediction: donor is supplied at an observed level and then
    # excluded. Returns predict.gam's list with se = TRUE, a vector otherwise.
    predict_population <- function(model, x_seq, fit_data, se = TRUE) {
        nd <- data.frame(GOI = x_seq)
        if (model_has_re(model)) {
            lev <- if (!is.null(fit_data) && "donor" %in% names(fit_data))
                as_donor_factor(fit_data$donor)[1] else NULL
            if (is.null(lev)) {
                # Fall back to a donor level stored in the model's var.summary.
                vs  <- model$var.summary[["donor"]]
                lev <- if (!is.null(vs) && is.factor(vs)) vs[1] else NULL
            }
            if (is.null(lev))
                return(predict(model, newdata = data.frame(GOI = x_seq),
                               se.fit = se))
            nd$donor <- lev
            return(predict(model, newdata = nd, se.fit = se,
                           exclude = DONOR_RE_LABEL))
        }
        predict(model, newdata = nd, se.fit = se)
    }

    # Null deviance against a donor-aware null, cached by frame size, response
    # sum, donor structure and weights.
    .p10b_null_dev_cache <- new.env(hash = TRUE, parent = emptyenv())
    get_null_deviance <- function(data) {
        if ("donor" %in% names(data)) data$donor <- as_donor_factor(data$donor)
        data   <- ensure_weights(data)
        use_re <- donor_re_ok(data)
        key <- paste0(nrow(data), "_", round(sum(data$Expression), 6),
                      "_", n_donor_levels(data), "_", as.integer(use_re),
                      "_", round(sum(data$w^2), 9))
        if (exists(key, envir = .p10b_null_dev_cache))
            return(get(key, envir = .p10b_null_dev_cache))
        set.seed(123)
        # The null carries the same weights as the fitted model.
        nd <- if (use_re)
            deviance(gam(Expression ~ 1 + s(donor, bs = "re"),
                         data = data, weights = w, method = GAM_METHOD))
        else
            deviance(gam(Expression ~ 1, data = data, weights = w))
        assign(key, nd, envir = .p10b_null_dev_cache)
        nd
    }

    # Deviance explained against the donor-aware null.
    de_explained <- function(model, data) {
        if (is.null(model) || is.null(data)) return(NA_real_)
        # Without a w column, the model's prior weights are used.
        if (!("w" %in% names(data)) && !is.null(model$prior.weights) &&
            length(model$prior.weights) == nrow(data))
            data$w <- as.numeric(model$prior.weights)
        nd <- tryCatch(get_null_deviance(data), error = function(e) NA_real_)
        if (!is.finite(nd) || nd <= 0) return(NA_real_)
        as.numeric((nd - deviance(model)) / nd)
    }

    # summary()$dev.expl, for the DE_vs_flat_null columns.
    de_flat_null <- function(model) {
        if (is.null(model)) return(NA_real_)
        tryCatch(as.numeric(summary(model)$dev.expl), error = function(e) NA_real_)
    }

    # Weights for the EP ratio and the TREP count: the frame's w, else the model's
    # prior weights, else unit weights.
    ep_weight_vec <- function(model, data) {
        if ("w" %in% names(data)) return(as.numeric(data$w))
        if (!is.null(model) && !is.null(model$prior.weights) &&
            length(model$prior.weights) == nrow(data))
            return(as.numeric(model$prior.weights))
        rep(1, nrow(data))
    }

    # Weighted population-level deviance explained, 1 - RSS_w / TSS_w, from the
    # same prediction vector as ep_values. It sizes the TREP set.
    pop_dev_explained <- function(model, data, mod_fitted) {
        if (is.null(model) || is.null(data) || is.null(mod_fitted))
            return(NA_real_)
        if (!("Expression" %in% names(data))) return(NA_real_)
        if (length(mod_fitted) != nrow(data)) return(NA_real_)
        w  <- ep_weight_vec(model, data)
        y  <- as.numeric(data$Expression)
        mf <- as.numeric(mod_fitted)
        ok <- is.finite(y) & is.finite(mf) & is.finite(w) & w > 0
        if (sum(ok) < 2L) return(NA_real_)
        ybar  <- sum(w[ok] * y[ok]) / sum(w[ok])
        denom <- sum(w[ok] * (y[ok] - ybar)^2)
        if (!is.finite(denom) || denom <= 0) return(NA_real_)
        as.numeric(1 - sum(w[ok] * (y[ok] - mf[ok])^2) / denom)
    }

    # TREP sizing fraction, clamped to [0, 1].
    trep_count_fraction <- function(de) {
        if (!is.finite(de)) return(NA_real_)
        max(0, min(1, de))
    }

    # The data argument is unused; kept for the existing call sites.
    calculate_prss <- function(model, data = NULL) p10b_calculate_prss(model)

    # A record of every fit in this arm, for the donor random-intercept check.
    .p10b_fit_registry <- list()
    # Called as register_fit(c(res, list(dataset = ds))); a NULL fit is skipped.
    register_fit <- function(res) {
        if (is.null(res) || is.null(res$label)) return(invisible(NULL))
        .p10b_fit_registry[[length(.p10b_fit_registry) + 1L]] <<- list(
            dataset  = if (is.null(res$dataset)) NA_character_ else res$dataset,
            label    = res$label,
            donor_re = isTRUE(res$donor_re),
            n_donors = res$n_donors)
        invisible(NULL)
    }


    # =============================================================================
    # SECTION A0c: FULL TWO-PHASE PRSS-REML GAM FITTER
    # =============================================================================
    # cl defaults to the arm's cluster; cl = NULL runs every iteration here.
    .fit_full_gam_p10b <- function(data, label,
                                   n_iter     = PRSS_N_ITER_P10,
                                   early_stop = PRSS_EARLY_STOP_P10,
                                   gamma      = GAM_GAMMA_P10,
                                   verbose    = TRUE,
                                   cl         = .p10b_cl) {
        set.seed(123)

        # Level donor and recompute the weights once, for every iteration and the null.
        data     <- prep_fit_frame(data)
        unique_x <- length(unique(data$GOI))
        max_k    <- min(10L, unique_x - 1L)
        n_don    <- n_donor_levels(data)
        use_re   <- donor_re_ok(data)

        if (verbose)
            cat(sprintf(paste0("    [%s] n=%d | unique_x=%d | max_k=%d | ",
                               "donors=%d | donor RE=%s | sum(w)=%.4f\n"),
                        label, nrow(data), unique_x, max_k, n_don,
                        if (use_re) "YES" else "NO", sum(data$w)))
        if (verbose && isTRUE(GAM_DONOR_RE) && !use_re)
            cat(sprintf(paste0("    Donor random intercept not fitted: %d ",
                               "donor level(s), minimum is %d. This fit is ",
                               "cell level.\n"),
                        n_don, DONOR_RE_MIN_LEVELS))

        if (unique_x < 4L) {
            if (verbose) cat("    Skipping: < 4 unique GOI values\n")
            return(NULL)
        }

        best_prss <- Inf; best_model <- NULL; best_k <- NULL
        best_iter <- 0L; prss_values <- numeric(n_iter)
        early_stop_counter <- 0L; early_stop_triggered <- FALSE
        # Indices of the iterations that ran.
        iters_run     <- integer(0)
        # Set when the winning iteration ran on a worker and must be refitted here.
        best_par_iter <- NULL

        # The k window around the current best.
        k_window <- function(bk)
            if (bk == max_k) c(max_k - 1L, max_k) else
            if (bk == 3L)    c(3L, 4L) else
            c(bk - 1L, bk, bk + 1L)
        k_range <- NULL

        # weights = w is the bare column name.
        run_local <- function(i, k) {
            tryCatch({
                iters_run <<- c(iters_run, i)
                set.seed(123L + i)
                m  <- gam(make_gam_formula(k, use_re),
                          data = data, weights = w,
                          method = GAM_METHOD,
                          select = GAM_SELECT, gamma = gamma)
                pc <- calculate_prss(m, data)
                list(ok = TRUE, i = i, k = k, prss = pc$PRSS,
                     model = m, local = TRUE)
            }, error = function(e)
                list(ok = FALSE, i = i, k = k, local = TRUE,
                     message = conditionMessage(e)))
        }

        # Applies one iteration's result in index order, for both the sequential and
        # the parallel paths. Only phase 2 updates k_range and counts errors toward
        # the early stop.
        apply_result <- function(i, r, phase) {
            if (!isTRUE(r$local)) iters_run <<- c(iters_run, i)
            if (!isTRUE(r$ok)) {
                if (verbose)
                    cat("      error iter", i, ":", r$message, "\n")
                prss_values[i] <<- NA
                if (phase == 2L)
                    early_stop_counter <<- early_stop_counter + 1L
                return(FALSE)
            }
            prss_values[i] <<- r$prss
            if (verbose)
                cat(sprintf("      %s iter %2d | k=%d | PRSS=%.4f\n",
                            if (phase == 1L) "init  " else "refine",
                            i, r$k, r$prss))
            if (r$prss < best_prss) {
                best_prss <<- r$prss; best_k <<- r$k; best_iter <<- i
                early_stop_counter <<- 0L
                if (!is.null(r$model)) {
                    best_model    <<- r$model
                    best_par_iter <<- NULL
                } else {
                    # Fitted on a worker, which returns summaries rather than
                    # a gam object. The master refits the winner once.
                    best_par_iter <<- list(i = i, k = r$k)
                }
                if (phase == 2L) k_range <<- k_window(best_k)
                TRUE
            } else {
                early_stop_counter <<- early_stop_counter + 1L
                FALSE
            }
        }

        use_cluster <- !is.null(cl) && length(cl) >= 2L

        if (use_cluster) {
            # The fitting frame and constants are exported once per fit.
            p10b_fit_ctx <- list(data   = data,
                                 bs     = GAM_BS,
                                 method = GAM_METHOD,
                                 select = GAM_SELECT,
                                 gamma  = gamma,
                                 use_re = use_re)
            # Only p10b_fit_ctx is exported from this environment.
            parallel::clusterExport(cl, "p10b_fit_ctx", envir = environment())
        }

        dispatch <- function(idx, ks) {
            if (!use_cluster)
                return(lapply(seq_along(idx),
                              function(j) run_local(idx[j], ks[j])))
            jobs <- Map(function(j, k) list(i = j, k = as.integer(k)), idx, ks)
            parallel::parLapply(cl, jobs, p10b_prss_worker)
        }

        # == Phase 1: exploration, k = 3..max_k ===============================
        # The eight iterations are independent, so they run as one batch and are
        # applied in index order.
        n_phase1 <- min(8L, n_iter)
        if (n_phase1 > 0L) {
            idx1 <- seq_len(n_phase1)
            ks1  <- 3L + ((idx1 - 1L) %% (max_k - 2L))
            out1 <- dispatch(idx1, ks1)
            for (j in seq_along(idx1)) apply_result(idx1[j], out1[[j]], 1L)
        }

        # Stop when every phase-1 iteration failed.
        if (best_iter == 0L || is.null(best_k)) {
            if (verbose)
                cat("    WARNING: all phase-1 iterations failed for",
                    label, "- aborting this fit\n")
            return(NULL)
        }
        k_range <- k_window(best_k)
        if (verbose) cat(sprintf("    Phase 1 best k=%d PRSS=%.4f\n",
                                 best_k, best_prss))

        # == Phase 2: adaptive refinement =====================================
        # At most 8 + PRSS_PHASE2_MAX_P10B iterations; skipped when n_iter < 9.
        phase2_end <- min(as.integer(n_iter),
                          8L + as.integer(PRSS_PHASE2_MAX_P10B))

        phase2_serial <- function(from_i) {
            if (phase2_end < 9L || from_i > phase2_end)
                return(invisible(NULL))
            for (i in from_i:phase2_end) {
                if (early_stop_counter >= early_stop) {
                    if (verbose) cat("      early stopping at iter", i, "\n")
                    early_stop_triggered <<- TRUE; break
                }
                if (i %% 3L == 0L) set.seed(200L + i)
                k <- if (i %% 3L == 0L) sample(k_range, 1L) else best_k
                apply_result(i, run_local(i, k), 2L)
            }
            invisible(NULL)
        }

        # Phase-2 iterations depend on earlier results. The parallel path runs a chunk
        # assuming no improvement within it; if one improves, the rest of the chunk is
        # discarded and run sequentially.
        phase2_parallel <- function() {
            n_workers <- max(1L, length(cl))
            i <- 9L
            while (i <= phase2_end) {
                if (early_stop_counter >= early_stop) {
                    if (verbose) cat("      early stopping at iter", i, "\n")
                    early_stop_triggered <<- TRUE
                    return(invisible(NULL))
                }
                idx <- i:min(i + n_workers - 1L, phase2_end)
                ks  <- vapply(idx, function(j) {
                    if (j %% 3L == 0L) {
                        set.seed(200L + j)
                        as.integer(sample(k_range, 1L))
                    } else as.integer(best_k)
                }, integer(1))
                out <- dispatch(idx, ks)

                for (jj in seq_along(idx)) {
                    j <- idx[jj]
                    if (early_stop_counter >= early_stop) {
                        if (verbose)
                            cat("      early stopping at iter", j, "\n")
                        early_stop_triggered <<- TRUE
                        return(invisible(NULL))
                    }
                    if (apply_result(j, out[[jj]], 2L) && j < phase2_end) {
                        if (verbose)
                            cat("      phase-2 speculation invalidated at iter",
                                j, "- remaining iterations run sequentially\n")
                        phase2_serial(j + 1L)
                        return(invisible(NULL))
                    }
                }
                i <- max(idx) + 1L
            }
            invisible(NULL)
        }

        if (phase2_end >= 9L) {
            if (use_cluster) phase2_parallel() else phase2_serial(9L)
        }

        # Refit the winning iteration here when it ran on a worker (deterministic).
        if (!is.null(best_par_iter) && identical(best_iter, best_par_iter$i)) {
            if (verbose)
                cat("      best model is from iter", best_iter,
                    "(k =", best_par_iter$k, ") on a worker - refitting here\n")
            set.seed(123L + best_iter)
            best_model <- gam(make_gam_formula(best_par_iter$k, use_re),
                              data = data, weights = w,
                              method = GAM_METHOD,
                              select = GAM_SELECT, gamma = gamma)
        }

        if (is.null(best_model)) {
            if (verbose)
                cat("    WARNING: no model retained for", label, "\n")
            return(NULL)
        }

        actual_iters <- if (length(iters_run) > 0L) max(iters_run) else 0L

        if (verbose)
            cat(sprintf(paste0("    Final best k=%d PRSS=%.4f  iter=%d  ",
                               "early_stop=%s\n"),
                        best_k, best_prss, best_iter, early_stop_triggered))

        list(best_model    = best_model,
             best_k        = best_k,
             best_prss     = best_prss,
             best_iter     = best_iter,
             donor_re      = use_re,
             n_donors      = n_don,
             prss_values   = prss_values[seq_len(max(1L, actual_iters))],
             actual_iters  = actual_iters,
             early_stopped = early_stop_triggered,
             data          = data,
             label         = label)
    }


    # =============================================================================
    # SECTION A0d: SELF-TESTS (bodies in Section Z)
    # =============================================================================
    # Both switches are FALSE by default; neither test touches the arm's files.
    if (isTRUE(RUN_SMOKE_TEST)) {
        p10b_smoke_test(environment())
        p10b_parallel_equivalence(environment())
    }
    if (isTRUE(RUN_STRESS_TEST)) p10b_stress_test(environment())


    # =============================================================================
    # SECTION A: SHARED HELPERS
    # =============================================================================

    # == A1. Plot theme =======================================================
    theme_p10b <- function(base_size = 10) {
        theme_minimal(base_size = base_size) +
            theme(
                plot.title       = element_text(face = "bold", size = base_size + 1),
                plot.subtitle    = element_text(size = base_size - 1),
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

    # == A5. show_plot: print to the Plots pane ===============================
    show_plot <- function(p) {
        tryCatch(print(p), error = function(e) invisible(NULL))
    }

    # == A6. Save PDF and PNG =================================================
    save_pp <- function(p, base, pdf_w = 10, pdf_h = 8,
                        png_w = 2000, png_h = 1600, res = 200) {
        # On error, close only the devices opened here.
        .dev_entry <- grDevices::dev.cur()
        .unwind <- function() {
            while (grDevices::dev.cur() != .dev_entry &&
                   grDevices::dev.cur() != 1L)
                tryCatch(grDevices::dev.off(),
                         error = function(e2) return(invisible(NULL)))
        }
        tryCatch({
            pdf(paste0(base, ".pdf"), width = pdf_w, height = pdf_h)
            print(p)
            dev.off()
        }, error = function(e) {
            .unwind()
            cat("  [save_pp] PDF error:", conditionMessage(e), "\n")
        })
        tryCatch({
            png(paste0(base, ".png"), width = png_w, height = png_h, res = res)
            print(p)
            dev.off()
        }, error = function(e) {
            .unwind()
            cat("  [save_pp] PNG error:", conditionMessage(e), "\n")
        })
    }


    # =============================================================================
    # SECTION B: MAIN DIAGNOSTIC FUNCTION (runs per dataset)
    # =============================================================================

    run_p10b_diagnostic <- function(dataset_tag) {

        meta      <- DATASET_META[[dataset_tag]]
        regions   <- meta$regions
        gam_dir   <- meta$gam_dir
        dfg_genes <- meta$dfg_genes
        IP_thresh <- meta$IP_thresh

        # IP_thresh is the reported IP, printed in titles; IP_ACTIVE is the threshold
        # the Pre/Post split uses, on the IP_CLASS_AXIS axis.
        IP_y_thresh <- meta$IP_y_thresh
        if (is.null(IP_y_thresh)) IP_y_thresh <- NA_real_
        # Without IP_y_thresh the cDFG split stops; there is no fallback to the GOI axis.
        IP_AXIS_USED <- IP_CLASS_AXIS
        if (identical(IP_CLASS_AXIS, "cDFG") && !is.finite(IP_y_thresh))
            stop("IP_CLASS_AXIS is \"cDFG\" but the IPRS RDS carries no ",
                 "IP_y_thresh for ", dataset_tag, ". Rerun Part10b_CP10K_IPRS_code.r for this arm.",
                 call. = FALSE)
        IP_ACTIVE <- if (identical(IP_AXIS_USED, "cDFG")) IP_y_thresh
                     else IP_thresh
        IP_ACTIVE_LABEL <- if (identical(IP_AXIS_USED, "cDFG"))
            "IP_y = f(IP_x) on cDFG" else "IP_x on GOI"
        if (is.finite(IP_ACTIVE))
            cat(sprintf(paste0("  Pre/Post split on %s | reported ",
                               "IP_x = %s | threshold used = %.6f\n"),
                        IP_ACTIVE_LABEL,
                        if (is.finite(IP_thresh)) sprintf("%.6f", IP_thresh)
                        else "NA",
                        IP_ACTIVE))


        cat("\n", strrep("=", 72), "\n")
        cat("  Part10b diagnostic:", dataset_tag, " [", ARM_TAG, "]\n")
        cat(strrep("=", 72), "\n\n")

        # == B1. Load combined SKCS scores ====================================
        scores_csv <- p10b_scores_csv(dataset_tag)
        if (!file.exists(scores_csv)) {
            cat("  ERROR: scores CSV not found. Looked in:\n")
            for (.p in p10b_scores_candidates(dataset_tag))
                cat("    ", .p, "\n", sep = "")
            return(invisible(NULL))
        }

        scores_all <- read.csv(scores_csv, stringsAsFactors = FALSE,
                               check.names = FALSE)
        cat("  Loaded", nrow(scores_all), "cells from", basename(scores_csv), "\n")

        # Split: GOI-positive scored cells vs GOI-zero excluded cells
        scores_pos  <- scores_all[!is.na(scores_all$spear_contrib), ]
        scores_goi0 <- scores_all[ is.na(scores_all$spear_contrib), ]
        cat("  GOI-positive (scored):", nrow(scores_pos),
            "| GOI-zero (excluded):", nrow(scores_goi0), "\n\n")

        spear  <- scores_pos$spear_norm
        kend   <- scores_pos$kend_contrib
        is_skcs <- scores_pos$is_skcs
        n_total <- nrow(scores_pos)

        # == B2. Load combined GOI-positive GAM data for threshold testing ====
        cat("  Loading combined GAM data for GAM fitting ...\n")
        region_dfs <- list()
        for (reg in regions) {
            gp <- file.path(gam_dir, "gam_data",
                            paste0(reg, "_",
                                   gsub("-", "_", gsub("/", "_", PRIMARY_STYPE)),
                                   "_selected.csv"))
            if (!file.exists(gp)) next
            df        <- read.csv(gp, stringsAsFactors = FALSE, check.names = FALSE)
            df$region <- reg

            # Donor is required in every frame that reaches gam().
            if (!("donor" %in% colnames(df)))
                stop(sprintf(paste0("Region %s: gam_data has no donor column, ",
                                    "so the donor random intercept cannot be ",
                                    "fitted."), reg), call. = FALSE)
            if (anyNA(df$donor))
                stop(sprintf("Region %s: %d cell(s) have a missing donor label.",
                             reg, sum(is.na(df$donor))), call. = FALSE)
            df$donor  <- as_donor_factor(df$donor)

            # which(), so an NA condition drops the row instead of adding an NA row.
            df        <- df[which(df$GOI > 0), ]
            if (nrow(df) >= MIN_CELLS) {
                cat(sprintf("    Region %-10s %6d cells | %3d donors\n",
                            reg, nrow(df), n_donor_levels(df$donor)))
                region_dfs[[reg]] <- df
            }
        }
        if (length(region_dfs) == 0) {
            cat("  ERROR: no GAM data CSVs found - GAM comparisons will be skipped\n\n")
        }
        df_combined <- if (length(region_dfs) > 0) {
            out <- do.call(rbind, region_dfs)
            rownames(out) <- NULL
            # Rebuild the pooled donor factor from character values, and its weights.
            prep_fit_frame(out)
        } else NULL

        # Align scores_pos with df_combined by barcode so GAM subsets are correct
        if (!is.null(df_combined)) {
            df_combined$barcode_match <- df_combined$barcode
            scores_pos$barcode_match  <- scores_pos$barcode
            common_bc <- intersect(df_combined$barcode, scores_pos$barcode)
            df_combined <- df_combined[df_combined$barcode %in% common_bc, ]
            scores_pos  <- scores_pos[ scores_pos$barcode  %in% common_bc, ]
            df_combined <- df_combined[order(df_combined$barcode), ]
            scores_pos  <- scores_pos[ order(scores_pos$barcode),  ]
            spear   <- scores_pos$spear_norm
            kend    <- scores_pos$kend_contrib
            is_skcs <- scores_pos$is_skcs
            n_total <- nrow(scores_pos)
            # Rebuild donor levels and weights after the alignment.
            df_combined <- prep_fit_frame(df_combined)

            # Non-finite GOI or Expression would be dropped by gam() and misalign the EP
            # ranking, so stop.
            .bad_goi <- !is.finite(df_combined$GOI)
            .bad_exp <- !is.finite(df_combined$Expression)
            if (any(.bad_goi) || any(.bad_exp))
                stop(sprintf(paste0("%s: %d of %d aligned cell(s) have a ",
                                    "non-finite GOI and %d a non-finite ",
                                    "Expression. Fix the Part 8 gam_data ",
                                    "export."),
                            dataset_tag, sum(.bad_goi), nrow(df_combined),
                            sum(.bad_exp)),
                     call. = FALSE)

            cat("  Aligned", n_total, "cells between scores CSV and GAM data\n")
            cat(sprintf("  Pooled frame: %d cells | %d donors | sum(w)=%.4f\n\n",
                        nrow(df_combined), n_donor_levels(df_combined),
                        sum(df_combined$w)))
        }

        # Median thresholds of the aligned cell set (unnamed).
        cur_S_thresh <- unname(quantile(spear, 0.50))
        cur_K_thresh <- unname(quantile(kend,  0.50))


        # =============================================================================
        # OBJECTIVE 1: Non-SKCS cell stratification
        # =============================================================================
        cat(strrep("-", 60), "\n")
        cat("  OBJECTIVE 1: Non-SKCS cell stratification\n")
        cat(strrep("-", 60), "\n\n")

        non_skcs_idx <- which(!is_skcs)
        s_ns  <- spear[non_skcs_idx]
        k_ns  <- kend[non_skcs_idx]

        # Three-way stratification
        inflection_zone <- s_ns >= 0   & k_ns < cur_K_thresh
        low_expressor   <- s_ns < 0    & k_ns < 0
        discordant      <- !(inflection_zone | low_expressor)

        strat_df <- data.frame(
            barcode       = scores_pos$barcode[non_skcs_idx],
            region        = scores_pos$region[ non_skcs_idx],
            spear_norm    = s_ns,
            kend_contrib  = k_ns,
            stratum       = ifelse(inflection_zone, "Inflection-zone",
                                   ifelse(low_expressor,   "Low-expressor",
                                          "Discordant")),
            stringsAsFactors = FALSE
        )

        strat_summary <- strat_df %>%
            group_by(stratum) %>%
            summarise(
                n              = n(),
                pct_of_nonskcs = round(100 * n() / nrow(strat_df), 1),
                pct_of_total   = round(100 * n() / n_total, 1),
                mean_spear     = round(mean(spear_norm), 3),
                mean_kend      = round(mean(kend_contrib),  3),
                .groups = "drop"
            )

        # Sheets of the per-dataset workbook, written at the end of the dataset.
        combined_tabs <- list()

        cat("  Non-SKCS stratification (dataset:", dataset_tag, "):\n")
        print(strat_summary, n = Inf)
        cat("\n")

        strat_csv <- file.path(OUT_DIR,
                               paste0(FILE_PREFIX, "_", dataset_tag,
                                      "_nonSKCS_stratification.csv"))
        write.csv(strat_df, strat_csv, row.names = FALSE)
        cat("  Saved:", basename(strat_csv), "\n")

        # Quadrant scatter
        quad_df <- data.frame(
            spear_norm    = spear,
            kend_contrib  = kend,
            Cell_group    = "Low-expressor",
            stringsAsFactors = FALSE
        )
        quad_df$Cell_group[ is_skcs]                              <- "SKCS-core"
        quad_df$Cell_group[!is_skcs & spear >= 0]                <- "Inflection-zone"
        quad_df$Cell_group[!is_skcs & spear <  0 & kend >= 0]   <- "Discordant"

        group_cols <- c("SKCS-core"       = "#4B0082",
                        "Inflection-zone" = "#E69F00",
                        "Low-expressor"   = "#B0B0B0",
                        "Discordant"      = "#56B4E9")

        n_by_group <- table(quad_df$Cell_group)
        group_labs <- setNames(
            paste0(names(n_by_group), "\nn=",
                   as.integer(n_by_group), " (",
                   round(100 * as.integer(n_by_group) / n_total, 1), "%)"),
            names(n_by_group)
        )


        # =============================================================================
        # OBJECTIVE 2: Fixed-threshold selection (S_norm > 0.05 and K > 0.5)
        # =============================================================================
        cat(strrep("-", 60), "\n")
        cat("  OBJECTIVE 2: Fixed-threshold selection (S_norm>0.05, K>0.5)\n")
        cat(strrep("-", 60), "\n\n")

        # Fixed thresholds
        FIXED_S_THRESH <- 0.05
        FIXED_K_THRESH <- 0.5

        sel_fixed    <- spear > FIXED_S_THRESH & kend > FIXED_K_THRESH
        n_fixed      <- sum(sel_fixed)
        pct_fixed    <- round(100 * n_fixed / n_total, 1)

        cat(sprintf("  Fixed thresholds: S > %.2f  AND  K > %.2f\n",
                    FIXED_S_THRESH, FIXED_K_THRESH))
        cat(sprintf("  Cells selected  : n = %d  (%.1f%% of GOI-positive cells)\n\n",
                    n_fixed, pct_fixed))

        # Reference arm: cells above the 50th percentiles of S and K.
        sel_50    <- spear > cur_S_thresh & kend > cur_K_thresh
        n_50pct   <- sum(sel_50)
        pct_50pct <- round(100 * n_50pct / n_total, 1)
        cat(sprintf(paste0("  50%% quantile thresholds: S > %.4f  AND  K > %.4f\n",
                           "  Cells selected  : n = %d  (%.1f%% of GOI-positive",
                           " cells)\n\n"),
                    cur_S_thresh, cur_K_thresh, n_50pct, pct_50pct))

        # Checks that this file's fixed cut and Part 10's is_skcs select the same cells.
        n_disagree <- sum(sel_fixed != is_skcs)
        gr10b(paste0("GR-11-", dataset_tag),
              if (n_disagree == 0L) "PASS" else "WARN",
              sprintf(paste0("this file's fixed cut (S>%.2f, K>%.2f) vs Part 10's",
                             " is_skcs: %d cell(s) differ%s"),
                      FIXED_S_THRESH, FIXED_K_THRESH, n_disagree,
                      if (n_disagree == 0L)
                          " (same thresholds)"
                      else
                          paste0(" -- Part 10's SKCS_S_THRESH/SKCS_K_THRESH and",
                                 " this file's FIXED_S_THRESH/FIXED_K_THRESH",
                                 " differ")))

        # Fit full PRSS GAM at fixed threshold
        de_fixed <- NA_real_
        res_fixed <- NULL
        res_50    <- NULL
        if (!is.null(df_combined) && n_fixed >= MIN_CELLS) {
            bc_fixed  <- scores_pos$barcode[sel_fixed]
            # Weights are recomputed from this subset's own donors.
            df_fixed  <- prep_fit_frame(
                df_combined[df_combined$barcode %in% bc_fixed, ])
            cat("  Fitting cDFG GAM at fixed threshold (full PRSS pipeline) ...\n")
            cat(sprintf("    frame: %d cells | %d donors | sum(w)=%.4f\n",
                        nrow(df_fixed), n_donor_levels(df_fixed),
                        sum(df_fixed$w)))
            res_fixed <- .fit_full_gam_p10b(df_fixed, label = "Fixed_Snorm05_K05",
                                            verbose = FALSE)
            register_fit(c(res_fixed, list(dataset = dataset_tag)))
            # DE against the donor-aware null.
            de_fixed  <- if (!is.null(res_fixed) && !is.null(res_fixed$best_model))
                de_explained(res_fixed$best_model, res_fixed$data) * 100
            else NA_real_
            cat(sprintf("  Fixed threshold DE: %.2f%%  (donor RE: %s, donors: %s)\n\n",
                        ifelse(is.na(de_fixed), NA, de_fixed),
                        if (!is.null(res_fixed) && isTRUE(res_fixed$donor_re))
                            "YES" else "NO",
                        if (is.null(res_fixed)) "NA" else
                            format(res_fixed$n_donors)))
        } else {
            cat("  GAM skipped (n < MIN_CELLS or no GAM data)\n\n")
        }

        # Reference DE at 50% threshold
        de_50pct <- NA_real_
        if (!is.null(df_combined) && n_50pct >= MIN_CELLS) {
            bc_50    <- scores_pos$barcode[sel_50]
            df_50    <- prep_fit_frame(
                df_combined[df_combined$barcode %in% bc_50, ])
            cat("  Fitting reference cDFG GAM at 50% threshold ...\n")
            cat(sprintf("    frame: %d cells | %d donors | sum(w)=%.4f\n",
                        nrow(df_50), n_donor_levels(df_50), sum(df_50$w)))
            res_50   <- .fit_full_gam_p10b(df_50, label = "Ref_50pct",
                                           verbose = FALSE)
            register_fit(c(res_50, list(dataset = dataset_tag)))
            de_50pct <- if (!is.null(res_50) && !is.null(res_50$best_model))
                de_explained(res_50$best_model, res_50$data) * 100
            else NA_real_
            cat(sprintf("  50%% reference DE : %.2f%%\n\n",
                        ifelse(is.na(de_50pct), NA, de_50pct)))
        }

        # Comparison summary; DE_vs_flat_null is summary()$dev.expl, kept for comparison.
        .re_of  <- function(r) if (is.null(r)) NA else isTRUE(r$donor_re)
        .don_of <- function(r) if (is.null(r)) NA_integer_ else
            as.integer(r$n_donors)
        .flat_of <- function(r) if (is.null(r) || is.null(r$best_model))
            NA_real_ else de_flat_null(r$best_model) * 100

        thresh_all <- data.frame(
            Method          = c("Fixed_Snorm05_K05", "Reference_50pct"),
            dataset         = dataset_tag,
            S_thresh        = c(FIXED_S_THRESH, round(cur_S_thresh, 4)),
            K_thresh        = c(FIXED_K_THRESH, round(cur_K_thresh, 4)),
            n_cells         = c(n_fixed,  n_50pct),
            pct_total       = c(pct_fixed, pct_50pct),
            DE_pct          = round(c(de_fixed, de_50pct), 2),
            Donor_RE        = c(.re_of(res_fixed),  .re_of(res_50)),
            N_Donors        = c(.don_of(res_fixed), .don_of(res_50)),
            DE_vs_flat_null = round(c(.flat_of(res_fixed), .flat_of(res_50)), 2),
            stringsAsFactors = FALSE
        )

        cat("  Threshold comparison:\n")
        print(thresh_all, row.names = FALSE)
        cat("\n")

        # == Threshold scatter ================================================
        scatter_df <- data.frame(
            S = spear, K = kend,
            # Grouped by the fixed and 50th-percentile selections.
            group = ifelse(sel_fixed, "Fixed (S_norm>0.05, K>0.5)",
                           ifelse(sel_50,    "50% quantile",
                                  "Non-selected")),
            stringsAsFactors = FALSE
        )
        scatter_df$group[sel_fixed & sel_50]  <- "Both (fixed + 50%)"
        scatter_df$group[sel_fixed & !sel_50] <- "Fixed only (added)"
        scatter_df$group[!sel_fixed & sel_50] <- "50% only (removed)"

        group_cols2 <- c(
            "Both (fixed + 50%)"  = "#4B0082",
            "Fixed only (added)"  = "#009E73",
            "50% only (removed)"  = "firebrick",
            "Non-selected"        = "#D0D0D0"
        )

        p_thresh_scatter <- ggplot(
            scatter_df[order(scatter_df$group == "Non-selected",
                             decreasing = TRUE), ],
            aes(x = S, y = K, colour = group)) +
            geom_point(alpha = 0.28, size = 0.7) +
            geom_vline(xintercept = FIXED_S_THRESH, colour = "#009E73",
                       linetype = "solid",  linewidth = 0.8) +
            geom_hline(yintercept = FIXED_K_THRESH, colour = "#009E73",
                       linetype = "solid",  linewidth = 0.8) +
            geom_vline(xintercept = cur_S_thresh, colour = "#4B0082",
                       linetype = "dashed", linewidth = 0.6) +
            geom_hline(yintercept = cur_K_thresh, colour = "#4B0082",
                       linetype = "dashed", linewidth = 0.6) +
            scale_colour_manual(values = group_cols2, name = "Cell group") +
            annotate("text",
                     x = FIXED_S_THRESH + diff(range(spear)) * 0.02,
                     y = max(kend) * 0.98,
                     label = sprintf("S_norm>0.05 K>0.5\nn=%d (%.1f%%)%s",
                                     n_fixed, pct_fixed,
                                     ifelse(is.na(de_fixed), "",
                                            sprintf("  DE=%.1f%%", de_fixed))),
                     hjust = 0, size = 3.0, colour = "#009E73", fontface = "bold") +
            annotate("text",
                     x = min(spear) + diff(range(spear)) * 0.02,
                     y = cur_K_thresh - diff(range(kend)) * 0.05,
                     label = sprintf("50%% q\nn=%d (%.1f%%)%s",
                                     n_50pct, pct_50pct,
                                     ifelse(is.na(de_50pct), "",
                                            sprintf("  DE=%.1f%%", de_50pct))),
                     hjust = 0, size = 3.0, colour = "#4B0082") +
            labs(title    = paste0(dataset_tag,
                                   " \u2014 Fixed threshold (S_norm>0.05, K>0.5) vs 50% quantile"),
                 subtitle = sprintf(
                     "Green solid = fixed (S_norm>0.05, K>0.5) | Purple dashed = 50%% quantile\n%s",
                     paste0("Fixed: n=", n_fixed, " (", pct_fixed, "%)  ",
                            ifelse(is.na(de_fixed), "",
                                   sprintf("DE=%.1f%%  ", de_fixed)),
                            "| 50%: n=", n_50pct, " (", pct_50pct, "%)  ",
                            ifelse(is.na(de_50pct), "",
                                   sprintf("DE=%.1f%%", de_50pct)))),
                 x = "Normalised Spearman contribution (Si / Si_max)",
                 y = "Kendall concordance contribution") +
            theme_p10b() +
            theme(legend.position = "right", aspect.ratio = 1)

        show_plot(p_thresh_scatter)
        save_pp(p_thresh_scatter,
                file.path(OUT_DIR, paste0(FILE_PREFIX, "_", dataset_tag,
                                          "_fixed_threshold_scatter")),
                pdf_w = 9, pdf_h = 9, png_w = 1800, png_h = 1800)
        cat("  Fixed threshold scatter saved\n\n")

        combined_tabs[["Threshold_Comparison"]] <- thresh_all


        # =============================================================================
        # SECTION D2: GAM CURVES AT FIXED THRESHOLD (S_norm>0.05, K>0.5)
        # =============================================================================
        cat(strrep("-", 60), "\n")
        cat("  SECTION D2: GAM curves at fixed threshold (S_norm>0.05, K>0.5)\n")
        cat(strrep("-", 60), "\n\n")

        # Initialise all Section D2 objects to NULL so the invisible() return
        # is safe even if the conditional block below is skipped.
        res_dfg_d2       <- NULL
        res_hkg_d2       <- NULL
        df_fixed_dfg     <- NULL
        ep_df            <- NULL
        trep_barcodes    <- NULL
        nontrep_barcodes <- NULL
        ip_status_bc     <- NULL
        trep_flag_bc     <- NULL

        if (!is.null(df_combined) && n_fixed >= MIN_CELLS) {

            bc_fixed_d2  <- scores_pos$barcode[sel_fixed]
            # data_ep below is this frame, so the EP null and the per-cell EP values use
            # the same donors and weights as the fit.
            df_fixed_dfg <- prep_fit_frame(
                df_combined[df_combined$barcode %in% bc_fixed_d2, ])

            cat("  Fitting cDFG GAM (verbose) ...\n")
            res_dfg_d2 <- .fit_full_gam_p10b(df_fixed_dfg, label = "cDFG_fixed",
                                             verbose = TRUE)
            register_fit(c(res_dfg_d2, list(dataset = dataset_tag)))

            # == Load the cHKG composite from norm.mtx (skipped for the HKG8 arm) ====
            hkg_data_ok   <- FALSE
            df_fixed_hkg  <- NULL
            res_hkg_d2    <- NULL

            if (!is.null(HKG_10)) {
                cat("  Loading cHKG data ...\n")
                hkg_genes_use <- HKG_10[HKG_10 != "18S rRNA"]

                mtx_path_b  <- file.path(gam_dir, "deg_support", "norm.mtx")
                bar_path_b  <- file.path(gam_dir, "deg_support", "norm_barcodes.tsv")
                feat_path_b <- file.path(gam_dir, "deg_support", "norm_features.tsv")

                if (all(file.exists(c(mtx_path_b, bar_path_b, feat_path_b)))) {
                    tryCatch({
                        X_hkg   <- readMM(mtx_path_b)
                        hkg_bc  <- readLines(bar_path_b, warn = FALSE)
                        hkg_ft  <- read.table(feat_path_b, sep = "\t", header = FALSE,
                                              stringsAsFactors = FALSE)[, 2]
                        rownames(X_hkg) <- hkg_bc
                        colnames(X_hkg) <- hkg_ft
                        hkg_avail <- hkg_genes_use[hkg_genes_use %in% colnames(X_hkg)]
                        if (length(hkg_avail) > 0) {
                            bc_in_mtx <- bc_fixed_d2[bc_fixed_d2 %in% hkg_bc]
                            if (length(bc_in_mtx) >= MIN_CELLS) {
                                sub_mat  <- as.matrix(
                                    X_hkg[bc_in_mtx, hkg_avail, drop = FALSE])
                                hkg_expr <- if (length(hkg_avail) == 1)
                                    as.numeric(sub_mat[, 1])
                                else rowMeans(sub_mat)
                                goi_vals <- df_fixed_dfg$GOI[
                                    match(bc_in_mtx, df_fixed_dfg$barcode)]
                                keep <- !is.na(goi_vals) & goi_vals > 0
                                # Donor pulled across by barcode; weights from this cell set.
                                don_vals <- df_fixed_dfg$donor[
                                    match(bc_in_mtx, df_fixed_dfg$barcode)]
                                df_fixed_hkg <- prep_fit_frame(data.frame(
                                    GOI        = goi_vals[keep],
                                    Expression = hkg_expr[keep],
                                    barcode    = bc_in_mtx[keep],
                                    donor      = don_vals[keep],
                                    stringsAsFactors = FALSE))
                                hkg_data_ok <- nrow(df_fixed_hkg) >= MIN_CELLS
                                cat(sprintf(paste0("    cHKG: %d cells, %d ",
                                                   "genes, %d donors, ",
                                                   "sum(w)=%.4f\n"),
                                            nrow(df_fixed_hkg), length(hkg_avail),
                                            n_donor_levels(df_fixed_hkg),
                                            sum(df_fixed_hkg$w)))
                            }
                        }
                        rm(X_hkg); gc(verbose = FALSE)
                    }, error = function(e)
                        cat("    norm.mtx error:", conditionMessage(e), "\n"))
                } else {
                    cat("    norm.mtx not found - cHKG skipped\n")
                }

                if (hkg_data_ok && !is.null(df_fixed_hkg)) {
                    cat("  Fitting cHKG GAM (verbose) ...\n")
                    res_hkg_d2 <- .fit_full_gam_p10b(df_fixed_hkg, label = "cHKG_fixed",
                                                     verbose = TRUE)
                    register_fit(c(res_hkg_d2, list(dataset = dataset_tag)))
                }
            } else {
                cat("  cHKG comparison skipped (GOI is an HKG gene)\n")
            }

            # == Build GAM curve plots ========================================
            # Curve labels are arm-aware: "cDFG" for Module 1, "cHKG8" for
            # the control arm, so no plot shows "cDFG" over housekeeping genes.
            dfg_curve_label <- paste0(DFG_DISPLAY_LABEL, " (SKCS-core)")
            hkg_curve_label <- "cHKG (SKCS-core)"
            gs_colors_d2 <- setNames(c("#4B0082", "#4F75DE"),
                                     c(dfg_curve_label, hkg_curve_label))
            gs_fills_d2  <- gs_colors_d2
            plot_data_d2 <- data.frame()

            for (item in list(
                list(res = res_dfg_d2, label = dfg_curve_label),
                list(res = res_hkg_d2, label = hkg_curve_label))) {

                if (is.null(item$res) || is.null(item$res$best_model)) next
                m   <- item$res$best_model
                dat <- item$res$data

                pts <- data.frame(GOI = dat$GOI, Expression = dat$Expression,
                                  Gene_Set = item$label, se.fit = NA_real_,
                                  type = "data", stringsAsFactors = FALSE)
                x_seq  <- seq(min(dat$GOI), max(dat$GOI), length.out = 1000)
                fit_df <- tryCatch({
                    set.seed(123)
                    # Population-level curve (donor excluded).
                    pred <- predict_population(m, x_seq, dat, se = TRUE)
                    data.frame(GOI = x_seq,
                               Expression = as.numeric(pred$fit),
                               Gene_Set = item$label,
                               se.fit = as.numeric(pred$se.fit),
                               type = "fit", stringsAsFactors = FALSE)
                }, error = function(e) NULL)

                plot_data_d2 <- rbind(plot_data_d2, pts)
                if (!is.null(fit_df)) plot_data_d2 <- rbind(plot_data_d2, fit_df)
            }

            if (nrow(plot_data_d2) > 0) {
                de_dfg_d2 <- if (!is.null(res_dfg_d2) &&
                                 !is.null(res_dfg_d2$best_model))
                    sprintf("%s DE=%.1f%%",
                            dfg_curve_label,
                            de_explained(res_dfg_d2$best_model,
                                         res_dfg_d2$data) * 100)
                else NA_character_
                de_hkg_d2 <- if (!is.null(res_hkg_d2) &&
                                 !is.null(res_hkg_d2$best_model))
                    sprintf("cHKG (SKCS-core) DE=%.1f%%",
                            de_explained(res_hkg_d2$best_model,
                                         res_hkg_d2$data) * 100)
                else NA_character_
                de_labs_d2 <- c(de_dfg_d2, de_hkg_d2)[!is.na(c(de_dfg_d2, de_hkg_d2))]
                n_pts_d2   <- sum(plot_data_d2$type == "data" &
                                      plot_data_d2$Gene_Set == dfg_curve_label)

                ABS_PAD  <- 0.08
                x_d2 <- plot_data_d2$GOI[plot_data_d2$type == "data"]
                y_d2 <- plot_data_d2$Expression[plot_data_d2$type == "data"]
                x_lim_d2 <- c(min(x_d2) - ABS_PAD, max(x_d2) + ABS_PAD)
                y_lim_d2 <- c(max(0, min(y_d2) - ABS_PAD), max(y_d2) + ABS_PAD)

            } else {
                cat("  No plot data - GAM curve skipped\n\n")
            }

            # == EP / TREP calculation ========================================
            cat("  Computing per-cell EP and TREP assignment ...\n")

            if (!is.null(res_dfg_d2) && !is.null(res_dfg_d2$best_model)) {

                gam_model_ep  <- res_dfg_d2$best_model
                # The frame the model was fitted on (same rows as df_fixed_dfg).
                data_ep       <- res_dfg_d2$data

                # EP null: flat intercept model, as in the published CEP-IP method. Its fitted
                # value is mean(Expression), the same null used for the non-selected cells.
                null_model_ep <- gam(Expression ~ 1, data = data_ep)

                null_fitted_ep  <- fitted(null_model_ep)
                # Population-level prediction, so neither side of the ratio carries donor
                # offsets.
                model_fitted_ep <- as.numeric(predict(gam_model_ep,
                                                      newdata = data_ep,
                                                      exclude = DONOR_RE_LABEL,
                                                      type    = "response"))

                # One value per row, or the subtractions below would recycle.
                if (length(null_fitted_ep)  != nrow(data_ep) ||
                    length(model_fitted_ep) != nrow(data_ep))
                    stop(sprintf(paste0("%s: the EP vectors are the wrong ",
                                        "length (null %d, model %d, rows %d). ",
                                        "Rows were dropped by the fit, so the ",
                                        "TREP ranking would be misaligned."),
                                dataset_tag, length(null_fitted_ep),
                                length(model_fitted_ep), nrow(data_ep)),
                         call. = FALSE)

                null_resid_ep   <- data_ep$Expression - null_fitted_ep
                model_resid_ep  <- data_ep$Expression - model_fitted_ep

                null_sq_ep      <- null_resid_ep^2
                model_sq_ep     <- model_resid_ep^2

                # EP is 0 when the null residual is zero, as for the non-selected cells.
                ep_values <- ifelse(
                    null_sq_ep > 0,
                    1 - (model_sq_ep / null_sq_ep),
                    0
                )

                sorted_idx   <- order(ep_values, decreasing = TRUE, na.last = TRUE)
                # The TREP count is sized by the weighted population-level DE of the same
                # prediction. Three DE values are printed: de_trep_count (sets the split),
                # de_trep_flat_null (summary()$dev.expl) and de_trep_donor_aware.
                de_trep_flat_null   <- de_flat_null(gam_model_ep)
                de_trep_donor_aware <- de_explained(gam_model_ep, data_ep)
                de_trep_count       <- pop_dev_explained(gam_model_ep, data_ep,
                                                         model_fitted_ep)
                if (!is.finite(de_trep_count)) de_trep_count <- de_trep_flat_null
                model_de     <- de_trep_count
                # target_cells is clamped to [0, n]; a non-finite DE sizes the set to 0.
                n_ep      <- nrow(data_ep)
                .trep_frac <- trep_count_fraction(model_de)
                if (!is.finite(.trep_frac)) {
                    cat(paste0("  WARN: TREP sizing DE is not finite ",
                               "(pop_dev_explained and summary()$dev.expl failed); ",
                               "target_cells set to 0.\n"))
                    .trep_frac <- 0
                }
                target_cells <- max(0L, min(
                    as.integer(round(n_ep * .trep_frac)), n_ep))

                trep_barcodes    <- if (target_cells > 0L)
                    data_ep$barcode[sorted_idx[seq_len(target_cells)]]
                else character(0)
                nontrep_barcodes <- if (target_cells < n_ep)
                    data_ep$barcode[sorted_idx[seq.int(target_cells + 1L, n_ep)]]
                else character(0)

                n_trep    <- length(trep_barcodes)
                n_nontrep <- length(nontrep_barcodes)

                cat(sprintf(paste0("  EP computed: n_TREP=%d  n_nonTREP=%d  ",
                                   "(DE=%.1f%% population level, sets the ",
                                   "split | %.1f%% summary()$dev.expl | ",
                                   "%.1f%% donor-aware)\n"),
                            n_trep, n_nontrep, de_trep_count * 100,
                            de_trep_flat_null * 100,
                            de_trep_donor_aware * 100))

                # Per-donor TREP composition (diagnostic only).
                if ("donor" %in% names(data_ep) && target_cells > 0L) {
                    .d_tot  <- table(droplevels(as_donor_factor(data_ep$donor)))
                    .d_trep <- table(droplevels(as_donor_factor(
                        data_ep$donor[sorted_idx[seq_len(target_cells)]])))
                    .d_frac <- as.numeric(.d_trep[names(.d_tot)]) /
                               as.numeric(.d_tot)
                    .d_frac[is.na(.d_frac)] <- 0
                    cat(sprintf(paste0("  TREP donor coverage: %d/%d donors ",
                                       "represented | per-donor TREP fraction ",
                                       "[%.3f, %.3f] median %.3f\n"),
                                sum(.d_frac > 0), length(.d_tot),
                                min(.d_frac), max(.d_frac),
                                stats::median(.d_frac)))
                }

                # == IP-dependent analysis (finite IP, non-degenerate TREP split) ====
                if (is.finite(IP_ACTIVE) && n_trep > 0L && n_nontrep > 0L) {

                    # Four-way group. On "cDFG" a cell is Pre-IP when its composite expression
                    # is below IP_y_thresh; on "GOI" when GOI <= IP_thresh.
                    ip_status_bc <- setNames(
                        if (identical(IP_AXIS_USED, "cDFG"))
                            ifelse(data_ep$Expression <  IP_ACTIVE,
                                   "Pre-IP", "Post-IP")
                        else
                            ifelse(data_ep$GOI        <= IP_ACTIVE,
                                   "Pre-IP", "Post-IP"),
                        data_ep$barcode
                    )
                    trep_flag_bc <- setNames(
                        ifelse(data_ep$barcode %in% trep_barcodes, "TREP", "non-TREP"),
                        data_ep$barcode
                    )

                    four_cols <- c(
                        "Pre-IP TREP"      = "#E69F00",
                        "Pre-IP non-TREP"  = "#1B3A8C",
                        "Post-IP TREP"     = "#009E73",
                        "Post-IP non-TREP" = "#D94F4F"
                    )

                    trep_groups    <- c("Pre-IP TREP",     "Post-IP TREP")
                    nontrep_groups <- c("Pre-IP non-TREP", "Post-IP non-TREP")


                    # == Plot 3: Four-colour GAM curve ========================
                    cat("  Building four-colour GAM curve plot ...\n")

                    plot_data_4c <- plot_data_d2
                    plot_data_4c$four_group <- NA_character_

                    is_dfg_data <- plot_data_4c$type == "data" &
                        plot_data_4c$Gene_Set == dfg_curve_label

                    dfg_rows <- which(is_dfg_data)
                    stopifnot(
                        "dfg_rows length must equal nrow(data_ep)" =
                            length(dfg_rows) == nrow(data_ep)
                    )
                    plot_data_4c$four_group[dfg_rows] <- paste(
                        ip_status_bc[data_ep$barcode],
                        trep_flag_bc[data_ep$barcode]
                    )

                    draw_ord_4c <- c(
                        "Pre-IP non-TREP"  = 1L,
                        "Post-IP non-TREP" = 2L,
                        "Pre-IP TREP"      = 3L,
                        "Post-IP TREP"     = 4L
                    )
                    plot_data_4c$draw_ord <- draw_ord_4c[
                        match(plot_data_4c$four_group, names(draw_ord_4c))]
                    plot_data_4c$draw_ord[is.na(plot_data_4c$draw_ord)] <- 0L

                    plot_data_4c_sorted <- plot_data_4c[
                        order(plot_data_4c$draw_ord), ]

                    nonsel_barcodes_gam <- setdiff(df_combined$barcode, data_ep$barcode)
                    df_nonsel_gam <- df_combined[df_combined$barcode %in%
                                                     nonsel_barcodes_gam, ]
                    df_nonsel_gam$four_group <- "Non SKCS-core"
                    cat(sprintf("  Non-SKCS-core cells for GAM bg layer: %d\n",
                                nrow(df_nonsel_gam)))

                    # Non-SKCS-core cells in grey.
                    four_cols_curve <- c(four_cols, "Non SKCS-core" = "#A9A9A9")
                    n_curve_4c      <- table(plot_data_4c$four_group[dfg_rows])
                    n_nonsel_gam    <- nrow(df_nonsel_gam)

                    set.seed(123)
                    p_curve_4c <- ggplot(plot_data_4c_sorted,
                                         aes(x = GOI, y = Expression)) +
                        geom_point(
                            data = df_nonsel_gam,
                            aes(x = GOI, y = Expression, colour = four_group),
                            alpha = 0.25, size = 0.6,
                            inherit.aes = FALSE) +
                        geom_point(
                            data = subset(plot_data_4c_sorted,
                                          type == "data" &
                                              Gene_Set == dfg_curve_label),
                            aes(colour = four_group),
                            position = position_jitter(width = 0.01,
                                                       height = 0, seed = 123),
                            alpha = 0.50, size = 0.9) +
                        geom_line(
                            data = subset(plot_data_4c_sorted,
                                          type == "fit" &
                                              Gene_Set == dfg_curve_label),
                            colour = "#4B0082", linewidth = 1.3) +
                        geom_ribbon(
                            data = subset(plot_data_4c_sorted,
                                          type == "fit" &
                                              Gene_Set == dfg_curve_label),
                            aes(ymin = Expression - 1.96 * se.fit,
                                ymax = Expression + 1.96 * se.fit),
                            fill = "#4B0082", alpha = 0.12, colour = NA) +
                        # The IP line is drawn on the axis the split was made on.
                        (if (identical(IP_AXIS_USED, "cDFG"))
                             geom_hline(yintercept = IP_ACTIVE,
                                        linetype = "dashed", colour = "black",
                                        linewidth = 0.7)
                         else
                             geom_vline(xintercept = IP_ACTIVE,
                                        linetype = "dashed", colour = "black",
                                        linewidth = 0.7)) +
                        (if (identical(IP_AXIS_USED, "cDFG"))
                             annotate("text",
                                      x = min(plot_data_4c_sorted$GOI,
                                              na.rm = TRUE),
                                      y = IP_ACTIVE,
                                      label = sprintf("IP_y=%.3f", IP_ACTIVE),
                                      hjust = 0, vjust = -0.5, size = 2.8,
                                      colour = "black")
                         else
                             annotate("text",
                                      x = IP_ACTIVE + 0.05,
                                      y = y_lim_d2[2] * 0.97,
                                      label = sprintf("IP=%.3f", IP_ACTIVE),
                                      hjust = 0, size = 2.8,
                                      colour = "black")) +
                        scale_colour_manual(
                            values = four_cols_curve,
                            labels = sapply(names(four_cols_curve), function(g) {
                                n_g <- if (g == "Non SKCS-core") n_nonsel_gam
                                       else as.integer(n_curve_4c[g])
                                if (is.na(n_g)) n_g <- 0L
                                sprintf("%s  n=%d", g, n_g)
                            }),
                            name   = "Cell group",
                            breaks = c("Pre-IP TREP", "Pre-IP non-TREP",
                                       "Post-IP TREP", "Post-IP non-TREP",
                                       "Non SKCS-core")) +
                        scale_x_continuous(
                            labels = scales::number_format(accuracy = 0.1)) +
                        scale_y_continuous(
                            labels = scales::number_format(accuracy = 0.1)) +
                        coord_cartesian(xlim = x_lim_d2, ylim = y_lim_d2) +
                        labs(title    = paste0("GAM (S_norm>0.05, K>0.5) 4-group: ",
                                               dataset_tag, " \u2014 Combined",
                                               "  [n=", n_pts_d2,
                                               ", gamma=", GAM_GAMMA_P10, "]"),
                             subtitle = paste0(
                                 de_dfg_d2,
                                 sprintf(paste0("  |  IP=%.3f  TREP=%d  ",
                                                "non-TREP=%d  |  split on %s"),
                                         IP_thresh, n_trep, n_nontrep,
                                         IP_ACTIVE_LABEL)),
                             x        = paste0(GOI, " Expression (", NORM_TAG,
                                               " log1p, GOI>0)"),
                             y        = paste0("Gene Set Mean Expression (",
                                               NORM_TAG, " log1p)")) +
                        theme_p10b(base_size = 10) +
                        theme(legend.position = "right",
                              plot.title      = element_text(face = "bold", size = 11),
                              aspect.ratio    = 1)

                    show_plot(p_curve_4c)
                    save_pp(p_curve_4c,
                            file.path(OUT_DIR, paste0(FILE_PREFIX, "_", dataset_tag,
                                                      "_Combined_GAM_curves_4group")),
                            pdf_w = 10, pdf_h = 9, png_w = 2000, png_h = 1800)
                    cat("  GAM curve (4-group) saved\n")

                    # == Plot 4: Four-colour SKCS scatter =====================
                    cat("  Building four-colour SKCS scatter ...\n")

                    skcs_4c <- data.frame(
                        S       = spear,
                        K       = kend,
                        barcode = scores_pos$barcode,
                        stringsAsFactors = FALSE
                    )
                    skcs_4c$group_4c <- "Non-selected"

                    in_fixed_mask <- skcs_4c$barcode %in% data_ep$barcode

                    skcs_4c$group_4c[in_fixed_mask] <- paste(
                        ip_status_bc[skcs_4c$barcode[in_fixed_mask]],
                        trep_flag_bc[skcs_4c$barcode[in_fixed_mask]]
                    )

                    draw_ord_skcs <- c(
                        "Non-selected"     = 0L,
                        "Pre-IP non-TREP"  = 1L,
                        "Post-IP non-TREP" = 2L,
                        "Pre-IP TREP"      = 3L,
                        "Post-IP TREP"     = 4L
                    )
                    skcs_4c$draw_ord <- draw_ord_skcs[
                        match(skcs_4c$group_4c, names(draw_ord_skcs))]
                    skcs_4c$draw_ord[is.na(skcs_4c$draw_ord)] <- 0L
                    skcs_4c <- skcs_4c[order(skcs_4c$draw_ord), ]

                    # Same grey as the curve plot.
                    four_scatter_cols <- c(four_cols, "Non-selected" = "#A9A9A9")

                    n_4c <- table(skcs_4c$group_4c)

                    p_skcs_4c <- ggplot(skcs_4c,
                                        aes(x = S, y = K)) +
                        geom_point(data = subset(skcs_4c, group_4c == "Non-selected"),
                                   aes(colour = group_4c), alpha = 0.25, size = 0.6) +
                        geom_point(data = subset(skcs_4c,
                                                 group_4c %in% c(trep_groups,
                                                                 nontrep_groups)),
                                   aes(colour = group_4c),
                                   alpha = 0.35, size = 0.65) +
                        geom_vline(xintercept = FIXED_S_THRESH, colour = "#009E73",
                                   linetype = "solid", linewidth = 0.7) +
                        geom_hline(yintercept = 0.5, colour = "#009E73",
                                   linetype = "solid", linewidth = 0.7) +
                        scale_colour_manual(
                            values = four_scatter_cols,
                            labels = sapply(names(four_scatter_cols), function(g) {
                                n_g <- as.integer(n_4c[g])
                                if (is.na(n_g)) n_g <- 0L
                                disp <- if (g == "Non-selected") "Non SKCS-core" else g
                                sprintf("%s  n=%d", disp, n_g)
                            }),
                            name   = "Cell group",
                            breaks = c("Pre-IP TREP", "Pre-IP non-TREP",
                                       "Post-IP TREP", "Post-IP non-TREP",
                                       "Non-selected")) +
                        labs(title    = paste0(dataset_tag,
                                               " \u2014 SKCS scatter: 4-group",
                                               " (IP \u00d7 TREP)"),
                             subtitle = sprintf(
                                 paste0("Fixed threshold: S_norm>0.05, K>0.5",
                                        "  |  IP=%.3f  |  TREP=%d  non-TREP=%d",
                                        "  |  split on %s"),
                                 IP_thresh, n_trep, n_nontrep,
                                 IP_ACTIVE_LABEL),
                             x = "Normalised Spearman contribution (Si / Si_max)",
                             y = "Kendall concordance contribution") +
                        theme_p10b() +
                        theme(legend.position = "right", aspect.ratio = 1)

                    show_plot(p_skcs_4c)
                    save_pp(p_skcs_4c,
                            file.path(OUT_DIR, paste0(FILE_PREFIX, "_", dataset_tag,
                                                      "_SKCS_scatter_4group")),
                            pdf_w = 10, pdf_h = 9, png_w = 2000, png_h = 1800)
                    cat("  SKCS scatter (4-group) saved\n\n")

                    # == Plot 5: 3D SKCS scatter (S x K x EP) =================
                    cat("  Building 3D SKCS-EP scatter (plotly) ...\n")
                    cat("  Computing EP for non-selected cells (same null as SKCS-core) ...\n")

                    # plotly and htmlwidgets are installed if missing; if unavailable, only the
                    # 3D HTML is skipped.
                    .have_plotly_3d <- TRUE
                    for (.pkg in c("plotly", "htmlwidgets")) {
                        if (requireNamespace(.pkg, quietly = TRUE)) next
                        cat("    ", .pkg, " not installed; attempting install ...\n",
                            sep = "")
                        tryCatch(
                            install.packages(.pkg,
                                             repos = "https://cloud.r-project.org",
                                             quiet = TRUE),
                            error   = function(e) invisible(NULL),
                            warning = function(w) invisible(NULL))
                        if (!requireNamespace(.pkg, quietly = TRUE)) {
                            cat("    [WARN] ", .pkg, " is unavailable; the 3D SKCS ",
                                "scatter HTML is skipped",
                                ".\n", sep = "")
                            .have_plotly_3d <- FALSE
                        }
                    }
                    if (.have_plotly_3d) {
                        suppressPackageStartupMessages({
                            library(plotly); library(htmlwidgets)
                        })
                    }

                    symlog10 <- function(x) sign(x) * log10(1 + abs(x))

                    gam_model_3d <- res_dfg_d2$best_model

                    y_bar_skcs <- mean(df_fixed_dfg$Expression)
                    cat(sprintf("    y_bar (SKCS-core null): %.6f\n", y_bar_skcs))

                    all_goi_df <- df_combined[df_combined$barcode %in% scores_pos$barcode, ]
                    all_goi_df <- all_goi_df[order(all_goi_df$barcode), ]
                    sp_aligned <- scores_pos[match(all_goi_df$barcode, scores_pos$barcode), ]

                    ep_skcs_lookup <- setNames(ep_values, data_ep$barcode)

                    nonsel_barcodes <- setdiff(all_goi_df$barcode, data_ep$barcode)
                    nonsel_df <- all_goi_df[all_goi_df$barcode %in% nonsel_barcodes, ]

                    if (nrow(nonsel_df) > 0) {
                        # Population-level prediction; these cells are outside the fitting frame.
                        gam_pred_nonsel <- as.numeric(
                            predict_population(gam_model_3d, nonsel_df$GOI,
                                               df_fixed_dfg, se = FALSE)
                        )
                        null_sq_nonsel  <- (nonsel_df$Expression - y_bar_skcs)^2
                        model_sq_nonsel <- (nonsel_df$Expression - gam_pred_nonsel)^2
                        ep_nonsel <- ifelse(
                            null_sq_nonsel > 0,
                            1 - (model_sq_nonsel / null_sq_nonsel),
                            0
                        )
                        ep_nonsel_lookup <- setNames(ep_nonsel, nonsel_df$barcode)
                    } else {
                        ep_nonsel_lookup <- setNames(numeric(0), character(0))
                    }

                    ep_all_consistent <- numeric(nrow(all_goi_df))
                    for (i in seq_len(nrow(all_goi_df))) {
                        bc <- all_goi_df$barcode[i]
                        if (bc %in% names(ep_skcs_lookup)) {
                            ep_all_consistent[i] <- ep_skcs_lookup[bc]
                        } else if (bc %in% names(ep_nonsel_lookup)) {
                            ep_all_consistent[i] <- ep_nonsel_lookup[bc]
                        } else {
                            ep_all_consistent[i] <- 0
                        }
                    }

                    cat(sprintf("    SKCS-core cells: %d (EP reused from ep_values)\n",
                                sum(all_goi_df$barcode %in% names(ep_skcs_lookup))))
                    cat(sprintf("    Non-selected cells: %d (EP computed with same null)\n",
                                length(nonsel_barcodes)))
                    cat(sprintf("    EP range: [%.1f, %.1f]\n",
                                min(ep_all_consistent), max(ep_all_consistent)))

                    # Verify TREP consistency
                    trep_ep_check <- ep_all_consistent[
                        all_goi_df$barcode %in% trep_barcodes]
                    nontrep_ep_check <- ep_all_consistent[
                        all_goi_df$barcode %in% nontrep_barcodes]
                    trep_cutoff_3d <- min(trep_ep_check, na.rm = TRUE)
                    n_nontrep_above <- sum(nontrep_ep_check > trep_cutoff_3d, na.rm = TRUE)

                    cat(sprintf("    TREP cutoff EP: %.4f\n", trep_cutoff_3d))
                    cat(sprintf("    non-TREP cells above TREP cutoff: %d (should be 0)\n",
                                n_nontrep_above))
                    if (n_nontrep_above > 0) {
                        cat("    *** WARNING: TREP/non-TREP EP inconsistency detected ***\n")
                    }

                    # == 5B. Assign 5-group labels ============================
                    skcs_3d <- data.frame(
                        S       = sp_aligned$spear_norm,
                        K       = sp_aligned$kend_contrib,
                        barcode = all_goi_df$barcode,
                        EP      = ep_all_consistent,
                        GOI     = all_goi_df$GOI,
                        # Carried so the IP marker can be placed on the cDFG axis.
                        Expression = all_goi_df$Expression,
                        stringsAsFactors = FALSE
                    )

                    skcs_3d$in_fixed <- skcs_3d$barcode %in% data_ep$barcode
                    skcs_3d$group_5c <- "Non-selected"
                    skcs_3d$group_5c[skcs_3d$in_fixed] <- paste(
                        ip_status_bc[skcs_3d$barcode[skcs_3d$in_fixed]],
                        trep_flag_bc[skcs_3d$barcode[skcs_3d$in_fixed]]
                    )

                    cat("    Group counts:\n")
                    print(table(skcs_3d$group_5c))

                    # == 5C. Diagnostic table =================================
                    diag_groups <- sort(unique(skcs_3d$group_5c))
                    diag_summary_3d <- do.call(rbind, lapply(diag_groups, function(grp) {
                        g <- skcs_3d[skcs_3d$group_5c == grp, ]
                        data.frame(
                            dataset    = dataset_tag,
                            group      = grp,
                            n          = nrow(g),
                            EP_min     = round(min(g$EP), 4),
                            EP_q25     = round(quantile(g$EP, 0.25), 4),
                            EP_median  = round(median(g$EP), 4),
                            EP_mean    = round(mean(g$EP), 4),
                            EP_q75     = round(quantile(g$EP, 0.75), 4),
                            EP_max     = round(max(g$EP), 4),
                            n_EP_pos   = sum(g$EP > 0),
                            n_EP_neg   = sum(g$EP <= 0),
                            pct_EP_pos = round(100 * sum(g$EP > 0) / nrow(g), 1),
                            stringsAsFactors = FALSE, row.names = NULL
                        )
                    }))

                    overall_row_3d <- data.frame(
                        dataset = dataset_tag, group = "ALL_GOI_positive",
                        n = nrow(skcs_3d),
                        EP_min = round(min(skcs_3d$EP), 4),
                        EP_q25 = round(quantile(skcs_3d$EP, 0.25), 4),
                        EP_median = round(median(skcs_3d$EP), 4),
                        EP_mean = round(mean(skcs_3d$EP), 4),
                        EP_q75 = round(quantile(skcs_3d$EP, 0.75), 4),
                        EP_max = round(max(skcs_3d$EP), 4),
                        n_EP_pos = sum(skcs_3d$EP > 0),
                        n_EP_neg = sum(skcs_3d$EP <= 0),
                        pct_EP_pos = round(100 * sum(skcs_3d$EP > 0) / nrow(skcs_3d), 1),
                        stringsAsFactors = FALSE, row.names = NULL
                    )

                    cutoff_row_3d <- data.frame(
                        dataset = dataset_tag, group = "TREP_CUTOFF_INFO",
                        n = sum(skcs_3d$barcode %in% trep_barcodes),
                        EP_min = round(trep_cutoff_3d, 4),
                        EP_q25 = NA_real_, EP_median = NA_real_, EP_mean = NA_real_,
                        EP_q75 = NA_real_, EP_max = NA_real_,
                        n_EP_pos = NA_integer_, n_EP_neg = NA_integer_,
                        # In the TREP_CUTOFF_INFO row this column holds the model's DE (percent).
                        pct_EP_pos = round(de_explained(gam_model_3d,
                                                        df_fixed_dfg) * 100, 1),
                        stringsAsFactors = FALSE, row.names = NULL
                    )

                    diag_summary_3d <- rbind(diag_summary_3d, overall_row_3d, cutoff_row_3d)

                    combined_tabs[["EP_Boundary_Diagnostic"]] <- diag_summary_3d
                    cat("  EP diagnostic -> Excel sheet",
                        "EP_Boundary_Diagnostic\n\n")

                    cat("  EP distribution by 5-group:\n")
                    print(diag_summary_3d[, c("group", "n", "EP_median",
                                              "EP_mean", "pct_EP_pos")],
                          row.names = FALSE)
                    cat("\n")

                    # == 5D. Build 3D scatter =================================
                    plotly_cols_5g <- c(
                        "Pre-IP TREP"      = "#E69F00",
                        "Pre-IP non-TREP"  = "#1B3A8C",
                        "Post-IP TREP"     = "#009E73",
                        "Post-IP non-TREP" = "#D94F4F",
                        "Non-selected"     = "#787878"
                    )

                    group_order_5g <- c(
                        "Non-selected"     = 1L,
                        "Pre-IP non-TREP"  = 2L,
                        "Post-IP non-TREP" = 3L,
                        "Pre-IP TREP"      = 4L,
                        "Post-IP TREP"     = 5L
                    )
                    skcs_3d$draw_ord <- group_order_5g[skcs_3d$group_5c]
                    skcs_3d <- skcs_3d[order(skcs_3d$draw_ord), ]

                    n_5g <- table(skcs_3d$group_5c)

                    skcs_3d$EP_plot <- symlog10(skcs_3d$EP)

                    cat(sprintf("    EP range (raw):    [%.1f, %.1f]\n",
                                min(skcs_3d$EP), max(skcs_3d$EP)))
                    cat(sprintf("    EP range (symlog): [%.2f, %.2f]\n",
                                min(skcs_3d$EP_plot), max(skcs_3d$EP_plot)))

                    skcs_3d$hover <- sprintf(
                        "Group: %s<br>S: %.3f<br>K: %.3f<br>EP: %.4f<br>symlog10(EP): %.3f",
                        skcs_3d$group_5c, skcs_3d$S, skcs_3d$K,
                        skcs_3d$EP, skcs_3d$EP_plot
                    )

                    skcs_3d$group_5c <- factor(
                        skcs_3d$group_5c, levels = names(plotly_cols_5g)
                    )

                    # Only the rendering is guarded; later sections do not use fig_3d.
                    if (isTRUE(.have_plotly_3d)) tryCatch({
                        fig_3d <- plot_ly()

                        for (grp in levels(skcs_3d$group_5c)) {
                            sub <- skcs_3d[skcs_3d$group_5c == grp, ]
                            if (nrow(sub) == 0) next
                            n_grp <- as.integer(n_5g[grp])
                            if (is.na(n_grp)) n_grp <- 0L
                            mk_size    <- if (grp == "Non-selected") 2.5 else 3
                            mk_opacity <- if (grp == "Non-selected") 0.2 else 0.55
                            fig_3d <- fig_3d %>%
                                add_trace(
                                    data = sub,
                                    x = ~S, y = ~K, z = ~EP_plot,
                                    type = "scatter3d", mode = "markers",
                                    marker = list(
                                        size    = mk_size,
                                        color   = plotly_cols_5g[grp],
                                        opacity = mk_opacity,
                                        line    = list(width = 0)
                                    ),
                                    text = ~hover, hoverinfo = "text",
                                    name = sprintf("%s (n=%d)", grp, n_grp)
                                )
                        }

                        trep_cutoff_symlog <- symlog10(trep_cutoff_3d)
                        s_range_3d <- range(skcs_3d$S)
                        k_range_3d <- range(skcs_3d$K)

                        boundary_x <- c(s_range_3d[1], s_range_3d[2], s_range_3d[2],
                                        s_range_3d[1], s_range_3d[1])
                        boundary_y <- c(k_range_3d[1], k_range_3d[1], k_range_3d[2],
                                        k_range_3d[2], k_range_3d[1])
                        boundary_z <- rep(trep_cutoff_symlog, 5)

                        fig_3d <- fig_3d %>%
                            add_trace(
                                type = "scatter3d", mode = "lines",
                                x = boundary_x, y = boundary_y, z = boundary_z,
                                line = list(
                                    color = "rgba(150, 100, 200, 0.50)",
                                    width = 3, dash  = "dash"
                                ),
                                hoverinfo = "text",
                                text = sprintf("TREP cutoff: EP = %.4f", trep_cutoff_3d),
                                name = sprintf("TREP cutoff (EP=%.4f)", trep_cutoff_3d),
                                showlegend = TRUE
                            )

                        symlog_ticks_val <- c(-1000, -100, -10, -1, 0, 0.5, 1)
                        symlog_ticks_pos <- symlog10(symlog_ticks_val)
                        symlog_ticks_lab <- as.character(symlog_ticks_val)

                        fig_3d <- fig_3d %>%
                            layout(
                                title = list(
                                    text = sprintf(
                                        "%s \u2014 3D SKCS scatter (S \u00d7 K \u00d7 EP)<br><sup>All GOI-positive cells | EP uses SKCS-core null (y_bar=%.4f) | TREP cutoff EP=%.4f | z = symlog10(EP)</sup>",
                                        dataset_tag, y_bar_skcs, trep_cutoff_3d),
                                    font = list(size = 14)
                                ),
                                scene = list(
                                    xaxis = list(
                                        title = list(text = "Spearman contribution (S)",
                                                     font = list(size = 11)),
                                        gridcolor = "rgb(230,230,230)",
                                        zerolinecolor = "rgb(100,100,100)"
                                    ),
                                    yaxis = list(
                                        title = list(text = "Kendall contribution (K)",
                                                     font = list(size = 11)),
                                        gridcolor = "rgb(230,230,230)",
                                        zerolinecolor = "rgb(100,100,100)"
                                    ),
                                    zaxis = list(
                                        title = list(text = "Explanatory Power  [symlog10(EP)]",
                                                     font = list(size = 11)),
                                        gridcolor = "rgb(230,230,230)",
                                        zerolinecolor = "rgb(100,100,100)",
                                        tickvals = symlog_ticks_pos,
                                        ticktext = symlog_ticks_lab
                                    ),
                                    camera = list(eye = list(x = 1.5, y = 1.5, z = 1.0)),
                                    bgcolor = "white"
                                ),
                                legend = list(font = list(size = 10), itemsizing = "constant"),
                                paper_bgcolor = "white", plot_bgcolor = "white",
                                margin = list(t = 80, b = 20, l = 20, r = 20)
                            )

                        # == 5E. IP location in S x K x EP space ==============
                        # Cells near the threshold, on the axis the split was made on.
                        .ip3d_axis <- if (identical(IP_AXIS_USED, "cDFG") &&
                                          "Expression" %in% names(skcs_3d))
                            skcs_3d$Expression else skcs_3d$GOI
                        ip_window <- 0.05 * diff(range(.ip3d_axis, na.rm = TRUE))
                        ip_cells <- skcs_3d[is.finite(.ip3d_axis) &
                                            abs(.ip3d_axis - IP_ACTIVE) <= ip_window, ]

                        if (nrow(ip_cells) >= 5) {
                            ip_S  <- median(ip_cells$S, na.rm = TRUE)
                            ip_K  <- median(ip_cells$K, na.rm = TRUE)
                            ip_EP <- median(ip_cells$EP, na.rm = TRUE)
                            ip_EP_sym <- symlog10(ip_EP)

                            cat(sprintf("    IP in 3D space (n=%d cells within +/-%.3f of IP=%.3f):\n",
                                        nrow(ip_cells), ip_window, IP_ACTIVE))
                            cat(sprintf("      median S=%.2f, median K=%.4f, median EP=%.4f\n",
                                        ip_S, ip_K, ip_EP))

                            fig_3d <- fig_3d %>%
                                add_trace(
                                    type = "scatter3d", mode = "markers",
                                    x = ip_S, y = ip_K, z = ip_EP_sym,
                                    marker = list(
                                        size = 6, color = "rgba(48, 0, 80, 0.95)",
                                        symbol = "cross",
                                        line = list(width = 1, color = "black")),
                                    hoverinfo = "text",
                                    text = sprintf(
                                        "Inflection Point<br>GOI=%.3f<br>median S=%.2f<br>median K=%.4f<br>median EP=%.4f",
                                        IP_ACTIVE, ip_S, ip_K, ip_EP),
                                    name = sprintf("IP (%s=%.3f)",
                                                   if (identical(IP_AXIS_USED, "cDFG"))
                                                       "cDFG" else "GOI",
                                                   IP_ACTIVE),
                                    showlegend = TRUE)

                            fig_3d <- fig_3d %>%
                                add_trace(
                                    type = "scatter3d", mode = "lines",
                                    x = rep(ip_S, 2), y = rep(ip_K, 2),
                                    z = range(skcs_3d$EP_plot, na.rm = TRUE),
                                    line = list(color = "rgba(48, 0, 80, 0.6)",
                                                width = 1.5),
                                    hoverinfo = "skip", showlegend = FALSE)

                            fig_3d <- fig_3d %>%
                                layout(
                                    title = list(
                                        text = sprintf(
                                            "%s \u2014 3D SKCS scatter (S \u00d7 K \u00d7 EP)<br><sup>IP at GOI=%.3f \u2192 S=%.1f, K=%.3f, EP=%.3f | TREP cutoff EP=%.4f | z = symlog10(EP)</sup>",
                                            dataset_tag, IP_ACTIVE, ip_S, ip_K, ip_EP,
                                            trep_cutoff_3d),
                                        font = list(size = 14)))
                        } else {
                            cat(sprintf("    IP window: only %d cells near IP=%.3f - marker skipped\n",
                                        nrow(ip_cells), IP_ACTIVE))
                        }

                        html_3d_path <- file.path(OUT_DIR, paste0(
                            FILE_PREFIX, "_", dataset_tag, "_SKCS_3D_scatter_EP.html"))
                        tryCatch({
                            htmlwidgets::saveWidget(
                                fig_3d, file = html_3d_path, selfcontained = TRUE,
                                title = paste0(dataset_tag,
                                               " - 3D SKCS scatter (S x K x EP)"))
                            cat("    (self-contained HTML)\n")
                        }, error = function(e) {
                            cat("    selfcontained failed (pandoc?), saving with libs/ ...\n")
                            htmlwidgets::saveWidget(
                                fig_3d, file = html_3d_path, selfcontained = FALSE,
                                title = paste0(dataset_tag,
                                               " - 3D SKCS scatter (S x K x EP)"))
                        })

                        cat(sprintf("  3D SKCS-EP scatter saved: %s\n\n",
                                    basename(html_3d_path)))
                    }, error = function(e) {
                        cat("  [WARN] the 3D SKCS scatter failed to render: ",
                            conditionMessage(e), "\n",
                            sep = "")
                    }) else
                        cat("  3D SKCS scatter skipped (plotly unavailable).\n\n")

                    # =============================================================================
                    # Analysis A: Bulbous edge EP analysis
                    # =============================================================================
                    cat("\n  Bulbous edge EP analysis ...\n")

                    trep_5c <- skcs_3d[grepl("TREP$", skcs_3d$group_5c) &
                                           !grepl("non-TREP", skcs_3d$group_5c), ]

                    bulb_rows <- list()
                    for (ip_grp in c("Pre-IP TREP", "Post-IP TREP")) {
                        sub_trep <- trep_5c[trep_5c$group_5c == ip_grp, ]
                        if (nrow(sub_trep) < 20) next

                        ep_q75 <- quantile(sub_trep$EP, 0.75)
                        ep_q25 <- quantile(sub_trep$EP, 0.25)

                        edge_cells <- sub_trep[sub_trep$EP >= ep_q75, ]
                        fall_cells <- sub_trep[sub_trep$EP <= ep_q25, ]

                        for (metric in c("EP", "K")) {
                            edge_v <- edge_cells[[metric]]
                            fall_v <- fall_cells[[metric]]

                            sw_e <- if (length(edge_v) >= 3 & length(edge_v) <= 5000)
                                suppressWarnings(shapiro.test(edge_v)$p.value)
                            else NA_real_
                            sw_f <- if (length(fall_v) >= 3 & length(fall_v) <= 5000)
                                suppressWarnings(shapiro.test(fall_v)$p.value)
                            else NA_real_

                            both_norm <- !is.na(sw_e) && !is.na(sw_f) &&
                                sw_e > 0.05 && sw_f > 0.05

                            if (both_norm) {
                                tt <- t.test(edge_v, fall_v, alternative = "greater")
                                test_name <- "Welch_t"
                                p_raw     <- tt$p.value
                            } else {
                                wt <- wilcox.test(edge_v, fall_v,
                                                  alternative = "greater",
                                                  exact = FALSE)
                                test_name <- "MannWhitneyU"
                                p_raw     <- wt$p.value
                            }

                            bulb_rows[[length(bulb_rows) + 1]] <- data.frame(
                                dataset = dataset_tag, IP_group = ip_grp,
                                metric = metric,
                                n_edge = length(edge_v), n_falling = length(fall_v),
                                edge_median = round(median(edge_v), 6),
                                fall_median = round(median(fall_v), 6),
                                edge_mean = round(mean(edge_v), 6),
                                fall_mean = round(mean(fall_v), 6),
                                sw_edge_p = round(sw_e, 4),
                                sw_fall_p = round(sw_f, 4),
                                both_normal = both_norm,
                                test_used = test_name, p_raw = p_raw,
                                stringsAsFactors = FALSE, row.names = NULL
                            )
                        }
                    }

                    bulb_df <- do.call(rbind, bulb_rows)
                    if (!is.null(bulb_df) && nrow(bulb_df) > 0) {
                        bulb_df$p_BH <- p.adjust(bulb_df$p_raw, method = "BH")
                        cat("  Bulbous edge analysis:\n")
                        print(bulb_df[, c("IP_group", "metric", "n_edge",
                                          "n_falling", "edge_median",
                                          "fall_median", "test_used", "p_BH")])
                        combined_tabs[["Bulbous_Edge"]] <- bulb_df
                        cat("  Bulbous edge -> Excel sheet Bulbous_Edge\n")
                    }

                    # =============================================================================
                    # Analysis B: Ki dispersion pre-IP vs post-IP TREP
                    # =============================================================================
                    cat("\n  Ki dispersion analysis (pre-IP vs post-IP TREP) ...\n")

                    trep_skcs <- data.frame(
                        barcode   = data_ep$barcode,
                        GOI_expr  = data_ep$GOI,
                        EP_value  = ep_values,
                        is_TREP   = data_ep$barcode %in% trep_barcodes,
                        IP_status = ip_status_bc[data_ep$barcode],
                        stringsAsFactors = FALSE
                    )
                    trep_skcs <- trep_skcs[trep_skcs$is_TREP == TRUE, ]
                    trep_skcs$S <- scores_pos$spear_norm[
                        match(trep_skcs$barcode, scores_pos$barcode)]
                    trep_skcs$K <- scores_pos$kend_contrib[
                        match(trep_skcs$barcode, scores_pos$barcode)]

                    # Residual against the population curve.
                    trep_skcs$gam_pred <- as.numeric(predict_population(
                        gam_model_3d, trep_skcs$GOI_expr,
                        df_fixed_dfg, se = FALSE))
                    trep_skcs$resid_abs <- abs(
                        df_fixed_dfg$Expression[
                            match(trep_skcs$barcode, df_fixed_dfg$barcode)] -
                            trep_skcs$gam_pred)

                    pre_trep  <- trep_skcs[trep_skcs$IP_status == "Pre-IP", ]
                    post_trep <- trep_skcs[trep_skcs$IP_status == "Post-IP", ]

                    # Each group needs at least 3 cells.
                    disp_df <- NULL
                    if (nrow(pre_trep) >= 3L && nrow(post_trep) >= 3L) {

                    disp_rows <- list()
                    for (metric in c("K", "resid_abs")) {
                        pre_v  <- pre_trep[[metric]]
                        post_v <- post_trep[[metric]]

                        sw_pre <- if (length(pre_v) >= 3 & length(pre_v) <= 5000)
                            suppressWarnings(shapiro.test(pre_v)$p.value)
                        else NA_real_
                        sw_post <- if (length(post_v) >= 3 & length(post_v) <= 5000)
                            suppressWarnings(shapiro.test(post_v)$p.value)
                        else NA_real_

                        both_norm <- !is.na(sw_pre) && !is.na(sw_post) &&
                            sw_pre > 0.05 && sw_post > 0.05

                        if (metric == "K") {
                            alt <- "greater"; grp1 <- post_v; grp2 <- pre_v
                        } else {
                            alt <- "greater"; grp1 <- pre_v; grp2 <- post_v
                        }

                        if (both_norm) {
                            tt <- t.test(grp1, grp2, alternative = alt)
                            test_name <- "Welch_t"; p_raw <- tt$p.value
                        } else {
                            wt <- wilcox.test(grp1, grp2,
                                              alternative = alt, exact = FALSE)
                            test_name <- "MannWhitneyU"; p_raw <- wt$p.value
                        }

                        pooled_sd <- sqrt(
                            ((length(pre_v) - 1) * var(pre_v) +
                                 (length(post_v) - 1) * var(post_v)) /
                                (length(pre_v) + length(post_v) - 2))
                        cohens_d <- abs(mean(post_v) - mean(pre_v)) /
                            (pooled_sd + 1e-12)

                        disp_rows[[length(disp_rows) + 1]] <- data.frame(
                            dataset = dataset_tag, metric = metric,
                            n_pre_IP_TREP = length(pre_v),
                            n_post_IP_TREP = length(post_v),
                            pre_median = round(median(pre_v), 6),
                            post_median = round(median(post_v), 6),
                            pre_mean = round(mean(pre_v), 6),
                            post_mean = round(mean(post_v), 6),
                            pre_sd = round(sd(pre_v), 6),
                            post_sd = round(sd(post_v), 6),
                            delta_median = round(median(post_v) - median(pre_v), 6),
                            sw_pre_p = round(sw_pre, 4),
                            sw_post_p = round(sw_post, 4),
                            both_normal = both_norm,
                            test_used = test_name,
                            hypothesis = if (metric == "K") "post > pre" else "pre > post",
                            p_raw = p_raw, cohens_d = round(cohens_d, 4),
                            stringsAsFactors = FALSE, row.names = NULL
                        )
                    }

                    disp_df <- do.call(rbind, disp_rows)
                    disp_df$p_BH <- p.adjust(disp_df$p_raw, method = "BH")

                    cat("  Ki dispersion results:\n")
                    print(disp_df[, c("metric", "pre_median", "post_median",
                                      "delta_median", "cohens_d",
                                      "test_used", "p_BH")])

                    trep_box_cols <- c("Pre-IP TREP" = "#E69F00",
                                       "Post-IP TREP" = "#009E73")

                    # == Boxplot 1: Kendall K =================================
                    k_pval <- disp_df$p_BH[disp_df$metric == "K"]
                    k_d    <- disp_df$cohens_d[disp_df$metric == "K"]

                    box_k_df <- data.frame(
                        K = c(pre_trep$K, post_trep$K),
                        IP_status = factor(
                            c(rep("Pre-IP TREP", nrow(pre_trep)),
                              rep("Post-IP TREP", nrow(post_trep))),
                            levels = c("Pre-IP TREP", "Post-IP TREP")),
                        stringsAsFactors = FALSE
                    )

                    set.seed(42)
                    box_k_jitter <- do.call(rbind, lapply(
                        levels(box_k_df$IP_status), function(grp) {
                            sub <- box_k_df[box_k_df$IP_status == grp, ]
                            if (nrow(sub) > 300L) sub[sample(nrow(sub), 300L), ]
                            else sub
                        }))

                    p_Ki_box <- ggplot(box_k_df, aes(x = IP_status, y = K,
                                                      fill = IP_status)) +
                        geom_jitter(data = box_k_jitter, aes(colour = IP_status),
                                    width = 0.20, height = 0, size = 1.8,
                                    alpha = 0.35, shape = 16) +
                        geom_boxplot(colour = "grey20", linewidth = 0.50,
                                     width = 0.45, outlier.shape = NA,
                                     alpha = 0.82) +
                        scale_fill_manual(values = trep_box_cols, guide = "none") +
                        scale_colour_manual(values = trep_box_cols, guide = "none") +
                        scale_x_discrete(expand = expansion(mult = 0, add = 0.55)) +
                        labs(title = paste0(dataset_tag,
                                            " \u2014 Kendall K: Pre-IP vs Post-IP TREP"),
                             subtitle = sprintf(
                                 "BH FDR p = %.2e | Cohen's d = %.3f | Pre median = %.4f | Post median = %.4f",
                                 k_pval, k_d, median(pre_trep$K), median(post_trep$K)),
                             x = "", y = "Kendall contribution (K)") +
                        theme_minimal(base_size = 10) +
                        theme(plot.background = element_rect(fill = "white", colour = NA),
                              panel.background = element_rect(fill = "white", colour = NA),
                              panel.grid.major = element_line(colour = "gray90"),
                              panel.grid.minor = element_blank(),
                              panel.border = element_rect(colour = "black", fill = NA,
                                                          linewidth = 0.35),
                              strip.text = element_text(face = "bold", size = 10),
                              plot.subtitle = element_text(colour = "gray30"),
                              legend.position = "none")

                    show_plot(p_Ki_box)
                    save_pp(p_Ki_box,
                            file.path(OUT_DIR, paste0(FILE_PREFIX, "_", dataset_tag,
                                "_Ki_prepost_TREP_boxplot")),
                            pdf_w = 8.4, pdf_h = 4.9,
                            png_w = 1820L, png_h = 1120L, res = 200)
                    cat("  Ki boxplot saved\n")

                    # == Boxplot 2: |residual| ================================
                    r_pval <- disp_df$p_BH[disp_df$metric == "resid_abs"]
                    r_d    <- disp_df$cohens_d[disp_df$metric == "resid_abs"]

                    box_r_df <- data.frame(
                        resid_abs = c(pre_trep$resid_abs, post_trep$resid_abs),
                        IP_status = factor(
                            c(rep("Pre-IP TREP", nrow(pre_trep)),
                              rep("Post-IP TREP", nrow(post_trep))),
                            levels = c("Pre-IP TREP", "Post-IP TREP")),
                        stringsAsFactors = FALSE
                    )

                    set.seed(42)
                    box_r_jitter <- do.call(rbind, lapply(
                        levels(box_r_df$IP_status), function(grp) {
                            sub <- box_r_df[box_r_df$IP_status == grp, ]
                            if (nrow(sub) > 300L) sub[sample(nrow(sub), 300L), ]
                            else sub
                        }))

                    p_resid_box <- ggplot(box_r_df,
                                          aes(x = IP_status, y = resid_abs,
                                              fill = IP_status)) +
                        geom_jitter(data = box_r_jitter, aes(colour = IP_status),
                                    width = 0.20, height = 0, size = 1.8,
                                    alpha = 0.35, shape = 16) +
                        geom_boxplot(colour = "grey20", linewidth = 0.50,
                                     width = 0.45, outlier.shape = NA,
                                     alpha = 0.82) +
                        scale_fill_manual(values = trep_box_cols, guide = "none") +
                        scale_colour_manual(values = trep_box_cols, guide = "none") +
                        scale_x_discrete(expand = expansion(mult = 0, add = 0.55)) +
                        labs(title = paste0(dataset_tag,
                                            " \u2014 GAM residual |y\u2212f(x)|: Pre vs Post-IP TREP"),
                             subtitle = sprintf(
                                 "BH FDR p = %.2e | Cohen's d = %.3f | Pre med = %.4f | Post med = %.4f",
                                 r_pval, r_d,
                                 median(pre_trep$resid_abs),
                                 median(post_trep$resid_abs)),
                             x = "", y = "|Residual| from GAM curve") +
                        theme_minimal(base_size = 10) +
                        theme(plot.background = element_rect(fill = "white", colour = NA),
                              panel.background = element_rect(fill = "white", colour = NA),
                              panel.grid.major = element_line(colour = "gray90"),
                              panel.grid.minor = element_blank(),
                              panel.border = element_rect(colour = "black", fill = NA,
                                                          linewidth = 0.35),
                              strip.text = element_text(face = "bold", size = 10),
                              plot.subtitle = element_text(colour = "gray30"),
                              legend.position = "none")

                    show_plot(p_resid_box)
                    save_pp(p_resid_box,
                            file.path(OUT_DIR, paste0(FILE_PREFIX, "_", dataset_tag,
                                "_resid_dispersion_prepost_TREP_boxplot")),
                            pdf_w = 8.4, pdf_h = 4.9,
                            png_w = 1820L, png_h = 1120L, res = 200)
                    cat("  Residual dispersion boxplot saved\n")

                    } else {
                        cat(sprintf(
                            "  Ki dispersion skipped: Pre-IP TREP n=%d, Post-IP TREP n=%d (min 3 each)\n",
                            nrow(pre_trep), nrow(post_trep)))
                    }   # end dispersion size guard

                    # Bulbous_Edge gets an empty sheet when it has no rows.
                    combined_tabs[["Ki_Dispersion"]] <-
                        if (!is.null(disp_df)) disp_df else data.frame()
                    if (is.null(combined_tabs[["Bulbous_Edge"]]))
                        combined_tabs[["Bulbous_Edge"]] <- data.frame()
                    cat("  TREP analysis -> Excel sheets Ki_Dispersion,",
                        "Bulbous_Edge\n\n")

                } else {
                    # Name the reason and the outputs not written.
                    cat("  IP-dependent analysis skipped: ",
                        if (!is.finite(IP_ACTIVE))
                            sprintf("threshold is NA on the %s axis (IPRS not available)",
                                    IP_AXIS_USED)
                        else
                            sprintf("degenerate TREP split (TREP=%d, non-TREP=%d)",
                                    n_trep, n_nontrep),
                        "\n", sep = "")
                    cat("  Not written for ", dataset_tag, ":\n", sep = "")
                    for (.miss in c(
                        paste0(FILE_PREFIX, "_", dataset_tag,
                               "_Combined_GAM_curves_4group.pdf/.png"),
                        paste0(FILE_PREFIX, "_", dataset_tag,
                               "_SKCS_scatter_4group.pdf/.png"),
                        paste0(FILE_PREFIX, "_", dataset_tag,
                               "_SKCS_3D_scatter_EP.html"),
                        paste0(FILE_PREFIX, "_", dataset_tag,
                               "_Ki_prepost_TREP_boxplot.pdf/.png"),
                        paste0(FILE_PREFIX, "_", dataset_tag,
                               "_resid_dispersion_prepost_TREP_boxplot.pdf/.png"),
                        paste0(FILE_PREFIX, "_", dataset_tag,
                               "_diagnostics.xlsx: the Ki_Dispersion and ",
                               "Bulbous_Edge sheets (the workbook is still ",
                               "written, with its other sheets)")))
                        cat("      ", .miss, "\n", sep = "")
                    if (!is.finite(IP_ACTIVE))
                        cat("  IPRS RDS expected at:",
                            "\n    ", iprs_rds_path,
                            "\n",
                            sep = "")
                    cat("\n")
                    # Assign trep_flag_bc for EP CSV even without IP
                    trep_flag_bc <- setNames(
                        ifelse(data_ep$barcode %in% trep_barcodes, "TREP", "non-TREP"),
                        data_ep$barcode
                    )
                }   # end if (is.finite(IP_thresh) && TREP split usable)

                # Export EP values to CSV (always, regardless of IP)
                ep_df <- data.frame(
                    barcode   = data_ep$barcode,
                    region    = data_ep$region,
                    GOI_expr  = data_ep$GOI,
                    S_value   = scores_pos$spear_norm[match(data_ep$barcode, scores_pos$barcode)],
                    K_value   = scores_pos$kend_contrib[match(data_ep$barcode, scores_pos$barcode)],
                    EP_value  = round(ep_values, 6),
                    EP_rank   = rank(-ep_values, ties.method = "first"),
                    is_TREP   = data_ep$barcode %in% trep_barcodes,
                    IP_status = if (is.finite(IP_ACTIVE))
                        ip_status_bc[data_ep$barcode] else NA_character_,
                    group_4   = if (is.finite(IP_ACTIVE))
                        paste(ip_status_bc[data_ep$barcode],
                              trep_flag_bc[data_ep$barcode])
                    else NA_character_,
                    stringsAsFactors = FALSE
                )
                ep_df <- ep_df[order(ep_df$EP_rank), ]

                ep_csv_path <- file.path(OUT_DIR, paste0(FILE_PREFIX, "_",
                                         dataset_tag, "_EP_TREP_values.csv"))
                tryCatch({
                    write.csv(ep_df, ep_csv_path, row.names = FALSE)
                    cat(sprintf("  EP/TREP CSV saved  (%d rows)\n\n", nrow(ep_df)))
                }, error = function(e) {
                    alt_path <- file.path(OUT_DIR, paste0(
                        FILE_PREFIX, "_", dataset_tag, "_EP_TREP_values_",
                        format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"))
                    cat(sprintf("  WARNING: could not write %s (%s)\n",
                                basename(ep_csv_path), conditionMessage(e)))
                    tryCatch({
                        write.csv(ep_df, alt_path, row.names = FALSE)
                        cat(sprintf("  EP/TREP CSV saved to fallback: %s  (%d rows)\n\n",
                                    basename(alt_path), nrow(ep_df)))
                    }, error = function(e2) {
                        cat(sprintf("  ERROR: fallback write also failed (%s)\n\n",
                                    conditionMessage(e2)))
                    })
                })

            } else {
                cat("  EP skipped: res_dfg_d2 not available\n\n")
            }   # end if(!is.null(res_dfg_d2))

        } else {
            cat("  Skipped: n_fixed < MIN_CELLS or no GAM data\n\n")
        }


        # == Workbook export, one per dataset =================================
        # The last three sheets are present only when their sections ran.
        xlsx_path <- file.path(OUT_DIR,
                               paste0(FILE_PREFIX, "_", dataset_tag,
                                      "_diagnostics.xlsx"))
        combined_tabs[["Threshold_Comparison"]]   <- thresh_all
        combined_tabs[["NonSKCS_Stratification"]] <- strat_summary
        .sheet_order <- c("Threshold_Comparison", "NonSKCS_Stratification",
                          "EP_Boundary_Diagnostic", "Bulbous_Edge",
                          "Ki_Dispersion")
        combined_tabs <- combined_tabs[intersect(.sheet_order,
                                                 names(combined_tabs))]
        write_xlsx(combined_tabs, path = xlsx_path)
        cat("  Excel saved:", basename(xlsx_path), "--",
            length(combined_tabs), "sheet(s):",
            paste(names(combined_tabs), collapse = ", "), "\n")
        cat("  Files from older runs, no longer written:",
            "\n")
        for (.old in c("_threshold_comparison.csv", "_threshold_comparison.xlsx",
                       "_TREP_analysis.xlsx", "_bulbous_edge_analysis.csv",
                       "_EP_all_GOI_boundary_diagnostic.csv"))
            cat("      ", paste0(FILE_PREFIX, "_", dataset_tag, .old), "\n",
                sep = "")
        cat("\n")

        invisible(list(
            thresh_all         = thresh_all,
            strat_summary      = strat_summary,
            n_fixed            = n_fixed,
            de_fixed           = de_fixed,
            donor_re_fixed     = if (is.null(res_fixed)) NA else
                                     isTRUE(res_fixed$donor_re),
            n_donors_fixed     = if (is.null(res_fixed)) NA_integer_ else
                                     as.integer(res_fixed$n_donors),
            de_fixed_flat_null = if (is.null(res_fixed) ||
                                     is.null(res_fixed$best_model)) NA_real_
                                 else de_flat_null(res_fixed$best_model) * 100,
            dataset_tag        = dataset_tag,
            sel_fixed          = sel_fixed,
            FIXED_S_THRESH     = FIXED_S_THRESH,
            FIXED_K_THRESH     = FIXED_K_THRESH,
            res_dfg_d2         = res_dfg_d2,
            res_hkg_d2         = res_hkg_d2,
            df_fixed_dfg       = df_fixed_dfg,
            ep_df              = ep_df,
            trep_barcodes      = trep_barcodes,
            nontrep_barcodes   = nontrep_barcodes,
            ip_status_bc       = ip_status_bc,
            trep_flag_bc       = trep_flag_bc,
            IP_thresh          = IP_thresh,
            IP_y_thresh        = IP_y_thresh,
            IP_active          = IP_ACTIVE,
            IP_axis_used       = IP_AXIS_USED,
            scores_pos_aligned = scores_pos,
            spear              = spear,
            kend               = kend,
            is_skcs            = is_skcs
        ))
    }


    # =============================================================================
    # SECTION D: EXECUTE
    # =============================================================================

    # vector() keeps a NULL slot for a failed dataset.
    all_results <- vector("list", length(DATASETS))
    names(all_results) <- DATASETS

    for (ds in DATASETS) {
        all_results[ds] <- list(tryCatch(
            run_p10b_diagnostic(ds),
            error = function(e) {
                cat("  ERROR in dataset", ds, ":", conditionMessage(e), "\n")
                NULL
            }
        ))
    }


    # =============================================================================
    # SECTION E1: POST-ANALYSIS QUALITY GUARDRAILS
    # =============================================================================
    cat("\n", strrep("#", 72), "\n")
    cat("  SECTION E1 - POST-ANALYSIS QUALITY GUARDRAILS\n")
    cat(strrep("#", 72), "\n\n")

    for (ds in names(all_results)) {
        r <- all_results[[ds]]
        if (is.null(r)) {
            gr10b(paste0("GR-06-", ds), "FAIL",
                  "Dataset processing returned NULL")
            next
        }

        n_goi_pos <- r$n_fixed + sum(!r$sel_fixed)
        gr10b(paste0("GR-06-", ds),
              if (n_goi_pos >= MIN_CELLS) "PASS" else "FAIL",
              sprintf("GOI-positive cells (total) = %d  (min = %d)",
                      n_goi_pos, MIN_CELLS))

        gr10b(paste0("GR-07-", ds),
              if (r$n_fixed >= MIN_CELLS) "PASS" else "WARN",
              sprintf("Fixed-threshold cells (S_norm>0.05, K>0.5) = %d  (min = %d)",
                      r$n_fixed, MIN_CELLS))

        de_val <- r$de_fixed
        gr10b(paste0("GR-08-", ds),
              if (!is.na(de_val) && is.finite(de_val)) "PASS" else "WARN",
              sprintf("Fixed-threshold DE = %s",
                      if (is.na(de_val)) "NA (GAM not fitted)"
                      else sprintf("%.2f%%", de_val)))

        # The split axis and threshold go in the ledger; anything other than cDFG fails.
        gr10b(paste0("GR-14-", ds),
              if (identical(r$IP_axis_used, "cDFG")) "PASS" else "FAIL",
              sprintf(paste0("Pre/Post split on %s (IP_CLASS_AXIS = %s) | ",
                             "threshold used = %s | IP_y from IPRS = %s"),
                      if (is.null(r$IP_axis_used)) "NA" else r$IP_axis_used,
                      IP_CLASS_AXIS,
                      if (is.null(r$IP_active) || !is.finite(r$IP_active)) "NA"
                      else sprintf("%.6f", r$IP_active),
                      if (is.null(r$IP_y_thresh) || !is.finite(r$IP_y_thresh)) "NA"
                      else sprintf("%.6f", r$IP_y_thresh)))

        ip_val <- r$IP_thresh
        gr10b(paste0("GR-09-", ds),
              if (!is.na(ip_val) && is.finite(ip_val)) "PASS" else "WARN",
              sprintf("IP_thresh = %s%s",
                      if (is.na(ip_val)) "NA (IPRS not available)"
                      else sprintf("%.6f", ip_val),
                      # The IP definition is part of what the number means.
                      if (is.na(ip_val)) "" else
                          paste0("  (", ip_from, ")")))

        # == Donor random intercept fitted on every model =====================
        .fits_ds <- Filter(function(f) identical(f$dataset, ds),
                           .p10b_fit_registry)
        if (length(.fits_ds) == 0L) {
            gr10b(paste0("GR-10-", ds), "WARN",
                  "no GAM fitted, so donor-RE coverage cannot be assessed")
        } else if (!isTRUE(GAM_DONOR_RE)) {
            gr10b(paste0("GR-10-", ds), "WARN",
                  sprintf(paste0("GAM_DONOR_RE is FALSE: all %d fit(s) are ",
                                 "cell level."),
                          length(.fits_ds)))
        } else {
            .re_fitted <- vapply(.fits_ds, function(f) isTRUE(f$donor_re),
                                 logical(1))
            .n_don_all <- vapply(.fits_ds, function(f)
                if (is.null(f$n_donors)) NA_integer_ else as.integer(f$n_donors),
                integer(1))
            .fallbacks <- vapply(.fits_ds, function(f) f$label,
                                 character(1))[!.re_fitted]
            gr10b(paste0("GR-10-", ds),
                  if (all(.re_fitted)) "PASS" else "WARN",
                  sprintf(paste0("donor RE on %d/%d fit(s); donors per fit ",
                                 "[%s, %s]; %d fallback(s)%s"),
                          sum(.re_fitted), length(.re_fitted),
                          if (all(is.na(.n_don_all))) "NA" else
                              format(min(.n_don_all, na.rm = TRUE)),
                          if (all(is.na(.n_don_all))) "NA" else
                              format(max(.n_don_all, na.rm = TRUE)),
                          length(.fallbacks),
                          if (length(.fallbacks) == 0L) "" else
                              paste0(": ", paste(.fallbacks, collapse = ", "))))
        }
    }

    # == Write final guardrail log ============================================
    gr10b_write()

    cat(sprintf("\n  Final guardrail tally: %d PASS / %d WARN / %d FAIL\n",
                .gr10b_pass, .gr10b_warn, .gr10b_fail))
    if (.gr10b_fail > 0) {
        cat("  Failures present -- review the ledger before interpreting.\n")
    } else if (.gr10b_warn > 0) {
        cat(sprintf("  %d passed, %d warning(s), no failures.\n",
                    .gr10b_pass, .gr10b_warn))
    } else {
        cat(sprintf("  All %d checks passed.\n", .gr10b_pass))
    }


    # =============================================================================
    # SECTION E: CROSS-DATASET SUMMARY
    # =============================================================================
    cat("\n", strrep("=", 72), "\n")
    cat("  Part10b COMPLETE - Cross-dataset summary [", ARM_TAG, "]\n")
    cat(strrep("=", 72), "\n\n")

    summary_rows <- lapply(all_results, function(r) {
        if (is.null(r)) return(NULL)
        r$thresh_all
    })
    summary_all <- do.call(rbind, Filter(Negate(is.null), summary_rows))

    if (!is.null(summary_all) && nrow(summary_all) > 0) {
        write.csv(summary_all,
                  file.path(OUT_DIR, paste0(FILE_PREFIX,
                                            "_summary_all_datasets.csv")),
                  row.names = FALSE)
        cat("  Fixed vs reference comparison across datasets:\n\n")
        print(summary_all, row.names = FALSE)
    }

    cat("\n  Output directory:", OUT_DIR, "\n")
    cat("  Files generated:\n")
    # Files removed at the end of the arm are not listed.
    .transient <- c("_RUN_INCOMPLETE.txt", basename(MSG_PATH_P10B))
    for (f in setdiff(sort(list.files(OUT_DIR)), .transient)) {
        sz <- round(file.size(file.path(OUT_DIR, f)) / 1e3, 0)
        cat(sprintf("    %-58s  %5d KB\n", f, sz))
    }

    # Save all results as one RDS (plots dropped, model environments slimmed).
    cat("\n")
    rds_path <- file.path(OUT_DIR, paste0(FILE_PREFIX, "_all_results.rds"))
    rds_save(all_results, rds_path, label = "Results object")
    cat(strrep("=", 72), "\n")

    cat("\n  Log written to:", LOG_PATH_P10B,
        if (isTRUE(LOG_MESSAGES_P10B)) " (console output and warnings)"
        else " (console output only)", "\n")
    cat(sprintf("\nPart10b SKCS TREP - log ended: %s\n",
                Sys.time()))
    cat("  Part 10b all done:", format(Sys.time()), "\n\n")

    # Stop the cluster; the on.exit handler covers an earlier abort.
    if (!is.null(.p10b_cl)) {
        try(parallel::stopCluster(.p10b_cl), silent = TRUE)
        .p10b_cl <- NULL
        cat("  PSOCK cluster stopped.\n")
    }
    .p10b_completed <- TRUE

    if (sink.number(type = "message") != 2L) sink(type = "message")
    while (sink.number() > 0) sink()

    # The arm finished, so remove the incomplete-run marker.
    suppressWarnings(file.remove(INCOMPLETE_MARKER_P10B))

    invisible(list(arm  = arm_tag, dir = OUT_DIR,
                   pass = .gr10b_pass, warn = .gr10b_warn, fail = .gr10b_fail))
}


# =============================================================================
# SECTION Z: SELF-TESTS
# =============================================================================
# Synthetic data built in memory; no file is read or written. Each function
# takes the arm environment to reach its helpers and cluster.

# == Synthetic frame ==========================================================
# n cells over n_donors donors with unbalanced cell counts, three regions,
# lognormal GOI and a logistic Expression with a per-donor offset (donor_sd).
# S and K are built so that SKCS-core is about 55-60% of cells.
p10b_make_synth <- function(n = 302L, n_donors = 15L,
                            region_n = c(100L, 102L, 100L),
                            min_cells = 5L, max_cells = 40L,
                            donor_sd = 0.45, seed = 42L) {
    set.seed(seed)

    cnt <- as.integer(round(seq(min_cells, max_cells, length.out = n_donors)))
    d   <- as.integer(n) - sum(cnt)
    i   <- 1L
    while (d != 0L) {
        j <- ((i - 1L) %% n_donors) + 1L
        if (d > 0L && cnt[j] < max_cells) { cnt[j] <- cnt[j] + 1L; d <- d - 1L }
        else if (d < 0L && cnt[j] > min_cells) { cnt[j] <- cnt[j] - 1L; d <- d + 1L }
        i <- i + 1L
    }
    donor <- factor(rep(sprintf("D%02d", seq_len(n_donors)), cnt))

    if (sum(region_n) != n)
        region_n <- c(rep(n %/% 3L, 2L), n - 2L * (n %/% 3L))
    # Shuffled so a donor's cells land in more than one region; the region
    # totals are exact because the labels are permuted, not resampled.
    region <- sample(rep(paste0("R", seq_along(region_n)), region_n))

    GOI <- rlnorm(n, 0.8, 0.6)
    Expression <- 3 / (1 + exp(-1.5 * (GOI - median(GOI)))) + rnorm(n, 0, 0.3)
    if (donor_sd > 0)
        Expression <- Expression +
            rnorm(n_donors, 0, donor_sd)[as.integer(donor)]

    dat <- data.frame(barcode    = paste0("SMOKE_", seq_len(n)),
                      GOI        = GOI,
                      Expression = Expression,
                      donor      = donor,
                      region     = region,
                      stringsAsFactors = FALSE)

    sig    <- (dat$Expression - mean(dat$Expression)) *
              (dat$GOI - median(dat$GOI))
    pos_s  <- integer(n); pos_s[order(sig, decreasing = TRUE)] <- seq_len(n)
    sig_k  <- sig + rnorm(n, 0, stats::sd(sig) / 6)
    pos_k  <- integer(n); pos_k[order(sig_k, decreasing = TRUE)] <- seq_len(n)
    n_s    <- round(0.60 * n)
    n_k    <- round(0.65 * n)

    scores_pos <- data.frame(
        barcode      = dat$barcode,
        spear_norm   = 0.05 + 0.30 * (n_s + 0.5 - pos_s) / n_s,
        kend_contrib = 0.50 + 0.40 * (n_k + 0.5 - pos_k) / n_k,
        stringsAsFactors = FALSE)
    scores_pos$is_skcs <- scores_pos$spear_norm > 0.05 &
                          scores_pos$kend_contrib > 0.5

    list(data = dat, scores = scores_pos, donor_counts = cnt)
}

# == Smoke test ===============================================================
p10b_smoke_test <- function(ctx) {

    g   <- function(nm) get(nm, envir = ctx)
    fit <- g(".fit_full_gam_p10b")
    de_explained       <- g("de_explained")
    predict_population <- g("predict_population")
    prep_fit_frame     <- g("prep_fit_frame")
    donor_re_ok        <- g("donor_re_ok")
    get_null_deviance  <- g("get_null_deviance")
    GAM_METHOD_        <- g("GAM_METHOD")

    cat("\n", strrep("=", 72), "\n", sep = "")
    cat("  PART 10b SMOKE TEST (synthetic, no external files)\n")
    cat(strrep("=", 72), "\n", sep = "")

    pass <- 0L; fail <- 0L
    report <- function(ok, label, detail = "") {
        if (isTRUE(ok)) { pass <<- pass + 1L; tag <- "PASS" }
        else            { fail <<- fail + 1L; tag <- "FAIL" }
        cat(sprintf("  [%s] %-24s %s\n", tag, label, detail))
        isTRUE(ok)
    }

    syn  <- p10b_make_synth()
    core <- syn$scores$barcode[syn$scores$is_skcs]
    dat  <- prep_fit_frame(syn$data[syn$data$barcode %in% core, ])
    cat(sprintf("  synthetic: %d cells | %d donors | %d regions | SKCS-core %d (%.1f%%)\n",
                nrow(syn$data), nlevels(syn$data$donor),
                length(unique(syn$data$region)), nrow(dat),
                100 * nrow(dat) / nrow(syn$data)))

    res <- fit(dat, label = "SMOKE", verbose = FALSE)
    if (is.null(res) || is.null(res$best_model)) {
        cat("  [FAIL] fit returned NULL - remaining checks skipped\n")
        cat("SMOKE TEST: 0/5 PASS - 5 FAILED\n")
        return(invisible(FALSE))
    }
    m  <- res$best_model
    fd <- res$data

    # == 1. Formula ===========================================================
    f_txt <- paste(deparse(stats::formula(m)), collapse = " ")
    ok1a  <- grepl("s(donor", f_txt, fixed = TRUE)
    # Collapse the factor to two levels: the RE must be dropped and the
    # donor-aware DE must then coincide with summary()$dev.expl.
    dat2  <- fd
    dat2$donor <- factor(ifelse(as.integer(factor(dat2$donor)) %% 2L == 0L,
                                "A", "B"))
    dat2  <- prep_fit_frame(dat2)
    res2  <- fit(dat2, label = "SMOKE_2donor", verbose = FALSE)
    f2    <- if (is.null(res2)) "" else
             paste(deparse(stats::formula(res2$best_model)), collapse = " ")
    ok1b  <- !is.null(res2) && !grepl("s(donor", f2, fixed = TRUE)
    de2a  <- if (is.null(res2)) NA_real_ else de_explained(res2$best_model, res2$data)
    de2b  <- if (is.null(res2)) NA_real_ else summary(res2$best_model)$dev.expl
    ok1c  <- isTRUE(abs(de2a - de2b) < 1e-8)
    report(ok1a && ok1b && ok1c, "1 formula",
           sprintf("15 donors -> s(donor) %s | 2 donors -> RE dropped %s, DE equal %s",
                   ok1a, ok1b, ok1c))

    # == 2. Deviance explained ================================================
    de_re   <- de_explained(m, fd)
    de_flat <- summary(m)$dev.expl
    report(is.finite(de_re) && abs(de_re - de_flat) >= 0.01, "2 DE estimator",
           sprintf("donor-aware %.4f vs summary()$dev.expl %.4f (diff %.4f)",
                   de_re, de_flat, abs(de_re - de_flat)))

    # == 3. Weights ===========================================================
    per_donor <- tapply(fd$w, fd$donor, sum)
    okw <- isTRUE(all.equal(sum(fd$w), as.numeric(nlevels(droplevels(fd$donor))))) &&
           isTRUE(all.equal(as.numeric(per_donor),
                            rep(1, length(per_donor)))) &&
           isTRUE(all.equal(as.numeric(m$prior.weights), as.numeric(fd$w)))
    old_w <- get("W_DONOR_WEIGHTS", envir = ctx)
    assign("W_DONOR_WEIGHTS", FALSE, envir = ctx)
    unw   <- prep_fit_frame(dat)
    assign("W_DONOR_WEIGHTS", old_w, envir = ctx)
    okw2  <- all(unw$w == 1)
    report(okw && okw2, "3 weights",
           sprintf("sum(w)=%.4f donors=%d, per-donor all 1 | prior.weights match | W=FALSE gives w==1 %s",
                   sum(fd$w), nlevels(droplevels(fd$donor)), okw2))

    # == 4. Population vs conditional prediction ==============================
    # Maximum over all donor levels.
    x_seq <- seq(min(fd$GOI), max(fd$GOI), length.out = 200)
    p_pop <- as.numeric(predict_population(m, x_seq, fd, se = FALSE))
    d_lev <- vapply(levels(fd$donor), function(lv) {
        nd <- data.frame(GOI = x_seq,
                         donor = factor(lv, levels = levels(fd$donor)))
        max(abs(p_pop - as.numeric(predict(m, newdata = nd))))
    }, numeric(1))
    dmax  <- max(d_lev)
    report(is.finite(dmax) && dmax > 0.05, "4 prediction",
           sprintf("max |population - donor-conditional| = %.4f over %d donors (first level %.4f)",
                   dmax, length(d_lev), d_lev[1]))

    # == 5. Null model ========================================================
    use_re    <- donor_re_ok(fd)
    null_re   <- mgcv::gam(Expression ~ 1 + s(donor, bs = "re"),
                           data = fd, weights = w, method = GAM_METHOD_)
    null_flat <- mgcv::gam(Expression ~ 1, data = fd)
    # The weighted flat null isolates the random effect.
    null_flat_w <- mgcv::gam(Expression ~ 1, data = fd, weights = w)
    f_null    <- paste(deparse(stats::formula(null_re)), collapse = " ")
    ok5 <- use_re &&
           grepl("s(donor", f_null, fixed = TRUE) &&
           abs(deviance(null_re) - deviance(null_flat))   > 1e-8 &&
           abs(deviance(null_re) - deviance(null_flat_w)) > 1e-8 &&
           isTRUE(all.equal(get_null_deviance(fd), deviance(null_re)))
    report(ok5, "5 null model",
           sprintf("RE null %.4f | weighted flat null %.4f | unweighted flat null %.4f",
                   deviance(null_re), deviance(null_flat_w),
                   deviance(null_flat)))

    # == EP / TREP block on the same frame, exercised for errors ==============
    null_ep <- if (use_re) null_re else null_flat
    ep <- ifelse((fd$Expression - fitted(null_ep))^2 == 0, NA_real_,
                 1 - ((fd$Expression - fitted(m))^2 /
                      (fd$Expression - fitted(null_ep))^2))
    n_ep <- nrow(fd)
    tgt  <- max(0L, min(as.integer(round(n_ep * de_re)), n_ep))
    cat(sprintf("  EP/TREP: DE=%.4f  target=%d/%d  TREP EP range [%.3f, %.3f]\n",
                de_re, tgt, n_ep,
                if (tgt > 0L) min(sort(ep, decreasing = TRUE)[seq_len(tgt)],
                                  na.rm = TRUE) else NA_real_,
                if (tgt > 0L) max(ep, na.rm = TRUE) else NA_real_))

    if (fail == 0L) cat("SMOKE TEST: 5/5 PASS\n")
    else cat(sprintf("SMOKE TEST: %d/5 PASS - %d FAILED\n", pass, fail))
    invisible(fail == 0L)
}

# == Parallel equivalence =====================================================
p10b_parallel_equivalence <- function(ctx) {
    cl <- get(".p10b_cl", envir = ctx)
    if (is.null(cl) || length(cl) < 2L) {
        cat("PARALLEL EQUIVALENCE: skipped (PAR_CORES_P10B < 2)\n")
        return(invisible(NULL))
    }
    fit            <- get(".fit_full_gam_p10b", envir = ctx)
    de_explained   <- get("de_explained", envir = ctx)
    prep_fit_frame <- get("prep_fit_frame", envir = ctx)

    syn  <- p10b_make_synth()
    core <- syn$scores$barcode[syn$scores$is_skcs]
    dat  <- prep_fit_frame(syn$data[syn$data$barcode %in% core, ])

    r_par <- fit(dat, label = "EQUIV_par", verbose = FALSE, cl = cl)
    r_seq <- fit(dat, label = "EQUIV_seq", verbose = FALSE, cl = NULL)

    if (is.null(r_par) || is.null(r_seq)) {
        cat("PARALLEL EQUIVALENCE: FAILED - one of the two fits returned NULL\n")
        return(invisible(FALSE))
    }

    checks <- list(
        best_k        = identical(r_par$best_k, r_seq$best_k),
        best_prss     = isTRUE(all.equal(r_par$best_prss, r_seq$best_prss,
                                         tolerance = 0)),
        best_iter     = identical(r_par$best_iter, r_seq$best_iter),
        coefficients  = isTRUE(all.equal(coef(r_par$best_model),
                                         coef(r_seq$best_model),
                                         tolerance = 0)),
        de_explained  = isTRUE(all.equal(de_explained(r_par$best_model, r_par$data),
                                         de_explained(r_seq$best_model, r_seq$data),
                                         tolerance = 0)),
        prior_weights = isTRUE(all.equal(as.numeric(r_par$best_model$prior.weights),
                                         as.numeric(r_seq$best_model$prior.weights),
                                         tolerance = 0)))
    # best_k, best_prss and best_iter are reported as one field.
    grouped <- c(selection    = checks$best_k && checks$best_prss &&
                                checks$best_iter,
                 coefficients = checks$coefficients,
                 de_explained = checks$de_explained,
                 prior_weights = checks$prior_weights)
    if (all(grouped)) {
        cat("PARALLEL EQUIVALENCE: 4/4 IDENTICAL\n")
    } else {
        cat(sprintf("PARALLEL EQUIVALENCE: %d/4 IDENTICAL - differ: %s\n",
                    sum(grouped), paste(names(grouped)[!grouped],
                                        collapse = ", ")))
    }
    invisible(all(grouped))
}

# == Stress test ==============================================================
# Sized to the PVM3_all pooled SKCS-core frame: 8,407 cells over 84 donors.
p10b_stress_test <- function(ctx) {
    cl <- get(".p10b_cl", envir = ctx)
    fit            <- get(".fit_full_gam_p10b", envir = ctx)
    de_explained   <- get("de_explained", envir = ctx)
    prep_fit_frame <- get("prep_fit_frame", envir = ctx)

    n <- 8407L; nd <- 84L
    set.seed(42)
    cnt <- sample(20:300, nd, replace = TRUE)
    d   <- n - sum(cnt)
    i   <- 1L
    while (d != 0L) {
        j <- ((i - 1L) %% nd) + 1L
        if (d > 0L && cnt[j] < 400L) { cnt[j] <- cnt[j] + 1L; d <- d - 1L }
        else if (d < 0L && cnt[j] > 5L) { cnt[j] <- cnt[j] - 1L; d <- d + 1L }
        i <- i + 1L
    }
    GOI <- rlnorm(n, 0.8, 0.6)
    don <- factor(rep(sprintf("S%02d", seq_len(nd)), cnt))
    dat <- prep_fit_frame(data.frame(
        barcode    = paste0("STRESS_", seq_len(n)),
        GOI        = GOI,
        Expression = 3 / (1 + exp(-1.5 * (GOI - median(GOI)))) +
                     rnorm(n, 0, 0.3) + rnorm(nd, 0, 0.45)[as.integer(don)],
        donor      = don,
        stringsAsFactors = FALSE))

    cat("\n", strrep("=", 72), "\n", sep = "")
    cat(sprintf("  PART 10b STRESS TEST - %d cells, %d donors\n", n, nd))
    cat(strrep("=", 72), "\n", sep = "")

    t_par <- NA_real_; r_par <- NULL
    if (!is.null(cl) && length(cl) >= 2L) {
        t0 <- proc.time()[["elapsed"]]
        r_par <- try(fit(dat, label = "STRESS_par", verbose = FALSE, cl = cl),
                     silent = TRUE)
        t_par <- proc.time()[["elapsed"]] - t0
    }
    t0 <- proc.time()[["elapsed"]]
    r_seq <- try(fit(dat, label = "STRESS_seq", verbose = FALSE, cl = NULL),
                 silent = TRUE)
    t_seq <- proc.time()[["elapsed"]] - t0

    ok_fit <- !inherits(r_seq, "try-error") && !is.null(r_seq) &&
              (is.null(cl) || length(cl) < 2L ||
               (!inherits(r_par, "try-error") && !is.null(r_par)))
    ok_time <- is.finite(t_par) && t_par > 0 && (t_seq / t_par) > 1.5
    de_s <- if (ok_fit) de_explained(r_seq$best_model, r_seq$data) else NA_real_
    de_p <- if (ok_fit && !is.null(r_par))
        de_explained(r_par$best_model, r_par$data) else de_s
    ok_de <- isTRUE(all.equal(de_p, de_s, tolerance = 1e-10))
    ok_w  <- ok_fit && isTRUE(all.equal(sum(r_seq$data$w), as.numeric(nd)))

    checks <- c(fit_completes = ok_fit, speedup = ok_time,
                de_agrees = ok_de, sum_w = ok_w)
    cat(sprintf("  wall time: parallel %.1f s | sequential %.1f s | ratio %s\n",
                t_par, t_seq,
                if (is.finite(t_par) && t_par > 0)
                    sprintf("%.2fx", t_seq / t_par) else "n/a"))
    cat(sprintf("  DE %.6f | sum(w) %.4f (donors %d)\n",
                de_s, if (ok_fit) sum(r_seq$data$w) else NA_real_, nd))
    if (all(checks)) cat(sprintf("STRESS TEST: 4/4 PASS  (parallel %.1f s, sequential %.1f s)\n",
                                 t_par, t_seq))
    else cat(sprintf("STRESS TEST: %d/4 PASS - failed: %s  (parallel %.1f s, sequential %.1f s)\n",
                     sum(checks), paste(names(checks)[!checks], collapse = ", "),
                     t_par, t_seq))
    invisible(all(checks))
}


# =============================================================================
# SECTION 0c: Driver
# =============================================================================
# The loop, folder creation, failure handling and summary live in CEP_arms.R.

cep_run_arms(PART_DIR_TAG, run_one_arm)
