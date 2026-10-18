# =============================================================================
# Part10_CP10K_SKCS_code.r
# SEA-AD CEP-IP/CMRF, Part 10, SKCS decomposition and TREP classification
#   (Figures 2 and 3, Supplementary Figures S5 and S7, Table S7)
#   For each arm: splits GOI-positive cells into SKCS-core and non-SKCS cells
#   (Spearman-Kendall cell-level score), refits the composite GAM on the
#   SKCS-core cells with a donor random intercept and PRSS selection of k,
#   and labels cells TREP or non-TREP by explanatory power.
#
# Usage (from the SEA-AD folder, or set SEAAD_ROOT), after Part 9
#   source("Part10_CP10K_SKCS_code.r")                          # all arms
#   CEP_ARM <- "PVM1_all"; source("Part10_CP10K_SKCS_code.r")   # one arm
#   Needs CEP_arms.R, R_config.R and R_config_HKG8.R in the same folder.
#
# Outputs (Part10/Part10_<arm>/): log, guardrail log, combined RDS and summary;
#   per dataset (_NeoC and _NonNeoC folders): merged PDF and PNG figures,
#   workbook, SKCS scores CSVs, DE comparison CSV and results RDS.
# =============================================================================


# =============================================================================
# SECTION 0: Libraries and shared arm definitions
# =============================================================================

.p10_rstudio_device <- grDevices::dev.cur()
if (.p10_rstudio_device == 1L) {
    tryCatch({ grDevices::dev.new()
        .p10_rstudio_device <- grDevices::dev.cur() },
        error = function(e) invisible(NULL))
}


library(mgcv)
library(ComplexHeatmap)
library(Matrix)
library(viridis)
library(circlize)
library(grid)
library(ggplot2)
library(gridExtra)
library(scales)
library(writexl)
library(dplyr)
library(qpdf)
library(parallel)

set.seed(42)

# CEP_arms.R holds the arm table, the folder layout and the driver.
# Set CEP_HOME before sourcing to use another folder.
if (!exists("CEP_HOME")) CEP_HOME <- Sys.getenv("SEAAD_ROOT", getwd())
.p10_harness <- file.path(CEP_HOME, "CEP_arms.R")
if (!file.exists(.p10_harness))
    stop("CEP_arms.R not found at: ", .p10_harness,
         call. = FALSE)
source(.p10_harness)

PART_TAG <- "Part10"


# =============================================================================
# SECTION 0d: Pure model helpers and PSOCK workers
# =============================================================================
# Defined at the top level so that PSOCK workers receive them without the arm's
# data; the configuration constants are passed in as arguments.

DONOR_RE_LABEL <- "s(donor)"

p10_gam_formula <- function(k, use_re, bs) {
    rhs <- sprintf("GOI + s(GOI, bs = \"%s\", k = %d)", bs, as.integer(k))
    if (use_re) rhs <- paste0(rhs, " + s(donor, bs = \"re\")")
    stats::as.formula(paste("Expression ~", rhs))
}

p10_goi_smooth_idx <- function(model) {
    if (is.null(model$smooth) || length(model$smooth) == 0)
        return(NA_integer_)
    labs <- vapply(model$smooth, function(s) s$label, character(1))
    i    <- which(labs == "s(GOI)")
    if (length(i) == 1L) return(i)
    i <- which(labs != DONOR_RE_LABEL)
    if (length(i) >= 1L) return(i[1])
    NA_integer_
}

p10_smooth_lambda <- function(model, i) {
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

p10_smooth_stat <- function(sm_summary, label, column) {
    st <- sm_summary$s.table
    if (is.null(st) || nrow(st) == 0) return(NA_real_)
    i <- match(label, rownames(st))
    if (is.na(i)) i <- 1L
    if (!(column %in% colnames(st))) return(NA_real_)
    unname(st[i, column])
}

p10_model_has_re <- function(model)
    !is.null(model$smooth) &&
        any(vapply(model$smooth, function(s) s$label, character(1)) ==
                DONOR_RE_LABEL)

p10_calculate_prss <- function(model) {
    RSS <- sum(stats::residuals(model)^2)
    i   <- p10_goi_smooth_idx(model)
    if (is.na(i))
        return(list(PRSS = RSS, RSS = RSS,
                    f_double_prime_integral = NA_real_, lambda = NA_real_))
    sm     <- model$smooth[[i]]
    lambda <- p10_smooth_lambda(model, i)
    S      <- sm$S[[1]]
    beta   <- stats::coef(model)[sm$first.para:sm$last.para]
    f_dpi  <- as.numeric(t(beta) %*% S %*% beta)
    list(PRSS = RSS + lambda * f_dpi, RSS = RSS,
         f_double_prime_integral = f_dpi, lambda = lambda)
}

p10_extract_reml_iters <- function(model) {
    if (is.null(model) || is.null(model$outer.info))
        return(data.frame(iteration = numeric(0), lambda = character(0),
                          score = numeric(0), stringsAsFactors = FALSE))
    tryCatch({
        n  <- model$outer.info$iter
        df <- data.frame(iteration = seq_len(n), lambda = NA_character_,
                         score = NA_real_, stringsAsFactors = FALSE)
        if (!is.null(model$outer.info$sp)) {
            sp <- model$outer.info$sp
            if (is.matrix(sp)) {
                for (i in seq_len(min(nrow(sp), n)))
                    df$lambda[i] <- paste(format(sp[i, ], digits = 6),
                                          collapse = ", ")
            } else if (is.vector(sp) && length(sp) >= n) {
                df$lambda <- format(sp[seq_len(n)], digits = 6)
            } else if (length(sp) > 0) {
                df$lambda <- format(sp[1], digits = 6)
            }
        }
        if (all(is.na(df$lambda)) && !is.null(model$sp))
            df$lambda <- paste(format(model$sp, digits = 6), collapse = ", ")
        if (!is.null(model$outer.info$score))
            df$score[seq_len(min(length(model$outer.info$score), n))] <-
                model$outer.info$score[
                    seq_len(min(length(model$outer.info$score), n))]
        df
    }, error = function(e) {
        cat("  Warning - extract_reml_iterations:", conditionMessage(e), "\n")
        data.frame(iteration = numeric(0), lambda = character(0),
                   score = numeric(0), stringsAsFactors = FALSE)
    })
}

# == One PRSS iteration, run on a worker ======================================
# Returns summaries only; the master refits the winning model.
p10_prss_worker <- function(job) {
    ctx <- p10_fit_ctx
    res <- tryCatch({
        set.seed(123L + job$i)
        m  <- mgcv::gam(p10_gam_formula(job$k, ctx$use_re, ctx$bs),
                        data = ctx$data, weights = w,
                        method = ctx$method,
                        select = ctx$select, gamma = ctx$gamma)
        pc <- p10_calculate_prss(m)
        ri <- p10_extract_reml_iters(m)
        if (nrow(ri) > 0) {
            ri$prss_iteration <- job$i
            ri$k_value        <- job$k
            ri$prss_score     <- pc$PRSS
            ri$is_best_model  <- FALSE
        }
        sm <- summary(m)
        rs <- if (!is.null(m$outer.info$score)) m$outer.info$score
              else m$gcv.ubre
        list(ok    = TRUE,
             prss  = pc$PRSS,
             reml  = if (nrow(ri) > 0) ri else NULL,
             params = list(
                 iteration     = job$i,
                 k             = job$k,
                 lambda        = pc$lambda,
                 PRSS          = pc$PRSS,
                 RSS           = pc$RSS,
                 f_dpi         = pc$f_double_prime_integral,
                 dev_explained =
                     if (is.finite(ctx$null_dev) && ctx$null_dev > 0)
                         as.numeric((ctx$null_dev - stats::deviance(m)) /
                                        ctx$null_dev)
                     else NA_real_,
                 reml_n        = if (!is.null(m$outer.info$iter))
                                     m$outer.info$iter else 1L,
                 reml_min      = min(rs, na.rm = TRUE),
                 reml_max      = max(rs, na.rm = TRUE),
                 smooth_edf    = p10_smooth_stat(sm, "s(GOI)", "edf"),
                 p_value       = p10_smooth_stat(sm, "s(GOI)", "p-value")))
    }, error = function(e) list(ok = FALSE, message = conditionMessage(e)))
    res$i <- job$i
    res$k <- job$k
    res
}


# =============================================================================
# SECTION 0e: RDS size control
# =============================================================================
# Fitted models are saved with their environments reset to baseenv(), which is
# lossless. Plots are dropped unless RDS_KEEP_PLOTS = TRUE (they are on disk).
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

# Strips every fitted model and drops plots unless RDS_KEEP_PLOTS; names and
# the shape of the object are unchanged.
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
    .p10_preamble <- character(0)

    assign("CEP_ARM", arm_tag, envir = globalenv())
    on.exit(suppressWarnings(rm("CEP_ARM", envir = globalenv())), add = TRUE)

    .p10_preamble <- tryCatch(
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

    .p10_needed <- c("ARM_TAG", "ARM_DIR", "R_EXPORT_ROOT", "GOI", "DFG",
                     "HKG_GENES", "DFG_COMPOSITE_LABEL", "RUN_SUFFIX",
                     "PRIMARY_STYPE", "NORM_TAG", "NORM_TARGET_SUM",
                     "NORM_CEILING")
    .p10_absent <- .p10_needed[!vapply(.p10_needed, exists, logical(1))]
    if (length(.p10_absent) > 0)
        stop("The configuration loaded from\n    ", config_path,
             "\n  does not define: ", paste(.p10_absent, collapse = ", "),
             call. = FALSE)

    if (!identical(ARM_TAG, arm_tag))
        stop("The configuration reports ARM_TAG '", ARM_TAG,
             "' but this arm is '", arm_tag, "'.", call. = FALSE)

    cat("  Config OK: ", ARM_TAG, " | ", GOI, " + ", length(DFG), " DFG | ",
        NORM_TAG, "\n\n", sep = "")

    set.seed(42)

    # == User-editable parameters =============================================
    SKCS_S_THRESH         <- 0.05   # applied to normalised Si (Si_norm)
    SKCS_K_THRESH         <- 0.5
    MIN_CELLS             <- 50L
    GAM_GAMMA_P10       <- 1.5
    PRSS_N_ITER_P10     <- 100L
    PRSS_EARLY_STOP_P10 <- 20L

    # == Estimator ============================================================
    # GAM_DONOR_RE: TRUE fits a donor random intercept; FALSE gives cell-level fits.
    GAM_DONOR_RE        <- TRUE
    DONOR_RE_MIN_LEVELS <- 3L
    # DONOR_RE_LABEL is defined at the top level (Section 0d) for the workers.

    # Inverse-donor-count weights; FALSE gives every cell w = 1.
    W_DONOR_WEIGHTS     <- TRUE

    # HKG8 arm only: extra pooled SKCS scatters with each gene as GOI and the mean
    # of the other eight as the composite. character(0) turns them off.
    ALT_GOI_GENES <- c("EEF1A1", "YWHAZ", "RAB7A")
    ALT_GOI_ARM   <- "HKG8_matched"

    # == Parallelism ==========================================================
    # PAR_CORES_P10 PSOCK workers; 1L runs everything in this process.
    .p10_detected <- suppressWarnings(parallel::detectCores())
    PAR_CORES_P10 <- if (is.na(.p10_detected)) 1L
                     else max(1L, as.integer(.p10_detected) - 1L)
    .p10_cl <- NULL

    # Cap on phase-2 PRSS iterations; the loop stops at 8 + PRSS_PHASE2_MAX_P10.
    PRSS_PHASE2_MAX_P10 <- 5L

    # == Gene definitions =====================================================
    GOI_10 <- GOI
    DFG_10 <- DFG            # unified panel; DFG_NEOC and DFG_NONNC are the same
    # Arm-aware title for the summary text.
    MODULE1_LBL <- if (GOI %in% HKG_GENES) "HKG control panel" else "Module 1"

    # == Labels and output routing (from R_config) ============================
    if (!exists("DFG_COMPOSITE_LABEL")) DFG_COMPOSITE_LABEL <- "cDFG"
    if (!exists("RUN_SUFFIX"))          RUN_SUFFIX <- ""

    # == HKG comparison gene set ==============================================
    # HKG arm: the DFG panel is the HKG module, so the HKG comparison is skipped.
    if (GOI %in% HKG_GENES) {
        HKG_10 <- NULL
        cat("  HKG comparison: SKIPPED (GOI is an HKG gene; DFG = cHKG)\n")
    } else {
        HKG_10 <- HKG_GENES
    }

    # == Display label for the DFG composite in plots =========================
    DFG_DISPLAY_LABEL <- DFG_COMPOSITE_LABEL   # "cDFG" or "cHKG8"

    # == Plot dimensions ======================================================
    PLOT_W_IN <- 4.86
    PLOT_H_IN <- 4.00
    GRID_W_IN <- 14
    GRID_H_IN <- 7

    # == Part 9 inputs ========================================================
    # Part 9 folder and file prefix for this arm, from CEP_arms.R.
    PART9_DIR       <- cep_arm_dir_for("Part9", arm_tag)
    .p9_prefix      <- cep_out_prefix("Part9", NORM_TAG, file_tag)
    PART9_NONNC_RDS <- file.path(PART9_DIR,
                                 paste0(.p9_prefix, "_NonNeoC_GAM_models.rds"))
    PART9_NEOC_RDS  <- file.path(PART9_DIR,
                                 paste0(.p9_prefix, "_NeoC_GAM_models.rds"))

    # == Output directories ===================================================
    # The driver has created out_dir at Part10/Part10_<subdir>/.
    OUT_DIR     <- out_dir
    # Outputs are split by dataset into <arm folder>_NeoC and _NonNeoC folders;
    # arm-level outputs (logs, combined RDS, summary) stay in the arm folder.
    ARM_OUT_DIR <- OUT_DIR
    DS_DIRS     <- character(0)

    # Every output name is built from this prefix.
    OUT_PREFIX  <- cep_out_prefix(PART_TAG, NORM_TAG, file_tag)
    FILE_PREFIX <- paste0(OUT_PREFIX, "_Module1_", gsub("-", "", GOI_10), "_",
                          DFG_COMPOSITE_LABEL)

    NONNC_REGIONS_P10 <- NONNC_REGIONS
    NEOC_REGIONS_P10  <- NEOC_REGIONS

    # == Log file, one per arm ================================================
    # Printed output and the message stream both go to the log; LOG_MESSAGES_P10 =
    # FALSE keeps warnings on the console instead.
    LOG_MESSAGES_P10 <- TRUE
    LOG_PATH_P10 <- file.path(OUT_DIR, paste0(OUT_PREFIX, "_GAM_Log.txt"))
    MSG_PATH_P10 <- file.path(OUT_DIR,
                              paste0(OUT_PREFIX, "_GAM_Log_messages.tmp"))

    while (sink.number() > 0) sink()
    # Only one message diversion can be active, so one conditional unwind is enough.
    if (sink.number(type = "message") != 2L) sink(type = "message")

    .p10_log_con <- file(LOG_PATH_P10, open = "wt")
    .p10_msg_con <- if (isTRUE(LOG_MESSAGES_P10))
        file(MSG_PATH_P10, open = "wt") else NULL
    sink(.p10_log_con, append = FALSE, split = TRUE)
    if (!is.null(.p10_msg_con)) sink(.p10_msg_con, type = "message")
    on.exit({
        tryCatch({
            if (!isTRUE(.p10_completed))
                cat("\n", strrep("!", 72),
                    "\n  ARM DID NOT COMPLETE. See \"*** ARM FAILED\" in the console.\n",
                    strrep("!", 72), "\n", sep = "")
            if (sink.number(type = "message") != 2L) sink(type = "message")
            while (sink.number() > 0) sink()
            if (!is.null(.p10_msg_con)) close(.p10_msg_con)
            close(.p10_log_con)
            # Append the message stream to the log and remove the temporary file.
            if (!is.null(.p10_msg_con) && file.exists(MSG_PATH_P10)) {
                .msg_txt <- readLines(MSG_PATH_P10, warn = FALSE)
                if (length(.msg_txt) > 0)
                    cat(c("", strrep("=", 72),
                          "  MESSAGE STREAM (warnings, messages, errors)",
                          strrep("=", 72), .msg_txt, ""),
                        sep = "\n", file = LOG_PATH_P10, append = TRUE)
                unlink(MSG_PATH_P10)
            }
        }, error = function(e) invisible(NULL))
    }, add = TRUE)

    # warn = 1 prints warnings as they happen, so they reach the log. Restored on exit.
    .p10_old_warn <- options(warn = 1)
    on.exit(options(.p10_old_warn), add = TRUE)

    # Lets the exit handler note in the log that the arm did not complete.
    .p10_completed <- FALSE

    # Replay what was printed before the log existed.
    if (length(.p10_preamble) > 0)
        cat(paste(.p10_preamble, collapse = "\n"), "\n", sep = "")

    cat("\n", strrep("=", 72), "\n")
    cat("  Part 10 -- SKCS Decomposition + Iterative GAM + TREP labels (Module 1)\n")
    cat("  GOI               :", GOI_10, "\n")
    cat("  DFG (", DFG_COMPOSITE_LABEL, ")   :", paste(DFG_10, collapse = ", "), "\n")
    cat("  HKG               :", paste(HKG_10, collapse = ", "), "\n")
    cat("  SKCS thresholds   : S_norm >", SKCS_S_THRESH, "(normalised) AND K >", SKCS_K_THRESH, "\n")
    cat("  Min cells         :", MIN_CELLS, "\n")
    cat("  GAM gamma         :", GAM_GAMMA_P10, "\n")
    cat("  Estimator         :",
        if (GAM_DONOR_RE)
            sprintf("donor random intercept, s(donor, bs = \"re\"), minimum %d donors",
                    DONOR_RE_MIN_LEVELS)
        else "cell level (donor random intercept DISABLED)", "\n")
    cat("  Weights           :",
        if (W_DONOR_WEIGHTS)
            "inverse donor cell count, w = 1/n_i (each donor sums to 1)"
        else "uniform, w = 1 (inverse-donor weighting DISABLED)", "\n")
    cat("  PRSS              : phase 1 = 8 iterations, phase 2 = up to",
        PRSS_PHASE2_MAX_P10, "(early stop", PRSS_EARLY_STOP_P10, ")\n")
    if (PRSS_EARLY_STOP_P10 > PRSS_PHASE2_MAX_P10)
        cat("                      early stop is unreachable at this cap:",
            PRSS_EARLY_STOP_P10, "consecutive non-improving iterations are",
            "needed and phase 2 runs at most", PRSS_PHASE2_MAX_P10, "\n")
    cat("  Alt-GOI scatters  :",
        if (identical(ARM_TAG, ALT_GOI_ARM) && length(ALT_GOI_GENES) > 0)
            paste(ALT_GOI_GENES, collapse = ", ")
        else "none (not the alternate-GOI arm)", "\n")

    # One PSOCK cluster per arm, released at the end or on exit.
    if (PAR_CORES_P10 > 1L) {
        .p10_cl <- parallel::makeCluster(PAR_CORES_P10)  # PSOCK: no fork
        on.exit({
            if (!is.null(.p10_cl))
                try(parallel::stopCluster(.p10_cl), silent = TRUE)
        }, add = TRUE)
        parallel::clusterEvalQ(.p10_cl, { library(mgcv); NULL })
        parallel::clusterExport(.p10_cl,
            c("DONOR_RE_LABEL", "p10_gam_formula", "p10_goi_smooth_idx",
              "p10_smooth_lambda", "p10_smooth_stat", "p10_model_has_re",
              "p10_calculate_prss", "p10_extract_reml_iters",
              "p10_prss_worker"))
        cat("  Parallel          :", PAR_CORES_P10,
            "PSOCK worker(s) for the PRSS iterations\n")
    } else {
        cat("  Parallel          : disabled (PAR_CORES_P10 = 1); sequential\n")
    }
    cat("  Cell source       : per-region gam_data CSVs (full M4-top, GOI > 0)\n")
    cat("  Part 9 RDS NeoC   :", PART9_NEOC_RDS, "\n")
    cat("  Part 9 RDS NonNeoC:", PART9_NONNC_RDS, "\n")
    cat("  Output dir        :", OUT_DIR, "\n")
    cat("  Dataset dirs      :",
        paste0(basename(ARM_OUT_DIR), "_NonNeoC"), "and",
        paste0(basename(ARM_OUT_DIR), "_NeoC"), "\n")
    cat(strrep("=", 72), "\n\n")


    # =============================================================================
    # SECTION A0: Pre-analysis guardrails
    # =============================================================================
    # Results go to the guardrail ledger and <prefix>_Guardrail_Log.txt.
    #   GR-01  norm.mtx max within the CP10K ceiling
    #   GR-02  expm1(norm) row sums pin NORM_TARGET_SUM
    #   GR-04  Raw count matrix integer-like
    #   GR-05  GAM CSV GOI/Expression ranges within ceiling
    #   GR-06  Cell counts match Part 8 expected
    #   GR-07  Part 9 RDS files exist

    GR_LOG_P10 <- file.path(OUT_DIR, paste0(OUT_PREFIX, "_Guardrail_Log.txt"))

    # Marker file written on entry and removed only when the arm completes.
    INCOMPLETE_MARKER <- file.path(OUT_DIR, "_RUN_INCOMPLETE.txt")
    writeLines(c(
        "This arm did not finish.",
        "",
        paste("Arm     :", arm_tag),
        paste("Part    :", PART_TAG),
        paste("Started :", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
        "",
        "Any files in this folder are from a partial run and must not be used.",
        "This file is removed automatically when the arm completes."
    ), INCOMPLETE_MARKER)

    .gr10_log   <- list()
    .gr10_pass  <- 0L
    .gr10_fail  <- 0L
    .gr10_warn  <- 0L

    gr10 <- function(check_id, status, detail) {
        tag <- switch(status, PASS = "[PASS]", FAIL = "[FAIL]", WARN = "[WARN]")
        cat(sprintf("  %s  %-14s  %s\n", tag, check_id, detail))
        .gr10_log[[length(.gr10_log) + 1L]] <<- list(id = check_id, status = status,
                                                       detail = detail)
        if (status == "PASS") .gr10_pass <<- .gr10_pass + 1L
        if (status == "FAIL") .gr10_fail <<- .gr10_fail + 1L
        if (status == "WARN") .gr10_warn <<- .gr10_warn + 1L
    }

    gr10_write <- function() {
        con <- file(GR_LOG_P10, "wt")
        writeLines(c(
            strrep("=", 72),
            paste0("Part 10 -- Guardrail Log [", NORM_TAG, "]"),
            paste("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
            paste("Arm :", ARM_TAG),
            paste0("Normalisation: adata.X = log1p(UMIs/total x ",
                   format(NORM_TARGET_SUM, scientific = FALSE),
                   "), exported unchanged"),
            strrep("=", 72), ""
        ), con)
        for (item in .gr10_log)
            writeLines(sprintf("  [%s]  %-14s  %s", item$status, item$id,
                               item$detail), con)
        writeLines(c("", strrep("=", 72),
                     sprintf("  Total : %d", .gr10_pass + .gr10_fail + .gr10_warn),
                     sprintf("  PASS  : %d", .gr10_pass),
                     sprintf("  WARN  : %d", .gr10_warn),
                     sprintf("  FAIL  : %d", .gr10_fail),
                     if (.gr10_fail == 0 && .gr10_warn == 0)
                         "\n  Result: no failures, no warnings"
                     else if (.gr10_fail == 0)
                         "\n  Result: no failures; review the warnings above"
                     else "\n  Result: FAILURES PRESENT -- review before interpreting",
                     strrep("=", 72)), con)
        close(con)
        cat("  Guardrail log written:", GR_LOG_P10, "\n")
    }

    cat("\n", strrep("=", 72), "\n")
    cat("  SECTION A0 -- PRE-ANALYSIS GUARDRAILS\n")
    cat(strrep("=", 72), "\n\n")

    # Normalisation checks in positive form; Part 9 runs the full battery.

    # Expected per-region counts from the exporter's summary, filtered to the
    # primary supertype; the count check is skipped if regions are still duplicated.
    p10_expected_counts <- function(ds_dir) {
        f <- file.path(ds_dir, "gam_data", "summary_gam_inputs.csv")
        if (!file.exists(f)) return(NULL)
        sm <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
        if (!all(c("region", "n_selected") %in% colnames(sm))) return(NULL)
        stype_col <- intersect(c("supertype", "Supertype", "stype", "cell_type"),
                               colnames(sm))
        if (length(stype_col) > 0) {
            keep <- as.character(sm[[stype_col[1]]]) == PRIMARY_STYPE
            if (!any(keep)) {
                cat("  GR-06 skipped:", basename(f), "has a", stype_col[1],
                    "column with no row for", PRIMARY_STYPE, "\n")
                return(NULL)
            }
            sm <- sm[keep, , drop = FALSE]
        }
        if (anyDuplicated(sm$region) > 0L) {
            cat("  GR-06 skipped:", basename(f),
                "has more than one row per region after filtering\n")
            return(NULL)
        }
        setNames(as.integer(sm$n_selected), sm$region)
    }

    p10_gam_file <- function(ds_dir, region) {
        file.path(ds_dir, "gam_data",
                  paste0(region, "_",
                         gsub("-", "_", gsub("/", "_", PRIMARY_STYPE)),
                         "_selected.csv"))
    }

    # == GR-07: Part 9 RDS files exist ========================================
    gr10("GR-07-NeoC",
         if (file.exists(PART9_NEOC_RDS)) "PASS" else "FAIL",
         sprintf("Part 9 NeoC RDS: %s",
                 if (file.exists(PART9_NEOC_RDS)) basename(PART9_NEOC_RDS)
                 else paste("NOT FOUND at", PART9_NEOC_RDS)))
    gr10("GR-07-NonNeoC",
         if (file.exists(PART9_NONNC_RDS)) "PASS" else "FAIL",
         sprintf("Part 9 NonNeoC RDS: %s",
                 if (file.exists(PART9_NONNC_RDS)) basename(PART9_NONNC_RDS)
                 else paste("NOT FOUND at", PART9_NONNC_RDS)))

    # Part 10 needs the Part 9 RDS; the arm stops here if either is missing.
    if (!file.exists(PART9_NEOC_RDS) || !file.exists(PART9_NONNC_RDS)) {
        gr10_write()
        stop("Part 9 output is missing for arm '", arm_tag, "':\n    ",
             paste(c(if (!file.exists(PART9_NEOC_RDS))  PART9_NEOC_RDS,
                     if (!file.exists(PART9_NONNC_RDS)) PART9_NONNC_RDS),
                   collapse = "\n    "),
             call. = FALSE)
    }

    for (.ds in list(
        list(tag = "NeoC",    dir = NEOC_DIR,  regions = NEOC_REGIONS_P10,
             dfg = DFG_NEOC),
        list(tag = "NonNeoC", dir = NONNC_DIR, regions = NONNC_REGIONS_P10,
             dfg = DFG_NONNC)
    )) {
        cat("  -- Dataset:", .ds$tag, "--------------------------------------\n")

        norm_mtx_path <- file.path(.ds$dir, "deg_support", "norm.mtx")
        raw_mtx_path  <- file.path(.ds$dir, "counts", "raw_selected", "matrix.mtx")

        if (file.exists(norm_mtx_path)) {
            X_n   <- as(readMM(norm_mtx_path), "CsparseMatrix")
            n_max <- max(X_n)

            # == GR-01: maximum within the CP10K ceiling ======================
            gr10(paste0("GR-01-", .ds$tag),
                 if (n_max > 0 && n_max <= NORM_CEILING + 1e-4) "PASS" else "FAIL",
                 sprintf("norm.mtx max = %.4f  (ceiling log1p(%s) = %.4f)",
                         n_max, format(NORM_TARGET_SUM, scientific = FALSE),
                         NORM_CEILING))

            # == GR-02: expm1 row sums pin the target =========================
            # A second normalisation, another log base or a gene-subset export breaks this.
            n_rows <- min(50L, nrow(X_n))
            rs     <- rowSums(expm1(X_n[seq_len(n_rows), , drop = FALSE]))
            rs_dev <- max(abs(rs - NORM_TARGET_SUM))
            gr10(paste0("GR-02-", .ds$tag),
                 if (rs_dev < NORM_TARGET_SUM * 1e-6) "PASS" else "FAIL",
                 sprintf("expm1(norm) row sums: mean = %.6f, max |dev| = %.2e, n = %d",
                         mean(rs), rs_dev, n_rows))

            rm(X_n); gc(verbose = FALSE)
        } else {
            gr10(paste0("GR-01-", .ds$tag), "FAIL",
                 paste("norm.mtx NOT FOUND:", norm_mtx_path))
        }

        # == GR-04: raw counts integer-valued =================================
        if (file.exists(raw_mtx_path)) {
            X_r    <- as(readMM(raw_mtx_path), "CsparseMatrix")
            r_samp <- as.numeric(X_r@x[seq_len(min(10000, length(X_r@x)))])
            is_int <- all(abs(r_samp - round(r_samp)) < 0.01)
            r_max  <- max(X_r)
            gr10(paste0("GR-04-", .ds$tag),
                 if (is_int && r_max > 20) "PASS" else "FAIL",
                 sprintf("raw max = %.0f, integer-like = %s",
                         r_max, if (is_int) "TRUE" else "FALSE"))
            rm(X_r); gc(verbose = FALSE)
        } else {
            gr10(paste0("GR-04-", .ds$tag), "FAIL",
                 paste("raw_selected matrix NOT FOUND:", raw_mtx_path))
        }

        .exp_counts <- p10_expected_counts(.ds$dir)

        # == GR-05: GAM CSV GOI/Expression ranges =============================
        for (reg in .ds$regions) {
            fpath <- p10_gam_file(.ds$dir, reg)
            if (!file.exists(fpath)) next
            df <- read.csv(fpath, stringsAsFactors = FALSE, check.names = FALSE)
            g_max <- max(df$GOI,        na.rm = TRUE)
            e_max <- max(df$Expression, na.rm = TRUE)
            ok    <- g_max <= NORM_CEILING + 1e-4 &&
                     e_max <= NORM_CEILING + 1e-4 && g_max > 0

            gr10(sprintf("GR-05-%s-%s", .ds$tag, reg),
                 if (ok) "PASS" else "FAIL",
                 sprintf("GOI max=%.3f  %s max=%.3f  (ceiling=%.4f)",
                         g_max, DFG_COMPOSITE_LABEL, e_max, NORM_CEILING))

            # == GR-06: cell count vs the exporter's summary ==================
            if (!is.null(.exp_counts) && reg %in% names(.exp_counts)) {
                exp_n <- .exp_counts[[reg]]
                obs_n <- nrow(df)
                gr10(sprintf("GR-06-%s-%s", .ds$tag, reg),
                     if (obs_n == exp_n) "PASS" else "FAIL",
                     sprintf("expected=%d  observed=%d", exp_n, obs_n))
            }
        }
    }

    # == Write pre-analysis guardrail log =====================================
    gr10_write()

    cat(sprintf("\n  Pre-analysis guardrails: %d PASS / %d WARN / %d FAIL\n",
                .gr10_pass, .gr10_warn, .gr10_fail))
    if (.gr10_fail > 0) {
        cat("  *** GUARDRAIL FAILURES -- review before interpreting results ***\n\n")
        warning("Part 10 pre-analysis guardrails: ", .gr10_fail, " FAIL(s)",
                call. = FALSE)
    } else {
        cat("  *** ALL PRE-ANALYSIS GUARDRAILS PASSED ***\n\n")
    }

    # == A0b. Donor random intercept: data, formula, term selection ===========
    # Local copies of the Part 9 helpers.

    # Donor enters every fitting frame as a factor with unused levels dropped.
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

    # Frames with too few donors fall back to the cell-level formula.
    donor_re_ok <- function(x)
        isTRUE(GAM_DONOR_RE) && n_donor_levels(x) >= DONOR_RE_MIN_LEVELS

    # w = 1/(cells this donor contributes to this frame), so each donor's weights
    # sum to 1. Recomputed for every fitted subset.
    donor_weights <- function(donor_fac) {
        if (!isTRUE(W_DONOR_WEIGHTS)) return(rep(1, length(donor_fac)))
        f <- as_donor_factor(donor_fac)
        as.numeric(1 / as.numeric(table(f)[f]))
    }

    # Ensures every fitting frame carries w.
    ensure_weights <- function(data) {
        if (!("w" %in% names(data))) {
            data$w <- if ("donor" %in% names(data))
                donor_weights(data$donor) else rep(1, nrow(data))
        }
        data
    }

    # Wrappers around the Section 0d helpers that supply the configuration constants.
    make_gam_formula <- function(k, use_re) p10_gam_formula(k, use_re, GAM_BS)
    goi_smooth_idx   <- p10_goi_smooth_idx
    smooth_lambda    <- p10_smooth_lambda
    smooth_stat      <- p10_smooth_stat
    model_has_re     <- p10_model_has_re

    # Population-level prediction and SE: donor is set to an observed level and
    # then excluded.
    predict_population <- function(model, x_seq, fit_data, se = TRUE) {
        nd <- data.frame(GOI = x_seq)
        if (model_has_re(model)) {
            lev <- if (!is.null(fit_data) && "donor" %in% names(fit_data))
                as_donor_factor(fit_data$donor)[1] else NULL
            if (is.null(lev)) {
                # Fall back to a donor level stored in the model's var.summary.
                vs <- model$var.summary[["donor"]]
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

    # Null deviance against a donor-aware null, cached by frame.
    .p10_null_dev_cache <- new.env(hash = TRUE, parent = emptyenv())
    get_null_deviance <- function(data) {
        if ("donor" %in% names(data)) data$donor <- as_donor_factor(data$donor)
        data   <- ensure_weights(data)
        use_re <- donor_re_ok(data)
        # sum(w^2) is in the key because the null deviance depends on the weights.
        key <- paste0(nrow(data), "_", round(sum(data$Expression), 6),
                      "_", n_donor_levels(data), "_", as.integer(use_re),
                      "_", round(sum(data$w^2), 9))
        if (exists(key, envir = .p10_null_dev_cache))
            return(get(key, envir = .p10_null_dev_cache))
        set.seed(123)
        # The null carries the same weights as the fitted model.
        nd <- if (use_re)
            deviance(gam(Expression ~ 1 + s(donor, bs = "re"),
                         data = data, weights = w, method = GAM_METHOD))
        else
            deviance(gam(Expression ~ 1, data = data, weights = w))
        assign(key, nd, envir = .p10_null_dev_cache)
        nd
    }

    # Deviance explained against the donor-aware null.
    de_explained <- function(model, data) {
        if (is.null(model) || is.null(data)) return(NA_real_)
        # Without a w column, the model's own prior weights are used for the null.
        if (!("w" %in% names(data)) && !is.null(model$prior.weights) &&
            length(model$prior.weights) == nrow(data))
            data$w <- as.numeric(model$prior.weights)
        nd <- tryCatch(get_null_deviance(data), error = function(e) NA_real_)
        if (!is.finite(nd) || nd <= 0) return(NA_real_)
        as.numeric((nd - deviance(model)) / nd)
    }

    # Weights for the EP ratio and the TREP count: the frame's w, else the model's
    # prior weights, else 1.
    ep_weight_vec <- function(model, data) {
        if ("w" %in% names(data)) return(as.numeric(data$w))
        if (!is.null(model) && !is.null(model$prior.weights) &&
            length(model$prior.weights) == nrow(data))
            return(as.numeric(model$prior.weights))
        rep(1, nrow(data))
    }

    # Weighted population-level deviance explained (1 - RSS_w / TSS_w), from the
    # same donor-free prediction that the EP ranking uses; it sizes the TREP set.
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

    # Fraction that sizes the TREP count, clamped to [0, 1].
    trep_count_fraction <- function(de) {
        if (!is.finite(de)) return(NA_real_)
        max(0, min(1, de))
    }

    # == A1. PRSS =============================================================
    # Penalty, coefficients and smoothing parameter come from the s(GOI) smooth.
    calculate_prss <- function(model, data = NULL) p10_calculate_prss(model)

    # == A2. REML iteration extraction ========================================
    extract_reml_iterations_enhanced <- p10_extract_reml_iters

    # == A3. Two-phase PRSS GAM ===============================================
    # cl defaults to the arm's cluster.
    fit_gam_prss_p10 <- function(data, region_id, gene_set_label,
                                 num_iterations = PRSS_N_ITER_P10,
                                 gamma          = GAM_GAMMA_P10,
                                 cl             = .p10_cl) {
        set.seed(123)
        unique_x <- length(unique(data$GOI))
        max_k    <- min(10L, unique_x - 1L)

        # Donor is levelled once so every iteration and the null see the same factor.
        if ("donor" %in% names(data)) data$donor <- as_donor_factor(data$donor)
        # Weights are recomputed from this frame's own donors.
        data$w <- if ("donor" %in% names(data)) donor_weights(data$donor)
                  else rep(1, nrow(data))
        n_don  <- n_donor_levels(data$donor)
        use_re <- donor_re_ok(data$donor)

        cat("  [", region_id, "|", gene_set_label, "]",
            "cells:", nrow(data), "| unique GOI:", unique_x,
            "| max k:", max_k, "| donors:", n_don,
            "| donor RE:", if (use_re) "YES" else "NO",
            "| sum(w):", round(sum(data$w), 4), "\n")

        if (isTRUE(GAM_DONOR_RE) && !use_re)
            cat("  Donor random intercept not fitted: ", n_don,
                " donor level(s), minimum is ", DONOR_RE_MIN_LEVELS,
                ". This fit is cell level.\n", sep = "")

        if (unique_x < 4L) {
            cat("  Skipping: fewer than 4 unique GOI values\n")
            return(NULL)
        }
    
        best_prss <- Inf; best_model <- NULL; best_k <- NULL
        best_iter <- 0L; prss_values <- numeric(num_iterations)
        model_params <- list(); reml_iters_list <- list()
        early_stop_counter <- 0L; early_stop_triggered <- FALSE
        # Iterations that actually ran, successfully or not.
        iters_run <- integer(0)
        # Set when the winning iteration was fitted on a worker.
        best_par_iter <- NULL

        # The k window around the incumbent best.
        k_window <- function(bk)
            if (bk == max_k) c(max_k - 1L, max_k) else
            if (bk == 3L)    c(3L, 4L) else
            c(bk - 1L, bk, bk + 1L)

        run_one <- function(i, k) {
            iters_run <<- c(iters_run, i)
            set.seed(123 + i)
            m <- gam(make_gam_formula(k, use_re),
                     data = data, weights = w,
                     method = GAM_METHOD,
                     select = GAM_SELECT, gamma = gamma)
            pc <- calculate_prss(m, data)
            ri <- extract_reml_iterations_enhanced(m)
            if (nrow(ri) > 0) {
                ri$prss_iteration <- i; ri$k_value <- k
                ri$prss_score <- pc$PRSS; ri$is_best_model <- FALSE
                reml_iters_list[[i]] <<- ri
            }
            prss_values[i] <<- pc$PRSS
            rs_i <- if (!is.null(m$outer.info$score)) m$outer.info$score
            else m$gcv.ubre
            rn_i <- if (!is.null(m$outer.info$iter))  m$outer.info$iter else 1L
            sm_i <- summary(m)
            model_params[[i]] <<- list(
                iteration     = i, k = k, lambda = pc$lambda,
                PRSS          = pc$PRSS, RSS = pc$RSS,
                f_dpi         = pc$f_double_prime_integral,
                dev_explained = de_explained(m, data),
                reml_n        = rn_i,
                reml_min      = min(rs_i, na.rm = TRUE),
                reml_max      = max(rs_i, na.rm = TRUE),
                # s(GOI) row of the smooth table.
                smooth_edf    = smooth_stat(sm_i, "s(GOI)", "edf"),
                p_value       = smooth_stat(sm_i, "s(GOI)", "p-value"))
            list(model = m, prss = pc$PRSS, k = k)
        }
    
        # == Iteration bookkeeping shared by both phases ======================
        # Every result is applied in index order, so sequential and parallel runs agree.

        run_local <- function(i, k) {
            tryCatch({
                res <- run_one(i, k)
                list(ok = TRUE, i = i, k = k, prss = res$prss,
                     model = res$model, local = TRUE)
            }, error = function(e)
                list(ok = FALSE, i = i, k = k, local = TRUE,
                     message = conditionMessage(e)))
        }

        apply_result <- function(i, r, phase) {
            if (!isTRUE(r$local)) iters_run <<- c(iters_run, i)
            if (!isTRUE(r$ok)) {
                cat("    error iter", i, ":", r$message, "\n")
                prss_values[i]    <<- NA
                model_params[[i]] <<- NULL
                if (phase == 2L)
                    early_stop_counter <<- early_stop_counter + 1L
                return(FALSE)
            }
            # run_one() has already written these for a local result.
            if (!isTRUE(r$local)) {
                prss_values[i]    <<- r$prss
                model_params[[i]] <<- r$params
                if (!is.null(r$reml)) reml_iters_list[[i]] <<- r$reml
            }
            cat(if (phase == 1L) "    init iter" else "    refine iter",
                i, "| k =", r$k, "| PRSS:", round(r$prss, 4), "\n")
            if (r$prss < best_prss) {
                best_prss <<- r$prss; best_k <<- r$k; best_iter <<- i
                early_stop_counter <<- 0L
                if (!is.null(r$model)) {
                    best_model    <<- r$model
                    best_par_iter <<- NULL
                } else {
                    # Fitted on a worker; the master refits the winner at the end.
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
            ctx <- list(data     = data,
                        bs       = GAM_BS,
                        method   = GAM_METHOD,
                        select   = GAM_SELECT,
                        gamma    = gamma,
                        use_re   = use_re,
                        null_dev = tryCatch(get_null_deviance(data),
                                            error = function(e) NA_real_))
            p10_fit_ctx <- ctx
            parallel::clusterExport(cl, "p10_fit_ctx", envir = environment())
        }

        dispatch <- function(idx, ks) {
            if (!use_cluster)
                return(lapply(seq_along(idx),
                              function(j) run_local(idx[j], ks[j])))
            jobs <- Map(function(j, k) list(i = j, k = as.integer(k)), idx, ks)
            parallel::parLapply(cl, jobs, p10_prss_worker)
        }

        # == Phase 1: exploration, k = 3..max_k ===============================
        # Independent iterations, dispatched as one batch and applied in index order.
        n_phase1 <- min(8L, num_iterations)
        if (n_phase1 > 0L) {
            idx1 <- seq_len(n_phase1)
            ks1  <- 3L + ((idx1 - 1L) %% (max_k - 2L))
            out1 <- dispatch(idx1, ks1)
            for (j in seq_along(idx1)) apply_result(idx1[j], out1[[j]], 1L)
        }

        if (best_iter == 0L || is.null(best_k)) {
            cat("  WARNING: all phase-1 iterations failed for",
                region_id, gene_set_label, "-- aborting\n")
            return(NULL)
        }
        k_range <- k_window(best_k)
        cat("    best k after phase 1:", best_k,
            "| PRSS:", round(best_prss, 4), "\n")

        # == Phase 2: adaptive refinement =====================================
        # Runs to 8 + PRSS_PHASE2_MAX_P10. In parallel, a chunk is run assuming no
        # improvement; if one improves, the remaining iterations run sequentially.
        phase2_end <- min(num_iterations, 8L + as.integer(PRSS_PHASE2_MAX_P10))

        phase2_serial <- function(from_i) {
            if (phase2_end < 9L || from_i > phase2_end)
                return(invisible(NULL))
            for (i in from_i:phase2_end) {
                if (early_stop_counter >= PRSS_EARLY_STOP_P10) {
                    cat("    early stopping at iter", i, "\n")
                    early_stop_triggered <<- TRUE; break
                }
                if (i %% 3L == 0L) set.seed(200L + i)
                k <- if (i %% 3L == 0L) sample(k_range, 1) else best_k
                apply_result(i, run_local(i, k), 2L)
            }
            invisible(NULL)
        }

        phase2_parallel <- function() {
            n_workers <- max(1L, length(cl))
            i <- 9L
            while (i <= phase2_end) {
                if (early_stop_counter >= PRSS_EARLY_STOP_P10) {
                    cat("    early stopping at iter", i, "\n")
                    early_stop_triggered <<- TRUE
                    return(invisible(NULL))
                }
                idx <- i:min(i + n_workers - 1L, phase2_end)
                ks  <- vapply(idx, function(j) {
                    if (j %% 3L == 0L) {
                        set.seed(200L + j)
                        as.integer(sample(k_range, 1))
                    } else as.integer(best_k)
                }, integer(1))
                out <- dispatch(idx, ks)

                for (jj in seq_along(idx)) {
                    j <- idx[jj]
                    if (early_stop_counter >= PRSS_EARLY_STOP_P10) {
                        cat("    early stopping at iter", j, "\n")
                        early_stop_triggered <<- TRUE
                        return(invisible(NULL))
                    }
                    if (apply_result(j, out[[jj]], 2L) && j < phase2_end) {
                        cat("    phase-2 speculation invalidated at iter", j,
                            "-- remaining iterations run sequentially\n")
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

        # Refit the winning model here if it was fitted on a worker.
        if (!is.null(best_par_iter) && identical(best_iter, best_par_iter$i)) {
            cat("    best model is from iter", best_iter,
                "(k =", best_par_iter$k, ") on a worker -- refitting here\n")
            set.seed(123L + best_iter)
            best_model <- gam(make_gam_formula(best_par_iter$k, use_re),
                              data = data, weights = w,
                              method = GAM_METHOD,
                              select = GAM_SELECT, gamma = gamma)
        }

        if (length(iters_run) == 0L ||
            all(is.na(prss_values[iters_run]))) {
            cat("  WARNING: all PRSS iterations failed for",
                region_id, gene_set_label, "\n")
            return(NULL)
        }
        actual_iters <- max(iters_run)
        valid_reml   <- Filter(Negate(is.null), reml_iters_list)
        all_reml     <- if (length(valid_reml) > 0)
            do.call(rbind, valid_reml) else NULL
        if (!is.null(all_reml) && nrow(all_reml) > 0)
            all_reml$is_best_model[all_reml$prss_iteration == best_iter] <- TRUE
    
        list(best_model     = best_model,
             best_iter      = best_iter,
             best_k         = best_k,
             best_prss      = best_prss,
             donor_re       = use_re,
             n_donors       = n_don,
             prss_values    = prss_values[seq_len(actual_iters)],
             model_params   = model_params[seq_len(actual_iters)],
             all_reml_iters = all_reml,
             early_stopped  = early_stop_triggered,
             actual_iters   = actual_iters,
             region         = region_id,
             gene_set       = gene_set_label,
             data           = data)
    }

    # == A5. Performance row ==================================================
    build_perf_row_p10 <- function(result,
                                   group_label = "SEAAD_P10_Module1_SKCS") {
        m       <- result$best_model
        pc      <- calculate_prss(m, result$data)
        sm      <- summary(m)
        has_re  <- model_has_re(m)
        data.frame(
            Region          = result$region,
            Group           = group_label,
            Gene_Set        = result$gene_set,
            PRSS            = pc$PRSS,
            RSS             = pc$RSS,
            Lambda_REML     = pc$lambda,
            k               = result$best_k,
            DE              = de_explained(m, result$data),
            DE_vs_flat_null = sm$dev.expl,
            Donor_RE        = has_re,
            N_Donors        = if (is.null(result$n_donors))
                                  n_donor_levels(result$data)
                              else result$n_donors,
            Donor_RE_EDF    = if (has_re)
                                  smooth_stat(sm, DONOR_RE_LABEL, "edf")
                              else NA_real_,
            Smooth_EDF      = smooth_stat(sm, "s(GOI)", "edf"),
            Smooth_p_value  = smooth_stat(sm, "s(GOI)", "p-value"),
            Actual_iters    = result$actual_iters,
            stringsAsFactors = FALSE)
    }


    # =============================================================================
    # SECTION A5: Plot helpers
    # =============================================================================

    theme_p10 <- function(base_size = 9) {
        theme_minimal(base_size = base_size) +
            theme(
                plot.title       = element_text(face = "bold", size = base_size + 1),
                plot.subtitle    = element_text(size = base_size - 1),
                axis.title       = element_text(face = "bold", size = base_size - 1),
                axis.text        = element_text(size = base_size - 2),
                legend.title     = element_text(size = base_size - 1),
                legend.text      = element_text(size = base_size - 2),
                panel.grid.major = element_line(color = "gray90"),
                panel.grid.minor = element_line(color = "gray95"),
                panel.border     = element_rect(color = "black", fill = NA,
                                                linewidth = 0.35),
                plot.background  = element_rect(fill = "white", color = NA),
                panel.background = element_rect(fill = "white", color = NA))
    }

    PRSS_PLOT_MAX_ITER_P10 <- 28L

    make_prss_plot_p10 <- function(result, compact = FALSE) {
        bs    <- if (compact) 7 else 9
        n_cap <- min(PRSS_PLOT_MAX_ITER_P10, result$actual_iters)
        iters <- seq_len(n_cap)
        prss_v <- result$prss_values[iters]
        df <- data.frame(Iteration = iters, PRSS = prss_v,
                         Is_Best   = iters == result$best_iter)
        df <- df[!is.na(df$PRSS), ]
    
        # Guard against empty df after zero-threshold truncation
        if (nrow(df) > 0) {
            prss_max    <- max(df$PRSS, na.rm = TRUE)
            zero_thresh <- max(prss_max * 0.01, 1e-6)
            first_zero  <- which(df$PRSS < zero_thresh)[1]
            if (!is.na(first_zero) && first_zero > 1L)
                df <- df[seq_len(first_zero - 1L), ]
        }
        if (nrow(df) == 0)
            return(ggplot() +
                       annotate("text", x = 0.5, y = 0.5, size = 3,
                                label = paste(result$region, result$gene_set,
                                              "PRSS N/A")) +
                       theme_void())
    
        best_row  <- df[df$Is_Best, ]
        beyond_28 <- result$best_iter > PRSS_PLOT_MAX_ITER_P10
        if (nrow(best_row) == 0) {
            best_row  <- df[which.min(df$PRSS), ]
            beyond_28 <- TRUE
        }
        x_max <- max(df$Iteration, na.rm = TRUE)
    
        ggplot(df, aes(x = Iteration, y = PRSS)) +
            geom_line(linewidth = 0.6, color = "gray30") +
            geom_point(aes(color = Is_Best, size = Is_Best)) +
            geom_vline(xintercept = best_row$Iteration[1],
                       linetype  = "dashed",
                       color     = if (beyond_28) "gray50" else "firebrick",
                       linewidth = 0.5) +
            geom_point(data = best_row, aes(x = Iteration, y = PRSS),
                       color = if (beyond_28) "gray40" else "firebrick",
                       size = 3.5, shape = 18, inherit.aes = FALSE) +
            scale_color_manual(values = c("FALSE" = "gray20", "TRUE" = "firebrick"),
                               guide = "none") +
            scale_size_manual( values = c("FALSE" = 1.2,     "TRUE" = 3.5),
                               guide = "none") +
            scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
            coord_cartesian(xlim = c(1, x_max)) +
            labs(title    = paste0(result$region, " | ", result$gene_set,
                                   " (Module1 SKCS-core) -- PRSS"),
                 subtitle = paste0("Best iter: ", result$best_iter,
                                   "  k=", result$best_k,
                                   "  PRSS=", round(result$best_prss, 3)),
                 x = "PRSS Iteration", y = "PRSS") +
            theme_p10(base_size = bs)
    }

    make_reml_plot_p10 <- function(result, compact = FALSE) {
        bs     <- if (compact) 7 else 9
        all_ri <- result$all_reml_iters
        if (is.null(all_ri) || nrow(all_ri) == 0)
            return(ggplot() +
                       annotate("text", x = 0.5, y = 0.5, size = 3,
                                label = "REML data unavailable") +
                       labs(title = paste0(result$region, " | ", result$gene_set,
                                           " (Module1 SKCS-core) -- REML")) +
                       theme_void())
        best_ri <- all_ri[all_ri$is_best_model == TRUE, ]
        if (nrow(best_ri) == 0)
            best_ri <- all_ri[all_ri$prss_iteration == result$best_iter, ]
        if (nrow(best_ri) == 0) best_ri <- all_ri
        best_ri   <- best_ri[order(best_ri$iteration), ]
        final_row <- best_ri[nrow(best_ri), ]
        ggplot(best_ri, aes(x = iteration, y = score)) +
            geom_line(linewidth = 0.6, color = "gray30") +
            geom_point(size = 2, color = "gray20") +
            geom_point(data = final_row, aes(x = iteration, y = score),
                       color = "firebrick", size = 3.5, shape = 18,
                       inherit.aes = FALSE) +
            labs(title    = paste0(result$region, " | ", result$gene_set,
                                   " (Module1 SKCS-core) -- REML"),
                 subtitle = paste0("PRSS iter ", result$best_iter,
                                   "  |  Final REML = ",
                                   round(final_row$score, 3)),
                 x = "REML Iteration", y = "REML Score") +
            theme_p10(base_size = bs)
    }

    # == SKCS scatter =========================================================
    make_skcs_scatter_p10 <- function(spear_norm, kend_contrib, is_skcs,
                                      label_title, global_rho, global_tau,
                                      n_skcs, n_nskcs,
                                      panel_lbl = MODULE1_LBL) {
        scat_df <- data.frame(
            Spearman_norm   = spear_norm,
            Kendall_contrib = kend_contrib,
            Group = ifelse(is_skcs, "SKCS-core", "Non-SKCS"),
            stringsAsFactors = FALSE)
        scat_df$draw_order <- ifelse(is_skcs, 2L, 1L)
        scat_df <- scat_df[order(scat_df$draw_order), ]
        n_total <- n_skcs + n_nskcs
    
        p <- ggplot(scat_df, aes(x = Spearman_norm, y = Kendall_contrib)) +
            geom_vline(xintercept = SKCS_S_THRESH, colour = "#00B050",
                       linetype = "solid", linewidth = 0.4) +
            geom_hline(yintercept = SKCS_K_THRESH, colour = "#00B050",
                       linetype = "solid", linewidth = 0.4) +
            geom_point(data   = subset(scat_df, Group == "Non-SKCS"),
                       colour = "#707070", alpha = 0.3, size = 0.9) +
            geom_point(data   = subset(scat_df, Group == "SKCS-core"),
                       colour = "#4B0082", alpha = 0.35, size = 1.2) +
            labs(title    = paste0(label_title, " -- ", panel_lbl,
                                   " SKCS score decomposition"),
                 subtitle = sprintf(
                     "rho=%.3f  tau=%.3f  |  SKCS-core: n=%d (S_norm>%.2f, K>%.1f)  |  Non-SKCS: n=%d",
                     global_rho, global_tau, n_skcs,
                     SKCS_S_THRESH, SKCS_K_THRESH, n_nskcs),
                 x = "Normalised Spearman contribution (Si / Si_max)",
                 y = "Kendall concordance contribution") +
            theme_p10(base_size = 10) +
            theme(legend.position  = "none",
                  plot.title        = element_text(face = "bold", size = 11),
                  plot.subtitle     = element_text(size = 8))
        p
    }

    `%||%` <- function(a, b) if (!is.null(a)) a else b

    # == Device-aware pane helper =============================================
    safe_dev_set_p10 <- function() {
        tryCatch(grDevices::dev.set(.p10_rstudio_device),
                 error = function(e) invisible(NULL))
    }

    emit_plot_pane <- function(p) {
        if (is.null(p)) return(invisible(NULL))
        safe_dev_set_p10()
        tryCatch(print(p), error = function(e) invisible(NULL))
    }

    # == CEP-IP scatter builder ===============================================
    # Cells coloured TREP (top round(n * DE) by explanatory power) or non-TREP,
    # with the population-level GAM fit and 95% band.
    make_cep_ip_scatter_p10 <- function(gam_model, sample_data,
                                        dataset_tag, panel_tag,
                                        goi_name = GOI_10) {
        if (is.null(gam_model) || is.null(sample_data) ||
            !all(c("GOI", "Expression") %in% names(sample_data)))
            return(NULL)
    
        set.seed(123)
        # EP per cell: 1 - model squared deviation / flat-null squared deviation, as in
        # the published CEP-IP method. Predictions exclude s(donor).
        null_model   <- gam(Expression ~ 1, data = sample_data)
        null_fitted  <- fitted(null_model)
        model_fitted <- as.numeric(predict(gam_model, newdata = sample_data,
                                           exclude = DONOR_RE_LABEL,
                                           type = "response"))
        null_sq_diff  <- (sample_data$Expression - null_fitted)^2
        model_sq_diff <- (sample_data$Expression - model_fitted)^2
        ep_value      <- ifelse(null_sq_diff > 0,
                                1 - model_sq_diff / null_sq_diff, 0)
    
        # The split size comes from the same population-level prediction.
        model_dev_explained <- pop_dev_explained(gam_model, sample_data,
                                                 model_fitted)
        if (!is.finite(model_dev_explained))
            model_dev_explained <- summary(gam_model)$dev.expl
        target_cells <- max(1L, min(
            round(nrow(sample_data) *
                      trep_count_fraction(model_dev_explained)),
            nrow(sample_data)))
        sorted_idx   <- order(ep_value, decreasing = TRUE)
        trep_idx     <- sorted_idx[seq_len(target_cells)]
    
        x_range  <- range(sample_data$GOI,        na.rm = TRUE)
        y_range  <- range(sample_data$Expression, na.rm = TRUE)
        ABS_PAD  <- 0.08
    
        pd_cep <- data.frame(
            GOI        = sample_data$GOI,
            Expression = sample_data$Expression,
            Group      = ifelse(seq_len(nrow(sample_data)) %in% trep_idx,
                                "TREP", "non-TREP"),
            stringsAsFactors = FALSE)
        pd_cep$draw_order <- ifelse(pd_cep$Group == "TREP", 2L, 1L)
        pd_cep <- pd_cep[order(pd_cep$draw_order), ]
    
        set.seed(123)
        # Population-level curve, random intercept excluded.
        x_seq_cep <- seq(x_range[1], x_range[2], length.out = 1000L)
        pred_x    <- data.frame(GOI = x_seq_cep)
        pg_pred   <- predict_population(gam_model, x_seq_cep, sample_data)
        pred_x$fit    <- as.numeric(pg_pred$fit)
        pred_x$se.fit <- as.numeric(pg_pred$se.fit)
    
        dev_pct <- round(model_dev_explained * 100, 2)
    
        ggplot() +
            geom_point(data = pd_cep,
                       aes(x = GOI, y = Expression, color = Group),
                       size = 1.8, alpha = 0.45) +
            geom_line(data = pred_x, aes(x = GOI, y = fit),
                      color = "#FFCC99", linewidth = 1.2, na.rm = TRUE) +
            geom_ribbon(data = pred_x,
                        aes(x    = GOI,
                            ymin = fit - 1.96 * se.fit,
                            ymax = fit + 1.96 * se.fit),
                        fill = "#FFCC99", alpha = 0.2, na.rm = TRUE) +
            scale_color_manual(
                values = c("TREP"     = "#4B0082",
                           "non-TREP" = "#C0C0C0"),
                name   = "CEP-IP") +
            theme_p10(base_size = 8) +
            theme(legend.position  = "right",
                  legend.text      = element_text(size = 7),
                  legend.title     = element_text(size = 8, face = "bold"),
                  legend.key.size  = unit(0.4, "cm"),
                  plot.title       = element_text(face = "bold", size = 9),
                  plot.subtitle    = element_text(size = 7)) +
            labs(title    = paste0(goi_name, " vs ", DFG_COMPOSITE_LABEL,
                                   " Expression -- ",
                                   dataset_tag, " Pooled  [", panel_tag, "]"),
                 subtitle = paste0("TREP sizing DE (flat null, population fit): ", dev_pct, "%  |  n = ",
                                   nrow(sample_data),
                                   " (TREP = ", target_cells,
                                   ", non-TREP = ",
                                   nrow(sample_data) - target_cells, ")"),
                 x = paste0(goi_name, " Expression (", NORM_TAG, " log1p)"),
                 y = paste0(DFG_COMPOSITE_LABEL, " Expression (", NORM_TAG,
                            " log1p)")) +
            scale_x_continuous(
                labels = number_format(accuracy=0.1),
                limits = c(x_range[1] - ABS_PAD, x_range[2] + ABS_PAD)) +
            scale_y_continuous(
                labels = number_format(accuracy=0.1),
                limits = c(max(0, y_range[1] - ABS_PAD),
                           y_range[2] + ABS_PAD))
    }

    # == File-save helper (PDF + PNG) =========================================
    save_plot_pdf_png <- function(p, base_path,
                                  pdf_w = 14, pdf_h = 7,
                                  png_w = 2800, png_h = 1400, png_res = 200) {
        tryCatch({
            pdf(paste0(base_path, ".pdf"), width = pdf_w, height = pdf_h)
            print(p); dev.off()
        }, error = function(e) {
            cat("    PDF error:", basename(base_path), "--",
                conditionMessage(e), "\n")
            try(dev.off(), silent = TRUE)
        })
        tryCatch({
            png(paste0(base_path, ".png"), width = png_w, height = png_h,
                res = png_res)
            print(p); dev.off()
        }, error = function(e) {
            cat("    PNG error:", basename(base_path), "--",
                conditionMessage(e), "\n")
            try(dev.off(), silent = TRUE)
        })
        safe_dev_set_p10()
    }

    safe_write_csv <- function(x, path, ...) {
        tryCatch(
            write.csv(x, path, ...),
            error = function(e)
                cat("    CSV write error:", basename(path), "--",
                    conditionMessage(e), "\n"))
    }

    # == Cell count audit sheet ===============================================
    build_cell_count_sheet_p10 <- function(tag, cell_audit_list,
                                           pooled_region_name) {
        rows <- list()
        for (nm in names(cell_audit_list)) {
            ca  <- cell_audit_list[[nm]]
            # Skipped regions return NA (neither pass nor fail).
            ok <- if (ca$n_used == 0L && !isTRUE(ca$n_goi_pos == 0L)) {
                NA   # skipped -- not a mismatch, not a pass
            } else {
                isTRUE(ca$n_goi_pos == ca$n_used)
            }
            rows[[length(rows) + 1L]] <- data.frame(
                Tag                 = tag,
                Region              = nm,
                Is_Pooled           = ca$is_pooled,
                N_Raw               = ca$n_raw,
                N_GOI_Neg           = ca$n_goi_neg,
                N_GOI_Pos           = ca$n_goi_pos,
                N_Used_for_Analysis = ca$n_used,
                N_Donors            = if (is.null(ca$n_donors)) NA_integer_
                                      else ca$n_donors,
                Sum_Weights         = if (is.null(ca$sum_w)) NA_real_
                                      else ca$sum_w,
                Tally_OK            = ok,
                Note                = ca$note,
                stringsAsFactors    = FALSE)
        }
        if (length(rows) == 0)
            return(data.frame(Note = "No audit data", stringsAsFactors = FALSE))
        do.call(rbind, rows)
    }


    # =============================================================================
    # SECTION B: SKCS decomposition
    # =============================================================================

    compute_skcs <- function(df, region_label) {
        goi_vec  <- df$GOI
        cdfg_vec <- df$Expression
        n_cells  <- nrow(df)
    
        rk_goi  <- rank(goi_vec,  ties.method = "average")
        rk_cdfg <- rank(cdfg_vec, ties.method = "average")
        mr_goi  <- mean(rk_goi)
        mr_cdfg <- mean(rk_cdfg)
        spear_contrib <- (rk_goi - mr_goi) * (rk_cdfg - mr_cdfg) / (n_cells - 1)
    
        # == Normalised Spearman contribution (Si_norm) =======================
        # Si_norm = Si / Si_max, with Si_max = (n-1)/4, maps Si to about [-1, 1].
        # Si_norm > 0.05 excludes the IP transition zone, where contributions are ambiguous.
        si_max       <- (n_cells - 1) / 4
        spear_norm   <- spear_contrib / si_max     # same sign as spear_contrib
    
        global_rho <- cor(goi_vec, cdfg_vec, method = "spearman")
        cat(sprintf("  [%s] Global Spearman rho = %.4f  |  Si_max = %.2f\n",
                    region_label, global_rho, si_max))
    
        KENDALL_MAX_N <- 2000L
        KENDALL_CHUNK <- 500L
    
        # Seed before the Kendall computation.
        set.seed(42)
        if (n_cells <= KENDALL_MAX_N) {
            sign_goi  <- outer(goi_vec,  goi_vec,  "-")
            sign_cdfg <- outer(cdfg_vec, cdfg_vec, "-")
            sg <- sign(sign_goi)
            sc <- sign(sign_cdfg)
            prod_sc <- sg * sc
            concordant  <- prod_sc > 0
            discordant  <- prod_sc < 0
            diag(concordant) <- FALSE
            diag(discordant) <- FALSE
            # tau_b: per-cell tie-adjusted denominator
            tied_goi  <- sg == 0; diag(tied_goi)  <- FALSE
            tied_cdfg <- sc == 0; diag(tied_cdfg) <- FALSE
            n0_goi  <- rowSums(tied_goi)
            n0_cdfg <- rowSums(tied_cdfg)
            Ci_minus_Di <- rowSums(concordant) - rowSums(discordant)
            denom_b <- sqrt((n_cells - 1 - n0_goi) * (n_cells - 1 - n0_cdfg))
            denom_b[denom_b == 0] <- 1  # guard: avoid division by zero
            kend_contrib <- Ci_minus_Di / denom_b
            rm(sign_goi, sign_cdfg, sg, sc, prod_sc,
               concordant, discordant, tied_goi, tied_cdfg)
            gc(verbose = FALSE)
        } else {
            cat(sprintf(
                "  [%s] n=%d > %d -- using chunked Kendall (chunk=%d)\n",
                region_label, n_cells, KENDALL_MAX_N, KENDALL_CHUNK))
            kend_contrib <- numeric(n_cells)
            n0_goi_vec   <- numeric(n_cells)
            n0_cdfg_vec  <- numeric(n_cells)
            chunk_starts <- seq(1L, n_cells, by = KENDALL_CHUNK)
            for (cs in chunk_starts) {
                ce     <- min(cs + KENDALL_CHUNK - 1L, n_cells)
                idx    <- cs:ce
                sg_blk <- sign(outer(goi_vec[idx],  goi_vec,  "-"))
                sc_blk <- sign(outer(cdfg_vec[idx], cdfg_vec, "-"))
                prod   <- sg_blk * sc_blk
                for (ii in seq_along(idx)) prod[ii, idx[ii]] <- 0L
                kend_contrib[idx] <- rowSums(prod > 0L) - rowSums(prod < 0L)
                # tau_b: count ties per cell (exclude self via zeroed diagonal)
                tied_g <- sg_blk == 0L; for (ii in seq_along(idx)) tied_g[ii, idx[ii]] <- FALSE
                tied_c <- sc_blk == 0L; for (ii in seq_along(idx)) tied_c[ii, idx[ii]] <- FALSE
                n0_goi_vec[idx]  <- rowSums(tied_g)
                n0_cdfg_vec[idx] <- rowSums(tied_c)
                rm(sg_blk, sc_blk, prod, tied_g, tied_c); gc(verbose = FALSE)
            }
            Ci_minus_Di <- kend_contrib   # still raw integer counts at this point
            denom_b <- sqrt((n_cells - 1 - n0_goi_vec) * (n_cells - 1 - n0_cdfg_vec))
            denom_b[denom_b == 0] <- 1  # guard: avoid division by zero
            kend_contrib <- Ci_minus_Di / denom_b
        }
    
        # tau_b sanity check (confirm this code path was executed)
        n_tied_cells <- sum(denom_b != (n_cells - 1))
        cat(sprintf("  [%s] Using tau_b denominator | cells with ties: %d / %d\n",
                    region_label, n_tied_cells, n_cells))
    
        global_tau <- cor(goi_vec, cdfg_vec, method = "kendall")
        cat(sprintf("  [%s] Global Kendall tau  = %.4f\n",
                    region_label, global_tau))
    
        is_skcs <- (spear_norm > SKCS_S_THRESH) &
            (kend_contrib  > SKCS_K_THRESH)
        n_skcs  <- sum(is_skcs)
        n_nskcs <- sum(!is_skcs)
    
        cat(sprintf(
            "  [%s] SKCS-core: %d (%.1f%%)  |  Non-SKCS: %d (%.1f%%)  [S_norm>%.2f, K>%.1f]\n",
            region_label,
            n_skcs,  100 * n_skcs  / n_cells,
            n_nskcs, 100 * n_nskcs / n_cells,
            SKCS_S_THRESH, SKCS_K_THRESH))
    
        ks_stat <- if (n_skcs >= 2 && n_nskcs >= 2) {
            suppressWarnings(
                ks.test(goi_vec[is_skcs], goi_vec[!is_skcs])$statistic)
        } else NA_real_
        if (!is.na(ks_stat) && ks_stat > 0.60)
            cat(sprintf(
                "  [%s] WARNING: KS = %.4f > 0.60 -- SKCS-core skewed toward high GOI.\n",
                region_label, ks_stat))
    
        list(spear_contrib = spear_contrib,   # raw Si
             spear_norm    = spear_norm,      # normalised Si = Si / Si_max
             si_max        = si_max,          # (n-1)/4
             kend_contrib  = kend_contrib,
             is_skcs       = is_skcs,
             n_skcs        = n_skcs,
             n_nskcs       = n_nskcs,
             global_rho    = global_rho,
             global_tau    = global_tau,
             ks_stat       = ks_stat)
    }


    # =============================================================================
    # SECTION C: Heatmap helper
    # =============================================================================

    zscore_cap <- function(x)
        pmin(pmax((x - mean(x, na.rm = TRUE)) /
                      (sd(x, na.rm = TRUE) + 1e-9), -1), 2.0)

    col_fun_mako <- colorRamp2(c(-1, 0, 1.0, 2.0),
                               viridis(4, option = "mako"))

    resolve_gene <- function(g, gene_syms) {
        if (g %in% gene_syms) return(g)
        NA_character_
    }

    # Returns draw functions for the merged PDF and saves a standalone PNG.
    draw_skcs_heatmap_p10 <- function(skcs_bc, expr_mat,
                                      dfg_pres_hm, hkg_pres_hm,
                                      region_label, dataset_tag,
                                      out_dir_heatmap,
                                      return_draw_fn = FALSE) {
        common_bc  <- intersect(skcs_bc, rownames(expr_mat))
        n_cells_hm <- length(common_bc)
    
        if (n_cells_hm < MIN_CELLS) {
            cat("  Heatmap skipped for", region_label,
                ": only", n_cells_hm, "barcodes\n")
            if (return_draw_fn) return(NULL)
            return(invisible(NULL))
        }
    
        cat("  Building SKCS-core heatmap:", region_label,
            "| n =", n_cells_hm, "\n")
    
        sub_expr <- expr_mat[common_bc, , drop = FALSE]
        goi_vals <- sub_expr[, GOI_10]
        cell_ord <- order(goi_vals)
        sub_expr <- sub_expr[cell_ord, ]
    
        cdfg_row <- if (length(dfg_pres_hm) == 0) NULL else
            if (length(dfg_pres_hm) == 1) sub_expr[, dfg_pres_hm] else
                rowMeans(sub_expr[, dfg_pres_hm, drop = FALSE])
        chkg_row <- if (length(hkg_pres_hm) == 0) NULL else
            if (length(hkg_pres_hm) == 1) sub_expr[, hkg_pres_hm] else
                rowMeans(sub_expr[, hkg_pres_hm, drop = FALSE])
    
        row_list <- list(
            matrix(sub_expr[, GOI_10], nrow = 1,
                   dimnames = list(GOI_10, NULL)))
        if (!is.null(cdfg_row))
            row_list[[length(row_list) + 1]] <-
            matrix(cdfg_row, nrow = 1, dimnames = list(DFG_DISPLAY_LABEL, NULL))
        for (g in dfg_pres_hm)
            row_list[[length(row_list) + 1]] <-
            matrix(sub_expr[, g], nrow = 1, dimnames = list(g, NULL))
        if (!is.null(chkg_row)) {
            row_list[[length(row_list) + 1]] <-
                matrix(chkg_row, nrow = 1, dimnames = list("cHKG", NULL))
            for (g in hkg_pres_hm)
                row_list[[length(row_list) + 1]] <-
                    matrix(sub_expr[, g], nrow = 1, dimnames = list(g, NULL))
        }
    
        expr_final <- do.call(rbind, row_list)
        mat_scaled <- t(apply(expr_final, 1, zscore_cap))
        dimnames(mat_scaled) <- dimnames(expr_final)
        rn <- rownames(mat_scaled)
    
        row_cat <- setNames(character(length(rn)), rn)
        row_cat[GOI_10]            <- "GOI"
        if (DFG_DISPLAY_LABEL %in% rn) row_cat[DFG_DISPLAY_LABEL] <- "cDFG_composite"
        row_cat[rn %in% dfg_pres_hm] <- "DFG"
        if ("cHKG" %in% rn) row_cat["cHKG"] <- "cHKG"
        row_cat[rn %in% hkg_pres_hm] <- "HKG"
    
        cat_cols <- c(GOI = "#CC2200", cDFG_composite = "#220044",
                      DFG = "#4B0082", cHKG = "#003355", HKG = "#006699")
    
        row_ann <- rowAnnotation(
            Category = row_cat,
            col      = list(Category = cat_cols),
            annotation_legend_param = list(
                Category = list(title     = "Gene category",
                                labels_gp = gpar(fontsize = 8),
                                title_gp  = gpar(fontsize = 9,
                                                 fontface = "bold"))),
            show_annotation_name = FALSE,
            width = unit(3, "mm"))
    
        row_name_gp <- gpar(
            fontsize = 8,
            fontface = ifelse(rn %in% c(DFG_DISPLAY_LABEL, "cHKG"), "bold", "plain"),
            col = ifelse(rn == GOI_10,          "#CC2200",
                         ifelse(rn == DFG_DISPLAY_LABEL,  "#220044",
                                ifelse(rn %in% dfg_pres_hm,     "#4B0082",
                                       ifelse(rn == "cHKG",            "#003355",
                                              ifelse(rn %in% hkg_pres_hm,     "#006699", "black"))))))
    
        block1    <- c(GOI_10, if (DFG_DISPLAY_LABEL %in% rn) DFG_DISPLAY_LABEL, dfg_pres_hm)
        block_label <- if (GOI %in% HKG_GENES) "GOI / cHKG" else "GOI / DFG"
        # The row split is built only over blocks that have rows.
        row_block   <- ifelse(rn %in% block1, block_label, "HKG")
        .blk_order  <- c(block_label, "HKG")
        row_split   <- factor(row_block,
                              levels = .blk_order[.blk_order %in% row_block])
    
        ht <- Heatmap(
            mat_scaled,
            name                 = "z-score",
            col                  = col_fun_mako,
            cluster_rows         = FALSE,
            show_row_names       = TRUE,
            row_names_side       = "right",
            row_names_gp         = row_name_gp,
            row_split            = row_split,
            row_title_gp         = gpar(fontsize = 9, fontface = "bold"),
            row_gap              = unit(1, "mm"),
            cluster_columns      = FALSE,
            show_column_names    = FALSE,
            column_title         = NULL,
            border               = TRUE,
            use_raster           = TRUE,
            raster_quality       = 3,
            right_annotation     = row_ann,
            heatmap_legend_param = list(
                title         = "z-score",
                title_gp      = gpar(fontsize = 9, fontface = "bold"),
                labels_gp     = gpar(fontsize = 8),
                legend_height = unit(3.5, "cm"),
                at            = c(-1, 0, 1, 2),
                labels        = c("-1", "0", "1", "2 (cap)")))
    
        cell_label  <- if (grepl("AllCells", region_label)) "All cells"
        else "SKCS-core"
        title_str   <- paste0(GOI_10, " vs DFG & HKG -- ", region_label,
                              " ", cell_label, " (n=", n_cells_hm,
                              " cells) [Module1; z-score - 1 to +2.0]")
        pdf_w <- min(50, max(14, round(n_cells_hm * 0.008 + 6)))
    
        draw_fn <- function() {
            draw(ht,
                 padding                = unit(c(10, 4, 12, 4), "mm"),
                 annotation_legend_side = "right",
                 heatmap_legend_side    = "right",
                 merge_legends          = TRUE)
            grid.text(title_str,
                      x = unit(0.01, "npc"), y = unit(0.99, "npc"),
                      just = c("left", "top"),
                      gp   = gpar(fontsize = 10, fontface = "bold",
                                  col = "#1A1A1A"))
        }
    
        # Standalone PNG
        hm_png <- file.path(out_dir_heatmap,
                            paste0(FILE_PREFIX, "_", dataset_tag, "_",
                                   region_label, "_SKCScore_Heatmap.png"))
        # The PNG keeps the PDF aspect; the resolution is reduced instead of the width.
        hm_png_res <- min(150, floor(4000 / pdf_w))
        tryCatch({
            png(hm_png, width = pdf_w, height = 700 / 150, units = "in",
                res = hm_png_res)
            draw_fn(); dev.off()
        }, error = function(e) {
            cat("    Heatmap PNG error:", conditionMessage(e), "\n")
            try(dev.off(), silent = TRUE)
        })
        safe_dev_set_p10()
    
        if (return_draw_fn)
            return(list(draw_fn = draw_fn, pdf_w = pdf_w))
        invisible(NULL)
    }


    # =============================================================================
    # SECTION D: Master pipeline
    # =============================================================================

    run_p10_pipeline <- function(dataset_tag,
                                 dataset_dir,
                                 regions,
                                 p9_rds_path,
                                 dfg_genes = DFG_10,
                                 hkg_genes = HKG_10,
                                 group_label_prefix = "SEAAD") {
    
        group_label        <- paste0(group_label_prefix, "_", dataset_tag,
                                     "_Module1_SKCS")
        pooled_key         <- "Pooled"   # internal key used in lists/CSVs
    
        set.seed(42)

        # Dataset output folders; local bindings, so the arm-level paths are unchanged.
        OUT_DIR     <- file.path(ARM_OUT_DIR,
                                 paste0(basename(ARM_OUT_DIR), "_", dataset_tag))
        OUT_DIR_SUB <- file.path(OUT_DIR, "Part10_subregions")
        dir.create(OUT_DIR_SUB, showWarnings = FALSE, recursive = TRUE)
        DS_DIRS <<- unique(c(DS_DIRS, OUT_DIR))

        cat("\n", strrep("=", 72), "\n")
        cat("  Part 10 Module 1 pipeline:", dataset_tag, "\n")
        cat("  Output dir:", OUT_DIR, "\n")
        cat("  Regions:", paste(regions, collapse = ", "), "\n")
        cat("  GOI:", GOI_10, "  DFG:",
            paste(dfg_genes, collapse = ", "), "\n")
        cat("  HKG:", paste(hkg_genes, collapse = ", "), "\n")
        cat(strrep("=", 72), "\n\n")
    
        # == D1. Load Part 9 baseline DE ======================================
        p9_baseline      <- setNames(rep(NA_real_, length(regions)), regions)
        p9_baseline_pool <- NA_real_
        p9_pooled_model  <- NULL   # full-population pooled GAM (for CEP-IP scatter)
        p9_pooled_data   <- NULL   # data the Part 9 pooled GAM was fit on
    
        if (file.exists(p9_rds_path)) {
            p9_obj         <- readRDS(p9_rds_path)
            p9_results     <- p9_obj$all_results
            p9_pooled_name <- p9_obj$pooled_region_name %||%
                paste0(dataset_tag, "_Pooled")

            # Report the estimator behind the Part 9 fits.
            .p9_re_flag <- p9_obj$donor_re
            .p9_re_seen <- vapply(p9_results, function(r) {
                if (is.null(r) || is.null(r$best_model)) return(NA)
                model_has_re(r$best_model)
            }, logical(1))
            .p9_re_n <- sum(.p9_re_seen, na.rm = TRUE)
            if (is.null(.p9_re_flag)) {
                cat("  Part 9 RDS carries no donor_re field.",
                    "Models with a donor random intercept:",
                    .p9_re_n, "/", sum(!is.na(.p9_re_seen)), "\n")
                if (.p9_re_n == 0L)
                    cat("  WARNING: the Part 9 fits are cell level; the baseline DE",
                        "is not on the same scale as the SKCS-core DE.\n")
            } else {
                cat("  Part 9 estimator: donor random intercept =",
                    if (isTRUE(.p9_re_flag)) "TRUE" else "FALSE",
                    "| models carrying it:", .p9_re_n, "/",
                    sum(!is.na(.p9_re_seen)), "\n")
            }

            # Baseline DE recomputed against the donor-aware null, as for the fits here.
            for (reg in regions) {
                res <- p9_results[[paste0(reg, "|DFG")]]
                if (!is.null(res) && !is.null(res$best_model))
                    p9_baseline[reg] <- de_explained(res$best_model, res$data)
            }
            res_pool <- p9_results[[paste0(p9_pooled_name, "|DFG")]]
            if (!is.null(res_pool) && !is.null(res_pool$best_model)) {
                p9_baseline_pool <- de_explained(res_pool$best_model,
                                                 res_pool$data)
                p9_pooled_model  <- res_pool$best_model
                p9_pooled_data   <- res_pool$data
            }
        
            cat("  Part 9 baseline DE (full population, cDFG composite,",
                "donor-aware null):\n")
            for (reg in regions)
                cat(sprintf("    %-12s : %s\n", reg,
                            ifelse(is.na(p9_baseline[reg]), "NA",
                                   sprintf("%.1f%%", p9_baseline[reg] * 100))))
            cat(sprintf("    %-12s : %s\n", "Pooled",
                        ifelse(is.na(p9_baseline_pool), "NA",
                               sprintf("%.1f%%", p9_baseline_pool * 100))))
            cat("\n")
            rm(p9_obj); gc(verbose = FALSE)
        } else {
            cat("  WARNING: Part 9 RDS not found at", p9_rds_path,
                "-- DE comparison skipped.\n\n")
        }
    
        # == D2. Load norm.mtx ================================================
        mtx_path  <- file.path(dataset_dir, "deg_support", "norm.mtx")
        bar_path  <- file.path(dataset_dir, "deg_support", "norm_barcodes.tsv")
        feat_path <- file.path(dataset_dir, "deg_support", "norm_features.tsv")
        X_norm <- NULL; norm_bc <- NULL; norm_genes <- NULL
        has_norm_mtx <- all(file.exists(c(mtx_path, bar_path, feat_path)))
    
        if (has_norm_mtx) {
            cat("  Loading normalised MTX for heatmaps + HKG ...\n")
            X_norm     <- readMM(mtx_path)
            norm_bc    <- readLines(bar_path, warn = FALSE)
            # Quoting and comment characters disabled for gene symbols.
            feat_df    <- read.table(feat_path, sep = "\t", header = FALSE,
                                     quote = "", comment.char = "",
                                     stringsAsFactors = FALSE)
            norm_genes <- feat_df[, 2]
            rownames(X_norm) <- norm_bc
            colnames(X_norm) <- norm_genes
            cat(sprintf("    Matrix: %d cells x %d genes\n",
                        nrow(X_norm), ncol(X_norm)))
        } else {
            cat("  WARNING: norm.mtx not found -- heatmaps and HKG skipped\n\n")
        }
    
        dfg_pres_hm <- character(0); hkg_pres_hm <- character(0)
        expr_mat_full <- NULL
    
        if (has_norm_mtx) {
            hm_genes_needed <- c(GOI_10, dfg_genes, hkg_genes)
            resolved_genes  <- setNames(
                sapply(hm_genes_needed, resolve_gene, gene_syms = norm_genes),
                hm_genes_needed)
            valid_genes   <- resolved_genes[!is.na(resolved_genes)]
            expr_mat_full <- as.matrix(X_norm[, valid_genes, drop = FALSE])
            colnames(expr_mat_full) <- names(valid_genes)
            # The heatmap indexes this matrix by the GOI, so a missing GOI skips heatmaps.
            if (!(GOI_10 %in% colnames(expr_mat_full))) {
                cat("  WARNING: GOI", GOI_10,
                    "is absent from norm_features.tsv --",
                    "heatmaps skipped for this dataset\n")
                has_norm_mtx  <- FALSE
                expr_mat_full <- NULL
            }
            dfg_pres_hm <- dfg_genes[dfg_genes %in% colnames(expr_mat_full)]
            hkg_pres_hm <- hkg_genes[hkg_genes %in% colnames(expr_mat_full)]
            cat("  DFG resolved (heatmap):", paste(dfg_pres_hm, collapse=", "),"\n")
            cat("  HKG resolved (heatmap):", paste(hkg_pres_hm, collapse=", "),"\n\n")
        }
    
        # == D3. Per-region data ==============================================
        cat("  Building per-region data (GOI > 0 filter) ...\n")
        region_dfs      <- list()
        region_dfs_all  <- list()   # full CSV (GOI-pos + GOI-neg), for CSV export
        cell_audit_list <- list()
    
        for (reg in regions) {
            gam_path <- p10_gam_file(dataset_dir, reg)
            if (!file.exists(gam_path)) {
                cat("  Region", reg, ": file not found --",
                    basename(gam_path), "\n")
                cell_audit_list[[reg]] <- list(
                    is_pooled = FALSE, n_raw = NA_integer_,
                    n_goi_neg = NA_integer_, n_goi_pos = NA_integer_,
                    n_used    = 0L, n_donors = NA_integer_,
                    sum_w     = NA_real_,
                    note      = "File not found -- region skipped")
                next
            }
            df_all    <- read.csv(gam_path, stringsAsFactors = FALSE,
                                  check.names = FALSE)
            df_all$region <- reg

            # Donor is carried through every frame that reaches a gam() call.
            if (!("donor" %in% colnames(df_all)))
                stop(sprintf(
                    "Region %s: gam_data has no donor column, so the donor random intercept cannot be fitted.",
                    reg))
            if (anyNA(df_all$donor))
                stop(sprintf("Region %s: %d cell(s) have a missing donor label.",
                             reg, sum(is.na(df_all$donor))))
            df_all$donor <- as_donor_factor(df_all$donor)
        
        
            # Composite recomputed from the individual DFG columns.
            dfg_in_all <- dfg_genes[dfg_genes %in% colnames(df_all)]
            if (length(dfg_in_all) < length(dfg_genes))
                stop(sprintf("Region %s: gam_data missing DFG column(s): %s",
                             reg, paste(setdiff(dfg_genes, dfg_in_all),
                                        collapse = ", ")))
            df_all$Expression <- rowMeans(
                as.matrix(df_all[, dfg_in_all, drop = FALSE]))
        
            n_raw     <- nrow(df_all)
            n_goi_neg <- sum(df_all$GOI <= 0)   # counts zero AND negative
            n_goi_pos <- n_raw - n_goi_neg
            df_pos    <- df_all[df_all$GOI > 0, ]
        
            n_don_reg <- n_donor_levels(df_pos$donor)
            # The region's own weights (for the audit); fitted subsets recompute theirs.
            df_pos$w  <- donor_weights(df_pos$donor)

            cat(sprintf(
                "  Region %s: %d raw | %d GOI=0 excluded | %d GOI>0 retained | %d donors | sum(w)=%.4f\n",
                reg, n_raw, n_goi_neg, n_goi_pos, n_don_reg, sum(df_pos$w)))
        
            if (n_goi_pos < MIN_CELLS) {
                cat("  Region", reg, ": n_goi_pos =", n_goi_pos,
                    "< MIN_CELLS -- skipped\n")
                cell_audit_list[[reg]] <- list(
                    is_pooled = FALSE, n_raw = n_raw,
                    n_goi_neg = n_goi_neg, n_goi_pos = n_goi_pos, n_used = 0L,
                    n_donors  = n_don_reg, sum_w = sum(df_pos$w),
                    note      = paste0("Skipped: n_goi_pos (", n_goi_pos,
                                       ") < MIN_CELLS (", MIN_CELLS, ")"))
                next
            }
            region_dfs[[reg]]     <- df_pos
            region_dfs_all[[reg]] <- df_all   # all cells (GOI-pos + GOI-neg)
            cell_audit_list[[reg]] <- list(
                is_pooled = FALSE, n_raw = n_raw,
                n_goi_neg = n_goi_neg, n_goi_pos = n_goi_pos,
                n_used    = n_goi_pos, n_donors = n_don_reg,
                sum_w     = sum(df_pos$w), note = "OK")
        }
    
        if (length(region_dfs) == 0) {
            cat("  No regions have sufficient data. Aborting for", dataset_tag,"\n")
            return(invisible(NULL))
        }
    
        # == D3b. Early cell count summary ====================================
        cat("\n  Cell count check (per-region):\n")
        for (nm in names(cell_audit_list)) {
            ca <- cell_audit_list[[nm]]
            ok <- if (ca$n_used == 0L && !isTRUE(ca$n_goi_pos == 0L)) NA
            else isTRUE(ca$n_goi_pos == ca$n_used)
            cat(sprintf(
                "    %-12s  n_raw=%-5s  n_goi_neg=%-5s  n_goi_pos=%-5s  n_used=%-5s  %s\n",
                nm,
                ifelse(is.na(ca$n_raw),     "NA", as.character(ca$n_raw)),
                ifelse(is.na(ca$n_goi_neg), "NA", as.character(ca$n_goi_neg)),
                ifelse(is.na(ca$n_goi_pos), "NA", as.character(ca$n_goi_pos)),
                as.character(ca$n_used),
                if (is.na(ca$n_raw)) "MISSING"
                else if (is.na(ok))  "SKIPPED"
                else if (ok)         "OK"
                else                 "MISMATCH"))
        }
        cat("\n")
    
        # == D4. SKCS decomposition, per region ===============================
        cat(strrep("-", 60), "\n")
        cat("  SECTION D4: Per-region SKCS decomposition\n")
        cat(strrep("-", 60), "\n\n")
    
        region_skcs      <- list()
        region_score_dfs <- list()
    
        for (reg in names(region_dfs)) {
            df <- region_dfs[[reg]]
            cat(strrep("-", 50), "\n")
            cat("  Region:", reg, "\n")
            set.seed(42 + which(names(region_dfs) == reg))
            skcs <- compute_skcs(df, reg)
            region_skcs[[reg]] <- skcs
        
            score_df <- data.frame(
                barcode       = df$barcode,
                region        = reg,
                dataset       = dataset_tag,
                GOI_expr      = df$GOI,
                cDFG_expr     = df$Expression,
                spear_contrib = skcs$spear_contrib,
                spear_norm    = skcs$spear_norm,
                kend_contrib  = skcs$kend_contrib,
                is_GOI_neg    = FALSE,
                is_skcs       = skcs$is_skcs,
                is_trep       = FALSE,
                stringsAsFactors = FALSE)
            region_score_dfs[[reg]] <- score_df
        
            # SKCS scatter -> subregion PNG + PDF (goes into merged PDF later)
            p_scatter <- make_skcs_scatter_p10(
                spear_norm    = skcs$spear_norm,
                kend_contrib  = skcs$kend_contrib,
                is_skcs       = skcs$is_skcs,
                label_title   = paste0(reg, " (", dataset_tag, ")"),
                global_rho    = skcs$global_rho,
                global_tau    = skcs$global_tau,
                n_skcs        = skcs$n_skcs,
                n_nskcs       = skcs$n_nskcs)
            save_plot_pdf_png(
                p_scatter,
                file.path(OUT_DIR_SUB, paste0(FILE_PREFIX, "_", dataset_tag,
                                              "_", reg, "_SKCS_scatter")),
                pdf_w = 7, pdf_h = 7, png_w = 1400, png_h = 1400, png_res = 200)
        
            # GOI distribution diagnostic -> subregion PNG + PDF
            diag_df <- data.frame(
                GOI_expr = df$GOI,
                Group    = ifelse(skcs$is_skcs, "SKCS-core", "Non-SKCS"),
                stringsAsFactors = FALSE)
            p_density <- ggplot(diag_df,
                                aes(x = GOI_expr, fill = Group,
                                    colour = Group)) +
                geom_density(alpha = 0.40, linewidth = 0.8) +
                scale_fill_manual(
                    values = c("SKCS-core" = "#4B0082", "Non-SKCS" = "#4F75DE")) +
                scale_colour_manual(
                    values = c("SKCS-core" = "#4B0082", "Non-SKCS" = "#4F75DE")) +
                labs(title    = paste0(reg, " -- ", GOI_10,
                                       " density: SKCS-core vs Non-SKCS"),
                     subtitle = "KS > 0.60 would indicate circularity",
                     x = paste0(GOI_10, " (", NORM_TAG, " log1p)"), y = "Density",
                     fill = "Group", colour = "Group") +
                theme_p10(base_size = 11) +
                theme(legend.position = "right")
        
            diag_base <- file.path(OUT_DIR_SUB,
                                   paste0(FILE_PREFIX, "_", dataset_tag,
                                          "_", reg, "_GOI_dist_diagnostic"))
            tryCatch({
                pdf(paste0(diag_base, ".pdf"), width = 7, height = 5)
                print(p_density)
                dev.off()
            }, error = function(e) {
                cat("    Diag PDF error:", conditionMessage(e), "\n")
                try(dev.off(), silent = TRUE)
            })
            tryCatch({
                png(paste0(diag_base, ".png"), width = 1400, height = 1000,
                    res = 200)
                print(p_density)
                dev.off()
            }, error = function(e) {
                cat("    Diag PNG error:", conditionMessage(e), "\n")
                try(dev.off(), silent = TRUE)
            })
            safe_dev_set_p10()
        }
    
        # == D5. SKCS decomposition, pooled ===================================
        cat(strrep("-", 60), "\n")
        cat("  SECTION D5: Pooled SKCS decomposition\n")
        cat(strrep("-", 60), "\n\n")
    
        df_pooled <- do.call(rbind, region_dfs)
        rownames(df_pooled) <- NULL
        # Pooled donor factor rebuilt: a donor in several regions becomes one level.
        df_pooled$donor <- as_donor_factor(df_pooled$donor)
        # Weights recomputed on the pooled frame.
        df_pooled$w     <- donor_weights(df_pooled$donor)
        cat("  Pooled dataset: n =", nrow(df_pooled), "cells from",
            length(region_dfs), "regions\n")
    
        set.seed(142)
        skcs_pooled <- compute_skcs(df_pooled, paste0(dataset_tag, "_Pooled"))
    
        {
            n_raw_pool <- sum(sapply(cell_audit_list[names(region_dfs)],
                                     function(ca) ca$n_raw), na.rm = TRUE)
            n_neg_pool <- sum(sapply(cell_audit_list[names(region_dfs)],
                                     function(ca) ca$n_goi_neg), na.rm = TRUE)
            n_pos_pool <- nrow(df_pooled)
            ok_pool    <- isTRUE(n_pos_pool == (n_raw_pool - n_neg_pool))
            cell_audit_list[[pooled_key]] <- list(
                is_pooled = TRUE,
                n_raw     = n_raw_pool,
                n_goi_neg = n_neg_pool,
                n_goi_pos = n_pos_pool,
                n_used    = n_pos_pool,
                n_donors  = n_donor_levels(df_pooled$donor),
                sum_w     = sum(df_pooled$w),
                note      = if (ok_pool) "OK (sum of constituent regions)"
                else paste0("MISMATCH: expected ",
                            n_raw_pool - n_neg_pool,
                            " but n_used = ", n_pos_pool))
            cat(sprintf(
                "  Pooled audit: n_raw=%d | n_goi_neg=%d | n_goi_pos=%d | %s\n",
                n_raw_pool, n_neg_pool, n_pos_pool,
                if (ok_pool) "OK" else "MISMATCH"))
        }
    
        # Pooled score data frame (is_trep filled in D7b)
        score_pooled_df <- data.frame(
            barcode       = df_pooled$barcode,
            region        = df_pooled$region,
            dataset       = dataset_tag,
            GOI_expr      = df_pooled$GOI,
            cDFG_expr     = df_pooled$Expression,
            spear_contrib = skcs_pooled$spear_contrib,
            spear_norm    = skcs_pooled$spear_norm,
            kend_contrib  = skcs_pooled$kend_contrib,
            is_GOI_neg    = FALSE,
            is_skcs       = skcs_pooled$is_skcs,
            is_trep       = FALSE,
            stringsAsFactors = FALSE)

        # Append individual DFG gene expression columns (for Part11a permutation)
        dfg_in_pooled <- dfg_genes[dfg_genes %in% colnames(df_pooled)]
        for (g in dfg_in_pooled)
            score_pooled_df[[g]] <- df_pooled[[g]]
        if (length(dfg_in_pooled) > 0)
            cat(sprintf("  Added %d/%d DFG expression columns to pooled scores\n",
                        length(dfg_in_pooled), length(dfg_genes)))

        # Append individual HKG gene expression columns (for Part11a HKG control)
        hkg_in_pooled <- hkg_genes[hkg_genes %in% colnames(df_pooled)]
        for (g in hkg_in_pooled)
            score_pooled_df[[g]] <- df_pooled[[g]]
        if (length(hkg_in_pooled) > 0)
            cat(sprintf("  Added %d/%d HKG expression columns to pooled scores\n",
                        length(hkg_in_pooled), length(hkg_genes)))

        # Append biological validation gene expression columns (GLYC + DAM + INFLAM)
        # for Part11d DAM biological validation
        bio_genes_var <- if (exists("BIO_GENES")) BIO_GENES
                         else if (exists("DAM_GENES")) DAM_GENES
                         else character(0)
        if (length(bio_genes_var) > 0) {
            bio_in_pooled <- bio_genes_var[bio_genes_var %in% colnames(df_pooled)]
            for (g in bio_in_pooled)
                score_pooled_df[[g]] <- df_pooled[[g]]
            if (length(bio_in_pooled) > 0)
                cat(sprintf("  Added %d/%d BIO expression columns to pooled scores\n",
                            length(bio_in_pooled), length(bio_genes_var)))
        } else {
            bio_in_pooled <- character(0)
        }
        score_pooled_csv <- file.path(OUT_DIR,
                                      paste0(FILE_PREFIX, "_", dataset_tag,
                                             "_Pooled_SKCS_scores.csv"))
    
        # Pooled SKCS scatter -> pooled PNG + PDF
        p_scatter_pool <- make_skcs_scatter_p10(
            spear_norm    = skcs_pooled$spear_norm,
            kend_contrib  = skcs_pooled$kend_contrib,
            is_skcs       = skcs_pooled$is_skcs,
            label_title   = paste0(dataset_tag, " Pooled"),
            global_rho    = skcs_pooled$global_rho,
            global_tau    = skcs_pooled$global_tau,
            n_skcs        = skcs_pooled$n_skcs,
            n_nskcs       = skcs_pooled$n_nskcs)
        save_plot_pdf_png(
            p_scatter_pool,
            file.path(OUT_DIR, paste0(FILE_PREFIX, "_", dataset_tag,
                                      "_Pooled_SKCS_scatter")),
            pdf_w = 7, pdf_h = 7, png_w = 1400, png_h = 1400, png_res = 200)
        cat("  Pooled SKCS scatter PDF+PNG saved\n")
        emit_plot_pane(p_scatter_pool)

        # == Pooled GOI distribution diagnostic ===============================
        # Density comparison faceted by region, with free y scales.
        diag_pool <- data.frame(
            GOI_expr = df_pooled$GOI,
            Group    = ifelse(skcs_pooled$is_skcs, "SKCS-core", "Non-SKCS"),
            Region   = df_pooled$region,
            stringsAsFactors = FALSE)
        p_dens_pool <- ggplot(diag_pool,
                              aes(x = GOI_expr, fill = Group, colour = Group)) +
            geom_density(alpha = 0.40, linewidth = 0.8) +
            scale_fill_manual(
                values = c("SKCS-core" = "#4B0082", "Non-SKCS" = "#4F75DE")) +
            scale_colour_manual(
                values = c("SKCS-core" = "#4B0082", "Non-SKCS" = "#4F75DE")) +
            facet_wrap(~ Region, nrow = 1, scales = "free_y") +
            labs(title    = paste0(dataset_tag, " Pooled -- ", GOI_10,
                                   " density: SKCS-core vs Non-SKCS"),
                 subtitle = sprintf(
                     "Pooled n=%d | SKCS-core n=%d | rho=%.3f | tau=%.3f | KS=%s  (KS > 0.60 would indicate circularity)",
                     nrow(df_pooled), skcs_pooled$n_skcs,
                     skcs_pooled$global_rho, skcs_pooled$global_tau,
                     ifelse(is.na(skcs_pooled$ks_stat), "NA",
                            sprintf("%.3f", skcs_pooled$ks_stat))),
                 x = paste0(GOI_10, " (", NORM_TAG, " log1p)"), y = "Density",
                 fill = "Group", colour = "Group") +
            theme_p10(base_size = 10) +
            theme(legend.position = "right")
        save_plot_pdf_png(
            p_dens_pool,
            file.path(OUT_DIR, paste0(FILE_PREFIX, "_", dataset_tag,
                                      "_Pooled_GOI_dist_diagnostic")),
            pdf_w = 14, pdf_h = 5, png_w = 2800, png_h = 1000, png_res = 200)
        cat("  Pooled GOI diagnostic PDF+PNG saved\n")
        emit_plot_pane(p_dens_pool)

        # == D5b. Alternate-GOI SKCS scatters =================================
        # HKG8 arm only: each ALT_GOI_GENES gene as GOI against the mean of the other
        # eight, on the pooled frame. Figures only.
        alt_scatter_plots <- list()
        alt_goi_req <- if (identical(ARM_TAG, ALT_GOI_ARM)) ALT_GOI_GENES
                       else character(0)
        alt_goi_use <- alt_goi_req[alt_goi_req %in% dfg_genes &
                                       alt_goi_req %in% colnames(df_pooled)]
        if (length(alt_goi_req) > length(alt_goi_use))
            cat("  Alternate GOI skipped (not a DFG column):",
                paste(setdiff(alt_goi_req, alt_goi_use), collapse = ", "), "\n")

        for (g_alt in alt_goi_use) {
            comp_cols <- setdiff(dfg_genes, g_alt)
            goi_alt   <- df_pooled[[g_alt]]
            comp_alt  <- rowMeans(cbind(
                df_pooled$GOI,
                as.matrix(df_pooled[, comp_cols, drop = FALSE])))
            keep_alt  <- is.finite(goi_alt) & goi_alt > 0 & is.finite(comp_alt)
            n_drop    <- sum(!keep_alt)
            df_alt    <- data.frame(GOI        = goi_alt[keep_alt],
                                    Expression = comp_alt[keep_alt])
            cat(sprintf(
                "  Alternate GOI %s: %d cells (%d dropped: %s = 0 or non-finite)\n",
                g_alt, nrow(df_alt), n_drop, g_alt))
            if (nrow(df_alt) < MIN_CELLS) {
                cat("  Alternate GOI", g_alt, ": below MIN_CELLS -- skipped\n")
                next
            }
            set.seed(142)
            skcs_alt <- compute_skcs(df_alt,
                                     paste0(dataset_tag, "_Pooled_", g_alt))
            p_alt <- make_skcs_scatter_p10(
                spear_norm   = skcs_alt$spear_norm,
                kend_contrib = skcs_alt$kend_contrib,
                is_skcs      = skcs_alt$is_skcs,
                label_title  = paste0(dataset_tag, " Pooled"),
                global_rho   = skcs_alt$global_rho,
                global_tau   = skcs_alt$global_tau,
                n_skcs       = skcs_alt$n_skcs,
                n_nskcs      = skcs_alt$n_nskcs,
                panel_lbl    = paste0("GOI ", g_alt, " / c",
                                      length(comp_cols) + 1L, "HKG"))
            save_plot_pdf_png(
                p_alt,
                file.path(OUT_DIR, paste0(FILE_PREFIX, "_", dataset_tag,
                                          "_Pooled_SKCS_scatter_GOI_", g_alt)),
                pdf_w = 7, pdf_h = 7, png_w = 1400, png_h = 1400, png_res = 200)
            cat("  Alternate-GOI SKCS scatter PDF+PNG saved:", g_alt, "\n")
            alt_scatter_plots[[g_alt]] <- p_alt
            emit_plot_pane(p_alt)
        }
        cat("\n")
    
        # == D6. Heatmaps =====================================================
        cat(strrep("-", 60), "\n")
        cat("  SECTION D6: Heatmaps\n")
        cat(strrep("-", 60), "\n\n")
    
        # Collect heatmap draw functions for inclusion in merged PDF
        heatmap_draw_fns <- list()
    
        if (!has_norm_mtx) {
            cat("  Heatmaps skipped (norm.mtx not available)\n\n")
        } else {
            # Per-region SKCS-core heatmaps -> OUT_DIR_SUB PNG only (PDF in merged)
            for (reg in names(region_skcs)) {
                skcs_r  <- region_skcs[[reg]]
                skcs_bc <- region_dfs[[reg]]$barcode[skcs_r$is_skcs]
                if (skcs_r$n_skcs >= MIN_CELLS) {
                    fn <- draw_skcs_heatmap_p10(
                        skcs_bc, expr_mat_full, dfg_pres_hm, hkg_pres_hm,
                        reg, dataset_tag,
                        out_dir_heatmap = OUT_DIR_SUB,
                        return_draw_fn  = TRUE)
                    if (!is.null(fn))
                        heatmap_draw_fns[[paste0("HM_", reg)]] <- fn
                }
            }
            # SKCS-core pooled heatmap -> OUT_DIR PNG only
            skcs_bc_pool <- df_pooled$barcode[skcs_pooled$is_skcs]
            if (skcs_pooled$n_skcs >= MIN_CELLS) {
                fn <- draw_skcs_heatmap_p10(
                    skcs_bc_pool, expr_mat_full, dfg_pres_hm, hkg_pres_hm,
                    "Pooled", dataset_tag,
                    out_dir_heatmap = OUT_DIR,
                    return_draw_fn  = TRUE)
                if (!is.null(fn))
                    heatmap_draw_fns[["HM_Pooled"]] <- fn
            }
            # All-cells pooled heatmap -> OUT_DIR PNG only
            all_bc_pool <- df_pooled$barcode
            if (length(all_bc_pool) >= MIN_CELLS) {
                fn <- draw_skcs_heatmap_p10(
                    all_bc_pool, expr_mat_full, dfg_pres_hm, hkg_pres_hm,
                    "AllCells_Pooled", dataset_tag,
                    out_dir_heatmap = OUT_DIR,
                    return_draw_fn  = TRUE)
                if (!is.null(fn))
                    heatmap_draw_fns[["HM_AllCells_Pooled"]] <- fn
            }
        }
    
        # == D7. GAM-REML-PRSS on SKCS-core cells =============================
        cat(strrep("-", 60), "\n")
        cat("  SECTION D7: GAM-REML-PRSS on SKCS-core cells\n")
        cat(strrep("-", 60), "\n\n")
    
        all_skcs_results   <- list()
        de_comparison_rows <- list()
        trep_map           <- list()
    
        run_gam_for_subset <- function(df_sub, region_id, skcs_info, p9_de) {
            n_sub <- nrow(df_sub)
            if (n_sub < MIN_CELLS) {
                cat("  ", region_id, ": n =", n_sub,
                    "< MIN_CELLS -- GAM skipped\n\n")
                return(list(results = list(), de_row = NULL,
                            trep = list(trep_barcodes = character(0),
                                        ep_df         = data.frame(),
                                        dev_explained = NA_real_,
                                        dev_explained_legacy = NA_real_)))
            }
            cat(strrep("-", 45), "\n")
            cat("  Region/set:", region_id, "| SKCS-core n =", n_sub, "\n")
            results_out <- list()
        
            # Donor reaches every gam() call in this block.
            donor_sub <- as_donor_factor(df_sub$donor)
            n_don_sub <- nlevels(donor_sub)
            # Weights recomputed on the subset.
            w_sub     <- donor_weights(donor_sub)
            dfg_data  <- data.frame(GOI        = df_sub$GOI,
                                    Expression = df_sub$Expression,
                                    donor      = donor_sub,
                                    w          = w_sub)
            cat("  Fitting SKCS-core DFG (cDFG) GAM ...\n")
            results_out[[paste0(region_id, "|DFG|SKCS")]] <- tryCatch(
                fit_gam_prss_p10(dfg_data, region_id, "DFG_SKCS"),
                error = function(e) {
                    cat("  FAILED:", conditionMessage(e), "\n"); NULL })
        
            # HKG GAM
            chkg_sub <- NULL
            if (has_norm_mtx && length(hkg_pres_hm) > 0) {
                bc_sub    <- df_sub$barcode
                idx       <- match(bc_sub, norm_bc)
                n_matched <- sum(!is.na(idx))
                if (n_matched >= MIN_CELLS) {
                    hkg_sub_mat <- matrix(0, nrow = n_sub,
                                          ncol = length(hkg_pres_hm),
                                          dimnames = list(bc_sub, hkg_pres_hm))
                    valid <- !is.na(idx)
                    hkg_sub_mat[valid, ] <- as.matrix(
                        expr_mat_full[idx[valid], hkg_pres_hm, drop = FALSE])
                    chkg_sub <- if (length(hkg_pres_hm) == 1)
                        as.numeric(hkg_sub_mat[, 1])
                    else rowMeans(hkg_sub_mat)
                    cat("  HKG joined:", n_matched, "/", n_sub,
                        "barcodes matched\n")
                } else {
                    cat("  HKG join: only", n_matched,
                        "matches -- HKG GAM skipped\n")
                }
            } else if (!has_norm_mtx) {
                cat("  HKG GAM skipped (norm.mtx unavailable)\n")
            }
            if (!is.null(chkg_sub)) {
                hkg_data <- data.frame(GOI        = df_sub$GOI,
                                       Expression = chkg_sub,
                                       donor      = donor_sub,
                                       w          = w_sub)
                cat("  Fitting SKCS-core HKG GAM ...\n")
                results_out[[paste0(region_id, "|HKG|SKCS")]] <- tryCatch(
                    fit_gam_prss_p10(hkg_data, region_id, "HKG_SKCS"),
                    error = function(e) {
                        cat("  FAILED:", conditionMessage(e), "\n"); NULL })
            }
        
            res_dfg <- results_out[[paste0(region_id, "|DFG|SKCS")]]
            # DE against the donor-aware null, as for the Part 9 baseline.
            de_skcs <- if (!is.null(res_dfg) && !is.null(res_dfg$best_model))
                de_explained(res_dfg$best_model, res_dfg$data) else NA_real_
            de_skcs_flat <- if (!is.null(res_dfg) && !is.null(res_dfg$best_model))
                summary(res_dfg$best_model)$dev.expl else NA_real_
            de_full <- if (!is.null(p9_de)) p9_de else NA_real_

            de_row <- data.frame(
                region       = region_id,
                dataset      = dataset_tag,
                n_full       = NA_integer_,
                n_skcs       = n_sub,
                n_nonskcs    = skcs_info$n_nskcs,
                n_donors     = n_don_sub,
                donor_RE     = if (!is.null(res_dfg))
                                   isTRUE(res_dfg$donor_re) else NA,
                rho          = skcs_info$global_rho,
                tau          = skcs_info$global_tau,
                ks_stat      = skcs_info$ks_stat,
                DE_full_pct  = round(de_full  * 100, 2),
                DE_skcs_pct  = round(de_skcs  * 100, 2),
                DE_skcs_flat_null_pct = round(de_skcs_flat * 100, 2),
                DE_delta_pct = round((de_skcs - de_full) * 100, 2),
                stringsAsFactors = FALSE)
        
            fmt_pct <- function(x) if (is.na(x)) "NA" else sprintf("%.1f%%", x)
            cat(sprintf(
                "  DE Part9=%s  |  DE SKCS-core=%s  |  delta=%s\n\n",
                fmt_pct(de_full  * 100),
                fmt_pct(de_skcs  * 100),
                if (is.na(de_skcs) || is.na(de_full)) "NA"
                else sprintf("%+.1f%%", (de_skcs - de_full) * 100)))
        
            # TREP classification: top round(n * DE) cells by EP, with DE from the same
            # population-level prediction; summary()$dev.expl is kept for reference.
            trep_barcodes <- character(0); ep_df <- data.frame()
            de_for_trep        <- NA_real_
            de_for_trep_legacy <- if (!is.null(res_dfg) &&
                                      !is.null(res_dfg$best_model))
                summary(res_dfg$best_model)$dev.expl else NA_real_

            if (!is.null(res_dfg) && !is.null(res_dfg$best_model) &&
                !is.na(de_skcs) && de_skcs > 0) {
                set.seed(123)
                # Flat null and population-level model prediction (s(donor) excluded).
                null_mdl    <- gam(Expression ~ 1, data = dfg_data)
                null_fitted <- fitted(null_mdl)
                mod_fitted  <- as.numeric(predict(res_dfg$best_model,
                                                  newdata = dfg_data,
                                                  exclude = DONOR_RE_LABEL,
                                                  type = "response"))
                null_sq     <- (dfg_data$Expression - null_fitted)^2
                mod_sq      <- (dfg_data$Expression - mod_fitted)^2
                ep_value    <- ifelse(null_sq > 0, 1 - mod_sq / null_sq, 0)

                # Count and ranking from the same model.
                de_for_trep_count <- pop_dev_explained(res_dfg$best_model,
                                                       dfg_data, mod_fitted)
                if (!is.finite(de_for_trep_count))
                    de_for_trep_count <- de_for_trep_legacy

                n_trep     <- max(1L, min(
                    round(n_sub * trep_count_fraction(de_for_trep_count)),
                    n_sub))
                sorted_idx <- order(ep_value, decreasing = TRUE)
                trep_idx   <- sorted_idx[seq_len(n_trep)]
                trep_barcodes  <- df_sub$barcode[trep_idx]
                de_for_trep    <- de_for_trep_count
            
                ep_df <- data.frame(
                    barcode      = df_sub$barcode,
                    EP_value     = ep_value,
                    Sorted_Index = match(seq_len(n_sub), sorted_idx),
                    is_trep      = seq_len(n_sub) %in% trep_idx,
                    stringsAsFactors = FALSE)
            
                cat(sprintf(
                    "  TREP: top %d / %d (%.1f%%) flagged | EP [%.4f, %.4f]\n",
                    n_trep, n_sub, 100 * de_for_trep_count,
                    min(ep_value, na.rm = TRUE), max(ep_value, na.rm = TRUE)))
                cat(sprintf(
                    "  TREP sizing DE: population-level %.4f | summary()$dev.expl %s\n",
                    de_for_trep_count,
                    if (is.finite(de_for_trep_legacy))
                        sprintf("%.4f", de_for_trep_legacy) else "NA"))
                # Per-donor TREP composition (diagnostic).
                if ("donor" %in% names(dfg_data)) {
                    .d_tot  <- table(droplevels(dfg_data$donor))
                    .d_trep <- table(droplevels(dfg_data$donor[trep_idx]))
                    .d_frac <- as.numeric(.d_trep[names(.d_tot)]) /
                               as.numeric(.d_tot)
                    .d_frac[is.na(.d_frac)] <- 0
                    cat(sprintf(
                        "  TREP donor coverage: %d/%d donors represented | per-donor TREP fraction [%.3f, %.3f] median %.3f\n",
                        sum(.d_frac > 0), length(.d_tot),
                        min(.d_frac), max(.d_frac), stats::median(.d_frac)))
                }
            } else {
                cat("  TREP skipped (no valid DFG GAM or DE <= 0)\n")
            }
        
            list(results = results_out, de_row = de_row,
                 trep    = list(trep_barcodes = trep_barcodes,
                                ep_df         = ep_df,
                                dev_explained = de_for_trep,
                                dev_explained_legacy = de_for_trep_legacy))
        }
    
        # Per-region GAMs
        for (reg in names(region_skcs)) {
            df_skcs <- region_dfs[[reg]][region_skcs[[reg]]$is_skcs, ]
            out <- run_gam_for_subset(df_skcs, reg,
                                      region_skcs[[reg]],
                                      p9_baseline[reg])
            all_skcs_results <- c(all_skcs_results, out$results)
            if (!is.null(out$de_row)) {
                out$de_row$n_full <- nrow(region_dfs[[reg]])
                de_comparison_rows[[reg]] <- out$de_row
            }
            trep_map[[reg]] <- out$trep
        }
    
        # Pooled GAM
        df_skcs_pool   <- df_pooled[skcs_pooled$is_skcs, ]
        skcs_info_pool <- list(global_rho = skcs_pooled$global_rho,
                               global_tau = skcs_pooled$global_tau,
                               ks_stat    = skcs_pooled$ks_stat,
                               n_nskcs    = skcs_pooled$n_nskcs)
        out_pool <- run_gam_for_subset(df_skcs_pool, pooled_key,
                                       skcs_info_pool, p9_baseline_pool)
        all_skcs_results <- c(all_skcs_results, out_pool$results)
        if (!is.null(out_pool$de_row)) {
            out_pool$de_row$n_full <- nrow(df_pooled)
            de_comparison_rows[[pooled_key]] <- out_pool$de_row
        }
        trep_map[[pooled_key]] <- out_pool$trep
    
        # == D7b. TREP labelling and CSV export ===============================
        cat(strrep("-", 60), "\n")
        cat("  SECTION D7b: TREP labelling and CSV export\n")
        cat(strrep("-", 60), "\n\n")
    
        for (reg in names(region_score_dfs)) {
            sd <- region_score_dfs[[reg]]
            tm <- trep_map[[reg]]
            if (!is.null(tm) && length(tm$trep_barcodes) > 0L)
                sd$is_trep <- sd$barcode %in% tm$trep_barcodes
        
            # Append GOI-negative cells with is_GOI_neg = TRUE, analysis cols = NA
            df_all_reg <- region_dfs_all[[reg]]
            goi_neg_bc <- df_all_reg$barcode[df_all_reg$GOI <= 0]
            if (length(goi_neg_bc) > 0) {
                df_neg_reg <- df_all_reg[df_all_reg$GOI <= 0, ]
                sd_neg <- data.frame(
                    barcode       = df_neg_reg$barcode,
                    region        = reg,
                    dataset       = dataset_tag,
                    GOI_expr      = df_neg_reg$GOI,
                    cDFG_expr     = df_neg_reg$Expression,
                    spear_contrib = NA_real_,
                    spear_norm    = NA_real_,
                    kend_contrib  = NA_real_,
                    is_GOI_neg    = TRUE,
                    is_skcs       = FALSE,
                    is_trep       = FALSE,
                    stringsAsFactors = FALSE)
                sd <- rbind(sd, sd_neg)
            }
        
            cat(sprintf(
                "  %-12s  total=%d | GOI_neg=%d | SKCS=%d (TREP=%d, non-TREP=%d) | non-SKCS=%d\n",
                reg, nrow(sd), sum(sd$is_GOI_neg),
                sum(sd$is_skcs), sum(sd$is_trep),
                sum(sd$is_skcs & !sd$is_trep), sum(!sd$is_skcs & !sd$is_GOI_neg)))
        
            safe_write_csv(
                sd,
                file.path(OUT_DIR_SUB,
                          paste0(FILE_PREFIX, "_", dataset_tag,
                                 "_", reg, "_SKCS_scores.csv")),
                row.names = FALSE)
            region_score_dfs[[reg]] <- sd
        }
    
        # Pooled TREP labels
        tm_pool <- trep_map[[pooled_key]]
        if (!is.null(tm_pool) && length(tm_pool$trep_barcodes) > 0L)
            score_pooled_df$is_trep <- score_pooled_df$barcode %in%
            tm_pool$trep_barcodes
    
        # Append GOI-negative cells from all regions into pooled CSV
        goi_neg_rows <- lapply(names(region_dfs_all), function(reg) {
            df_a   <- region_dfs_all[[reg]]
            df_neg <- df_a[df_a$GOI <= 0, ]
            if (nrow(df_neg) == 0) return(NULL)
            neg_df <- data.frame(
                barcode       = df_neg$barcode,
                region        = reg,
                dataset       = dataset_tag,
                GOI_expr      = df_neg$GOI,
                cDFG_expr     = df_neg$Expression,
                spear_contrib = NA_real_,
                spear_norm    = NA_real_,
                kend_contrib  = NA_real_,
                is_GOI_neg    = TRUE,
                is_skcs       = FALSE,
                is_trep       = FALSE,
                stringsAsFactors = FALSE)
            # Add DFG/HKG expression columns (NA for GOI-negative, or actual
            # expression if the columns exist in the source data)
            for (g in dfg_in_pooled)
                neg_df[[g]] <- if (g %in% colnames(df_neg)) df_neg[[g]]
                               else NA_real_
            for (g in hkg_in_pooled)
                neg_df[[g]] <- if (g %in% colnames(df_neg)) df_neg[[g]]
                               else NA_real_
            for (g in bio_in_pooled)
                neg_df[[g]] <- if (g %in% colnames(df_neg)) df_neg[[g]]
                               else NA_real_
            neg_df
        })
        goi_neg_rows <- Filter(Negate(is.null), goi_neg_rows)
        if (length(goi_neg_rows) > 0) {
            score_pooled_df <- rbind(score_pooled_df,
                                     do.call(rbind, goi_neg_rows))
        }
        cat(sprintf(
            "  %-12s  total=%d | GOI_neg=%d | SKCS=%d (TREP=%d, non-TREP=%d) | non-SKCS=%d\n",
            "Pooled", nrow(score_pooled_df),
            sum(score_pooled_df$is_GOI_neg),
            sum(score_pooled_df$is_skcs), sum(score_pooled_df$is_trep),
            sum(score_pooled_df$is_skcs & !score_pooled_df$is_trep),
            sum(!score_pooled_df$is_skcs & !score_pooled_df$is_GOI_neg)))
        safe_write_csv(score_pooled_df, score_pooled_csv, row.names = FALSE)
        cat("  Pooled SKCS scores CSV saved:", basename(score_pooled_csv), "\n\n")
    
        # == D9. Workbook export ==============================================
        sheet_cell_counts_p10 <- build_cell_count_sheet_p10(
            dataset_tag, cell_audit_list, pooled_key)
    
        n_audit_fail <- sum(!is.na(sheet_cell_counts_p10$Tally_OK) &
                                !sheet_cell_counts_p10$Tally_OK)
        if (n_audit_fail > 0)
            cat("  *** CELL COUNT AUDIT:", n_audit_fail, "MISMATCH(ES) ***\n")
        else
            cat("  Cell count audit: all regions PASS\n")
    
        valid_skcs   <- Filter(Negate(is.null), all_skcs_results)
        sheet_perf   <- if (length(valid_skcs) > 0)
            do.call(rbind, lapply(valid_skcs, build_perf_row_p10, group_label))
        else
            data.frame(Note = "No SKCS GAM fits completed",
                       stringsAsFactors = FALSE)
        if (!is.null(sheet_perf) && "Smooth_p_value" %in% names(sheet_perf))
            sheet_perf$Smooth_p_BH_FDR <- p.adjust(sheet_perf$Smooth_p_value,
                                                   "BH")
    
        de_comp_df <- if (length(de_comparison_rows) > 0)
            do.call(rbind, de_comparison_rows)
        else
            data.frame(Note = "No DE comparison data", stringsAsFactors = FALSE)
    
        # TREP summary
        trep_summary_rows <- list()
        for (key in c(names(region_score_dfs), pooled_key)) {
            sd <- if (key == pooled_key) score_pooled_df
            else                   region_score_dfs[[key]]
            if (is.null(sd) || nrow(sd) == 0) next
            tm <- trep_map[[key]]
            trep_summary_rows[[key]] <- data.frame(
                Region           = key,
                Dataset          = dataset_tag,
                Is_Pooled        = (key == pooled_key),
                N_Total          = nrow(sd),
                N_SKCS_core      = sum(sd$is_skcs),
                N_TREP           = sum(sd$is_trep),
                N_SKCS_nonTREP   = sum(sd$is_skcs & !sd$is_trep),
                N_non_SKCS       = sum(!sd$is_skcs),
                DE_used_for_TREP = if (!is.null(tm) && !is.na(tm$dev_explained))
                    tm$dev_explained else NA_real_,
                # summary()$dev.expl, for reference.
                DE_summary_dev_expl =
                    if (!is.null(tm) && !is.null(tm$dev_explained_legacy) &&
                        !is.na(tm$dev_explained_legacy))
                        tm$dev_explained_legacy else NA_real_,
                stringsAsFactors = FALSE)
        }
        trep_summary_df <- if (length(trep_summary_rows) > 0)
            do.call(rbind, trep_summary_rows) else data.frame()
    
        ep_detail_rows <- list()
        for (key in names(trep_map)) {
            tm <- trep_map[[key]]
            if (is.null(tm$ep_df) || nrow(tm$ep_df) == 0) next
            df_ep <- tm$ep_df
            df_ep$Region  <- key; df_ep$Dataset <- dataset_tag
            ep_detail_rows[[key]] <- df_ep[, c("Region", "Dataset", "barcode",
                                               "EP_value", "Sorted_Index",
                                               "is_trep")]
        }
        ep_detail_df <- if (length(ep_detail_rows) > 0)
            do.call(rbind, ep_detail_rows) else data.frame()
    
        safe_sheet <- function(x, note) {
            if (is.null(x) || (is.data.frame(x) && nrow(x) == 0))
                data.frame(Note = paste0("No data (", note, ")"),
                           stringsAsFactors = FALSE)
            else x
        }
    
        xlsx_path <- file.path(OUT_DIR,
                               paste0(FILE_PREFIX, "_", dataset_tag,
                                      "_GAM_Results.xlsx"))
        write_xlsx(list(
            Cell_Count_Audit      = safe_sheet(sheet_cell_counts_p10, "counts"),
            GAM_Performance_SKCS  = safe_sheet(sheet_perf,              "perf"),
            DE_Comparison         = safe_sheet(de_comp_df,              "DE"),
            TREP_Summary          = safe_sheet(trep_summary_df,         "TREP"),
            EP_Detail             = safe_sheet(ep_detail_df,            "EP")),
            path = xlsx_path)
        cat("  Excel saved:", basename(xlsx_path), "\n")
    
        de_csv <- file.path(OUT_DIR,
                            paste0(FILE_PREFIX, "_", dataset_tag,
                                   "_summary_DE_comparison.csv"))
        safe_write_csv(de_comp_df, de_csv, row.names = FALSE)
        cat("  DE comparison CSV saved:", basename(de_csv), "\n\n")
    
        # == D10. RDS =========================================================
        rds_path <- file.path(OUT_DIR,
                              paste0(FILE_PREFIX, "_", dataset_tag,
                                     "_GAM_models.rds"))
        rds_save(list(all_skcs_results  = all_skcs_results,
                      region_skcs       = region_skcs,
                      skcs_pooled       = skcs_pooled,
                      de_comparison     = de_comp_df,
                      sheet_cell_counts = sheet_cell_counts_p10,
                      trep_map          = trep_map,
                      region_score_dfs  = region_score_dfs,
                      score_pooled_df   = score_pooled_df,
                      p9_baseline       = p9_baseline,
                      p9_baseline_pool  = p9_baseline_pool,
                      skcs_s_thresh     = SKCS_S_THRESH,
                      skcs_k_thresh     = SKCS_K_THRESH,
                      min_cells         = MIN_CELLS,
                      dataset_tag       = dataset_tag,
                      regions           = regions,
                      goi               = GOI_10,
                      dfg_genes         = DFG_10,
                      hkg_genes         = HKG_10,
                      module             = paste0("Module 1 (", GOI_10, " / ",
                                                 DFG_COMPOSITE_LABEL, " composite)")),
                 rds_path)
        cat("\n")
    
        # == D11. Plots, single merged PDF per dataset ========================
        # Every plot is written to a temporary PDF page; the pages are merged with qpdf.
        cat(strrep("-", 60), "\n")
        cat("  SECTION D11: Plots -- single merged PDF\n")
        cat(strrep("-", 60), "\n\n")
    
        final_pdf_path <- file.path(OUT_DIR,
                                    paste0(FILE_PREFIX, "_", dataset_tag,
                                           "_GAM_Results_2.pdf"))
    
        tmp_dir <- tempfile("p10_pdf_")
        dir.create(tmp_dir)
        on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)
        pdf_parts <- character(0)
        part_idx  <- 0L
        next_part <- function(suffix) {
            part_idx <<- part_idx + 1L
            file.path(tmp_dir, sprintf("%03d_%s.pdf", part_idx, suffix))
        }
    
        # helper: add a single ggplot page to the temp PDF queue
        add_ggplot_page <- function(p, suffix, pdf_w, pdf_h) {
            if (is.null(p)) return(invisible(NULL))
            path <- next_part(suffix)
            tryCatch({
                pdf(path, width = pdf_w, height = pdf_h)
                print(p); dev.off()
                pdf_parts <<- c(pdf_parts, path)
            }, error = function(e) {
                cat("    PDF page error (", suffix, "):", conditionMessage(e), "\n")
                try(dev.off(), silent = TRUE)
            })
            safe_dev_set_p10()
        }
    
        # helper: add a grid grob page
        add_grob_page <- function(g, suffix, pdf_w, pdf_h) {
            if (is.null(g)) return(invisible(NULL))
            path <- next_part(suffix)
            tryCatch({
                pdf(path, width = pdf_w, height = pdf_h)
                grid::grid.newpage(); grid::grid.draw(g); dev.off()
                pdf_parts <<- c(pdf_parts, path)
            }, error = function(e) {
                cat("    PDF grob error (", suffix, "):", conditionMessage(e), "\n")
                try(dev.off(), silent = TRUE)
            })
            safe_dev_set_p10()
        }
    
        # helper: add a heatmap page
        add_heatmap_page <- function(draw_fn, pdf_w, pdf_h, suffix) {
            if (is.null(draw_fn)) return(invisible(NULL))
            path <- next_part(suffix)
            tryCatch({
                pdf(path, width = pdf_w, height = pdf_h)
                draw_fn(); dev.off()
                pdf_parts <<- c(pdf_parts, path)
            }, error = function(e) {
                cat("    Heatmap PDF error (", suffix, "):",
                    conditionMessage(e), "\n")
                try(dev.off(), silent = TRUE)
            })
            safe_dev_set_p10()
        }
    
        # == (a) SKCS scatter, per region and pooled ==========================
        for (reg in names(region_dfs)) {
            skcs_r    <- region_skcs[[reg]]
            p_scatter <- make_skcs_scatter_p10(
                spear_norm    = skcs_r$spear_norm,
                kend_contrib  = skcs_r$kend_contrib,
                is_skcs       = skcs_r$is_skcs,
                label_title   = paste0(reg, " (", dataset_tag, ")"),
                global_rho    = skcs_r$global_rho,
                global_tau    = skcs_r$global_tau,
                n_skcs        = skcs_r$n_skcs,
                n_nskcs       = skcs_r$n_nskcs)
            add_ggplot_page(p_scatter, paste0("scatter_", reg), 7, 7)
        }
        add_ggplot_page(p_scatter_pool, "scatter_Pooled", 7, 7)
        # Alternate-GOI pooled scatters (HKG8 arm only; empty otherwise)
        for (g_alt in names(alt_scatter_plots))
            add_ggplot_page(alt_scatter_plots[[g_alt]],
                            paste0("scatter_Pooled_GOI_", g_alt), 7, 7)
    
        # == (b) GOI diagnostic, per region (grid) and pooled =================
        for (reg in names(region_dfs)) {
            df_r   <- region_dfs[[reg]]
            skcs_r <- region_skcs[[reg]]
            diag_df_r <- data.frame(
                GOI_expr = df_r$GOI,
                Group    = ifelse(skcs_r$is_skcs, "SKCS-core", "Non-SKCS"),
                stringsAsFactors = FALSE)
            pd <- ggplot(diag_df_r, aes(x = GOI_expr, fill = Group,
                                        colour = Group)) +
                geom_density(alpha = 0.40, linewidth = 0.8) +
                scale_fill_manual(
                    values = c("SKCS-core" = "#4B0082", "Non-SKCS" = "#4F75DE")) +
                scale_colour_manual(
                    values = c("SKCS-core" = "#4B0082", "Non-SKCS" = "#4F75DE")) +
                labs(title = paste0(reg, " -- density"),
                     x = paste0(GOI_10, " (", NORM_TAG, " log1p)"), y = "Density",
                     fill = "Group", colour = "Group") +
                theme_p10(base_size = 10) +
                theme(legend.position = "right")
            add_ggplot_page(pd, paste0("goi_diag_", reg), 7, 5)
        }
        add_ggplot_page(p_dens_pool, "goi_diag_Pooled", 14, 5)
    
        # == (c) Heatmaps =====================================================
        for (hm_key in names(heatmap_draw_fns)) {
            hm_info <- heatmap_draw_fns[[hm_key]]
            add_heatmap_page(hm_info$draw_fn, hm_info$pdf_w, 7, hm_key)
        }
    
        # == (d) PRSS and REML grids (14 x 7) =================================
        keys_vec <- names(valid_skcs)
        fits_pp  <- 4L
        n_pages  <- ceiling(length(keys_vec) / fits_pp)
        for (pg in seq_len(n_pages)) {
            idx_s     <- (pg - 1L) * fits_pp + 1L
            idx_e     <- min(pg * fits_pp, length(keys_vec))
            page_keys <- keys_vec[idx_s:idx_e]
            panel_list <- list()
            for (k in page_keys) {
                res    <- valid_skcs[[k]]
                p_prss <- tryCatch(make_prss_plot_p10(res, compact = TRUE),
                                   error = function(e) NULL)
                p_reml <- tryCatch(make_reml_plot_p10(res, compact = TRUE),
                                   error = function(e) NULL)
                panel_list <- c(panel_list,
                                list(if (!is.null(p_prss)) p_prss else
                                    ggplot()+theme_void()+annotate("text",0.5,0.5,
                                                                   label=paste(k,"PRSS N/A"))),
                                list(if (!is.null(p_reml)) p_reml else
                                    ggplot()+theme_void()+annotate("text",0.5,0.5,
                                                                   label=paste(k,"REML N/A"))))
            }
            g_diag <- tryCatch(
                gridExtra::arrangeGrob(
                    grobs = panel_list, ncol = length(page_keys), nrow = 2L,
                    top   = grid::textGrob(
                        paste0("Module1 SKCS-core PRSS & REML -- ",
                               dataset_tag,
                               "  (page ", pg, " / ", n_pages, ")"),
                        gp = grid::gpar(fontface = "bold", fontsize = 11))),
                error = function(e) NULL)
            add_grob_page(g_diag, sprintf("prss_reml_pg%02d", pg),
                          GRID_W_IN, GRID_H_IN)
        }
    
        # == (e) GAM curve plots (4.86 x 4.00), per region and pooled =========
        cat("  Building GAM curve plots ...\n")
        gs_map    <- c("DFG_SKCS" = DFG_DISPLAY_LABEL,
                       "HKG_SKCS" = "cHKG")
        gs_colors <- setNames(c("#4B0082", "#4F75DE"),
                              c(DFG_DISPLAY_LABEL, "cHKG"))
        gs_fills  <- gs_colors
    
        plot_regions_all <- c(names(region_skcs), pooled_key)
    
        for (reg in plot_regions_all) {
            plot_data <- data.frame(); has_any <- FALSE
        
            for (gs in c("DFG_SKCS", "HKG_SKCS")) {
                key_base <- sub("_SKCS$", "", gs)
                res      <- valid_skcs[[paste0(reg, "|", key_base, "|SKCS")]]
                if (is.null(res) || is.null(res$best_model)) next
                has_any <- TRUE
                m   <- res$best_model; dat <- res$data
                label <- gs_map[gs]
                pts <- data.frame(GOI = dat$GOI, Expression = dat$Expression,
                                  Gene_Set = label, se.fit = NA_real_,
                                  type = "data", row.names = NULL,
                                  stringsAsFactors = FALSE)
                x_seq  <- seq(min(dat$GOI), max(dat$GOI), length.out = 1000L)
                fit_df <- tryCatch({
                    set.seed(123)
                    # Population-level curve, random intercept excluded.
                    pred <- predict_population(m, x_seq, dat)
                    data.frame(GOI = x_seq, Expression = as.numeric(pred$fit),
                               Gene_Set = label, se.fit = as.numeric(pred$se.fit),
                               type = "fit", row.names = NULL,
                               stringsAsFactors = FALSE)
                }, error = function(e) NULL)
                plot_data <- rbind(plot_data, pts)
                if (!is.null(fit_df)) plot_data <- rbind(plot_data, fit_df)
            }
            if (!has_any || nrow(plot_data) == 0) next
        
            de_labs <- sapply(c("DFG_SKCS", "HKG_SKCS"), function(gs) {
                res <- valid_skcs[[paste0(reg, "|", sub("_SKCS$","",gs),
                                          "|SKCS")]]
                if (is.null(res) || is.null(res$best_model))
                    return(NA_character_)
                sprintf("%s DE=%.1f%%", gs_map[gs],
                        de_explained(res$best_model, res$data) * 100)
            })
            de_labs <- de_labs[!is.na(de_labs)]
        
            de_full_val <- if (reg == pooled_key) p9_baseline_pool
            else p9_baseline[reg]
            de_full_lab <- if (!is.na(de_full_val))
                sprintf("Part9 DE=%.1f%%", de_full_val * 100) else ""
        
            n_pts <- max(
                sum(plot_data$type=="data" & plot_data$Gene_Set==DFG_DISPLAY_LABEL),
                sum(plot_data$type=="data" & plot_data$Gene_Set=="cHKG"))
        
            x_data <- plot_data$GOI[plot_data$type == "data"]
            y_data <- plot_data$Expression[plot_data$type == "data"]
            ABS_PAD <- 0.08
            x_lim <- c(min(x_data, na.rm=TRUE) - ABS_PAD,
                       max(x_data, na.rm=TRUE) + ABS_PAD)
            y_lim <- c(max(0, min(y_data, na.rm=TRUE) - ABS_PAD),
                       max(y_data, na.rm=TRUE) + ABS_PAD)
        
            reg_display <- ifelse(reg == pooled_key, "Pooled", reg)
            p_curve <- ggplot(plot_data,
                              aes(x = GOI, y = Expression, colour = Gene_Set)) +
                geom_point(data = subset(plot_data, type == "data"),
                           position = position_jitter(width=0.01, height=0,
                                                      seed=123),
                           alpha = 0.25, size = 0.8) +
                geom_line(data = subset(plot_data, type == "fit"),
                          linewidth = 1.2, na.rm = TRUE) +
                geom_ribbon(data = subset(plot_data, type == "fit"),
                            aes(ymin  = Expression - 1.96 * se.fit,
                                ymax  = Expression + 1.96 * se.fit,
                                fill  = Gene_Set),
                            alpha = 0.20, colour = NA, na.rm = TRUE) +
                scale_colour_manual(values = gs_colors, name = "Gene Set") +
                scale_fill_manual(  values = gs_fills,  name = "Gene Set") +
                scale_x_continuous(labels = number_format(accuracy=0.1),
                                   limits = x_lim) +
                scale_y_continuous(labels = number_format(accuracy=0.1),
                                   limits = y_lim) +
                labs(title    = paste0("GAM (Module1 SKCS-core): ", dataset_tag,
                                       " -- ", reg_display,
                                       "  [n=", n_pts,
                                       ", gamma=", GAM_GAMMA_P10, "]"),
                     subtitle = paste(c(de_labs, de_full_lab),
                                      collapse = "  |  "),
                     x = paste0(GOI_10, " (", NORM_TAG, " log1p)"),
                     y = paste0("Gene Set Expression (", NORM_TAG, " log1p)")) +
                theme_p10(base_size = 8) +
                theme(legend.position = "right",
                      legend.text     = element_text(size = 7),
                      legend.title    = element_text(size = 8, face = "bold"),
                      legend.key.size = unit(0.4, "cm"),
                      plot.title      = element_text(face = "bold", size = 9),
                      plot.subtitle   = element_text(size = 7))
        
            # Add to merged PDF
            add_ggplot_page(p_curve, paste0("gam_curve_", reg),
                            PLOT_W_IN, PLOT_H_IN)
        
            # Individual PNG
            reg_safe <- gsub("[^A-Za-z0-9_]", "_", reg_display)
            png_out  <- if (reg == pooled_key) OUT_DIR else OUT_DIR_SUB
            tryCatch({
                png(file.path(png_out,
                              paste0(FILE_PREFIX, "_", dataset_tag,
                                     "_", reg_safe, "_GAM_curves.png")),
                    width  = round(PLOT_W_IN * 300),
                    height = round(PLOT_H_IN * 300), res = 300)
                print(p_curve); dev.off()
            }, error = function(e) {
                cat("    GAM curve PNG error:", reg_safe, "--",
                    conditionMessage(e), "\n")
                try(dev.off(), silent = TRUE)
            })
            safe_dev_set_p10()
        
            if (reg == pooled_key) emit_plot_pane(p_curve)
        }
    
        # == (e2) GAM curve plots, all GOI-positive cells (4.86 x 4.00) =======
        # As (e), fitted on all GOI-positive cells rather than the SKCS-core subset.
        cat("  Building GAM curve plots (All GOI+) ...\n")
    
        for (reg in plot_regions_all) {
            plot_data_gp <- data.frame(); has_any_gp <- FALSE
            de_labs_gp   <- character(0)
        
            # Select data: pooled or per-region (all GOI-positive cells)
            if (reg == pooled_key) {
                df_gp <- df_pooled
            } else {
                df_gp <- region_dfs[[reg]]
                if (is.null(df_gp)) next
            }
            n_gp <- nrow(df_gp)
            if (n_gp < MIN_CELLS) next
        
            for (gs in c("DFG_SKCS", "HKG_SKCS")) {
                label <- gs_map[gs]
            
                if (gs == "DFG_SKCS") {
                    # Both branches recompute the composite as rowMeans over the DFG columns.
                    if (reg == pooled_key) {
                        dfg_present_gp <- DFG_10[DFG_10 %in% colnames(df_gp)]
                        dfg_missing_gp <- setdiff(DFG_10, dfg_present_gp)
                        dfg_mat_gp     <- as.matrix(
                            df_gp[, dfg_present_gp, drop = FALSE])
                        if (length(dfg_missing_gp) > 0) {
                            pad_gp <- matrix(0, nrow = nrow(dfg_mat_gp),
                                             ncol = length(dfg_missing_gp),
                                             dimnames = list(NULL, dfg_missing_gp))
                            dfg_mat_gp <- cbind(dfg_mat_gp, pad_gp)[, DFG_10,
                                                                    drop = FALSE]
                        }
                        .don_gp <- as_donor_factor(df_gp$donor)
                        gs_data_gp <- data.frame(
                            GOI        = df_gp$GOI,
                            Expression = rowMeans(dfg_mat_gp),
                            donor      = .don_gp,
                            w          = donor_weights(.don_gp))
                    } else {
                        .don_gp <- as_donor_factor(df_gp$donor)
                        gs_data_gp <- data.frame(
                            GOI        = df_gp$GOI,
                            Expression = df_gp$Expression,
                            donor      = .don_gp,
                            w          = donor_weights(.don_gp))
                    }
                } else {
                    # HKG: join from norm matrix
                    if (!has_norm_mtx || length(hkg_pres_hm) == 0) next
                    bc_gp      <- df_gp$barcode
                    idx_gp     <- match(bc_gp, norm_bc)
                    n_match_gp <- sum(!is.na(idx_gp))
                    if (n_match_gp < MIN_CELLS) next
                    hkg_mat_gp <- matrix(0, nrow = n_gp,
                                         ncol = length(hkg_pres_hm),
                                         dimnames = list(bc_gp, hkg_pres_hm))
                    valid_gp <- !is.na(idx_gp)
                    hkg_mat_gp[valid_gp, ] <- as.matrix(
                        expr_mat_full[idx_gp[valid_gp], hkg_pres_hm, drop = FALSE])
                    chkg_gp <- if (length(hkg_pres_hm) == 1)
                        as.numeric(hkg_mat_gp[, 1])
                    else rowMeans(hkg_mat_gp)
                    .don_gp <- as_donor_factor(df_gp$donor)
                    gs_data_gp <- data.frame(
                        GOI        = df_gp$GOI,
                        Expression = chkg_gp,
                        donor      = .don_gp,
                        w          = donor_weights(.don_gp))
                }
            
                # Same two-phase PRSS search as the SKCS-core fits and Part 9.
                res_gp <- tryCatch(
                    fit_gam_prss_p10(gs_data_gp, reg, gs),
                    error = function(e) {
                        cat("    All GOI+ PRSS fit FAILED:", reg, "|", gs,
                            "--", conditionMessage(e), "\n")
                        NULL
                    })
                if (is.null(res_gp) || is.null(res_gp$best_model)) next
                m_gp <- res_gp$best_model
            
                has_any_gp <- TRUE
                de_labs_gp <- c(de_labs_gp,
                                sprintf("%s DE=%.1f%%", label,
                                        de_explained(m_gp, gs_data_gp) * 100))
            
                pts_gp <- data.frame(GOI = gs_data_gp$GOI,
                                     Expression = gs_data_gp$Expression,
                                     Gene_Set = label, se.fit = NA_real_,
                                     type = "data", row.names = NULL,
                                     stringsAsFactors = FALSE)
                x_seq_gp <- seq(min(gs_data_gp$GOI), max(gs_data_gp$GOI),
                                length.out = 1000L)
                fit_df_gp <- tryCatch({
                    set.seed(123)
                    # Population-level curve, random intercept excluded.
                    pred_gp <- predict_population(m_gp, x_seq_gp, gs_data_gp)
                    data.frame(GOI = x_seq_gp,
                               Expression = as.numeric(pred_gp$fit),
                               Gene_Set = label,
                               se.fit = as.numeric(pred_gp$se.fit),
                               type = "fit", row.names = NULL,
                               stringsAsFactors = FALSE)
                }, error = function(e) NULL)
                plot_data_gp <- rbind(plot_data_gp, pts_gp)
                if (!is.null(fit_df_gp))
                    plot_data_gp <- rbind(plot_data_gp, fit_df_gp)
            }
            if (!has_any_gp || nrow(plot_data_gp) == 0) next
        
            # Part 9 baseline DE for subtitle
            de_full_val_gp <- if (reg == pooled_key) p9_baseline_pool
            else p9_baseline[reg]
            de_full_lab_gp <- if (!is.na(de_full_val_gp))
                sprintf("Part9 DE=%.1f%%", de_full_val_gp * 100) else ""
        
            n_pts_gp <- max(
                sum(plot_data_gp$type=="data" & plot_data_gp$Gene_Set==DFG_DISPLAY_LABEL),
                sum(plot_data_gp$type=="data" & plot_data_gp$Gene_Set=="cHKG"))
        
            x_data_gp <- plot_data_gp$GOI[plot_data_gp$type == "data"]
            y_data_gp <- plot_data_gp$Expression[plot_data_gp$type == "data"]
            ABS_PAD_GP <- 0.08
            x_lim_gp <- c(min(x_data_gp, na.rm=TRUE) - ABS_PAD_GP,
                           max(x_data_gp, na.rm=TRUE) + ABS_PAD_GP)
            y_lim_gp <- c(max(0, min(y_data_gp, na.rm=TRUE) - ABS_PAD_GP),
                           max(y_data_gp, na.rm=TRUE) + ABS_PAD_GP)
        
            reg_display_gp <- ifelse(reg == pooled_key, "Pooled", reg)
            p_curve_gp <- ggplot(plot_data_gp,
                                 aes(x = GOI, y = Expression,
                                     colour = Gene_Set)) +
                geom_point(data = subset(plot_data_gp, type == "data"),
                           position = position_jitter(width=0.01, height=0,
                                                      seed=123),
                           alpha = 0.25, size = 0.8) +
                geom_line(data = subset(plot_data_gp, type == "fit"),
                          linewidth = 1.2, na.rm = TRUE) +
                geom_ribbon(data = subset(plot_data_gp, type == "fit"),
                            aes(ymin  = Expression - 1.96 * se.fit,
                                ymax  = Expression + 1.96 * se.fit,
                                fill  = Gene_Set),
                            alpha = 0.20, colour = NA, na.rm = TRUE) +
                scale_colour_manual(values = gs_colors, name = "Gene Set") +
                scale_fill_manual(  values = gs_fills,  name = "Gene Set") +
                scale_x_continuous(labels = number_format(accuracy=0.1),
                                   limits = x_lim_gp) +
                scale_y_continuous(labels = number_format(accuracy=0.1),
                                   limits = y_lim_gp) +
                labs(title    = paste0("GAM (All GOI+): ", dataset_tag,
                                       " -- ", reg_display_gp,
                                       "  [n=", n_pts_gp,
                                       ", gamma=", GAM_GAMMA_P10, "]"),
                     subtitle = paste(c(de_labs_gp, de_full_lab_gp),
                                      collapse = "  |  "),
                     x = paste0(GOI_10, " (", NORM_TAG, " log1p)"),
                     y = paste0("Gene Set Expression (", NORM_TAG, " log1p)")) +
                theme_p10(base_size = 8) +
                theme(legend.position = "right",
                      legend.text     = element_text(size = 7),
                      legend.title    = element_text(size = 8, face = "bold"),
                      legend.key.size = unit(0.4, "cm"),
                      plot.title      = element_text(face = "bold", size = 9),
                      plot.subtitle   = element_text(size = 7))
        
            # Add to merged PDF
            add_ggplot_page(p_curve_gp, paste0("gam_curve_allgoipos_", reg),
                            PLOT_W_IN, PLOT_H_IN)
        
            # Individual PNG
            reg_safe_gp <- gsub("[^A-Za-z0-9_]", "_", reg_display_gp)
            png_out_gp  <- if (reg == pooled_key) OUT_DIR else OUT_DIR_SUB
            tryCatch({
                png(file.path(png_out_gp,
                              paste0(FILE_PREFIX, "_", dataset_tag,
                                     "_", reg_safe_gp,
                                     "_AllGOIpos_GAM_curves.png")),
                    width  = round(PLOT_W_IN * 300),
                    height = round(PLOT_H_IN * 300), res = 300)
                print(p_curve_gp); dev.off()
            }, error = function(e) {
                cat("    All GOI+ GAM curve PNG error:", reg_safe_gp, "--",
                    conditionMessage(e), "\n")
                try(dev.off(), silent = TRUE)
            })
            safe_dev_set_p10()
        
            if (reg == pooled_key) emit_plot_pane(p_curve_gp)
        }
        cat("  All GOI+ GAM curve plots done.\n")
    
        # == (g) DE comparison bar (14 x 7) ===================================
        # Part 9 bars are omitted where the baseline is NA.
        de_bar_plot <- NULL
        if (nrow(de_comp_df) > 0 && "DE_skcs_pct" %in% names(de_comp_df) &&
            any(!is.na(de_comp_df$DE_skcs_pct))) {
        
            # Long format; a missing Part 9 baseline does not drop that region's SKCS bar.
            de_p9   <- data.frame(region = de_comp_df$region,
                                  DE_pct = de_comp_df$DE_full_pct,
                                  type   = "Part 9 baseline",
                                  stringsAsFactors = FALSE)
            lbl_skcs <- paste0(MODULE1_LBL, " SKCS-core")
            de_skcs <- data.frame(region = de_comp_df$region,
                                  DE_pct = de_comp_df$DE_skcs_pct,
                                  type   = lbl_skcs,
                                  stringsAsFactors = FALSE)
            # Drop individual NA rows (missing baseline bars are simply absent)
            de_long <- rbind(de_p9[!is.na(de_p9$DE_pct), ],
                             de_skcs[!is.na(de_skcs$DE_pct), ])
        
            all_regs <- c(sort(regions), pooled_key)
            de_long$region <- factor(de_long$region,
                                     levels = all_regs[all_regs %in%
                                                           de_long$region])
            de_long$type <- factor(de_long$type,
                                   levels = c("Part 9 baseline", lbl_skcs))

            # The y range includes negative DE values.
            de_y_lo   <- min(0, min(de_long$DE_pct, na.rm = TRUE) - 10)
            de_y_hi   <- max(110, max(de_long$DE_pct, na.rm = TRUE) + 10)
            de_breaks <- seq(floor(de_y_lo / 25) * 25,
                             ceiling(de_y_hi / 25) * 25, by = 25)
            de_breaks <- de_breaks[de_breaks >= de_y_lo & de_breaks <= de_y_hi]
            # Labels sit above a positive bar and below a negative one.
            de_long$lab_vjust <- ifelse(de_long$DE_pct >= 0, -0.4, 1.3)
        
            de_bar_plot <- ggplot(de_long,
                                  aes(x = region, y = DE_pct, fill = type)) +
                geom_col(position = position_dodge(width = 0.7), width = 0.6) +
                geom_text(aes(label = sprintf("%.1f%%", DE_pct),
                              vjust = lab_vjust),
                          position = position_dodge(width = 0.7),
                          size = 3, colour = "gray20") +
                geom_hline(yintercept = 50, linetype = "dashed",
                           colour = "gray50", linewidth = 0.5) +
                scale_fill_manual(
                    values = setNames(c("#B0B0B0", "#4B0082"),
                                      c("Part 9 baseline", lbl_skcs)),
                    name = NULL) +
                scale_y_continuous(
                    labels = function(x) paste0(x, "%"),
                    breaks = de_breaks,
                    limits = c(de_y_lo, de_y_hi)) +
                labs(title    = paste0(MODULE1_LBL, " GAM deviance explained -- ",
                                       "Part 9 baseline vs SKCS-core"),
                     subtitle = paste0(dataset_tag, " | ", PRIMARY_STYPE,
                                       " | GOI=", GOI_10,
                                       " | DFG=", DFG_DISPLAY_LABEL,
                                       " composite"),
                     x = "Brain region", y = "Deviance explained (%)") +
                theme_minimal(base_size = 11) +
                theme(
                    plot.title       = element_text(face = "bold", size = 12),
                    plot.subtitle    = element_text(size = 10),
                    axis.title       = element_text(face = "bold", size = 10),
                    axis.title.y     = element_text(face = "bold", size = 10),
                    axis.text        = element_text(size = 9),
                    panel.grid.major = element_line(color = "gray90"),
                    panel.grid.minor = element_blank(),
                    panel.border     = element_rect(color = "black",
                                                    fill = NA, linewidth = 0.35),
                    plot.background  = element_rect(fill = "white", color = NA),
                    panel.background = element_rect(fill = "white", color = NA),
                    legend.position  = "right",
                    legend.text      = element_text(size = 10),
                    legend.key.size  = unit(0.5, "cm"))
            add_ggplot_page(de_bar_plot, "de_bar", 10, 6)   # 10 x 6 in in merged PDF
        
            # DE bar PNG, 10 x 6 in at 200 dpi
            de_bar_png <- file.path(OUT_DIR,
                                    paste0(FILE_PREFIX, "_", dataset_tag,
                                           "_DE_comparison_bar.png"))
            tryCatch({
                png(de_bar_png, width = 2000, height = 1200, res = 200)
                print(de_bar_plot); dev.off()
                cat("  DE comparison bar PNG saved:",
                    basename(de_bar_png), "\n")
            }, error = function(e) {
                cat("    DE bar PNG error:", conditionMessage(e), "\n")
                try(dev.off(), silent = TRUE)
            })
            safe_dev_set_p10()
            emit_plot_pane(de_bar_plot)
        }
        # == CEP-IP scatter, pooled SKCS-core =================================
        cep_ip_scatter_plot <- NULL
    
        res_dfg_pool <- valid_skcs[[paste0(pooled_key, "|DFG|SKCS")]]
        if (!is.null(res_dfg_pool) && !is.null(res_dfg_pool$best_model)) {
            cep_ip_scatter_plot <- make_cep_ip_scatter_p10(
                gam_model   = res_dfg_pool$best_model,
                sample_data = res_dfg_pool$data,
                dataset_tag = dataset_tag,
                panel_tag   = "Module1 SKCS-core")
        
            # Into merged PDF + standalone PNG
            add_ggplot_page(cep_ip_scatter_plot, "cep_ip_scatter",
                            PLOT_W_IN, PLOT_H_IN)
        
            # CEP-IP scatter PNG
            tryCatch({
                png(file.path(OUT_DIR,
                              paste0(FILE_PREFIX, "_", dataset_tag,
                                     "_Pooled_CEP_IP_scatter.png")),
                    width  = round(PLOT_W_IN * 300),
                    height = round(PLOT_H_IN * 300), res = 300)
                print(cep_ip_scatter_plot); dev.off()
                cat("  CEP-IP scatter PNG saved (Module1 SKCS-core)\n")
            }, error = function(e) {
                cat("    CEP-IP PNG error:", conditionMessage(e), "\n")
                try(dev.off(), silent = TRUE)
            })
            safe_dev_set_p10()
            emit_plot_pane(cep_ip_scatter_plot)
        } else {
            cat("  CEP-IP scatter skipped (no Pooled DFG GAM)\n")
        }
    
        # == CEP-IP scatter, pooled Part 9 baseline ===========================
        # Same plot on the Part 9 pooled full-population GAM.
        cep_ip_scatter_p9_plot <- NULL
    
        if (!is.null(p9_pooled_model) && !is.null(p9_pooled_data) &&
            all(c("GOI", "Expression") %in% names(p9_pooled_data))) {
            cep_ip_scatter_p9_plot <- make_cep_ip_scatter_p10(
                gam_model   = p9_pooled_model,
                sample_data = p9_pooled_data,
                dataset_tag = dataset_tag,
                panel_tag   = "Part 9 baseline")
        
            # Into merged PDF + standalone PNG
            add_ggplot_page(cep_ip_scatter_p9_plot, "cep_ip_scatter_Part9",
                            PLOT_W_IN, PLOT_H_IN)
        
            tryCatch({
                png(file.path(OUT_DIR,
                              paste0(FILE_PREFIX, "_", dataset_tag,
                                     "_Pooled_CEP_IP_scatter_Part9_baseline.png")),
                    width  = round(PLOT_W_IN * 300),
                    height = round(PLOT_H_IN * 300), res = 300)
                print(cep_ip_scatter_p9_plot); dev.off()
                cat("  CEP-IP scatter PNG saved (Part 9 baseline)\n")
            }, error = function(e) {
                cat("    CEP-IP (Part 9) PNG error:", conditionMessage(e), "\n")
                try(dev.off(), silent = TRUE)
            })
            safe_dev_set_p10()
            emit_plot_pane(cep_ip_scatter_p9_plot)
        } else {
            cat("  CEP-IP scatter (Part 9 baseline) skipped",
                "(no Part 9 pooled GAM / data)\n")
        }
    
        # == Merge all PDF parts ==============================================
        if (length(pdf_parts) == 0) {
            cat("  No PDF parts -- final PDF will not be produced.\n")
        } else {
            tryCatch({
                qpdf::pdf_combine(input = pdf_parts, output = final_pdf_path)
                cat("  Merged PDF saved:", basename(final_pdf_path),
                    "(", length(pdf_parts), "parts,",
                    "mixed page sizes)\n")
            }, error = function(e) {
                cat("  qpdf merge failed:", conditionMessage(e), "\n")
                cat("  Copying first part as fallback.\n")
                file.copy(pdf_parts[1], final_pdf_path, overwrite = TRUE)
            })
        }
    
        # == D12. Summary =====================================================
        cat("\n", strrep("=", 72), "\n")
        cat("  Part 10", dataset_tag, "COMPLETE\n")
        cat(strrep("=", 72), "\n\n")
        cat("  DE summary (Part9 baseline vs Module1 SKCS-core):\n")
        if (nrow(de_comp_df) > 0 && "DE_full_pct" %in% names(de_comp_df))
            for (i in seq_len(nrow(de_comp_df))) {
                r <- de_comp_df[i, ]
                cat(sprintf(
                    "    %-12s  n_full=%-5s  n_skcs=%-5d  KS=%s  DE_p9=%s  DE_skcs=%s  delta=%s\n",
                    r$region,
                    ifelse(is.na(r$n_full), "NA", as.character(r$n_full)),
                    r$n_skcs,
                    ifelse(is.na(r$ks_stat), "NA", sprintf("%.3f", r$ks_stat)),
                    ifelse(is.na(r$DE_full_pct),  "NA%",
                           sprintf("%.1f%%", r$DE_full_pct)),
                    ifelse(is.na(r$DE_skcs_pct),  "NA%",
                           sprintf("%.1f%%", r$DE_skcs_pct)),
                    ifelse(is.na(r$DE_delta_pct), "NA%",
                           sprintf("%+.1f%%", r$DE_delta_pct))))
            }
    
        invisible(list(all_skcs_results  = all_skcs_results,
                       region_skcs       = region_skcs,
                       skcs_pooled       = skcs_pooled,
                       de_comparison     = de_comp_df,
                       sheet_cell_counts = sheet_cell_counts_p10,
                       trep_map          = trep_map,
                       region_score_dfs  = region_score_dfs,
                       score_pooled_df   = score_pooled_df,
                       plots_pane        = list(
                           pooled_scatter  = p_scatter_pool,
                           pooled_goi_diag = p_dens_pool,
                           de_bar          = de_bar_plot,
                           cep_ip_scatter  = cep_ip_scatter_plot,
                           cep_ip_scatter_p9 = cep_ip_scatter_p9_plot)))
    }


    # =============================================================================
    # SECTION E: Execute
    # =============================================================================

    cat("\n\n", strrep("#", 72), "\n")
    cat(sprintf("  RUNNING: NonNeoC  [%s / %s composite]\n", GOI_10, DFG_COMPOSITE_LABEL))
    cat(strrep("#", 72), "\n")
    set.seed(42)
    p10_nonnc <- run_p10_pipeline(
        dataset_tag = "NonNeoC",
        dataset_dir = NONNC_DIR,
        regions     = NONNC_REGIONS_P10,
        p9_rds_path = PART9_NONNC_RDS,
        dfg_genes   = DFG_10,
        hkg_genes   = HKG_10
    )

    cat("\n\n", strrep("#", 72), "\n")
    cat(sprintf("  RUNNING: NeoC  [%s / %s composite]\n", GOI_10, DFG_COMPOSITE_LABEL))
    cat(strrep("#", 72), "\n")
    set.seed(42)
    p10_neoc <- run_p10_pipeline(
        dataset_tag = "NeoC",
        dataset_dir = NEOC_DIR,
        regions     = NEOC_REGIONS_P10,
        p9_rds_path = PART9_NEOC_RDS,
        dfg_genes   = DFG_10,
        hkg_genes   = HKG_10
    )


    # =============================================================================
    # SECTION F: Non-SKCS second pass (placeholder)
    # =============================================================================
    cat("\n", strrep("=", 72), "\n")
    cat("SECTION F: Non-SKCS / non-TREP second-pass GAM (optional)\n")
    cat(strrep("=", 72), "\n\n")
    cat("  Not run in this script. Minimum cells:", MIN_CELLS, "\n\n")


    # =============================================================================
    # SECTION G: Combined RDS and summary text
    # =============================================================================

    cat("\n", strrep("=", 72), "\n")
    cat(sprintf("  Part 10 OVERALL COMPLETE (%s / %s)\n", GOI_10, DFG_COMPOSITE_LABEL))
    cat(strrep("=", 72), "\n\n")
    cat("  Arm dir (log, guardrails, combined RDS, summary):", OUT_DIR, "\n")
    for (d in DS_DIRS) cat("  Dataset dir                                   :",
                           d, "\n")
    cat("\n")

    p10_combined_rds <- file.path(OUT_DIR,
                                  paste0(OUT_PREFIX, "_combined_results.rds"))
    rds_save(
        list(NonNeoC       = p10_nonnc,
             NeoC          = p10_neoc,
             skcs_s_thresh = SKCS_S_THRESH,
             skcs_k_thresh = SKCS_K_THRESH,
             min_cells     = MIN_CELLS,
             gam_gamma     = GAM_GAMMA_P10,
             goi           = GOI_10,
             dfg_genes     = DFG_10,
             hkg_genes     = HKG_10,
             dfg_label     = DFG_COMPOSITE_LABEL,
             module         = paste0(GOI_10, " / ", DFG_COMPOSITE_LABEL, " composite"),
             timestamp     = Sys.time()),
        p10_combined_rds, label = "Overall RDS")

    txt_path <- file.path(OUT_DIR,
                          paste0(OUT_PREFIX, "_summary_rationale.txt"))

    collect_counts <- function(result_obj, dataset_tag, regions_vec) {
        if (is.null(result_obj)) return(character(0))
        lines <- character(0)
        for (reg in regions_vec) {
            skcs_r <- result_obj$region_skcs[[reg]]
            if (is.null(skcs_r)) next
            n_total <- skcs_r$n_skcs + skcs_r$n_nskcs
            lines <- c(lines, sprintf(
                "    %-12s  total=%-5d  SKCS-core=%-5d (%.1f%%)  non-SKCS=%-5d (%.1f%%)",
                reg, n_total,
                skcs_r$n_skcs,  100 * skcs_r$n_skcs  / n_total,
                skcs_r$n_nskcs, 100 * skcs_r$n_nskcs / n_total))
        }
        sc <- result_obj$skcs_pooled
        if (!is.null(sc)) {
            n_tot_c <- sc$n_skcs + sc$n_nskcs
            lines <- c(lines, sprintf(
                "    %-12s  total=%-5d  SKCS-core=%-5d (%.1f%%)  non-SKCS=%-5d (%.1f%%)",
                "Pooled", n_tot_c,
                sc$n_skcs,  100 * sc$n_skcs  / n_tot_c,
                sc$n_nskcs, 100 * sc$n_nskcs / n_tot_c))
        }
        lines
    }

    lines_nonnc <- collect_counts(p10_nonnc, "NonNeoC", NONNC_REGIONS_P10)
    lines_neoc  <- collect_counts(p10_neoc,  "NeoC",    NEOC_REGIONS_P10)

    txt_lines <- c(
        strrep("=", 72),
        paste0("Part 10 - ", MODULE1_LBL,
               " SKCS Decomposition: Rationale & Cell Count Summary"),
        strrep("=", 72), "",
        toupper(MODULE1_LBL),
        strrep("-", nchar(MODULE1_LBL)),
        paste0("  GOI  : ", GOI_10),
        paste0("  DFG  : ", DFG_COMPOSITE_LABEL, " composite of  ",
               paste(DFG_10, collapse = ", ")),
        paste0("  HKG  : ", paste(HKG_10, collapse = ", ")),
        paste0("  Cell source : per-region gam_data CSVs (full M4-top)"),
        paste0("  GOI filter  : ", GOI_10, " > 0"),
        paste0("  Baseline DE : Part 9 RDS (full-population GAM)"), "",
        "SKCS THRESHOLDS",
        "---------------",
        paste0("  S_norm > ", SKCS_S_THRESH,
               " (calibration-purity floor: >= 5% of max single-cell rank-product)"),
        paste0("  K > ", SKCS_K_THRESH,
               " (concordant with >= 75% of all pairwise comparisons)"),
        "  Normalisation: Si_norm = Si / Si_max, where Si_max = (n-1)/4",
        "  Si_norm > 0.05 excludes ambiguous IP transition zone for CMRF calibration", "",
        "TREP CLASSIFICATION (EP-based)",
        "------------------------------",
        "  EP = 1 - model_residual^2 / null_residual^2  (per cell)",
        "  Model residual uses the population-level prediction (s(donor) excluded).",
        "  TREP = top round(n x DE) cells by EP descending, where DE is the",
        "  weighted population-level deviance explained of that same prediction.",
        "  Non-SKCS cells: is_trep = FALSE by default.", "",
        "CELL COUNTS",
        "-----------",
        "NON-NEOCORTICAL REGIONS (NonNeoC)",
        strrep("-", 40))
    txt_lines <- c(txt_lines, lines_nonnc, "",
                   "NEOCORTICAL REGIONS (NeoC)",
                   strrep("-", 40),
                   lines_neoc, "",
                   paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
                   strrep("=", 72))

    writeLines(txt_lines, txt_path)
    cat("  Summary txt saved:", basename(txt_path), "\n\n")

    # File listing
    cat("  Files generated:\n")
    # The arm folder plus each dataset folder and its subregions folder.
    .listing_dirs <- unique(c(OUT_DIR, DS_DIRS,
                              file.path(DS_DIRS, "Part10_subregions")))
    .listing_dirs <- .listing_dirs[dir.exists(.listing_dirs)]
    all_files <- unlist(lapply(.listing_dirs, function(d)
        file.path(d, list.files(d, recursive = FALSE))))
    for (f in sort(all_files)) {
        if (isTRUE(file.info(f)$isdir)) next
        sz  <- round(file.size(f) / 1e3, 0)
        rel <- sub(paste0(R_EXPORT_ROOT, .Platform$file.sep), "", f, fixed = TRUE)
        cat(sprintf("    %-80s  %6d KB\n", rel, sz))
    }


    # =============================================================================
    # SECTION H: Global cell count guardrail
    # =============================================================================

    cat("\n", strrep("#", 72), "\n")
    cat("  SECTION H -- GLOBAL CELL COUNT GUARDRAIL\n")
    cat(strrep("#", 72), "\n\n")

    all_counts_p10 <- rbind(
        if (!is.null(p10_nonnc$sheet_cell_counts))
            p10_nonnc$sheet_cell_counts else NULL,
        if (!is.null(p10_neoc$sheet_cell_counts))
            p10_neoc$sheet_cell_counts  else NULL)

    if (is.null(all_counts_p10) || nrow(all_counts_p10) == 0) {
        cat("  WARNING: no cell count audit data found.\n\n")
    } else {
        cat("  Full cell count audit:\n\n")
        print(all_counts_p10, row.names = FALSE)
    
        failed_p10 <- all_counts_p10[
            !is.na(all_counts_p10$Tally_OK) &
                !all_counts_p10$Tally_OK, ]
        cat("\n")
        if (nrow(failed_p10) == 0) {
            cat(strrep("*", 72), "\n")
            cat("  GUARDRAIL PASSED: N_Used_for_Analysis == N_GOI_Pos\n")
            cat("  for ALL regions in BOTH NeoC and NonNeoC datasets.\n")
            cat(strrep("*", 72), "\n\n")
        } else {
            cat(strrep("!", 72), "\n")
            cat("  *** GUARDRAIL FAILED ***\n")
            print(failed_p10, row.names = FALSE)
            cat(strrep("!", 72), "\n\n")
            warning("Part10 guardrail FAILED for ",
                    nrow(failed_p10), " region(s).", call. = FALSE)
        }
    }

    cat("  Part 10 all done.\n\n")


    # =============================================================================
    # SECTION I: Post-analysis quality guardrails
    # =============================================================================
    cat("\n", strrep("#", 72), "\n")
    cat("  SECTION I -- POST-ANALYSIS QUALITY GUARDRAILS\n")
    cat(strrep("#", 72), "\n\n")

    # Collect SKCS statistics across datasets
    for (.ds_res in list(
        list(tag = "NeoC",    result = p10_neoc,  regions = NEOC_REGIONS_P10),
        list(tag = "NonNeoC", result = p10_nonnc, regions = NONNC_REGIONS_P10)
    )) {
        if (is.null(.ds_res$result)) next
        .res <- .ds_res$result

        # == GR-08: SKCS-core cells exist in every region =====================
        for (reg in .ds_res$regions) {
            skcs_r <- .res$region_skcs[[reg]]
            if (is.null(skcs_r)) next
            n_skcs  <- skcs_r$n_skcs
            n_total <- n_skcs + skcs_r$n_nskcs
            gr10(sprintf("GR-08-%s-%s", .ds_res$tag, reg),
                 if (n_skcs >= MIN_CELLS) "PASS"
                 else if (n_skcs > 0) "WARN" else "FAIL",
                 sprintf("SKCS-core=%d/%d (%.1f%%)",
                         n_skcs, n_total,
                         if (n_total > 0) 100 * n_skcs / n_total else 0))
        }

        # == GR-09: SKCS-core GAM DE > 0 for pooled ===========================
        .skcs_gam <- Filter(Negate(is.null), .res$all_skcs_results %||% list())
        if (length(.skcs_gam) > 0) {
            .de_skcs <- sapply(.skcs_gam, function(r) {
                if (is.null(r$best_model)) NA_real_
                else tryCatch(de_explained(r$best_model, r$data),
                              error = function(e) NA_real_)
            })
            .de_valid <- .de_skcs[!is.na(.de_skcs)]
            if (length(.de_valid) > 0) {
                gr10(sprintf("GR-09-%s", .ds_res$tag),
                     if (all(.de_valid >= 0 & .de_valid <= 1)) "PASS" else "WARN",
                     sprintf("SKCS-core GAM DE range: [%.1f%%, %.1f%%]  (n=%d, donor-aware null)",
                             min(.de_valid) * 100, max(.de_valid) * 100,
                             length(.de_valid)))
            }

            # == GR-11: donor random intercept fitted on every model ==========
            .re_fitted <- vapply(.skcs_gam,
                                 function(r) isTRUE(r$donor_re), logical(1))
            .n_don_all <- vapply(.skcs_gam, function(r)
                if (is.null(r$n_donors)) NA_integer_
                else as.integer(r$n_donors), integer(1))
            .fallbacks <- names(.skcs_gam)[!.re_fitted]
            if (!isTRUE(GAM_DONOR_RE)) {
                gr10(sprintf("GR-11-%s", .ds_res$tag), "WARN",
                     paste0("GAM_DONOR_RE is FALSE: all ", length(.skcs_gam),
                            " fit(s) are cell level."))
            } else {
                gr10(sprintf("GR-11-%s", .ds_res$tag),
                     if (length(.fallbacks) == 0L) "PASS" else "WARN",
                     sprintf(
                         "donor RE on %d/%d fit(s); donors per fit [%s, %s]; %d fallback(s)%s",
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

        # == GR-10: TREP cell counts are reasonable ===========================
        if (!is.null(.res$trep_map)) {
            for (reg in names(.res$trep_map)) {
                .tm <- .res$trep_map[[reg]]
                if (is.null(.tm) || is.null(.tm$ep_df) ||
                    nrow(.tm$ep_df) == 0L) next
                n_trep <- sum(.tm$ep_df$is_trep)
                n_tot  <- nrow(.tm$ep_df)
                gr10(sprintf("GR-10-%s-%s", .ds_res$tag, reg),
                     if (n_trep > 0 && n_trep <= n_tot) "PASS" else "WARN",
                     sprintf("TREP=%d/%d (%.1f%%) | sizing DE=%s",
                             n_trep, n_tot,
                             if (n_tot > 0) 100 * n_trep / n_tot else 0,
                             if (!is.null(.tm$dev_explained) &&
                                 is.finite(.tm$dev_explained))
                                 sprintf("%.4f", .tm$dev_explained) else "NA"))
            }
        }
    }

    # == Update guardrail log with post-analysis checks =======================
    gr10_write()
    cat(sprintf("\n  Final guardrail tally: %d PASS / %d WARN / %d FAIL\n",
                .gr10_pass, .gr10_warn, .gr10_fail))
    if (.gr10_fail > 0) {
        cat("  Failures present -- review the ledger before interpreting.\n")
    } else if (.gr10_warn > 0) {
        cat(sprintf("  %d passed, %d warning(s), no failures.\n",
                    .gr10_pass, .gr10_warn))
    } else {
        cat(sprintf("  All %d checks passed.\n", .gr10_pass))
    }

    cat("\n  Log written to:", LOG_PATH_P10,
        if (isTRUE(LOG_MESSAGES_P10)) " (console output and warnings)"
        else " (console output only)", "\n")
    cat("  Part 10 all done:", format(Sys.time()), "\n\n")

    # Released here on the normal path; on.exit() covers an early abort.
    if (!is.null(.p10_cl)) {
        try(parallel::stopCluster(.p10_cl), silent = TRUE)
        .p10_cl <- NULL
        cat("  PSOCK cluster stopped.\n")
    }
    .p10_completed <- TRUE

    if (sink.number(type = "message") != 2L) sink(type = "message")
    while (sink.number() > 0) sink()

    # The arm reached the end, so the folder holds a complete run.
    suppressWarnings(file.remove(INCOMPLETE_MARKER))

    invisible(list(arm  = arm_tag, dir = OUT_DIR,
                   pass = .gr10_pass, warn = .gr10_warn, fail = .gr10_fail))
}


# =============================================================================
# SECTION 0c: Driver
# =============================================================================
# The arm loop and its summary live in CEP_arms.R.

cep_run_arms(PART_TAG, run_one_arm)
