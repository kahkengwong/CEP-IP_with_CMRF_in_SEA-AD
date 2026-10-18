# =============================================================================
# Part10b_CP10K_IPRS_code.r
# SEA-AD CEP-IP/CMRF, Part 10b IPRS, inflection point detection
#   (Figure 3, Supplementary Figure S6, Table S8)
#   For each arm: fits a GAM of the partner composite on the GOI over all
#   GOI-positive cells pooled from the Part 10 RDS (donor random intercept,
#   inverse-donor weights), classifies purple cells by explanatory power, fits
#   the residual GAM R(x), takes its first neg->pos zero-crossing as the IP with
#   a Brent 95% CI, and scores it with IPRS = (C1 + C2 + C3 + C4 + C5) / 5.
#
# Usage (from the SEA-AD folder, or set SEAAD_ROOT), after Part 10
#   source("Part10b_CP10K_IPRS_code.r")                         # all arms
#   CEP_ARM <- "PVM1_all"; source("Part10b_CP10K_IPRS_code.r")  # one arm
#   Needs CEP_arms.R, R_config.R and R_config_HKG8.R in the same folder.
#
# Outputs (Part10b_IPRS/Part10b_IPRS_<arm>/): results RDS (and a copy for
#   Part 10b SKCS_TREP), workbook, summary CSV, scatter, R(x) and pitcher
#   figures (PDF, PNG), RDS verification log, guardrail log and run log.
# =============================================================================


# =============================================================================
# SECTION 0: Libraries and shared arm definitions
# =============================================================================

library(mgcv)
library(Matrix)        # readMM for norm.mtx
library(ggplot2)
library(dplyr)
library(writexl)
library(scales)
library(MASS)          # kde2d for the expression-space density background
library(viridis)       # turbo palette for the Ki/Si colouring
library(ggnewscale)    # dual fill scale: density background + point fill

set.seed(123)

# CEP_arms.R holds the arm table, the folder layout and the driver, shared with
# Parts 9 onwards. Set CEP_HOME before sourcing this file to point elsewhere.
if (!exists("CEP_HOME")) CEP_HOME <- Sys.getenv("SEAAD_ROOT", getwd())
.iprs_harness <- file.path(CEP_HOME, "CEP_arms.R")
if (!file.exists(.iprs_harness))
    stop("CEP_arms.R not found at: ", .iprs_harness,
         "\n  It holds the arm table shared by Parts 9 onwards.",
         call. = FALSE)
source(.iprs_harness)

PART_TAG <- "Part10b_IPRS"

# Self-tests (Section Z). Leave FALSE for a production run.
RUN_SMOKE_TEST  <- FALSE
RUN_STRESS_TEST <- FALSE

# Label mgcv gives the donor random intercept.
DONOR_RE_LABEL <- "s(donor)"

# RStudio Plots pane. No dev.new(): RStudio routes print() to its own pane.
emit_plot_iprs <- function(p) {
    print(p)
}


# =============================================================================
# SECTION 0b: The pipeline for one arm
# =============================================================================

run_one_arm <- function(arm_tag, config_file, out_dir, file_tag) {

    config_path <- file.path(CEP_HOME, config_file)
    if (!file.exists(config_path))
        stop("Configuration not found: ", config_path, call. = FALSE)

    # Output printed before the log exists is buffered and replayed into it.
    .iprs_preamble <- character(0)

    assign("CEP_ARM", arm_tag, envir = globalenv())
    on.exit(suppressWarnings(rm("CEP_ARM", envir = globalenv())), add = TRUE)

    .iprs_preamble <- tryCatch(
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

    .iprs_needed <- c("ARM_TAG", "ARM_DIR", "R_EXPORT_ROOT", "GOI", "DFG",
                      "HKG_GENES", "DFG_COMPOSITE_LABEL", "RUN_SUFFIX",
                      "PRIMARY_STYPE", "NORM_TAG", "NORM_TARGET_SUM",
                      "NORM_CEILING", "NEOC_DIR", "NONNC_DIR",
                      "NEOC_REGIONS", "NONNC_REGIONS", "DFG_NEOC", "DFG_NONNC",
                      "GAM_BS", "GAM_METHOD", "GAM_SELECT")
    .iprs_absent <- .iprs_needed[!vapply(.iprs_needed, exists, logical(1))]
    if (length(.iprs_absent) > 0)
        stop("The configuration loaded from\n    ", config_path,
             "\n  does not define: ", paste(.iprs_absent, collapse = ", "),
             call. = FALSE)

    if (!identical(ARM_TAG, arm_tag))
        stop("The configuration reports ARM_TAG '", ARM_TAG,
             "' but this arm is '", arm_tag, "'.", call. = FALSE)

    .iprs_preamble <- c(.iprs_preamble, utils::capture.output(
        cat("  Config OK: ", ARM_TAG, " | ", GOI, " + ", length(DFG), " DFG | ",
            NORM_TAG, "\n", sep = "")))

    set.seed(123)

    # == User parameters =====================================================
    DATASETS <- c("NonNeoC", "NeoC")

    # GAM_BS, GAM_METHOD and GAM_SELECT come from R_config.R.
    GAM_GAMMA_IPRS <- 1.5

    # == Estimator ===========================================================
    # The donor is the experimental unit; FALSE gives the cell-level fit.
    GAM_DONOR_RE        <- TRUE
    DONOR_RE_MIN_LEVELS <- 3L
    # DONOR_RE_LABEL is defined at the top level in Section 0.

    # Donor random intercept on the residual GAM, switchable separately.
    RESID_DONOR_RE      <- TRUE

    # IP_y = f(IP_x): the Pre/Post classification threshold on the response axis.
    # FALSE: IP_y is still computed and exported, but no line is drawn.
    IP_Y_ENABLE         <- TRUE

    # IP_y line style for the scatter and pitcher plots.
    IP_Y_LINE_COL       <- "grey10"
    IP_Y_LINE_TYPE      <- "dashed"

    # Inverse-donor-count weights. FALSE gives every cell w = 1.
    W_DONOR_WEIGHTS     <- TRUE

    # == IP definition =======================================================
    # IPRS R(x)=0 zero-crossing with Brent CI, on the full GOI-positive population.


    STAT_ALPHA        <- 0.05
    SHAPIRO_ALPHA     <- 0.05
    POST_IP_MIN_CELLS <- 5L
    MIN_PURPLE_CELLS  <- 10L

    IPRS_VERY_STRONG <- 0.80
    IPRS_STRONG      <- 0.60
    IPRS_MODERATE    <- 0.40
    IPRS_WEAK        <- 0.20

    # == Part 10 input location ==============================================
    PART10_DIR <- cep_arm_dir_for("Part10", arm_tag)
    P10_PREFIX <- cep_out_prefix("Part10", NORM_TAG,
                                 cep_file_tag_for(arm_tag))
    P10_FILE_PREFIX <- paste0(P10_PREFIX, "_Module1_", gsub("-", "", GOI), "_",
                              DFG_COMPOSITE_LABEL)

    # The per-dataset folder is tried first, then the arm root.
    p10_rds_candidates <- function(ds) {
        fn <- paste0(P10_FILE_PREFIX, "_", ds, "_GAM_models.rds")
        c(file.path(PART10_DIR, paste0(basename(PART10_DIR), "_", ds), fn),
          file.path(PART10_DIR, fn))
    }
    p10_rds_path <- function(ds) {
        cand <- p10_rds_candidates(ds)
        hit  <- cand[file.exists(cand)]
        if (length(hit) > 0L) hit[1] else cand[1]
    }

    # == Output directory and file naming ====================================
    OUT_DIR <- out_dir
    dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

    OUT_PREFIX <- cep_out_prefix(PART_TAG, NORM_TAG,
                                 cep_file_tag_for(arm_tag))

    # Where Part 10b reads the RDS this part writes. A non-canonical PART_TAG
    # redirects the copy to a matching side folder.
    PART_TAG_CANON  <- "Part10b_IPRS"
    .p10b_suffix    <- if (identical(PART_TAG, PART_TAG_CANON)) ""
                       else if (startsWith(PART_TAG, PART_TAG_CANON))
                           substring(PART_TAG, nchar(PART_TAG_CANON) + 1L)
                       else paste0("_", PART_TAG)
    P10B_DIR_TAG    <- paste0("Part10b_SKCS", .p10b_suffix)  # Part 10b's PART_DIR_TAG
    P10B_FILE_TAG   <- paste0("Part10b", .p10b_suffix)       # Part 10b's PART_TAG
    if (nzchar(.p10b_suffix))
        cat("  PART_TAG is not \"", PART_TAG_CANON,
            "\", so the Part 10b copy goes to ",
            P10B_DIR_TAG, "/.\n", sep = "")
    P10B_DIR        <- cep_arm_dir_for(P10B_DIR_TAG, arm_tag)
    P10B_PREFIX     <- cep_out_prefix(P10B_FILE_TAG, NORM_TAG,
                                      cep_file_tag_for(arm_tag))
    P10B_IPRS_RDS   <- file.path(P10B_DIR,
                                 paste0(P10B_PREFIX, "_IPRS_results.rds"))

    # == Per-dataset input metadata ==========================================
    DATASET_DIRS <- list(NonNeoC = NONNC_DIR, NeoC = NEOC_DIR)
    DATASET_REGS <- list(NonNeoC = NONNC_REGIONS, NeoC = NEOC_REGIONS)
    DATASET_DFG  <- list(NonNeoC = DFG_NONNC, NeoC = DFG_NEOC)

    # == Log file ============================================================
    # Messages go to a temporary file that is folded into the log at the end.
    LOG_MESSAGES_IPRS <- TRUE
    LOG_PATH_IPRS <- file.path(OUT_DIR,
                               paste0(OUT_PREFIX, "_Log_",
                                      format(Sys.time(), "%Y%m%d_%H%M%S"),
                                      ".txt"))
    MSG_PATH_IPRS <- file.path(OUT_DIR,
                               paste0(OUT_PREFIX, "_Log_messages.tmp"))

    while (sink.number() > 0) sink()
    # sink.number(type = "message") is 2 when nothing is diverted.
    if (sink.number(type = "message") != 2L) sink(type = "message")

    .iprs_log_con <- file(LOG_PATH_IPRS, open = "wt")
    .iprs_msg_con <- if (isTRUE(LOG_MESSAGES_IPRS))
        file(MSG_PATH_IPRS, open = "wt") else NULL
    sink(.iprs_log_con, append = FALSE, split = TRUE)
    if (!is.null(.iprs_msg_con)) sink(.iprs_msg_con, type = "message")

    # Flag an incomplete run in the log; the error itself prints to the console.
    .iprs_completed <- FALSE

    on.exit({
        tryCatch({
            if (!isTRUE(.iprs_completed))
                cat("\n", strrep("!", 72),
                    "\n  Arm did not complete; the error is printed in the console.\n",
                    strrep("!", 72), "\n", sep = "")
            if (sink.number(type = "message") != 2L) sink(type = "message")
            while (sink.number() > 0) sink()
            if (!is.null(.iprs_msg_con)) close(.iprs_msg_con)
            close(.iprs_log_con)
            if (!is.null(.iprs_msg_con) && file.exists(MSG_PATH_IPRS)) {
                .msg_txt <- readLines(MSG_PATH_IPRS, warn = FALSE)
                if (length(.msg_txt) > 0)
                    cat(c("", strrep("=", 72),
                          "  MESSAGE STREAM (warnings, messages, errors)",
                          strrep("=", 72), .msg_txt, ""),
                        sep = "\n", file = LOG_PATH_IPRS, append = TRUE)
                unlink(MSG_PATH_IPRS)
            }
        }, error = function(e) invisible(NULL))
    }, add = TRUE)

    # warn = 1 prints warnings immediately so they reach the log.
    .iprs_old_warn <- options(warn = 1)
    on.exit(options(.iprs_old_warn), add = TRUE)

    # Replay what was printed before the log existed.
    if (length(.iprs_preamble) > 0)
        cat(paste(.iprs_preamble, collapse = "\n"), "\n\n", sep = "")

    cat("\n", strrep("=", 72), "\n")
    cat("  Part10b IPRS -- IP detection on the full GOI-positive population\n")
    cat("  Arm          :", ARM_TAG, "\n")
    cat("  GOI          :", GOI, "\n")
    cat("  Composite    :", DFG_COMPOSITE_LABEL, "\n")
    cat("  Normalisation:", NORM_TAG, " log1p(UMIs/total x ",
        format(NORM_TARGET_SUM, scientific = FALSE), ")\n")
    cat("  Population   : all GOI>0 cells\n")
    cat("  GAM spec     : Expression ~ GOI + s(GOI,", GAM_BS, ") gamma =",
        GAM_GAMMA_IPRS, "\n")
    cat("  Estimator    :",
        if (GAM_DONOR_RE)
            sprintf("donor random intercept, s(donor, bs = \"re\"), minimum %d donors",
                    DONOR_RE_MIN_LEVELS)
        else "cell level (donor random intercept off)", "\n")
    cat("  Weights      :",
        if (W_DONOR_WEIGHTS)
            "inverse donor cell count, w = 1/n_i (each donor sums to 1)"
        else "uniform, w = 1 (inverse-donor weighting off)", "\n")
    cat("  R(x) donor RE:",
        if (RESID_DONOR_RE) "ON"
        else "OFF (residual GAM fitted at cell level)", "\n")
    cat("  IP method    :", "IPRS R(x)=0 zero-crossing + Brent CI", "\n")
    cat("  Datasets     :", paste(DATASETS, collapse = ", "), "\n")
    cat("  Part 10 dir  :", PART10_DIR, "\n")
    cat("  Output dir   :", OUT_DIR, "\n")
    cat("  Part 10b RDS :", P10B_IPRS_RDS, "\n")
    cat(strrep("=", 72), "\n\n")

    # Marker written on entry and removed on a clean finish.
    INCOMPLETE_MARKER_IPRS <- file.path(OUT_DIR, "_RUN_INCOMPLETE.txt")
    writeLines(c(
        "This arm did not finish.",
        "",
        paste("Arm     :", arm_tag),
        paste("Part    :", PART_TAG),
        paste("Started :", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
        "",
        "Any files in this folder are from a partial run and must not be used.",
        "This file is removed automatically when the arm completes."
    ), INCOMPLETE_MARKER_IPRS)


    # =========================================================================
    # SECTION A0b: DONOR RANDOM INTERCEPT, WEIGHTS, PREDICTION, DE
    # =========================================================================
    # Local copies of the Part 9 and Part 10 helpers.

    # Donor as a factor with unused levels dropped.
    as_donor_factor <- function(x) droplevels(factor(as.character(x)))

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

    # Below the minimum donor count the cell-level formula is used.
    donor_re_ok <- function(x)
        isTRUE(GAM_DONOR_RE) && n_donor_levels(x) >= DONOR_RE_MIN_LEVELS

    # w = 1 / (cells this donor contributes to this frame); each donor sums to 1.
    donor_weights <- function(donor_fac) {
        if (!isTRUE(W_DONOR_WEIGHTS)) return(rep(1, length(donor_fac)))
        f <- as_donor_factor(donor_fac)
        as.numeric(1 / as.numeric(table(f)[f]))
    }

    # Level the donor factor and recompute the weights from this frame's own
    # donors. Idempotent, so applying it twice gives the same column.
    prep_fit_frame <- function(data) {
        if ("donor" %in% names(data)) data$donor <- as_donor_factor(data$donor)
        data$w <- if ("donor" %in% names(data)) donor_weights(data$donor)
                  else rep(1, nrow(data))
        data
    }

    ensure_weights <- function(data) {
        if (!("w" %in% names(data))) {
            data$w <- if ("donor" %in% names(data))
                donor_weights(data$donor) else rep(1, nrow(data))
        }
        data
    }

    # Full-population GAM formula.
    make_gam_formula <- function(k, use_re) {
        rhs <- sprintf("GOI + s(GOI, bs = \"%s\", k = %d)", GAM_BS, as.integer(k))
        if (use_re) rhs <- paste0(rhs, " + s(donor, bs = \"re\")")
        stats::as.formula(paste("Expression ~", rhs))
    }

    # Residual GAM formula.
    make_resid_formula <- function(k, use_re) {
        rhs <- sprintf("s(GOI, bs = \"tp\", k = %d)", as.integer(k))
        if (use_re) rhs <- paste0(rhs, " + s(donor, bs = \"re\")")
        stats::as.formula(paste("residual ~", rhs))
    }

    model_has_re <- function(model)
        !is.null(model$smooth) &&
            any(vapply(model$smooth, function(s) s$label, character(1)) ==
                    DONOR_RE_LABEL)

    # Population-level prediction: donor supplied at an observed level, then
    # excluded. resp is the predictor column in newdata.
    predict_population <- function(model, x_seq, fit_data, se = TRUE) {
        nd <- data.frame(GOI = x_seq)
        if (model_has_re(model)) {
            lev <- if (!is.null(fit_data) && "donor" %in% names(fit_data))
                as_donor_factor(fit_data$donor)[1] else NULL
            if (is.null(lev)) {
                # The representative donor level comes from $var.summary.
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

    # Null deviance against a donor-aware, identically weighted null (cached).
    .iprs_null_dev_cache <- new.env(hash = TRUE, parent = emptyenv())
    get_null_deviance <- function(data) {
        if ("donor" %in% names(data)) data$donor <- as_donor_factor(data$donor)
        data   <- ensure_weights(data)
        use_re <- donor_re_ok(data)
        key <- paste0(nrow(data), "_", round(sum(data$Expression), 6),
                      "_", n_donor_levels(data), "_", as.integer(use_re),
                      "_", round(sum(data$w^2), 9))
        if (exists(key, envir = .iprs_null_dev_cache))
            return(get(key, envir = .iprs_null_dev_cache))
        set.seed(123)
        nd <- if (use_re)
            deviance(gam(Expression ~ 1 + s(donor, bs = "re"),
                         data = data, weights = w, method = GAM_METHOD))
        else
            deviance(gam(Expression ~ 1, data = data, weights = w))
        assign(key, nd, envir = .iprs_null_dev_cache)
        nd
    }

    # The intercept-only null model, used by the EP computation and the scatter.
    fit_null_model <- function(data) {
        data <- ensure_weights(data)
        set.seed(123)
        if (donor_re_ok(data))
            gam(Expression ~ 1 + s(donor, bs = "re"),
                data = data, weights = w, method = GAM_METHOD)
        else
            gam(Expression ~ 1, data = data, weights = w)
    }

    # Deviance explained against the donor-aware null.
    de_explained <- function(model, data) {
        if (is.null(model) || is.null(data)) return(NA_real_)
        if (!("w" %in% names(data)) && !is.null(model$prior.weights) &&
            length(model$prior.weights) == nrow(data))
            data$w <- as.numeric(model$prior.weights)
        nd <- tryCatch(get_null_deviance(data), error = function(e) NA_real_)
        if (!is.finite(nd) || nd <= 0) return(NA_real_)
        as.numeric((nd - deviance(model)) / nd)
    }

    # Flat-null DE, kept for the DE_vs_flat_null column.
    de_flat_null <- function(model) {
        if (is.null(model)) return(NA_real_)
        tryCatch(as.numeric(summary(model)$dev.expl),
                 error = function(e) NA_real_)
    }

    # Weights for the EP ratio and the purple count: frame w, model weights, or 1.
    ep_weight_vec <- function(model, data) {
        if ("w" %in% names(data)) return(as.numeric(data$w))
        if (!is.null(model) && !is.null(model$prior.weights) &&
            length(model$prior.weights) == nrow(data))
            return(as.numeric(model$prior.weights))
        rep(1, nrow(data))
    }

    # Weighted population-level deviance explained, 1 - RSS_w / TSS_w, from the
    # same prediction the EP ranking uses.
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

    # Population-level fitted response at given GOI values (donor term excluded).
    fitted_response_at <- function(model, x) {
        if (is.null(model) || !length(x)) return(rep(NA_real_, length(x)))
        ok <- is.finite(x)
        out <- rep(NA_real_, length(x))
        if (!any(ok)) return(out)
        nd <- data.frame(GOI = x[ok])
        dl <- tryCatch(model$var.summary[["donor"]], error = function(e) NULL)
        if (!is.null(dl)) {
            lv <- if (is.factor(dl)) levels(dl) else as.character(dl)
            nd$donor <- factor(lv[1], levels = lv)
        }
        out[ok] <- tryCatch(
            as.numeric(predict(model, newdata = nd,
                               exclude = DONOR_RE_LABEL, type = "response")),
            error = function(e) rep(NA_real_, sum(ok)))
        out
    }

    # Population-level prediction and its DE in one call.
    pop_dev_explained_model <- function(model, data) {
        if (is.null(model) || is.null(data)) return(NA_real_)
        mf <- tryCatch(as.numeric(predict(model, newdata = data,
                                          exclude = DONOR_RE_LABEL,
                                          type = "response")),
                       error = function(e) NULL)
        pop_dev_explained(model, data, mf)
    }

    # Fraction that sizes the purple set, clamped to [0, 1].
    purple_count_fraction <- function(de) {
        if (!is.finite(de)) return(NA_real_)
        max(0, min(1, de))
    }

    # Fits the full-population GAM.
    fit_full_pop_gam <- function(full_data, label = "full population",
                                 verbose = TRUE) {
        full_data <- prep_fit_frame(full_data)
        unique_x  <- length(unique(full_data$GOI))
        k_gp      <- min(10L, unique_x - 1L)
        n_don     <- n_donor_levels(full_data)
        use_re    <- donor_re_ok(full_data)

        if (unique_x < 4L) {
            if (verbose)
                cat("  ERROR [", label,
                    "]: fewer than 4 unique GOI values, cannot fit GAM\n",
                    sep = "")
            return(NULL)
        }
        if (verbose && isTRUE(GAM_DONOR_RE) && !use_re)
            cat(sprintf(paste0("  [%s] donor random intercept not fitted: %d ",
                               "donor level(s), minimum is %d; cell-level fit.\n"),
                        label, n_don, DONOR_RE_MIN_LEVELS))

        set.seed(123)
        m <- tryCatch(
            gam(make_gam_formula(k_gp, use_re),
                data    = full_data,
                weights = w,                      # bare column name
                method  = GAM_METHOD,
                select  = GAM_SELECT,
                gamma   = GAM_GAMMA_IPRS),
            error = function(e) {
                if (verbose)
                    cat("  ERROR [", label, "]: GAM fitting failed: ",
                        conditionMessage(e), "\n", sep = "")
                NULL
            })
        if (is.null(m)) return(NULL)

        list(model = m, data = full_data, k = k_gp,
             donor_re = use_re, n_donors = n_don)
    }


    # =========================================================================
    # SECTION A0: PRE-ANALYSIS GUARDRAILS
    # =========================================================================
    # Input checks before any IPRS fitting.

    GR_LOG_IPRS <- file.path(OUT_DIR, paste0(OUT_PREFIX, "_Guardrail_Log.txt"))

    .griprs_log   <- list()
    .griprs_pass  <- 0L
    .griprs_fail  <- 0L
    .griprs_warn  <- 0L

    griprs <- function(check_id, status, detail) {
        # A zero-length detail would vanish from the log, so substitute a note.
        if (length(detail) != 1L || is.na(detail))
            detail <- paste0("[detail unavailable: the message was built from a",
                             " missing value and came out zero-length]")
        tag <- switch(status, PASS = "[PASS]", FAIL = "[FAIL]", WARN = "[WARN]")
        cat(sprintf("  %s  %-16s  %s\n", tag, check_id, detail))
        .griprs_log[[length(.griprs_log) + 1L]] <<- list(
            id = check_id, status = status, detail = detail)
        if (status == "PASS") .griprs_pass <<- .griprs_pass + 1L
        if (status == "FAIL") .griprs_fail <<- .griprs_fail + 1L
        if (status == "WARN") .griprs_warn <<- .griprs_warn + 1L
    }

    griprs_write <- function() {
        con <- file(GR_LOG_IPRS, "wt")
        # Close the connection even if a write fails.
        on.exit(try(close(con), silent = TRUE), add = TRUE)
        writeLines(c(
            strrep("=", 72),
            "Part 10b IPRS -- Guardrail Log",
            paste("Arm :", ARM_TAG),
            paste("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
            paste0("Normalisation: ", NORM_TAG, " log1p(UMIs/total x ",
                   format(NORM_TARGET_SUM, scientific = FALSE), ")"),
            paste0("Estimator: ",
                   if (GAM_DONOR_RE) "donor random intercept" else "cell level",
                   " | weights: ",
                   if (W_DONOR_WEIGHTS) "inverse donor cell count" else "uniform"),
            strrep("=", 72), ""
        ), con)
        for (item in .griprs_log)
            writeLines(sprintf("  [%s]  %-16s  %s", item$status, item$id,
                               item$detail), con)
        writeLines(c("", strrep("=", 72),
                     sprintf("  Total : %d",
                             .griprs_pass + .griprs_fail + .griprs_warn),
                     sprintf("  PASS  : %d", .griprs_pass),
                     sprintf("  WARN  : %d", .griprs_warn),
                     sprintf("  FAIL  : %d", .griprs_fail),
                     if (.griprs_fail == 0 && .griprs_warn == 0)
                         "\n  All guardrails passed"
                     else if (.griprs_fail == 0)
                         "\n  Passed with warnings"
                     else "\n  Guardrail failures detected",
                     strrep("=", 72)), con)
        close(con)
        cat("  Guardrail log written:", GR_LOG_IPRS, "\n")
    }

    cat("\n", strrep("=", 72), "\n")
    cat("  SECTION A0 -- PRE-ANALYSIS GUARDRAILS\n")
    cat(strrep("=", 72), "\n\n")

    # == Normalisation checks per dataset ====================================
    for (.ds in list(list(tag = "NeoC",    dir = NEOC_DIR),
                     list(tag = "NonNeoC", dir = NONNC_DIR))) {
        cat("  -- Dataset:", .ds$tag, "--------------------------------\n")

        norm_mtx_path <- file.path(.ds$dir, "deg_support", "norm.mtx")
        raw_mtx_path  <- file.path(.ds$dir, "counts", "raw_selected",
                                   "matrix.mtx")

        if (file.exists(norm_mtx_path)) {
            X_n   <- readMM(norm_mtx_path)
            n_max <- max(X_n)

            griprs(paste0("GR-01-", .ds$tag),
                   if (n_max > 0 && n_max <= NORM_CEILING + 1e-6) "PASS"
                   else "FAIL",
                   sprintf("norm.mtx max = %.4f  (ceiling log1p(%s) = %.4f)",
                           n_max, format(NORM_TARGET_SUM, scientific = FALSE),
                           NORM_CEILING))

            # expm1 of a log1p(CP10K) row sums to NORM_TARGET_SUM.
            n_r  <- nrow(X_n)
            idx  <- seq_len(min(50L, n_r))
            rs   <- rowSums(expm1(as.matrix(X_n[idx, , drop = FALSE])))
            dev  <- max(abs(rs - NORM_TARGET_SUM))
            griprs(paste0("GR-02-", .ds$tag),
                   if (dev < 1e-6 * NORM_TARGET_SUM) "PASS" else "FAIL",
                   sprintf(paste0("expm1(norm) row sums: mean = %.6f, ",
                                  "max |dev| = %.2e, n = %d"),
                           mean(rs), dev, length(idx)))
            rm(X_n); gc(verbose = FALSE)
        } else {
            griprs(paste0("GR-01-", .ds$tag), "FAIL",
                   paste("norm.mtx not found:", norm_mtx_path))
            griprs(paste0("GR-02-", .ds$tag), "FAIL",
                   paste("norm.mtx not found, so the CP10K row-sum identity",
                         "could not be checked:", norm_mtx_path))
        }

        if (file.exists(raw_mtx_path)) {
            X_r    <- readMM(raw_mtx_path)
            r_samp <- as.numeric(X_r@x[seq_len(min(10000L, length(X_r@x)))])
            is_int <- all(abs(r_samp - round(r_samp)) < 0.01)
            r_max  <- max(X_r)
            griprs(paste0("GR-04-", .ds$tag),
                   if (is_int && r_max > 20) "PASS" else "FAIL",
                   sprintf("raw max = %.0f, integer-like = %s",
                           r_max, if (is_int) "TRUE" else "FALSE"))
            rm(X_r); gc(verbose = FALSE)
        } else {
            griprs(paste0("GR-04-", .ds$tag), "FAIL",
                   paste("raw_selected matrix not found:", raw_mtx_path))
        }
    }

    # == Part 10 RDS files exist and are loadable ============================
    p10_data <- list()
    for (ds in DATASETS) {
        rds_path <- p10_rds_path(ds)
        rds_ok   <- file.exists(rds_path)
        loadable <- FALSE
        if (rds_ok) {
            obj <- tryCatch(readRDS(rds_path), error = function(e) NULL)
            loadable <- !is.null(obj)
            if (loadable) p10_data[[ds]] <- obj
        }
        griprs(paste0("GR-05-", ds),
               if (rds_ok && loadable) "PASS" else "FAIL",
               if (rds_ok)
                   sprintf("Part 10 RDS: %s  loadable=%s",
                           file.path(basename(dirname(rds_path)),
                                     basename(rds_path)),
                           if (loadable) "TRUE" else "FALSE")
               else
                   paste("Part 10 RDS not found in either location:",
                         paste(p10_rds_candidates(ds), collapse = "  |  ")))
    }

    if (length(p10_data) == 0L) {
        griprs_write()
        stop("No Part 10 RDS could be loaded for arm ", arm_tag,
             ". Run Part 10 for this arm first.", call. = FALSE)
    }

    # == Donor lookup from the per-region gam_data CSVs ======================
    # Taken from the RDS if present, otherwise mapped from the CSVs by barcode.
    build_donor_map <- function(ds) {
        gam_dir <- DATASET_DIRS[[ds]]
        regs    <- DATASET_REGS[[ds]]
        stype   <- gsub("-", "_", gsub("/", "_", PRIMARY_STYPE))
        bc  <- character(0)
        don <- character(0)
        for (reg in regs) {
            gp <- file.path(gam_dir, "gam_data",
                            paste0(reg, "_", stype, "_selected.csv"))
            if (!file.exists(gp)) next
            df <- read.csv(gp, stringsAsFactors = FALSE, check.names = FALSE)
            if (!("donor" %in% colnames(df)))
                stop(sprintf(paste0("Region %s: gam_data has no donor column, ",
                                    "so the donor random intercept cannot be ",
                                    "fitted."), reg), call. = FALSE)
            if (!("barcode" %in% colnames(df)))
                stop(sprintf("Region %s: gam_data has no barcode column.", reg),
                     call. = FALSE)
            bc  <- c(bc,  df$barcode)
            don <- c(don, as.character(df$donor))
        }
        if (length(bc) == 0L)
            stop("No gam_data CSVs found for ", ds,
                 " under ", file.path(gam_dir, "gam_data"), call. = FALSE)
        dup <- duplicated(bc)
        if (any(dup)) {
            # A barcode must not map to two different donors.
            chk <- tapply(don, bc, function(v) length(unique(v)))
            if (any(chk > 1L))
                stop(sprintf(paste0("%s: %d barcode(s) map to more than one ",
                                    "donor across the region CSVs."),
                             ds, sum(chk > 1L)), call. = FALSE)
            bc  <- bc[!dup]
            don <- don[!dup]
        }
        setNames(don, bc)
    }


    # =========================================================================
    # SECTION A0c: FULL-POPULATION GAM EXTRACTION
    # =========================================================================
    # Fits a GAM on all GOI > 0 cells pooled from the Part 10 RDS.

    extract_pooled_gam <- function(ds) {
        obj <- p10_data[[ds]]
        if (is.null(obj)) return(NULL)

        score_df <- obj[["score_pooled_df"]]
        if (is.null(score_df)) {
            cat("  ERROR [", ds, "]: score_pooled_df not found in Part 10 RDS.\n",
                "  Top-level names: ", paste(names(obj), collapse = ", "), "\n",
                sep = "")
            return(NULL)
        }

        req_cols <- c("GOI_expr", "cDFG_expr", "kend_contrib", "spear_norm",
                      "is_skcs", "barcode")
        missing  <- setdiff(req_cols, names(score_df))
        if (length(missing) > 0) {
            cat("  ERROR [", ds, "]: score_pooled_df missing columns: ",
                paste(missing, collapse = ", "), "\n", sep = "")
            return(NULL)
        }

        # Keep GOI > 0 (Part 10 appends GOI-negative cells with NA scores).
        n_pre    <- nrow(score_df)
        # which() drops NA subscripts.
        score_df <- score_df[which(score_df$GOI_expr > 0), ]
        n_post   <- nrow(score_df)
        if (n_pre != n_post)
            cat(sprintf("  [%s] GOI>0 filter: %d -> %d (removed %d non-positive cells)\n",
                        ds, n_pre, n_post, n_pre - n_post))

        n_full <- nrow(score_df)
        cat("  [", ds, "] Full GOI-positive population: n = ", n_full,
            " cells\n", sep = "")

        # == Donor ===========================================================
        if ("donor" %in% names(score_df)) {
            don_vals <- as.character(score_df$donor)
            cat("  [", ds, "] donor taken from the Part 10 RDS\n", sep = "")
        } else {
            dmap     <- build_donor_map(ds)
            don_vals <- unname(dmap[score_df$barcode])
            n_miss   <- sum(is.na(don_vals))
            if (n_miss > 0)
                stop(sprintf(paste0("%s: %d of %d barcode(s) in ",
                                    "score_pooled_df have no donor in the ",
                                    "gam_data CSVs."),
                             ds, n_miss, n_full), call. = FALSE)
            cat("  [", ds, "] donor matched from the gam_data CSVs by barcode (",
                length(dmap), " barcodes in the map)\n", sep = "")
        }

        full_data <- prep_fit_frame(data.frame(
            GOI        = score_df$GOI_expr,
            Expression = score_df$cDFG_expr,
            barcode    = score_df$barcode,
            region     = score_df$region,
            donor      = don_vals,
            stringsAsFactors = FALSE))

        # A non-finite GOI or Expression would misalign the EP ranking: hard stop.
        .bad_goi <- !is.finite(full_data$GOI)
        .bad_exp <- !is.finite(full_data$Expression)
        if (any(.bad_goi) || any(.bad_exp))
            stop(sprintf(paste0("%s: %d of %d cell(s) have a non-finite GOI ",
                                "and %d a non-finite %s."),
                        ds, sum(.bad_goi), nrow(full_data), sum(.bad_exp),
                        DFG_COMPOSITE_LABEL),
                 call. = FALSE)

        cat(sprintf("  [%s] donors = %d | sum(w) = %.4f\n",
                    ds, n_donor_levels(full_data), sum(full_data$w)))

        # == Full-population GAM =============================================
        fit <- fit_full_pop_gam(full_data, label = ds)
        if (is.null(fit)) return(NULL)
        gam_full  <- fit$model
        full_data <- fit$data

        # == Donor-aware deviance explained ==================================
        de      <- de_explained(gam_full, full_data)
        de_flat <- de_flat_null(gam_full)
        # Population-level DE that sizes the purple set.
        de_pop  <- pop_dev_explained_model(gam_full, full_data)
        if (!is.finite(de_pop)) de_pop <- de_flat
        cat(sprintf(paste0("  [%s] Full-population GAM: n=%d | k=%d | ",
                           "donor RE=%s | DE=%.1f%% | DE_vs_flat_null=%.1f%% | ",
                           "gamma=%.1f\n"),
                    ds, n_full, fit$k, if (fit$donor_re) "YES" else "NO",
                    de * 100, de_flat * 100, GAM_GAMMA_IPRS))

        skcs_pooled <- obj[["skcs_pooled"]]
        if (!is.null(skcs_pooled))
            cat(sprintf("  [%s] SKCS reference: rho=%.4f | tau=%.4f | SKCS-core=%d/%d\n",
                        ds, skcs_pooled$global_rho, skcs_pooled$global_tau,
                        skcs_pooled$n_skcs,
                        skcs_pooled$n_skcs + skcs_pooled$n_nskcs))

        list(gam_model    = gam_full,
             data         = full_data,
             dev_expl     = de,
             dev_expl_flat = de_flat,
             # The quantity that sizes the purple set.
             dev_expl_pop = de_pop,
             donor_re     = fit$donor_re,
             n_donors     = fit$n_donors,
             k_used       = fit$k,
             key_used     = "Pooled|DFG|FullPop",
             score_df     = score_df,
             skcs_info    = skcs_pooled)
    }

    gam_inputs <- list()
    for (ds in DATASETS) {
        gam_inputs[[ds]] <- tryCatch(extract_pooled_gam(ds),
                                     error = function(e) {
                                         cat("  ERROR [", ds, "]: ",
                                             conditionMessage(e), "\n", sep = "")
                                         NULL
                                     })
        if (!is.null(gam_inputs[[ds]]))
            cat(sprintf("  [%s] Pooled GAM: n=%d cells | DE=%.1f%% | key=%s\n",
                        ds, nrow(gam_inputs[[ds]]$data),
                        gam_inputs[[ds]]$dev_expl * 100,
                        gam_inputs[[ds]]$key_used))
    }
    cat("\n")

    # == Pooled full GOI-positive cell counts ================================
    for (ds in DATASETS) {
        gi <- gam_inputs[[ds]]
        if (is.null(gi)) {
            griprs(paste0("GR-06-", ds), "FAIL",
                   "No valid pooled GAM input extracted")
        } else {
            n_cells <- nrow(gi$data)
            griprs(paste0("GR-06-", ds),
                   if (n_cells >= MIN_PURPLE_CELLS) "PASS" else "FAIL",
                   sprintf("Pooled full GOI-positive cells = %d  (min = %d)",
                           n_cells, MIN_PURPLE_CELLS))
        }
    }

    # == Full-population pooled GAM DE > 0 ===================================
    for (ds in DATASETS) {
        gi <- gam_inputs[[ds]]
        if (is.null(gi)) {
            griprs(paste0("GR-07-", ds), "FAIL",
                   "No valid pooled GAM, cannot check DE")
        } else {
            de_pct <- gi$dev_expl * 100
            griprs(paste0("GR-07-", ds),
                   if (!is.na(de_pct) && de_pct > 0) "PASS" else "FAIL",
                   sprintf(paste0("Full-population pooled GAM DE = %.1f%% ",
                                  "(donor-aware null) | vs flat null = %.1f%%"),
                           de_pct, gi$dev_expl_flat * 100))
        }
    }

    # == Donor RE on the full-population GAM =================================
    for (ds in DATASETS) {
        gi <- gam_inputs[[ds]]
        if (is.null(gi)) {
            griprs(paste0("GR-11-", ds), "WARN",
                   "no GAM fitted, so donor-RE coverage cannot be assessed")
        } else if (!isTRUE(GAM_DONOR_RE)) {
            griprs(paste0("GR-11-", ds), "WARN",
                   "GAM_DONOR_RE is FALSE: the full-population fit is cell level")
        } else {
            griprs(paste0("GR-11-", ds),
                   if (isTRUE(gi$donor_re)) "PASS" else "WARN",
                   sprintf(paste0("full-population GAM donor RE = %s | ",
                                  "donors = %d | sum(w) = %.4f"),
                           if (isTRUE(gi$donor_re)) "YES" else
                               sprintf("NO (below %d)", DONOR_RE_MIN_LEVELS),
                           gi$n_donors, sum(gi$data$w)))
        }
    }

    # == Write pre-analysis guardrail log ====================================
    # Written here and again at the end.
    griprs_write()
    cat(sprintf("\n  Pre-analysis guardrails: %d PASS / %d WARN / %d FAIL\n",
                .griprs_pass, .griprs_warn, .griprs_fail))
    if (.griprs_fail > 0) {
        cat("  Guardrail failures detected\n\n")
        warning("Part 10b IPRS pre-analysis guardrails: ", .griprs_fail,
                " FAIL(s)", call. = FALSE)
    } else {
        cat("  All pre-analysis guardrails passed\n\n")
    }


    # =========================================================================
    # SECTION 1: Gene set labels and members
    # =========================================================================
    # The composite label is arm-aware.
    GS_KEY    <- DFG_COMPOSITE_LABEL
    GS_LABELS <- list(
        cDFG = paste0(DFG_COMPOSITE_LABEL, " composite (full GOI+, Pooled)")
    )
    # Safe default; overwritten per dataset in Section 4.
    GS_MEMBERS <- list(cDFG = "unknown")


    # =========================================================================
    # SECTION 2: Tuning parameters
    # =========================================================================
    # Thresholds are set in the user parameter block.


    # =========================================================================
    # SECTION 3: Helper functions
    # =========================================================================

    # == 3a. SW-gated two-group stat test ====================================
    stat_test_zc <- function(pre_r, post_r) {
        out <- list(stat_test   = NA_character_,
                    stat_p_raw  = NA_real_,
                    stat_p_adj  = NA_real_,
                    both_normal = NA)

        if (length(pre_r) < 3L || length(post_r) < POST_IP_MIN_CELLS)
            return(out)

        # shapiro.test() accepts at most 5,000 observations; Wilcoxon above that.
        SHAPIRO_MAX <- 5000L
        if (length(pre_r) > SHAPIRO_MAX || length(post_r) > SHAPIRO_MAX) {
            sw_pre <- NA_real_; sw_post <- NA_real_
            both_normal <- FALSE
        } else {
            sw_pre  <- if (length(pre_r)  >= 3L)
                suppressWarnings(shapiro.test(pre_r))$p.value  else NA_real_
            sw_post <- if (length(post_r) >= 3L)
                suppressWarnings(shapiro.test(post_r))$p.value else NA_real_
            both_normal <- !is.na(sw_pre) && !is.na(sw_post) &&
                sw_pre > SHAPIRO_ALPHA && sw_post > SHAPIRO_ALPHA
        }
        out$both_normal <- both_normal

        if (isTRUE(both_normal)) {
            tt <- suppressWarnings(
                t.test(post_r, pre_r, alternative = "greater", var.equal = FALSE))
            out$stat_test  <- "Welch_t"
            out$stat_p_raw <- tt$p.value
        } else {
            mw <- suppressWarnings(
                wilcox.test(post_r, pre_r, alternative = "greater", exact = FALSE))
            out$stat_test  <- "MannWhitneyU"
            out$stat_p_raw <- mw$p.value
        }
        # One test per dataset, so stat_p_adj equals stat_p_raw.
        out$stat_p_adj <- out$stat_p_raw
        out
    }

    # == 3b. Brent CI on R(x) GAM confidence bands ===========================
    find_ci_brent <- function(residual_gam, x_seq, r_hat, r_se,
                              multiplier, near_x, fit_data = NULL) {
        band     <- r_hat + multiplier * r_se
        sign_chg <- which(diff(sign(band)) != 0L)
        if (length(sign_chg) == 0L) return(NA_real_)

        if (length(sign_chg) > 1L)
            sign_chg <- sign_chg[which.min(abs(x_seq[sign_chg] - near_x))]

        bracket <- c(x_seq[sign_chg], x_seq[sign_chg + 1L])

        tryCatch({
            uniroot(
                # newdata must carry every model term, including donor.
                f        = function(x) {
                    p <- predict_population(residual_gam, x, fit_data, se = TRUE)
                    as.numeric(p$fit) + multiplier * as.numeric(p$se.fit)
                },
                interval = bracket,
                tol      = .Machine$double.eps^0.5
            )$root
        }, error = function(e) {
            x1 <- bracket[1L]; x2 <- bracket[2L]
            y1 <- band[sign_chg]; y2 <- band[sign_chg + 1L]
            x1 + (0 - y1) * (x2 - x1) / (y2 - y1)
        })
    }

    # == 3b2. Proximity fallback CI bound ====================================
    find_ci_proximity_fallback <- function(x_seq, r_hat, r_se, multiplier,
                                           ip_value, side = c("lower", "upper")) {
        side <- match.arg(side)
        band <- r_hat + multiplier * r_se

        idx <- if (side == "lower") which(x_seq <= ip_value)
               else                 which(x_seq >= ip_value)
        if (length(idx) < 3L) return(NA_real_)

        band_region <- band[idx]
        if (any(diff(sign(band_region)) != 0L)) return(NA_real_)

        x_seq[idx[which.min(abs(band_region))]]
    }

    # == 3c. IPRS component calculation ======================================
    compute_iprs <- function(ip_value, ip_ci_lower, ip_ci_upper,
                             stat_p_adj, x_seq, r_hat, r_upper, r_lower,
                             residual_data) {

        # == C1: Precision ===================================================
        # C1 requires an IP; without one it is 0.
        goi_range <- max(x_seq) - min(x_seq)
        have_ip   <- !is.na(ip_value)
        both_ci   <- have_ip && !is.na(ip_ci_lower) && !is.na(ip_ci_upper)
        lower_ok  <- have_ip && !is.na(ip_ci_lower)
        upper_ok  <- have_ip && !is.na(ip_ci_upper)
        # ci_width depends only on the bounds; it is exported as CI_Width.
        ci_width  <- if (!is.na(ip_ci_lower) && !is.na(ip_ci_upper))
            abs(ip_ci_upper - ip_ci_lower) else NA_real_

        C1 <- if (both_ci && goi_range > 0) {
            max(0, 1 - ci_width / goi_range)
        } else if ((lower_ok || upper_ok) && goi_range > 0) {
            half_width <- if (lower_ok) abs(ip_value - ip_ci_lower)
                          else          abs(ip_ci_upper - ip_value)
            0.5 * max(0, 1 - half_width / goi_range)
        } else {
            0
        }

        C1_note <- if (!have_ip)
            sprintf(paste0("No IP: R(x) has no neg-to-pos crossing, so there ",
                           "is no precision to grade -> C1=0. Band-crossing ",
                           "width was %s and is reported as CI_Width only."),
                    if (is.na(ci_width)) "NA" else sprintf("%.4f", ci_width))
        else if (both_ci)
            sprintf("Both CI bounds: CI_width=%.4f / GOI_range=%.4f → C1=%.4f",
                    ci_width, goi_range, C1)
        else if (lower_ok || upper_ok) {
            known_bound <- if (lower_ok) ip_ci_lower else ip_ci_upper
            bound_label <- if (lower_ok) "lower" else "upper"
            half_w      <- abs(ip_value - known_bound)
            sprintf("Partial CI (%s only): half_width=%.4f / GOI_range=%.4f → 0.5×max(0,1-ratio)=%.4f",
                    bound_label, half_w, goi_range, C1)
        } else
            "Both CI bounds NA → C1=0"

        # == C2: Statistical =================================================
        GW_ALPHA <- 0.05
        scale_c2 <- atanh(0.5) / (-log10(GW_ALPHA))
        C2 <- if (!is.na(stat_p_adj))
            tanh(scale_c2 * (-log10(max(stat_p_adj, 1e-300))))
        else
            0

        C2_note <- if (is.na(stat_p_adj))
            "stat_p_adj not available → C2=0"
        else
            sprintf("tanh(%.4f * -log10(%.4e)) → C2=%.4f",
                    scale_c2, stat_p_adj, C2)

        # == C3: Abruptness ==================================================
        deriv_at_ip <- if (!is.na(ip_value) && length(x_seq) >= 3L) {
            i_ip <- which.min(abs(x_seq - ip_value))
            n    <- length(x_seq)
            if (i_ip > 1L && i_ip < n)
                abs((r_hat[i_ip + 1L] - r_hat[i_ip - 1L]) /
                        (x_seq[i_ip + 1L] - x_seq[i_ip - 1L]))
            else if (i_ip == 1L)
                abs((r_hat[2L] - r_hat[1L]) / (x_seq[2L] - x_seq[1L]))
            else
                abs((r_hat[n] - r_hat[n - 1L]) / (x_seq[n] - x_seq[n - 1L]))
        } else NA_real_

        C3 <- if (!is.na(deriv_at_ip)) tanh(deriv_at_ip) else 0

        C3_note <- if (is.na(deriv_at_ip))
            "IP not available → C3=0"
        else
            sprintf("|R'(IP)|=%.4f → tanh(%.4f) → C3=%.4f",
                    deriv_at_ip, deriv_at_ip, C3)

        # == C4: GAM ribbon width penalty ====================================
        amplitude   <- diff(range(r_hat, na.rm = TRUE))
        mean_ribbon <- mean(r_upper - r_lower, na.rm = TRUE)
        C4          <- exp(-mean_ribbon / (amplitude + 1e-6))

        C4_note <- sprintf(
            "mean_ribbon=%.4f / amplitude=%.4f → exp(-ratio=%.4f) → C4=%.4f",
            mean_ribbon, amplitude, mean_ribbon / (amplitude + 1e-6), C4)

        # == C5: R(x) smoothness + post-IP sign fidelity =====================
        tv <- sum(abs(diff(r_hat)), na.rm = TRUE)

        if (!is.na(ip_value) && nrow(residual_data) > 0L) {
            post_cells    <- residual_data$residual[residual_data$GOI >= ip_value]
            n_post        <- length(post_cells)
            neg_post_frac <- if (n_post > 0L) sum(post_cells < 0) / n_post else 0
        } else {
            neg_post_frac <- 0
        }
        sign_fidelity <- (1 - neg_post_frac)^2

        if (amplitude < 1e-6) {
            C5      <- 0
            C5_note <- sprintf("amplitude=%.6f < 1e-6 (flat R(x)) → C5=0", amplitude)
        } else {
            excess  <- tv / (amplitude + 1e-6) - 1
            wave_sc <- exp(-max(0, excess))
            C5      <- wave_sc * sign_fidelity
            C5_note <- sprintf(
                "TV=%.4f/amp=%.4f → wave=%.4f | neg_post=%.3f → sign=(1-%.3f)^2=%.4f | C5=%.4f*%.4f=%.4f",
                tv, amplitude, wave_sc, neg_post_frac, neg_post_frac,
                sign_fidelity, wave_sc, sign_fidelity, C5)
        }

        IPRS     <- (C1 + C2 + C3 + C4 + C5) / 5
        any_zero <- any(c(C1, C2, C3, C4, C5) == 0)

        IPRS_tier <- if (IPRS >= IPRS_VERY_STRONG) "Very Strong" else
            if (IPRS >= IPRS_STRONG)      "Strong"      else
                if (IPRS >= IPRS_MODERATE)    "Moderate"    else
                    if (IPRS >= IPRS_WEAK)        "Weak"        else
                        "Failed"

        IPRS_flag <- if (any_zero) "Visual review needed [C=0]" else
            if (IPRS_tier %in% c("Very Strong", "Strong")) "None" else
                if (IPRS_tier %in% c("Moderate", "Weak"))  "Visual review needed" else
                    "No IP - do not proceed"

        # No IP: the tier is Failed; the components are kept as computed.
        if (!have_ip) {
            IPRS_tier <- "Failed"
            IPRS_flag <- "No IP - do not proceed"
        }

        cat(sprintf("  IPRS components: C1=%.4f [%s]\n", C1, C1_note))
        cat(sprintf("                   C2=%.4f [%s]\n", C2, C2_note))
        cat(sprintf("                   C3=%.4f [%s]\n", C3, C3_note))
        cat(sprintf("                   C4=%.4f [%s]\n", C4, C4_note))
        cat(sprintf("                   C5=%.4f [%s]\n", C5, C5_note))
        cat(sprintf("  IPRS = (%.4f+%.4f+%.4f+%.4f+%.4f)/5 = %.4f  [%s]  %s\n",
                    C1, C2, C3, C4, C5, IPRS, IPRS_tier,
                    if (IPRS_flag == "None") "" else paste0("*** ", IPRS_flag, " ***")))

        list(C1 = C1, C2 = C2, C3 = C3, C4 = C4, C5 = C5,
             IPRS = IPRS, IPRS_tier = IPRS_tier, IPRS_flag = IPRS_flag,
             C1_note = C1_note, C2_note = C2_note, C3_note = C3_note,
             C4_note = C4_note, C5_note = C5_note,
             deriv_at_ip = deriv_at_ip,
             ci_width    = ci_width,
             goi_range   = goi_range,
             amplitude   = amplitude,
             mean_ribbon = mean_ribbon,
             tv          = tv)
    }

    # == 3d. Core IP detection function ======================================
    # sample_data must have columns: GOI, Expression
    detect_ip <- function(sample_data, gam_model, dev_explained_pct,
                          dataset_tag) {

        # Level the donor factor and recompute weights for this frame.
        sample_data <- prep_fit_frame(sample_data)

        # Step 1: Purple-cell classification (top EP cells)
        # Flat intercept null, ratio of squared deviations, top round(n * DE) cells.
        null_fitted  <- fitted(gam(Expression ~ 1, data = sample_data))
        model_fitted <- as.numeric(predict(gam_model, newdata = sample_data,
                                           exclude = DONOR_RE_LABEL,
                                           type = "response"))
        # Both vectors must have one value per row.
        if (length(null_fitted) != nrow(sample_data) ||
            length(model_fitted) != nrow(sample_data))
            stop(sprintf(paste0("%s: the EP vectors are the wrong length ",
                                "(null %d, model %d, rows %d)."),
                        dataset_tag, length(null_fitted), length(model_fitted),
                        nrow(sample_data)),
                 call. = FALSE)
        null_sq      <- (sample_data$Expression - null_fitted)^2
        model_sq     <- (sample_data$Expression - model_fitted)^2
        expl_power   <- ifelse(null_sq > 0, 1 - model_sq / null_sq, 0)

        # DE fraction clamped to [0, 1] before it multiplies n.
        n_purple  <- max(1L, min(
            round(nrow(sample_data) * purple_count_fraction(dev_explained_pct)),
            nrow(sample_data)))
        is_purple <- rep(FALSE, nrow(sample_data))
        is_purple[order(expl_power, decreasing = TRUE)[seq_len(n_purple)]] <- TRUE

        # Step 2: Signed residuals for purple cells
        # Against the population-level curve.
        residual      <- sample_data$Expression -
            as.numeric(predict_population(gam_model, sample_data$GOI,
                                          sample_data, se = FALSE))
        purple_idx    <- which(is_purple)
        # Donor travels into the residual frame; weights recomputed for the subset.
        residual_data <- prep_fit_frame(data.frame(
            GOI      = sample_data$GOI[purple_idx],
            residual = residual[purple_idx],
            donor    = sample_data$donor[purple_idx],
            stringsAsFactors = FALSE
        ))

        cat("  Purple cells:", nrow(residual_data),
            "| Above:", sum(residual_data$residual > 0),
            "| Below:", sum(residual_data$residual < 0),
            "| donors:", n_donor_levels(residual_data),
            "| sum(w):", round(sum(residual_data$w), 4), "\n")

        # Per-donor purple composition (diagnostic only).
        if ("donor" %in% names(sample_data)) {
            .d_tot  <- table(droplevels(as_donor_factor(sample_data$donor)))
            .d_pur  <- table(droplevels(as_donor_factor(
                sample_data$donor[purple_idx])))
            .d_frac <- as.numeric(.d_pur[names(.d_tot)]) / as.numeric(.d_tot)
            .d_frac[is.na(.d_frac)] <- 0
            cat(sprintf(paste0("  Purple donor coverage: %d/%d donors ",
                               "represented | per-donor purple fraction ",
                               "[%.3f, %.3f] median %.3f\n"),
                        sum(.d_frac > 0), length(.d_tot),
                        min(.d_frac), max(.d_frac), stats::median(.d_frac)))
        }

        FAILED_list <- function(status) {
            iprs_fail <- list(C1 = 0, C2 = 0, C3 = 0, C4 = 0, C5 = 0,
                              IPRS = 0, IPRS_tier = "Failed",
                              IPRS_flag = "No IP - do not proceed",
                              C1_note = "IP detection failed",
                              C2_note = "IP detection failed",
                              C3_note = "IP detection failed",
                              C4_note = "IP detection failed",
                              C5_note = "IP detection failed",
                              deriv_at_ip = NA_real_, ci_width  = NA_real_,
                              goi_range   = NA_real_, amplitude = NA_real_,
                              mean_ribbon = NA_real_, tv        = NA_real_)
            list(ip_value = NA_real_, manual_flag = TRUE, status = status,
                 resid_donor_re = NA, resid_n_donors = n_donor_levels(residual_data),
                 stat_test = NA_character_, stat_p_raw = NA_real_,
                 stat_p_adj = NA_real_, both_normal = NA,
                 ip_ci_lower = NA_real_, ip_ci_upper = NA_real_,
                 ci_lower_method = NA_character_, ci_upper_method = NA_character_,
                 neg_pos_count = 0L,
                 ip_y_value = NA_real_, ip_y_ci_lower = NA_real_,
                 ip_y_ci_upper = NA_real_,
                 iprs = iprs_fail,
                 residual_data = residual_data,
                 residual_pred = data.frame(GOI     = numeric(0),
                                            r_hat   = numeric(0),
                                            se      = numeric(0),
                                            r_upper = numeric(0),
                                            r_lower = numeric(0)))
        }

        if (nrow(residual_data) < MIN_PURPLE_CELLS ||
            length(unique(residual_data$GOI)) < 3L) {
            cat("  Insufficient purple cells, skipping IP detection\n")
            return(FAILED_list("insufficient_purple_cells"))
        }

        # Step 3: Residual GAM R(x)
        # Donor random intercept included; the k=8 fallback keeps the same formula.
        set.seed(123)
        n_unique     <- length(unique(residual_data$GOI))
        k_resid      <- min(20L, n_unique - 1L)
        resid_re     <- isTRUE(RESID_DONOR_RE) && donor_re_ok(residual_data)
        if (isTRUE(GAM_DONOR_RE) && !isTRUE(RESID_DONOR_RE))
            cat("  R(x) GAM: donor random intercept off for the residual",
                "fit (RESID_DONOR_RE = FALSE)\n")
        if (isTRUE(GAM_DONOR_RE) && isTRUE(RESID_DONOR_RE) && !resid_re)
            cat(sprintf(paste0("  R(x) GAM: donor random intercept not fitted: %d ",
                               "donor level(s), minimum is %d\n"),
                        n_donor_levels(residual_data), DONOR_RE_MIN_LEVELS))
        residual_gam <- tryCatch(
            gam(make_resid_formula(k_resid, resid_re),
                data = residual_data, weights = w, method = "REML"),
            error = function(e) {
                cat("  Fallback k=8:", conditionMessage(e), "\n")
                gam(make_resid_formula(min(8L, n_unique - 1L), resid_re),
                    data = residual_data, weights = w, method = "REML")
            }
        )
        # The s(GOI) row of s.table, matched by name.
        gam_sum  <- summary(residual_gam)
        .rx_row  <- match("s(GOI)", rownames(gam_sum$s.table))
        if (is.na(.rx_row)) .rx_row <- 1L
        cat("  R(x) GAM: donor RE =", if (resid_re) "YES" else "NO",
            "| edf =", sprintf("%.2f", gam_sum$s.table[.rx_row, "edf"]),
            "| Dev =", sprintf("%.1f%%", gam_sum$dev.expl * 100),
            "| p =", sprintf("%.2e", gam_sum$s.table[.rx_row, "p-value"]), "\n")

        # Step 4: 2000-pt grid prediction
        x_seq <- seq(min(residual_data$GOI), max(residual_data$GOI),
                     length.out = 2000L)
        # Population-level R(x).
        pred  <- predict_population(residual_gam, x_seq, residual_data, se = TRUE)
        r_hat   <- as.numeric(pred$fit)
        r_se    <- as.numeric(pred$se.fit)
        r_upper <- r_hat + 1.96 * r_se
        r_lower <- r_hat - 1.96 * r_se

        # Step 5: First neg->pos zero-crossing
        sc_all <- which(diff(sign(r_hat)) != 0L)
        sc_np  <- sc_all[r_hat[sc_all] < 0 & r_hat[sc_all + 1L] > 0]
        ip_value    <- NA_real_
        manual_flag <- FALSE

        if (length(sc_np) > 0L) {
            i_cross  <- sc_np[1L]
            x0 <- x_seq[i_cross];      y0 <- r_hat[i_cross]
            x1 <- x_seq[i_cross + 1L]; y1 <- r_hat[i_cross + 1L]
            ip_value <- x0 - y0 * (x1 - x0) / (y1 - y0)
            cat("  IP (R(x) zero-crossing):", sprintf("%.6f", ip_value),
                "| neg\u2192pos crossings:", length(sc_np), "\n")
        } else {
            manual_flag <- TRUE
            cat("  R(x) has no neg\u2192pos zero-crossing; R(x) is entirely",
                ifelse(all(r_hat >= 0), "positive", "negative"), "\n")
        }

        # Step 6: Brent 95% CI with proximity fallback
        ci_anchor <- if (!is.na(ip_value)) ip_value else median(x_seq)

        ip_ci_lower <- tryCatch(
            find_ci_brent(residual_gam, x_seq, r_hat, r_se, +1.96, ci_anchor,
                          fit_data = residual_data),
            error = function(e) NA_real_)
        ci_lower_method <- if (!is.na(ip_ci_lower)) "brent" else NA_character_

        ip_ci_upper <- tryCatch(
            find_ci_brent(residual_gam, x_seq, r_hat, r_se, -1.96, ci_anchor,
                          fit_data = residual_data),
            error = function(e) NA_real_)
        ci_upper_method <- if (!is.na(ip_ci_upper)) "brent" else NA_character_

        if (is.na(ip_ci_lower) && !is.na(ip_value)) {
            fb <- find_ci_proximity_fallback(x_seq, r_hat, r_se,
                                             +1.96, ip_value, "lower")
            if (!is.na(fb)) {
                ip_ci_lower     <- fb
                ci_lower_method <- "proximity_fallback"
                cat("  CI lower: Brent NA → proximity fallback =",
                    sprintf("%.4f", ip_ci_lower), "\n")
            }
        }
        if (is.na(ip_ci_upper) && !is.na(ip_value)) {
            fb <- find_ci_proximity_fallback(x_seq, r_hat, r_se,
                                             -1.96, ip_value, "upper")
            if (!is.na(fb)) {
                ip_ci_upper     <- fb
                ci_upper_method <- "proximity_fallback"
                cat("  CI upper: Brent NA → proximity fallback =",
                    sprintf("%.4f", ip_ci_upper), "\n")
            }
        }

        # Swap if inverted
        if (!is.na(ip_ci_lower) && !is.na(ip_ci_upper) &&
            ip_ci_lower > ip_ci_upper) {
            tmp <- ip_ci_lower; ip_ci_lower <- ip_ci_upper; ip_ci_upper <- tmp
            tmp <- ci_lower_method; ci_lower_method <- ci_upper_method
            ci_upper_method <- tmp
        }

        cat("  95% CI: [",
            if (is.na(ip_ci_lower)) "NA"
            else sprintf("%.4f [%s]", ip_ci_lower, ci_lower_method), ",",
            if (is.na(ip_ci_upper)) "NA"
            else sprintf("%.4f [%s]", ip_ci_upper, ci_upper_method),
            "] | width:",
            if (is.na(ip_ci_lower) || is.na(ip_ci_upper)) "NA"
            else sprintf("%.4f", abs(ip_ci_upper - ip_ci_lower)), "\n")

        # Step 7: Stat test pre vs post
        if (!is.na(ip_value)) {
            pre_r  <- residual_data$residual[residual_data$GOI <  ip_value]
            post_r <- residual_data$residual[residual_data$GOI >= ip_value]
            st     <- stat_test_zc(pre_r, post_r)
            cat("  Stat:", st$stat_test,
                sprintf("| p (unadjusted): %s",
                        if (is.na(st$stat_p_adj)) "NA"
                        else sprintf("%.3e", st$stat_p_adj)), "\n")
            # Stat test not run when either side is too short; C2 is then 0.
            if (is.na(st$stat_p_adj))
                cat(sprintf(paste0("  Stat test not run: pre=%d (need >= 3), ",
                                   "post=%d (need >= %d). C2 will be 0.\n"),
                            length(pre_r), length(post_r), POST_IP_MIN_CELLS))
        } else {
            st <- list(stat_test = NA_character_, stat_p_raw = NA_real_,
                       stat_p_adj = NA_real_, both_normal = NA)
        }

        # Step 8: IPRS
        iprs <- compute_iprs(ip_value, ip_ci_lower, ip_ci_upper,
                             st$stat_p_adj, x_seq, r_hat, r_upper, r_lower,
                             residual_data)

        # Step 9: response-axis threshold IP_y = f(IP_x) and its mapped CI.
        ip_y_value <- unname(fitted_response_at(gam_model, ip_value))
        ip_y_bnds  <- unname(fitted_response_at(gam_model,
                                                c(ip_ci_lower, ip_ci_upper)))
        if (all(is.finite(ip_y_bnds))) ip_y_bnds <- sort(ip_y_bnds)
        cat(sprintf(paste0("  IP_y = f(IP_x) = %s | mapped CI [%s, %s]",
                           "  (classification threshold)\n"),
                    if (is.finite(ip_y_value)) sprintf("%.6f", ip_y_value) else "NA",
                    if (is.finite(ip_y_bnds[1])) sprintf("%.4f", ip_y_bnds[1]) else "NA",
                    if (is.finite(ip_y_bnds[2])) sprintf("%.4f", ip_y_bnds[2]) else "NA"))

        list(ip_value        = ip_value,
             ip_y_value      = ip_y_value,
             ip_y_ci_lower   = ip_y_bnds[1],
             ip_y_ci_upper   = ip_y_bnds[2],
             manual_flag     = manual_flag,
             status          = "ok",
             stat_test       = st$stat_test,
             stat_p_raw      = st$stat_p_raw,
             stat_p_adj      = st$stat_p_adj,
             both_normal     = st$both_normal,
             ip_ci_lower     = ip_ci_lower,
             ip_ci_upper     = ip_ci_upper,
             ci_lower_method = ci_lower_method,
             ci_upper_method = ci_upper_method,
             neg_pos_count   = length(sc_np),
             resid_donor_re  = resid_re,
             resid_n_donors  = n_donor_levels(residual_data),
             iprs            = iprs,
             residual_data   = residual_data,
             residual_pred   = data.frame(GOI     = x_seq,
                                          r_hat   = r_hat,
                                          se      = r_se,
                                          r_upper = r_upper,
                                          r_lower = r_lower))
    }



    # =========================================================================
    # SECTION A0d: SELF-TESTS (bodies in Section Z)
    # =========================================================================
    # Off by default; neither test touches the file system.
    if (isTRUE(RUN_SMOKE_TEST))  iprs_smoke_test(environment())
    if (isTRUE(RUN_STRESS_TEST)) iprs_stress_test(environment())


    # =========================================================================
    # SECTION 4: Run IP detection across datasets
    # =========================================================================
    set.seed(123)
    cat("\n", strrep("=", 72), "\n")
    cat(sprintf("IP DETECTION | Zero-crossing + Brent CI + IPRS | STAT_ALPHA=%.2f\n",
                STAT_ALPHA))
    cat("Input     : Full", GOI, "-positive population pooled",
        DFG_COMPOSITE_LABEL, "GAM\n")
    cat("Datasets  :", paste(DATASETS, collapse = ", "), "\n")
    cat(strrep("=", 72), "\n\n")

    ip_results <- list()

    for (ds in DATASETS) {
        gi <- gam_inputs[[ds]]
        if (is.null(gi)) {
            cat("  Skipping", ds, ": no valid GAM input\n\n")
            next
        }

        gs_genes <- if (ds == "NonNeoC") DFG_NONNC else DFG_NEOC
        GS_MEMBERS[["cDFG"]] <- paste(gs_genes, collapse = ", ")

        key <- paste0(ds, "|", GS_KEY)
        cat("\n--- Processing:", ds, "| cDFG composite ---\n")
        cat("  n cells (full GOI-positive population):", nrow(gi$data), "\n")
        # Population-level DE sizes the purple set; the other two are for reference.
        cat("  Dev explained (population level, sets n_purple):",
            sprintf("%.1f%%", gi$dev_expl_pop * 100),
            sprintf("| summary()$dev.expl: %.1f%% | donor-aware DE: %.1f%%\n",
                    gi$dev_expl_flat * 100, gi$dev_expl * 100))

        ip_results[[key]] <- tryCatch(
            detect_ip(gi$data, gi$gam_model, gi$dev_expl_pop, dataset_tag = ds),
            error = function(e) {
                cat("  ERROR in detect_ip:", conditionMessage(e), "\n")
                NULL
            }
        )
    }


    # =========================================================================
    # SECTION 5: Console summary
    # =========================================================================
    cat("\n\n", strrep("=", 72), "\n")
    cat("IP SUMMARY - IP Reliability Score (IPRS)\n")
    cat(strrep("=", 72), "\n\n")
    cat(sprintf("%-28s %8s %8s %8s %6s %6s %6s %6s %6s %6s  %-14s  %s\n",
                "Fit (Dataset|GS)", "IP", "CI_Lo", "CI_Hi",
                "C1", "C2", "C3", "C4", "C5", "IPRS", "Tier", "Flag"))
    cat(strrep("-", 140), "\n")
    for (key in names(ip_results)) {
        r <- ip_results[[key]]; if (is.null(r)) next
        iprs <- r$iprs
        cat(sprintf("%-28s %8s %8s %8s %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f  %-14s  %s\n",
                    key,
                    if (is.na(r$ip_value))    "NA" else sprintf("%.4f", r$ip_value),
                    if (is.na(r$ip_ci_lower)) "NA" else sprintf("%.4f", r$ip_ci_lower),
                    if (is.na(r$ip_ci_upper)) "NA" else sprintf("%.4f", r$ip_ci_upper),
                    iprs$C1, iprs$C2, iprs$C3, iprs$C4, iprs$C5, iprs$IPRS,
                    iprs$IPRS_tier, iprs$IPRS_flag))
    }
    cat("\n")
    for (key in names(ip_results)) {
        r <- ip_results[[key]]; if (is.null(r)) next
        if (r$iprs$IPRS_flag != "None")
            cat(sprintf("  [FLAG] %s | IPRS=%.3f (%s) → %s\n",
                        key, r$iprs$IPRS, r$iprs$IPRS_tier, r$iprs$IPRS_flag))
    }


    # =========================================================================
    # SECTION 5b: IP position diagnostics (console only)
    # =========================================================================
    cat("\n\n", strrep("=", 72), "\n")
    cat("SECTION 5b: IP position diagnostics\n")
    cat(strrep("=", 72), "\n\n")

    for (ds in DATASETS) {
        gi <- gam_inputs[[ds]]
        if (is.null(gi)) next

        key <- paste0(ds, "|", GS_KEY)
        r   <- ip_results[[key]]

        goi_vals  <- gi$data$GOI
        goi_range <- diff(range(goi_vals))

        ip_iprs   <- if (!is.null(r)) r$ip_value       else NA_real_
        ip_ci_lo  <- if (!is.null(r)) r$ip_ci_lower    else NA_real_
        ip_ci_hi  <- if (!is.null(r)) r$ip_ci_upper    else NA_real_
        iprs_tier <- if (!is.null(r)) r$iprs$IPRS_tier else "N/A"

        # u_IP = F_GOI(x_IP): where the IP sits on the copula u-axis.
        u_iprs <- if (!is.na(ip_iprs)) mean(goi_vals <= ip_iprs) else NA_real_

        # Mean Si_norm over the ten cells nearest the IP in GOI.
        si_at_iprs_ip <- if (!is.na(ip_iprs) && !is.null(gi$score_df)) {
            nearest_10 <- order(abs(gi$score_df$GOI_expr - ip_iprs))[
                1:min(10L, nrow(gi$score_df))]
            mean(gi$score_df$spear_norm[nearest_10])
        } else NA_real_

        cat(sprintf("  ── %s ──────────────────\n", ds))
        cat(sprintf("  n cells             : %d\n", nrow(gi$data)))
        cat(sprintf("  GOI range           : [%.4f, %.4f]  (range = %.4f)\n",
                    min(goi_vals), max(goi_vals), goi_range))
        cat(sprintf("  IP (IPRS R(x)=0)    : %s    [u = %s]\n",
                    if (is.na(ip_iprs)) "NA" else sprintf("%.4f", ip_iprs),
                    if (is.na(u_iprs))  "NA" else sprintf("%.4f", u_iprs)))
        if (!is.na(ip_ci_lo) && !is.na(ip_ci_hi))
            cat(sprintf("  IPRS 95%% CI         : [%.4f, %.4f]  width = %.4f\n",
                        ip_ci_lo, ip_ci_hi, ip_ci_hi - ip_ci_lo))
        cat(sprintf("  IPRS tier           : %s\n", iprs_tier))
        cat(sprintf("  Mean Si_norm near IP (10 closest cells) : %s\n",
                    if (is.na(si_at_iprs_ip)) "NA"
                    else sprintf("%.6f", si_at_iprs_ip)))
        cat("\n")
    }


    # =========================================================================
    # SECTION 6: Build ip_table and environment vectors
    # =========================================================================
    ip_table <- do.call(rbind, lapply(names(ip_results), function(key) {
        r     <- ip_results[[key]]
        parts <- strsplit(key, "\\|")[[1]]
        if (is.null(r)) return(data.frame(
            Dataset = parts[1], Gene_Set = parts[2],
            IP_Value = NA_real_, IP_CI_Lower = NA_real_, IP_CI_Upper = NA_real_,
            CI_Width = NA_real_,
            IP_Y_Value = NA_real_, IP_Y_CI_Lower = NA_real_,
            IP_Y_CI_Upper = NA_real_,
            C1 = NA_real_, C2 = NA_real_, C3 = NA_real_,
            C4 = NA_real_, C5 = NA_real_, IPRS = NA_real_,
            IPRS_Tier = NA_character_, IPRS_Flag = NA_character_,
            stringsAsFactors = FALSE))
        iprs <- r$iprs
        data.frame(Dataset = parts[1], Gene_Set = parts[2],
                   IP_Value = r$ip_value, IP_CI_Lower = r$ip_ci_lower,
                   IP_CI_Upper = r$ip_ci_upper, CI_Width = iprs$ci_width,
                   # Derived, not estimated. IPRS grades IP_x.
                   IP_Y_Value = r$ip_y_value, IP_Y_CI_Lower = r$ip_y_ci_lower,
                   IP_Y_CI_Upper = r$ip_y_ci_upper,
                   C1 = iprs$C1, C2 = iprs$C2, C3 = iprs$C3,
                   C4 = iprs$C4, C5 = iprs$C5, IPRS = iprs$IPRS,
                   IPRS_Tier = iprs$IPRS_tier, IPRS_Flag = iprs$IPRS_flag,
                   stringsAsFactors = FALSE)
    }))
    if (is.null(ip_table)) {
        cat("  WARNING: ip_table is NULL, no IP results to tabulate.\n")
        ip_table <- data.frame(Dataset = character(0), Gene_Set = character(0),
                               IP_Value = numeric(0), IP_CI_Lower = numeric(0),
                               IP_CI_Upper = numeric(0), CI_Width = numeric(0),
                               IP_Y_Value = numeric(0),
                               IP_Y_CI_Lower = numeric(0),
                               IP_Y_CI_Upper = numeric(0),
                               C1 = numeric(0), C2 = numeric(0), C3 = numeric(0),
                               C4 = numeric(0), C5 = numeric(0), IPRS = numeric(0),
                               IPRS_Tier = character(0), IPRS_Flag = character(0),
                               stringsAsFactors = FALSE)
    }
    rownames(ip_table) <- NULL

    inflection_points <- setNames(
        ip_table$IP_Value,
        paste0(ip_table$Dataset, "|", ip_table$Gene_Set))

    inflection_points_DFG <- setNames(
        ip_table$IP_Value[ip_table$Gene_Set == GS_KEY],
        ip_table$Dataset[ip_table$Gene_Set == GS_KEY])

    # Response-axis classification threshold, keyed the same way.
    inflection_points_DFG_y <- setNames(
        ip_table$IP_Y_Value[ip_table$Gene_Set == GS_KEY],
        ip_table$Dataset[ip_table$Gene_Set == GS_KEY])

    # One population; the Population column is kept.
    ip_table_all <- cbind(Population = "full_GOI_positive", ip_table)
    rownames(ip_table_all) <- NULL

    cat("\n=== Environment vectors created ===\n")
    cat("inflection_points     (", length(inflection_points),
        "entries, keyed 'dataset|cDFG')\n")
    cat("inflection_points_DFG (", length(inflection_points_DFG),
        "entries, keyed by dataset):\n")
    for (nm in names(inflection_points_DFG))
        cat(sprintf("  %-10s  IP = %s\n", nm,
                    if (is.na(inflection_points_DFG[nm])) "NA"
                    else sprintf("%.6f", inflection_points_DFG[nm])))

    cat("\nFinal IP values with IPRS (full GOI-positive population):\n")
    print(ip_table_all[, c("Population", "Dataset", "Gene_Set", "IP_Value",
                           "IP_CI_Lower", "IP_CI_Upper", "C1", "C2", "C3",
                           "C4", "C5", "IPRS", "IPRS_Tier", "IPRS_Flag")],
          row.names = FALSE)



    # =========================================================================
    # SECTION 7: Workbook sheets
    # =========================================================================
    cat("\nExporting results...\n")

    # One builder, called once per population.
    .build_ip_sheet <- function(res_list, pop_label) do.call(rbind,
                                lapply(names(res_list), function(key) {
        r     <- res_list[[key]]; if (is.null(r)) return(NULL)
        parts <- strsplit(key, "\\|")[[1]]
        ds    <- parts[1]
        iprs  <- r$iprs
        goi_label <- GOI
        data.frame(
            Population       = pop_label,
            Dataset          = ds,
            Group            = "SEA-AD_Microglia_Pooled",
            Gene_Set         = parts[2],
            Gene_Set_Label   = GS_LABELS[["cDFG"]],
            Gene_Members     = GS_MEMBERS[["cDFG"]],
            GOI_Label        = goi_label,
            Source_GAM       = basename(p10_rds_path(ds)),
            Donor_RE         = if (is.null(gam_inputs[[ds]])) NA
                               else isTRUE(gam_inputs[[ds]]$donor_re),
            N_Donors         = if (is.null(gam_inputs[[ds]])) NA_integer_
                               else as.integer(gam_inputs[[ds]]$n_donors),
            Resid_Donor_RE   = if (is.null(r$resid_donor_re)) NA
                               else isTRUE(r$resid_donor_re),
            DE_Full_Pop      = if (is.null(gam_inputs[[ds]])) NA_real_
                               else gam_inputs[[ds]]$dev_expl,
            DE_vs_flat_null  = if (is.null(gam_inputs[[ds]])) NA_real_
                               else gam_inputs[[ds]]$dev_expl_flat,
            # The DE that sized the purple set, beside the two reporting DEs.
            DE_Purple_Sizing = if (is.null(gam_inputs[[ds]])) NA_real_
                               else gam_inputs[[ds]]$dev_expl_pop,
            IP_Value         = r$ip_value,
            CI_Lower_95      = r$ip_ci_lower,
            CI_Upper_95      = r$ip_ci_upper,
            CI_Width         = iprs$ci_width,
            GOI_Range        = iprs$goi_range,
            CI_Lower_Method  = if (is.null(r$ci_lower_method)) NA_character_
                               else r$ci_lower_method,
            CI_Upper_Method  = if (is.null(r$ci_upper_method)) NA_character_
                               else r$ci_upper_method,
            # Pre/Post threshold and the IP_x CI mapped through the curve.
            IP_Y_Value       = if (is.null(r$ip_y_value)) NA_real_
                               else r$ip_y_value,
            IP_Y_CI_Lower    = if (is.null(r$ip_y_ci_lower)) NA_real_
                               else r$ip_y_ci_lower,
            IP_Y_CI_Upper    = if (is.null(r$ip_y_ci_upper)) NA_real_
                               else r$ip_y_ci_upper,
            NegPos_Crossings = r$neg_pos_count,
            N_Purple_Cells   = nrow(r$residual_data),
            Stat_Test        = r$stat_test,
            Stat_P_Raw       = r$stat_p_raw,
            Stat_P_Adj       = r$stat_p_adj,
            C1_Precision     = iprs$C1,  C1_Note = iprs$C1_note,
            C2_Statistical   = iprs$C2,  C2_Note = iprs$C2_note,
            C3_Abruptness    = iprs$C3,  C3_Note = iprs$C3_note,
            C4_GAM_Ribbon    = iprs$C4,  C4_Note = iprs$C4_note,
            C5_Smoothness    = iprs$C5,  C5_Note = iprs$C5_note,
            Rx_Amplitude     = iprs$amplitude,
            Rx_MeanRibbon    = iprs$mean_ribbon,
            Rx_TotalVar      = iprs$tv,
            Deriv_at_IP      = iprs$deriv_at_ip,
            IPRS             = iprs$IPRS,
            IPRS_Tier        = iprs$IPRS_tier,
            IPRS_Flag        = iprs$IPRS_flag,
            Method           = "ZeroCrossing|BrentCI|IPRS=(C1+C2+C3+C4+C5)/5",
            stringsAsFactors = FALSE)
    }))

    ip_sheet <- .build_ip_sheet(ip_results, "full_GOI_positive")
    rownames(ip_sheet) <- NULL

    resid_pred_df <- do.call(rbind, lapply(names(ip_results), function(key) {
        r     <- ip_results[[key]]; if (is.null(r)) return(NULL)
        if (nrow(r$residual_pred) == 0L) return(NULL)
        parts <- strsplit(key, "\\|")[[1]]
        p     <- r$residual_pred[seq(1L, nrow(r$residual_pred), by = 10L), ]
        data.frame(Dataset = parts[1], Gene_Set = parts[2], p,
                   stringsAsFactors = FALSE)
    }))
    if (is.null(resid_pred_df) || length(resid_pred_df) == 0)
        resid_pred_df <- data.frame()
    rownames(resid_pred_df) <- NULL

    # =========================================================================
    # SECTION 7b: Workbook export (3 sheets)
    # =========================================================================

    method_notes <- data.frame(
        Step = paste0("Step ", 1:8),
        Description = c(
            paste0("Cell classification: EP = 1 - model_residual\u00b2 / ",
                   "null_residual\u00b2; top (dev_expl \u00d7 n) cells are 'purple'. ",
                   "Input data = Part 10 pooled full GOI-positive population ",
                   "(all regions combined, all GOI > 0 cells). The null model ",
                   "carries the donor random intercept and the same inverse-donor ",
                   "weights as the fitted model."),
            paste0("Signed residuals: residual = observed Expression - ",
                   "population-level GAM fitted value, for each purple cell. ",
                   "Population level means the donor random intercept is excluded ",
                   "from the prediction."),
            paste0("Residual GAM R(x): thin-plate spline (REML) fitted to ",
                   "(GOI, residual) for purple cells, plus s(donor, bs = \"re\"); ",
                   "k = min(20, n_unique-1). Weights recomputed from the purple ",
                   "subset's own donors."),
            "Fine grid: R(x) and SE evaluated at 2000 equally-spaced GOI values.",
            paste0("IP = first neg\u2192pos zero-crossing of the R(x) GAM ",
                   "(linear interpolation). If none, IP=NA."),
            paste0("95% CI: Brent primary; proximity fallback if one bound is NA. ",
                   "CI_Lower/Upper_Method records 'brent' or 'proximity_fallback'."),
            paste0("Stat test: Shapiro-Wilk (\u03b1=", SHAPIRO_ALPHA,
                   ") per group; both normal \u2192 Welch t; either non-normal ",
                   "\u2192 Mann-Whitney U. One test per dataset, so no ",
                   "multiplicity correction is applied and Stat_P_Adj ",
                   "equals Stat_P_Raw."),
            paste0("IPRS = (C1+C2+C3+C4+C5)/5. Tiers: \u22650.80 Very Strong | ",
                   "0.60 Strong | 0.40 Moderate | 0.20 Weak | <0.20 Failed.")),
        stringsAsFactors = FALSE)

    xlsx_path <- file.path(OUT_DIR, paste0(OUT_PREFIX, "_results.xlsx"))
    if (is.null(ip_sheet) || nrow(ip_sheet) == 0)
        cat("  WARNING: ip_sheet empty, IP_Results sheet will be blank\n")

    write_xlsx(list(
        IP_Results           = if (!is.null(ip_sheet) && nrow(ip_sheet) > 0)
                                   ip_sheet else data.frame(),
        Residual_Predictions = if (nrow(resid_pred_df) > 0) resid_pred_df
                               else data.frame(),
        Method_Steps         = method_notes
    ), path = xlsx_path)
    cat("Exported:", basename(xlsx_path), "\n")


    # =========================================================================
    # SECTION 8: Plot colour constants + save helper
    # =========================================================================
    ip_line_color <- "#2E8B57"
    manual_color  <- "#CC0000"

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
            cat("  save_pp PDF error:", conditionMessage(e), "\n")
        })
        tryCatch({
            png(paste0(base, ".png"), width = png_w, height = png_h, res = res)
            print(p)
            dev.off()
        }, error = function(e) {
            .unwind()
            cat("  save_pp PNG error:", conditionMessage(e), "\n")
        })
    }

    iprs_subtitle <- function(r) {
        iprs <- r$iprs
        sprintf("IPRS=%.3f (%s)  |  C1=%.3f  C2=%.3f  C3=%.3f  C4=%.3f  C5=%.3f  |  %s",
                iprs$IPRS, iprs$IPRS_tier,
                iprs$C1, iprs$C2, iprs$C3, iprs$C4, iprs$C5,
                if (iprs$IPRS_flag == "None") "OK" else iprs$IPRS_flag)
    }


    # =========================================================================
    # SECTION 9: Scatter plot (GOI vs cDFG Expression)
    # =========================================================================
    # Purple = TREP (EP-ranked top set), grey = Non-TREP; IP line and 95% CI band.
    create_scatter <- function(key, res_list = ip_results) {
        r <- res_list[[key]]; if (is.null(r)) return(NULL)
        parts     <- strsplit(key, "\\|")[[1]]
        ds        <- parts[1]
        gi        <- gam_inputs[[ds]]
        gam_model <- gi$gam_model
        data_ep   <- gi$data
        # Same population-level DE detect_ip() uses.
        dev_expl  <- gi$dev_expl_pop
        if (is.null(dev_expl) || !is.finite(dev_expl))
            dev_expl <- gi$dev_expl_flat
        if (is.null(gam_model) || is.null(data_ep)) return(NULL)
        pop_lbl   <- "full GOI+"
        goi_label <- GOI

        # Flat null and population-level prediction, as in detect_ip().
        null_fitted  <- fitted(gam(Expression ~ 1, data = data_ep))
        model_fitted <- as.numeric(predict(gam_model, newdata = data_ep,
                                           exclude = DONOR_RE_LABEL,
                                           type = "response"))
        # Same length check as detect_ip().
        if (length(null_fitted) != nrow(data_ep) ||
            length(model_fitted) != nrow(data_ep)) {
            cat("  Scatter skipped for", ds, "- the EP vectors are the wrong",
                "length (null", length(null_fitted), ", model",
                length(model_fitted), ", rows", nrow(data_ep), ")\n")
            return(NULL)
        }
        null_sq      <- (data_ep$Expression - null_fitted)^2
        model_sq     <- (data_ep$Expression - model_fitted)^2
        expl_power   <- ifelse(null_sq > 0, 1 - model_sq / null_sq, 0)

        # Integer index; purple count clamped to [1, n] as in detect_ip().
        row_idx   <- seq_len(nrow(data_ep))
        n_dev     <- max(1L, min(as.integer(round(nrow(data_ep) *
                                                  purple_count_fraction(dev_expl))),
                                 nrow(data_ep)))
        dev_idx   <- sort(order(expl_power, decreasing = TRUE)[seq_len(n_dev)])

        # TREP = EP-ranked top set; factor keeps TREP first in the legend.
        plot_data <- data.frame(
            GOI        = data_ep$GOI,
            Expression = data_ep$Expression,
            Group      = factor(ifelse(row_idx %in% dev_idx,
                                       "TREP", "Non-TREP"),
                                levels = c("TREP", "Non-TREP")))
        plot_data <- plot_data[order(plot_data$Group == "TREP"), ]

        pred_data     <- data.frame(GOI = seq(min(plot_data$GOI),
                                              max(plot_data$GOI), length.out = 200L))
        pg            <- predict_population(gam_model, pred_data$GOI, data_ep,
                                            se = TRUE)
        pred_data$fit    <- as.numeric(pg$fit)
        pred_data$se.fit <- as.numeric(pg$se.fit)

        lc        <- ip_line_color
        title_col <- if (r$manual_flag) manual_color else "black"

        title <- if (r$manual_flag || is.na(r$ip_value)) {
            paste0(ds, " | ", DFG_COMPOSITE_LABEL, " Pooled (", pop_lbl,
                   ")  *** NO IP DETECTED ***")
        } else {
            paste0(ds, " | ", DFG_COMPOSITE_LABEL, " Pooled (", pop_lbl, ")",
                   "  IP(IPRS): ", sprintf("%.4f", r$ip_value),
                   "  95% CI [",
                   ifelse(is.na(r$ip_ci_lower), "NA", sprintf("%.4f", r$ip_ci_lower)),
                   ", ",
                   ifelse(is.na(r$ip_ci_upper), "NA", sprintf("%.4f", r$ip_ci_upper)),
                   "]")
        }

        p <- ggplot()
        if (!is.na(r$ip_ci_lower) && !is.na(r$ip_ci_upper))
            p <- p + annotate("rect",
                              xmin = r$ip_ci_lower, xmax = r$ip_ci_upper,
                              ymin = -Inf, ymax = Inf, fill = lc, alpha = 0.08)

        p <- p +
            geom_point(data = plot_data,
                       aes(x = GOI, y = Expression, color = Group),
                       size = 1.8, alpha = 0.45) +
            geom_ribbon(data = pred_data,
                        aes(x = GOI, ymin = fit - 1.96 * se.fit,
                            ymax = fit + 1.96 * se.fit),
                        fill = "#FFCC99", alpha = 0.20, na.rm = TRUE) +
            geom_line(data = pred_data, aes(x = GOI, y = fit),
                      color = "#FFCC99", linewidth = 1.2, na.rm = TRUE)

        if (!is.na(r$ip_value))
            p <- p + geom_vline(xintercept = r$ip_value,
                                linetype = "dashed", color = lc, linewidth = 0.8)
        if (!is.na(r$ip_ci_lower))
            p <- p + geom_vline(xintercept = r$ip_ci_lower,
                                linetype = "dotdash", color = lc, linewidth = 0.6)
        if (!is.na(r$ip_ci_upper))
            p <- p + geom_vline(xintercept = r$ip_ci_upper,
                                linetype = "dotdash", color = lc, linewidth = 0.6)

        # Horizontal line: Pre/Post classification threshold IP_y = f(IP_x).
        if (isTRUE(IP_Y_ENABLE) && !is.null(r$ip_y_value) &&
            is.finite(r$ip_y_value))
            p <- p + geom_hline(yintercept = r$ip_y_value,
                                linetype = IP_Y_LINE_TYPE,
                                color = IP_Y_LINE_COL,
                                linewidth = 0.8)

        # Caption built next to the line so the two agree.
        ip_y_val_sub  <- if (is.null(r$ip_y_value)) NA_real_ else r$ip_y_value
        ip_y_drawn    <- isTRUE(IP_Y_ENABLE) && is.finite(ip_y_val_sub)
        ip_y_subtitle <- if (!is.finite(ip_y_val_sub)) {
            "IP_y unavailable for this dataset"
        } else {
            paste0("IP_y = ", sprintf("%.4f", ip_y_val_sub),
                   if (ip_y_drawn) " (horizontal dashed line)" else "",
                   " = f(IP_x): the fitted population curve read at the IP")
        }

        p + scale_color_manual(
                values = c("TREP"     = "#4B0082",
                           "Non-TREP" = "#C0C0C0")) +
            theme_minimal() +
            theme(
                panel.grid.major = element_line(color = "#EEEEEE"),
                panel.grid.minor = element_line(color = "#F5F5F5"),
                legend.position  = "bottom",
                plot.title       = element_text(size = 10, face = "bold",
                                                hjust = 0.5, color = title_col),
                plot.subtitle    = element_text(size = 8, hjust = 0.5),
                panel.border     = element_rect(color = "black", fill = NA,
                                                linewidth = 0.5),
                axis.line        = element_blank()) +
            labs(title    = title,
                 subtitle = paste0(ip_y_subtitle, "\n", iprs_subtitle(r)),
                 x        = paste0(goi_label, " Expression (", NORM_TAG,
                                    " log1p, GOI>0)"),
                 y        = paste0(DFG_COMPOSITE_LABEL,
                                    " Gene Set Mean Expression (", NORM_TAG,
                                    " log1p)"),
                 color    = "Cell group") +
            scale_x_continuous(
                limits = range(data_ep$GOI, na.rm = TRUE) +
                    c(-1, 1) * 0.05 * diff(range(data_ep$GOI, na.rm = TRUE)))
    }


    # =========================================================================
    # SECTION 10: R(x) diagnostic plot
    # =========================================================================
    create_rx <- function(key, res_list = ip_results, pop_lbl = "full GOI+") {
        r <- res_list[[key]]; if (is.null(r)) return(NULL)
        parts     <- strsplit(key, "\\|")[[1]]
        ds        <- parts[1]
        goi_label <- GOI

        pred      <- r$residual_pred
        rd        <- r$residual_data
        iprs      <- r$iprs
        lc        <- ip_line_color
        title_col <- if (r$manual_flag) manual_color else "black"

        n_bins_rx <- 100L
        p_goi     <- rd$GOI

        if (nrow(rd) > 0L) {
            if (length(unique(p_goi)) > n_bins_rx) {
                breaks_rx  <- seq(min(p_goi), max(p_goi),
                                  length.out = n_bins_rx + 1L)
                bin_ids_rx <- cut(p_goi, breaks = breaks_rx,
                                  include.lowest = TRUE, labels = FALSE)
                bin_ctr_rx <- (breaks_rx[-length(breaks_rx)] + breaks_rx[-1L]) / 2
                bin_w_rx   <- breaks_rx[2L] - breaks_rx[1L]
            } else {
                uq_rx      <- sort(unique(p_goi))
                bin_ids_rx <- match(p_goi, uq_rx)
                bin_ctr_rx <- uq_rx
                bin_w_rx   <- if (length(uq_rx) > 1L) min(diff(uq_rx)) else 0.1
            }
            bs_rx <- data.frame(
                goi_center    = bin_ctr_rx,
                mean_residual = sapply(seq_along(bin_ctr_rx), function(b) {
                    idx <- which(bin_ids_rx == b)
                    if (length(idx) > 0L) mean(rd$residual[idx]) else NA_real_
                }))
            bs_rx    <- bs_rx[!is.na(bs_rx$mean_residual), ]
            # Bar width from the spacing between surviving bin centres.
            bar_w_rx <- if (nrow(bs_rx) > 1L)
                min(diff(bs_rx$goi_center)) * 0.8
            else bin_w_rx * 0.8
        } else {
            bs_rx    <- data.frame(goi_center = numeric(0),
                                   mean_residual = numeric(0))
            bar_w_rx <- 0.1
        }

        stat_txt <- if (!is.na(r$stat_test))
            sprintf("%s | p = %.3e (unadjusted)", r$stat_test,
                    ifelse(is.na(r$stat_p_adj), NA, r$stat_p_adj))
        else "Stat: NA"

        subtitle <- if (r$status == "insufficient_purple_cells") {
            paste0("*** FAILED - insufficient purple cells (",
                   nrow(rd), " cells, ", length(unique(rd$GOI)),
                   " unique GOI values) ***")
        } else if (r$manual_flag) {
            paste0("*** NO IP - R(x) has no neg\u2192pos crossing ***\n",
                   iprs_subtitle(r))
        } else {
            paste0("IP = ", sprintf("%.4f", r$ip_value),
                   "  95% CI [",
                   ifelse(is.na(r$ip_ci_lower), "NA", sprintf("%.4f", r$ip_ci_lower)),
                   ", ",
                   ifelse(is.na(r$ip_ci_upper), "NA", sprintf("%.4f", r$ip_ci_upper)),
                   "]  |  ", stat_txt, "\n", iprs_subtitle(r))
        }

        p <- ggplot()
        if (!is.na(r$ip_ci_lower) && !is.na(r$ip_ci_upper))
            p <- p + annotate("rect",
                              xmin = r$ip_ci_lower, xmax = r$ip_ci_upper,
                              ymin = -Inf, ymax = Inf, fill = lc, alpha = 0.08)

        if (nrow(bs_rx) > 0L)
            p <- p + geom_bar(data = bs_rx,
                              aes(x = goi_center, y = mean_residual),
                              stat = "identity", fill = "#4B0082", color = "#3D0065",
                              alpha = 0.40, width = bar_w_rx, linewidth = 0.4)

        if (nrow(pred) > 0L)
            p <- p +
                geom_ribbon(data = pred,
                            aes(x = GOI, ymin = r_lower, ymax = r_upper),
                            fill = "#FFCC99", alpha = 0.35, na.rm = TRUE) +
                geom_line(data = pred, aes(x = GOI, y = r_hat),
                          color = "#FFCC99", linewidth = 1.5, na.rm = TRUE)

        p <- p + geom_hline(yintercept = 0, linetype = "dotted",
                            color = "gray40", linewidth = 0.5)
        if (!is.na(r$ip_ci_lower))
            p <- p + geom_vline(xintercept = r$ip_ci_lower,
                                linetype = "dotdash", color = lc, linewidth = 0.7)
        if (!is.na(r$ip_ci_upper))
            p <- p + geom_vline(xintercept = r$ip_ci_upper,
                                linetype = "dotdash", color = lc, linewidth = 0.7)
        if (!is.na(r$ip_value))
            p <- p + geom_vline(xintercept = r$ip_value,
                                linetype = "dashed", color = lc, linewidth = 1.0)

        p + theme_minimal() +
            theme(panel.border  = element_rect(color = "black", fill = NA,
                                               linewidth = 0.5),
                  axis.line     = element_blank(),
                  plot.title    = element_text(size = 10, face = "bold",
                                               hjust = 0.5, color = title_col),
                  plot.subtitle = element_text(size = 7.5, hjust = 0.5),
                  axis.text     = element_text(size = 7)) +
            labs(title    = paste0(ds, " | ", DFG_COMPOSITE_LABEL,
                                   " Pooled (", pop_lbl, ")",
                                   " -- R(x) = E[residual | purple, GOI=x]",
                                   "  [Zero-crossing | Brent CI | IPRS]"),
                 subtitle = subtitle,
                 x        = paste0(goi_label, " Expression (", NORM_TAG,
                                    " log1p, GOI>0)"),
                 y        = "Expected Signed Residual R(x)")
    }


    # =========================================================================
    # SECTION 11: Generate, emit to Plots pane, and save all plots
    # =========================================================================
    cat("\n\nGenerating scatter plots...\n")
    for (key in names(ip_results)) {
        tryCatch({
            p <- create_scatter(key)
            if (!is.null(p)) {
                emit_plot_iprs(p)
                ds <- strsplit(key, "\\|")[[1]][1]
                save_pp(p,
                        file.path(OUT_DIR, paste0(OUT_PREFIX, "_", ds, "_scatter")),
                        pdf_w = 10, pdf_h = 9, png_w = 2000, png_h = 1800)
                cat("  Scatter saved:", ds, "\n")
            }
        }, error = function(e) cat("Scatter error:", key, conditionMessage(e), "\n"))
    }

    cat("\nGenerating R(x) diagnostic plots...\n")
    for (key in names(ip_results)) {
        tryCatch({
            p <- create_rx(key)
            if (!is.null(p)) {
                emit_plot_iprs(p)
                ds <- strsplit(key, "\\|")[[1]][1]
                save_pp(p,
                        file.path(OUT_DIR, paste0(OUT_PREFIX, "_", ds, "_Rx")),
                        pdf_w = 10, pdf_h = 9, png_w = 2000, png_h = 1800)
                cat("  R(x) plot saved:", ds, "\n")
            }
        }, error = function(e) cat("R(x) error:", key, conditionMessage(e), "\n"))
    }


    # =========================================================================
    # SECTION 11b: Expression-space pitcher plots
    # =========================================================================
    # x = GOI, y = composite expression (log1p), all GOI > 0 cells.
    # Fill: K_i or S_i_norm (turbo); 2D density background; GAM curve and IP.

    cat("\nGenerating expression-space pitcher plots...\n")

    # Theme matching Part10d's theme_copula
    theme_pitcher_expr <- function(base_size = 10) {
        theme_minimal(base_size = base_size) +
            theme(
                plot.title       = element_text(face = "bold", size = base_size + 1),
                plot.subtitle    = element_text(size = base_size - 1,
                                                color = "gray30"),
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

    # == Helper: build one expression-space pitcher plot =====================
    # fill_vec   : numeric vector (same length as nrow(cell_df)) mapped to fill
    # fill_name  : legend title (expression or string)
    # fill_limits: length-2 numeric for turbo palette limits
    # fill_label : short string for subtitle (e.g. "K_i", "S_i_norm")
    # ip_y_val   : IP_y; NA draws no horizontal line
    build_pitcher_expr <- function(cell_df, dens_grid, dens_max,
                                    pred_x,
                                    med_goi, med_cdfg,
                                    ip_val, ip_lo, ip_hi,
                                    fill_col, fill_name, fill_limits, fill_label,
                                    ds, n_full, dev_pct, goi_label,
                                    ip_y_val = NA_real_) {

        p <- ggplot() +
            # Layer 1: 2D density background, same purple gradient as copula
            geom_raster(data = dens_grid,
                        aes(x = x, y = y, fill = log1p_dens),
                        interpolate = TRUE, na.rm = TRUE) +
            scale_fill_gradientn(
                colours = c("#FFFFFF00", "#FFFFFF00", "#E8E0FFA0", "#C8B8F8C8",
                            "#9E7AE8E0", "#7B2FBEEE", "#4B0082F8", "#2A0050FF"),
                values  = scales::rescale(c(0, 0.20, 0.35, 0.50,
                                            0.65, 0.78, 0.90, 1.00)),
                limits  = c(0, dens_max),
                oob     = scales::squish,
                name    = "Density\n(log1p)",
                na.value = "#FFFFFF00"
            ) +
            # Layer 2: Cell points, fill colour only, uniform size (no EP)
            ggnewscale::new_scale_fill() +
            geom_point(data = cell_df,
                       aes(x = GOI, y = cDFG, fill = .data[[fill_col]]),
                       shape = 21, size = 2.2, alpha = 0.72,
                       colour = "gray15", stroke = 0.15) +
            {
                # Quantile-based colour placement: map turbo colour stops to
                # data quantiles so each colour band contains equal cell counts.
                # This prevents skewed distributions from producing monotone fills.
                fill_data <- cell_df[[fill_col]]
                n_col     <- 9L   # number of turbo colours
                q_probs   <- seq(0, 1, length.out = n_col)
                q_vals    <- quantile(fill_data, probs = q_probs, na.rm = TRUE)
                q_mapped  <- (q_vals - fill_limits[1]) / diff(fill_limits)
                q_mapped  <- pmin(pmax(q_mapped, 0), 1)
                # Enforce strictly increasing (ties possible in quantiles)
                for (ii in 2:length(q_mapped))
                    if (q_mapped[ii] <= q_mapped[ii - 1])
                        q_mapped[ii] <- q_mapped[ii - 1] + 1e-6
                scale_fill_gradientn(
                    colours = viridis::turbo(256)[round(seq(0.05 * 255 + 1,
                                                           0.95 * 255 + 1,
                                                           length.out = n_col))],
                    values  = q_mapped,
                    limits  = fill_limits,
                    oob     = scales::squish,
                    name    = fill_name
                )
            }

        # Layer 4: GAM curve, dashed dark teal with SE ribbon
        p <- p +
            geom_ribbon(data = pred_x,
                        aes(x = GOI, ymin = lower, ymax = upper),
                        inherit.aes = FALSE,
                        fill = "#006D6D", alpha = 0.15, na.rm = TRUE) +
            geom_line(data = pred_x,
                      aes(x = GOI, y = fit),
                      inherit.aes = FALSE,
                      colour = "#006D6D", linewidth = 1.0,
                      linetype = "dashed", na.rm = TRUE)

        # Layer 5: IP + CI + reference lines (use ip_line_color = green)
        if (!is.na(ip_lo) && !is.na(ip_hi))
            p <- p +
                annotate("rect",
                         xmin = ip_lo, xmax = ip_hi,
                         ymin = -Inf,  ymax = Inf,
                         fill = "#2E8B57", alpha = 0.06)
        if (!is.na(ip_val))
            p <- p +
                geom_vline(xintercept = ip_val, linetype = "dashed",
                           colour = "#2E8B57", linewidth = 0.55)
        p <- p +
            geom_vline(xintercept = med_goi,  linetype = "dashed",
                       colour = "gray50", linewidth = 0.4) +
            geom_hline(yintercept = med_cdfg, linetype = "dashed",
                       colour = "gray50", linewidth = 0.4)

        # IP_y classification threshold, near-black dashed.
        if (isTRUE(IP_Y_ENABLE) && is.finite(ip_y_val))
            p <- p + geom_hline(yintercept = ip_y_val,
                                linetype = IP_Y_LINE_TYPE,
                                colour = IP_Y_LINE_COL, linewidth = 0.7)

        # Title, subtitle, axes, theme
        ip_str <- if (is.na(ip_val)) "NA"
                  else sprintf("%.4f  95%% CI [%s, %s]",
                               ip_val,
                               ifelse(is.na(ip_lo), "NA", sprintf("%.4f", ip_lo)),
                               ifelse(is.na(ip_hi), "NA", sprintf("%.4f", ip_hi)))

        lim_str <- sprintf("%.2f-%.2f", fill_limits[1], fill_limits[2])

        p + labs(title    = paste0(ds,
                                   " -- Expression-space pitcher geometry",
                                   "  [n=", n_full, "]"),
                 subtitle = paste0(
                     "Background: 2D density (log1p).  ",
                     "Point fill: ", fill_label, " (turbo, ", lim_str, ").  ",
                     "Donor-aware DE: ", dev_pct, "%.\n",
                     "Dashed dark teal: GAM fit.  ",
                     "Green vertical: x_IP = ", ip_str, ".\n",
                     "Black dashed horizontal: y_IP = f(x_IP) = ",
                     if (is.finite(ip_y_val)) sprintf("%.4f", ip_y_val) else "NA",
                     " -- the Pre/Post classification threshold.  ",
                     "Grey dashed: marginal medians."),
                 x = paste0(goi_label, " Expression (", NORM_TAG, " log1p)"),
                 y = paste0(DFG_COMPOSITE_LABEL, " Expression (", NORM_TAG,
                             " log1p)")) +
            theme_pitcher_expr() +
            theme(legend.position = "right")
    }


    # == Generate Ki and Si pitcher plots for each dataset ===================
    for (key in names(ip_results)) {
        tryCatch({
            r     <- ip_results[[key]]; if (is.null(r)) next
            parts <- strsplit(key, "\\|")[[1]]
            ds    <- parts[1]
            gi    <- gam_inputs[[ds]]
            if (is.null(gi) || is.null(gi$score_df)) {
                cat("  Pitcher plot skipped:", ds, "(no score_df)\n")
                next
            }

            gam_model <- gi$gam_model
            score_df  <- gi$score_df
            data_full <- gi$data        # GOI, Expression
            goi_label <- GOI
            n_full    <- nrow(data_full)

            cat("  Building expression-space pitcher plots:", ds,
                "| n =", n_full, "\n")

            # == Score vectors from SKCS decomposition (all cells) ===========
            Ki_all      <- score_df$kend_contrib
            Si_norm_all <- score_df$spear_norm

            cat(sprintf("    Ki range:      [%.4f, %.4f]\n",
                        min(Ki_all), max(Ki_all)))
            cat(sprintf("    Si_norm range: [%.4f, %.4f]\n",
                        min(Si_norm_all), max(Si_norm_all)))

            # == 2D density background (computed once, shared) ===============
            dens_2d <- kde2d(data_full$GOI, data_full$Expression,
                             n = 200,
                             lims = c(range(data_full$GOI)        + c(-0.3, 0.3),
                                      range(data_full$Expression) + c(-0.3, 0.3)))
            dens_grid <- expand.grid(x = dens_2d$x, y = dens_2d$y)
            dens_grid$density     <- as.vector(dens_2d$z)
            dens_grid$log1p_dens  <- log1p(dens_grid$density)
            dens_max              <- max(dens_grid$log1p_dens)

            # == Cell-level data frame (both Ki and Si_norm present) =========
            cell_df <- data.frame(
                GOI     = data_full$GOI,
                cDFG    = data_full$Expression,
                Ki      = Ki_all,
                Si_norm = Si_norm_all,
                stringsAsFactors = FALSE
            )

            # == GAM prediction curve + SE ribbon ============================
            set.seed(123)
            pred_x  <- data.frame(GOI = seq(min(data_full$GOI),
                                            max(data_full$GOI),
                                            length.out = 1000L))
            pg_pred <- predict_population(gam_model, pred_x$GOI, data_full,
                                          se = TRUE)
            pred_x$fit     <- as.numeric(pg_pred$fit)
            pred_x$se.fit  <- as.numeric(pg_pred$se.fit)
            pred_x$lower   <- pred_x$fit - 1.96 * pred_x$se.fit
            pred_x$upper   <- pred_x$fit + 1.96 * pred_x$se.fit

            # == Reference lines and IP ======================================
            med_goi  <- median(data_full$GOI)
            med_cdfg <- median(data_full$Expression)
            ip_val   <- r$ip_value
            ip_lo    <- r$ip_ci_lower
            ip_hi    <- r$ip_ci_upper
            # IP_y classification threshold on both pitchers.
            ip_y_val <- if (is.null(r$ip_y_value)) NA_real_ else r$ip_y_value
            dev_pct  <- round(gi$dev_expl * 100, 2)

            # == (a) Ki pitcher plot =========================================
            # Colour limits from the data range.
            ki_range <- range(Ki_all, na.rm = TRUE)

            p_ki <- build_pitcher_expr(
                cell_df    = cell_df,
                dens_grid  = dens_grid,
                dens_max   = dens_max,
                pred_x     = pred_x,
                med_goi    = med_goi,
                med_cdfg   = med_cdfg,
                ip_val     = ip_val,
                ip_lo      = ip_lo,
                ip_hi      = ip_hi,
                fill_col   = "Ki",
                fill_name  = "Kendall\nconcordance\ncontribution",
                fill_limits = ki_range,
                fill_label = "K_i",
                ds         = ds,
                n_full     = n_full,
                dev_pct    = dev_pct,
                goi_label  = goi_label,
                ip_y_val   = ip_y_val)

            emit_plot_iprs(p_ki)
            save_pp(p_ki,
                    file.path(OUT_DIR, paste0(OUT_PREFIX, "_", ds,
                                              "_pitcher_Ki_expr")),
                    pdf_w = 11, pdf_h = 10, png_w = 2400, png_h = 2200)
            cat("    Ki pitcher plot saved:", ds, "\n")

            # == (b) Si_norm pitcher plot ====================================
            # Colour limits from the data range.
            si_range <- range(Si_norm_all, na.rm = TRUE)

            p_si <- build_pitcher_expr(
                cell_df    = cell_df,
                dens_grid  = dens_grid,
                dens_max   = dens_max,
                pred_x     = pred_x,
                med_goi    = med_goi,
                med_cdfg   = med_cdfg,
                ip_val     = ip_val,
                ip_lo      = ip_lo,
                ip_hi      = ip_hi,
                fill_col   = "Si_norm",
                fill_name  = "Normalised\nSpearman\ncontribution",
                fill_limits = si_range,
                fill_label = "S_i_norm",
                ds         = ds,
                n_full     = n_full,
                dev_pct    = dev_pct,
                goi_label  = goi_label,
                ip_y_val   = ip_y_val)

            emit_plot_iprs(p_si)
            save_pp(p_si,
                    file.path(OUT_DIR, paste0(OUT_PREFIX, "_", ds,
                                              "_pitcher_Si_expr")),
                    pdf_w = 11, pdf_h = 10, png_w = 2400, png_h = 2200)
            cat("    Si pitcher plot saved:", ds, "\n")

        }, error = function(e) {
            cat("  Pitcher plot error:", key, conditionMessage(e), "\n")
        })
    }


    # =========================================================================
    # SECTION 12: IPRS tier summary + CSV
    # =========================================================================
    cat("\nIPRS tier summary -- cDFG (", length(ip_results), "datasets):\n")
    tier_counts <- table(sapply(ip_results,
                                function(r) if (!is.null(r)) r$iprs$IPRS_tier
                                            else "NULL"))
    for (tier in c("Very Strong", "Strong", "Moderate", "Weak", "Failed"))
        cat(sprintf("  %-12s : %d\n", tier,
                    if (tier %in% names(tier_counts)) tier_counts[[tier]] else 0L))

    # One population; the Population column is kept.
    write.csv(ip_table_all,
              file.path(OUT_DIR, paste0(OUT_PREFIX, "_summary.csv")),
              row.names = FALSE)
    cat("Summary CSV saved:", paste0(OUT_PREFIX, "_summary.csv"), "\n")


    # =========================================================================
    # SECTION 13: Save standalone RDS for Part 10b + verify + log
    # =========================================================================

    # == Helper: dual-sink writer ============================================
    LOG_LINES <- character(0)
    log_cat <- function(...) {
        msg <- paste0(...)
        cat(msg)
        LOG_LINES <<- c(LOG_LINES, msg)
    }

    log_cat("\n", strrep("=", 72), "\n")
    log_cat("  SECTION 13: Saving ", OUT_PREFIX, "_results.rds + Verification\n")
    log_cat("  Run timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
    log_cat(strrep("=", 72), "\n\n")

    # == 13a. Build and save the RDS =========================================
    iprs_rds_path <- file.path(OUT_DIR, paste0(OUT_PREFIX, "_results.rds"))

    iprs_save_obj <- list(
        inflection_points     = inflection_points,
        # inflection_points_DFG is the field Part 10b reads (full GOI-positive IP).
        inflection_points_DFG      = inflection_points_DFG,
        inflection_points_DFG_full = inflection_points_DFG,
        ip_population_for_part10b  = "full",
        ip_definition              = "IPRS_full",
        # IP_y = f(IP_x): the value Part 10b splits Pre/Post on.
        inflection_points_DFG_y      = inflection_points_DFG_y,
        inflection_points_DFG_y_full = inflection_points_DFG_y,
        ip_y_definition            = "f(IP_x): population-level fitted cDFG at the reported IP",
        ip_y_enabled               = IP_Y_ENABLE,
        ip_table_all               = ip_table_all,
        ip_method_used        = "IPRS",
        ip_table              = ip_table,
        ip_results            = ip_results,
        DATASETS              = DATASETS,
        arm                   = ARM_TAG,
        norm_tag              = NORM_TAG,
        goi                   = GOI,
        composite             = DFG_COMPOSITE_LABEL,
        # The estimator the IP was produced under.
        gam_donor_re          = GAM_DONOR_RE,
        resid_donor_re_switch = RESID_DONOR_RE,
        donor_re_min_levels   = DONOR_RE_MIN_LEVELS,
        w_donor_weights       = W_DONOR_WEIGHTS,
        full_pop_donor_re     = setNames(
            sapply(DATASETS, function(d)
                if (is.null(gam_inputs[[d]])) NA else
                    isTRUE(gam_inputs[[d]]$donor_re)), DATASETS),
        full_pop_n_donors     = setNames(
            sapply(DATASETS, function(d)
                if (is.null(gam_inputs[[d]])) NA_integer_ else
                    as.integer(gam_inputs[[d]]$n_donors)), DATASETS),
        dev_expl              = setNames(
            sapply(DATASETS, function(d)
                if (is.null(gam_inputs[[d]])) NA_real_ else
                    gam_inputs[[d]]$dev_expl), DATASETS),
        dev_expl_flat_null    = setNames(
            sapply(DATASETS, function(d)
                if (is.null(gam_inputs[[d]])) NA_real_ else
                    gam_inputs[[d]]$dev_expl_flat), DATASETS),
        source_rds            = setNames(lapply(DATASETS, p10_rds_path), DATASETS),
        timestamp             = Sys.time()
    )
    saveRDS(iprs_save_obj, file = iprs_rds_path)

    # Copy into Part 10b's arm folder, where Part 10b reads it.
    tryCatch({
        dir.create(P10B_DIR, showWarnings = FALSE, recursive = TRUE)
        saveRDS(iprs_save_obj, file = P10B_IPRS_RDS)
        log_cat("  Part 10b copy : ", P10B_IPRS_RDS, "\n")
    }, error = function(e)
        log_cat("  WARNING: could not write the Part 10b copy: ",
                conditionMessage(e), "\n"))

    log_cat("  RDS saved : ", iprs_rds_path, "\n")
    log_cat("  File size : ",
            sprintf("%.1f KB", file.size(iprs_rds_path) / 1e3), "\n")
    log_cat("  mtime     : ",
            format(file.mtime(iprs_rds_path), "%Y-%m-%d %H:%M:%S"), "\n\n")

    # == 13b. Verification on fresh read-back ================================
    #
    #  CHECK 1  inflection_points_DFG names match DATASETS
    #  CHECK 2  IP values match ip_table$IP_Value exactly (no serialisation loss)
    #  CHECK 3  ip_table dimensions match number of processed datasets
    #  CHECK 4  ip_results list is intact (same keys, non-NULL entries retained)
    #  CHECK 5  inflection_points_DFG_y survives the round trip and
    #           ip_y_definition is present
    log_cat(strrep("-", 72), "\n")
    log_cat("  VERIFICATION (independent read-back from disk)\n")
    log_cat(strrep("-", 72), "\n\n")

    verify <- readRDS(iprs_rds_path)

    verify_results <- list()

    # CHECK 1: inflection_points_DFG names
    exp_names <- sort(names(inflection_points_DFG))
    obs_names <- sort(names(verify$inflection_points_DFG))
    ck1_pass  <- identical(obs_names, exp_names)
    ck1_detail <- sprintf("expected names=%s  read_back names=%s",
                           paste(exp_names, collapse=","),
                           paste(obs_names,  collapse=","))
    verify_results[["check1"]] <- list(
        status = if (ck1_pass) "PASS" else "FAIL", detail = ck1_detail)
    log_cat(sprintf("  CHECK 1 [inflection_points_DFG names]  %s  %s\n",
                    if (ck1_pass) "PASS" else "FAIL", ck1_detail))

    # CHECK 2: IP values, exact numeric match via all.equal
    ck2_pass   <- TRUE
    ck2_detail <- character(0)
    for (ds in names(verify$inflection_points_DFG)) {
        ip_orig  <- inflection_points_DFG[ds]
        ip_readb <- verify$inflection_points_DFG[ds]
        both_na  <- is.na(ip_orig) && is.na(ip_readb)
        num_ok   <- both_na || isTRUE(all.equal(as.numeric(ip_orig),
                                                 as.numeric(ip_readb),
                                                 tolerance = 1e-12))
        if (!num_ok) {
            ck2_pass   <- FALSE
            ck2_detail <- c(ck2_detail,
                            sprintf("%s: orig=%.10f  readback=%.10f",
                                    ds, ip_orig, ip_readb))
        } else {
            ck2_detail <- c(ck2_detail,
                            sprintf("%s: %.6f OK", ds,
                                    ifelse(is.na(ip_orig), NA, ip_orig)))
        }
    }
    ck2_str <- paste(ck2_detail, collapse = " | ")
    verify_results[["check2"]] <- list(
        status = if (ck2_pass) "PASS" else "FAIL", detail = ck2_str)
    log_cat(sprintf("  CHECK 2 [IP value fidelity]            %s  %s\n",
                    if (ck2_pass) "PASS" else "FAIL", ck2_str))

    # CHECK 3: ip_table row count
    exp_rows <- length(ip_results)
    obs_rows <- nrow(verify$ip_table)
    ck3_pass  <- obs_rows == exp_rows
    ck3_detail <- sprintf("expected=%d rows  read_back=%d rows", exp_rows, obs_rows)
    verify_results[["check3"]] <- list(
        status = if (ck3_pass) "PASS" else "FAIL", detail = ck3_detail)
    log_cat(sprintf("  CHECK 3 [ip_table row count]           %s  %s\n",
                    if (ck3_pass) "PASS" else "FAIL", ck3_detail))

    # CHECK 4: ip_results keys and non-NULL structure
    exp_keys <- sort(names(ip_results))
    obs_keys <- sort(names(verify$ip_results))
    ck4_pass  <- identical(obs_keys, exp_keys)
    ck4_detail <- sprintf("keys match=%s  n_keys=%d", ck4_pass, length(obs_keys))
    verify_results[["check4"]] <- list(
        status = if (ck4_pass) "PASS" else "FAIL", detail = ck4_detail)
    log_cat(sprintf("  CHECK 4 [ip_results keys]              %s  %s\n",
                    if (ck4_pass) "PASS" else "FAIL", ck4_detail))

    # CHECK 5: IP_y round trip and ip_y_definition.
    exp_y_names <- sort(names(inflection_points_DFG_y))
    obs_y_names <- sort(names(verify$inflection_points_DFG_y))
    ck5_pass    <- identical(obs_y_names, exp_y_names)
    ck5_detail  <- character(0)
    if (!ck5_pass) {
        ck5_detail <- sprintf("names differ: expected=%s  read_back=%s",
                              paste(exp_y_names, collapse = ","),
                              paste(obs_y_names,  collapse = ","))
    } else {
        if (length(obs_y_names) == 0L)
            ck5_detail <- "no IP_y entries (no dataset produced an IP)"
        for (ds in obs_y_names) {
            ipy_orig  <- inflection_points_DFG_y[ds]
            ipy_readb <- verify$inflection_points_DFG_y[ds]
            both_na   <- is.na(ipy_orig) && is.na(ipy_readb)
            num_ok    <- both_na || isTRUE(all.equal(as.numeric(ipy_orig),
                                                      as.numeric(ipy_readb),
                                                      tolerance = 1e-12))
            if (!num_ok) {
                ck5_pass   <- FALSE
                ck5_detail <- c(ck5_detail,
                                sprintf("%s: orig=%.12f  readback=%.12f", ds,
                                        as.numeric(ipy_orig),
                                        as.numeric(ipy_readb)))
            } else {
                ck5_detail <- c(ck5_detail,
                                sprintf("%s: %s OK", ds,
                                        if (is.na(ipy_orig)) "NA"
                                        else sprintf("%.6f", ipy_orig)))
            }
        }
        # An absent ip_y_definition is a FAIL.
        if (is.null(verify$ip_y_definition) ||
            !nzchar(as.character(verify$ip_y_definition)[1])) {
            ck5_pass   <- FALSE
            ck5_detail <- c(ck5_detail, "ip_y_definition absent or empty")
        }
    }
    ck5_str <- paste(ck5_detail, collapse = " | ")
    verify_results[["check5"]] <- list(
        status = if (ck5_pass) "PASS" else "FAIL", detail = ck5_str)
    log_cat(sprintf("  CHECK 5 [IP_y value fidelity]          %s  %s\n",
                    if (ck5_pass) "PASS" else "FAIL", ck5_str))

    # Summary
    overall_fail <- any(sapply(verify_results, function(v) v$status == "FAIL"))
    log_cat("\n")
    log_cat(strrep("-", 72), "\n")
    log_cat(sprintf("  OVERALL RESULT: %s\n",
                    if (!overall_fail) "all checks passed"
                    else "one or more checks failed"))
    log_cat(strrep("-", 72), "\n\n")

    # == 13c. Usage reminder =================================================
    log_cat("  HOW PART 10b PICKS THIS UP\n")
    log_cat("  Nothing to add: Part 10b already reads\n")
    log_cat("    ", P10B_IPRS_RDS, "\n")
    log_cat("  and sets DATASET_META$<ds>$IP_thresh from inflection_points_DFG\n")
    log_cat("  and DATASET_META$<ds>$IP_y_thresh from inflection_points_DFG_y.\n")
    log_cat("\n")
    log_cat("  IP values to be applied (full GOI-positive population):\n")
    for (ds in names(inflection_points_DFG)) {
        ipy_ds <- if (ds %in% names(inflection_points_DFG_y))
                      as.numeric(inflection_points_DFG_y[ds]) else NA_real_
        log_cat(sprintf("    %-10s  IP_thresh = %-18s  IP_y_thresh = %s\n", ds,
                        if (is.na(inflection_points_DFG[ds]))
                            "NA (IPRS failed)"
                        else sprintf("%.6f", inflection_points_DFG[ds]),
                        if (is.na(ipy_ds)) "NA (no IP_y)"
                        else sprintf("%.6f", ipy_ds)))
    }

    log_cat("\n", strrep("=", 72), "\n")
    log_cat("  SECTION 13 COMPLETE\n")
    log_cat("  IPRS RDS  : ", basename(iprs_rds_path), "\n")
    log_cat("  Excel     : ", OUT_PREFIX, "_results.xlsx\n")
    log_cat("  Summary   : ", OUT_PREFIX, "_summary.csv\n")
    log_cat("  Plots     : ", OUT_PREFIX, "_{dataset}_scatter.pdf/png\n")
    log_cat("              ", OUT_PREFIX, "_{dataset}_Rx.pdf/png\n")
    log_cat("              ", OUT_PREFIX, "_{dataset}_pitcher_Ki_expr.pdf/png\n")
    log_cat("              ", OUT_PREFIX, "_{dataset}_pitcher_Si_expr.pdf/png\n")
    log_cat(strrep("=", 72), "\n")

    # == 13d. Write log file =================================================
    VERIFY_LOG_PATH <- file.path(OUT_DIR, paste0(OUT_PREFIX, "_RDS_verify.log"))
    writeLines(LOG_LINES, con = VERIFY_LOG_PATH)
    cat("\n  Verification log saved:", VERIFY_LOG_PATH, "\n")


    # =========================================================================
    # SECTION 14: POST-ANALYSIS QUALITY GUARDRAILS
    # =========================================================================
    cat("\n", strrep("#", 72), "\n")
    cat("  SECTION 14 -- POST-ANALYSIS QUALITY GUARDRAILS\n")
    cat(strrep("#", 72), "\n\n")

    # == IP values finite where detected =====================================
    for (ds in DATASETS) {
        key <- paste0(ds, "|", GS_KEY)
        r   <- ip_results[[key]]
        if (is.null(r)) {
            griprs(paste0("GR-08-", ds), "FAIL",
                   paste0("no IP result for this dataset"))
            next
        }
        if (r$status == "ok" && !r$manual_flag) {
            griprs(paste0("GR-08-", ds),
                   if (is.finite(r$ip_value)) "PASS" else "FAIL",
                   sprintf("IP = %s  finite = %s",
                           if (is.na(r$ip_value)) "NA"
                           else sprintf("%.6f", r$ip_value),
                           if (is.finite(r$ip_value)) "TRUE" else "FALSE"))
        } else {
            griprs(paste0("GR-08-", ds), "WARN",
                   sprintf("IP detection status = %s (manual_flag = %s)",
                           r$status, r$manual_flag))
        }
    }

    # == IPRS values in [0, 1] ===============================================
    for (ds in DATASETS) {
        key <- paste0(ds, "|", GS_KEY)
        r   <- ip_results[[key]]
        if (is.null(r)) {
            griprs(paste0("GR-09-", ds), "FAIL",
                   "no IP result for this dataset, so there is no IPRS")
            next
        }
        iprs_val <- r$iprs$IPRS
        griprs(paste0("GR-09-", ds),
               if (!is.na(iprs_val) && iprs_val >= 0 && iprs_val <= 1) "PASS"
               else "FAIL",
               sprintf("IPRS = %.4f  tier = %s  (expected [0, 1])",
                       iprs_val, r$iprs$IPRS_tier))
    }

    # == IP_y finite and inside the observed cDFG range ======================
    for (ds in DATASETS) {
        gi <- gam_inputs[[ds]]
        # Index by name membership so a missing entry gives a zero-length value.
        ipy <- if (ds %in% names(inflection_points_DFG_y))
            suppressWarnings(as.numeric(inflection_points_DFG_y[[ds]]))
        else numeric(0)
        # A missing IP_y is recorded as a FAIL.
        if (is.null(gi) || length(ipy) == 0L) {
            griprs(paste0("GR-15-", ds), "FAIL",
                   sprintf(paste0("IP_y could not be checked for %s (%s)."),
                           ds,
                           if (is.null(gi)) "no GAM inputs for this dataset"
                           else "inflection_points_DFG_y has no entry"))
            next
        }
        yr <- range(gi$data$Expression, na.rm = TRUE)
        ok <- is.finite(ipy) && ipy >= yr[1] && ipy <= yr[2]
        griprs(paste0("GR-15-", ds), if (ok) "PASS" else "FAIL",
               sprintf(paste0("IP_y = %s | cDFG range [%.4f, %.4f] | ",
                              "classification threshold"),
                       if (is.finite(ipy)) sprintf("%.6f", ipy) else "NA",
                       yr[1], yr[2]))
    }

    # == Donor RE on the residual GAM ========================================
    for (ds in DATASETS) {
        key <- paste0(ds, "|", GS_KEY)
        r   <- ip_results[[key]]
        if (is.null(r)) {
            griprs(paste0("GR-12-", ds), "FAIL",
                   paste0("no IP result for this dataset, so no residual GAM"))
            next
        }
        if (!isTRUE(GAM_DONOR_RE) || !isTRUE(RESID_DONOR_RE)) {
            griprs(paste0("GR-12-", ds), "WARN",
                   sprintf("residual GAM fitted at cell level (%s)",
                           if (!isTRUE(GAM_DONOR_RE)) "GAM_DONOR_RE = FALSE"
                           else "RESID_DONOR_RE = FALSE"))
        } else {
            griprs(paste0("GR-12-", ds),
                   if (isTRUE(r$resid_donor_re)) "PASS" else "WARN",
                   sprintf("residual GAM R(x) donor RE = %s | donors = %s",
                           if (isTRUE(r$resid_donor_re)) "YES"
                           else sprintf("NO (below %d)", DONOR_RE_MIN_LEVELS),
                           if (is.null(r$resid_n_donors)) "NA"
                           else format(r$resid_n_donors)))
        }
    }

    # == RDS readback verification ===========================================
    rds_verify_ok <- !overall_fail
    griprs("GR-10",
           if (rds_verify_ok) "PASS" else "FAIL",
           sprintf("RDS readback verification: %d/%d checks passed",
                   sum(sapply(verify_results, function(v) v$status == "PASS")),
                   length(verify_results)))

    # == Write final guardrail log (overwrites pre-analysis version) =========
    griprs_write()

    cat(sprintf("\n  Final guardrail tally: %d PASS / %d WARN / %d FAIL\n",
                .griprs_pass, .griprs_warn, .griprs_fail))
    if (.griprs_fail > 0) {
        cat("  Guardrail failures detected\n")
    } else {
        cat(sprintf("  All %d guardrails passed\n",
                    .griprs_pass + .griprs_warn))
    }


    # =========================================================================
    # SECTION 15: Final summary
    # =========================================================================
    cat("\n", strrep("=", 72), "\n")
    cat("  Part10b IPRS COMPLETE -- [", ARM_TAG, "]\n")
    cat("  Output directory:", OUT_DIR, "\n")
    cat("  Files generated:\n")
    # The marker and the temporary message file are not listed.
    .transient <- c("_RUN_INCOMPLETE.txt", basename(MSG_PATH_IPRS))
    for (f in setdiff(sort(list.files(OUT_DIR)), .transient)) {
        sz <- round(file.size(file.path(OUT_DIR, f)) / 1e3, 0)
        cat(sprintf("    %-55s  %5d KB\n", f, sz))
    }
    cat("\n  IP values written for Part 10b:\n")
    for (nm in names(inflection_points_DFG))
        cat(sprintf("    %-10s  IP_thresh = %s\n", nm,
                    if (is.na(inflection_points_DFG[nm])) "NA (no crossing)"
                    else sprintf("%.6f", inflection_points_DFG[nm])))
    cat(strrep("=", 72), "\n")

    cat("\n  Log written to:", LOG_PATH_IPRS,
        if (isTRUE(LOG_MESSAGES_IPRS)) " (console output and warnings)"
        else " (console output only)", "\n")
    cat("  Part 10b IPRS all done:", format(Sys.time()), "\n\n")

    .iprs_completed <- TRUE

    if (sink.number(type = "message") != 2L) sink(type = "message")
    while (sink.number() > 0) sink()

    # Run complete: remove the incomplete-run marker.
    suppressWarnings(file.remove(INCOMPLETE_MARKER_IPRS))

    invisible(list(arm  = arm_tag, dir = OUT_DIR,
                   pass = .griprs_pass, warn = .griprs_warn,
                   fail = .griprs_fail))
}


# =============================================================================
# SECTION Z: SELF-TESTS
# =============================================================================
# Synthetic data built in memory; no file is read or written.

# == Synthetic frame =========================================================
# n cells over n_donors unbalanced donors; Expression is a planted logistic
# in GOI. donor_sd adds a per-donor offset (0 removes it).
iprs_make_synth <- function(n = 500L, n_donors = 20L,
                            min_cells = 10L, max_cells = 50L,
                            cross_q = 0.55, donor_sd = 0.4, seed = 42L) {
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

    GOI      <- rlnorm(n, 0.8, 0.6)
    x_cross  <- as.numeric(quantile(GOI, cross_q))
    Expression <- 3 / (1 + exp(-1.5 * (GOI - x_cross))) + rnorm(n, 0, 0.3)
    if (donor_sd > 0)
        Expression <- Expression +
            rnorm(n_donors, 0, donor_sd)[as.integer(donor)]

    list(data = data.frame(barcode    = paste0("SMOKE_", seq_len(n)),
                           GOI        = GOI,
                           Expression = Expression,
                           donor      = donor,
                           stringsAsFactors = FALSE),
         x_cross = x_cross,
         donor_counts = cnt)
}

# == Smoke test ==============================================================
iprs_smoke_test <- function(ctx) {

    g <- function(nm) get(nm, envir = ctx)
    fit_full_pop_gam   <- g("fit_full_pop_gam")
    de_explained       <- g("de_explained")
    predict_population <- g("predict_population")
    prep_fit_frame     <- g("prep_fit_frame")
    detect_ip          <- g("detect_ip")
    donor_re_ok        <- g("donor_re_ok")
    # Size the purple set as the run does.
    pop_dev_explained_model <- g("pop_dev_explained_model")

    cat("\n", strrep("=", 72), "\n", sep = "")
    cat("  PART 10b IPRS SMOKE TEST (synthetic, no external files)\n")
    cat(strrep("=", 72), "\n", sep = "")

    pass <- 0L; fail <- 0L
    report <- function(ok, label, detail = "") {
        if (isTRUE(ok)) { pass <<- pass + 1L; tag <- "PASS" }
        else            { fail <<- fail + 1L; tag <- "FAIL" }
        cat(sprintf("  [%s] %-22s %s\n", tag, label, detail))
        isTRUE(ok)
    }

    syn <- iprs_make_synth()
    dat <- prep_fit_frame(syn$data)
    cat(sprintf("  synthetic: %d cells | %d donors | planted crossing at GOI = %.4f\n",
                nrow(dat), nlevels(dat$donor), syn$x_cross))

    fit <- fit_full_pop_gam(dat, label = "SMOKE", verbose = FALSE)
    if (is.null(fit)) {
        cat("  [FAIL] full-population GAM returned NULL, remaining checks skipped\n")
        cat("SMOKE TEST: 0/5 PASS -- 5 FAILED\n")
        return(invisible(FALSE))
    }
    m  <- fit$model
    fd <- fit$data

    # == 1. Formula and DE estimator =========================================
    f_txt   <- paste(deparse(stats::formula(m)), collapse = " ")
    ok1a    <- grepl("s(donor", f_txt, fixed = TRUE)
    de_re   <- de_explained(m, fd)
    de_flat <- summary(m)$dev.expl
    ok1b    <- is.finite(de_re) && abs(de_re - de_flat) > 1e-6
    report(ok1a && ok1b, "1 formula + DE",
           sprintf("s(donor) at %d donors: %s | DE %.4f vs summary()$dev.expl %.4f",
                   fit$n_donors, ok1a, de_re, de_flat))

    # == 2. Weights ==========================================================
    per_donor <- tapply(fd$w, fd$donor, sum)
    ok2 <- isTRUE(all.equal(sum(fd$w),
                            as.numeric(nlevels(droplevels(fd$donor))))) &&
           isTRUE(all.equal(as.numeric(per_donor), rep(1, length(per_donor)))) &&
           isTRUE(all.equal(as.numeric(m$prior.weights), as.numeric(fd$w)))
    report(ok2, "2 weights",
           sprintf("sum(w)=%.4f over %d donors; per-donor totals all 1; prior.weights match",
                   sum(fd$w), nlevels(droplevels(fd$donor))))

    # == 3. Population vs conditional prediction =============================
    # Checked over every donor level.
    x_seq <- seq(min(fd$GOI), max(fd$GOI), length.out = 200)
    p_pop <- as.numeric(predict_population(m, x_seq, fd, se = FALSE))
    d_lev <- vapply(levels(fd$donor), function(lv) {
        nd <- data.frame(GOI = x_seq,
                         donor = factor(lv, levels = levels(fd$donor)))
        max(abs(p_pop - as.numeric(predict(m, newdata = nd))))
    }, numeric(1))
    report(max(d_lev) > 0.05, "3 prediction",
           sprintf("max |population - donor-conditional| = %.4f over %d donors",
                   max(d_lev), length(d_lev)))

    # == 4. IP detection =====================================================
    # Sized by the population-level DE.
    de_pop_sm <- pop_dev_explained_model(m, fd)
    if (!is.finite(de_pop_sm)) de_pop_sm <- de_flat
    cat(sprintf(paste0("  [ -- ] %-22s population %.4f | summary()$dev.expl ",
                       "%.4f | donor-aware %.4f\n"),
                "sizing DE", de_pop_sm, de_flat, de_re))
    res <- tryCatch(detect_ip(fd, m, de_pop_sm, dataset_tag = "SMOKE"),
                    error = function(e) {
                        cat("    detect_ip error:", conditionMessage(e), "\n")
                        NULL
                    })
    goi_rng <- diff(range(fd$GOI))
    ok4 <- !is.null(res) && is.finite(res$ip_value) &&
           abs(res$ip_value - syn$x_cross) < 0.25 * goi_rng
    report(ok4, "4 IP detection",
           if (is.null(res)) "detect_ip returned NULL"
           else sprintf("IP = %s | planted = %.4f | |diff| = %s of GOI range",
                        if (is.na(res$ip_value)) "NA"
                        else sprintf("%.4f", res$ip_value),
                        syn$x_cross,
                        if (is.na(res$ip_value)) "NA"
                        else sprintf("%.1f%%",
                                     100 * abs(res$ip_value - syn$x_cross) / goi_rng)))

    # == 5. Residual GAM =====================================================
    ok5 <- FALSE
    if (!is.null(res) && nrow(res$residual_data) > 0) {
        rd     <- prep_fit_frame(res$residual_data)
        use_re <- donor_re_ok(rd)
        k_r    <- min(20L, length(unique(rd$GOI)) - 1L)
        rg     <- mgcv::gam(get("make_resid_formula", envir = ctx)(k_r, use_re),
                            data = rd, weights = w, method = "REML")
        f_rg   <- paste(deparse(stats::formula(rg)), collapse = " ")
        xs     <- seq(min(rd$GOI), max(rd$GOI), length.out = 100)
        pp     <- as.numeric(predict_population(rg, xs, rd, se = FALSE))
        pc     <- as.numeric(predict(rg, newdata = data.frame(
                      GOI = xs, donor = rd$donor[which.max(abs(
                          tapply(rd$residual, rd$donor, mean)[rd$donor]))])))
        ok5 <- grepl("s(donor", f_rg, fixed = TRUE) && isTRUE(res$resid_donor_re) &&
               max(abs(pp - pc)) > 1e-8
        report(ok5, "5 residual GAM",
               sprintf("formula has s(donor): %s | population vs conditional differ by %.4f",
                       grepl("s(donor", f_rg, fixed = TRUE), max(abs(pp - pc))))
    } else {
        report(FALSE, "5 residual GAM", "no residual data to refit")
    }

    # Diagnostic: the same IP with the residual donor RE switched.
    old_rre <- get("RESID_DONOR_RE", envir = ctx)
    assign("RESID_DONOR_RE", !isTRUE(old_rre), envir = ctx)
    alt <- tryCatch(detect_ip(fd, m, de_pop_sm, dataset_tag = "SMOKE_ALT"),
                    error = function(e) NULL)
    assign("RESID_DONOR_RE", old_rre, envir = ctx)
    cat(sprintf(paste0("  [ -- ] %-22s RESID_DONOR_RE=%s gives IP %s; ",
                       "=%s gives IP %s (planted %.4f)\n"),
                "resid RE sensitivity", old_rre,
                if (is.null(res) || is.na(res$ip_value)) "NA"
                else sprintf("%.4f", res$ip_value),
                !isTRUE(old_rre),
                if (is.null(alt) || is.na(alt$ip_value)) "NA"
                else sprintf("%.4f", alt$ip_value),
                syn$x_cross))

    if (fail == 0L) cat("SMOKE TEST: 5/5 PASS\n")
    else cat(sprintf("SMOKE TEST: %d/5 PASS -- %d FAILED\n", pass, fail))
    invisible(fail == 0L)
}

# == Stress test =============================================================
# Sized to the PVM3_all pooled frame: 8,407 cells over 84 donors.
iprs_stress_test <- function(ctx) {
    g <- function(nm) get(nm, envir = ctx)
    fit_full_pop_gam <- g("fit_full_pop_gam")
    de_explained     <- g("de_explained")
    prep_fit_frame   <- g("prep_fit_frame")
    detect_ip        <- g("detect_ip")
    pop_dev_explained_model <- g("pop_dev_explained_model")

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
    don     <- factor(rep(sprintf("S%02d", seq_len(nd)), cnt))
    GOI     <- rlnorm(n, 0.8, 0.6)
    x_cross <- as.numeric(quantile(GOI, 0.55))
    dat <- prep_fit_frame(data.frame(
        barcode    = paste0("STRESS_", seq_len(n)),
        GOI        = GOI,
        Expression = 3 / (1 + exp(-1.5 * (GOI - x_cross))) +
                     rnorm(n, 0, 0.3) + rnorm(nd, 0, 0.4)[as.integer(don)],
        donor      = don,
        stringsAsFactors = FALSE))

    cat("\n", strrep("=", 72), "\n", sep = "")
    cat(sprintf("  PART 10b IPRS STRESS TEST -- %d cells, %d donors\n", n, nd))
    cat(strrep("=", 72), "\n", sep = "")

    t0  <- proc.time()[["elapsed"]]
    fit <- try(fit_full_pop_gam(dat, label = "STRESS", verbose = FALSE),
               silent = TRUE)
    t_gam <- proc.time()[["elapsed"]] - t0

    ok_gam <- !inherits(fit, "try-error") && !is.null(fit)
    de     <- if (ok_gam) de_explained(fit$model, fit$data) else NA_real_
    # Size the purple set as the run does.
    de_pop_st <- if (ok_gam) pop_dev_explained_model(fit$model, fit$data)
                 else NA_real_
    if (ok_gam && !is.finite(de_pop_st))
        de_pop_st <- summary(fit$model)$dev.expl

    t1  <- proc.time()[["elapsed"]]
    res <- if (ok_gam)
        try(detect_ip(fit$data, fit$model, de_pop_st, dataset_tag = "STRESS"),
            silent = TRUE) else NULL
    t_ip <- proc.time()[["elapsed"]] - t1

    ok_ip <- !is.null(res) && !inherits(res, "try-error")
    ok_w  <- ok_gam && isTRUE(all.equal(sum(fit$data$w), as.numeric(nd)))

    checks <- c(gam_completes = ok_gam, ip_runs = ok_ip, sum_w = ok_w)
    cat(sprintf("  wall time: full-population GAM %.1f s | IP detection %.1f s | total %.1f s\n",
                t_gam, t_ip, t_gam + t_ip))
    cat(sprintf("  DE = %s | sum(w) = %s (donors %d) | IP = %s\n",
                if (is.na(de)) "NA" else sprintf("%.6f", de),
                if (ok_gam) sprintf("%.4f", sum(fit$data$w)) else "NA", nd,
                if (!ok_ip) "NA" else
                    if (is.na(res$ip_value)) "NA (no crossing)"
                    else sprintf("%.4f", res$ip_value)))
    if (all(checks))
        cat(sprintf("STRESS TEST: 3/3 PASS  (GAM %.1f s, IP %.1f s)\n",
                    t_gam, t_ip))
    else
        cat(sprintf("STRESS TEST: %d/3 PASS -- failed: %s  (GAM %.1f s, IP %.1f s)\n",
                    sum(checks), paste(names(checks)[!checks], collapse = ", "),
                    t_gam, t_ip))
    invisible(all(checks))
}


# =============================================================================
# SECTION 0c: Driver
# =============================================================================
# The loop and the summary table live in CEP_arms.R.

cep_run_arms(PART_TAG, run_one_arm)
