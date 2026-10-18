# =============================================================================
# Part10h_CP10K_Prognosis_4_Multivariate_code.r
# SEA-AD CEP-IP/CMRF, Part 10h, multivariate (Figure 10C and 10D, Table S12)
#   Asks whether the CMRF module adds to a gene score plus a cell count. The
#   baselines enter as covariates: a nested F (1 df) of
#   y ~ cov + score + cluster proportion + module against the same model
#   without the module, per dataset x block x outcome. AUCell is the primary
#   score arm; UCell and AddModuleScore are substituted as robustness arms.
#
# Usage (from the SEA-AD folder, or set SEAAD_ROOT), after Part 10h Prognosis_1
#   source("Part10h_CP10K_Prognosis_4_Multivariate_code.r")
#
# Outputs (Part10h_Prognosis_4_Multivariate/<arm>/): workbook, results RDS,
#   figures (PDF, PNG), guardrail log and console log.
# =============================================================================

# =============================================================================
# SECTION 0: Libraries and shared arm definitions
# =============================================================================

suppressPackageStartupMessages({
    library(writexl)
    library(ggplot2)
})

set.seed(42)

if (!exists("CEP_HOME")) CEP_HOME <- Sys.getenv("SEAAD_ROOT", getwd())
.p10mv_harness <- file.path(CEP_HOME, "CEP_arms.R")
if (!file.exists(.p10mv_harness))
    stop("CEP_arms.R not found at: ", .p10mv_harness, call. = FALSE)
source(.p10mv_harness)

PART_TAG     <- "Part10h"
PART_DIR_TAG <- "Part10h_Prognosis_4_Multivariate"
CMP_DIR_TAG  <- "Part10h_Prognosis_1_Comparator"

MIN_DONORS <- 20L
ALPHA      <- 0.05

MV_OUTCOMES <- c("MMSE score", "CASI score")
MV_PATH_ADD <- c("Braak stage", "Thal phase", "CERAD score")   # the +BTC block

# The predictors and what each one is.
MV_MODULE <- "module"
MV_ABUND  <- "B1"
# Score arms, substituted one at a time; the first is the primary arm.
MV_SCORES <- c("B2", "B4", "B5")
MV_SCORE  <- MV_SCORES[1]
MV_LABEL  <- c(module = "CMRF module Post-IP (state x abundance)",
               B1 = "Cluster proportion (supertype) -- pure abundance",
               B2 = "Module score (AUCell) -- pure state",
               B4 = "Module score (UCell) -- pure state",
               B5 = "Module score (AddModuleScore) -- pure state")
MV_SHORT  <- c(module = "Module", B1 = "Cluster %", B2 = "AUCell",
               B4 = "UCell", B5 = "AddModuleScore")
MV_FILL   <- c(Module = "#2D2D2D", `Cluster %` = "#B2182B",
               AUCell = "#1B3A6B", UCell = "#1B7837",
               AddModuleScore = "#E08214")

# Above this joint VIF the coefficients are not separately readable.
MV_VIF_MAX <- 5

# MMSE and CASI are higher-is-better, so the expected module sign is negative.
MV_EXPECT <- c("MMSE score" = -1, "CASI score" = -1)

# Every outcome needs a declared sign.
local({
    .no_sign <- setdiff(MV_OUTCOMES, names(MV_EXPECT))
    if (length(.no_sign) > 0)
        stop("MV_EXPECT declares no expected sign for: ",
             paste(.no_sign, collapse = ", "), call. = FALSE)
    if (!all(MV_EXPECT[MV_OUTCOMES] %in% c(-1, 1)))
        stop("MV_EXPECT values must be -1 or +1.", call. = FALSE)
    .all_pred <- c(MV_MODULE, MV_ABUND, MV_SCORES)
    if (!all(.all_pred %in% names(MV_LABEL)) ||
        !all(.all_pred %in% names(MV_SHORT)))
        stop("MV_LABEL and MV_SHORT must name every predictor and score.",
             call. = FALSE)
    if (!all(unname(MV_SHORT[.all_pred]) %in% names(MV_FILL)))
        stop("MV_FILL has no colour for one of the predictors.", call. = FALSE)
    if (anyDuplicated(MV_SCORES) > 0L || length(MV_SCORES) == 0L)
        stop("MV_SCORES must be a non-empty set of distinct keys.",
             call. = FALSE)
    invisible(NULL)
})

# =============================================================================
# Pre-specification, declared before any model is fitted
# =============================================================================
MV_PRESPEC <- c(
    "QUESTION. Is the CMRF module percentage more than a gene score plus a",
    "cell count -- that is, does it survive when AUCell and the supertype",
    "cluster proportion are both already in the model?",
    "",
    "PRIMARY TEST. Nested F, 1 df, of",
    "    y ~ covariates + AUCell + cluster proportion + module",
    "against",
    "    y ~ covariates + AUCell + cluster proportion",
    "fitted on one common complete-case set per dataset x block x outcome.",
    "",
    "FAMILY. Benjamini-Hochberg within dataset x block, over the outcomes.",
    "",
    "VERDICT RULE, fixed before any model was fitted:",
    "  HOLDS       BH q < 0.05 AND the module coefficient carries its",
    "              expected sign (negative for MMSE and CASI).",
    "  WRONG SIGN  BH q < 0.05 but the sign is inverted -- a suppression",
    "              artefact, NOT support.",
    "  DOES NOT HOLD   BH q >= 0.05 with joint VIF below 5.",
    "  NOT READABLE    joint VIF at or above 5 on any of the three terms.",
    "",
    "INTERPRETATION, fixed in advance:",
    "  If the module HOLDS on a majority of readable rows, it is not",
    "  decomposable into its two rivals and no rearrangement of a gene-set",
    "  score reproduces it.",
    "  If it DOES NOT HOLD, the module is a score times an abundance. That is",
    "  a clean negative result and is to be reported as the answer, not",
    "  retested with a different predictor.",
    "",
    "SCOPE. The module is unchanged and the question is the one Part 10h has",
    "asked throughout: does the module offer anything the baselines do not.",
    "The baselines are entered as covariates, as +BTC enters pathology. Every",
    "model fitted is reported, not only those that favour the module.",
    "",
    "SCORE ARMS. The same specification is run three times with the",
    "score term SWAPPED, never stacked: AUCell, UCell, AddModuleScore. They",
    "correlate 0.986-0.999, so two of them in one model is the same variable",
    "twice. AUCell is the PRIMARY arm and carries the verdict above; the other",
    "two are ROBUSTNESS arms whose only job is to say whether the verdict",
    "depends on which score is used. A disagreement is reported, not resolved.")


run_one_arm <- function(arm_tag, config_file, out_dir, file_tag) {

    config_path <- file.path(CEP_HOME, config_file)
    if (!file.exists(config_path))
        stop("Configuration not found: ", config_path, call. = FALSE)

    .p10mv_preamble <- character(0)
    assign("CEP_ARM", arm_tag, envir = globalenv())
    on.exit(suppressWarnings(rm("CEP_ARM", envir = globalenv())), add = TRUE)

    .p10mv_preamble <- tryCatch(
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

    # One RDS, linear models only: small config surface.
    .p10mv_needed <- c("ARM_TAG", "GOI", "DFG", "NORM_TAG", "PRIMARY_STYPE")
    .p10mv_absent <- .p10mv_needed[!vapply(.p10mv_needed, exists, logical(1))]
    if (length(.p10mv_absent) > 0)
        stop("The configuration loaded from\n    ", config_path,
             "\n  does not define: ", paste(.p10mv_absent, collapse = ", "),
             call. = FALSE)
    if (!identical(ARM_TAG, arm_tag))
        stop("The configuration reports ARM_TAG '", ARM_TAG,
             "' but this arm is '", arm_tag, "'.", call. = FALSE)

    .p10mv_preamble <- c(.p10mv_preamble, utils::capture.output(
        cat("  Config OK: ", ARM_TAG, " | ", GOI, " + ", length(DFG),
            " DFG | ", NORM_TAG, "\n", sep = "")))

    # == Paths ================================================================
    OUT_DIR <- out_dir
    dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
    FILE_PREFIX <- cep_out_prefix(PART_TAG, NORM_TAG, cep_file_tag_for(arm_tag))
    CMP_RDS <- file.path(cep_arm_dir_for(CMP_DIR_TAG, arm_tag),
                         paste0(FILE_PREFIX, "_Comparator_results.rds"))

    # == Logging ==============================================================
    LOG_PATH_10MV <- file.path(OUT_DIR, paste0(FILE_PREFIX, "_console_log_",
                                               format(Sys.time(),
                                                      "%Y%m%d_%H%M%S"), ".txt"))
    MSG_PATH_10MV <- file.path(OUT_DIR, paste0(FILE_PREFIX,
                                               "_console_log_messages.tmp"))
    while (sink.number() > 0) sink()
    if (sink.number(type = "message") != 2L) sink(type = "message")
    .p10mv_log_con <- file(LOG_PATH_10MV, open = "wt")
    .p10mv_msg_con <- file(MSG_PATH_10MV, open = "wt")
    sink(.p10mv_log_con, append = FALSE, split = TRUE)
    sink(.p10mv_msg_con, type = "message")
    on.exit({
        tryCatch({
            if (!isTRUE(.p10mv_completed))
                cat("\n", strrep("!", 72),
                    "\n  Arm did not complete; the error is printed to the console.\n",
                    strrep("!", 72), "\n", sep = "")
            if (sink.number(type = "message") != 2L) sink(type = "message")
            while (sink.number() > 0) sink()
            close(.p10mv_msg_con); close(.p10mv_log_con)
            if (file.exists(MSG_PATH_10MV)) {
                .msg <- readLines(MSG_PATH_10MV, warn = FALSE)
                if (length(.msg) > 0)
                    cat(c("", strrep("=", 72),
                          "  MESSAGE STREAM (warnings, messages, errors)",
                          strrep("=", 72), .msg, ""),
                        sep = "\n", file = LOG_PATH_10MV, append = TRUE)
                unlink(MSG_PATH_10MV)
            }
        }, error = function(e) invisible(NULL))
    }, add = TRUE)

    .p10mv_completed <- FALSE
    .p10mv_warn_opt  <- getOption("warn")
    .p10mv_pm_opt    <- getOption("warnPartialMatchDollar")
    options(warn = 1)
    options(warnPartialMatchDollar = TRUE)
    on.exit({
        options(warn = .p10mv_warn_opt)
        options(warnPartialMatchDollar = .p10mv_pm_opt)
    }, add = TRUE)

    cat(sprintf("Part10h multivariate - log started: %s\n\n", Sys.time()))
    if (length(.p10mv_preamble) > 0)
        cat(paste(.p10mv_preamble, collapse = "\n"), "\n\n", sep = "")

    cat("\n", strrep("=", 72), "\n")
    cat("  Part10h multivariate: does the module add to the baselines?\n")
    cat("  Arm          :", ARM_TAG, "\n")
    cat("  GOI          :", GOI, "|", NORM_TAG, "\n")
    cat("  Primary input:", basename(CMP_RDS), "\n")
    cat("  Unit         : donor\n")
    cat("  Output dir   :", OUT_DIR, "\n")
    cat(strrep("=", 72), "\n\n")

    INCOMPLETE_MARKER_10MV <- file.path(OUT_DIR, "_RUN_INCOMPLETE.txt")
    writeLines(c("This arm did not finish.", "",
                 paste("Arm     :", arm_tag),
                 paste("Part    :", PART_DIR_TAG),
                 paste("Started :", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), "",
                 "Any files in this folder are from a partial run.",
                 "This file is removed automatically when the arm completes."),
               INCOMPLETE_MARKER_10MV)

    # =========================================================================
    # GUARDRAIL LEDGER
    # =========================================================================
    GR_LOG_10MV <- file.path(OUT_DIR, paste0(FILE_PREFIX, "_Guardrail_Log.txt"))
    .gr_log <- list(); .gr_pass <- 0L; .gr_fail <- 0L; .gr_warn <- 0L
    gr10mv <- function(id, status, detail) {
        # A zero-length or NA detail still prints a line.
        if (length(detail) == 0L || (length(detail) == 1L && is.na(detail)))
            detail <- paste("[detail unavailable: the value passed was NA",
                            "or zero-length]")
        detail <- paste(as.character(detail), collapse = " ")
        .gr_log[[length(.gr_log) + 1L]] <<- list(id = id, status = status,
                                                 detail = detail)
        if (status == "PASS") .gr_pass <<- .gr_pass + 1L
        if (status == "FAIL") .gr_fail <<- .gr_fail + 1L
        if (status == "WARN") .gr_warn <<- .gr_warn + 1L
        cat(sprintf("  %s  %-18s  %s\n",
                    switch(status, PASS = "[PASS]", FAIL = "[FAIL]",
                           WARN = "[WARN]"), id, detail))
    }

    cat("PRE-SPECIFICATION\n")
    cat(paste0("  ", MV_PRESPEC, collapse = "\n"), "\n\n", sep = "")

    # =========================================================================
    # STEP 1: the comparator's donor frame
    # =========================================================================
    cat("STEP 1 -- comparator donor frame\n")
    if (!file.exists(CMP_RDS))
        stop("The comparator's results RDS was not found:\n    ", CMP_RDS,
             call. = FALSE)
    RCMP <- readRDS(CMP_RDS)
    DF   <- RCMP[["donor_frame"]]
    MCMP <- RCMP[["meta"]]
    if (!is.data.frame(DF))
        stop("`donor_frame` is absent or not a data frame.",
             call. = FALSE)

    .req <- c(MV_MODULE, MV_ABUND, MV_SCORE)
    .miss <- .req[!(.req %in% names(DF))]
    if (length(.miss) > 0)
        stop("The comparator's donor frame does not carry: ",
             paste(vapply(.miss, function(k) MV_LABEL[[k]], character(1)),
                   collapse = ", "),
             call. = FALSE)
    # A robustness arm whose score is absent is dropped and named.
    SCORES <- MV_SCORES[MV_SCORES %in% names(DF)]
    .score_missing <- setdiff(MV_SCORES, SCORES)

    # The covariate and outcome columns travel with the comparator's meta.
    as_map <- function(x) {
        if (is.null(x)) return(NULL)
        if (is.list(x) && all(vapply(x, function(v)
            is.character(v) && length(v) == 1L, logical(1)))) x <- unlist(x)
        if (!is.character(x) || is.null(names(x))) return(NULL)
        x <- x[nzchar(names(x))]
        if (sum(x %in% names(DF)) < 2L) return(NULL)
        x
    }
    BUILTIN_OUT_COL <- c("MMSE score" = "MMSE_num", "CASI score" = "CASI_num",
                         "Braak stage" = "Braak_num", "Thal phase" = "Thal_num",
                         "CERAD score" = "CERAD_num")
    .out_map <- as_map(MCMP[["outcomes"]])
    OUT_ROUTE <- "comparator outcome map"
    if (is.null(.out_map)) { .out_map <- BUILTIN_OUT_COL
                             OUT_ROUTE <- "built-in" }
    col_of <- function(labs) { v <- unname(.out_map[labs]); names(v) <- labs
                               v[!is.na(v) & v %in% names(DF)] }
    OUT <- col_of(MV_OUTCOMES)
    BTC <- col_of(MV_PATH_ADD)

    BASE_COVS <- MCMP[["covariate_columns"]]
    if (is.null(BASE_COVS) || !is.character(BASE_COVS))
        BASE_COVS <- c("Age_num", "APOE_e4_dose", "Sex_num", "PMI_num")
    BASE_COVS <- BASE_COVS[BASE_COVS %in% names(DF)]
    if (length(BASE_COVS) == 0L || length(OUT) == 0L)
        stop("The base covariates or the outcomes did not resolve to columns ",
             "of the comparator's donor frame.", call. = FALSE)
    BLOCKS <- list(cov4 = character(0), `+BTC` = unname(BTC))

    DF$dataset <- as.character(DF[["dataset"]])
    DATASETS <- unique(DF$dataset[!is.na(DF$dataset)])

    gr10mv("GR-10MV-01", if (length(BTC) < 3L || length(.score_missing) > 0L)
                             "WARN" else "PASS",
           sprintf(paste("%s (schema %s) | %d row(s) over %s | base = %s |",
                         "outcomes %s | +BTC = %s | outcome map via %s |",
                         "score arm(s): %s.%s"),
                   basename(CMP_RDS),
                   if (is.null(MCMP[["schema"]])) "not recorded"
                   else MCMP[["schema"]],
                   nrow(DF), paste(DATASETS, collapse = " + "),
                   paste(BASE_COVS, collapse = " + "),
                   paste(names(OUT), collapse = ", "),
                   if (length(BTC) == 0L) "not available, cov4 only"
                   else paste(names(BTC), collapse = " + "),
                   OUT_ROUTE,
                   paste(sprintf("%s%s", vapply(SCORES, function(k)
                       MV_SHORT[[k]], character(1)),
                       ifelse(SCORES == MV_SCORE, " (primary)", "")),
                       collapse = ", "),
                   if (length(.score_missing) == 0L) "" else
                       sprintf(paste(" Arm(s) not built: %s (absent from the",
                                     "comparator donor frame)."),
                               paste(vapply(.score_missing, function(k)
                                   MV_SHORT[[k]], character(1)),
                                   collapse = ", "))))

    # =========================================================================
    # STEP 2: correlations among the three predictors
    # =========================================================================
    cat("\nSTEP 2 -- how the predictors relate\n")
    cor_rows <- list()
    for (ds in DATASETS) {
        d <- DF[DF$dataset == ds, , drop = FALSE]
        pr <- c(MV_MODULE, MV_ABUND, SCORES)
        for (i in seq_along(pr)) for (j in seq_along(pr)) {
            if (j <= i) next
            a <- suppressWarnings(as.numeric(d[[pr[i]]]))
            b <- suppressWarnings(as.numeric(d[[pr[j]]]))
            ok <- is.finite(a) & is.finite(b)
            if (sum(ok) < MIN_DONORS) next
            cor_rows[[length(cor_rows) + 1L]] <- data.frame(
                dataset = ds, x = MV_LABEL[[pr[i]]], y = MV_LABEL[[pr[j]]],
                n = sum(ok),
                spearman = round(suppressWarnings(
                    stats::cor(a[ok], b[ok], method = "spearman")), 3),
                pearson = round(suppressWarnings(
                    stats::cor(a[ok], b[ok])), 3),
                stringsAsFactors = FALSE)
        }
    }
    COR_DF <- if (length(cor_rows)) do.call(rbind, cor_rows) else NULL
    if (!is.null(COR_DF)) {
        rownames(COR_DF) <- NULL
        for (k in seq_len(nrow(COR_DF)))
            cat(sprintf("  %-8s %-42s vs %-42s Spearman %+.3f (n = %d)\n",
                        COR_DF$dataset[k], COR_DF$x[k], COR_DF$y[k],
                        COR_DF$spearman[k], COR_DF$n[k]))
    }

    # =========================================================================
    # STEP 3: six models per row
    # =========================================================================
    cat(sprintf(paste("\nSTEP 3 -- models per dataset x block x outcome,",
                      "%d score arm(s)\n"), length(SCORES)))

    .zs <- function(x) {
        sdv <- sd(x, na.rm = TRUE)
        if (!is.finite(sdv) || sdv < 1e-12) return(x - mean(x, na.rm = TRUE))
        (x - mean(x, na.rm = TRUE)) / sdv
    }
    bq <- function(v) if (length(v)) paste(sprintf("`%s`", v), collapse = " + ")
                      else ""

    # Specifications declared as data so the sheet and the fits stay in step.
    # Those with the score term are refitted once per arm, the rest once per row.
    SPEC_FREE <- list(
        list(id = "M0", terms = character(0), label = "covariates only"),
        list(id = "M1", terms = MV_MODULE,    label = "+ module"),
        list(id = "M3", terms = MV_ABUND,     label = "+ cluster proportion"))
    # Built per arm so the label names the score fitted.
    spec_score <- function(sk) list(
        list(id = "M2", terms = sk,
             label = sprintf("+ %s", MV_SHORT[[sk]])),
        list(id = "M4", terms = c(sk, MV_ABUND),
             label = sprintf("+ %s + cluster proportion  (THE RIVAL)",
                             MV_SHORT[[sk]])),
        list(id = "M5", terms = c(sk, MV_ABUND, MV_MODULE),
             label = sprintf("+ %s + cluster proportion + module  (THE TEST)",
                             MV_SHORT[[sk]])),
        # Descriptive only.
        list(id = "M6", terms = c(sk, MV_MODULE),
             label = sprintf("+ %s + module  (DESCRIPTIVE)", MV_SHORT[[sk]])))

    fit_spec <- function(sp, covs, D) {
        rhs <- c(sp$terms, covs)
        tryCatch(stats::lm(stats::as.formula(
            if (length(rhs)) paste("y ~", bq(rhs)) else "y ~ 1"),
            data = D), error = function(e) NULL)
    }

    mod_rows <- list(); dec_rows <- list(); .n_fit <- 0L; .n_rowset <- 0L
    for (ds in DATASETS) for (bn in names(BLOCKS)) for (lab in names(OUT)) {
        ycol <- OUT[[lab]]
        covs <- setdiff(c(BASE_COVS, BLOCKS[[bn]]), ycol)
        covs <- covs[covs %in% names(DF)]
        # One complete-case set spanning every score arm, so all arms share donors.
        need <- unique(c(ycol, covs, MV_MODULE, MV_ABUND, SCORES))
        d0 <- DF[DF$dataset == ds, , drop = FALSE]
        if (!all(need %in% names(d0))) next
        D <- data.frame(y = suppressWarnings(as.numeric(d0[[ycol]])))
        for (v in unique(c(covs, MV_MODULE, MV_ABUND, SCORES)))
            D[[v]] <- suppressWarnings(as.numeric(d0[[v]]))
        D <- D[stats::complete.cases(D), , drop = FALSE]
        if (nrow(D) < MIN_DONORS) next
        for (cn in names(D)) D[[cn]] <- .zs(D[[cn]])
        .n_rowset <- .n_rowset + 1L

        add_mod <- function(sp, f, sk) {
            if (is.null(f)) return(invisible(NULL))
            .n_fit <<- .n_fit + 1L
            s2 <- summary(f); cf <- s2$coefficients
            emit <- function(tm) {
                rn <- paste0("`", tm, "`")
                if (!(rn %in% rownames(cf))) rn <- tm
                if (!(rn %in% rownames(cf))) return(NULL)
                data.frame(dataset = ds, block = bn, outcome = lab,
                           score_arm = sk, model = sp$id,
                           model_terms = sp$label,
                           predictor = MV_LABEL[[tm]], predictor_key = tm,
                           n = nrow(D),
                           beta = round(unname(cf[rn, 1]), 4),
                           se = round(unname(cf[rn, 2]), 4),
                           p_value = unname(cf[rn, 4]),
                           r2 = round(s2$r.squared, 4),
                           stringsAsFactors = FALSE)
            }
            if (length(sp$terms) == 0L)
                mod_rows[[length(mod_rows) + 1L]] <<- data.frame(
                    dataset = ds, block = bn, outcome = lab, score_arm = sk,
                    model = sp$id, model_terms = sp$label,
                    predictor = "(covariates only)", predictor_key = "",
                    n = nrow(D), beta = NA_real_, se = NA_real_,
                    p_value = NA_real_, r2 = round(s2$r.squared, 4),
                    stringsAsFactors = FALSE)
            else for (tm in sp$terms) {
                r <- emit(tm)
                if (!is.null(r)) mod_rows[[length(mod_rows) + 1L]] <<- r
            }
            invisible(NULL)
        }

        # == Score-free models, fitted once per row =========================
        free <- list()
        for (sp in SPEC_FREE) {
            f <- fit_spec(sp, covs, D)
            free[[sp$id]] <- f
            add_mod(sp, f, "")
        }
        if (is.null(free[["M0"]])) next

        # == One specification set per score arm ============================
        for (sk in SCORES) {
            sps <- spec_score(sk)
            fits <- list()
            for (sp in sps) {
                f <- fit_spec(sp, covs, D)
                fits[[sp$id]] <- f
                add_mod(sp, f, MV_SHORT[[sk]])
            }
            if (is.null(fits[["M4"]]) || is.null(fits[["M5"]])) next

            # Verdict: nested F of M5 against M4.
            av <- tryCatch(stats::anova(fits[["M4"]], fits[["M5"]]),
                           error = function(e) NULL)
            f_stat <- if (is.null(av)) NA_real_ else unname(av[["F"]][2])
            f_p    <- if (is.null(av)) NA_real_ else unname(av[["Pr(>F)"]][2])

            s5 <- summary(fits[["M5"]]); c5 <- s5$coefficients
            getc <- function(cf, tm, k) {
                if (is.null(cf)) return(NA_real_)
                rn <- paste0("`", tm, "`")
                if (!(rn %in% rownames(cf))) rn <- tm
                if (!(rn %in% rownames(cf))) return(NA_real_)
                unname(cf[rn, k])
            }
            c6 <- if (is.null(fits[["M6"]])) NULL else
                summary(fits[["M6"]])$coefficients

            # VIF by its definition: 1 / (1 - R2 of that term on the others).
            vif_on <- function(tm, oth) {
                f <- tryCatch(stats::lm(stats::as.formula(
                    paste0("`", tm, "` ~ ", bq(c(oth, covs)))), data = D),
                    error = function(e) NULL)
                if (is.null(f)) return(NA_real_)
                r2 <- summary(f)$r.squared
                if (!is.finite(r2) || r2 >= 1) return(Inf)
                1 / (1 - r2)
            }
            .three <- c(sk, MV_ABUND, MV_MODULE)
            v_mod <- vif_on(MV_MODULE, setdiff(.three, MV_MODULE))
            v_sc  <- vif_on(sk,        setdiff(.three, sk))
            v_ab  <- vif_on(MV_ABUND,  setdiff(.three, MV_ABUND))
            # M6 collinearity: module against the score alone.
            v6 <- max(vif_on(MV_MODULE, sk), vif_on(sk, MV_MODULE))

            dec_rows[[length(dec_rows) + 1L]] <- data.frame(
                dataset = ds, block = bn, outcome = lab,
                score_arm = MV_SHORT[[sk]], score_key = sk,
                is_primary = identical(sk, MV_SCORE), n = nrow(D),
                r2_covariates = round(summary(free[["M0"]])$r.squared, 4),
                r2_rival = round(summary(fits[["M4"]])$r.squared, 4),
                r2_full = round(s5$r.squared, 4),
                dR2_module_over_rival = round(s5$r.squared -
                    summary(fits[["M4"]])$r.squared, 4),
                F_1df = round(f_stat, 3), p_F = f_p,
                beta_module = round(getc(c5, MV_MODULE, 1), 4),
                p_module_t = getc(c5, MV_MODULE, 4),
                beta_score = round(getc(c5, sk, 1), 4),
                p_score = getc(c5, sk, 4),
                beta_cluster = round(getc(c5, MV_ABUND, 1), 4),
                p_cluster = getc(c5, MV_ABUND, 4),
                # Descriptive: module against the score alone.
                r2_module_plus_score = if (is.null(fits[["M6"]])) NA_real_
                    else round(summary(fits[["M6"]])$r.squared, 4),
                beta_module_vs_score = round(getc(c6, MV_MODULE, 1), 4),
                p_module_vs_score = getc(c6, MV_MODULE, 4),
                vif_module_vs_score = round(v6, 2),
                vif_module = round(v_mod, 2), vif_score = round(v_sc, 2),
                vif_cluster = round(v_ab, 2),
                expected_sign = unname(MV_EXPECT[lab]),
                stringsAsFactors = FALSE)
        }
    }

    MOD <- if (length(mod_rows)) do.call(rbind, mod_rows) else NULL
    DEC <- if (length(dec_rows)) do.call(rbind, dec_rows) else NULL
    if (!is.null(MOD)) rownames(MOD) <- NULL
    if (is.null(DEC))
        stop("No row produced both the rival and the full model, so the ",
             "pre-specified test could not be run on this arm.", call. = FALSE)
    rownames(DEC) <- NULL

    # =========================================================================
    # STEP 4: BH, readability and the verdict
    # =========================================================================
    cat("\nSTEP 4 -- the pre-specified verdict\n")

    # BH within dataset x block x score arm, over the outcomes.
    .fam <- paste(DEC$dataset, DEC$block, DEC$score_key, sep = "|")
    DEC$q_F <- NA_real_; DEC$q_module_t <- NA_real_
    DEC$q_module_vs_score <- NA_real_
    for (k in unique(.fam)) {
        i <- .fam == k
        DEC$q_F[i] <- stats::p.adjust(DEC$p_F[i], method = "BH")
        DEC$q_module_t[i] <- stats::p.adjust(DEC$p_module_t[i], method = "BH")
        DEC$q_module_vs_score[i] <-
            stats::p.adjust(DEC$p_module_vs_score[i], method = "BH")
    }
    # `Inf` is perfect collinearity, so the test is !is.na(), not is.finite().
    .unread <- (!is.na(DEC$vif_module) & DEC$vif_module >= MV_VIF_MAX) |
               (!is.na(DEC$vif_score) & DEC$vif_score >= MV_VIF_MAX) |
               (!is.na(DEC$vif_cluster) & DEC$vif_cluster >= MV_VIF_MAX) |
               is.na(DEC$vif_module) | is.na(DEC$vif_score) |
               is.na(DEC$vif_cluster)
    DEC$max_vif <- pmax(DEC$vif_module, DEC$vif_score, DEC$vif_cluster,
                        na.rm = FALSE)
    .sig  <- !is.na(DEC$q_F) & DEC$q_F < ALPHA
    # A surviving coefficient with the wrong sign is a suppression artefact.
    .right <- !is.na(DEC$beta_module) &
              sign(DEC$beta_module) == sign(DEC$expected_sign)
    DEC$verdict <- ifelse(.unread, "not readable",
                   ifelse(.sig & .right, "HOLDS",
                   ifelse(.sig & !.right, "WRONG SIGN",
                          "does not hold")))
    DEC$readable <- !.unread

    # The verdict is the primary arm's; robustness arms are reported separately.
    PRI <- DEC[DEC$is_primary, , drop = FALSE]
    .rd <- PRI[PRI$readable, , drop = FALSE]
    .n_hold <- sum(.rd$verdict == "HOLDS")
    .n_wrong <- sum(.rd$verdict == "WRONG SIGN")
    for (k in seq_len(nrow(DEC)))
        cat(sprintf(paste0("  %-8s %-5s %-12s %-14s n=%3d | dR2 %+.4f",
                           " | F %6s q %s | VIF %4s | %s\n"),
                    DEC$dataset[k], DEC$block[k], DEC$outcome[k],
                    DEC$score_arm[k], DEC$n[k],
                    DEC$dR2_module_over_rival[k],
                    if (is.finite(DEC$F_1df[k])) sprintf("%.2f", DEC$F_1df[k])
                    else "n/a",
                    if (is.finite(DEC$q_F[k])) sprintf("%.4f", DEC$q_F[k])
                    else "n/a",
                    if (is.finite(DEC$max_vif[k]))
                        sprintf("%.1f", DEC$max_vif[k]) else "n/a",
                    DEC$verdict[k]))

    gr10mv("GR-10MV-02", if (nrow(DEC) == 0L) "FAIL" else "PASS",
           sprintf(paste("%d decomposition row(s) over %d dataset x block x",
                         "outcome set(s) and %d score arm(s), %d model fit(s).",
                         "Donor counts %s"),
                   nrow(DEC), .n_rowset, length(SCORES), .n_fit,
                   paste(range(DEC$n), collapse = "-")))

    gr10mv("GR-10MV-03",
           if (nrow(.rd) == 0L) "WARN"
           else if (.n_wrong > 0L) "WARN"
           else if (.n_hold > nrow(.rd) / 2) "PASS" else "WARN",
           sprintf(paste("Pre-specified verdict, %s arm only. Nested F,",
                         "1 df, of `cov + score + cluster + module` against",
                         "`cov + score + cluster`: %s. Of %d readable row(s)",
                         "of %d, HOLDS %d, does not hold %d, WRONG SIGN %d.",
                         "%s%s"),
                   MV_SHORT[[MV_SCORE]],
                   paste(sprintf("%s/%s/%s %s", PRI$dataset, PRI$block,
                                 PRI$outcome, PRI$verdict), collapse = "; "),
                   nrow(.rd), nrow(PRI), .n_hold,
                   sum(.rd$verdict == "does not hold"), .n_wrong,
                   if (nrow(.rd) == 0L)
                       paste("No row is readable on this arm; the question",
                             "is not answered.")
                   else if (.n_wrong > 0L)
                       paste("Rows with the wrong sign are not support.")
                   else if (.n_hold > nrow(.rd) / 2)
                       paste("The module is not decomposable into a gene score",
                             "plus a cell count on this arm.")
                   else if (.n_hold > 0L)
                       paste("Mixed; the pre-specified majority is not met.")
                   else
                       paste("The module is a score times an abundance on this",
                             "arm."),
                   if (nrow(.rd) < nrow(PRI))
                       sprintf(paste(" %d row(s) not readable at joint VIF",
                                     ">= %g."),
                               nrow(PRI) - nrow(.rd), MV_VIF_MAX) else ""))

    # range()/max() on an all-NA column return +/-Inf with a warning.
    # MV_LABEL and MV_SHORT share key order, so a label maps back by position.
    .short_of <- function(lab) unname(MV_SHORT[match(lab, MV_LABEL)])
    .rng <- function(v) if (all(is.na(v))) "n/a" else
        paste(range(v, na.rm = TRUE), collapse = "-")
    .mx <- function(v) if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE)
    gr10mv("GR-10MV-04", if (any(!DEC$readable)) "WARN" else "PASS",
           sprintf(paste("Collinearity. Joint VIF ranges:",
                         "module %s, score %s, cluster %s (refusal at %g).%s"),
                   .rng(DEC$vif_module), .rng(DEC$vif_score),
                   .rng(DEC$vif_cluster),
                   MV_VIF_MAX,
                   if (is.null(COR_DF)) "" else
                       sprintf(" Predictor correlations: %s",
                               paste(sprintf("%s %s/%s rho %+.2f",
                                             COR_DF$dataset,
                                             .short_of(COR_DF$x),
                                             .short_of(COR_DF$y),
                                             COR_DF$spearman),
                                     collapse = "; "))))

    # Descriptive: M6 is judged on its own collinearity.
    DEC$readable_vs_score <- !is.na(DEC$vif_module_vs_score) &
                             DEC$vif_module_vs_score < MV_VIF_MAX
    .m6 <- DEC[DEC$is_primary, , drop = FALSE]
    .m6_ok <- sum(.m6$readable_vs_score & !is.na(.m6$q_module_vs_score))
    .m6_hold <- sum(.m6$readable_vs_score & !is.na(.m6$q_module_vs_score) &
                    .m6$q_module_vs_score < ALPHA)
    gr10mv("GR-10MV-09", if (.m6_ok == 0L) "WARN" else "PASS",
           sprintf(paste("Descriptive (%s arm). M6, module plus score:",
                         "%d of %d row(s) readable, VIF %s. Module",
                         "significant with the score alone on %d; with",
                         "abundance also present, %d. Per row: %s"),
                   MV_SHORT[[MV_SCORE]], .m6_ok, nrow(.m6),
                   .rng(.m6$vif_module_vs_score), .m6_hold, .n_hold,
                   paste(sprintf("%s/%s/%s %s -> %s", .m6$dataset, .m6$block,
                                 .m6$outcome,
                                 ifelse(!.m6$readable_vs_score,
                                        "not readable",
                                 ifelse(!is.na(.m6$q_module_vs_score) &
                                        .m6$q_module_vs_score < ALPHA,
                                        "sig", "ns")),
                                 .m6$verdict), collapse = "; ")))

    # Does the verdict depend on which score is in the model?
    if (length(SCORES) < 2L)
        gr10mv("GR-10MV-10", "WARN",
               sprintf(paste("only the %s arm was built; score dependence",
                             "not tested"), MV_SHORT[[MV_SCORE]]))
    else {
        .vk <- paste(DEC$dataset, DEC$block, DEC$outcome, sep = "|")
        .agree_rows <- vapply(unique(.vk), function(k)
            length(unique(DEC$verdict[.vk == k])) == 1L, logical(1))
        .tal <- vapply(SCORES, function(sk) {
            z <- DEC[DEC$score_key == sk, , drop = FALSE]
            sprintf(paste("%s HOLDS %d / does not hold %d / wrong sign %d /",
                          "not readable %d"),
                    MV_SHORT[[sk]], sum(z$verdict == "HOLDS"),
                    sum(z$verdict == "does not hold"),
                    sum(z$verdict == "WRONG SIGN"),
                    sum(z$verdict == "not readable"))
        }, character(1))
        gr10mv("GR-10MV-10", if (all(.agree_rows)) "PASS" else "WARN",
               sprintf(paste("Score-arm robustness. %s. The scores agree on",
                             "%d of %d row(s). %s"),
                       paste(.tal, collapse = "; "),
                       sum(.agree_rows), length(.agree_rows),
                       if (all(.agree_rows))
                           paste("The verdict does not depend on the score.")
                       else paste("The verdict differs by score on",
                                  sum(!.agree_rows),
                                  "row(s).")))
    }

    # The nested F and the module's t-test in M5 are the same test; they must
    # agree. Rows missing a p on either side are counted separately.
    .cmp <- is.finite(DEC$p_F) & is.finite(DEC$p_module_t)
    .agree <- .cmp & abs(DEC$p_F - DEC$p_module_t) < 1e-6
    .dmax <- if (!any(.cmp)) NA_real_ else
        max(abs(DEC$p_F[.cmp] - DEC$p_module_t[.cmp]))
    gr10mv("GR-10MV-06",
           if (any(.cmp) && all(.agree[.cmp])) "PASS"
           else if (!any(.cmp)) "WARN" else "FAIL",
           sprintf(paste("the nested F (1 df) and the module's t-test in the",
                         "full model agree on %d of %d comparable row(s), max",
                         "|dp| = %.2e%s"),
                   sum(.agree), sum(.cmp), .dmax,
                   if (any(!.cmp))
                       sprintf(paste("; %d row(s) not comparable"), sum(!.cmp))
                   else ""))

    # =========================================================================
    # STEP 5: figure
    # =========================================================================
    cat("\nSTEP 5 -- figure\n")
    .n_fig <- 0L
    save_fig <- function(p, base, w = 10, h = 6, dpi = 300) {
        render <- function(open_fn, path) tryCatch({
            open_fn(); on.exit(grDevices::dev.off(), add = TRUE); print(p); TRUE
        }, error = function(e) {
            cat(sprintf("    not written (%s): %s\n", basename(path),
                        conditionMessage(e))); FALSE })
        pp <- paste0(base, ".pdf"); gp <- paste0(base, ".png")
        ok1 <- render(function() grDevices::pdf(pp, width = w, height = h), pp)
        ok2 <- render(function() grDevices::png(
            gp, width = w * dpi, height = h * dpi, res = dpi), gp)
        if (ok1 && ok2) .n_fig <<- .n_fig + 1L
        invisible(NULL)
    }
    theme_mv <- function() theme_bw(base_size = 11) +
        theme(panel.grid.minor = element_blank(),
              strip.background = element_rect(fill = "grey92", colour = NA),
              plot.title = element_text(face = "bold"),
              plot.subtitle = element_text(size = 9, colour = "grey25"))

    # All three coefficients from the full model, primary arm.
    FL <- do.call(rbind, lapply(
        list(c(MV_MODULE, "beta_module", "p_module_t", "vif_module"),
             c(MV_SCORE,  "beta_score",  "p_score",    "vif_score"),
             c(MV_ABUND,  "beta_cluster", "p_cluster", "vif_cluster")),
        function(z) data.frame(
            PRI[, c("dataset", "block", "outcome", "n", "verdict")],
            predictor = MV_SHORT[[z[1]]],
            beta = PRI[[z[2]]],
            p = PRI[[z[3]]], vif = PRI[[z[4]]],
            stringsAsFactors = FALSE)))
    # CI rebuilt from the Models sheet, which holds se for every M5 term.
    .m5 <- MOD[MOD$model == "M5" & MOD$score_arm == MV_SHORT[[MV_SCORE]], ,
               drop = FALSE]
    # FL$predictor already holds the SHORT label, so it keys directly.
    .key <- paste(FL$dataset, FL$block, FL$outcome, FL$predictor, sep = "|")
    .m5key <- paste(.m5$dataset, .m5$block, .m5$outcome,
                    MV_SHORT[.m5$predictor_key], sep = "|")
    FL$se <- .m5$se[match(.key, .m5key)]
    FL$lo <- FL$beta - 1.96 * FL$se
    FL$hi <- FL$beta + 1.96 * FL$se
    FL$sig <- !is.na(FL$p) & FL$p < ALPHA
    FL$pred_f <- factor(FL$predictor, levels = rev(unname(MV_SHORT)))
    FL$panel <- factor(paste(FL$dataset, FL$block, sep = "  "),
                       levels = unique(paste(rep(DATASETS,
                                                 each = length(BLOCKS)),
                                             names(BLOCKS), sep = "  ")))
    FL$panel <- droplevels(FL$panel)

    p1 <- ggplot(FL, aes(x = beta, y = pred_f, colour = predictor)) +
        geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
        geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.18,
                       linewidth = 0.6, na.rm = TRUE) +
        geom_point(aes(shape = sig), size = 2.8, na.rm = TRUE) +
        scale_shape_manual(values = c(`TRUE` = 19, `FALSE` = 1),
                           name = "p < 0.05", na.translate = FALSE) +
        scale_colour_manual(values = MV_FILL, guide = "none") +
        facet_grid(panel ~ outcome) +
        labs(title = sprintf("All three predictors in one model [%s, %s arm]",
                             ARM_TAG, MV_SHORT[[MV_SCORE]]),
             subtitle = paste(sprintf(paste("y ~ covariates + %s + cluster",
                                            "proportion + module,"),
                                      MV_SHORT[[MV_SCORE]]),
                              "standardised, one common donor set per panel.",
                              "\nThe question is whether the module keeps a",
                              "coefficient when a pure state score and a pure",
                              "abundance are both already present.",
                              "\nThe module is unchanged; the baselines are",
                              "entered as covariates."),
             x = "standardised coefficient in the full model (95% CI)",
             y = NULL) +
        theme_mv()
    print(p1)
    save_fig(p1, file.path(OUT_DIR, paste0(FILE_PREFIX,
                                           "_Multivariate_Coefficients")),
             w = 9, h = min(24, max(5, 2.1 * length(levels(FL$panel)))))

    # Verdict per row per score, drawn only with more than one score arm.
    if (length(SCORES) > 1L) {
        VG <- DEC
        VG$row_lab <- paste(VG$dataset, VG$block, VG$outcome, sep = " | ")
        VG$row_f <- factor(VG$row_lab, levels = rev(unique(VG$row_lab)))
        VG$arm_f <- factor(VG$score_arm,
                           levels = unname(MV_SHORT[SCORES]))
        VG$verdict_f <- factor(VG$verdict,
                               levels = c("HOLDS", "does not hold",
                                          "WRONG SIGN", "not readable"))
        VG$tile <- sprintf("%s\nq %s  VIF %s",
                           VG$verdict,
                           ifelse(is.finite(VG$q_F),
                                  sprintf("%.3f", VG$q_F), "n/a"),
                           ifelse(is.finite(VG$max_vif),
                                  sprintf("%.1f", VG$max_vif), "n/a"))
        p2 <- ggplot(VG, aes(x = arm_f, y = row_f, fill = verdict_f)) +
            geom_tile(colour = "white", linewidth = 0.6) +
            geom_text(aes(label = tile), size = 2.2, colour = "black",
                      lineheight = 0.95) +
            scale_fill_manual(values = c(HOLDS = "#BFE3BF",
                                         `does not hold` = "#E8E8E8",
                                         `WRONG SIGN` = "#F2C29B",
                                         `not readable` = "#CFCFCF"),
                              name = NULL, drop = FALSE) +
            labs(title = sprintf("Does the verdict depend on the score? [%s]",
                                 ARM_TAG),
                 subtitle = paste("The same model with the score term",
                                  "SWAPPED, never stacked: the three scores",
                                  "correlate 0.99 and two in one model is the",
                                  "same variable twice.",
                                  sprintf(paste("\n%s is the pre-specified",
                                                "arm and carries the verdict;",
                                                "the others are robustness."),
                                          MV_SHORT[[MV_SCORE]])),
                 x = NULL, y = NULL) +
            theme_mv() +
            theme(panel.grid = element_blank(),
                  axis.text.y = element_text(size = 7.5))
        print(p2)
        save_fig(p2, file.path(OUT_DIR, paste0(FILE_PREFIX,
                                               "_Multivariate_ScoreArms")),
                 w = 8, h = min(24, max(5, 0.75 * length(unique(VG$row_lab))
                                        + 3)))
    }

    # =========================================================================
    # STEP 6: output
    # =========================================================================
    cat("\nSTEP 6 -- output\n")
    wb <- list(Prespecification = data.frame(line = MV_PRESPEC,
                                             stringsAsFactors = FALSE),
               Decomposition = DEC)
    if (!is.null(MOD)) wb[["Models"]] <- MOD
    if (!is.null(COR_DF)) wb[["Correlations"]] <- COR_DF
    xlsx_path <- file.path(OUT_DIR, paste0(FILE_PREFIX, "_Multivariate.xlsx"))
    tryCatch({
        writexl::write_xlsx(wb, xlsx_path)
        cat(sprintf("  %d sheet(s) written to %s\n", length(wb),
                    basename(xlsx_path)))
    }, error = function(e)
        cat("  workbook not written:", conditionMessage(e), "\n"))

    rds_out <- list(decomposition = DEC, models = MOD, correlations = COR_DF,
                    prespecification = MV_PRESPEC,
                    meta = list(arm_tag = ARM_TAG, config_file = config_file,
                                norm_tag = NORM_TAG, GOI = GOI,
                                supertype = PRIMARY_STYPE,
                                comparator_rds = CMP_RDS,
                                base_covariates = BASE_COVS, blocks = BLOCKS,
                                outcomes = OUT, vif_max = MV_VIF_MAX,
                                written = format(Sys.time(),
                                                 "%Y-%m-%d %H:%M:%S"),
                                R_version = R.version.string,
                                schema = "M10MV v1"))
    rds_path <- file.path(OUT_DIR,
                          paste0(FILE_PREFIX, "_Multivariate_results.rds"))
    tryCatch({
        saveRDS(rds_out, rds_path)
        cat(sprintf("  %s written (%.2f MB)\n", basename(rds_path),
                    file.size(rds_path) / 1024^2))
    }, error = function(e)
        cat("  RDS not written:", conditionMessage(e), "\n"))

    # =========================================================================
    # STEP 7: remaining guardrails
    # =========================================================================
    cat("\nSTEP 7 -- guardrails\n")

    # Every fitted model must appear in the Models sheet.
    .fitted   <- .n_fit
    .in_sheet <- if (is.null(MOD)) 0L else
        nrow(unique(MOD[, c("dataset", "block", "outcome", "score_arm",
                            "model")]))
    # Nominal: score-free specifications once per row set, score ones per arm.
    .nominal <- .n_rowset * (length(SPEC_FREE) +
                             length(spec_score(MV_SCORE)) * length(SCORES))
    gr10mv("GR-10MV-05", if (.in_sheet >= .fitted) "PASS" else "FAIL",
           sprintf(paste("%d model(s) fitted, %d in the Models sheet",
                         "(nominal %d). %d specification(s): %s.%s"),
                   .fitted, .in_sheet, .nominal,
                   length(SPEC_FREE) + length(spec_score(MV_SCORE)),
                   paste(vapply(c(SPEC_FREE, spec_score(MV_SCORE)), function(z)
                       sprintf("%s %s", z$id, z$label), character(1)),
                       collapse = "; "),
                   if (.fitted < .nominal)
                       sprintf(paste(" %d specification-row(s) did not",
                                     "fit."), .nominal - .fitted) else ""))

    gr10mv("GR-10MV-07", if (.n_fig == 0L) "WARN" else "PASS",
           sprintf(paste("%d figure(s) written as PDF and PNG"), .n_fig))

    .files <- c(rds_path, xlsx_path); .ok <- file.exists(.files)
    gr10mv("GR-10MV-08", if (all(.ok)) "PASS" else "FAIL",
           sprintf("outputs: %s",
                   paste(sprintf("%s %s", basename(.files),
                                 ifelse(.ok, "written", "MISSING")),
                         collapse = " | ")))

    .gr_lines <- c(
        strrep("=", 72),
        "Part 10h multivariate - Guardrail Log [CP10K]",
        paste("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
        paste("Arm :", ARM_TAG),
        strrep("=", 72), "",
        "PRE-SPECIFICATION:",
        paste0("  ", MV_PRESPEC), "",
        strrep("=", 72), "",
        vapply(.gr_log, function(g)
            sprintf("  %s  %-18s  %s",
                    switch(g$status, PASS = "[PASS]", FAIL = "[FAIL]",
                           WARN = "[WARN]"), g$id, g$detail), character(1)),
        "", strrep("=", 72),
        sprintf("  Total : %d", length(.gr_log)),
        sprintf("  PASS  : %d", .gr_pass),
        sprintf("  WARN  : %d", .gr_warn),
        sprintf("  FAIL  : %d", .gr_fail), "",
        if (.gr_fail > 0) "  Result: failures present"
        else "  Result: no failures",
        strrep("=", 72))
    writeLines(.gr_lines, GR_LOG_10MV)
    cat(sprintf("\n  Final guardrail tally: %d PASS / %d WARN / %d FAIL\n",
                .gr_pass, .gr_warn, .gr_fail))

    if (file.exists(INCOMPLETE_MARKER_10MV)) unlink(INCOMPLETE_MARKER_10MV)
    .p10mv_completed <- TRUE
    cat("\n", strrep("=", 72), "\n")
    cat("  Part10h stage 4 multivariate complete [", ARM_TAG, "]\n", sep = "")
    cat("  Output directory:", OUT_DIR, "\n")
    cat(strrep("=", 72), "\n")

    invisible(NULL)
}

cep_run_arms(PART_DIR_TAG, run_one_arm)
