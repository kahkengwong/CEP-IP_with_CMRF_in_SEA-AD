# =============================================================================
# Part10h_CP10K_Prognosis_2_CogVsPath_code.r
# SEA-AD CEP-IP/CMRF, Part 10h, cognition against pathology (Figure 9B and 9C, Table S12)
#   Tests whether the CMRF module still predicts cognition (MMSE, CASI, MOCA,
#   cognitive status) once plaque and tangle burden are in the model: nested
#   blocks cov4, +ADNC, +BTC (Braak, Thal, CERAD) and +ALL, the attenuation of
#   each coefficient, the mirror (pathology adjusted for MMSE and CASI), and a
#   joint module plus non-module model. Baselines are read from the comparator.
#
# Usage (from the SEA-AD folder, or set SEAAD_ROOT), after Part 10h Prognosis
#   and the comparator
#   source("Part10h_CP10K_Prognosis_2_CogVsPath_code.r")
#
# Outputs (Part10h_Prognosis_2_CogVsPath/<arm>/): workbook, results RDS,
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
.p10cp_harness <- file.path(CEP_HOME, "CEP_arms.R")
if (!file.exists(.p10cp_harness))
    stop("CEP_arms.R not found at: ", .p10cp_harness,
         "\n  It holds the arm table shared by Parts 9 onwards.", call. = FALSE)
source(.p10cp_harness)

PART_TAG     <- "Part10h"
PART_DIR_TAG <- "Part10h_Prognosis_2_CogVsPath"
P10H_DIR_TAG <- "Part10h_Prognosis"

BASE_BLOCK <- "cov4"
MIN_DONORS <- 20L
# Below this |beta| floor pct_retained is NA and abs_change_beta is reported.
MIN_BETA_FOR_PCT <- 0.10
ALPHA      <- 0.05
VIF_READABLE <- 5

# Cognition, most to least objectively measured; reported separately, never
# pooled. MOCA is not among Part 10h's outcomes, so it has no cov4 cross-check.
COG_OUTCOMES  <- c("MMSE score", "CASI score", "MOCA score", "Cognitive status")
# The mirror: pathology outcomes.
PATH_OUTCOMES <- c("ADNC rank", "CERAD score", "Braak stage", "Thal phase")

# Pathology additions, by Part 10h outcome labels.
PATH_ADD <- list(
    "+ADNC" = "ADNC rank",
    "+BTC"  = c("Braak stage", "Thal phase", "CERAD score"),
    "+ALL"  = c("ADNC rank", "Braak stage", "Thal phase", "CERAD score"))
HEADLINE_BLOCK <- "+BTC"

# cov4 takes the emphatic colour and +ALL the muted one; order comes from breaks.
BLOCK_FILL <- c(cov4 = "#7B2D8E", `+ADNC` = "#5BA4CF",
                `+BTC` = "#1B3A6B", `+ALL` = "#B0B0B0")

# == Baselines, read from the comparator =====================================
# Donor-level columns only; no score is recomputed.
COMPARATOR_DIR_TAG <- "Part10h_Prognosis_1_Comparator"
BASE_KEYS <- c("B1", "B2", "B3", "B4", "B5")
# == The rest of the 2x2, from Part 10h prop_df =============================
# Keys are prop_df column names; labels are Part 10h PRED_LABEL strings.
MODFAM_KEYS <- c(pct_excl_post = "Non-module Post-IP",
                 pct_pre       = "Module Pre-IP",
                 pct_excl_pre  = "Non-module Pre-IP")

# The mirror runs on the cell populations only, not on the gene-score baselines.
MIRROR_KEYS <- c("module", names(MODFAM_KEYS))

# == Joint model =============================================================
JOINT_A      <- "module"          # must be a key of PRED_KEYS
JOINT_B      <- "pct_excl_post"   # must be a key of MODFAM_KEYS
JOINT_BLOCKS <- c("cov4", "+BTC")
# A row at or above this VIF gets no verdict.
JOINT_VIF_MAX <- 5

# Figure height cap (inches).
MAX_FIG_H <- 22

# == Publication panels ======================================================
# Fixed x limits, shared across arms so assembled panels are comparable.
# A guardrail fails if any interval falls outside them.
PUB_ATT_XLIM  <- c(-0.90, 0.30)   # cognition
PUB_MIR_XLIM  <- c(-0.40, 0.80)   # mirror, pathology
PUB_OUTCOMES  <- c("MMSE score", "CASI score")   # cognitive status is demoted
PUB_KEY       <- "module"

# Figure heights (inches) and dodge width for the publication panels.
PUB_ATT_H     <- 3.20
PUB_MIR_H     <- 3.80
PUB_DODGE     <- 0.65

PRED_KEYS <- c(module = "CMRF module Post-IP",
               MODFAM_KEYS,
               B1 = "Cluster proportion (supertype)",
               B2 = "Module score (AUCell)",
               B3 = "Changepoint split (segmented)",
               B4 = "Module score (UCell)",
               B5 = "Module score (AddModuleScore)")
local({
    if (!(JOINT_A %in% names(PRED_KEYS)) ||
        !(JOINT_B %in% names(MODFAM_KEYS)))
        stop("JOINT_A must be a PRED_KEYS key and JOINT_B a MODFAM_KEYS key.",
             call. = FALSE)
    if (anyDuplicated(names(PRED_KEYS)) > 0L)
        stop("PRED_KEYS has a duplicated key; a predictor would be fitted ",
             "twice and counted twice in every tally.", call. = FALSE)
})
PRED_FILL <- c("CMRF module Post-IP" = "#2D2D2D",
               # B5 already uses #E08214, so Non-module Post-IP differs from Part 10h.
               "Non-module Post-IP" = "#8C510A",
               "Module Pre-IP" = "#2166AC",
               "Non-module Pre-IP" = "#8073AC",
               "Cluster proportion (supertype)" = "#B2182B",
               "Module score (AUCell)" = "#1B3A6B",
               "Changepoint split (segmented)" = "#7B2D8E",
               "Module score (UCell)" = "#1B7837",
               "Module score (AddModuleScore)" = "#E08214")
# Short facet labels; an unmapped name would drop a facet silently, so it is checked.
PRED_SHORT <- c("CMRF module Post-IP" = "Module",
                "Non-module Post-IP" = "NonMod Post",
                "Module Pre-IP" = "Module Pre",
                "Non-module Pre-IP" = "NonMod Pre",
                "Cluster proportion (supertype)" = "Cluster %",
                "Module score (AUCell)" = "AUCell",
                "Changepoint split (segmented)" = "Changepoint",
                "Module score (UCell)" = "UCell",
                "Module score (AddModuleScore)" = "AddModuleScore")
local({
    if (!all(names(PRED_FILL) %in% unname(PRED_KEYS)) ||
        !setequal(names(PRED_SHORT), unname(PRED_KEYS)))
        stop("PRED_KEYS, PRED_FILL and PRED_SHORT are out of step.", call. = FALSE)
})


run_one_arm <- function(arm_tag, config_file, out_dir, file_tag) {

    config_path <- file.path(CEP_HOME, config_file)
    if (!file.exists(config_path))
        stop("Configuration not found: ", config_path, call. = FALSE)

    .p10cp_preamble <- character(0)
    assign("CEP_ARM", arm_tag, envir = globalenv())
    on.exit(suppressWarnings(rm("CEP_ARM", envir = globalenv())), add = TRUE)

    .p10cp_preamble <- tryCatch(
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

    # One RDS in, linear models out: no matrix, no obs table, no dataset dir.
    .p10cp_needed <- c("ARM_TAG", "GOI", "DFG", "NORM_TAG", "PRIMARY_STYPE")
    .p10cp_absent <- .p10cp_needed[!vapply(.p10cp_needed, exists, logical(1))]
    if (length(.p10cp_absent) > 0)
        stop("The configuration loaded from\n    ", config_path,
             "\n  does not define: ", paste(.p10cp_absent, collapse = ", "),
             call. = FALSE)
    if (!identical(ARM_TAG, arm_tag))
        stop("The configuration reports ARM_TAG '", ARM_TAG,
             "' but this arm is '", arm_tag, "'.", call. = FALSE)

    .p10cp_preamble <- c(.p10cp_preamble, utils::capture.output(
        cat("  Config OK: ", ARM_TAG, " | ", GOI, " + ", length(DFG),
            " DFG | ", NORM_TAG, "\n", sep = "")))

    DATASETS <- c("NeoC", "NonNeoC")

    # == Paths ===============================================================
    OUT_DIR <- out_dir
    dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
    FILE_PREFIX <- cep_out_prefix(PART_TAG, NORM_TAG, cep_file_tag_for(arm_tag))

    # Part 10h results object, written beside its workbook.
    P10H_DIR <- cep_arm_dir_for(P10H_DIR_TAG, arm_tag)
    P10H_RDS <- file.path(P10H_DIR,
                          paste0(FILE_PREFIX, "_Prognosis_results.rds"))

    # == Logging =============================================================
    LOG_PATH_10CP <- file.path(OUT_DIR, paste0(FILE_PREFIX, "_console_log_",
                                               format(Sys.time(),
                                                      "%Y%m%d_%H%M%S"), ".txt"))
    MSG_PATH_10CP <- file.path(OUT_DIR, paste0(FILE_PREFIX,
                                               "_console_log_messages.tmp"))
    while (sink.number() > 0) sink()
    if (sink.number(type = "message") != 2L) sink(type = "message")
    .p10cp_log_con <- file(LOG_PATH_10CP, open = "wt")
    .p10cp_msg_con <- file(MSG_PATH_10CP, open = "wt")
    sink(.p10cp_log_con, append = FALSE, split = TRUE)
    sink(.p10cp_msg_con, type = "message")
    on.exit({
        tryCatch({
            if (!isTRUE(.p10cp_completed))
                cat("\n", strrep("!", 72),
                    "\n  Arm did not complete; the error is printed to the console.\n",
                    strrep("!", 72), "\n", sep = "")
            if (sink.number(type = "message") != 2L) sink(type = "message")
            while (sink.number() > 0) sink()
            close(.p10cp_msg_con); close(.p10cp_log_con)
            if (file.exists(MSG_PATH_10CP)) {
                .msg <- readLines(MSG_PATH_10CP, warn = FALSE)
                if (length(.msg) > 0)
                    cat(c("", strrep("=", 72),
                          "  MESSAGE STREAM (warnings, messages, errors)",
                          strrep("=", 72), .msg, ""),
                        sep = "\n", file = LOG_PATH_10CP, append = TRUE)
                unlink(MSG_PATH_10CP)
            }
        }, error = function(e) invisible(NULL))
    }, add = TRUE)

    .p10cp_completed <- FALSE
    .p10cp_warn_opt  <- getOption("warn")
    .p10cp_pm_opt    <- getOption("warnPartialMatchDollar")
    options(warn = 1)
    options(warnPartialMatchDollar = TRUE)
    on.exit({
        options(warn = .p10cp_warn_opt)
        options(warnPartialMatchDollar = .p10cp_pm_opt)
    }, add = TRUE)

    cat(sprintf("Part10h stage 2 Cog-vs-Path - log started: %s\n\n",
                Sys.time()))
    if (length(.p10cp_preamble) > 0)
        cat(paste(.p10cp_preamble, collapse = "\n"), "\n\n", sep = "")

    cat("\n", strrep("=", 72), "\n")
    cat("  Part10h stage 2 - Cognition vs Pathology\n")
    cat("  Arm    :", ARM_TAG, "|", GOI, "|", NORM_TAG, "\n")
    cat("  Input  :", basename(P10H_RDS), "(+ the comparator RDS, optional)\n")
    cat("  Unit   : donor\n")
    cat("  Output :", OUT_DIR, "\n")
    cat(strrep("=", 72), "\n\n")

    INCOMPLETE_MARKER_10CP <- file.path(OUT_DIR, "_RUN_INCOMPLETE.txt")
    writeLines(c("This arm did not finish.", "",
                 paste("Arm     :", arm_tag),
                 paste("Part    :", PART_DIR_TAG),
                 paste("Started :", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), "",
                 "Any files in this folder are from a partial run.",
                 "This file is removed automatically when the arm completes."),
               INCOMPLETE_MARKER_10CP)

    # =============================================================================
    # GUARDRAIL LEDGER
    # =============================================================================
    GR_LOG_10CP <- file.path(OUT_DIR, paste0(FILE_PREFIX, "_Guardrail_Log.txt"))
    .gr_log <- list(); .gr_pass <- 0L; .gr_fail <- 0L; .gr_warn <- 0L
    gr10cp <- function(id, status, detail) {
        # A zero-length detail would vanish from console and ledger.
        if (length(detail) != 1L || is.na(detail))
            detail <- "[detail unavailable]"
        tag <- switch(status, PASS = "[PASS]", FAIL = "[FAIL]", WARN = "[WARN]")
        cat(sprintf("  %s  %-18s  %s\n", tag, id, detail))
        .gr_log[[length(.gr_log) + 1L]] <<- list(id = id, status = status,
                                                 detail = detail)
        if (identical(status, "PASS")) .gr_pass <<- .gr_pass + 1L
        else if (identical(status, "WARN")) .gr_warn <<- .gr_warn + 1L
        else .gr_fail <<- .gr_fail + 1L
        invisible(NULL)
    }


    # =============================================================================
    # HELPERS
    # =============================================================================
    .z <- function(x) {
        x <- suppressWarnings(as.numeric(x))
        s <- stats::sd(x, na.rm = TRUE)
        if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
        (x - mean(x, na.rm = TRUE)) / s
    }
    .stars <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "***",
                         ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", ""))))
    # max/min that return NA on an empty or all-NA vector.
    .rng <- function(x, f) {
        x <- x[is.finite(x)]
        if (length(x) == 0L) NA_real_ else f(x)
    }
    .safe_vif <- function(fit) {
        X <- try(stats::model.matrix(fit)[, -1, drop = FALSE], silent = TRUE)
        if (inherits(X, "try-error") || ncol(X) < 2L) return(c(x = NA_real_))
        v <- vapply(seq_len(ncol(X)), function(j) {
            r2 <- suppressWarnings(summary(stats::lm(X[, j] ~ X[, -j,
                                                                drop = FALSE]))$r.squared)
            if (!is.finite(r2) || r2 >= 1) NA_real_ else 1 / (1 - r2)
        }, numeric(1))
        stats::setNames(v, colnames(X))
    }

    # =============================================================================
    # STEP 1: Part 10h's results object
    # =============================================================================
    cat("STEP 1 -- Part 10h results object\n")
    if (!file.exists(P10H_RDS))
        stop("Part 10h's results RDS was not found:\n    ", P10H_RDS,
             call. = FALSE)
    R10H   <- readRDS(P10H_RDS)
    DECL   <- R10H[["declarations"]]
    DSETS  <- R10H[["datasets"]]
    SHEETS <- R10H[["sheets"]]
    M10H   <- R10H[["meta"]]
    if (is.null(DECL) || is.null(DSETS))
        stop("Part 10h's RDS carries no `declarations` or no `datasets`.",
             call. = FALSE)

    # == Donor frame =========================================================
    DONOR <- list()
    for (ds in names(DSETS)) {
        pd <- DSETS[[ds]][["prop_df"]]
        if (!is.data.frame(pd) || is.null(pd[["donor_id"]]) ||
            is.null(pd[["pct_post"]])) next
        pd$.donor <- as.character(pd[["donor_id"]])
        pd$dataset <- ds
        pd$module <- suppressWarnings(as.numeric(pd[["pct_post"]]))
        DONOR[[ds]] <- pd
    }
    if (length(DONOR) == 0L)
        stop("No dataset carries a `prop_df` with `donor_id` and `pct_post`.",
             call. = FALSE)
    .common <- Reduce(intersect, lapply(DONOR, names))
    DF <- do.call(rbind, lapply(DONOR, function(d) d[, .common, drop = FALSE]))
    rownames(DF) <- NULL

    # == Key to column and label to column ===================================
    # Part 10h declarations first; built-ins are fallbacks.
    BUILTIN_COV_COL <- c(age = "Age_num", apoe = "APOE_e4_dose", sex = "Sex_num",
                         pmi = "PMI_num", rin = "RIN_num", ph = "BrainPH_num",
                         weight = "BrainWeight_num", mito = "MitoFrac_num",
                         numi = "nUMI_num", ngene = "nGene_num",
                         intronexon = "IntronExon_num", doublet = "Doublet_num",
                         satur = "Satur_num", late = "LATE_num",
                         lewyfull = "LewyFull_num", athero = "Athero_num",
                         arteriolo = "Arteriolo_num", microinf = "Microinf_num",
                         eduyears = "EduYears_num", mmseint = "MMSEint_num")
    BUILTIN_OUT_COL <- c(
        "Cognitive status" = "CogStatus_num", "MMSE score" = "MMSE_num",
        "CASI score" = "CASI_num", "MOCA score" = "MOCA_num",
    "CERAD score" = "CERAD_num",
        "Braak stage" = "Braak_num", "Thal phase" = "Thal_num",
        "ADNC rank" = "ADNC_rank", "CAA score" = "CAA_num",
        "APOE e4 dose" = "APOE_e4_dose", "Age at death" = "Age_num",
        "Sex" = "Sex_num", "Brain pH" = "BrainPH_num", "PMI" = "PMI_num",
        "RIN" = "RIN_num", "Fresh brain weight" = "BrainWeight_num")
    as_map <- function(x) {
        if (is.null(x)) return(NULL)
        if (is.list(x) && all(vapply(x, function(v)
            is.character(v) && length(v) == 1L, logical(1)))) x <- unlist(x)
        if (!is.character(x) || is.null(names(x))) return(NULL)
        x <- x[nzchar(names(x))]
        if (sum(x %in% names(DF)) < 2L) return(NULL)
        x
    }
    .cov_map <- as_map(DECL[["COV_MAP"]]); COV_ROUTE <- "Part 10h COV_MAP"
    if (is.null(.cov_map)) { .cov_map <- BUILTIN_COV_COL; COV_ROUTE <- "built-in" }
    .out_map <- as_map(DECL[["CLINICAL_PARAMS"]]); OUT_ROUTE <- "Part 10h CLINICAL_PARAMS"
    if (is.null(.out_map)) { .out_map <- BUILTIN_OUT_COL; OUT_ROUTE <- "built-in" }

    COV_SETS <- DECL[["COV_SETS"]]
    if (is.null(COV_SETS[[BASE_BLOCK]]))
        stop("Part 10h's COV_SETS has no '", BASE_BLOCK, "' block.", call. = FALSE)
    BASE_COVS <- unname(.cov_map[COV_SETS[[BASE_BLOCK]]])
    BASE_COVS <- BASE_COVS[!is.na(BASE_COVS) & BASE_COVS %in% names(DF)]
    col_of <- function(labs) {
        v <- unname(.out_map[labs]); names(v) <- labs
        v[!is.na(v) & v %in% names(DF)]
    }
    COG  <- col_of(COG_OUTCOMES)
    PATH <- col_of(PATH_OUTCOMES)
    PATH_COLS <- lapply(PATH_ADD, col_of)

    # A covariate counts only if its mapped column exists in the donor frame.
    .cov_ok <- COV_SETS[[BASE_BLOCK]][
        !is.na(.cov_map[COV_SETS[[BASE_BLOCK]]]) &
        .cov_map[COV_SETS[[BASE_BLOCK]]] %in% names(DF)]
    .lost <- c(setdiff(COG_OUTCOMES, names(COG)),
               setdiff(PATH_OUTCOMES, names(PATH)),
               setdiff(COV_SETS[[BASE_BLOCK]], .cov_ok))
    if (length(BASE_COVS) == 0L || length(COG) == 0L)
        stop("The base block or the cognitive outcomes did not resolve to columns ",
             "of the donor frame. Unresolved: ", paste(.lost, collapse = ", "),
             "\n  Columns available: ", paste(utils::head(names(DF), 60),
                                              collapse = ", "), call. = FALSE)

    gr10cp("GR-10CP-01", if (length(.lost) > 0L) "WARN" else "PASS",
           sprintf(paste("%s (written %s) | %d row(s), %d dataset(s) |",
                         "%s = %s | blocks %s | %d cognitive, %d pathology",
                         "outcome(s) | maps via %s, %s | donors usable: %s%s"),
                   basename(P10H_RDS),
                   if (is.null(M10H[["written"]])) "not recorded"
                   else M10H[["written"]],
                   nrow(DF), length(DONOR), BASE_BLOCK,
                   paste(BASE_COVS, collapse = " + "),
                   paste(sprintf("%s(+%d)", names(PATH_COLS),
                                 vapply(PATH_COLS, length, integer(1))),
                         collapse = " "),
                   length(COG), length(PATH), COV_ROUTE, OUT_ROUTE,
                   paste(sprintf("%s %d", names(COG),
                                 vapply(unname(COG), function(c1)
                                     sum(is.finite(suppressWarnings(
                                         as.numeric(DF[[c1]])))),
                                     integer(1))), collapse = ", "),
                   if (length(.lost) > 0L)
                       paste0(" | not resolved: ",
                              paste(unique(.lost), collapse = ", ")) else ""))

    # =============================================================================
    # STEP 1b: baselines, read from the comparator
    # =============================================================================
    PREDS <- "module"
    CMP_RDS <- file.path(cep_arm_dir_for(COMPARATOR_DIR_TAG, arm_tag),
                         paste0(FILE_PREFIX, "_Comparator_results.rds"))
    .cmp_note <- ""
    if (!file.exists(CMP_RDS))
        .cmp_note <- paste0("not found: ", basename(CMP_RDS))
    else {
        RC <- tryCatch(readRDS(CMP_RDS), error = function(e) NULL)
        cd <- if (is.null(RC)) NULL else RC[["donor_frame"]]
        if (is.null(cd) || !all(c(".donor", "dataset", "module") %in% names(cd)))
            .cmp_note <- "the comparator RDS carries no usable `donor_frame`"
        else {
            have <- intersect(BASE_KEYS, names(cd))
            if (length(have) == 0L)
                .cmp_note <- "no baseline column in the comparator donor frame"
            else {
                DF <- merge(DF, cd[, c(".donor", "dataset", "module", have)],
                            by = c(".donor", "dataset"), all.x = TRUE,
                            suffixes = c("", ".cmp"))
                # The join is verified on the module column both frames carry.
                .jd <- suppressWarnings(max(abs(DF$module - DF$module.cmp),
                                            na.rm = TRUE))
                DF$module.cmp <- NULL
                if (!is.finite(.jd) || .jd > 1e-8)
                    .cmp_note <- sprintf(
                        "the join does not reproduce the module column (max |dev| %s)",
                        if (is.finite(.jd)) sprintf("%.3e", .jd) else "NA")
                else PREDS <- c("module", have)
            }
        }
    }
    # == The rest of the 2x2 =================================================
    # Columns of prop_df; an absent, non-numeric or constant column is dropped.
    .mf_have <- character(0); .mf_drop <- character(0)
    for (k in names(MODFAM_KEYS)) {
        if (!(k %in% names(DF))) {
            .mf_drop <- c(.mf_drop, paste0(k, " (absent from prop_df)"))
            next
        }
        v <- suppressWarnings(as.numeric(DF[[k]]))
        if (all(is.na(v))) {
            .mf_drop <- c(.mf_drop, paste0(k, " (all NA)"))
        } else if (isTRUE(stats::sd(v, na.rm = TRUE) == 0)) {
            .mf_drop <- c(.mf_drop, paste0(k, " (constant)"))
        } else {
            DF[[k]] <- v
            .mf_have <- c(.mf_have, k)
        }
    }
    # Order sets the facet order: the module family sits next to the module.
    PREDS <- c("module", .mf_have, setdiff(PREDS, "module"))
    gr10cp("GR-10CP-11", if (length(.mf_have) == 0L) "WARN"
                         else if (length(.mf_drop) > 0L) "WARN" else "PASS",
           sprintf("module family: %d of %d present (%s)%s",
                   length(.mf_have), length(MODFAM_KEYS),
                   if (length(.mf_have))
                       paste(vapply(.mf_have, function(k) MODFAM_KEYS[[k]],
                                    character(1)), collapse = ", ")
                   else "none",
                   if (length(.mf_drop))
                       paste0(" | dropped: ", paste(.mf_drop, collapse = ", "))
                   else ""))

    gr10cp("GR-10CP-08", if (length(PREDS) == 1L) "WARN" else "PASS",
           if (length(PREDS) == 1L)
               sprintf("baselines not added (%s)", .cmp_note)
           else sprintf(paste("baselines from %s, joined on donor + dataset;",
                              "join reproduces the module exactly (max |dev|",
                              "%.1e). Predictors: %s"),
                        basename(CMP_RDS), .jd,
                        paste(vapply(PREDS, function(k) PRED_KEYS[[k]],
                                     character(1)), collapse = ", ")))

    # =============================================================================
    # STEP 1c: how alike are the cognitive measures?
    # =============================================================================
    # Is cognitive status close to a coarsened MMSE?
    cor_rows <- list()
    for (ds in unique(DF$dataset)) {
        d <- DF[DF$dataset == ds, , drop = FALSE]
        for (i in seq_along(COG)) for (j in seq_along(COG)) if (i < j) {
            x <- suppressWarnings(as.numeric(d[[COG[[i]]]]))
            y <- suppressWarnings(as.numeric(d[[COG[[j]]]]))
            ok <- is.finite(x) & is.finite(y)
            cor_rows[[length(cor_rows) + 1L]] <- data.frame(
                dataset = ds, a = names(COG)[i], b = names(COG)[j],
                n = sum(ok),
                n_levels_a = length(unique(x[is.finite(x)])),
                n_levels_b = length(unique(y[is.finite(y)])),
                pearson = if (sum(ok) >= MIN_DONORS) round(suppressWarnings(
                    stats::cor(x[ok], y[ok])), 3) else NA_real_,
                spearman = if (sum(ok) >= MIN_DONORS) round(suppressWarnings(
                    stats::cor(x[ok], y[ok], method = "spearman")), 3) else NA_real_,
                kendall = if (sum(ok) >= MIN_DONORS) round(suppressWarnings(
                    stats::cor(x[ok], y[ok], method = "kendall")), 3) else NA_real_,
                stringsAsFactors = FALSE)
        }
    }
    COR <- do.call(rbind, cor_rows)
    if (!is.null(COR)) {
        rownames(COR) <- NULL
        cat("\nSTEP 1c -- agreement among the cognitive measures\n")
        print(COR)
        .cs <- COR[COR$a == "Cognitive status" | COR$b == "Cognitive status", ,
                   drop = FALSE]
        gr10cp("GR-10CP-09",
               if (nrow(.cs) == 0L) "WARN" else "PASS",
               if (nrow(.cs) == 0L)
                   "cognitive status is not among the resolved outcomes"
               else sprintf(paste("cognitive status vs the administered",
                                  "scores: %s. %d level(s) against MMSE's %d"),
                            paste(sprintf("%s/%s rho %.3f (n %d)",
                                          .cs$dataset,
                                          ifelse(.cs$a == "Cognitive status",
                                                 .cs$b, .cs$a),
                                          .cs$spearman, .cs$n),
                                  collapse = "; "),
                            max(c(.cs$n_levels_a[.cs$a == "Cognitive status"],
                                  .cs$n_levels_b[.cs$b == "Cognitive status"])),
                            max(c(.cs$n_levels_a[.cs$a == "MMSE score"],
                                  .cs$n_levels_b[.cs$b == "MMSE score"], 0))))
    }

    # =============================================================================
    # STEP 2: the models
    # =============================================================================
    cat("\nSTEP 2 -- module against cognition, adjusted for pathology\n")

    # Filter first, standardise second, so the attenuation ratio is not a rescaling.
    fit_block <- function(d, out_col, covs, pred_col = "module") {
        covs <- setdiff(covs, out_col)          # never regress a variable on itself
        D <- data.frame(y = suppressWarnings(as.numeric(d[[out_col]])),
                        x = suppressWarnings(as.numeric(d[[pred_col]])))
        for (c1 in covs) D[[c1]] <- suppressWarnings(as.numeric(d[[c1]]))
        keep <- stats::complete.cases(D)
        D <- D[keep, , drop = FALSE]
        if (nrow(D) < MIN_DONORS) return(NULL)
        D$y <- .z(D$y); D$x <- .z(D$x)
        if (anyNA(D$y) || anyNA(D$x)) return(NULL)
        f  <- stats::as.formula(paste("y ~ x",
            if (length(covs)) paste("+", paste(sprintf("`%s`", covs),
                                               collapse = " + ")) else ""))
        f0 <- stats::as.formula(paste("y ~",
            if (length(covs)) paste(sprintf("`%s`", covs), collapse = " + ")
            else "1"))
        fit  <- try(stats::lm(f,  data = D), silent = TRUE)
        fit0 <- try(stats::lm(f0, data = D), silent = TRUE)
        if (inherits(fit, "try-error") || inherits(fit0, "try-error")) return(NULL)
        s <- summary(fit)
        if (!("x" %in% rownames(s$coefficients))) return(NULL)
        vf <- .safe_vif(fit)
        list(n = nrow(D), donors = d$.donor[keep],
             beta = unname(s$coefficients["x", 1]),
             se = unname(s$coefficients["x", 2]),
             p = unname(s$coefficients["x", 4]),
             delta_r2 = s$r.squared - summary(fit0)$r.squared,
             vif_module = unname(vf["x"]),
             max_vif = .rng(vf, max))
    }

    BLOCKS <- c(list(cov4 = character(0)), PATH_COLS)
    # Keys joined with '|' because outcome labels contain spaces.
    .dkey <- function(...) paste(..., sep = "|")
    rows <- list(); donor_sets <- list()
    for (ds in unique(DF$dataset)) {
        d <- DF[DF$dataset == ds, , drop = FALSE]
        for (pk in PREDS) for (bn in names(BLOCKS)) for (lab in names(COG)) {
            r <- fit_block(d, COG[[lab]], c(BASE_COVS, unname(BLOCKS[[bn]])),
                           pred_col = pk)
            if (is.null(r)) next
            donor_sets[[.dkey(ds, pk, bn, lab)]] <- r$donors
            rows[[length(rows) + 1L]] <- data.frame(
                dataset = ds, predictor = PRED_KEYS[[pk]], key = pk,
                block = bn, outcome = lab, outcome_class = "cognitive",
                n_extra_covariates = length(BLOCKS[[bn]]), n = r$n,
                beta = round(r$beta, 4), se = round(r$se, 4), p_value = r$p,
                delta_r2 = round(r$delta_r2, 4),
                vif_module = round(r$vif_module, 2),
                max_vif = round(r$max_vif, 2), stringsAsFactors = FALSE)
        }
    }
    MOD <- do.call(rbind, rows)
    if (is.null(MOD)) stop("No cognitive model was fitted.", call. = FALSE)
    rownames(MOD) <- NULL
    # BH within dataset x predictor x block, over the cognitive outcomes.
    MOD$p_adj <- NA_real_
    for (k in unique(paste(MOD$dataset, MOD$key, MOD$block))) {
        i <- paste(MOD$dataset, MOD$key, MOD$block) == k
        MOD$p_adj[i] <- stats::p.adjust(MOD$p_value[i], method = "BH")
    }
    MOD$sig <- .stars(MOD$p_adj)

    for (ds in unique(MOD$dataset)) {
        cat(sprintf("\n  %s\n", ds))
        for (lab in names(COG)) {
            s <- MOD[MOD$dataset == ds & MOD$key == "module" &
                     MOD$outcome == lab, , drop = FALSE]
            if (nrow(s) == 0L) next
            cat(sprintf("    %-18s %s\n", lab,
                paste(sprintf("%s beta %+.3f%s (dR2 %.3f, n %d)", s$block, s$beta,
                              s$sig, s$delta_r2, s$n), collapse = "  |  ")))
        }
    }

    # =============================================================================
    # STEP 3: attenuation
    # =============================================================================
    cat("\nSTEP 3 -- how much of the coefficient survives\n")
    att <- list()
    for (ds in unique(MOD$dataset)) for (pk in PREDS) for (lab in names(COG)) {
        b0 <- MOD[MOD$dataset == ds & MOD$key == pk &
                  MOD$block == BASE_BLOCK & MOD$outcome == lab, , drop = FALSE]
        if (nrow(b0) == 0L) next
        for (bn in setdiff(names(BLOCKS), BASE_BLOCK)) {
            b1 <- MOD[MOD$dataset == ds & MOD$key == pk & MOD$block == bn &
                      MOD$outcome == lab, , drop = FALSE]
            if (nrow(b1) == 0L) next
            # Flag an attenuation across two different donor sets.
            same <- setequal(donor_sets[[.dkey(ds, pk, BASE_BLOCK, lab)]],
                             donor_sets[[.dkey(ds, pk, bn, lab)]])
            att[[length(att) + 1L]] <- data.frame(
                dataset = ds, predictor = PRED_KEYS[[pk]], key = pk,
                outcome = lab, block = bn,
                n_base = b0$n, n_adj = b1$n, same_donors = same,
                beta_base = b0$beta, beta_adj = b1$beta,
                # Below the beta floor only the absolute change is reported.
                pct_retained = if (abs(b0$beta) >= MIN_BETA_FOR_PCT)
                    round(100 * b1$beta / b0$beta, 1) else NA_real_,
                abs_change_beta = round(b1$beta - b0$beta, 4),
                below_beta_floor = abs(b0$beta) < MIN_BETA_FOR_PCT,
                padj_base = b0$p_adj, padj_adj = b1$p_adj,
                dR2_base = b0$delta_r2, dR2_adj = b1$delta_r2,
                survives = !is.na(b1$p_adj) & b1$p_adj < ALPHA,
                vif_module = b1$vif_module, max_vif = b1$max_vif,
                stringsAsFactors = FALSE)
        }
    }
    ATT <- do.call(rbind, att)
    if (!is.null(ATT)) {
        rownames(ATT) <- NULL
        h <- ATT[ATT$block == HEADLINE_BLOCK & ATT$key == "module", ,
                 drop = FALSE]
        cat(sprintf("  headline block %s (cov4 + Braak + Thal + CERAD), module:\n",
                    HEADLINE_BLOCK))
        for (i in seq_len(nrow(h)))
            cat(sprintf(paste("    %-8s %-18s beta %+.3f -> %+.3f (%s),",
                              "BH %.4f -> %.4f  %s\n"),
                        h$dataset[i], h$outcome[i], h$beta_base[i],
                        h$beta_adj[i],
                        if (is.finite(h$pct_retained[i]))
                            sprintf("%s%% retained", format(h$pct_retained[i]))
                        else sprintf(paste("change %+.3f; |beta| below the",
                                           "%.2f floor"),
                                     h$abs_change_beta[i], MIN_BETA_FOR_PCT),
                        h$padj_base[i], h$padj_adj[i],
                        if (h$survives[i]) "survives" else "lost"))
    }

    # =============================================================================
    # STEP 4: the mirror
    # =============================================================================
    cat("\nSTEP 4 -- the mirror: pathology adjusted for cognition\n")
    # The mirror adjusts for MMSE and CASI only; adding MOCA would drop donors.
    COG_AS_COV <- unname(COG[intersect(c("MMSE score", "CASI score"),
                                       names(COG))])
    # Cell populations only; the gene-score baselines are excluded.
    .mirror <- intersect(MIRROR_KEYS, PREDS)
    rev_rows <- list()
    for (ds in unique(DF$dataset)) {
        d <- DF[DF$dataset == ds, , drop = FALSE]
        for (pk in .mirror) for (lab in names(PATH))
            for (bn in c("cov4", "+COG")) {
            cv <- if (bn == "cov4") BASE_COVS else c(BASE_COVS, COG_AS_COV)
            r <- fit_block(d, PATH[[lab]], cv, pred_col = pk)
            if (is.null(r)) next
            rev_rows[[length(rev_rows) + 1L]] <- data.frame(
                dataset = ds, predictor = PRED_KEYS[[pk]], key = pk,
                block = bn, outcome = lab, outcome_class = "pathology",
                n = r$n, beta = round(r$beta, 4), se = round(r$se, 4),
                p_value = r$p, delta_r2 = round(r$delta_r2, 4),
                vif_module = round(r$vif_module, 2), max_vif = round(r$max_vif, 2),
                stringsAsFactors = FALSE)
        }
    }
    REV <- do.call(rbind, rev_rows)
    if (!is.null(REV)) {
        rownames(REV) <- NULL
        REV$p_adj <- NA_real_
        # The BH family includes the predictor key.
        for (k in unique(paste(REV$dataset, REV$key, REV$block))) {
            i <- paste(REV$dataset, REV$key, REV$block) == k
            REV$p_adj[i] <- stats::p.adjust(REV$p_value[i], method = "BH")
        }
        REV$sig <- .stars(REV$p_adj)
        for (ds in unique(REV$dataset)) for (pk in unique(REV$key)) {
            s0 <- REV[REV$dataset == ds & REV$key == pk &
                      REV$block == "cov4", , drop = FALSE]
            s1 <- REV[REV$dataset == ds & REV$key == pk &
                      REV$block == "+COG", , drop = FALSE]
            if (nrow(s0) == 0L && nrow(s1) == 0L) next
            cat(sprintf(paste("  %-8s %-20s pathology surviving BH: %d of",
                              "%d at cov4, %d of %d with MMSE + CASI",
                              "added\n"),
                        ds, PRED_KEYS[[pk]],
                        sum(s0$p_adj < ALPHA, na.rm = TRUE), nrow(s0),
                        sum(s1$p_adj < ALPHA, na.rm = TRUE), nrow(s1)))
        }
    }

    # =============================================================================
    # STEP 4b: the joint model
    # =============================================================================
    cat("\nSTEP 4b -- module and non-module in one model\n")

    # Both predictors standardised on the same complete-case set.
    fit_joint <- function(d, out_col, covs, a_col, b_col) {
        covs <- setdiff(covs, out_col)
        D <- data.frame(y = suppressWarnings(as.numeric(d[[out_col]])),
                        a = suppressWarnings(as.numeric(d[[a_col]])),
                        b = suppressWarnings(as.numeric(d[[b_col]])))
        for (c1 in covs) D[[c1]] <- suppressWarnings(as.numeric(d[[c1]]))
        keep <- stats::complete.cases(D)
        D <- D[keep, , drop = FALSE]
        if (nrow(D) < MIN_DONORS) return(NULL)
        D$y <- .z(D$y); D$a <- .z(D$a); D$b <- .z(D$b)
        if (anyNA(D$y) || anyNA(D$a) || anyNA(D$b)) return(NULL)
        f <- stats::as.formula(paste("y ~ a + b",
            if (length(covs)) paste("+", paste(sprintf("`%s`", covs),
                                               collapse = " + ")) else ""))
        f0 <- stats::as.formula(paste("y ~",
            if (length(covs)) paste(sprintf("`%s`", covs), collapse = " + ")
            else "1"))
        fit  <- try(stats::lm(f,  data = D), silent = TRUE)
        fit0 <- try(stats::lm(f0, data = D), silent = TRUE)
        if (inherits(fit, "try-error") || inherits(fit0, "try-error"))
            return(NULL)
        s <- summary(fit)
        if (!all(c("a", "b") %in% rownames(s$coefficients))) return(NULL)
        vf <- .safe_vif(fit)
        list(n = nrow(D),
             beta_a = unname(s$coefficients["a", 1]),
             se_a   = unname(s$coefficients["a", 2]),
             p_a    = unname(s$coefficients["a", 4]),
             beta_b = unname(s$coefficients["b", 1]),
             se_b   = unname(s$coefficients["b", 2]),
             p_b    = unname(s$coefficients["b", 4]),
             vif_a  = unname(vf["a"]), vif_b = unname(vf["b"]),
             delta_r2 = s$r.squared - summary(fit0)$r.squared,
             cor_ab = suppressWarnings(stats::cor(D$a, D$b,
                                                  use = "complete.obs")))
    }

    JNT <- NULL
    if (!(JOINT_B %in% PREDS)) {
        gr10cp("GR-10CP-12", "WARN",
               sprintf("the joint model did not run: `%s` is not among the predictors", JOINT_B))
    } else {
        jr <- list()
        for (ds in unique(DF$dataset)) {
            d <- DF[DF$dataset == ds, , drop = FALSE]
            for (bn in intersect(JOINT_BLOCKS, names(BLOCKS)))
                for (lab in names(COG)) {
                    r <- fit_joint(d, COG[[lab]],
                                   c(BASE_COVS, unname(BLOCKS[[bn]])),
                                   JOINT_A, JOINT_B)
                    if (is.null(r)) next
                    jr[[length(jr) + 1L]] <- data.frame(
                        dataset = ds, block = bn, outcome = lab, n = r$n,
                        term_a = PRED_KEYS[[JOINT_A]],
                        beta_a = round(r$beta_a, 4), se_a = round(r$se_a, 4),
                        p_a = r$p_a, vif_a = round(r$vif_a, 2),
                        term_b = PRED_KEYS[[JOINT_B]],
                        beta_b = round(r$beta_b, 4), se_b = round(r$se_b, 4),
                        p_b = r$p_b, vif_b = round(r$vif_b, 2),
                        cor_ab = round(r$cor_ab, 3),
                        delta_r2_joint = round(r$delta_r2, 4),
                        stringsAsFactors = FALSE)
                }
        }
        JNT <- do.call(rbind, jr)
        if (is.null(JNT))
            gr10cp("GR-10CP-12", "WARN",
                   sprintf(paste("the joint model produced no row: every",
                                 "dataset-block-outcome fell below MIN_DONORS =",
                                 "%d or gave a degenerate fit"), MIN_DONORS))
    }

    if (!is.null(JNT)) {
        rownames(JNT) <- NULL
        # One BH family per term.
        JNT$p_adj_a <- NA_real_; JNT$p_adj_b <- NA_real_
        for (k in unique(paste(JNT$dataset, JNT$block))) {
            i <- paste(JNT$dataset, JNT$block) == k
            JNT$p_adj_a[i] <- stats::p.adjust(JNT$p_a[i], method = "BH")
            JNT$p_adj_b[i] <- stats::p.adjust(JNT$p_b[i], method = "BH")
        }
        .vmax <- pmax(JNT$vif_a, JNT$vif_b)
        JNT$readable <- !is.na(.vmax) & .vmax < JOINT_VIF_MAX
        .sa <- !is.na(JNT$p_adj_a) & JNT$p_adj_a < ALPHA
        .sb <- !is.na(JNT$p_adj_b) & JNT$p_adj_b < ALPHA
        JNT$verdict <- ifelse(!JNT$readable, "not readable",
                       ifelse(.sa & .sb, "both",
                       ifelse(.sa, "module only",
                       ifelse(.sb, "non-module only", "neither"))))
        for (i in seq_len(nrow(JNT)))
            cat(sprintf(paste("  %-8s %-6s %-18s n %3d | module %+.3f (BH",
                              "%.4f) | non-module %+.3f (BH %.4f) | VIF",
                              "%.1f/%.1f | %s\n"),
                        JNT$dataset[i], JNT$block[i], JNT$outcome[i], JNT$n[i],
                        JNT$beta_a[i], JNT$p_adj_a[i], JNT$beta_b[i],
                        JNT$p_adj_b[i], JNT$vif_a[i], JNT$vif_b[i],
                        JNT$verdict[i]))
    }

    # =============================================================================
    # STEP 5: figures
    # =============================================================================
    cat("\nSTEP 5 -- figures\n")
    .n_fig <- 0L
    save_fig <- function(p, base, w = 9, h = 5.5, dpi = 300) {
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
    theme_cp <- function() theme_bw(base_size = 11) +
        theme(plot.title = element_text(face = "bold", size = 12),
              plot.subtitle = element_text(size = 8, colour = "grey30"),
              strip.background = element_rect(fill = "#1F3864"),
              strip.text = element_text(colour = "white", face = "bold", size = 9),
              panel.grid.minor = element_blank(), legend.position = "bottom",
              plot.background = element_rect(fill = "white", colour = NA))

    # Reversed factor so blocks read downwards in legend order; breaks fix the legend.
    MOD$block_f <- factor(MOD$block, levels = rev(names(BLOCKS)))
    MOD$out_f   <- factor(MOD$outcome, levels = rev(names(COG)))
    # Short facet labels; an unmapped label would drop a facet silently.
    .ps <- unname(PRED_SHORT[MOD$predictor])
    if (anyNA(.ps))
        stop("No short facet label for: ",
             paste(unique(MOD$predictor[is.na(.ps)]), collapse = ", "),
             call. = FALSE)
    MOD$pred_f <- factor(.ps, levels = unname(PRED_SHORT[
        vapply(PREDS, function(k) PRED_KEYS[[k]], character(1))]))
    MOD$lo <- MOD$beta - 1.96 * MOD$se
    MOD$hi <- MOD$beta + 1.96 * MOD$se
        # group is declared so the dodge order depends on the block alone.
    p1 <- ggplot(MOD, aes(x = beta, y = out_f, colour = block_f,
                          group = block_f)) +
        geom_vline(xintercept = 0, colour = "grey50", linetype = "dashed") +
        # One layer draws the point and its interval, so they stay paired.
        geom_pointrange(aes(xmin = lo, xmax = hi,
                            shape = factor(p_adj < ALPHA,
                                           levels = c(TRUE, FALSE))),
                        position = position_dodge(width = 0.65),
                        linewidth = 0.4, size = 0.6, fatten = 4,
                        orientation = "y") +
        scale_shape_manual(values = c(`TRUE` = 19, `FALSE` = 1), name = "BH < 0.05",
                           drop = FALSE, na.translate = FALSE) +
        scale_colour_manual(values = BLOCK_FILL, name = NULL,
                            breaks = names(BLOCKS)) +
        facet_grid(pred_f ~ dataset) +
        labs(title = sprintf("Coefficient on cognition, before and after adjusting for pathology [%s]",
                             ARM_TAG),
             subtitle = paste("Filled = clears BH within block; blocks read",
                              "downwards in legend order. +BTC is the headline",
                              "(cov4 + Braak + Thal + CERAD); +ALL adds ADNC on",
                              "top and is collinear by construction.",
                              "\nDIRECTION: MMSE and CASI score higher when",
                              "cognition is better, cognitive status higher",
                              "when it is worse, so a NEGATIVE beta on the",
                              "first two and a POSITIVE beta on the third both",
                              "mean more module = worse cognition."),
             x = "standardised beta (95% CI)", y = NULL) +
        theme_cp()
    print(p1)
    save_fig(p1, file.path(OUT_DIR, paste0(FILE_PREFIX, "_CogVsPath_Attenuation")),
             w = 9, h = min(MAX_FIG_H, max(5.5, 2.2 * length(PREDS))))

    p2 <- ggplot(MOD, aes(x = out_f, y = delta_r2, fill = block_f)) +
        geom_col(position = position_dodge(width = 0.75), width = 0.7,
                 colour = "black", linewidth = 0.25) +
        scale_fill_manual(values = BLOCK_FILL, name = NULL,
                          breaks = names(BLOCKS)) +
        facet_grid(pred_f ~ dataset) + coord_flip() +
        labs(title = sprintf("Variance each predictor explains beyond its covariates [%s]",
                             ARM_TAG),
             subtitle = "delta R-squared over the covariates-only model, same donors within a block",
             x = NULL, y = expression(paste(Delta, R^2))) +
        theme_cp()
    print(p2)
    save_fig(p2, file.path(OUT_DIR, paste0(FILE_PREFIX, "_CogVsPath_DeltaR2")),
             w = 9, h = min(MAX_FIG_H, max(5, 2 * length(PREDS))))

    # == Publication panels ==================================================
    # Same styling as p1; only the subset, axis limits and titles differ.
    .clip <- character(0)
    PUB <- MOD[MOD$key == PUB_KEY & MOD$outcome %in% PUB_OUTCOMES, , drop = FALSE]
    if (nrow(PUB) > 0L) {
        PUB$out_f <- factor(PUB$outcome, levels = rev(PUB_OUTCOMES))
        .o <- PUB$lo < PUB_ATT_XLIM[1] | PUB$hi > PUB_ATT_XLIM[2]
        if (any(.o)) .clip <- c(.clip, sprintf(
            "attenuation %s/%s/%s [%.3f, %.3f]", PUB$dataset[.o], PUB$outcome[.o],
            PUB$block[.o], PUB$lo[.o], PUB$hi[.o]))
            # group is declared so the dodge order depends on the block alone.
        p3 <- ggplot(PUB, aes(x = beta, y = out_f, colour = block_f,
                              group = block_f)) +
            geom_vline(xintercept = 0, colour = "grey50", linetype = "dashed") +
            # Both shape levels are declared so the legend is the same on every arm.
            # One layer draws the point and its interval, so they stay paired.
            geom_pointrange(aes(xmin = lo, xmax = hi,
                                shape = factor(p_adj < ALPHA,
                                               levels = c(TRUE, FALSE))),
                            position = position_dodge(width = PUB_DODGE),
                            linewidth = 0.4, size = 0.6, fatten = 4,
                            orientation = "y") +
            scale_shape_manual(values = c(`TRUE` = 19, `FALSE` = 1),
                               name = "BH < 0.05", drop = FALSE,
                               na.translate = FALSE) +
            scale_colour_manual(values = BLOCK_FILL, name = NULL,
                                breaks = names(BLOCKS)) +
            scale_x_continuous(limits = PUB_ATT_XLIM) +
            facet_grid(. ~ dataset) +
            labs(title = sprintf("Cognition before and after adjusting for pathology [%s]",
                                 ARM_TAG),
                 subtitle = paste("Module Post-IP only. Filled = clears BH within",
                                  "block. +BTC is the headline (cov4 + Braak +",
                                  "Thal + CERAD); +ALL adds ADNC and is",
                                  "collinear, shown as sensitivity.",
                                  "\nMMSE and CASI are higher when cognition is",
                                  "better, so a NEGATIVE beta means more module",
                                  "cells with worse cognition.",
                                  sprintf("x axis FIXED at [%.2f, %.2f] across arms.",
                                          PUB_ATT_XLIM[1], PUB_ATT_XLIM[2])),
                 x = "standardised beta (95% CI)", y = NULL) +
            theme_cp()
        print(p3)
        save_fig(p3, file.path(OUT_DIR,
                 paste0(FILE_PREFIX, "_CogVsPath_Pub_Attenuation")),
                 w = 9, h = PUB_ATT_H)
    }

    if (!is.null(REV)) {
        MIR <- REV[REV$key == PUB_KEY, , drop = FALSE]
    } else MIR <- NULL
    if (!is.null(MIR) && nrow(MIR) > 0L) {
        MIR$lo <- MIR$beta - 1.96 * MIR$se
        MIR$hi <- MIR$beta + 1.96 * MIR$se
        MIR$block_f <- factor(MIR$block, levels = c("+COG", "cov4"))
        MIR_FILL <- c(cov4 = BLOCK_FILL[["cov4"]],
                      `+COG` = BLOCK_FILL[["+BTC"]])
        MIR$out_f   <- factor(MIR$outcome, levels = rev(names(PATH)))
        .o <- MIR$lo < PUB_MIR_XLIM[1] | MIR$hi > PUB_MIR_XLIM[2]
        if (any(.o)) .clip <- c(.clip, sprintf(
            "mirror %s/%s/%s [%.3f, %.3f]", MIR$dataset[.o], MIR$outcome[.o],
            MIR$block[.o], MIR$lo[.o], MIR$hi[.o]))
            # group is declared so the dodge order depends on the block alone.
        p4 <- ggplot(MIR, aes(x = beta, y = out_f, colour = block_f,
                              group = block_f)) +
            geom_vline(xintercept = 0, colour = "grey50", linetype = "dashed") +
            # Both shape levels are declared so the legend is the same on every arm.
            # One layer draws the point and its interval, so they stay paired.
            geom_pointrange(aes(xmin = lo, xmax = hi,
                                shape = factor(p_adj < ALPHA,
                                               levels = c(TRUE, FALSE))),
                            position = position_dodge(width = PUB_DODGE),
                            linewidth = 0.4, size = 0.6, fatten = 4,
                            orientation = "y") +
            scale_shape_manual(values = c(`TRUE` = 19, `FALSE` = 1),
                               name = "BH < 0.05", drop = FALSE,
                               na.translate = FALSE) +
            # The mirror blocks are cov4 and +COG; +COG takes +BTC's colour.
            scale_colour_manual(values = MIR_FILL, name = NULL,
                                breaks = c("cov4", "+COG"),
                                labels = c("cov4", "+MMSE+CASI")) +
            scale_x_continuous(limits = PUB_MIR_XLIM) +
            facet_grid(. ~ dataset) +
            labs(title = sprintf("The mirror: pathology before and after adjusting for cognition [%s]",
                                 ARM_TAG),
                 subtitle = paste("Module Post-IP only. Filled = clears BH within block.",
                                  "\nA population that predicts pathology but stops",
                                  "doing so once MMSE and CASI are in the model",
                                  "carries a one-directional association.",
                                  sprintf("x axis FIXED at [%.2f, %.2f] across arms.",
                                          PUB_MIR_XLIM[1], PUB_MIR_XLIM[2])),
                 x = "standardised beta (95% CI)", y = NULL) +
            theme_cp()
        print(p4)
        save_fig(p4, file.path(OUT_DIR,
                 paste0(FILE_PREFIX, "_CogVsPath_Pub_Mirror")),
                 w = 9, h = PUB_MIR_H)
    }

    # =============================================================================
    # STEP 6: output
    # =============================================================================
    cat("\nSTEP 6 -- output\n")
    readme <- data.frame(
        item = c("question", "predictor", "base block", "headline block",
                 "cognitive outcomes", "pathology outcomes", "mirror", "unit"),
        value = c(
            "does the module still predict cognition once plaque and tangle burden are in the model",
            "CMRF module Post-IP %, the rest of the 2x2 and the five comparator baselines",
            paste(BASE_COVS, collapse = " + "),
            paste0(HEADLINE_BLOCK, " = cov4 + Braak + Thal + CERAD"),
            paste(names(COG), collapse = ", "),
            paste(names(PATH), collapse = ", "),
            "pathology outcomes adjusted for cov4 + MMSE + CASI (Reverse sheet)",
            "donor"), stringsAsFactors = FALSE)
    # Plotting columns are dropped from the sheet and the RDS.
    .plot_cols <- c("block_f", "out_f", "pred_f", "lo", "hi")
    MOD_OUT <- MOD[, setdiff(names(MOD), .plot_cols), drop = FALSE]
    wb <- list(Read_Me = readme, Models = MOD_OUT)
    if (!is.null(ATT)) wb[["Attenuation"]] <- ATT
    if (!is.null(COR)) wb[["Cognitive_Agreement"]] <- COR
    if (!is.null(REV)) wb[["Reverse"]] <- REV
    if (!is.null(JNT)) wb[["Joint_Mod_vs_NonMod"]] <- JNT
    xlsx_path <- file.path(OUT_DIR, paste0(FILE_PREFIX, "_CogVsPath.xlsx"))
    tryCatch({
        write_xlsx(wb, xlsx_path)
        cat(sprintf("  %d sheet(s) written to %s\n", length(wb),
                    basename(xlsx_path)))
    }, error = function(e)
        cat(sprintf("  workbook not written: %s\n", conditionMessage(e))))

    rds_out <- list(models = MOD_OUT, attenuation = ATT, reverse = REV,
                    joint = JNT,
                    correlations = COR, predictors = PRED_KEYS[PREDS],
                    meta = list(arm_tag = ARM_TAG, config_file = config_file,
                                norm_tag = NORM_TAG, GOI = GOI,
                                supertype = PRIMARY_STYPE,
                                base_covariates = BASE_COVS,
                                cognitive = COG, pathology = PATH,
                                blocks = BLOCKS, headline_block = HEADLINE_BLOCK,
                                part10h_rds = P10H_RDS,
                                written = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                                R_version = R.version.string,
                                joint_terms = c(a = JOINT_A, b = JOINT_B),
                                schema = "M10CP v2"))
    rds_path <- file.path(OUT_DIR, paste0(FILE_PREFIX, "_CogVsPath_results.rds"))
    tryCatch({
        saveRDS(rds_out, rds_path)
        cat(sprintf("  %s written (%.2f MB)\n", basename(rds_path),
                    file.size(rds_path) / 1e6))
    }, error = function(e)
        cat(sprintf("  results RDS not written: %s\n", conditionMessage(e))))

    # =============================================================================
    # STEP 7: guardrails
    # =============================================================================
    cat("\nSTEP 7 -- guardrails\n")

    # The cov4 baseline must reproduce Part 10h.
    MVSH <- SHEETS[["Multivariate"]]
    .xc <- NULL; .xc_note <- ""
    .pl <- DECL[["PRED_LABEL"]]
    .mod_names <- unique(c("Module Post-IP", "pct_post",
                           if (!is.null(.pl)) names(.pl)[.pl == "Module Post-IP"],
                           if (!is.null(.pl)) unname(.pl["pct_post"])))
    .mod_names <- .mod_names[!is.na(.mod_names) & nzchar(.mod_names)]
    if (is.data.frame(MVSH) && all(c("outcome", "covariate_set", "dataset") %in%
                                   names(MVSH))) {
        .pcol <- NA_character_
        for (cn in intersect(c("predictor", "arm"), names(MVSH)))
            if (any(as.character(MVSH[[cn]]) %in% .mod_names)) { .pcol <- cn; break }
        if (!is.na(.pcol)) {
            ref <- MVSH[as.character(MVSH[[.pcol]]) %in% .mod_names &
                        as.character(MVSH$covariate_set) == BASE_BLOCK, ,
                        drop = FALSE]
            if (!is.null(ref[["tier"]]) && nrow(ref) > 0L) {
                tt <- stats::aggregate(list(n = suppressWarnings(
                    as.numeric(ref$n_donors))), by = list(tier = ref$tier),
                    FUN = function(v) max(v, na.rm = TRUE))
                ref <- ref[ref$tier == tt$tier[which.max(tt$n)], , drop = FALSE]
            }
            # Module rows only.
            mine <- MOD[MOD$block == BASE_BLOCK & MOD$key == "module", ,
                        drop = FALSE]
            cmp <- merge(mine[, c("dataset", "outcome", "n", "beta", "delta_r2")],
                         ref[, intersect(c("dataset", "outcome", "n_donors",
                                           "beta_std", "delta_R2"), names(ref))],
                         by = c("dataset", "outcome"))
            if (nrow(cmp) > 0L) {
                cmp$d_beta <- abs(suppressWarnings(as.numeric(cmp$beta_std)) - cmp$beta)
                cmp$d_n <- abs(suppressWarnings(as.numeric(cmp$n_donors)) - cmp$n)
                cmp$d_r2 <- abs(suppressWarnings(as.numeric(cmp$delta_R2)) - cmp$delta_r2)
                .xc <- cmp
            } else .xc_note <- "no row merged on dataset + outcome"
        } else .xc_note <- "no module name matched the predictor or arm column"
    } else .xc_note <- "no usable Multivariate sheet"
    if (is.null(.xc))
        gr10cp("GR-10CP-04", "WARN",
               sprintf("cov4 baseline not checked against Part 10h (%s)", .xc_note))
    else {
        .bad <- sum(.xc$d_beta > 0.01 | .xc$d_n > 0 | .xc$d_r2 > 0.01, na.rm = TRUE)
        gr10cp("GR-10CP-04", if (.bad > 0L) "FAIL" else "PASS",
               sprintf(paste("cov4 vs Part 10h over %d row(s): max |d beta|",
                             "%.4f, |d n| %d, |d delta-R2| %.4f.%s"),
                       nrow(.xc), max(.xc$d_beta, na.rm = TRUE),
                       as.integer(max(.xc$d_n, na.rm = TRUE)),
                       max(.xc$d_r2, na.rm = TRUE),
                       if (.bad > 0L)
                           sprintf(" %d row(s) disagree", .bad)
                       else sprintf(" Unverified here (no Part 10h row): %s",
                                    paste(setdiff(names(COG),
                                                  unique(.xc$outcome)),
                                          collapse = ", "))))
    }

    # Donor sets must match across blocks; the ratio is read only above the beta floor.
    if (is.null(ATT)) gr10cp("GR-10CP-02", "WARN", "no attenuation row was produced")
    else {
        .diff <- sum(!ATT$same_donors)
        .flr <- sum(ATT$below_beta_floor)
        gr10cp("GR-10CP-02", if (.diff > 0L) "WARN" else "PASS",
               sprintf(paste("realised n per block: %s | %d of %d attenuation",
                             "row(s) compare different donor sets%s | %d",
                             "row(s) below the |beta| %.2f floor"),
                       paste(vapply(names(BLOCKS), function(b) {
                           v <- MOD$n[MOD$block == b]
                           sprintf("%s %s-%s", b,
                                   format(.rng(v, min)), format(.rng(v, max)))
                       }, character(1)), collapse = ", "),
                       .diff, nrow(ATT),
                       if (.diff > 0L) " (sample differs)"
                       else "",
                       .flr, MIN_BETA_FOR_PCT))
    }

    # Collinearity of the pathology blocks.
    .vhi <- MOD[is.finite(MOD$max_vif) & MOD$max_vif >= VIF_READABLE, , drop = FALSE]
    gr10cp("GR-10CP-03", if (nrow(.vhi) > 0L) "WARN" else "PASS",
           sprintf("worst VIF per block: %s | %d model(s) at or above %g%s",
                   paste(vapply(names(BLOCKS), function(b)
                       sprintf("%s %.1f", b,
                               .rng(MOD$max_vif[MOD$block == b], max)),
                       character(1)), collapse = ", "),
                   nrow(.vhi), VIF_READABLE,
                   if (nrow(.vhi) > 0L)
                       paste0(". Expected in +ALL; the headline is ", HEADLINE_BLOCK)
                   else ""))

    # Headline result.
    if (is.null(ATT)) gr10cp("GR-10CP-05", "WARN", "no headline row to report")
    else {
        h <- ATT[ATT$block == HEADLINE_BLOCK & ATT$key == "module", ,
                 drop = FALSE]
        surv <- h[h$survives, , drop = FALSE]
        gr10cp("GR-10CP-05", "PASS",
               sprintf(paste("headline at %s: %d of %d dataset-outcome",
                             "row(s) survive. %s | %s"),
                       HEADLINE_BLOCK, nrow(surv), nrow(h),
                       if (nrow(surv) == 0L)
                           "nothing survives"
                       else paste("surviving:",
                                  paste(sprintf("%s %s (%s, BH %.4f)",
                                                surv$dataset, surv$outcome,
                                                ifelse(is.finite(
                                                    surv$pct_retained),
                                                    sprintf("%.0f%% of beta",
                                                            surv$pct_retained),
                                                    sprintf(paste(
                                                        "beta %+.3f, below",
                                                        "the %.2f floor"),
                                                        surv$abs_change_beta,
                                                        MIN_BETA_FOR_PCT)),
                                                surv$padj_adj),
                                        collapse = "; ")),
                       if (is.null(REV)) "the mirror did not run"
                       else sprintf("mirror: %d of %d pathology row(s) survive adding MMSE + CASI",
                                    sum(REV$p_adj[REV$block == "+COG"] < ALPHA,
                                        na.rm = TRUE),
                                    sum(REV$block == "+COG"))))
    }

    # Do the external baselines also survive the headline block?
    .bh <- if (is.null(ATT)) NULL else
        ATT[ATT$block == HEADLINE_BLOCK & ATT$key %in% BASE_KEYS &
            ATT$outcome %in% c("MMSE score", "CASI score"), , drop = FALSE]
    if (is.null(.bh) || nrow(.bh) == 0L)
        gr10cp("GR-10CP-10", "WARN",
               "no baseline row at the headline block")
    else {
        .mh <- ATT[ATT$block == HEADLINE_BLOCK & ATT$key == "module" &
                   ATT$outcome %in% c("MMSE score", "CASI score"), ,
                   drop = FALSE]
        gr10cp("GR-10CP-10", "PASS",
               sprintf(paste("baselines at %s, MMSE/CASI rows: module %d of",
                             "%d; %s"),
                       HEADLINE_BLOCK, sum(.mh$survives), nrow(.mh),
                       paste(vapply(unique(.bh$predictor), function(p) {
                           z <- .bh[.bh$predictor == p, , drop = FALSE]
                           sprintf("%s %d of %d", p, sum(z$survives), nrow(z))
                       }, character(1)), collapse = "; ")))
    }

    # Joint model and retention, module against non-module.
    if (!is.null(JNT)) {
        .jh <- JNT[JNT$block == HEADLINE_BLOCK &
                   JNT$outcome %in% c("MMSE score", "CASI score"), ,
                   drop = FALSE]
        .jr <- .jh[.jh$readable, , drop = FALSE]
        .att_pair <- function(pk) {
            if (is.null(ATT)) return("no attenuation frame")
            z <- ATT[ATT$block == HEADLINE_BLOCK & ATT$key == pk &
                     ATT$outcome %in% c("MMSE score", "CASI score"), ,
                     drop = FALSE]
            if (nrow(z) == 0L) return("no row")
            paste(vapply(seq_len(nrow(z)), function(i)
                sprintf("%s %s %s", z$dataset[i], z$outcome[i],
                        if (isTRUE(is.finite(z$pct_retained[i])))
                            sprintf("%.0f%%", z$pct_retained[i])
                        else "below the beta floor"),
                character(1)), collapse = "; ")
        }
        gr10cp("GR-10CP-12",
               if (nrow(.jh) == 0L) "WARN"
               else if (nrow(.jr) == 0L) "WARN" else "PASS",
               sprintf(paste("joint model at %s over MMSE/CASI: %d of %d",
                             "row(s) readable at VIF < %.0f (module-non-module",
                             "r %s). Verdicts: %s.%s | retention at %s:",
                             "module: %s || non-module: %s"),
                       HEADLINE_BLOCK, nrow(.jr), nrow(.jh), JOINT_VIF_MAX,
                       if (nrow(.jh) && any(is.finite(.jh$cor_ab)))
                           sprintf("%.2f to %.2f",
                                   min(.jh$cor_ab, na.rm = TRUE),
                                   max(.jh$cor_ab, na.rm = TRUE))
                       else "NA",
                       if (nrow(.jh))
                           paste(vapply(unique(.jh$verdict), function(v)
                               sprintf("%s %d", v, sum(.jh$verdict == v)),
                               character(1)), collapse = ", ")
                       else "none",
                       if (nrow(.jr) == 0L && nrow(.jh) > 0L)
                           " No row is readable."
                       else "",
                       HEADLINE_BLOCK, .att_pair(JOINT_A), .att_pair(JOINT_B)))
    }

    # Fixed axes are safe only if no interval falls outside them.
    gr10cp("GR-10CP-13", if (length(.clip) > 0L) "FAIL" else "PASS",
           sprintf(paste("publication panels on fixed axes: attenuation [%.2f,",
                         "%.2f], mirror [%.2f, %.2f]. %s"),
                   PUB_ATT_XLIM[1], PUB_ATT_XLIM[2],
                   PUB_MIR_XLIM[1], PUB_MIR_XLIM[2],
                   if (length(.clip) > 0L)
                       paste0(length(.clip), " interval(s) clipped: ",
                              paste(.clip, collapse = "; "))
                   else "0 interval(s) clipped."))

    gr10cp("GR-10CP-06", if (.n_fig == 0L) "WARN" else "PASS",
           sprintf("%d figure(s) written as PDF and PNG",
                   .n_fig))
    .files <- c(rds_path, xlsx_path)
    .have <- file.exists(.files)
    gr10cp("GR-10CP-07", if (all(.have)) "PASS" else "FAIL",
           sprintf("outputs: %s", paste(sprintf("%s %s", basename(.files),
                                                ifelse(.have, "written", "MISSING")),
                                        collapse = " | ")))

    # =============================================================================
    # LEDGER
    # =============================================================================
    .gr_lines <- c(
        strrep("=", 72),
        "Part 10h Cognition vs Pathology - Guardrail Log [CP10K]",
        paste("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
        paste("Arm :", ARM_TAG),
        "Unit: donor.",
        strrep("=", 72), "",
        vapply(.gr_log, function(g)
            sprintf("  %s  %-18s  %s",
                    switch(g$status, PASS = "[PASS]", FAIL = "[FAIL]",
                           WARN = "[WARN]"), g$id, g$detail),
            character(1)),
        "", strrep("=", 72),
        sprintf("  Total : %d", length(.gr_log)),
        sprintf("  PASS  : %d", .gr_pass),
        sprintf("  WARN  : %d", .gr_warn),
        sprintf("  FAIL  : %d", .gr_fail), "",
        if (.gr_fail > 0) "  Result: failures present"
        else "  Result: no failures",
        strrep("=", 72))
    writeLines(.gr_lines, GR_LOG_10CP)
    cat(sprintf("\n  Final guardrail tally: %d PASS / %d WARN / %d FAIL\n",
                .gr_pass, .gr_warn, .gr_fail))

    if (file.exists(INCOMPLETE_MARKER_10CP)) unlink(INCOMPLETE_MARKER_10CP)
    .p10cp_completed <- TRUE

    cat("\n", strrep("=", 72), "\n")
    cat("  Part10h stage 2 Cog-vs-Path complete [", ARM_TAG, "]\n", sep = "")
    cat("  Output directory:", OUT_DIR, "\n")
    cat(strrep("=", 72), "\n")

    invisible(NULL)
}

cep_run_arms(PART_DIR_TAG, run_one_arm)
