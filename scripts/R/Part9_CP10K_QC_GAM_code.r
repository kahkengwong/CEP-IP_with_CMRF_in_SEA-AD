# =============================================================================
# Part9_CP10K_QC_GAM_code.r
# SEA-AD CEP-IP/CMRF, Part 9, export QC, heatmaps and GAM fits
#   (Figure 2, Supplementary Figures S2, S3 and S5, Table S7)
#   For each arm (Mic13_PVM3, PVM3_all, PVM1_all, HKG8_matched): checks the
#   Part 8 export (cell counts, CP10K normalisation), adds missing gene
#   columns to the per-region GAM CSVs, draws the GOI and partner-gene
#   heatmaps, and fits GAMs of each partner composite on the GOI per region
#   and pooled, with a donor random intercept and PRSS selection of k.
#
# Usage (from the SEA-AD folder, or set SEAAD_ROOT), after Part 8
#   source("Part9_CP10K_QC_GAM_code.r")                         # all arms
#   CEP_ARM <- "PVM1_all"; source("Part9_CP10K_QC_GAM_code.r")  # one arm
#   Needs CEP_arms.R, R_config.R and R_config_HKG8.R in the same folder.
#
# Outputs (Part9/Part9_<arm>/): check ledger (text, workbook), heatmaps (PDF,
#   PNG), GAM results workbook, GAM models RDS, GAM curve figures and log.
# =============================================================================


# =============================================================================
# SECTION 0: Libraries
# =============================================================================

library(Matrix)
library(writexl)
library(mgcv)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(scales)
library(viridis)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(parallel)

set.seed(123)

# =============================================================================
# SECTION 0a: Shared arm definitions
# =============================================================================
# CEP_arms.R holds the arm table, the folder layout and the driver.
if (!exists("CEP_HOME")) CEP_HOME <- Sys.getenv("SEAAD_ROOT", getwd())
.p9_harness <- file.path(CEP_HOME, "CEP_arms.R")
if (!file.exists(.p9_harness))
    stop("CEP_arms.R not found at: ", .p9_harness,
         "\n  It holds the arm table shared by Parts 9 onwards.",
         call. = FALSE)
source(.p9_harness)

PART_TAG <- "Part9"


# =============================================================================
# SECTION 0d: Pure model helpers and PSOCK workers
# =============================================================================
# Defined at the top level so that PSOCK workers receive them without the arm's
# data; the configuration constants are passed in as arguments.

DONOR_RE_LABEL <- "s(donor)"

p9_gam_formula <- function(k, use_re, bs) {
    rhs <- sprintf("GOI + s(GOI, bs = \"%s\", k = %d)", bs, as.integer(k))
    if (use_re) rhs <- paste0(rhs, " + s(donor, bs = \"re\")")
    stats::as.formula(paste("Expression ~", rhs))
}

p9_goi_smooth_idx <- function(model) {
    if (is.null(model$smooth) || length(model$smooth) == 0) return(NA_integer_)
    labs <- vapply(model$smooth, function(s) s$label, character(1))
    i    <- which(labs == "s(GOI)")
    if (length(i) == 1L) return(i)
    i <- which(labs != DONOR_RE_LABEL)
    if (length(i) >= 1L) return(i[1])
    NA_integer_
}

p9_smooth_lambda <- function(model, i) {
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

p9_smooth_stat <- function(sm_summary, label, column) {
    st <- sm_summary$s.table
    if (is.null(st) || nrow(st) == 0) return(NA_real_)
    i <- match(label, rownames(st))
    if (is.na(i)) i <- 1L
    if (!(column %in% colnames(st))) return(NA_real_)
    unname(st[i, column])
}

p9_model_has_re <- function(model) {
    !is.null(model$smooth) &&
        any(vapply(model$smooth, function(s) s$label, character(1)) ==
                DONOR_RE_LABEL)
}

p9_calculate_prss <- function(model) {
    RSS <- sum(stats::residuals(model)^2)
    i   <- p9_goi_smooth_idx(model)
    if (is.na(i))
        return(list(PRSS = RSS, RSS = RSS,
                    f_double_prime_integral = NA_real_, lambda = NA_real_))
    sm     <- model$smooth[[i]]
    lambda <- p9_smooth_lambda(model, i)
    S      <- sm$S[[1]]
    beta   <- stats::coef(model)[sm$first.para:sm$last.para]
    f_dpi  <- as.numeric(t(beta) %*% S %*% beta)
    list(PRSS = RSS + lambda * f_dpi, RSS = RSS,
         f_double_prime_integral = f_dpi, lambda = lambda)
}

p9_extract_reml_iters <- function(model) {
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

# == One phase-2 PRSS iteration, run on a worker ==============================
# Returns summaries only; the master refits the winning model.
p9_prss_worker <- function(job) {
    ctx <- p9_fit_ctx
    res <- tryCatch({
        set.seed(123L + job$i)
        m  <- mgcv::gam(p9_gam_formula(job$k, ctx$use_re, ctx$bs),
                        data = ctx$data, weights = w,
                        method = ctx$method,
                        select = ctx$select, gamma = ctx$gamma)
        pc <- p9_calculate_prss(m)
        ri <- p9_extract_reml_iters(m)
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
                 smooth_edf    = p9_smooth_stat(sm, "s(GOI)", "edf"),
                 p_value       = p9_smooth_stat(sm, "s(GOI)", "p-value")))
    }, error = function(e) list(ok = FALSE, message = conditionMessage(e)))
    res$i <- job$i
    res$k <- job$k
    res
}

# == One leave-one-gene-out refit, run on a worker ============================
p9_loo_worker <- function(job) {
    ctx <- p9_loo_ctx
    de  <- tryCatch({
        loo_avg <- if (length(job$genes_in) == 1L)
            ctx$mat[, job$genes_in]
        else
            rowMeans(ctx$mat[, job$genes_in, drop = FALSE])
        d <- data.frame(GOI = ctx$goi_vec, Expression = loo_avg,
                        donor = ctx$donor_vec,
                        w = if (is.null(ctx$w)) rep(1, length(ctx$goi_vec))
                            else ctx$w)
        set.seed(123L)
        m  <- mgcv::gam(p9_gam_formula(ctx$k_use, ctx$use_re, ctx$bs),
                        data = d, weights = w,
                        method = ctx$method,
                        select = ctx$select, gamma = ctx$gamma)
        set.seed(123L)
        nd <- if (ctx$use_re)
            stats::deviance(mgcv::gam(Expression ~ 1 + s(donor, bs = "re"),
                                      data = d, weights = w,
                                      method = ctx$method))
        else
            stats::deviance(mgcv::gam(Expression ~ 1, data = d,
                                      weights = w))
        if (is.finite(nd) && nd > 0)
            as.numeric((nd - stats::deviance(m)) / nd)
        else NA_real_
    }, error = function(e) NA_real_)
    list(gene = job$gene_out, de = de)
}


# =============================================================================
# SECTION 0b: The pipeline for one arm
# =============================================================================
# Everything below runs once per arm; the arm's configuration is loaded first.

run_one_arm <- function(arm_tag, config_file, out_dir, file_tag) {

    config_path <- file.path(CEP_HOME, config_file)
    if (!file.exists(config_path))
        stop("Configuration not found: ", config_path, call. = FALSE)

    # Output printed before the log exists is buffered and replayed into it.
    .p9_preamble <- character(0)

    assign("CEP_ARM", arm_tag, envir = globalenv())
    on.exit(suppressWarnings(rm("CEP_ARM", envir = globalenv())), add = TRUE)

    .p9_preamble <- tryCatch(
        utils::capture.output({
            cat("\n", strrep("=", 76), "\n", sep = "")
            cat("  ARM: ", arm_tag, "\n", sep = "")
            cat("  Config: ", config_path, "\n", sep = "")
            cat(strrep("=", 76), "\n", sep = "")
            # The configuration reads CEP_ARM; it is removed when this function exits.
            source(config_path)
        }),
        error = function(e) {
            cat("\n", strrep("=", 76), "\n", sep = "")
            cat("  ARM: ", arm_tag, "\n", sep = "")
            cat("  Config: ", config_path, "\n", sep = "")
            cat(strrep("=", 76), "\n", sep = "")
            stop(e)
        })

    # == Confirm the configuration ============================================
    .p9_needed <- c("ARM_TAG", "ARM_DIR", "R_EXPORT_ROOT", "SHARED_DIR",
                    "GOI", "DFG", "SELECT_GENES", "MOD1_GENES", "HKG_GENES",
                    "DFG_COMPOSITE_LABEL", "RUN_SUFFIX", "PRIMARY_STYPE",
                    "NORM_TAG", "NORM_TARGET_SUM", "NORM_CEILING")
    .p9_absent <- .p9_needed[!vapply(.p9_needed, exists, logical(1))]
    if (length(.p9_absent) > 0)
        stop("The configuration loaded from\n    ", config_path,
             "\n  does not define: ", paste(.p9_absent, collapse = ", "),
             call. = FALSE)

    if (!identical(ARM_TAG, arm_tag))
        stop("The configuration reports ARM_TAG '", ARM_TAG,
             "' but this arm is '", arm_tag,
             "'.", call. = FALSE)

    if (!dir.exists(R_EXPORT_ROOT))
        stop("R_EXPORT_ROOT does not exist: ", R_EXPORT_ROOT,
             "\n  Config: ", config_path, call. = FALSE)

    cat("  Config OK: ", ARM_TAG, " | ", GOI, " + ", length(DFG), " DFG | ",
        NORM_TAG, "\n\n", sep = "")

    # == Arm mode =============================================================
    # The HKG comparator is HKG_GENES minus the panel under test (empty in HKG8).
    IS_HKG_ARM <- GOI %in% HKG_GENES

    # == Fallback for PRSS_EARLY_STOP =========================================
    if (!exists("PRSS_EARLY_STOP")) {
        PRSS_EARLY_STOP <- 20L
        cat("  NOTE: PRSS_EARLY_STOP not in R_config.R -- using default 20\n")
    }

    # == Normalisation tolerances =============================================
    # NORM_TOL: tolerance for the reconstruction and row-sum checks (float32 CSVs).
    NORM_TOL     <- 1e-4
    NORM_CSV_TOL <- 1e-4

    # == HKG comparator set ===================================================
    HKG_USE <- setdiff(HKG_GENES, c(GOI, DFG))
    if (IS_HKG_ARM && length(HKG_USE) > 0) {
        stop("The GOI is a housekeeping gene but ", length(HKG_USE),
             " housekeeping gene(s) remain outside the panel: ",
             paste(HKG_USE, collapse = ", "),
             call. = FALSE)
    }
    if (length(HKG_USE) == 0) {
        cat("  Control arm: no separate HKG comparator; GR-11 does not apply.\n")
    } else if (length(HKG_USE) < length(HKG_GENES)) {
        cat("  NOTE: removed", length(HKG_GENES) - length(HKG_USE),
            "gene(s) from the HKG comparator because they are in the",
            "GOI/DFG panel\n")
    }

    # == Output paths =========================================================
    # out_dir is <RESULTS_ROOT>/Part9/Part9_<subdir>/, created by the driver.
    PART9_DIR     <- out_dir
    SUBREGION_DIR <- file.path(PART9_DIR, "subregions")

    if (!dir.exists(PART9_DIR))     dir.create(PART9_DIR,     recursive = TRUE)
    if (!dir.exists(SUBREGION_DIR)) dir.create(SUBREGION_DIR, recursive = TRUE)

    # Every output name is built from this prefix.
    OUT_PREFIX <- cep_out_prefix(PART_TAG, NORM_TAG, file_tag)
    LOG_PATH   <- file.path(PART9_DIR, paste0(OUT_PREFIX, "_QC_GAM_Log.txt"))
    GR_TXT_PATH <- file.path(PART9_DIR, paste0(OUT_PREFIX, "_Guardrail_Log.txt"))
    GR_XLSX_PATH <- file.path(PART9_DIR,
                              paste0(OUT_PREFIX, "_Guardrail_Summary.xlsx"))

    # == Log file, one per arm ================================================
    # Printed output and the message stream both go to the log; LOG_MESSAGES = FALSE
    # keeps warnings on the console instead.
    LOG_MESSAGES <- TRUE
    MSG_PATH <- file.path(PART9_DIR,
                          paste0(OUT_PREFIX, "_QC_GAM_Log_messages.tmp"))

    while (sink.number() > 0) sink()
    # Only one message diversion can be active, so one conditional unwind is enough.
    if (sink.number(type = "message") != 2L) sink(type = "message")
    log_con <- file(LOG_PATH, open = "wt")
    msg_con <- if (isTRUE(LOG_MESSAGES)) file(MSG_PATH, open = "wt") else NULL
    sink(log_con, append = FALSE, split = TRUE)
    if (!is.null(msg_con)) sink(msg_con, type = "message")
    on.exit({
        tryCatch({
            if (sink.number(type = "message") != 2L) sink(type = "message")
            while (sink.number() > 0) sink()
            if (!is.null(msg_con)) close(msg_con)
            close(log_con)
            # Append the message stream to the log and remove the temporary file.
            if (!is.null(msg_con) && file.exists(MSG_PATH)) {
                .msg_txt <- readLines(MSG_PATH, warn = FALSE)
                if (length(.msg_txt) > 0)
                    cat(c("", strrep("=", 76),
                          "  MESSAGE STREAM (warnings, messages, errors)",
                          strrep("=", 76), .msg_txt, ""),
                        sep = "\n", file = LOG_PATH, append = TRUE)
                unlink(MSG_PATH)
            }
        }, error = function(e) invisible(NULL))
    }, add = TRUE)

    # warn = 1 prints warnings as they happen, so they reach the log. Restored on exit.
    .p9_old_warn <- options(warn = 1)
    on.exit(options(.p9_old_warn), add = TRUE)

    # Replay what was printed before the log existed.
    if (length(.p9_preamble) > 0)
        cat(paste(.p9_preamble, collapse = "\n"), "\n", sep = "")

    cat("\n", strrep("=", 76), "\n", sep = "")
    cat("  Part 9 -- QC, heatmaps and GAM-REML-PRSS  [", NORM_TAG, "]\n", sep = "")
    cat("  Started    : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", sep = "")
    cat("  Arm        : ", ARM_TAG, "\n", sep = "")
    cat("  Supertype  : ", PRIMARY_STYPE, "\n", sep = "")
    cat("  GOI        : ", GOI, "\n", sep = "")
    cat("  DFG (", length(DFG), ")    : ", paste(DFG, collapse = ", "), "\n", sep = "")
    cat("  Composite  : ", DFG_COMPOSITE_LABEL, "\n", sep = "")
    cat("  HKG (", length(HKG_USE), ")    : ",
        if (length(HKG_USE)) paste(HKG_USE, collapse = ", ")
        else "none -- control arm, panel only", "\n", sep = "")
    cat("  Norm       : ", NORM_TAG, " -- log1p(UMIs/total x ",
        format(NORM_TARGET_SUM, scientific = FALSE), "), ceiling ",
        sprintf("%.4f", NORM_CEILING), "\n", sep = "")
    cat("  Export tree: ", R_EXPORT_ROOT, "  (read; Part 9.2 also writes)\n",
        sep = "")
    cat("  Output dir : ", PART9_DIR, "\n", sep = "")
    cat("  Log file   : ", LOG_PATH,
        if (isTRUE(LOG_MESSAGES)) "  (console output and warnings)"
        else "  (console output only)", "\n", sep = "")
    cat(strrep("=", 76), "\n\n", sep = "")


    # =============================================================================
    # SECTION 0b: Shared check ledger
    # =============================================================================
    # One ledger for the QC (QC-nn) and GAM (GR-nn) checks.

    .chk_log  <- list()
    .chk_pass <- 0L
    .chk_warn <- 0L
    .chk_fail <- 0L

    log_result <- function(check_id, description, status,
                           detail = "", expected = NA, observed = NA) {
        symbol <- switch(status, PASS = "[PASS]", WARN = "[WARN]", FAIL = "[FAIL]")
        cat(sprintf("  %s  %-10s %s\n", symbol, check_id, description))
        if (nchar(detail) > 0)
            cat(sprintf("             detail   : %s\n", detail))
        if (!is.na(expected))
            cat(sprintf("             expected = %s   observed = %s\n",
                        expected, observed))
        .chk_log[[length(.chk_log) + 1L]] <<- data.frame(
            Check_ID    = check_id,
            Description = description,
            Status      = status,
            Detail      = detail,
            Expected    = as.character(expected),
            Observed    = as.character(observed),
            stringsAsFactors = FALSE
        )
        if (status == "PASS") .chk_pass <<- .chk_pass + 1L
        if (status == "WARN") .chk_warn <<- .chk_warn + 1L
        if (status == "FAIL") .chk_fail <<- .chk_fail + 1L
        invisible(NULL)
    }

    # Writes the ledger to text and to a workbook; called after QC and at the end.
    write_guardrail_outputs <- function() {
        con <- file(GR_TXT_PATH, "wt")
        writeLines(c(
            strrep("=", 76),
            paste0("Part 9 check ledger -- ", NORM_TAG),
            paste("Date        :", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
            paste("Arm         :", ARM_TAG),
            paste("Export root :", R_EXPORT_ROOT),
            paste0("Normalisation: adata.X = log1p(UMIs/total x ",
                   format(NORM_TARGET_SUM, scientific = FALSE),
                   "), exported unchanged"),
            strrep("=", 76), ""
        ), con)
        for (item in .chk_log)
            writeLines(sprintf("  [%s]  %-10s  %s%s",
                               item$Status, item$Check_ID, item$Description,
                               if (nchar(item$Detail) > 0)
                                   paste0("\n                          ", item$Detail)
                               else ""), con)
        writeLines(c("", strrep("=", 76),
                     sprintf("  Total : %d", .chk_pass + .chk_warn + .chk_fail),
                     sprintf("  PASS  : %d", .chk_pass),
                     sprintf("  WARN  : %d", .chk_warn),
                     sprintf("  FAIL  : %d", .chk_fail),
                     if (.chk_fail == 0) "  Result: no failures"
                     else "  Result: failures present",
                     strrep("=", 76)), con)
        close(con)

        chk_df <- do.call(rbind, .chk_log)
        rownames(chk_df) <- NULL
        write_xlsx(list(Check_Ledger = chk_df), path = GR_XLSX_PATH)

        cat("  Ledger written : ", GR_TXT_PATH,  "\n", sep = "")
        cat("  Ledger workbook: ", GR_XLSX_PATH, "\n", sep = "")
        invisible(chk_df)
    }


    # =============================================================================
    # SECTION 0c: Shared path and file helpers
    # =============================================================================

    # Dataset descriptors. The label is the sub-directory name used by the Part 8
    # exporter both inside the arm tree and under shared/.
    DATASETS <- list(
        list(tag = "NeoC",    label = "neocortex",
             dir = NEOC_DIR,  regions = NEOC_REGIONS,  dfg = DFG_NEOC,
             rec_field = "n_selected_neoc",
             plot_label = "Neocortex (DFC, MTG, AnG, STG, V1C, ITG)",
             file_label = "Neocortex"),
        list(tag = "NonNeoC", label = "nonneocortex",
             dir = NONNC_DIR, regions = NONNC_REGIONS, dfg = DFG_NONNC,
             rec_field = "n_selected_nonneoc",
             plot_label = "Non-neocortex (HIP, MEC, LEC, FI)",
             file_label = "NonNeocortex")
    )

    # Per-region GAM input file, matching the exporter's naming rule
    # TARGET_STYPE.replace('-','_').replace('/','_') + "_selected.csv".
    gam_file_path <- function(ds_dir, region, supertype = PRIMARY_STYPE) {
        file.path(ds_dir, "gam_data",
                  paste0(region, "_",
                         gsub("-", "_", gsub("/", "_", supertype)),
                         "_selected.csv"))
    }

    # First region with a per-region CSV; regions under 5 cells write none.
    first_available_region <- function(ds) {
        for (r in ds$regions)
            if (file.exists(gam_file_path(ds$dir, r))) return(r)
        NA_character_
    }

    read_gam_summary <- function(ds_dir) {
        p <- file.path(ds_dir, "gam_data", "summary_gam_inputs.csv")
        if (!file.exists(p)) return(NULL)
        read.csv(p, stringsAsFactors = FALSE, check.names = FALSE)
    }

    # Expected per-region selected counts, from the exporter's summary.
    expected_region_counts <- function(ds_dir) {
        summ <- read_gam_summary(ds_dir)
        if (is.null(summ) || !all(c("region", "n_selected") %in% colnames(summ)))
            return(NULL)
        setNames(as.integer(summ$n_selected), summ$region)
    }

    # Part 8 export record (already validated by R_config.R).
    EXPORT_RECORD <- local({
        p <- file.path(ARM_DIR, "Part8_export_record.csv")
        if (!file.exists(p)) NULL
        else read.csv(p, stringsAsFactors = FALSE, check.names = FALSE)
    })

    rec_field <- function(field, default = NA_character_) {
        if (is.null(EXPORT_RECORD)) return(default)
        v <- EXPORT_RECORD$value[EXPORT_RECORD$field == field]
        if (!length(v)) return(default)
        as.character(v[1])
    }

    USE_M4_TOP <- identical(toupper(rec_field("use_m4_top", "TRUE")), "TRUE")

    # Safeguard against a file with lone CR line endings.
    read_lines_clean <- function(path)
        sub("\r$", "", readLines(path, warn = FALSE))

    count_lines <- function(path) length(read_lines_clean(path))

    # Feature files hold gene symbols; quoting and comment characters are disabled.
    read_features <- function(path)
        read.table(path, sep = "\t", header = FALSE, quote = "",
                   comment.char = "", stringsAsFactors = FALSE)

    # Row count of a wide CSV without parsing every column.
    count_csv_rows <- function(path) {
        hdr <- read.csv(path, nrows = 1, stringsAsFactors = FALSE,
                        check.names = FALSE)
        cc  <- c("character", rep("NULL", ncol(hdr) - 1L))
        nrow(read.csv(path, colClasses = cc, stringsAsFactors = FALSE,
                      check.names = FALSE))
    }

    # Dimension line of a MatrixMarket file; returns c(nrow, ncol, nnz).
    read_mtx_header <- function(path) {
        con  <- file(path, "r")
        on.exit(close(con), add = TRUE)
        dims <- NULL
        repeat {
            ln <- readLines(con, n = 1, warn = FALSE)
            if (length(ln) == 0) break
            ln <- trimws(ln)
            if (nchar(ln) == 0) next
            if (startsWith(ln, "%")) next
            dims <- as.integer(strsplit(ln, "\\s+")[[1]])
            break
        }
        dims
    }

    # Loads a normalised matrix directory into a cell x gene sparse matrix with
    # dimnames set. Used by the QC, patch, heatmap and GAM sections.
    load_norm_matrix <- function(ds_dir, as_csparse = TRUE) {
        mtx_path  <- file.path(ds_dir, "deg_support", "norm.mtx")
        bar_path  <- file.path(ds_dir, "deg_support", "norm_barcodes.tsv")
        feat_path <- file.path(ds_dir, "deg_support", "norm_features.tsv")
        for (p in c(mtx_path, bar_path, feat_path))
            if (!file.exists(p)) stop("File not found: ", p)

        X        <- readMM(mtx_path)
        if (as_csparse) X <- as(X, "CsparseMatrix")
        barcodes <- read_lines_clean(bar_path)
        features <- read_features(feat_path)
        genes    <- features[, 2]
        rownames(X) <- barcodes
        colnames(X) <- genes
        list(X = X, barcodes = barcodes, genes = genes)
    }


    # =============================================================================
    # PART 9.1 -- EXPORT QC
    # =============================================================================
    # QC-01  R_config.R integrity; arm tree and shared tree present
    # QC-02  Shared files; donor_clinical carries the tissue-quality covariates
    # QC-03  GAM CSV structure: base columns, DFG columns, HKG columns
    # QC-04  GOI and composite value ranges within the CP10K ceiling
    # QC-05  Per-region selected counts vs summary_gam_inputs.csv
    # QC-06  MTX dimensions vs barcode and feature line counts
    # QC-07  GAM barcodes are a subset of the selected-cell matrix barcodes
    # QC-08  norm_GOI_DFG_selected.csv agrees with the GAM CSV GOI column
    # QC-09  Clinical metadata completeness in obs_selected.csv
    # QC-10  Dataset and selected cell counts vs the Part 8 export record
    # QC-11  Scorer columns present and finite (M4-top arms only)
    # QC-12  Normalisation guardrails, positive form
    # QC-13  QC summary
    # =============================================================================

    cat("\n", strrep("#", 76), "\n", sep = "")
    cat("  PART 9.1 -- EXPORT QC\n")
    cat(strrep("#", 76), "\n\n", sep = "")


    # =============================================================================
    # QC-01  R_config.R integrity
    # =============================================================================
    cat("-- QC-01  R_config.R integrity ------------------------------------------\n")

    required_vars <- c("R_EXPORT_ROOT", "ARM_TAG", "ARM_DIR",
                       "NEOC_DIR", "NONNC_DIR", "SHARED_DIR",
                       "GOI", "DFG", "DFG_NEOC", "DFG_NONNC", "SELECT_GENES",
                       "MOD1_GENES", "HKG_GENES",
                       "DFG_COMPOSITE_LABEL", "RUN_SUFFIX",
                       "NEOC_REGIONS", "NONNC_REGIONS", "PRIMARY_STYPE",
                       "NORM_TARGET_SUM", "NORM_CEILING", "NORM_TAG",
                       "REGION_COL", "STYPE_COL", "DONOR_COL",
                       "ADNC_COL", "BRAAK_COL",
                       "GAM_GAMMA", "GAM_METHOD", "GAM_SELECT", "GAM_BS",
                       "GAM_MAX_K", "PRSS_N_ITER")

    missing_vars <- required_vars[!sapply(required_vars, exists)]
    if (length(missing_vars) == 0) {
        log_result("QC-01a", "All required config variables defined", "PASS",
                   detail = paste(length(required_vars), "variables checked"))
    } else {
        log_result("QC-01a", "Required config variables missing", "FAIL",
                   detail = paste(missing_vars, collapse = ", "))
    }

    required_dirs <- c(R_EXPORT_ROOT, ARM_DIR, SHARED_DIR)
    for (ds in DATASETS) {
        required_dirs <- c(required_dirs,
            ds$dir,
            file.path(ds$dir, "gam_data"),
            file.path(ds$dir, "counts", "raw_selected"),
            file.path(ds$dir, "deg_support"),
            file.path(ds$dir, "metadata"),
            file.path(ds$dir, "embeddings"),
            file.path(SHARED_DIR, ds$label, "counts", "raw_dataset"))
    }

    missing_dirs <- required_dirs[!dir.exists(required_dirs)]
    if (length(missing_dirs) == 0) {
        log_result("QC-01b", "All required directories exist", "PASS",
                   detail = paste(length(required_dirs), "directories checked"))
    } else {
        log_result("QC-01b", "Missing directories", "FAIL",
                   detail = paste(missing_dirs, collapse = " | "))
    }

    # The export record is required for the count checks below.
    if (is.null(EXPORT_RECORD)) {
        log_result("QC-01c", "Part8_export_record.csv present", "FAIL",
                   detail = paste("not found in", ARM_DIR))
    } else {
        log_result("QC-01c", "Part8_export_record.csv present", "PASS",
                   detail = sprintf("h5ad = %s | release = %s | exported %s",
                                    rec_field("h5ad"), rec_field("release"),
                                    rec_field("exported_utc")))
        log_result("QC-01d", "Export record normalisation matches config",
                   if (identical(rec_field("norm_tag"), NORM_TAG) &&
                       identical(rec_field("norm_target_sum"),
                                 format(NORM_TARGET_SUM, scientific = FALSE)))
                       "PASS" else "FAIL",
                   detail = sprintf("record: %s at target_sum %s | config: %s at %s",
                                    rec_field("norm_tag"),
                                    rec_field("norm_target_sum"), NORM_TAG,
                                    format(NORM_TARGET_SUM, scientific = FALSE)))
    }


    # =============================================================================
    # QC-02  Shared files
    # =============================================================================
    cat("\n-- QC-02  Shared files --------------------------------------------------\n")

    shared_files <- c("gene_metadata.csv", "donor_clinical.csv", "obs_all_cells.csv",
                      "embedding_scVI.csv", "embedding_umap.csv")

    for (fn in shared_files) {
        fp <- file.path(SHARED_DIR, fn)
        if (!file.exists(fp)) {
            log_result("QC-02", paste("Missing:", fn), "FAIL")
        } else {
            log_result("QC-02", paste("File present:", fn), "PASS",
                       detail = paste0(round(file.size(fp) / 1e6, 1), " MB"))
        }
    }

    if (file.exists(file.path(SHARED_DIR, "gene_metadata.csv"))) {
        gene_meta <- read.csv(file.path(SHARED_DIR, "gene_metadata.csv"),
                              stringsAsFactors = FALSE, check.names = FALSE)
        log_result("QC-02", "gene_metadata row count",
                   if (nrow(gene_meta) == 36601) "PASS" else "WARN",
                   expected = 36601, observed = nrow(gene_meta))
    }

    if (file.exists(file.path(SHARED_DIR, "donor_clinical.csv"))) {
        donor_clin <- read.csv(file.path(SHARED_DIR, "donor_clinical.csv"),
                               stringsAsFactors = FALSE, check.names = FALSE)
        log_result("QC-02", "donor_clinical row count",
                   if (nrow(donor_clin) == 84) "PASS" else "WARN",
                   expected = 84, observed = nrow(donor_clin))

        # The table carries the tissue-quality covariates and no region column.
        quality_cols <- c("PMI", "RIN", "Brain pH", "Fresh Brain Weight",
                          "Severely Affected Donor")
        q_missing <- setdiff(quality_cols, colnames(donor_clin))
        log_result("QC-02", "donor_clinical tissue-quality covariates",
                   if (length(q_missing) == 0) "PASS" else "WARN",
                   detail = if (length(q_missing) == 0)
                       paste("all present:", paste(quality_cols, collapse = ", "))
                   else paste("missing:", paste(q_missing, collapse = ", ")))

        log_result("QC-02", "donor_clinical has no region column",
                   if (!("region" %in% colnames(donor_clin))) "PASS" else "WARN",
                   detail = if ("region" %in% colnames(donor_clin))
                       "a region column is present; one arbitrary region per donor"
                   else "region correctly absent")
    }

    if (file.exists(file.path(SHARED_DIR, "obs_all_cells.csv"))) {
        obs_all_hdr <- read.csv(file.path(SHARED_DIR, "obs_all_cells.csv"),
                                nrows = 5, stringsAsFactors = FALSE,
                                check.names = FALSE)
        log_result("QC-02", "obs_all_cells has a barcode column",
                   if ("barcode" %in% colnames(obs_all_hdr)) "PASS" else "FAIL",
                   detail = paste("first 10 cols:",
                                  paste(head(colnames(obs_all_hdr), 10),
                                        collapse = ", ")))
    }


    # =============================================================================
    # QC-03  GAM CSV structure
    # =============================================================================
    cat("\n-- QC-03  GAM CSV column structure --------------------------------------\n")

    required_gam_cols <- c("barcode", "region", "supertype", "donor",
                           "GOI", "Expression")

    check_gam_structure <- function(ds) {
        summ <- read_gam_summary(ds$dir)
        if (is.null(summ)) {
            log_result("QC-03", paste(ds$tag, "summary_gam_inputs.csv missing"),
                       "FAIL")
            return(invisible(NULL))
        }
        log_result("QC-03", paste(ds$tag, "summary_gam_inputs.csv rows"),
                   if (nrow(summ) > 0) "PASS" else "FAIL",
                   detail = paste(nrow(summ), "region x supertype combinations"))

        # The CP10K exporter writes "file" and "n_selected".
        if (!("file" %in% colnames(summ))) {
            log_result("QC-03", paste(ds$tag, "summary_gam_inputs.csv schema"),
                       "FAIL",
                       detail = paste("expected a 'file' column; found:",
                                      paste(colnames(summ), collapse = ", ")))
            return(invisible(NULL))
        }

        primary_rows <- summ[summ$supertype == PRIMARY_STYPE, , drop = FALSE]
        if (nrow(primary_rows) == 0) {
            log_result("QC-03", paste(ds$tag, "rows for", PRIMARY_STYPE), "FAIL",
                       detail = paste("supertypes present:",
                                      paste(unique(summ$supertype), collapse = ", ")))
            return(invisible(NULL))
        }

        for (i in seq_len(nrow(primary_rows))) {
            fpath <- file.path(ds$dir, "gam_data", primary_rows$file[i])
            if (!file.exists(fpath)) {
                log_result("QC-03",
                           paste(ds$tag, primary_rows$region[i], "selected CSV"),
                           "FAIL",
                           detail = paste("file not found:", basename(fpath)))
                next
            }
            df       <- read.csv(fpath, stringsAsFactors = FALSE,
                                 check.names = FALSE)
            miss     <- setdiff(required_gam_cols, colnames(df))
            dfg_miss <- setdiff(ds$dfg,  colnames(df))
            hkg_miss <- setdiff(HKG_USE, colnames(df))

            status <- if (length(miss) > 0) "FAIL"
                      else if (length(dfg_miss) > 0 || length(hkg_miss) > 0) "WARN"
                      else "PASS"

            detail_parts <- character(0)
            if (length(miss) > 0) {
                detail_parts <- paste("missing base cols:",
                                      paste(miss, collapse = ","))
            } else {
                if (length(dfg_miss) > 0)
                    detail_parts <- c(detail_parts,
                        paste("DFG missing:", paste(dfg_miss, collapse = ",")))
                if (length(hkg_miss) > 0)
                    detail_parts <- c(detail_parts,
                        paste("HKG missing:", paste(hkg_miss, collapse = ",")))
                if (length(detail_parts) > 0)
                    detail_parts <- c(detail_parts,
                                      "(Part 9.2 adds these from norm.mtx)")
                else
                    detail_parts <- paste("all required columns present; total cols =",
                                          ncol(df))
            }
            log_result("QC-03",
                       paste(ds$tag, primary_rows$region[i],
                             "column set, n =", nrow(df)),
                       status, detail = paste(detail_parts, collapse = " | "))
        }
    }

    for (ds in DATASETS) check_gam_structure(ds)


    # =============================================================================
    # QC-04  GOI and composite value ranges
    # =============================================================================
    cat("\n-- QC-04  GOI and composite value ranges --------------------------------\n")

    check_value_ranges <- function(ds) {
        goi_selected <- GOI %in% SELECT_GENES
        for (reg in ds$regions) {
            fpath <- gam_file_path(ds$dir, reg)
            if (!file.exists(fpath)) next
            df <- read.csv(fpath, stringsAsFactors = FALSE, check.names = FALSE)

            goi_range  <- range(df$GOI,        na.rm = TRUE)
            expr_range <- range(df$Expression, na.rm = TRUE)
            n_na_goi   <- sum(is.na(df$GOI))
            n_na_expr  <- sum(is.na(df$Expression))
            n_neg_goi  <- sum(df$GOI        < 0, na.rm = TRUE)
            n_neg_expr <- sum(df$Expression < 0, na.rm = TRUE)

            # Upper bound log1p(NORM_TARGET_SUM); no lower bound.
            goi_ok  <- goi_range[1]  >= 0 &&
                       goi_range[2]  <= NORM_CEILING + NORM_TOL && n_na_goi  == 0
            expr_ok <- expr_range[1] >= 0 &&
                       expr_range[2] <= NORM_CEILING + NORM_TOL && n_na_expr == 0

            status <- if (goi_ok && expr_ok) "PASS"
                      else if (n_neg_goi > 0 || n_neg_expr > 0) "FAIL"
                      else "WARN"

            log_result("QC-04", paste(ds$tag, reg, "GOI / composite ranges"), status,
                       detail = sprintf(
                           "GOI [%.3f, %.3f] NA=%d neg=%d | %s [%.3f, %.3f] NA=%d neg=%d | ceiling %.4f",
                           goi_range[1], goi_range[2], n_na_goi, n_neg_goi,
                           DFG_COMPOSITE_LABEL,
                           expr_range[1], expr_range[2], n_na_expr, n_neg_expr,
                           NORM_CEILING))

            # Part 8 selected on GOI > 0, so a non-positive GOI value is flagged.
            if (goi_selected) {
                n_nonpos <- sum(df$GOI <= 0, na.rm = TRUE)
                log_result("QC-04b", paste(ds$tag, reg, "GOI > 0 as selected"),
                           if (n_nonpos == 0) "PASS" else "FAIL",
                           detail = sprintf("%d of %d cells have GOI <= 0 (%s in SELECT_GENES)",
                                            n_nonpos, nrow(df), GOI))
            }
        }
    }

    for (ds in DATASETS) check_value_ranges(ds)


    # =============================================================================
    # QC-05  Per-region selected counts vs the exporter's summary
    # =============================================================================
    cat("\n-- QC-05  Per-region selected counts ------------------------------------\n")

    check_region_counts <- function(ds) {
        exp_counts <- expected_region_counts(ds$dir)
        if (is.null(exp_counts)) {
            log_result("QC-05", paste(ds$tag, "expected counts unavailable"), "FAIL",
                       detail = "summary_gam_inputs.csv missing or lacks n_selected")
            return(invisible(NULL))
        }
        for (reg in ds$regions) {
            fpath <- gam_file_path(ds$dir, reg)
            if (!file.exists(fpath)) {
                # Regions below the exporter's 5-cell floor have no CSV and no summary row.
                status <- if (reg %in% names(exp_counts)) "FAIL" else "WARN"
                log_result("QC-05", paste(ds$tag, reg, "selected CSV"), status,
                           detail = if (status == "FAIL")
                               "summary lists this region but the file is absent"
                           else "no CSV and no summary row; region below the export floor")
                next
            }
            obs <- nrow(read.csv(fpath, stringsAsFactors = FALSE))
            if (!(reg %in% names(exp_counts))) {
                log_result("QC-05", paste(ds$tag, reg, "selected count"), "WARN",
                           detail = paste("CSV present with", obs,
                                          "rows but no summary row"))
                next
            }
            exp <- exp_counts[[reg]]
            log_result("QC-05", paste(ds$tag, reg, "selected count"),
                       if (obs == exp) "PASS" else "FAIL",
                       expected = exp, observed = obs)
        }
    }

    for (ds in DATASETS) check_region_counts(ds)


    # =============================================================================
    # QC-06  MTX dimensions
    # =============================================================================
    cat("\n-- QC-06  MTX dimensions ------------------------------------------------\n")

    check_mtx <- function(label, mtx_dir, prefix = "matrix") {
        mtx_path  <- file.path(mtx_dir, paste0(prefix, ".mtx"))
        bar_path  <- file.path(mtx_dir, paste0(prefix, "_barcodes.tsv"))
        feat_path <- file.path(mtx_dir, paste0(prefix, "_features.tsv"))

        if (!all(file.exists(c(mtx_path, bar_path, feat_path)))) {
            log_result("QC-06", paste(label, "MTX file set"), "FAIL",
                       detail = paste0("expected prefix '", prefix, "'; missing: ",
                                       paste(basename(c(mtx_path, bar_path,
                                                        feat_path))[
                                           !file.exists(c(mtx_path, bar_path,
                                                          feat_path))],
                                           collapse = ", ")))
            return(invisible(NULL))
        }

        dims <- read_mtx_header(mtx_path)
        if (is.null(dims) || length(dims) < 2L || anyNA(dims[1:2])) {
            # An unreadable header is caught before the comparison.
            log_result("QC-06", paste(label, "MTX header"), "FAIL",
                       detail = paste("no dimension line found in",
                                      basename(mtx_path)))
            return(invisible(NULL))
        }
        n_rows_mtx <- dims[1]
        n_cols_mtx <- dims[2]
        n_barcodes <- count_lines(bar_path)
        n_features <- count_lines(feat_path)

        status <- if (n_rows_mtx == n_barcodes && n_cols_mtx == n_features)
            "PASS" else "FAIL"

        log_result("QC-06", paste(label, "MTX dims vs barcodes and features"),
                   status,
                   detail = sprintf("MTX %d x %d | barcodes %d | features %d",
                                    n_rows_mtx, n_cols_mtx, n_barcodes, n_features))
    }

    for (ds in DATASETS) {
        check_mtx(paste(ds$tag, "raw_selected"),
                  file.path(ds$dir, "counts", "raw_selected"), prefix = "matrix")
        check_mtx(paste(ds$tag, "raw_dataset (shared)"),
                  file.path(SHARED_DIR, ds$label, "counts", "raw_dataset"),
                  prefix = "matrix")
        check_mtx(paste(ds$tag, "norm"),
                  file.path(ds$dir, "deg_support"), prefix = "norm")
    }


    # =============================================================================
    # QC-07  Barcode consistency
    # =============================================================================
    cat("\n-- QC-07  Barcode consistency -------------------------------------------\n")

    check_barcode_subset <- function(ds) {
        bar_path <- file.path(ds$dir, "counts", "raw_selected",
                              "matrix_barcodes.tsv")
        if (!file.exists(bar_path)) {
            log_result("QC-07", paste(ds$tag, "raw_selected barcodes"), "FAIL",
                       detail = paste("not found:", bar_path))
            return(invisible(NULL))
        }
        mtx_set <- read_lines_clean(bar_path)

        for (reg in ds$regions) {
            fpath <- gam_file_path(ds$dir, reg)
            if (!file.exists(fpath)) next
            gam_bc   <- read.csv(fpath, stringsAsFactors = FALSE,
                                 check.names = FALSE)$barcode
            n_total  <- length(gam_bc)
            n_in_mtx <- sum(gam_bc %in% mtx_set)
            log_result("QC-07", paste(ds$tag, reg, "GAM barcodes in count matrix"),
                       if (n_in_mtx == n_total) "PASS" else "FAIL",
                       detail = sprintf("%d / %d found", n_in_mtx, n_total))
        }
    }

    for (ds in DATASETS) check_barcode_subset(ds)


    # =============================================================================
    # QC-08  Dense GOI/DFG CSV against the GAM CSV
    # =============================================================================
    cat("\n-- QC-08  norm_GOI_DFG_selected.csv spot-check --------------------------\n")

    check_norm_csv <- function(ds) {
        norm_path <- file.path(ds$dir, "deg_support", "norm_GOI_DFG_selected.csv")
        spot_reg  <- first_available_region(ds)
        if (is.na(spot_reg)) {
            log_result("QC-08", paste(ds$tag, "no per-region CSV to check against"),
                       "FAIL")
            return(invisible(NULL))
        }
        gam_path <- gam_file_path(ds$dir, spot_reg)
        if (!file.exists(norm_path) || !file.exists(gam_path)) {
            log_result("QC-08", paste(ds$tag, "norm and GAM files present"), "FAIL",
                       detail = paste("checked:", basename(norm_path), "and",
                                      basename(gam_path)))
            return(invisible(NULL))
        }

        norm_df   <- read.csv(norm_path, stringsAsFactors = FALSE,
                              check.names = FALSE)
        gam_df    <- read.csv(gam_path,  stringsAsFactors = FALSE,
                              check.names = FALSE)
        common_bc <- intersect(norm_df$barcode, gam_df$barcode)

        if (length(common_bc) < 10) {
            log_result("QC-08", paste(ds$tag, "barcode overlap"), "WARN",
                       detail = paste("only", length(common_bc),
                                      "common barcodes"))
            return(invisible(NULL))
        }

        norm_goi <- as.numeric(norm_df[match(common_bc, norm_df$barcode), GOI])
        gam_goi  <- as.numeric(gam_df[ match(common_bc, gam_df$barcode),  "GOI"])
        r        <- suppressWarnings(cor(norm_goi, gam_goi, use = "complete.obs"))
        max_diff <- max(abs(norm_goi - gam_goi), na.rm = TRUE)

        status <- if (!is.na(r) && r > 0.999 && max_diff < NORM_CSV_TOL) "PASS"
                  else if (!is.na(r) && r > 0.99)                        "WARN"
                  else                                                   "FAIL"

        log_result("QC-08", paste(ds$tag, "GOI in dense CSV vs GAM CSV"), status,
                   detail = sprintf("r = %.6f | max_diff = %.2e | n = %d | region %s",
                                    r, max_diff, length(common_bc), spot_reg))
    }

    for (ds in DATASETS) check_norm_csv(ds)


    # =============================================================================
    # QC-09  Clinical metadata completeness
    # =============================================================================
    cat("\n-- QC-09  Clinical metadata completeness --------------------------------\n")

    check_clinical <- function(ds) {
        meta_path <- file.path(ds$dir, "metadata", "obs_selected.csv")
        if (!file.exists(meta_path)) {
            log_result("QC-09", paste(ds$tag, "obs_selected.csv missing"), "FAIL")
            return(invisible(NULL))
        }
        meta <- read.csv(meta_path, stringsAsFactors = FALSE, check.names = FALSE)

        clinical_cols <- c(ADNC_COL, BRAAK_COL, "CERAD score", "Thal",
                           "Cognitive Status", "APOE Genotype", DONOR_COL)

        for (col in clinical_cols) {
            if (!(col %in% colnames(meta))) {
                similar <- grep(gsub(" ", ".", col, fixed = TRUE),
                                colnames(meta), value = TRUE, ignore.case = TRUE)
                log_result("QC-09", paste(ds$tag, shQuote(col), "column absent"),
                           "WARN",
                           detail = if (length(similar) > 0)
                               paste("closest match:", paste(similar, collapse = ", "))
                           else paste("total cols in file:", ncol(meta)))
                next
            }
            n_na <- sum(is.na(meta[[col]]) | meta[[col]] == "")
            pct  <- round(100 * n_na / nrow(meta), 1)
            log_result("QC-09", paste(ds$tag, shQuote(col), "NA rate"),
                       if (pct == 0) "PASS" else if (pct < 5) "WARN" else "FAIL",
                       detail = sprintf("%d NA or empty out of %d cells (%.1f%%)",
                                        n_na, nrow(meta), pct))
        }
    }

    for (ds in DATASETS) check_clinical(ds)


    # =============================================================================
    # QC-10  Cell count audit against the export record
    # =============================================================================
    cat("\n-- QC-10  Cell count audit ----------------------------------------------\n")

    check_cell_counts <- function(ds) {
        sel_path <- file.path(ds$dir, "metadata", "obs_selected.csv")
        ds_path  <- file.path(ds$dir, "metadata", "obs_dataset.csv")

        exp_sel <- suppressWarnings(as.integer(rec_field(ds$rec_field)))

        if (file.exists(sel_path)) {
            n_sel <- count_csv_rows(sel_path)
            if (is.na(exp_sel)) {
                log_result("QC-10", paste(ds$tag, "selected cell count"), "WARN",
                           detail = paste("observed", n_sel,
                                          "; export record has no",
                                          ds$rec_field))
            } else {
                log_result("QC-10", paste(ds$tag, "selected cell count"),
                           if (n_sel == exp_sel) "PASS" else "FAIL",
                           expected = exp_sel, observed = n_sel)
            }
        } else {
            log_result("QC-10", paste(ds$tag, "obs_selected.csv missing"), "FAIL")
        }

        # Dataset-wide obs rows are compared with the shared raw_dataset barcodes.
        if (file.exists(ds_path)) {
            n_ds <- count_csv_rows(ds_path)
            shared_bar <- file.path(SHARED_DIR, ds$label, "counts", "raw_dataset",
                                    "matrix_barcodes.tsv")
            if (file.exists(shared_bar)) {
                n_shared <- count_lines(shared_bar)
                log_result("QC-10", paste(ds$tag, "dataset cell count"),
                           if (n_ds == n_shared) "PASS" else "FAIL",
                           expected = n_shared, observed = n_ds)
            } else {
                log_result("QC-10", paste(ds$tag, "dataset cell count"), "WARN",
                           detail = paste("observed", n_ds,
                                          "; shared raw_dataset barcodes not found"))
            }
        } else {
            log_result("QC-10", paste(ds$tag, "obs_dataset.csv missing"), "FAIL")
        }
    }

    for (ds in DATASETS) check_cell_counts(ds)


    # =============================================================================
    # QC-11  Scorer columns
    # =============================================================================
    cat("\n-- QC-11  Scorer columns ------------------------------------------------\n")

    check_scores <- function(ds) {
        sc_path <- file.path(ds$dir, "metadata", "scores_selected_all_methods.csv")

        # Part 8 writes the scorer table only for the M4-top arms.
        if (!file.exists(sc_path)) {
            log_result("QC-11", paste(ds$tag, "scores_selected_all_methods.csv"),
                       if (USE_M4_TOP) "FAIL" else "PASS",
                       detail = if (USE_M4_TOP)
                           "M4-top arm: the scorer table should be present"
                       else "all-cells arm: no scorer table expected")
            return(invisible(NULL))
        }

        sc_df      <- read.csv(sc_path, stringsAsFactors = FALSE,
                               check.names = FALSE)
        score_cols <- grep("^score_", colnames(sc_df), value = TRUE)

        if (length(score_cols) == 0) {
            log_result("QC-11", paste(ds$tag, "no score_* columns found"), "FAIL")
            return(invisible(NULL))
        }

        for (col in score_cols) {
            n_na  <- sum(is.na(sc_df[[col]]))
            n_inf <- sum(!is.finite(sc_df[[col]]), na.rm = TRUE)
            sc_r  <- range(sc_df[[col]], na.rm = TRUE)
            log_result("QC-11", paste(ds$tag, col),
                       if (n_na == 0 && n_inf == 0) "PASS" else "WARN",
                       detail = sprintf("range [%.4f, %.4f] NA=%d Inf=%d",
                                        sc_r[1], sc_r[2], n_na, n_inf))
        }
    }

    for (ds in DATASETS) check_scores(ds)


    # =============================================================================
    # QC-12  Normalisation guardrails, positive form
    # =============================================================================
    # norm.mtx spans the full gene space, so these tests are exact.
    cat("\n-- QC-12  Normalisation guardrails --------------------------------------\n")

    check_normalisation <- function(ds) {
        norm_dir <- file.path(ds$dir, "deg_support")
        raw_dir  <- file.path(ds$dir, "counts", "raw_selected")

        required <- c(file.path(norm_dir, c("norm.mtx", "norm_barcodes.tsv",
                                            "norm_features.tsv")),
                      file.path(raw_dir,  c("matrix.mtx", "matrix_barcodes.tsv",
                                            "matrix_features.tsv")))
        if (!all(file.exists(required))) {
            log_result("QC-12", paste(ds$tag, "norm and raw files present"), "FAIL",
                       detail = paste("missing:",
                                      paste(basename(required[!file.exists(required)]),
                                            collapse = ", ")))
            return(invisible(NULL))
        }

        nm       <- load_norm_matrix(ds$dir, as_csparse = TRUE)
        X_norm   <- nm$X
        norm_bc  <- nm$barcodes
        norm_gen <- nm$genes

        # == Maximum within the CP10K ceiling =================================
        norm_max <- max(X_norm)
        log_result("QC-12a", paste(ds$tag, "norm.mtx maximum within the ceiling"),
                   if (norm_max > 0 && norm_max <= NORM_CEILING + NORM_TOL)
                       "PASS" else "FAIL",
                   detail = sprintf("max = %.4f (ceiling log1p(%s) = %.4f)",
                                    norm_max,
                                    format(NORM_TARGET_SUM, scientific = FALSE),
                                    NORM_CEILING))

        # == expm1 row sums pin the target ====================================
        # A second normalisation, another log base or a gene-subset export breaks this.
        # One row-block subset, then a sparse expm1.
        n_rows_check  <- min(50L, nrow(X_norm))
        row_sums_norm <- rowSums(expm1(X_norm[seq_len(n_rows_check), ,
                                              drop = FALSE]))
        rs_mean <- mean(row_sums_norm)
        rs_dev  <- max(abs(row_sums_norm - NORM_TARGET_SUM))
        log_result("QC-12b", paste(ds$tag, "expm1(norm) row sums pin the target"),
                   if (rs_dev < NORM_TARGET_SUM * 1e-6) "PASS" else "FAIL",
                   detail = sprintf("mean = %.6f, max |deviation| = %.2e, n = %d cells",
                                    rs_mean, rs_dev, n_rows_check))

        X_raw    <- readMM(file.path(raw_dir, "matrix.mtx"))
        X_raw    <- as(X_raw, "CsparseMatrix")
        raw_bc   <- read_lines_clean(file.path(raw_dir, "matrix_barcodes.tsv"))
        raw_feat <- read_features(file.path(raw_dir, "matrix_features.tsv"))
        raw_gen  <- raw_feat[, 2]
        rownames(X_raw) <- raw_bc
        colnames(X_raw) <- raw_gen

        # == Reconstruction from raw counts ===================================
        common_bc <- intersect(norm_bc, raw_bc)
        n_sample  <- min(50L, length(common_bc))

        if (n_sample < 10) {
            log_result("QC-12c", paste(ds$tag, "reconstruction from raw counts"),
                       "WARN",
                       detail = paste("only", length(common_bc),
                                      "common barcodes"))
        } else {
            sample_bc <- common_bc[seq_len(n_sample)]
            norm_idx  <- match(sample_bc, norm_bc)
            raw_idx   <- match(sample_bc, raw_bc)

            # The DFG set comes from the dataset descriptor.
            test_genes  <- unique(c(GOI, ds$dfg))
            random_pool <- setdiff(intersect(norm_gen, raw_gen), test_genes)
            set.seed(42)
            if (length(random_pool) >= 5)
                test_genes <- c(test_genes, sample(random_pool, 5))
            test_genes <- test_genes[test_genes %in% norm_gen &
                                     test_genes %in% raw_gen]

            if (length(test_genes) < 5) {
                log_result("QC-12c",
                           paste(ds$tag, "reconstruction from raw counts"), "FAIL",
                           detail = sprintf("only %d gene(s) in both feature spaces",
                                            length(test_genes)))
            } else {
                norm_gi <- match(test_genes, norm_gen)
                raw_gi  <- match(test_genes, raw_gen)

                norm_sub <- as.matrix(X_norm[norm_idx, norm_gi, drop = FALSE])
                raw_sub  <- as.matrix(X_raw[ raw_idx,  raw_gi,  drop = FALSE])

                # Totals over the full gene space, as in the release normalisation.
                row_totals <- rowSums(X_raw[raw_idx, , drop = FALSE])
                row_totals[row_totals == 0] <- 1
                recon <- log1p(sweep(raw_sub, 1, row_totals, "/") * NORM_TARGET_SUM)

                max_diff <- max(abs(norm_sub - recon))
                log_result("QC-12c",
                           paste(ds$tag,
                                 "norm.mtx reconstructs from raw counts"),
                           if (max_diff < NORM_TOL) "PASS" else "FAIL",
                           detail = sprintf(
                               "max_diff = %.2e vs log1p(UMIs/total x %s), %d genes x %d cells",
                               max_diff,
                               format(NORM_TARGET_SUM, scientific = FALSE),
                               length(test_genes), n_sample))
            }
        }

        # == norm.mtx GOI against the GAM CSV GOI =============================
        # Spot-check region from the dataset's own region list.
        spot_reg <- first_available_region(ds)
        gam_path <- if (is.na(spot_reg)) "" else gam_file_path(ds$dir, spot_reg)
        if (nzchar(gam_path) && file.exists(gam_path)) {
            gam_df  <- read.csv(gam_path, stringsAsFactors = FALSE,
                                check.names = FALSE)
            goi_col <- match(GOI, norm_gen)
            if (is.na(goi_col)) {
                log_result("QC-12d", paste(ds$tag, "GOI present in norm.mtx"),
                           "FAIL", detail = paste(GOI, "not in norm_features.tsv"))
            } else {
                bc_in_norm <- match(gam_df$barcode, norm_bc)
                valid      <- !is.na(bc_in_norm)
                if (sum(valid) >= 10) {
                    norm_goi <- as.numeric(X_norm[bc_in_norm[valid], goi_col])
                    gam_goi  <- gam_df$GOI[valid]
                    ic_diff  <- max(abs(norm_goi - gam_goi))
                    log_result("QC-12d",
                               paste(ds$tag, "norm.mtx GOI vs GAM CSV GOI"),
                               if (ic_diff < NORM_CSV_TOL) "PASS" else "FAIL",
                               detail = sprintf("max_diff = %.2e over %d cells, region %s",
                                                ic_diff, sum(valid), spot_reg))
                } else {
                    log_result("QC-12d",
                               paste(ds$tag, "norm.mtx GOI vs GAM CSV GOI"), "WARN",
                               detail = sprintf("only %d matched barcodes",
                                                sum(valid)))
                }
            }
        } else {
            log_result("QC-12d", paste(ds$tag, "GAM CSV for the spot-check region"),
                       "WARN", detail = paste("not found:", basename(gam_path)))
        }

        # == Raw counts are integers ==========================================
        raw_sample <- as.numeric(X_raw@x[seq_len(min(10000L, length(X_raw@x)))])
        is_int     <- all(abs(raw_sample - round(raw_sample)) < 0.01)
        raw_max    <- max(X_raw)
        log_result("QC-12e", paste(ds$tag, "raw counts are integer-valued"),
                   if (is_int && raw_max > 20) "PASS" else "FAIL",
                   detail = sprintf("max = %.0f, integer-like = %s",
                                    raw_max, if (is_int) "TRUE" else "FALSE"))

        rm(X_norm, X_raw, nm); gc(verbose = FALSE)
    }

    for (ds in DATASETS) check_normalisation(ds)


    # =============================================================================
    # QC-13  QC summary
    # =============================================================================
    cat("\n", strrep("=", 76), "\n", sep = "")
    cat("  QC SUMMARY\n")
    cat(strrep("=", 76), "\n\n", sep = "")

    cat(sprintf("  Total checks : %d\n",   .chk_pass + .chk_warn + .chk_fail))
    cat(sprintf("  PASS         : %d\n",   .chk_pass))
    cat(sprintf("  WARN         : %d\n",   .chk_warn))
    cat(sprintf("  FAIL         : %d\n\n", .chk_fail))

    .qc_df       <- do.call(rbind, .chk_log)
    .qc_non_pass <- .qc_df[.qc_df$Status != "PASS", , drop = FALSE]

    if (nrow(.qc_non_pass) > 0) {
        cat("  Non-PASS items:\n")
        for (i in seq_len(nrow(.qc_non_pass))) {
            cat(sprintf("    [%s] %s -- %s\n           %s\n",
                        .qc_non_pass$Status[i], .qc_non_pass$Check_ID[i],
                        .qc_non_pass$Description[i], .qc_non_pass$Detail[i]))
        }
        cat("\n")
    } else {
        cat("  No non-PASS items.\n\n")
    }

    QC_FAIL_COUNT <- .chk_fail
    write_guardrail_outputs()

    if (QC_FAIL_COUNT > 0) {
        cat("\n  QC failures are present.\n\n")
        warning("Part 9 QC: ", QC_FAIL_COUNT, " FAIL(s)", call. = FALSE)
    }


    # =============================================================================
    # PART 9.2 -- GENE COLUMN PATCH
    # =============================================================================
    # Appends any missing GOI, DFG or HKG column from deg_support/norm.mtx to the
    # per-region GAM CSVs, which Part 10 reads; present columns are left alone.
    # =============================================================================

    cat("\n", strrep("#", 76), "\n", sep = "")
    cat("  PART 9.2 -- GENE COLUMN PATCH\n")
    cat(strrep("#", 76), "\n\n", sep = "")

    cat("  Writes HKG columns into the GAM CSVs under ", R_EXPORT_ROOT, "\n\n", sep = "")

    PATCH_GENES <- unique(c(GOI, DFG, HKG_USE))
    cat("  Required columns:", paste(PATCH_GENES, collapse = ", "), "\n")
    cat(sprintf("  (1 GOI + %d DFG + %d HKG = %d unique)\n\n",
                length(DFG), length(HKG_USE), length(PATCH_GENES)))

    patch_gam_csvs <- function(ds) {
        cat(sprintf("  === %s ===\n", ds$tag))

        norm_paths <- file.path(ds$dir, "deg_support",
                                c("norm.mtx", "norm_barcodes.tsv",
                                  "norm_features.tsv"))
        if (!all(file.exists(norm_paths))) {
            cat("    SKIP: norm.mtx file set not found\n\n")
            return(invisible(NULL))
        }

        nm       <- load_norm_matrix(ds$dir, as_csparse = TRUE)
        X_norm   <- nm$X
        norm_bc  <- nm$barcodes
        norm_gen <- nm$genes

        cat(sprintf("    norm.mtx: %d cells x %d genes\n",
                    nrow(X_norm), ncol(X_norm)))

        genes_in_matrix <- PATCH_GENES[PATCH_GENES %in% norm_gen]
        genes_missing   <- setdiff(PATCH_GENES, norm_gen)
        if (length(genes_missing) > 0)
            cat(sprintf("    WARNING: %d gene(s) not in norm.mtx: %s\n",
                        length(genes_missing),
                        paste(genes_missing, collapse = ", ")))
        cat(sprintf("    Genes available in norm.mtx: %d / %d\n",
                    length(genes_in_matrix), length(PATCH_GENES)))

        if (length(genes_in_matrix) == 0) {
            rm(X_norm, nm); gc(verbose = FALSE)
            cat("    Nothing to patch.\n\n")
            return(invisible(NULL))
        }

        # One dense extraction for the whole panel rather than one per gene.
        panel_mat <- as.matrix(X_norm[, genes_in_matrix, drop = FALSE])
        rownames(panel_mat) <- norm_bc

        n_patched <- 0L
        n_skipped <- 0L

        for (reg in ds$regions) {
            gam_path <- gam_file_path(ds$dir, reg)
            if (!file.exists(gam_path)) {
                cat(sprintf("    %s: file not found -- skipped\n", reg))
                next
            }

            gam_df       <- read.csv(gam_path, stringsAsFactors = FALSE,
                                     check.names = FALSE)
            genes_to_add <- setdiff(genes_in_matrix, colnames(gam_df))

            if (length(genes_to_add) == 0) {
                cat(sprintf("    %s: all %d genes already present -- no patch\n",
                            reg, length(genes_in_matrix)))
                n_skipped <- n_skipped + 1L
                next
            }

            bc_idx    <- match(gam_df$barcode, norm_bc)
            n_matched <- sum(!is.na(bc_idx))
            if (n_matched < nrow(gam_df)) {
                cat(sprintf("    %s: WARNING -- only %d/%d barcodes matched norm.mtx\n",
                            reg, n_matched, nrow(gam_df)))
            }

            valid <- !is.na(bc_idx)
            for (g in genes_to_add) {
                expr_vals        <- rep(NA_real_, nrow(gam_df))
                expr_vals[valid] <- panel_mat[bc_idx[valid], g]
                gam_df[[g]]      <- expr_vals
            }

            write.csv(gam_df, gam_path, row.names = FALSE)
            n_patched <- n_patched + 1L

            cat(sprintf("    %s: +%d column(s) (%s) | %d cells | %d barcodes matched\n",
                        reg, length(genes_to_add),
                        paste(genes_to_add, collapse = ", "),
                        nrow(gam_df), n_matched))
        }

        rm(X_norm, nm, panel_mat); gc(verbose = FALSE)
        cat(sprintf("  %s summary: %d region(s) patched, %d already complete\n\n",
                    ds$tag, n_patched, n_skipped))
    }

    for (ds in DATASETS) patch_gam_csvs(ds)

    cat("  Gene column patch complete.\n\n")


    # =============================================================================
    # PART 9.3 -- HEATMAPS
    # =============================================================================
    # Row order: block 1 (OXPHOS) GOI, composite, individual DFGs;
    # block 2 (HKG) cHKG, individual HKGs.
    # =============================================================================

    cat("\n", strrep("#", 76), "\n", sep = "")
    cat("  PART 9.3 -- HEATMAPS\n")
    cat(strrep("#", 76), "\n\n", sep = "")

    HM_GOI <- GOI
    HM_DFG <- DFG
    HM_HKG <- HKG_USE
    HM_ALL <- unique(c(HM_GOI, HM_DFG, HM_HKG))

    CDFG_LBL <- DFG_COMPOSITE_LABEL
    CHKG_LBL <- "cHKG"

    cat("  GOI       :", HM_GOI, "\n")
    cat("  DFG       :", paste(HM_DFG, collapse = ", "), "\n")
    cat("  HKG       :", paste(HM_HKG, collapse = ", "), "\n")
    cat("  Composite :", CDFG_LBL, "\n\n")


    # == 9.3.1  Load the normalised expression for the panel ==================
    load_expression <- function(dataset_dir, genes_needed) {

        cat("  Loading normalised matrix from:\n    ", dataset_dir, "\n")
        nm <- load_norm_matrix(dataset_dir, as_csparse = TRUE)
        X            <- nm$X
        gene_symbols <- nm$genes

        cat("    Matrix:", nrow(X), "cells x", ncol(X), "genes\n")

        present   <- genes_needed[genes_needed %in% gene_symbols]
        not_found <- setdiff(genes_needed, gene_symbols)

        if (length(not_found) > 0)
            cat("    WARNING: not in the matrix, excluded:",
                paste(not_found, collapse = ", "), "\n")

        # The GOI is indexed by name below, so its absence stops the run here.
        if (!(HM_GOI %in% present))
            stop("GOI '", HM_GOI, "' not found in ", dataset_dir,
                 "/deg_support/norm_features.tsv")

        expr_mat <- as.matrix(X[, present, drop = FALSE])

        cat("    Genes resolved:", length(present), "/", length(genes_needed), "\n")
        cat("    Columns:", paste(colnames(expr_mat), collapse = ", "), "\n\n")

        rm(X, nm); gc(verbose = FALSE)
        expr_mat   # cells x genes, CP10K log1p, no further transform
    }


    # == 9.3.2  Build the per-dataset heatmap input ===========================
    build_heatmap_input <- function(dataset_dir, regions, expr_mat) {

        meta_path <- file.path(dataset_dir, "metadata", "obs_selected.csv")
        if (!file.exists(meta_path))
            stop("obs_selected.csv not found: ", meta_path)
        meta <- read.csv(meta_path, stringsAsFactors = FALSE, check.names = FALSE)
        rownames(meta) <- meta$barcode
        cat("  Metadata loaded:", nrow(meta), "selected cells\n")

        cat("  Region column   :", REGION_COL,
            "| present:", REGION_COL %in% colnames(meta), "\n")
        cat("  Supertype column:", STYPE_COL,
            "| present:", STYPE_COL %in% colnames(meta), "\n")

        keep      <- (meta[[REGION_COL]] %in% regions) &
                     (meta[[STYPE_COL]] == PRIMARY_STYPE)
        meta_filt <- meta[keep, , drop = FALSE]
        cat("  Cells passing the region and supertype filter:",
            nrow(meta_filt), "\n")

        if (nrow(meta_filt) == 0) {
            cat("  Regions in metadata   :",
                paste(sort(unique(meta[[REGION_COL]])), collapse = ", "), "\n")
            cat("  Supertypes in metadata:",
                paste(sort(unique(meta[[STYPE_COL]])), collapse = ", "), "\n")
            stop("Filter returned 0 cells. Check REGION_COL, STYPE_COL and PRIMARY_STYPE.")
        }

        common_bc <- intersect(rownames(meta_filt), rownames(expr_mat))
        if (length(common_bc) == 0)
            stop("No common barcodes between the metadata and the expression matrix")

        meta_filt <- meta_filt[common_bc, , drop = FALSE]
        expr_filt <- expr_mat[common_bc, , drop = FALSE]

        cat("  Common barcodes:", nrow(meta_filt), "\n")
        cat("  Region counts:\n")
        print(table(meta_filt[[REGION_COL]]))

        # Drop regions with no selected cells from the factor levels.
        regions_present <- regions[regions %in% unique(as.character(
            meta_filt[[REGION_COL]]))]
        if (length(regions_present) < length(regions))
            cat("  Regions with no cells, omitted from the panel split:",
                paste(setdiff(regions, regions_present), collapse = ", "), "\n")

        region_fac <- factor(meta_filt[[REGION_COL]], levels = regions_present)
        goi_expr   <- expr_filt[, HM_GOI]
        cell_order <- order(region_fac, goi_expr)

        expr_s <- expr_filt[cell_order, , drop = FALSE]
        goi_s  <- goi_expr[cell_order]
        reg_s  <- region_fac[cell_order]

        dfg_pres <- HM_DFG[HM_DFG %in% colnames(expr_s)]
        hkg_pres <- HM_HKG[HM_HKG %in% colnames(expr_s)]

        cat("  DFG present:", paste(dfg_pres, collapse = ", "), "\n")
        cat("  HKG present:", paste(hkg_pres, collapse = ", "), "\n")

        if (length(dfg_pres) == 0)
            stop("No DFG genes found in the expression matrix")

        make_composite <- function(genes) {
            if (length(genes) == 0) return(NULL)
            if (length(genes) == 1) return(expr_s[, genes])
            rowMeans(expr_s[, genes, drop = FALSE])
        }

        cdfg_vec <- make_composite(dfg_pres)
        chkg_vec <- make_composite(hkg_pres)

        row_list <- list(
            matrix(expr_s[, HM_GOI], nrow = 1, dimnames = list(HM_GOI,   NULL)),
            matrix(cdfg_vec,         nrow = 1, dimnames = list(CDFG_LBL, NULL))
        )
        for (g in dfg_pres)
            row_list[[length(row_list) + 1L]] <-
                matrix(expr_s[, g], nrow = 1, dimnames = list(g, NULL))

        if (!is.null(chkg_vec)) {
            row_list[[length(row_list) + 1L]] <-
                matrix(chkg_vec, nrow = 1, dimnames = list(CHKG_LBL, NULL))
            for (g in hkg_pres)
                row_list[[length(row_list) + 1L]] <-
                    matrix(expr_s[, g], nrow = 1, dimnames = list(g, NULL))
        }

        expr_final <- do.call(rbind, row_list)   # genes x cells
        cat("  Matrix assembled:", nrow(expr_final), "rows x",
            ncol(expr_final), "cols\n")
        cat("  Row order:", paste(rownames(expr_final), collapse = " | "), "\n\n")

        zscore_cap <- function(x)
            pmin(pmax((x - mean(x, na.rm = TRUE)) /
                      (sd(x, na.rm = TRUE) + 1e-9), -1), 2.0)
        mat_scaled <- t(apply(expr_final, 1, zscore_cap))
        dimnames(mat_scaled) <- dimnames(expr_final)

        n_per_region <- table(reg_s)[regions_present]
        n_per_region[is.na(n_per_region)] <- 0L

        list(
            mat_scaled      = mat_scaled,
            region_s        = reg_s,
            goi_s           = goi_s,
            regions_present = regions_present,
            n_per_region    = n_per_region,
            dfg_pres        = dfg_pres,
            hkg_pres        = hkg_pres,
            n_cells         = ncol(mat_scaled)
        )
    }


    # == 9.3.3  Draw and save =================================================
    draw_heatmap <- function(hm_input, dataset_label, output_pdf, output_png) {

        mat       <- hm_input$mat_scaled
        reg_s     <- hm_input$region_s
        regions   <- hm_input$regions_present
        n_per_reg <- hm_input$n_per_region
        n_cells   <- hm_input$n_cells
        dfg_pres  <- hm_input$dfg_pres
        hkg_pres  <- hm_input$hkg_pres
        rn        <- rownames(mat)

        cat("  Drawing heatmap:", dataset_label,
            "| rows:", nrow(mat), "| cols:", ncol(mat), "\n")

        reg_labels <- sapply(regions, function(r) {
            n   <- as.integer(n_per_reg[r])
            pct <- round(100 * n / n_cells, 1)
            paste0(r, "\nn=", n, "\n(", pct, "%)")
        })

        n_reg_levels    <- length(levels(reg_s))
        reg_block_fills <- setNames(rep("#000000", n_reg_levels), levels(reg_s))

        col_ann <- HeatmapAnnotation(
            Region = anno_block(
                gp        = gpar(fill = reg_block_fills[levels(reg_s)]),
                labels    = reg_labels[levels(reg_s)],
                labels_gp = gpar(col = "white", fontsize = 7, fontface = "bold")
            ),
            show_legend          = FALSE,
            show_annotation_name = FALSE
        )

        row_cat <- setNames(character(length(rn)), rn)
        row_cat[HM_GOI]            <- "GOI"
        row_cat[CDFG_LBL]          <- "cDFG"
        row_cat[rn %in% dfg_pres]  <- "DFG"
        if (CHKG_LBL %in% rn) row_cat[CHKG_LBL] <- "cHKG"
        row_cat[rn %in% hkg_pres]  <- "HKG"

        cat_cols <- c(GOI  = "#CC2200", cDFG = "#220044",
                      DFG  = "#4B0082", cHKG = "#003355",
                      HKG  = "#006699")

        row_ann <- rowAnnotation(
            Category = row_cat,
            col      = list(Category = cat_cols),
            annotation_legend_param = list(
                Category = list(title     = "Gene category",
                                labels_gp = gpar(fontsize = 8),
                                title_gp  = gpar(fontsize = 9, fontface = "bold"))
            ),
            show_annotation_name = FALSE,
            width = unit(3, "mm")
        )

        row_name_gp <- gpar(
            fontsize = 8,
            fontface = ifelse(rn %in% c(CDFG_LBL, CHKG_LBL), "bold", "plain"),
            col      = ifelse(rn == HM_GOI,          "#CC2200",
                       ifelse(rn == CDFG_LBL,        "#220044",
                       ifelse(rn %in% dfg_pres,      "#4B0082",
                       ifelse(rn == CHKG_LBL,        "#003355",
                       ifelse(rn %in% hkg_pres,      "#006699",
                              "black")))))
        )

        col_fun <- colorRamp2(c(-1, 0, 1.0, 2.0), viridis(4, option = "mako"))

        block1 <- c(HM_GOI, CDFG_LBL, dfg_pres)
        block2 <- c(if (CHKG_LBL %in% rn) CHKG_LBL, hkg_pres)

        # The row split is built only over blocks that have rows.
        row_block   <- ifelse(rn %in% block2, "2: HKG", "1: OXPHOS")
        block_order <- c("1: OXPHOS", "2: HKG")
        row_split   <- factor(row_block,
                              levels = block_order[block_order %in% row_block])

        col_split <- factor(reg_s, levels = regions)

        pdf_w <- min(50, max(18, round(n_cells * 0.008 + 6)))
        pdf_h <- 7

        ht <- Heatmap(
            mat,
            name                 = "z-score",
            col                  = col_fun,
            cluster_rows         = FALSE,
            show_row_names       = TRUE,
            row_names_side       = "right",
            row_names_gp         = row_name_gp,
            row_split            = row_split,
            row_title_gp         = gpar(fontsize = 9, fontface = "bold"),
            row_gap              = unit(1, "mm"),
            cluster_columns      = FALSE,
            show_column_names    = FALSE,
            column_split         = col_split,
            column_title         = NULL,
            column_gap           = unit(0.8, "mm"),
            top_annotation       = col_ann,
            right_annotation     = row_ann,
            border               = TRUE,
            use_raster           = TRUE,
            raster_quality       = 3,
            heatmap_legend_param = list(
                title         = "z-score",
                title_gp      = gpar(fontsize = 9, fontface = "bold"),
                labels_gp     = gpar(fontsize = 8),
                legend_height = unit(3.5, "cm"),
                at            = c(-1, 0, 1, 2),
                labels        = c("-1", "0", "1", "2 (cap)")
            )
        )

        title_str <- paste0(
            HM_GOI, " vs ", CDFG_LBL, " and HKG -- ", dataset_label,
            "  [", ARM_TAG, " | ", PRIMARY_STYPE,
            " | z-score capped [-1, 2.0]]"
        )
        draw_args <- list(
            ht,
            padding                = unit(c(10, 4, 4, 4), "mm"),
            annotation_legend_side = "right",
            heatmap_legend_side    = "right",
            merge_legends          = TRUE
        )

        # Qualified so the generic is not taken from the search path.
        emit <- function() {
            do.call(ComplexHeatmap::draw, draw_args)
            grid.text(title_str,
                      x    = unit(0.01, "npc"), y = unit(0.99, "npc"),
                      just = c("left", "top"),
                      gp   = gpar(fontsize = 10, fontface = "bold",
                                  col = "#1A1A1A"))
        }

        cat("  Saving PDF:", output_pdf, "\n")
        pdf(output_pdf, width = pdf_w, height = pdf_h)
        emit()
        dev.off()
        cat("  PDF saved.\n")

        # PNG keeps the PDF geometry; resolution drops above the pixel-width cap.
        png_max_px <- 4000
        png_res    <- min(150, floor(png_max_px / pdf_w))
        cat("  Saving PNG:", output_png,
            sprintf("(%.0f x %.0f in at %d dpi)\n", pdf_w, pdf_h, png_res))
        png(output_png, width = pdf_w, height = pdf_h, units = "in", res = png_res)
        emit()
        dev.off()
        cat("  PNG saved.\n\n")

        invisible(ht)
    }


    # == 9.3.4  Run both datasets =============================================
    hm_summary <- list()

    for (ds in DATASETS) {
        cat(strrep("=", 70), "\n", sep = "")
        cat("  ", toupper(ds$file_label), "\n", sep = "")
        cat(strrep("=", 70), "\n\n", sep = "")

        ds_expr  <- load_expression(ds$dir, HM_ALL)
        ds_input <- build_heatmap_input(ds$dir, ds$regions, ds_expr)

        draw_heatmap(
            hm_input      = ds_input,
            dataset_label = ds$plot_label,
            output_pdf    = file.path(PART9_DIR,
                                paste0(OUT_PREFIX, "_Heatmap_", ds$file_label,
                                       "_", HM_GOI, ".pdf")),
            output_png    = file.path(PART9_DIR,
                                paste0(OUT_PREFIX, "_Heatmap_", ds$file_label,
                                       "_", HM_GOI, ".png"))
        )

        hm_summary[[ds$tag]] <- ds_input$n_cells
        rm(ds_expr, ds_input); gc(verbose = FALSE)
    }

    cat(strrep("=", 70), "\n", sep = "")
    cat("  Part 9.3 complete. Outputs in:", PART9_DIR, "\n\n")
    cat("  Row order, both heatmaps:\n")
    cat("    1: OXPHOS --", HM_GOI, "(GOI) |", CDFG_LBL, "|",
        paste(HM_DFG, collapse = " "), "\n")
    cat("    2: HKG    --", CHKG_LBL, "|", paste(HM_HKG, collapse = " "), "\n")
    cat("  Cells:",
        paste(names(hm_summary), unlist(hm_summary), sep = " = ",
              collapse = " | "), "\n")
    cat(strrep("=", 70), "\n", sep = "")


    # =============================================================================
    # PART 9.4 -- GAM-REML-PRSS
    # =============================================================================
    # =============================================================================

    cat("\n\n", strrep("#", 76), "\n", sep = "")
    cat("  PART 9.4 -- GAM-REML-PRSS\n")
    cat(strrep("#", 76), "\n\n", sep = "")

    # Re-seeded so the GAM stream does not depend on the QC and heatmap sections.
    set.seed(123)

    # == Estimator ============================================================
    # GAM_DONOR_RE: TRUE fits a donor random intercept; FALSE gives cell-level fits.
    GAM_DONOR_RE        <- TRUE
    DONOR_RE_MIN_LEVELS <- 3L

    # Inverse-donor-count weights; FALSE gives every cell w = 1.
    W_DONOR_WEIGHTS <- TRUE

    # Leave-one-gene-out contribution analysis (Gene_Contribution_LOO sheet); slow.
    RUN_LOO <- FALSE

    # Cap on phase-2 PRSS iterations.
    PRSS_PHASE2_MAX <- 5L
    # DONOR_RE_LABEL is defined at the top level (Section 0d) for the workers.

    cat("  Estimator: ",
        if (GAM_DONOR_RE)
            sprintf("donor random intercept, s(donor, bs = \"re\"), minimum %d donors",
                    DONOR_RE_MIN_LEVELS)
        else "cell level (donor random intercept disabled)", "\n", sep = "")
    cat("  Weights  : ",
        if (W_DONOR_WEIGHTS)
            "inverse donor cell count, w = 1/n_i (each donor sums to 1)"
        else "uniform, w = 1 (inverse-donor weighting disabled)",
        "\n", sep = "")
    cat("  PRSS     : phase 1 = 8 iterations, phase 2 = up to ",
        PRSS_PHASE2_MAX, " (early stop ", PRSS_EARLY_STOP, ")\n", sep = "")
    cat("  LOO      : ", if (isTRUE(RUN_LOO)) "enabled" else
        "skipped (RUN_LOO = FALSE)", "\n", sep = "")

    # == Parallelism ==========================================================
    # PAR_CORES PSOCK workers; 1L runs everything in this process.
    .p9_detected <- suppressWarnings(parallel::detectCores())
    PAR_CORES <- if (is.na(.p9_detected)) 1L
                 else max(1L, as.integer(.p9_detected) - 1L)

    .p9_cl <- NULL
    if (PAR_CORES > 1L) {
        .p9_cl <- parallel::makeCluster(PAR_CORES)   # PSOCK
        # Released whether the arm finishes or fails, alongside the log sink.
        on.exit({
            if (!is.null(.p9_cl)) try(parallel::stopCluster(.p9_cl), silent = TRUE)
        }, add = TRUE)
        parallel::clusterEvalQ(.p9_cl, { library(mgcv); NULL })
        parallel::clusterExport(.p9_cl,
            c("DONOR_RE_LABEL", "p9_gam_formula", "p9_goi_smooth_idx",
              "p9_smooth_lambda", "p9_smooth_stat", "p9_model_has_re",
              "p9_calculate_prss", "p9_extract_reml_iters",
              "p9_prss_worker", "p9_loo_worker"))
        cat("  Parallel : ", PAR_CORES, " PSOCK worker(s) for phase-2 PRSS ",
            "iterations and LOO refits\n", sep = "")
    } else {
        cat("  Parallel : disabled (PAR_CORES = 1); sequential fitting\n")
    }

    # Pre-analysis checks, recorded in the ledger.
    log_result("GR-00", "Pre-analysis input checks performed in Part 9.1",
               if (QC_FAIL_COUNT == 0) "PASS" else "FAIL",
               detail = sprintf(
                   "QC-01 to QC-13 ran against the same files: %d PASS / %d WARN / %d FAIL",
                   .chk_pass, .chk_warn, QC_FAIL_COUNT))


    # =============================================================================
    # SECTION A: SHARED HELPER FUNCTIONS
    # =============================================================================

    # == A0. Donor random intercept: data, formula, term selection ============
    # Donor enters every fitting frame as a factor with unused levels dropped.
    as_donor_factor <- function(x) droplevels(factor(as.character(x)))

    n_donor_levels <- function(data) {
        if (is.null(data) || !("donor" %in% names(data))) return(0L)
        nlevels(as_donor_factor(data$donor))
    }

    # Frames with too few donors fall back to the cell-level formula.
    donor_re_ok <- function(data) {
        isTRUE(GAM_DONOR_RE) && n_donor_levels(data) >= DONOR_RE_MIN_LEVELS
    }

    # w = 1/(cells this donor contributes to this frame), so each donor's weights
    # sum to 1.
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
    make_gam_formula <- function(k, use_re) p9_gam_formula(k, use_re, GAM_BS)
    goi_smooth_idx   <- p9_goi_smooth_idx
    smooth_lambda    <- p9_smooth_lambda
    smooth_stat      <- p9_smooth_stat
    model_has_re     <- p9_model_has_re

    # Terms to drop for a population-level prediction.
    re_exclude <- function(model) {
        if (model_has_re(model)) DONOR_RE_LABEL else character(0)
    }

    # == A1. Per-region data from the gam_data CSVs ===========================
    build_region_data <- function(dataset_dir, regions, dfg_genes, hkg_genes,
                                  supertype = PRIMARY_STYPE,
                                  min_cells = 10) {

        hkg_expr_df  <- NULL
        hkg_pres_all <- character(0)

        norm_paths <- file.path(dataset_dir, "deg_support",
                                c("norm.mtx", "norm_barcodes.tsv",
                                  "norm_features.tsv"))

        if (all(file.exists(norm_paths))) {
            cat("  Loading HKG expression from norm.mtx ...\n")
            nm       <- load_norm_matrix(dataset_dir, as_csparse = TRUE)
            gene_sym <- nm$genes

            hkg_pres_all <- hkg_genes[hkg_genes %in% gene_sym]
            if (length(hkg_genes) == 0) {
                cat("  No HKG comparator for this arm; panel only.\n")
            } else if (length(hkg_pres_all) == 0) {
                cat("  WARNING: no HKG gene found in norm.mtx\n")
            } else {
                cat("  HKG genes found in norm.mtx:",
                    paste(hkg_pres_all, collapse = ", "), "\n")
                hkg_sub     <- as.matrix(nm$X[, hkg_pres_all, drop = FALSE])
                hkg_expr_df <- as.data.frame(hkg_sub)
                hkg_expr_df$barcode <- nm$barcodes
            }
            rm(nm); gc(verbose = FALSE)
        } else {
            cat("  WARNING: norm.mtx not found -- HKG fits will be skipped\n")
        }

        region_data <- lapply(regions, function(reg) {
            fpath <- gam_file_path(dataset_dir, reg, supertype)
            if (!file.exists(fpath)) {
                cat("  Region", reg, ": file not found --", basename(fpath), "\n")
                return(NULL)
            }
            df    <- read.csv(fpath, stringsAsFactors = FALSE, check.names = FALSE)
            n_raw <- nrow(df)

            # Normally a no-op: Part 8 already selected on GOI > 0.
            keep_goi  <- df[["GOI"]] > 0
            n_goi_neg <- sum(!keep_goi)
            df        <- df[keep_goi, , drop = FALSE]
            n_goi_pos <- nrow(df)

            cat("  Region", reg, ": n_raw =", n_raw,
                "| GOI<=0 excluded:", n_goi_neg,
                "| GOI>0 retained:", n_goi_pos, "\n")

            if (n_goi_pos < min_cells) {
                cat("  Region", reg, ": only", n_goi_pos,
                    "GOI-positive cells -- skipping (min", min_cells, ")\n")
                return(NULL)
            }

            n <- n_goi_pos

            # Donor is carried into every fitting frame.
            if (!("donor" %in% colnames(df)))
                stop(sprintf(
                    "Region %s: gam_data has no donor column, so the donor random intercept cannot be fitted.",
                    reg))
            donor_vec <- as_donor_factor(df[["donor"]])
            if (anyNA(donor_vec))
                stop(sprintf("Region %s: %d cell(s) have a missing donor label.",
                             reg, sum(is.na(donor_vec))))
            # Weights from this region's own cell counts.
            weight_vec <- donor_weights(donor_vec)

            dfg_pres <- dfg_genes[dfg_genes %in% colnames(df)]
            goi_vec  <- df[["GOI"]]
            dfg_mat  <- as.matrix(df[, dfg_pres, drop = FALSE])

            # The composite is recomputed from the individual DFG columns.
            if (length(dfg_pres) < length(dfg_genes))
                stop(sprintf("Region %s: gam_data is missing DFG column(s): %s",
                             reg, paste(setdiff(dfg_genes, dfg_pres),
                                        collapse = ", ")))
            cdfg <- rowMeans(dfg_mat)

            chkg     <- NULL
            hkg_pres <- character(0)
            hkg_mat  <- NULL

            if (!is.null(hkg_expr_df) && length(hkg_pres_all) > 0) {
                bc_gam    <- df[["barcode"]]
                idx       <- match(bc_gam, hkg_expr_df$barcode)
                n_matched <- sum(!is.na(idx))

                # Every GAM barcode must be in norm.mtx; a shortfall stops the run.
                if (n_matched < n) {
                    stop(sprintf(
                        paste0("Region %s: %d of %d GAM barcodes are absent from ",
                               "norm_barcodes.tsv. The gam_data CSVs and ",
                               "deg_support/norm.mtx are not from the same Part 8 ",
                               "run for this arm."),
                        reg, n - n_matched, n))
                }

                hkg_mat <- as.matrix(
                    hkg_expr_df[idx, hkg_pres_all, drop = FALSE])
                rownames(hkg_mat) <- bc_gam
                hkg_pres <- hkg_pres_all
                chkg     <- rowMeans(hkg_mat)
                cat("  Region", reg, ": HKG matched", n_matched, "/", n,
                    "barcodes\n")
            }

            cat("  Region", reg,
                "| n_raw:", n_raw,
                "| n_goi_neg:", n_goi_neg,
                "| n_goi_pos:", n_goi_pos,
                "| n_trimmed: 0",
                "| n_used:", n,
                "| donors:", nlevels(donor_vec),
                "| sum(w):", round(sum(weight_vec), 4),
                "| unique GOI:", length(unique(goi_vec)),
                "| DFG:", length(dfg_pres),
                "| HKG:", length(hkg_pres), "\n")

            list(
                dfg         = data.frame(GOI = goi_vec, Expression = cdfg,
                                         donor = donor_vec, w = weight_vec),
                hkg         = if (!is.null(chkg))
                                  data.frame(GOI = goi_vec, Expression = chkg,
                                             donor = donor_vec, w = weight_vec)
                              else NULL,
                goi_vec     = goi_vec,
                donor_vec   = donor_vec,
                weight_vec  = weight_vec,
                dfg_mat     = dfg_mat,
                hkg_mat     = hkg_mat,
                dfg_present = dfg_pres,
                hkg_present = hkg_pres,
                n_raw       = n_raw,
                n_goi_neg   = n_goi_neg,
                n_goi_pos   = n_goi_pos,
                n_trimmed   = 0L,
                n_used      = n,
                n_donors    = nlevels(donor_vec)
            )
        })
        names(region_data) <- regions
        region_data
    }


    # == A1b. Pooled region ===================================================
    build_pooled_region <- function(region_data, tag) {
        valid_rds <- Filter(Negate(is.null), region_data)
        if (length(valid_rds) == 0) {
            cat("  build_pooled_region: no valid regions to pool\n")
            return(NULL)
        }

        pooled_name <- paste0(tag, "_Pooled")
        cat("  Building pooled region '", pooled_name, "' from ",
            length(valid_rds), " region(s): ",
            paste(names(valid_rds), collapse = ", "), "\n", sep = "")

        all_dfg_present <- unique(unlist(lapply(valid_rds, `[[`, "dfg_present")))

        dfg_list <- lapply(valid_rds, function(rd) {
            mat       <- rd$dfg_mat
            missing_g <- setdiff(all_dfg_present, colnames(mat))
            if (length(missing_g) > 0) {
                pad <- matrix(0, nrow = nrow(mat), ncol = length(missing_g),
                              dimnames = list(rownames(mat), missing_g))
                mat <- cbind(mat, pad)
            }
            mat[, all_dfg_present, drop = FALSE]
        })
        pooled_dfg_mat <- do.call(rbind, dfg_list)
        pooled_goi     <- unlist(lapply(valid_rds, `[[`, "goi_vec"))
        pooled_cdfg    <- rowMeans(pooled_dfg_mat)

        # Donor joined as character; a donor in several regions becomes one level.
        pooled_donor <- as_donor_factor(
            unlist(lapply(valid_rds,
                          function(rd) as.character(rd$donor_vec)),
                   use.names = FALSE))
        if (length(pooled_donor) != length(pooled_goi))
            stop("build_pooled_region: donor and GOI vectors have different ",
                 "lengths (", length(pooled_donor), " vs ", length(pooled_goi),
                 ").")

        # Weights recomputed on the pooled frame.
        pooled_w <- donor_weights(pooled_donor)

        has_hkg         <- sapply(valid_rds, function(rd) !is.null(rd$hkg_mat))
        pooled_hkg_mat  <- NULL
        pooled_hkg_pres <- character(0)
        pooled_chkg     <- NULL
        pooled_hkg_df   <- NULL

        # The pooled HKG fit runs only if every contributing region has HKG data.
        if (all(has_hkg)) {
            all_hkg_present <- unique(unlist(lapply(valid_rds, `[[`, "hkg_present")))
            hkg_list <- lapply(valid_rds, function(rd) {
                mat       <- rd$hkg_mat
                missing_g <- setdiff(all_hkg_present, colnames(mat))
                if (length(missing_g) > 0) {
                    pad <- matrix(0, nrow = nrow(mat), ncol = length(missing_g),
                                  dimnames = list(rownames(mat), missing_g))
                    mat <- cbind(mat, pad)
                }
                mat[, all_hkg_present, drop = FALSE]
            })
            pooled_hkg_mat  <- do.call(rbind, hkg_list)
            pooled_hkg_pres <- all_hkg_present
            pooled_chkg     <- rowMeans(pooled_hkg_mat)
            pooled_hkg_df   <- data.frame(GOI = pooled_goi,
                                          Expression = pooled_chkg,
                                          donor = pooled_donor,
                                          w = pooled_w)
        } else if (any(has_hkg)) {
            cat("  Pooled HKG fit skipped:", sum(!has_hkg), "of",
                length(has_hkg), "region(s) have no HKG data.\n")
        }

        n_pooled      <- length(pooled_goi)
        n_raw_sum     <- sum(sapply(valid_rds, function(rd) rd$n_raw))
        n_goi_neg_sum <- sum(sapply(valid_rds, function(rd) rd$n_goi_neg))
        n_goi_pos_sum <- sum(sapply(valid_rds, function(rd) rd$n_goi_pos))

        cat("  Pooled cells:", n_pooled,
            "| n_raw (sum):", n_raw_sum,
            "| n_goi_neg (sum):", n_goi_neg_sum,
            "| n_goi_pos (sum):", n_goi_pos_sum,
            "| n_trimmed: 0",
            "| donors:", nlevels(pooled_donor),
            "| sum(w):", round(sum(pooled_w), 4),
            "| unique GOI:", length(unique(pooled_goi)),
            "| DFG genes:", length(all_dfg_present),
            "| HKG genes:", length(pooled_hkg_pres), "\n")

        list(
            dfg         = data.frame(GOI = pooled_goi, Expression = pooled_cdfg,
                                     donor = pooled_donor, w = pooled_w),
            hkg         = pooled_hkg_df,
            goi_vec     = pooled_goi,
            donor_vec   = pooled_donor,
            weight_vec  = pooled_w,
            dfg_mat     = pooled_dfg_mat,
            hkg_mat     = pooled_hkg_mat,
            dfg_present = all_dfg_present,
            hkg_present = pooled_hkg_pres,
            n_raw       = n_raw_sum,
            n_goi_neg   = n_goi_neg_sum,
            n_goi_pos   = n_goi_pos_sum,
            n_trimmed   = 0L,
            n_used      = n_pooled,
            n_donors    = nlevels(pooled_donor)
        )
    }


    # == A1c. Cell count audit sheet ==========================================
    build_cell_count_sheet <- function(tag, region_data, regions,
                                       pooled_region_name) {
        rows <- list()

        for (reg in regions) {
            rd <- region_data[[reg]]
            if (is.null(rd)) {
                rows[[length(rows) + 1L]] <- data.frame(
                    Tag            = tag,
                    Region         = reg,
                    Is_Pooled      = FALSE,
                    N_Raw          = NA_integer_,
                    N_GOI_Neg      = NA_integer_,
                    N_GOI_Pos      = NA_integer_,
                    N_Trimmed      = NA_integer_,
                    N_Used_for_GAM = NA_integer_,
                    N_Donors       = NA_integer_,
                    Tally_OK       = NA,
                    Note = "Region skipped (insufficient cells or file missing)",
                    stringsAsFactors = FALSE)
                next
            }
            ok <- isTRUE(rd$n_used == rd$n_goi_pos)
            rows[[length(rows) + 1L]] <- data.frame(
                Tag            = tag,
                Region         = reg,
                Is_Pooled      = FALSE,
                N_Raw          = rd$n_raw,
                N_GOI_Neg      = rd$n_goi_neg,
                N_GOI_Pos      = rd$n_goi_pos,
                N_Trimmed      = 0L,
                N_Used_for_GAM = rd$n_used,
                N_Donors       = rd$n_donors,
                Tally_OK       = ok,
                Note           = if (ok) "OK" else
                    paste0("MISMATCH: n_goi_pos=", rd$n_goi_pos,
                           " but n_used=", rd$n_used),
                stringsAsFactors = FALSE)
        }

        prd <- region_data[[pooled_region_name]]
        if (!is.null(prd)) {
            ok_p <- isTRUE(prd$n_used == prd$n_goi_pos)
            rows[[length(rows) + 1L]] <- data.frame(
                Tag            = tag,
                Region         = pooled_region_name,
                Is_Pooled      = TRUE,
                N_Raw          = prd$n_raw,
                N_GOI_Neg      = prd$n_goi_neg,
                N_GOI_Pos      = prd$n_goi_pos,
                N_Trimmed      = 0L,
                N_Used_for_GAM = prd$n_used,
                N_Donors       = prd$n_donors,
                Tally_OK       = ok_p,
                Note           = if (ok_p) "OK (sum of constituent regions)" else
                    paste0("MISMATCH: n_goi_pos=", prd$n_goi_pos,
                           " but n_used=", prd$n_used),
                stringsAsFactors = FALSE)
        }

        do.call(rbind, rows)
    }


    # == A2. PRSS =============================================================
    # PRSS = RSS of the fitted model plus the s(GOI) roughness penalty;
    # body in p9_calculate_prss() (Section 0d).
    calculate_prss <- function(model, data = NULL) p9_calculate_prss(model)


    # == A3. Equation string ==================================================
    # Donor BLUPs are left out of the GOI-response equation.
    extract_equation <- function(model) {
        cf <- coef(model)
        nm <- names(cf)
        is_goi_smooth <- grepl("^s\\(GOI\\)", nm)
        is_any_smooth <- grepl("^s\\(", nm)
        pt <- cf[!is_any_smooth]
        st <- cf[ is_goi_smooth]
        lhs <- sprintf("%.9f", pt[1])
        if (length(pt) > 1) lhs <- paste0(lhs, " + ", sprintf("%.9f*x", pt[2]))
        rhs <- paste(mapply(function(v, i) sprintf("%.9f*\u03c6%d(x)", v, i),
                            st, seq_along(st)), collapse = " + ")
        paste0("f(x) = ", lhs, if (length(st) > 0) " + " else "", rhs)
    }


    # == A4. REML iteration extraction ========================================
    # Body in Section 0d, shared with the workers.
    extract_reml_iterations_enhanced <- p9_extract_reml_iters


    # == A5. Convergence details ==============================================
    extract_convergence_details <- function(model) {
        blank <- list(convergence_criterion = "Unknown",
                      gradient_value = NA_character_,
                      gradient_norm = NA_real_, relative_score_change = NA_real_,
                      iterations = NA_integer_, max_iterations_reached = NA,
                      final_reml_score = NA_real_, optimizer = NA_character_,
                      monotonic_decrease = "NA")
        if (is.null(model) || is.null(model$outer.info)) return(blank)
        tryCatch({
            iters <- model$outer.info$iter
            gv <- NA_character_; gn <- NA_real_
            if (!is.null(model$outer.info$grad)) {
                g <- model$outer.info$grad
                if (is.vector(g) && length(g) > 0) {
                    fg <- g[length(g)]; gv <- format(fg, digits = 8); gn <- abs(fg)
                } else if (is.matrix(g) && nrow(g) > 0) {
                    fg <- g[nrow(g), ]
                    gv <- paste(format(fg, digits = 8), collapse = ", ")
                    gn <- sqrt(sum(fg^2))
                }
            }
            rsc <- NA_real_
            if (!is.null(model$outer.info$score) &&
                length(model$outer.info$score) > 1) {
                sc  <- model$outer.info$score
                rsc <- tail(abs(diff(sc)) / abs(sc[-length(sc)]), 1)
            }
            mono <- "NO"
            if (!is.null(model$outer.info$score) &&
                length(model$outer.info$score) > 1)
                mono <- if (all(diff(model$outer.info$score) <= 1e-10)) "YES"
                        else "NO"
            crit <- "Unknown"
            if (!is.na(iters) && !is.null(model$control$maxit) &&
                iters >= model$control$maxit) {
                crit <- "Maximum number of iterations reached"
            } else if (!is.na(gn) && gn < 1e-5) {
                crit <- "Gradient-based convergence"
            } else if (!is.na(rsc) && rsc < 1e-6) {
                crit <- "Score-based convergence"
            } else if (isTRUE(model$converged)) {
                crit <- "Fallback convergence"
            }
            list(convergence_criterion  = crit,
                 gradient_value         = gv,
                 gradient_norm          = gn,
                 relative_score_change  = rsc,
                 iterations             = iters,
                 max_iterations_reached =
                     if (!is.na(iters) && !is.null(model$control$maxit))
                         iters >= model$control$maxit else NA,
                 final_reml_score       = model$gcv.ubre,
                 optimizer              = paste(model$optimizer, collapse = ", "),
                 monotonic_decrease     = mono)
        }, error = function(e) {
            cat("  Warning - extract_convergence:", conditionMessage(e), "\n")
            blank
        })
    }


    # == A6. GAM fitting with two-phase PRSS optimisation =====================
    fit_gam_prss <- function(data, region_id, gene_set_label,
                             num_iterations = PRSS_N_ITER,
                             gamma = GAM_GAMMA,
                             cl = NULL) {
        set.seed(123)
        unique_x <- length(unique(data$GOI))
        max_k    <- min(10L, unique_x - 1L)

        # Donor is levelled once for all iterations, the null and LOO refits.
        if ("donor" %in% names(data)) data$donor <- as_donor_factor(data$donor)
        # The frame reaching gam() always carries w.
        data   <- ensure_weights(data)
        n_don  <- n_donor_levels(data)
        use_re <- donor_re_ok(data)

        cat("  [", region_id, "|", gene_set_label, "]",
            "cells:", nrow(data), "| unique GOI:", unique_x,
            "| max k:", max_k, "| donors:", n_don,
            "| donor RE:", if (use_re) "YES" else "NO",
            "| sum(w):", round(sum(data$w), 4), "\n")

        if (isTRUE(GAM_DONOR_RE) && !use_re)
            cat("  Donor random intercept not fitted: ", n_don,
                " donor level(s), minimum is ", DONOR_RE_MIN_LEVELS,
                ". This fit is cell level (D-08).\n", sep = "")

        if (unique_x < 4) {
            cat("  Skipping: fewer than 4 unique GOI values\n")
            return(NULL)
        }

        best_prss <- Inf; best_model <- NULL; best_k <- NULL
        best_iter <- 0;   prss_values <- numeric(num_iterations)
        model_params <- list(); reml_iters_list <- list()
        early_stop_counter <- 0; early_stop_triggered <- FALSE
        # Set when the winner was fitted on a worker, so the master refits it.
        best_par_iter <- NULL

        # The k window around the incumbent best.
        k_window <- function(bk)
            if (bk == max_k) c(max_k - 1L, max_k) else
            if (bk == 3L)    c(3L, 4L) else
            c(bk - 1L, bk, bk + 1L)

        run_one <- function(i, k) {
            set.seed(123 + i)
            gam_model <- gam(make_gam_formula(k, use_re),
                             data = data, weights = w,
                             method = GAM_METHOD,
                             select = GAM_SELECT, gamma = gamma)
            pc <- calculate_prss(gam_model, data)
            ri <- extract_reml_iterations_enhanced(gam_model)
            if (nrow(ri) > 0) {
                ri$prss_iteration <- i; ri$k_value <- k
                ri$prss_score <- pc$PRSS; ri$is_best_model <- FALSE
                reml_iters_list[[i]] <<- ri
            }
            prss_values[i] <<- pc$PRSS
            reml_scores_i <- if (!is.null(gam_model$outer.info$score))
                gam_model$outer.info$score else gam_model$gcv.ubre
            reml_n_i   <- if (!is.null(gam_model$outer.info$iter))
                gam_model$outer.info$iter else 1L
            reml_min_i <- min(reml_scores_i, na.rm = TRUE)
            reml_max_i <- max(reml_scores_i, na.rm = TRUE)
            sm_i <- summary(gam_model)
            model_params[[i]] <<- list(
                iteration     = i,
                k             = k,
                lambda        = pc$lambda,
                PRSS          = pc$PRSS,
                RSS           = pc$RSS,
                f_dpi         = pc$f_double_prime_integral,
                dev_explained = de_explained(gam_model, data),
                reml_n        = reml_n_i,
                reml_min      = reml_min_i,
                reml_max      = reml_max_i,
                # The s(GOI) row of s.table, not the last row.
                smooth_edf    = smooth_stat(sm_i, "s(GOI)", "edf"),
                p_value       = smooth_stat(sm_i, "s(GOI)", "p-value"))
            list(model = gam_model, prss = pc$PRSS, k = k)
        }

        # == Iteration bookkeeping shared by both phases ======================
        # Applied in index order, so the sequential and parallel paths match.
        iters_run <- integer(0)

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
            iters_run <<- c(iters_run, i)
            if (!isTRUE(r$ok)) {
                cat("    error iter", i, ":", r$message, "\n")
                prss_values[i]    <<- NA
                model_params[[i]] <<- NULL
                if (phase == 2L)
                    early_stop_counter <<- early_stop_counter + 1
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
                early_stop_counter <<- 0
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
                early_stop_counter <<- early_stop_counter + 1
                FALSE
            }
        }

        use_cluster <- !is.null(cl) && length(cl) >= 2L

        if (use_cluster) {
            # Exported once per fit; the null deviance is sent as a number.
            ctx <- list(data     = data,
                        bs       = GAM_BS,
                        method   = GAM_METHOD,
                        select   = GAM_SELECT,
                        gamma    = gamma,
                        use_re   = use_re,
                        null_dev = tryCatch(get_null_deviance(data),
                                            error = function(e) NA_real_))
            p9_fit_ctx <- ctx
            # Only ctx and p9_fit_ctx are exported.
            parallel::clusterExport(cl, "p9_fit_ctx", envir = environment())
        }

        dispatch <- function(idx, ks) {
            if (!use_cluster)
                return(lapply(seq_along(idx),
                              function(j) run_local(idx[j], ks[j])))
            jobs <- Map(function(j, k) list(i = j, k = as.integer(k)), idx, ks)
            parallel::parLapply(cl, jobs, p9_prss_worker)
        }

        # == Phase 1: exploration, k = 3..max_k ===============================
        # Independent iterations, run as one batch.
        n_phase1 <- min(8L, num_iterations)
        if (n_phase1 > 0L) {
            idx1 <- seq_len(n_phase1)
            ks1  <- 3L + ((idx1 - 1L) %% (max_k - 2L))
            out1 <- dispatch(idx1, ks1)
            for (j in seq_along(idx1)) apply_result(idx1[j], out1[[j]], 1L)
        }

        # Returns NULL if every phase-1 iteration failed.
        if (best_iter == 0L || is.null(best_k)) {
            cat("  WARNING: no phase-1 iteration converged for",
                region_id, gene_set_label, "-- skipping\n")
            return(NULL)
        }

        k_range <- k_window(best_k)
        cat("    best k after phase 1:", best_k,
            "| PRSS:", round(best_prss, 4), "\n")

        # == Phase 2: adaptive refinement =====================================
        # Run in chunks assuming no improvement; after an improvement the rest of
        # the chunk is discarded and finished by phase2_serial().

        phase2_end <- min(num_iterations, 8L + as.integer(PRSS_PHASE2_MAX))

        phase2_serial <- function(from_i) {
            if (phase2_end < 9L || from_i > phase2_end)
                return(invisible(NULL))
            for (i in from_i:phase2_end) {
                if (early_stop_counter >= PRSS_EARLY_STOP) {
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
                if (early_stop_counter >= PRSS_EARLY_STOP) {
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
                    if (early_stop_counter >= PRSS_EARLY_STOP) {
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

        # If the winner was fitted on a worker, refit it here (deterministic).
        if (!is.null(best_par_iter) && identical(best_iter, best_par_iter$i)) {
            cat("    best model is from iter", best_iter,
                "(k =", best_par_iter$k, ") on a worker -- refitting here\n")
            set.seed(123L + best_iter)
            best_model <- gam(make_gam_formula(best_par_iter$k, use_re),
                              data = data, weights = w,
                              method = GAM_METHOD,
                              select = GAM_SELECT, gamma = gamma)
        }

        # Only the iterations that ran are reported.
        if (length(iters_run) == 0L || all(is.na(prss_values[iters_run]))) {
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


    # == A7. Sheet-building helpers ===========================================
    build_convergence_row <- function(result, group_label) {
        cv <- extract_convergence_details(result$best_model)
        data.frame(
            Region                  = result$region,
            Group                   = group_label,
            Gene_Set                = result$gene_set,
            REML_monotonic_decrease = cv$monotonic_decrease,
            Convergence_criterion   = cv$convergence_criterion,
            Gradient_value          = cv$gradient_value,
            Gradient_norm           = cv$gradient_norm,
            Relative_score_change   = cv$relative_score_change,
            Iterations              = cv$iterations,
            Max_iterations_reached  = cv$max_iterations_reached,
            Final_REML_score        = cv$final_reml_score,
            Optimizer               = cv$optimizer,
            stringsAsFactors = FALSE)
    }

    build_prss_rows <- function(result, group_label) {
        n     <- result$actual_iters
        iters <- seq_len(n)
        get_p <- function(i, key) {
            x <- result$model_params[[i]]
            if (is.null(x)) return(NA)
            val <- x[[key]]
            if (is.null(val)) NA else val
        }
        data.frame(
            Region              = result$region,
            Group               = group_label,
            Gene_Set            = result$gene_set,
            Best_Model_Tag      = iters == result$best_iter,
            k_Value             = sapply(iters, function(i) get_p(i, "k")),
            Lambda_Value        = sapply(iters, function(i) get_p(i, "lambda")),
            PRSS_Iteration      = iters,
            PRSS_Value          = result$prss_values,
            REML_num_iterations = sapply(iters, function(i) get_p(i, "reml_n")),
            REML_min_score      = sapply(iters, function(i) get_p(i, "reml_min")),
            REML_max_score      = sapply(iters, function(i) get_p(i, "reml_max")),
            stringsAsFactors = FALSE)
    }

    build_reml_best_rows <- function(result, group_label) {
        all_ri <- result$all_reml_iters
        if (is.null(all_ri) || nrow(all_ri) == 0) return(NULL)
        best_ri <- all_ri[all_ri$is_best_model == TRUE, ]
        if (nrow(best_ri) == 0) return(NULL)
        m <- result$best_model
        hess_str <- NA_character_
        if (!is.null(m$outer.info$hess)) {
            h <- m$outer.info$hess
            if (is.matrix(h))
                hess_str <- paste(apply(h, 1, function(r)
                    paste(format(r, digits = 4, scientific = TRUE),
                          collapse = " ")), collapse = "; ")
        }
        nr <- nrow(best_ri)
        grad_vec <- rep(NA_real_, nr)
        if (!is.null(m$outer.info$grad)) {
            g <- m$outer.info$grad
            if      (is.matrix(g) && nrow(g) >= nr)
                grad_vec <- g[seq_len(nr), 1]
            else if (is.vector(g) && length(g) >= nr)
                grad_vec <- g[seq_len(nr)]
            else if (length(g) > 0)
                grad_vec[seq_len(min(length(g), nr))] <-
                    as.numeric(g)[seq_len(min(length(g), nr))]
        }
        data.frame(
            Region         = result$region,
            Group          = group_label,
            Gene_Set       = result$gene_set,
            Best_Model     = TRUE,
            REML_Iteration = best_ri$iteration,
            REML_Score     = best_ri$score,
            Lambda_Value   = best_ri$lambda,
            Gradient       = grad_vec,
            Hessian        = c(rep(NA_character_, nr - 1), hess_str),
            Hessian_Notes  = c(rep("Not calculated for intermediate steps", nr - 1),
                               "Final Hessian at convergence"),
            stringsAsFactors = FALSE)
    }

    # Null-deviance cache; keeps the null gam() out of the seeded stream.
    # Null includes the donor intercept; key includes the donor structure.
    .null_dev_cache <- new.env(hash = TRUE, parent = emptyenv())
    get_null_deviance <- function(data) {
        if ("donor" %in% names(data)) data$donor <- as_donor_factor(data$donor)
        data   <- ensure_weights(data)
        use_re <- donor_re_ok(data)
        # sum(w^2) is in the key: the null deviance depends on the weights.
        key <- paste0(nrow(data), "_", round(sum(data$Expression), 6),
                      "_", n_donor_levels(data), "_", as.integer(use_re),
                      "_", round(sum(data$w^2), 9))
        if (exists(key, envir = .null_dev_cache))
            return(get(key, envir = .null_dev_cache))
        set.seed(123)
        # The null carries the same weights as the fitted model.
        nd <- if (use_re)
            deviance(gam(Expression ~ 1 + s(donor, bs = "re"),
                         data = data, weights = w, method = GAM_METHOD))
        else
            deviance(gam(Expression ~ 1, data = data, weights = w))
        assign(key, nd, envir = .null_dev_cache)
        nd
    }

    # Deviance explained against the donor-aware null.
    de_explained <- function(model, data) {
        if (is.null(model)) return(NA_real_)
        # Without a w column, the model's own prior weights are used for the null.
        if (!("w" %in% names(data)) && !is.null(model$prior.weights) &&
            length(model$prior.weights) == nrow(data))
            data$w <- as.numeric(model$prior.weights)
        nd <- tryCatch(get_null_deviance(data), error = function(e) NA_real_)
        if (!is.finite(nd) || nd <= 0) return(NA_real_)
        as.numeric((nd - deviance(model)) / nd)
    }

    build_gam_row <- function(result, group_label) {
        m  <- result$best_model
        pc <- calculate_prss(m, result$data)
        sm <- summary(m)

        # GOI smooth statistics, by name; donor BLUPs counted separately.
        cf           <- coef(m)
        nm           <- names(cf)
        is_goi_sm    <- grepl("^s\\(GOI\\)", nm)
        is_any_sm    <- grepl("^s\\(", nm)
        smooth_coefs <- cf[is_goi_sm]

        has_re  <- model_has_re(m)
        null_dev <- get_null_deviance(result$data)

        data.frame(
            Region                     = result$region,
            Group                      = group_label,
            Gene_Set                   = result$gene_set,
            PRSS                       = pc$PRSS,
            RSS                        = pc$RSS,
            Lambda_REML                = pc$lambda,
            f_double_prime_int         = pc$f_double_prime_integral,
            k                          = result$best_k,
            Model_deviance             = deviance(m),
            Null_deviance              = null_dev,
            DE                         = de_explained(m, result$data),
            DE_vs_flat_null            = sm$dev.expl,
            Donor_RE                   = has_re,
            N_Donors                   = if (is.null(result$n_donors))
                                             n_donor_levels(result$data)
                                         else result$n_donors,
            Donor_RE_EDF               = if (has_re)
                                             smooth_stat(sm, DONOR_RE_LABEL, "edf")
                                         else NA_real_,
            Equation                   = extract_equation(m),
            Smooth_EDF                 = smooth_stat(sm, "s(GOI)", "edf"),
            Smooth_p_value             = smooth_stat(sm, "s(GOI)", "p-value"),
            Smooth_p_BH_FDR            = NA_real_,
            Total_EDF                  = sum(sm$edf),
            Parametric_DF              = sum(!is_any_sm),
            Is_Significantly_Nonlinear = NA,
            Nonzero_Basis_Functions    = sum(abs(smooth_coefs) > 1e-8),
            Total_Basis_Functions      = length(smooth_coefs),
            Smooth_Coef_Count          = length(smooth_coefs),
            stringsAsFactors = FALSE)
    }

    # Predicted_value includes the donor BLUP; Predicted_pop_value excludes it
    # (the plotted curve).
    build_obs_pred_rows <- function(result) {
        m   <- result$best_model
        ex  <- re_exclude(m)
        fit_cond <- as.numeric(predict(m, newdata = result$data))
        fit_pop  <- if (length(ex) > 0)
            as.numeric(predict(m, newdata = result$data, exclude = ex))
        else fit_cond
        data.frame(
            Region              = result$region,
            Gene_Set            = result$gene_set,
            GOI                 = result$data$GOI,
            Donor               = if ("donor" %in% names(result$data))
                                      as.character(result$data$donor)
                                  else NA_character_,
            Observed_value      = result$data$Expression,
            Predicted_value     = fit_cond,
            Predicted_pop_value = fit_pop,
            stringsAsFactors = FALSE)
    }

    build_basis_rows <- function(result) {
        m            <- result$best_model
        X            <- predict(m, type = "lpmatrix")
        # The s(donor) block of the lpmatrix is disjoint from this one.
        smooth_cols  <- grep("^s\\(GOI\\)", colnames(X))
        bm           <- X[, smooth_cols, drop = FALSE]
        colnames(bm) <- paste0("phi_", seq_len(ncol(bm)))
        data.frame(Region   = result$region,
                   Gene_Set = result$gene_set,
                   GOI      = result$data$GOI,
                   as.data.frame(bm), stringsAsFactors = FALSE)
    }


    # == A8. LOO gene contribution ============================================
    # LOO refits use the same formula and donor-aware null as the full fit.
    fit_gam_loo <- function(data, best_k, gamma = GAM_GAMMA) {
        set.seed(123)
        k_use <- max(3L, min(best_k, length(unique(data$GOI)) - 1L))
        k_use <- min(k_use, 10L)
        if ("donor" %in% names(data)) data$donor <- as_donor_factor(data$donor)
        data   <- ensure_weights(data)
        use_re <- donor_re_ok(data)
        tryCatch({
            set.seed(123)
            m <- gam(make_gam_formula(k_use, use_re),
                     data = data, weights = w,
                     method = GAM_METHOD,
                     select = GAM_SELECT, gamma = gamma)
            de_explained(m, data)
        }, error = function(e) NA_real_)
    }

    build_gene_contribution_sheet <- function(region, rd, full_results,
                                              group_label, gamma = GAM_GAMMA,
                                              cl = NULL) {
        rows <- list()
        for (gs in c("DFG", "HKG")) {
            key      <- paste0(region, "|", gs)
            full_res <- full_results[[key]]
            if (is.null(full_res) || is.null(full_res$best_model)) next

            de_full   <- de_explained(full_res$best_model, full_res$data)
            best_k    <- full_res$best_k
            goi_vec    <- rd$goi_vec
            donor_vec  <- rd$donor_vec
            # Weights recomputed if the region record lacks them.
            weight_vec <- if (is.null(rd$weight_vec))
                donor_weights(donor_vec) else rd$weight_vec

            if (gs == "DFG") {
                genes <- rd$dfg_present; mat <- rd$dfg_mat
            } else {
                genes <- rd$hkg_present; mat <- rd$hkg_mat
            }
            if (is.null(mat) || length(genes) == 0) next

            if (length(genes) < 2) {
                rows[[length(rows) + 1L]] <- data.frame(
                    Region = region, Group = group_label, Gene_Set = gs,
                    Removed_Gene = "(only 1 gene -- LOO not applicable)",
                    DE_Full = de_full, DE_LOO = NA_real_,
                    Delta_DE = NA_real_, Pct_Contribution = NA_real_,
                    stringsAsFactors = FALSE)
                next
            }

            # One independent refit per removed gene; k and the random-effect choice
            # are resolved once.
            jobs <- lapply(genes, function(g)
                list(gene_out = g, genes_in = setdiff(genes, g)))

            de_list <- if (!is.null(cl) && length(cl) >= 2L) {
                k_use <- max(3L, min(best_k, length(unique(goi_vec)) - 1L))
                k_use <- min(k_use, 10L)
                loo_probe <- data.frame(GOI = goi_vec,
                                        Expression = goi_vec,
                                        donor = donor_vec)
                p9_loo_ctx <- list(goi_vec   = goi_vec,
                                   donor_vec = as_donor_factor(donor_vec),
                                   w         = weight_vec,
                                   mat       = mat[, genes, drop = FALSE],
                                   k_use     = k_use,
                                   use_re    = donor_re_ok(loo_probe),
                                   bs        = GAM_BS,
                                   method    = GAM_METHOD,
                                   select    = GAM_SELECT,
                                   gamma     = gamma)
                parallel::clusterExport(cl, "p9_loo_ctx",
                                        envir = environment())
                parallel::parLapply(cl, jobs, p9_loo_worker)
            } else {
                lapply(jobs, function(job) {
                    loo_avg <- if (length(job$genes_in) == 1)
                        mat[, job$genes_in]
                    else
                        rowMeans(mat[, job$genes_in, drop = FALSE])
                    loo_data <- data.frame(GOI = goi_vec,
                                           Expression = loo_avg,
                                           donor = donor_vec,
                                           w = weight_vec)
                    list(gene = job$gene_out,
                         de   = fit_gam_loo(loo_data, best_k, gamma = gamma))
                })
            }

            # Rows follow the gene order.
            for (res_loo in de_list) {
                gene_out <- res_loo$gene
                de_loo   <- res_loo$de
                delta    <- de_full - de_loo
                pct_c    <- if (!is.na(de_full) && de_full > 0)
                    100 * delta / de_full else NA_real_
                cat("    LOO [", region, "|", gs, "] remove", gene_out,
                    "| DE_full:", round(de_full, 4),
                    "| DE_loo:",  round(de_loo, 4),
                    "| delta:",   round(delta, 4), "\n")
                rows[[length(rows) + 1L]] <- data.frame(
                    Region = region, Group = group_label, Gene_Set = gs,
                    Removed_Gene = gene_out,
                    DE_Full = de_full, DE_LOO = de_loo,
                    Delta_DE = delta, Pct_Contribution = pct_c,
                    stringsAsFactors = FALSE)
            }
        }
        if (length(rows) == 0) return(NULL)
        do.call(rbind, rows)
    }


    # =============================================================================
    # SECTION B: PLOT HELPER FUNCTIONS
    # =============================================================================
    # Strings drawn on graphics devices are ASCII.

    PRSS_PLOT_MAX_ITER <- 28L

    PNG_CURVE_W_IN <- 4.86
    PNG_CURVE_H_IN <- 4.00
    PNG_CURVE_RES  <- 300L

    # Axis label for the normalised expression scale.
    EXPR_UNIT_LBL <- paste0(NORM_TAG, " log1p")

    theme_gam <- function(base_size = 9) {
        theme_minimal(base_size = base_size) +
            theme(
                plot.title       = element_text(face = "bold",
                                                size = base_size + 1),
                plot.subtitle    = element_text(size = base_size - 1),
                axis.title       = element_text(face = "bold",
                                                size = base_size - 1),
                axis.text        = element_text(size = base_size - 2),
                panel.grid.major = element_line(color = "gray90"),
                panel.grid.minor = element_line(color = "gray95"),
                panel.border     = element_rect(color = "black", fill = NA,
                                                linewidth = 0.35),
                plot.background  = element_rect(fill = "white", color = NA),
                panel.background = element_rect(fill = "white", color = NA)
            )
    }

    make_prss_plot <- function(result, compact = FALSE) {
        bs     <- if (compact) 7 else 9
        n_cap  <- min(PRSS_PLOT_MAX_ITER, result$actual_iters)
        iters  <- seq_len(n_cap)
        prss_v <- result$prss_values[iters]
        df <- data.frame(Iteration = iters, PRSS = prss_v,
                         Is_Best = iters == result$best_iter)
        df <- df[!is.na(df$PRSS), ]

        if (nrow(df) > 0) {
            prss_max    <- max(df$PRSS, na.rm = TRUE)
            zero_thresh <- max(prss_max * 0.01, 1e-6)
            first_zero  <- which(df$PRSS < zero_thresh)[1]
            if (!is.na(first_zero) && first_zero > 1)
                df <- df[seq_len(first_zero - 1), ]
        }

        if (nrow(df) == 0)
            return(ggplot() +
                       annotate("text", x = 0.5, y = 0.5, size = 3,
                                label = paste(result$region, result$gene_set,
                                              "PRSS N/A")) +
                       theme_void())

        best_row  <- df[df$Is_Best, ]
        beyond_28 <- result$best_iter > PRSS_PLOT_MAX_ITER
        if (nrow(best_row) == 0) {
            best_row  <- df[which.min(df$PRSS), ]
            beyond_28 <- TRUE
        }
        x_max <- max(df$Iteration, na.rm = TRUE)

        ggplot(df, aes(x = Iteration, y = PRSS)) +
            geom_line(linewidth = 0.6, color = "gray30") +
            geom_point(aes(color = Is_Best, size = Is_Best)) +
            geom_vline(xintercept = best_row$Iteration[1],
                       linetype = "dashed",
                       color = if (beyond_28) "gray50" else "firebrick",
                       linewidth = 0.5) +
            geom_point(data = best_row, aes(x = Iteration, y = PRSS),
                       color = if (beyond_28) "gray40" else "firebrick",
                       size = 3.5, shape = 18, inherit.aes = FALSE) +
            scale_color_manual(values = c("FALSE" = "gray20", "TRUE" = "firebrick"),
                               guide = "none") +
            scale_size_manual(values  = c("FALSE" = 1.2,     "TRUE" = 3.5),
                              guide = "none") +
            scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
            coord_cartesian(xlim = c(1, x_max)) +
            labs(title    = paste0(result$region, " | ", result$gene_set,
                                   " - PRSS"),
                 subtitle = paste0("Best iter: ", result$best_iter,
                                   "  k=", result$best_k,
                                   "  PRSS=", round(result$best_prss, 3)),
                 x = "PRSS Iteration", y = "PRSS") +
            theme_gam(base_size = bs)
    }

    make_reml_plot <- function(result, compact = FALSE) {
        bs     <- if (compact) 7 else 9
        all_ri <- result$all_reml_iters
        if (is.null(all_ri) || nrow(all_ri) == 0)
            return(ggplot() +
                       annotate("text", x = 0.5, y = 0.5, size = 3,
                                label = "REML data unavailable") +
                       labs(title = paste0(result$region, " | ", result$gene_set,
                                           " - REML")) +
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
                                   " - REML"),
                 subtitle = paste0("PRSS iter ", result$best_iter,
                                   "  |  Final REML = ",
                                   round(final_row$score, 3)),
                 x = "REML Iteration", y = "REML Score") +
            theme_gam(base_size = bs)
    }

    make_gam_curve_plot <- function(region, vres, goi_label, gamma_used,
                                    region_gammas = NULL) {
        if (!is.null(region_gammas) && region %in% names(region_gammas))
            gamma_used <- region_gammas[[region]]

        gs_map          <- c(DFG = DFG_COMPOSITE_LABEL, HKG = "cHKG")
        gene_set_colors <- setNames(c("#4B0082", "#4F75DE"),
                                    c(DFG_COMPOSITE_LABEL, "cHKG"))
        gene_set_fills  <- gene_set_colors
        plot_data <- data.frame()

        for (gs in c("DFG", "HKG")) {
            res <- vres[[paste0(region, "|", gs)]]
            if (is.null(res) || is.null(res$best_model)) next
            m     <- res$best_model
            dat   <- res$data
            label <- unname(gs_map[gs])

            pts <- data.frame(GOI        = dat$GOI,
                              Expression = dat$Expression,
                              Gene_Set   = label,
                              se.fit     = NA_real_,
                              type       = "data",
                              row.names  = NULL,
                              stringsAsFactors = FALSE)

            x_seq <- seq(min(dat$GOI), max(dat$GOI), length.out = 1000)

            # Population-level curve: donor set to an observed level, then excluded.
            newd <- data.frame(GOI = x_seq)
            if ("donor" %in% names(dat)) newd$donor <- dat$donor[1]
            ex_terms <- re_exclude(m)
            pred_pop <- function() {
                if (length(ex_terms) > 0)
                    predict(m, newdata = newd, se.fit = TRUE, exclude = ex_terms)
                else
                    predict(m, newdata = newd, se.fit = TRUE)
            }

            fit_df <- tryCatch({
                pred <- pred_pop()
                data.frame(GOI        = x_seq,
                           Expression = as.numeric(pred$fit),
                           Gene_Set   = label,
                           se.fit     = as.numeric(pred$se.fit),
                           type       = "fit",
                           row.names  = NULL,
                           stringsAsFactors = FALSE)
            }, error = function(e) {
                cat("  [make_gam_curve_plot] predict() failed for",
                    region, gs, ":", conditionMessage(e), "\n")
                NULL
            }, warning = function(w) {
                pred <- suppressWarnings(pred_pop())
                data.frame(GOI        = x_seq,
                           Expression = as.numeric(pred$fit),
                           Gene_Set   = label,
                           se.fit     = as.numeric(pred$se.fit),
                           type       = "fit",
                           row.names  = NULL,
                           stringsAsFactors = FALSE)
            })

            plot_data <- rbind(plot_data, pts)
            if (!is.null(fit_df)) plot_data <- rbind(plot_data, fit_df)
        }

        if (nrow(plot_data) == 0)
            return(ggplot() +
                       annotate("text", x = 0.5, y = 0.5,
                                label = paste("No data / predict failed:", region),
                                size = 4) +
                       theme_void())

        # DE against the donor-aware null, as in the sheets.
        de_labels <- sapply(c("DFG", "HKG"), function(gs) {
            res <- vres[[paste0(region, "|", gs)]]
            if (is.null(res) || is.null(res$best_model)) return(NA_character_)
            sprintf("%s DE=%.1f%%", unname(gs_map[gs]),
                    de_explained(res$best_model, res$data) * 100)
        })
        de_labels   <- de_labels[!is.na(de_labels)]
        n_cells_dfg <- sum(plot_data$type == "data" &
                           plot_data$Gene_Set == DFG_COMPOSITE_LABEL)
        n_cells_hkg <- sum(plot_data$type == "data" &
                           plot_data$Gene_Set == "cHKG")
        n_cells     <- max(n_cells_dfg, n_cells_hkg)

        ggplot(plot_data, aes(x = GOI, y = Expression, color = Gene_Set)) +
            geom_point(data     = subset(plot_data, type == "data"),
                       position = position_jitter(width = 0.01, height = 0,
                                                  seed = 123),
                       alpha = 0.25, size = 0.8) +
            geom_line(data      = subset(plot_data, type == "fit"),
                      linewidth = 1.2) +
            geom_ribbon(data = subset(plot_data, type == "fit"),
                        aes(ymin = Expression - 1.96 * se.fit,
                            ymax = Expression + 1.96 * se.fit,
                            fill = Gene_Set),
                        alpha = 0.20, color = NA) +
            scale_color_manual(values = gene_set_colors, name = "Gene Set") +
            scale_fill_manual( values = gene_set_fills,  name = "Gene Set") +
            labs(
                title    = paste0("GAM: Gene Set vs ", goi_label,
                                  " -- ", region,
                                  "  [n=", n_cells, ", GOI>0, gamma=",
                                  gamma_used, "]"),
                subtitle = paste(de_labels, collapse = "  |  "),
                x        = paste0(goi_label, " Expression (", EXPR_UNIT_LBL,
                                  ", GOI>0)"),
                y        = paste0("Gene Set Mean Expression (", EXPR_UNIT_LBL, ")")
            ) +
            scale_x_continuous(labels = number_format(accuracy = 0.1)) +
            scale_y_continuous(labels = number_format(accuracy = 0.1)) +
            theme_gam() +
            theme(plot.title      = element_text(face = "bold", size = 11),
                  plot.subtitle   = element_text(size = 9),
                  legend.position = "right")
    }


    # =============================================================================
    # SECTION C: MASTER RUN FUNCTION
    # =============================================================================

    run_gam_pipeline <- function(tag, dataset_dir, regions, dfg_genes, hkg_genes,
                                 supertype      = PRIMARY_STYPE,
                                 gamma          = GAM_GAMMA,
                                 region_gammas  = NULL,
                                 include_pooled = TRUE,
                                 cl             = NULL) {

        group_label        <- paste0("SEAAD_", tag)
        pooled_region_name <- paste0(tag, "_Pooled")
        pooled_group_label <- paste0("SEAAD_", tag, "_Pooled")

        cat("\n", strrep("=", 70), "\n", sep = "")
        cat("  PART 9.4 -- GAM-REML-PRSS |", tag,
            "| supertype:", supertype, "| gamma:", gamma, "\n")
        cat("  Pipeline: max_k=10 | Phase1=8 iters modulo",
            "| Phase2=iters 9-100 | EarlyStop=", PRSS_EARLY_STOP, "\n")
        cat("  GOI filter: GOI > 0\n")
        cat("  Tail trim: NONE\n")
        if (!is.null(region_gammas) && length(region_gammas) > 0)
            cat("  Region gamma overrides:",
                paste(names(region_gammas), region_gammas,
                      sep = "=", collapse = ", "), "\n")
        cat("  Pooled region:", if (include_pooled) "YES" else "NO", "\n")
        cat("  GOI :", GOI, "\n")
        cat("  DFG :", paste(dfg_genes, collapse = ", "), "\n")
        cat("  HKG :", paste(hkg_genes, collapse = ", "), "\n")
        cat("  Regions:", paste(regions, collapse = ", "), "\n")
        cat(strrep("=", 70), "\n\n", sep = "")

        # == C1. Per-region data ==============================================
        cat("Building per-region data (GOI>0 filter applied) ...\n")
        region_data <- build_region_data(dataset_dir, regions, dfg_genes,
                                         hkg_genes, supertype = supertype)
        n_valid_regions <- sum(!sapply(region_data, is.null))
        cat("Regions with sufficient cells:", n_valid_regions, "/",
            length(regions), "\n\n")

        # == C1b. Pooled region ===============================================
        if (include_pooled && n_valid_regions > 0) {
            cat("Building pooled region ...\n")
            pooled_rd <- build_pooled_region(region_data, tag)
            if (!is.null(pooled_rd)) {
                region_data[[pooled_region_name]] <- pooled_rd
                cat("  Pooled region '", pooled_region_name,
                    "' added to region_data.\n\n", sep = "")
            }
        }

        # == C1c. Cell count audit ============================================
        cat("Building cell count audit sheet ...\n")
        sheet_cell_counts <- build_cell_count_sheet(
            tag, region_data, regions, pooled_region_name)

        cat("\n  Cell count audit -- GOI filter summary:\n")
        cat(sprintf("  %-22s  %8s  %8s  %8s  %8s  %8s  %s\n",
                    "Region", "N_Raw", "GOI_Neg", "GOI_Pos",
                    "Trimmed", "N_Used", "Tally"))
        cat("  ", strrep("-", 78), "\n", sep = "")
        for (i in seq_len(nrow(sheet_cell_counts))) {
            r <- sheet_cell_counts[i, ]
            cat(sprintf("  %-22s  %8s  %8s  %8s  %8s  %8s  %s\n",
                        r$Region,
                        ifelse(is.na(r$N_Raw),          "-", r$N_Raw),
                        ifelse(is.na(r$N_GOI_Neg),      "-", r$N_GOI_Neg),
                        ifelse(is.na(r$N_GOI_Pos),      "-", r$N_GOI_Pos),
                        ifelse(is.na(r$N_Trimmed),      "-", r$N_Trimmed),
                        ifelse(is.na(r$N_Used_for_GAM), "-", r$N_Used_for_GAM),
                        ifelse(is.na(r$Tally_OK), "SKIP",
                               ifelse(r$Tally_OK, "OK", "*** MISMATCH ***"))))
        }
        cat("\n")

        n_mismatches <- sum(!is.na(sheet_cell_counts$Tally_OK) &
                                !sheet_cell_counts$Tally_OK)
        if (n_mismatches > 0) {
            cat("  Cell count guardrail: ", n_mismatches,
                " mismatch(es)\n", sep = "")
            print(sheet_cell_counts[!sheet_cell_counts$Tally_OK, ])
        } else {
            cat("  Cell count audit: all regions pass",
                "(N_Used_for_GAM == N_GOI_Pos)\n")
        }
        cat("\n")

        all_fit_regions <- names(Filter(Negate(is.null), region_data))

        # == C2. GAM fits =====================================================
        cat("========== Running GAM fits:", length(all_fit_regions),
            "regions x 2 gene sets ==========\n")
        all_results <- list()

        for (reg in all_fit_regions) {
            rd        <- region_data[[reg]]
            is_pooled <- reg == pooled_region_name
            reg_gamma <- if (!is_pooled && !is.null(region_gammas) &&
                             reg %in% names(region_gammas))
                region_gammas[[reg]] else gamma

            cat("\n--- Region:", reg, if (is_pooled) "[POOLED]" else "",
                "| gamma =", reg_gamma, "---\n")

            cat(" Fitting DFG ...\n")
            key_dfg <- paste0(reg, "|DFG")
            all_results[[key_dfg]] <- tryCatch(
                fit_gam_prss(rd$dfg, reg, "DFG", gamma = reg_gamma, cl = cl),
                error = function(e) {
                    cat("  FAILED:", conditionMessage(e), "\n"); NULL })

            if (!is.null(rd$hkg)) {
                cat(" Fitting HKG ...\n")
                key_hkg <- paste0(reg, "|HKG")
                all_results[[key_hkg]] <- tryCatch(
                    fit_gam_prss(rd$hkg, reg, "HKG", gamma = reg_gamma, cl = cl),
                    error = function(e) {
                        cat("  FAILED:", conditionMessage(e), "\n"); NULL })
            } else {
                cat(" HKG skipped (no HKG data)\n")
            }
        }

        n_success <- sum(!sapply(all_results, is.null))
        cat("\nCompleted:", n_success, "fits (of max",
            length(all_fit_regions) * 2, ")\n")
        valid_results <- Filter(Negate(is.null), all_results)

        # == C3. BH FDR =======================================================
        dfg_results <- Filter(function(r) r$gene_set == "DFG", valid_results)
        hkg_results <- Filter(function(r) r$gene_set == "HKG", valid_results)

        sheet_perf_dfg <- if (length(dfg_results) > 0)
            do.call(rbind, lapply(dfg_results, function(r)
                build_gam_row(r, if (r$region == pooled_region_name)
                    pooled_group_label else group_label))) else NULL
        sheet_perf_hkg <- if (length(hkg_results) > 0)
            do.call(rbind, lapply(hkg_results, function(r)
                build_gam_row(r, if (r$region == pooled_region_name)
                    pooled_group_label else group_label))) else NULL

        if (!is.null(sheet_perf_dfg) && nrow(sheet_perf_dfg) > 0) {
            sheet_perf_dfg$Smooth_p_BH_FDR <- p.adjust(
                sheet_perf_dfg$Smooth_p_value, "BH")
            sheet_perf_dfg$Is_Significantly_Nonlinear <-
                sheet_perf_dfg$Smooth_EDF > 1 &
                sheet_perf_dfg$Smooth_p_BH_FDR < 0.05
            sheet_perf_dfg$Is_Pooled <-
                sheet_perf_dfg$Region == pooled_region_name
        }
        if (!is.null(sheet_perf_hkg) && nrow(sheet_perf_hkg) > 0) {
            sheet_perf_hkg$Smooth_p_BH_FDR <- p.adjust(
                sheet_perf_hkg$Smooth_p_value, "BH")
            sheet_perf_hkg$Is_Significantly_Nonlinear <-
                sheet_perf_hkg$Smooth_EDF > 1 &
                sheet_perf_hkg$Smooth_p_BH_FDR < 0.05
            sheet_perf_hkg$Is_Pooled <-
                sheet_perf_hkg$Region == pooled_region_name
        }

        # == C4. Remaining sheets =============================================
        sheet_conv <- do.call(rbind, lapply(valid_results, function(r)
            build_convergence_row(r, if (r$region == pooled_region_name)
                pooled_group_label else group_label)))
        sheet_prss <- do.call(rbind, lapply(valid_results, function(r)
            build_prss_rows(r, if (r$region == pooled_region_name)
                pooled_group_label else group_label)))

        reml_best_list <- Filter(Negate(is.null),
            lapply(valid_results, function(r)
                build_reml_best_rows(r,
                    if (r$region == pooled_region_name)
                        pooled_group_label else group_label)))
        sheet_reml_best <- if (length(reml_best_list) > 0)
            do.call(rbind, reml_best_list)
        else
            data.frame(Region = character(), Group = character(),
                       Gene_Set = character(), stringsAsFactors = FALSE)

        sheet_obs_pred <- do.call(rbind,
            lapply(valid_results, build_obs_pred_rows))

        basis_list <- lapply(valid_results, build_basis_rows)
        phi_counts <- sapply(basis_list, function(df)
            sum(grepl("^phi_", names(df))))
        max_phi <- if (length(phi_counts) > 0 && max(phi_counts) > 0)
            max(phi_counts) else 0L
        pad_phi <- function(df) {
            ex <- sum(grepl("^phi_", names(df)))
            if (max_phi > 0 && ex < max_phi)
                for (j in (ex + 1):max_phi)
                    df[[paste0("phi_", j)]] <- NA_real_
            df
        }
        sheet_basis <- if (length(basis_list) > 0)
            do.call(rbind, lapply(basis_list, pad_phi))
        else
            data.frame(Region = character(), Gene_Set = character(),
                       stringsAsFactors = FALSE)

        # == C5. LOO ==========================================================
        # Run only when RUN_LOO is TRUE.
        sheet_loo <- NULL
        if (!isTRUE(RUN_LOO)) {
            cat("\n========== LOO skipped (RUN_LOO = FALSE) ==========\n")
        } else {
        cat("\n========== LOO gene contribution analysis ==========\n")
        loo_list <- lapply(all_fit_regions, function(reg) {
            rd <- region_data[[reg]]
            if (is.null(rd)) return(NULL)
            is_pooled     <- reg == pooled_region_name
            reg_gamma_loo <- if (!is_pooled && !is.null(region_gammas) &&
                                 reg %in% names(region_gammas))
                region_gammas[[reg]] else gamma
            reg_group_label <- if (is_pooled) pooled_group_label else group_label
            cat("  LOO region:", reg,
                if (is_pooled) "[POOLED]" else "",
                "| gamma =", reg_gamma_loo, "\n")
            tryCatch(
                build_gene_contribution_sheet(reg, rd, all_results,
                                              reg_group_label,
                                              gamma = reg_gamma_loo,
                                              cl = cl),
                error = function(e) {
                    cat("  LOO FAILED for", reg, ":", conditionMessage(e), "\n")
                    NULL })
        })
        loo_valid <- Filter(Negate(is.null), loo_list)
        sheet_loo <- if (length(loo_valid) > 0) do.call(rbind, loo_valid) else NULL
        cat("LOO sheet rows:",
            if (is.null(sheet_loo)) 0 else nrow(sheet_loo), "\n")
        }

        # == C6. Workbook =====================================================
        safe_df <- function(x, ...) {
            if (is.null(x) || (is.data.frame(x) && nrow(x) == 0))
                data.frame(Note = paste("No data (", ..., ")"),
                           stringsAsFactors = FALSE)
            else x
        }

        xlsx_path <- file.path(PART9_DIR,
                               paste0(OUT_PREFIX, "_", tag, "_GAM_Results_1.xlsx"))
        write_xlsx(list(
            Cell_Count_Audit          = safe_df(sheet_cell_counts, "cell counts"),
            REML_Convergence          = safe_df(sheet_conv,        "convergence"),
            All_Models_PRSS           = safe_df(sheet_prss,        "PRSS"),
            REML_Iterations_BestModel = safe_df(sheet_reml_best,   "REML iters"),
            GAM_Performance_DFG       = safe_df(sheet_perf_dfg,    "DFG perf"),
            GAM_Performance_HKG       = safe_df(sheet_perf_hkg,    "HKG perf"),
            Observed_Predicted        = safe_df(sheet_obs_pred,    "obs/pred"),
            Basis_Functions           = safe_df(sheet_basis,       "basis"),
            Gene_Contribution_LOO     = safe_df(sheet_loo,         "LOO")
        ), path = xlsx_path)
        cat("\nSaved:", xlsx_path, "\n")

        # == C7. RDS ==========================================================
        rds_path <- file.path(PART9_DIR,
                              paste0(OUT_PREFIX, "_", tag, "_GAM_models.rds"))
        # The estimator is recorded so that later parts can check it.
        .re_flags <- vapply(valid_results,
                            function(r) isTRUE(r$donor_re), logical(1))
        saveRDS(list(all_results        = all_results,
                     region_data        = region_data,
                     regions            = regions,
                     donor_re           = GAM_DONOR_RE,
                     donor_re_min_levels = DONOR_RE_MIN_LEVELS,
                     donor_re_fitted    = .re_flags,
                     n_donors           = vapply(valid_results, function(r)
                         if (is.null(r$n_donors)) NA_integer_
                         else as.integer(r$n_donors), integer(1)),
                     all_fit_regions    = all_fit_regions,
                     pooled_region_name = pooled_region_name,
                     sheet_loo          = sheet_loo,
                     sheet_cell_counts  = sheet_cell_counts,
                     tag                = tag,
                     group_label        = group_label,
                     supertype          = supertype,
                     gamma              = gamma,
                     region_gammas      = region_gammas,
                     include_pooled     = include_pooled,
                     arm_tag            = ARM_TAG,
                     norm_tag           = NORM_TAG,
                     goi                = GOI,
                     dfg_genes          = dfg_genes,
                     hkg_genes          = hkg_genes),
                rds_path)
        cat("Models saved:", rds_path, "\n")

        # == C8. Plots ========================================================
        pdf_path <- file.path(PART9_DIR,
                              paste0(OUT_PREFIX, "_", tag, "_GAM_Results_2.pdf"))
        cat("\nGenerating plots -> PDF:", pdf_path, "\n")
        cat("  GAM curve PNGs -> subregion dir:", SUBREGION_DIR, "\n\n")

        safe_dev_set <- function(device_id) {
            tryCatch(grDevices::dev.set(device_id),
                     error = function(e) invisible(NULL))
        }

        emit_plot <- function(p) {
            if (is.null(p)) return(invisible(NULL))
            tryCatch(print(p), error = function(e) invisible(NULL))
            safe_dev_set(pdf_device)
            tryCatch(print(p), error = function(e) invisible(NULL))
            safe_dev_set(rstudio_device)
        }

        emit_grid <- function(g) {
            if (is.null(g)) return(invisible(NULL))
            tryCatch({ grid::grid.newpage(); grid::grid.draw(g) },
                     error = function(e) invisible(NULL))
            safe_dev_set(pdf_device)
            tryCatch({ grid::grid.newpage(); grid::grid.draw(g) },
                     error = function(e) invisible(NULL))
            safe_dev_set(rstudio_device)
        }

        rstudio_device <- grDevices::dev.cur()
        if (rstudio_device == 1L) {
            tryCatch({ grDevices::dev.new()
                       rstudio_device <- grDevices::dev.cur() },
                     error = function(e) invisible(NULL))
        }

        pdf(pdf_path, width = 14, height = 7)
        pdf_device <- grDevices::dev.cur()
        safe_dev_set(rstudio_device)

        # PRSS and REML grid pages, PDF only
        all_panel_pairs <- list()
        for (key in names(valid_results)) {
            res    <- valid_results[[key]]
            p_prss <- tryCatch(make_prss_plot(res, compact = TRUE),
                               error = function(e) NULL)
            p_reml <- tryCatch(make_reml_plot(res, compact = TRUE),
                               error = function(e) NULL)
            all_panel_pairs[[key]] <- list(prss = p_prss, reml = p_reml, key = key)
        }
        fits_pp  <- 4L
        n_pages  <- ceiling(length(all_panel_pairs) / fits_pp)
        keys_vec <- names(all_panel_pairs)

        for (pg in seq_len(n_pages)) {
            idx_start <- (pg - 1L) * fits_pp + 1L
            idx_end   <- min(pg * fits_pp, length(all_panel_pairs))
            page_keys <- keys_vec[idx_start:idx_end]
            panel_list <- list()
            for (k in page_keys) {
                pp <- all_panel_pairs[[k]]
                panel_list <- c(panel_list,
                    list(if (!is.null(pp$prss)) pp$prss else
                         ggplot() + theme_void() +
                             annotate("text", 0.5, 0.5,
                                      label = paste(k, "PRSS N/A"))),
                    list(if (!is.null(pp$reml)) pp$reml else
                         ggplot() + theme_void() +
                             annotate("text", 0.5, 0.5,
                                      label = paste(k, "REML N/A"))))
            }
            g <- tryCatch(
                gridExtra::arrangeGrob(
                    grobs = panel_list,
                    ncol  = length(page_keys), nrow = 2L,
                    top   = grid::textGrob(
                        paste0("PRSS (first 28 iters) & REML -- ", tag,
                               "  (page ", pg, " of ", n_pages, ")"),
                        gp = grid::gpar(fontface = "bold", fontsize = 11))),
                error = function(e) {
                    cat("  grid error page", pg, ":", conditionMessage(e), "\n")
                    NULL })
            emit_grid(g)
        }

        # GAM curve plots, PDF plus one PNG per region
        cat("  GAM curve plots (PDF and PNG per region) ...\n")
        for (reg in all_fit_regions) {
            has_result <- any(sapply(c("DFG", "HKG"), function(gs)
                !is.null(valid_results[[paste0(reg, "|", gs)]])))
            if (!has_result) next

            p_curve <- tryCatch(
                make_gam_curve_plot(reg, valid_results, GOI,
                                    gamma_used    = gamma,
                                    region_gammas = region_gammas),
                error = function(e) {
                    cat("  Curve plot error", reg, ":", conditionMessage(e), "\n")
                    NULL })
            emit_plot(p_curve)

            if (!is.null(p_curve)) {
                reg_safe <- gsub("[^A-Za-z0-9_]", "_", reg)
                png_path <- file.path(SUBREGION_DIR,
                    paste0(OUT_PREFIX, "_", tag, "_GAMcurve_", reg_safe, ".png"))
                tryCatch({
                    png(png_path,
                        width  = round(PNG_CURVE_W_IN * PNG_CURVE_RES),
                        height = round(PNG_CURVE_H_IN * PNG_CURVE_RES),
                        res    = PNG_CURVE_RES)
                    print(p_curve)
                    dev.off()
                    cat("  PNG saved:", basename(png_path), "\n")
                }, error = function(e) {
                    cat("  WARNING: PNG failed for", reg,
                        ":", conditionMessage(e), "\n")
                    try(dev.off(), silent = TRUE)
                })
                safe_dev_set(rstudio_device)
            }
        }

        # LOO heatmap and bar plots, PDF only
        cat("  LOO plots (PDF only) ...\n")
        if (!is.null(sheet_loo) && nrow(sheet_loo) > 0) {
            loo_clean <- sheet_loo[!is.na(sheet_loo$Delta_DE), ]
            if (nrow(loo_clean) > 0) {
                loo_clean$Region_short <- loo_clean$Region
                p_heat <- ggplot(loo_clean,
                    aes(x = Removed_Gene, y = Region_short, fill = Delta_DE * 100)) +
                    geom_tile(color = "white", linewidth = 0.5) +
                    geom_text(aes(label = sprintf("%.1f", Delta_DE * 100)),
                              size = 2.6, color = "black") +
                    facet_wrap(~Gene_Set, scales = "free_x") +
                    scale_fill_distiller(palette = "Spectral",
                        name   = expression(Delta * "DE (%)"),
                        labels = number_format(accuracy = 0.1)) +
                    labs(title    = paste("LOO Gene Contribution --", tag),
                         subtitle = paste0("dDE = DE(full) - DE(LOO)",
                                           " | red = gene helps,",
                                           " blue = gene suppresses"),
                         x = "Removed Gene", y = "Region") +
                    theme_gam() +
                    theme(axis.text.x = element_text(angle = 35, hjust = 1,
                                                     size = 7),
                          strip.text  = element_text(face = "bold", size = 10))
                emit_plot(p_heat)
            }

            bar_plots <- list()
            for (reg in unique(sheet_loo$Region)) {
                dd_loo <- sheet_loo[sheet_loo$Region == reg &
                                        !is.na(sheet_loo$Delta_DE), ]
                if (nrow(dd_loo) == 0) next
                gene_order <- dd_loo %>%
                    group_by(Removed_Gene) %>%
                    summarise(mean_delta = mean(Delta_DE, na.rm = TRUE),
                              .groups = "drop") %>%
                    arrange(desc(mean_delta)) %>%
                    pull(Removed_Gene)
                dd_loo$Removed_Gene <- factor(dd_loo$Removed_Gene,
                                              levels = gene_order)
                bar_plots[[reg]] <- ggplot(dd_loo,
                    aes(x = Removed_Gene, y = Delta_DE * 100, fill = Gene_Set)) +
                    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
                    geom_hline(yintercept = 0, linewidth = 0.4, color = "black") +
                    geom_text(aes(label = sprintf("%.1f%%", Delta_DE * 100),
                                  vjust = ifelse(Delta_DE >= 0, -0.4, 1.2)),
                              position = position_dodge(width = 0.7),
                              size = 2.2, color = "gray20") +
                    scale_fill_manual(
                        values = c("DFG" = "#4B0082", "HKG" = "#4F75DE"),
                        name = "Gene Set") +
                    labs(title = reg, x = NULL, y = expression(Delta * "DE (%)")) +
                    theme_gam(base_size = 7) +
                    theme(legend.position = "none",
                          axis.text.x = element_text(angle = 30, hjust = 1))
            }
            if (length(bar_plots) > 0) {
                g_bars <- tryCatch(
                    gridExtra::arrangeGrob(
                        grobs = bar_plots, ncol = 2,
                        top   = grid::textGrob(
                            paste0("LOO Gene Contribution -- ", tag,
                                   "\ndDE = DE(full) - DE(LOO);",
                                   " positive = gene contributes"),
                            gp = grid::gpar(fontface = "bold", fontsize = 10))),
                    error = function(e) {
                        cat("  LOO grid error:", conditionMessage(e), "\n")
                        NULL })
                emit_grid(g_bars)
            }
        }

        safe_dev_set(pdf_device)
        dev.off()
        safe_dev_set(rstudio_device)
        cat("PDF saved:", pdf_path, "\n")

        # == C9. Console summary ==============================================
        cat("\n========== Summary statistics --", tag, "==========\n")
        for (gs in c("DFG", "HKG")) {
            gs_res <- Filter(function(r) r$gene_set == gs, valid_results)
            if (length(gs_res) == 0) next

            gs_regional <- Filter(function(r) r$region != pooled_region_name,
                                  gs_res)
            gs_pooled   <- Filter(function(r) r$region == pooled_region_name,
                                  gs_res)

            summarise_group <- function(rlist, label) {
                if (length(rlist) == 0) return(invisible(NULL))
                k_vals    <- sapply(rlist, function(r) r$best_k)
                prss_vals <- sapply(rlist, function(r) r$best_prss)
                de_vals   <- sapply(rlist, function(r)
                    if (!is.null(r$best_model))
                        de_explained(r$best_model, r$data)
                    else NA_real_)
                # The s(GOI) row, not row 1 by position.
                edf_vals  <- sapply(rlist, function(r) {
                    if (!is.null(r$best_model))
                        smooth_stat(summary(r$best_model), "s(GOI)", "edf")
                    else NA_real_ })
                don_vals  <- sapply(rlist, function(r)
                    if (is.null(r$n_donors)) NA_integer_ else r$n_donors)
                n_re      <- sum(vapply(rlist,
                                        function(r) isTRUE(r$donor_re),
                                        logical(1)))
                cat("\nGene Set:", gs, "(", label, ":",
                    length(rlist), "region(s))\n")
                cat("  Mean best k:     ",
                    round(mean(k_vals,    na.rm = TRUE), 2), "\n")
                cat("  Mean best PRSS:  ",
                    round(mean(prss_vals, na.rm = TRUE), 4), "\n")
                cat("  Mean DE:         ",
                    round(mean(de_vals,   na.rm = TRUE) * 100, 2),
                    "%  (vs donor-aware null)\n")
                cat("  Mean smooth EDF: ",
                    round(mean(edf_vals,  na.rm = TRUE), 4), "\n")
                cat("  Donor RE fitted: ", n_re, "/", length(rlist),
                    " | mean donors: ",
                    round(mean(don_vals, na.rm = TRUE), 1), "\n", sep = "")
            }

            summarise_group(gs_regional, "per-region")
            summarise_group(gs_pooled,   "pooled")
        }

        cat("\n========== Part 9.4", tag, "complete ==========\n")
        cat("  Excel :", xlsx_path, "\n")
        cat("  RDS   :", rds_path,  "\n")
        cat("  PDF   :", pdf_path,  "\n")
        cat("  PNGs  :", SUBREGION_DIR, "\n\n")

        invisible(list(all_results        = all_results,
                       region_data        = region_data,
                       valid_results      = valid_results,
                       pooled_region_name = pooled_region_name,
                       sheet_loo          = sheet_loo,
                       sheet_cell_counts  = sheet_cell_counts))
    }


    # =============================================================================
    # SECTION D: EXECUTE
    # =============================================================================

    cat("\n\n", strrep("#", 70), "\n", sep = "")
    cat("  RUNNING: NeoC  |  gamma = 1.5  |  no tail trimming\n")
    cat(strrep("#", 70), "\n", sep = "")
    neoc_final <- run_gam_pipeline(
        tag            = "NeoC",
        dataset_dir    = NEOC_DIR,
        regions        = NEOC_REGIONS,
        dfg_genes      = DFG_NEOC,
        hkg_genes      = HKG_USE,
        supertype      = PRIMARY_STYPE,
        gamma          = 1.5,
        region_gammas  = NULL,
        include_pooled = TRUE,
        cl             = .p9_cl
    )

    cat("\n\n", strrep("#", 70), "\n", sep = "")
    cat("  RUNNING: NonNeoC  |  gamma = 1.5  |  no tail trimming\n")
    cat(strrep("#", 70), "\n", sep = "")
    nonnc_final <- run_gam_pipeline(
        tag            = "NonNeoC",
        dataset_dir    = NONNC_DIR,
        regions        = NONNC_REGIONS,
        dfg_genes      = DFG_NONNC,
        hkg_genes      = HKG_USE,
        supertype      = PRIMARY_STYPE,
        gamma          = 1.5,
        region_gammas  = NULL,
        include_pooled = TRUE,
        cl             = .p9_cl
    )


    # =============================================================================
    # SECTION E: GLOBAL CELL COUNT GUARDRAIL
    # =============================================================================

    cat("\n", strrep("#", 70), "\n", sep = "")
    cat("  SECTION E -- GLOBAL CELL COUNT GUARDRAIL\n")
    cat(strrep("#", 70), "\n\n", sep = "")

    all_counts <- rbind(
        if (!is.null(neoc_final$sheet_cell_counts))
            neoc_final$sheet_cell_counts  else NULL,
        if (!is.null(nonnc_final$sheet_cell_counts))
            nonnc_final$sheet_cell_counts else NULL
    )

    if (is.null(all_counts) || nrow(all_counts) == 0) {
        cat("  WARNING: no cell count audit data found.\n\n")
        log_result("GR-06", "Global cell count audit", "WARN",
                   detail = "no audit rows produced")
    } else {
        cat("  Full cell count audit across both datasets:\n\n")
        cat(sprintf("  %-8s  %-22s  %-7s  %8s  %8s  %8s  %8s  %8s  %s\n",
                    "Tag", "Region", "Pooled",
                    "N_Raw", "GOI_Neg", "GOI_Pos",
                    "Trimmed", "N_Used", "Tally"))
        cat("  ", strrep("-", 96), "\n", sep = "")
        for (i in seq_len(nrow(all_counts))) {
            r <- all_counts[i, ]
            cat(sprintf("  %-8s  %-22s  %-7s  %8s  %8s  %8s  %8s  %8s  %s\n",
                        r$Tag,
                        r$Region,
                        ifelse(isTRUE(r$Is_Pooled), "YES", "no"),
                        ifelse(is.na(r$N_Raw),          "-", r$N_Raw),
                        ifelse(is.na(r$N_GOI_Neg),      "-", r$N_GOI_Neg),
                        ifelse(is.na(r$N_GOI_Pos),      "-", r$N_GOI_Pos),
                        ifelse(is.na(r$N_Trimmed),      "-", r$N_Trimmed),
                        ifelse(is.na(r$N_Used_for_GAM), "-", r$N_Used_for_GAM),
                        ifelse(is.na(r$Tally_OK), "SKIP",
                               ifelse(r$Tally_OK, "OK",
                                      "*** MISMATCH ***"))))
        }

        failed_rows <- all_counts[!is.na(all_counts$Tally_OK) &
                                      !all_counts$Tally_OK, ]
        cat("\n")

        if (nrow(failed_rows) == 0) {
            cat(strrep("*", 70), "\n", sep = "")
            cat("  Guardrail passed: N_Used_for_GAM == N_GOI_Pos for every region\n")
            cat(strrep("*", 70), "\n\n", sep = "")
        } else {
            cat(strrep("!", 70), "\n", sep = "")
            cat("  Guardrail failed\n")
            print(failed_rows, row.names = FALSE)
            cat(strrep("!", 70), "\n\n", sep = "")
            warning("Cell count guardrail failed for ",
                    nrow(failed_rows), " region(s).", call. = FALSE)
        }

        log_result("GR-06", "Global cell count audit",
                   if (nrow(failed_rows) == 0) "PASS" else "FAIL",
                   detail = sprintf("%d row(s) audited, %d mismatch(es)",
                                    nrow(all_counts), nrow(failed_rows)))
    }

    cat("  Log written to :", LOG_PATH, "\n")
    cat("  Subregion PNGs :", SUBREGION_DIR, "\n")


    # =============================================================================
    # SECTION F: POST-GAM QUALITY GUARDRAILS
    # =============================================================================
    cat("\n", strrep("#", 70), "\n", sep = "")
    cat("  SECTION F -- POST-GAM QUALITY GUARDRAILS\n")
    cat(strrep("#", 70), "\n\n", sep = "")

    .all_valid <- c(
        if (!is.null(neoc_final$valid_results))  neoc_final$valid_results  else list(),
        if (!is.null(nonnc_final$valid_results)) nonnc_final$valid_results else list()
    )

    if (length(.all_valid) > 0) {

        # == All GAMs converged ===============================================
        .converged   <- sapply(.all_valid, function(r) {
            if (is.null(r$best_model)) return(NA)
            isTRUE(r$best_model$converged)
        })
        .n_conv      <- sum(.converged, na.rm = TRUE)
        .n_fail_conv <- sum(!.converged, na.rm = TRUE)
        .n_na_conv   <- sum(is.na(.converged))

        log_result("GR-07", "GAM convergence",
                   if (.n_fail_conv == 0) "PASS" else "WARN",
                   detail = sprintf("%d/%d converged, %d failed, %d skipped",
                                    .n_conv, length(.converged),
                                    .n_fail_conv, .n_na_conv))

        # == Deviance explained within [0, 1] =================================
        .de_vals  <- sapply(.all_valid, function(r) {
            if (is.null(r$best_model)) return(NA_real_)
            tryCatch(de_explained(r$best_model, r$data),
                     error = function(e) NA_real_)
        })
        .de_valid <- .de_vals[!is.na(.de_vals)]
        if (length(.de_valid) > 0) {
            .de_range <- range(.de_valid)
            log_result("GR-08", "Deviance explained within range",
                       if (all(.de_valid >= 0 & .de_valid <= 1)) "PASS" else "WARN",
                       detail = sprintf(
                           "DE range [%.2f%%, %.2f%%] over %d models, vs donor-aware null",
                           .de_range[1] * 100, .de_range[2] * 100,
                           length(.de_valid)))
        }

        # == PRSS finite and positive =========================================
        .prss_vals  <- sapply(.all_valid, function(r) r$best_prss)
        .prss_valid <- .prss_vals[!is.na(.prss_vals)]
        if (length(.prss_valid) > 0) {
            .prss_ok <- all(is.finite(.prss_valid) & .prss_valid > 0)
            log_result("GR-09", "PRSS finite and positive",
                       if (.prss_ok) "PASS" else "FAIL",
                       detail = sprintf("range [%.4f, %.4f] over %d models",
                                        min(.prss_valid), max(.prss_valid),
                                        length(.prss_valid)))
        }

        # == Best k within the configured range ===============================
        .k_vals  <- sapply(.all_valid, function(r) r$best_k)
        .k_valid <- .k_vals[!is.na(.k_vals)]
        if (length(.k_valid) > 0) {
            .k_range <- range(.k_valid)
            log_result("GR-10", "Best k within the configured range",
                       if (all(.k_valid >= 3 & .k_valid <= GAM_MAX_K))
                           "PASS" else "WARN",
                       detail = sprintf("best k range [%d, %d], max_k = %d, %d models",
                                        .k_range[1], .k_range[2], GAM_MAX_K,
                                        length(.k_valid)))
        }

        # == DFG DE exceeds HKG DE on average =================================
        .dfg_de <- sapply(Filter(function(r) r$gene_set == "DFG", .all_valid),
                          function(r) {
            if (is.null(r$best_model)) NA_real_
            else tryCatch(de_explained(r$best_model, r$data),
                          error = function(e) NA_real_)
        })
        .hkg_de <- sapply(Filter(function(r) r$gene_set == "HKG", .all_valid),
                          function(r) {
            if (is.null(r$best_model)) NA_real_
            else tryCatch(de_explained(r$best_model, r$data),
                          error = function(e) NA_real_)
        })
        if (length(.dfg_de) > 0 && length(.hkg_de) > 0) {
            .dfg_mean <- mean(.dfg_de, na.rm = TRUE)
            .hkg_mean <- mean(.hkg_de, na.rm = TRUE)
            log_result("GR-11", "Mean DFG deviance explained exceeds HKG",
                       if (.dfg_mean > .hkg_mean) "PASS" else "WARN",
                       detail = sprintf(
                           "DFG %.2f%% vs HKG %.2f%%, vs donor-aware null",
                           .dfg_mean * 100, .hkg_mean * 100))
        }

        # == Donor random intercept fitted on every model =====================
        # Fits with too few donors fall back to the cell-level formula (WARN).
        .re_fitted <- vapply(.all_valid, function(r) isTRUE(r$donor_re),
                             logical(1))
        .n_don_all <- vapply(.all_valid, function(r)
            if (is.null(r$n_donors)) NA_integer_ else as.integer(r$n_donors),
            integer(1))
        .fallbacks <- names(.all_valid)[!.re_fitted]
        if (!isTRUE(GAM_DONOR_RE)) {
            log_result("GR-12", "Donor random intercept", "WARN",
                       detail = paste0(
                           "GAM_DONOR_RE is FALSE: all ", length(.all_valid),
                           " fit(s) are cell level."))
        } else {
            log_result("GR-12", "Donor random intercept fitted on every model",
                       if (length(.fallbacks) == 0) "PASS" else "WARN",
                       expected = length(.all_valid),
                       observed = sum(.re_fitted),
                       detail = sprintf(
                           "donors per fit [%s, %s]; %d fallback(s) to cell level%s",
                           if (all(is.na(.n_don_all))) "NA" else
                               format(min(.n_don_all, na.rm = TRUE)),
                           if (all(is.na(.n_don_all))) "NA" else
                               format(max(.n_don_all, na.rm = TRUE)),
                           length(.fallbacks),
                           if (length(.fallbacks) == 0) "" else
                               paste0(": ", paste(.fallbacks, collapse = ", "))))
        }

    } else {
        log_result("GR-07", "No valid GAM results to check", "WARN")
    }


    # =============================================================================
    # CLOSE
    # =============================================================================

    write_guardrail_outputs()

    cat(sprintf("\n  Final ledger: %d PASS / %d WARN / %d FAIL\n",
                .chk_pass, .chk_warn, .chk_fail))
    if (.chk_fail > 0) {
        cat("  Failures are present.\n")
    } else {
        cat("  No failures recorded.\n")
    }

    # The cluster is released here; on.exit covers an earlier abort.
    if (!is.null(.p9_cl)) {
        try(parallel::stopCluster(.p9_cl), silent = TRUE)
        .p9_cl <- NULL
        cat("  PSOCK cluster stopped.\n")
    }

    cat("\n  Outputs in:", PART9_DIR, "\n")
    cat("  Arm complete:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

    if (sink.number(type = "message") != 2L) sink(type = "message")
    while (sink.number() > 0) sink()

    invisible(list(arm = arm_tag, dir = PART9_DIR,
                   pass = .chk_pass, warn = .chk_warn, fail = .chk_fail))
}


# =============================================================================
# SECTION 0c: Driver
# =============================================================================
# Sourcing this file runs every arm in CEP_ARM_TABLE; set CEP_ARM first to
# run one. The loop lives in CEP_arms.R.

cep_run_arms(PART_TAG, run_one_arm)
