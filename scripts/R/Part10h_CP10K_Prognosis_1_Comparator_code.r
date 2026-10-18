# =============================================================================
# Part10h_CP10K_Prognosis_1_Comparator_code.r
# SEA-AD CEP-IP/CMRF, Part 10h, comparator (Figure 10A and 10B, Table S12)
#   Compares the donor-level CMRF module Post-IP percentage with five external
#   baselines on the same donors, covariates and outcomes: B1 cluster
#   proportion, B2 AUCell, B3 segmented changepoint, B4 UCell and
#   B5 AddModuleScore. Marginal models, encompassing-test duels, and donor-level
#   stress tests (block sweep, bootstrap, leave-one-donor-out).
#
# Usage (from the SEA-AD folder, or set SEAAD_ROOT), after Part 10h Prognosis
#   source("Part10h_CP10K_Prognosis_1_Comparator_code.r")
#
# Outputs (Part10h_Prognosis_1_Comparator/<arm>/): workbook, results RDS,
#   figures (PDF, PNG), guardrail log and console log.
# =============================================================================

# =============================================================================
# SECTION 0: Libraries and shared arm definitions
# =============================================================================

suppressPackageStartupMessages({
    library(Matrix)
    library(writexl)
    library(ggplot2)
    library(parallel)
})
.HAVE_SEGMENTED <- requireNamespace("segmented", quietly = TRUE)
.HAVE_AUCELL    <- requireNamespace("AUCell", quietly = TRUE)
.HAVE_UCELL     <- requireNamespace("UCell", quietly = TRUE)

set.seed(42)

if (!exists("CEP_HOME")) CEP_HOME <- Sys.getenv("SEAAD_ROOT", getwd())
.p10hc_harness <- file.path(CEP_HOME, "CEP_arms.R")
if (!file.exists(.p10hc_harness))
    stop("CEP_arms.R not found at: ", .p10hc_harness,
         "\n  It holds the arm table shared by Parts 9 onwards.", call. = FALSE)
source(.p10hc_harness)

PART_TAG      <- "Part10h"
PART_DIR_TAG  <- "Part10h_Prognosis_1_Comparator"
P10H_DIR_TAG  <- "Part10h_Prognosis"

# == B1 denominator ===========================================================
# Microglial and PVM supertypes; Monocyte and Lymphocyte are excluded.
MICROGLIAL_SUPERTYPES <- c("Micro-PVM_1", "Micro-PVM_2",
                           "Micro-PVM_2_1-SEAAD", "Micro-PVM_2_3-SEAAD",
                           "Micro-PVM_3-SEAAD", "Micro-PVM_4-SEAAD")

# == AUCell ===================================================================
AUCELL_TOP_FRAC <- 0.05   # AUCell's default rank threshold
AUCELL_CHECK_N  <- 500L   # cells sampled for the package cross-check
AUCELL_CELL_BLOCK <- 2000L # cells per parallel chunk

# == UCell ====================================================================
# UCell default; a gene ranked beyond this adds a fixed penalty.
UCELL_MAXRANK   <- 1500L

# == AddModuleScore ===========================================================
# Seurat defaults (Scanpy score_genes uses 25 bins and 50 controls).
AMS_NBIN        <- 24L
AMS_CTRL        <- 100L
AMS_SEED        <- 42L

# == Parallelism ==============================================================
# PSOCK workers; two cores are left for the session.
PAR_ENABLE  <- TRUE
# detectCores() can return NA.
.p10hc_cores <- suppressWarnings(parallel::detectCores(logical = TRUE))
if (!is.finite(.p10hc_cores)) .p10hc_cores <- 2L
PAR_WORKERS <- as.integer(max(1L, min(.p10hc_cores - 2L, 16L)))

# == Changepoint ==============================================================
SEG_MAX_CELLS <- 40000L   # cells subsampled for segmented()
SEG_SEED      <- 42L

# == Reading the dataset matrix ===============================================
RAW_MAX_CELLS <- 150000L
NORM_TOL      <- 1e-6

# == Modelling ================================================================
# Outcome classes: control, technical (donor means) and clinical.
# Headline counts are clinical only.
TECH_PATTERN    <- "\\(donor mean\\)$"

PRIMARY_BLOCK   <- "cov4"   # the block every headline row is fitted on
MIN_DONORS      <- 20L      # a model below this is not fitted
BOOT_N          <- 10000L   # donor bootstrap replicates
# The bootstrap runs on the clinical outcomes with the largest module beta.
BOOT_MAX_OUTCOMES <- 8L
VIF_READABLE    <- 5        # joint VIF above this: coefficients not separable
ALPHA           <- 0.05

# == Predictor labels and colours =============================================
PRED_KEYS <- c(module = "CMRF module Post-IP",
               B1 = "Cluster proportion (supertype)",
               B2 = "Module score (AUCell)",
               B3 = "Changepoint split (segmented)",
               B4 = "Module score (UCell)",
               B5 = "Module score (AddModuleScore)")
PRED_FILL <- c("CMRF module Post-IP" = "#2D2D2D",
               "Cluster proportion (supertype)" = "#B2182B",
               "Module score (AUCell)" = "#1B3A6B",
               "Changepoint split (segmented)" = "#7B2D8E",
               "Module score (UCell)" = "#1B7837",
               "Module score (AddModuleScore)" = "#E08214")
# The three module scores.
SCORE_KEYS <- c("B2", "B4", "B5")


run_one_arm <- function(arm_tag, config_file, out_dir, file_tag) {

    config_path <- file.path(CEP_HOME, config_file)
    if (!file.exists(config_path))
        stop("Configuration not found: ", config_path, call. = FALSE)

    .p10hc_preamble <- character(0)
    assign("CEP_ARM", arm_tag, envir = globalenv())
    on.exit(suppressWarnings(rm("CEP_ARM", envir = globalenv())), add = TRUE)

    .p10hc_preamble <- tryCatch(
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

    .p10hc_needed <- c("ARM_TAG", "R_EXPORT_ROOT", "GOI", "DFG", "NORM_TAG",
                         "PRIMARY_STYPE", "NEOC_DIR", "NONNC_DIR",
                        "SHARED_DIR", "NORM_TARGET_SUM", "STYPE_COL",
                        "DONOR_COL")
    .p10hc_absent <- .p10hc_needed[!vapply(.p10hc_needed, exists, logical(1))]
    if (length(.p10hc_absent) > 0)
        stop("The configuration loaded from\n    ", config_path,
             "\n  does not define: ", paste(.p10hc_absent, collapse = ", "),
             call. = FALSE)
    if (!identical(ARM_TAG, arm_tag))
        stop("The configuration reports ARM_TAG '", ARM_TAG,
             "' but this arm is '", arm_tag, "'.", call. = FALSE)

    .p10hc_preamble <- c(.p10hc_preamble, utils::capture.output(
        cat("  Config OK: ", ARM_TAG, " | ", GOI, " + ", length(DFG),
            " DFG | ", NORM_TAG, "\n", sep = "")))

    DATASETS <- c("NeoC", "NonNeoC")

    # == Paths ================================================================
    OUT_DIR <- out_dir
    dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
    FILE_PREFIX <- cep_out_prefix(PART_TAG, NORM_TAG, cep_file_tag_for(arm_tag))

    # Part 10h's results object, written beside its workbook.
    P10H_DIR <- cep_arm_dir_for(P10H_DIR_TAG, arm_tag)
    P10H_RDS <- file.path(P10H_DIR,
                          paste0(FILE_PREFIX, "_Prognosis_results.rds"))
    ds_dir <- function(ds) if (identical(ds, "NeoC")) NEOC_DIR else NONNC_DIR

    # == Logging ==============================================================
    LOG_PATH_10HC <- file.path(OUT_DIR, paste0(FILE_PREFIX, "_console_log_",
                                               format(Sys.time(),
                                                      "%Y%m%d_%H%M%S"), ".txt"))
    MSG_PATH_10HC <- file.path(OUT_DIR, paste0(FILE_PREFIX,
                                               "_console_log_messages.tmp"))
    while (sink.number() > 0) sink()
    if (sink.number(type = "message") != 2L) sink(type = "message")
    .p10hc_log_con <- file(LOG_PATH_10HC, open = "wt")
    .p10hc_msg_con <- file(MSG_PATH_10HC, open = "wt")
    sink(.p10hc_log_con, append = FALSE, split = TRUE)
    sink(.p10hc_msg_con, type = "message")
    on.exit({
        tryCatch({
            if (!isTRUE(.p10hc_completed))
                cat("\n", strrep("!", 72),
                    "\n  Arm did not complete; the error is printed to the console.\n",
                    strrep("!", 72), "\n", sep = "")
            if (sink.number(type = "message") != 2L) sink(type = "message")
            while (sink.number() > 0) sink()
            close(.p10hc_msg_con); close(.p10hc_log_con)
            if (file.exists(MSG_PATH_10HC)) {
                .msg <- readLines(MSG_PATH_10HC, warn = FALSE)
                if (length(.msg) > 0)
                    cat(c("", strrep("=", 72),
                          "  MESSAGE STREAM (warnings, messages, errors)",
                          strrep("=", 72), .msg, ""),
                        sep = "\n", file = LOG_PATH_10HC, append = TRUE)
                unlink(MSG_PATH_10HC)
            }
        }, error = function(e) invisible(NULL))
    }, add = TRUE)

    .p10hc_completed <- FALSE
    .p10hc_warn_opt  <- getOption("warn")
    .p10hc_pm_opt    <- getOption("warnPartialMatchDollar")
    options(warn = 1)
    options(warnPartialMatchDollar = TRUE)
    on.exit({
        options(warn = .p10hc_warn_opt)
        options(warnPartialMatchDollar = .p10hc_pm_opt)
    }, add = TRUE)

    cat(sprintf("Part10h Prognosis Comparator - log started: %s\n\n", Sys.time()))
    if (length(.p10hc_preamble) > 0)
        cat(paste(.p10hc_preamble, collapse = "\n"), "\n\n", sep = "")

    cat("\n", strrep("=", 72), "\n")
    cat("  Part10h Comparator - external baselines against the CMRF module\n")
    cat("  Arm          :", ARM_TAG, "\n")
    cat("  GOI          :", GOI, "|", NORM_TAG, "\n")
    cat("  Primary input:", basename(P10H_RDS), "\n")
    cat("  Also read    : obs_dataset.csv and the dataset raw matrix\n")
    cat("  Unit         : donor\n")
    cat("  Output dir   :", OUT_DIR, "\n")
    cat(strrep("=", 72), "\n\n")

    INCOMPLETE_MARKER_10HC <- file.path(OUT_DIR, "_RUN_INCOMPLETE.txt")
    writeLines(c("This arm did not finish.", "",
                 paste("Arm     :", arm_tag),
                 paste("Part    :", PART_DIR_TAG),
                 paste("Started :", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), "",
                 "Any files in this folder are from a partial run.",
                 "This file is removed automatically when the arm completes."),
               INCOMPLETE_MARKER_10HC)

    # =============================================================================
    # GUARDRAIL LEDGER
    # =============================================================================
    GR_LOG_10HC <- file.path(OUT_DIR, paste0(FILE_PREFIX, "_Guardrail_Log.txt"))
    .gr_log <- list(); .gr_pass <- 0L; .gr_fail <- 0L; .gr_warn <- 0L
    gr10hc <- function(id, status, detail) {
        # A zero-length detail would vanish from console and ledger.
        if (length(detail) != 1L || is.na(detail))
            detail <- paste0("[detail unavailable: the message was built from",
                             " a missing value and came out zero-length]")
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
    .safe_vif <- function(fit) {
        X <- try(stats::model.matrix(fit)[, -1, drop = FALSE], silent = TRUE)
        if (inherits(X, "try-error") || ncol(X) < 2L) return(rep(NA_real_, 1L))
        v <- vapply(seq_len(ncol(X)), function(j) {
            r2 <- suppressWarnings(summary(stats::lm(X[, j] ~ X[, -j,
                                                                drop = FALSE]))$r.squared)
            if (!is.finite(r2) || r2 >= 1) NA_real_ else 1 / (1 - r2)
        }, numeric(1))
        stats::setNames(v, colnames(X))
    }

    # CP10K log1p normalisation; X is cells x genes raw counts, full gene set.
    p10hc_cp10k <- function(X, target = NORM_TARGET_SUM) {
        X <- as(X, "CsparseMatrix")
        tot <- Matrix::rowSums(X); tot[tot == 0] <- 1
        Y <- X
        Y@x <- log1p(Y@x / tot[Y@i + 1L] * target)
        Y
    }

    # =============================================================================
    # PARALLELISM
    # =============================================================================
    # One PSOCK cluster per arm, stopped on exit.
    .cl <- NULL
    p10hc_cluster <- function() {
        if (!isTRUE(PAR_ENABLE) || !isTRUE(PAR_WORKERS >= 2L)) return(NULL)
        if (!is.null(.cl)) return(.cl)
        cl <- tryCatch(parallel::makePSOCKcluster(PAR_WORKERS),
                       error = function(e) {
                           cat(sprintf("  parallel disabled: %s\n",
                                       conditionMessage(e))); NULL })
        if (is.null(cl)) return(NULL)
        parallel::clusterEvalQ(cl, suppressPackageStartupMessages(
            library(Matrix)))
        .cl <<- cl
        cat(sprintf("  parallel: %d PSOCK worker(s)\n", PAR_WORKERS))
        cl
    }
    on.exit({
        if (!is.null(.cl)) tryCatch(parallel::stopCluster(.cl),
                                    error = function(e) invisible(NULL))
    }, add = TRUE)
    # Rebuild a task function in a clean environment, so the large
    # run_one_arm frame is not sent to the workers.
    p10hc_fn <- function(f, ...) {
        vals <- list(...)
        vals <- lapply(vals, function(v) {
            if (is.function(v)) environment(v) <- globalenv()
            v
        })
        environment(f) <- list2env(vals, parent = globalenv())
        f
    }

    # Falls back to lapply when there is no cluster.
    p10hc_map <- function(X, FUN, ...) {
        cl <- p10hc_cluster()
        if (is.null(cl)) lapply(X, FUN, ...)
        else tryCatch(parallel::parLapply(cl, X, FUN, ...),
                      error = function(e) {
                          cat(sprintf("  parallel call failed (%s); running sequentially\n",
                                      conditionMessage(e)))
                          lapply(X, FUN, ...)
                      })
    }

    # =============================================================================
    # AUCELL AND UCELL ON THE SPARSE COLUMNS
    # =============================================================================
    # Only the nonzero entries of each cell are ranked; a zero-valued
    # signature gene is placed among the zeros.
    # UCell in the same pass: average ties, cap at UCELL_MAXRANK, and a fixed
    # penalty beyond the cap. AUCell breaks ties at random.
    p10hc_auc_chunk <- function(Ec, si, k, max_area, ng, seed,
                                maxrank = UCELL_MAXRANK) {
        set.seed(seed)
        Ec <- as(Ec, "CsparseMatrix")
        p <- Ec@p; ii <- Ec@i; xx <- Ec@x
        nc <- ncol(Ec); nsig <- length(si)
        out <- numeric(nc); u_out <- numeric(nc); n_thin <- 0L
        for (j in seq_len(nc)) {
            lo <- p[j] + 1L; hi <- p[j + 1L]
            nnz <- hi - lo + 1L
            r <- numeric(0); ur <- numeric(0); n_hit <- 0L
            if (nnz > 0L) {
                rows <- ii[lo:hi] + 1L
                vals <- xx[lo:hi]
                # descending rank among the nonzeros, ties at random
                ord <- order(-vals, stats::runif(nnz))
                rk <- integer(nnz); rk[ord] <- seq_len(nnz)
                hit <- match(si, rows)
                r <- rk[hit[!is.na(hit)]]
                # descending rank among the nonzeros, ties averaged (UCell)
                asc <- vals[rev(ord)]
                uv <- vals[hit[!is.na(hit)]]
                n_hit <- length(uv)
                if (n_hit > 0L) {
                    n_le <- findInterval(uv, asc)
                    n_lt <- findInterval(uv, asc, left.open = TRUE)
                    ur <- (nnz - n_le) + ((n_le - n_lt) + 1) / 2
                }
            } else nnz <- 0L
            n_zero_sig <- nsig - n_hit
            if (n_zero_sig > 0L)
                ur <- c(ur, rep((nnz + (ng - nnz + 1) / 2), n_zero_sig))
            insig <- ur > maxrank
            u_out[j] <- if (all(insig)) 0 else {
                ur[insig] <- maxrank + 1
                1 - (sum(ur) - nsig * (nsig + 1) / 2) / (nsig * maxrank)
            }
            n_zero_auc <- nsig - length(r)
            if (n_zero_auc > 0L && nnz < k) {
                # Fewer than k nonzeros: a zero signature gene gets a random
                # rank in the zero block.
                n_zero <- ng - nnz
                r <- c(r, nnz + sample.int(n_zero, n_zero_auc))
                n_thin <- n_thin + 1L
            }
            r <- sort(r[r <= k])
            out[j] <- if (length(r) == 0L) 0 else
                (length(r) * k - sum(r) + length(r)) / max_area
        }
        list(auc = out, ucell = u_out, n_thin = n_thin)
    }

    p10hc_aucell <- function(E, sig, top_frac = AUCELL_TOP_FRAC,
                             block = AUCELL_CELL_BLOCK, seed = 42L) {
        sig <- intersect(sig, rownames(E))
        if (length(sig) == 0L) return(NULL)
        ng <- nrow(E); nc <- ncol(E)
        k  <- max(1L, floor(ng * top_frac))
        # AUC = sum(k - r_i + 1) over hits at r_i <= k; max_area puts the
        # signature at ranks 1..M.
        max_area <- k * length(sig) - length(sig) * (length(sig) - 1) / 2
        si <- match(sig, rownames(E))
        idx <- split(seq_len(nc), ceiling(seq_len(nc) / block))
        t0 <- proc.time()
        # Chunks are cut here so a worker never needs the whole of E.
        chunks <- lapply(seq_along(idx), function(b)
            list(m = E[, idx[[b]], drop = FALSE], seed = seed + b))
        worker <- p10hc_fn(function(ch)
            p10hc_auc_chunk(ch$m, si, k, max_area, ng, ch$seed, maxrank),
            p10hc_auc_chunk = p10hc_auc_chunk, si = si, k = k,
            max_area = max_area, ng = ng, maxrank = UCELL_MAXRANK)
        res <- p10hc_map(chunks, worker)
        rm(chunks)
        out <- numeric(nc); names(out) <- colnames(E)
        uout <- numeric(nc); names(uout) <- colnames(E)
        for (b in seq_along(idx)) {
            out[idx[[b]]]  <- res[[b]]$auc
            uout[idx[[b]]] <- res[[b]]$ucell
        }
        n_thin <- sum(vapply(res, function(r) r$n_thin, integer(1)))
        cat(sprintf("    AUCell + UCell: %d cell(s) in %d chunk(s), %.1f sec\n",
                    nc, length(idx), (proc.time() - t0)[["elapsed"]]))
        # Both are bounded in [0, 1]; values outside are clamped and counted.
        n_out <- sum(out < 0 | out > 1, na.rm = TRUE) +
                 sum(uout < 0 | uout > 1, na.rm = TRUE)
        list(auc = pmin(pmax(out, 0), 1), ucell = pmin(pmax(uout, 0), 1),
             n_clamped = n_out, k = k, maxrank = UCELL_MAXRANK,
             n_sig = length(sig), n_thin = n_thin, n_cells = nc)
    }

    # =============================================================================
    # ADDMODULESCORE
    # =============================================================================
    # Signature mean minus the mean of control genes drawn from the same
    # expression bins. Bins use rank(ties.method = "first"), so no RNG.
    p10hc_addmodulescore <- function(E, sig, nbin = AMS_NBIN, ctrl = AMS_CTRL,
                                     seed = AMS_SEED) {
        sig <- intersect(sig, rownames(E))
        if (length(sig) == 0L || ncol(E) == 0L || nrow(E) < 2L) return(NULL)
        avg <- Matrix::rowMeans(E)
        nbin <- max(1L, min(as.integer(nbin), nrow(E)))
        bin <- cut(rank(avg, ties.method = "first"), breaks = nbin,
                   labels = FALSE)
        names(bin) <- rownames(E)
        set.seed(seed)
        cu <- character(0); n_short <- 0L
        for (g in sig) {
            pool <- names(bin)[bin == bin[[g]]]
            n <- min(as.integer(ctrl), length(pool))
            if (n < ctrl) n_short <- n_short + 1L
            cu <- c(cu, pool[sample.int(length(pool), n)])
        }
        cu <- unique(cu)
        if (length(cu) == 0L) return(NULL)
        sc <- Matrix::colMeans(E[sig, , drop = FALSE]) -
              Matrix::colMeans(E[cu, , drop = FALSE])
        names(sc) <- colnames(E)
        list(score = as.numeric(sc), n_sig = length(sig), n_ctrl = length(cu),
             n_bins_used = length(unique(bin[sig])), n_short = n_short,
             n_ctrl_in_sig = length(intersect(cu, sig)),
             n_ctrl_mito = sum(grepl("^MT-", cu)),
             mean_sig = mean(avg[sig]), mean_ctrl = mean(avg[cu]),
             cells = colnames(E))
    }

    # =============================================================================
    # STEP 1: Part 10h's results object
    # =============================================================================
    # prop_df is the donor frame: pct_post plus the clinical columns.
    cat("STEP 1 -- Part 10h results object\n")
    if (!file.exists(P10H_RDS))
        stop("Part 10h's results RDS was not found:\n    ", P10H_RDS,
             call. = FALSE)
    R10H <- readRDS(P10H_RDS)
    DECL   <- R10H[["declarations"]]
    DSETS  <- R10H[["datasets"]]
    SHEETS <- R10H[["sheets"]]
    M10H   <- R10H[["meta"]]
    if (is.null(DECL) || is.null(DSETS))
        stop("Part 10h's RDS carries no `declarations` or no `datasets`.",
             call. = FALSE)

    COV_SETS      <- DECL[["COV_SETS"]]
    OUTCOME_ORDER <- DECL[["MV_OUTCOME_ORDER"]]
    CONTROL_OUT   <- DECL[["MV_CONTROL_OUTCOMES"]]
    if (is.null(COV_SETS) || is.null(OUTCOME_ORDER))
        stop("Part 10h's declarations carry no `COV_SETS` or no ",
             "`MV_OUTCOME_ORDER`.",
             call. = FALSE)
    if (!PRIMARY_BLOCK %in% names(COV_SETS))
        stop("The primary block '", PRIMARY_BLOCK, "' is not among Part 10h's ",
             "covariate sets (", paste(names(COV_SETS), collapse = ", "), ").",
             call. = FALSE)

    # == Donor frame: prop_df per dataset =====================================
    DONOR <- list()
    for (ds in names(DSETS)) {
        pd <- DSETS[[ds]][["prop_df"]]
        if (!is.data.frame(pd) || nrow(pd) == 0L) next
        if (is.null(pd[["donor_id"]]) || is.null(pd[["pct_post"]])) next
        pd$.donor <- as.character(pd[["donor_id"]])
        pd$dataset <- ds
        pd$module <- suppressWarnings(as.numeric(pd[["pct_post"]]))
        DONOR[[ds]] <- pd
    }
    if (length(DONOR) == 0L)
        stop("No dataset in Part 10h's `datasets` carries a `prop_df` with both ",
             "`donor_id` and `pct_post`.\n  Frames present: ",
             paste(unique(unlist(lapply(DSETS, names))), collapse = ", "),
             call. = FALSE)
    # rbind on the common columns only.
    .common <- Reduce(intersect, lapply(DONOR, names))
    .dropped <- setdiff(unique(unlist(lapply(DONOR, names))), .common)
    DONOR_DF <- do.call(rbind, lapply(DONOR, function(d)
        d[, .common, drop = FALSE]))
    rownames(DONOR_DF) <- NULL
    if (length(.dropped) > 0L)
        cat(sprintf("  NOTE: %d column(s) present in one prop_df and not the other, dropped: %s\n",
                    length(.dropped), paste(.dropped, collapse = ", ")))
    gr10hc("GR-10HC-02", "PASS",
           sprintf(paste("module predictor = prop_df$pct_post keyed on donor_id,",
                         "over %s | %d donor-dataset row(s), %d column(s)"),
                   paste(sprintf("%s %d", names(DONOR),
                                 vapply(DONOR, nrow, integer(1))),
                         collapse = " and "),
                   nrow(DONOR_DF), ncol(DONOR_DF)))

    # == Key to column, label to column =======================================
    # Fallbacks; Part 10h's COV_MAP and CLINICAL_PARAMS are preferred.
    BUILTIN_COV_COL <- c(age = "Age_num", apoe = "APOE_e4_dose", sex = "Sex_num",
                         pmi = "PMI_num", rin = "RIN_num", ph = "BrainPH_num",
                         weight = "BrainWeight_num", mito = "MitoFrac_num",
                         numi = "nUMI_num", ngene = "nGene_num",
                         intronexon = "IntronExon_num", doublet = "Doublet_num",
                         satur = "Satur_num", study = "Study_num",
                         late = "LATE_num", lewyfull = "LewyFull_num",
                         athero = "Athero_num", arteriolo = "Arteriolo_num",
                         microinf = "Microinf_num", eduyears = "EduYears_num",
                         mmseint = "MMSEint_num", casiint = "CASIint_num")
    BUILTIN_OUT_COL <- c(
        "Cognitive status" = "CogStatus_num", "MMSE score" = "MMSE_num",
        "CASI score" = "CASI_num", "CERAD score" = "CERAD_num",
        "Braak stage" = "Braak_num", "Thal phase" = "Thal_num",
        "ADNC rank" = "ADNC_rank", "CAA score" = "CAA_num",
        "APOE e4 dose" = "APOE_e4_dose", "LATE stage" = "LATE_num",
        "Lewy body stage" = "LewyFull_num", "Atherosclerosis" = "Athero_num",
        "Arteriolosclerosis" = "Arteriolo_num", "Microinfarcts" = "Microinf_num",
        "Age at death" = "Age_num", "Sex" = "Sex_num", "Brain pH" = "BrainPH_num",
        "PMI" = "PMI_num", "RIN" = "RIN_num",
        "Fresh brain weight" = "BrainWeight_num",
        "Years of education" = "EduYears_num",
        "Months from last MMSE" = "MMSEint_num",
        "Mito fraction (donor mean)" = "MitoFrac_num",
        "UMIs per cell (donor mean)" = "nUMI_num",
        "Genes detected (donor mean)" = "nGene_num",
        "Doublet score (donor mean)" = "Doublet_num",
        "Intronic:exonic ratio (donor mean)" = "IntronExon_num",
        "Sequencing saturation (donor mean)" = "Satur_num")

    # A usable map: named strings whose values are donor-frame columns.
    as_map <- function(x) {
        if (is.null(x)) return(NULL)
        if (is.list(x) && all(vapply(x, function(v)
            is.character(v) && length(v) == 1L, logical(1))))
            x <- unlist(x)
        if (!is.character(x) || is.null(names(x))) return(NULL)
        x <- x[nzchar(names(x))]
        if (sum(x %in% names(DONOR_DF)) < 2L) return(NULL)
        x
    }
    .cov_map <- as_map(DECL[["COV_MAP"]])
    .out_map <- as_map(DECL[["CLINICAL_PARAMS"]])
    COV_ROUTE <- if (is.null(.cov_map)) "built-in" else "Part 10h COV_MAP"
    OUT_ROUTE <- if (is.null(.out_map)) "built-in" else "Part 10h CLINICAL_PARAMS"
    if (is.null(.cov_map)) .cov_map <- BUILTIN_COV_COL
    if (is.null(.out_map)) .out_map <- BUILTIN_OUT_COL

    resolve <- function(keys, map) {
        out <- unname(map[keys])
        out[!is.na(out) & out %in% names(DONOR_DF)]
    }
    COVS <- lapply(COV_SETS, resolve, map = .cov_map)
    .cov_lost <- lapply(names(COV_SETS), function(b)
        COV_SETS[[b]][is.na(.cov_map[COV_SETS[[b]]]) |
                      !(.cov_map[COV_SETS[[b]]] %in% names(DONOR_DF))])
    names(.cov_lost) <- names(COV_SETS)
    .cov_lost <- .cov_lost[vapply(.cov_lost, length, integer(1)) > 0L]

    OUT_MAP <- .out_map[OUTCOME_ORDER]
    names(OUT_MAP) <- OUTCOME_ORDER
    OUT_MAP <- OUT_MAP[!is.na(OUT_MAP) & OUT_MAP %in% names(DONOR_DF)]
    OUT_MAP <- OUT_MAP[vapply(OUT_MAP, function(c1)
        sum(is.finite(suppressWarnings(as.numeric(DONOR_DF[[c1]])))) >= MIN_DONORS,
        logical(1))]
    OUT_LOST <- setdiff(OUTCOME_ORDER, names(OUT_MAP))

    # Keep Part 10h's modelled outcome set, so the BH family matches.
    MVSH <- SHEETS[["Multivariate"]]
    .p10h_out <- if (is.data.frame(MVSH) && !is.null(MVSH[["outcome"]]))
        unique(as.character(MVSH$outcome)) else character(0)
    .out_extra <- character(0)
    if (length(.p10h_out) > 0L) {
        .keep <- intersect(names(OUT_MAP), .p10h_out)
        if (length(.keep) >= 4L) {
            .out_extra <- setdiff(names(OUT_MAP), .p10h_out)
            OUT_MAP <- OUT_MAP[.keep]
        }
    }

    # Outcome classes.
    OUT_CLASS <- ifelse(names(OUT_MAP) %in% CONTROL_OUT, "control",
                 ifelse(grepl(TECH_PATTERN, names(OUT_MAP)), "technical",
                        "clinical"))
    names(OUT_CLASS) <- names(OUT_MAP)
    CLIN_LABS <- names(OUT_CLASS)[OUT_CLASS == "clinical"]
    if (length(CLIN_LABS) < 4L)
        stop("Only ", length(CLIN_LABS), " clinical outcome(s) after removing ",
             "the controls and the technical donor means.", call. = FALSE)

    if (length(COVS[[PRIMARY_BLOCK]]) == 0L)
        stop("Not one covariate of the primary block '", PRIMARY_BLOCK,
             "' resolved to a column of the donor frame. Keys: ",
             paste(COV_SETS[[PRIMARY_BLOCK]], collapse = ", "),
             "\n  Columns available: ",
             paste(utils::head(names(DONOR_DF), 80), collapse = ", "),
             call. = FALSE)
    if (length(OUT_MAP) < 4L)
        stop("Only ", length(OUT_MAP), " outcome(s) resolved to a column of the ",
             "donor frame; at least four are needed for a comparison. Unresolved: ",
             paste(OUT_LOST, collapse = ", "), call. = FALSE)

    gr10hc("GR-10HC-01",
           if (length(OUT_LOST) > 0L || length(.cov_lost) > 0L) "WARN" else "PASS",
           sprintf(paste("%s (schema %s, written %s) | %d donor-dataset row(s), %d",
                         "outcome(s) of %d declared, %d covariate block(s): %s |",
                         "covariate map via %s, outcome map via %s%s%s%s%s"),
                   basename(P10H_RDS),
                   if (is.null(M10H[["schema"]])) "not recorded" else M10H[["schema"]],
                   if (is.null(M10H[["written"]])) "not recorded" else M10H[["written"]],
                   nrow(DONOR_DF), length(OUT_MAP), length(OUTCOME_ORDER),
                   length(COVS),
                   paste(sprintf("%s(%d)", names(COVS),
                                 vapply(COVS, length, integer(1))), collapse = " "),
                   COV_ROUTE, OUT_ROUTE,
                   sprintf(" | classes: %s",
                           paste(sprintf("%s %d", names(table(OUT_CLASS)),
                                         as.integer(table(OUT_CLASS))),
                                 collapse = ", ")),
                   if (length(.out_extra) > 0L)
                       paste0(" | not modelled by Part 10h, dropped: ",
                              paste(.out_extra, collapse = ", "))
                   else "",
                   if (length(OUT_LOST) > 0L)
                       paste0(" | not resolved: ",
                              paste(OUT_LOST, collapse = ", ")) else "",
                   if (length(.cov_lost) > 0L)
                       paste0(" | covariate keys not resolved: ",
                              paste(sprintf("%s[%s]", names(.cov_lost),
                                            vapply(.cov_lost, paste, character(1),
                                                   collapse = ",")),
                                    collapse = " ")) else ""))

    # == Part 10h IP, for the changepoint comparison ==========================
    # psi is on the GOI axis, so IP_x_reported is the comparable quantity.
    P10H_IP <- list(); .ip_src <- NA_character_; .ip_col <- NA_character_
    for (.sh in names(SHEETS)) {
        fr <- SHEETS[[.sh]]
        if (!is.data.frame(fr) || is.null(fr[["dataset"]])) next
        ipc <- intersect(c("IP_x_reported", "IP_x", "ip_x"), names(fr))
        if (length(ipc) == 0L) next
        for (dsx in unique(as.character(fr$dataset)))
            P10H_IP[[dsx]] <- suppressWarnings(as.numeric(
                fr[[ipc[1]]][as.character(fr$dataset) == dsx][1]))
        .ip_src <- .sh; .ip_col <- ipc[1]
        break
    }

    # =============================================================================
    # STEP 2: the five baselines
    # =============================================================================
    cat("\nSTEP 2 -- external baselines\n")
    BASE <- list(); .b_report <- list(); .seg_report <- list()
    .u_report <- list(); .ams_report <- list()
    # Per-cell scores per dataset: barcode, auc, ucell, ams, past_psi
    .percell <- list()
    for (ds in DATASETS) {
        cat(sprintf("\n  %s\n", ds))
        obs_p <- file.path(ds_dir(ds), "metadata", "obs_dataset.csv")
        if (!file.exists(obs_p)) {
            cat("    SKIP: obs_dataset.csv not found\n"); next
        }
        obs <- read.csv(obs_p, stringsAsFactors = FALSE, check.names = FALSE)
        if (!("barcode" %in% names(obs))) obs$barcode <- obs[[1]]
        if (!(STYPE_COL %in% names(obs)) || !(DONOR_COL %in% names(obs))) {
            cat("    SKIP: obs_dataset.csv lacks the supertype or donor column\n")
            next
        }
        obs$.st <- as.character(obs[[STYPE_COL]])
        obs$.dn <- trimws(as.character(obs[[DONOR_COL]]))
        obs$.dn[obs$.dn == "" | tolower(obs$.dn) %in% c("na", "nan")] <- NA
        obs <- obs[!is.na(obs$.dn) & !is.na(obs$.st), , drop = FALSE]

        # == B1: cluster proportion ===========================================
        micro <- obs[obs$.st %in% MICROGLIAL_SUPERTYPES, , drop = FALSE]
        if (nrow(micro) == 0L) {
            cat("    SKIP: no cell carries a declared microglial supertype\n")
            next
        }
        den <- table(micro$.dn)
        num <- table(micro$.dn[micro$.st == PRIMARY_STYPE])
        b1 <- data.frame(.donor = names(den),
                         B1 = 100 * as.numeric(num[names(den)]) / as.numeric(den),
                         n_micro = as.integer(den), stringsAsFactors = FALSE)
        b1$B1[is.na(b1$B1)] <- 0
        cat(sprintf("    B1: %s is %.2f%% of %s microglial cell(s) over %d donor(s)\n",
                    PRIMARY_STYPE,
                    100 * sum(micro$.st == PRIMARY_STYPE) / nrow(micro),
                    format(nrow(micro), big.mark = ","), nrow(b1)))

        # == One read of the dataset matrix ===================================
        st_bc <- obs$barcode[obs$.st == PRIMARY_STYPE]
        rd <- file.path(SHARED_DIR, basename(ds_dir(ds)), "counts", "raw_dataset")
        mp <- file.path(rd, "matrix.mtx"); bp <- file.path(rd, "matrix_barcodes.tsv")
        fp <- file.path(rd, "matrix_features.tsv")
        E <- NULL
        if (all(file.exists(c(mp, bp, fp)))) {
            bc <- readLines(bp, warn = FALSE)
            if (length(bc) > RAW_MAX_CELLS)
                cat(sprintf("    B2/B3 SKIP: raw_dataset holds %s cell(s), above RAW_MAX_CELLS\n",
                            format(length(bc), big.mark = ",")))
            else {
                ft <- read.table(fp, sep = "\t", header = FALSE,
                                 stringsAsFactors = FALSE)
                X <- as(readMM(mp), "CsparseMatrix")
                if (nrow(X) != length(bc) || ncol(X) != nrow(ft))
                    cat("    B2/B3 SKIP: raw_dataset dimensions disagree with its index files\n")
                else {
                    rownames(X) <- bc
                    colnames(X) <- make.unique(as.character(ft[, 2]))
                    keep <- intersect(st_bc, bc)
                    X <- X[keep, , drop = FALSE]
                    E <- Matrix::t(p10hc_cp10k(X))
                    rm(X); gc(verbose = FALSE)
                    cat(sprintf("    matrix: %d gene(s) x %d %s cell(s)\n",
                                nrow(E), ncol(E), PRIMARY_STYPE))
                    # Zero columns: barcode mismatch with obs_dataset.csv.
                    if (ncol(E) == 0L) {
                        cat("    B2/B3/B4/B5 SKIP: no supertype barcode of ",
                            "obs_dataset.csv is present in the raw matrix\n",
                            sep = "")
                        E <- NULL
                    }
                }
            }
        } else cat(sprintf("    B2/B3 SKIP: raw_dataset not found under %s\n", rd))

        b2 <- NULL; b3 <- NULL; b4 <- NULL; b5 <- NULL; seg_rep <- NULL
        if (!is.null(E)) {
            dn_cell <- obs$.dn[match(colnames(E), obs$barcode)]
            # Per-cell scores are kept for the leave-one-region-out stage.
            .pc <- data.frame(barcode = colnames(E), donor = dn_cell,
                              auc = NA_real_, ucell = NA_real_,
                              ams = NA_real_, past_psi = NA,
                              stringsAsFactors = FALSE)

            # == B2: AUCell ===================================================
            sig <- unique(c(GOI, DFG))
            au <- p10hc_aucell(E, sig)
            if (is.null(au)) cat("    B2 SKIP: no signature gene in the matrix\n")
            else {
                .pc$auc <- as.numeric(au$auc)
                b2 <- stats::aggregate(list(B2 = au$auc),
                                       by = list(.donor = dn_cell), FUN = mean)
                cat(sprintf("    B2: AUCell over %d signature gene(s), top %d rank(s); donor means %.4f-%.4f\n",
                            au$n_sig, au$k, min(b2$B2), max(b2$B2)))
                # Cross-check against the AUCell package when installed.
                au$check_rho <- NA_real_; au$check_n <- 0L
                au$check_note <- if (.HAVE_AUCELL) "" else
                    "the AUCell package is not installed"
                if (.HAVE_AUCELL) {
                    set.seed(SEG_SEED)
                    ci <- if (ncol(E) > AUCELL_CHECK_N)
                        sample.int(ncol(E), AUCELL_CHECK_N) else seq_len(ncol(E))
                    ref <- tryCatch({
                        rk <- AUCell::AUCell_buildRankings(
                            as.matrix(E[, ci, drop = FALSE]), plotStats = FALSE,
                            verbose = FALSE)
                        a2 <- AUCell::AUCell_calcAUC(
                            list(sig = intersect(sig, rownames(E))), rk,
                            aucMaxRank = au$k, verbose = FALSE)
                        as.numeric(AUCell::getAUC(a2)["sig", ])
                    }, error = function(e) {
                        au$check_note <<- conditionMessage(e); NULL })
                    if (!is.null(ref) && length(ref) == length(ci)) {
                        au$check_rho <- suppressWarnings(stats::cor(
                            au$auc[ci], ref, method = "spearman"))
                        au$check_n <- length(ci)
                        cat(sprintf("    B2 cross-check vs the AUCell package on %d cell(s): Spearman %.4f\n",
                                    au$check_n, au$check_rho))
                    }
                }
                .b_report[[paste0(ds, "_B2")]] <- au

                # == B4: UCell ================================================
                # Computed in the same scan as B2.
                .pc$ucell <- as.numeric(au$ucell)
                b4 <- stats::aggregate(list(B4 = au$ucell),
                                       by = list(.donor = dn_cell), FUN = mean)
                cat(sprintf("    B4: UCell over %d signature gene(s), maxRank %d; donor means %.4f-%.4f\n",
                            au$n_sig, au$maxrank, min(b4$B4), max(b4$B4)))
                uc <- list(check_rho = NA_real_, check_n = 0L,
                           maxrank = au$maxrank, n_sig = au$n_sig,
                           check_note = if (.HAVE_UCELL) ""
                                        else "the UCell package is absent")
                if (.HAVE_UCELL) {
                    set.seed(SEG_SEED)
                    ci <- if (ncol(E) > AUCELL_CHECK_N)
                        sample.int(ncol(E), AUCELL_CHECK_N) else seq_len(ncol(E))
                    uref <- tryCatch({
                        um <- UCell::ScoreSignatures_UCell(
                            E[, ci, drop = FALSE],
                            features = list(sig = intersect(sig, rownames(E))),
                            maxRank = au$maxrank, ncores = 1)
                        as.numeric(um[, 1])
                    }, error = function(e) {
                        uc$check_note <<- conditionMessage(e); NULL })
                    if (!is.null(uref) && length(uref) == length(ci)) {
                        uc$check_rho <- suppressWarnings(stats::cor(
                            au$ucell[ci], uref, method = "spearman"))
                        uc$check_n <- length(ci)
                        cat(sprintf("    B4 cross-check vs the UCell package on %d cell(s): Spearman %.4f\n",
                                    uc$check_n, uc$check_rho))
                    }
                }
                .u_report[[paste0(ds, "_B4")]] <- uc

                # == B5: AddModuleScore =======================================
                am <- p10hc_addmodulescore(E, sig)
                if (is.null(am))
                    cat("    B5 SKIP: no signature gene in the matrix, or",
                        " fewer than two genes to bin\n", sep = "")
                else {
                    .pc$ams <- am$score
                    b5 <- stats::aggregate(list(B5 = am$score),
                                           by = list(.donor = dn_cell),
                                           FUN = mean)
                    cat(sprintf("    B5: AddModuleScore, %d control gene(s) from %d bin(s); donor means %.4f-%.4f\n",
                                am$n_ctrl, am$n_bins_used, min(b5$B5),
                                max(b5$B5)))
                    cat(sprintf("        control set: %d of %d gene(s) match '^MT-', mean expression %.4f against the signature's %.4f\n",
                                am$n_ctrl_mito, am$n_ctrl, am$mean_ctrl,
                                am$mean_sig))
                    am$ds <- ds
                    .ams_report[[ds]] <- am
                }
            }

            # == B3: segmented ================================================
            gi <- match(GOI, rownames(E))
            dfi <- match(intersect(DFG, rownames(E)), rownames(E))
            if (is.na(gi) || length(dfi) == 0L)
                cat("    B3 SKIP: the GOI or every DFG is absent from the matrix\n")
            else if (!.HAVE_SEGMENTED)
                cat("    B3 SKIP: the `segmented` package is not installed\n")
            else {
                gx <- as.numeric(E[gi, ])
                cy <- if (length(dfi) == 1L) as.numeric(E[dfi, ])
                      else Matrix::colMeans(E[dfi, , drop = FALSE])
                set.seed(SEG_SEED)
                idx <- if (length(gx) > SEG_MAX_CELLS)
                    sample.int(length(gx), SEG_MAX_CELLS) else seq_along(gx)
                # segmented() needs a named data frame for seg.Z.
                segD <- data.frame(x = gx[idx], y = cy[idx])
                segD <- segD[is.finite(segD$x) & is.finite(segD$y), ,
                             drop = FALSE]
                .seg_err <- ""
                sg <- NULL
                if (nrow(segD) >= 50L) {
                    fit0 <- stats::lm(y ~ x, data = segD)
                    # Start at the median of the GOI axis.
                    sg <- tryCatch(
                        segmented::segmented(fit0, seg.Z = ~x,
                            psi = stats::median(segD$x, na.rm = TRUE)),
                        error = function(e) { .seg_err <<- conditionMessage(e)
                                              NULL })
                    if (is.null(sg))
                        sg <- tryCatch(segmented::segmented(fit0, seg.Z = ~x),
                                       error = function(e) {
                                           .seg_err <<- conditionMessage(e)
                                           NULL })
                } else .seg_err <- sprintf("only %d finite (x, y) pair(s)",
                                           nrow(segD))
                psi <- if (is.null(sg) || is.null(sg$psi)) NA_real_
                       else as.numeric(sg$psi[1, "Est."])
                psi_se <- if (is.null(sg) || is.null(sg$psi)) NA_real_
                          else as.numeric(sg$psi[1, "St.Err"])
                if (!is.finite(psi))
                    cat(sprintf("    B3 SKIP: segmented() returned no breakpoint%s\n",
                                if (nzchar(.seg_err))
                                    paste0(" -- ", .seg_err) else ""))
                else {
                    past <- gx > psi
                    .pc$past_psi <- past
                    b3 <- stats::aggregate(list(B3 = past),
                                           by = list(.donor = dn_cell),
                                           FUN = function(v) 100 * mean(v))
                    seg_rep <- list(ds = ds, psi = psi, se = psi_se,
                                    n_used = length(idx), n_total = length(gx),
                                    pct_past = 100 * mean(past))
                    cat(sprintf("    B3: breakpoint on the %s axis at %.4f (SE %.4f), %.1f%% of cells past it\n",
                                GOI, psi, psi_se, 100 * mean(past)))
                }
            }
            .percell[[ds]] <- .pc
            rm(E); gc(verbose = FALSE)
        }
        .seg_report[[ds]] <- seg_rep

        d <- b1
        d$dataset <- ds
        if (!is.null(b2)) d <- merge(d, b2, by = ".donor", all.x = TRUE)
        else d$B2 <- NA_real_
        if (!is.null(b3)) d <- merge(d, b3, by = ".donor", all.x = TRUE)
        else d$B3 <- NA_real_
        if (!is.null(b4)) d <- merge(d, b4, by = ".donor", all.x = TRUE)
        else d$B4 <- NA_real_
        if (!is.null(b5)) d <- merge(d, b5, by = ".donor", all.x = TRUE)
        else d$B5 <- NA_real_
        BASE[[ds]] <- d
    }
    BASE_DF <- if (length(BASE) > 0) do.call(rbind, BASE) else NULL
    if (is.null(BASE_DF)) stop("No dataset produced a baseline.", call. = FALSE)
    rownames(BASE_DF) <- NULL

    .avail <- vapply(c("B1", "B2", "B3", "B4", "B5"), function(k)
        sum(is.finite(BASE_DF[[k]])), integer(1))
    gr10hc("GR-10HC-03",
           if (.avail[["B1"]] == 0L) "FAIL"
           else if (any(.avail == 0L)) "WARN" else "PASS",
           sprintf(paste("baselines built: B1 cluster proportion %d row(s), B2",
                         "AUCell %d, B3 changepoint %d, B4 UCell %d,",
                         "B5 AddModuleScore %d, over %d dataset(s).%s"),
                   .avail[["B1"]], .avail[["B2"]], .avail[["B3"]],
                   .avail[["B4"]], .avail[["B5"]], length(BASE),
                   if (any(.avail == 0L))
                       paste(" absent:",
                             paste(names(.avail)[.avail == 0L], collapse = ", "))
                   else ""))

    .seg_ok <- Filter(Negate(is.null), .seg_report)
    gr10hc("GR-10HC-04",
           if (length(.seg_ok) == 0L) "WARN" else "PASS",
           if (length(.seg_ok) == 0L)
               paste("no changepoint was fitted")
           else sprintf(paste("changepoint on the %s axis: %s"),
                        GOI,
                        paste(vapply(.seg_ok, function(r)
                            sprintf("%s psi=%.4f (SE %.4f, %s of %s cell(s) used), %.1f%% past; Part 10h IP %s",
                                    r$ds, r$psi, r$se,
                                    format(r$n_used, big.mark = ","),
                                    format(r$n_total, big.mark = ","),
                                    r$pct_past,
                                    if (is.null(P10H_IP[[r$ds]]) ||
                                        !is.finite(P10H_IP[[r$ds]]))
                                        "not found in the RDS sheets"
                                    else sprintf("%.4f (%s$%s), difference %.4f",
                                                 P10H_IP[[r$ds]], .ip_src,
                                                 .ip_col,
                                                 r$psi - P10H_IP[[r$ds]])),
                            character(1)), collapse = " | ")))

    # Built-in AUCell against the AUCell package.
    .chk <- Filter(function(r) is.finite(r$check_rho), .b_report)
    if (length(.chk) > 0)
        gr10hc("GR-10HC-05",
               if (min(vapply(.chk, function(r) r$check_rho, numeric(1))) < 0.99)
                   "WARN" else "PASS",
               sprintf(paste("built-in AUCell against the AUCell package: %s"),
                       paste(vapply(names(.chk), function(k)
                           sprintf("%s rho %.4f on %d cell(s)", k,
                                   .chk[[k]]$check_rho, .chk[[k]]$check_n),
                           character(1)), collapse = " | ")))
    else if (.HAVE_AUCELL && length(.b_report) > 0)
        gr10hc("GR-10HC-05", "WARN",
               sprintf(paste("the AUCell package is installed but the",
                             "cross-check did not complete (%s)"),
                       paste(unique(vapply(.b_report, function(r)
                           if (nzchar(r$check_note)) r$check_note else "no reason recorded",
                           character(1))), collapse = "; ")))
    else if (length(.b_report) == 0L)
        gr10hc("GR-10HC-05", "WARN",
               sprintf(paste("B2 was not computed on any dataset; the AUCell",
                             "package is %s"),
                       if (.HAVE_AUCELL) "installed" else "not installed"))
    else
        gr10hc("GR-10HC-05", "WARN",
               paste("the AUCell package is not installed; B2 is the built-in",
                     "statistic, top", sprintf("%.0f%%", 100 * AUCELL_TOP_FRAC),
                     "of ranks, not cross-checked"))

    # Built-in UCell against the UCell package.
    .uchk <- Filter(function(r) is.finite(r$check_rho), .u_report)
    if (length(.uchk) > 0)
        gr10hc("GR-10HC-14",
               if (min(vapply(.uchk, function(r) r$check_rho,
                              numeric(1))) < 0.99)
                   "WARN" else "PASS",
               sprintf(paste("built-in UCell against the UCell package: %s"),
                       paste(vapply(names(.uchk), function(k)
                           sprintf("%s rho %.4f on %d cell(s)", k,
                                   .uchk[[k]]$check_rho, .uchk[[k]]$check_n),
                           character(1)), collapse = " | ")))
    else if (length(.u_report) > 0)
        gr10hc("GR-10HC-14", "WARN",
               sprintf(paste("B4 rests on the built-in UCell statistic alone",
                             "(%s), maxRank %d, not cross-checked"),
                       paste(unique(vapply(.u_report, function(r)
                           if (nzchar(r$check_note)) r$check_note
                           else "no reason recorded", character(1))),
                           collapse = "; "),
                       UCELL_MAXRANK))
    else
        gr10hc("GR-10HC-14", "WARN",
               paste("B4 was not computed on any dataset"))

    # AddModuleScore control set: mitochondrial share and short bins.
    if (length(.ams_report) > 0) {
        .mito_frac <- max(vapply(.ams_report, function(r)
            r$n_ctrl_mito / max(1L, r$n_ctrl), numeric(1)))
        .short <- sum(vapply(.ams_report, function(r) r$n_short, integer(1)))
        gr10hc("GR-10HC-16",
               if (.mito_frac >= 0.25 || .short > 0L) "WARN" else "PASS",
               sprintf(paste("B5 control set: %s.%s%s"),
                       paste(vapply(names(.ams_report), function(k) {
                           r <- .ams_report[[k]]
                           sprintf(paste("%s %d control gene(s) from %d",
                                         "bin(s), %d matching '^MT-', %d",
                                         "already in the signature; mean",
                                         "expression %.4f against the",
                                         "signature's %.4f"),
                                   k, r$n_ctrl, r$n_bins_used, r$n_ctrl_mito,
                                   r$n_ctrl_in_sig, r$mean_ctrl, r$mean_sig)
                       }, character(1)), collapse = " | "),
                       if (.mito_frac >= 0.25)
                           sprintf(paste(" At least %.0f%% of the control set",
                                         "is mitochondrial."),
                                   100 * .mito_frac)
                       else "",
                       if (.short > 0L)
                           sprintf(paste(" %d signature gene(s) sit in a bin",
                                         "smaller than AMS_CTRL = %d."),
                                   .short, AMS_CTRL)
                       else ""))
    } else
        gr10hc("GR-10HC-16", "WARN",
               paste("B5 was not computed on any dataset"))

    # =============================================================================
    # STEP 3: one donor frame
    # =============================================================================
    # The five baselines are joined to DONOR_DF on donor and dataset.
    cat("\nSTEP 3 -- donor frame\n")
    # Measure the donor-key overlap; try `Donor ID` if low; stop on zero.
    .ov <- function(a, b) length(intersect(unique(a), unique(b)))
    .key_used <- "prop_df$donor_id"
    .n_ov <- .ov(DONOR_DF$.donor, BASE_DF$.donor)
    if (.n_ov < 0.5 * length(unique(BASE_DF$.donor)) &&
        !is.null(DONOR_DF[["Donor ID"]])) {
        .alt <- .ov(as.character(DONOR_DF[["Donor ID"]]), BASE_DF$.donor)
        if (.alt > .n_ov) {
            DONOR_DF$.donor <- as.character(DONOR_DF[["Donor ID"]])
            .key_used <- "prop_df$`Donor ID`"; .n_ov <- .alt
        }
    }
    if (.n_ov == 0L)
        stop("No donor identifier in Part 10h's prop_df matches the donor ",
             "column of obs_dataset.csv, so every baseline would join to NA.",
             "\n  prop_df keys look like: ",
             paste(utils::head(unique(DONOR_DF$.donor), 3), collapse = ", "),
             "\n  obs_dataset keys look like: ",
             paste(utils::head(unique(BASE_DF$.donor), 3), collapse = ", "),
             call. = FALSE)
    DF <- merge(DONOR_DF, BASE_DF, by = c(".donor", "dataset"), all.x = TRUE)
    cat(sprintf("  %d donor-dataset row(s) over %d donor(s), %d column(s) | join on %s, %d donor(s) matched\n",
                nrow(DF), length(unique(DF$.donor)), ncol(DF), .key_used, .n_ov))
    gr10hc("GR-10HC-11",
           if (.n_ov < 0.9 * length(unique(BASE_DF$.donor))) "WARN" else "PASS",
           sprintf(paste("donor join: %d of %d obs_dataset donor(s) matched a",
                         "prop_df row, on %s"),
                   .n_ov, length(unique(BASE_DF$.donor)), .key_used))
    PRED <- c("module", "B1", "B2", "B3", "B4", "B5")
    PRED <- PRED[vapply(PRED, function(k)
        !is.null(DF[[k]]) && sum(is.finite(DF[[k]])) >= MIN_DONORS, logical(1))]
    cat(sprintf("  predictors with at least %d donor(s): %s\n", MIN_DONORS,
                paste(vapply(PRED, function(k) PRED_KEYS[[k]], character(1)),
                      collapse = ", ")))
    if (!"module" %in% PRED)
        stop("The module predictor has fewer than ", MIN_DONORS,
             " finite donor value(s).",
             call. = FALSE)

    # Spearman correlations among the predictors, per dataset.
    cor_rows <- list()
    for (ds in unique(DF$dataset)) {
        d <- DF[DF$dataset == ds, , drop = FALSE]
        for (i in seq_along(PRED)) for (j in seq_along(PRED)) if (i < j) {
            x <- suppressWarnings(as.numeric(d[[PRED[i]]]))
            y <- suppressWarnings(as.numeric(d[[PRED[j]]]))
            ok <- is.finite(x) & is.finite(y)
            cor_rows[[length(cor_rows) + 1L]] <- data.frame(
                dataset = ds, a = PRED_KEYS[[PRED[i]]], b = PRED_KEYS[[PRED[j]]],
                n = sum(ok),
                spearman = if (sum(ok) >= MIN_DONORS)
                    round(suppressWarnings(stats::cor(x[ok], y[ok],
                                                      method = "spearman")), 3)
                else NA_real_, stringsAsFactors = FALSE)
        }
    }
    COR_DF <- if (length(cor_rows) > 0) do.call(rbind, cor_rows) else NULL
    if (!is.null(COR_DF)) { rownames(COR_DF) <- NULL; print(COR_DF) }

    # =============================================================================
    # STEP 4: marginal models -- every predictor, every block, every outcome
    # =============================================================================
    cat("\nSTEP 4 -- marginal models\n")
    # Standardise after the complete-case filter.
    fit_one <- function(d, out, pred, covs) {
        # An outcome that is also a covariate is dropped from the covariates.
        cv <- setdiff(covs[covs %in% names(d)], out)
        D <- data.frame(y = suppressWarnings(as.numeric(d[[out]])),
                        x = suppressWarnings(as.numeric(d[[pred]])))
        for (c1 in cv) D[[c1]] <- suppressWarnings(as.numeric(d[[c1]]))
        D <- D[stats::complete.cases(D), , drop = FALSE]
        if (nrow(D) < MIN_DONORS) return(NULL)
        D$y <- .z(D$y); D$x <- .z(D$x)
        if (all(is.na(D$y)) || all(is.na(D$x))) return(NULL)
        # Covariate names are back-quoted.
        f <- stats::as.formula(paste("y ~ x",
            if (length(cv) > 0)
                paste("+", paste(sprintf("`%s`", cv), collapse = " + "))
            else ""))
        fit <- try(stats::lm(f, data = D), silent = TRUE)
        if (inherits(fit, "try-error")) return(NULL)
        s <- summary(fit)
        if (!("x" %in% rownames(s$coefficients))) return(NULL)
        vf <- .safe_vif(fit)
        # delta_R2 over the covariates-only model, as Part 10h reports it.
        f0 <- stats::as.formula(paste("y ~",
            if (length(cv) > 0) paste(sprintf("`%s`", cv), collapse = " + ")
            else "1"))
        fit0 <- try(stats::lm(f0, data = D), silent = TRUE)
        r2_cov <- if (inherits(fit0, "try-error")) NA_real_
                  else summary(fit0)$r.squared
        list(n = nrow(D), beta = unname(s$coefficients["x", 1]),
             se = unname(s$coefficients["x", 2]),
             p = unname(s$coefficients["x", 4]),
             r2 = s$adj.r.squared, r2_full = s$r.squared,
             r2_cov = r2_cov, delta_r2 = s$r.squared - r2_cov,
             vif_x = unname(vf["x"]),
             max_vif = suppressWarnings(max(vf, na.rm = TRUE)))
    }
    mrows <- list()
    for (ds in unique(DF$dataset)) {
        d <- DF[DF$dataset == ds, , drop = FALSE]
        for (bl in names(COVS)) for (pk in PRED)
            for (lab in names(OUT_MAP)) {
            r <- fit_one(d, OUT_MAP[[lab]], pk, COVS[[bl]])
            if (is.null(r)) next
            mrows[[length(mrows) + 1L]] <- data.frame(
                dataset = ds, block = bl, predictor = PRED_KEYS[[pk]], key = pk,
                outcome = lab, outcome_class = unname(OUT_CLASS[[lab]]),
                is_control = lab %in% CONTROL_OUT, n = r$n,
                beta = round(r$beta, 4), se = round(r$se, 4), p_value = r$p,
                delta_r2 = round(r$delta_r2, 4),
                r2_covariates_only = round(r$r2_cov, 4),
                adj_r2 = round(r$r2, 4), vif_pred = round(r$vif_x, 2),
                max_vif = round(r$max_vif, 2), stringsAsFactors = FALSE)
        }
    }
    MARG <- do.call(rbind, mrows)
    if (is.null(MARG)) stop("No marginal model was fitted.", call. = FALSE)
    rownames(MARG) <- NULL
    # BH within dataset x block x predictor, the same family as Part 10h.
    MARG$p_adj <- NA_real_
    for (k in unique(paste(MARG$dataset, MARG$block, MARG$key))) {
        i <- paste(MARG$dataset, MARG$block, MARG$key) == k
        MARG$p_adj[i] <- stats::p.adjust(MARG$p_value[i], method = "BH")
    }
    MARG$sig <- .stars(MARG$p_adj)
    prim <- MARG[MARG$block == PRIMARY_BLOCK &
                 MARG$outcome_class == "clinical", , drop = FALSE]
    for (ds in unique(prim$dataset)) {
        cat(sprintf("\n  %s [%s], clinical outcomes significant after BH:\n",
                    ds, PRIMARY_BLOCK))
        for (pk in PRED) {
            s <- prim[prim$dataset == ds & prim$key == pk, , drop = FALSE]
            cat(sprintf("    %-32s %d of %d | median |beta| %.3f | median dR2 %.3f\n",
                        PRED_KEYS[[pk]], sum(s$p_adj < ALPHA, na.rm = TRUE),
                        nrow(s), stats::median(abs(s$beta), na.rm = TRUE),
                        stats::median(s$delta_r2, na.rm = TRUE)))
        }
    }

    # Module row at the primary block against Part 10h's Multivariate sheet.
    .xc <- NULL; .tier <- NA_character_; .xc_note <- ""
    .pl <- DECL[["PRED_LABEL"]]
    .mod_names <- unique(c("Module Post-IP", "pct_post",
                           if (!is.null(.pl)) names(.pl)[.pl == "Module Post-IP"],
                           if (!is.null(.pl)) unname(.pl["pct_post"])))
    .mod_names <- .mod_names[!is.na(.mod_names) & nzchar(.mod_names)]
    if (is.data.frame(MVSH) && all(c("outcome", "covariate_set",
                                     "dataset") %in% names(MVSH))) {
        .pcol <- NA_character_
        for (cn in intersect(c("predictor", "arm"), names(MVSH)))
            if (any(as.character(MVSH[[cn]]) %in% .mod_names)) { .pcol <- cn; break }
        .cset <- if (PRIMARY_BLOCK %in% as.character(MVSH$covariate_set))
            PRIMARY_BLOCK else NA_character_
        if (is.na(.pcol) || is.na(.cset)) {
            .xc_note <- sprintf(
                paste("no row matched: module name(s) tried [%s] against",
                      "column(s) [%s] holding [%s]; covariate_set '%s' against",
                      "[%s]"),
                paste(.mod_names, collapse = " | "),
                paste(intersect(c("predictor", "arm"), names(MVSH)),
                      collapse = ", "),
                paste(utils::head(unique(unlist(lapply(
                    intersect(c("predictor", "arm"), names(MVSH)),
                    function(cn) as.character(MVSH[[cn]])))), 20),
                    collapse = " | "),
                PRIMARY_BLOCK,
                paste(unique(as.character(MVSH$covariate_set)), collapse = " | "))
            ref <- MVSH[0, , drop = FALSE]
        } else
            ref <- MVSH[as.character(MVSH[[.pcol]]) %in% .mod_names &
                        as.character(MVSH$covariate_set) == .cset, ,
                        drop = FALSE]
        if (!is.null(ref[["tier"]]) && nrow(ref) > 0L) {
            # the all-donor tier is the one with the most donors
            .tt <- stats::aggregate(list(n = suppressWarnings(
                as.numeric(ref$n_donors))), by = list(tier = ref$tier),
                FUN = function(v) max(v, na.rm = TRUE))
            .tier <- .tt$tier[which.max(.tt$n)]
            ref <- ref[ref$tier == .tier, , drop = FALSE]
        } else .tier <- NA_character_
        mine <- MARG[MARG$key == "module" & MARG$block == PRIMARY_BLOCK, ,
                     drop = FALSE]
        cmp <- merge(mine[, c("dataset", "outcome", "n", "beta", "delta_r2")],
                     ref[, intersect(c("dataset", "outcome", "n_donors",
                                       "beta_std", "delta_R2"), names(ref))],
                     by = c("dataset", "outcome"))
        if (nrow(cmp) > 0L) {
            cmp$d_beta <- abs(suppressWarnings(as.numeric(cmp$beta_std)) -
                              cmp$beta)
            cmp$d_n <- abs(suppressWarnings(as.numeric(cmp$n_donors)) - cmp$n)
            cmp$d_r2 <- abs(suppressWarnings(as.numeric(cmp$delta_R2)) -
                            cmp$delta_r2)
            .xc <- cmp
        } else if (!nzchar(.xc_note))
            .xc_note <- sprintf(
                paste("%d Part 10h row(s) selected but none merged on",
                      "dataset+outcome: Part 10h dataset(s) [%s] outcome(s)",
                      "[%s]; comparator dataset(s) [%s] outcome(s) [%s]"),
                nrow(ref), paste(unique(ref$dataset), collapse = " | "),
                paste(utils::head(unique(ref$outcome), 20), collapse = " | "),
                paste(unique(mine$dataset), collapse = " | "),
                paste(utils::head(unique(mine$outcome), 20), collapse = " | "))
    }

    # =============================================================================
    # STEP 5: the duels
    # =============================================================================
    cat("\nSTEP 5 -- duels: module against each baseline\n")
    duel_one <- function(d, out, covs) {
        # Standardise after the complete-case filter.
        cv <- setdiff(covs[covs %in% names(d)], out)
        function(rival) {
            D <- data.frame(y = suppressWarnings(as.numeric(d[[out]])),
                            module = suppressWarnings(as.numeric(d[["module"]])),
                            rival = suppressWarnings(as.numeric(d[[rival]])))
            for (c1 in cv) D[[c1]] <- suppressWarnings(as.numeric(d[[c1]]))
            D <- D[stats::complete.cases(D), , drop = FALSE]
            if (nrow(D) < MIN_DONORS) return(NULL)
            D$y <- .z(D$y); D$module <- .z(D$module); D$rival <- .z(D$rival)
            if (anyNA(D$y) || anyNA(D$module) || anyNA(D$rival)) return(NULL)
            base <- paste(if (length(cv) > 0)
                paste(sprintf("`%s`", cv), collapse = " + ") else "1")
            f_j <- stats::as.formula(paste("y ~ module + rival +", base))
            f_m <- stats::as.formula(paste("y ~ module +", base))
            f_r <- stats::as.formula(paste("y ~ rival +", base))
            jt <- try(stats::lm(f_j, data = D), silent = TRUE)
            fm <- try(stats::lm(f_m, data = D), silent = TRUE)
            fr <- try(stats::lm(f_r, data = D), silent = TRUE)
            if (any(vapply(list(jt, fm, fr), function(z) inherits(z, "try-error"),
                           logical(1)))) return(NULL)
            s <- summary(jt)$coefficients
            vf <- .safe_vif(jt)
            list(n = nrow(D),
                 p_module = if ("module" %in% rownames(s)) s["module", 4] else NA_real_,
                 p_rival = if ("rival" %in% rownames(s)) s["rival", 4] else NA_real_,
                 beta_module = if ("module" %in% rownames(s)) s["module", 1] else NA_real_,
                 beta_rival = if ("rival" %in% rownames(s)) s["rival", 1] else NA_real_,
                 r2_module = summary(fm)$adj.r.squared,
                 r2_rival = summary(fr)$adj.r.squared,
                 r2_joint = summary(jt)$adj.r.squared,
                 vif_joint = suppressWarnings(max(vf[c("module", "rival")],
                                                  na.rm = TRUE)))
        }
    }
    drows <- list()
    for (ds in unique(DF$dataset)) {
        d <- DF[DF$dataset == ds, , drop = FALSE]
        for (rival in setdiff(PRED, "module")) {
            for (lab in names(OUT_MAP)) {
                r <- duel_one(d, OUT_MAP[[lab]],
                              COVS[[PRIMARY_BLOCK]])(rival)
                if (is.null(r)) next
                drows[[length(drows) + 1L]] <- data.frame(
                    dataset = ds, duel = sprintf("module_vs_%s", rival),
                    rival_label = PRED_KEYS[[rival]], outcome = lab,
                    outcome_class = unname(OUT_CLASS[[lab]]),
                    is_control = lab %in% CONTROL_OUT, n = r$n,
                    beta_module = round(r$beta_module, 4),
                    beta_rival = round(r$beta_rival, 4),
                    p_module = r$p_module, p_rival = r$p_rival,
                    r2_module_only = round(r$r2_module, 4),
                    r2_rival_only = round(r$r2_rival, 4),
                    r2_joint = round(r$r2_joint, 4),
                    dR2_module_minus_rival = round(r$r2_module - r$r2_rival, 4),
                    vif_joint = round(r$vif_joint, 2),
                    stringsAsFactors = FALSE)
            }
        }
    }
    DUEL <- do.call(rbind, drows)
    if (!is.null(DUEL)) {
        rownames(DUEL) <- NULL
        # BH within a duel, over the outcomes, on each side.
        DUEL$p_module_adj <- NA_real_; DUEL$p_rival_adj <- NA_real_
        for (k in unique(paste(DUEL$dataset, DUEL$duel))) {
            i <- paste(DUEL$dataset, DUEL$duel) == k
            DUEL$p_module_adj[i] <- stats::p.adjust(DUEL$p_module[i], method = "BH")
            DUEL$p_rival_adj[i] <- stats::p.adjust(DUEL$p_rival[i], method = "BH")
        }
        DUEL$verdict <- ifelse(
            DUEL$p_module_adj < ALPHA & DUEL$p_rival_adj >= ALPHA, "module only",
            ifelse(DUEL$p_module_adj >= ALPHA & DUEL$p_rival_adj < ALPHA,
                   "rival only",
            ifelse(DUEL$p_module_adj < ALPHA & DUEL$p_rival_adj < ALPHA, "both",
                   "neither")))
        DUEL$readable <- ifelse(is.finite(DUEL$vif_joint) &
                                DUEL$vif_joint >= VIF_READABLE,
                                sprintf(paste0("joint VIF >= %g: not ",
                                               "separately readable"),
                                        VIF_READABLE), "")
        for (k in unique(paste(DUEL$dataset, DUEL$duel))) {
            s <- DUEL[paste(DUEL$dataset, DUEL$duel) == k & !DUEL$is_control, ,
                      drop = FALSE]
            cat(sprintf("  %-46s module only %d, rival only %d, both %d, neither %d of %d%s\n",
                        k, sum(s$verdict == "module only"),
                        sum(s$verdict == "rival only"), sum(s$verdict == "both"),
                        sum(s$verdict == "neither"), nrow(s),
                        if (any(nzchar(s$readable)))
                            sprintf(" | %d row(s) not separately readable",
                                    sum(nzchar(s$readable))) else ""))
        }
    }


    # =============================================================================
    # STEP 6: stress tests
    # =============================================================================
    cat("\nSTEP 6 -- stress tests\n")

    # == 6a: covariate-block sweep ============================================
    sweep_rows <- list()
    for (ds in unique(MARG$dataset)) for (pk in PRED) for (bl in names(COVS)) {
        s <- MARG[MARG$dataset == ds & MARG$key == pk & MARG$block == bl &
                  !MARG$is_control, , drop = FALSE]
        if (nrow(s) == 0L) next
        sweep_rows[[length(sweep_rows) + 1L]] <- data.frame(
            dataset = ds, predictor = PRED_KEYS[[pk]], block = bl,
            n_outcomes = nrow(s), n_sig = sum(s$p_adj < ALPHA, na.rm = TRUE),
            median_abs_beta = round(stats::median(abs(s$beta), na.rm = TRUE), 4),
            median_n = stats::median(s$n),
            worst_vif = round(suppressWarnings(max(s$vif_pred, na.rm = TRUE)), 2),
            stringsAsFactors = FALSE)
    }
    SWEEP <- do.call(rbind, sweep_rows)
    if (!is.null(SWEEP)) {
        rownames(SWEEP) <- NULL
        cat("  block sweep (non-control outcomes, significant after BH):\n")
        for (ds in unique(SWEEP$dataset)) {
            s <- SWEEP[SWEEP$dataset == ds, , drop = FALSE]
            for (pk in unique(s$predictor))
                cat(sprintf("    %-10s %-32s %s\n", ds, pk,
                            paste(sprintf("%s %d/%d", s$block[s$predictor == pk],
                                          s$n_sig[s$predictor == pk],
                                          s$n_outcomes[s$predictor == pk]),
                                  collapse = "  ")))
        }
    }

    # Bootstrap outcomes: clinical, ranked by |module beta|, capped.
    .bo <- MARG[MARG$block == PRIMARY_BLOCK & MARG$key == "module" &
                MARG$outcome_class == "clinical", , drop = FALSE]
    .bo <- stats::aggregate(list(b = abs(.bo$beta)),
                            by = list(outcome = .bo$outcome), FUN = max)
    BOOT_LABS <- utils::head(.bo$outcome[order(-.bo$b)], BOOT_MAX_OUTCOMES)
    cat(sprintf("  bootstrap outcomes (top %d clinical by |module beta| at %s): %s\n",
                length(BOOT_LABS), PRIMARY_BLOCK,
                paste(BOOT_LABS, collapse = ", ")))

    # == 6b: donor bootstrap of the R2 difference =============================
    # Resampling is by donor.
    boot_rows <- list()
    if (!is.null(DUEL)) for (ds in unique(DF$dataset)) {
        d <- DF[DF$dataset == ds, , drop = FALSE]
        donors <- unique(d$.donor[is.finite(d$module)])
        if (length(donors) < MIN_DONORS) next
        cv0 <- COVS[[PRIMARY_BLOCK]][COVS[[PRIMARY_BLOCK]] %in% names(d)]
        # One task per rival x outcome, seeded from SEG_SEED + task index.
        tasks <- list()
        for (rival in setdiff(PRED, "module")) for (lab in BOOT_LABS)
            tasks[[length(tasks) + 1L]] <- list(rival = rival, lab = lab)
        bt_worker <- p10hc_fn(function(ti) {
            tk <- tasks[[ti]]
            out <- OUT_MAP[[tk$lab]]
            cv <- setdiff(cv0, out)
            dd <- d[, intersect(c(".donor", out, "module", tk$rival, cv),
                                names(d)), drop = FALSE]
            set.seed(SEG_SEED + ti)
            base <- if (length(cv) > 0)
                paste("+", paste(sprintf("`%s`", cv), collapse = " + ")) else ""
            f_m <- stats::as.formula(paste("y ~ module", base))
            f_r <- stats::as.formula(paste("y ~ rival", base))
            deltas <- rep(NA_real_, BOOT_N)
            # Both replicate fits must be at this rank; others are skipped and counted.
            rank_want <- 2L + length(cv)
            n_rank <- 0L      # replicates dropped for rank deficiency
            n_small <- 0L     # replicates dropped for too few complete donors
            n_err <- 0L       # replicates dropped for a fit error or NA input
            for (b in seq_len(BOOT_N)) {
                take <- sample(donors, length(donors), replace = TRUE)
                db <- dd[match(take, dd$.donor), , drop = FALSE]
                # Standardise after the complete-case filter.
                D <- data.frame(y = suppressWarnings(as.numeric(db[[out]])),
                                module = suppressWarnings(as.numeric(db[["module"]])),
                                rival = suppressWarnings(as.numeric(db[[tk$rival]])),
                                check.names = FALSE)
                for (c1 in cv) D[[c1]] <- suppressWarnings(as.numeric(db[[c1]]))
                D <- D[stats::complete.cases(D), , drop = FALSE]
                if (nrow(D) < MIN_DONORS) { n_small <- n_small + 1L; next }
                D$y <- .z(D$y); D$module <- .z(D$module); D$rival <- .z(D$rival)
                if (anyNA(D$y) || anyNA(D$module) || anyNA(D$rival)) {
                    n_err <- n_err + 1L; next
                }
                fm <- try(stats::lm(f_m, data = D), silent = TRUE)
                fr <- try(stats::lm(f_r, data = D), silent = TRUE)
                if (inherits(fm, "try-error") || inherits(fr, "try-error")) {
                    n_err <- n_err + 1L; next
                }
                # Skip unless both fits are at full rank.
                if (!identical(as.integer(fm$rank), rank_want) ||
                    !identical(as.integer(fr$rank), rank_want)) {
                    n_rank <- n_rank + 1L; next
                }
                # r.squared: the covariates-only R2 cancels, as in delta_r2.
                deltas[b] <- summary(fm)$r.squared - summary(fr)$r.squared
            }
            dv <- deltas[is.finite(deltas)]
            if (length(dv) < BOOT_N / 10) return(NULL)
            ci <- stats::quantile(dv, c(0.025, 0.975), names = FALSE)
            data.frame(dataset = ds, duel = sprintf("module_vs_%s", tk$rival),
                       rival_label = PRED_KEYS[[tk$rival]], outcome = tk$lab,
                       outcome_class = "clinical",
                       is_control = tk$lab %in% CONTROL_OUT,
                       n_boot = length(dv),
                       n_requested = BOOT_N,
                       n_skip_rank = n_rank,
                       n_skip_small = n_small,
                       n_skip_error = n_err,
                       median_dR2 = round(stats::median(dv), 4),
                       ci_lo = round(ci[1], 4), ci_hi = round(ci[2], 4),
                       pct_module_ahead = round(100 * mean(dv > 0), 1),
                       ci_excludes_zero = ci[1] > 0 | ci[2] < 0,
                       stringsAsFactors = FALSE)
        },
        tasks = tasks, OUT_MAP = OUT_MAP, cv0 = cv0, d = d, donors = donors,
        .z = .z, MIN_DONORS = MIN_DONORS, BOOT_N = BOOT_N, SEG_SEED = SEG_SEED,
        PRED_KEYS = PRED_KEYS, CONTROL_OUT = CONTROL_OUT, ds = ds)
        bt <- p10hc_map(seq_along(tasks), bt_worker)
        for (z in bt) if (!is.null(z))
            boot_rows[[length(boot_rows) + 1L]] <- z
    }
    BOOT <- do.call(rbind, boot_rows)
    if (!is.null(BOOT)) {
        rownames(BOOT) <- NULL
        cat(sprintf(paste("  donor bootstrap: %d row(s) at %d replicate(s),",
                          "%d with a CI excluding zero\n"),
                    nrow(BOOT), BOOT_N, sum(BOOT$ci_excludes_zero)))
        cat(sprintf(paste("  replicates dropped: %d rank-deficient, %d below",
                          "MIN_DONORS, %d fit error; %d row(s) below the",
                          "requested count\n"),
                    sum(BOOT$n_skip_rank), sum(BOOT$n_skip_small),
                    sum(BOOT$n_skip_error), sum(BOOT$n_boot < BOOT$n_requested)))
    }

    # == 6c: leave-one-donor-out ==============================================
    loo_rows <- list()
    for (ds in unique(DF$dataset)) {
        d <- DF[DF$dataset == ds, , drop = FALSE]
        cv <- COVS[[PRIMARY_BLOCK]][COVS[[PRIMARY_BLOCK]] %in% names(d)]
        for (pk in PRED) for (lab in names(OUT_MAP)) {
            out <- OUT_MAP[[lab]]
            full <- fit_one(d, out, pk, cv)
            if (is.null(full)) next
            donors <- unique(d$.donor)
            bs <- rep(NA_real_, length(donors)); ps <- bs
            for (i in seq_along(donors)) {
                r <- fit_one(d[d$.donor != donors[i], , drop = FALSE], out, pk, cv)
                if (is.null(r)) next
                bs[i] <- r$beta; ps[i] <- r$p
            }
            ok <- is.finite(bs)
            if (!any(ok)) next
            shift <- abs(bs[ok] - full$beta)
            loo_rows[[length(loo_rows) + 1L]] <- data.frame(
                dataset = ds, predictor = PRED_KEYS[[pk]], outcome = lab,
                outcome_class = unname(OUT_CLASS[[lab]]),
                is_control = lab %in% CONTROL_OUT, n_donors = sum(ok),
                beta_full = round(full$beta, 4),
                max_abs_shift = round(max(shift), 4),
                max_pct_shift = if (abs(full$beta) > 1e-8)
                    round(100 * max(shift) / abs(full$beta), 1) else NA_real_,
                worst_donor = donors[ok][which.max(shift)],
                n_flip_sign = sum(sign(bs[ok]) != sign(full$beta)),
                n_lose_sig = if (full$p < ALPHA)
                    sum(ps[ok] >= ALPHA, na.rm = TRUE) else NA_integer_,
                stringsAsFactors = FALSE)
        }
    }
    LOO <- do.call(rbind, loo_rows)
    if (!is.null(LOO)) {
        rownames(LOO) <- NULL
        l <- LOO[!LOO$is_control, , drop = FALSE]
        cat(sprintf("  leave-one-donor-out: %d row(s); %d lose significance to a single donor, %d flip sign\n",
                    nrow(l), sum(l$n_lose_sig > 0, na.rm = TRUE),
                    sum(l$n_flip_sign > 0, na.rm = TRUE)))
        for (pk in unique(l$predictor)) {
            s <- l[l$predictor == pk, , drop = FALSE]
            i <- which.max(s$max_pct_shift)
            if (length(i) == 1L)
                cat(sprintf("    %-32s largest shift %.1f%% (%s on %s)\n", pk,
                            s$max_pct_shift[i], s$worst_donor[i], s$outcome[i]))
        }
    }

    # =============================================================================
    # STEP 7: figures
    # =============================================================================
    cat("\nSTEP 7 -- figures\n")
    .n_fig <- 0L
    save_fig <- function(p, base, w = 11, h = 8, dpi = 300) {
        render <- function(open_fn, path) {
            tryCatch({
                open_fn(); on.exit(grDevices::dev.off(), add = TRUE); print(p)
                TRUE
            }, error = function(e) {
                cat(sprintf("    figure not written (%s): %s\n", basename(path),
                            conditionMessage(e))); FALSE })
        }
        pp <- paste0(base, ".pdf"); gp <- paste0(base, ".png")
        ok1 <- render(function() grDevices::pdf(pp, width = w, height = h), pp)
        ok2 <- render(function() grDevices::png(gp, width = w * dpi,
                                                height = h * dpi, res = dpi), gp)
        if (ok1 && ok2) .n_fig <<- .n_fig + 1L
        invisible(NULL)
    }
    theme_hc <- function() theme_bw(base_size = 11) +
        theme(plot.title = element_text(face = "bold", size = 12),
              plot.subtitle = element_text(size = 8, colour = "grey30"),
              strip.background = element_rect(fill = "#1F3864"),
              strip.text = element_text(colour = "white", face = "bold", size = 9),
              panel.grid.minor = element_blank(),
              legend.position = "bottom",
              plot.background = element_rect(fill = "white", colour = NA))

    # == 7a: coefficient forest ===============================================
    fd <- MARG[MARG$block == PRIMARY_BLOCK, , drop = FALSE]
    if (nrow(fd) > 0) {
        # position_dodge() puts the first level at the bottom; breaks fix the legend.
        fd$pred_f <- factor(fd$predictor, levels = rev(unname(PRED_KEYS)))
        fd$out_f <- factor(fd$outcome, levels = rev(names(OUT_MAP)))
        fd$lo <- fd$beta - 1.96 * fd$se; fd$hi <- fd$beta + 1.96 * fd$se
        p1 <- ggplot(fd, aes(x = beta, y = out_f, colour = pred_f)) +
            geom_vline(xintercept = 0, colour = "grey50", linetype = "dashed") +
            geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0,
                           position = position_dodge(width = 0.7),
                           linewidth = 0.4) +
            geom_point(aes(shape = p_adj < ALPHA), size = 2.2,
                       position = position_dodge(width = 0.7)) +
            scale_shape_manual(values = c(`TRUE` = 19, `FALSE` = 1),
                               name = "BH < 0.05", na.translate = FALSE) +
            scale_colour_manual(values = PRED_FILL, name = NULL,
                                breaks = unname(PRED_KEYS)) +
            facet_wrap(~dataset, ncol = 2) +
            labs(title = sprintf("Standardised coefficient per predictor [%s, %s]",
                                 ARM_TAG, PRIMARY_BLOCK),
                 subtitle = paste("Same donors, same covariate block, same",
                                  "outcomes; predictors read downwards in",
                                  "legend order. Open points do not clear BH.",
                                  "Control outcomes are included and labelled",
                                  "in the Marginal sheet.",
                                  "\nDIRECTION: higher-is-better outcomes",
                                  "(MMSE, CASI) and higher-is-worse outcomes",
                                  "(staging, cognitive status) sit on opposite",
                                  "sides of zero for the same biological",
                                  "direction."),
                 x = "beta (outcome SD per predictor SD)", y = NULL) +
            theme_hc()
        print(p1)
        # Height scales with outcomes and predictors.
        save_fig(p1, file.path(OUT_DIR, paste0(FILE_PREFIX, "_Comparator_Forest")),
                 w = 11,
                 h = max(5, length(OUT_MAP) * 0.11 * length(PRED) + 3))
    }

    # == 7b: paired delta-R2 ==================================================
    if (!is.null(DUEL) && nrow(DUEL) > 0) {
        dd <- DUEL[DUEL$outcome_class == "clinical", , drop = FALSE]
        p2 <- ggplot(dd, aes(x = r2_rival_only, y = r2_module_only,
                             colour = rival_label)) +
            geom_abline(slope = 1, intercept = 0, colour = "grey50",
                        linetype = "dashed") +
            geom_point(size = 2.4, alpha = 0.85) +
            scale_colour_manual(values = PRED_FILL, name = NULL) +
            facet_wrap(~dataset, ncol = 2) +
            coord_equal() +
            labs(title = sprintf("Adjusted R2: the module against each baseline [%s, %s]",
                                 ARM_TAG, PRIMARY_BLOCK),
                 subtitle = paste("One point per outcome. Above the line the",
                                  "module explains more; below it the baseline",
                                  "does. The bootstrap CI for each point is in",
                                  "the Stress_Bootstrap sheet."),
                 x = "adjusted R2, baseline alone",
                 y = "adjusted R2, module alone") +
            theme_hc()
        print(p2)
        save_fig(p2, file.path(OUT_DIR, paste0(FILE_PREFIX, "_Comparator_DeltaR2")),
                 w = 10, h = 6)

        tally <- as.data.frame(table(dataset = dd$dataset, duel = dd$rival_label,
                                     verdict = dd$verdict),
                               stringsAsFactors = FALSE)
        tally$verdict <- factor(tally$verdict,
                                levels = c("module only", "both", "neither",
                                           "rival only"))
        p3 <- ggplot(tally, aes(x = duel, y = Freq, fill = verdict)) +
            geom_col(colour = "black", linewidth = 0.3, width = 0.65) +
            scale_fill_manual(values = c("module only" = "#2D2D2D",
                                         "both" = "#7B7B7B",
                                         "neither" = "#D9D9D9",
                                         "rival only" = "#B2182B"), name = NULL) +
            facet_wrap(~dataset, ncol = 2) +
            coord_flip() +
            labs(title = sprintf("Encompassing-test verdicts [%s, %s]", ARM_TAG,
                                 PRIMARY_BLOCK),
                 subtitle = paste("y ~ covariates + module + rival, BH-corrected",
                                  "within a duel over the non-control outcomes.",
                                  "Red is the rival winning outright and is",
                                  "reported as it falls."),
                 x = NULL, y = "outcomes") +
            theme_hc()
        print(p3)
        save_fig(p3, file.path(OUT_DIR, paste0(FILE_PREFIX, "_Comparator_Duels")),
                 w = 10, h = max(5, 2 + 0.7 * length(setdiff(PRED, "module"))))
    }

    # =============================================================================
    # STEP 8: output
    # =============================================================================
    cat("\nSTEP 8 -- output\n")
    wb <- list(Donor_Frame = DF[, intersect(c(".donor", "dataset", "module", "B1",
                                              "B2", "B3", "B4", "B5", "n_micro",
                                              unname(OUT_MAP)),
                                            names(DF)), drop = FALSE],
               Marginal = MARG)
    if (length(.ams_report) > 0)
        wb[["Score_Controls"]] <- do.call(rbind, lapply(.ams_report, function(r)
            data.frame(dataset = r$ds, n_signature = r$n_sig,
                       n_control = r$n_ctrl, n_bins_used = r$n_bins_used,
                       n_control_mito = r$n_ctrl_mito,
                       n_control_in_signature = r$n_ctrl_in_sig,
                       n_bins_short_of_ctrl = r$n_short,
                       mean_expr_signature = round(r$mean_sig, 5),
                       mean_expr_control = round(r$mean_ctrl, 5),
                       stringsAsFactors = FALSE)))
    if (!is.null(COR_DF))  wb[["Predictor_Correlations"]] <- COR_DF
    if (!is.null(.xc))     wb[["Part10h_Crosscheck"]] <- .xc
    if (!is.null(DUEL))    wb[["Duels"]] <- DUEL
    if (!is.null(SWEEP))   wb[["Stress_Block_Sweep"]] <- SWEEP
    if (!is.null(BOOT))    wb[["Stress_Bootstrap"]] <- BOOT
    if (!is.null(LOO))     wb[["Stress_LeaveOneDonor"]] <- LOO
    if (length(.seg_ok) > 0)
        wb[["Changepoint"]] <- do.call(rbind, lapply(.seg_ok, function(r)
            data.frame(dataset = r$ds, breakpoint = round(r$psi, 5),
                       se = round(r$se, 5), n_used = r$n_used,
                       n_total = r$n_total, pct_past = round(r$pct_past, 2),
                       stringsAsFactors = FALSE)))
    wb <- wb[!vapply(wb, is.null, logical(1))]
    xlsx_path <- file.path(OUT_DIR, paste0(FILE_PREFIX, "_Comparator.xlsx"))
    tryCatch({
        write_xlsx(wb, xlsx_path)
        cat(sprintf("  %d sheet(s) written to %s\n", length(wb),
                    basename(xlsx_path)))
    }, error = function(e)
        cat(sprintf("  workbook not written: %s\n", conditionMessage(e))))

    rds_out <- list(donor_frame = DF, per_cell = .percell,
                    marginal = MARG, duels = DUEL,
                    correlations = COR_DF, block_sweep = SWEEP, bootstrap = BOOT,
                    leave_one_out = LOO, changepoint = .seg_ok,
                    meta = list(arm_tag = ARM_TAG, arm_requested = arm_tag,
                                config_file = config_file, norm_tag = NORM_TAG,
                                GOI = GOI, DFG = DFG, supertype = PRIMARY_STYPE,
                                microglial_supertypes = MICROGLIAL_SUPERTYPES,
                                predictors = PRED_KEYS[PRED],
                                primary_block = PRIMARY_BLOCK,
                                outcomes = OUT_MAP, control_outcomes = CONTROL_OUT,
                            covariate_columns = COVS,
                            cov_map_route = COV_ROUTE,
                            outcome_map_route = OUT_ROUTE,
                                part10h_rds = P10H_RDS,
                                written = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                                R_version = R.version.string,
                                schema = "M10HC v1"))
    rds_path <- file.path(OUT_DIR, paste0(FILE_PREFIX, "_Comparator_results.rds"))
    tryCatch({
        saveRDS(rds_out, rds_path)
        cat(sprintf("  %s written (%.2f MB)\n", basename(rds_path),
                    file.size(rds_path) / 1e6))
    }, error = function(e)
        cat(sprintf("  results RDS not written: %s\n", conditionMessage(e))))

    # =============================================================================
    # STEP 9: GUARDRAILS
    # =============================================================================
    cat("\nSTEP 9 -- guardrails\n")

    if (is.null(.xc))
        gr10hc("GR-10HC-13", "WARN",
               sprintf(paste("Part 10h's Multivariate sheet could not be",
                             "matched to the module rows. %s"),
                       if (nzchar(.xc_note)) .xc_note
                       else "No Multivariate sheet in the RDS."))
    else {
        .bad <- sum(.xc$d_beta > 0.01 | .xc$d_n > 0 | .xc$d_r2 > 0.01,
                    na.rm = TRUE)
        gr10hc("GR-10HC-13", if (.bad > 0L) "FAIL" else "PASS",
               sprintf(paste("module row at %s against Part 10h's own",
                             "Multivariate sheet%s over %d matched",
                             "dataset-outcome row(s): max |d beta| %.4f, max",
                             "|d n| %d, max |d delta-R2| %.4f. %s"),
                       PRIMARY_BLOCK,
                       if (is.na(.tier)) "" else sprintf(" (tier '%s')", .tier),
                       nrow(.xc), max(.xc$d_beta, na.rm = TRUE),
                       as.integer(max(.xc$d_n, na.rm = TRUE)),
                       max(.xc$d_r2, na.rm = TRUE),
                       if (.bad > 0L)
                           paste(.bad, "row(s) disagree beyond tolerance")
                       else "Matches Part 10h on the shared rows."))
    }

    # Module against the cluster proportion (B1), reported either way.
    if (!is.null(DUEL)) {
        b1 <- DUEL[DUEL$duel == "module_vs_B1" &
                   DUEL$outcome_class == "clinical", , drop = FALSE]
        if (nrow(b1) == 0L)
            gr10hc("GR-10HC-06", "WARN",
                   "the cluster-proportion duel produced no row")
        else {
            won <- sum(b1$verdict == "module only")
            lost <- sum(b1$verdict == "rival only")
            gr10hc("GR-10HC-06",
                   if (lost > won) "WARN" else "PASS",
                   sprintf(paste("pre-registered row: module against the cluster",
                                 "proportion, over %d clinical outcome(s) and",
                                 "%d dataset(s): module only %d, rival only %d,",
                                 "both %d, neither %d. %s"),
                           length(unique(b1$outcome)), length(unique(b1$dataset)),
                           won, lost, sum(b1$verdict == "both"),
                           sum(b1$verdict == "neither"),
                           if (lost > won)
                               "The cluster proportion wins more outcomes."
                           else if (won == 0L && lost == 0L)
                               "Neither predictor clears BH in this duel."
                           else "The module wins or ties."))
        }
    }

    # Do the three module scores give the same duel tally?
    .sk <- intersect(SCORE_KEYS, PRED)
    if (!is.null(DUEL) && length(.sk) >= 2L) {
        .sd <- DUEL[DUEL$duel %in% sprintf("module_vs_%s", .sk) &
                    DUEL$outcome_class == "clinical", , drop = FALSE]
        .tal <- lapply(.sk, function(k) {
            s <- .sd[.sd$duel == sprintf("module_vs_%s", k), , drop = FALSE]
            c(mod = sum(s$verdict == "module only"),
              riv = sum(s$verdict == "rival only"),
              both = sum(s$verdict == "both"),
              nei = sum(s$verdict == "neither"))
        })
        names(.tal) <- .sk
        .agree <- length(unique(vapply(.tal, function(v)
            paste(v, collapse = "/"), character(1)))) == 1L
        .rho <- if (is.null(COR_DF)) NA_real_ else {
            .lab <- unname(PRED_KEYS[.sk])
            .p <- COR_DF[COR_DF$a %in% .lab & COR_DF$b %in% .lab, ,
                         drop = FALSE]
            if (nrow(.p) == 0L) NA_real_ else min(.p$spearman, na.rm = TRUE)
        }
        gr10hc("GR-10HC-15", if (.agree) "PASS" else "WARN",
               sprintf(paste("score family (%s), clinical rows, module only /",
                             "rival only / both / neither: %s. Lowest pairwise",
                             "Spearman among the scores %s. %s"),
                       paste(unname(PRED_KEYS[.sk]), collapse = ", "),
                       paste(vapply(.sk, function(k)
                           sprintf("%s %d/%d/%d/%d", k, .tal[[k]][["mod"]],
                                   .tal[[k]][["riv"]], .tal[[k]][["both"]],
                                   .tal[[k]][["nei"]]), character(1)),
                           collapse = " | "),
                       if (is.finite(.rho)) sprintf("%.3f", .rho)
                       else "not computed",
                       if (.agree)
                           paste("The three scores give the same verdict",
                                 "tally.")
                       else
                           paste("The tallies differ between scoring",
                                 "methods.")))
    }

    # Collinearity between the module and each rival.
    if (!is.null(DUEL)) {
        nr <- sum(is.finite(DUEL$vif_joint) & DUEL$vif_joint >= VIF_READABLE)
        gr10hc("GR-10HC-07", if (nr > 0L) "WARN" else "PASS",
               sprintf(paste("joint VIF on the module term: %d of %d duel row(s)",
                             "at or above %g%s"),
                       nr, nrow(DUEL), VIF_READABLE,
                       if (nr > 0L)
                           "; neither coefficient is separately readable there"
                       else ""))
    }

    # Stress-test summary.
    .fragile <- if (is.null(LOO)) NA_integer_ else
        sum(LOO$n_lose_sig > 0, na.rm = TRUE)
    gr10hc("GR-10HC-08",
           if (is.null(LOO) || is.null(BOOT)) "WARN" else "PASS",
           sprintf(paste("stress tests: block sweep over %d block(s); donor",
                         "bootstrap %s; leave-one-donor-out %s"),
                   length(COVS),
                   if (is.null(BOOT)) "not run"
                   else sprintf(paste("%d row(s) over the %d outcome(s) with",
                                      "the largest module coefficient (of %d",
                                      "non-control), %d with a CI excluding",
                                      "zero"),
                                nrow(BOOT), length(BOOT_LABS),
                                sum(!(names(OUT_MAP) %in% CONTROL_OUT)),
                                sum(BOOT$ci_excludes_zero)),
                   if (is.null(LOO)) "not run"
                   else sprintf("%d row(s), %d losing significance to one donor",
                                nrow(LOO), .fragile)))

    # Bootstrap replicate accounting.
    if (is.null(BOOT))
        gr10hc("GR-10HC-17", "WARN",
               "donor bootstrap not run")
    else {
        .sk_r <- sum(BOOT$n_skip_rank)
        .sk_s <- sum(BOOT$n_skip_small)
        .sk_e <- sum(BOOT$n_skip_error)
        .short <- sum(BOOT$n_boot < BOOT$n_requested)
        .worst <- if (nrow(BOOT)) min(BOOT$n_boot) else NA_integer_
        gr10hc("GR-10HC-17",
               if (.short == 0L) "PASS" else "WARN",
               sprintf(paste("bootstrap replicate accounting over %d row(s) at",
                             "%d requested replicate(s): %d dropped for rank",
                             "deficiency, %d for fewer than MIN_DONORS",
                             "complete donors, %d for a fit error. %s"),
                       nrow(BOOT), BOOT_N, .sk_r, .sk_s, .sk_e,
                       if (.short == 0L)
                           "Every row used its full requested count."
                       else sprintf(paste("%d row(s) fell short of the",
                                          "requested count, the smallest at",
                                          "%d replicate(s)."),
                                    .short, as.integer(.worst))))
    }

    gr10hc("GR-10HC-12", if (is.null(.cl)) "WARN" else "PASS",
           if (is.null(.cl))
               sprintf(paste("no parallel cluster was created (PAR_ENABLE %s,",
                             "PAR_WORKERS %d), so every stage ran",
                             "sequentially"), PAR_ENABLE, PAR_WORKERS)
           else sprintf(paste("%d PSOCK worker(s). Parallelised: the AUCell",
                              "scan over cell chunks of %d, and the donor",
                              "bootstrap over %d rival x outcome task(s) per",
                              "dataset"),
                        PAR_WORKERS, AUCELL_CELL_BLOCK,
                        length(setdiff(PRED, "module")) * length(BOOT_LABS)))

    gr10hc("GR-10HC-09", if (.n_fig == 0L) "WARN" else "PASS",
           sprintf(paste("%d figure(s) written as PDF and PNG"), .n_fig))

    .files <- c(rds_path, xlsx_path)
    .have <- file.exists(.files)
    gr10hc("GR-10HC-10", if (all(.have)) "PASS" else "FAIL",
           sprintf("outputs: %s", paste(sprintf("%s %s", basename(.files),
                                                ifelse(.have, "written",
                                                       "MISSING")),
                                        collapse = " | ")))

    # =============================================================================
    # LEDGER
    # =============================================================================
    .gr_lines <- c(
        strrep("=", 72),
        "Part 10h Prognosis Comparator - Guardrail Log [CP10K]",
        paste("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
        paste("Arm :", ARM_TAG),
        "Unit: donor",
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
    writeLines(.gr_lines, GR_LOG_10HC)
    cat(sprintf("\n  Final guardrail tally: %d PASS / %d WARN / %d FAIL\n",
                .gr_pass, .gr_warn, .gr_fail))

    if (file.exists(INCOMPLETE_MARKER_10HC)) unlink(INCOMPLETE_MARKER_10HC)
    .p10hc_completed <- TRUE

    cat("\n", strrep("=", 72), "\n")
    cat("  Part10h Comparator complete [", ARM_TAG, "]\n", sep = "")
    cat("  Output directory:", OUT_DIR, "\n")
    cat(strrep("=", 72), "\n")

    invisible(NULL)
}

cep_run_arms(PART_DIR_TAG, run_one_arm)
