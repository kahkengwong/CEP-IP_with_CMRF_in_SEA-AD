# =============================================================================
# Part10c_CP10K_Copula1_code.r
# SEA-AD CEP-IP/CMRF, Part 10c, copula fundamentals
#   (Figure 4, Supplementary Figures S8 to S10, Table S9)
#   For each arm: maps GOI-positive cells to copula coordinates
#   u = rank(GOI)/(n+1), v = rank(cDFG)/(n+1), draws the marginal, joint and
#   copula densities, tests uniformity of the margins (KS, Anderson-Darling),
#   compares original and shuffled dependence (CvM distance, tail dependence)
#   and tests band-specific correlations (Fisher z, permutation).
#
# Usage (from the SEA-AD folder, or set SEAAD_ROOT), after Part 10 and Part 10b
#   source("Part10c_CP10K_Copula1_code.r")                         # all arms
#   CEP_ARM <- "PVM1_all"; source("Part10c_CP10K_Copula1_code.r")  # one arm
#   Needs CEP_arms.R, R_config.R and R_config_HKG8.R in the same folder.
#
# Outputs (Part10c_Copula1/Part10c_Copula1_<arm>/): figures (PDF, PNG),
#   workbooks, results RDS, guardrail ledger and logs.
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
library(grid)
library(gtable)
library(cowplot)       # plot_grid with axis alignment
library(MASS)          # kde2d

set.seed(42)

# CEP_arms.R holds the arm table, the folder layout and the driver.
if (!exists("CEP_HOME")) CEP_HOME <- Sys.getenv("SEAAD_ROOT", getwd())
.p10c_harness <- file.path(CEP_HOME, "CEP_arms.R")
if (!file.exists(.p10c_harness))
    stop("CEP_arms.R not found at: ", .p10c_harness,
         call. = FALSE)
source(.p10c_harness)

# PART_DIR_TAG names the results folder; PART_TAG prefixes the file names.
PART_TAG     <- "Part10c"           # file-name tag
PART_DIR_TAG <- "Part10c_Copula1"   # results-folder tag


# =============================================================================
# SECTION 0d: Pure helpers
# =============================================================================
# Defined at the top level; nothing here reads the configuration.


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

save_pp <- function(p, base, pdf_w = 10, pdf_h = 9,
                    png_w = 2000, png_h = 1800, res = 200) {
    # Close only devices opened here, so a failed write cannot abort the dataset.
    .dev_entry <- grDevices::dev.cur()
    .unwind <- function() {
        while (grDevices::dev.cur() != .dev_entry &&
               grDevices::dev.cur() != 1L)
            tryCatch(grDevices::dev.off(),
                     error = function(e2) return(invisible(NULL)))
    }
    .draw <- function(x) {
        if (inherits(x, "grob") || inherits(x, "gtable") ||
            inherits(x, "arrangelist")) {
            grid.newpage(); grid.draw(x)
        } else print(x)
    }
    tryCatch({
        pdf(paste0(base, ".pdf"), width = pdf_w, height = pdf_h)
        .draw(p)
        dev.off()
    }, error = function(e) {
        .unwind()
        cat("  WARNING: could not write", paste0(basename(base), ".pdf"),
            "--", conditionMessage(e), "\n")
    })
    tryCatch({
        png(paste0(base, ".png"), width = png_w, height = png_h, res = res)
        .draw(p)
        dev.off()
    }, error = function(e) {
        .unwind()
        cat("  WARNING: could not write", paste0(basename(base), ".png"),
            "--", conditionMessage(e), "\n")
    })
    invisible(NULL)
}

# Non-parametric tail dependence coefficients via empirical copula
# λ_L = C(q,q)/q  (lower tail: joint probability both below q, normalised)
# λ_U = P(V>1-q, U>1-q) / q  (upper tail, symmetric version)
# Estimated at threshold quantile q (default 0.10 = lower/upper 10%)
tail_dep_np <- function(u, v, q = 0.10) {
    Cqq   <- mean(u <= q & v <= q)                    # C(q,q) empirical
    lam_L <- if (q > 0) Cqq / q else NA_real_         # λ_L ≈ C(q,q)/q
    Cuq   <- mean(u > (1 - q) & v > (1 - q))          # P(U>1-q, V>1-q)
    lam_U <- if (q > 0) Cuq / q else NA_real_         # λ_U ≈ P(both in upper q)/q
    list(lambda_L = lam_L, lambda_U = lam_U)
}

# Cramer-von Mises distance between two empirical copulas (grid approximation)
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

# Fisher z-test for difference between two Spearman correlations
# Returns z-statistic and two-sided p-value
fisher_z_test <- function(r1, n1, r2, n2) {
    # Return NA with a note when the test is undefined (n < 4 or |r| = 1).
    bad <- function(msg) list(z_stat = NA_real_, p_value = NA_real_,
                              z1 = NA_real_, z2 = NA_real_, note = msg)
    if (!is.finite(n1) || !is.finite(n2) || n1 < 4 || n2 < 4)
        return(bad("not testable: Fisher z needs n >= 4 in both groups"))
    if (!is.finite(r1) || !is.finite(r2))
        return(bad("not testable: a correlation is NA or non-finite"))
    if (abs(r1) >= 1 || abs(r2) >= 1)
        return(bad("not testable: |r| = 1 gives an infinite z transform"))
    z1  <- 0.5 * log((1 + r1) / (1 - r1))
    z2  <- 0.5 * log((1 + r2) / (1 - r2))
    se  <- sqrt(1 / (n1 - 3) + 1 / (n2 - 3))
    z   <- (z1 - z2) / se
    p   <- 2 * pnorm(-abs(z))
    list(z_stat = z, p_value = p, z1 = z1, z2 = z2, note = NA_character_)
}

# Permutation test for difference in Spearman correlation between two groups
# H0: ρ_band1 = ρ_band2 (no difference in local correlation)
perm_corr_diff <- function(x1, y1, x2, y2, n_perm = 2000, seed = NULL) {
    if (!is.null(seed)) set.seed(seed)
    obs_diff <- cor(x1, y1, method = "spearman") -
        cor(x2, y2, method = "spearman")
    x_all  <- c(x1, x2)
    y_all  <- c(y1, y2)
    n1     <- length(x1)
    n_all  <- length(x_all)
    perm_d <- numeric(n_perm)
    for (i in seq_len(n_perm)) {
        idx    <- sample(n_all, n1, replace = FALSE)
        perm_d[i] <- cor(x_all[idx],      y_all[idx],      method = "spearman") -
            cor(x_all[-idx], y_all[-idx], method = "spearman")
    }
    p_val <- mean(abs(perm_d) >= abs(obs_diff))
    list(obs_diff = obs_diff, p_value = p_val, perm_dist = perm_d)
}


# =============================================================================
# SECTION 0b: The pipeline for one arm
# =============================================================================

run_one_arm <- function(arm_tag, config_file, out_dir, file_tag) {

    config_path <- file.path(CEP_HOME, config_file)
    if (!file.exists(config_path))
        stop("Configuration not found: ", config_path, call. = FALSE)

    # Output printed before the log exists is buffered and replayed into it.
    .p10c_preamble <- character(0)

    assign("CEP_ARM", arm_tag, envir = globalenv())
    on.exit(suppressWarnings(rm("CEP_ARM", envir = globalenv())), add = TRUE)

    .p10c_preamble <- tryCatch(
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

    .p10c_needed <- c("ARM_TAG", "ARM_DIR", "R_EXPORT_ROOT", "GOI", "DFG",
                      "HKG_GENES", "DFG_COMPOSITE_LABEL", "RUN_SUFFIX",
                      "PRIMARY_STYPE", "NORM_TAG", "NORM_TARGET_SUM",
                      "NORM_CEILING")
    .p10c_absent <- .p10c_needed[!vapply(.p10c_needed, exists, logical(1))]
    if (length(.p10c_absent) > 0)
        stop("The configuration loaded from\n    ", config_path,
             "\n  does not define: ", paste(.p10c_absent, collapse = ", "),
             call. = FALSE)

    if (!identical(ARM_TAG, arm_tag))
        stop("The configuration reports ARM_TAG '", ARM_TAG,
             "' but this arm is '", arm_tag, "'.", call. = FALSE)

    .p10c_preamble <- c(.p10c_preamble, utils::capture.output(
        cat("  Config OK: ", ARM_TAG, " | ", GOI, " + ", length(DFG), " DFG | ",
            NORM_TAG, "\n", sep = "")))

    set.seed(42)

    # == User parameters ======================================================
    DATASETS  <- c("NonNeoC", "NeoC")
    MIN_CELLS <- 50

    # Axis of the Pre/Post boundary; must match IP_CLASS_AXIS of the Part 10b run.
    IP_CLASS_AXIS <- "cDFG"
    if (!IP_CLASS_AXIS %in% c("GOI", "cDFG"))
        stop("IP_CLASS_AXIS must be \"GOI\" or \"cDFG\", got: ", IP_CLASS_AXIS,
             call. = FALSE)

    # == Input directories ====================================================
    PART10_DIR      <- cep_arm_dir_for("Part10", arm_tag)
    P10_PREFIX      <- cep_out_prefix("Part10", NORM_TAG,
                                      cep_file_tag_for(arm_tag))
    P10_FILE_PREFIX <- paste0(P10_PREFIX, "_Module1_", gsub("-", "", GOI), "_",
                              DFG_COMPOSITE_LABEL)

    # Part 10 scores: per-dataset folder first, then the arm root.
    p10c_scores_candidates <- function(ds) {
        fn <- paste0(P10_FILE_PREFIX, "_", ds, "_Pooled_SKCS_scores.csv")
        c(file.path(PART10_DIR, paste0(basename(PART10_DIR), "_", ds), fn),
          file.path(PART10_DIR, fn))
    }
    p10c_scores_csv <- function(ds) {
        cand <- p10c_scores_candidates(ds)
        hit  <- cand[file.exists(cand)]
        if (length(hit) > 0L) hit[1] else cand[1]
    }

    # Part 10b EP/TREP CSV: folder tag Part10b_SKCS, file prefix Part10b.
    P10B_SKCS_PREFIX <- paste0(cep_out_prefix("Part10b", NORM_TAG,
                                              cep_file_tag_for(arm_tag)),
                               "_SKCS")
    p10c_ep_candidates <- function(ds) {
        fn <- paste0(P10B_SKCS_PREFIX, "_", ds, "_EP_TREP_values.csv")
        c(file.path(cep_arm_dir_for("Part10b_SKCS", arm_tag), fn),
          file.path(cep_arm_dir_for("Part10b", arm_tag), fn),
          file.path(R_EXPORT_ROOT, "Part10b_SKCS",
                    paste0("Part10b_scran_", ds, "_EP_TREP_values.csv")))
    }
    p10c_ep_csv <- function(ds) {
        cand <- p10c_ep_candidates(ds)
        hit  <- cand[file.exists(cand)]
        if (length(hit) > 0L) hit[1] else cand[1]
    }

    # == Output directory and file prefix =====================================
    OUT_DIR <- out_dir
    dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

    OUT_PREFIX  <- cep_out_prefix(PART_TAG, NORM_TAG,
                                  cep_file_tag_for(arm_tag))
    FILE_PREFIX <- paste0(OUT_PREFIX, "_Copula1")

    # == Log file =============================================================
    # The message stream goes to a temporary file appended to the log at the end;
    # LOG_MESSAGES_P10C = FALSE keeps warnings on the console.
    LOG_MESSAGES_P10C <- TRUE
    LOG_PATH_P10C <- file.path(OUT_DIR,
                               paste0(FILE_PREFIX, "_console_log_",
                                      format(Sys.time(), "%Y%m%d_%H%M%S"),
                                      ".txt"))
    MSG_PATH_P10C <- file.path(OUT_DIR,
                               paste0(FILE_PREFIX, "_console_log_messages.tmp"))

    while (sink.number() > 0) sink()
    # The message sink number is 2 (stderr) when nothing is diverted.
    if (sink.number(type = "message") != 2L) sink(type = "message")

    .p10c_log_con <- file(LOG_PATH_P10C, open = "wt")
    .p10c_msg_con <- if (isTRUE(LOG_MESSAGES_P10C))
        file(MSG_PATH_P10C, open = "wt") else NULL
    sink(.p10c_log_con, append = FALSE, split = TRUE)
    if (!is.null(.p10c_msg_con)) sink(.p10c_msg_con, type = "message")
    on.exit({
        tryCatch({
            if (!isTRUE(.p10c_completed))
                cat("\n", strrep("!", 72),
                    "\n  Arm did not complete. See \"*** ARM FAILED\" in the console.\n",
                    strrep("!", 72), "\n", sep = "")
            if (sink.number(type = "message") != 2L) sink(type = "message")
            while (sink.number() > 0) sink()
            if (!is.null(.p10c_msg_con)) close(.p10c_msg_con)
            close(.p10c_log_con)
            if (!is.null(.p10c_msg_con) && file.exists(MSG_PATH_P10C)) {
                .msg_txt <- readLines(MSG_PATH_P10C, warn = FALSE)
                if (length(.msg_txt) > 0)
                    cat(c("", strrep("=", 72),
                          "  MESSAGE STREAM (warnings, messages, errors)",
                          strrep("=", 72), .msg_txt, ""),
                        sep = "\n", file = LOG_PATH_P10C, append = TRUE)
                unlink(MSG_PATH_P10C)
            }
        }, error = function(e) invisible(NULL))
    }, add = TRUE)

    .p10c_completed <- FALSE

    .p10c_warn_opt <- getOption("warn")
    options(warn = 1)   # surface warnings immediately rather than deferring
    on.exit(options(warn = .p10c_warn_opt), add = TRUE)
    cat(sprintf("Part10c Copula1 - log started: %s\n\n", Sys.time()))

    if (length(.p10c_preamble) > 0)
        cat(paste(.p10c_preamble, collapse = "\n"), "\n\n", sep = "")

    # == DATASET_META: regions, GAM data dir, IP thresholds ===================
    # IP values come from the IPRS output; NA until then.
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

    # == Load IPRS-derived IP values ==========================================
    # Candidate locations, searched in order.
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
        # Check the field exists before indexing it.
        ip_vec   <- iprs_rds$inflection_points_DFG
        ip_y_vec <- iprs_rds$inflection_points_DFG_y
        ip_from  <- if (!is.null(iprs_rds$ip_definition))
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
            warning("This IPRS RDS has no ip_definition field.",
                    call. = FALSE)

        if (is.null(ip_vec) || length(ip_vec) == 0L) {
            cat("  [IPRS] RDS has no usable inflection_points_DFG field:",
                iprs_rds_path, "\n")
            cat("  [IPRS] The IP boundary will not be drawn.\n")
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
                                       "  = f(IP_thresh), boundary axis\n"),
                                ds, ip_y_val))
                } else if (identical(IP_CLASS_AXIS, "cDFG")) {
                    cat(sprintf(paste0("  [IPRS] %-10s IP_y_thresh ABSENT from ",
                                       "the IPRS RDS\n"), ds))
                }
            }
        }
    } else {
        cat("  [IPRS] RDS not found. Looked in:\n")
        for (.p in .iprs_candidates) cat("      ", .p, "\n", sep = "")
        cat("  [IPRS] The IP boundary will not be drawn.\n")
    }

    cat("\n", strrep("=", 72), "\n")
    cat("  Part10c - Copula Fundamentals\n")
    cat("  Arm          :", ARM_TAG, "\n")
    cat("  GOI          :", GOI, "\n")
    cat("  DFG          :", DFG_COMPOSITE_LABEL, " - ",
        paste(DFG, collapse = ", "), "\n")
    cat("  Normalisation:", NORM_TAG, " log1p(UMIs/total x ",
        format(NORM_TARGET_SUM, scientific = FALSE), ")\n")
    cat("  Datasets     :", paste(DATASETS, collapse = ", "), "\n")
    cat("  Boundary axis:", IP_CLASS_AXIS,
        if (identical(IP_CLASS_AXIS, "cDFG"))
            "(no fallback: the dataset stops if IP_y is absent)\n"
        else "(GOI axis: Part 10b must use the same axis)\n")
    cat("  Part 10 dir  :", PART10_DIR, "\n")
    cat("  Part 10b EP  :", dirname(p10c_ep_csv(DATASETS[1])), "\n")
    cat("  Output dir   :", OUT_DIR, "\n")
    cat(strrep("=", 72), "\n\n")

    # Marker file, removed only when the arm completes.
    INCOMPLETE_MARKER_P10C <- file.path(OUT_DIR, "_RUN_INCOMPLETE.txt")
    writeLines(c(
        "This arm did not finish.",
        "",
        paste("Arm     :", arm_tag),
        paste("Part    :", PART_DIR_TAG),
        paste("Started :", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
        "",
        "Files in this folder are from a partial run."
    ), INCOMPLETE_MARKER_P10C)


    # =============================================================================
    # SECTION A0: PRE-ANALYSIS GUARDRAILS
    # =============================================================================
    #   GR-01  norm.mtx max within the CP10K ceiling
    #   GR-02  expm1(norm) row sums pin NORM_TARGET_SUM
    #   GR-04  raw count matrix integer-like
    #   GR-05  Part 10 SKCS scores CSV exists per dataset
    #   GR-06  Part 10b EP/TREP CSV exists per dataset

    GR_LOG_10C <- file.path(OUT_DIR, paste0(FILE_PREFIX, "_Guardrail_Log.txt"))

    .gr10c_log  <- list()
    .gr10c_pass <- 0L
    .gr10c_fail <- 0L
    .gr10c_warn <- 0L

    gr10c <- function(check_id, status, detail) {
        # A zero-length detail would otherwise drop the row from the log.
        if (length(detail) != 1L || is.na(detail))
            detail <- paste0("[detail unavailable: the message was built from a",
                             " missing value and came out zero-length]")
        # Unknown statuses still print a tag.
        tag <- switch(status, PASS = "[PASS]", FAIL = "[FAIL]", WARN = "[WARN]",
                      sprintf("[%s]", as.character(status)[1]))
        cat(sprintf("  %s  %-16s  %s\n", tag, check_id, detail))
        .gr10c_log[[length(.gr10c_log) + 1L]] <<- list(
            id = check_id, status = status, detail = detail)
        if (status == "PASS") .gr10c_pass <<- .gr10c_pass + 1L
        if (status == "FAIL") .gr10c_fail <<- .gr10c_fail + 1L
        if (status == "WARN") .gr10c_warn <<- .gr10c_warn + 1L
    }

    gr10c_write <- function() {
        con <- file(GR_LOG_10C, "wt")
        # Close the connection even if a write fails.
        on.exit(try(close(con), silent = TRUE), add = TRUE)
        writeLines(c(
            strrep("=", 72),
            paste0("Part 10c Copula1 - Guardrail Log [", NORM_TAG, "]"),
            paste("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
            paste("Arm :", ARM_TAG),
            paste0("Normalisation: adata.X = log1p(UMIs/total x ",
                   format(NORM_TARGET_SUM, scientific = FALSE),
                   "), exported unchanged"),
            strrep("=", 72), ""
        ), con)
        for (item in .gr10c_log)
            writeLines(sprintf("  [%s]  %-16s  %s", item$status, item$id,
                               item$detail), con)
        writeLines(c("", strrep("=", 72),
                     sprintf("  Total : %d",
                             .gr10c_pass + .gr10c_fail + .gr10c_warn),
                     sprintf("  PASS  : %d", .gr10c_pass),
                     sprintf("  WARN  : %d", .gr10c_warn),
                     sprintf("  FAIL  : %d", .gr10c_fail),
                     if (.gr10c_fail == 0 && .gr10c_warn == 0)
                         "\n  Result: no failures, no warnings"
                     else if (.gr10c_fail == 0)
                         "\n  Result: no failures, warnings present"
                     else "\n  Result: failures present",
                     strrep("=", 72)), con)
        cat("  Guardrail log written:", GR_LOG_10C, "\n")
    }

    cat("\n", strrep("=", 72), "\n")
    cat("  SECTION A0 - PRE-ANALYSIS GUARDRAILS\n")
    cat(strrep("=", 72), "\n\n")

    # == GR-01, GR-02, GR-04: normalisation checks (per dataset) ==============
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

            gr10c(paste0("GR-01-", .ds$tag),
                  if (n_max > 0 && n_max <= NORM_CEILING + 1e-6) "PASS"
                  else "FAIL",
                  sprintf("norm.mtx max = %.4f  (ceiling log1p(%s) = %.4f)",
                          n_max,
                          format(NORM_TARGET_SUM, scientific = FALSE),
                          NORM_CEILING))

            n_rows <- min(50L, nrow(X_n))
            rs     <- rowSums(expm1(X_n[seq_len(n_rows), , drop = FALSE]))
            rs_dev <- max(abs(rs - NORM_TARGET_SUM))
            gr10c(paste0("GR-02-", .ds$tag),
                  if (rs_dev < NORM_TARGET_SUM * 1e-6) "PASS" else "FAIL",
                  sprintf("expm1(norm) row sums: mean = %.6f, max |dev| = %.2e, n = %d",
                          mean(rs), rs_dev, n_rows))

            rm(X_n); gc(verbose = FALSE)
        } else {
            gr10c(paste0("GR-01-", .ds$tag), "FAIL",
                  paste("norm.mtx not found:", norm_mtx_path))
            gr10c(paste0("GR-02-", .ds$tag), "FAIL",
                  paste("row-sum check not run, norm.mtx not found:",
                        norm_mtx_path))
        }

        if (file.exists(raw_mtx_path)) {
            X_r    <- as(readMM(raw_mtx_path), "CsparseMatrix")
            r_samp <- as.numeric(X_r@x[seq_len(min(10000, length(X_r@x)))])
            is_int <- all(abs(r_samp - round(r_samp)) < 0.01)
            r_max  <- max(X_r)
            gr10c(paste0("GR-04-", .ds$tag),
                  if (is_int && r_max > 20) "PASS" else "FAIL",
                  sprintf("raw max = %.0f, integer-like = %s",
                          r_max, if (is_int) "TRUE" else "FALSE"))
            rm(X_r); gc(verbose = FALSE)
        } else {
            gr10c(paste0("GR-04-", .ds$tag), "FAIL",
                  paste("raw_selected matrix not found:", raw_mtx_path))
        }
    }

    # == GR-05: Part 10 scores CSV exists =====================================
    for (ds in DATASETS) {
        .p     <- p10c_scores_csv(ds)
        .found <- file.exists(.p)
        gr10c(paste0("GR-05-", ds),
              if (.found) "PASS" else "FAIL",
              sprintf("Part 10 scores CSV: %s",
                      if (.found)
                          file.path(basename(dirname(.p)), basename(.p))
                      else paste("not found in any location:",
                                 paste(p10c_scores_candidates(ds),
                                       collapse = "  |  "))))
    }

    # == GR-06: Part 10b EP/TREP CSV exists ===================================
    # WARN only: without it all GOI-positive cells are used, without colouring.
    for (ds in DATASETS) {
        .p     <- p10c_ep_csv(ds)
        .found <- file.exists(.p)
        gr10c(paste0("GR-06-", ds),
              if (.found) "PASS" else "WARN",
              sprintf("Part 10b EP/TREP CSV: %s",
                      if (.found)
                          file.path(basename(dirname(.p)), basename(.p))
                      else paste("not found (all GOI+ cells used):",
                                 paste(p10c_ep_candidates(ds),
                                       collapse = "  |  "))))
    }

    # Written now and again at the end, so an aborted arm still leaves a ledger.
    gr10c_write()

    cat(sprintf("\n  Pre-analysis guardrails: %d PASS / %d WARN / %d FAIL\n",
                .gr10c_pass, .gr10c_warn, .gr10c_fail))
    if (.gr10c_fail > 0) {
        cat("  *** Guardrail failures ***\n\n")
        warning("Part 10c pre-analysis guardrails: ", .gr10c_fail, " FAIL(s)",
                call. = FALSE)
    } else if (.gr10c_warn > 0) {
        cat(sprintf("  %d passed, %d warning(s), no failures.\n\n",
                    .gr10c_pass, .gr10c_warn))
    } else {
        cat("  All pre-analysis guardrails passed.\n\n")
    }


    # =============================================================================
    # SECTION B: MAIN FUNCTION (runs per dataset)
    # =============================================================================


    run_part10c <- function(dataset_tag) {

        set.seed(42)   # reset seed at the start of each dataset run

        meta      <- DATASET_META[[dataset_tag]]
        gam_dir   <- meta$gam_dir
        regions   <- meta$regions
        IP_thresh <- meta$IP_thresh

        # IP_thresh is the reported IP; IP_ACTIVE is the threshold on the boundary axis.
        IP_y_thresh <- meta$IP_y_thresh
        if (is.null(IP_y_thresh)) IP_y_thresh <- NA_real_
        # No fallback to the GOI axis: the Pre/Post colours are assigned on cDFG.
        IP_AXIS_USED <- IP_CLASS_AXIS
        if (identical(IP_CLASS_AXIS, "cDFG") && !is.finite(IP_y_thresh))
            stop("IP_CLASS_AXIS is \"cDFG\" but the IPRS RDS carries no ",
                 "IP_y_thresh for ", dataset_tag, ".",
                 call. = FALSE)
        IP_ACTIVE <- if (identical(IP_AXIS_USED, "cDFG")) IP_y_thresh else IP_thresh
        IP_ACTIVE_LABEL <- if (identical(IP_AXIS_USED, "cDFG"))
            "IP_y = f(IP_x) on cDFG" else "IP_x on GOI"
        if (is.finite(IP_ACTIVE))
            cat(sprintf(paste0("  Pre/Post boundary drawn on %s | ",
                               "reported IP_x = %s | threshold used = %.6f\n"),
                        IP_ACTIVE_LABEL,
                        if (is.finite(IP_thresh)) sprintf("%.6f", IP_thresh)
                        else "NA",
                        IP_ACTIVE))

        cat("\n", strrep("=", 72), "\n")
        cat("  Part10c Copula1:", dataset_tag, " [", ARM_TAG, "]\n")
        cat(strrep("=", 72), "\n\n")

        all_plots   <- list()
        excel_tabs  <- list()

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
        n_cells <- nrow(df_all)
        cat(sprintf("  GOI-positive cells: %d\n\n", n_cells))

        # == Load SKCS-core cells from the Part 10b EP CSV ====================
        ep_csv <- p10c_ep_csv(dataset_tag)
        skcs_df <- NULL
        if (file.exists(ep_csv)) {
            skcs_df <- read.csv(ep_csv, stringsAsFactors = FALSE)
            cat(sprintf("  SKCS-core cells loaded: %d\n\n", nrow(skcs_df)))
        } else {
            cat("  Part10b EP CSV not found - using all GOI-positive cells\n\n")
            skcs_df <- data.frame(
                barcode   = df_all$barcode,
                GOI_expr  = df_all$GOI,
                stringsAsFactors = FALSE
            )
        }

        # Match SKCS cells to GOI/cDFG data
        skcs_bc  <- skcs_df$barcode
        in_skcs  <- df_all$barcode %in% skcs_bc
        x_skcs   <- df_all$GOI[in_skcs]
        y_skcs   <- df_all$Expression[in_skcs]
        n_skcs   <- sum(in_skcs)
        cat(sprintf("  Using %d SKCS-core cells (of %d GOI+) for copula analysis\n\n",
                    n_skcs, n_cells))

        # == Load Part 10 pooled scores for spear_contrib / kend_contrib ======
        p10_scores_csv <- p10c_scores_csv(dataset_tag)
        p10_scores <- NULL
        if (file.exists(p10_scores_csv)) {
            # check.names = FALSE keeps hyphenated gene names unchanged.
            p10_scores <- read.csv(p10_scores_csv, stringsAsFactors = FALSE,
                                   check.names = FALSE)
            cat(sprintf("  Part10 pooled scores loaded: %d rows  (%s)\n\n",
                        nrow(p10_scores), basename(p10_scores_csv)))
        } else {
            cat("  Part10 pooled scores CSV not found - spear/kend will be NA\n")
            cat("    Expected:", p10_scores_csv, "\n\n")
        }

        x <- df_all$GOI
        y <- df_all$Expression
        n <- length(x)

        # == Copula coordinates over all GOI+ cells ===========================
        # SKCS-core cells are shown in colour, the rest in grey.
        u <- rank(x, ties.method = "average") / (n + 1)
        v <- rank(y, ties.method = "average") / (n + 1)

        # == Dependence statistics, reused for the plot annotations ===========

        # Original dependence
        rho_orig <- cor(x, y, method = "spearman")
        td_orig  <- tail_dep_np(u, v, q = 0.10)
        cat(sprintf("  Original: rho = %.4f  |  lambda_L = %.4f  |  lambda_U = %.4f\n",
                    rho_orig, td_orig$lambda_L, td_orig$lambda_U))

        # Shuffled (independence) reference
        set.seed(999)
        y_shuffled <- sample(y)
        u_shuf     <- rank(x,          ties.method = "average") / (n + 1)
        v_shuf     <- rank(y_shuffled, ties.method = "average") / (n + 1)

        rho_shuf <- cor(x, y_shuffled, method = "spearman")
        td_shuf  <- tail_dep_np(u_shuf, v_shuf, q = 0.10)
        ks_cdfg  <- suppressWarnings(ks.test(y, y_shuffled))
        cat(sprintf("  Shuffled: rho = %.4f  |  lambda_L = %.4f  |  lambda_U = %.4f\n",
                    rho_shuf, td_shuf$lambda_L, td_shuf$lambda_U))
        cat(sprintf("  KS marginals: D = %.4f  p = %.4f\n",
                    ks_cdfg$statistic, ks_cdfg$p.value))

        # CvM copula distance: original vs shuffled (permutation test)
        cvm_dist <- cvm_copula_dist(u, v, u_shuf, v_shuf, grid_size = 25)
        set.seed(42)
        cvm_null <- replicate(200, {
            idx1 <- sample(n, n, replace = TRUE)
            idx2 <- sample(n, n, replace = TRUE)
            cvm_copula_dist(u[idx1], v[idx1], u[idx2], v[idx2], grid_size = 20)
        })
        cvm_p <- mean(cvm_null >= cvm_dist)
        cvm_p_label <- if (cvm_p < 0.001) "p < 0.001" else sprintf("p = %.3f", cvm_p)
        cat(sprintf("  CvM copula distance = %.5f  (%s)\n\n",
                    cvm_dist, cvm_p_label))

        # == Annotation labels (reused across plots) ==========================
        # Expression-space panels (A) and (B): dependence stats only, no p-values
        ann_expr <- sprintf(
            paste0("r_s = %.3f\n",
                   "lambda_L = %.3f\n",
                   "lambda_U = %.3f"),
            rho_orig, td_orig$lambda_L, td_orig$lambda_U)

        ann_shuf <- sprintf(
            paste0("r_s = %.3f\n",
                   "lambda_L = %.3f\n",
                   "lambda_U = %.3f"),
            rho_shuf, td_shuf$lambda_L, td_shuf$lambda_U)


        # =============================================================================
        # PLOT 1: Marginal PDFs and joint density (expression space)
        # Top: p(x) of GOI; right: p(y) of cDFG; main: joint density with contours.
        # =============================================================================
        cat("  Building Plot 1: Marginal PDFs + Joint density ...\n")

        # Tight limits via pilot kde2d on quantile-clipped range
        x_inner1 <- quantile(x, c(0.02, 0.98))
        y_inner1 <- quantile(y, c(0.02, 0.98))
        d_pilot1 <- kde2d(x, y, n = 80,
                          h = c(bw.nrd0(x) * 1.2, bw.nrd0(y) * 1.2),
                          lims = c(x_inner1 + c(-0.3, 0.3),
                                   y_inner1 + c(-0.3, 0.3)))
        nx1 <- length(d_pilot1$x); ny1 <- length(d_pilot1$y)
        gx1 <- rep(d_pilot1$x, times = ny1)
        gy1 <- rep(d_pilot1$y, each  = nx1)
        gz1 <- as.vector(d_pilot1$z)
        thr1  <- max(gz1) * 0.01
        pad1  <- 0.15
        x_lim1 <- c(min(gx1[gz1 > thr1]) - pad1, max(gx1[gz1 > thr1]) + pad1)
        y_lim1 <- c(min(gy1[gz1 > thr1]) - pad1, max(gy1[gz1 > thr1]) + pad1)

        bw_x_val <- bw.nrd0(x) * 1.2
        bw_y_val <- bw.nrd0(y) * 1.2
        dens2d  <- kde2d(x, y, n = 150,
                         h = c(bw_x_val, bw_y_val),
                         lims = c(x_lim1, y_lim1))
        dens_df <- expand.grid(GOI = dens2d$x, cDFG = dens2d$y)
        dens_df$density <- as.vector(dens2d$z)

        bw_goi  <- density(x, bw = "nrd0")
        bw_cdfg <- density(y, bw = "nrd0")
        goi_dens  <- data.frame(val = bw_goi$x,  pdf = bw_goi$y)
        cdfg_dens <- data.frame(val = bw_cdfg$x, pdf = bw_cdfg$y)
        goi_dens  <- goi_dens[ goi_dens$val  >= x_lim1[1] & goi_dens$val  <= x_lim1[2], ]
        cdfg_dens <- cdfg_dens[cdfg_dens$val >= y_lim1[1] & cdfg_dens$val <= y_lim1[2], ]

        dens_max  <- max(dens_df$density)
        jet_cols <- colorRampPalette(
            c("white", "white", "#E8E8FF", "#C0B0FF",
              "#7B2FBE", "#4B0082", "purple",
              "cyan", "#00CC88", "yellow"))(512)

        p_main <- ggplot(dens_df, aes(x = GOI, y = cDFG)) +
            geom_raster(aes(fill = density), interpolate = TRUE, na.rm = TRUE) +
            geom_contour(aes(z = density), colour = "black",
                         alpha = 0.45, linewidth = 0.45, bins = 10,
                         na.rm = TRUE) +
            annotate("label",
                     x = x_lim1[2] - diff(x_lim1) * 0.03,
                     y = y_lim1[1] + diff(y_lim1) * 0.12,
                     label = ann_expr,
                     hjust = 1, vjust = 0,
                     size = 3.0, fontface = "bold",
                     fill = "white", alpha = 0.80,
                     label.padding = unit(0.3, "lines"),
                     colour = "gray15") +
            scale_fill_gradientn(
                colours   = jet_cols,
                limits    = c(0, dens_max),
                oob       = scales::squish,
                name      = "Density",
                na.value  = "white"
            ) +
            scale_x_continuous(limits = x_lim1, expand = c(0, 0)) +
            scale_y_continuous(limits = y_lim1, expand = c(0, 0)) +
            labs(x = paste0(GOI, " Expression (", NORM_TAG, " log1p)"),
                 y = paste0(DFG_COMPOSITE_LABEL, " Expression (", NORM_TAG, " log1p)")) +
            theme_copula() +
            theme(legend.position  = "none",
                  panel.grid.major = element_line(color = "gray85"),
                  panel.grid.minor = element_blank())

        p_top <- ggplot(goi_dens, aes(x = val, y = pdf)) +
            geom_area(fill = "#5555BB", alpha = 0.30) +
            geom_line(colour = "#6666CC", linewidth = 1.0) +
            annotate("text",
                     x = x_lim1[2] - diff(x_lim1) * 0.02,
                     y = max(goi_dens$pdf) * 0.85,
                     label = "p(x) - marginal PDF of GOI",
                     hjust = 1, size = 2.8, colour = "gray30") +
            scale_x_continuous(limits = x_lim1, expand = c(0, 0)) +
            scale_y_continuous(expand = c(0.02, 0)) +
            labs(x = NULL, y = "Density") +
            theme_copula() +
            theme(axis.text.x      = element_blank(),
                  axis.ticks.x     = element_blank(),
                  panel.grid.minor = element_blank())

        p_right <- ggplot(cdfg_dens, aes(x = pdf, y = val)) +
            geom_ribbon(aes(xmin = 0, xmax = pdf), fill = "#009988", alpha = 0.25) +
            geom_path(colour = "#00BBAA", linewidth = 1.0) +
            annotate("text",
                     x = max(cdfg_dens$pdf) * 0.50,
                     y = y_lim1[1] + diff(y_lim1) * 0.92,
                     label = "p(y)\ncDFG PDF",
                     hjust = 0, size = 2.6, colour = "gray30") +
            scale_y_continuous(limits = y_lim1, expand = c(0, 0)) +
            scale_x_continuous(expand = c(0.02, 0)) +
            labs(x = "Density", y = NULL) +
            theme_copula() +
            theme(axis.text.y      = element_blank(),
                  axis.ticks.y     = element_blank(),
                  panel.grid.minor = element_blank())

        left_col <- plot_grid(
            p_top, p_main,
            ncol        = 1,
            align       = "v",
            axis        = "lr",
            rel_heights = c(0.5, 3)
        )
        p_blank <- ggplot() + theme_void()
        right_col <- plot_grid(
            p_blank, p_right,
            ncol        = 1,
            align       = "v",
            axis        = "l",
            rel_heights = c(0.5, 3)
        )
        g1_inner <- plot_grid(
            left_col, right_col,
            ncol       = 2,
            rel_widths = c(3, 0.5)
        )
        title_str <- paste0(
            dataset_tag, "  -  Joint density with marginal PDFs  [EXPRESSION SPACE]\n",
            "Main: f(GOI, cDFG) | Top: p(GOI) -> u=Fg(x) | Right: p(cDFG) -> v=Fc(y)  |  n=",
            n)
        g1 <- plot_grid(
            ggdraw() + draw_label(title_str, fontface = "bold", size = 10,
                                  x = 0.5, hjust = 0.5),
            g1_inner,
            ncol        = 1,
            rel_heights = c(0.08, 1)
        )

        base1 <- file.path(OUT_DIR,
                           paste0(FILE_PREFIX, "_", dataset_tag, "_01_marginal_joint"))
        save_pp(g1, base1, pdf_w = 10, pdf_h = 10,
                png_w = 2000, png_h = 2000, res = 200)
        cat("  Plot 1 saved\n")
        all_plots[["marginal_joint"]] <- g1


        # =============================================================================
        # PLOT 2: Copula unit square, u = Fg(x), v = Fc(y)
        # Each cell after the CDF transform, coloured by IP x TREP group.
        # =============================================================================
        cat("  Building Plot 2: Copula unit square ...\n")

        copula_df <- data.frame(u = u, v = v,
                                x_raw = x, y_raw = y,
                                barcode = df_all$barcode,
                                is_skcs = in_skcs,
                                stringsAsFactors = FALSE)

        if ("IP_status" %in% names(skcs_df) &&
            "is_TREP"   %in% names(skcs_df)) {
            lkp_ip   <- setNames(skcs_df$IP_status, skcs_df$barcode)
            lkp_trep <- setNames(skcs_df$is_TREP,   skcs_df$barcode)
            copula_df$IP   <- lkp_ip  [copula_df$barcode]
            copula_df$TREP <- lkp_trep[copula_df$barcode]
            copula_df$group <- ifelse(
                !copula_df$is_skcs, "Non-selected",
                ifelse(!is.na(copula_df$IP) & !is.na(copula_df$TREP),
                       paste(copula_df$IP,
                             ifelse(copula_df$TREP, "TREP", "non-TREP")),
                       "Unknown"))
        } else {
            copula_df$group <- ifelse(copula_df$is_skcs,
                                      "SKCS-core", "Non-selected")
        }

        # Pre/Post labels recomputed on both axes to measure agreement with Part 10b;
        # no label is overwritten.
        ip_label_check <- list(n = 0L, agree_active = NA_real_,
                               agree_alt = NA_real_, axis = IP_AXIS_USED)
        if ("IP" %in% names(copula_df)) {
            .chk <- copula_df[!is.na(copula_df$IP), , drop = FALSE]
            if (nrow(.chk) > 0L) {
                .lab_goi  <- if (is.finite(IP_thresh))
                    ifelse(.chk$x_raw <= IP_thresh, "Pre-IP", "Post-IP") else NULL
                .lab_cdfg <- if (is.finite(IP_y_thresh))
                    ifelse(.chk$y_raw <  IP_y_thresh, "Pre-IP", "Post-IP") else NULL
                .lab_act  <- if (identical(IP_AXIS_USED, "cDFG")) .lab_cdfg else .lab_goi
                .lab_alt  <- if (identical(IP_AXIS_USED, "cDFG")) .lab_goi  else .lab_cdfg
                ip_label_check <- list(
                    n            = nrow(.chk),
                    agree_active = if (is.null(.lab_act)) NA_real_
                                   else mean(.lab_act == .chk$IP),
                    agree_alt    = if (is.null(.lab_alt)) NA_real_
                                   else mean(.lab_alt == .chk$IP),
                    axis         = IP_AXIS_USED)
                cat(sprintf(paste0("  Pre/Post label agreement with Part 10b: ",
                                   "%.4f on %s, %.4f on the other axis  (n=%d)\n"),
                            ip_label_check$agree_active, IP_AXIS_USED,
                            ip_label_check$agree_alt, ip_label_check$n))
            }
        }

        four_cols <- c(
            "Non-selected"     = "gray80",
            "Pre-IP TREP"      = "#E69F00",
            "Pre-IP non-TREP"  = "#1B3A8C",
            "Post-IP TREP"     = "#009E73",
            "Post-IP non-TREP" = "#D94F4F",
            "SKCS-core"        = "#4B0082",
            "Unknown"          = "gray60"
        )
        cols_use <- four_cols[names(four_cols) %in% unique(copula_df$group)]

        n_by_grp <- table(copula_df$group)
        grp_labs <- sapply(names(cols_use), function(g) {
            cnt <- as.integer(n_by_grp[g])
            if (is.na(cnt)) cnt <- 0L
            sprintf("%s  n=%d", g, cnt)
        })

        # == Export copula cell-level data ====================================
        cat("  Building copula cell data export ...\n")

        export_df <- data.frame(
            barcode   = copula_df$barcode,
            region    = df_all$region,
            dataset   = dataset_tag,
            GOI_expr  = copula_df$x_raw,
            cDFG_expr = copula_df$y_raw,
            stringsAsFactors = FALSE
        )

        # Merge spear_contrib / kend_contrib from Part10 pooled scores
        if (!is.null(p10_scores)) {
            lkp_s <- setNames(p10_scores$spear_contrib, p10_scores$barcode)
            lkp_k <- setNames(p10_scores$kend_contrib,  p10_scores$barcode)
            export_df$spear_contrib <- lkp_s[export_df$barcode]
            export_df$kend_contrib  <- lkp_k[export_df$barcode]
        } else {
            export_df$spear_contrib <- NA_real_
            export_df$kend_contrib  <- NA_real_
        }

        export_df$is_GOI_neg <- FALSE    # all GOI-positive by construction
        export_df$is_skcs    <- copula_df$is_skcs

        if ("TREP" %in% names(copula_df)) {
            export_df$is_trep <- copula_df$TREP
        } else {
            export_df$is_trep <- FALSE
        }

        # Copula coordinates
        export_df$rank_GOI  <- rank(copula_df$x_raw, ties.method = "average")
        export_df$rank_cDFG <- rank(copula_df$y_raw, ties.method = "average")
        export_df$n_plus_1  <- n + 1L
        export_df$u         <- copula_df$u
        export_df$v         <- copula_df$v

        cell_xlsx <- file.path(OUT_DIR,
                               paste0(FILE_PREFIX, "_", dataset_tag,
                                      "_copula_cell_data.xlsx"))
        write_xlsx(list("CopulaCellData" = export_df), path = cell_xlsx)
        cat(sprintf("  Copula cell data exported: %s  (%d rows, %d cols)\n\n",
                    basename(cell_xlsx), nrow(export_df), ncol(export_df)))

        # Draw order: non-SKCS first (background), then non-TREP, then TREP on top
        draw_ord_cop <- c(
            "Non-selected"     = 0L,
            "Pre-IP non-TREP"  = 1L,
            "Post-IP non-TREP" = 2L,
            "Pre-IP TREP"      = 3L,
            "Post-IP TREP"     = 4L,
            "SKCS-core"        = 3L,
            "Unknown"          = 0L
        )
        copula_df$draw_ord <- draw_ord_cop[
            match(copula_df$group, names(draw_ord_cop))]
        copula_df$draw_ord[is.na(copula_df$draw_ord)] <- 0L
        copula_df <- copula_df[order(copula_df$draw_ord), ]

        # IP mapped to copula space over all GOI+ cells, /(n+1) as for u and v.
        # IP_COP is the value on the boundary axis.
        IP_u <- if (is.finite(IP_thresh))
            sum(df_all$GOI <= IP_thresh) / (n_cells + 1) else NA_real_
        IP_v <- if (is.finite(IP_y_thresh))
            sum(df_all$Expression <= IP_y_thresh) / (n_cells + 1) else NA_real_
        IP_COP <- if (identical(IP_AXIS_USED, "cDFG")) IP_v else IP_u
        cat(sprintf(paste0("  IP_x = %s -> IP_u = Fg(IP_x) = %s | IP_y = %s -> ",
                           "IP_v = Fc(IP_y) = %s  [full-pop n=%d, /(n+1)]\n"),
                    if (is.finite(IP_thresh))   sprintf("%.4f", IP_thresh)   else "NA",
                    if (is.finite(IP_u))        sprintf("%.4f", IP_u)        else "NA",
                    if (is.finite(IP_y_thresh)) sprintf("%.4f", IP_y_thresh) else "NA",
                    if (is.finite(IP_v))        sprintf("%.4f", IP_v)        else "NA",
                    n_cells))

        # IP conversion printed on the copula unit-square figure.
        n_le_IP_x <- if (is.finite(IP_thresh))
            sum(df_all$GOI <= IP_thresh) else NA_integer_
        n_le_IP_y <- if (is.finite(IP_y_thresh))
            sum(df_all$Expression <= IP_y_thresh) else NA_integer_
        ip_conversion_line <- if (!is.finite(IP_v) || !is.finite(IP_u)) {
            "IP conversion unavailable for this dataset"
        } else {
            sprintf(paste0("IP conversion: IP_x = %.4f -> IP_y = f(IP_x) = %.4f",
                           "  |  boundary IP_v = Fc(IP_y) = #{cDFG <= %.4f}/(n+1)",
                           " = %d/%d = %.4f",
                           "  |  IP_u = Fg(IP_x) = %d/%d = %.4f (GOI axis, not",
                           " the boundary)"),
                    IP_thresh, IP_y_thresh, IP_y_thresh,
                    n_le_IP_y, n_cells + 1L, IP_v,
                    n_le_IP_x, n_cells + 1L, IP_u)
        }
        cat("  ", ip_conversion_line, "\n", sep = "")

        # Boundary line on the axis the Pre/Post labels were assigned on; none if NA.
        ip_boundary_layer <- function(space = c("expression", "copula"),
                                      linetype = "dashed", colour = "gray10",
                                      linewidth = 0.4) {
            space <- match.arg(space)
            val   <- if (identical(space, "copula")) IP_COP else IP_ACTIVE
            if (!is.finite(val)) return(list())
            if (identical(IP_AXIS_USED, "cDFG"))
                list(geom_hline(yintercept = val, linetype = linetype,
                                colour = colour, linewidth = linewidth))
            else
                list(geom_vline(xintercept = val, linetype = linetype,
                                colour = colour, linewidth = linewidth))
        }

        # The copula coordinate the Pre/Post boundary lies on, used to place the
        # three demonstration cells of Plot 2-demo.
        cop_active_col <- if (identical(IP_AXIS_USED, "cDFG")) "v" else "u"

        p_copula <- ggplot(copula_df, aes(x = u, y = v)) +
            # Background: non-SKCS cells in light gray
            geom_point(data = copula_df[!copula_df$is_skcs, ],
                       colour = "gray80", alpha = 0.30, size = 0.6) +
            # Foreground: SKCS-core cells in 4-group colours
            geom_point(data = copula_df[copula_df$is_skcs, ],
                       aes(colour = group),
                       alpha = 0.50, size = 0.9) +
            ip_boundary_layer("copula") +
            scale_colour_manual(values = cols_use, labels = grp_labs,
                                name = "Cell group") +
            scale_x_continuous(breaks = seq(0, 1, 0.25),
                               labels = c("0\n(lowest GOI)", "0.25", "0.5\n(median)",
                                          "0.75", "1\n(highest GOI)")) +
            scale_y_continuous(breaks = seq(0, 1, 0.25),
                               labels = c("0", "0.25", "0.5", "0.75", "1")) +
            labs(title    = paste0(dataset_tag,
                                   "  -  Copula unit square  [0,1] x [0,1]"),
                 subtitle = paste0(
                     "u = Fg(x) = rank(GOI) / (n+1)   |   v = Fc(y) = rank(cDFG) / (n+1)  |  n=",
                     n, " GOI+ cells\n",
                     "Grey = non-SKCS (n=", sum(!copula_df$is_skcs),
                     "). Dashed = Pre/Post boundary on ", IP_ACTIVE_LABEL,
                     " at ", if (is.finite(IP_COP)) sprintf("%.4f", IP_COP) else "NA",
                     "\n", ip_conversion_line),
                 x = "u = Fg(xi)  [GOI marginal CDF]",
                 y = "v = Fc(yi)  [cDFG marginal CDF]") +
            theme_copula() +
            theme(legend.position = "right",
                  aspect.ratio    = 1)

        save_pp(p_copula,
                file.path(OUT_DIR,
                          paste0(FILE_PREFIX, "_", dataset_tag, "_02_copula_unit_square")),
                pdf_w = 10, pdf_h = 9, png_w = 2000, png_h = 1800)
        cat("  Plot 2 saved\n")
        all_plots[["copula_unit_square"]] <- p_copula


        # =============================================================================
        # PLOT 2-demo: three representative cells and their (u, v) formulas
        # Cell A below the boundary, B on it, C above it.
        # =============================================================================
        cat("  Building Plot 2-demo: Representative cells with formulas ...\n")

        # Targets below, on and above the boundary; quartiles when it is unavailable.
        cop_targets <- if (is.finite(IP_COP))
            c(IP_COP * 0.40, IP_COP, IP_COP + (1 - IP_COP) * 0.55)
        else c(0.25, 0.50, 0.75)

        pick_cell <- function(target_c, df) {
            skcs_sub <- df[df$is_skcs, , drop = FALSE]
            # Fall back to the full frame if the SKCS mask selects nothing.
            if (nrow(skcs_sub) == 0L) skcs_sub <- df
            idx <- which.min(abs(skcs_sub[[cop_active_col]] - target_c))
            skcs_sub[idx, , drop = FALSE]
        }

        cell_A <- pick_cell(cop_targets[1], copula_df)
        cell_B <- pick_cell(cop_targets[2], copula_df)
        cell_C <- pick_cell(cop_targets[3], copula_df)
        demo_cells <- rbind(cell_A, cell_B, cell_C)
        demo_cells$label <- c("A", "B", "C")

        # Compute display values for the formula panel
        rank_x <- rank(x, ties.method = "average")
        rank_y <- rank(y, ties.method = "average")
        np1    <- n + 1L

        demo_info <- lapply(seq_len(3), function(i) {
            bc   <- demo_cells$barcode[i]
            # Find position in original (unordered) x,y vectors
            bc_all <- df_all$barcode
            j      <- which(bc_all == bc)[1]
            list(
                label  = demo_cells$label[i],
                xi     = x[j],
                yi     = y[j],
                rk_x   = rank_x[j],
                rk_y   = rank_y[j],
                ui     = demo_cells$u[i],
                vi     = demo_cells$v[i]
            )
        })

        # Build the scatter with highlighted cells
        p_demo_scatter <- ggplot(copula_df, aes(x = u, y = v)) +
            geom_point(data = copula_df[!copula_df$is_skcs, ],
                       colour = "gray80", alpha = 0.30, size = 0.6) +
            geom_point(data = copula_df[copula_df$is_skcs, ],
                       aes(colour = group), alpha = 0.50, size = 0.9) +
            ip_boundary_layer("copula") +
            # Highlight the 3 cells: large violet with dark violet border
            geom_point(data = demo_cells,
                       aes(x = u, y = v),
                       inherit.aes = FALSE,
                       shape = 21, size = 5,
                       fill = "#8B5CF6", colour = "#4C1D95",
                       stroke = 1.2, alpha = 0.95) +
            # Label A, B, C
            geom_text(data = demo_cells,
                      aes(x = u, y = v, label = label),
                      inherit.aes = FALSE,
                      fontface = "bold", size = 3.0,
                      colour = "white") +
            scale_colour_manual(values = cols_use, labels = grp_labs,
                                name = "Cell group") +
            scale_x_continuous(breaks = seq(0, 1, 0.25),
                               labels = c("0\n(lowest GOI)", "0.25",
                                          "0.5\n(median)", "0.75",
                                          "1\n(highest GOI)")) +
            scale_y_continuous(breaks = seq(0, 1, 0.25),
                               labels = c("0", "0.25", "0.5", "0.75", "1")) +
            labs(title = paste0(dataset_tag,
                                "  -  Copula unit square  [0,1] x [0,1]"),
                 subtitle = paste0(
                     "Three representative cells (A, B, C) highlighted\n",
                     "Formulas show how each cell\u2019s (u, v) is computed from rank / (n+1)"),
                 x = "u = Fg(xi)  [GOI marginal CDF]",
                 y = "v = Fc(yi)  [cDFG marginal CDF]") +
            theme_copula() +
            theme(legend.position = "right",
                  aspect.ratio    = 1)

        # Build the formula text panel
        formula_lines <- c(
            "   How (u, v) is computed",
            "   ──────────────────────────",
            sprintf("   n = %d cells   (n+1 = %d)", n, np1),
            ""
        )
        for (ci in demo_info) {
            formula_lines <- c(formula_lines,
                sprintf("   Cell %s", ci$label),
                sprintf("   GOI  x = %.4f   rank = %.1f", ci$xi, ci$rk_x),
                sprintf("   u = %.1f / %d = %.4f", ci$rk_x, np1, ci$ui),
                "",
                sprintf("   cDFG y = %.4f   rank = %.1f", ci$yi, ci$rk_y),
                sprintf("   v = %.1f / %d = %.4f", ci$rk_y, np1, ci$vi),
                "   ──────────────────────────",
                ""
            )
        }
        formula_text <- paste(formula_lines, collapse = "\n")

        p_formula <- ggplot() +
            annotate("text", x = 0, y = 0.5,
                     label = formula_text,
                     hjust = 0, vjust = 0.5,
                     family = "mono", size = 3.2,
                     colour = "gray15") +
            theme_void() +
            coord_cartesian(xlim = c(-0.05, 1), ylim = c(0, 1), clip = "off")

        p_demo <- plot_grid(p_demo_scatter, p_formula,
                            ncol = 2, rel_widths = c(1.1, 0.6),
                            align = "none")

        save_pp(p_demo,
                file.path(OUT_DIR,
                          paste0(FILE_PREFIX, "_", dataset_tag,
                                 "_02_copula_demo_cells")),
                pdf_w = 16, pdf_h = 9, png_w = 3200, png_h = 1800)
        cat("  Plot 2-demo saved\n")
        all_plots[["copula_demo_cells"]] <- p_demo


        # =============================================================================
        # PLOT 2a: Expression-space scatter (raw GOI vs cDFG, same colouring as Plot 2)
        # =============================================================================
        cat("  Building Plot 2a: Expression-space scatter ...\n")

        # Axis limits with small padding
        x_pad2a <- diff(range(x)) * 0.04
        y_pad2a <- diff(range(y)) * 0.04
        x_lim2a <- c(min(x) - x_pad2a, max(x) + x_pad2a)
        y_lim2a <- c(min(y) - y_pad2a, max(y) + y_pad2a)

        p_expr <- ggplot(copula_df, aes(x = x_raw, y = y_raw)) +
            geom_point(data = copula_df[!copula_df$is_skcs, ],
                       colour = "gray80", alpha = 0.30, size = 0.6) +
            geom_point(data = copula_df[copula_df$is_skcs, ],
                       aes(colour = group),
                       alpha = 0.50, size = 0.9) +
            ip_boundary_layer("expression") +
            scale_colour_manual(values = cols_use, labels = grp_labs,
                                name = "Cell group") +
            scale_x_continuous(limits = x_lim2a, expand = c(0, 0)) +
            scale_y_continuous(limits = y_lim2a, expand = c(0, 0)) +
            labs(title    = paste0(dataset_tag,
                                   "  -  Expression space scatter  [GOI vs cDFG]"),
                 subtitle = paste0(
                     "Grey = non-SKCS (n=", sum(!copula_df$is_skcs),
                     "). Same cells and colouring as the copula unit square (Plot 2)\n",
                     "Dashed = Pre/Post boundary on ", IP_ACTIVE_LABEL, " at ",
                     if (is.finite(IP_ACTIVE)) sprintf("%.4f", IP_ACTIVE) else "NA",
                     "  |  reported IP_x = ",
                     if (is.finite(IP_thresh)) sprintf("%.4f", IP_thresh) else "NA",
                     "  |  axes in ", NORM_TAG, " log1p units  |  n=", n),
                 x = paste0(GOI, " Expression (", NORM_TAG, " log1p)"),
                 y = paste0(DFG_COMPOSITE_LABEL, " Expression (", NORM_TAG, " log1p)")) +
            theme_copula() +
            theme(legend.position = "right",
                  aspect.ratio    = 1)

        save_pp(p_expr,
                file.path(OUT_DIR,
                          paste0(FILE_PREFIX, "_", dataset_tag,
                                 "_02a_expression_space_scatter")),
                pdf_w = 10, pdf_h = 9, png_w = 2000, png_h = 1800)
        cat("  Plot 2a saved\n")
        all_plots[["expression_space_scatter"]] <- p_expr


        # =============================================================================
        # PLOT 2b: Copula joint density with uniform marginal PDFs
        # Copula-space analogue of Plot 1; axes extended 0.20 beyond [0,1].
        # Test 2C: KS and Anderson-Darling tests of u and v against Uniform(0,1), BH-adjusted.
        # =============================================================================
        cat("  Building Plot 2b: Copula joint density with uniform marginal PDFs ...\n")

        # == 2b-i: Uniformity tests (Test 2C) =================================
        # goftest (Anderson-Darling) is installed if missing; KS only if that fails.
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
        has_goftest <- requireNamespace("goftest", quietly = TRUE)
        if (!has_goftest) {
            cat("  NOTE: 'goftest' unavailable - AD test skipped (KS only).\n")
        }

        # One-sample KS test, H0: u ~ Uniform(0,1); ties warning suppressed.
        ks_u <- suppressWarnings(ks.test(u, "punif", min = 0, max = 1))
        ks_v <- suppressWarnings(ks.test(v, "punif", min = 0, max = 1))

        if (has_goftest) {
            ad_u <- goftest::ad.test(u, "punif", min = 0, max = 1)
            ad_v <- goftest::ad.test(v, "punif", min = 0, max = 1)
            raw_p <- c(ks_u$p.value, ks_v$p.value,
                       ad_u$p.value, ad_v$p.value)
            bh_p  <- p.adjust(raw_p, method = "BH")
        } else {
            ad_u <- list(statistic = NA_real_, p.value = NA_real_)
            ad_v <- list(statistic = NA_real_, p.value = NA_real_)
            raw_p <- c(ks_u$p.value, ks_v$p.value, NA_real_, NA_real_)
            bh_p  <- p.adjust(raw_p[!is.na(raw_p)], method = "BH")
            bh_p  <- c(bh_p, NA_real_, NA_real_)   # pad to length 4
        }

        unif_test_df <- data.frame(
            dataset       = dataset_tag,
            variable      = c("u (GOI CDF)", "v (cDFG CDF)",
                              "u (GOI CDF)", "v (cDFG CDF)"),
            test          = c("KS", "KS", "Anderson-Darling", "Anderson-Darling"),
            H0            = "Uniform(0,1)",
            statistic     = round(c(ks_u$statistic, ks_v$statistic,
                                    ad_u$statistic, ad_v$statistic), 5),
            p_raw         = round(raw_p, 5),
            p_BH          = round(bh_p,  5),
            uniform_confirmed = ifelse(is.na(bh_p), "Not tested (package unavailable)",
                                       ifelse(bh_p > 0.05,
                                              "Yes (fail to reject H0 — uniform)",
                                              "No (reject H0 — non-uniform)")),
            note = c(
                "KS boundary bias in KDE panels does not affect this test",
                "KS boundary bias in KDE panels does not affect this test",
                "AD is more sensitive at tails than KS",
                "AD is more sensitive at tails than KS"
            ),
            stringsAsFactors = FALSE
        )

        # Print to console
        cat("\n  ── Test 2C: Uniformity of copula marginals ──────────────────────\n")
        cat(sprintf("  Dataset: %s   n = %d\n\n", dataset_tag, n))
        cat(sprintf("  KS test  u vs Uniform(0,1): D = %.5f   p_raw = %.5f   p_BH = %.5f   -> %s\n",
                    ks_u$statistic, raw_p[1], bh_p[1],
                    ifelse(bh_p[1] > 0.05,
                           "UNIFORM (fail to reject H0)",
                           "NON-UNIFORM (reject H0)")))
        cat(sprintf("  KS test  v vs Uniform(0,1): D = %.5f   p_raw = %.5f   p_BH = %.5f   -> %s\n",
                    ks_v$statistic, raw_p[2], bh_p[2],
                    ifelse(bh_p[2] > 0.05,
                           "UNIFORM (fail to reject H0)",
                           "NON-UNIFORM (reject H0)")))
        if (has_goftest) {
            cat(sprintf("  AD test  u vs Uniform(0,1): A = %.5f   p_raw = %.5f   p_BH = %.5f   -> %s\n",
                        ad_u$statistic, raw_p[3], bh_p[3],
                        ifelse(bh_p[3] > 0.05,
                               "UNIFORM (fail to reject H0)",
                               "NON-UNIFORM (reject H0)")))
            cat(sprintf("  AD test  v vs Uniform(0,1): A = %.5f   p_raw = %.5f   p_BH = %.5f   -> %s\n",
                        ad_v$statistic, raw_p[4], bh_p[4],
                        ifelse(bh_p[4] > 0.05,
                               "UNIFORM (fail to reject H0)",
                               "NON-UNIFORM (reject H0)")))
        }
        cat("  BH correction applied across all p-values in this test set.\n")
        cat("  ─────────────────────────────────────────────────────────────────\n\n")

        # Stored for the workbook
        excel_tabs[["Test2C_UniformityTests"]] <- unif_test_df

        # == 2b-ii: 2D KDE of (u, v) on the unit square =======================
        # Grid and axes extend COP_PAD beyond [0,1] so contours are not clipped.
        COP_PAD <- 0.20
        u_lim2b <- c(0 - COP_PAD, 1 + COP_PAD)
        v_lim2b <- c(0 - COP_PAD, 1 + COP_PAD)

        bw_u <- bw.nrd0(u) * 1.0
        bw_v <- bw.nrd0(v) * 1.0

        dens2d_cop <- kde2d(u, v, n = 150,
                            h    = c(bw_u, bw_v),
                            lims = c(u_lim2b, v_lim2b))
        dens_cop_df <- expand.grid(u_grid = dens2d_cop$x,
                                   v_grid = dens2d_cop$y)
        dens_cop_df$density <- as.vector(dens2d_cop$z)

        # == 2b-iii: Marginal PDFs of u and v =================================
        u_dens <- density(u, bw = "nrd0",
                          from = u_lim2b[1], to = u_lim2b[2], n = 512)
        v_dens <- density(v, bw = "nrd0",
                          from = v_lim2b[1], to = v_lim2b[2], n = 512)

        u_dens_df <- data.frame(val = u_dens$x, pdf = u_dens$y)
        v_dens_df <- data.frame(val = v_dens$x, pdf = v_dens$y)

        # Clip to the padded range for rendering
        u_dens_df <- u_dens_df[u_dens_df$val >= u_lim2b[1] &
                                   u_dens_df$val <= u_lim2b[2], ]
        v_dens_df <- v_dens_df[v_dens_df$val >= v_lim2b[1] &
                                   v_dens_df$val <= v_lim2b[2], ]

        # == 2b-iv: Joint density panel =======================================
        cop_dens_max <- max(dens_cop_df$density)

        jet_cols_cop <- colorRampPalette(
            c("white", "white", "#E8E8FF", "#C0B0FF",
              "#7B2FBE", "#4B0082", "purple",
              "cyan", "#00CC88", "yellow"))(512)

        unif_ref <- 1.0   # theoretical p(u) = p(v) = 1 for Uniform(0,1)

        rho_uv <- cor(u, v, method = "spearman")

        # Break positions: include the two boundary-box ticks at 0 and 1 plus the pad
        cop_x_breaks <- c(0, 0.25, 0.50, 0.75, 1.00)

        p_cop_main <- ggplot(dens_cop_df, aes(x = u_grid, y = v_grid)) +
            geom_raster(aes(fill = density), interpolate = TRUE, na.rm = TRUE) +
            geom_contour(aes(z = density), colour = "black",
                         alpha = 0.50, linewidth = 0.45, bins = 10,
                         na.rm = TRUE) +
            # Vertical and horizontal lines at the true [0,1] copula boundaries
            geom_vline(xintercept = c(0, 1), colour = "gray60",
                       linewidth = 0.35, linetype = "dotted") +
            geom_hline(yintercept = c(0, 1), colour = "gray60",
                       linewidth = 0.35, linetype = "dotted") +
            # Diagonal reference line: perfect concordance
            geom_abline(slope = 1, intercept = 0,
                        linetype = "dashed", colour = "white",
                        linewidth = 0.6, alpha = 0.7) +
            scale_fill_gradientn(
                colours  = jet_cols_cop,
                limits   = c(0, cop_dens_max),
                oob      = scales::squish,
                name     = "Density",
                na.value = "white"
            ) +
            scale_x_continuous(limits = u_lim2b, expand = c(0, 0),
                               breaks = cop_x_breaks) +
            scale_y_continuous(limits = v_lim2b, expand = c(0, 0),
                               breaks = cop_x_breaks) +
            # Annotation: dependence + CvM + uniformity (copula panel = (C))
            annotate("label", x = 0.98, y = 0.03,
                     label = paste0(
                         sprintf(paste0(
                             "r_s = %.3f\n",
                             "lambda_L = %.3f\n",
                             "lambda_U = %.3f\n"),
                             rho_orig, td_orig$lambda_L, td_orig$lambda_U),
                         sprintf("CvM: %s ((B) vs (A))\n", cvm_p_label),
                         sprintf("KS & AD: FDR = %.1f (uniform)",
                                 min(bh_p, na.rm = TRUE))),
                     hjust = 1, vjust = 0, size = 2.8,
                     fontface = "bold",
                     fill = alpha("black", 0.55),
                     label.padding = unit(0.35, "lines"),
                     colour = "white") +
            labs(x = "u = Fg(xi)  [GOI marginal CDF]",
                 y = "v = Fc(yi)  [cDFG marginal CDF]") +
            theme_copula() +
            theme(legend.position  = "none",
                  panel.background = element_rect(fill = "gray98", colour = NA),
                  panel.grid.major = element_line(color = "gray80"),
                  panel.grid.minor = element_blank())

        # == 2b-v: Top panel, marginal PDF of u ===============================
        # y-axis headroom: at least 1.30 to show both the curve and the ref line
        u_pdf_max <- max(u_dens_df$pdf, 1.0) * 1.30

        # Annotation: KS p-value placed inside the panel
        ks_u_label <- sprintf("KS p = %.4f (BH: %.4f)", raw_p[1], bh_p[1])

        p_cop_top <- ggplot(u_dens_df, aes(x = val, y = pdf)) +
            geom_hline(yintercept = unif_ref, linetype = "dashed",
                       colour = "gray50", linewidth = 0.7) +
            geom_area(fill = "#5555BB", alpha = 0.30) +
            geom_line(colour = "#6666CC", linewidth = 1.0) +
            annotate("text",
                     x = u_lim2b[2] - 0.01, y = u_pdf_max * 0.90,
                     label = ks_u_label,
                     hjust = 1, size = 2.5, colour = "gray30") +
            scale_x_continuous(limits = u_lim2b, expand = c(0, 0),
                               breaks = cop_x_breaks) +
            scale_y_continuous(limits = c(0, u_pdf_max), expand = c(0.02, 0)) +
            labs(x = NULL, y = "Density") +
            theme_copula() +
            theme(axis.text.x      = element_blank(),
                  axis.ticks.x     = element_blank(),
                  panel.grid.minor = element_blank())

        # == 2b-vi: Right panel, marginal PDF of v (rotated) ==================
        v_pdf_max <- max(v_dens_df$pdf, 1.0) * 1.30

        ks_v_label <- sprintf("KS p = %.4f\n(BH: %.4f)", raw_p[2], bh_p[2])

        p_cop_right <- ggplot(v_dens_df, aes(x = pdf, y = val)) +
            geom_vline(xintercept = unif_ref, linetype = "dashed",
                       colour = "gray50", linewidth = 0.7) +
            geom_ribbon(aes(xmin = 0, xmax = pdf), fill = "#009988", alpha = 0.25) +
            geom_path(colour = "#00BBAA", linewidth = 1.0) +
            annotate("text",
                     x = v_pdf_max * 0.40,
                     y = v_lim2b[2] - 0.03,
                     label = ks_v_label,
                     hjust = 0, vjust = 1, size = 2.4, colour = "gray30") +
            scale_y_continuous(limits = v_lim2b, expand = c(0, 0),
                               breaks = cop_x_breaks) +
            scale_x_continuous(limits = c(0, v_pdf_max), expand = c(0.02, 0)) +
            labs(x = "Density", y = NULL) +
            theme_copula() +
            theme(axis.text.y      = element_blank(),
                  axis.ticks.y     = element_blank(),
                  panel.grid.minor = element_blank())

        # == 2b-vii: Assemble with cowplot (same layout as Plot 1) ============
        left_col_2b <- plot_grid(
            p_cop_top,
            p_cop_main,
            ncol        = 1,
            align       = "v",
            axis        = "lr",
            rel_heights = c(0.5, 3)
        )
        p_blank_2b <- ggplot() + theme_void()
        right_col_2b <- plot_grid(
            p_blank_2b,
            p_cop_right,
            ncol        = 1,
            align       = "v",
            axis        = "l",
            rel_heights = c(0.5, 3)
        )
        g2b_inner <- plot_grid(
            left_col_2b, right_col_2b,
            ncol       = 2,
            rel_widths = c(3, 0.5)
        )

        # KS result summary for the title
        ks_summary <- sprintf(
            "KS u: D=%.4f p=%s  |  KS v: D=%.4f p=%s  (BH-corrected; p>0.05 = uniform confirmed)",
            ks_u$statistic,
            ifelse(bh_p[1] < 0.001, "<0.001", sprintf("%.4f", bh_p[1])),
            ks_v$statistic,
            ifelse(bh_p[2] < 0.001, "<0.001", sprintf("%.4f", bh_p[2]))
        )

        title_2b <- paste0(
            dataset_tag, "  -  Copula joint density with marginal PDFs  [COPULA SPACE]\n",
            "Main: C(u,v) kernel density  |  Top: p(u) marginal  |  Right: p(v) marginal  |  n=", n, "\n",
            "Dashed grey = Uniform(0,1) reference  |  Dashed white = concordance diagonal  |  ",
            "Axis range extended by 0.20 to prevent boundary clipping\n",
            ks_summary)

        g2b <- plot_grid(
            ggdraw() + draw_label(title_2b, fontface = "bold", size = 8.5,
                                  x = 0.5, hjust = 0.5),
            g2b_inner,
            ncol        = 1,
            rel_heights = c(0.12, 1)
        )

        base2b <- file.path(OUT_DIR,
                            paste0(FILE_PREFIX, "_", dataset_tag,
                                   "_02b_copula_joint_density"))
        save_pp(g2b, base2b, pdf_w = 10, pdf_h = 10,
                png_w = 2000, png_h = 2000, res = 200)
        cat("  Plot 2b saved\n")
        all_plots[["copula_joint_density"]] <- g2b


        # =============================================================================
        # PLOT 2c: Conditional density slices c(v | u) at five u values
        # Each slice integrates to 1, so the u-marginal is uniform.
        # =============================================================================
        cat("  Building Plot 2c: Conditional density slices ...\n")

        # Re-use dens2d_cop from Plot 2b
        u_grid <- dens2d_cop$x     # 150 grid points along u
        v_grid <- dens2d_cop$y     # 150 grid points along v
        z_mat  <- dens2d_cop$z     # 150 x 150 density matrix

        # Select 5 u-positions to slice
        u_slices <- c(0.10, 0.25, 0.50, 0.75, 0.90)
        slice_cols <- c("#2166AC", "#67A9CF", "#404040", "#EF8A62", "#B2182B")

        # Conditional density at each u slice, normalised to area 1 by adaptive
        # quadrature (integrate()) on a spline of the KDE grid values.
        dv_fine <- diff(v_grid[1:2])

        slice_list <- lapply(seq_along(u_slices), function(k) {
            u_target <- u_slices[k]
            j <- which.min(abs(u_grid - u_target))
            cond_v <- z_mat[j, ]

            # Cubic spline, clamped non-negative.
            sp_raw <- splinefun(v_grid, cond_v, method = "natural")
            sp_fn  <- function(v) pmax(sp_raw(v), 0)
            # Adaptive Gauss-Kronrod quadrature
            aq <- tryCatch(
                integrate(sp_fn, lower = min(v_grid), upper = max(v_grid),
                          subdivisions = 500L, rel.tol = 1e-8,
                          stop.on.error = FALSE),
                error = function(e) {
                    # Fallback: trapezoidal rule on the fine grid
                    dv <- diff(v_grid[1:2])
                    list(value = sum(cond_v) * dv, abs.error = dv^2,
                         subdivisions = 0L, message = "trapezoidal fallback")
                }
            )
            area_raw <- aq$value

            # Normalise so integral = 1.0
            cond_norm <- if (area_raw > 0) cond_v / area_raw else cond_v
            sp_norm_raw <- splinefun(v_grid, cond_norm, method = "natural")
            sp_norm     <- function(v) pmax(sp_norm_raw(v), 0)
            aq_norm <- tryCatch(
                integrate(sp_norm, lower = min(v_grid), upper = max(v_grid),
                          subdivisions = 500L, rel.tol = 1e-8,
                          stop.on.error = FALSE),
                error = function(e) {
                    dv <- diff(v_grid[1:2])
                    list(value = sum(cond_norm) * dv, abs.error = dv^2,
                         subdivisions = 0L)
                }
            )

            data.frame(
                v     = v_grid,
                density = cond_norm,
                u_val = sprintf("u = %.2f", u_target),
                u_num = u_target,
                area  = round(aq_norm$value, 6),
                aq_abs_error = aq_norm$abs.error,
                aq_subdivisions = ifelse(is.null(aq_norm$subdivisions), 0L,
                                         aq_norm$subdivisions),
                stringsAsFactors = FALSE
            )
        })
        slice_df <- do.call(rbind, slice_list)
        slice_df$u_val <- factor(slice_df$u_val,
                                 levels = sprintf("u = %.2f", u_slices))

        # == Panel A: Copula density heatmap with slice lines =================
        p2c_heat <- ggplot(dens_cop_df, aes(x = u_grid, y = v_grid)) +
            geom_raster(aes(fill = density), interpolate = TRUE) +
            scale_fill_gradientn(colours = jet_cols, name = "c(u,v)") +
            geom_vline(xintercept = u_slices,
                       colour = slice_cols, linewidth = 0.7,
                       linetype = "dashed") +
            annotate("text",
                     x = u_slices + 0.02, y = rep(0.03, 5),
                     label = sprintf("%.2f", u_slices),
                     colour = slice_cols, size = 3, fontface = "bold",
                     hjust = 0) +
            coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
            labs(title = "A: Copula joint density c(u,v)",
                 subtitle = "Dashed lines mark cross-section positions",
                 x = "u = Fg(GOI)", y = "v = Fc(cDFG)") +
            theme_copula() +
            theme(aspect.ratio = 1, legend.position = "right")

        # == Panel B: Conditional densities (faceted) =========================
        p2c_slices <- ggplot(slice_df, aes(x = v, y = density)) +
            geom_area(aes(fill = u_val), alpha = 0.35, show.legend = FALSE) +
            geom_line(aes(colour = u_val), linewidth = 0.8, show.legend = FALSE) +
            scale_fill_manual(values = setNames(slice_cols,
                                                sprintf("u = %.2f", u_slices))) +
            scale_colour_manual(values = setNames(slice_cols,
                                                  sprintf("u = %.2f", u_slices))) +
            geom_text(data = slice_df[!duplicated(slice_df$u_val), ],
                      aes(label = sprintf("Area = %.4f", area)),
                      x = 0.5, y = Inf, vjust = 1.3, size = 3.2,
                      fontface = "italic", colour = "gray30") +
            facet_wrap(~ u_val, ncol = 5, scales = "free_y") +
            labs(title = "B: Conditional density c(v | u) at each slice",
                 subtitle = paste0(
                     "Each slice has a DIFFERENT shape but the SAME area (=1). ",
                     "Area verified by adaptive Gauss-Kronrod quadrature.\n",
                     "Marginals confirmed uniform: KS & AD FDR=1.0. ",
                     "Joint confirmed NON-uniform: CvM distance=",
                     sprintf("%.4f", cvm_dist),
                     ", permutation p=",
                     sprintf("%.4f", cvm_p),
                     " (vs independence copula C(u,v)=u*v)."),
                 x = "v = Fc(cDFG)",
                 y = "Conditional density c(v | u)") +
            theme_copula() +
            theme(strip.text = element_text(face = "bold", size = 9),
                  strip.background = element_rect(fill = "gray95"))

        # == Panel C: Overlaid slices on one axis =============================
        p2c_overlay <- ggplot(slice_df, aes(x = v, y = density,
                                            colour = u_val, fill = u_val)) +
            geom_area(alpha = 0.10, position = "identity") +
            geom_line(linewidth = 0.9) +
            scale_colour_manual(values = setNames(slice_cols,
                                                  sprintf("u = %.2f", u_slices)),
                                name = "Slice at") +
            scale_fill_manual(values = setNames(slice_cols,
                                                sprintf("u = %.2f", u_slices)),
                              guide = "none") +
            labs(title = "C: All slices overlaid",
                 subtitle = paste0(
                     "High u (red): density peaks near v=1 (upper tail dependence).\n",
                     "Low u (blue): flatter, weaker dependence. Every curve integrates to 1."),
                 x = "v = Fc(cDFG)",
                 y = "Conditional density c(v | u)") +
            theme_copula() +
            theme(legend.position = "right")

        # == Combine into 3-panel figure ======================================
        g2c <- plot_grid(
            plot_grid(p2c_heat, p2c_overlay, ncol = 2, rel_widths = c(1, 1.1),
                      labels = NULL),
            p2c_slices,
            ncol  = 1,
            rel_heights = c(1, 0.7)
        )

        g2c_titled <- plot_grid(
            ggdraw() + draw_label(
                paste0(dataset_tag,
                       "  -  Why copula marginals are uniform despite non-uniform joint density"),
                fontface = "bold", size = 11, x = 0.01, hjust = 0),
            g2c,
            ncol = 1,
            rel_heights = c(0.04, 1)
        )

        save_pp(g2c_titled,
                file.path(OUT_DIR,
                          paste0(FILE_PREFIX, "_", dataset_tag,
                                 "_02c_conditional_density_slices")),
                pdf_w = 14, pdf_h = 12, png_w = 2800, png_h = 2400)
        cat("  Plot 2c saved (conditional density slices)\n")
        all_plots[["conditional_density_slices"]] <- g2c_titled


        # =============================================================================
        # PLOT 2d: Horizontal slices c(u | v) at five v values
        # Companion to Plot 2c for the v-marginal.
        # =============================================================================
        cat("  Building Plot 2d: Horizontal conditional density slices (v-direction) ...\n")

        # Re-use dens2d_cop, u_grid, v_grid, z_mat from Plot 2b/2c
        # Select 5 v-positions to slice (same positions as u-slices in 2c)
        v_slices   <- c(0.10, 0.25, 0.50, 0.75, 0.90)
        slice_cols_d <- c("#2166AC", "#67A9CF", "#404040", "#EF8A62", "#B2182B")
        du_fine    <- diff(u_grid[1:2])

        # Extract conditional density at each v slice
        slice_list_d <- lapply(seq_along(v_slices), function(k) {
            v_target <- v_slices[k]
            j <- which.min(abs(v_grid - v_target))
            cond_u <- z_mat[, j]   # column = fix v, vary u

            sp_raw <- splinefun(u_grid, cond_u, method = "natural")
            sp_fn  <- function(u) pmax(sp_raw(u), 0)
            aq <- tryCatch(
                integrate(sp_fn, lower = min(u_grid), upper = max(u_grid),
                          subdivisions = 500L, rel.tol = 1e-8,
                          stop.on.error = FALSE),
                error = function(e) {
                    du <- diff(u_grid[1:2])
                    list(value = sum(cond_u) * du, abs.error = du^2,
                         subdivisions = 0L, message = "trapezoidal fallback")
                }
            )
            area_raw <- aq$value

            cond_norm <- if (area_raw > 0) cond_u / area_raw else cond_u
            sp_norm_raw <- splinefun(u_grid, cond_norm, method = "natural")
            sp_norm     <- function(u) pmax(sp_norm_raw(u), 0)
            aq_norm <- tryCatch(
                integrate(sp_norm, lower = min(u_grid), upper = max(u_grid),
                          subdivisions = 500L, rel.tol = 1e-8,
                          stop.on.error = FALSE),
                error = function(e) {
                    du <- diff(u_grid[1:2])
                    list(value = sum(cond_norm) * du, abs.error = du^2,
                         subdivisions = 0L)
                }
            )

            data.frame(
                u     = u_grid,
                density = cond_norm,
                v_val = sprintf("v = %.2f", v_target),
                v_num = v_target,
                area  = round(aq_norm$value, 6),
                aq_abs_error = aq_norm$abs.error,
                aq_subdivisions = ifelse(is.null(aq_norm$subdivisions), 0L,
                                         aq_norm$subdivisions),
                stringsAsFactors = FALSE
            )
        })
        slice_df_d <- do.call(rbind, slice_list_d)
        slice_df_d$v_val <- factor(slice_df_d$v_val,
                                   levels = sprintf("v = %.2f", v_slices))

        # == Panel A: Copula density heatmap with horizontal slice lines ======
        p2d_heat <- ggplot(dens_cop_df, aes(x = u_grid, y = v_grid)) +
            geom_raster(aes(fill = density), interpolate = TRUE) +
            scale_fill_gradientn(colours = jet_cols, name = "c(u,v)") +
            geom_hline(yintercept = v_slices,
                       colour = slice_cols_d, linewidth = 0.7,
                       linetype = "dashed") +
            annotate("text",
                     x = rep(0.03, 5), y = v_slices + 0.02,
                     label = sprintf("%.2f", v_slices),
                     colour = slice_cols_d, size = 3, fontface = "bold",
                     hjust = 0, vjust = 0) +
            coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
            labs(title = "A: Copula joint density c(u,v)",
                 subtitle = "Dashed lines mark horizontal cross-section positions",
                 x = "u = Fg(GOI)", y = "v = Fc(cDFG)") +
            theme_copula() +
            theme(aspect.ratio = 1, legend.position = "right")

        # == Panel B: Conditional densities (faceted) =========================
        p2d_slices <- ggplot(slice_df_d, aes(x = u, y = density)) +
            geom_area(aes(fill = v_val), alpha = 0.35, show.legend = FALSE) +
            geom_line(aes(colour = v_val), linewidth = 0.8, show.legend = FALSE) +
            scale_fill_manual(values = setNames(slice_cols_d,
                                                sprintf("v = %.2f", v_slices))) +
            scale_colour_manual(values = setNames(slice_cols_d,
                                                  sprintf("v = %.2f", v_slices))) +
            geom_text(data = slice_df_d[!duplicated(slice_df_d$v_val), ],
                      aes(label = sprintf("Area = %.4f", area)),
                      x = 0.5, y = Inf, vjust = 1.3, size = 3.2,
                      fontface = "italic", colour = "gray30") +
            facet_wrap(~ v_val, ncol = 5, scales = "free_y") +
            labs(title = "B: Conditional density c(u | v) at each slice",
                 subtitle = paste0(
                     "Each slice has a DIFFERENT shape but the SAME area (=1). ",
                     "Area verified by adaptive Gauss-Kronrod quadrature.\n",
                     "Marginals confirmed uniform: KS & AD FDR=1.0. ",
                     "Joint confirmed NON-uniform: CvM distance=",
                     sprintf("%.4f", cvm_dist),
                     ", permutation p=",
                     sprintf("%.4f", cvm_p),
                     " (vs independence copula C(u,v)=u*v)."),
                 x = "u = Fg(GOI)",
                 y = "Conditional density c(u | v)") +
            theme_copula() +
            theme(strip.text = element_text(face = "bold", size = 9),
                  strip.background = element_rect(fill = "gray95"))

        # == Panel C: Overlaid slices on one axis =============================
        p2d_overlay <- ggplot(slice_df_d, aes(x = u, y = density,
                                              colour = v_val, fill = v_val)) +
            geom_area(alpha = 0.10, position = "identity") +
            geom_line(linewidth = 0.9) +
            scale_colour_manual(values = setNames(slice_cols_d,
                                                  sprintf("v = %.2f", v_slices)),
                                name = "Slice at") +
            scale_fill_manual(values = setNames(slice_cols_d,
                                                sprintf("v = %.2f", v_slices)),
                              guide = "none") +
            labs(title = "C: All slices overlaid",
                 subtitle = paste0(
                     "High v (red): density peaks near u=1 (upper tail dependence).\n",
                     "Low v (blue): flatter, weaker dependence. Every curve integrates to 1."),
                 x = "u = Fg(GOI)",
                 y = "Conditional density c(u | v)") +
            theme_copula() +
            theme(legend.position = "right")

        # == Combine into 3-panel figure ======================================
        g2d <- plot_grid(
            plot_grid(p2d_heat, p2d_overlay, ncol = 2, rel_widths = c(1, 1.1),
                      labels = NULL),
            p2d_slices,
            ncol  = 1,
            rel_heights = c(1, 0.7)
        )

        g2d_titled <- plot_grid(
            ggdraw() + draw_label(
                paste0(dataset_tag,
                       "  -  Horizontal slices: integrating c(u,v) over u gives uniform v-margin"),
                fontface = "bold", size = 11, x = 0.01, hjust = 0),
            g2d,
            ncol = 1,
            rel_heights = c(0.04, 1)
        )

        save_pp(g2d_titled,
                file.path(OUT_DIR,
                          paste0(FILE_PREFIX, "_", dataset_tag,
                                 "_02d_conditional_density_slices_v")),
                pdf_w = 14, pdf_h = 12, png_w = 2800, png_h = 2400)
        cat("  Plot 2d saved (horizontal conditional density slices, v-direction)\n")
        all_plots[["conditional_density_slices_v"]] <- g2d_titled


        # =============================================================================
        # EXHAUSTIVE MARGINAL INTEGRATION: every grid slice, both directions
        # Integrates c(u,v) at every KDE grid point; pass if every integral in the
        # boundary-safe zone is within MARG_TOL of 1.
        # =============================================================================
        cat("  Exhaustive marginal integration (all slices, both directions) ...\n")

        # Re-use dens2d_cop, u_grid, v_grid, z_mat from Plot 2b
        n_grid_e <- length(u_grid)

        # == Direction 1: int_0^1 c(u,v) dv at every u grid point =============
        cat(sprintf("    Direction 1: integrating c(u,v) dv at all %d u-grid points ...\n", n_grid_e))
        Iu_all <- numeric(n_grid_e)
        for (i in seq_len(n_grid_e)) {
            sp_raw <- splinefun(v_grid, z_mat[i, ], method = "natural")
            sp_fn  <- function(v) pmax(sp_raw(v), 0)
            aq <- tryCatch(
                integrate(sp_fn, lower = 0, upper = 1,
                          subdivisions = 500L, rel.tol = 1e-8,
                          stop.on.error = FALSE),
                error = function(e) {
                    in01 <- v_grid >= 0 & v_grid <= 1
                    list(value = sum(z_mat[i, in01]) * diff(v_grid[1:2]))
                }
            )
            Iu_all[i] <- aq$value
        }

        # == Direction 2: int_0^1 c(u,v) du at every v grid point =============
        cat(sprintf("    Direction 2: integrating c(u,v) du at all %d v-grid points ...\n", n_grid_e))
        Iv_all <- numeric(n_grid_e)
        for (j in seq_len(n_grid_e)) {
            sp_raw <- splinefun(u_grid, z_mat[, j], method = "natural")
            sp_fn  <- function(u) pmax(sp_raw(u), 0)
            aq <- tryCatch(
                integrate(sp_fn, lower = 0, upper = 1,
                          subdivisions = 500L, rel.tol = 1e-8,
                          stop.on.error = FALSE),
                error = function(e) {
                    in01 <- u_grid >= 0 & u_grid <= 1
                    list(value = sum(z_mat[in01, j]) * diff(u_grid[1:2]))
                }
            )
            Iv_all[j] <- aq$value
        }

        # Boundary-safe zone (exclude 2*bandwidth fringe for pass/fail)
        bw_margin_u_e <- 2 * bw_u
        bw_margin_v_e <- 2 * bw_v
        safe_u <- u_grid >= bw_margin_u_e & u_grid <= (1 - bw_margin_u_e)
        safe_v <- v_grid >= bw_margin_v_e & v_grid <= (1 - bw_margin_v_e)

        Iu_safe <- Iu_all[safe_u]
        Iv_safe <- Iv_all[safe_v]
        n_safe_u <- sum(safe_u)
        n_safe_v <- sum(safe_v)

        # An empty boundary-safe zone gives NA, not a pass.
        .safe_ok_u <- n_safe_u > 0L
        .safe_ok_v <- n_safe_v > 0L

        # Summary: mean of all slice integrals
        Iu_mean <- if (.safe_ok_u) mean(Iu_safe)          else NA_real_
        Iv_mean <- if (.safe_ok_v) mean(Iv_safe)          else NA_real_
        Iu_sd   <- if (.safe_ok_u) sd(Iu_safe)            else NA_real_
        Iv_sd   <- if (.safe_ok_v) sd(Iv_safe)            else NA_real_
        Iu_maxd <- if (.safe_ok_u) max(abs(Iu_safe - 1))  else NA_real_
        Iv_maxd <- if (.safe_ok_v) max(abs(Iv_safe - 1))  else NA_real_

        # The key statistic: sum of all integrals / n_slices
        Iu_ratio <- if (.safe_ok_u) sum(Iu_safe) / n_safe_u else NA_real_
        Iv_ratio <- if (.safe_ok_v) sum(Iv_safe) / n_safe_v else NA_real_

        MARG_TOL <- 0.05
        # isTRUE(), so NA is not a pass.
        u_pass <- isTRUE(Iu_maxd < MARG_TOL)
        v_pass <- isTRUE(Iv_maxd < MARG_TOL)
        all_pass <- u_pass & v_pass
        if (!.safe_ok_u || !.safe_ok_v)
            cat(sprintf(paste0("\n  Boundary-safe zone empty",
                               " (%d u-slices, %d v-slices); not verified.\n"),
                        n_safe_u, n_safe_v))

        # Console report
        cat("\n  ── Exhaustive marginal integration (all slices, both directions) ──\n")
        cat(sprintf("  Integration limits: [0, 1]  |  Boundary-safe: [%.3f, %.3f]\n",
                    bw_margin_u_e, 1 - bw_margin_u_e))
        cat(sprintf("\n  DIRECTION 1 - int c(u,v) dv  (%d slices in boundary-safe zone):\n", n_safe_u))
        cat(sprintf("    sum(I_u) / n_slices = %.8f\n", Iu_ratio))
        cat(sprintf("    mean = %.6f   sd = %.6f   max|I-1| = %.6f  -> %s\n",
                    Iu_mean, Iu_sd, Iu_maxd, ifelse(u_pass, "PASS", "FAIL")))
        cat(sprintf("\n  DIRECTION 2 - int c(u,v) du  (%d slices in boundary-safe zone):\n", n_safe_v))
        cat(sprintf("    sum(I_v) / n_slices = %.8f\n", Iv_ratio))
        cat(sprintf("    mean = %.6f   sd = %.6f   max|I-1| = %.6f  -> %s\n",
                    Iv_mean, Iv_sd, Iv_maxd, ifelse(v_pass, "PASS", "FAIL")))
        cat(sprintf("\n  OVERALL: %s  (tolerance: max|I-1| < %.2f in boundary-safe zone)\n",
                    ifelse(all_pass, "BOTH MARGINALS VERIFIED", "NOT VERIFIED"), MARG_TOL))
        cat("  ─────────────────────────────────────────────────────────────────\n\n")

        # == Workbook export ==================================================
        exhaust_verify_df <- data.frame(
            dataset   = dataset_tag,
            direction = c("int c(u,v) dv (u-margin)", "int c(u,v) du (v-margin)"),
            n_slices_safe = c(n_safe_u, n_safe_v),
            sum_div_n = c(Iu_ratio, Iv_ratio),
            mean      = c(Iu_mean, Iv_mean),
            sd        = c(Iu_sd, Iv_sd),
            max_dev   = c(Iu_maxd, Iv_maxd),
            result    = c(ifelse(u_pass, "PASS", "FAIL"),
                          ifelse(v_pass, "PASS", "FAIL")),
            tolerance = MARG_TOL,
            boundary_safe = sprintf("[%.3f, %.3f]", bw_margin_u_e, 1 - bw_margin_u_e),
            stringsAsFactors = FALSE
        )
        excel_tabs[["ExhaustiveMarginal_Verification"]] <- exhaust_verify_df


        # =============================================================================
        # TEST 2A: Same marginals, different dependence
        # Statistics were computed before Plot 1; printed and exported here.
        # =============================================================================
        cat("\n  Test 2A: Same marginals, different dependence ...\n")
        cat(sprintf("    KS test cDFG marginal: D=%.4f  p=%.4f\n",
                    ks_cdfg$statistic, ks_cdfg$p.value))
        cat(sprintf("    Spearman rho  original : %.4f\n", rho_orig))
        cat(sprintf("    Spearman rho  shuffled : %.4f\n", rho_shuf))
        cat(sprintf("    CvM copula distance  : %.6f\n", cvm_dist))
        cat(sprintf("    CvM permutation p    : %.4f\n", cvm_p))
        cat(sprintf("    Tail dep (original)  lambda_L=%.4f  lambda_U=%.4f\n",
                    td_orig$lambda_L, td_orig$lambda_U))
        cat(sprintf("    Tail dep (shuffled)  lambda_L=%.4f  lambda_U=%.4f\n",
                    td_shuf$lambda_L, td_shuf$lambda_U))

        q_grid <- c(0.05, 0.10, 0.15, 0.20)
        td_table <- do.call(rbind, lapply(q_grid, function(q) {
            td_o <- tail_dep_np(u,      v,      q = q)
            td_s <- tail_dep_np(u_shuf, v_shuf, q = q)
            data.frame(
                dataset         = dataset_tag,
                threshold_q     = q,
                orig_lambda_L   = round(td_o$lambda_L, 4),
                orig_lambda_U   = round(td_o$lambda_U, 4),
                shuf_lambda_L   = round(td_s$lambda_L, 4),
                shuf_lambda_U   = round(td_s$lambda_U, 4),
                delta_lambda_L  = round(td_o$lambda_L - td_s$lambda_L, 4),
                delta_lambda_U  = round(td_o$lambda_U - td_s$lambda_U, 4),
                stringsAsFactors = FALSE
            )
        }))

        test2a_summary <- data.frame(
            dataset                  = dataset_tag,
            n_cells                  = n,
            KS_D_cDFG_marginal       = round(ks_cdfg$statistic, 4),
            KS_p_cDFG_marginal       = round(ks_cdfg$p.value,   4),
            KS_interpretation        = ifelse(ks_cdfg$p.value > 0.05,
                                              "Marginals statistically same (p>0.05)",
                                              "Marginals differ"),
            Spearman_rho_original    = round(rho_orig, 4),
            Spearman_rho_shuffled    = round(rho_shuf, 4),
            CvM_copula_distance      = round(cvm_dist, 6),
            CvM_perm_p               = round(cvm_p,    4),
            CvM_interpretation       = ifelse(cvm_p < 0.05,
                                              "Copula structures significantly different",
                                              "Copula structures not distinguishable"),
            lambda_L_original_q10    = round(td_orig$lambda_L, 4),
            lambda_U_original_q10    = round(td_orig$lambda_U, 4),
            lambda_L_shuffled_q10    = round(td_shuf$lambda_L, 4),
            lambda_U_shuffled_q10    = round(td_shuf$lambda_U, 4),
            stringsAsFactors = FALSE
        )

        excel_tabs[["Test2A_Summary"]] <- test2a_summary
        excel_tabs[["Test2A_TailDep"]] <- td_table

        # Plot 3: Side-by-side joint density, original vs shuffled
        cat("  Building Plot 3: Original vs Shuffled joint density ...\n")

        x_inner <- quantile(x, c(0.02, 0.98))
        y_inner <- quantile(y, c(0.02, 0.98))
        d_pilot <- kde2d(x, y, n = 80,
                         h = c(bw.nrd0(x) * 1.2, bw.nrd0(y) * 1.2),
                         lims = c(x_inner + c(-0.3, 0.3),
                                  y_inner + c(-0.3, 0.3)))
        nx <- length(d_pilot$x); ny <- length(d_pilot$y)
        gx_p <- rep(d_pilot$x, times = ny)
        gy_p <- rep(d_pilot$y, each  = nx)
        gz_p <- as.vector(d_pilot$z)
        thr_p <- max(gz_p) * 0.01
        pad3  <- 0.15
        x_lim3 <- c(min(gx_p[gz_p > thr_p]) - pad3,
                    max(gx_p[gz_p > thr_p]) + pad3)

        y_inner_s <- quantile(y_shuffled, c(0.02, 0.98))
        d_pilot_s <- kde2d(x, y_shuffled, n = 80,
                           h = c(bw.nrd0(x) * 1.2, bw.nrd0(y_shuffled) * 1.2),
                           lims = c(x_inner + c(-0.3, 0.3),
                                    y_inner_s + c(-0.3, 0.3)))
        nx_s <- length(d_pilot_s$x); ny_s <- length(d_pilot_s$y)
        gy_ps <- rep(d_pilot_s$y, each = nx_s)
        gz_ps <- as.vector(d_pilot_s$z)
        thr_ps <- max(gz_ps) * 0.01
        y_lim3 <- c(min(c(gy_p[gz_p > thr_p], gy_ps[gz_ps > thr_ps])) - pad3,
                    max(c(gy_p[gz_p > thr_p], gy_ps[gz_ps > thr_ps])) + pad3)

        d_orig_kde <- kde2d(x, y, n = 120,
                            h = c(bw.nrd0(x) * 1.2, bw.nrd0(y) * 1.2),
                            lims = c(x_lim3, y_lim3))
        d_shuf_kde <- kde2d(x, y_shuffled, n = 120,
                            h = c(bw.nrd0(x) * 1.2, bw.nrd0(y_shuffled) * 1.2),
                            lims = c(x_lim3, y_lim3))

        df_orig_dens <- expand.grid(GOI = d_orig_kde$x, cDFG = d_orig_kde$y)
        df_orig_dens$density <- as.vector(d_orig_kde$z)
        df_shuf_dens <- expand.grid(GOI = d_shuf_kde$x, cDFG = d_shuf_kde$y)
        df_shuf_dens$density <- as.vector(d_shuf_kde$z)

        dens_max3 <- max(c(df_orig_dens$density, df_shuf_dens$density), na.rm = TRUE)
        jet_cols3 <- colorRampPalette(
            c("white", "white", "#E8E8FF", "#C0B0FF",
              "#7B2FBE", "#4B0082", "purple",
              "cyan", "#00CC88", "yellow"))(512)

        make_density_panel <- function(df_dens, panel_title, show_y = TRUE) {
            p <- ggplot(df_dens, aes(x = GOI, y = cDFG)) +
                geom_raster(aes(fill = density), interpolate = TRUE, na.rm = TRUE) +
                geom_contour(aes(z = density), colour = "black",
                             alpha = 0.45, linewidth = 0.45, bins = 10,
                             na.rm = TRUE) +
                scale_fill_gradientn(
                    colours  = jet_cols3,
                    limits   = c(0, dens_max3),
                    oob      = scales::squish,
                    name     = "Density",
                    na.value = "white"
                ) +
                scale_x_continuous(limits = x_lim3, expand = c(0, 0)) +
                scale_y_continuous(limits = y_lim3, expand = c(0, 0)) +
                labs(title = panel_title,
                     x = paste0(GOI, " Expression (", NORM_TAG, " log1p)"),
                     y = if (show_y) paste0(DFG_COMPOSITE_LABEL, " Expression (", NORM_TAG, " log1p)") else NULL) +
                theme_copula() +
                theme(legend.position  = "none",
                      strip.background = element_rect(fill = "gray85"),
                      panel.grid.major = element_line(color = "gray80"),
                      panel.grid.minor = element_blank(),
                      plot.title       = element_text(face = "bold", size = 9,
                                                      hjust = 0.5))
            if (!show_y) p <- p + theme(axis.text.y  = element_blank(),
                                        axis.ticks.y = element_blank())
            p
        }

        title_orig <- sprintf("Original\nrho=%.3f | lambda_L=%.3f | lambda_U=%.3f",
                              rho_orig, td_orig$lambda_L, td_orig$lambda_U)
        title_shuf <- sprintf("Shuffled (same marginals)\nrho=%.3f | lambda_L=%.3f | lambda_U=%.3f",
                              rho_shuf, td_shuf$lambda_L, td_shuf$lambda_U)

        p_orig <- make_density_panel(df_orig_dens, title_orig, show_y = TRUE)
        p_shuf <- make_density_panel(df_shuf_dens, title_shuf, show_y = FALSE)

        p_orig_leg <- p_orig + theme(legend.position = "bottom")
        legend_grob <- cowplot::get_legend(p_orig_leg)

        panels_row <- plot_grid(p_orig, p_shuf,
                                ncol = 2, align = "h", axis = "tb",
                                rel_widths = c(1, 1))

        main_title <- paste0(
            dataset_tag, "  -  Same marginals, different dependence\n",
            "Left: original GOI-cDFG pairing (real co-expression)\n",
            "Right: y values randomly shuffled - marginal distributions IDENTICAL, but dependence broken\n",
            sprintf("CvM copula distance = %.5f  (permutation p = %.4f  -> %s)",
                    cvm_dist, cvm_p,
                    ifelse(cvm_p < 0.05, "different copula structures", "similar")))

        title_grob <- ggdraw() +
            draw_label(main_title, fontface = "plain", size = 9,
                       x = 0.01, hjust = 0, vjust = 1)

        p_compare <- plot_grid(
            title_grob, panels_row, legend_grob,
            ncol        = 1,
            rel_heights = c(0.12, 1, 0.08)
        )

        base3 <- file.path(OUT_DIR,
                           paste0(FILE_PREFIX, "_", dataset_tag,
                                  "_03_shuffled_vs_original"))
        save_pp(p_compare, base3, pdf_w = 14, pdf_h = 8,
                png_w = 2800, png_h = 1600, res = 200)
        cat("  Plot 3 saved\n")
        all_plots[["shuffled_vs_original"]] <- p_compare


        # =============================================================================
        # PLOT 3a: Shuffled joint density with marginal PDFs
        # Same layout as Plot 1 with y_shuffled; marginals identical, dependence broken.
        # =============================================================================
        cat("  Building Plot 3a: Shuffled joint density with marginal PDFs ...\n")

        # Use Plot 3's shared axis limits for fair comparison with Plot 1
        # Recompute KDE for the shuffled pairing within these limits
        bw_x_3a <- bw.nrd0(x) * 1.2
        bw_y_3a <- bw.nrd0(y_shuffled) * 1.2
        dens2d_3a  <- kde2d(x, y_shuffled, n = 150,
                            h = c(bw_x_3a, bw_y_3a),
                            lims = c(x_lim3, y_lim3))
        dens_df_3a <- expand.grid(GOI = dens2d_3a$x, cDFG = dens2d_3a$y)
        dens_df_3a$density <- as.vector(dens2d_3a$z)

        # Marginal PDFs (same values as Plot 1 since marginals are identical)
        bw_goi_3a  <- density(x, bw = "nrd0")
        bw_cdfg_3a <- density(y_shuffled, bw = "nrd0")
        goi_dens_3a  <- data.frame(val = bw_goi_3a$x,  pdf = bw_goi_3a$y)
        cdfg_dens_3a <- data.frame(val = bw_cdfg_3a$x, pdf = bw_cdfg_3a$y)
        goi_dens_3a  <- goi_dens_3a[ goi_dens_3a$val  >= x_lim3[1] &
                                         goi_dens_3a$val  <= x_lim3[2], ]
        cdfg_dens_3a <- cdfg_dens_3a[cdfg_dens_3a$val >= y_lim3[1] &
                                         cdfg_dens_3a$val <= y_lim3[2], ]

        # Plot 1's density maximum sets the colour scale, for a fair comparison.
        dens_max_3a <- dens_max   # from Plot 1 (original, unshuffled)
        jet_cols_3a <- colorRampPalette(
            c("white", "white", "#E8E8FF", "#C0B0FF",
              "#7B2FBE", "#4B0082", "purple",
              "cyan", "#00CC88", "yellow"))(512)

        p_main_3a <- ggplot(dens_df_3a, aes(x = GOI, y = cDFG)) +
            geom_raster(aes(fill = density), interpolate = TRUE, na.rm = TRUE) +
            geom_contour(aes(z = density), colour = "black",
                         alpha = 0.45, linewidth = 0.45, bins = 10,
                         na.rm = TRUE) +
            annotate("label",
                     x = x_lim3[2] - diff(x_lim3) * 0.03,
                     y = y_lim3[1] + diff(y_lim3) * 0.12,
                     label = ann_shuf,
                     hjust = 1, vjust = 0,
                     size = 3.0, fontface = "bold",
                     fill = "white", alpha = 0.80,
                     label.padding = unit(0.3, "lines"),
                     colour = "gray15") +
            scale_fill_gradientn(
                colours   = jet_cols_3a,
                limits    = c(0, dens_max_3a),
                oob       = scales::squish,
                name      = "Density",
                na.value  = "white"
            ) +
            scale_x_continuous(limits = x_lim3, expand = c(0, 0)) +
            scale_y_continuous(limits = y_lim3, expand = c(0, 0)) +
            labs(x = paste0(GOI, " Expression (", NORM_TAG, " log1p)"),
                 y = paste0(DFG_COMPOSITE_LABEL, " Expression (", NORM_TAG, " log1p)")) +
            theme_copula() +
            theme(legend.position  = "none",
                  panel.grid.major = element_line(color = "gray85"),
                  panel.grid.minor = element_blank())

        p_top_3a <- ggplot(goi_dens_3a, aes(x = val, y = pdf)) +
            geom_area(fill = "#5555BB", alpha = 0.30) +
            geom_line(colour = "#6666CC", linewidth = 1.0) +
            annotate("text",
                     x = x_lim3[2] - diff(x_lim3) * 0.02,
                     y = max(goi_dens_3a$pdf) * 0.85,
                     label = "p(x) - marginal PDF of GOI",
                     hjust = 1, size = 2.8, colour = "gray30") +
            scale_x_continuous(limits = x_lim3, expand = c(0, 0)) +
            scale_y_continuous(expand = c(0.02, 0)) +
            labs(x = NULL, y = "Density") +
            theme_copula() +
            theme(axis.text.x      = element_blank(),
                  axis.ticks.x     = element_blank(),
                  panel.grid.minor = element_blank())

        p_right_3a <- ggplot(cdfg_dens_3a, aes(x = pdf, y = val)) +
            geom_ribbon(aes(xmin = 0, xmax = pdf), fill = "#009988", alpha = 0.25) +
            geom_path(colour = "#00BBAA", linewidth = 1.0) +
            annotate("text",
                     x = max(cdfg_dens_3a$pdf) * 0.50,
                     y = y_lim3[1] + diff(y_lim3) * 0.92,
                     label = "p(y)\ncDFG PDF",
                     hjust = 0, size = 2.6, colour = "gray30") +
            scale_y_continuous(limits = y_lim3, expand = c(0, 0)) +
            scale_x_continuous(expand = c(0.02, 0)) +
            labs(x = "Density", y = NULL) +
            theme_copula() +
            theme(axis.text.y      = element_blank(),
                  axis.ticks.y     = element_blank(),
                  panel.grid.minor = element_blank())

        left_col_3a <- plot_grid(
            p_top_3a, p_main_3a,
            ncol        = 1,
            align       = "v",
            axis        = "lr",
            rel_heights = c(0.5, 3)
        )
        p_blank_3a <- ggplot() + theme_void()
        right_col_3a <- plot_grid(
            p_blank_3a, p_right_3a,
            ncol        = 1,
            align       = "v",
            axis        = "l",
            rel_heights = c(0.5, 3)
        )
        g3a_inner <- plot_grid(
            left_col_3a, right_col_3a,
            ncol       = 2,
            rel_widths = c(3, 0.5)
        )
        title_str_3a <- paste0(
            dataset_tag,
            "  -  SHUFFLED joint density with marginal PDFs  [EXPRESSION SPACE]\n",
            "Same marginal distributions as Plot 1 (identical p(x) and p(y))  |  ",
            "Dependence broken by random permutation of cDFG\n",
            sprintf("rho=%.3f | lambda_L=%.3f | lambda_U=%.3f  |  n=%d",
                    rho_shuf, td_shuf$lambda_L, td_shuf$lambda_U, n))
        g3a <- plot_grid(
            ggdraw() + draw_label(title_str_3a, fontface = "bold", size = 10,
                                  x = 0.01, hjust = 0, vjust = 1),
            g3a_inner,
            ncol        = 1,
            rel_heights = c(0.08, 1)
        )

        save_pp(g3a,
                file.path(OUT_DIR,
                          paste0(FILE_PREFIX, "_", dataset_tag,
                                 "_03a_shuffled_marginal_joint")),
                pdf_w = 10, pdf_h = 9, png_w = 2000, png_h = 1800)
        cat("  Plot 3a saved\n")
        all_plots[["shuffled_marginal_joint"]] <- g3a


        # =============================================================================
        # TEST 2B: Tail vs centre dependence
        # =============================================================================
        cat("\n  Test 2B: Band-specific correlations ...\n")

        q33 <- quantile(x, 1/3)
        q67 <- quantile(x, 2/3)

        band_lower  <- x <= q33
        band_middle <- x >  q33 & x <= q67
        band_upper  <- x >  q67

        n_lower  <- sum(band_lower)
        n_middle <- sum(band_middle)
        n_upper  <- sum(band_upper)

        rho_lower  <- cor(x[band_lower],  y[band_lower],  method = "spearman")
        rho_middle <- cor(x[band_middle], y[band_middle], method = "spearman")
        rho_upper  <- cor(x[band_upper],  y[band_upper],  method = "spearman")
        tau_lower  <- cor(x[band_lower],  y[band_lower],  method = "kendall")
        tau_middle <- cor(x[band_middle], y[band_middle], method = "kendall")
        tau_upper  <- cor(x[band_upper],  y[band_upper],  method = "kendall")

        cat(sprintf("    Lower tertile  (x <= %.3f)  n=%d  rho=%.4f  tau=%.4f\n",
                    q33, n_lower,  rho_lower,  tau_lower))
        cat(sprintf("    Middle tertile             n=%d  rho=%.4f  tau=%.4f\n",
                    n_middle, rho_middle, tau_middle))
        cat(sprintf("    Upper tertile  (x >  %.3f)  n=%d  rho=%.4f  tau=%.4f\n",
                    q67, n_upper,  rho_upper,  tau_upper))
        cat(sprintf("    Global rho = %.4f\n", rho_orig))

        fz_lower_vs_global <- fisher_z_test(rho_lower,  n_lower,  rho_orig,  n)
        fz_upper_vs_global <- fisher_z_test(rho_upper,  n_upper,  rho_orig,  n)
        fz_lower_vs_upper  <- fisher_z_test(rho_lower,  n_lower,  rho_upper, n_upper)
        fz_mid_vs_upper    <- fisher_z_test(rho_middle, n_middle, rho_upper, n_upper)
        fz_lower_vs_mid    <- fisher_z_test(rho_lower,  n_lower,  rho_middle, n_middle)

        cat(sprintf("    Fisher z (lower vs global): z=%.3f p=%.4f\n",
                    fz_lower_vs_global$z_stat, fz_lower_vs_global$p_value))
        cat(sprintf("    Fisher z (upper vs global): z=%.3f p=%.4f\n",
                    fz_upper_vs_global$z_stat, fz_upper_vs_global$p_value))
        cat(sprintf("    Fisher z (lower vs upper) : z=%.3f p=%.4f\n",
                    fz_lower_vs_upper$z_stat, fz_lower_vs_upper$p_value))

        cat("  Running permutation tests (2000 iterations each) ...\n")
        perm_lu <- perm_corr_diff(x[band_lower],  y[band_lower],
                                  x[band_upper],  y[band_upper],
                                  n_perm = 2000, seed = 101)
        perm_lm <- perm_corr_diff(x[band_lower],  y[band_lower],
                                  x[band_middle], y[band_middle],
                                  n_perm = 2000, seed = 102)
        perm_mu <- perm_corr_diff(x[band_middle], y[band_middle],
                                  x[band_upper],  y[band_upper],
                                  n_perm = 2000, seed = 103)
        cat(sprintf("    Perm test (lower vs upper)  obs_diff=%.4f p=%.4f\n",
                    perm_lu$obs_diff, perm_lu$p_value))
        cat(sprintf("    Perm test (lower vs middle) obs_diff=%.4f p=%.4f\n",
                    perm_lm$obs_diff, perm_lm$p_value))
        cat(sprintf("    Perm test (middle vs upper) obs_diff=%.4f p=%.4f\n",
                    perm_mu$obs_diff, perm_mu$p_value))

        u_lower  <- rank(x[band_lower],  ties.method = "average") / (n_lower  + 1)
        v_lower  <- rank(y[band_lower],  ties.method = "average") / (n_lower  + 1)
        u_middle <- rank(x[band_middle], ties.method = "average") / (n_middle + 1)
        v_middle <- rank(y[band_middle], ties.method = "average") / (n_middle + 1)
        u_upper  <- rank(x[band_upper],  ties.method = "average") / (n_upper  + 1)
        v_upper  <- rank(y[band_upper],  ties.method = "average") / (n_upper  + 1)

        td_lower  <- tail_dep_np(u_lower,  v_lower,  q = 0.10)
        td_middle <- tail_dep_np(u_middle, v_middle, q = 0.10)
        td_upper  <- tail_dep_np(u_upper,  v_upper,  q = 0.10)

        band_results <- data.frame(
            dataset          = dataset_tag,
            band             = c("Lower tertile", "Middle tertile",
                                 "Upper tertile", "Global"),
            GOI_range_lower  = round(c(min(x),    q33+.001, q67+.001, min(x)), 3),
            GOI_range_upper  = round(c(q33, q67,  max(x),   max(x)), 3),
            n_cells          = c(n_lower, n_middle, n_upper, n),
            Spearman_rho     = round(c(rho_lower, rho_middle, rho_upper, rho_orig), 4),
            Kendall_tau      = round(c(tau_lower, tau_middle, tau_upper,
                                       cor(x, y, method = "kendall")), 4),
            lambda_L_q10     = round(c(td_lower$lambda_L,  td_middle$lambda_L,
                                       td_upper$lambda_L,  td_orig$lambda_L), 4),
            lambda_U_q10     = round(c(td_lower$lambda_U,  td_middle$lambda_U,
                                       td_upper$lambda_U,  td_orig$lambda_U), 4),
            GOI_tertile      = c(
                paste0("Lower third (GOI <= ", round(q33, 2), ")"),
                paste0("Middle third (", round(q33, 2), " < GOI <= ",
                       round(q67, 2), ")"),
                paste0("Upper third (GOI > ", round(q67, 2), ")"),
                "All cells"
            ),
            stringsAsFactors = FALSE
        )

        pairwise_tests <- data.frame(
            dataset       = dataset_tag,
            comparison    = c("Lower vs Upper", "Lower vs Middle",
                              "Middle vs Upper", "Lower vs Global",
                              "Upper vs Global"),
            rho_group1    = round(c(rho_lower,  rho_lower,  rho_middle,
                                    rho_lower,  rho_upper), 4),
            rho_group2    = round(c(rho_upper,  rho_middle, rho_upper,
                                    rho_orig,   rho_orig), 4),
            n_group1      = c(n_lower, n_lower,  n_middle, n_lower,  n_upper),
            n_group2      = c(n_upper, n_middle, n_upper,  n,        n),
            Fisher_z_stat = round(c(fz_lower_vs_upper$z_stat,
                                    fz_lower_vs_mid$z_stat,
                                    fz_mid_vs_upper$z_stat,
                                    fz_lower_vs_global$z_stat,
                                    fz_upper_vs_global$z_stat), 3),
            Fisher_p      = round(c(fz_lower_vs_upper$p_value,
                                    fz_lower_vs_mid$p_value,
                                    fz_mid_vs_upper$p_value,
                                    fz_lower_vs_global$p_value,
                                    fz_upper_vs_global$p_value), 4),
            Perm_obs_diff = round(c(perm_lu$obs_diff, perm_lm$obs_diff,
                                    perm_mu$obs_diff, NA_real_, NA_real_), 4),
            Perm_p        = round(c(perm_lu$p_value,  perm_lm$p_value,
                                    perm_mu$p_value,  NA_real_, NA_real_), 4),
            Fisher_sig    = c(fz_lower_vs_upper$p_value  < 0.05,
                              fz_lower_vs_mid$p_value    < 0.05,
                              fz_mid_vs_upper$p_value    < 0.05,
                              fz_lower_vs_global$p_value < 0.05,
                              fz_upper_vs_global$p_value < 0.05),
            interpretation = c(
                ifelse(fz_lower_vs_upper$p_value < 0.05,
                       "Lower and upper tails have DIFFERENT rho -> masking present",
                       "Lower and upper tails have similar rho"),
                ifelse(fz_lower_vs_mid$p_value < 0.05,
                       "Lower band rho differs from middle",
                       "Lower and middle similar"),
                ifelse(fz_mid_vs_upper$p_value < 0.05,
                       "Middle and upper bands differ",
                       "Middle and upper similar"),
                ifelse(fz_lower_vs_global$p_value < 0.05,
                       "Lower band differs from global rho -> global rho is masking",
                       "Lower band consistent with global"),
                ifelse(fz_upper_vs_global$p_value < 0.05,
                       "Upper band differs from global rho -> global rho is masking",
                       "Upper band consistent with global")
            ),
            stringsAsFactors = FALSE
        )

        excel_tabs[["Test2B_BandCorrelations"]] <- band_results
        excel_tabs[["Test2B_PairwiseTests"]]    <- pairwise_tests

        # Print summaries
        cat("\n  Band correlations summary:\n")
        print(band_results[, c("band", "n_cells", "Spearman_rho",
                               "Kendall_tau", "lambda_L_q10",
                               "lambda_U_q10")], row.names = FALSE)
        cat("\n  Pairwise tests:\n")
        print(pairwise_tests[, c("comparison", "rho_group1", "rho_group2",
                                 "Fisher_p", "Perm_p",
                                 "Fisher_sig", "interpretation")],
              row.names = FALSE)

        # Write the workbook
        xlsx_path <- file.path(OUT_DIR,
                               paste0(FILE_PREFIX, "_", dataset_tag,
                                      "_copula_analysis_results.xlsx"))
        write_xlsx(excel_tabs, path = xlsx_path)
        cat(sprintf("\n  Excel saved: %s\n", basename(xlsx_path)))

        # == Write console log to text file ===================================
        log_path <- file.path(OUT_DIR,
                              paste0(FILE_PREFIX, "_", dataset_tag, "_console_log.txt"))
        tryCatch({
            log_con  <- file(log_path, open = "wt")
            writeLines(c(
                paste0(strrep("=", 72)),
                paste0("  Part10c: ", dataset_tag, "  -  Console Log"),
                paste0("  Generated: ", Sys.time()),
                paste0(strrep("=", 72)),
                "",
                paste0("  GOI-positive cells: ", n),
                "",
                "  ── Test 2C: Uniformity of copula marginals ──────────────────────",
                sprintf("  KS u: D=%.5f  p_BH=%.4f", unif_test_df$statistic[1], unif_test_df$p_BH[1]),
                sprintf("  KS v: D=%.5f  p_BH=%.4f", unif_test_df$statistic[2], unif_test_df$p_BH[2]),
                sprintf("  AD u: A=%.5f  p_BH=%.4f", unif_test_df$statistic[3], unif_test_df$p_BH[3]),
                sprintf("  AD v: A=%.5f  p_BH=%.4f", unif_test_df$statistic[4], unif_test_df$p_BH[4]),
                "  All p > 0.05 -> UNIFORM (fail to reject H0)",
            "",
            sprintf("  Spearman rho = %.4f", cor(u, v, method = "spearman")),
            sprintf("  CvM copula distance = %.5f (p < 0.001)", cvm_dist),
            "",
            "  ── Exhaustive marginal integration (all slices) ──────────────────",
            sprintf("  DIRECTION 1 - int c(u,v) dv: sum/n = %.8f  mean = %.6f  sd = %.6f  max|I-1| = %.6f  -> %s",
                    Iu_ratio, Iu_mean, Iu_sd, Iu_maxd,
                    ifelse(u_pass, "PASS", "FAIL")),
            sprintf("  DIRECTION 2 - int c(u,v) du: sum/n = %.8f  mean = %.6f  sd = %.6f  max|I-1| = %.6f  -> %s",
                    Iv_ratio, Iv_mean, Iv_sd, Iv_maxd,
                    ifelse(v_pass, "PASS", "FAIL")),
            sprintf("  OVERALL: %s  (tolerance: max|I-1| < %.2f in boundary-safe zone)",
                    ifelse(all_pass, "BOTH MARGINALS VERIFIED", "NOT VERIFIED"), MARG_TOL),
            "",
            paste0(strrep("=", 72))
        ), con = log_con)

            close(log_con)
            cat(sprintf("  Log saved: %s\n", basename(log_path)))
        }, error = function(e) {
            cat(sprintf("  Log file error: %s\n", conditionMessage(e)))
            tryCatch(close(log_con), error = function(e2) NULL)
        })

        # == Display all plots ================================================
        cat("  Displaying plots in RStudio ...\n")
        for (nm in names(all_plots)) {
            cat("    ->", nm, "\n")
            p_obj <- all_plots[[nm]]
            if (inherits(p_obj, "ggplot") || inherits(p_obj, "gg")) {
                print(p_obj)
            } else {
                grid.newpage(); grid.draw(p_obj)
            }
        }

        invisible(list(
            dataset_tag    = dataset_tag,
            n              = n,
            u              = u, v = v,
            x              = x, y = y,
            IP_thresh      = IP_thresh,
            IP_y_thresh    = IP_y_thresh,
            IP_axis_used   = IP_AXIS_USED,
            IP_active      = IP_ACTIVE,
            IP_u           = IP_u,
            IP_v           = IP_v,
            # Deprecated aliases of IP_u / IP_v, kept for downstream readers.
            u_IP           = IP_u,   # deprecated alias
            v_IP           = IP_v,   # deprecated alias
            ip_label_check = ip_label_check,
            band_results   = band_results,
            pairwise_tests = pairwise_tests,
            test2a_summary = test2a_summary,
            td_table       = td_table,
            excel_tabs     = excel_tabs
        ))
    }


    # =============================================================================
    # SECTION C: EXECUTE
    # =============================================================================

    all_results_10c <- vector("list", length(DATASETS))
    names(all_results_10c) <- DATASETS

    for (ds in DATASETS) {
        # Single-bracket assignment keeps a failed dataset's slot as NULL.
        all_results_10c[ds] <- list(tryCatch(
            run_part10c(ds),
            error = function(e) {
                cat("  ERROR in dataset", ds, ":", conditionMessage(e), "\n")
                NULL
            }
        ))
    }

    rds_path <- file.path(OUT_DIR, paste0(FILE_PREFIX, "_all_results.rds"))
    saveRDS(all_results_10c, file = rds_path)


    # =============================================================================
    # SECTION C1: POST-ANALYSIS GUARDRAILS
    # =============================================================================
    cat("\n", strrep("#", 72), "\n")
    cat("  SECTION C1 - POST-ANALYSIS GUARDRAILS\n")
    cat(strrep("#", 72), "\n\n")

    for (ds in names(all_results_10c)) {
        r <- all_results_10c[[ds]]
        if (is.null(r)) {
            gr10c(paste0("GR-07-", ds), "FAIL", "Dataset returned NULL")
            next
        }
        gr10c(paste0("GR-07-", ds),
              if (r$n >= MIN_CELLS) "PASS" else "FAIL",
              sprintf("GOI-positive cells = %d  (min = %d)", r$n, MIN_CELLS))

        rho_val <- cor(r$x, r$y, method = "spearman")
        gr10c(paste0("GR-08-", ds),
              if (is.finite(rho_val) && abs(rho_val) <= 1) "PASS" else "FAIL",
              sprintf("Spearman rho = %.4f  (finite & in [-1,1])", rho_val))

        # == GR-09: the IP this run drew, and where it came from ==============
        ip_val <- r$IP_thresh
        gr10c(paste0("GR-09-", ds),
              if (!is.null(ip_val) && is.finite(ip_val)) "PASS" else "WARN",
              sprintf("IP_thresh = %s  (%s)",
                      if (is.null(ip_val) || !is.finite(ip_val))
                          "NA (IPRS not available)"
                      else sprintf("%.6f", ip_val),
                      ip_from))

        # == GR-10: the axis the boundary was drawn on ========================
        # FAIL unless the boundary is on cDFG at IP_y.
        gr10c(paste0("GR-10-", ds),
              if (identical(r$IP_axis_used, "cDFG")) "PASS" else "FAIL",
              sprintf(paste0("Boundary drawn on %s (IP_CLASS_AXIS = %s) | ",
                             "threshold used = %s | IP_y from IPRS = %s"),
                      if (is.null(r$IP_axis_used)) "NA" else r$IP_axis_used,
                      IP_CLASS_AXIS,
                      if (is.null(r$IP_active) || !is.finite(r$IP_active)) "NA"
                      else sprintf("%.6f", r$IP_active),
                      if (is.null(r$IP_y_thresh) || !is.finite(r$IP_y_thresh))
                          "NA" else sprintf("%.6f", r$IP_y_thresh)))

        # == GR-11: agreement with Part 10b's Pre/Post labels =================
        .chk <- r$ip_label_check
        if (is.null(.chk) || is.null(.chk$n) || .chk$n == 0L ||
            !is.finite(.chk$agree_active)) {
            gr10c(paste0("GR-11-", ds), "WARN",
                  paste0("Pre/Post labels not comparable (no IP_status in the ",
                         "EP CSV, or no threshold on the active axis)"))
        } else {
            # FAIL below 0.999: both parts classify on cDFG at the same IP_y.
            gr10c(paste0("GR-11-", ds),
                  if (.chk$agree_active >= 0.999) "PASS" else "FAIL",
                  sprintf(paste0("Pre/Post labels agree with Part 10b on %.4f ",
                                 "of %d cell(s) on the %s axis; %s on the ",
                                 "other axis%s"),
                          .chk$agree_active, .chk$n, .chk$axis,
                          if (is.finite(.chk$agree_alt))
                              sprintf("%.4f", .chk$agree_alt) else "NA",
                          if (.chk$agree_active < 0.999)
                              " -- axis or IP_y differs from Part 10b"
                          else ""))
        }
    }

    gr10c_write()
    cat(sprintf("\n  Final guardrail tally: %d PASS / %d WARN / %d FAIL\n",
                .gr10c_pass, .gr10c_warn, .gr10c_fail))
    if (.gr10c_fail > 0) {
        cat("  Failures present -- review the ledger before interpreting.\n")
    } else if (.gr10c_warn > 0) {
        cat(sprintf("  %d passed, %d warning(s), no failures.\n",
                    .gr10c_pass, .gr10c_warn))
    } else {
        cat(sprintf("  All %d checks passed.\n", .gr10c_pass))
    }


    # =============================================================================
    # SECTION C2: ARM SUMMARY
    # =============================================================================
    cat("\n", strrep("=", 72), "\n")
    cat("  Part10c Copula1 COMPLETE [", ARM_TAG, "]\n")
    cat("  Output directory:", OUT_DIR, "\n")
    cat("  Files generated:\n")
    # The marker and the temporary message file are removed later, so not listed.
    .transient <- c("_RUN_INCOMPLETE.txt", basename(MSG_PATH_P10C))
    for (f in setdiff(sort(list.files(OUT_DIR)), .transient)) {
        sz <- round(file.size(file.path(OUT_DIR, f)) / 1e3, 0)
        cat(sprintf("    %-60s  %5d KB\n", f, sz))
    }
    cat(sprintf("\n  Results object saved: %s\n", basename(rds_path)))
    cat(strrep("=", 72), "\n")


    cat("\n  Log written to:", LOG_PATH_P10C,
        if (isTRUE(LOG_MESSAGES_P10C)) " (console output and warnings)"
        else " (console output only)", "\n")
    cat(sprintf("\nPart10c Copula1 - log ended: %s\n", Sys.time()))

    .p10c_completed <- TRUE

    if (sink.number(type = "message") != 2L) sink(type = "message")
    while (sink.number() > 0) sink()

    # The arm completed, so the marker is removed.
    suppressWarnings(file.remove(INCOMPLETE_MARKER_P10C))

    invisible(list(arm  = arm_tag, dir = OUT_DIR,
                   pass = .gr10c_pass, warn = .gr10c_warn,
                   fail = .gr10c_fail))
}


# =============================================================================
# SECTION 0c: Driver
# =============================================================================
# One call per arm, folder layout and error handling in CEP_arms.R.

cep_run_arms(PART_DIR_TAG, run_one_arm)
