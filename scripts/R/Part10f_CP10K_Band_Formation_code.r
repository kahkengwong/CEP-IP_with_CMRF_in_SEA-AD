# =============================================================================
# Part10f_CP10K_Band_Formation_code.r
# SEA-AD CEP-IP/CMRF, Part 10f, band-formation figures
#   (Figure 6, Supplementary Figures S11 and S12)
#   For each arm and dataset: eleven figures (BF01 to BF11) showing, layer by
#   layer, how the CMRF adaptive tolerance band is built in (u, v) copula space
#   and which cells each layer admits, from the spine alone to the final
#   classification. Counts and their donor spread are written to a workbook.
#
# Usage (from the SEA-AD folder, or set SEAAD_ROOT), after
#   Part10e_CP10K_CMRF_code.r (Part10e_CP10K_Validation1_code.r optional)
#   source("Part10f_CP10K_Band_Formation_code.r")                         # all arms
#   CEP_ARM <- "PVM1_all"; source("Part10f_CP10K_Band_Formation_code.r")  # one arm
#   Needs CEP_arms.R, R_config.R and R_config_HKG8.R in the same folder.
#
# Outputs (Part10f_Band_Formation/Part10f_Band_Formation_<arm>/): figures
#   (PDF, PNG), counts workbook, guardrail ledger and log.
# =============================================================================


# =============================================================================
# SECTION 0: Libraries and shared arm definitions
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(grid)
library(gtable)
library(scales)
library(writexl)

set.seed(42)

# CEP_arms.R holds the arm table, the folder layout and the driver.
if (!exists("CEP_HOME")) CEP_HOME <- Sys.getenv("SEAAD_ROOT", getwd())
.p10f_harness <- file.path(CEP_HOME, "CEP_arms.R")
if (!file.exists(.p10f_harness))
    stop("CEP_arms.R not found at: ", .p10f_harness,
         "\n  It holds the arm table shared by Parts 9 onwards.", call. = FALSE)
source(.p10f_harness)

PART_TAG     <- "Part10f"              # file-name tag
PART_DIR_TAG <- "Part10f_Band_Formation"   # results-folder tag

# == CMRF band parameters =====================================================
# Must match Part 10e CMRF_KAPPA / CMRF_BW / CMRF_COVERAGE; verified per cell.
KAPPA          <- 0.50
BW             <- 0.10
CMRF_COVERAGE  <- 0.95
V10F_BAND_TOL  <- 1e-8   # max |recomputed - stored| delta_adapt before a stop

# == Donor reporting ==========================================================
V10F_TVD_PERM  <- 200L   # label permutations for the TVD null
V10F_SEED      <- 42L

# == Subgroup labelling loaded ================================================
# Loaded for the label-axis check; every panel is coloured from group_orig.
P10F_GROUP_COL <- "group_refined"

U_GRID <- seq(0, 1, length.out = 1001)


# =============================================================================
# SECTION 0d: Pure helpers
# =============================================================================

# Coerces a logical, integer or character flag column to logical.
as_flag <- function(x) {
    if (is.logical(x)) return(x)
    if (is.numeric(x)) return(x != 0)
    y <- toupper(trimws(as.character(x)))
    out <- rep(NA, length(y))
    out[y %in% c("TRUE", "T", "1", "YES")] <- TRUE
    out[y %in% c("FALSE", "F", "0", "NO")] <- FALSE
    out
}

# A flag used in a subset must have no NA (NA gives an all-NA phantom row).
require_flag <- function(x, nm) {
    f <- as_flag(x)
    if (anyNA(f))
        stop(sum(is.na(f)), " missing value(s) in `", nm,
             "`.",
             call. = FALSE)
    f
}

gauss_K <- function(u, IP_u, bw) exp(-(u - IP_u)^2 / (2 * bw^2))

# Reads IP_u / IP_v / IP_x / IP_y, falling back to the u_IP-style names.
ip_field_10f <- function(r, new_nm, old_nm, ds) {
    if (!is.null(r[[new_nm]])) return(as.numeric(r[[new_nm]])[1])
    if (!is.null(r[[old_nm]])) {
        cat(sprintf(paste0("    NOTE: %s read from the deprecated `%s`",
                           " for %s.\n"),
                    new_nm, old_nm, ds))
        return(as.numeric(r[[old_nm]])[1])
    }
    NA_real_
}

# geom_label's border-width argument is `linewidth` from ggplot2 3.5.0 and
# `label.size` before; probe which name the installed version accepts.
.LABEL_BORDER_PARAM <- local({
    accepts <- function(nm) {
        ok <- TRUE
        withCallingHandlers(
            tryCatch({
                lay <- do.call(annotate,
                               c(list("label", x = 1, y = 1, label = "t"),
                                 stats::setNames(list(0), nm)))
                p <- ggplot(data.frame(x = 1, y = 1), aes(x, y)) +
                     geom_blank() + lay
                invisible(ggplot_build(p))
            }, error = function(e) ok <<- FALSE),
            warning = function(w) {
                if (grepl("unknown parameter", conditionMessage(w),
                          ignore.case = TRUE)) ok <<- FALSE
                invokeRestart("muffleWarning")
            })
        ok
    }
    if (accepts("linewidth")) "linewidth"
    else if (accepts("label.size")) "label.size"
    else NA_character_
})

# Reported once if neither name is accepted; labels keep the default border.
if (is.na(.LABEL_BORDER_PARAM))
    cat("NOTE: this ggplot2 accepts neither `linewidth` nor `label.size` on\n",
        "      geom_label; the BF07 QC labels will keep their default border.\n",
        sep = "")

ann_label_10f <- function(...) {
    args <- list("label", ...)
    if (!is.na(.LABEL_BORDER_PARAM)) args[[.LABEL_BORDER_PARAM]] <- 0
    do.call(annotate, args)
}

delta_adapt_fn <- function(u, delta_base, IP_u, kappa = KAPPA, bw = BW)
    delta_base * (1 + kappa * gauss_K(u, IP_u, bw))

show_plot <- function(p) {
    tryCatch({
        if (inherits(p, "ggplot")) print(p) else { grid.newpage(); grid.draw(p) }
    }, error = function(e) invisible(NULL))
}

save_pp <- function(p, base, pdf_w = 10, pdf_h = 7,
                    png_w = 2200, png_h = 1540, res = 200) {
    is_grid <- inherits(p, "grob") || inherits(p, "gtable") ||
        inherits(p, "arrangelist")
    draw <- function() if (is_grid) { grid.newpage(); grid.draw(p) } else print(p)
    tryCatch({
        pdf(paste0(base, ".pdf"), width = pdf_w, height = pdf_h)
        draw(); dev.off()
    }, error = function(e) {
        tryCatch(dev.off(), error = function(e2) invisible(NULL))
        cat("  [save_pp] PDF error:", conditionMessage(e), "\n")
    })
    tryCatch({
        png(paste0(base, ".png"), width = png_w, height = png_h, res = res)
        draw(); dev.off()
    }, error = function(e) {
        tryCatch(dev.off(), error = function(e2) invisible(NULL))
        cat("  [save_pp] PNG error:", conditionMessage(e), "\n")
    })
    invisible(NULL)
}

th_band <- function(bs = 10) {
    theme_minimal(base_size = bs) +
        theme(
            plot.title       = element_text(face = "bold", size = bs + 1),
            plot.subtitle    = element_text(size = bs - 1, colour = "gray30"),
            axis.title       = element_text(face = "bold", size = bs - 1),
            axis.text        = element_text(size = bs - 2),
            panel.grid.major = element_line(colour = "gray90"),
            panel.grid.minor = element_blank(),
            panel.border     = element_rect(colour = "black", fill = NA,
                                            linewidth = 0.35),
            plot.background  = element_rect(fill = "white", colour = NA),
            panel.background = element_rect(fill = "white", colour = NA),
            legend.background = element_rect(fill = alpha("white", 0.8),
                                             colour = NA)
        )
}

# == Palette (Part 10e Plot 05 scheme) ========================================
FOUR_COLS <- c(
    "Pre-IP TREP"      = "#E69F00",
    "Pre-IP non-TREP"  = "#1B3A8C",
    "Post-IP TREP"     = "#009E73",
    "Post-IP non-TREP" = "#D94F4F"
)
GRAY_COL   <- "gray65"
ALL_COLS   <- c("Non-selected" = GRAY_COL, FOUR_COLS)
ALL_LEVS   <- names(ALL_COLS)

COL_SPINE        <- "#00695C"
COL_IP           <- "#1A237E"
COL_ADAPT        <- "#D55E00"
COL_ADAPT_BORDER <- "#B84A00"
COL_FLAT         <- "#0072B2"

class_levels <- c("Original + Refined", "Recovered",
                  "Original only", "Excluded")

# ===========================================================================
# Donor spread of a set of cells, with its own permutation null
# ===========================================================================
# TVD is biased upward at small n, so each TVD is compared with the 95th
# percentile of a null built by permuting set membership with donor fixed.
v10f_tvd <- function(pa, pb) {
    lev <- union(names(pa), names(pb))
    a <- as.numeric(pa[lev]); a[is.na(a)] <- 0
    b <- as.numeric(pb[lev]); b[is.na(b)] <- 0
    0.5 * sum(abs(a - b))
}

v10f_donor_row <- function(label, donor_set, donor_ref, donor_pool,
                           n_perm = V10F_TVD_PERM, seed = V10F_SEED) {
    d  <- as.character(donor_set); d <- d[!is.na(d) & nzchar(d)]
    dr <- as.character(donor_ref); dr <- dr[!is.na(dr) & nzchar(dr)]
    n_set <- length(donor_set)
    if (length(d) == 0L)
        return(data.frame(set = label, n_cells = n_set, n_donors = 0L,
                          top_donor_share = NA_real_, TVD_vs_OrigRef = NA_real_,
                          TVD_null_p95 = NA_real_, stringsAsFactors = FALSE))
    p <- table(d) / length(d)
    tvd <- if (length(dr) == 0L) NA_real_ else
        round(v10f_tvd(p, table(dr) / length(dr)), 4)

    null_p95 <- NA_real_
    dp <- as.character(donor_pool); dp <- dp[!is.na(dp) & nzchar(dp)]
    if (length(dr) >= 2L && length(d) >= 2L &&
        length(dp) >= length(d) + length(dr)) {
        set.seed(seed)
        # membership permuted with donor held fixed: the set and the reference
        # are drawn as disjoint blocks of their observed sizes from the pool
        v <- vapply(seq_len(n_perm), function(i) {
            idx <- sample.int(length(dp))
            a <- table(dp[idx[seq_len(length(d))]])
            b <- table(dp[idx[length(d) + seq_len(length(dr))]])
            v10f_tvd(a / sum(a), b / sum(b))
        }, numeric(1L))
        null_p95 <- round(as.numeric(quantile(v, 0.95, names = FALSE)), 4)
    }
    data.frame(set = label, n_cells = n_set, n_donors = length(p),
               top_donor_share = round(max(p), 4), TVD_vs_OrigRef = tvd,
               TVD_null_p95 = null_p95, stringsAsFactors = FALSE)
}


# =============================================================================
# SECTION 0b: The pipeline for one arm
# =============================================================================

run_one_arm <- function(arm_tag, config_file, out_dir, file_tag) {

    config_path <- file.path(CEP_HOME, config_file)
    if (!file.exists(config_path))
        stop("Configuration not found: ", config_path, call. = FALSE)

    .p10f_preamble <- character(0)
    assign("CEP_ARM", arm_tag, envir = globalenv())
    on.exit(suppressWarnings(rm("CEP_ARM", envir = globalenv())), add = TRUE)

    .p10f_preamble <- tryCatch(
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

    .p10f_needed <- c("ARM_TAG", "R_EXPORT_ROOT", "GOI", "DFG",
                      "DFG_COMPOSITE_LABEL", "PRIMARY_STYPE", "NORM_TAG",
                      "NEOC_DIR", "NONNC_DIR", "NEOC_REGIONS", "NONNC_REGIONS")
    .p10f_absent <- .p10f_needed[!vapply(.p10f_needed, exists, logical(1))]
    if (length(.p10f_absent) > 0)
        stop("The configuration loaded from\n    ", config_path,
             "\n  does not define: ", paste(.p10f_absent, collapse = ", "), call. = FALSE)
    if (!identical(ARM_TAG, arm_tag))
        stop("The configuration reports ARM_TAG '", ARM_TAG,
             "' but this arm is '", arm_tag, "'.", call. = FALSE)

    .p10f_preamble <- c(.p10f_preamble, utils::capture.output(
        cat("  Config OK: ", ARM_TAG, " | ", GOI, " + ", length(DFG), " DFG | ",
            NORM_TAG, "\n", sep = "")))

    set.seed(V10F_SEED)
    DATASETS <- c("NeoC", "NonNeoC")

    # == Input paths ==========================================================
    CMRF_DIR  <- cep_arm_dir_for("Part10e_CMRF", arm_tag)
    CMRF_RDS  <- file.path(CMRF_DIR, paste0(
        cep_out_prefix("Part10e", NORM_TAG, cep_file_tag_for(arm_tag)),
        "_CMRF_all_results.rds"))
    V1_DIR    <- cep_arm_dir_for("Part10e_Validation1", arm_tag)
    V1_PREFIX <- paste0(cep_out_prefix("Part10e", NORM_TAG,
                                       cep_file_tag_for(arm_tag)), "_V1")
    p10f_v1_rds <- function(ds)
        file.path(V1_DIR, paste0(V1_PREFIX, "_", ds, "_results.rds"))

    # == Output ===============================================================
    OUT_DIR <- out_dir
    dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
    FILE_PREFIX <- cep_out_prefix(PART_TAG, NORM_TAG, cep_file_tag_for(arm_tag))

    # == Logging ==============================================================
    LOG_MESSAGES_10F <- TRUE
    LOG_PATH_10F <- file.path(OUT_DIR, paste0(FILE_PREFIX, "_console_log_",
                                              format(Sys.time(), "%Y%m%d_%H%M%S"),
                                              ".txt"))
    MSG_PATH_10F <- file.path(OUT_DIR, paste0(FILE_PREFIX,
                                              "_console_log_messages.tmp"))
    while (sink.number() > 0) sink()
    if (sink.number(type = "message") != 2L) sink(type = "message")
    # Set before on.exit() registers so the handler never meets a missing object.
    .p10f_completed <- FALSE
    .p10f_log_con <- file(LOG_PATH_10F, open = "wt")
    .p10f_msg_con <- if (isTRUE(LOG_MESSAGES_10F))
        file(MSG_PATH_10F, open = "wt") else NULL
    sink(.p10f_log_con, append = FALSE, split = TRUE)
    if (!is.null(.p10f_msg_con)) sink(.p10f_msg_con, type = "message")
    # Each clean-up stage has its own tryCatch so the sinks are always released.
    on.exit({
        tryCatch({
            if (!isTRUE(.p10f_completed))
                cat("\n", strrep("!", 72),
                    "\n  Arm did not complete; the error is printed to the console.\n",
                    strrep("!", 72), "\n", sep = "")
        }, error = function(e) invisible(NULL))
        tryCatch({
            if (sink.number(type = "message") != 2L) sink(type = "message")
            while (sink.number() > 0) sink()
        }, error = function(e) invisible(NULL))
        tryCatch(if (!is.null(.p10f_msg_con)) close(.p10f_msg_con),
                 error = function(e) invisible(NULL))
        tryCatch(close(.p10f_log_con), error = function(e) invisible(NULL))
        tryCatch({
            if (!is.null(.p10f_msg_con) && file.exists(MSG_PATH_10F)) {
                .msg_txt <- readLines(MSG_PATH_10F, warn = FALSE)
                if (length(.msg_txt) > 0)
                    cat(c("", strrep("=", 72),
                          "  MESSAGE STREAM (warnings, messages, errors)",
                          strrep("=", 72), .msg_txt, ""),
                        sep = "\n", file = LOG_PATH_10F, append = TRUE)
                unlink(MSG_PATH_10F)
            }
        }, error = function(e) invisible(NULL))
    }, add = TRUE)

    .p10f_warn_opt  <- getOption("warn")
    options(warn = 1)
    on.exit(options(warn = .p10f_warn_opt), add = TRUE)
    cat(sprintf("Part10f Band Formation - log started: %s\n\n", Sys.time()))
    if (length(.p10f_preamble) > 0)
        cat(paste(.p10f_preamble, collapse = "\n"), "\n\n", sep = "")

    cat("\n", strrep("=", 72), "\n")
    cat("  Part10f - progressive band formation in copula space\n")
    cat("  Arm          :", ARM_TAG, "\n")
    cat("  GOI          :", GOI, "\n")
    cat("  Normalisation:", NORM_TAG, "\n")
    cat("  Band         : kappa =", KAPPA, "| bw =", BW,
        "| coverage =", CMRF_COVERAGE, "\n")
    cat("  Colour by    : group_orig (SKCS core, all 11 panels)\n")
    cat("  Labels loaded:", P10F_GROUP_COL, "\n")
    cat("  Part 10e RDS :", CMRF_RDS, "\n")
    cat("  Output dir   :", OUT_DIR, "\n")
    cat(strrep("=", 72), "\n\n")

    INCOMPLETE_MARKER_10F <- file.path(OUT_DIR, "_RUN_INCOMPLETE.txt")
    writeLines(c("This arm did not finish.", "",
                 paste("Arm     :", arm_tag),
                 paste("Part    :", PART_DIR_TAG),
                 paste("Started :", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), "",
                 "Any files in this folder are from a partial run.",
                 "This file is removed automatically when the arm completes."),
               INCOMPLETE_MARKER_10F)


    # =======================================================================
    # SECTION A0: Guardrail ledger
    # =======================================================================
    GR_LOG_10F <- file.path(OUT_DIR, paste0(FILE_PREFIX, "_Guardrail_Log.txt"))
    .gr10f_log  <- list()
    .gr10f_pass <- 0L; .gr10f_fail <- 0L; .gr10f_warn <- 0L

    gr10f <- function(check_id, status, detail) {
        # A zero-length detail would vanish from console and ledger.
        if (length(detail) != 1L || is.na(detail))
            detail <- paste0("[detail unavailable: the message was built from a",
                             " missing value and came out zero-length]")
        # An unknown status is shown and counted as a failure.
        if (!status %in% c("PASS", "FAIL", "WARN")) {
            detail <- paste0("[status '", as.character(status)[1],
                             "' is not PASS/FAIL/WARN] ", detail)
            status <- "FAIL"
        }
        tag <- switch(status, PASS = "[PASS]", FAIL = "[FAIL]", WARN = "[WARN]")
        cat(sprintf("  %s  %-18s  %s\n", tag, check_id, detail))
        .gr10f_log[[length(.gr10f_log) + 1L]] <<- list(
            id = check_id, status = status, detail = detail)
        if (status == "PASS") .gr10f_pass <<- .gr10f_pass + 1L
        if (status == "FAIL") .gr10f_fail <<- .gr10f_fail + 1L
        if (status == "WARN") .gr10f_warn <<- .gr10f_warn + 1L
    }

    gr10f_write <- function() {
        con <- file(GR_LOG_10F, "wt")
        # Close the connection on any exit.
        on.exit(try(close(con), silent = TRUE), add = TRUE)
        writeLines(c(strrep("=", 72),
                     paste0("Part 10f Band Formation - Guardrail Log [",
                            NORM_TAG, "]"),
                     paste("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
                     paste("Arm :", ARM_TAG),
                     paste("Band: kappa =", KAPPA, "| bw =", BW,
                           "| coverage =", CMRF_COVERAGE),
                     strrep("=", 72), ""), con)
        for (item in .gr10f_log)
            writeLines(sprintf("  [%s]  %-18s  %s", item$status, item$id,
                               item$detail), con)
        writeLines(c("", strrep("=", 72),
                     sprintf("  Total : %d",
                             .gr10f_pass + .gr10f_fail + .gr10f_warn),
                     sprintf("  PASS  : %d", .gr10f_pass),
                     sprintf("  WARN  : %d", .gr10f_warn),
                     sprintf("  FAIL  : %d", .gr10f_fail),
                     if (.gr10f_fail == 0 && .gr10f_warn == 0)
                         "\n  Result: no failures, no warnings"
                     else if (.gr10f_fail == 0)
                         "\n  Result: no failures; review the warnings above"
                     else "\n  Result: failures present",
                     strrep("=", 72)), con)
        cat("  Guardrail log written:", GR_LOG_10F, "\n")
    }

    cat("  SECTION A0 - PRE-ANALYSIS GUARDRAILS\n\n")
    gr10f("GR-10F-01-CMRF",
          if (file.exists(CMRF_RDS)) "PASS" else "FAIL",
          if (file.exists(CMRF_RDS))
              sprintf("Part 10e results: %s (%.1f MB)",
                      file.path(basename(dirname(CMRF_RDS)), basename(CMRF_RDS)),
                      file.size(CMRF_RDS) / 1e6)
          else sprintf("not found: %s", CMRF_RDS))
    if (!file.exists(CMRF_RDS))
        stop("Part 10e results RDS not found: ", CMRF_RDS,
             call. = FALSE)
    all_10e <- readRDS(CMRF_RDS)

    for (ds in DATASETS) {
        .v <- p10f_v1_rds(ds)
        gr10f(paste0("GR-10F-02-", ds),
              if (file.exists(.v)) "PASS" else "WARN",
              if (file.exists(.v))
                  sprintf(paste0("Validation1 RDS present; %s loaded for the",
                                 " axis check"),
                          P10F_GROUP_COL)
              else sprintf(paste0("Validation1 RDS absent; the loaded labels",
                                  " fall back to group_orig (%s)"),
                           basename(.v)))
    }
    cat(sprintf("\n  Pre-analysis guardrails: %d PASS / %d WARN / %d FAIL\n\n",
                .gr10f_pass, .gr10f_warn, .gr10f_fail))


    # =======================================================================
    # SECTION A1: Donor lookup
    # =======================================================================
    # Donor matched from the per-region gam_data CSVs by barcode.
    .p10f_donor_cache <- list()
    p10f_donor_map <- function(ds) {
        if (!is.null(.p10f_donor_cache[[ds]])) return(.p10f_donor_cache[[ds]])
        gam_dir <- if (identical(ds, "NeoC")) NEOC_DIR else NONNC_DIR
        regs    <- if (identical(ds, "NeoC")) NEOC_REGIONS else NONNC_REGIONS
        stype   <- gsub("-", "_", gsub("/", "_", PRIMARY_STYPE))
        bc <- character(0); don <- character(0)
        for (reg in regs) {
            gp <- file.path(gam_dir, "gam_data",
                            paste0(reg, "_", stype, "_selected.csv"))
            if (!file.exists(gp)) next
            d <- read.csv(gp, stringsAsFactors = FALSE, check.names = FALSE)
            if (!all(c("barcode", "donor") %in% colnames(d)))
                stop(sprintf("Region %s: gam_data lacks barcode or donor.", reg),
                     call. = FALSE)
            bc  <- c(bc,  d$barcode)
            don <- c(don, as.character(d$donor))
        }
        if (length(bc) == 0L)
            stop("No gam_data CSVs found for ", ds, call. = FALSE)
        # Duplicated barcodes must agree on the donor; a conflict stops the run.
        dup <- duplicated(bc)
        if (any(dup)) {
            .first <- setNames(don[!dup], bc[!dup])
            .clash <- bc[dup][don[dup] != .first[bc[dup]]]
            .clash <- unique(.clash[!is.na(.clash)])
            if (length(.clash) > 0)
                stop(sprintf(paste0("%s: %d barcode(s) appear in more than one",
                                    " region with different donors, e.g. %s."),
                             ds, length(.clash),
                             paste(utils::head(.clash, 3), collapse = ", ")),
                     call. = FALSE)
            cat(sprintf(paste0("    %d duplicated barcode(s) across regions,",
                               " all agreeing on the donor; keeping one each\n"),
                        sum(dup)))
            bc <- bc[!dup]; don <- don[!dup]
        }
        m <- setNames(don, bc)
        .p10f_donor_cache[[ds]] <<- m
        m
    }


    # =======================================================================
    # SECTION B: Main function (per dataset)
    # =======================================================================

    run_V10F_dataset <- function(dataset_tag) {

        ds <- dataset_tag
        cat("\n", strrep("=", 72), "\n")
        cat("  Part10f Band Formation:", ds, " [", ARM_TAG, "]\n")
        cat(strrep("=", 72), "\n\n")

        r <- all_10e[[ds]]
        if (is.null(r))
            stop("Part 10e returned no results for ", ds, call. = FALSE)
        need <- c("delta", "df_refined")
        absent <- need[!need %in% names(r) |
                           vapply(need, function(k) is.null(r[[k]]), logical(1))]
        if (is.null(r$IP_u) && is.null(r$u_IP)) absent <- c("IP_u", absent)
        if (length(absent) > 0)
            stop("Part 10e results for ", ds, " lack: ",
                 paste(absent, collapse = ", "), call. = FALSE)

        df         <- as.data.frame(r$df_refined, stringsAsFactors = FALSE)
        IP_u       <- ip_field_10f(r, "IP_u", "u_IP", ds)
        delta_base <- as.numeric(r$delta)[1]
        if (!is.finite(IP_u) || !is.finite(delta_base) || delta_base <= 0)
            stop("Part 10e gave a non-usable IP_u (", IP_u, ") or delta (",
                 delta_base, ") for ", ds, call. = FALSE)

        # stage2_res is NULL when Stage 2 removed nothing.
        stage2_p <- tryCatch({
            v <- r$stage2_res$wilcox_p_final
            if (is.null(v)) NA_real_ else as.numeric(v)[1]
        }, error = function(e) NA_real_)
        # The RSk-fCV gate is dual (p < alpha and relative delta > noise floor);
        # both terms are read.
        .s2 <- r$stage2_res
        stage2_delta <- tryCatch(
            if (is.null(.s2$relative_delta_final)) NA_real_
            else as.numeric(.s2$relative_delta_final)[1],
            error = function(e) NA_real_)
        stage2_nf <- tryCatch(
            if (is.null(.s2$noise_floor)) NA_real_
            else as.numeric(.s2$noise_floor)[1],
            error = function(e) NA_real_)
        stage2_alpha <- tryCatch(
            if (is.null(.s2$alpha)) NA_real_ else as.numeric(.s2$alpha)[1],
            error = function(e) NA_real_)
        # Small p values in scientific notation.
        fmt_p_10f <- function(x) {
            if (!is.finite(x)) "NA"
            else if (x > 0 && x < 1e-4) format(x, digits = 3, scientific = TRUE)
            else sprintf("%.4f", x)
        }

        # Measured spine monotonicity from Part 10e.
        mono <- r$spine_monotone
        mono_txt <- if (is.null(mono)) {
            "monotonicity not recorded by Part 10e"
        } else if (is.list(mono)) {
            ok <- isTRUE(mono$monotone) || isTRUE(mono$is_monotone)
            md <- suppressWarnings(as.numeric(mono$min_diff))
            sprintf("spine %s on the grid%s",
                    if (ok) "non-decreasing" else "DECREASING somewhere",
                    if (length(md) == 1L && is.finite(md))
                        sprintf(" (min diff %.3e)", md) else "")
        } else {
            sprintf("spine %s on the grid",
                    if (isTRUE(as.logical(mono)[1])) "non-decreasing"
                    else "DECREASING somewhere")
        }
        cat("    Spine:", mono_txt, "\n")

        required_cols <- c("barcode", "u_all", "v_all", "g_u", "resid_abs",
                           "delta_adapt", "pass_stage1", "is_orig_skcs",
                           "is_recovered", "cell_class", "group_orig",
                           "GOI_expr")
        miss <- setdiff(required_cols, names(df))
        if (length(miss) > 0)
            stop("Part 10e df_refined lacks: ", paste(miss, collapse = ", "),
                 call. = FALSE)
        df$barcode <- as.character(df$barcode)

        # ggplot drops non-finite points while counts use nrow(), so stop on any.
        .n_uv_bad <- sum(!is.finite(as.numeric(df$u_all)) |
                         !is.finite(as.numeric(df$v_all)))
        if (.n_uv_bad > 0)
            stop(.n_uv_bad, " cell(s) have a non-finite u_all or v_all in ",
                 "df_refined for ", ds, ".",
                 call. = FALSE)
        # approx() needs at least two distinct u values.
        .n_u_distinct <- length(unique(as.numeric(df$u_all)))
        if (.n_u_distinct < 2L)
            stop("df_refined for ", ds, " has ", .n_u_distinct,
                 " distinct u_all value(s); no band can be drawn.", call. = FALSE)
        .n_bc_dup <- sum(duplicated(df$barcode))
        if (.n_bc_dup > 0)
            stop(.n_bc_dup, " duplicated barcode(s) in df_refined for ", ds,
                 ".", call. = FALSE)

        # == Donor ============================================================
        donor_map <- p10f_donor_map(ds)
        df$donor  <- unname(donor_map[df$barcode])
        n_no_donor <- sum(is.na(df$donor))
        cat(sprintf("    Cells %d | donors %d distinct | %d without a donor\n",
                    nrow(df),
                    length(unique(df$donor[!is.na(df$donor)])), n_no_donor))

        # == Subgroup labels: Validation 1's refit when available =============
        group_src <- "Part 10e group_orig"
        df$group_plot <- as.character(df$group_orig)
        if (identical(P10F_GROUP_COL, "group_refined")) {
            v1p <- p10f_v1_rds(ds)
            if (file.exists(v1p)) {
                v1 <- tryCatch(readRDS(v1p), error = function(e) NULL)
                rc <- if (!is.null(v1)) v1$rc_v3 else NULL
                if (!is.null(rc) && all(c("barcode", "group_refined") %in% names(rc))) {
                    # Duplicates in rc_v3: the first label is used.
                    .v1_dup <- sum(duplicated(as.character(rc$barcode)))
                    if (.v1_dup > 0)
                        cat(sprintf(paste0("    WARNING: Validation 1 rc_v3 has",
                                           " %d duplicated barcode(s); the",
                                           " first label for each is used\n"),
                                    .v1_dup))
                    m <- setNames(as.character(rc$group_refined),
                                  as.character(rc$barcode))
                    hit <- df$barcode %in% names(m)
                    df$group_plot[hit] <- unname(m[df$barcode[hit]])
                    group_src <- sprintf("Validation1 group_refined (%d of %d cells matched)",
                                         sum(hit), nrow(df))
                }
            }
        }
        cat("    Subgroup labels from:", group_src, "\n")

        # Label-axis check inputs: IP_y from Part 10e and cDFG_expr in df_refined.
        IP_y_10f <- suppressWarnings(ip_field_10f(r, "IP_y", "y_IP", ds))
        # IP_v = Fc(IP_y), the Pre/Post boundary in copula space; recomputed from
        # df as sum(cDFG_expr <= IP_y)/(n+1) when the RDS lacks it.
        IP_v_10f <- suppressWarnings(ip_field_10f(r, "IP_v", "v_IP", ds))
        IP_v_src <- "Part 10e RDS"
        if (!is.finite(IP_v_10f) && is.finite(IP_y_10f) &&
            all(c("cDFG_expr", "v_all") %in% names(df))) {
            .vok <- !is.na(df$v_all) & !is.na(df$cDFG_expr)
            if (any(.vok)) {
                IP_v_10f <- sum(df$cDFG_expr[.vok] <= IP_y_10f) / (sum(.vok) + 1)
                IP_v_src <- "recomputed here from df_refined"
            }
        }
        cat(sprintf("    Pre/Post boundary in copula space: IP_v = %s  (%s)\n",
                    if (is.finite(IP_v_10f)) sprintf("%.6f", IP_v_10f) else "NA",
                    if (is.finite(IP_v_10f)) IP_v_src
                    else "no IP_y or no cDFG_expr; not drawn"))

        # Horizontal IP_v line (the Pre/Post boundary the colours follow), drawn
        # beside the IP_u band-kernel centre; list() when IP_v is unavailable.
        COL_VIP <- "gray10"
        ip_hline_v <- function(linewidth = 0.5, label = FALSE, x_lab = 0.015) {
            if (!is.finite(IP_v_10f)) return(list())
            out <- list(geom_hline(yintercept = IP_v_10f, colour = COL_VIP,
                                   linetype = "dashed", linewidth = linewidth))
            if (label)
                out <- c(out, list(annotate("text", x = x_lab,
                                            y = IP_v_10f + 0.022,
                                            label = sprintf("IP_v=%.3f", IP_v_10f),
                                            colour = COL_VIP, hjust = 0,
                                            size = 2.8, fontface = "bold")))
            out
        }
        # Caption line naming the two IP lines.
        vip_note <- if (is.finite(IP_v_10f))
            sprintf(paste0("Navy vertical = IP_u (band-kernel centre); dark red",
                           " dashed = IP_v = Fc(IP_y) = %.4f, the Pre/Post",
                           " boundary the four-group colours are assigned on",
                           " (cDFG at IP_y)."), IP_v_10f)
        else paste0("Navy vertical = IP_u (band-kernel centre). The Pre/Post",
                    " boundary IP_v is not drawn: Part 10e supplied no IP_y.")

        axis_chk <- list(n = 0L, agree = NA_real_, IP_y = IP_y_10f,
                         have = FALSE)
        # axis_chk carries the reason when the check cannot run.
        .why <- if (!is.finite(IP_y_10f) && !("cDFG_expr" %in% names(df)))
            "Part 10e supplied no IP_y AND df_refined has no cDFG_expr"
        else if (!is.finite(IP_y_10f))
            "Part 10e supplied no IP_y"
        else if (!("cDFG_expr" %in% names(df)))
            "df_refined has no cDFG_expr"
        else "no cell has both a Pre/Post label and a non-missing cDFG_expr"
        axis_chk$why <- .why
        # The column held in group_plot.
        .post_col <- if (grepl("group_refined", group_src, fixed = TRUE))
            "group_refined" else "group_orig"
        # Label-axis agreement for group_plot and group_orig.
        .axis_agree <- function(labcol) {
            .lab <- ifelse(grepl("Post-IP", labcol, fixed = TRUE), "Post-IP",
                    ifelse(grepl("Pre-IP",  labcol, fixed = TRUE), "Pre-IP",
                           NA_character_))
            .rule <- ifelse(is.na(df$cDFG_expr), NA_character_,
                     ifelse(df$cDFG_expr < IP_y_10f, "Pre-IP", "Post-IP"))
            .ok <- !is.na(.lab) & !is.na(.rule)
            if (!any(.ok)) return(list(n = 0L, agree = NA_real_))
            list(n = sum(.ok), agree = mean(.lab[.ok] == .rule[.ok]))
        }
        if (is.finite(IP_y_10f) && "cDFG_expr" %in% names(df)) {
            .a_plot <- .axis_agree(df$group_plot)
            .a_orig <- .axis_agree(as.character(df$group_orig))
            if (.a_plot$n > 0L)
                axis_chk <- list(n = .a_plot$n, agree = .a_plot$agree,
                                 IP_y = IP_y_10f, have = TRUE, why = "",
                                 n_orig = .a_orig$n, agree_orig = .a_orig$agree,
                                 have_orig = .a_orig$n > 0L)
        }
        # Present on every path, so the guardrail never reaches a NULL field.
        if (is.null(axis_chk$n_orig))     axis_chk$n_orig     <- 0L
        if (is.null(axis_chk$agree_orig)) axis_chk$agree_orig <- NA_real_
        if (is.null(axis_chk$have_orig))  axis_chk$have_orig  <- FALSE
        cat(sprintf("    Label axis check: %s\n",
                    if (!isTRUE(axis_chk$have))
                        sprintf("not run (IP_y %s, cDFG_expr %s)",
                                if (is.finite(IP_y_10f)) "present" else "absent",
                                if ("cDFG_expr" %in% names(df)) "present"
                                else "absent")
                    else sprintf(paste0("%s = %.2f%% of %d cell(s)",
                                        " [loaded, not drawn]",
                                        " | group_orig = %s [all 11 panels]",
                                        " | rule `cDFG_expr < IP_y`,",
                                        " IP_y = %.6f"),
                                 .post_col,
                                 100 * axis_chk$agree, axis_chk$n,
                                 if (isTRUE(axis_chk$have_orig))
                                     sprintf("%.2f%% of %d cell(s)",
                                             100 * axis_chk$agree_orig,
                                             axis_chk$n_orig)
                                 else "not run (no labelled cell)",
                                 axis_chk$IP_y)))

        # == Flags must be NA-free before any subset ==========================
        df$is_orig_skcs <- require_flag(df$is_orig_skcs, "is_orig_skcs")
        df$is_recovered <- require_flag(df$is_recovered, "is_recovered")
        pass_stage1     <- require_flag(df$pass_stage1, "pass_stage1")
        df$pass_stage1  <- pass_stage1
        resid_abs       <- as.numeric(df$resid_abs)
        delta_adapt     <- as.numeric(df$delta_adapt)
        if (anyNA(resid_abs) || anyNA(delta_adapt))
            stop("resid_abs or delta_adapt holds NA in df_refined for ", ds,
                 call. = FALSE)

        # == The drawn band must be the applied band ==========================
        da_check  <- delta_adapt_fn(as.numeric(df$u_all), delta_base, IP_u)
        band_dev  <- max(abs(da_check - delta_adapt))
        if (!is.finite(band_dev) || band_dev > V10F_BAND_TOL)
            stop(sprintf(paste0("Drawn band does not match the Part 10e band:",
                                " max |recomputed - stored| delta_adapt = %.3e",
                                " (tolerance %.1e); KAPPA = %.4f, BW = %.4f."),
                        band_dev, V10F_BAND_TOL, KAPPA, BW), call. = FALSE)
        cat(sprintf("    Band verified against Part 10e: max |dev| = %.3e\n",
                    band_dev))

        # == Ribbon grid ======================================================
        # approx() needs finite g_u; tied u values are reported.
        n_gu_na <- sum(!is.finite(as.numeric(df$g_u)))
        if (n_gu_na > 0)
            stop(n_gu_na, " cell(s) have a non-finite g_u in df_refined for ",
                 ds, "; the spine cannot be interpolated onto the grid.",
                 call. = FALSE)
        n_u_tied <- sum(duplicated(df$u_all))
        if (n_u_tied > 0)
            cat(sprintf("    %d tied u value(s); approx() collapses them\n",
                        n_u_tied))

        u_full  <- as.numeric(U_GRID)
        g_full  <- as.numeric(suppressWarnings(
            approx(x = df$u_all, y = df$g_u, xout = u_full, rule = 2)$y))
        Kw_full <- as.numeric(gauss_K(u_full, IP_u, BW))
        da_full <- delta_base * (1 + KAPPA * Kw_full)

        df_rf <- data.frame(
            u    = u_full,
            g    = g_full,
            a_up = pmin(g_full + da_full, 1),
            a_lo = pmax(g_full - da_full, 0),
            f_up = pmin(g_full + delta_base, 1),
            f_lo = pmax(g_full - delta_base, 0)
        )

        adapt_half <- BW * sqrt(-2 * log(0.02))
        u_adapt_lo <- max(0, IP_u - adapt_half)
        u_adapt_hi <- min(1, IP_u + adapt_half)

        # == Gap zone =========================================================
        # Stage 1 is band proximity only: |v - g(u)| <= delta_adapt.
        pass_adapt <- pass_stage1
        pass_flat  <- pass_stage1 & (resid_abs <= delta_base)
        in_gap     <- pass_adapt & !pass_flat

        df$cell_class <- factor(as.character(df$cell_class),
                                levels = class_levels)
        if (anyNA(df$cell_class))
            stop(sum(is.na(df$cell_class)), " cell(s) carry a cell_class ",
                 "outside the four CMRF classes.", call. = FALSE)

        df$gap_fate <- NA_character_
        df$gap_fate[in_gap] <- ifelse(
            df$cell_class[in_gap] == "Recovered",  "Gap -> Recovered",
            ifelse(df$cell_class[in_gap] == "Excluded",
                   "Gap -> Excluded (Stage 2)", "Gap -> Other"))

        flat_cand    <- !df$is_orig_skcs & pass_stage1 & (resid_abs <= delta_base)
        df_flat_cand <- df[flat_cand, ]
        n_flat_cand  <- nrow(df_flat_cand)
        BLUE_FILL    <- "#2171B5"

        df$plot_col <- ifelse(df$group_plot %in% names(FOUR_COLS),
                              df$group_plot, "Non-selected")
        df$plot_col <- factor(df$plot_col, levels = ALL_LEVS)

        # Every panel is coloured from group_orig (the Part 10b SKCS-core labels).
        df$plot_col_orig <- ifelse(
            as.character(df$group_orig) %in% names(FOUR_COLS),
            as.character(df$group_orig), "Non-selected")
        df$plot_col_orig <- factor(df$plot_col_orig, levels = ALL_LEVS)
        n_core_col <- sum(df$plot_col_orig != "Non-selected")
        n_ref_col  <- sum(df$plot_col      != "Non-selected")
        n_col_diff <- n_ref_col - n_core_col
        cat(sprintf(paste0("    Coloured cells: %d under group_orig (the SKCS",
                           " core, used on all 11 panels) | %d under %s",
                           " (loaded, not drawn) | difference %d\n"),
                    n_core_col, n_ref_col, .post_col, n_col_diff))

        # Spine training set (TREP SKCS-core cells, both sides of the IP), from
        # group_orig; which() avoids NA phantom rows.
        df$spine_train_grp <- ifelse(
            as.character(df$group_orig) %in% c("Pre-IP TREP", "Post-IP TREP"),
            as.character(df$group_orig), NA_character_)
        n_spine_pre   <- sum(df$spine_train_grp == "Pre-IP TREP",  na.rm = TRUE)
        n_spine_post  <- sum(df$spine_train_grp == "Post-IP TREP", na.rm = TRUE)
        n_spine_train <- n_spine_pre + n_spine_post
        spine_n_10e   <- suppressWarnings(as.numeric(r$spine_n_cells)[1])
        if (!is.finite(spine_n_10e)) spine_n_10e <- NA_real_
        cat(sprintf(paste0("    Spine training set: %d TREP SKCS-core cells",
                           " (Pre-IP %d, Post-IP %d) | Part 10e spine_n_cells:",
                           " %s\n"),
                    n_spine_train, n_spine_pre, n_spine_post,
                    if (is.finite(spine_n_10e))
                        format(as.integer(spine_n_10e)) else "not recorded"))

        # df_all_*: every cell. df_std_*: every cell except the non-SKCS gap
        # cells, which panels with a gap symbol draw as dark circles.
        df_all_gray <- df[df$plot_col_orig == "Non-selected", ]
        df_all_col  <- df[df$plot_col_orig != "Non-selected", ]
        .std_keep   <- !(in_gap & !df$is_orig_skcs)
        df_std_gray <- df[.std_keep & df$plot_col_orig == "Non-selected", ]
        df_std_col  <- df[.std_keep & df$plot_col_orig != "Non-selected", ]
        df_gap         <- df[in_gap, ]

        # The Stage 1 candidate pool is non-SKCS by definition.
        df_gap_cand   <- df[in_gap & !df$is_orig_skcs, ]
        n_gap_cand    <- nrow(df_gap_cand)

        df_gap_skcs   <- df[in_gap & df$is_orig_skcs, ]
        n_gap_skcs    <- nrow(df_gap_skcs)
        GAP_FILL      <- "#2C2C2C"   # non-SKCS gap: the Stage 1 candidates
        # SKCS-core gap cells are drawn as ordinary four-group cells.

        n_gap_all <- nrow(df_gap)
        gap_fate_table <- if (n_gap_all > 0) table(df_gap$gap_fate) else integer(0)
        fate_subtitle <- if (n_gap_all > 0)
            paste(paste0(names(gap_fate_table), ": n=",
                         as.integer(gap_fate_table)), collapse = "   |   ")
        else "No gap-zone cells"

        ZOOM_HALF <- 0.25
        u_lo_zoom <- max(0, IP_u - ZOOM_HALF)
        u_hi_zoom <- min(1, IP_u + ZOOM_HALF)

        # == Five-class legend (explicit limits) ==============================
        legend_scale <- scale_colour_manual(
            values = ALL_COLS, name = "Cell group",
            limits = ALL_LEVS, drop = FALSE,
            guide  = guide_legend(override.aes = list(
                size  = c(1.8, rep(2.2, 4)),
                shape = rep(16L, 5),
                alpha = c(0.65, rep(0.75, 4)))))

        # == Output path helper ===============================================
        bf_path <- function(tag)
            file.path(OUT_DIR, paste0(FILE_PREFIX, "_", tag, "_", ds))

        # == Caption carried by every figure ==================================
        band_note_fn <- function(lab_src) sprintf(paste0(
            "delta_base = %.5f, a donor-weighted %.0f%% quantile of |v - g(u)|",
            " over the original SKCS-core.\n",
            "Spine = mgcv penalised thin-plate spline, s(u, bs='tp') +",
            " s(donor, bs='re'), inverse-donor weights, REML; drawn as the",
            " population curve (exclude = s(donor)). No shape-constrained",
            " basis is used; monotonicity is verified per fit: %s.\n",
            "Subgroup labels: %s.  Arm: %s.\n%s"),
            delta_base, 100 * CMRF_COVERAGE, mono_txt, lab_src, ARM_TAG,
            vip_note)
        group_src_orig <- sprintf(paste0("Part 10e group_orig, the SKCS-core",
                                         " labels Part 10b assigned (%d core",
                                         " cell(s); pre-CMRF)"), n_core_col)
        band_note      <- band_note_fn(group_src_orig)   # all eleven panels
        cat("    ", band_note, "\n", sep = "")

        cat(sprintf(paste0("    Gap zone %d = non-SKCS %d (dark fill) +",
                           " SKCS-core %d (drawn as ordinary core cells) | Stage 1 flat",
                           " candidates %d\n"),
                    n_gap_all, n_gap_cand, n_gap_skcs, n_flat_cand))


        # ===================================================================
        # BF01: SPINE ONLY
        # ===================================================================
        cat("  BF01: Spine only ...\n")

        # g(IP_u) is evaluated on the fitted spine; v = 0.5 is the cDFG median.
        v_spine_at_ip <- as.numeric(approx(u_full, g_full, xout = IP_u,
                                           rule = 2)$y)
        goi_at_ip     <- as.numeric(approx(df$u_all, df$GOI_expr, xout = IP_u,
                                           rule = 2)$y)
        if (!is.finite(v_spine_at_ip))
            stop("BF01: g(IP_u) is not finite for ", ds, call. = FALSE)
        cat(sprintf(paste0("    g(IP_u) = %.6f | cDFG median line v = 0.5",
                           " | deviation %+.6f\n"),
                    v_spine_at_ip, v_spine_at_ip - 0.5))
    
        p_bf01 <- ggplot() +
            # All cells as gray background
            geom_point(data = df,
                       aes(x = u_all, y = v_all),
                       colour = "gray80", size = 0.5, alpha = 0.35) +
            # Pre-IP and Post-IP TREP: the spine's training cells
            geom_point(data = df[which(df$spine_train_grp == "Pre-IP TREP"), ],
                       aes(x = u_all, y = v_all),
                       colour = unname(FOUR_COLS["Pre-IP TREP"]),
                       size = 1.25, alpha = 0.45) +
            geom_point(data = df[which(df$spine_train_grp == "Post-IP TREP"), ],
                       aes(x = u_all, y = v_all),
                       colour = unname(FOUR_COLS["Post-IP TREP"]),
                       size = 1.25, alpha = 0.45) +
            # spine
            geom_line(data = df_rf, aes(x = u, y = g),
                      colour = COL_SPINE, linewidth = 1.2) +
            # IP_v line, IP crosshairs at (IP_u, g(IP_u)), cDFG median in gray
            ip_hline_v(label = TRUE) +
            geom_vline(xintercept = IP_u, colour = COL_IP,
                       linewidth = 0.6, linetype = "dashed") +
            geom_hline(yintercept = v_spine_at_ip, colour = COL_IP,
                       linewidth = 0.5, linetype = "dashed") +
            geom_hline(yintercept = 0.5, colour = "gray55",
                       linewidth = 0.4, linetype = "dotted") +
            # IP point and label
            annotate("point",
                     x = IP_u, y = v_spine_at_ip,
                     shape = 21, size = 4.5, fill = "white",
                     colour = COL_IP, stroke = 1.2) +
            annotate("text",
                     x = IP_u + 0.02, y = max(0.06, v_spine_at_ip - 0.12),
                     label = sprintf(
                         paste0("IP_u = %.4f\ng(IP_u) = %.4f\n",
                                "cDFG median v = 0.5 (offset %+.4f)\n",
                                "GOI value = %.4f"),
                         IP_u, v_spine_at_ip, v_spine_at_ip - 0.5, goi_at_ip),
                     size = 3.0, hjust = 0, fontface = "bold",
                     colour = "gray15", lineheight = 0.9) +
            # Spine label
            annotate("text",
                     x = 0.72, y = g_full[which.min(abs(u_full - 0.72))] + 0.05,
                     label = "g(u) = E[v | u, TREP SKCS-core]",
                     size = 2.8, colour = COL_SPINE, fontface = "italic") +
            # v = 0.5 label
            annotate("text", x = 0.02, y = 0.53,
                     label = "v = 0.5\n[cDFG median]",
                     size = 2.5, colour = "gray45", fontface = "italic",
                     hjust = 0, lineheight = 0.85) +
            coord_cartesian(xlim = c(0, 1), ylim = c(-0.04, 1.04)) +
            scale_y_continuous(expand = expansion(mult = 0.01)) +
            labs(
                title = sprintf(
                    "BF01: Spine g(u) -- %s", ds),
                subtitle = paste0(
                    sprintf(paste0("Gray = all cells.  Amber = Pre-IP TREP (%d),",
                                   " green = Post-IP TREP (%d): together the %d",
                                   " TREP SKCS-core cells the spine is fitted",
                                   " on%s.  "),
                            n_spine_pre, n_spine_post, n_spine_train,
                            if (is.finite(spine_n_10e) &&
                                spine_n_10e != n_spine_train)
                                sprintf(" [Part 10e reports %d]", spine_n_10e)
                            else ""),
                    "Teal = penalised spline g(u); monotonicity verdict in the caption.\n",
                    sprintf(paste0("Navy crosshairs: (IP_u, g(IP_u)).  ",
                                   "IP_u=%.6f, g(IP_u)=%.6f.  ",
                                   "Gray dotted line is the cDFG marginal ",
                                   "median v=0.5 (offset %+.6f); it is not ",
                                   "the IP."),
                            IP_u, v_spine_at_ip, v_spine_at_ip - 0.5)),
                x = "u = Fx(GOI)  [GOI marginal CDF]",
                y = "v = Fc(cDFG)  [cDFG marginal CDF]"
            ) +
            th_band(10) +
            theme(legend.position = "none",
                  plot.subtitle = element_text(size = 7.5, lineheight = 1.2))
    
        # Caption: band provenance and spine monotonicity.
        p_bf01 <- p_bf01 +
            labs(caption = band_note) +
            theme(plot.caption = element_text(size = 6.2, colour = "gray35",
                                              hjust = 0, lineheight = 1.1))
    
        save_pp(p_bf01,
                bf_path("BF01_spine_only"),
                pdf_w = 11, pdf_h = 9, png_w = 2200, png_h = 1800)
        show_plot(p_bf01)
        cat("  BF01 saved\n")
    
    
        # ===================================================================
        # BF02: SPINE + 5-CLASS CELLS (NO BAND)
        # ===================================================================
        cat("  BF02: Spine + 5-class cells (no band) ...\n")
    
        p_bf02 <- ggplot() +
            # == Gray non-selected cells (background) =========================
        geom_point(data = df_all_gray,
                   aes(x = u_all, y = v_all, colour = plot_col_orig),
                   size = 0.55, alpha = 0.15) +
            # == 4-group coloured cells (all, no gap distinction) =============
        geom_point(data = df_all_col,
                   aes(x = u_all, y = v_all, colour = plot_col_orig),
                   size = 1.25, alpha = 0.55) +
            # == spine ========================================================
        geom_line(data = df_rf, aes(x = u, y = g),
                  colour = COL_SPINE, linewidth = 1.1) +
            # == IP_u reference line ==========================================
        # IP_v line (the Pre/Post boundary)
        ip_hline_v(label = FALSE) +
        geom_vline(xintercept = IP_u, colour = COL_IP,
                   linewidth = 0.6, linetype = "dashed") +
            annotate("text", x = IP_u, y = -0.01,
                     label = "IP_u", size = 2.8, colour = COL_IP,
                     fontface = "bold") +
            # == Legend =======================================================
        legend_scale +
            coord_cartesian(xlim = c(0, 1), ylim = c(-0.04, 1.04)) +
            scale_y_continuous(expand = expansion(mult = 0.01)) +
            labs(
                title = sprintf(
                    "BF02: Spine + 5-class cells (no tolerance band) -- %s", ds),
                subtitle = paste0(
                    "All GOI+ cells assigned to their IP x TREP subgroup ",
                    "(Pre/Post-IP x TREP/non-TREP) or Non-selected (gray).\n",
                    "Teal = penalised spline g(u); monotonicity verdict in the caption.  ",
                    "No tolerance band shown."),
                x = "u = Fx(GOI)",
                y = "v = Fc(cDFG)"
            ) +
            th_band(10) +
            theme(legend.position = "right",
                  legend.text = element_text(size = 8),
                  plot.subtitle = element_text(size = 7.5, lineheight = 1.2))
    
        # Caption: band provenance and spine monotonicity.
        p_bf02 <- p_bf02 +
            labs(caption = band_note) +
            theme(plot.caption = element_text(size = 6.2, colour = "gray35",
                                              hjust = 0, lineheight = 1.1))
    
        save_pp(p_bf02,
                bf_path("BF02_cells_only"),
                pdf_w = 11, pdf_h = 9, png_w = 2200, png_h = 1800)
        show_plot(p_bf02)
        cat("  BF02 saved\n")
    
    
        # ===================================================================
        # BF03: SPINE + FLAT BAND + 5-CLASS CELLS
        # ===================================================================
        cat("  BF03: Flat band + 5-class cells ...\n")
    
        p_bf03 <- ggplot() +
            # == Flat (delta_base) band ribbon (blue) =========================
        geom_ribbon(data = df_rf, aes(x = u, ymin = f_lo, ymax = f_up),
                    fill = COL_FLAT, alpha = 0.18) +
            # == Gray non-selected cells (background) =========================
        geom_point(data = df_all_gray,
                   aes(x = u_all, y = v_all, colour = plot_col_orig),
                   size = 0.55, alpha = 0.15) +
            # == 4-group coloured cells (all, no gap distinction) =============
        geom_point(data = df_all_col,
                   aes(x = u_all, y = v_all, colour = plot_col_orig),
                   size = 1.25, alpha = 0.55) +
            # == spine ========================================================
        geom_line(data = df_rf, aes(x = u, y = g),
                  colour = COL_SPINE, linewidth = 1.1) +
            # == Flat band boundary lines (blue dashed) =======================
        geom_line(data = df_rf, aes(x = u, y = f_up),
                  colour = COL_FLAT, linewidth = 0.5, linetype = "dashed") +
            geom_line(data = df_rf, aes(x = u, y = f_lo),
                      colour = COL_FLAT, linewidth = 0.5, linetype = "dashed") +
            # == IP_u reference line ==========================================
        # IP_v line (the Pre/Post boundary)
        ip_hline_v(label = FALSE) +
        geom_vline(xintercept = IP_u, colour = COL_IP,
                   linewidth = 0.6, linetype = "dashed") +
            annotate("text", x = IP_u, y = -0.01,
                     label = "IP_u", size = 2.8, colour = COL_IP,
                     fontface = "bold") +
            # == Band annotation ==============================================
        annotate("text",
                 x = 0.78, y = 0.15,
                 label = sprintf("delta_base = %.5f\n(flat band, kappa = 0)",
                                 delta_base),
                 size = 2.6, colour = COL_FLAT, hjust = 0,
                 lineheight = 0.85) +
            # == Legend =======================================================
        legend_scale +
            coord_cartesian(xlim = c(0, 1), ylim = c(-0.04, 1.04)) +
            scale_y_continuous(expand = expansion(mult = 0.01)) +
            labs(
                title = sprintf(
                    "BF03: Flat tolerance band (delta_base) + 5-class cells -- %s",
                    ds),
                subtitle = sprintf(
                    paste0("Blue ribbon = flat band: g(u) +/- delta_base (kappa=0).  ",
                           "delta_base = %.5f: a donor-weighted 95th-percentile of |v - g(u)| over the original SKCS-core.\n",
                           "SKCS-core cells (pre-CMRF) in their IP x TREP subgroup ",
                           "(Pre/Post-IP x TREP/non-TREP); all other GOI+ cells gray."),
                    delta_base),
                x = "u = Fx(GOI)",
                y = "v = Fc(cDFG)"
            ) +
            th_band(10) +
            theme(legend.position = "right",
                  legend.text = element_text(size = 8),
                  plot.subtitle = element_text(size = 7.5, lineheight = 1.2))
    
        # Caption: band provenance and spine monotonicity.
        p_bf03 <- p_bf03 +
            labs(caption = band_note) +
            theme(plot.caption = element_text(size = 6.2, colour = "gray35",
                                              hjust = 0, lineheight = 1.1))
    
        save_pp(p_bf03,
                bf_path("BF03_flat_band"),
                pdf_w = 11, pdf_h = 9, png_w = 2200, png_h = 1800)
        show_plot(p_bf03)
        cat("  BF03 saved\n")
    
    
        # ===================================================================
        # BF04: ADD ADAPTIVE BAND (NO GAP CIRCLES)
        # ===================================================================
        cat("  BF04: Adaptive band (no gap circles) ...\n")
    
        p_bf04 <- ggplot() +
            # == Flat band ribbon (blue, behind) ==============================
        geom_ribbon(data = df_rf, aes(x = u, ymin = f_lo, ymax = f_up),
                    fill = COL_FLAT, alpha = 0.18) +
            # == Adaptive band ribbon (vermillion, in front) ==================
        geom_ribbon(data = df_rf, aes(x = u, ymin = a_lo, ymax = a_up),
                    fill = COL_ADAPT, alpha = 0.18) +
            # == Gray non-selected cells (background) =========================
        geom_point(data = df_all_gray,
                   aes(x = u_all, y = v_all, colour = plot_col_orig),
                   size = 0.55, alpha = 0.15) +
            # == 4-group coloured cells (foreground, ALL cells) ===============
        geom_point(data = df_all_col, 
                   aes(x = u_all, y = v_all, colour = plot_col_orig),
                   size = 1.25, alpha = 0.55) +
            # == spine ========================================================
        geom_line(data = df_rf, aes(x = u, y = g),
                  colour = COL_SPINE, linewidth = 1.1) +
            # == Adaptive band boundary lines (vermillion solid, darker) ======
            geom_line(data = df_rf, aes(x = u, y = a_up),
                      colour = COL_ADAPT_BORDER, linewidth = 0.6,
                      linetype = "solid") +
            geom_line(data = df_rf, aes(x = u, y = a_lo),
                      colour = COL_ADAPT_BORDER, linewidth = 0.6,
                      linetype = "solid") +
            # == Flat band boundary lines (blue dashed) =======================
        geom_line(data = df_rf, aes(x = u, y = f_up),
                  colour = COL_FLAT, linewidth = 0.45, linetype = "dashed") +
            geom_line(data = df_rf, aes(x = u, y = f_lo),
                      colour = COL_FLAT, linewidth = 0.45, linetype = "dashed") +
            # == Adaptive band start/end vertical dashed lines ================
        geom_vline(xintercept = c(u_adapt_lo, u_adapt_hi),
                   colour = "gray15", linewidth = 0.45, linetype = "dashed") +
            # == IP_u reference line ==========================================
        # IP_v line (the Pre/Post boundary)
        ip_hline_v(label = FALSE) +
        geom_vline(xintercept = IP_u, colour = COL_IP,
                   linewidth = 0.6, linetype = "dashed") +
            annotate("text", x = IP_u, y = -0.01,
                     label = "IP_u", size = 2.8, colour = COL_IP,
                     fontface = "bold") +
            # == Band annotation labels =======================================
        annotate("text",
                 x = IP_u + 0.02, y = 0.12,
                 label = sprintf(
                     "Adaptive band\nkappa = %.2f, bw = %.2f\ndelta_max = %.5f",
                     KAPPA, BW, delta_base * (1 + KAPPA)),
                 size = 2.5, colour = COL_ADAPT_BORDER, hjust = 0,
                 lineheight = 0.85) +
            annotate("label",
                     x = 0.80, y = 0.18,
                     label = sprintf("Adaptive bound\n(kappa=%.1f)", KAPPA),
                     size = 2.3, colour = COL_ADAPT_BORDER, fill = "white", alpha = 0.8, lineheight = 0.85) +
            annotate("label",
                     x = 0.80, y = 0.08,
                     label = sprintf("Flat bound\n(kappa=0)"),
                     size = 2.3, colour = COL_FLAT, fill = "white", alpha = 0.8, lineheight = 0.85) +
            # == Legend =======================================================
        legend_scale +
            coord_cartesian(xlim = c(0, 1), ylim = c(-0.04, 1.04)) +
            scale_y_continuous(expand = expansion(mult = 0.01)) +
            labs(
                title = sprintf(
                    "BF04: Adaptive band overlay (no gap-cell markers) -- %s", ds),
                subtitle = sprintf(
                    paste0("Vermillion = adaptive band: delta_adapt(u) = delta_base * ",
                           "(1 + kappa*K(u)).  Blue = flat band (kappa=0).\n",
                           "Coloured cells are the SKCS core (pre-CMRF); all ",
                           "other GOI+ cells gray.  ",
                           "kappa = %.2f, bw = %.2f.  ",
                           "Adaptive band active zone (K > 0.02): u in [%.3f, %.3f]  ",
                           "(= IP_u +/- %.3f)."),
                    KAPPA, BW, u_adapt_lo, u_adapt_hi, adapt_half),
                x = "u = Fx(GOI)",
                y = "v = Fc(cDFG)"
            ) +
            th_band(10) +
            theme(legend.position = "right",
                  legend.text = element_text(size = 8),
                  plot.subtitle = element_text(size = 7.5, lineheight = 1.2))
    
        # Caption: band provenance and spine monotonicity.
        p_bf04 <- p_bf04 +
            labs(caption = band_note) +
            theme(plot.caption = element_text(size = 6.2, colour = "gray35",
                                              hjust = 0, lineheight = 1.1))
    
        save_pp(p_bf04,
                bf_path("BF04_adaptive_band"),
                pdf_w = 11, pdf_h = 9, png_w = 2200, png_h = 1800)
        show_plot(p_bf04)
        cat("  BF04 saved\n")
    
    
        # ===================================================================
        # BF05: FULL BAND PLOT WITH GAP CIRCLES
        # (matches Part10e Plot 05 exactly)
        # ===================================================================
        cat("  BF05: Full band + gap circles ...\n")
    
        p_bf05 <- ggplot() +
            # == Ribbons (full range) =========================================
        # flat band (blue, behind)
        geom_ribbon(data = df_rf, aes(x = u, ymin = f_lo, ymax = f_up),
                    fill = COL_FLAT, alpha = 0.18) +
            # adaptive band (vermillion, in front)
            geom_ribbon(data = df_rf, aes(x = u, ymin = a_lo, ymax = a_up),
                        fill = COL_ADAPT, alpha = 0.18) +
            # == Gray non-gap cells (background layer) ========================
        geom_point(data = df_std_gray,
                   aes(x = u_all, y = v_all, colour = plot_col_orig),
                   size = 0.55, alpha = 0.15) +
            # == 4-group coloured non-gap cells (foreground layer) ============
        geom_point(data = df_std_col, 
                   aes(x = u_all, y = v_all, colour = plot_col_orig),
                   size = 1.25, alpha = 0.55) +
            # == spine (full range) ===========================================
        geom_line(data = df_rf, aes(x = u, y = g),
                  colour = COL_SPINE, linewidth = 1.1) +
            # == Flat-band Stage 1 candidates: blue circles ===================
        {
            if (n_flat_cand > 0)
                geom_point(data = df_flat_cand,
                           aes(x = u_all, y = v_all),
                           shape = 21, size = 3.0,
                           fill = BLUE_FILL, colour = "white",
                           stroke = 0.5, alpha = 0.88,
                           inherit.aes = FALSE)
            else
                geom_blank()
        } +
            # == Non-SKCS gap cells: dark circles =============================
            {
                if (n_gap_cand > 0)
                    geom_point(data = df_gap_cand,
                               aes(x = u_all, y = v_all),
                               shape = 21, size = 3.0,
                               fill = GAP_FILL, colour = "white",
                               stroke = 0.5, alpha = 0.88,
                               inherit.aes = FALSE)
                else
                    geom_blank()
            } +
            # == Band boundary lines (full range) =============================
        geom_line(data = df_rf, aes(x = u, y = a_up),
                  colour = COL_ADAPT, linewidth = 0.55, linetype = "solid") +
            geom_line(data = df_rf, aes(x = u, y = a_lo),
                      colour = COL_ADAPT, linewidth = 0.55, linetype = "solid") +
            geom_line(data = df_rf, aes(x = u, y = f_up),
                      colour = COL_FLAT, linewidth = 0.45, linetype = "dashed") +
            geom_line(data = df_rf, aes(x = u, y = f_lo),
                      colour = COL_FLAT, linewidth = 0.45, linetype = "dashed") +
            # == Adaptive band start/end vertical dashed lines ================
        geom_vline(xintercept = c(u_adapt_lo, u_adapt_hi),
                   colour = "gray15", linewidth = 0.45, linetype = "dashed") +
            # == IP_u reference line ==========================================
        # IP_v line (the Pre/Post boundary)
        ip_hline_v(label = FALSE) +
        geom_vline(xintercept = IP_u, colour = COL_IP,
                   linewidth = 0.6, linetype = "dashed") +
            # == Annotations ==================================================
        annotate("text", x = IP_u, y = -0.01,
                 label = "IP_u", size = 2.8, colour = COL_IP,
                 fontface = "bold") +
            # == Zoom window boundary markers (subtle) ========================
        geom_vline(xintercept = c(u_lo_zoom, u_hi_zoom),
                   colour = "gray60", linewidth = 0.3, linetype = "dotted") +
            # == Gap zone annotation ==========================================
        {
            if (n_gap_all > 0)
                annotate("label",
                         x = IP_u + 0.04,
                         y = g_full[which.min(abs(u_full - IP_u))] + 0.15,
                         label = sprintf(paste0("Gap zone (n=%d)\ndark %d non-SKCS",
                                                " | %d SKCS-core, as ordinary cells\nMissed by flat,",
                                                " captured by adaptive"),
                                         n_gap_all, n_gap_cand, n_gap_skcs),
                         size = 2.4, fill = alpha("white", 0.85), colour = "gray20",
                         lineheight = 0.85)
            else
                geom_blank()
        } +
            # == Legend =======================================================
        legend_scale +
            coord_cartesian(xlim = c(0, 1), ylim = c(-0.04, 1.04)) +
            scale_y_continuous(expand = expansion(mult = 0.01)) +
            labs(
                title = sprintf(
                    "BF05: Full range + gap-zone cells  u in [0, 1]  -- %s", ds),
                subtitle = sprintf(
                    paste0("Vermillion = adaptive band; blue = flat band.  ",
                           "Blue circles (n=%d) = Stage 1 candidates within delta_base.  ",
                           "Gap cells (n=%d) = missed by flat, captured by adaptive: ",
                           "dark fill (n=%d) non-SKCS, the Stage 1 candidates BF08 draws; ",
                           "the other %d are already SKCS-core and keep the ordinary four-group dot.\n",
                           "These cells sit in the orange-only zone between ",
                           "flat and adaptive boundaries.  ",
                           "Their Stage 2 downstream fates: %s\n",
                           "Adaptive band active zone (K > 0.02): u in [%.3f, %.3f]  ",
                           "(= IP_u +/- 2.80*bw = IP_u +/- %.3f).  ",
                           "Zoom window [%.2f, %.2f] shown by gray dotted lines."),
                    n_flat_cand, n_gap_all, n_gap_cand, n_gap_skcs, fate_subtitle,
                    u_adapt_lo, u_adapt_hi,
                    adapt_half,
                    u_lo_zoom, u_hi_zoom),
                x = "u = Fx(GOI)",
                y = "v = Fc(cDFG)"
            ) +
            th_band(10) +
            theme(legend.position = "right",
                  legend.text = element_text(size = 8),
                  plot.subtitle = element_text(size = 7.5, lineheight = 1.2))
    
        # Caption: band provenance and spine monotonicity.
        p_bf05 <- p_bf05 +
            labs(caption = band_note) +
            theme(plot.caption = element_text(size = 6.2, colour = "gray35",
                                              hjust = 0, lineheight = 1.1))
    
        save_pp(p_bf05,
                bf_path("BF05_full_band"),
                pdf_w = 11, pdf_h = 9, png_w = 2200, png_h = 1800)
        show_plot(p_bf05)
        cat("  BF05 saved\n")
    
    
        # ===================================================================
        # BF06: WAIST ZOOM
        # Zoomed to [IP_u - ZOOM_HALF, IP_u + ZOOM_HALF]
        # Same colour scheme as BF05 (Plot 05 colours)
        # ===================================================================
        cat("  BF06: Waist zoom ...\n")
    
        # == bw footprint markers for zoom view ===============================
        bw_marks <- c(IP_u - 2*BW, IP_u - BW, IP_u + BW, IP_u + 2*BW)
        bw_labels <- c("-2bw", "-bw", "+bw", "+2bw")
        # Only include markers within zoom window
        in_zoom_bw <- bw_marks >= u_lo_zoom & bw_marks <= u_hi_zoom
        bw_marks_z <- bw_marks[in_zoom_bw]
        bw_labels_z <- bw_labels[in_zoom_bw]
    
        # == Subset cells visible in zoom window ==============================
        df_zoom <- df[df$u_all >= (u_lo_zoom - 0.02) & df$u_all <= (u_hi_zoom + 0.02), ]
        # Zoom membership is carried on the barcode.
        df_zoom$is_gap <- df_zoom$barcode %in% df_gap$barcode
        # Zoom-window form of df_std_*.
        .zoom_std_keep <- !(df_zoom$is_gap & !df_zoom$is_orig_skcs)
        df_zoom_std_gray <- df_zoom[.zoom_std_keep &
                                        df_zoom$plot_col_orig == "Non-selected", ]
        df_zoom_std_col  <- df_zoom[.zoom_std_keep &
                                        df_zoom$plot_col_orig != "Non-selected", ]
        df_zoom_gap <- df_zoom[df_zoom$is_gap, ]
        # Split kept for the subtitle counts.
        df_zoom_gap_skcs <- df_zoom_gap[df_zoom_gap$is_orig_skcs, ]
        df_zoom_gap_cand <- df_zoom_gap[!df_zoom_gap$is_orig_skcs, ]
        n_gap_zoom_skcs  <- nrow(df_zoom_gap_skcs)
        n_gap_zoom_cand  <- nrow(df_zoom_gap_cand)
    
        n_gap_zoom <- nrow(df_zoom_gap)
        gap_fate_zoom_table <- if (n_gap_zoom > 0)
            table(df_zoom_gap$gap_fate) else integer(0)
        fate_zoom_subtitle <- if (n_gap_zoom > 0) {
            paste(paste0(names(gap_fate_zoom_table), ": n=",
                         as.integer(gap_fate_zoom_table)),
                  collapse = "   |   ")
        } else {
            "No gap-zone cells in zoom window"
        }
    
        p_bf06 <- ggplot() +
            # == Ribbons (full range, clipped by coord_cartesian) =============
        # flat band (blue, behind)
        geom_ribbon(data = df_rf, aes(x = u, ymin = f_lo, ymax = f_up),
                    fill = COL_FLAT, alpha = 0.18) +
            # adaptive band (vermillion, in front)
            geom_ribbon(data = df_rf, aes(x = u, ymin = a_lo, ymax = a_up),
                        fill = COL_ADAPT, alpha = 0.18) +
            # == Gray non-gap cells (background) ==============================
        geom_point(data = df_zoom_std_gray,
                   aes(x = u_all, y = v_all, colour = plot_col_orig),
                   size = 0.8, alpha = 0.15) +
            # == 4-group coloured non-gap cells (foreground) ==================
        geom_point(data = df_zoom_std_col, 
                   aes(x = u_all, y = v_all, colour = plot_col_orig),
                   size = 1.55, alpha = 0.55) +
            # == spine ========================================================
        geom_line(data = df_rf, aes(x = u, y = g),
                  colour = COL_SPINE, linewidth = 1.1) +
            # == Non-SKCS gap cells: dark circles =============================
        {
            if (n_gap_zoom_cand > 0)
                geom_point(data = df_zoom_gap_cand,
                           aes(x = u_all, y = v_all),
                           shape = 21, size = 3.5,
                           fill = GAP_FILL, colour = "white",
                           stroke = 0.5, alpha = 0.88,
                           inherit.aes = FALSE)
            else
                geom_blank()
        } +
            # == Adaptive band boundary lines (vermillion solid) ==============
        geom_line(data = df_rf, aes(x = u, y = a_up),
                  colour = COL_ADAPT, linewidth = 0.55, linetype = "solid") +
            geom_line(data = df_rf, aes(x = u, y = a_lo),
                      colour = COL_ADAPT, linewidth = 0.55, linetype = "solid") +
            # == Flat band boundary lines (blue dashed) =======================
        geom_line(data = df_rf, aes(x = u, y = f_up),
                  colour = COL_FLAT, linewidth = 0.45, linetype = "dashed") +
            geom_line(data = df_rf, aes(x = u, y = f_lo),
                      colour = COL_FLAT, linewidth = 0.45, linetype = "dashed") +
            # == IP_u reference line ==========================================
        # IP_v line (the Pre/Post boundary)
        ip_hline_v(label = FALSE) +
        geom_vline(xintercept = IP_u, colour = COL_IP,
                   linewidth = 0.6, linetype = "dashed") +
            # == v = 0.5 reference line =======================================
        geom_hline(yintercept = 0.5, colour = "gray40",
                   linewidth = 0.35, linetype = "dashed") +
            # == bw footprint markers (dotdash verticals) =====================
        {
            if (length(bw_marks_z) > 0)
                geom_vline(
                    data = data.frame(xint = bw_marks_z),
                    aes(xintercept = xint),
                    colour = "gray50", linewidth = 0.35,
                    linetype = "dotdash",
                    inherit.aes = FALSE)
            else
                geom_blank()
        } +
            # == bw footprint labels on x-axis ================================
        {
            if (length(bw_marks_z) > 0)
                annotate("text",
                         x = bw_marks_z, y = rep(-0.03, length(bw_marks_z)),
                         label = bw_labels_z,
                         size = 2.3, colour = COL_ADAPT, fontface = "italic")
            else
                geom_blank()
        } +
            # == IP_u label ===================================================
        annotate("text", x = IP_u, y = -0.05,
                 label = "IP_u", size = 2.8, colour = COL_IP,
                 fontface = "bold") +
            # == Gap zone annotation label ====================================
        {
            if (n_gap_zoom > 0)
                annotate("label",
                         x = IP_u + 0.04,
                         y = g_full[which.min(abs(u_full - IP_u))] + 0.18,
                         label = sprintf(
                             paste0("Gap zone (n=%d)\ndark %d non-SKCS |",
                                    " %d SKCS-core, as ordinary cells\nMissed by flat,",
                                    " captured by adaptive"),
                             n_gap_zoom, n_gap_zoom_cand, n_gap_zoom_skcs),
                         size = 2.4, fill = alpha("white", 0.85), colour = "gray20",
                         lineheight = 0.85)
            else
                geom_blank()
        } +
            # == Adaptive/flat bound labels ===================================
        annotate("label",
                 x = u_hi_zoom - 0.02,
                 y = approx(u_full, df_rf$a_lo,
                            xout = u_hi_zoom - 0.02, rule = 2)$y + 0.04,
                 label = sprintf("Adaptive bound\n(kappa=%.1f)", KAPPA),
                 size = 2.3, colour = COL_ADAPT, fill = "white", alpha = 0.8, lineheight = 0.85) +
            annotate("label",
                     x = u_hi_zoom - 0.02,
                     y = approx(u_full, df_rf$f_lo,
                                xout = u_hi_zoom - 0.02, rule = 2)$y - 0.04,
                     label = "Flat bound\n(kappa=0)",
                     size = 2.3, colour = COL_FLAT, fill = "white", alpha = 0.8, lineheight = 0.85) +
            # == Legend =======================================================
        legend_scale +
            coord_cartesian(xlim = c(u_lo_zoom, u_hi_zoom),
                            ylim = c(-0.06, 1.06)) +
            scale_y_continuous(expand = expansion(mult = 0.01)) +
            labs(
                title = sprintf(
                    "BF06: Waist Zoom  u in [%.2f, %.2f]  -- %s",
                    u_lo_zoom, u_hi_zoom, ds),
                subtitle = sprintf(
                    paste0("Vermillion = adaptive band; blue = flat band.  ",
                           "Gap cells in window (n=%d) = missed by flat, captured by adaptive: ",
                           "dark fill (n=%d) non-SKCS; the other %d are SKCS-core and keep the ordinary four-group dot.\n",
                           "These cells sit in the orange-only zone between ",
                           "flat and adaptive boundaries.  ",
                           "Their Stage 2 downstream fates: %s\n",
                           "Adaptive band active zone (K > 0.02): u in [%.3f, %.3f]  ",
                           "(= IP_u +/- 2.80*bw = IP_u +/- %.3f)"),
                    n_gap_zoom, n_gap_zoom_cand, n_gap_zoom_skcs,
                    fate_zoom_subtitle,
                    u_adapt_lo, u_adapt_hi, adapt_half),
                x = "u  =  Fx(GOI)",
                y = "v = Fc(cDFG)"
            ) +
            th_band(10) +
            theme(legend.position = "right",
                  legend.text = element_text(size = 8),
                  plot.subtitle = element_text(size = 7.5, lineheight = 1.2))
    
        # Caption: band provenance and spine monotonicity.
        p_bf06 <- p_bf06 +
            labs(caption = band_note) +
            theme(plot.caption = element_text(size = 6.2, colour = "gray35",
                                              hjust = 0, lineheight = 1.1))
    
        save_pp(p_bf06,
                bf_path("BF06_waist_zoom"),
                pdf_w = 11, pdf_h = 9, png_w = 2200, png_h = 1800)
        show_plot(p_bf06)
        cat("  BF06 saved\n")
    
    
        # ===================================================================
        # BF07: WAIST ZOOM WITH QC1-QC3 ANNOTATIONS
        # Gap cells rendered as gray (no special treatment).
        # QC1: max band width at IP_u (K=1)
        # QC2: band width at IP_u +/- bw (K=0.6065)
        # QC3: fastest widening/narrowing rate at IP_u +/- bw (|dK/du|=6.065)
        # ===================================================================
        cat("  BF07: Waist zoom with QC annotations ...\n")
    
        # == Compute QC values at key u positions =============================
        g_at_uIP    <- approx(u_full, g_full, xout = IP_u, rule = 2)$y
        g_at_uIP_m  <- approx(u_full, g_full, xout = IP_u - BW, rule = 2)$y
        g_at_uIP_p  <- approx(u_full, g_full, xout = IP_u + BW, rule = 2)$y
    
        # Band boundaries at IP_u (K=1, max width)
        da_at_uIP   <- delta_base * (1 + KAPPA * 1.0)   # delta_max
        a_up_uIP    <- min(g_at_uIP + da_at_uIP, 1)
        a_lo_uIP    <- max(g_at_uIP - da_at_uIP, 0)
    
        # Band boundaries at IP_u - bw (K = e^{-0.5} = 0.6065)
        K_at_1bw    <- exp(-0.5)   # 0.6065
        da_at_m1bw  <- delta_base * (1 + KAPPA * K_at_1bw)
        a_up_m1bw   <- min(g_at_uIP_m + da_at_m1bw, 1)
        a_lo_m1bw   <- max(g_at_uIP_m - da_at_m1bw, 0)
    
        # Band boundaries at IP_u + bw
        da_at_p1bw  <- delta_base * (1 + KAPPA * K_at_1bw)
        a_up_p1bw   <- min(g_at_uIP_p + da_at_p1bw, 1)
        a_lo_p1bw   <- max(g_at_uIP_p - da_at_p1bw, 0)
    
        # |dK/du| at +/-bw = 1/(bw*sqrt(e)) = 6.065
        dKdu_peak   <- 1 / (BW * sqrt(exp(1)))   # 6.065
    
        # == All cells as regular (gap cells become gray) =====================
        df_qc <- df
        df_qc$plot_col_qc <- df_qc$plot_col_orig
        df_qc_gray <- df_qc[df_qc$plot_col_qc == "Non-selected", ]
        df_qc_col  <- df_qc[df_qc$plot_col_qc != "Non-selected", ]
    
        # == Arrow specifications =============================================
        arrow_both <- grid::arrow(angle = 20, length = unit(0.08, "inches"),
                                  ends = "both", type = "closed")
    
        # == QC3: Tangent slopes on adaptive boundary at +/-bw ================
        eps <- 0.005
        a_up_slope_m <- (approx(df_rf$u, df_rf$a_up,
                                xout = IP_u - BW + eps, rule = 2)$y -
                             approx(df_rf$u, df_rf$a_up,
                                    xout = IP_u - BW - eps, rule = 2)$y) / (2*eps)
        a_up_slope_p <- (approx(df_rf$u, df_rf$a_up,
                                xout = IP_u + BW + eps, rule = 2)$y -
                             approx(df_rf$u, df_rf$a_up,
                                    xout = IP_u + BW - eps, rule = 2)$y) / (2*eps)
        a_lo_slope_m <- (approx(df_rf$u, df_rf$a_lo,
                                xout = IP_u - BW + eps, rule = 2)$y -
                             approx(df_rf$u, df_rf$a_lo,
                                    xout = IP_u - BW - eps, rule = 2)$y) / (2*eps)
        a_lo_slope_p <- (approx(df_rf$u, df_rf$a_lo,
                                xout = IP_u + BW + eps, rule = 2)$y -
                             approx(df_rf$u, df_rf$a_lo,
                                    xout = IP_u + BW - eps, rule = 2)$y) / (2*eps)
    
        tang_hu <- 0.04   # tangent half-length in u-units
    
        df_tang <- data.frame(
            x_start = c(IP_u - BW - tang_hu, IP_u + BW - tang_hu),
            x_end   = c(IP_u - BW + tang_hu, IP_u + BW + tang_hu),
            y_start = c(a_up_m1bw - a_up_slope_m * tang_hu,
                        a_up_p1bw - a_up_slope_p * tang_hu),
            y_end   = c(a_up_m1bw + a_up_slope_m * tang_hu,
                        a_up_p1bw + a_up_slope_p * tang_hu),
            loc     = c("-bw (upper)", "+bw (upper)")
        )
        df_tang_lo <- data.frame(
            x_start = c(IP_u - BW - tang_hu, IP_u + BW - tang_hu),
            x_end   = c(IP_u - BW + tang_hu, IP_u + BW + tang_hu),
            y_start = c(a_lo_m1bw - a_lo_slope_m * tang_hu,
                        a_lo_p1bw - a_lo_slope_p * tang_hu),
            y_end   = c(a_lo_m1bw + a_lo_slope_m * tang_hu,
                        a_lo_p1bw + a_lo_slope_p * tang_hu),
            loc     = c("-bw (lower)", "+bw (lower)")
        )
    
        # == QC annotation colours ============================================
        QC_COL1 <- "#8B0000"   # dark red for QC1
        QC_COL2 <- "#4A148C"   # deep purple for QC2
        QC_COL3 <- "#202020"   # dark gray for QC3 tangents
        qc2_xoff <- 0.008
    
        p_bf07 <- ggplot() +
            # == Ribbons ======================================================
        geom_ribbon(data = df_rf, aes(x = u, ymin = f_lo, ymax = f_up),
                    fill = COL_FLAT, alpha = 0.18) +
            geom_ribbon(data = df_rf, aes(x = u, ymin = a_lo, ymax = a_up),
                        fill = COL_ADAPT, alpha = 0.18) +
            # == Gray cells (very faint background) ===========================
        geom_point(data = df_qc_gray,
                   aes(x = u_all, y = v_all, colour = plot_col_qc),
                   size = 0.6, alpha = 0.15) +
            # == 4-group coloured cells (foreground) ==========================
        geom_point(data = df_qc_col,
                   aes(x = u_all, y = v_all, colour = plot_col_qc),
                   size = 1.55, alpha = 0.45) +
            # == spine ========================================================
        geom_line(data = df_rf, aes(x = u, y = g),
                  colour = COL_SPINE, linewidth = 1.1) +
            # == Band boundaries ==============================================
        geom_line(data = df_rf, aes(x = u, y = a_up),
                  colour = COL_ADAPT_BORDER, linewidth = 0.55,
                  linetype = "solid") +
            geom_line(data = df_rf, aes(x = u, y = a_lo),
                      colour = COL_ADAPT_BORDER, linewidth = 0.55,
                      linetype = "solid") +
            geom_line(data = df_rf, aes(x = u, y = f_up),
                      colour = COL_FLAT, linewidth = 0.45, linetype = "dashed") +
            geom_line(data = df_rf, aes(x = u, y = f_lo),
                      colour = COL_FLAT, linewidth = 0.45, linetype = "dashed") +
            # == Reference lines: IP_u and +/-bw only =========================
        # IP_v line (the Pre/Post boundary)
        ip_hline_v(label = FALSE) +
        geom_vline(xintercept = IP_u, colour = COL_IP,
                   linewidth = 0.5, linetype = "dashed") +
            geom_vline(xintercept = c(IP_u - BW, IP_u + BW),
                       colour = "gray50", linewidth = 0.35,
                       linetype = "dotdash") +
        
            # ==============================================================
        # QC1: maximum band width at IP_u (K = 1)
        # ==============================================================
        annotate("segment", x = IP_u, xend = IP_u,
                 y = a_lo_uIP, yend = a_up_uIP,
                 colour = QC_COL1, linewidth = 0.6, arrow = arrow_both) +
            ann_label_10f(x = IP_u, y = a_up_uIP,
                          label = sprintf("QC1  2*delta_max = %.4f", 2 * da_at_uIP),
                          size = 2.3, colour = QC_COL1, vjust = -0.15,
                          fill = alpha("white", 0.85)) +
            # ==============================================================
        # QC2: band width at IP_u +/- bw (K = e^-0.5 = 0.6065)
        # ==============================================================
        annotate("segment", x = IP_u - BW - qc2_xoff, xend = IP_u - BW - qc2_xoff,
                 y = a_lo_m1bw, yend = a_up_m1bw,
                 colour = QC_COL2, linewidth = 0.5, arrow = arrow_both) +
            annotate("segment", x = IP_u + BW + qc2_xoff, xend = IP_u + BW + qc2_xoff,
                     y = a_lo_p1bw, yend = a_up_p1bw,
                     colour = QC_COL2, linewidth = 0.5, arrow = arrow_both) +
            ann_label_10f(x = IP_u + BW + qc2_xoff, y = a_lo_p1bw,
                          label = sprintf("QC2  2*delta(IP_u+bw) = %.4f\nK = %.4f",
                                          2 * da_at_p1bw, K_at_1bw),
                          size = 2.2, colour = QC_COL2, vjust = 1.15, hjust = 0,
                          fill = alpha("white", 0.85), lineheight = 0.9) +
            # ==============================================================
        # QC3: Tangent lines at +/-bw on adaptive band boundaries
        # (steepest band slope = fastest widening/narrowing)
        # ==============================================================
        geom_segment(data = df_tang,
                     aes(x = x_start, xend = x_end,
                         y = y_start, yend = y_end),
                     colour = QC_COL3, linewidth = 0.4,
                     linetype = "solid",
                     inherit.aes = FALSE) +
            geom_segment(data = df_tang_lo,
                         aes(x = x_start, xend = x_end,
                             y = y_start, yend = y_end),
                         colour = QC_COL3, linewidth = 0.4,
                         linetype = "solid",
                         inherit.aes = FALSE) +
            # Contact dots at tangent points (upper and lower boundaries)
            annotate("point",
                     x = c(IP_u - BW, IP_u + BW),
                     y = c(a_up_m1bw, a_up_p1bw),
                     shape = 16, size = 2.5, colour = QC_COL3) +
            annotate("point",
                     x = c(IP_u - BW, IP_u + BW),
                     y = c(a_lo_m1bw, a_lo_p1bw),
                     shape = 16, size = 2.5, colour = QC_COL3) +
        
            # == Legend (explicit limits, as legend_scale) ====================
        scale_colour_manual(
            values = ALL_COLS,
            name   = "Cell group",
            limits = ALL_LEVS,
            drop   = FALSE,
            guide  = guide_legend(
                override.aes = list(
                    size  = c(1.8, rep(2.2, 4)),
                    shape = rep(16L, 5),
                    alpha = c(0.65, rep(0.75, 4))
                )
            )
        ) +
            coord_cartesian(xlim = c(u_lo_zoom, u_hi_zoom),
                            ylim = c(-0.04, 1.04)) +
            scale_y_continuous(expand = expansion(mult = 0.01)) +
            labs(
                title = sprintf(
                    "BF07: Waist Zoom with QC Tangents  u in [%.2f, %.2f]  -- %s",
                    u_lo_zoom, u_hi_zoom, ds),
                subtitle = sprintf(
                    paste0("Tangent lines at IP_u +/- bw: steepest adaptive band ",
                           "boundary slope (|dK/du| = %.3f).\n",
                           "QC1: K(IP_u) = 1 (max width at IP_u).  ",
                           "QC2: K(IP_u +/- bw) = 0.6065 (60.65%% of max widening).  ",
                           "QC3: |dK/du| peaks at +/-bw = %.3f.\n",
                           "delta_base = %.5f, delta_max = %.5f.  ",
                           "kappa = %.2f, bw = %.2f."),
                    dKdu_peak, dKdu_peak, delta_base, da_at_uIP, KAPPA, BW),
                x = "u  =  Fx(GOI)",
                y = "v = Fc(cDFG)"
            ) +
            th_band(10) +
            theme(legend.position = "right",
                  legend.text = element_text(size = 8),
                  plot.subtitle = element_text(size = 7.5, lineheight = 1.2))
    
        # Caption: band provenance and spine monotonicity.
        p_bf07 <- p_bf07 +
            labs(caption = band_note) +
            theme(plot.caption = element_text(size = 6.2, colour = "gray35",
                                              hjust = 0, lineheight = 1.1))
    
        save_pp(p_bf07,
                bf_path("BF07_QC_zoom"),
                pdf_w = 12, pdf_h = 9, png_w = 2400, png_h = 1800)
        show_plot(p_bf07)
        cat("  BF07 saved\n")
    
    
        # ===================================================================
        # BF08: PLOT #6 -- STAGE 1 CANDIDATE POOL
        # Blue + dark circles showing the full non-SKCS candidate pool
        # entering RSk-fCV.  Tan triangles deferred to BF09.
        # (Tan aesthetics defined here, used by BF09+.)
        # ===================================================================
        cat("  BF08: Plot #6 -- Stage 1 candidate pool ...\n")
    
        # == Tan triangle aesthetics (Part10e consistency) ====================
        TAN_FILL   <- "#8B6914"
        TAN_BORDER <- "white"
        TAN_ALPHA  <- 0.85
        TAN_STROKE <- 0.4
        TAN_SIZE   <- 3.0
        TAN_SHAPE  <- 25L     # inverted triangle
    
        # == Identify Original only cells (SKCS failing C2, outside band) =====
        df_all_orig_only <- df[df$is_orig_skcs &
                                   df$cell_class == "Original only", ]
        n_orig_only      <- nrow(df_all_orig_only)
        cat(sprintf("    Original only (C2-fail, outside band): n=%d\n",
                    n_orig_only))

        # BF09 legend: seven classes with counts that partition the dataset.
        # The four subgroups are the in-band core (Original + Refined).
        BF09_LEVS <- c(names(FOUR_COLS), "Original only",
                       "Stage 1 within delta_base", "Gap cells (non-SKCS)",
                       "Non-selected")

        .c9 <- rep("Non-selected", nrow(df))
        .c9[flat_cand]                              <- "Stage 1 within delta_base"
        .c9[in_gap & !df$is_orig_skcs]              <- "Gap cells (non-SKCS)"
        .is_or9 <- df$cell_class == "Original + Refined"
        .c9[.is_or9] <- as.character(df$plot_col_orig[.is_or9])
        .c9[df$cell_class == "Original only"]       <- "Original only"
        if (any(.c9[.is_or9] == "Non-selected"))
            stop(sum(.c9[.is_or9] == "Non-selected"),
                 " Original+Refined cell(s) carry no four-group label.",
                 call. = FALSE)
        df$bf09_class <- factor(.c9, levels = BF09_LEVS)

        bf09_n   <- table(df$bf09_class)
        bf09_lab <- sprintf("%s (n=%d)", BF09_LEVS, as.integer(bf09_n[BF09_LEVS]))
        names(bf09_lab) <- BF09_LEVS
        bf09_col <- c(FOUR_COLS, "Original only" = TAN_FILL,
                      "Stage 1 within delta_base" = BLUE_FILL,
                      "Gap cells (non-SKCS)" = GAP_FILL,
                      "Non-selected" = GRAY_COL)[BF09_LEVS]
        bf09_shp <- c(rep(16L, length(FOUR_COLS)), TAN_SHAPE, 21L, 21L, 16L)
        names(bf09_shp) <- BF09_LEVS

        cat("    BF09 legend classes:\n")
        for (.l in BF09_LEVS)
            cat(sprintf("      %-26s %d\n", .l, as.integer(bf09_n[.l])))
        cat(sprintf("      %-26s %d (must equal %d)\n", "TOTAL",
                    sum(as.integer(bf09_n)), nrow(df)))

        # Off-panel key layer, clipped away by coord_cartesian, so a class
        # with zero cells still shows as (n=0).
        df_bf09_key <- data.frame(
            u_all = rep(-1, length(BF09_LEVS)),
            v_all = rep(-1, length(BF09_LEVS)),
            bf09_class = factor(BF09_LEVS, levels = BF09_LEVS))
    
        p_bf08 <- ggplot() +
            # == Ribbons ======================================================
        geom_ribbon(data = df_rf, aes(x = u, ymin = f_lo, ymax = f_up),
                    fill = COL_FLAT, alpha = 0.18) +
            geom_ribbon(data = df_rf, aes(x = u, ymin = a_lo, ymax = a_up),
                        fill = COL_ADAPT, alpha = 0.18) +
            # == Gray non-gap cells (background) ==============================
        geom_point(data = df_std_gray,
                   aes(x = u_all, y = v_all, colour = plot_col_orig),
                   size = 0.55, alpha = 0.15) +
            # == 4-group coloured non-gap cells (foreground) ==================
        geom_point(data = df_std_col, 
                   aes(x = u_all, y = v_all, colour = plot_col_orig),
                   size = 1.25, alpha = 0.55) +
            # == spine ========================================================
        geom_line(data = df_rf, aes(x = u, y = g),
                  colour = COL_SPINE, linewidth = 1.1) +
            # == Flat-band Stage 1 candidates: blue circles ===================
        {
            if (n_flat_cand > 0)
                geom_point(data = df_flat_cand,
                           aes(x = u_all, y = v_all),
                           shape = 21, size = 3.0,
                           fill = BLUE_FILL, colour = "white",
                           stroke = 0.5, alpha = 0.88,
                           inherit.aes = FALSE)
            else
                geom_blank()
        } +
            # == Gap-zone cells: dark circles =================================
        {
            if (n_gap_cand > 0)
                geom_point(data = df_gap_cand,
                           aes(x = u_all, y = v_all),
                           shape = 21, size = 3.0,
                           fill = GAP_FILL, colour = "white",
                           stroke = 0.5, alpha = 0.88,
                           inherit.aes = FALSE)
            else
                geom_blank()
        } +
            # == Band boundary lines ==========================================
        geom_line(data = df_rf, aes(x = u, y = a_up),
                  colour = COL_ADAPT, linewidth = 0.55, linetype = "solid") +
            geom_line(data = df_rf, aes(x = u, y = a_lo),
                      colour = COL_ADAPT, linewidth = 0.55, linetype = "solid") +
            geom_line(data = df_rf, aes(x = u, y = f_up),
                      colour = COL_FLAT, linewidth = 0.45, linetype = "dashed") +
            geom_line(data = df_rf, aes(x = u, y = f_lo),
                      colour = COL_FLAT, linewidth = 0.45, linetype = "dashed") +
            # == Adaptive band extent verticals ===============================
        geom_vline(xintercept = c(u_adapt_lo, u_adapt_hi),
                   colour = "gray15", linewidth = 0.45, linetype = "dashed") +
            # == IP_u reference line ==========================================
        # IP_v line (the Pre/Post boundary)
        ip_hline_v(label = FALSE) +
        geom_vline(xintercept = IP_u, colour = COL_IP,
                   linewidth = 0.6, linetype = "dashed") +
            annotate("text", x = IP_u, y = -0.01,
                     label = "IP_u", size = 2.8, colour = COL_IP,
                     fontface = "bold") +
            # == Gap annotation ===============================================
        {
            if (n_gap_cand > 0)
                annotate("label",
                         x = IP_u + 0.04,
                         y = g_full[which.min(abs(u_full - IP_u))] + 0.15,
                         label = sprintf(
                             "Dark cells (n=%d): non-SKCS, inside\ndelta_adapt but outside delta_base",
                             n_gap_cand),
                         size = 2.4, fill = alpha("white", 0.85), colour = "gray20",
                         lineheight = 0.85)
            else
                geom_blank()
        } +
            # == Legend =======================================================
        legend_scale +
            coord_cartesian(xlim = c(0, 1), ylim = c(-0.04, 1.04)) +
            scale_y_continuous(expand = expansion(mult = 0.01)) +
            labs(
                title = sprintf(
                    "BF08: Plot #6 -- Stage 1 candidate pool: %s",
                    ds),
                subtitle = sprintf(
                    paste0("Blue circles (n=%d) = Stage 1 candidates within delta_base.  ",
                           "Dark circles (n=%d) = NON-SKCS gap cells (inside delta_adapt but outside delta_base).\n",
                           "Together blue + dark = all %d non-SKCS Stage 1 candidates considered for RSk-fCV recovery.  ",
                           "The gap zone also holds %d SKCS-core cell(s), drawn as ordinary four-group cells here and everywhere: ",
                           "the Stage 1 pool is non-SKCS by definition."),
                    n_flat_cand, n_gap_cand, n_flat_cand + n_gap_cand, n_gap_skcs),
                x = "u = Fx(GOI)",
                y = "v = Fc(cDFG)"
            ) +
            th_band(10) +
            theme(legend.position = "right",
                  legend.text = element_text(size = 8),
                  plot.subtitle = element_text(size = 7.5, lineheight = 1.2))
    
        # Caption: band provenance and spine monotonicity.
        p_bf08 <- p_bf08 +
            labs(caption = band_note) +
            theme(plot.caption = element_text(size = 6.2, colour = "gray35",
                                              hjust = 0, lineheight = 1.1))
    
        save_pp(p_bf08,
                bf_path("BF08_stage1_pool"),
                pdf_w = 11, pdf_h = 9, png_w = 2200, png_h = 1800)
        show_plot(p_bf08)
        cat("  BF08 saved\n")
    
    
        # ===================================================================
        # BF09: PLOT #7 -- GAP CELLS + ALL ORIGINAL ONLY (OUTSIDE BAND)
        # Original only = SKCS cells failing C2 (outside adaptive band).
        # Band proximity is the only Stage 1 criterion.
        # ===================================================================
        cat("  BF09: Plot #7 -- Gap cells + all Original only ...\n")
        cat(sprintf("    Original only (outside band): n=%d\n", n_orig_only))
    
        p_bf09 <- ggplot() +
            # == Ribbons ======================================================
        geom_ribbon(data = df_rf, aes(x = u, ymin = f_lo, ymax = f_up),
                    fill = COL_FLAT, alpha = 0.18) +
            geom_ribbon(data = df_rf, aes(x = u, ymin = a_lo, ymax = a_up),
                        fill = COL_ADAPT, alpha = 0.18) +
            # == Gray non-gap cells (background) ==============================
        geom_point(data = df_std_gray,
                   aes(x = u_all, y = v_all, colour = plot_col_orig),
                   # the key comes from df_bf09_key below
                   size = 0.55, alpha = 0.15, show.legend = FALSE) +
            # == 4-group coloured non-gap cells (foreground) ==================
        geom_point(data = df_std_col, 
                   aes(x = u_all, y = v_all, colour = plot_col_orig),
                   size = 1.25, alpha = 0.55, show.legend = FALSE) +
            # == spine ========================================================
        geom_line(data = df_rf, aes(x = u, y = g),
                  colour = COL_SPINE, linewidth = 1.1) +
            # == Flat-band Stage 1 candidates: blue circles ===================
        {
            if (n_flat_cand > 0)
                geom_point(data = df_flat_cand,
                           aes(x = u_all, y = v_all),
                           shape = 21, size = 3.0,
                           fill = BLUE_FILL, colour = "white",
                           stroke = 0.5, alpha = 0.88,
                           inherit.aes = FALSE)
            else
                geom_blank()
        } +
            # == Non-SKCS gap cells: dark circles =============================
        {
            if (n_gap_cand > 0)
                geom_point(data = df_gap_cand,
                           aes(x = u_all, y = v_all),
                           shape = 21, size = 3.0,
                           fill = GAP_FILL, colour = "white",
                           stroke = 0.5, alpha = 0.88,
                           inherit.aes = FALSE)
            else
                geom_blank()
        } +
            # == All Original only: tan inverted triangles (outside band) =====
        geom_point(data = df_all_orig_only,
                   aes(x = u_all, y = v_all),
                   shape = TAN_SHAPE, size = TAN_SIZE,
                   fill = TAN_FILL, colour = TAN_BORDER,
                   stroke = TAN_STROKE, alpha = TAN_ALPHA,
                   inherit.aes = FALSE) +
            # == Band boundary lines ==========================================
        geom_line(data = df_rf, aes(x = u, y = a_up),
                  colour = COL_ADAPT, linewidth = 0.55, linetype = "solid") +
            geom_line(data = df_rf, aes(x = u, y = a_lo),
                      colour = COL_ADAPT, linewidth = 0.55, linetype = "solid") +
            geom_line(data = df_rf, aes(x = u, y = f_up),
                      colour = COL_FLAT, linewidth = 0.45, linetype = "dashed") +
            geom_line(data = df_rf, aes(x = u, y = f_lo),
                      colour = COL_FLAT, linewidth = 0.45, linetype = "dashed") +
            # == Adaptive band extent verticals ===============================
        geom_vline(xintercept = c(u_adapt_lo, u_adapt_hi),
                   colour = "gray15", linewidth = 0.45, linetype = "dashed") +
            # == IP_u reference line ==========================================
        # IP_v line (the Pre/Post boundary)
        ip_hline_v(label = FALSE) +
        geom_vline(xintercept = IP_u, colour = COL_IP,
                   linewidth = 0.6, linetype = "dashed") +
            annotate("text", x = IP_u, y = -0.01,
                     label = "IP_u", size = 2.8, colour = COL_IP,
                     fontface = "bold") +
            # == Gap annotation ===============================================
        {
            if (n_gap_all > 0)
                annotate("label",
                         x = IP_u + 0.04,
                         y = g_full[which.min(abs(u_full - IP_u))] + 0.15,
                         label = sprintf(
                             paste0("Gap cells (n=%d): inside\ndelta_adapt but",
                                    " outside delta_base\ndark %d non-SKCS |",
                                    " %d SKCS-core, as ordinary cells"),
                             n_gap_all, n_gap_cand, n_gap_skcs),
                         size = 2.4, fill = alpha("white", 0.85), colour = "gray20",
                         lineheight = 0.85)
            else
                geom_blank()
        } +
            # == Original only annotation =====================================
        annotate("label",
                 x = 0.75, y = 0.08,
                 label = sprintf(
                     "Tan (v) = Original only (n=%d)\nSKCS outside band (C2 fail)",
                     n_orig_only),
                 size = 2.4, fill = alpha("white", 0.85), colour = TAN_FILL,
                 fontface = "bold", lineheight = 0.85) +
            # == Off-panel layer that keys the seven-class legend =============
        geom_point(data = df_bf09_key,
                   aes(x = u_all, y = v_all, colour = bf09_class,
                       shape = bf09_class, fill = bf09_class),
                   size = 2.2, show.legend = TRUE) +
            # == Legend: seven classes with counts ============================
        scale_colour_manual(values = bf09_col, labels = bf09_lab,
                            limits = BF09_LEVS, drop = FALSE,
                            name = "Cell group") +
            scale_shape_manual(values = bf09_shp, labels = bf09_lab,
                               limits = BF09_LEVS, drop = FALSE,
                               name = "Cell group") +
            # Shapes 21 and 25 are the filled forms, so fill is mapped too or
            # the blue, dark and tan keys draw hollow.
            scale_fill_manual(values = bf09_col, labels = bf09_lab,
                              limits = BF09_LEVS, drop = FALSE,
                              name = "Cell group") +
            guides(colour = guide_legend(order = 1,
                                         override.aes = list(size = 2.4, alpha = 1)),
                   shape  = guide_legend(order = 1,
                                         override.aes = list(size = 2.4, alpha = 1)),
                   fill   = guide_legend(order = 1,
                                         override.aes = list(size = 2.4, alpha = 1))) +
            coord_cartesian(xlim = c(0, 1), ylim = c(-0.04, 1.04)) +
            scale_y_continuous(expand = expansion(mult = 0.01)) +
            labs(
                title = sprintf(
                    "BF09: Plot #7 -- All Original only (outside band): %s",
                    ds),
                subtitle = sprintf(
                    paste0("Blue circles (n=%d) = Stage 1 candidates within delta_base.  ",
                           "Gap cells (n=%d): dark fill (n=%d) non-SKCS; the other %d are SKCS-core and keep the ordinary four-group dot.  ",
                           "Tan (v) (n=%d) = Original only: SKCS-core cells outside the adaptive band.\n",
                           "Band proximity is the only Stage 1 criterion; the former arm-concordance condition was dropped from CMRF.  ",
                           "All Original only cells have |v - g(u)| > delta_adapt(u)."),
                    n_flat_cand, n_gap_all, n_gap_cand, n_gap_skcs, n_orig_only),
                x = "u = Fx(GOI)",
                y = "v = Fc(cDFG)"
            ) +
            th_band(10) +
            theme(legend.position = "right",
                  legend.text = element_text(size = 8),
                  plot.subtitle = element_text(size = 7.5, lineheight = 1.2))
    
        # Caption: band provenance and spine monotonicity.
        p_bf09 <- p_bf09 +
            labs(caption = band_note) +
            theme(plot.caption = element_text(size = 6.2, colour = "gray35",
                                              hjust = 0, lineheight = 1.1))
    
        save_pp(p_bf09,
                bf_path("BF09_all_orig_only"),
                pdf_w = 11, pdf_h = 9, png_w = 2200, png_h = 1800)
        show_plot(p_bf09)
        cat("  BF09 saved\n")
    
    
        # ===================================================================
        # BF10: PLOT #8 -- RSk-fCV OUTCOME
        # Recovered (purple diamonds), Stage 2 removed (dark gray crosses),
        # Original only (tan triangles).  No blue/dark circles.
        # Shows the RSk-fCV gate outcome: who passed vs who was removed.
        # ===================================================================
        cat("  BF10: Plot #8 -- RSk-fCV outcome ...\n")
    
        # == Identify Recovered cells =========================================
        df_recovered <- df[df$cell_class == "Recovered", ]
        n_recovered  <- nrow(df_recovered)
        cat(sprintf("    Recovered: n=%d\n", n_recovered))
    
        # Purple diamond aesthetics (Part10e consistency)
        RECOV_FILL   <- "#7C3AED"
        RECOV_BORDER <- "white"
        RECOV_ALPHA  <- 0.88
        RECOV_STROKE <- 0.5
        RECOV_SIZE   <- 3.0
        RECOV_SHAPE  <- 23L     # filled diamond
    
        # == Identify Stage 2 removed cells ===================================
        # Non-SKCS cells that passed Stage 1 (band proximity to delta_adapt)
        # but were removed by the RSk-fCV iterative gate
        stage1_all    <- !df$is_orig_skcs & pass_stage1
        s2_removed    <- stage1_all & df$cell_class == "Excluded"
        # Part 10e's removed_by_stage2 is used when present; a disagreement with
        # the derived set is reported, not fatal.
        n_s2_disagree <- 0L
        if ("removed_by_stage2" %in% names(df)) {
            s2_col <- require_flag(df$removed_by_stage2, "removed_by_stage2")
            n_s2_disagree <- sum(s2_col != s2_removed)
            s2_removed <- s2_col
            s2_source  <- if (n_s2_disagree == 0L)
                "Part 10e removed_by_stage2 (agrees with the derived set)"
            else sprintf(paste0("Part 10e removed_by_stage2; it DISAGREES with",
                                " the set derived from pass_stage1 and",
                                " cell_class on %d cell(s) -- the column wins"),
                         n_s2_disagree)
        } else {
            s2_source  <- "derived from pass_stage1 and cell_class (Part 10e wrote no removed_by_stage2)"
        }
        df_s2_removed <- df[s2_removed, ]
        n_s2_removed  <- nrow(df_s2_removed)
        cat("    Stage 2 set from:", s2_source, "\n")
        cat(sprintf("    Stage 2 removed: n=%d\n", n_s2_removed))
        # Recovered and Stage-2-removed partition the non-SKCS Stage 1 pool.
        n_stage1_pool <- sum(stage1_all)
        s1_reconciles <- (n_recovered + n_s2_removed) == n_stage1_pool
        cat(sprintf(paste0("    Check: Recovered(%d) + Removed(%d) = %d",
                           " against the non-SKCS Stage 1 pool (%d) -- %s\n"),
                    n_recovered, n_s2_removed, n_recovered + n_s2_removed,
                    n_stage1_pool,
                    if (s1_reconciles) "agree" else "disagree"))
    
        # Dark gray cross aesthetics
        S2R_COL    <- "#404040"    # dark gray
        S2R_SIZE   <- 2.8
        S2R_STROKE <- 0.7
        S2R_ALPHA  <- 0.75
        S2R_SHAPE  <- 4L           # cross (x)
    
        # All cells as background, coloured by the SKCS core; the CMRF outcome
        # is carried by shape.
        df_gray_bf10 <- df_all_gray
        df_col_bf10  <- df_all_col
        n_orig_ref   <- nrow(df[df$cell_class == "Original + Refined", ])
    
        p_bf10 <- ggplot() +
            # == Ribbons ======================================================
        geom_ribbon(data = df_rf, aes(x = u, ymin = f_lo, ymax = f_up),
                    fill = COL_FLAT, alpha = 0.18) +
            geom_ribbon(data = df_rf, aes(x = u, ymin = a_lo, ymax = a_up),
                        fill = COL_ADAPT, alpha = 0.18) +
            # == Gray cells (background) ======================================
        geom_point(data = df_gray_bf10,
                   aes(x = u_all, y = v_all, colour = plot_col_orig),
                   size = 0.55, alpha = 0.15) +
            # == 4-group coloured cells (foreground) ==========================
        geom_point(data = df_col_bf10,
                   aes(x = u_all, y = v_all, colour = plot_col_orig),
                   size = 1.25, alpha = 0.55) +
            # == spine ========================================================
        geom_line(data = df_rf, aes(x = u, y = g),
                  colour = COL_SPINE, linewidth = 1.1) +
            # == Stage 2 removed: dark gray crosses ===========================
        {
            if (n_s2_removed > 0)
                geom_point(data = df_s2_removed,
                           aes(x = u_all, y = v_all),
                           shape = S2R_SHAPE, size = S2R_SIZE,
                           colour = S2R_COL,
                           stroke = S2R_STROKE, alpha = S2R_ALPHA,
                           inherit.aes = FALSE)
            else
                geom_blank()
        } +
            # == All Original only: tan inverted triangles ====================
        geom_point(data = df_all_orig_only,
                   aes(x = u_all, y = v_all),
                   shape = TAN_SHAPE, size = TAN_SIZE,
                   fill = TAN_FILL, colour = TAN_BORDER,
                   stroke = TAN_STROKE, alpha = TAN_ALPHA,
                   inherit.aes = FALSE) +
            # == Recovered: purple diamonds (top layer) =======================
        {
            if (n_recovered > 0)
                geom_point(data = df_recovered,
                           aes(x = u_all, y = v_all),
                           shape = RECOV_SHAPE, size = RECOV_SIZE,
                           fill = RECOV_FILL, colour = RECOV_BORDER,
                           stroke = RECOV_STROKE, alpha = RECOV_ALPHA,
                           inherit.aes = FALSE)
            else
                geom_blank()
        } +
            # == Band boundary lines ==========================================
        geom_line(data = df_rf, aes(x = u, y = a_up),
                  colour = COL_ADAPT, linewidth = 0.55, linetype = "solid") +
            geom_line(data = df_rf, aes(x = u, y = a_lo),
                      colour = COL_ADAPT, linewidth = 0.55, linetype = "solid") +
            geom_line(data = df_rf, aes(x = u, y = f_up),
                      colour = COL_FLAT, linewidth = 0.45, linetype = "dashed") +
            geom_line(data = df_rf, aes(x = u, y = f_lo),
                      colour = COL_FLAT, linewidth = 0.45, linetype = "dashed") +
            # == IP_u reference line ==========================================
        # IP_v line (the Pre/Post boundary)
        ip_hline_v(label = FALSE) +
        geom_vline(xintercept = IP_u, colour = COL_IP,
                   linewidth = 0.6, linetype = "dashed") +
            # v = 0.5 (cDFG marginal median) in gray
            geom_hline(yintercept = 0.5, colour = "gray55",
                       linewidth = 0.4, linetype = "dotted") +
            annotate("text", x = IP_u + 0.015, y = -0.01,
                     label = sprintf("IP_u=%.3f", IP_u),
                     size = 2.5, colour = COL_IP, fontface = "bold") +
            # == Recovered annotation =========================================
        annotate("label",
                 x = 0.18, y = 0.88,
                 label = sprintf("Purple <> = Recovered (n=%d)",
                                 n_recovered),
                 size = 2.4, fill = alpha("white", 0.85),
                 colour = RECOV_FILL,
                 fontface = "bold", lineheight = 0.85) +
            # == Stage 2 removed annotation ===================================
        annotate("label",
                 x = 0.18, y = 0.78,
                 label = sprintf("Gray x = Removed by RSk-fCV (n=%d)",
                                 n_s2_removed),
                 size = 2.4, fill = alpha("white", 0.85),
                 colour = S2R_COL,
                 fontface = "bold", lineheight = 0.85) +
            # == Original only annotation =====================================
        annotate("label",
                 x = 0.78, y = 0.08,
                 label = sprintf("Tan (v) = Original only (n=%d)",
                                 n_orig_only),
                 size = 2.4, fill = alpha("white", 0.85),
                 colour = TAN_FILL,
                 fontface = "bold", lineheight = 0.85) +
            # == Legend =======================================================
        legend_scale +
            coord_cartesian(xlim = c(0, 1), ylim = c(-0.04, 1.04)) +
            scale_y_continuous(expand = expansion(mult = 0.01)) +
            labs(
                title = sprintf(
                    "BF10: Plot #8 -- RSk-fCV outcome: %s", ds),
                subtitle = sprintf(
                    paste0("Purple <> (n=%d) = Recovered (passed Stage 1 + Stage 2).  ",
                           "Gray x (n=%d) = removed by RSk-fCV iterative gate.\n",
                           "Tan (v) (n=%d) = Original only (SKCS outside band).  ",
                           "CMRF-core = Orig+Ref (%d) + Recovered (%d) = %d cells."),
                    n_recovered, n_s2_removed, n_orig_only,
                    n_orig_ref, n_recovered, n_orig_ref + n_recovered),
                x = "u = Fx(GOI)",
                y = "v = Fc(cDFG)"
            ) +
            th_band(10) +
            theme(legend.position = "right",
                  legend.text = element_text(size = 8),
                  plot.subtitle = element_text(size = 7.5, lineheight = 1.2))
    
        # Caption: band provenance and spine monotonicity.
        p_bf10 <- p_bf10 +
            labs(caption = band_note) +
            theme(plot.caption = element_text(size = 6.2, colour = "gray35",
                                              hjust = 0, lineheight = 1.1))
    
        save_pp(p_bf10,
                bf_path("BF10_rskfcv_outcome"),
                pdf_w = 11, pdf_h = 9, png_w = 2200, png_h = 1800)
        show_plot(p_bf10)
        cat("  BF10 saved\n")
    
    
        # ===================================================================
        # BF11: PLOT #9 (FINAL) -- COMPLETE CMRF CLASSIFICATION
        # Clean final view: Orig+Ref as coloured background,
        # Recovered (purple diamonds).  Original only cells rendered
        # as regular gray background (not highlighted).
        # No gap circles, no tan triangles -- clean result only.
        # ===================================================================
        cat("  BF11: Plot #9 (final) -- Complete CMRF classification ...\n")
    
        n_orig_ref <- nrow(df[df$cell_class == "Original + Refined", ])
        n_cmrf_core <- n_orig_ref + n_recovered
    
        cat(sprintf("    Original+Refined: %d\n", n_orig_ref))
        cat(sprintf("    Recovered: %d\n", n_recovered))
        cat(sprintf("    CMRF-core: %d\n", n_cmrf_core))
        cat(sprintf("    Original only: %d\n", n_orig_only))
    
        # Original only cells are drawn as gray background in the final plot.
        df_bf11 <- df
        df_bf11$plot_col_final <- as.character(df_bf11$plot_col_orig)
        df_bf11$plot_col_final[df_bf11$cell_class == "Original only"] <- "Non-selected"
        df_bf11$plot_col_final <- factor(df_bf11$plot_col_final,
                                         levels = c("Non-selected", names(FOUR_COLS)))
    
        df_gray_all <- df_bf11[df_bf11$plot_col_final == "Non-selected", ]
        df_col_all  <- df_bf11[df_bf11$plot_col_final != "Non-selected", ]

        # BF11 legend: eight classes with counts that partition the dataset.
        # "Excluded (RSk-fCV)" reached Stage 2 and was rejected by the gate;
        # "Non-selected" never reached the gate.
        BF11_LEVS <- c(names(FOUR_COLS), "Recovered", "Original only",
                       "Excluded (RSk-fCV)", "Non-selected")

        .cls <- rep("Non-selected", nrow(df_bf11))
        # Order matters: each assignment overwrites the previous, so the most
        # specific class is written last.
        .cls[s2_removed]                            <- "Excluded (RSk-fCV)"
        .cls[df_bf11$cell_class == "Original only"] <- "Original only"
        .is_or <- df_bf11$cell_class == "Original + Refined"
        .cls[.is_or] <- as.character(df_bf11$plot_col_orig[.is_or])
        .cls[df_bf11$cell_class == "Recovered"]     <- "Recovered"
        # An Original+Refined cell is SKCS-core, so group_orig gives it one of
        # the four subgroups and never "Non-selected". Stop if that breaks.
        if (any(.cls[.is_or] == "Non-selected"))
            stop(sum(.cls[.is_or] == "Non-selected"),
                 " Original+Refined cell(s) carry no four-group label.",
                 call. = FALSE)
        df_bf11$bf11_class <- factor(.cls, levels = BF11_LEVS)

        bf11_n   <- table(df_bf11$bf11_class)
        bf11_lab <- sprintf("%s (n=%d)", BF11_LEVS, as.integer(bf11_n[BF11_LEVS]))
        names(bf11_lab) <- BF11_LEVS
        bf11_col <- c(FOUR_COLS, "Recovered" = RECOV_FILL,
                      "Original only" = TAN_FILL,
                      "Excluded (RSk-fCV)" = S2R_COL,
                      "Non-selected" = GRAY_COL)[BF11_LEVS]
        bf11_shp <- c(rep(16L, length(FOUR_COLS)), RECOV_SHAPE, TAN_SHAPE,
                      S2R_SHAPE, 16L)
        names(bf11_shp) <- BF11_LEVS

        cat("    BF11 legend classes:\n")
        for (.l in BF11_LEVS)
            cat(sprintf("      %-20s %d\n", .l, as.integer(bf11_n[.l])))
        cat(sprintf("      %-20s %d (must equal %d)\n", "TOTAL",
                    sum(as.integer(bf11_n)), nrow(df_bf11)))

        # One off-panel point per class keys the legend, so a class with zero
        # cells still shows as (n=0) instead of silently vanishing.
        df_bf11_key <- data.frame(
            u_all = rep(-1, length(BF11_LEVS)),
            v_all = rep(-1, length(BF11_LEVS)),
            bf11_class = factor(BF11_LEVS, levels = BF11_LEVS))
    
        p_bf11 <- ggplot() +
            # == Adaptive band ribbon only (steel gray, like Part10e 04b) =====
        geom_ribbon(data = df_rf, aes(x = u, ymin = a_lo, ymax = a_up),
                    fill = "gray70", alpha = 0.25) +
            # == Gray cells (Non-selected + Original only, background) ========
            # Legend keys come from df_bf11_key below.
        geom_point(data = df_gray_all,
                   aes(x = u_all, y = v_all, colour = plot_col_final),
                   size = 0.55, alpha = 0.15, show.legend = FALSE) +
            # == 4-group coloured cells (Orig+Ref only, foreground) ===========
        geom_point(data = df_col_all,
                   aes(x = u_all, y = v_all, colour = plot_col_final),
                   size = 1.25, alpha = 0.55, show.legend = FALSE) +
            # == spine ========================================================
        geom_line(data = df_rf, aes(x = u, y = g),
                  colour = COL_SPINE, linewidth = 1.1) +
            # == Recovered: purple diamonds (top layer) =======================
        {
            if (n_recovered > 0)
                geom_point(data = df_recovered,
                           aes(x = u_all, y = v_all),
                           shape = RECOV_SHAPE, size = RECOV_SIZE,
                           fill = RECOV_FILL, colour = RECOV_BORDER,
                           stroke = RECOV_STROKE, alpha = RECOV_ALPHA,
                           inherit.aes = FALSE)
            else
                geom_blank()
        } +
            # == Off-panel layer that keys the eight-class legend =============
            # Drawn at (-1, -1) and clipped away by coord_cartesian.
        geom_point(data = df_bf11_key,
                   aes(x = u_all, y = v_all, colour = bf11_class,
                       shape = bf11_class, fill = bf11_class),
                   size = 2.2, show.legend = TRUE) +
            # == Adaptive band boundary lines (darkened for final plot) =======
        geom_line(data = df_rf, aes(x = u, y = a_up),
                  colour = COL_ADAPT_BORDER, linewidth = 0.75, linetype = "solid") +
            geom_line(data = df_rf, aes(x = u, y = a_lo),
                      colour = COL_ADAPT_BORDER, linewidth = 0.75, linetype = "solid") +
            # == IP_u and v=0.5 reference lines ===============================
        # IP_v line (the Pre/Post boundary)
        ip_hline_v(label = TRUE) +
        geom_vline(xintercept = IP_u, colour = COL_IP,
                   linewidth = 0.6, linetype = "dashed") +
            # v = 0.5 (cDFG marginal median) in gray
            geom_hline(yintercept = 0.5, colour = "gray55",
                       linewidth = 0.4, linetype = "dotted") +
            # == IP label =====================================================
        annotate("text", x = IP_u + 0.015, y = -0.01,
                 label = sprintf("IP_u=%.3f", IP_u),
                 size = 2.5, colour = COL_IP, fontface = "bold") +
            annotate("text", x = 0.02, y = 0.53,
                     label = "v=0.5 [cDFG median]", size = 2.3,
                     colour = "gray45", fontface = "italic", hjust = 0) +
            # == Recovered annotation =========================================
        annotate("label",
                 x = 0.18, y = 0.88,
                 label = sprintf("Purple <> = Recovered (n=%d)",
                                 n_recovered),
                 size = 2.6, fill = alpha("white", 0.85),
                 colour = RECOV_FILL,
                 fontface = "bold", lineheight = 0.85) +
            # == CMRF-core count annotation ===================================
        annotate("label",
                 x = 0.02, y = 0.97,
                 label = sprintf(
                     "CMRF-core: %d = Orig+Ref (%d) + Recovered (%d)",
                     n_cmrf_core, n_orig_ref, n_recovered),
                 size = 2.4, fill = alpha("white", 0.85),
                 colour = "gray20", hjust = 0,
                 fontface = "bold", lineheight = 0.85) +
            # == Legend: eight classes with counts ============================
        scale_colour_manual(values = bf11_col, labels = bf11_lab,
                            limits = BF11_LEVS, drop = FALSE,
                            name = "Cell group") +
            scale_shape_manual(values = bf11_shp, labels = bf11_lab,
                               limits = BF11_LEVS, drop = FALSE,
                               name = "Cell group") +
            # Shapes 23 and 25 (the diamond and the triangle) are the filled
            # forms, so fill is mapped as well or their keys draw hollow.
            scale_fill_manual(values = bf11_col, labels = bf11_lab,
                              limits = BF11_LEVS, drop = FALSE,
                              name = "Cell group") +
            # All three scales carry the same name, labels and limits, so
            # ggplot merges them into one key showing each class's own symbol.
            guides(colour = guide_legend(order = 1,
                                         override.aes = list(size = 2.4, alpha = 1)),
                   shape  = guide_legend(order = 1,
                                         override.aes = list(size = 2.4, alpha = 1)),
                   fill   = guide_legend(order = 1,
                                         override.aes = list(size = 2.4, alpha = 1))) +
            coord_cartesian(xlim = c(0, 1), ylim = c(-0.04, 1.04)) +
            scale_y_continuous(expand = expansion(mult = 0.01)) +
            labs(
                title = sprintf(
                    "BF11: Plot #9 (final) -- Complete CMRF classification: %s",
                    ds),
                subtitle = sprintf(
                    paste0("CMRF-core = %d cells: Original+Refined (%d) + Recovered (%d).  ",
                           "Original only (%d) shown as background.\n",
                           "Spine g(u) fitted on TREP SKCS-core cells from both ",
                           "sides of the IP (n=%d).  ",
                           "Adaptive band: delta_base=%.5f, kappa=%.2f, bw=%.2f.\n",
                           "RSk-fCV gate: Wilcoxon p=%s (alpha=%s), ",
                           "rel. RMSE increase=%s vs background variability=%s."),
                    n_cmrf_core, n_orig_ref, n_recovered, n_orig_only,
                    n_spine_train, delta_base, KAPPA, BW,
                    fmt_p_10f(stage2_p),
                    if (is.finite(stage2_alpha)) sprintf("%.2f", stage2_alpha)
                    else "not recorded",
                    if (is.finite(stage2_delta)) sprintf("%.4f", stage2_delta)
                    else "NA",
                    if (is.finite(stage2_nf)) sprintf("%.4f", stage2_nf)
                    else "NA"),
                x = "u = Fx(GOI)  [GOI marginal CDF]",
                y = "v = Fc(cDFG)  [cDFG marginal CDF]"
            ) +
            th_band(10) +
            theme(legend.position = "right",
                  legend.text = element_text(size = 8),
                  plot.subtitle = element_text(size = 7.5, lineheight = 1.2))
    
        # Caption: band provenance and spine monotonicity.
        p_bf11 <- p_bf11 +
            labs(caption = band_note) +
            theme(plot.caption = element_text(size = 6.2, colour = "gray35",
                                              hjust = 0, lineheight = 1.1))
    
        save_pp(p_bf11,
                bf_path("BF11_final_classification"),
                pdf_w = 11, pdf_h = 9, png_w = 2200, png_h = 1800)
        show_plot(p_bf11)
        cat("  BF11 saved\n")

        # == Counts with their donor spread ===================================
        cat("\n  Counts and donor spread\n")
        ref_don <- df$donor[df$cell_class == "Original + Refined"]
        pool    <- df$donor
        sets <- list(
            "Original + Refined"      = df$donor[df$cell_class == "Original + Refined"],
            "Recovered"               = df$donor[df$cell_class == "Recovered"],
            "Original only"           = df$donor[df$cell_class == "Original only"],
            "Excluded"                = df$donor[df$cell_class == "Excluded"],
            "Gap zone (all)"          = df_gap$donor,
            "Gap zone (non-SKCS)"     = df_gap_cand$donor,
            "Stage 1 flat candidates" = df_flat_cand$donor,
            "Stage 1 pool (blue+dark)" = c(df_flat_cand$donor, df_gap_cand$donor))
        don_tab <- do.call(rbind, lapply(names(sets), function(nm)
            v10f_donor_row(nm, sets[[nm]], ref_don, pool)))
        don_tab <- cbind(dataset = ds, don_tab, stringsAsFactors = FALSE)
        # TVD_above_null flags sets above their own null; Original + Refined is
        # the reference set.
        don_tab$is_reference   <- don_tab$set == "Original + Refined"
        don_tab$TVD_above_null <- !don_tab$is_reference &
            is.finite(don_tab$TVD_vs_OrigRef) & is.finite(don_tab$TVD_null_p95) &
            don_tab$TVD_vs_OrigRef > don_tab$TVD_null_p95
        for (i in seq_len(nrow(don_tab)))
            cat(sprintf("    %-26s n=%6d | donors=%3d | top share=%s | TVD=%s (null p95 %s)%s\n",
                        don_tab$set[i], don_tab$n_cells[i], don_tab$n_donors[i],
                        if (is.na(don_tab$top_donor_share[i])) "NA"
                        else sprintf("%.3f", don_tab$top_donor_share[i]),
                        if (is.na(don_tab$TVD_vs_OrigRef[i])) "NA"
                        else sprintf("%.3f", don_tab$TVD_vs_OrigRef[i]),
                        if (is.na(don_tab$TVD_null_p95[i])) "NA"
                        else sprintf("%.3f", don_tab$TVD_null_p95[i]),
                        if (isTRUE(don_tab$is_reference[i])) "  [reference set]"
                        else if (isTRUE(don_tab$TVD_above_null[i])) "  [above null]"
                        else ""))

        class_tab <- data.frame(
            dataset = ds, arm = ARM_TAG,
            IP_u = round(IP_u, 6), delta_base = round(delta_base, 6),
            IP_y = round(IP_y_10f, 6), IP_v = round(IP_v_10f, 6),
            IP_v_source = if (is.finite(IP_v_10f)) IP_v_src else NA_character_,
            delta_max = round(delta_base * (1 + KAPPA), 6),
            kappa = KAPPA, bw = BW, coverage = CMRF_COVERAGE,
            band_max_dev = signif(band_dev, 3),
            n_spine_train = n_spine_train, n_spine_pre = n_spine_pre,
            n_spine_post = n_spine_post, spine_n_cells_10e = spine_n_10e,
            rskfcv_rel_delta = stage2_delta, rskfcv_noise_floor = stage2_nf,
            rskfcv_alpha = stage2_alpha,
            spine_monotone = mono_txt,
            group_label_source = group_src,
            group_label_source_drawn = "Part 10e group_orig (SKCS core), all 11 panels",
            group_col_post = .post_col,
            n_coloured_drawn = n_core_col,
            n_coloured_if_refined = n_ref_col,
            n_gap_skcs_core = n_gap_skcs,
            # BF11 legend counts
            bf11_pre_trep      = as.integer(bf11_n[["Pre-IP TREP"]]),
            bf11_pre_nontrep   = as.integer(bf11_n[["Pre-IP non-TREP"]]),
            bf11_post_trep     = as.integer(bf11_n[["Post-IP TREP"]]),
            bf11_post_nontrep  = as.integer(bf11_n[["Post-IP non-TREP"]]),
            bf11_recovered     = as.integer(bf11_n[["Recovered"]]),
            bf11_original_only = as.integer(bf11_n[["Original only"]]),
            bf11_excluded_gate = as.integer(bf11_n[["Excluded (RSk-fCV)"]]),
            bf11_non_selected  = as.integer(bf11_n[["Non-selected"]]),
            # BF09 legend counts
            bf09_pre_trep      = as.integer(bf09_n[["Pre-IP TREP"]]),
            bf09_pre_nontrep   = as.integer(bf09_n[["Pre-IP non-TREP"]]),
            bf09_post_trep     = as.integer(bf09_n[["Post-IP TREP"]]),
            bf09_post_nontrep  = as.integer(bf09_n[["Post-IP non-TREP"]]),
            bf09_original_only = as.integer(bf09_n[["Original only"]]),
            bf09_flat_cand     = as.integer(bf09_n[["Stage 1 within delta_base"]]),
            bf09_gap_nonskcs   = as.integer(bf09_n[["Gap cells (non-SKCS)"]]),
            bf09_non_selected  = as.integer(bf09_n[["Non-selected"]]),
            # Axis agreement of both label columns
            axis_agree_refined = round(axis_chk$agree, 6),
            axis_agree_orig = round(axis_chk$agree_orig, 6),
            axis_n_orig = axis_chk$n_orig,
            n_cells = nrow(df),
            n_donors = length(unique(df$donor[!is.na(df$donor)])),
            n_cells_no_donor = n_no_donor,
            n_orig_ref_class = sum(df$cell_class == "Original + Refined"),
            n_recovered_class = sum(df$cell_class == "Recovered"),
            n_orig_only_class = sum(df$cell_class == "Original only"),
            n_excluded_class = sum(df$cell_class == "Excluded"),
            n_gap_all = n_gap_all, n_gap_non_skcs = n_gap_cand,
            n_flat_candidates = n_flat_cand,
            n_stage2_removed = n_s2_removed, stage2_set_source = s2_source,
            n_stage2_disagree = n_s2_disagree,
            u_adapt_lo = round(u_adapt_lo, 6), u_adapt_hi = round(u_adapt_hi, 6),
            rskfcv_wilcox_p = stage2_p,
            stringsAsFactors = FALSE)

        invisible(list(dataset_tag = ds, IP_u = IP_u, delta_base = delta_base,
            n_spine_train = n_spine_train, n_spine_pre = n_spine_pre,
            n_spine_post = n_spine_post, spine_n_10e = spine_n_10e,
            stage2_delta = stage2_delta, stage2_nf = stage2_nf,
            stage2_alpha = stage2_alpha,
            IP_y = IP_y_10f, IP_v = IP_v_10f,
            IP_v_source = if (is.finite(IP_v_10f)) IP_v_src else NA_character_,
                       s2_source = s2_source, n_s2_removed = n_s2_removed,
                       n_s2_disagree = n_s2_disagree,
                       n_recovered = n_recovered,
                       n_stage1_pool = n_stage1_pool,
                       s1_reconciles = s1_reconciles,
                       band_dev = band_dev, mono_txt = mono_txt,
                       group_src = group_src, axis_chk = axis_chk,
                       group_src_orig = group_src_orig,
                       group_col_post = .post_col,
                       n_core_col = n_core_col, n_ref_col = n_ref_col,
                       n_col_diff = n_col_diff,
                       n_no_donor = n_no_donor,
                       n_gap_all = n_gap_all, n_gap_cand = n_gap_cand,
                       n_gap_skcs = n_gap_skcs,
                       bf11_n = bf11_n, bf11_levs = BF11_LEVS,
                       bf09_n = bf09_n, bf09_levs = BF09_LEVS,
                       n_flat_cand = n_flat_cand, stage2_p = stage2_p,
                       donor_tab = don_tab, class_tab = class_tab))
    }


    # =======================================================================
    # SECTION C: Execute
    # =======================================================================

    all_results <- vector("list", length(DATASETS))
    names(all_results) <- DATASETS
    for (ds in DATASETS) {
        # x[[ds]] <- NULL deletes the element; single-bracket with list()
        # stores it, so a failed dataset still reaches the ledger.
        all_results[ds] <- list(tryCatch(
            run_V10F_dataset(ds),
            error = function(e) {
                cat("  ERROR in dataset", ds, ":", conditionMessage(e), "\n")
                NULL
            }))
    }


    # =======================================================================
    # SECTION D: Workbook
    # =======================================================================
    cat("\n", strrep("=", 72), "\n")
    cat("  SAVING COUNTS WORKBOOK\n")
    cat(strrep("=", 72), "\n\n")

    ok_ds <- DATASETS[!vapply(all_results, is.null, logical(1L))]
    if (length(ok_ds) > 0) {
        tabs <- list(
            Band_Summary = do.call(rbind, lapply(ok_ds, function(d)
                all_results[[d]]$class_tab)),
            Donor_Spread = do.call(rbind, lapply(ok_ds, function(d)
                all_results[[d]]$donor_tab)))
        xlsx_out <- file.path(OUT_DIR,
                              paste0(FILE_PREFIX, "_Band_Formation_Counts.xlsx"))
        tryCatch({
            write_xlsx(tabs, path = xlsx_out)
            cat(sprintf("  Excel saved: %s (%d sheet(s))\n",
                        basename(xlsx_out), length(tabs)))
        }, error = function(e)
            cat("  Excel error:", conditionMessage(e), "\n"))
    } else {
        cat("  No dataset completed; nothing to write.\n")
    }


    # =======================================================================
    # SECTION E: Post-analysis guardrails
    # =======================================================================
    cat("\n  POST-ANALYSIS GUARDRAILS\n\n")

    for (ds in DATASETS) {
        r <- all_results[[ds]]
        if (is.null(r)) {
            gr10f(paste0("GR-10F-03-", ds), "FAIL", "Dataset returned NULL")
            next
        }
        ct <- r$class_tab
        # The gap-zone split must add back to the total.
        .gap_adds <- (r$n_gap_cand + r$n_gap_skcs) == r$n_gap_all
        gr10f(paste0("GR-10F-03-", ds), if (.gap_adds) "PASS" else "FAIL",
              sprintf(paste0("%d cell(s): OrigRef %d, Recovered %d, Orig-only",
                             " %d, Excluded %d | gap zone %d = non-SKCS %d",
                             " (dark fill) + SKCS-core %d%s"),
                      ct$n_cells, ct$n_orig_ref_class, ct$n_recovered_class,
                      ct$n_orig_only_class, ct$n_excluded_class,
                      r$n_gap_all, r$n_gap_cand, r$n_gap_skcs,
                      if (.gap_adds) ""
                      else " -- split does not add up"))

        # The ribbon drawn is the band Part 10e applied.
        gr10f(paste0("GR-10F-04-", ds),
              if (isTRUE(r$band_dev <= V10F_BAND_TOL)) "PASS" else "FAIL",
              sprintf("band recomputed vs stored delta_adapt: max |dev| = %.3e (tolerance %.1e, kappa %.2f, bw %.2f)",
                      r$band_dev, V10F_BAND_TOL, KAPPA, BW))

        # Spine monotonicity.
        gr10f(paste0("GR-10F-05-", ds),
              if (grepl("non-decreasing", r$mono_txt, fixed = TRUE)) "PASS"
              else "WARN",
              sprintf("%s",
                      r$mono_txt))

        # Axis the labels were made on (cDFG at IP_y).
        .ax <- r$axis_chk
        # BF11 legend classes must sum to n_cells.
        .b11    <- as.integer(r$bf11_n[r$bf11_levs])
        .b11_ok <- sum(.b11) == ct$n_cells
        gr10f(paste0("GR-10F-14-", ds), if (.b11_ok) "PASS" else "FAIL",
              sprintf("BF11 legend: %s | total %d against n_cells %d%s",
                      paste(sprintf("%s %d", r$bf11_levs, .b11), collapse = ", "),
                      sum(.b11), ct$n_cells,
                      if (.b11_ok) "" else " -- classes do not partition"))

        # BF09 legend classes must sum to n_cells.
        .b09    <- as.integer(r$bf09_n[r$bf09_levs])
        .b09_ok <- sum(.b09) == ct$n_cells
        gr10f(paste0("GR-10F-15-", ds), if (.b09_ok) "PASS" else "FAIL",
              sprintf("BF09 legend: %s | total %d against n_cells %d%s",
                      paste(sprintf("%s %d", r$bf09_levs, .b09), collapse = ", "),
                      sum(.b09), ct$n_cells,
                      if (.b09_ok) "" else " -- classes do not partition"))

        # Both label columns are judged.
        .ax_orig_ok <- isTRUE(.ax$have_orig) && isTRUE(.ax$agree_orig >= 0.999)
        .ax_plot_ok <- isTRUE(.ax$have) && isTRUE(.ax$agree >= 0.999)
        gr10f(paste0("GR-10F-10-", ds),
              if (.ax_plot_ok && .ax_orig_ok) "PASS"
              else if (!isTRUE(.ax$have)) "WARN"
              else if (!.ax_plot_ok || (isTRUE(.ax$have_orig) && !.ax_orig_ok))
                  "FAIL"
              else "WARN",
              if (!isTRUE(.ax$have))
                  sprintf("label axis not checked: %s",
                          if (length(.ax$why) == 1L && !is.na(.ax$why) &&
                              nzchar(.ax$why)) .ax$why
                          else "Part 10e supplied no IP_y, or df_refined has no cDFG_expr")
              else sprintf(paste0("rule `cDFG_expr < IP_y` (IP_y = %.6f, drawn",
                                  " at IP_v = %s): %s reproduces %.2f%% over",
                                  " %d cell(s) [loaded, not drawn];",
                                  " group_orig %s [all 11 panels]%s"),
                           .ax$IP_y,
                           if (length(r$IP_v) == 1L && is.finite(r$IP_v))
                               sprintf("%.6f", r$IP_v) else "not drawn",
                           r$group_col_post,
                           100 * .ax$agree, .ax$n,
                           if (isTRUE(.ax$have_orig))
                               sprintf("reproduces %.2f%% over %d cell(s)",
                                       100 * .ax$agree_orig, .ax$n_orig)
                           else "not checked: no labelled cell",
                           paste0(
                               if (!.ax_plot_ok)
                                   " -- refined labels not on the cDFG axis"
                               else "",
                               if (isTRUE(.ax$have_orig) && !.ax_orig_ok)
                                   " -- group_orig labels not on the cDFG axis"
                               else "")))

        # Which labelling was loaded.
        gr10f(paste0("GR-10F-06-", ds),
              if (grepl("Validation1", r$group_src, fixed = TRUE)) "PASS"
              else "WARN",
              sprintf("subgroup labels loaded: %s",
                      r$group_src))

        # Colours must be the SKCS core: coloured = OrigRef + Orig-only, and
        # the group_orig / group_refined difference = Recovered.
        .core_adds <- (ct$n_orig_ref_class + ct$n_orig_only_class) == r$n_core_col
        .diff_is_rec <- r$n_col_diff == ct$n_recovered_class
        gr10f(paste0("GR-10F-13-", ds),
              if (is.finite(r$n_core_col) && r$n_core_col > 0 &&
                  r$n_core_col <= r$n_ref_col && .core_adds && .diff_is_rec)
                  "PASS"
              else if (!.core_adds || !.diff_is_rec) "FAIL" else "WARN",
              sprintf(paste0("all 11 panels colour %d SKCS-core cell(s) from",
                             " group_orig = OrigRef %d + Orig-only %d%s;",
                             " %s would colour %d, a difference of %d against",
                             " Recovered %d%s"),
                      r$n_core_col, ct$n_orig_ref_class, ct$n_orig_only_class,
                      if (.core_adds) "" else " -- do not add up",
                      r$group_col_post, r$n_ref_col, r$n_col_diff,
                      ct$n_recovered_class,
                      if (.diff_is_rec) "" else " -- disagree"))

        # Donor spread of the Recovered set, against its own null.
        dt  <- r$donor_tab
        rec <- dt[dt$set == "Recovered", , drop = FALSE]
        # Verdict on Recovered only; gate-defined sets are expected to differ.
        .above <- if (!is.null(dt$TVD_above_null))
            dt$set[!is.na(dt$TVD_above_null) & dt$TVD_above_null]
        else character(0)
        .above_txt <- if (length(.above) > 0L)
            sprintf(" | also above own null p95: %s",
                    paste(.above, collapse = ", "))
        else " | other sets at or below own null p95"
        if (nrow(rec) == 1L) {
            ok <- isTRUE(rec$top_donor_share[1] < 0.5) &&
                isTRUE(rec$TVD_vs_OrigRef[1] <= rec$TVD_null_p95[1])
            gr10f(paste0("GR-10F-07-", ds), if (ok) "PASS" else "WARN",
                  sprintf(paste0("Recovered: %d cell(s) over %d donor(s), largest",
                                 " share %s, TVD vs OrigRef %s against a null p95",
                                 " of %s%s%s"),
                          rec$n_cells[1], rec$n_donors[1],
                          if (is.na(rec$top_donor_share[1])) "NA"
                          else sprintf("%.3f", rec$top_donor_share[1]),
                          if (is.na(rec$TVD_vs_OrigRef[1])) "NA"
                          else sprintf("%.3f", rec$TVD_vs_OrigRef[1]),
                          if (is.na(rec$TVD_null_p95[1])) "NA"
                          else sprintf("%.3f", rec$TVD_null_p95[1]),
                          if (r$n_no_donor > 0)
                              sprintf(" | %d cell(s) without a donor",
                                      r$n_no_donor) else "",
                          .above_txt))
        } else {
            gr10f(paste0("GR-10F-07-", ds), "WARN",
                  paste0("no Recovered row in the donor table", .above_txt))
        }

        # BF10 Stage 2 removals and their reconciliation with the Stage 1 pool.
        gr10f(paste0("GR-10F-09-", ds),
              if (!isTRUE(r$s1_reconciles)) "FAIL"
              else if (grepl("removed_by_stage2", r$s2_source, fixed = TRUE) &&
                       isTRUE(r$n_s2_disagree == 0L)) "PASS" else "WARN",
              sprintf(paste0("BF10 Stage 2 removals: %d cell(s) | Recovered %d",
                             " + Removed %d = %d against the non-SKCS Stage 1",
                             " pool %d (%s) | %s"),
                      r$n_s2_removed, r$n_recovered, r$n_s2_removed,
                      r$n_recovered + r$n_s2_removed, r$n_stage1_pool,
                      if (isTRUE(r$s1_reconciles)) "agree" else "disagree",
                      r$s2_source))

        # BF01 spine training set against Part 10e spine_n_cells.
        gr10f(paste0("GR-10F-11-", ds),
              if (!is.finite(r$spine_n_10e)) "WARN"
              else if (identical(as.integer(r$spine_n_10e),
                                 as.integer(r$n_spine_train))) "PASS" else "FAIL",
              if (!is.finite(r$spine_n_10e))
                  sprintf(paste0("Part 10e recorded no spine_n_cells; BF01 draws",
                                 " %d cell(s) (Pre-IP %d,",
                                 " Post-IP %d), not checked"),
                          r$n_spine_train, r$n_spine_pre, r$n_spine_post)
              else sprintf(paste0("BF01 spine training set: %d cell(s)",
                                  " (Pre-IP %d, Post-IP %d) against Part 10e",
                                  " spine_n_cells = %d"),
                           r$n_spine_train, r$n_spine_pre, r$n_spine_post,
                           as.integer(r$spine_n_10e)))

        # The dual RSk-fCV gate, both terms.
        gr10f(paste0("GR-10F-12-", ds),
              if (is.finite(r$stage2_delta) && is.finite(r$stage2_nf))
                  "PASS" else "WARN",
              if (is.finite(r$stage2_delta) && is.finite(r$stage2_nf))
                  sprintf(paste0("Gate closed with rel. RMSE increase %.4f",
                                 " against background variability %.4f (%s);",
                                 " Wilcoxon alpha = %s"),
                          r$stage2_delta, r$stage2_nf,
                          if (r$stage2_delta <= r$stage2_nf)
                              "effect-size criterion"
                          else "Wilcoxon criterion",
                          if (is.finite(r$stage2_alpha))
                              sprintf("%.2f", r$stage2_alpha) else "not recorded")
              else "no Stage 2 effect size recorded")

        # The RSk-fCV p quoted in BF11.
        gr10f(paste0("GR-10F-08-", ds),
              if (is.finite(r$stage2_p)) "PASS" else "WARN",
              if (is.finite(r$stage2_p))
                  sprintf("RSk-fCV Wilcoxon p = %s (quoted in BF11)",
                          if (r$stage2_p > 0 && r$stage2_p < 1e-4)
                              format(r$stage2_p, digits = 3, scientific = TRUE)
                          else sprintf("%.4f", r$stage2_p))
              else "no Stage 2 result recorded")
    }

    gr10f_write()
    cat(sprintf("\n  Final guardrail tally: %d PASS / %d WARN / %d FAIL\n",
                .gr10f_pass, .gr10f_warn, .gr10f_fail))
    if (.gr10f_fail > 0) {
        cat("  Failures present; see the ledger.\n")
    } else if (.gr10f_warn > 0) {
        cat(sprintf("  %d passed, %d warning(s), no failures.\n",
                    .gr10f_pass, .gr10f_warn))
    } else {
        cat(sprintf("  All %d checks passed.\n", .gr10f_pass))
    }

    cat("\n", strrep("=", 72), "\n")
    cat("  Part10f Band Formation COMPLETE [", ARM_TAG, "]\n")
    cat("  Output directory:", OUT_DIR, "\n")
    cat("  Files generated:\n")
    .transient <- c("_RUN_INCOMPLETE.txt", basename(MSG_PATH_10F))
    for (f in setdiff(sort(list.files(OUT_DIR)), .transient)) {
        sz <- round(file.size(file.path(OUT_DIR, f)) / 1e3, 0)
        cat(sprintf("    %-65s  %5d KB\n", f, sz))
    }
    cat(strrep("=", 72), "\n")

    cat("\n  Log written to:", LOG_PATH_10F, "\n")
    cat(sprintf("\nPart10f Band Formation - log ended: %s\n", Sys.time()))

    .p10f_completed <- TRUE
    if (sink.number(type = "message") != 2L) sink(type = "message")
    while (sink.number() > 0) sink()
    suppressWarnings(file.remove(INCOMPLETE_MARKER_10F))

    invisible(list(arm = arm_tag, dir = OUT_DIR, pass = .gr10f_pass,
                   warn = .gr10f_warn, fail = .gr10f_fail))
}


# =============================================================================
# SECTION 0c: Driver
# =============================================================================

cep_run_arms(PART_DIR_TAG, run_one_arm)
