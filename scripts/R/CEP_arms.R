# =============================================================================
# CEP_arms.R
# SEA-AD CEP-IP/CMRF, shared arm table and driver (Parts 9 onwards)
#   Defines the four arms (Mic13_PVM3, PVM3_all, PVM1_all, HKG8_matched), their
#   output folders <part>/<part>_<subdir>/ and file prefixes
#   <part>_<norm>_<file_tag>_, and runs a part over every arm, or over the one
#   arm named in CEP_ARM.
#
# Usage: sourced by the analysis scripts, not run directly.
#   source(file.path(CEP_HOME, "CEP_arms.R")); cep_run_arms(PART_TAG, run_one_arm)
#
# Outputs: none of its own; creates each part's arm folders.
# =============================================================================

# == Installation paths =======================================================
# CEP_HOME is the SEA-AD folder; set it before sourcing to override.
if (!exists("CEP_HOME")) CEP_HOME <- Sys.getenv("SEAAD_ROOT", getwd())
if (!dir.exists(CEP_HOME)) {
    stop("CEP_HOME does not exist: ", CEP_HOME,
         "\n  Set CEP_HOME to the SEA-AD folder before sourcing.",
         call. = FALSE)
}
if (!exists("RESULTS_ROOT")) RESULTS_ROOT <- CEP_HOME


# == The arms =================================================================
#   arm      ARM_TAG reported by the configuration; also the Part 8 subtree
#   config   configuration file that defines the arm
#   subdir   folder name within the part folder, prefixed with the part tag
#   file_tag inserted into every output file name
CEP_ARM_TABLE <- data.frame(
    arm      = c("Mic13_PVM3",  "PVM3_all",  "PVM1_all",  "HKG8_matched"),
    config   = c("R_config.R",  "R_config.R", "R_config.R",
                 "R_config_HKG8.R"),
    subdir   = c("Mic13_PVM3",  "PVM_3",     "PVM_1",     "HKG8"),
    file_tag = c("M4",          "PVM_3",     "PVM_1",     "HKG8"),
    label    = c("Mic13 PVM_3 (M4-top, primary)",
                 "PVM_3 all cells (mask specificity)",
                 "PVM_1 all cells (supertype specificity)",
                 "HKG8 control (panel specificity)"),
    stringsAsFactors = FALSE)


# == Path and name helpers ====================================================
cep_part_root <- function(part_tag) file.path(RESULTS_ROOT, part_tag)

cep_arm_dir <- function(part_tag, subdir)
    file.path(cep_part_root(part_tag), paste0(part_tag, "_", subdir))

cep_out_prefix <- function(part_tag, norm_tag, file_tag)
    paste0(part_tag, "_", norm_tag,
           if (nzchar(file_tag)) paste0("_", file_tag) else "")

# Folder of an arm within a part, looked up by ARM_TAG.
cep_arm_dir_for <- function(part_tag, arm_tag) {
    i <- match(arm_tag, CEP_ARM_TABLE$arm)
    if (is.na(i))
        stop("Unknown arm '", arm_tag, "'. Known arms: ",
             paste(CEP_ARM_TABLE$arm, collapse = ", "), call. = FALSE)
    cep_arm_dir(part_tag, CEP_ARM_TABLE$subdir[i])
}

cep_file_tag_for <- function(arm_tag) {
    i <- match(arm_tag, CEP_ARM_TABLE$arm)
    if (is.na(i))
        stop("Unknown arm '", arm_tag, "'.", call. = FALSE)
    CEP_ARM_TABLE$file_tag[i]
}


# == Arm selection ============================================================
cep_select_arms <- function() {
    if (exists("CEP_ARM", envir = globalenv()) &&
        !is.null(get("CEP_ARM", envir = globalenv()))) {
        a <- get("CEP_ARM", envir = globalenv())
        if (!a %in% CEP_ARM_TABLE$arm) {
            stop("CEP_ARM is '", a, "', which is not in CEP_ARM_TABLE. ",
                 "Known arms: ", paste(CEP_ARM_TABLE$arm, collapse = ", "),
                 call. = FALSE)
        }
        CEP_ARM_TABLE[CEP_ARM_TABLE$arm == a, , drop = FALSE]
    } else {
        CEP_ARM_TABLE
    }
}


# == Driver ===================================================================
# Runs arm_fun once per selected arm; a failed arm is reported and the rest run.
cep_run_arms <- function(part_tag, arm_fun) {

    selected  <- cep_select_arms()
    part_root <- cep_part_root(part_tag)
    if (!dir.exists(part_root)) dir.create(part_root, recursive = TRUE)

    # Files left directly in the part folder, outside the arm folders.
    loose <- setdiff(list.files(part_root),
                     list.dirs(part_root, recursive = FALSE, full.names = FALSE))
    if (length(loose) > 0) {
        cat("\n  NOTE: ", length(loose), " file(s) outside the arm folders in ",
            part_root, ":\n", sep = "")
        for (f in head(loose, 8)) cat("          ", f, "\n", sep = "")
        if (length(loose) > 8)
            cat("          ... and ", length(loose) - 8, " more\n", sep = "")
        cat("\n")
    }

    cat("\n", strrep("#", 76), "\n", sep = "")
    cat("  ", part_tag, "\n", sep = "")
    cat("  Home        : ", CEP_HOME, "\n", sep = "")
    cat("  Results     : ", part_root, "\n", sep = "")
    cat("  Arms to run : ", nrow(selected), " of ", nrow(CEP_ARM_TABLE),
        "\n", sep = "")
    for (i in seq_len(nrow(selected)))
        cat("      ", format(selected$arm[i], width = 14), " -> ",
            format(basename(cep_arm_dir(part_tag, selected$subdir[i])),
                   width = 20),
            selected$label[i], "\n", sep = "")
    cat(strrep("#", 76), "\n", sep = "")

    t0      <- Sys.time()
    outcome <- vector("list", nrow(selected))

    for (i in seq_len(nrow(selected))) {
        row     <- selected[i, ]
        out_dir <- cep_arm_dir(part_tag, row$subdir)
        if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

        outcome[[i]] <- tryCatch(
            arm_fun(row$arm, row$config, out_dir, row$file_tag),
            error = function(e) {
                # Restore the console if the arm failed while its log was open.
                while (sink.number() > 0) sink()
                cat("\n  ARM FAILED: ", row$arm, "\n", sep = "")
                cat("  ", conditionMessage(e), "\n\n", sep = "")
                list(pass = NA_integer_, warn = NA_integer_,
                     fail = NA_integer_, error = conditionMessage(e))
            })
    }

    cat("\n", strrep("#", 76), "\n", sep = "")
    cat("  ", part_tag, " COMPLETE -- ",
        format(round(difftime(Sys.time(), t0, units = "mins"), 1)), "\n",
        sep = "")
    cat(strrep("#", 76), "\n\n", sep = "")
    cat(sprintf("  %-14s  %-22s  %6s  %6s  %6s  %s\n",
                "Arm", "Folder", "PASS", "WARN", "FAIL", "Status"))
    cat("  ", strrep("-", 82), "\n", sep = "")
    for (i in seq_along(outcome)) {
        o <- outcome[[i]]
        cat(sprintf("  %-14s  %-22s  %6s  %6s  %6s  %s\n",
                    selected$arm[i],
                    basename(cep_arm_dir(part_tag, selected$subdir[i])),
                    ifelse(is.na(o$pass), "-", o$pass),
                    ifelse(is.na(o$warn), "-", o$warn),
                    ifelse(is.na(o$fail), "-", o$fail),
                    if (!is.null(o$error)) "FAILED TO RUN"
                    else if (isTRUE(o$fail > 0)) "check failures"
                    else "ok"))
    }
    cat("\n  Results under: ", part_root, "\n", sep = "")

    # Leave the session holding the primary arm's configuration.
    suppressWarnings(rm(list = intersect("CEP_ARM",
                                         ls(envir = globalenv(),
                                            all.names = TRUE)),
                        envir = globalenv()))
    source(file.path(CEP_HOME, "R_config.R"))
    cat("  Session configuration reset to the primary arm: ",
        get("ARM_TAG", envir = globalenv()), "\n\n", sep = "")

    invisible(outcome)
}

cat("CEP_arms.R loaded.\n")
cat("  CEP_HOME    :", CEP_HOME, "\n")
cat("  RESULTS_ROOT:", RESULTS_ROOT, "\n")
cat("  Arms        :", paste(CEP_ARM_TABLE$arm, collapse = ", "), "\n")
