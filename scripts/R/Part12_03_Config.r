# =============================================================================
# Part12_03_Config.r
# SEA-AD CEP-IP/CMRF, Part 12 rerun, shared settings and functions
#   Sourced first by every Part 12 R step; not run on its own. Reads
#   Part12_01_config.json and provides paths and inputs, logging, the donor
#   bootstrap, BH, the number format, the stats-file writer and the
#   pre-specified readings.
#
# Usage (from the SEA-AD folder, or set SEAAD_ROOT)
#   sourced by each step: source("Part12_rerun/Part12_03_Config.r")
#
# Outputs: none of its own.
# =============================================================================

for (.p in c("jsonlite", "readxl", "Matrix"))
    if (!requireNamespace(.p, quietly = TRUE))
        stop(sprintf("package '%s' is missing: install.packages(\"%s\")", .p, .p),
             call. = FALSE)
rm(.p)

# == Root, config, session settings ===========================================
P12_ROOT <- Sys.getenv("SEAAD_ROOT", getwd())
P12_CFG_PATH <- file.path(P12_ROOT, "Part12_rerun", "Part12_01_config.json")
if (!file.exists(P12_CFG_PATH))
    stop("config not found: ", P12_CFG_PATH, call. = FALSE)
P12_CFG <- jsonlite::fromJSON(P12_CFG_PATH, simplifyVector = TRUE,
                              simplifyDataFrame = FALSE, simplifyMatrix = FALSE)
P12_NEED_CONFIG <- 4L
if (!isTRUE(as.integer(P12_CFG$config_version) >= P12_NEED_CONFIG))
    stop(sprintf("config is version %s; version %d or later is needed",
                 P12_CFG$config_version, P12_NEED_CONFIG), call. = FALSE)
P12_RERUN <- file.path(P12_ROOT, P12_CFG$rerun_dir)
P12_ARMS <- P12_CFG$arm_order
P12_DATASETS <- names(P12_CFG$datasets)

invisible(Sys.setlocale("LC_COLLATE", P12_CFG$inference$collate_locale))
RNGkind(P12_CFG$inference$rng_kind)
options(warn = 1, stringsAsFactors = FALSE)

p12_path <- function(...) file.path(P12_ROOT, ...)
p12_rel <- function(path) {
    p <- normalizePath(path, winslash = "/", mustWork = FALSE)
    pre <- paste0(P12_ROOT, "/")
    ifelse(startsWith(p, pre), substring(p, nchar(pre) + 1L), p)
}

# == Folders for one step =====================================================
p12_step_dirs <- function(step, name, parent = NULL) {
    # One folder per step; a figure's sections sit in the figure folder
    # `parent` (e.g. 07_Fig15/15AB_ContamDoublets).
    out <- if (is.null(parent)) file.path(P12_RERUN, sprintf("%s_%s", step, name)) else
        file.path(P12_RERUN, parent, name)
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    list(out = out, tables = out, figures = out, logs = out)
}

# Removes stale outputs before a step runs; a file that cannot be removed
# stops the step.
p12_clean <- function(dirs) {
    for (d in unique(dirs)) {
        f <- list.files(d, full.names = TRUE)
        f <- f[!dir.exists(f)]
        if (!length(f)) next
        suppressWarnings(file.remove(f))
        left <- f[file.exists(f)]
        if (length(left))
            stop(sprintf("cannot remove old output(s), probably open in another program: %s",
                         paste(basename(left), collapse = ", ")), call. = FALSE)
    }
}

# == Logging ==================================================================
P12 <- new.env()

p12_say <- function(...) {
    s <- paste0(...)
    cat(s, "\n", sep = "")
    if (!is.null(P12$con) && isOpen(P12$con)) {
        writeLines(iconv(s, to = "ASCII", sub = "?"), P12$con)
        flush(P12$con)
    }
    invisible(s)
}

p12_chk <- function(ok, desc, detail = "", warn_only = FALSE) {
    P12$id <- P12$id + 1L
    ok <- isTRUE(ok)
    tag <- if (ok) "PASS" else if (warn_only) "WARN" else "FAIL"
    P12$n[[tag]] <- P12$n[[tag]] + 1L
    p12_say(sprintf("[%s] %s-%02d %s | %s", tag, P12$step, P12$id, desc, detail))
    invisible(ok)
}

p12_note <- function(desc, detail = "") p12_say(sprintf("[NOTE] %s | %s", desc, detail))

p12_abort <- function(msg) stop(msg, call. = FALSE)

p12_output <- function(path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    P12$outputs <- c(P12$outputs, path)
    invisible(path)
}

p12_log_pkgs <- function(pkgs) {
    v <- vapply(pkgs, function(p) as.character(utils::packageVersion(p)), "")
    p12_say("Packages: ", paste(sprintf("%s %s", pkgs, v), collapse = " | "))
}

# Opens the log and writes the header. `script` is the file name in Part12_rerun.
p12_start <- function(step, script, dirs, pkgs = character(0)) {
    # Close a log left open by an interrupted run.
    if (!is.null(P12$con)) try(close(P12$con), silent = TRUE)
    P12$con <- NULL
    P12$tree <- NULL
    P12$step <- step
    P12$script <- script
    P12$dirs <- dirs
    P12$id <- 0L
    P12$n <- c(PASS = 0L, FAIL = 0L, WARN = 0L)
    P12$outputs <- character(0)
    P12$t0 <- Sys.time()
    P12$done <- FALSE
    P12$run_id <- paste0(format(P12$t0, "%Y%m%d-%H%M"), "-",
                         if (.Platform$OS.type == "windows") "windows" else "fedora")
    # Log named after the step folder: Part12_<folder>_log.txt.
    P12$log <- file.path(dirs$logs, sprintf("Part12_%s_log.txt", basename(dirs$out)))
    P12$con <- file(P12$log, open = "wt")
    sp <- file.path(P12_RERUN, script)
    p12_say("Part 12 rerun -- step ", step)
    p12_say("Script: ", script,
            if (file.exists(sp)) "" else " | NOT FOUND in Part12_rerun")
    p12_say("Shared file: Part12_03_Config.r")
    p12_say("Run ID: ", P12$run_id)
    p12_say("Start: ", format(P12$t0, "%Y-%m-%d %H:%M:%S"))
    p12_say("Machine: ", Sys.info()[["nodename"]], " | ", R.version.string)
    p12_log_pkgs(unique(c("jsonlite", "readxl", "Matrix", pkgs)))
    p12_say(sprintf("Config: %s | version %s", p12_rel(P12_CFG_PATH),
                    P12_CFG$config_version))
    p12_say("Labels version: ", P12_CFG$labels$labels_version,
            " | cut token: ", P12_CFG$labels$cut_token)
    p12_chk(file.exists(sp), "script saved in Part12_rerun", script)
    p12_say("")
}

# Writes the outputs list, the verdict and the ledger line. On FAIL the
# outputs are renamed *_REJECTED.
p12_finish <- function() {
    if (isTRUE(P12$done)) return(invisible(NULL))
    P12$done <- TRUE
    verdict <- if (P12$n[["FAIL"]] == 0L) "PASS" else "FAIL"
    dur <- as.numeric(difftime(Sys.time(), P12$t0, units = "mins"))
    p12_say("")
    outs <- unique(P12$outputs)
    if (verdict == "PASS") {
        p12_say("Outputs:")
        for (f in outs)
            p12_say(sprintf("  %s | %.1f MB", p12_rel(f),
                            file.size(f) / 2^20))
    } else {
        # A file that cannot be renamed is deleted, or else named in the log.
        for (f in outs[file.exists(outs)]) {
            to <- sub("(\\.[A-Za-z]+)(\\.gz)?$", "_REJECTED\\1\\2", f)
            if (!suppressWarnings(file.rename(f, to))) suppressWarnings(file.remove(f))
            if (file.exists(f))
                p12_say("  COULD NOT REJECT ", p12_rel(f), ": delete it by hand")
        }
    }
    p12_say(sprintf("Checks: %d PASS, %d FAIL, %d WARN | duration %.1f min",
                    P12$n[["PASS"]], P12$n[["FAIL"]], P12$n[["WARN"]], dur))
    p12_say("VERDICT: ", verdict)
    close(P12$con)
    cat(sprintf("%s\t%s\t%s\t%s\t%s\tPASS=%d FAIL=%d WARN=%d\t%.1f min\n",
                format(Sys.time(), "%Y-%m-%d %H:%M"), P12$step, P12$script,
                P12$run_id, verdict, P12$n[["PASS"]], P12$n[["FAIL"]],
                P12$n[["WARN"]], dur),
        file = p12_path(P12_CFG$ledger_file), append = TRUE)
    if (verdict == "FAIL")
        message("Step ", P12$step, " FAILED: see ", P12$log)
    invisible(verdict)
}

# Runs a step's main function; any error is logged as a FAIL.
p12_run <- function(main) {
    tryCatch(main(),
             error = function(e) {
                 p12_say("")
                 p12_chk(FALSE, "run stopped", conditionMessage(e))
             },
             interrupt = function(e) {
                 p12_say("")
                 p12_chk(FALSE, "run stopped", "interrupted by the user")
             })
    p12_finish()
}

# == Inputs ===================================================================
.p12_fill <- function(x, arm = NULL, dataset = NULL) {
    if (!is.null(arm)) {
        x <- gsub("{arm}", arm, x, fixed = TRUE)
        x <- gsub("{export_dir}", P12_CFG$arms[[arm]]$export_dir, x, fixed = TRUE)
    }
    if (!is.null(dataset)) {
        x <- gsub("{dataset}", dataset, x, fixed = TRUE)
        x <- gsub("{dataset_dir}", P12_CFG$datasets[[dataset]], x, fixed = TRUE)
    }
    x
}

# Every file under the root outside the excluded folders, listed once per session.
.p12_tree <- function() {
    if (is.null(P12$tree)) {
        top <- list.dirs(P12_ROOT, recursive = FALSE, full.names = FALSE)
        keep <- setdiff(top, P12_CFG$inputs$exclude_from_search)
        f <- c(list.files(P12_ROOT, full.names = TRUE),
               unlist(lapply(file.path(P12_ROOT, keep), list.files,
                             recursive = TRUE, full.names = TRUE)))
        P12$tree <- f[!dir.exists(f)]
    }
    P12$tree
}

# Resolves one input from the config. A pattern must match exactly one file;
# archive_fallback searches the archive only when no other copy exists.
p12_input <- function(key, arm = NULL, dataset = NULL) {
    spec <- P12_CFG$inputs[[key]]
    if (is.null(spec)) p12_abort(sprintf("no input '%s' in the config", key))
    label <- paste(c(key, arm, dataset), collapse = " ")
    if (!is.null(spec$path)) {
        p <- p12_path(.p12_fill(spec$path, arm, dataset))
        if (!file.exists(p)) p12_abort(sprintf("input %s not found: %s", label, p))
        return(p)
    }
    if (is.null(spec$pattern)) {
        p <- p12_path(spec$folder)
        if (!dir.exists(p)) p12_abort(sprintf("input %s folder not found: %s", label, p))
        return(p)
    }
    pat <- .p12_fill(spec$pattern, arm, dataset)
    skip <- P12_CFG$inputs$skip_file_prefix
    pick <- function(files) {
        b <- basename(files)
        files[grepl(pat, b, ignore.case = TRUE) & !startsWith(b, skip)]
    }
    if (isTRUE(spec$recursive) && spec$folder == ".") {
        hits <- pick(.p12_tree())
    } else {
        base <- p12_path(spec$folder)
        if (!dir.exists(base))
            p12_abort(sprintf("input %s folder not found: %s", label, base))
        hits <- pick(list.files(base, recursive = isTRUE(spec$recursive),
                                full.names = TRUE))
    }
    if (!length(hits) && isTRUE(spec$archive_fallback)) {
        arch <- p12_path(P12_CFG$archive_dir)
        hits <- if (dir.exists(arch))
            pick(list.files(arch, recursive = TRUE, full.names = TRUE)) else character(0)
        if (length(hits) == 1L)
            p12_note(sprintf("input %s taken from the archive", label), p12_rel(hits))
    }
    hits <- unique(normalizePath(hits, winslash = "/"))
    if (length(hits) != 1L)
        p12_abort(sprintf("input %s: %d files match '%s'%s", label, length(hits), pat,
                          if (length(hits)) paste0(": ", paste(p12_rel(hits),
                                                               collapse = ", ")) else ""))
    hits
}

# Records an input in the log.
p12_log_input <- function(path, label) {
    p12_say(sprintf("  input %-28s %s | %.1f MB", label, p12_rel(path),
                    file.size(path) / 2^20))
    invisible(path)
}

# Exact column name, ignoring case.
p12_col <- function(df, name, what) {
    hit <- names(df)[tolower(names(df)) == tolower(name)]
    if (length(hit) != 1L)
        p12_abort(sprintf("%s: column '%s' %s. Columns present: %s", what, name,
                          if (length(hit)) "is duplicated" else "not found",
                          paste(names(df), collapse = ", ")))
    hit
}

# Excel sheet reader; refuses a workbook that Excel has open.
p12_read_sheet <- function(path, sheet) {
    lock <- file.path(dirname(path), paste0("~$", basename(path)))
    if (file.exists(lock))
        p12_abort(sprintf("%s is open in Excel; close it and rerun", basename(path)))
    sh <- readxl::excel_sheets(path)
    if (!(sheet %in% sh))
        p12_abort(sprintf("%s has no sheet '%s'. Sheets: %s", basename(path), sheet,
                          paste(sh, collapse = ", ")))
    as.data.frame(readxl::read_excel(path, sheet = sheet, guess_max = 1e6))
}

# Reads only the named columns of a CSV, all as text (names matched ignoring case).
p12_read_csv_cols <- function(path, cols, what) {
    hdr <- utils::read.csv(path, nrows = 1L, check.names = FALSE, colClasses = "character")
    hit <- vapply(cols, function(cn) p12_col(hdr, cn, what), "")
    cc <- rep("NULL", ncol(hdr))
    cc[match(hit, names(hdr))] <- "character"
    out <- utils::read.csv(path, check.names = FALSE, colClasses = cc)
    out <- out[, hit, drop = FALSE]
    names(out) <- names(cols)
    out
}

p12_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

p12_donor_key <- function(x) gsub("[^A-Z0-9]", "", toupper(as.character(x)))

# == Statistics ===============================================================
p12_seed <- function(step) {
    code <- P12_CFG$inference$step_codes[[step]]
    if (is.null(code)) p12_abort(sprintf("no step code for '%s' in the config", step))
    as.integer(P12_CFG$inference$master_seed) + as.integer(code)
}

# Donor resample indices from the step seed: a B x n matrix of positions into `donors`.
p12_boot_index <- function(donors, step, B = P12_CFG$inference$bootstrap_B) {
    set.seed(p12_seed(step))
    n <- length(donors)
    matrix(sample.int(n, n * B, replace = TRUE), nrow = B, ncol = n)
}

p12_boot_ci <- function(v, level = P12_CFG$inference$ci_level) {
    v <- v[is.finite(v)]
    if (length(v) < 20L) return(c(lo = NA_real_, hi = NA_real_))
    a <- (1 - level) / 2
    q <- stats::quantile(v, c(a, 1 - a), names = FALSE, type = 7)
    c(lo = q[1], hi = q[2])
}

# Two-sided bootstrap p: min(1, 2(k+1)/(B+1)), k = resamples on the far side of zero.
p12_boot_p <- function(v) {
    v <- v[is.finite(v)]
    if (length(v) < 20L) return(NA_real_)
    k <- min(sum(v <= 0), sum(v >= 0))
    min(1, 2 * (k + 1) / (length(v) + 1))
}

p12_boot_floor <- function(B = P12_CFG$inference$bootstrap_B) 2 / (B + 1)

p12_bh <- function(p) {
    out <- rep(NA_real_, length(p))
    ok <- is.finite(p)
    if (any(ok)) out[ok] <- stats::p.adjust(p[ok], method = "BH")
    out
}

p12_stars <- function(fdr) {
    s <- P12_CFG$inference$stars
    ifelse(!is.finite(fdr), "",
           ifelse(fdr < s[["**"]], "**", ifelse(fdr < s[["*"]], "*", "")))
}

# Text format for p and FDR; a value at the resampling floor is written with "<=".
p12_fmt <- function(x, floor = NA_real_) {
    out <- sprintf(P12_CFG$inference$text_number_format, x)
    at <- is.finite(x) & is.finite(floor) & x <= floor * (1 + 1e-9)
    out[at] <- paste0("<=", sprintf(P12_CFG$inference$text_number_format, floor))
    out[!is.finite(x)] <- "NA"
    out
}

p12_vif <- function(M) {
    out <- stats::setNames(rep(NA_real_, ncol(M)), colnames(M))
    C <- tryCatch(stats::cor(M, use = "complete.obs"), error = function(e) NULL)
    if (is.null(C) || any(!is.finite(C))) return(out)
    iv <- tryCatch(solve(C), error = function(e) NULL)
    if (is.null(iv)) return(stats::setNames(rep(Inf, ncol(M)), colnames(M)))
    d <- diag(iv)
    d[!is.finite(d)] <- Inf
    d[d < 1] <- 1
    stats::setNames(d, colnames(M))
}

# Within-donor centring for donors with at least `min_n` finite values.
p12_centre_within <- function(v, donor,
                              min_n = P12_CFG$floors$finite_units_within_donor) {
    out <- rep(NA_real_, length(v))
    for (u in unique(donor)) {
        i <- which(donor == u & is.finite(v))
        if (length(i) >= min_n) out[i] <- v[i] - mean(v[i])
    }
    out
}

p12_ascii_ok <- function(x) {
    x <- as.character(unlist(x))
    all(is.na(x) | !grepl("[^\\x01-\\x7F]", x, perl = TRUE))
}

# == Stats files ==============================================================
# Writes one panel's stats file: configured columns, BH within family, ASCII text.
p12_write_stats <- function(df, path) {
    want <- P12_CFG$stats_columns
    if (!identical(names(df), want))
        p12_abort(sprintf("stats columns differ from the config: have %s",
                          paste(names(df), collapse = ", ")))
    if (anyNA(df$family) || any(!nzchar(df$family)))
        p12_abort("every stats row needs a family")
    re <- rep(NA_real_, nrow(df))
    for (f in unique(df$family)) {
        i <- which(df$family == f)
        re[i] <- p12_bh(df$p[i])
    }
    same <- (is.na(re) & is.na(df$fdr)) |
        (!is.na(re) & !is.na(df$fdr) & abs(re - df$fdr) < 1e-12)
    if (!all(same)) p12_abort("FDR column is not BH within family")
    if (!p12_ascii_ok(df)) p12_abort("a stats field is not ASCII")
    p12_output(path)
    utils::write.csv(df, path, row.names = FALSE)
    invisible(path)
}

# == Pre-specified readings ===================================================
# Names a panel's pre-specified reading in the log before any computation; the
# full text is in the readings file.
p12_reading <- function(panel) {
    f <- p12_path(P12_CFG$readings_file)
    if (!file.exists(f)) p12_abort(paste0("readings file not found: ", f))
    x <- readLines(f, warn = FALSE, encoding = "UTF-8")
    head <- grep(sprintf("^### %s[. ]", panel), x)
    if (length(head) != 1L)
        p12_abort(sprintf("readings file has %d sections for panel %s", length(head), panel))
    p12_say(sprintf("Pre-set reading: %s | %s", sub("^### ", "", x[head]), p12_rel(f)))
}

# == Outputs listed in an earlier step's log ==================================
# Paths listed in the "Outputs:" block of a log (named by path).
p12_log_outputs <- function(log_path) {
    x <- readLines(log_path, warn = FALSE)
    if (!any(x == "VERDICT: PASS"))
        p12_abort(sprintf("%s does not end with VERDICT: PASS", basename(log_path)))
    m <- regmatches(x, regexec("^  (\\S+) \\| [0-9.]+ MB$", x))
    m <- m[lengths(m) == 2L]
    stats::setNames(vapply(m, `[`, "", 2L), vapply(m, `[`, "", 2L))
}

# Path of an earlier step's output, checked as listed in that step's log. `folder`
# is the step folder in Part12_rerun (e.g. "05_Units", "07_Fig15/15CD_DecontXStrata").
p12_upstream_file <- function(folder, file) {
    # Outputs are matched by file name in Part12_<folder>_log.txt.
    lg <- file.path(P12_RERUN, folder, sprintf("Part12_%s_log.txt", basename(folder)))
    if (!file.exists(lg)) p12_abort(sprintf("log not found: %s", p12_rel(lg)))
    want <- p12_log_outputs(lg)
    k <- names(want)[basename(names(want)) == file]
    f <- file.path(P12_RERUN, folder, file)
    ok <- length(k) == 1L && file.exists(f)
    p12_chk(ok, sprintf("%s listed in the %s log", file, folder), p12_rel(f))
    if (!ok) p12_abort(sprintf("%s is missing or not listed in its log", file))
    f
}

# Workers: all logical cores but one, at most `n_jobs`; P12_WORKERS overrides.
p12_n_workers <- function(n_jobs) {
    n <- if (exists("P12_WORKERS", envir = globalenv()))
        as.integer(get("P12_WORKERS", envir = globalenv())) else
        max(1L, parallel::detectCores(logical = TRUE) - 1L)
    max(1L, min(n, as.integer(n_jobs)))
}

# Splits 1..n into `k` contiguous chunks.
p12_chunks <- function(n, k) {
    k <- max(1L, min(as.integer(k), n))
    split(seq_len(n), ceiling(seq_len(n) * k / n))
}

