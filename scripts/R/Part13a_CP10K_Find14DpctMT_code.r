# =============================================================================
# Part13a_CP10K_Find14DpctMT_code.r
# SEA-AD CEP-IP/CMRF, Part 13a, donor %MT for the Figure 14D model
#   Searches the project folder for donor %MT columns, fits the Figure 14D model
#   (outcome ~ cov4 + Module post-IP% + %MT) on each candidate as in Part 13,
#   and writes the one set of donor values that reproduces the published 14D
#   rows (n exactly, beta within TOL_BETA) as the MT_FILE that Part 13 reads.
#
# Usage (from the SEA-AD folder, or set SEAAD_ROOT), before Part 13
#   source("Part13a_CP10K_Find14DpctMT_code.r")
#
# Outputs (Part13_M4Share/): Part13_M4_pctMT_for14D.csv (only when exactly one
#   set passes), ranked candidates CSV, and log.
# =============================================================================

# == CONFIGURATION ===========================================================
if (!exists("CEP_HOME")) CEP_HOME <- Sys.getenv("SEAAD_ROOT", getwd())
OUT_DIR       <- file.path(CEP_HOME, "Part13_M4Share")
COMP_RDS      <- NA     # NA = search CEP_HOME for COMP_BASENAME
COMP_BASENAME <- "Part10h_CP10K_M4_Comparator_results.rds"
FILE_REGEX    <- "(Part ?12|12q|Battery|Fig ?14|SuppTable_Fig14)"   # matched on the full path
SEARCH_EXCLUDE <- c("Previous_version", "Archive", "archive", "Part13_M4Share")
MAX_MB        <- 300
MT_NAME_REGEX <- "mito|pctmt|mtpct|percentmt|fracmt|mtfrac|mtlog|logmt|pmt|^mt$|^mtlog1p$"
ARM_KEEP      <- c("m4", "mic13pvm3", "mic13")           # normalised arm values kept
ARM_OTHER     <- c("pvm3all", "pvm3", "pvm1all", "pvm1", "hkg8matched", "hkg8")
DATASETS      <- c("NeoC", "NonNeoC")
MIN_DONORS    <- 20L
TOL_BETA      <- 1.5e-3
# Published Figure 14D rows: dataset, outcome label, beta, n
PUB14D <- data.frame(dataset = c("NeoC", "NonNeoC", "NonNeoC"),
                     outcome = c("MMSE score", "Cognitive status", "MMSE score"),
                     beta = c(-0.529, 0.415, -0.337), n = c(74L, 80L, 76L),
                     fdr = c(0.003, 0.015, 0.028), stringsAsFactors = FALSE)
CLIN_LAB <- c("Cognitive status", "MMSE score", "CASI score", "CERAD score", "Braak stage",
              "Thal phase", "ADNC rank", "CAA score", "APOE e4 dose")
REGION_SET <- c(mtg = "NeoC", stg = "NeoC", itg = "NeoC", dfc = "NeoC", a9 = "NeoC", ba9 = "NeoC",
                ang = "NeoC", v1c = "NeoC", v1 = "NeoC", middletemporalgyrus = "NeoC",
                superiortemporalgyrus = "NeoC", inferiortemporalgyrus = "NeoC",
                dorsolateralprefrontalcortex = "NeoC", angulargyrus = "NeoC",
                primaryvisualcortex = "NeoC",
                hip = "NonNeoC", mec = "NonNeoC", lec = "NonNeoC", fi = "NonNeoC",
                hippocampus = "NonNeoC", medialentorhinalcortex = "NonNeoC",
                lateralentorhinalcortex = "NonNeoC", frontoinsula = "NonNeoC",
                frontalinsula = "NonNeoC")

# == LOG =====================================================================
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
STAMP   <- format(Sys.time(), "%Y%m%d_%H%M%S")
CONSOLE <- file.path(OUT_DIR, sprintf("Part13a_CP10K_Find14DpctMT_console_log_%s.txt", STAMP))
CAND_CSV <- file.path(OUT_DIR, "Part13a_CP10K_Find14DpctMT_candidates.csv")
OUT_MT  <- file.path(OUT_DIR, "Part13_M4_pctMT_for14D.csv")
.con <- file(CONSOLE, open = "wt"); sink(.con, split = TRUE)
cat(strrep("=", 78), "\n  Part 13a - find the Figure 14D %MT\n  ", format(Sys.time()), "\n  ",
    R.version.string, "\n  Root: ", CEP_HOME, "\n", strrep("=", 78), "\n", sep = "")

# == HELPERS =================================================================
nrm <- function(x) gsub("[^a-z0-9]", "", tolower(as.character(x)))
num <- function(x) suppressWarnings(as.numeric(x))
.z <- function(x) { x <- num(x); s <- stats::sd(x, na.rm = TRUE)
    if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x))); (x - mean(x, na.rm = TRUE)) / s }
excluded <- function(p) Reduce(`|`, lapply(SEARCH_EXCLUDE, function(s) grepl(s, p, fixed = TRUE)),
                               rep(FALSE, length(p)))

fit_14d <- function(D, y, covs) {           # Part 13 fit_ols, reduced to the module row
    covs <- setdiff(covs, y); vars <- unique(c(y, "module", "mt", covs))
    X <- as.data.frame(lapply(D[, vars, drop = FALSE], num)); names(X) <- vars
    X <- X[stats::complete.cases(X), , drop = FALSE]
    if (nrow(X) < MIN_DONORS) return(c(n = nrow(X), beta = NA, p = NA))
    X[[y]] <- .z(X[[y]]); X$module <- .z(X$module); X$mt <- .z(X$mt)
    if (anyNA(X[, c(y, "module", "mt")])) return(c(n = nrow(X), beta = NA, p = NA))
    f <- stats::as.formula(paste(sprintf("`%s`", y), "~ module + mt +",
                                 paste(sprintf("`%s`", covs), collapse = " + ")))
    cf <- summary(stats::lm(f, data = X))$coefficients
    if (!"module" %in% rownames(cf)) return(c(n = nrow(X), beta = NA, p = NA))
    c(n = nrow(X), beta = cf["module", 1], p = cf["module", 4])
}

read_tables <- function(path) {
    ext <- tolower(tools::file_ext(path)); out <- list()
    tryCatch({
        if (ext == "csv") {
            out[["csv"]] <- if (requireNamespace("data.table", quietly = TRUE))
                as.data.frame(data.table::fread(path, showProgress = FALSE)) else
                utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
        } else if (ext == "xlsx") {
            for (sh in openxlsx::getSheetNames(path))
                out[[paste0("sheet:", sh)]] <- openxlsx::read.xlsx(path, sheet = sh)
        } else if (ext == "rds") {
            walk <- function(x, nm, depth) {
                if (is.data.frame(x)) { out[[nm]] <<- as.data.frame(x); return(invisible()) }
                if (is.list(x) && depth < 3 && length(x) < 500)
                    for (k in seq_along(x)) walk(x[[k]], paste0(nm, "$", names(x)[k] %||% k), depth + 1)
            }
            walk(readRDS(path), "rds", 0)
        }
    }, error = function(e) cat(sprintf("  unreadable: %s (%s)\n", path, conditionMessage(e))))
    out
}
`%||%` <- function(a, b) if (is.null(a) || !length(a) || is.na(a) || !nzchar(a)) b else a

# == STEP 1: the comparator donor frame (the frame Part 13 fits on) ==========
cat("\nSTEP 1 -- comparator donor frame\n")
comp_path <- if (!is.na(COMP_RDS)) COMP_RDS else {
    h <- list.files(CEP_HOME, pattern = paste0("^", gsub(".", "\\.", COMP_BASENAME, fixed = TRUE), "$"),
                    recursive = TRUE, full.names = TRUE)
    h <- h[!excluded(h)]
    if (length(h) != 1L) stop(sprintf("%s: %d match(es) after exclusions; set COMP_RDS. %s", COMP_BASENAME,
                                      length(h), paste(h, collapse = "; ")), call. = FALSE)
    h
}
RC <- readRDS(comp_path)
DF <- as.data.frame(RC$donor_frame)
DF$.donor <- as.character(DF$.donor); DF$dataset <- as.character(DF$dataset); DF$module <- num(DF$module)
OUT_MAP <- unlist(RC$meta$outcomes); COV4 <- unlist(RC$meta$covariate_columns$cov4)
cat(sprintf("  %s | donors %s\n", comp_path,
            paste(table(DF$dataset)[DATASETS], collapse = "/")))
donor_key <- nrm(DF$.donor)

# == STEP 2: candidate files =================================================
cat("\nSTEP 2 -- candidate files\n")
files <- list.files(CEP_HOME, pattern = "\\.(csv|xlsx|rds)$", recursive = TRUE, full.names = TRUE,
                    ignore.case = TRUE)
files <- files[grepl(FILE_REGEX, files, ignore.case = TRUE) & !excluded(files)]
big <- file.size(files) / 1e6 > MAX_MB
if (any(big)) cat(sprintf("  %d file(s) above %d MB skipped: %s\n", sum(big), MAX_MB,
                          paste(basename(files[big]), collapse = "; ")))
files <- files[!big]
cat(sprintf("  %d file(s) to scan\n", length(files)))

# == STEP 3: build candidate donor-level %MT versions ========================
cat("\nSTEP 3 -- candidate %MT columns\n")
VERS <- list()
add_version <- function(src, vals) {         # vals: data.frame(.donor, dataset, mt)
    vals <- vals[is.finite(vals$mt), , drop = FALSE]
    if (nrow(vals) < MIN_DONORS) return(invisible())
    VERS[[length(VERS) + 1L]] <<- list(src = src, vals = vals)
}
for (f in files) {
    tabs <- read_tables(f)
    path_ds <- if (grepl("nonneoc", tolower(f))) "NonNeoC" else if (grepl("neoc", tolower(f))) "NeoC" else NA
    path_arm <- if (grepl("(^|[^a-z])(m4|mic13)", tolower(basename(f)))) "M4" else
                if (grepl("pvm_?3|pvm_?1|hkg8", tolower(basename(f)))) "other" else NA
    for (tn in names(tabs)) {
        t <- tabs[[tn]]
        if (!is.data.frame(t) || nrow(t) < MIN_DONORS || ncol(t) < 2) next
        nm <- names(t)
        # donor column: the column sharing the most values with the comparator donors
        ov <- vapply(nm, function(c) sum(unique(nrm(t[[c]])) %in% donor_key), numeric(1))
        if (max(ov) < MIN_DONORS) next
        dcol <- nm[which.max(ov)]
        # region set: a column holding NeoC/NonNeoC, else a region column, else the path
        ds <- rep(NA_character_, nrow(t)); ds_src <- NA
        for (c in nm) {
            v <- nrm(t[[c]]); okv <- v[!is.na(v) & nzchar(v)]
            if (length(okv) && mean(okv %in% c("neoc", "nonneoc")) > 0.95) {
                ds <- ifelse(v == "neoc", "NeoC", ifelse(v == "nonneoc", "NonNeoC", NA)); ds_src <- c; break }
        }
        if (is.na(ds_src)) for (c in nm) {
            v <- nrm(t[[c]]); okv <- v[!is.na(v) & nzchar(v)]
            if (length(okv) && mean(okv %in% names(REGION_SET)) > 0.9) {
                ds <- unname(REGION_SET[v]); ds_src <- paste0("region:", c); break }
        }
        if (is.na(ds_src) && !is.na(path_ds)) { ds <- rep(path_ds, nrow(t)); ds_src <- "path" }
        if (is.na(ds_src)) next
        # arm
        keep <- rep(TRUE, nrow(t)); arm_src <- "none"
        for (c in nm) {
            v <- nrm(t[[c]]); okv <- unique(v[!is.na(v) & nzchar(v)])
            if (length(okv) >= 2 && length(okv) <= 8 && any(okv %in% c(ARM_KEEP, ARM_OTHER))) {
                keep <- v %in% ARM_KEEP; arm_src <- c; break }
        }
        if (arm_src == "none" && identical(path_arm, "other")) next
        if (!any(keep)) next
        # cell-count column for weighted means
        wcol <- nm[nrm(nm) %in% c("ncells", "ncell", "cells", "narm", "ncellsarm", "cellcount", "nunit", "n")][1]
        mtcols <- nm[grepl(MT_NAME_REGEX, nrm(nm)) & vapply(t, function(x) sum(is.finite(num(x))) >= MIN_DONORS,
                                                            logical(1))]
        for (mc in mtcols) {
            base <- data.frame(.donor = DF$.donor[match(nrm(t[[dcol]]), donor_key)], dataset = ds,
                               mt = num(t[[mc]]), w = if (!is.na(wcol)) num(t[[wcol]]) else NA_real_,
                               stringsAsFactors = FALSE)[keep, , drop = FALSE]
            base <- base[!is.na(base$.donor) & !is.na(base$dataset) & is.finite(base$mt), , drop = FALSE]
            if (nrow(base) < MIN_DONORS) next
            src <- sprintf("%s | %s | col %s | donor %s | set %s | arm %s", f, tn, mc, dcol, ds_src, arm_src)
            key <- paste(base$.donor, base$dataset)
            if (!any(duplicated(key))) { add_version(paste(src, "| as stored"), base[, 1:3]); next }
            g <- split(base, key)
            agg <- function(fun) do.call(rbind, lapply(g, function(x)
                data.frame(.donor = x$.donor[1], dataset = x$dataset[1], mt = fun(x), stringsAsFactors = FALSE)))
            add_version(paste(src, "| donor mean"), agg(function(x) mean(x$mt)))
            add_version(paste(src, "| donor median"), agg(function(x) stats::median(x$mt)))
            if (!is.na(wcol)) add_version(paste(src, "| donor mean weighted by", wcol),
                                          agg(function(x) stats::weighted.mean(x$mt, x$w, na.rm = TRUE)))
        }
    }
}
if ("MitoFrac_num" %in% names(DF))
    add_version("comparator donor_frame | MitoFrac_num (Micro-PVM_3-wide donor mean) | as stored",
                data.frame(.donor = DF$.donor, dataset = DF$dataset, mt = num(DF$MitoFrac_num)))
cat(sprintf("  %d candidate column version(s) before transforms\n", length(VERS)))
if (!length(VERS)) stop("No candidate %MT column was found. Widen FILE_REGEX or MT_NAME_REGEX.", call. = FALSE)

# == STEP 4: fit the Figure 14D model on every version =======================
cat("\nSTEP 4 -- Figure 14D reproduction\n")
TRANS <- list(none = function(x) x, log1p = function(x) log1p(x),
              log1p_x100 = function(x) log1p(100 * x), log = function(x) ifelse(x > 0, log(x), NA))
rows <- list(); vals_store <- list()
for (i in seq_along(VERS)) for (tr in names(TRANS)) {
    v <- VERS[[i]]$vals
    if (tr == "log1p_x100" && max(v$mt, na.rm = TRUE) > 1) next
    v$mt <- TRANS[[tr]](v$mt)
    D <- merge(DF, v, by = c(".donor", "dataset"), all.x = TRUE)
    p9 <- lapply(DATASETS, function(ds) {
        d <- D[D$dataset == ds, ]
        r <- t(vapply(CLIN_LAB, function(lab) fit_14d(d, OUT_MAP[[lab]], COV4), numeric(3)))
        data.frame(dataset = ds, outcome = CLIN_LAB, n = r[, "n"], beta = r[, "beta"], p = r[, "p"],
                   fdr = stats::p.adjust(r[, "p"], "BH"), stringsAsFactors = FALSE)
    })
    p9 <- do.call(rbind, p9)
    m <- merge(PUB14D, p9, by = c("dataset", "outcome"), suffixes = c("_pub", "_fit"))
    dbeta <- max(abs(m$beta_fit - m$beta_pub)); n_ok <- all(m$n_fit == m$n_pub, na.rm = FALSE)
    id <- sprintf("V%03d_%s", i, tr)
    vals_store[[id]] <- v
    rows[[length(rows) + 1L]] <- data.frame(
        id = id, pass = isTRUE(n_ok && dbeta <= TOL_BETA), max_abs_dbeta = dbeta,
        n_fit = paste(m$n_fit, collapse = "/"), n_pub = paste(m$n_pub, collapse = "/"),
        beta_fit = paste(sprintf("%.3f", m$beta_fit), collapse = "/"),
        beta_pub = paste(sprintf("%.3f", m$beta_pub), collapse = "/"),
        fdr_fit = paste(sprintf("%.3f", m$fdr_fit), collapse = "/"),
        fdr_pub = paste(sprintf("%.3f", m$fdr_pub), collapse = "/"),
        transform = tr, source = VERS[[i]]$src, stringsAsFactors = FALSE)
}
CAND <- do.call(rbind, rows)
CAND <- CAND[order(!CAND$pass, CAND$max_abs_dbeta), ]
utils::write.csv(CAND, CAND_CSV, row.names = FALSE)
cat(sprintf("  %d version(s) tried, %d pass. Rows are published-row order: %s\n", nrow(CAND), sum(CAND$pass),
            paste(PUB14D$dataset, PUB14D$outcome, collapse = "; ")))
print(utils::head(CAND[, c("id", "pass", "max_abs_dbeta", "n_fit", "beta_fit", "transform", "source")], 10),
      row.names = FALSE)

# == STEP 5: write MT_FILE when exactly one set of values passes =============
cat("\nSTEP 5 -- result\n")
ok <- CAND[CAND$pass, , drop = FALSE]
sig <- function(id) { v <- vals_store[[id]]; v <- v[order(v$dataset, v$.donor), ]
    paste(v$.donor, v$dataset, signif(v$mt, 10), collapse = ";") }
uniq <- if (nrow(ok)) unique(vapply(ok$id, sig, "")) else character(0)
if (length(uniq) == 1L) {
    v <- vals_store[[ok$id[1]]]
    out <- data.frame(donor_id = v$.donor, dataset = v$dataset, pct_mt = v$mt)
    out <- out[order(out$dataset, out$donor_id), ]
    utils::write.csv(out, OUT_MT, row.names = FALSE)
    cat(sprintf("  One set of donor values reproduces Figure 14D (%d equivalent source(s)).\n", nrow(ok)))
    for (j in seq_len(nrow(ok))) cat(sprintf("    %s [%s]\n", ok$source[j], ok$transform[j]))
    cat(sprintf("  Written: %s | %d donor rows\n", OUT_MT, nrow(out)))
} else if (length(uniq) > 1L) {
    cat(sprintf("  %d different sets of donor values pass; nothing written. See %s.\n",
                length(uniq), basename(CAND_CSV)))
} else {
    cat(sprintf(paste("  No version reproduces Figure 14D; nothing written. The closest are in %s.",
                      "\n"), basename(CAND_CSV)))
}
cat("\nConsole log: ", CONSOLE, "\n", sep = "")
sink(); close(.con)
