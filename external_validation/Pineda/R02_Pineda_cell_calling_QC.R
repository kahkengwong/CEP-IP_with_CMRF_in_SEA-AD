# =============================================================================
# R02_Pineda_cell_calling_QC.R
# External validation 1 (Pineda et al. 2024, BA9 batch 210526): cell calling and QC of
# each Cell Ranger raw matrix. Runs in native R on Windows, as the SEA-AD R parts did.
#
# Per library:
#   1. Reads counts/<LIB>_raw_feature_bc_matrix.h5.
#   2. emptyDrops (lower = 100, niters = 10,000); cells are barcodes with FDR <= 0.01.
#      These are the settings of SEA-AD Part12_06 (Gerrits cell calling).
#   3. Gene floor: called nuclei with fewer than 500 detected genes are removed.
#   4. scDblFinder with default settings, on the nuclei that pass the floor.
#   5. Saves R02_<LIB>.rds with the counts of the kept nuclei (called, at least 500 genes,
#      and singlet), the table of all called nuclei with their QC values and flags, and the
#      library summary row.
# No %MT cutoff is applied here. The authors' published QC is applied in R12.
#
# Random numbers: each library has its own seed (42 + its row in the manifest). The seed is
# set before emptyDrops and again before scDblFinder, with RNGkind fixed and
# BiocParallel::SerialParam(). The results therefore do not depend on the number of
# workers, the order in which libraries finish, or the OS.
#
# Resuming: libraries that already have R02_<LIB>.rds are skipped. Libraries whose raw .h5
# is not in counts/ yet are listed and skipped. Run the script again after copying more
# libraries. To process a library again, delete its R02_<LIB>.rds first.
# The summary workbook is rebuilt from all R02_<LIB>.rds files on every run.
#
# Usage: set ROOT and WORKERS below, then run the whole script in RStudio, or
#   source("C:/Users/Wong/Desktop/SEA-AD/External_Dataset_1/R02_Pineda_cell_calling_QC.R")
# While it runs, open R02_cell_calling/R02_progress.txt in Notepad++ to see each library
# as it finishes (the RStudio console is quiet until all workers return).
#
# Packages: DropletUtils, scDblFinder, SingleCellExperiment, BiocParallel, Matrix, openxlsx
#   BiocManager::install(c("DropletUtils", "scDblFinder"), update = FALSE)
#
# Outputs, all in ROOT/R02_cell_calling/:
#   R02_<LIB>.rds                    one per library
#   R02_cell_calling_summary.xlsx    sheets Libraries, By_group, Not_yet_processed, Settings, Session
#   R02_cell_calling_log.txt         log of every run (appended), with the session info
#   R02_progress.txt                 one line per library as the workers finish it
# =============================================================================

# The data folder: PINEDA_ROOT when set (the container sets it to /data), otherwise the Windows folder.
ROOT    <- Sys.getenv("PINEDA_ROOT", "C:/Users/Wong/Desktop/SEA-AD/External_Dataset_1")
# Libraries processed at the same time; each holds one raw matrix in memory. PINEDA_WORKERS
# overrides it. The results do not depend on the number of workers.
WORKERS <- as.integer(Sys.getenv("PINEDA_WORKERS", "6"))

SET <- list(lower = 100L, niters = 10000L, fdr = 0.01, min_genes = 500L, seed_base = 42L,
            ambient_umi_range = c(1L, 100L),
            rng_kind = c("Mersenne-Twister", "Inversion", "Rejection"))

MANIFEST <- file.path(ROOT, "Pineda_BA9_210526_manifest.tsv")
COUNTS   <- file.path(ROOT, "counts")
CR_MET   <- file.path(ROOT, "05_cellranger_metrics.csv")
OUT      <- file.path(ROOT, "R02_cell_calling")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
LOG      <- file.path(OUT, "R02_cell_calling_log.txt")
PROGRESS <- file.path(OUT, "R02_progress.txt")
SUMMARY  <- file.path(OUT, "R02_cell_calling_summary.xlsx")

# Writes one line to the console and to the log.
say <- function(...) {
    s <- paste0(...)
    cat(s, "\n", sep = "")
    cat(s, "\n", sep = "", file = LOG, append = TRUE)
}

need <- c("DropletUtils", "scDblFinder", "SingleCellExperiment", "SummarizedExperiment",
          "BiocParallel", "Matrix", "openxlsx")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Install these packages first: ", paste(miss, collapse = ", "))
if (!file.exists(MANIFEST)) stop("Manifest not found: ", MANIFEST)
if (!dir.exists(COUNTS)) stop("counts folder not found: ", COUNTS)

t0 <- Sys.time()
say("")
say("== R02_Pineda_cell_calling_QC.R started ", format(t0, "%Y-%m-%d %H:%M:%S"))
say("ROOT ", ROOT, " | workers ", WORKERS)
say("Settings: emptyDrops lower ", SET$lower, ", niters ", SET$niters, ", cells FDR <= ",
    SET$fdr, " | gene floor ", SET$min_genes, " | scDblFinder defaults | seed ",
    SET$seed_base, " + manifest row | RNGkind ", paste(SET$rng_kind, collapse = ", "))
say(R.version.string)
for (p in need) say("  ", p, " ", as.character(utils::packageVersion(p)))

# ---- Libraries ---------------------------------------------------------------
man <- utils::read.delim(MANIFEST, stringsAsFactors = FALSE, check.names = FALSE)
libs <- man[[1]]
say("Manifest: ", length(libs), " libraries (first column: ", names(man)[1], ")")
if (length(libs) != 63L) say("NOTE: the manifest does not hold 63 libraries")
if (anyDuplicated(libs)) stop("The manifest has duplicated library names.")

jobs <- data.frame(lib = libs, row = seq_along(libs),
                   h5 = file.path(COUNTS, paste0(libs, "_raw_feature_bc_matrix.h5")),
                   rds = file.path(OUT, paste0("R02_", libs, ".rds")),
                   stringsAsFactors = FALSE)
jobs$seed <- SET$seed_base + jobs$row
# A raw .h5 of 0 bytes is a copy that did not complete; it is treated as not copied.
h5_size <- file.size(jobs$h5)
have_h5 <- !is.na(h5_size) & h5_size > 0
empty_h5 <- !is.na(h5_size) & h5_size == 0
done <- file.exists(jobs$rds)
todo <- jobs[have_h5 & !done, ]
say(sprintf("Raw .h5 present %d | already processed %d | to process now %d | raw .h5 not yet copied %d",
            sum(have_h5), sum(done), nrow(todo), sum(!have_h5 & !done)))
if (any(empty_h5 & !done))
    say("WARNING: ", sum(empty_h5 & !done), " raw .h5 files are 0 bytes (the copy did not ",
        "complete) and are skipped: ", paste(jobs$lib[empty_h5 & !done], collapse = ", "))

# ---- One library (runs in a worker; uses only its arguments and package functions) ----
process_library <- function(job, set, progress) {
    t1 <- Sys.time()
    note <- function(s) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), s),
                            file = progress, append = TRUE)
    out <- tryCatch({
        RNGkind(set$rng_kind[1], set$rng_kind[2], set$rng_kind[3])

        sce <- DropletUtils::read10xCounts(job$h5, col.names = TRUE, delayed = FALSE)
        M <- methods::as(SingleCellExperiment::counts(sce), "dgCMatrix")
        genes <- data.frame(id = rownames(sce),
                            symbol = SummarizedExperiment::rowData(sce)$Symbol,
                            stringsAsFactors = FALSE)
        rm(sce)
        i_mt <- grep("^MT-", genes$symbol)
        if (length(i_mt) != 13L) stop(sprintf("%d MT- genes found, expected 13", length(i_mt)))

        # Ambient pool: barcodes with 1 to 100 UMIs (recorded only).
        tot <- Matrix::colSums(M)
        amb <- which(tot >= set$ambient_umi_range[1] & tot <= set$ambient_umi_range[2])
        amb_umis <- sum(tot[amb])
        amb_mt <- if (amb_umis > 0) sum(M[i_mt, amb]) / amb_umis else NA_real_

        # Cell calling.
        set.seed(job$seed)
        ed <- DropletUtils::emptyDrops(M, lower = set$lower, niters = set$niters,
                                       BPPARAM = BiocParallel::SerialParam())
        ok_fdr <- !is.na(ed$FDR)
        called <- which(ok_fdr & ed$FDR <= set$fdr)
        # Barcodes not called whose p-value hit the iteration limit (more iterations
        # could have called them).
        n_limited <- sum(ok_fdr & ed$FDR > set$fdr & ed$Limited)
        n_barcodes <- ncol(M)
        Mc <- M[, called, drop = FALSE]
        rm(M); gc(verbose = FALSE)

        tot_c <- Matrix::colSums(Mc)
        n_genes <- Matrix::colSums(Mc > 0)
        mt <- Matrix::colSums(Mc[i_mt, , drop = FALSE])
        nuc <- data.frame(library = job$lib, barcode = colnames(Mc), total_umis = tot_c,
                          n_genes = n_genes, mt_umis = mt, pct_mt = 100 * mt / tot_c,
                          emptydrops_FDR = ed$FDR[called],
                          pass_gene_floor = n_genes >= set$min_genes,
                          scDblFinder_score = NA_real_, scDblFinder_class = NA_character_,
                          keep = FALSE, stringsAsFactors = FALSE, row.names = NULL)

        # Doublets, on the nuclei that pass the gene floor.
        fl <- which(nuc$pass_gene_floor)
        if (length(fl) < 100L) stop(sprintf("only %d nuclei pass the gene floor", length(fl)))
        s2 <- SingleCellExperiment::SingleCellExperiment(list(counts = Mc[, fl, drop = FALSE]))
        set.seed(job$seed)
        s2 <- scDblFinder::scDblFinder(s2, BPPARAM = BiocParallel::SerialParam())
        nuc$scDblFinder_score[fl] <- s2$scDblFinder.score
        nuc$scDblFinder_class[fl] <- as.character(s2$scDblFinder.class)
        rm(s2)
        nuc$keep <- nuc$pass_gene_floor & nuc$scDblFinder_class %in% "singlet"

        K <- Mc[, nuc$keep, drop = FALSE]
        colnames(K) <- paste0(job$lib, "_", colnames(K))
        kept <- nuc[nuc$keep, ]
        row <- data.frame(
            library = job$lib, group = sub("^[0-9]+_([A-Za-z]+)_.*$", "\\1", job$lib),
            manifest_row = job$row, seed = job$seed,
            n_barcodes = n_barcodes, n_called = length(called),
            n_limited_not_called = n_limited,
            n_below_gene_floor = sum(!nuc$pass_gene_floor),
            n_doublets = sum(nuc$scDblFinder_class %in% "doublet"),
            pct_doublets = 100 * mean(nuc$scDblFinder_class[fl] == "doublet"),
            n_kept = nrow(kept),
            median_umis_kept = stats::median(kept$total_umis),
            median_genes_kept = stats::median(kept$n_genes),
            median_mt_umis_kept = stats::median(kept$mt_umis),
            median_pct_mt_kept = stats::median(kept$pct_mt),
            n_ambient_barcodes = length(amb), ambient_pct_mt = 100 * amb_mt,
            minutes = round(as.numeric(difftime(Sys.time(), t1, units = "mins")), 1),
            stringsAsFactors = FALSE)

        tmp <- paste0(job$rds, ".part")
        saveRDS(list(counts = K, genes = genes, nuclei = nuc, library = row, settings = set), tmp)
        if (!file.rename(tmp, job$rds)) stop("could not rename ", tmp)
        note(sprintf("DONE %s | called %d | below floor %d | doublets %.1f%% | kept %d | %.1f min",
                     job$lib, row$n_called, row$n_below_gene_floor, row$pct_doublets,
                     row$n_kept, row$minutes))
        list(ok = TRUE, lib = job$lib, row = row)
    }, error = function(e) {
        note(sprintf("FAILED %s | %s", job$lib, conditionMessage(e)))
        list(ok = FALSE, lib = job$lib, msg = conditionMessage(e))
    })
    gc(verbose = FALSE)
    out
}

# ---- Run ---------------------------------------------------------------------
if (nrow(todo)) {
    job_list <- lapply(seq_len(nrow(todo)), function(i) as.list(todo[i, ]))
    n_work <- max(1L, min(as.integer(WORKERS), nrow(todo)))
    say("Processing ", nrow(todo), " libraries with ", n_work, " workers. Progress: ", PROGRESS)
    cat(sprintf("== run started %s, %d libraries, %d workers\n",
                format(t0, "%Y-%m-%d %H:%M:%S"), nrow(todo), n_work),
        file = PROGRESS, append = TRUE)
    if (n_work > 1L) {
        cl <- parallel::makePSOCKcluster(n_work)
        res <- tryCatch(parallel::parLapplyLB(cl, job_list, process_library, set = SET,
                                              progress = PROGRESS),
                        finally = parallel::stopCluster(cl))
    } else {
        res <- lapply(job_list, process_library, set = SET, progress = PROGRESS)
    }
    for (r in res) {
        if (isTRUE(r$ok)) {
            x <- r$row
            say(sprintf("  DONE   %-28s called %6d | below floor %5d | doublets %5.1f%% | kept %6d | median UMIs %6.0f | median %%MT %.2f | %.1f min",
                        x$library, x$n_called, x$n_below_gene_floor, x$pct_doublets,
                        x$n_kept, x$median_umis_kept, x$median_pct_mt_kept, x$minutes))
        } else {
            say("  FAILED ", r$lib, " | ", r$msg)
        }
    }
} else {
    say("No library to process in this run.")
}

# ---- Summary workbook (all processed libraries) -----------------------------------
say("Rebuilding the summary from all R02_<LIB>.rds files")
rds_all <- jobs$rds[file.exists(jobs$rds)]
rows <- do.call(rbind, lapply(rds_all, function(f) readRDS(f)$library))
if (!is.null(rows) && file.exists(CR_MET)) {
    cr <- utils::read.csv(CR_MET, check.names = FALSE, stringsAsFactors = FALSE)
    if ("Estimated Number of Cells" %in% names(cr)) {
        est <- as.numeric(gsub(",", "", cr[["Estimated Number of Cells"]]))
        rows$cellranger_estimated_cells <- est[match(rows$library, cr[[1]])]
    }
}
not_yet <- jobs[!file.exists(jobs$rds), c("lib", "row")]
not_yet$raw_h5_copied <- have_h5[not_yet$row]
names(not_yet) <- c("library", "manifest_row", "raw_h5_copied")

by_group <- NULL
if (!is.null(rows)) {
    by_group <- do.call(rbind, lapply(split(rows, rows$group), function(d) data.frame(
        group = d$group[1], libraries = nrow(d), nuclei_called = sum(d$n_called),
        nuclei_kept = sum(d$n_kept), median_kept_per_library = stats::median(d$n_kept),
        median_pct_doublets = stats::median(d$pct_doublets),
        median_umis_kept = stats::median(d$median_umis_kept),
        median_pct_mt_kept = stats::median(d$median_pct_mt_kept),
        libraries_with_limited_barcodes = sum(d$n_limited_not_called > 0))))
    say(sprintf("Processed %d of %d libraries | nuclei called %d | kept %d",
                nrow(rows), nrow(jobs), sum(rows$n_called), sum(rows$n_kept)))
    for (i in seq_len(nrow(by_group))) {
        g <- by_group[i, ]
        say(sprintf("  %-5s libraries %2d | kept %7d (median %6.0f per library) | median doublets %.1f%% | median UMIs %.0f | median %%MT %.2f",
                    g$group, g$libraries, g$nuclei_kept, g$median_kept_per_library,
                    g$median_pct_doublets, g$median_umis_kept, g$median_pct_mt_kept))
    }
    if (any(rows$n_limited_not_called > 0))
        say("NOTE: ", sum(rows$n_limited_not_called > 0), " libraries have barcodes that were ",
            "not called and whose p-value hit the 10,000-iteration limit (column n_limited_not_called).")
}

settings <- data.frame(
    setting = c("emptyDrops lower", "emptyDrops niters", "cells: emptyDrops FDR at most",
                "gene floor (detected genes, at least)", "doublets", "%MT cutoff",
                "seed", "RNGkind", "BiocParallel", "ambient pool (recorded only)",
                "kept nuclei"),
    value = c(SET$lower, SET$niters, SET$fdr, SET$min_genes,
              "scDblFinder, default settings, per library, on nuclei passing the gene floor",
              "none", paste(SET$seed_base, "+ manifest row, set before emptyDrops and before scDblFinder"),
              paste(SET$rng_kind, collapse = ", "), "SerialParam()",
              sprintf("barcodes with %d to %d UMIs", SET$ambient_umi_range[1], SET$ambient_umi_range[2]),
              "called, at least 500 genes, and scDblFinder singlet"),
    stringsAsFactors = FALSE)

wb <- openxlsx::createWorkbook()
bold <- openxlsx::createStyle(textDecoration = "bold")
add_sheet <- function(name, d) {
    openxlsx::addWorksheet(wb, name)
    if (is.null(d) || !nrow(d)) d <- data.frame(note = "none")
    openxlsx::writeData(wb, name, d, headerStyle = bold)
    openxlsx::freezePane(wb, name, firstRow = TRUE)
}
add_sheet("Libraries", rows)
add_sheet("By_group", by_group)
add_sheet("Not_yet_processed", not_yet)
add_sheet("Settings", settings)
add_sheet("Session", data.frame(line = utils::capture.output(utils::sessionInfo())))
openxlsx::saveWorkbook(wb, SUMMARY, overwrite = TRUE)
say("Summary: ", SUMMARY)

say("Session info:")
for (l in utils::capture.output(utils::sessionInfo())) say("  ", l)
say(sprintf("== Finished %s (%.1f min)", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
