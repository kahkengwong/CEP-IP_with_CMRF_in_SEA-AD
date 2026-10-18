# =============================================================================
# Part12_06_GerritsCells.r
# SEA-AD CEP-IP/CMRF, Part 12, step 06: Gerrits (GSE148822) cells
#   Builds the sample table (36 single-nucleus samples, 18 donors) from the
#   series matrix, matches each deposit matrix to one sample, calls cells with
#   emptyDrops (fixed seed per sample), and records per nucleus depth, %MT,
#   MT-CO3 and its partners, lineage markers and the soup index. Samples run in
#   parallel; set P12_WORKERS before sourcing to limit workers.
#
# Usage (from the SEA-AD folder, or set SEAAD_ROOT)
#   source("Part12_rerun/Part12_06_GerritsCells.r")
#
# Outputs (Part12_rerun/06_GerritsCells/): 06_gerrits_samples.csv, 06_gerrits_cells.csv.gz, log.
# =============================================================================

source(file.path(Sys.getenv("SEAAD_ROOT", getwd()), "Part12_rerun", "Part12_03_Config.r"))
if (!requireNamespace("DropletUtils", quietly = TRUE))
    stop("DropletUtils is required: BiocManager::install(\"DropletUtils\")", call. = FALSE)

STEP <- "06"
D <- p12_step_dirs(STEP, "GerritsCells")
p12_clean(c(D$tables, D$figures))
p12_start(STEP, "Part12_06_GerritsCells.r", D, pkgs = c("DropletUtils", "BiocParallel"))

F_SAMPLES <- file.path(D$tables, "06_gerrits_samples.csv")
F_CELLS <- file.path(D$tables, "06_gerrits_cells.csv.gz")

# Reads the !Sample_ rows of a GEO series matrix into a list of vectors.
read_series_matrix <- function(path) {
    con <- gzfile(path, "rt")
    on.exit(close(con))
    x <- readLines(con, warn = FALSE)
    x <- x[startsWith(x, "!Sample_")]
    out <- list()
    for (l in x) {
        f <- strsplit(l, "\t", fixed = TRUE)[[1]]
        key <- f[1]
        val <- gsub('^"|"$', "", f[-1])
        if (key == "!Sample_characteristics_ch1") {
            lab <- sub(":.*$", "", val[1])
            key <- paste0(key, "|", lab)
            val <- trimws(sub("^[^:]*:", "", val))
        }
        if (is.null(out[[key]])) out[[key]] <- val
    }
    out
}

# One sample, start to finish, in a worker: sets its own RNG kind and seed and
# returns an error as a value.
p12_gerrits_sample <- function(job) {
    out <- list(ok = FALSE, sample = job$sample, msg = "")
    res <- tryCatch({
        RNGkind(job$rng_kind)
        M <- Matrix::readMM(job$path)
        transposed <- FALSE
        if (nrow(M) != job$n_genes && ncol(M) == job$n_genes) {
            M <- Matrix::t(M)
            transposed <- TRUE
        }
        if (nrow(M) != job$n_genes)
            stop(sprintf("matrix is %d x %d for %d genes", nrow(M), ncol(M), job$n_genes))
        M <- methods::as(M, "CsparseMatrix")
        tot <- Matrix::colSums(M)
        amb <- which(tot >= job$amb[1] & tot <= job$amb[2])
        amb_tot <- sum(tot[amb])
        amb_mt <- if (amb_tot > 0) sum(M[job$i_mt, amb]) / amb_tot else NA_real_

        set.seed(job$seed)
        ed <- DropletUtils::emptyDrops(M, lower = job$lower, niters = job$niters,
                                       BPPARAM = BiocParallel::SerialParam())
        keep <- which(!is.na(ed$FDR) & ed$FDR <= job$fdr)
        srow <- data.frame(
            sample = job$sample, gsm = job$gsm, donor = job$donor, region = job$region,
            group = job$group, file_stem = job$file_stem, n_barcodes = ncol(M),
            n_called = length(keep), n_ambient_barcodes = length(amb),
            ambient_umis = amb_tot, ambient_mito_frac = amb_mt,
            median_umis = NA_real_, median_mito_frac = NA_real_,
            n_goi_positive = 0L, seed = job$seed)
        cells <- NULL
        if (length(keep)) {
            Mc <- M[, keep, drop = FALSE]
            totc <- tot[keep]
            mt <- Matrix::colSums(Mc[job$i_mt, , drop = FALSE])
            lin <- vapply(job$lin_i, function(g)
                if (length(g)) Matrix::colSums(Mc[g, , drop = FALSE]) else rep(0, ncol(Mc)),
                numeric(ncol(Mc)))
            if (!is.matrix(lin))
                lin <- matrix(lin, nrow = 1L, dimnames = list(NULL, names(job$lin_i)))
            best <- max.col(lin, ties.method = "first")
            lineage <- ifelse(rowSums(lin) > 0, colnames(lin)[best], NA_character_)
            own <- lin[cbind(seq_len(nrow(lin)), best)]
            own[is.na(lineage)] <- 0
            cells <- data.frame(sample = job$sample, gsm = job$gsm, donor = job$donor,
                                region = job$region, group = job$group,
                                cell = keep, total_umis = totc,
                                n_genes = Matrix::colSums(Mc > 0), mt_umis = mt,
                                mito_frac = mt / totc,
                                goi = as.numeric(Mc[job$i_goi, ]))
            for (g in names(job$i_dfg))
                cells[[paste0("dfg_", g)]] <- as.numeric(Mc[job$i_dfg[[g]], ])
            for (k in colnames(lin)) cells[[paste0("lin_", k)]] <- lin[, k]
            cells$lineage <- lineage
            cells$soup_index <- (rowSums(lin) - own) / totc
            srow$median_umis <- stats::median(totc)
            srow$median_mito_frac <- stats::median(mt / totc)
            srow$n_goi_positive <- sum(cells$goi > 0)
        }
        list(ok = TRUE, sample = job$sample, transposed = transposed,
             srow = srow, cells = cells)
    }, error = function(e) list(ok = FALSE, sample = job$sample, msg = conditionMessage(e)))
    gc(verbose = FALSE)
    res
}

main <- function() {
    G <- P12_CFG$gerrits
    CC <- G$cell_calling
    AMB <- as.integer(G$ambient_umi_range)

    # == A. sample table from the series matrix ==============================
    p12_say("SECTION A -- samples and donors (GSE148822 series matrix)")
    f_sm <- p12_log_input(p12_input("gerrits_series_matrix"), "series matrix")
    SM <- read_series_matrix(f_sm)
    need <- c("!Sample_title", "!Sample_geo_accession",
              "!Sample_characteristics_ch1|nuclei population",
              "!Sample_characteristics_ch1|region",
              "!Sample_characteristics_ch1|sample group")
    miss <- setdiff(need, names(SM))
    if (length(miss)) p12_abort(paste("series matrix lacks:", paste(miss, collapse = ", ")))
    S <- data.frame(title = SM[["!Sample_title"]], gsm = SM[["!Sample_geo_accession"]],
                    population = SM[["!Sample_characteristics_ch1|nuclei population"]],
                    region = SM[["!Sample_characteristics_ch1|region"]],
                    group = SM[["!Sample_characteristics_ch1|sample group"]])
    S <- S[grepl(G$sample_title_regex, S$title), ]
    rownames(S) <- NULL
    m <- regmatches(S$title, regexec(G$sample_title_regex, S$title))
    S$donor <- paste0("GerritsDonor", vapply(m, `[`, "", 2L))
    S$region_code <- vapply(m, `[`, "", 3L)
    p12_chk(nrow(S) == G$n_samples, "single-nucleus samples in the series matrix",
            sprintf("%d (expected %d)", nrow(S), G$n_samples))
    dpn <- table(S$donor)
    p12_chk(length(dpn) == G$n_donors && all(dpn == G$samples_per_donor),
            "donors and samples per donor",
            sprintf("%d donors | samples per donor: %s", length(dpn),
                    paste(sort(unique(as.integer(dpn))), collapse = ",")))
    p12_chk(all(S$population == G$sorted_fraction), "every sample is the sorted fraction",
            paste(unique(S$population), collapse = ", "))
    p12_say(sprintf("  groups: %s", paste(sprintf("%s %d", names(table(S$group)),
                                                   table(S$group)), collapse = " | ")))
    p12_say("")

    # == B. deposit files ====================================================
    p12_say("SECTION B -- deposit matrices")
    dep <- p12_input("gerrits_deposit")
    mtx <- sort(list.files(dep, pattern = paste0(gsub(".", "\\.", G$matrix_suffix,
                                                      fixed = TRUE), "$"),
                           full.names = TRUE))
    mtx <- mtx[!startsWith(basename(mtx), P12_CFG$inputs$skip_file_prefix)]
    stem <- substring(basename(mtx), 1L, nchar(basename(mtx)) - nchar(G$matrix_suffix))
    # A file belongs to one sample by its GSM accession, else by sample title.
    pick <- vapply(stem, function(st) {
        h <- which(vapply(S$gsm, function(g) grepl(paste0("(^|_)", g, "(_|$)"), st), TRUE))
        if (!length(h))
            # Deposit names drop the "SN_" prefix of the GEO title.
            h <- which(vapply(sub("^SN_", "", S$title), function(t)
                grepl(paste0("(^|_)", t, "(_|$)"), st), TRUE))
        if (length(h) == 1L) h else NA_integer_
    }, 0L)
    snrna <- !is.na(pick)
    p12_note("matrix files in the deposit",
             sprintf("%d | %d match a single-nucleus sample", length(mtx), sum(snrna)))
    mtx <- mtx[snrna]; stem <- stem[snrna]; pick <- pick[snrna]
    p12_chk(length(pick) == nrow(S) && !anyDuplicated(pick),
            "each single-nucleus sample has exactly one matrix file",
            sprintf("%d files for %d samples | %d duplicated", length(pick), nrow(S),
                    sum(duplicated(pick))))
    if (P12$n[["FAIL"]] > 0L) p12_abort("sample table or deposit files do not match")
    S$file_stem <- NA_character_
    S$file_stem[pick] <- stem
    feat <- file.path(dep, paste0(S$file_stem, G$features_suffix))
    p12_chk(all(file.exists(feat)), "a features file for every sample",
            sprintf("%d of %d", sum(file.exists(feat)), length(feat)))
    if (!all(file.exists(feat))) p12_abort("features files missing")

    # Every sample must carry the same feature list.
    sym <- function(f) as.character(utils::read.delim(gzfile(f), header = FALSE,
                                                      stringsAsFactors = FALSE)[[G$features_symbol_column]])
    genes <- sym(feat[1])
    same <- vapply(feat[-1], function(f) identical(sym(f), genes), TRUE)
    p12_chk(all(same), "every sample has the same feature list",
            sprintf("%d genes | %d samples differ", length(genes), sum(!same)))
    GN <- P12_CFG$genes
    MT13 <- GN$mt_coding
    i_mt <- match(MT13, genes)
    i_goi <- match("MT-CO3", genes)
    i_dfg <- stats::setNames(match(GN$cdfg, genes), GN$cdfg)
    p12_chk(!anyNA(i_mt) && !is.na(i_goi) && !anyNA(i_dfg),
            "13 MT genes, MT-CO3 and its eight partners present", "checked")
    if (anyNA(i_mt) || is.na(i_goi) || anyNA(i_dfg)) p12_abort("MT genes missing")
    LIN <- c(P12_CFG$panels$other_lineage, list(microglia = P12_CFG$panels$microglia))
    LIN_I <- lapply(LIN, function(g) stats::na.omit(match(g, genes)))
    for (k in names(LIN)) {
        ms <- setdiff(LIN[[k]], genes)
        p12_chk(!length(ms), sprintf("lineage panel %s genes present", k),
                if (length(ms)) paste("missing:", paste(ms, collapse = ", ")) else "all",
                warn_only = TRUE)
    }
    p12_say("")

    # == C. per sample =======================================================
    p12_say(sprintf("SECTION C -- cell calling (emptyDrops, lower %d, FDR %.2f, %d iterations)",
                    CC$lower, CC$fdr, CC$niters))
    seed0 <- p12_seed(STEP)
    jobs <- lapply(seq_len(nrow(S)), function(i) list(
        i = i, sample = S$title[i], gsm = S$gsm[i], donor = S$donor[i],
        region = S$region_code[i], group = S$group[i], file_stem = S$file_stem[i],
        path = file.path(dep, paste0(S$file_stem[i], G$matrix_suffix)),
        n_genes = length(genes), i_mt = i_mt, i_goi = i_goi, i_dfg = i_dfg,
        lin_i = LIN_I, amb = AMB, lower = CC$lower, niters = CC$niters, fdr = CC$fdr,
        seed = seed0 + i, rng_kind = P12_CFG$inference$rng_kind))

    # All logical cores but one (or P12_WORKERS). Seeds are per sample, so the
    # result does not depend on the number of workers.
    n_work <- if (exists("P12_WORKERS", envir = globalenv()))
        as.integer(get("P12_WORKERS", envir = globalenv())) else
        max(1L, parallel::detectCores(logical = TRUE) - 1L)
    n_work <- max(1L, min(n_work, length(jobs)))
    p12_note("workers", sprintf("%d of %d logical cores (PSOCK); each holds one raw matrix at a time",
                                n_work, parallel::detectCores(logical = TRUE)))
    t_c <- Sys.time()
    if (n_work > 1L) {
        cl <- parallel::makeCluster(n_work, type = "PSOCK")
        on.exit(parallel::stopCluster(cl), add = TRUE)
        res <- parallel::parLapplyLB(cl, jobs, p12_gerrits_sample)
    } else {
        res <- lapply(jobs, p12_gerrits_sample)
    }
    p12_note("cell calling time", sprintf("%.1f min",
                                          as.numeric(difftime(Sys.time(), t_c, units = "mins"))))

    srows <- list(); crows <- list()
    for (r in res) {
        if (!isTRUE(r$ok)) {
            p12_chk(FALSE, sprintf("%s processed", r$sample), r$msg)
            next
        }
        if (isTRUE(r$transposed))
            p12_note(sprintf("%s matrix stored barcodes x genes", r$sample), "transposed")
        p12_chk(r$srow$n_called > 0L, sprintf("%s emptyDrops called cells", r$sample),
                sprintf("%d called", r$srow$n_called))
        if (r$srow$n_called == 0L) next
        srows[[length(srows) + 1L]] <- r$srow
        crows[[length(crows) + 1L]] <- r$cells
        x <- r$srow
        p12_say(sprintf("  %-18s barcodes %8d | called %6d | median UMIs %6.0f | median %%MT %.3f | ambient %%MT %.3f | MT-CO3+ %d",
                        x$sample, x$n_barcodes, x$n_called, x$median_umis,
                        100 * x$median_mito_frac, 100 * x$ambient_mito_frac, x$n_goi_positive))
    }
    SS <- do.call(rbind, srows); CE <- do.call(rbind, crows)
    if (is.null(SS)) p12_abort("no sample produced called cells")
    p12_say("")

    # == D. totals and write =================================================
    p12_say("SECTION D -- totals")
    p12_chk(nrow(SS) == nrow(S), "every sample produced a cell set",
            sprintf("%d of %d", nrow(SS), nrow(S)))
    low <- SS$sample[SS$n_called < G$min_cells_per_sample]
    p12_chk(!length(low), sprintf("samples with at least %d called nuclei", G$min_cells_per_sample),
            if (length(low)) paste("below:", paste(low, collapse = ", ")) else "all",
            warn_only = TRUE)
    p12_say(sprintf("  %d nuclei called | median per sample %.0f | lineage share: %s",
                    nrow(CE), stats::median(SS$n_called),
                    paste(sprintf("%s %.2f", names(table(CE$lineage)),
                                  table(CE$lineage) / nrow(CE)), collapse = " | ")))
    SS$config_version <- CE$config_version <- P12_CFG$config_version
    SS$run_id <- CE$run_id <- P12$run_id
    if (P12$n[["FAIL"]] > 0L) p12_abort("checks failed; tables not written")
    p12_output(F_SAMPLES)
    p12_output(F_CELLS)
    utils::write.csv(SS, F_SAMPLES, row.names = FALSE)
    con <- gzfile(F_CELLS, "wb")
    utils::write.csv(CE, con, row.names = FALSE)
    close(con)
    back <- utils::read.csv(F_CELLS)
    p12_chk(nrow(back) == nrow(CE), "tables written and read back",
            sprintf("samples %d | nuclei %d", nrow(SS), nrow(back)))
}

p12_run(main)
