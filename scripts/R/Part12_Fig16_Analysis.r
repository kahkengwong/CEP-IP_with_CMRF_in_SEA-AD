# =============================================================================
# Part12_Fig16_Analysis.r
# SEA-AD CEP-IP/CMRF, Part 12, Figure 16 analysis (measurement checks)
#   Sections 13 to 17: donor share of variance and split-half (16A, 16B), label
#   agreement between arms (16C), within-donor co-expression (16D), offset
#   share (16E) and library fields (16F, 16G). Stops at the first section that
#   does not pass. The figures are drawn by Part12_Fig16_Plots.r.
#
# Usage (from the SEA-AD folder, or set SEAAD_ROOT), after steps 02 to 05
#   source("Part12_rerun/Part12_Fig16_Analysis.r")
#
# Outputs (Part12_rerun/08_Fig16/<section>/): stats, plot data and logs.
# =============================================================================

source(file.path(Sys.getenv("SEAAD_ROOT", getwd()), "Part12_rerun", "Part12_03_Config.r"))

SCRIPT <- "Part12_Fig16_Analysis.r"
FIG_DIR <- "08_Fig16"   # parent of every Figure 16 section folder

# == Helpers shared by the sections ===========================================
# Spearman from ranks; NA when either side has no spread.
p12_rs <- function(x, y) {
    if (length(x) < 3L) return(NA_real_)
    rx <- rank(x); ry <- rank(y)
    sx <- stats::sd(rx); sy <- stats::sd(ry)
    if (!is.finite(sx) || !is.finite(sy) || sx == 0 || sy == 0) return(NA_real_)
    stats::cor(rx, ry)
}

# Partial Spearman: Pearson correlation of the ranks of x and y after the
# ranks of the columns of Z are regressed out.
p12_partial_rs <- function(x, y, Z) {
    R <- cbind(1, apply(as.matrix(Z), 2L, rank))
    rx <- stats::lm.fit(R, rank(x))$residuals
    ry <- stats::lm.fit(R, rank(y))$residuals
    if (stats::sd(rx) == 0 || stats::sd(ry) == 0) return(NA_real_)
    stats::cor(rx, ry)
}

# Median of each row of a matrix.
p12_row_median <- function(M) {
    M <- t(apply(M, 1L, sort.int, method = "quick", na.last = TRUE))
    n <- ncol(M)
    if (n %% 2L) M[, (n + 1L) / 2L] else (M[, n / 2L] + M[, n / 2L + 1L]) / 2
}

# Adds the three bookkeeping columns and orders the stats columns.
p12_finalise <- function(X) {
    X$run_id <- rep(P12$run_id, nrow(X))
    X$config_version <- rep(P12_CFG$config_version, nrow(X))
    X$labels_version <- rep(P12_CFG$labels$labels_version, nrow(X))
    rownames(X) <- NULL
    X[, P12_CFG$stats_columns]
}

# One stats row.
p12_row <- function(panel, arm, level, quantity, measure, n_units, n_donors, est,
                    lo = NA_real_, hi = NA_real_, p = NA_real_, family) {
    data.frame(panel = panel, arm = arm, level = level, quantity = quantity, measure = measure,
               n_units = as.integer(n_units), n_donors = as.integer(n_donors),
               estimate = as.numeric(est), lo = as.numeric(lo), hi = as.numeric(hi),
               p = as.numeric(p), fdr = NA_real_, family = family)
}

p12_bh_by_family <- function(X) {
    for (f in unique(X$family)) { i <- X$family == f; X$fdr[i] <- p12_bh(X$p[i]) }
    X
}

# 16B bootstrap on one chunk of resamples: median over splits of rs between
# the half means (A, Bm: splits x donors) of the drawn donors.
p12_16b_chunk <- function(idx, A, Bm) {
    out <- numeric(nrow(idx))
    for (b in seq_len(nrow(idx))) {
        a <- A[, idx[b, ], drop = FALSE]; m <- Bm[, idx[b, ], drop = FALSE]
        ra <- apply(a, 1L, rank); rm <- apply(m, 1L, rank)
        ra <- sweep(ra, 2L, colMeans(ra)); rm <- sweep(rm, 2L, colMeans(rm))
        den <- sqrt(colSums(ra^2) * colSums(rm^2))
        r <- ifelse(den > 0, colSums(ra * rm) / den, NA_real_)
        out[b] <- stats::median(r, na.rm = TRUE)
    }
    out
}

# 16F/16G bootstrap on one chunk of resamples: each row is repeated as often
# as its donor was drawn. Returns a resamples x tasks matrix of rs.
p12_16fg_chunk <- function(idx, tasks, nd) {
    out <- matrix(NA_real_, nrow(idx), length(tasks))
    for (b in seq_len(nrow(idx))) {
        w <- tabulate(idx[b, ], nbins = nd)
        for (k in seq_along(tasks)) {
            t <- tasks[[k]]
            r <- rep.int(seq_along(t$d), w[t$d])
            if (length(unique(t$d[r])) < t$min_d) next
            rx <- rank(t$x[r]); ry <- rank(t$y[r])
            sx <- stats::sd(rx); sy <- stats::sd(ry)
            if (sx > 0 && sy > 0) out[b, k] <- stats::cor(rx, ry)
        }
    }
    out
}

# Runs a chunked bootstrap on the cluster when there is one.
p12_boot_run <- function(cl, idx, FUN, ...) {
    parts <- lapply(p12_chunks(nrow(idx), max(1L, if (is.null(cl)) 1L else length(cl))),
                    function(r) idx[r, , drop = FALSE])
    res <- if (is.null(cl)) lapply(parts, FUN, ...) else parallel::parLapply(cl, parts, FUN, ...)
    if (is.matrix(res[[1]])) do.call(rbind, res) else unlist(res, use.names = FALSE)
}

p12_cluster <- function(n_jobs) {
    n <- p12_n_workers(n_jobs)
    p12_note("workers", sprintf("%d of %d logical cores (PSOCK)", n,
                                parallel::detectCores(logical = TRUE)))
    if (n > 1L) parallel::makeCluster(n, type = "PSOCK") else NULL
}

# =============================================================================
# SECTION 13: 16A and 16B
# 16A: donor share of variance per arm and quantity, against within-region
# donor shuffles. 16B: split-half rs over 200 random splits, donor bootstrap.
# =============================================================================
p12_section_13 <- function() {
    STEP <- "13"
    D <- p12_step_dirs(STEP, "16AB_DonorStructure", parent = FIG_DIR)
    p12_clean(D$out)
    p12_start(STEP, SCRIPT, D)
    F_A <- file.path(D$out, "Stats_16A.csv")
    F_B <- file.path(D$out, "Stats_16B.csv")
    QA <- c("module_post_pct", "module_resid_tq4", "mito_log1p", "contamination", "ctrl_pct_null")
    QB <- c("module_post_pct", "module_resid_tq4", "mito_log1p")
    N_SPLITS <- 200L

    main <- function() {
        R_PERM <- as.integer(P12_CFG$inference$permutation_R)
        B <- as.integer(P12_CFG$inference$bootstrap_B)
        MIND <- as.integer(P12_CFG$floors$donors_per_test)
        p12_reading("16A")
        p12_reading("16B")
        p12_say("SECTION A -- inputs")
        U <- utils::read.csv(p12_upstream_file("05_Units", "05_units.csv"), check.names = FALSE,
                             colClasses = c(donor_key = "character", region = "character"))
        p12_chk(all(U$labels_version == P12_CFG$labels$labels_version),
                "unit table carries this labels version", sprintf("%d units", nrow(U)))
        miss <- setdiff(QA, names(U))
        if (length(miss)) p12_abort(paste("unit table lacks:", paste(miss, collapse = ", ")))
        p12_say("")

        # == 16A ==============================================================
        p12_say(sprintf("SECTION B -- 16A, donor share against %d within-region shuffles", R_PERM))
        set.seed(p12_seed(STEP))
        rows_a <- list()
        for (a in P12_ARMS) {
            for (q in QA) {
                u <- U[U$arm == a & is.finite(U[[q]]), ]
                y <- u[[q]]; don <- as.integer(factor(u$donor_key)); nd <- max(c(0L, don))
                n_d <- tabulate(don, nd)
                share <- function(dd) {
                    s <- rowsum(y, dd, reorder = TRUE)[, 1]
                    ssb <- sum(s^2 / n_d) - length(y) * mean(y)^2
                    ssb / sum((y - mean(y))^2)
                }
                est <- p <- lo <- hi <- NA_real_
                if (nd >= MIND && isTRUE(stats::var(y) > 0)) {
                    est <- share(don)
                    reg <- split(seq_along(y), u$region)
                    nul <- numeric(R_PERM)
                    for (r in seq_len(R_PERM)) {
                        dd <- don
                        for (g in reg) if (length(g) > 1L) dd[g] <- don[g][sample.int(length(g))]
                        nul[r] <- share(dd)
                    }
                    p <- (sum(nul >= est - 1e-12) + 1) / (R_PERM + 1)
                    lo <- stats::median(nul); hi <- stats::quantile(nul, 0.95, names = FALSE)
                }
                rows_a[[length(rows_a) + 1L]] <- p12_row("16A", a, "units", q,
                                                         "donor share of variance", nrow(u), nd,
                                                         est, lo, hi, p, paste("16A", a, sep = "|"))
                p12_say(sprintf("  %-6s %-17s units %4d | donors %3d | donor share %5.1f%% | null median %5.1f%%, 95th %5.1f%% | p %s",
                                a, q, nrow(u), nd, 100 * est, 100 * lo, 100 * hi,
                                p12_fmt(p, 1 / (R_PERM + 1))))
            }
        }
        A <- p12_bh_by_family(do.call(rbind, rows_a))
        p12_say("")

        # == 16B ==============================================================
        p12_say(sprintf("SECTION C -- 16B, split-half rs over %d splits (B = %d)", N_SPLITS, B))
        cl <- p12_cluster(B)
        if (!is.null(cl)) on.exit(parallel::stopCluster(cl), add = TRUE)
        # Splits use their own stream (step seed + 1), carried from task to task,
        # because p12_boot_index() re-seeds on every call.
        set.seed(p12_seed(STEP) + 1L)
        split_state <- get(".Random.seed", envir = globalenv())
        rows_b <- list()
        for (a in P12_ARMS) {
            for (q in QB) {
                u <- U[U$arm == a & is.finite(U[[q]]), ]
                k <- table(u$donor_key)
                dons <- names(k)[k >= 2L]
                u <- u[u$donor_key %in% dons, ]
                est <- lo <- hi <- p <- NA_real_
                if (length(dons) >= MIND) {
                    byd <- split(u[[q]], u$donor_key)[dons]
                    A1 <- A2 <- matrix(NA_real_, N_SPLITS, length(dons))
                    assign(".Random.seed", split_state, envir = globalenv())
                    for (s in seq_len(N_SPLITS)) for (j in seq_along(dons)) {
                        v <- byd[[j]][sample.int(length(byd[[j]]))]
                        h <- length(v) %/% 2L
                        A1[s, j] <- mean(v[seq_len(h)]); A2[s, j] <- mean(v[-seq_len(h)])
                    }
                    split_state <- get(".Random.seed", envir = globalenv())
                    est <- stats::median(vapply(seq_len(N_SPLITS), function(s)
                        p12_rs(A1[s, ], A2[s, ]), 0), na.rm = TRUE)
                    idx <- p12_boot_index(dons, STEP, B)
                    bs <- p12_boot_run(cl, idx, p12_16b_chunk, A = A1, Bm = A2)
                    ci <- p12_boot_ci(bs); lo <- ci[["lo"]]; hi <- ci[["hi"]]; p <- p12_boot_p(bs)
                }
                rows_b[[length(rows_b) + 1L]] <- p12_row("16B", a, "donors", q,
                                                         "split-half rs, median over splits",
                                                         nrow(u), length(dons), est, lo, hi, p,
                                                         paste("16B", a, sep = "|"))
                p12_say(sprintf("  %-6s %-17s donors %3d | split-half rs %+.3f [%+.3f, %+.3f] | p %s",
                                a, q, length(dons), est, lo, hi, p12_fmt(p, p12_boot_floor(B))))
            }
        }
        Bt <- p12_bh_by_family(do.call(rbind, rows_b))
        p12_say("")

        p12_say("SECTION D -- write")
        if (P12$n[["FAIL"]] > 0L) p12_abort("checks failed; stats not written")
        p12_write_stats(p12_finalise(A), F_A)
        p12_write_stats(p12_finalise(Bt), F_B)
        p12_chk(TRUE, "stats files written", sprintf("16A %d rows | 16B %d rows", nrow(A), nrow(Bt)))
    }

    p12_run(main)
}

# =============================================================================
# SECTION 14: 16C
# Cohen's kappa for Module membership on cells shared with the PVM_3 arm
# (M4 positive control, HKG8 negative control); donor-bootstrap interval.
# =============================================================================
p12_section_14 <- function() {
    STEP <- "14"
    D <- p12_step_dirs(STEP, "16C_LabelAgreement", parent = FIG_DIR)
    p12_clean(D$out)
    p12_start(STEP, SCRIPT, D)
    F_C <- file.path(D$out, "Stats_16C.csv")
    REF <- "PVM_3"; TESTS <- c("M4", "HKG8")

    kappa4 <- function(n) {
        N <- sum(n)
        if (N == 0) return(NA_real_)
        po <- (n[1] + n[4]) / N
        pe <- ((n[1] + n[2]) * (n[1] + n[3]) + (n[4] + n[3]) * (n[4] + n[2])) / N^2
        if (pe >= 1) NA_real_ else (po - pe) / (1 - pe)
    }

    main <- function() {
        B <- as.integer(P12_CFG$inference$bootstrap_B)
        MINC <- as.integer(P12_CFG$floors$cells_per_group_per_donor)
        MIND <- as.integer(P12_CFG$floors$donors_per_test)
        p12_reading("16C")
        p12_say("SECTION A -- inputs")
        L <- utils::read.csv(p12_upstream_file("04_Labels", "04_cell_labels.csv.gz"),
                             check.names = FALSE,
                             colClasses = c(barcode = "character", donor_key = "character",
                                            arm = "character", dataset = "character",
                                            region = "character"))
        p12_chk(all(L$labels_version == P12_CFG$labels$labels_version),
                "labels carry this labels version", "checked")
        ref <- L[L$arm == REF, c("barcode", "is_module")]
        names(ref)[2] <- "ref_module"
        p12_say("")

        p12_say(sprintf("SECTION B -- kappa against %s (regions with at least %d shared cells)", REF, MINC))
        rows <- list()
        for (a in TESTS) {
            t <- merge(L[L$arm == a, c("barcode", "is_module", "donor_key", "dataset", "region")],
                       ref, by = "barcode")
            nreg <- table(paste(t$dataset, t$region))
            t <- t[paste(t$dataset, t$region) %in% names(nreg)[nreg >= MINC], ]
            t$cell <- ifelse(t$is_module & t$ref_module, 1L, ifelse(t$is_module, 2L,
                             ifelse(t$ref_module, 3L, 4L)))
            for (lv in c("pooled", P12_DATASETS)) {
                s <- if (lv == "pooled") t else t[t$dataset == lv, ]
                dons <- sort(unique(s$donor_key))
                M <- matrix(0, length(dons), 4)
                if (nrow(s)) {
                    tb <- table(factor(s$donor_key, levels = dons), factor(s$cell, levels = 1:4))
                    M <- matrix(as.numeric(tb), nrow = length(dons))
                }
                n4 <- colSums(M)
                k <- kappa4(n4)
                lo <- hi <- NA_real_
                if (length(dons) >= MIND) {
                    idx <- p12_boot_index(dons, STEP, B)
                    W <- t(apply(idx, 1L, tabulate, nbins = length(dons)))
                    S <- W %*% M
                    bs <- apply(S, 1L, kappa4)
                    ci <- p12_boot_ci(bs); lo <- ci[["lo"]]; hi <- ci[["hi"]]
                }
                rows[[length(rows) + 1L]] <- p12_row("16C", a, lv, "kappa",
                                                     sprintf("Module membership vs %s", REF),
                                                     sum(n4), length(dons), k, lo, hi, NA,
                                                     "16C")
                for (j in 1:4) {
                    nm <- c("share_both", "share_test_only", "share_ref_only", "share_neither")[j]
                    rows[[length(rows) + 1L]] <- p12_row("16C", a, lv, nm,
                                                         sprintf("share of shared cells vs %s", REF),
                                                         sum(n4), length(dons),
                                                         if (sum(n4)) n4[j] / sum(n4) else NA,
                                                         family = "16C|shares")
                }
                p12_say(sprintf("  %-5s %-8s shared cells %6d | donors %3d | kappa %.3f [%.3f, %.3f] | both %.1f%% test only %.1f%% ref only %.1f%% neither %.1f%%",
                                a, lv, sum(n4), length(dons), k, lo, hi,
                                100 * n4[1] / max(1, sum(n4)), 100 * n4[2] / max(1, sum(n4)),
                                100 * n4[3] / max(1, sum(n4)), 100 * n4[4] / max(1, sum(n4))))
            }
        }
        X <- do.call(rbind, rows)
        kp <- X[X$quantity == "kappa" & X$level == "pooled", ]
        k_m4 <- kp$estimate[kp$arm == "M4"]; k_hk <- kp$estimate[kp$arm == "HKG8"]
        p12_note("pre-set thresholds, pooled kappa",
                 sprintf("M4 %.3f (> 0.60: %s) | HKG8 %.3f (|kappa| < 0.20: %s)", k_m4,
                         isTRUE(k_m4 > 0.60), k_hk, isTRUE(abs(k_hk) < 0.20)))
        p12_say("")
        p12_say("SECTION C -- write")
        if (P12$n[["FAIL"]] > 0L) p12_abort("checks failed; stats not written")
        p12_write_stats(p12_finalise(X), F_C)
        p12_chk(TRUE, "stats file written", sprintf("%d rows", nrow(X)))
    }

    p12_run(main)
}

# =============================================================================
# SECTION 15: 16D
# Within-donor rs of MT-CO3 with the cDFG composite on the PVM_3 arm, raw
# and depth-partialled; control pair ACTB with cHKG8.
# =============================================================================
p12_section_15 <- function() {
    STEP <- "15"
    D <- p12_step_dirs(STEP, "16D_WithinDonorRho", parent = FIG_DIR)
    p12_clean(D$out)
    p12_start(STEP, SCRIPT, D)
    F_D <- file.path(D$out, "Stats_16D.csv")
    F_DP <- file.path(D$out, "Plot_16D_donors.csv")

    main <- function() {
        B <- as.integer(P12_CFG$inference$bootstrap_B)
        MINC <- as.integer(P12_CFG$floors$cells_per_group_per_donor)
        MIND <- as.integer(P12_CFG$floors$donors_per_test)
        p12_reading("16D")
        p12_say("SECTION A -- inputs")
        L <- utils::read.csv(p12_upstream_file("04_Labels", "04_cell_labels.csv.gz"),
                             check.names = FALSE,
                             colClasses = c(barcode = "character", donor_key = "character",
                                            arm = "character"))
        p12_chk(all(L$labels_version == P12_CFG$labels$labels_version),
                "labels carry this labels version", "checked")
        X <- methods::as(Matrix::readMM(p12_upstream_file("02_Extract", "02_counts.mtx.gz")),
                         "CsparseMatrix")
        genes <- readLines(p12_upstream_file("02_Extract", "02_counts_genes.tsv"))
        cells <- readLines(p12_upstream_file("02_Extract", "02_counts_cells.tsv"))
        m <- L[L$arm == "PVM_3", ]
        ci <- match(m$barcode, cells)
        p12_chk(all(!is.na(ci)), "every PVM_3 arm cell is in the step 02 count matrix",
                sprintf("%d of %d", sum(!is.na(ci)), length(ci)))
        PAIRS <- list(
            list(name = "MT-CO3 vs cDFG", goi = "MT-CO3", genes = P12_CFG$genes$cdfg),
            list(name = "ACTB vs cHKG8", goi = P12_CFG$arms$HKG8$goi,
                 genes = P12_CFG$genes[[P12_CFG$arms$HKG8$partners]]))
        need <- unique(unlist(lapply(PAIRS, function(p) c(p$goi, p$genes))))
        gi <- match(need, genes)
        if (anyNA(gi)) p12_abort(paste("genes missing from step 02:",
                                       paste(need[is.na(gi)], collapse = ", ")))
        Y <- as.matrix(X[ci, gi, drop = FALSE]); colnames(Y) <- need
        tot <- m$total_umis
        E <- log1p(Y / tot * 1e4)
        p12_say(sprintf("  population: %d cells from %d donors", nrow(m),
                        length(unique(m$donor_key))))
        p12_say("")

        p12_say(sprintf("SECTION B -- per-donor rs (donors with at least %d cells, B = %d)", MINC, B))
        n_d <- table(m$donor_key)
        dons <- names(n_d)[n_d >= MINC]
        rows <- list(); plot_rows <- list()
        for (pr in PAIRS) {
            g <- E[, pr$goi]; cmp <- rowMeans(E[, pr$genes, drop = FALSE])
            per <- t(vapply(dons, function(d) {
                i <- which(m$donor_key == d)
                c(length(i), p12_rs(g[i], cmp[i]), p12_partial_rs(g[i], cmp[i], log(tot[i])))
            }, numeric(3)))
            plot_rows[[pr$name]] <- data.frame(donor_key = dons, pair = pr$name,
                                               n_cells = as.integer(per[, 1]), rs_raw = per[, 2],
                                               rs_partial = per[, 3])
            pooled <- p12_rs(g, cmp)
            dm_g <- tapply(g, m$donor_key, mean); dm_c <- tapply(cmp, m$donor_key, mean)
            between <- p12_rs(dm_g, dm_c[names(dm_g)])
            for (k in 2:3) {
                v <- per[, k]; ok <- is.finite(v)
                est <- stats::median(v[ok]); lo <- hi <- NA_real_
                if (sum(ok) >= MIND) {
                    idx <- p12_boot_index(dons[ok], STEP, B)
                    bs <- p12_row_median(matrix(v[ok][idx], nrow = B))
                    cc <- p12_boot_ci(bs); lo <- cc[["lo"]]; hi <- cc[["hi"]]
                }
                qn <- if (k == 2) "within_median_raw" else "within_median_partial"
                rows[[length(rows) + 1L]] <- p12_row("16D", "PVM_3", "within donor", qn, pr$name,
                                                     sum(per[ok, 1]), sum(ok), est, lo, hi,
                                                     family = "16D")
                if (k == 2) rows[[length(rows) + 1L]] <- p12_row(
                    "16D", "PVM_3", "within donor", "within_iqr_raw", pr$name, sum(per[ok, 1]),
                    sum(ok), est, stats::quantile(v[ok], 0.25, names = FALSE),
                    stats::quantile(v[ok], 0.75, names = FALSE), family = "16D|iqr")
            }
            rows[[length(rows) + 1L]] <- p12_row("16D", "PVM_3", "pooled", "pooled_raw", pr$name,
                                                 nrow(m), length(unique(m$donor_key)), pooled,
                                                 family = "16D")
            rows[[length(rows) + 1L]] <- p12_row("16D", "PVM_3", "between donors", "between_raw",
                                                 pr$name, length(dm_g), length(dm_g), between,
                                                 family = "16D")
            p12_say(sprintf("  %-15s donors %3d | within median raw %.4f | partial %.4f | pooled %.4f | between %.4f",
                            pr$name, sum(is.finite(per[, 2])), stats::median(per[, 2], na.rm = TRUE),
                            stats::median(per[, 3], na.rm = TRUE), pooled, between))
        }
        Xs <- do.call(rbind, rows)
        wm <- Xs$estimate[Xs$quantity == "within_median_raw" & Xs$measure == "MT-CO3 vs cDFG"]
        p12_note("pre-set marks on the raw within-donor median",
                 sprintf("%.4f -> %s", wm, if (!is.finite(wm)) "untestable" else if (wm >= 0.75)
                     "at or above 0.75" else if (wm >= 0.40) "0.40 to 0.75" else "below 0.40"))
        p12_say("")
        p12_say("SECTION C -- write")
        if (P12$n[["FAIL"]] > 0L) p12_abort("checks failed; stats not written")
        p12_write_stats(p12_finalise(Xs), F_D)
        PD <- do.call(rbind, plot_rows); rownames(PD) <- NULL
        p12_output(F_DP)
        utils::write.csv(PD, F_DP, row.names = FALSE)
        p12_chk(TRUE, "stats and plot files written", sprintf("%d stats rows | %d donor rows",
                                                              nrow(Xs), nrow(PD)))
    }

    p12_run(main)
}

# =============================================================================
# SECTION 16: 16E
# Offset share a / (a + b x median depth) from per-donor OLS of MT UMIs on
# total UMIs; post-IP minus pre-IP, paired donor bootstrap, 10 pp margin.
# =============================================================================
p12_section_16 <- function() {
    STEP <- "16"
    D <- p12_step_dirs(STEP, "16E_Offset", parent = FIG_DIR)
    p12_clean(D$out)
    p12_start(STEP, SCRIPT, D)
    F_E <- file.path(D$out, "Stats_16E.csv")
    MARGIN <- 10

    main <- function() {
        B <- as.integer(P12_CFG$inference$bootstrap_B)
        MINC <- as.integer(P12_CFG$floors$cells_per_group_per_donor)
        MIND <- as.integer(P12_CFG$floors$donors_per_test)
        EQ <- as.numeric(P12_CFG$inference$equivalence_ci_level)
        p12_reading("16E")
        p12_say("SECTION A -- inputs")
        L <- utils::read.csv(p12_upstream_file("04_Labels", "04_cell_labels.csv.gz"),
                             check.names = FALSE,
                             colClasses = c(barcode = "character", donor_key = "character",
                                            arm = "character", group = "character"))
        p12_chk(all(L$labels_version == P12_CFG$labels$labels_version),
                "labels carry this labels version", "checked")
        p12_chk(all(is.finite(L$mt_umis) & is.finite(L$total_umis)),
                "MT UMIs and total UMIs on every labelled cell", "checked")
        # Every Module cell must be in one of the two branches.
        gm <- unique(L$group[L$is_module])
        p12_chk(all(gm %in% c("Module post-IP", "Module pre-IP")),
                "every Module cell is Module post-IP or Module pre-IP", paste(sort(gm), collapse = ", "))
        p12_say("")

        p12_say(sprintf("SECTION B -- offset share per branch (%d cells per donor and branch, %d donors, B = %d)",
                        MINC, MIND, B))
        fit <- function(x, y) {
            vx <- stats::var(x)
            if (length(x) < 2L || !is.finite(vx) || vx <= 0) return(c(NA_real_, NA_real_))
            b <- stats::cov(x, y) / vx
            c(mean(y) - b * mean(x), b)
        }
        share <- function(a, b, d) { den <- a + b * d; if (is.finite(den) && den != 0) a / den else NA_real_ }
        rows <- list()
        for (arm in P12_ARMS) {
            m <- L[L$arm == arm & L$is_module, ]
            dref <- stats::median(m$total_umis)
            m$br <- ifelse(m$group == "Module post-IP", "post", "pre")
            FT <- list()
            for (br in c("post", "pre")) {
                s <- m[m$br == br, ]
                sp <- split(seq_len(nrow(s)), s$donor_key)
                sp <- sp[lengths(sp) >= MINC]
                ab <- t(vapply(sp, function(i) fit(s$total_umis[i], s$mt_umis[i]), numeric(2)))
                ok <- is.finite(ab[, 1]) & is.finite(ab[, 2])
                FT[[br]] <- data.frame(donor = names(sp)[ok], a = ab[ok, 1], b = ab[ok, 2])
            }
            est_one <- function(fp, fq) {
                if (nrow(fp) < MIND || nrow(fq) < MIND) return(rep(NA_real_, 7))
                ap <- stats::median(fp$a); bp <- stats::median(fp$b)
                aq <- stats::median(fq$a); bq <- stats::median(fq$b)
                sp <- share(ap, bp, dref); sq <- share(aq, bq, dref)
                c(100 * (sp - sq), sp, sq, ap, aq, bp, bq)
            }
            e <- est_one(FT$post, FT$pre)
            uni <- sort(unique(c(FT$post$donor, FT$pre$donor)))
            lo95 <- hi95 <- lo90 <- hi90 <- p <- NA_real_
            if (is.finite(e[1])) {
                idx <- p12_boot_index(uni, STEP, B)
                ip <- match(FT$post$donor, uni); iq <- match(FT$pre$donor, uni)
                bs <- vapply(seq_len(B), function(b) {
                    w <- tabulate(idx[b, ], nbins = length(uni))
                    rp <- rep.int(seq_along(ip), w[ip]); rq <- rep.int(seq_along(iq), w[iq])
                    est_one(FT$post[rp, ], FT$pre[rq, ])[1]
                }, 0)
                c95 <- p12_boot_ci(bs, 0.95); c90 <- p12_boot_ci(bs, EQ)
                lo95 <- c95[["lo"]]; hi95 <- c95[["hi"]]; lo90 <- c90[["lo"]]; hi90 <- c90[["hi"]]
                p <- p12_boot_p(bs)
            }
            # The verdict is set after BH over the arms, below.
            verdict <- if (!is.finite(e[1])) "untestable" else "pending"
            np <- nrow(FT$post); nq <- nrow(FT$pre)
            rows[[length(rows) + 1L]] <- p12_row("16E", arm, "branches", "d_offset_share_pp",
                                                 "post-IP minus pre-IP", NA, length(uni), e[1],
                                                 lo95, hi95, p, "16E")
            rows[[length(rows) + 1L]] <- p12_row("16E", arm, "branches", "d_offset_share_pp_90",
                                                 verdict, NA, length(uni), e[1], lo90, hi90,
                                                 family = "16E|90")
            rows[[length(rows) + 1L]] <- p12_row("16E", arm, "post-IP", "offset_share_pct",
                                                 "offset share", NA, np, 100 * e[2], family = "16E|branch")
            rows[[length(rows) + 1L]] <- p12_row("16E", arm, "pre-IP", "offset_share_pct",
                                                 "offset share", NA, nq, 100 * e[3], family = "16E|branch")
            rows[[length(rows) + 1L]] <- p12_row("16E", arm, "post-IP", "intercept",
                                                 "median intercept, MT UMIs", NA, np, e[4],
                                                 family = "16E|branch")
            rows[[length(rows) + 1L]] <- p12_row("16E", arm, "pre-IP", "intercept",
                                                 "median intercept, MT UMIs", NA, nq, e[5],
                                                 family = "16E|branch")
            rows[[length(rows) + 1L]] <- p12_row("16E", arm, "all Module cells", "reference_depth",
                                                 "median total UMIs", nrow(m), NA, dref,
                                                 family = "16E|branch")
        }
        X <- p12_bh_by_family(do.call(rbind, rows))
        # Verdicts: non-zero is FDR < 0.05 (BH over the arms); the margin uses the
        # 90% interval against +/-10 pp.
        for (arm in P12_ARMS) {
            d <- X[X$arm == arm & X$quantity == "d_offset_share_pp", ]
            i9 <- which(X$arm == arm & X$quantity == "d_offset_share_pp_90")
            if (!nrow(d) || !is.finite(d$estimate)) next
            within <- is.finite(X$lo[i9]) && X$lo[i9] > -MARGIN && X$hi[i9] < MARGIN
            nonzero <- is.finite(d$fdr) && d$fdr < 0.05
            X$measure[i9] <- if (within && !nonzero) "within the margin" else
                if (within) "non-zero but within the margin" else
                if (nonzero) "non-zero and not shown to be within the margin" else
                    "underpowered: neither non-zero nor shown to be within the margin"
        }
        for (arm in P12_ARMS) {
            d <- X[X$arm == arm & X$quantity == "d_offset_share_pp", ]
            d9 <- X[X$arm == arm & X$quantity == "d_offset_share_pp_90", ]
            sp <- X$estimate[X$arm == arm & X$quantity == "offset_share_pct"]
            nd <- X$n_donors[X$arm == arm & X$quantity == "offset_share_pct"]
            p12_say(sprintf("  %-6s donors post %3d, pre %3d | share post %5.1f%%, pre %5.1f%% | difference %+.1f pp, 95%% [%+.1f, %+.1f], 90%% [%+.1f, %+.1f] | p %s | FDR %s | %s",
                            arm, nd[1], nd[2], sp[1], sp[2], d$estimate, d$lo, d$hi, d9$lo, d9$hi,
                            p12_fmt(d$p, p12_boot_floor(B)), p12_fmt(d$fdr, p12_boot_floor(B)),
                            d9$measure))
        }
        p12_say("")
        p12_say("SECTION C -- write")
        if (P12$n[["FAIL"]] > 0L) p12_abort("checks failed; stats not written")
        p12_write_stats(p12_finalise(X), F_E)
        p12_chk(TRUE, "stats file written", sprintf("%d rows", nrow(X)))
    }

    p12_run(main)
}

# =============================================================================
# SECTION 17: 16F and 16G
# rs of each unit quantity with each library field, between and within
# donor; contrasts against %MT, depth and HKG8, one shared donor bootstrap.
# =============================================================================
p12_section_17 <- function() {
    STEP <- "17"
    D <- p12_step_dirs(STEP, "16FG_LibraryFields", parent = FIG_DIR)
    p12_clean(D$out)
    p12_start(STEP, SCRIPT, D)
    F_F <- file.path(D$out, "Stats_16F.csv")
    F_G <- file.path(D$out, "Stats_16G.csv")
    QS <- c(module = "module_post_pct", mito = "mito_log1p", ncells = "log_n_cells",
            depth = "log_median_umi")
    FIELDS <- c(intronic = "16F", antisense = "16G", neun_fans = "16G", cdna_400 = "16G",
                pcr_cycles = "16G", saturation = "16G")
    LEVELS <- c("between donors", "within donor")

    main <- function() {
        B <- as.integer(P12_CFG$inference$bootstrap_B)
        MIND <- as.integer(P12_CFG$floors$donors_per_test)
        MINU <- as.integer(P12_CFG$floors$finite_units_within_donor)
        MAXMISS <- as.numeric(P12_CFG$floors$max_missing_share)
        p12_reading("16F")
        p12_reading("16G")
        p12_say("SECTION A -- inputs")
        U <- utils::read.csv(p12_upstream_file("05_Units", "05_units.csv"), check.names = FALSE,
                             colClasses = c(donor_key = "character"))
        p12_chk(all(U$labels_version == P12_CFG$labels$labels_version),
                "unit table carries this labels version", sprintf("%d units", nrow(U)))
        miss <- setdiff(QS, names(U))
        if (length(miss)) p12_abort(paste("unit table lacks:", paste(miss, collapse = ", ")))
        use <- character(0)
        for (f in names(FIELDS)) {
            col <- paste0("lib_", f)
            share_na <- if (col %in% names(U)) mean(!is.finite(U[[col]])) else 1
            keep <- col %in% names(U) && share_na <= MAXMISS
            p12_note(sprintf("field %s", f),
                     if (!(col %in% names(U))) "absent from the unit table: dropped" else
                         sprintf("%.1f%% of units missing: %s", 100 * share_na,
                                 if (keep) "kept" else "dropped (above the missing-value limit)"))
            if (keep) use <- c(use, f)
        }
        p12_chk("intronic" %in% use, "the intronic rate is available for 16F", paste(use, collapse = ", "))
        dons <- sort(unique(U$donor_key)); ND <- length(dons)
        p12_say("")

        # Tasks: one per arm, field, level and quantity.
        tasks <- list(); meta <- list()
        for (a in P12_ARMS) for (f in use) for (lv in LEVELS) for (q in names(QS)) {
            u <- U[U$arm == a, ]
            x <- u[[QS[[q]]]]; y <- u[[paste0("lib_", f)]]
            if (lv == "between donors") {
                ok <- is.finite(x) & is.finite(y)
                if (!any(ok)) { x <- y <- numeric(0); d <- character(0) } else {
                    mx <- tapply(x[ok], u$donor_key[ok], mean); my <- tapply(y[ok], u$donor_key[ok], mean)
                    d <- names(mx); x <- as.numeric(mx); y <- as.numeric(my[d])
                }
            } else {
                ok <- is.finite(x) & is.finite(y)
                d <- u$donor_key[ok]
                x <- p12_centre_within(x[ok], d, MINU); y <- p12_centre_within(y[ok], d, MINU)
                k <- is.finite(x) & is.finite(y); x <- x[k]; y <- y[k]; d <- d[k]
            }
            tasks[[length(tasks) + 1L]] <- list(x = x, y = y, d = match(d, dons), min_d = MIND)
            meta[[length(meta) + 1L]] <- data.frame(arm = a, field = f, level = lv, q = q,
                                                    n_units = length(x),
                                                    n_donors = length(unique(d)))
        }
        MT <- do.call(rbind, meta)
        MT$obs <- vapply(tasks, function(t) if (length(unique(t$d)) >= MIND) p12_rs(t$x, t$y) else
            NA_real_, 0)

        p12_say(sprintf("SECTION B -- donor bootstrap over %d donors, %d tasks (B = %d)", ND,
                        length(tasks), B))
        cl <- p12_cluster(B)
        if (!is.null(cl)) on.exit(parallel::stopCluster(cl), add = TRUE)
        idx <- p12_boot_index(dons, STEP, B)
        t0 <- Sys.time()
        BS <- p12_boot_run(cl, idx, p12_16fg_chunk, tasks = tasks, nd = ND)
        p12_note("bootstrap time", sprintf("%.1f min", as.numeric(difftime(Sys.time(), t0,
                                                                          units = "mins"))))
        p12_say("")

        rows <- list()
        key <- function(a, f, lv, q) which(MT$arm == a & MT$field == f & MT$level == lv & MT$q == q)
        for (a in P12_ARMS) for (f in use) for (lv in LEVELS) {
            panel <- FIELDS[[f]]
            fam <- paste(panel, a, lv, "rs", sep = "|")
            famc <- paste(panel, a, lv, "contrast", sep = "|")
            for (q in names(QS)) {
                i <- key(a, f, lv, q)
                ok <- is.finite(MT$obs[i])
                ci <- if (ok) p12_boot_ci(BS[, i]) else c(lo = NA_real_, hi = NA_real_)
                rows[[length(rows) + 1L]] <- p12_row(panel, a, lv, q, f, MT$n_units[i],
                                                     MT$n_donors[i], MT$obs[i], ci[["lo"]],
                                                     ci[["hi"]], if (ok) p12_boot_p(BS[, i]) else NA,
                                                     fam)
            }
            im <- key(a, f, lv, "module")
            con <- list(module_minus_mito = list(est = abs(MT$obs[im]) - abs(MT$obs[key(a, f, lv, "mito")]),
                                                 bs = abs(BS[, im]) - abs(BS[, key(a, f, lv, "mito")])),
                        module_minus_depth = list(est = MT$obs[im] - MT$obs[key(a, f, lv, "depth")],
                                                  bs = BS[, im] - BS[, key(a, f, lv, "depth")]))
            if (a != "HKG8") {
                ih <- key("HKG8", f, lv, "module")
                con$arm_minus_hkg8 <- list(est = MT$obs[im] - MT$obs[ih], bs = BS[, im] - BS[, ih])
            }
            for (cn in names(con)) {
                ok <- is.finite(con[[cn]]$est)
                ci <- if (ok) p12_boot_ci(con[[cn]]$bs) else c(lo = NA_real_, hi = NA_real_)
                rows[[length(rows) + 1L]] <- p12_row(panel, a, lv, cn, f, MT$n_units[im],
                                                     MT$n_donors[im], con[[cn]]$est, ci[["lo"]],
                                                     ci[["hi"]], if (ok) p12_boot_p(con[[cn]]$bs) else NA,
                                                     famc)
            }
        }
        X <- p12_bh_by_family(do.call(rbind, rows))
        for (a in P12_ARMS) for (lv in LEVELS) {
            p12_say(sprintf("  %s, %s", a, lv))
            for (f in use) {
                g <- X[X$arm == a & X$level == lv & X$measure == f, ]
                v <- function(q) g[g$quantity == q, ]
                m <- v("module")
                p12_say(sprintf("    %-11s Module rs %+.3f [%+.3f, %+.3f] FDR %s | %%MT %+.3f | cells %+.3f | depth %+.3f | donors %d",
                                f, m$estimate, m$lo, m$hi, p12_fmt(m$fdr, p12_boot_floor(B)),
                                v("mito")$estimate, v("ncells")$estimate, v("depth")$estimate,
                                m$n_donors))
            }
        }
        p12_say("")
        p12_say("SECTION C -- the four conditions for 16G, within donor")
        excl <- function(r) nrow(r) == 1L && is.finite(r$lo) && (r$lo > 0 || r$hi < 0)
        for (a in setdiff(P12_ARMS, "HKG8")) for (f in setdiff(use, "intronic")) {
            g <- X[X$panel == "16G" & X$arm == a & X$level == "within donor" & X$measure == f, ]
            m <- g[g$quantity == "module", ]
            c1 <- nrow(m) == 1L && is.finite(m$fdr) && m$fdr < 0.05
            c2 <- excl(g[g$quantity == "module_minus_depth", ])
            c3 <- excl(g[g$quantity == "arm_minus_hkg8", ])
            mm <- g[g$quantity == "module_minus_mito", ]
            c4 <- nrow(mm) == 1L && is.finite(mm$hi) && mm$hi >= 0
            testable <- nrow(m) == 1L && is.finite(m$estimate)
            p12_say(sprintf("  %-6s %-11s non-zero %-5s | vs depth %-5s | vs HKG8 %-5s | not matched by %%MT %-5s -> %s",
                            a, f, c1, c2, c3, c4,
                            if (!testable) "untestable (below the donor floor)" else
                                if (c1 && c2 && c3 && c4) "orders Module post-IP%" else
                                    "does not order Module post-IP%"))
        }
        p12_say("")
        p12_say("SECTION D -- write")
        if (P12$n[["FAIL"]] > 0L) p12_abort("checks failed; stats not written")
        p12_write_stats(p12_finalise(X[X$panel == "16F", ]), F_F)
        p12_write_stats(p12_finalise(X[X$panel == "16G", ]), F_G)
        p12_chk(TRUE, "stats files written", sprintf("16F %d rows | 16G %d rows",
                                                     sum(X$panel == "16F"), sum(X$panel == "16G")))
    }

    p12_run(main)
}

# =============================================================================
# RUN
# Sections in order; stop at the first that does not PASS.
# =============================================================================
local({
    # Clear every Figure 16 section folder first, so no output from an earlier
    # run survives.
    for (d in list.dirs(file.path(P12_RERUN, FIG_DIR), recursive = FALSE)) p12_clean(d)
    for (code in c("13", "14", "15", "16", "17")) {
        v <- get(paste0("p12_section_", code))()
        if (!identical(v, "PASS"))
            stop(sprintf("section %s ended with verdict %s; later sections not run. See its log.",
                         code, v), call. = FALSE)
    }
    message("Figure 16 analysis: sections 13 to 17 all PASS. Next: run Part12_Fig16_Plots.r")
})
