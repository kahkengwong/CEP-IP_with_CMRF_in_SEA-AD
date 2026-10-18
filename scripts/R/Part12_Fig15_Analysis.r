# =============================================================================
# Part12_Fig15_Analysis.r
# SEA-AD CEP-IP/CMRF, Part 12, Figure 15 analysis (Figure 15, Supplementary Figure 10)
#   Sections 07-11: contamination covariate and doublet leave-out (15A-B),
#   contamination by cell group and tertile (15C-D), Gerrits ambient pool
#   (15E-G), depth thinning (15H) and own cytoplasm against other lineages
#   (S10). Stops at the first section that does not pass.
#
# Usage (from the SEA-AD folder, or set SEAAD_ROOT), after steps 02, 04, 05 and 06
#   source("Part12_rerun/Part12_Fig15_Analysis.r")
#   then Part12_Fig15_Plots.r draws the figures.
#
# Outputs (Part12_rerun/07_Fig15/<section>/): stats, plot data and logs.
# =============================================================================

source(file.path(Sys.getenv("SEAAD_ROOT", getwd()), "Part12_rerun", "Part12_03_Config.r"))

SCRIPT <- "Part12_Fig15_Analysis.r"
FIG_DIR <- "07_Fig15"   # every Figure 15 section folder sits inside this one

# =============================================================================
# SECTION 07: 15A and 15B
#   15A: cov4 clinical model refitted with donor DecontX contamination.
#   15B: the same model with Module post-IP nuclei above each doublet cut removed.
# Outputs (07_Fig15/15AB_ContamDoublets/): Stats_15A.csv, Stats_15B.csv.
# =============================================================================
p12_section_07 <- function() {
    STEP <- "07"
    D <- p12_step_dirs(STEP, "15AB_ContamDoublets", parent = FIG_DIR)
    p12_clean(c(D$tables, D$figures))
    p12_start(STEP, SCRIPT, D)

    F_15A <- file.path(D$tables, "Stats_15A.csv")
    F_15B <- file.path(D$tables, "Stats_15B.csv")

    # One cov4 OLS model: standardised Module beta, 95% CI, OLS p, dR2 and n.
    # `extra` is an added covariate column.
    p12_fit_clinical <- function(d, y, x, covs, extra = NULL, min_n) {
        cols <- unique(c(y, x, covs, extra))
        d <- d[stats::complete.cases(d[, cols, drop = FALSE]), , drop = FALSE]
        out <- c(n = nrow(d), beta = NA, lo = NA, hi = NA, p = NA, dr2 = NA)
        if (nrow(d) < min_n) return(out)
        if (!isTRUE(stats::sd(d[[y]]) > 0) || !isTRUE(stats::sd(d[[x]]) > 0)) return(out)
        d$.y <- as.numeric(scale(d[[y]])); d$.x <- as.numeric(scale(d[[x]]))
        rhs0 <- c(sprintf("`%s`", c(covs, extra)))
        f0 <- stats::as.formula(paste(".y ~", if (length(rhs0)) paste(rhs0, collapse = " + ") else "1"))
        f1 <- stats::update(f0, . ~ . + .x)
        m0 <- stats::lm(f0, data = d); m1 <- stats::lm(f1, data = d)
        cf <- summary(m1)$coefficients
        if (!(".x" %in% rownames(cf))) return(out)
        ci <- stats::confint(m1, ".x", level = 0.95)
        out[c("beta", "lo", "hi", "p")] <- c(cf[".x", "Estimate"], ci[1], ci[2],
                                             cf[".x", "Pr(>|t|)"])
        out["dr2"] <- summary(m1)$r.squared - summary(m0)$r.squared
        out
    }

    main <- function() {
        CL <- P12_CFG$clinical
        COVS <- unlist(CL$cov4)
        OUT <- vapply(CL$outcomes, function(o) o$column, "")
        MIN_N <- as.integer(P12_CFG$floors$donors_per_test)
        CUTS <- as.numeric(P12_CFG$doublet_cuts)

        p12_reading("15A")
        p12_reading("15B")

        # == Inputs ==========================================================
        p12_say("SECTION A -- inputs")
        DN <- utils::read.csv(p12_upstream_file("05_Units", "05_donors.csv"),
                              check.names = FALSE,
                              colClasses = c(donor_key = "character", donor_id = "character"))
        L <- utils::read.csv(p12_upstream_file("04_Labels", "04_cell_labels.csv.gz"),
                             check.names = FALSE,
                             colClasses = c(barcode = "character", donor_key = "character"))
        p12_chk(all(DN$labels_version == P12_CFG$labels$labels_version) &&
                    all(L$labels_version == P12_CFG$labels$labels_version),
                "donor table and labels carry this labels version", "checked")
        miss <- setdiff(c(COVS, OUT, "module_post_pct", "contam_super_nonmodule", "n_super",
                          "n_post"), names(DN))
        if (length(miss)) p12_abort(paste("donor table lacks:", paste(miss, collapse = ", ")))
        p12_note("Module post-IP% denominator", P12_CFG$quantities$module_post_pct_denominator)
        p12_say("")

        rows_a <- list(); rows_b <- list()
        add <- function(lst, panel, arm, ds, quantity, measure, fit, family) {
            lst[[length(lst) + 1L]] <- data.frame(
                panel = panel, arm = arm, level = ds, quantity = quantity, measure = measure,
                n_units = NA_integer_, n_donors = as.integer(fit[["n"]]),
                estimate = unname(fit[["est"]]), lo = unname(fit[["lo"]]),
                hi = unname(fit[["hi"]]), p = unname(fit[["p"]]), fdr = NA_real_,
                family = family)
            lst
        }

        # == 15A =============================================================
        p12_say("SECTION B -- 15A, Module post-IP% before and after the contamination covariate")
        for (a in P12_ARMS) for (ds in P12_DATASETS) {
            d <- DN[DN$arm == a & DN$dataset == ds, ]
            for (k in names(OUT)) {
                y <- OUT[[k]]
                cv <- if (y %in% COVS) setdiff(COVS, y) else COVS
                u <- p12_fit_clinical(d, y, "module_post_pct", cv, min_n = MIN_N)
                dd <- d[is.finite(d$contam_super_nonmodule), ]
                us <- p12_fit_clinical(dd, y, "module_post_pct", cv, min_n = MIN_N)
                ad <- p12_fit_clinical(d, y, "module_post_pct", cv,
                                       extra = "contam_super_nonmodule", min_n = MIN_N)
                fam_u <- paste("15A", a, ds, "unadjusted", sep = "|")
                fam_a <- paste("15A", a, ds, "contamination-adjusted", sep = "|")
                fam_d <- paste("15A", a, ds, "dR2", sep = "|")
                rows_a <- add(rows_a, "15A", a, ds, "module_beta_unadjusted", k,
                              c(n = u[["n"]], est = u[["beta"]], lo = u[["lo"]], hi = u[["hi"]],
                                p = u[["p"]]), fam_u)
                rows_a <- add(rows_a, "15A", a, ds, "module_beta_contam_adjusted", k,
                              c(n = ad[["n"]], est = ad[["beta"]], lo = ad[["lo"]], hi = ad[["hi"]],
                                p = ad[["p"]]), fam_a)
                rows_a <- add(rows_a, "15A", a, ds, "dR2_unadjusted", k,
                              c(n = u[["n"]], est = u[["dr2"]], lo = NA, hi = NA, p = NA), fam_d)
                rows_a <- add(rows_a, "15A", a, ds, "dR2_unadjusted_same_donors", k,
                              c(n = us[["n"]], est = us[["dr2"]], lo = NA, hi = NA, p = NA), fam_d)
                rows_a <- add(rows_a, "15A", a, ds, "dR2_contam_adjusted", k,
                              c(n = ad[["n"]], est = ad[["dr2"]], lo = NA, hi = NA, p = NA), fam_d)
                ret <- if (is.finite(us[["dr2"]]) && us[["dr2"]] > 0) ad[["dr2"]] / us[["dr2"]] else NA
                rows_a <- add(rows_a, "15A", a, ds, "retention", k,
                              c(n = ad[["n"]], est = ret, lo = NA, hi = NA, p = NA), fam_d)
            }
        }
        A <- do.call(rbind, rows_a)
        for (f in unique(A$family)) { i <- A$family == f; A$fdr[i] <- p12_bh(A$p[i]) }
        for (a in P12_ARMS) for (ds in P12_DATASETS) {
            u <- A[A$arm == a & A$level == ds & A$quantity == "module_beta_unadjusted", ]
            v <- A[A$arm == a & A$level == ds & A$quantity == "module_beta_contam_adjusted", ]
            p12_say(sprintf("  %-6s %-8s donors %3d/%3d | FDR<0.05 before %d, after %d | lost %d",
                            a, ds, max(u$n_donors, na.rm = TRUE), max(v$n_donors, na.rm = TRUE),
                            sum(u$fdr < 0.05, na.rm = TRUE), sum(v$fdr < 0.05, na.rm = TRUE),
                            sum(u$fdr < 0.05 & !(v$fdr < 0.05), na.rm = TRUE)))
        }
        p12_say("")

        # == 15B =============================================================
        p12_say("SECTION C -- 15B, Module post-IP nuclei removed at each doublet threshold")
        denom_sup <- identical(P12_CFG$quantities$module_post_pct_denominator, "supertype")
        for (cut in CUTS) {
            tag <- sprintf("%.2f", cut)
            rm_n <- with(L[L$group == "Module post-IP" & is.finite(L$doublet_score) &
                               L$doublet_score >= cut, ],
                         tapply(barcode, paste(arm, dataset, donor_key, sep = "|"), length))
            key <- paste(DN$arm, DN$dataset, DN$donor_key, sep = "|")
            r <- as.numeric(rm_n[key]); r[is.na(r)] <- 0
            den <- if (denom_sup) DN$n_super else DN$n_cells
            DN[[paste0("module_nodbl_", tag)]] <- 100 * (DN$n_post - r) / (den - r)
            for (a in P12_ARMS) for (ds in P12_DATASETS) {
                i <- DN$arm == a & DN$dataset == ds
                share <- sum(r[i]) / max(1, sum(DN$n_post[i]))
                fam <- paste("15B", a, ds, tag, sep = "|")
                rows_b <- add(rows_b, "15B", a, ds, paste0("removed_share_", tag), "all",
                              c(n = sum(i), est = share, lo = NA, hi = NA, p = NA),
                              paste("15B", a, ds, "removed", sep = "|"))
                d <- DN[i, ]
                for (k in names(OUT)) {
                    y <- OUT[[k]]
                    cv <- if (y %in% COVS) setdiff(COVS, y) else COVS
                    fit <- p12_fit_clinical(d, y, paste0("module_nodbl_", tag), cv, min_n = MIN_N)
                    rows_b <- add(rows_b, "15B", a, ds, paste0("module_beta_doublet_", tag), k,
                                  c(n = fit[["n"]], est = fit[["beta"]], lo = fit[["lo"]],
                                    hi = fit[["hi"]], p = fit[["p"]]), fam)
                    rows_b <- add(rows_b, "15B", a, ds, paste0("dR2_doublet_", tag), k,
                                  c(n = fit[["n"]], est = fit[["dr2"]], lo = NA, hi = NA, p = NA),
                                  paste("15B", a, ds, "dR2", sep = "|"))
                }
                p12_say(sprintf("  cut %s %-6s %-8s Module post-IP nuclei removed %5.1f%%", tag, a, ds,
                                100 * share))
            }
        }
        B <- do.call(rbind, rows_b)
        for (f in unique(B$family)) { i <- B$family == f; B$fdr[i] <- p12_bh(B$p[i]) }
        p12_say("")

        # == Write ===========================================================
        p12_say("SECTION D -- write")
        fin <- function(X) {
            X$run_id <- P12$run_id
            X$config_version <- P12_CFG$config_version
            X$labels_version <- P12_CFG$labels$labels_version
            rownames(X) <- NULL
            X[, P12_CFG$stats_columns]
        }
        A <- fin(A); B <- fin(B)
        p12_chk(all(is.finite(A$p[A$quantity == "module_beta_unadjusted"])),
                "15A unadjusted model fitted for every arm, region set and outcome",
                sprintf("%d of %d", sum(is.finite(A$p[A$quantity == "module_beta_unadjusted"])),
                        sum(A$quantity == "module_beta_unadjusted")), warn_only = TRUE)
        if (P12$n[["FAIL"]] > 0L) p12_abort("checks failed; stats not written")
        p12_write_stats(A, F_15A)
        p12_write_stats(B, F_15B)
        p12_chk(nrow(utils::read.csv(F_15A)) == nrow(A) && nrow(utils::read.csv(F_15B)) == nrow(B),
                "stats files written and read back",
                sprintf("15A %d rows | 15B %d rows", nrow(A), nrow(B)))
    }

    p12_run(main)
}

# =============================================================================
# SECTION 08: 15C and 15D
#   15C: DecontX contamination by cell group, donor-level contrasts.
#   15D: GOI-composite rs by within-donor contamination tertile, T3 minus T1.
# Outputs (07_Fig15/15CD_DecontXStrata/): Stats_15C.csv, Stats_15D.csv, Plot_15D_tertiles.csv.
# =============================================================================
p12_section_08 <- function() {
    STEP <- "08"
    D <- p12_step_dirs(STEP, "15CD_DecontXStrata", parent = FIG_DIR)
    p12_clean(c(D$tables, D$figures))
    p12_start(STEP, SCRIPT, D)

    F_15C <- file.path(D$tables, "Stats_15C.csv")
    F_15D <- file.path(D$tables, "Stats_15D.csv")
    F_15D_PLOT <- file.path(D$tables, "Plot_15D_tertiles.csv")

    # 15D bootstrap on one chunk of resamples: rs(T1) and rs(T3) per resample.
    # Argument names avoid x, which parLapply() passes on to clusterApply().
    p12_15d_chunk <- function(idx, dons, by_dt, gx, gy) {
        out <- matrix(NA_real_, nrow(idx), 2L)
        for (b in seq_len(nrow(idx))) {
            db <- dons[idx[b, ]]
            for (k in 1:2) {
                rows <- unlist(by_dt[paste(db, c(1L, 3L)[k], sep = ".")], use.names = FALSE)
                out[b, k] <- suppressWarnings(stats::cor(gx[rows], gy[rows], method = "spearman"))
            }
        }
        out
    }

    # Median of each row of a matrix, without apply() over 10,000 rows.
    p12_row_median <- function(M) {
        M <- t(apply(M, 1L, sort.int, method = "quick"))
        n <- ncol(M)
        if (n %% 2L) M[, (n + 1L) / 2L] else (M[, n / 2L] + M[, n / 2L + 1L]) / 2
    }

    main <- function() {
        MINC <- as.integer(P12_CFG$floors$cells_per_group_per_donor)
        MIND <- as.integer(P12_CFG$floors$donors_per_test)
        B <- as.integer(P12_CFG$inference$bootstrap_B)

        p12_reading("15C")
        p12_reading("15D")

        # == Inputs ==========================================================
        p12_say("SECTION A -- inputs")
        L <- utils::read.csv(p12_upstream_file("04_Labels", "04_cell_labels.csv.gz"),
                             check.names = FALSE,
                             colClasses = c(barcode = "character", donor_key = "character",
                                            group = "character", arm = "character"))
        p12_chk(all(L$labels_version == P12_CFG$labels$labels_version),
                "labels carry this labels version", "checked")
        f_mtx <- p12_upstream_file("02_Extract", "02_counts.mtx.gz")
        f_gen <- p12_upstream_file("02_Extract", "02_counts_genes.tsv")
        f_cel <- p12_upstream_file("02_Extract", "02_counts_cells.tsv")
        X <- methods::as(Matrix::readMM(f_mtx), "CsparseMatrix")
        genes <- readLines(f_gen); cells <- readLines(f_cel)
        p12_chk(nrow(X) == length(cells) && ncol(X) == length(genes),
                "step 02 count matrix matches its cell and gene lists", "checked")
        p12_chk(all(is.finite(L$contamination)), "every labelled cell has a contamination value",
                sprintf("%d missing", sum(!is.finite(L$contamination))))
        L$grp3 <- ifelse(L$group == "Module post-IP", "post",
                         ifelse(L$group == "Module pre-IP", "pre", "non"))
        p12_say("")

        # == 15C =============================================================
        p12_say(sprintf("SECTION B -- 15C, donor-level contrasts (%d cells per group, %d donors, B = %d)",
                        MINC, MIND, B))
        CONTR <- list(post_minus_pre = c("post", "pre"), post_minus_nonmodule = c("post", "non"),
                      pre_minus_nonmodule = c("pre", "non"))
        rows_c <- list()
        for (a in P12_ARMS) {
            la <- L[L$arm == a, ]
            med <- tapply(la$contamination, list(la$donor_key, la$grp3), stats::median)
            cnt <- tapply(la$contamination, list(la$donor_key, la$grp3), length)
            cnt[is.na(cnt)] <- 0L
            for (g in c("post", "pre", "non")) {
                ok <- cnt[, g] >= MINC
                v <- med[ok, g]
                rows_c[[length(rows_c) + 1L]] <- data.frame(
                    panel = "15C", arm = a, level = "pooled region sets",
                    quantity = paste0("donor_median_", g), measure = "DecontX contamination",
                    n_units = NA_integer_, n_donors = sum(ok),
                    estimate = if (length(v)) stats::median(v) else NA_real_, lo = NA_real_,
                    hi = NA_real_, p = NA_real_, fdr = NA_real_,
                    family = paste("15C", a, "group medians", sep = "|"))
            }
            for (k in names(CONTR)) {
                g1 <- CONTR[[k]][1]; g2 <- CONTR[[k]][2]
                ok <- cnt[, g1] >= MINC & cnt[, g2] >= MINC
                dd <- med[ok, g1] - med[ok, g2]
                est <- lo <- hi <- p <- NA_real_
                if (sum(ok) >= MIND) {
                    idx <- p12_boot_index(names(dd), STEP, B)
                    bs <- p12_row_median(matrix(dd[idx], nrow = B))
                    est <- stats::median(dd)
                    ci <- p12_boot_ci(bs); lo <- ci[["lo"]]; hi <- ci[["hi"]]
                    p <- p12_boot_p(bs)
                }
                rows_c[[length(rows_c) + 1L]] <- data.frame(
                    panel = "15C", arm = a, level = "pooled region sets", quantity = k,
                    measure = "DecontX contamination", n_units = NA_integer_,
                    n_donors = sum(ok), estimate = est, lo = lo, hi = hi, p = p, fdr = NA_real_,
                    family = paste("15C", a, sep = "|"))
                p12_say(sprintf("  %-6s %-22s donors %3d | median difference %+.4f [%+.4f, %+.4f] | p %s",
                                a, k, sum(ok), est, lo, hi, p12_fmt(p, p12_boot_floor(B))))
            }
        }
        C <- do.call(rbind, rows_c)
        for (f in unique(C$family)) { i <- C$family == f; C$fdr[i] <- p12_bh(C$p[i]) }
        p12_say("")

        # == 15D =============================================================
        p12_say(sprintf("SECTION C -- 15D, co-expression by within-donor contamination tertile (B = %d)", B))
        # Donor draws made once here from the step seed, so the result does not
        # depend on the number of workers.
        n_work <- p12_n_workers(B)
        p12_note("workers", sprintf("%d of %d logical cores (PSOCK)", n_work,
                                    parallel::detectCores(logical = TRUE)))
        cl <- NULL
        if (n_work > 1L) {
            cl <- parallel::makeCluster(n_work, type = "PSOCK")
            on.exit(parallel::stopCluster(cl), add = TRUE)
        }
        rows_d <- list(); plot_d <- list()
        jx <- match(L$barcode, cells)
        p12_chk(all(!is.na(jx)), "every labelled cell is in the step 02 count matrix",
                sprintf("%d missing", sum(is.na(jx))))
        tot <- utils::read.csv(p12_upstream_file("02_Extract", "02_cells.csv.gz"),
                               colClasses = c(barcode = "character"))
        tot <- tot$total_umis[match(L$barcode, tot$barcode)]
        for (a in P12_ARMS) {
            goi <- P12_CFG$arms[[a]]$goi
            gi <- match(goi, genes)
            if (is.na(gi)) p12_abort(sprintf("%s GOI %s missing from step 02", a, goi))
            i <- which(L$arm == a & L$is_module)
            m <- L[i, ]
            m$goi_expr <- log1p(as.numeric(X[jx[i], gi]) / tot[i] * 1e4)
            # Tertiles within donor; ties broken by barcode.
            m <- m[order(m$donor_key, m$contamination, m$barcode), ]
            n_d <- table(m$donor_key)
            keep_d <- names(n_d)[n_d >= MINC]
            m <- m[m$donor_key %in% keep_d, ]
            m$pos <- stats::ave(seq_len(nrow(m)), m$donor_key, FUN = seq_along)
            m$n_d <- as.integer(n_d[m$donor_key])
            m$tert <- pmin(3L, 1L + (3L * (m$pos - 1L)) %/% m$n_d)
            tie0 <- mean(m$contamination == 0)
            p12_say(sprintf("  %-6s GOI %-6s donors %3d | Module cells %5d | cells at zero contamination %.1f%%",
                            a, goi, length(keep_d), nrow(m), 100 * tie0))
            rs_t <- function(rows) suppressWarnings(stats::cor(m$goi_expr[rows], m$composite[rows],
                                                               method = "spearman"))
            by_dt <- split(seq_len(nrow(m)), list(m$donor_key, m$tert), drop = FALSE)
            dons <- keep_d
            obs <- vapply(1:3, function(t) rs_t(which(m$tert == t)), 0)
            # Plot data: median contamination of each tertile.
            plot_d[[a]] <- data.frame(arm = a, tertile = 1:3,
                                      n_cells = as.integer(tabulate(m$tert, 3L)),
                                      contamination_median = vapply(1:3, function(t)
                                          stats::median(m$contamination[m$tert == t]), 0),
                                      rs = obs)
            est_d <- obs[3] - obs[1]
            bt <- matrix(NA_real_, B, 3)
            if (length(dons) >= MIND) {
                idx <- p12_boot_index(dons, STEP, B)
                ch <- p12_chunks(B, n_work)
                parts <- lapply(ch, function(r) idx[r, , drop = FALSE])
                t0 <- Sys.time()
                res <- if (is.null(cl))
                    lapply(parts, p12_15d_chunk, dons = dons, by_dt = by_dt,
                           gx = m$goi_expr, gy = m$composite) else
                    parallel::parLapply(cl, parts, p12_15d_chunk, dons = dons, by_dt = by_dt,
                                        gx = m$goi_expr, gy = m$composite)
                R <- do.call(rbind, res)
                bt[, 1] <- R[, 1]; bt[, 3] <- R[, 2]
                p12_say(sprintf("         bootstrap %.1f min", as.numeric(difftime(Sys.time(), t0,
                                                                                  units = "mins"))))
            }
            dif <- bt[, 3] - bt[, 1]
            for (t in 1:3) {
                ci <- if (t == 2L) c(lo = NA_real_, hi = NA_real_) else p12_boot_ci(bt[, t])
                rows_d[[length(rows_d) + 1L]] <- data.frame(
                    panel = "15D", arm = a, level = "pooled region sets",
                    quantity = paste0("rs_T", t), measure = sprintf("%s vs composite", goi),
                    n_units = sum(m$tert == t), n_donors = length(dons), estimate = obs[t],
                    lo = ci[["lo"]], hi = ci[["hi"]], p = NA_real_, fdr = NA_real_,
                    family = paste("15D", a, "tertiles", sep = "|"))
            }
            ok <- length(dons) >= MIND
            if (!ok) est_d <- NA_real_
            ci <- if (ok) p12_boot_ci(dif) else c(lo = NA_real_, hi = NA_real_)
            pv <- if (ok) p12_boot_p(dif) else NA_real_
            rows_d[[length(rows_d) + 1L]] <- data.frame(
                panel = "15D", arm = a, level = "pooled region sets",
                quantity = "rs_T3_minus_T1", measure = sprintf("%s vs composite", goi),
                n_units = nrow(m), n_donors = length(dons), estimate = est_d,
                lo = ci[["lo"]], hi = ci[["hi"]], p = pv, fdr = NA_real_, family = "15D")
            p12_say(sprintf("         rs T1 %.3f | T2 %.3f | T3 %.3f | T3 - T1 %+.3f [%+.3f, %+.3f] | p %s",
                            obs[1], obs[2], obs[3], est_d, ci[["lo"]], ci[["hi"]],
                            p12_fmt(pv, p12_boot_floor(B))))
        }
        Dd <- do.call(rbind, rows_d)
        for (f in unique(Dd$family)) { i <- Dd$family == f; Dd$fdr[i] <- p12_bh(Dd$p[i]) }
        p12_say("")

        # == Write ===========================================================
        p12_say("SECTION D -- write")
        fin <- function(Z) {
            Z$run_id <- P12$run_id
            Z$config_version <- P12_CFG$config_version
            Z$labels_version <- P12_CFG$labels$labels_version
            rownames(Z) <- NULL
            Z[, P12_CFG$stats_columns]
        }
        C <- fin(C); Dd <- fin(Dd)
        if (P12$n[["FAIL"]] > 0L) p12_abort("checks failed; stats not written")
        p12_write_stats(C, F_15C)
        p12_write_stats(Dd, F_15D)
        PD <- do.call(rbind, plot_d); rownames(PD) <- NULL
        p12_output(F_15D_PLOT)
        utils::write.csv(PD, F_15D_PLOT, row.names = FALSE)
        p12_chk(nrow(utils::read.csv(F_15C)) == nrow(C) && nrow(utils::read.csv(F_15D)) == nrow(Dd),
                "stats files written and read back",
                sprintf("15C %d rows | 15D %d rows", nrow(C), nrow(Dd)))
    }

    p12_run(main)
}

# =============================================================================
# SECTION 09: 15E, 15F and 15G (Gerrits)
#   15E: ambient against nuclear %MT per sample. 15F: raw and partial rs of
#   MT-CO3 with its partners. 15G: top minus bottom soup-index tertile.
# Outputs (07_Fig15/15EG_GerritsAmbient/): Stats_15E/F/G.csv.
# =============================================================================
p12_section_09 <- function() {
    STEP <- "09"
    D <- p12_step_dirs(STEP, "15EG_GerritsAmbient", parent = FIG_DIR)
    p12_clean(c(D$tables, D$figures))
    p12_start(STEP, SCRIPT, D)

    F_15E <- file.path(D$tables, "Stats_15E.csv")
    F_15F <- file.path(D$tables, "Stats_15F.csv")
    F_15G <- file.path(D$tables, "Stats_15G.csv")
    F_PLOT_E <- file.path(D$tables, "Plot_15E_samples.csv")
    F_PLOT_F <- file.path(D$tables, "Plot_15F_samples.csv")
    F_PLOT_G <- file.path(D$tables, "Plot_15G_tertiles.csv")

    # Partial Spearman: correlation of the ranks of x and y after the ranks of Z
    # are regressed out.
    p12_partial_spearman <- function(x, y, Z) {
        R <- cbind(1, apply(Z, 2L, rank))
        rx <- stats::lm.fit(R, rank(x))$residuals
        ry <- stats::lm.fit(R, rank(y))$residuals
        if (stats::sd(rx) == 0 || stats::sd(ry) == 0) return(NA_real_)
        stats::cor(rx, ry)
    }

    p12_spearman <- function(x, y) {
        if (stats::sd(x) == 0 || stats::sd(y) == 0) return(NA_real_)
        stats::cor(x, y, method = "spearman")
    }

    # Median over samples with a donor bootstrap. `v`: one value per sample (NA
    # below its floor); `donor`: the sample's donor.
    p12_boot_samples <- function(v, donor, B) {
        ok <- is.finite(v)
        out <- list(est = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_,
                    n_samples = sum(ok), n_donors = length(unique(donor[ok])), boot = NULL,
                    donors = sort(unique(donor[ok])))
        if (out$n_donors < as.integer(P12_CFG$floors$donors_per_test)) return(out)
        dons <- sort(unique(donor[ok]))
        per <- lapply(dons, function(d) v[ok & donor == d])
        k <- max(lengths(per))
        V <- t(vapply(per, function(z) c(z, rep(NA_real_, k - length(z))), numeric(k)))
        if (k == 1L) V <- matrix(V, ncol = 1L)
        idx <- p12_boot_index(dons, STEP, B)
        M <- do.call(cbind, lapply(seq_len(k), function(j) matrix(V[idx, j], nrow = B)))
        bs <- apply(M, 1L, stats::median, na.rm = TRUE)
        ci <- p12_boot_ci(bs)
        out$est <- stats::median(v[ok]); out$lo <- ci[["lo"]]; out$hi <- ci[["hi"]]
        out$p <- p12_boot_p(bs); out$boot <- bs
        out
    }

    main <- function() {
        G <- P12_CFG$gerrits
        B <- as.integer(P12_CFG$inference$bootstrap_B)
        PART <- P12_CFG$genes$cdfg
        MIN_S <- as.integer(G$min_cells_per_sample)
        MIN_T <- as.integer(G$min_cells_per_tertile)

        p12_reading("15E")
        p12_reading("15F")
        p12_reading("15G")

        # == Inputs ==========================================================
        p12_say("SECTION A -- inputs")
        S <- utils::read.csv(p12_upstream_file("06_GerritsCells", "06_gerrits_samples.csv"),
                             check.names = FALSE)
        CE <- utils::read.csv(p12_upstream_file("06_GerritsCells", "06_gerrits_cells.csv.gz"),
                              check.names = FALSE)
        pc <- paste0("dfg_", PART)
        miss <- setdiff(c("sample", "donor", "total_umis", "mt_umis", "goi", "soup_index", pc),
                        names(CE))
        if (length(miss)) p12_abort(paste("Gerrits cell table lacks:", paste(miss, collapse = ", ")))
        p12_chk(nrow(S) == G$n_samples && length(unique(S$donor)) == G$n_donors,
                "36 samples from 18 donors", sprintf("%d samples | %d donors", nrow(S),
                                                      length(unique(S$donor))))
        rows <- list()
        add <- function(panel, quantity, measure, r, family) {
            rows[[length(rows) + 1L]] <<- data.frame(
                panel = panel, arm = "Gerrits", level = "sample", quantity = quantity,
                measure = measure, n_units = r$n_samples, n_donors = r$n_donors,
                estimate = r$est, lo = r$lo, hi = r$hi, p = r$p, fdr = NA_real_, family = family)
        }
        p12_say("")

        # == 15E =============================================================
        p12_say("SECTION B -- 15E, ambient against nuclear %MT")
        nuc <- tapply(CE$mt_umis, CE$sample, sum) / tapply(CE$total_umis, CE$sample, sum)
        S$nuclear_mito_frac <- as.numeric(nuc[S$sample])
        S$log_ratio <- log(S$ambient_mito_frac / S$nuclear_mito_frac)
        p12_chk(all(is.finite(S$log_ratio)), "log ratio finite on every sample",
                sprintf("%d not finite", sum(!is.finite(S$log_ratio))))
        r <- p12_boot_samples(S$log_ratio, S$donor, B)
        add("15E", "log_ambient_over_nuclear_mito", "pooled MT fraction", r, "15E")
        p12_say(sprintf("  median log ratio %.3f [%.3f, %.3f] (ratio %.2f) | samples above 0: %d of %d | p %s",
                        r$est, r$lo, r$hi, exp(r$est), sum(S$log_ratio > 0, na.rm = TRUE),
                        nrow(S), p12_fmt(r$p, p12_boot_floor(B))))
        p12_say("")

        # == Per-sample correlations (15F and 15G) ===========================
        p12_say("SECTION C -- per-sample correlations in the MT-CO3-positive nuclei")
        P <- CE[CE$goi > 0, ]
        ex <- function(cnt) log1p(cnt / P$total_umis * 1e4)
        P$x <- ex(P$goi)
        for (g in PART) P[[paste0("y_", g)]] <- ex(P[[paste0("dfg_", g)]])
        P$log_depth <- log(P$total_umis)
        raw <- part <- tdif <- matrix(NA_real_, nrow(S), length(PART), dimnames = list(S$sample, PART))
        tert_rs <- list()
        n_pop <- n_tmin <- stats::setNames(integer(nrow(S)), S$sample)
        for (s in S$sample) {
            q <- P[P$sample == s, ]
            n_pop[s] <- nrow(q)
            if (nrow(q) >= MIN_S) {
                Z <- cbind(q$soup_index, q$log_depth)
                for (g in PART) {
                    y <- q[[paste0("y_", g)]]
                    raw[s, g] <- p12_spearman(q$x, y)
                    part[s, g] <- p12_partial_spearman(q$x, y, Z)
                }
            }
            # Soup tertiles within sample, ties broken by cell index.
            o <- order(q$soup_index, q$cell)
            tt <- integer(nrow(q))
            tt[o] <- pmin(3L, 1L + (3L * (seq_along(o) - 1L)) %/% max(1L, nrow(q)))
            n_tmin[s] <- if (nrow(q)) min(tabulate(tt, 3L)) else 0L
            if (n_tmin[s] >= MIN_T) {
                for (g in PART) {
                    y <- q[[paste0("y_", g)]]
                    r3 <- vapply(1:3, function(t) p12_spearman(q$x[tt == t], y[tt == t]), 0)
                    tdif[s, g] <- r3[3] - r3[1]
                    tert_rs[[length(tert_rs) + 1L]] <- data.frame(sample = s, partner = g,
                                                                  tertile = 1:3, rs = r3)
                }
            }
        }
        zero_soup <- tapply(P$soup_index == 0, P$sample, mean)
        p12_note("share of population nuclei with a soup index of 0, per sample",
                 sprintf("median %.1f%% (range %.1f%% to %.1f%%); above two thirds, the top tertile also holds zeros",
                         100 * stats::median(zero_soup), 100 * min(zero_soup), 100 * max(zero_soup)))
        med_row <- function(M) apply(M, 1L, function(z) if (all(is.na(z))) NA_real_ else
            stats::median(z, na.rm = TRUE))
        S$rs_raw <- med_row(raw); S$rs_partial <- med_row(part); S$rs_tert_diff <- med_row(tdif)
        p12_say(sprintf("  population per sample: median %d nuclei (range %d to %d) | 15F samples %d | 15G samples %d",
                        as.integer(stats::median(n_pop)), min(n_pop), max(n_pop),
                        sum(is.finite(S$rs_raw)), sum(is.finite(S$rs_tert_diff))))
        p12_say("")

        # == 15F =============================================================
        p12_say("SECTION D -- 15F, raw and partial co-expression")
        rr <- p12_boot_samples(S$rs_raw, S$donor, B)
        rp <- p12_boot_samples(S$rs_partial, S$donor, B)
        add("15F", "rs_raw", "median over 8 partners", rr, "15F|primary")
        add("15F", "rs_partial", "median over 8 partners", rp, "15F|primary")
        ratio <- list(est = rp$est / rr$est, lo = NA_real_, hi = NA_real_, p = NA_real_,
                      n_samples = rp$n_samples, n_donors = rp$n_donors)
        # The ratio interval pairs partial and raw resamples of the same draw.
        if (!is.null(rp$boot) && !is.null(rr$boot) && identical(rp$donors, rr$donors)) {
            ci <- p12_boot_ci(rp$boot / rr$boot); ratio$lo <- ci[["lo"]]; ratio$hi <- ci[["hi"]]
        }
        add("15F", "partial_over_raw", "median over 8 partners", ratio, "15F|ratio")
        # Per-sample raw spread (median, quartiles): the Gerrits range read by 15H.
        qv <- stats::quantile(S$rs_raw, c(0.25, 0.5, 0.75), na.rm = TRUE, names = FALSE)
        add("15F", "rs_raw_per_sample_quartiles", "median over 8 partners",
            list(est = qv[2], lo = qv[1], hi = qv[3], p = NA_real_,
                 n_samples = sum(is.finite(S$rs_raw)),
                 n_donors = length(unique(S$donor[is.finite(S$rs_raw)]))), "15F|reference")
        p12_say(sprintf("  raw %.3f [%.3f, %.3f] p %s | partial %.3f [%.3f, %.3f] p %s | partial/raw %.3f",
                        rr$est, rr$lo, rr$hi, p12_fmt(rr$p, p12_boot_floor(B)),
                        rp$est, rp$lo, rp$hi, p12_fmt(rp$p, p12_boot_floor(B)), ratio$est))
        for (g in PART) {
            a <- p12_boot_samples(raw[, g], S$donor, B)
            b <- p12_boot_samples(part[, g], S$donor, B)
            add("15F", "rs_raw", g, a, "15F|partners raw")
            add("15F", "rs_partial", g, b, "15F|partners partial")
            p12_say(sprintf("    %-7s raw %.3f | partial %.3f", g, a$est, b$est))
        }
        p12_say("")

        # == 15G =============================================================
        p12_say("SECTION E -- 15G, top minus bottom soup-index tertile")
        rg <- p12_boot_samples(S$rs_tert_diff, S$donor, B)
        add("15G", "rs_top_minus_bottom", "median over 8 partners", rg, "15G|primary")
        p12_say(sprintf("  difference %+.4f [%+.4f, %+.4f] | p %s", rg$est, rg$lo, rg$hi,
                        p12_fmt(rg$p, p12_boot_floor(B))))
        for (g in PART) {
            a <- p12_boot_samples(tdif[, g], S$donor, B)
            add("15G", "rs_top_minus_bottom", g, a, "15G|partners")
            p12_say(sprintf("    %-7s %+.4f", g, a$est))
        }
        p12_say("")

        # == Write ===========================================================
        p12_say("SECTION F -- write")
        X <- do.call(rbind, rows)
        for (f in unique(X$family)) { i <- X$family == f; X$fdr[i] <- p12_bh(X$p[i]) }
        X$run_id <- P12$run_id
        X$config_version <- P12_CFG$config_version
        X$labels_version <- P12_CFG$labels$labels_version
        X <- X[, P12_CFG$stats_columns]
        if (P12$n[["FAIL"]] > 0L) p12_abort("checks failed; stats not written")
        p12_write_stats(X[X$panel == "15E", ], F_15E)
        p12_write_stats(X[X$panel == "15F", ], F_15F)
        p12_write_stats(X[X$panel == "15G", ], F_15G)
        # Plot data: per-sample values behind 15E to 15G.
        PE <- S[, c("sample", "donor", "ambient_mito_frac", "nuclear_mito_frac", "log_ratio")]
        PF <- data.frame(sample = rep(S$sample, length(PART)),
                         donor = rep(S$donor, length(PART)),
                         partner = rep(PART, each = nrow(S)),
                         rs_raw = as.vector(raw), rs_partial = as.vector(part))
        PG <- if (length(tert_rs)) do.call(rbind, tert_rs) else
            data.frame(sample = character(0), partner = character(0), tertile = integer(0),
                       rs = numeric(0))
        for (f in c(F_PLOT_E, F_PLOT_F, F_PLOT_G)) p12_output(f)
        utils::write.csv(PE, F_PLOT_E, row.names = FALSE)
        utils::write.csv(PF, F_PLOT_F, row.names = FALSE)
        utils::write.csv(PG, F_PLOT_G, row.names = FALSE)
        p12_chk(sum(vapply(c(F_15E, F_15F, F_15G), function(f) nrow(utils::read.csv(f)), 0L)) ==
                    nrow(X), "stats files written and read back", sprintf("%d rows", nrow(X)))
    }

    p12_run(main)
}

# =============================================================================
# SECTION 10: 15H
#   Depth thinning of the PVM_3 arm to target depths (and to Gerrits' MT
#   fraction), median rs over the 8 partners against the Gerrits IQR; no test.
# Outputs (07_Fig15/15H_DepthThinning/): Stats_15H.csv.
# =============================================================================
p12_section_10 <- function() {
    STEP <- "10"
    D <- p12_step_dirs(STEP, "15H_DepthThinning", parent = FIG_DIR)
    p12_clean(c(D$tables, D$figures))
    p12_start(STEP, SCRIPT, D)

    F_15H <- file.path(D$tables, "Stats_15H.csv")

    main <- function() {
        TH <- P12_CFG$thinning
        PART <- P12_CFG$genes$cdfg
        MT13 <- P12_CFG$genes$mt_coding
        GOI <- "MT-CO3"
        R_N <- as.integer(TH$repeats)

        p12_reading("15H")

        # == Inputs ==========================================================
        p12_say("SECTION A -- inputs")
        L <- utils::read.csv(p12_upstream_file("04_Labels", "04_cell_labels.csv.gz"),
                             check.names = FALSE,
                             colClasses = c(barcode = "character", arm = "character"))
        p12_chk(all(L$labels_version == P12_CFG$labels$labels_version),
                "labels carry this labels version", "checked")
        X <- methods::as(Matrix::readMM(p12_upstream_file("02_Extract", "02_counts.mtx.gz")),
                         "CsparseMatrix")
        genes <- readLines(p12_upstream_file("02_Extract", "02_counts_genes.tsv"))
        cells <- readLines(p12_upstream_file("02_Extract", "02_counts_cells.tsv"))
        C02 <- utils::read.csv(p12_upstream_file("02_Extract", "02_cells.csv.gz"),
                               colClasses = c(barcode = "character"))
        GS <- utils::read.csv(p12_upstream_file("06_GerritsCells", "06_gerrits_samples.csv"))
        F15 <- utils::read.csv(p12_upstream_file(file.path(FIG_DIR, "15EG_GerritsAmbient"), "Stats_15F.csv"))

        HK_GOI <- P12_CFG$arms$HKG8$goi
        HK_PART <- P12_CFG$genes[[P12_CFG$arms$HKG8$partners]]
        need <- unique(c(GOI, PART, MT13, HK_GOI, HK_PART))
        gi <- match(need, genes)
        if (anyNA(gi)) p12_abort(paste("genes missing from step 02:",
                                       paste(need[is.na(gi)], collapse = ", ")))
        bc <- L$barcode[L$arm == "PVM_3"]
        ci <- match(bc, cells)
        p12_chk(all(!is.na(ci)), "every PVM_3 arm cell is in the step 02 count matrix",
                sprintf("%d of %d", sum(!is.na(ci)), length(ci)))
        Y <- as.matrix(X[ci, gi, drop = FALSE]); colnames(Y) <- need
        tot <- as.numeric(C02$total_umis[match(bc, C02$barcode)])
        pos <- Y[, GOI] > 0
        p12_chk(mean(pos) == 1, "every PVM_3 arm cell has MT-CO3 above zero (the arm population)",
                sprintf("%d of %d", sum(pos), length(pos)), warn_only = TRUE)
        Y <- Y[pos, , drop = FALSE]; tot <- tot[pos]
        mt_cols <- match(MT13, need)

        g_depth <- round(stats::median(GS$median_umis))
        g_mito <- stats::median(GS$median_mito_frac)
        ref <- F15[F15$quantity == "rs_raw_per_sample_quartiles", ]
        p12_chk(nrow(ref) == 1L && all(is.finite(c(ref$lo, ref$hi))),
                "Gerrits reference range read from Stats_15F",
                if (nrow(ref) == 1L) sprintf("IQR %.3f to %.3f (median %.3f)", ref$lo, ref$hi,
                                             ref$estimate) else "row missing: rerun step 09")
        if (nrow(ref) != 1L) p12_abort("Stats_15F has no per-sample quartile row")
        targets <- sort(unique(c(as.numeric(TH$target_depths),
                                 if (isTRUE(TH$add_gerrits_depth)) g_depth)), decreasing = TRUE)
        p12_say(sprintf("  population: %d cells | median depth %.0f UMIs | median MT fraction %.4f",
                        nrow(Y), stats::median(tot), stats::median(rowSums(Y[, mt_cols]) / tot)))
        p12_say(sprintf("  Gerrits: depth %d UMIs | MT fraction %.5f | targets: %s | repeats %d",
                        g_depth, g_mito, paste(targets, collapse = ", "), R_N))
        p12_say("")

        # == Simulation ======================================================
        p12_say("SECTION B -- thinning")
        rs_med1 <- function(Z, dep, goi, part) {
            ok <- Z[, goi] > 0 & dep > 0
            if (sum(ok) < 3L) return(c(rs = NA_real_, n = sum(ok)))
            g <- log1p(Z[ok, goi] / dep[ok] * 1e4)
            r <- vapply(part, function(p) {
                v <- log1p(Z[ok, p] / dep[ok] * 1e4)
                if (stats::sd(v) == 0 || stats::sd(g) == 0) NA_real_ else
                    stats::cor(g, v, method = "spearman")
            }, 0)
            c(rs = stats::median(r, na.rm = TRUE), n = sum(ok))
        }
        # Module series (MT-CO3, 8 partners) and, as a descriptive control on the
        # same cells, the housekeeping series (HKG8 GOI with its partners).
        rs_med <- function(Z, dep) {
            a <- rs_med1(Z, dep, GOI, PART)
            h <- rs_med1(Z, dep, HK_GOI, HK_PART)
            c(rs = a[["rs"]], n = a[["n"]], rs_hk = h[["rs"]], n_hk = h[["n"]])
        }
        set.seed(p12_seed(STEP))
        rows <- list()
        obs <- rs_med(Y, tot)
        rows[[1]] <- data.frame(level = "observed", target = NA_real_, rep = 0L, rs = obs[["rs"]],
                                n = obs[["n"]], rs_hk = obs[["rs_hk"]], pm = 1, depth = stats::median(tot),
                                mito = stats::median(rowSums(Y[, mt_cols]) / tot),
                                mt_umis = stats::median(rowSums(Y[, mt_cols])))
        for (tg in targets) {
            keep <- pmin(1, tg / pmax(tot, 1))
            dep <- pmax(1, tot * keep)
            for (r in seq_len(R_N)) {
                Z <- matrix(stats::rbinom(length(Y), as.vector(Y), rep(keep, ncol(Y))),
                            nrow = nrow(Y), dimnames = dimnames(Y))
                for (cond in vapply(TH$match_mito_fraction, function(m)
                    if (isTRUE(m)) "depth+mito" else "depth", "")) {
                    W <- Z
                    pm <- 1
                    if (cond == "depth+mito") {
                        cur <- stats::median(rowSums(W[, mt_cols]) / dep)
                        if (cur > 0) {
                            pm <- min(1, g_mito / cur)
                            W[, mt_cols] <- stats::rbinom(nrow(W) * length(mt_cols),
                                                          as.vector(W[, mt_cols]), pm)
                        }
                    }
                    o <- rs_med(W, dep)
                    rows[[length(rows) + 1L]] <- data.frame(
                        level = cond, target = tg, rep = r, rs = o[["rs"]], n = o[["n"]],
                        rs_hk = o[["rs_hk"]], pm = pm,
                        depth = stats::median(dep), mito = stats::median(rowSums(W[, mt_cols]) / dep),
                        mt_umis = stats::median(rowSums(W[, mt_cols])))
                }
            }
        }
        R <- do.call(rbind, rows)
        # Median and range over repeats; a repeat with no finite rs is counted and left out.
        Rs <- R[R$level != "observed", ]
        grp <- split(Rs, list(Rs$level, Rs$target), drop = TRUE)
        S <- do.call(rbind, lapply(grp, function(g) {
            f <- is.finite(g$rs)
            data.frame(level = g$level[1], target = g$target[1],
                       rs = if (any(f)) stats::median(g$rs[f]) else NA_real_,
                       lo = if (any(f)) min(g$rs[f]) else NA_real_,
                       hi = if (any(f)) max(g$rs[f]) else NA_real_,
                       rs_hk = stats::median(g$rs_hk, na.rm = TRUE),
                       hk_lo = suppressWarnings(min(g$rs_hk, na.rm = TRUE)),
                       hk_hi = suppressWarnings(max(g$rs_hk, na.rm = TRUE)),
                       n = stats::median(g$n), depth = stats::median(g$depth),
                       mito = stats::median(g$mito), mt_umis = stats::median(g$mt_umis),
                       pm = stats::median(g$pm),
                       n_rep = sum(f))
        }))
        rownames(S) <- NULL
        S <- S[order(S$level, -S$target), ]
        p12_chk(all(S$n_rep == R_N), "every repeat gave a finite rs",
                sprintf("%d of %d condition-target cells complete", sum(S$n_rep == R_N), nrow(S)),
                warn_only = TRUE)
        p12_say(sprintf("  observed: rs %.3f | %d cells | median MT fraction %.4f", obs[["rs"]],
                        as.integer(obs[["n"]]), rows[[1]]$mito))
        for (i in seq_len(nrow(S)))
            p12_say(sprintf("  %-10s target %5.0f | rs %.3f [%.3f, %.3f] | housekeeping %.3f | cells %6.0f | MT fraction %.5f | MT UMIs %.1f | MT keep %.4f",
                            S$level[i], S$target[i], S$rs[i], S$lo[i], S$hi[i], S$rs_hk[i], S$n[i],
                            S$mito[i], S$mt_umis[i], S$pm[i]))
        if (any(S$level == "depth+mito" & S$pm >= 1))
            p12_note("MT genes not thinned at some targets",
                     "the thinned MT fraction was already at or below Gerrits'")
        at <- S[S$target == g_depth, ]
        for (i in seq_len(nrow(at)))
            p12_say(sprintf("  at Gerrits' depth, %-10s rs %.3f against the Gerrits IQR %.3f to %.3f: %s",
                            at$level[i], at$rs[i], ref$lo, ref$hi,
                            if (at$rs[i] > ref$hi) "above" else if (at$rs[i] < ref$lo) "below" else
                                "inside"))
        p12_say("")

        # == Write ===========================================================
        p12_say("SECTION C -- write")
        mk <- function(level, quantity, measure, n, est, lo, hi) data.frame(
            panel = "15H", arm = "PVM_3", level = level, quantity = quantity, measure = measure,
            n_units = as.integer(n), n_donors = NA_integer_, estimate = est, lo = lo, hi = hi,
            p = NA_real_, fdr = NA_real_, family = "15H")
        out <- list(mk("observed", "rs_median_over_partners", "observed", obs[["n"]], obs[["rs"]],
                       NA, NA),
                    mk("observed", "rs_housekeeping", "observed", obs[["n_hk"]], obs[["rs_hk"]],
                       NA, NA),
                    mk("observed", "median_depth", "observed", obs[["n"]], stats::median(tot),
                       NA, NA))
        for (i in seq_len(nrow(S))) {
            tgt <- sprintf("target %d", as.integer(S$target[i]))
            out[[length(out) + 1L]] <- mk(S$level[i], "rs_median_over_partners", tgt, S$n[i],
                                          S$rs[i], S$lo[i], S$hi[i])
            out[[length(out) + 1L]] <- mk(S$level[i], "rs_housekeeping", tgt, NA, S$rs_hk[i],
                                          S$hk_lo[i], S$hk_hi[i])
            out[[length(out) + 1L]] <- mk(S$level[i], "median_depth", tgt, S$n[i], S$depth[i], NA, NA)
            out[[length(out) + 1L]] <- mk(S$level[i], "median_mito_fraction", tgt, S$n[i], S$mito[i],
                                          NA, NA)
        }
        out[[length(out) + 1L]] <- mk("Gerrits", "reference_depth", "median over samples",
                                      nrow(GS), g_depth, NA, NA)
        out[[length(out) + 1L]] <- mk("Gerrits", "reference_mito_fraction", "median over samples",
                                      nrow(GS), g_mito, NA, NA)
        out[[length(out) + 1L]] <- mk("Gerrits", "reference_rs_raw_iqr", "per-sample raw rs",
                                      ref$n_units, ref$estimate, ref$lo, ref$hi)
        H <- do.call(rbind, out)
        H$run_id <- P12$run_id
        H$config_version <- P12_CFG$config_version
        H$labels_version <- P12_CFG$labels$labels_version
        H <- H[, P12_CFG$stats_columns]
        if (P12$n[["FAIL"]] > 0L) p12_abort("checks failed; stats not written")
        p12_write_stats(H, F_15H)
        p12_chk(nrow(utils::read.csv(F_15H)) == nrow(H), "stats file written and read back",
                sprintf("%d rows", nrow(H)))
    }

    p12_run(main)
}

# =============================================================================
# SECTION 11: Supplementary Figure 10
#   Own cytoplasm or other lineages: own-ER, ribosomal and NEAT1 shares and the
#   other-lineage index, Module post-IP minus pre-IP per donor, median over donors.
# Outputs (07_Fig15/S10_Cargo/): Stats_S10.csv.
# =============================================================================
p12_section_11 <- function() {
    STEP <- "11"
    D <- p12_step_dirs(STEP, "S10_Cargo", parent = FIG_DIR)
    p12_clean(c(D$tables, D$figures))
    p12_start(STEP, SCRIPT, D)

    F_S10 <- file.path(D$tables, "Stats_S10.csv")

    # Median of each row of a matrix.
    p12_row_median <- function(M) {
        M <- t(apply(M, 1L, sort.int, method = "quick"))
        n <- ncol(M)
        if (n %% 2L) M[, (n + 1L) / 2L] else (M[, n / 2L] + M[, n / 2L + 1L]) / 2
    }

    main <- function() {
        MINC <- as.integer(P12_CFG$floors$cells_per_group_per_donor)
        MIND <- as.integer(P12_CFG$floors$donors_per_test)
        B <- as.integer(P12_CFG$inference$bootstrap_B)

        p12_reading("Supplementary Figure 10")

        # == Inputs ==========================================================
        p12_say("SECTION A -- inputs")
        L <- utils::read.csv(p12_upstream_file("04_Labels", "04_cell_labels.csv.gz"),
                             check.names = FALSE,
                             colClasses = c(barcode = "character", donor_key = "character",
                                            group = "character", arm = "character",
                                            dataset = "character"))
        p12_chk(all(L$labels_version == P12_CFG$labels$labels_version),
                "labels carry this labels version", "checked")
        soup_col <- paste0("soup_combined__", L$dataset)
        need <- c("own_er_fracnonmt", "ribosome_fracnonmt", "nuclear_retained_fracnonmt",
                  unique(soup_col))
        miss <- setdiff(need, names(L))
        if (length(miss)) p12_abort(paste("labels lack:", paste(miss, collapse = ", ")))
        # Other-lineage index under each cell's own dataset guard, as per cent.
        sc <- as.matrix(L[, unique(soup_col), drop = FALSE])
        L$other_lineage <- 100 * sc[cbind(seq_len(nrow(L)), match(soup_col, colnames(sc)))]
        MET <- c(own_er = "own_er_fracnonmt", ribosome = "ribosome_fracnonmt",
                 neat1 = "nuclear_retained_fracnonmt", other_lineage = "other_lineage")
        for (k in names(MET))
            p12_chk(mean(is.finite(L[[MET[[k]]]][L$is_module])) == 1,
                    sprintf("%s finite on every Module cell", k),
                    sprintf("%d missing", sum(!is.finite(L[[MET[[k]]]][L$is_module]))))
        p12_say("")

        # == Per arm =========================================================
        p12_say(sprintf("SECTION B -- post-IP minus pre-IP, per donor (%d cells per group, %d donors, B = %d)",
                        MINC, MIND, B))
        rows <- list()
        for (a in P12_ARMS) {
            m <- L[L$arm == a & L$is_module, ]
            m$g <- ifelse(m$group == "Module post-IP", "post", "pre")
            m$g <- factor(m$g, levels = c("post", "pre"))
            cnt <- table(m$donor_key, m$g)
            ok <- rownames(cnt)[cnt[, "post"] >= MINC & cnt[, "pre"] >= MINC]
            p12_say(sprintf("  %-6s donors with %d or more cells in both groups: %d", a, MINC,
                            length(ok)))
            for (k in names(MET)) {
                v <- m[[MET[[k]]]]
                med <- tapply(v, list(m$donor_key, m$g), stats::median)
                med <- med[, c("post", "pre"), drop = FALSE]
                for (g in c("post", "pre"))
                    rows[[length(rows) + 1L]] <- data.frame(
                        panel = "S10", arm = a, level = "pooled region sets",
                        quantity = paste0("donor_median_", g), measure = k, n_units = NA_integer_,
                        n_donors = length(ok),
                        estimate = if (length(ok)) stats::median(med[ok, g]) else NA_real_,
                        lo = NA_real_, hi = NA_real_, p = NA_real_, fdr = NA_real_,
                        family = paste("S10", a, "group medians", sep = "|"))
                dd <- med[ok, "post"] - med[ok, "pre"]
                est <- lo <- hi <- p <- NA_real_
                if (length(ok) >= MIND) {
                    idx <- p12_boot_index(ok, STEP, B)
                    bs <- p12_row_median(matrix(dd[idx], nrow = B))
                    est <- stats::median(dd)
                    ci <- p12_boot_ci(bs); lo <- ci[["lo"]]; hi <- ci[["hi"]]
                    p <- p12_boot_p(bs)
                }
                rows[[length(rows) + 1L]] <- data.frame(
                    panel = "S10", arm = a, level = "pooled region sets",
                    quantity = "post_minus_pre", measure = k, n_units = NA_integer_,
                    n_donors = length(ok), estimate = est, lo = lo, hi = hi, p = p, fdr = NA_real_,
                    family = paste("S10", a, sep = "|"))
                p12_say(sprintf("    %-14s %+.4f [%+.4f, %+.4f] | p %s", k, est, lo, hi,
                                p12_fmt(p, p12_boot_floor(B))))
            }
        }
        X <- do.call(rbind, rows)
        for (f in unique(X$family)) { i <- X$family == f; X$fdr[i] <- p12_bh(X$p[i]) }
        p12_say("")

        # == Write ===========================================================
        p12_say("SECTION C -- write")
        X$run_id <- P12$run_id
        X$config_version <- P12_CFG$config_version
        X$labels_version <- P12_CFG$labels$labels_version
        rownames(X) <- NULL
        X <- X[, P12_CFG$stats_columns]
        if (P12$n[["FAIL"]] > 0L) p12_abort("checks failed; stats not written")
        p12_write_stats(X, F_S10)
        p12_chk(nrow(utils::read.csv(F_S10)) == nrow(X), "stats file written and read back",
                sprintf("%d rows", nrow(X)))
    }

    p12_run(main)
}

# =============================================================================
# RUN: sections in order; stops at the first that does not pass
# =============================================================================
local({
    # Clears every Figure 15 section folder first, so no stale output survives a failed run.
    for (d in list.dirs(file.path(P12_RERUN, FIG_DIR), recursive = FALSE)) p12_clean(d)
    for (code in c("07", "08", "09", "10", "11")) {
        v <- get(paste0("p12_section_", code))()
        if (!identical(v, "PASS"))
            stop(sprintf("section %s ended with verdict %s; later sections not run. See its log.",
                         code, v), call. = FALSE)
    }
    message("Figure 15 analysis: sections 07 to 11 all PASS. Next: run Part12_Fig15_Plots.r")
})
