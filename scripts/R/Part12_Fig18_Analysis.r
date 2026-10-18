# =============================================================================
# Part12_Fig18_Analysis.r
# SEA-AD CEP-IP/CMRF, Part 12, Figure 18 and Supplementary Figure 17 analysis
#   Sections 22 to 24: nuclear-encoded panels against Module post-IP% and %MT
#   (18A), %MT against Module membership in pooled cells and within donor
#   (18B, 18C), and donor and region requirements (18D, 18E). Stops at the
#   first section that does not pass. The figures are drawn by
#   Part12_Fig18_Plots.r.
#
# Usage (from the SEA-AD folder, or set SEAAD_ROOT), after steps 02, 04 and 05
# and Part12_Fig16_Analysis.r
#   source("Part12_rerun/Part12_Fig18_Analysis.r")
#
# Outputs (Part12_rerun/10_Fig18/<section>/): stats, plot data and logs.
# =============================================================================

source(file.path(Sys.getenv("SEAAD_ROOT", getwd()), "Part12_rerun", "Part12_03_Config.r"))

SCRIPT <- "Part12_Fig18_Analysis.r"
FIG_DIR <- "10_Fig18"   # parent of every Figure 18 section folder
if (!isTRUE(as.integer(P12_CFG$config_version) >= 16L))
    stop(sprintf("config is version %s; Figure 18 needs version 16 or later",
                 P12_CFG$config_version), call. = FALSE)
F18 <- P12_CFG$figure18
PW18 <- F18$power

# == Helpers shared by the sections ===========================================
# Stop unless this file is sourced from its copy in Part12_rerun.
local({
    src <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
    want <- file.path(P12_RERUN, SCRIPT)
    if (!is.null(src) && !identical(tolower(normalizePath(src, winslash = "/", mustWork = FALSE)),
                                    tolower(normalizePath(want, winslash = "/", mustWork = FALSE))))
        stop(sprintf("run this file from Part12_rerun: save it there, then source(\"%s\")", want),
             call. = FALSE)
})

p12_finalise <- function(X) {
    X$run_id <- rep(P12$run_id, nrow(X))
    X$config_version <- rep(P12_CFG$config_version, nrow(X))
    X$labels_version <- rep(P12_CFG$labels$labels_version, nrow(X))
    rownames(X) <- NULL
    X[, P12_CFG$stats_columns]
}
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
p12_cluster <- function(n_jobs) {
    n <- p12_n_workers(n_jobs)
    p12_note("workers", sprintf("%d of %d logical cores (PSOCK)", n,
                                parallel::detectCores(logical = TRUE)))
    if (n > 1L) parallel::makeCluster(n, type = "PSOCK") else NULL
}
p12_write_plot <- function(df, path) {
    df$run_id <- rep(P12$run_id, nrow(df))   # rep() so an empty table still writes
    p12_output(path)
    utils::write.csv(df, path, row.names = FALSE)
}
# AUC of %MT for Module membership by the rank-sum identity (midranks for ties).
p12_auc <- function(x, pos) {
    n1 <- sum(pos); n0 <- sum(!pos)
    if (n1 < 1L || n0 < 1L) return(NA_real_)
    r <- rank(x)
    (sum(r[pos]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
# Best cut by Youden's J (cells above it called Module) and kappa at that cut.
p12_best_cut <- function(x, pos) {
    n1 <- sum(pos); n0 <- sum(!pos)
    out <- c(cut = NA_real_, kappa = NA_real_, youden = NA_real_)
    if (n1 < 1L || n0 < 1L) return(out)
    o <- order(x); xs <- x[o]; ps <- pos[o]
    tp <- n1 - cumsum(ps); fp <- n0 - cumsum(!ps)
    last <- c(xs[-1] != xs[-length(xs)], TRUE)
    tp <- tp[last]; fp <- fp[last]; cuts <- xs[last]
    j <- tp / n1 + (1 - fp / n0) - 1
    i <- which.max(j)
    a <- tp[i]; b <- fp[i]; cc <- n0 - b; d <- n1 - a; n <- n1 + n0
    acc <- (a + cc) / n
    pe <- ((a + b) / n) * ((a + d) / n) + ((cc + d) / n) * ((b + cc) / n)
    c(cut = cuts[i], kappa = if (pe < 1) (acc - pe) / (1 - pe) else NA_real_, youden = j[i])
}
# Donor share of variance across units (the 16A statistic).
p12_share <- function(y, d) {
    m <- mean(y); sst <- sum((y - m)^2)
    if (!(sst > 0)) return(NA_real_)
    s <- rowsum(y, d, reorder = FALSE)[, 1]; n <- rowsum(rep(1, length(y)), d, reorder = FALSE)[, 1]
    (sum(s^2 / n) - length(y) * m^2) / sst
}
# Split-half rs by the 16B rule, median over S random splits; NA if fewer than
# `min_d` donors have 2 or more units. Base R only (runs on the workers).
p12_split_half <- function(y, d, S, min_d) {
    tb <- table(d); keep <- names(tb)[tb >= 2L]
    if (length(keep) < min_d) return(NA_real_)
    i <- d %in% keep; y <- y[i]; d <- match(d[i], keep)
    o <- order(d); y <- y[o]; d <- d[o]
    nd <- tabulate(d); k <- length(nd)
    start <- cumsum(c(0L, nd[-k])); pos <- seq_along(d) - start[d]; h <- nd %/% 2L
    rs <- rep(NA_real_, S)
    for (s in seq_len(S)) {
        oo <- order(d, stats::runif(length(y)))
        rk <- integer(length(y)); rk[oo] <- pos
        a <- rk <= h[d]
        ma <- rowsum(y[a], d[a], reorder = TRUE)[, 1] / h
        mb <- rowsum(y[!a], d[!a], reorder = TRUE)[, 1] / (nd - h)
        ra <- rank(ma); rb <- rank(mb)
        if (stats::sd(ra) > 0 && stats::sd(rb) > 0) rs[s] <- stats::cor(ra, rb)
    }
    if (all(is.na(rs))) NA_real_ else stats::median(rs, na.rm = TRUE)
}
# Rescales the donor deviation so the donor share becomes `target`.
p12_rescale <- function(y, d, target, observed) {
    if (!(is.finite(target) && is.finite(observed) && target > 0 && target < 1 &&
          observed > 0 && observed < 1)) return(y)
    k <- sqrt((target / (1 - target)) / (observed / (1 - observed)))
    g <- mean(y); dm <- stats::ave(y, d)
    g + k * (dm - g) + (y - dm)
}
# OLS p of the last column of X (intercept already in X).
p12_ols_p_last <- function(yv, X) {
    q <- qr(X)
    if (q$rank < ncol(X)) return(NA_real_)
    b <- qr.coef(q, yv); r <- yv - X %*% b
    df <- nrow(X) - ncol(X)
    if (df < 1L) return(NA_real_)
    s2 <- sum(r^2) / df
    XtXi <- chol2inv(qr.R(q))
    se <- sqrt(s2 * XtXi[ncol(X), ncol(X)])
    2 * stats::pt(-abs(b[ncol(X)] / se), df)
}

# 18A verdict; an identical copy is in Part12_Fig18_Plots.r.
VK18 <- c("second only", "Both", "NE panel only", "Neither", "Unreadable", "Untestable")
p12_18_verdict <- function(fdr_s, fdr_n, vif, fitted) {
    alpha <- as.numeric(P12_CFG$inference$alpha); vr <- as.numeric(P12_CFG$floors$vif_readable)
    ss <- is.finite(fdr_s) & fdr_s < alpha; sn <- is.finite(fdr_n) & fdr_n < alpha
    ifelse(!fitted, "Untestable",
           ifelse(!(is.finite(vif) & vif < vr), "Unreadable",
                  ifelse(ss & !sn, "second only",
                         ifelse(!ss & sn, "NE panel only",
                                ifelse(ss & sn, "Both", "Neither")))))
}

# =============================================================================
# SECTION 22: 18A
# Donor panel scores (median over cells of the mean log-normalised count of
# the panel's detected genes); outcome ~ cov4 + second predictor (Module
# post-IP% or %MT) + panel score, standardised, OLS; BH over the outcomes.
# =============================================================================
p12_section_22 <- function() {
    STEP <- "22"
    D <- p12_step_dirs(STEP, "18A_NuclearPanels", parent = FIG_DIR)
    p12_clean(D$out)
    p12_start(STEP, SCRIPT, D)
    F_S <- file.path(D$out, "Stats_18A.csv")
    F_V <- file.path(D$out, "Plot_18A_verdicts.csv")
    F_P <- file.path(D$out, "Plot_18A_panels.csv")

    fit_18a <- function(d, y, covs, sec, min_n) {
        cols <- unique(c(y, sec, ".panel", covs))
        d <- d[stats::complete.cases(d[, cols, drop = FALSE]), , drop = FALSE]
        out <- list(n = nrow(d), fitted = FALSE)
        if (nrow(d) < min_n) return(out)
        if (!isTRUE(stats::sd(d[[y]]) > 0) || !isTRUE(stats::sd(d[[sec]]) > 0) ||
            !isTRUE(stats::sd(d$.panel) > 0)) return(out)
        d$.y <- as.numeric(scale(d[[y]]))
        d$.s <- as.numeric(scale(d[[sec]])); d$.n <- as.numeric(scale(d$.panel))
        rhs <- sprintf("`%s`", covs)
        f <- stats::as.formula(paste(".y ~", paste(c(rhs, ".s", ".n"), collapse = " + ")))
        m <- stats::lm(f, data = d)
        cf <- summary(m)$coefficients
        if (!all(c(".s", ".n") %in% rownames(cf))) return(out)
        list(n = nrow(d), fitted = TRUE, b_s = cf[".s", 1], p_s = cf[".s", 4],
             b_n = cf[".n", 1], p_n = cf[".n", 4],
             vif = unname(p12_vif(cbind(d$.s, d$.n))[1]))
    }

    main <- function() {
        MIND <- as.integer(P12_CFG$floors$donors_per_test)
        MINC <- as.integer(P12_CFG$floors$cells_per_group_per_donor)
        DET <- as.numeric(P12_CFG$nuclear_panels$detect_min)
        MING <- as.integer(F18$min_panel_genes)
        CL <- P12_CFG$clinical
        COVS <- unlist(CL$cov4)
        OUT <- vapply(CL$outcomes, function(o) o$column, "")
        OCLASS <- vapply(CL$outcomes, function(o) o$class, "")
        PAN <- F18$panels
        PK <- vapply(PAN, function(p) p$key, "")
        SECOND <- c(module = "module_post_pct", mito = "mito_log1p")
        p12_reading("Figure 18 definitions,")
        p12_reading("18A")

        p12_say("SECTION A -- inputs")
        L <- utils::read.csv(p12_upstream_file("04_Labels", "04_cell_labels.csv.gz"), check.names = FALSE,
                             colClasses = c(barcode = "character", donor_key = "character"))
        X <- methods::as(Matrix::readMM(p12_upstream_file("02_Extract", "02_counts.mtx.gz")),
                         "CsparseMatrix")
        genes <- readLines(p12_upstream_file("02_Extract", "02_counts_genes.tsv"))
        cells <- readLines(p12_upstream_file("02_Extract", "02_counts_cells.tsv"))
        SETS <- utils::read.csv(p12_upstream_file("02_Extract", "02_matched_sets.csv"),
                                stringsAsFactors = FALSE)
        DN <- utils::read.csv(p12_upstream_file("05_Units", "05_donors.csv"), check.names = FALSE,
                              colClasses = c(donor_key = "character", donor_id = "character"))
        p12_chk(all(L$labels_version == P12_CFG$labels$labels_version) &&
                    all(DN$labels_version == P12_CFG$labels$labels_version),
                "labels and donor table carry this labels version", "checked")
        p12_chk(nrow(X) == length(cells) && ncol(X) == length(genes),
                "count matrix, gene list and cell list agree", sprintf("%d x %d", nrow(X), ncol(X)))
        j0 <- match(L$barcode, cells)
        p12_chk(!anyNA(j0), "every labelled cell is in the count matrix", sprintf("%d cells", nrow(L)))
        miss <- setdiff(c(COVS, OUT, SECOND), names(DN))
        if (length(miss)) p12_abort(paste("donor table lacks:", paste(miss, collapse = ", ")))
        sets_needed <- vapply(PAN, function(p) p$set, "")
        for (a in P12_ARMS) {
            st <- P12_CFG$arms[[a]]$supertype
            have <- unique(SETS$set[SETS$supertype == st])
            p12_chk(all(sets_needed %in% have), sprintf("%s: every panel set is in 02_matched_sets.csv", a),
                    if (all(sets_needed %in% have)) sprintf("%d sets for %s", length(sets_needed), st) else
                        paste("missing:", paste(setdiff(sets_needed, have), collapse = ", ")))
        }
        if (P12$n[["FAIL"]] > 0L) p12_abort("input checks failed; nothing computed")
        p12_say("")

        p12_say("SECTION B -- donor panel scores")
        pan_rows <- list(); scores <- list()
        for (a in P12_ARMS) for (ds in P12_DATASETS) {
            st <- P12_CFG$arms[[a]]$supertype
            r <- which(L$arm == a & L$dataset == ds)
            j <- j0[r]; tot <- L$total_umis[r]; don <- L$donor_key[r]
            if (!all(is.finite(tot) & tot > 0)) p12_abort(sprintf("%s %s: a cell has no total UMIs", a, ds))
            ncell <- table(don); okd <- names(ncell)[ncell >= MINC]
            for (k in seq_along(PAN)) {
                g <- unique(SETS$gene[SETS$supertype == st & SETS$set == PAN[[k]]$set])
                gi <- match(g, genes)
                M <- X[j, gi[!is.na(gi)], drop = FALSE]
                det <- Matrix::colMeans(M > 0) >= DET
                ng <- sum(det)
                pan_rows[[length(pan_rows) + 1L]] <- data.frame(
                    arm = a, dataset = ds, panel = PK[k], n_genes = length(g),
                    n_in_matrix = sum(!is.na(gi)), n_detected = ng, n_donors = length(okd))
                key <- paste(a, ds, PK[k], sep = "|")
                if (ng < MING || !length(okd)) { scores[[key]] <- NULL; next }
                Md <- M[, det, drop = FALSE]
                Md <- Matrix::Diagonal(x = 1e4 / tot) %*% Md
                Md@x <- log1p(Md@x)
                v <- Matrix::rowMeans(Md)
                sc <- tapply(v[don %in% okd], don[don %in% okd], stats::median)
                scores[[key]] <- sc
            }
            p12_say(sprintf("  %-6s %-8s cells %5d | donors with %d or more cells: %3d of %3d", a, ds,
                            length(r), MINC, length(okd), length(ncell)))
        }
        PR <- do.call(rbind, pan_rows)
        low <- PR[PR$n_detected < MING, ]
        p12_note("panels below the gene floor",
                 if (nrow(low)) paste(sprintf("%s %s %s (%d)", low$arm, low$dataset, low$panel, low$n_detected),
                                      collapse = "; ") else "none")
        p12_chk(all(PR$n_in_matrix == PR$n_genes), "every panel gene is in the count matrix",
                sprintf("%d panel rows", nrow(PR)))
        p12_say("")

        p12_say("SECTION C -- models")
        rows <- list(); vr <- list()
        for (a in P12_ARMS) for (ds in P12_DATASETS) {
            d0 <- DN[DN$arm == a & DN$dataset == ds, ]
            for (k in seq_along(PAN)) {
                sc <- scores[[paste(a, ds, PK[k], sep = "|")]]
                d0$.panel <- if (is.null(sc)) NA_real_ else as.numeric(sc[d0$donor_key])
                for (sp in names(SECOND)) for (o in names(OUT)) {
                    y <- OUT[[o]]
                    cv <- if (y %in% COVS) setdiff(COVS, y) else COVS
                    f <- if (is.null(sc)) list(n = 0L, fitted = FALSE) else fit_18a(d0, y, cv, SECOND[[sp]], MIND)
                    g <- function(z) if (isTRUE(f$fitted)) f[[z]] else NA_real_
                    fam <- paste("18A", a, ds, PK[k], sp, sep = "|")
                    ms <- paste(sp, PK[k], o, sep = "|")
                    rows[[length(rows) + 1L]] <- p12_row("18A", a, ds, "beta_second", ms, NA, f$n, g("b_s"),
                                                         p = g("p_s"), family = paste(fam, "second", sep = "|"))
                    rows[[length(rows) + 1L]] <- p12_row("18A", a, ds, "beta_panel", ms, NA, f$n, g("b_n"),
                                                         p = g("p_n"), family = paste(fam, "panel", sep = "|"))
                    rows[[length(rows) + 1L]] <- p12_row("18A", a, ds, "vif", ms, NA, f$n, g("vif"),
                                                         family = "18A|descriptive")
                    vr[[length(vr) + 1L]] <- data.frame(arm = a, dataset = ds, predictor = sp, panel = PK[k],
                                                        outcome = o, class = OCLASS[[o]], n_donors = f$n,
                                                        fitted = isTRUE(f$fitted), vif = g("vif"),
                                                        beta_second = g("b_s"), beta_panel = g("b_n"),
                                                        stringsAsFactors = FALSE)
                }
            }
        }
        S <- p12_bh_by_family(do.call(rbind, rows))
        V <- do.call(rbind, vr)
        key <- function(q) match(paste(V$arm, V$dataset, paste(V$predictor, V$panel, V$outcome, sep = "|"), q),
                                 paste(S$arm, S$level, S$measure, S$quantity))
        V$fdr_second <- S$fdr[key("beta_second")]; V$fdr_panel <- S$fdr[key("beta_panel")]
        V$verdict <- p12_18_verdict(V$fdr_second, V$fdr_panel, V$vif, V$fitted)
        p12_say("  M4 primary rows (cognition, both region sets), verdict counts per panel:")
        for (sp in names(SECOND)) {
            p12_say(sprintf("   second predictor: %s", sp))
            for (k in PK) {
                g <- V[V$arm == "M4" & V$predictor == sp & V$panel == k & V$class == "cognition", ]
                tb <- table(factor(g$verdict, VK18))
                p12_say(sprintf("    %-18s %s", k, paste(sprintf("%s %d", names(tb), tb), collapse = ", ")))
            }
        }
        for (a in P12_ARMS) for (sp in names(SECOND)) {
            g <- V[V$arm == a & V$predictor == sp, ]
            tb <- table(factor(g$verdict, VK18))
            p12_say(sprintf("  %-6s %-6s all rows: %s | max VIF %.2f", a, sp,
                            paste(sprintf("%s %d", names(tb), tb), collapse = ", "),
                            suppressWarnings(max(g$vif, na.rm = TRUE))))
        }
        p12_say("")

        p12_say("SECTION D -- write")
        if (P12$n[["FAIL"]] > 0L) p12_abort("checks failed; stats not written")
        p12_write_stats(p12_finalise(S), F_S)
        p12_write_plot(V, F_V)
        p12_write_plot(PR, F_P)
        p12_chk(TRUE, "stats and plot files written", sprintf("18A %d rows | %d verdicts", nrow(S), nrow(V)))
    }

    p12_run(main)
}

# =============================================================================
# SECTION 23: 18B and 18C
# 18B: AUC of %MT for Module membership, cells pooled, best cut and kappa.
# 18C: within-donor AUC, median over donors, donor bootstrap against 0.5.
# =============================================================================
p12_section_23 <- function() {
    STEP <- "23"
    D <- p12_step_dirs(STEP, "18BC_MTOverlap", parent = FIG_DIR)
    p12_clean(D$out)
    p12_start(STEP, SCRIPT, D)
    F_B <- file.path(D$out, "Stats_18B.csv")
    F_C <- file.path(D$out, "Stats_18C.csv")
    F_BC <- file.path(D$out, "Plot_18B_cells.csv")
    F_CD <- file.path(D$out, "Plot_18C_donors.csv")

    main <- function() {
        B <- as.integer(P12_CFG$inference$bootstrap_B)
        MIND <- as.integer(P12_CFG$floors$donors_per_test)
        MINC <- as.integer(P12_CFG$floors$cells_per_group_per_donor)
        BANDS <- as.numeric(F18$auc_bands)
        p12_reading("Figure 18 definitions,")
        p12_reading("18B")
        p12_reading("18C")
        p12_say("SECTION A -- inputs")
        L <- utils::read.csv(p12_upstream_file("04_Labels", "04_cell_labels.csv.gz"), check.names = FALSE,
                             colClasses = c(barcode = "character", donor_key = "character"))
        L <- L[, c("arm", "dataset", "donor_key", "is_module", "mito_pct", "labels_version")]
        p12_chk(all(L$labels_version == P12_CFG$labels$labels_version), "labels carry this labels version",
                sprintf("%d cells", nrow(L)))
        L$is_module <- as.logical(L$is_module)
        p12_chk(all(is.finite(L$mito_pct)) && !anyNA(L$is_module), "%MT and Module membership on every cell",
                "checked")
        p12_say("")

        p12_say("SECTION B -- 18B, cells pooled over donors")
        rb <- list()
        for (a in P12_ARMS) for (ds in P12_DATASETS) {
            d <- L[L$arm == a & L$dataset == ds, ]
            au <- p12_auc(d$mito_pct, d$is_module); bc <- p12_best_cut(d$mito_pct, d$is_module)
            band <- if (!is.finite(au)) "NA" else if (au < BANDS[1]) "does not recover" else
                if (au < BANDS[2]) "partial" else "close to a threshold"
            for (q in c("auc", "kappa", "best_cut", "youden"))
                rb[[length(rb) + 1L]] <- p12_row("18B", a, ds, q, "Module vs non-Module cells, %MT",
                                                 nrow(d), length(unique(d$donor_key)),
                                                 c(auc = au, kappa = bc[["kappa"]], best_cut = bc[["cut"]],
                                                   youden = bc[["youden"]])[[q]], family = "18B|descriptive")
            p12_say(sprintf("  %-6s %-8s Module %5d | non-Module %5d | AUC %.3f (%s) | kappa %.3f | best cut %.2f%%",
                            a, ds, sum(d$is_module), sum(!d$is_module), au, band, bc[["kappa"]], bc[["cut"]]))
        }
        SB <- do.call(rbind, rb)
        p12_say("")

        p12_say(sprintf("SECTION C -- 18C, within donor (%d or more cells per side)", MINC))
        pd <- list()
        for (a in P12_ARMS) for (ds in P12_DATASETS) {
            d <- L[L$arm == a & L$dataset == ds, ]
            for (u in unique(d$donor_key)) {
                z <- d[d$donor_key == u, ]
                n1 <- sum(z$is_module); n0 <- sum(!z$is_module)
                if (n1 < MINC || n0 < MINC) next
                pd[[length(pd) + 1L]] <- data.frame(arm = a, dataset = ds, donor_key = u, n_module = n1,
                                                    n_nonmodule = n0, auc = p12_auc(z$mito_pct, z$is_module))
            }
        }
        PD <- if (length(pd)) do.call(rbind, pd) else
            data.frame(arm = character(0), dataset = character(0), donor_key = character(0),
                       n_module = integer(0), n_nonmodule = integer(0), auc = numeric(0))
        rc <- list(); gi <- 0L
        for (ds in P12_DATASETS) for (a in P12_ARMS) {
            gi <- gi + 1L
            v <- PD$auc[PD$arm == a & PD$dataset == ds]
            est <- lo <- hi <- p <- NA_real_
            if (length(v) >= MIND) {
                est <- stats::median(v)
                set.seed(p12_seed(STEP) + 1000L * gi)
                idx <- matrix(sample.int(length(v), length(v) * B, replace = TRUE), nrow = B)
                bs <- apply(matrix(v[idx], nrow = B), 1L, stats::median) - 0.5
                ci <- p12_boot_ci(bs + 0.5); lo <- ci[["lo"]]; hi <- ci[["hi"]]; p <- p12_boot_p(bs)
            }
            rc[[length(rc) + 1L]] <- p12_row("18C", a, ds, "median_within_donor_auc", "AUC of %MT, within donor",
                                             sum(PD$n_module[PD$arm == a & PD$dataset == ds] +
                                                     PD$n_nonmodule[PD$arm == a & PD$dataset == ds]),
                                             length(v), est, lo, hi, p, paste("18C", ds, sep = "|"))
        }
        SC <- p12_bh_by_family(do.call(rbind, rc))
        for (i in seq_len(nrow(SC)))
            p12_say(sprintf("  %-6s %-8s donors %3d | median AUC %s | FDR %s", SC$arm[i], SC$level[i],
                            SC$n_donors[i],
                            if (is.finite(SC$estimate[i])) sprintf("%.3f [%.3f, %.3f]", SC$estimate[i], SC$lo[i],
                                                                   SC$hi[i]) else "untestable",
                            p12_fmt(SC$fdr[i], p12_boot_floor(B))))
        p12_say("")

        p12_say("SECTION D -- write")
        if (P12$n[["FAIL"]] > 0L) p12_abort("checks failed; stats not written")
        p12_write_stats(p12_finalise(SB), F_B)
        p12_write_stats(p12_finalise(SC), F_C)
        p12_write_plot(L[, c("arm", "dataset", "donor_key", "is_module", "mito_pct")], F_BC)
        p12_write_plot(PD, F_CD)
        p12_chk(TRUE, "stats and plot files written", sprintf("18B %d | 18C %d rows | %d donor AUCs",
                                                             nrow(SB), nrow(SC), nrow(PD)))
    }

    p12_run(main)
}

# =============================================================================
# SECTION 24: 18D and 18E
# 18D: donor-count simulation of rules R1 to R3 at rescaled donor shares;
# each simulation cell is seeded and run whole by one worker. 18E: donor
# share and split-half by minimum regions per donor. Closure: at 1.0x the
# donor share equals the 16A estimate.
# =============================================================================
p12_18d_cell <- function(spec, data, par) {
    RNGkind("L'Ecuyer-CMRG")
    set.seed(spec$seed)
    d <- data[[spec$arm]]
    y <- d$y_eff[[as.character(spec$effect)]]
    nd <- d$nd; n <- min(spec$n, nd)
    thr <- d$nullp[as.character(spec$n)]
    h1 <- h2 <- h3 <- t2 <- t3 <- 0L; nd_ok <- 0L
    for (b in seq_len(par$draws)) {
        pick <- sort(sample.int(nd, n))
        sel <- d$dcode %in% pick
        if (sum(sel) < 10L) next
        yy <- y[sel]; dd <- d$dcode[sel]
        if (length(unique(dd)) < 3L) next
        nd_ok <- nd_ok + 1L
        s <- p12_share(yy, dd)
        if (is.finite(s) && is.finite(thr) && s > thr) h1 <- h1 + 1L
        r2 <- p12_split_half(yy, dd, par$n_splits, par$min_d)
        if (is.finite(r2)) { t2 <- t2 + 1L; if (r2 >= par$split_rule) h2 <- h2 + 1L }
        if (!is.null(d$clin)) {
            dm <- tapply(yy, dd, mean)
            cl <- d$clin[match(as.integer(names(dm)), d$clin$dcode), , drop = FALSE]
            ok <- stats::complete.cases(cl[, c("out", par$covs), drop = FALSE])
            if (sum(ok) >= par$min_d) {
                t3 <- t3 + 1L
                z <- function(v) { s <- stats::sd(v); if (is.finite(s) && s > 0) (v - mean(v)) / s else v * NA }
                yv <- z(cl$out[ok]); xv <- z(as.numeric(dm)[ok])
                if (all(is.finite(yv)) && all(is.finite(xv))) {
                    # A covariate with no spread in this subsample is dropped.
                    cvm <- as.matrix(cl[ok, par$covs, drop = FALSE])
                    cvm <- cvm[, apply(cvm, 2L, function(v) isTRUE(stats::sd(v) > 0)), drop = FALSE]
                    Xm <- cbind(1, cvm, xv)
                    p <- p12_ols_p_last(yv, Xm)
                    if (is.finite(p) && p < par$alpha) h3 <- h3 + 1L
                }
            }
        }
    }
    data.frame(arm = spec$arm, effect = spec$effect, n_donors = n, target_share = d$obs * spec$effect,
               draws = nd_ok, r1_rate = if (nd_ok) h1 / nd_ok else NA_real_,
               r2_rate = if (nd_ok) h2 / nd_ok else NA_real_,
               r2_evaluable = if (nd_ok) t2 / nd_ok else NA_real_,
               r3_rate = if (nd_ok) h3 / nd_ok else NA_real_,
               r3_evaluable = if (nd_ok) t3 / nd_ok else NA_real_)
}

p12_section_24 <- function() {
    STEP <- "24"
    D <- p12_step_dirs(STEP, "18DE_DonorPower", parent = FIG_DIR)
    p12_clean(D$out)
    p12_start(STEP, SCRIPT, D, pkgs = "parallel")
    F_D <- file.path(D$out, "Stats_18D.csv")
    F_DC <- file.path(D$out, "Plot_18D_curves.csv")
    F_E <- file.path(D$out, "Stats_18E.csv")
    F_EC <- file.path(D$out, "Plot_18E_regions.csv")

    main <- function() {
        MIND <- as.integer(P12_CFG$floors$donors_per_test)
        CL <- P12_CFG$clinical
        COVS <- unlist(CL$cov4)
        OUTC <- CL$outcomes[[PW18$outcome]]$column
        EFF <- as.numeric(PW18$effects); GRID <- as.integer(PW18$n_grid)
        DETECT <- as.numeric(PW18$detect)
        p12_reading("Figure 18 definitions,")
        p12_reading("18D")
        p12_reading("18E")
        p12_reading("Supplementary Figure 17")

        p12_say("SECTION A -- inputs")
        U <- utils::read.csv(p12_upstream_file("05_Units", "05_units.csv"), check.names = FALSE,
                             colClasses = c(donor_key = "character", region = "character"))
        DN <- utils::read.csv(p12_upstream_file("05_Units", "05_donors.csv"), check.names = FALSE,
                              colClasses = c(donor_key = "character", donor_id = "character"))
        A16 <- utils::read.csv(p12_upstream_file(file.path("08_Fig16", "16AB_DonorStructure"), "Stats_16A.csv"),
                               stringsAsFactors = FALSE)
        p12_chk(all(U$labels_version == P12_CFG$labels$labels_version), "unit table carries this labels version",
                sprintf("%d units", nrow(U)))
        miss <- setdiff(c(COVS, OUTC), names(DN))
        if (length(miss)) p12_abort(paste("donor table lacks:", paste(miss, collapse = ", ")))
        p12_say("")

        p12_say("SECTION B -- arms, observed shares and the null per donor count")
        DATA <- list()
        for (a in P12_ARMS) {
            u <- U[U$arm == a & is.finite(U$module_post_pct), ]
            dons <- sort(unique(u$donor_key))
            dcode <- match(u$donor_key, dons); y <- u$module_post_pct
            obs <- p12_share(y, dcode)
            a16 <- A16$estimate[A16$arm == a & A16$quantity == "module_post_pct"]
            p12_chk(length(a16) == 1L && isTRUE(abs(a16 - obs) < 1e-9),
                    sprintf("%s donor share equals the 16A estimate", a),
                    sprintf("%.6f against 16A %.6f", obs, if (length(a16)) a16 else NA))
            # One row per donor; each field takes the donor's non-missing value,
            # and the two region-set rows must agree.
            fl <- c(OUTC, COVS)
            dv0 <- DN[DN$arm == a, c("donor_key", fl), drop = FALSE]
            conflict <- 0L
            dv <- do.call(rbind, lapply(split(dv0, dv0$donor_key), function(z) {
                out <- z[1L, , drop = FALSE]
                for (v in fl) {
                    w <- unique(z[[v]][!is.na(z[[v]])])
                    if (length(w) > 1L) conflict <<- conflict + 1L
                    out[[v]] <- if (length(w)) w[1L] else NA
                }
                out
            }))
            p12_chk(conflict == 0L, sprintf("%s cov4 and %s agree between a donor's region-set rows", a, OUTC),
                    sprintf("%d donors | %d conflicting fields", nrow(dv), conflict))
            clin <- data.frame(dcode = match(dv$donor_key, dons), out = dv[[OUTC]])
            for (cv in COVS) clin[[cv]] <- dv[[cv]]
            clin <- clin[!is.na(clin$dcode), ]
            # The null per donor count, from its own seed stream.
            set.seed(p12_seed(STEP) + 100L * match(a, P12_ARMS))
            nullp <- c()
            for (n in GRID[GRID <= length(dons)]) {
                vals <- c()
                for (i in seq_len(as.integer(PW18$null_draws))) {
                    pick <- sample.int(length(dons), n)
                    sel <- dcode %in% pick
                    yy <- y[sel]; dd <- dcode[sel]; rr <- u$region[sel]
                    if (length(unique(dd)) < 3L) next
                    nv <- vapply(seq_len(as.integer(PW18$null_perm)), function(k) {
                        lab <- dd
                        for (g in split(seq_along(rr), rr)) if (length(g) > 1L) lab[g] <- dd[g][sample.int(length(g))]
                        p12_share(yy, lab) }, 0)
                    nv <- nv[is.finite(nv)]
                    if (length(nv)) vals <- c(vals, stats::quantile(nv, 0.95, names = FALSE))
                }
                if (length(vals)) nullp[as.character(n)] <- stats::median(vals)
            }
            y_eff <- lapply(stats::setNames(EFF, as.character(EFF)), function(e) p12_rescale(y, dcode, obs * e, obs))
            sh <- vapply(EFF, function(e) p12_share(y_eff[[as.character(e)]], dcode), 0)
            p12_chk(all(abs(sh - obs * EFF) < 1e-9), sprintf("%s rescaled donor shares equal their targets", a),
                    paste(sprintf("%.2fx %.4f", EFF, sh), collapse = " | "))
            DATA[[a]] <- list(dcode = dcode, nd = length(dons), obs = obs, nullp = nullp, y_eff = y_eff,
                              clin = clin, region = u$region, y = y)
            p12_say(sprintf("  %-6s units %4d | donors %3d | observed share %.3f | null p95 at n: %s", a, nrow(u),
                            length(dons), obs, paste(sprintf("%s %.3f", names(nullp), nullp), collapse = ", ")))
        }
        if (P12$n[["FAIL"]] > 0L) p12_abort("closure against 16A failed; nothing simulated")
        p12_say("")

        p12_say(sprintf("SECTION C -- 18D, %d draws per point", as.integer(PW18$draws)))
        specs <- list(); ci <- 0L
        for (a in P12_ARMS) for (e in EFF) for (n in GRID[GRID <= DATA[[a]]$nd]) {
            ci <- ci + 1L
            specs[[ci]] <- list(arm = a, effect = e, n = n, seed = p12_seed(STEP) + 10000L + ci)
        }
        par18 <- list(draws = as.integer(PW18$draws), n_splits = as.integer(PW18$n_splits), min_d = MIND,
                      split_rule = as.numeric(PW18$split_rule), covs = COVS,
                      alpha = as.numeric(P12_CFG$inference$alpha))
        DATA_W <- lapply(DATA, function(d) d[c("dcode", "nd", "obs", "nullp", "y_eff", "clin")])
        cl <- p12_cluster(length(specs))
        if (!is.null(cl)) {
            on.exit(parallel::stopCluster(cl), add = TRUE)
            parallel::clusterExport(cl, c("p12_share", "p12_split_half", "p12_ols_p_last"), envir = globalenv())
        }
        t0 <- Sys.time()
        res <- if (is.null(cl)) lapply(specs, p12_18d_cell, data = DATA_W, par = par18) else
            parallel::parLapplyLB(cl, specs, p12_18d_cell, data = DATA_W, par = par18)
        CV <- do.call(rbind, res)
        p12_note("simulation time", sprintf("%.1f min | %d cells", as.numeric(difftime(Sys.time(), t0,
                                                                                    units = "mins")), nrow(CV)))
        cross_at <- function(d, col) {
            d <- d[order(d$n_donors), ]; ok <- is.finite(d[[col]]) & d[[col]] >= DETECT
            if (any(ok)) d$n_donors[which(ok)[1]] else NA_real_
        }
        rows <- list()
        for (a in P12_ARMS) for (e in EFF) {
            d <- CV[CV$arm == a & abs(CV$effect - e) < 1e-9, ]
            for (rr in c("r1_rate", "r2_rate", "r3_rate")) {
                tag <- c(r1_rate = "R1", r2_rate = "R2", r3_rate = "R3")[[rr]]
                rows[[length(rows) + 1L]] <- p12_row("18D", a, sprintf("%.2fx", e), "n_at_80", tag, NA,
                                                     max(d$n_donors), cross_at(d, rr), family = "18D|descriptive")
                rows[[length(rows) + 1L]] <- p12_row("18D", a, sprintf("%.2fx", e), "rate_at_max_n", tag, NA,
                                                     max(d$n_donors), d[[rr]][which.max(d$n_donors)],
                                                     family = "18D|descriptive")
            }
            x <- rows[length(rows) - c(5, 3, 1)]
            p12_say(sprintf("  %-6s %.2fx target share %.3f | n at %.0f%%: R1 %s, R2 %s, R3 %s", a, e,
                            d$target_share[1], 100 * DETECT,
                            ifelse(is.na(x[[1]]$estimate), "not reached", x[[1]]$estimate),
                            ifelse(is.na(x[[2]]$estimate), "not reached", x[[2]]$estimate),
                            ifelse(is.na(x[[3]]$estimate), "not reached", x[[3]]$estimate)))
        }
        SD <- do.call(rbind, rows)
        hk <- SD[SD$arm == "HKG8" & SD$quantity == "n_at_80" & SD$measure == "R1", ]
        p12_note("validity check: HKG8 R1 reaches 80% at no donor count",
                 if (all(is.na(hk$estimate))) "holds" else
                     sprintf("does not hold: %s", paste(sprintf("%s at %d", hk$level, hk$estimate)[!is.na(hk$estimate)],
                                                         collapse = ", ")))
        p12_say("")

        p12_say("SECTION D -- 18E, minimum regions per donor")
        er <- list()
        for (a in P12_ARMS) {
            d <- DATA[[a]]
            nreg <- tabulate(d$dcode, d$nd)
            top <- max(nreg)
            for (k in seq(3L, max(3L, top))) {
                keep <- which(nreg >= k)
                sel <- d$dcode %in% keep
                yy <- d$y[sel]; dd <- d$dcode[sel]
                share <- if (length(keep) >= as.integer(PW18$share_floor_donors)) p12_share(yy, dd) else NA_real_
                set.seed(p12_seed(STEP) + 50000L + 100L * match(a, P12_ARMS) + k)
                sh <- p12_split_half(yy, dd, as.integer(PW18$n_splits), MIND)
                er[[length(er) + 1L]] <- data.frame(arm = a, min_regions = k, donors = length(keep),
                                                    donors_all = d$nd, donors_frac = length(keep) / d$nd,
                                                    units = sum(sel), donor_share = share, split_half = sh)
            }
        }
        ER <- do.call(rbind, er)
        re <- list()
        for (i in seq_len(nrow(ER))) for (q in c("donors_frac", "donor_share", "split_half"))
            re[[length(re) + 1L]] <- p12_row("18E", ER$arm[i], sprintf("%d+ regions", ER$min_regions[i]), q,
                                             "minimum regions per donor", ER$units[i], ER$donors[i], ER[[q]][i],
                                             family = "18E|descriptive")
        SE <- do.call(rbind, re)
        MARG <- as.numeric(PW18$rise_margin)
        for (a in P12_ARMS) {
            z <- ER[ER$arm == a, ]
            p12_say(sprintf("  %-6s %s", a, paste(sprintf("%d+: kept %.2f share %s split %s", z$min_regions,
                                                          z$donors_frac, ifelse(is.finite(z$donor_share),
                                                                                sprintf("%.2f", z$donor_share), "NA"),
                                                          ifelse(is.finite(z$split_half),
                                                                 sprintf("%.2f", z$split_half), "NA")),
                                                  collapse = " | ")))
            gain <- function(v) { v <- v[is.finite(v)]; if (length(v) < 2L) 0 else max(v[-1] - v[1]) }
            p12_note(sprintf("%s rise over 3+ regions", a),
                     sprintf("donor share %+.3f, split-half %+.3f (margin %.2f)", gain(z$donor_share),
                             gain(z$split_half), MARG))
        }
        p12_say("")

        p12_say("SECTION E -- write")
        if (P12$n[["FAIL"]] > 0L) p12_abort("checks failed; stats not written")
        p12_write_stats(p12_finalise(SD), F_D)
        p12_write_stats(p12_finalise(SE), F_E)
        p12_write_plot(CV, F_DC)
        p12_write_plot(ER, F_EC)
        p12_chk(TRUE, "stats and plot files written", sprintf("18D %d | 18E %d rows", nrow(SD), nrow(SE)))
    }

    p12_run(main)
}

# =============================================================================
# RUN
# Sections in order; stop at the first that does not PASS.
# =============================================================================
## RUN, in order; stop at the first section that does not PASS
################################################################################
local({
    # Clear every Figure 18 section folder first, so no output from an earlier
    # run survives.
    for (d in list.dirs(file.path(P12_RERUN, FIG_DIR), recursive = FALSE)) p12_clean(d)
    for (code in c("22", "23", "24")) {
        v <- get(paste0("p12_section_", code))()
        if (!identical(v, "PASS"))
            stop(sprintf("section %s ended with verdict %s; later sections not run. See its log.",
                         code, v), call. = FALSE)
    }
    message("Figure 18 analysis: sections 22 to 24 all PASS. Next: run Part12_Fig18_Plots.r")
})
