# =============================================================================
# Part12_Fig17_Analysis.r
# SEA-AD CEP-IP/CMRF, Part 12, Figure 17 analysis (Figure 17, Supplementary Figure 12)
#   Sections 19 and 20: histology and CPS against five unit quantities (17A,
#   17B, S12), and Module post-IP% with %MT in one model plus their
#   collinearity (17C, 17D, 17E). Stops at the first section that does not
#   pass. The figures are drawn by Part12_Fig17_Plots.r.
#
# Usage (from the SEA-AD folder, or set SEAAD_ROOT), after steps 02 to 05
#   source("Part12_rerun/Part12_Fig17_Analysis.r")
#
# Outputs (Part12_rerun/09_Fig17/<section>/): stats, plot data and logs.
# =============================================================================

source(file.path(Sys.getenv("SEAAD_ROOT", getwd()), "Part12_rerun", "Part12_03_Config.r"))

SCRIPT <- "Part12_Fig17_Analysis.r"
FIG_DIR <- "09_Fig17"   # parent of every Figure 17 section folder
if (!isTRUE(as.integer(P12_CFG$config_version) >= 14L))
    stop(sprintf("config is version %s; Figure 17 needs version 14 or later",
                 P12_CFG$config_version), call. = FALSE)

# == Helpers shared by the sections ===========================================
# Spearman from ranks; NA when either side has no spread.
p12_rs <- function(x, y) {
    if (length(x) < 3L) return(NA_real_)
    rx <- rank(x); ry <- rank(y)
    sx <- stats::sd(rx); sy <- stats::sd(ry)
    if (!is.finite(sx) || !is.finite(sy) || sx == 0 || sy == 0) return(NA_real_)
    stats::cor(rx, ry)
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

# One statistic on rows r of a task (rs, abs_rs, abs_r, or duel: the
# standardised rank coefficients b_x, b_z). Same function for the observed
# value and each resample; NA when fewer than min_d donors. Base R only.
BOOT_MIN_D <- 3L
p12_17_stat <- function(t, r, min_d = t$min_d) {
    k <- if (t$kind == "duel") 2L else 1L
    na <- rep(NA_real_, k)
    if (length(r) < 3L || length(unique(t$d[r])) < min_d) return(na)
    if (t$kind == "abs_r") { x <- t$x[r]; y <- t$y[r] } else { x <- rank(t$x[r]); y <- rank(t$y[r]) }
    sx <- stats::sd(x); sy <- stats::sd(y)
    if (!is.finite(sx) || !is.finite(sy) || sx == 0 || sy == 0) return(na)
    if (t$kind == "rs") return(stats::cor(x, y))
    if (t$kind %in% c("abs_rs", "abs_r")) return(abs(stats::cor(x, y)))
    z <- rank(t$z[r])
    if (!isTRUE(stats::sd(z) > 0)) return(na)
    ryx <- stats::cor(y, x); ryz <- stats::cor(y, z); rxz <- stats::cor(x, z)
    den <- 1 - rxz^2
    if (!is.finite(den) || den < 1e-12) return(na)
    c((ryx - ryz * rxz) / den, (ryz - ryx * rxz) / den)
}

# Bootstrap on one chunk of resamples; each task owns column col (and col + 1 for a duel).
p12_17_chunk <- function(idx, tasks, nd, ncol_out, stat, min_d_boot) {
    out <- matrix(NA_real_, nrow(idx), ncol_out)
    for (b in seq_len(nrow(idx))) {
        w <- tabulate(idx[b, ], nbins = nd)
        for (t in tasks) {
            if (!t$testable) next
            r <- rep.int(seq_along(t$d), w[t$d])
            v <- stat(t, r, min_d_boot)
            out[b, t$col:(t$col + length(v) - 1L)] <- v
        }
    }
    out
}

# Donor bootstrap for every task, indices drawn once per arm and frame.
# Seed: step seed plus 1000 times the group's position. Returns B x ncol.
p12_17_boot <- function(tasks, group_of, dons_by_group, step, B, cl, ncol_out) {
    left <- setdiff(unique(group_of), names(dons_by_group))
    if (length(left)) p12_abort(paste("tasks without a donor group:", paste(left, collapse = ", ")))
    BS <- matrix(NA_real_, B, ncol_out)
    for (i in seq_along(dons_by_group)) {
        a <- names(dons_by_group)[i]
        k <- which(group_of == a)
        nd <- length(dons_by_group[[a]])
        if (!length(k) || nd < BOOT_MIN_D) next
        tk <- tasks[k]
        cols <- unlist(lapply(tk, function(t) t$col:(t$col + (if (t$kind == "duel") 1L else 0L))))
        loc <- 1L
        for (j in seq_along(tk)) {
            w <- if (tk[[j]]$kind == "duel") 2L else 1L
            tk[[j]]$col <- loc; loc <- loc + w
        }
        set.seed(p12_seed(step) + 1000L * i)
        idx <- matrix(sample.int(nd, nd * B, replace = TRUE), nrow = B, ncol = nd)
        res <- p12_boot_run(cl, idx, p12_17_chunk, tasks = tk, nd = nd, ncol_out = loc - 1L,
                            stat = p12_17_stat, min_d_boot = BOOT_MIN_D)
        BS[, cols] <- res
        p12_note(sprintf("%s bootstrap", a), sprintf("%d donors | %d columns", nd, length(cols)))
    }
    BS
}

# Builds one task at a level from unit-level vectors: donor means between
# donors, or donor-centred values within donor (donors with >= min_u units).
p12_17_task <- function(kind, vars, donor, level, dons, min_d, min_u) {
    ok <- Reduce(`&`, lapply(vars, is.finite))
    d <- donor[ok]
    v <- lapply(vars, function(z) z[ok])
    n_units <- sum(ok)
    if (level == "between donors") {
        if (length(d)) {
            v <- lapply(v, function(z) { m <- tapply(z, d, mean); m })
            dk <- names(v[[1]])
            v <- lapply(v, function(z) as.numeric(z[dk]))
            d <- dk
        }
    } else {
        v <- lapply(v, function(z) p12_centre_within(z, d, min_u))
        k <- Reduce(`&`, lapply(v, is.finite))
        v <- lapply(v, function(z) z[k]); d <- d[k]; n_units <- sum(k)
    }
    t <- c(list(kind = kind, d = match(d, dons), min_d = min_d), v)
    t$n_units <- n_units
    t$n_donors <- length(unique(d))
    t$testable <- t$n_donors >= min_d
    t
}

# Row median ignoring NA; NA where a row has no finite value.
p12_row_median_na <- function(M) {
    if (is.null(dim(M))) return(M)
    if (ncol(M) == 1L) return(M[, 1])
    apply(M, 1L, function(z) if (any(is.finite(z))) stats::median(z, na.rm = TRUE) else NA_real_)
}

# Verdict of a two-predictor row (17C, 17D). Copied in Part12_Fig17_Plots.r.
p12_17_verdict <- function(fdr_m, fdr_t, vif, fitted) {
    alpha <- as.numeric(P12_CFG$inference$alpha)
    vr <- as.numeric(P12_CFG$floors$vif_readable)
    sm <- is.finite(fdr_m) & fdr_m < alpha
    st <- is.finite(fdr_t) & fdr_t < alpha
    ifelse(!fitted, "Untestable",
           ifelse(!(is.finite(vif) & vif < vr), "Unreadable",
                  ifelse(sm & !st, "Module only",
                         ifelse(!sm & st, "%MT only",
                                ifelse(sm & st, "Both", "Neither")))))
}

# Histology region name to snRNA-seq region name, by the config map: NA for a
# region without sequencing, "" for an unmapped name.
p12_17_map_region <- function(x) {
    RM <- P12_CFG$region_map
    ns <- tolower(RM$no_sequencing)
    k <- tolower(trimws(as.character(x)))
    al <- unlist(RM$alias)
    out <- unname(al[k])
    out[is.na(out)] <- ""
    out[k %in% ns | tolower(out) %in% ns] <- NA_character_
    out
}

# The histology and CPS join onto the step 05 units, with input checks.
# Returns the joined unit table and the measure and CPS column names.
p12_17_join <- function(U) {
    HI <- P12_CFG$histology
    HC <- P12_CFG$inputs$neuropath_table
    CC <- P12_CFG$inputs$cps_workbook
    MEAS <- unlist(lapply(HI$families, function(f) f$measures), use.names = FALSE)
    CPS <- vapply(HI$cps, function(z) z$column, "")
    GLOBAL <- CPS[vapply(HI$cps, function(z) z$scope, "") == "global"]
    p12_chk(!anyDuplicated(MEAS) && length(MEAS) == 28L && !any(tolower(MEAS) %in% tolower(HI$excluded)),
            "the config lists 28 distinct stain measures, none also excluded",
            sprintf("%d listed", length(MEAS)))

    f_np <- p12_input("neuropath_table"); p12_log_input(f_np, "neuropathology table")
    f_cp <- p12_input("cps_workbook"); p12_log_input(f_cp, "CPS workbook")
    H <- p12_read_sheet(f_np, HC$sheet)
    C <- p12_read_sheet(f_cp, CC$sheet)
    p12_note("sheets", sprintf("neuropathology '%s' %d rows | CPS '%s' %d rows", HC$sheet,
                               nrow(H), CC$sheet, nrow(C)))

    # == Neuropathology ==
    hd <- p12_col(H, HC$donor, "neuropathology sheet")
    hr <- p12_col(H, HC$region, "neuropathology sheet")
    ha <- p12_col(H, HC$analysis_region, "neuropathology sheet")
    mcol <- vapply(MEAS, function(m) p12_col(H, m, "neuropathology sheet"), "")
    xcol <- vapply(HI$excluded, function(m) p12_col(H, m, "neuropathology sheet"), "")
    rest <- setdiff(names(H), c(hd, hr, ha, mcol, xcol))
    p12_chk(!length(rest), "every column of the sheet is a key, a used measure or an excluded measure",
            if (length(rest)) paste(rest, collapse = ", ") else
                sprintf("%d used, %d excluded", length(mcol), length(xcol)))
    av <- unique(as.character(H[[ha]]))
    p12_chk(length(av) == 1L, "the sheet has one analysis region",
            paste(av, collapse = ", "))
    reg <- p12_17_map_region(H[[hr]])
    bad <- unique(as.character(H[[hr]])[!is.na(reg) & reg == ""])
    p12_chk(!length(bad), "every neuropathology region is mapped or has no sequencing",
            if (length(bad)) paste("unmapped:", paste(bad, collapse = ", ")) else "all mapped")
    p12_note("neuropathology rows without sequencing, left out",
             sprintf("%d rows: %s", sum(is.na(reg)),
                     paste(sort(unique(as.character(H[[hr]])[is.na(reg)])), collapse = ", ")))
    keep <- !is.na(reg) & reg != ""
    HN <- data.frame(donor_key = p12_donor_key(H[[hd]][keep]), region = reg[keep],
                     stringsAsFactors = FALSE)
    nonnum <- 0L
    for (i in seq_along(MEAS)) {
        raw <- H[[mcol[i]]][keep]
        v <- p12_num(raw)
        txt <- toupper(trimws(as.character(raw)))
        nonnum <- nonnum + sum(is.na(v) & !is.nan(v) & !is.na(raw) & nzchar(txt) &
                                   !(txt %in% c("NA", "NAN")))
        HN[[MEAS[i]]] <- v
    }
    p12_chk(nonnum == 0L, "every stain value is numeric or empty",
            sprintf("%d non-numeric entries", nonnum))
    hk <- paste(HN$donor_key, HN$region, sep = "|")
    p12_chk(!anyDuplicated(hk), "one neuropathology row per donor and region",
            sprintf("%d rows | %d donors", nrow(HN), length(unique(HN$donor_key))))

    # == CPS ==
    cd <- p12_col(C, CC$donor, "CPS sheet")
    cr <- p12_col(C, CC$region, "CPS sheet")
    ccol <- vapply(CPS, function(m) p12_col(C, m, "CPS sheet"), "")
    creg <- p12_17_map_region(C[[cr]])
    bad <- unique(as.character(C[[cr]])[!is.na(creg) & creg == ""])
    p12_chk(!length(bad), "every CPS region is mapped or has no sequencing",
            if (length(bad)) paste("unmapped:", paste(bad, collapse = ", ")) else "all mapped")
    keep <- !is.na(creg) & creg != ""
    CN <- data.frame(donor_key = p12_donor_key(C[[cd]][keep]), region = creg[keep],
                     stringsAsFactors = FALSE)
    nonnum <- 0L
    for (i in seq_along(CPS)) {
        raw <- C[[ccol[i]]][keep]
        v <- p12_num(raw)
        txt <- toupper(trimws(as.character(raw)))
        nonnum <- nonnum + sum(is.na(v) & !is.nan(v) & !is.na(raw) & nzchar(txt) &
                                   !(txt %in% c("NA", "NAN")))
        CN[[CPS[i]]] <- v
    }
    p12_chk(nonnum == 0L, "every CPS value is numeric or empty",
            sprintf("%d non-numeric entries", nonnum))
    ck <- paste(CN$donor_key, CN$region, sep = "|")
    p12_chk(!anyDuplicated(ck), "one CPS row per donor and region",
            sprintf("%d rows | %d donors", nrow(CN), length(unique(CN$donor_key))))
    spread <- vapply(GLOBAL, function(g) {
        r <- tapply(CN[[g]], CN$donor_key, function(z) {
            z <- z[is.finite(z)]; if (length(z) > 1L) diff(range(z)) else 0 })
        max(c(0, r), na.rm = TRUE)
    }, 0)
    p12_chk(all(spread < 1e-9), "each global CPS score is constant within donor",
            paste(sprintf("%s max range %.2e", GLOBAL, spread), collapse = " | "))

    # == Join ==
    # Every unit-table donor should appear in both sheets.
    ud <- unique(U$donor_key)
    p12_chk(all(ud %in% HN$donor_key), "every unit-table donor has a neuropathology row",
            sprintf("%d of %d donors%s", sum(ud %in% HN$donor_key), length(ud),
                    if (all(ud %in% HN$donor_key)) "" else
                        paste0(" | missing: ", paste(utils::head(setdiff(ud, HN$donor_key), 8),
                                                     collapse = ", "))), warn_only = TRUE)
    p12_chk(all(ud %in% CN$donor_key), "every unit-table donor has a CPS row",
            sprintf("%d of %d donors%s", sum(ud %in% CN$donor_key), length(ud),
                    if (all(ud %in% CN$donor_key)) "" else
                        paste0(" | missing: ", paste(utils::head(setdiff(ud, CN$donor_key), 8),
                                                     collapse = ", "))), warn_only = TRUE)
    uk <- paste(U$donor_key, U$region, sep = "|")
    J <- U
    jh <- match(uk, hk); jc <- match(uk, ck)
    for (m in MEAS) J[[m]] <- HN[[m]][jh]
    for (m in CPS) J[[m]] <- CN[[m]][jc]
    J$joined_histology <- !is.na(jh)
    J$joined_cps <- !is.na(jc)
    for (a in P12_ARMS) {
        i <- J$arm == a
        miss <- sort(unique(J$region[i & !J$joined_histology]))
        p12_note(sprintf("%s join", a),
                 sprintf("%d units | histology %d | CPS %d | regions without histology: %s",
                         sum(i), sum(J$joined_histology[i]), sum(J$joined_cps[i]),
                         if (length(miss)) paste(miss, collapse = ", ") else "none"))
    }
    p12_chk(all(vapply(P12_ARMS, function(a) any(J$joined_histology[J$arm == a]), TRUE)) &&
                all(vapply(P12_ARMS, function(a) any(J$joined_cps[J$arm == a]), TRUE)),
            "every arm has joined units", "histology and CPS")
    list(J = J, meas = MEAS, cps = CPS, global = GLOBAL)
}

# Measure to stain family key, and family key to label and block. Copied in Part12_Fig17_Plots.r.
p12_17_families <- function() {
    F <- P12_CFG$histology$families
    list(of = unlist(lapply(F, function(f) stats::setNames(rep(f$key, length(f$measures)),
                                                           f$measures))),
         keys = vapply(F, function(f) f$key, ""),
         label = stats::setNames(vapply(F, function(f) f$label, ""), vapply(F, function(f) f$key, "")),
         block = stats::setNames(vapply(F, function(f) f$block, ""), vapply(F, function(f) f$key, "")),
         n = stats::setNames(vapply(F, function(f) length(f$measures), 0L),
                             vapply(F, function(f) f$key, "")))
}

LEVELS <- c("within donor", "between donors")
QS <- c("module_post_pct", "mito_log1p", "cyto_share", "contamination", "log_n_cells")

# =============================================================================
# SECTION 19: 17A, 17B and Supplementary Figure 12
# Per arm, level and quantity, Spearman rs of each unit quantity with the 28
# stain measures and the CPS scores; donor bootstrap. 17A: IBA1 measures.
# 17B: CPS scores. S12: per stain family median rs, and |Module rs| - |%MT rs|.
# =============================================================================
p12_section_19 <- function() {
    STEP <- "19"
    D <- p12_step_dirs(STEP, "17AB_Histology", parent = FIG_DIR)
    p12_clean(D$out)
    p12_start(STEP, SCRIPT, D, pkgs = "parallel")
    F_A <- file.path(D$out, "Stats_17A.csv")
    F_B <- file.path(D$out, "Stats_17B.csv")
    F_S <- file.path(D$out, "Stats_S12.csv")
    F_SM <- file.path(D$out, "Stats_S12_measures.csv")
    F_SC <- file.path(D$out, "Plot_S12_counts.csv")
    F_J <- file.path(D$out, "Joined_17AB_units.csv")

    main <- function() {
        B <- as.integer(P12_CFG$inference$bootstrap_B)
        MIND <- as.integer(P12_CFG$floors$donors_per_test)
        MINU <- as.integer(P12_CFG$floors$finite_units_within_donor)
        HI <- P12_CFG$histology
        IBA <- vapply(HI$iba1_17A, function(z) z$measure, "")
        SIGN <- stats::setNames(vapply(HI$iba1_17A, function(z) as.numeric(z$predicted_sign), 0), IBA)
        p12_reading("Definitions")
        p12_reading("17A")
        p12_reading("17B")
        p12_reading("Supplementary Figure 12")

        p12_say("SECTION A -- inputs and join")
        U <- utils::read.csv(p12_upstream_file("05_Units", "05_units.csv"), check.names = FALSE,
                             colClasses = c(donor_key = "character", region = "character"))
        p12_chk(all(U$labels_version == P12_CFG$labels$labels_version),
                "unit table carries this labels version", sprintf("%d units", nrow(U)))
        miss <- setdiff(QS, names(U))
        if (length(miss)) p12_abort(paste("unit table lacks:", paste(miss, collapse = ", ")))
        JN <- p12_17_join(U)
        if (P12$n[["FAIL"]] > 0L) p12_abort("input checks failed; nothing computed")
        J <- JN$J; MEAS <- JN$meas; CPS <- JN$cps; GLOBAL <- JN$global
        FAM <- p12_17_families()
        p12_chk(all(IBA %in% MEAS[FAM$of[MEAS] == "microglia_IBA1"]) && length(IBA) == 6L,
                "the 17A measures are the 6 IBA1 measures", paste(length(IBA), "listed"))
        # Donors of the arm's units joined to histology or CPS.
        DONS <- lapply(stats::setNames(P12_ARMS, P12_ARMS), function(a)
            sort(unique(J$donor_key[J$arm == a & (J$joined_histology | J$joined_cps)])))
        p12_say("")

        # == Tasks: one per arm, level, quantity and measure or score ==
        p12_say("SECTION B -- tasks")
        tasks <- list(); meta <- list(); col <- 1L
        add_task <- function(a, lv, q, m, set) {
            u <- J[J$arm == a, ]
            t <- p12_17_task("rs", list(x = u[[q]], y = u[[m]]), u$donor_key, lv, DONS[[a]], MIND, MINU)
            t$col <- col
            tasks[[length(tasks) + 1L]] <<- t
            meta[[length(meta) + 1L]] <<- data.frame(arm = a, level = lv, q = q, measure = m,
                                                     set = set, n_units = t$n_units,
                                                     n_donors = t$n_donors, col = col,
                                                     stringsAsFactors = FALSE)
            col <<- col + 1L
        }
        for (a in P12_ARMS) for (lv in LEVELS) for (q in QS) {
            for (m in MEAS) add_task(a, lv, q, m, "histology")
            for (m in if (lv == "within donor") setdiff(CPS, GLOBAL) else CPS) add_task(a, lv, q, m, "cps")
        }
        MT <- do.call(rbind, meta)
        MT$obs <- vapply(tasks, function(t) if (t$testable) p12_17_stat(t, seq_along(t$d)) else
            NA_real_, 0)
        p12_note("tasks", sprintf("%d (%d histology, %d CPS) | %d testable", nrow(MT),
                                  sum(MT$set == "histology"), sum(MT$set == "cps"),
                                  sum(is.finite(MT$obs))))
        p12_say("")

        p12_say(sprintf("SECTION C -- donor bootstrap, donors drawn per arm (B = %d)", B))
        cl <- p12_cluster(B)
        if (!is.null(cl)) on.exit(parallel::stopCluster(cl), add = TRUE)
        t0 <- Sys.time()
        BS <- p12_17_boot(tasks, MT$arm, DONS, STEP, B, cl, col - 1L)
        p12_note("bootstrap time", sprintf("%.1f min", as.numeric(difftime(Sys.time(), t0,
                                                                          units = "mins"))))
        p12_chk(nrow(BS) == B && ncol(BS) == nrow(MT), "bootstrap matrix complete",
                sprintf("%d x %d", nrow(BS), ncol(BS)))
        p12_say("")

        key <- function(a, lv, q, m) which(MT$arm == a & MT$level == lv & MT$q == q & MT$measure == m)
        one <- function(panel, a, lv, q, m, family) {
            i <- key(a, lv, q, m)
            ok <- is.finite(MT$obs[i])
            ci <- if (ok) p12_boot_ci(BS[, MT$col[i]]) else c(lo = NA_real_, hi = NA_real_)
            p12_row(panel, a, lv, q, m, MT$n_units[i], MT$n_donors[i], MT$obs[i], ci[["lo"]],
                    ci[["hi"]], if (ok && !is.na(family)) p12_boot_p(BS[, MT$col[i]]) else NA,
                    if (is.na(family)) "S12|per-measure|descriptive" else family)
        }

        # == 17A ==
        rows <- list()
        for (a in P12_ARMS) for (lv in LEVELS) for (q in QS) for (m in IBA)
            rows[[length(rows) + 1L]] <- one("17A", a, lv, q, m, paste("17A", a, lv, q, sep = "|"))
        A <- p12_bh_by_family(do.call(rbind, rows))

        # == 17B ==
        rows <- list()
        for (a in P12_ARMS) for (lv in LEVELS) for (q in QS)
            for (m in if (lv == "within donor") setdiff(CPS, GLOBAL) else CPS)
                rows[[length(rows) + 1L]] <- one("17B", a, lv, q, m, paste("17B", a, lv, q, sep = "|"))
        Bt <- p12_bh_by_family(do.call(rbind, rows))

        # == S12, per measure (no test) and per family ==
        rows <- list()
        for (a in P12_ARMS) for (lv in LEVELS) for (q in QS) for (m in MEAS)
            rows[[length(rows) + 1L]] <- one("S12", a, lv, q, m, NA)
        SM <- do.call(rbind, rows)
        rows <- list(); cnt <- list()
        for (a in P12_ARMS) for (lv in LEVELS) {
            for (q in QS) for (f in FAM$keys) {
                ms <- names(FAM$of)[FAM$of == f]
                ii <- vapply(ms, function(m) key(a, lv, q, m), 0L)
                tested <- ii[is.finite(MT$obs[ii])]
                cnt[[length(cnt) + 1L]] <- data.frame(arm = a, level = lv, quantity = q, family = f,
                                                      n_tested = length(tested), n_total = length(ms))
                est <- lo <- hi <- p <- NA_real_
                if (length(tested)) {
                    est <- stats::median(MT$obs[tested])
                    bs <- p12_row_median_na(BS[, MT$col[tested], drop = FALSE])
                    ci <- p12_boot_ci(bs); lo <- ci[["lo"]]; hi <- ci[["hi"]]; p <- p12_boot_p(bs)
                }
                rows[[length(rows) + 1L]] <- p12_row("S12", a, lv, q, f,
                                                     if (length(tested)) max(MT$n_units[tested]) else 0L,
                                                     if (length(tested)) max(MT$n_donors[tested]) else 0L,
                                                     est, lo, hi, p, paste("S12A", a, lv, q, sep = "|"))
            }
            for (f in FAM$keys) {
                ms <- names(FAM$of)[FAM$of == f]
                im <- vapply(ms, function(m) key(a, lv, "module_post_pct", m), 0L)
                it <- vapply(ms, function(m) key(a, lv, "mito_log1p", m), 0L)
                ok <- is.finite(MT$obs[im]) & is.finite(MT$obs[it])
                est <- lo <- hi <- p <- NA_real_
                if (any(ok)) {
                    est <- stats::median(abs(MT$obs[im[ok]]) - abs(MT$obs[it[ok]]))
                    bs <- p12_row_median_na(abs(BS[, MT$col[im[ok]], drop = FALSE]) -
                                                abs(BS[, MT$col[it[ok]], drop = FALSE]))
                    ci <- p12_boot_ci(bs); lo <- ci[["lo"]]; hi <- ci[["hi"]]; p <- p12_boot_p(bs)
                }
                rows[[length(rows) + 1L]] <- p12_row("S12", a, lv, "abs_module_minus_abs_mito", f,
                                                     if (any(ok)) max(MT$n_units[im[ok]]) else 0L,
                                                     if (any(ok)) max(MT$n_donors[im[ok]]) else 0L,
                                                     est, lo, hi, p, paste("S12B", a, lv, sep = "|"))
            }
        }
        S <- p12_bh_by_family(do.call(rbind, rows))
        SC <- do.call(rbind, cnt)

        # == Log: the readings' rows ==
        p12_say("SECTION D -- results against the readings")
        p12_say("  17A, within donor, Module post-IP% (predicted sign: + activated, - inactivated, - process area and length)")
        for (a in P12_ARMS) {
            g <- A[A$arm == a & A$level == "within donor" & A$quantity == "module_post_pct", ]
            hit <- 0L; opp <- 0L
            for (m in IBA) {
                r <- g[g$measure == m, ]; t <- A[A$arm == a & A$level == "within donor" & A$measure == m, ]
                sig <- is.finite(r$fdr) && r$fdr < 0.05
                dir <- if (!is.finite(r$estimate) || SIGN[[m]] == 0) "" else
                    if (sign(r$estimate) == SIGN[[m]]) "predicted" else "opposite"
                if (sig && dir == "predicted") hit <- hit + 1L
                if (sig && dir == "opposite") opp <- opp + 1L
                p12_say(sprintf("    %-6s %-48s Module %+.3f FDR %-10s %-9s | %%MT %+.3f | cyto %+.3f | donors %d",
                                a, m, r$estimate, p12_fmt(r$fdr, p12_boot_floor(B)), dir,
                                t$estimate[t$quantity == "mito_log1p"],
                                t$estimate[t$quantity == "cyto_share"], r$n_donors))
            }
            p12_say(sprintf("    %-6s shape measures significant in the predicted direction: %d of %d | opposite: %d",
                            a, hit, sum(SIGN != 0), opp))
        }
        p12_say("  17B, Module post-IP% against CPS")
        for (a in P12_ARMS) for (lv in LEVELS) {
            g <- Bt[Bt$arm == a & Bt$level == lv & Bt$quantity == "module_post_pct", ]
            p12_say(sprintf("    %-6s %-15s %s", a, lv,
                            paste(sprintf("%s %+.3f (%s)", g$measure, g$estimate,
                                          p12_fmt(g$fdr, p12_boot_floor(B))), collapse = " | ")))
        }
        p12_say("  S12B, |Module rs| - |%MT rs| per stain family")
        for (a in P12_ARMS) for (lv in LEVELS) {
            g <- S[S$arm == a & S$level == lv & S$quantity == "abs_module_minus_abs_mito", ]
            p12_say(sprintf("    %-6s %-15s %s", a, lv,
                            paste(sprintf("%s %+.3f%s", FAM$label[g$measure], g$estimate,
                                          p12_stars(g$fdr)), collapse = " | ")))
        }
        p12_say("")

        p12_say("SECTION E -- write")
        if (P12$n[["FAIL"]] > 0L) p12_abort("checks failed; stats not written")
        p12_write_stats(p12_finalise(A), F_A)
        p12_write_stats(p12_finalise(Bt), F_B)
        p12_write_stats(p12_finalise(S), F_S)
        p12_write_stats(p12_finalise(SM), F_SM)
        SC$run_id <- P12$run_id
        p12_output(F_SC); utils::write.csv(SC, F_SC, row.names = FALSE)
        keepc <- c("arm", "dataset", "donor_key", "region", QS, MEAS, CPS,
                   "joined_histology", "joined_cps")
        JW <- J[, keepc]; JW$run_id <- P12$run_id
        p12_chk(p12_ascii_ok(JW[, c("arm", "dataset", "donor_key", "region")]),
                "joined table text fields are ASCII", "checked")
        p12_output(F_J); utils::write.csv(JW, F_J, row.names = FALSE)
        p12_chk(TRUE, "stats and plot files written",
                sprintf("17A %d | 17B %d | S12 %d | S12 measures %d rows", nrow(A), nrow(Bt),
                        nrow(S), nrow(SM)))
    }

    p12_run(main)
}

# =============================================================================
# SECTION 20: 17C, 17D and 17E
# 17C: rank(measure) ~ rank(Module post-IP%) + rank(%MT) per arm, level and
# measure, donor bootstrap. 17D: the 15A cov4 clinical model with both
# predictors, OLS. 17E: |r| of Module post-IP% with %MT and its VIF.
# =============================================================================
p12_section_20 <- function() {
    STEP <- "20"
    D <- p12_step_dirs(STEP, "17CE_ModuleVsMT", parent = FIG_DIR)
    p12_clean(D$out)
    p12_start(STEP, SCRIPT, D, pkgs = "parallel")
    F_C <- file.path(D$out, "Stats_17C.csv")
    F_D <- file.path(D$out, "Stats_17D.csv")
    F_E <- file.path(D$out, "Stats_17E.csv")
    F_CV <- file.path(D$out, "Plot_17C_verdicts.csv")
    F_DV <- file.path(D$out, "Plot_17D_verdicts.csv")

    # One cov4 OLS model with both predictors: coefficients, CI, p, single-predictor fits and VIF.
    fit_duel <- function(d, y, covs, min_n) {
        cols <- unique(c(y, "module_post_pct", "mito_log1p", covs))
        d <- d[stats::complete.cases(d[, cols, drop = FALSE]), , drop = FALSE]
        out <- list(n = nrow(d), fitted = FALSE)
        if (nrow(d) < min_n) return(out)
        if (!isTRUE(stats::sd(d[[y]]) > 0) || !isTRUE(stats::sd(d$module_post_pct) > 0) ||
            !isTRUE(stats::sd(d$mito_log1p) > 0)) return(out)
        d$.y <- as.numeric(scale(d[[y]]))
        d$.m <- as.numeric(scale(d$module_post_pct)); d$.t <- as.numeric(scale(d$mito_log1p))
        rhs <- sprintf("`%s`", covs)
        base <- paste(".y ~", if (length(rhs)) paste(rhs, collapse = " + ") else "1")
        m2 <- stats::lm(stats::as.formula(paste(base, "+ .m + .t")), data = d)
        cf <- summary(m2)$coefficients
        if (!all(c(".m", ".t") %in% rownames(cf))) return(out)
        ci <- stats::confint(m2, c(".m", ".t"), level = 0.95)
        s_m <- stats::coef(stats::lm(stats::as.formula(paste(base, "+ .m")), data = d))[[".m"]]
        s_t <- stats::coef(stats::lm(stats::as.formula(paste(base, "+ .t")), data = d))[[".t"]]
        r12 <- stats::cor(d$.m, d$.t)
        list(n = nrow(d), fitted = TRUE,
             b_m = cf[".m", "Estimate"], lo_m = ci[".m", 1], hi_m = ci[".m", 2], p_m = cf[".m", "Pr(>|t|)"],
             b_t = cf[".t", "Estimate"], lo_t = ci[".t", 1], hi_t = ci[".t", 2], p_t = cf[".t", "Pr(>|t|)"],
             single_m = s_m, single_t = s_t, vif = unname(p12_vif(cbind(d$.m, d$.t))[1]))
    }

    main <- function() {
        B <- as.integer(P12_CFG$inference$bootstrap_B)
        MIND <- as.integer(P12_CFG$floors$donors_per_test)
        MINU <- as.integer(P12_CFG$floors$finite_units_within_donor)
        CL <- P12_CFG$clinical
        COVS <- unlist(CL$cov4)
        OUT <- vapply(CL$outcomes, function(o) o$column, "")
        OCLASS <- vapply(CL$outcomes, function(o) o$class, "")
        p12_reading("Definitions")
        p12_reading("17C")
        p12_reading("17D")
        p12_reading("17E")

        p12_say("SECTION A -- inputs and join")
        U <- utils::read.csv(p12_upstream_file("05_Units", "05_units.csv"), check.names = FALSE,
                             colClasses = c(donor_key = "character", region = "character"))
        DN <- utils::read.csv(p12_upstream_file("05_Units", "05_donors.csv"), check.names = FALSE,
                              colClasses = c(donor_key = "character", donor_id = "character"))
        p12_chk(all(U$labels_version == P12_CFG$labels$labels_version) &&
                    all(DN$labels_version == P12_CFG$labels$labels_version),
                "unit and donor tables carry this labels version", "checked")
        miss <- setdiff(c(COVS, OUT, "module_post_pct", "mito_log1p"), names(DN))
        if (length(miss)) p12_abort(paste("donor table lacks:", paste(miss, collapse = ", ")))
        JN <- p12_17_join(U)
        if (P12$n[["FAIL"]] > 0L) p12_abort("input checks failed; nothing computed")
        J <- JN$J; MEAS <- JN$meas
        FAM <- p12_17_families()
        # Donor groups: histology-joined units, and the donor table per region set.
        DONS <- list()
        for (a in P12_ARMS) {
            DONS[[paste(a, "HU", sep = "|")]] <-
                sort(unique(J$donor_key[J$arm == a & J$joined_histology]))
            for (ds in P12_DATASETS) {
                d <- DN[DN$arm == a & DN$dataset == ds, ]
                d <- d[stats::complete.cases(d[, COVS, drop = FALSE]) & is.finite(d$module_post_pct) &
                           is.finite(d$mito_log1p), , drop = FALSE]
                DONS[[paste(a, "CG", ds, sep = "|")]] <- sort(unique(d$donor_key))
            }
        }
        p12_say("")

        # == Tasks: 17C duels, 17E |r| ==
        p12_say("SECTION B -- tasks")
        tasks <- list(); meta <- list(); col <- 1L
        push <- function(t, panel, a, lv, m, grp) {
            t$col <- col
            tasks[[length(tasks) + 1L]] <<- t
            meta[[length(meta) + 1L]] <<- data.frame(panel = panel, arm = a, level = lv, measure = m,
                                                     group = grp,
                                                     n_units = t$n_units, n_donors = t$n_donors,
                                                     col = col, stringsAsFactors = FALSE)
            col <<- col + if (t$kind == "duel") 2L else 1L
        }
        for (a in P12_ARMS) for (lv in LEVELS) {
            u <- J[J$arm == a, ]
            for (m in MEAS)
                push(p12_17_task("duel", list(y = u[[m]], x = u$module_post_pct, z = u$mito_log1p),
                                 u$donor_key, lv, DONS[[paste(a, "HU", sep = "|")]], MIND, MINU),
                     "17C", a, lv, m, paste(a, "HU", sep = "|"))
        }
        for (a in P12_ARMS) {
            # Histology-unit rows use the units joined to a neuropathology row.
            u <- J[J$arm == a & J$joined_histology, ]
            push(p12_17_task("abs_rs", list(x = u$module_post_pct, y = u$mito_log1p), u$donor_key,
                             "within donor", DONS[[paste(a, "HU", sep = "|")]], MIND, MINU), "17E", a,
                 "within donor HU", "pooled", paste(a, "HU", sep = "|"))
            push(p12_17_task("abs_rs", list(x = u$module_post_pct, y = u$mito_log1p), u$donor_key,
                             "between donors", DONS[[paste(a, "HU", sep = "|")]], MIND, MINU), "17E", a,
                 "between donors HU", "pooled", paste(a, "HU", sep = "|"))
            for (ds in P12_DATASETS) {
                d <- DN[DN$arm == a & DN$dataset == ds, ]
                d <- d[stats::complete.cases(d[, COVS, drop = FALSE]), , drop = FALSE]
                ok <- is.finite(d$module_post_pct) & is.finite(d$mito_log1p)
                t <- list(kind = "abs_r", x = d$module_post_pct[ok], y = d$mito_log1p[ok],
                          d = match(d$donor_key[ok], DONS[[paste(a, "CG", ds, sep = "|")]]), min_d = MIND,
                          n_units = sum(ok),
                          n_donors = length(unique(d$donor_key[ok])))
                t$testable <- t$n_donors >= MIND
                push(t, "17E", a, "between donors CG", ds, paste(a, "CG", ds, sep = "|"))
            }
        }
        MT <- do.call(rbind, meta)
        OBS <- lapply(tasks, function(t) if (t$testable) p12_17_stat(t, seq_along(t$d)) else
            rep(NA_real_, if (t$kind == "duel") 2L else 1L))
        p12_note("tasks", sprintf("%d 17C models, %d 17E rows | testable %d",
                                  sum(MT$panel == "17C"), sum(MT$panel == "17E"),
                                  sum(vapply(OBS, function(v) all(is.finite(v)), TRUE))))
        p12_say("")

        p12_say(sprintf("SECTION C -- donor bootstrap, donors drawn per arm (B = %d)", B))
        cl <- p12_cluster(B)
        if (!is.null(cl)) on.exit(parallel::stopCluster(cl), add = TRUE)
        t0 <- Sys.time()
        BS <- p12_17_boot(tasks, MT$group, DONS, STEP, B, cl, col - 1L)
        p12_note("bootstrap time", sprintf("%.1f min", as.numeric(difftime(Sys.time(), t0,
                                                                          units = "mins"))))
        p12_chk(nrow(BS) == B && ncol(BS) == col - 1L, "bootstrap matrix complete",
                sprintf("%d x %d", nrow(BS), ncol(BS)))
        p12_say("")

        # == 17C ==
        p12_say("SECTION D -- 17C, Module post-IP% and %MT in one rank model, histology")
        rows <- list(); vr <- list()
        for (k in which(MT$panel == "17C")) {
            t <- tasks[[k]]; a <- MT$arm[k]; lv <- MT$level[k]; m <- MT$measure[k]
            ob <- OBS[[k]]; ok <- all(is.finite(ob))
            bm <- BS[, MT$col[k]]; bt <- BS[, MT$col[k] + 1L]
            cm <- if (ok) p12_boot_ci(bm) else c(lo = NA, hi = NA)
            ct <- if (ok) p12_boot_ci(bt) else c(lo = NA, hi = NA)
            rows[[length(rows) + 1L]] <- p12_row("17C", a, lv, "beta_module", m, MT$n_units[k],
                                                 MT$n_donors[k], ob[1], cm[["lo"]], cm[["hi"]],
                                                 if (ok) p12_boot_p(bm) else NA,
                                                 paste("17C", a, lv, "module", sep = "|"))
            rows[[length(rows) + 1L]] <- p12_row("17C", a, lv, "beta_mito", m, MT$n_units[k],
                                                 MT$n_donors[k], ob[2], ct[["lo"]], ct[["hi"]],
                                                 if (ok) p12_boot_p(bt) else NA,
                                                 paste("17C", a, lv, "mito", sep = "|"))
            r12 <- rsm <- rst <- NA_real_
            if (t$testable) {
                r12 <- p12_rs(t$x, t$z); rsm <- p12_rs(t$x, t$y); rst <- p12_rs(t$z, t$y)
            }
            vif <- if (is.finite(r12) && abs(r12) < 1) 1 / (1 - r12^2) else if (is.finite(r12)) Inf else NA
            for (qq in c("vif", "rs_module_marginal", "rs_mito_marginal"))
                rows[[length(rows) + 1L]] <- p12_row("17C", a, lv, qq, m, MT$n_units[k], MT$n_donors[k],
                                                     c(vif = vif, rs_module_marginal = rsm,
                                                       rs_mito_marginal = rst)[[qq]],
                                                     family = "17C|descriptive")
            vr[[length(vr) + 1L]] <- data.frame(arm = a, level = lv, measure = m,
                                                stain_family = FAM$of[[m]],
                                                block = FAM$block[[FAM$of[[m]]]],
                                                n_units = MT$n_units[k], n_donors = MT$n_donors[k],
                                                fitted = ok, vif = vif,
                                                beta_module = ob[1], beta_mito = ob[2],
                                                rs_module = rsm, rs_mito = rst,
                                                stringsAsFactors = FALSE)
        }
        C <- p12_bh_by_family(do.call(rbind, rows))
        CV <- do.call(rbind, vr)
        gf <- function(q) C$fdr[match(paste(CV$arm, CV$level, CV$measure, q),
                                      paste(C$arm, C$level, C$measure, C$quantity))]
        CV$fdr_module <- gf("beta_module"); CV$fdr_mito <- gf("beta_mito")
        CV$module_sign_flip <- is.finite(CV$beta_module) & is.finite(CV$rs_module) &
            CV$beta_module * CV$rs_module < 0
        CV$mito_sign_flip <- is.finite(CV$beta_mito) & is.finite(CV$rs_mito) &
            CV$beta_mito * CV$rs_mito < 0
        CV$verdict <- p12_17_verdict(CV$fdr_module, CV$fdr_mito, CV$vif, CV$fitted)
        VK <- c("Module only", "%MT only", "Both", "Neither", "Unreadable", "Untestable")
        for (a in P12_ARMS) for (lv in LEVELS) {
            g <- CV[CV$arm == a & CV$level == lv, ]
            tb <- table(factor(g$verdict, VK))
            p12_say(sprintf("  %-6s %-15s %s | %%MT only with a Module sign flip %d | donors %d",
                            a, lv, paste(sprintf("%s %d", VK, tb), collapse = ", "),
                            sum(g$verdict == "%MT only" & g$module_sign_flip),
                            max(g$n_donors)))
        }
        p12_say("")

        # == 17D ==
        p12_say("SECTION E -- 17D, Module post-IP% and %MT in one cov4 model, clinical")
        rows <- list(); vr <- list()
        for (a in P12_ARMS) for (ds in P12_DATASETS) {
            d <- DN[DN$arm == a & DN$dataset == ds, ]
            for (k in names(OUT)) {
                y <- OUT[[k]]
                cv <- if (y %in% COVS) setdiff(COVS, y) else COVS
                f <- fit_duel(d, y, cv, MIND)
                g <- function(z) if (isTRUE(f$fitted)) f[[z]] else NA_real_
                rows[[length(rows) + 1L]] <- p12_row("17D", a, ds, "beta_module", k, NA, f$n, g("b_m"),
                                                     g("lo_m"), g("hi_m"), g("p_m"),
                                                     paste("17D", a, ds, "module", sep = "|"))
                rows[[length(rows) + 1L]] <- p12_row("17D", a, ds, "beta_mito", k, NA, f$n, g("b_t"),
                                                     g("lo_t"), g("hi_t"), g("p_t"),
                                                     paste("17D", a, ds, "mito", sep = "|"))
                rows[[length(rows) + 1L]] <- p12_row("17D", a, ds, "vif", k, NA, f$n, g("vif"),
                                                     family = "17D|descriptive")
                vr[[length(vr) + 1L]] <- data.frame(arm = a, dataset = ds, outcome = k,
                                                    class = OCLASS[[k]], n_donors = f$n,
                                                    fitted = isTRUE(f$fitted), vif = g("vif"),
                                                    beta_module = g("b_m"), beta_mito = g("b_t"),
                                                    single_module = g("single_m"),
                                                    single_mito = g("single_t"),
                                                    stringsAsFactors = FALSE)
            }
        }
        Dd <- p12_bh_by_family(do.call(rbind, rows))
        DV <- do.call(rbind, vr)
        gd <- function(q) Dd$fdr[match(paste(DV$arm, DV$dataset, DV$outcome, q),
                                       paste(Dd$arm, Dd$level, Dd$measure, Dd$quantity))]
        DV$fdr_module <- gd("beta_module"); DV$fdr_mito <- gd("beta_mito")
        DV$module_sign_flip <- is.finite(DV$beta_module) & is.finite(DV$single_module) &
            DV$beta_module * DV$single_module < 0
        DV$mito_sign_flip <- is.finite(DV$beta_mito) & is.finite(DV$single_mito) &
            DV$beta_mito * DV$single_mito < 0
        DV$verdict <- p12_17_verdict(DV$fdr_module, DV$fdr_mito, DV$vif, DV$fitted)
        for (a in P12_ARMS) for (ds in P12_DATASETS) {
            g <- DV[DV$arm == a & DV$dataset == ds, ]
            p12_say(sprintf("  %-6s %-8s %s", a, ds,
                            paste(sprintf("%s: %s", g$outcome, g$verdict), collapse = " | ")))
        }
        p12_say("")

        # == 17E ==
        p12_say("SECTION F -- 17E, collinearity of Module post-IP% and %MT")
        rows <- list()
        for (k in which(MT$panel == "17E")) {
            ob <- OBS[[k]]; ok <- is.finite(ob)
            ci <- if (ok) p12_boot_ci(BS[, MT$col[k]]) else c(lo = NA, hi = NA)
            vif <- if (ok && ob < 1) 1 / (1 - ob^2) else if (ok) Inf else NA_real_
            rows[[length(rows) + 1L]] <- p12_row("17E", MT$arm[k], MT$level[k], "abs_r", MT$measure[k],
                                                 MT$n_units[k], MT$n_donors[k], ob, ci[["lo"]],
                                                 ci[["hi"]], family = "17E|descriptive")
            rows[[length(rows) + 1L]] <- p12_row("17E", MT$arm[k], MT$level[k], "vif", MT$measure[k],
                                                 MT$n_units[k], MT$n_donors[k], vif,
                                                 family = "17E|descriptive")
            p12_say(sprintf("  %-6s %-18s %-8s |r| %s | VIF %s | donors %d",
                            MT$arm[k], MT$level[k], MT$measure[k],
                            if (ok) sprintf("%.3f [%.3f, %.3f]", ob, ci[["lo"]], ci[["hi"]]) else
                                "untestable", if (ok) sprintf("%.1f", vif) else "NA", MT$n_donors[k]))
        }
        E <- do.call(rbind, rows)
        p12_say("")

        p12_say("SECTION G -- write")
        if (P12$n[["FAIL"]] > 0L) p12_abort("checks failed; stats not written")
        p12_write_stats(p12_finalise(C), F_C)
        p12_write_stats(p12_finalise(Dd), F_D)
        p12_write_stats(p12_finalise(E), F_E)
        CV$run_id <- P12$run_id; DV$run_id <- P12$run_id
        p12_output(F_CV); utils::write.csv(CV, F_CV, row.names = FALSE)
        p12_output(F_DV); utils::write.csv(DV, F_DV, row.names = FALSE)
        p12_chk(TRUE, "stats and plot files written",
                sprintf("17C %d | 17D %d | 17E %d rows", nrow(C), nrow(Dd), nrow(E)))
    }

    p12_run(main)
}
# =============================================================================
# Run in order; stop at the first section that does not pass
# =============================================================================
local({
    # Clear every Figure 17 section folder and the figure folder first.
    for (d in list.dirs(file.path(P12_RERUN, FIG_DIR), recursive = FALSE)) p12_clean(d)
    for (code in c("19", "20")) {
        v <- get(paste0("p12_section_", code))()
        if (!identical(v, "PASS"))
            stop(sprintf("section %s ended with verdict %s; later sections not run. See its log.",
                         code, v), call. = FALSE)
    }
    message("Figure 17 analysis: sections 19 and 20 all PASS. Next: run Part12_Fig17_Plots.r")
})
