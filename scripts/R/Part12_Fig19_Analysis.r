# =============================================================================
# Part12_Fig19_Analysis.r
# SEA-AD CEP-IP/CMRF, Part 12, Figure 19 analysis (severely affected donors)
#   Sections 26 to 30: the SA donors and random-removal draws (26), 16A, 16B
#   and 16F without the eleven (27, 19A), 17A to 17C (28) and 15A, 17D and 18A
#   (29) without the eleven (19B, 19C), and leave-one-donor-out influence (30).
#   Uses its own config (Part12_Fig19_config.json) and pre-set readings. Stops
#   at the first section that does not pass. Figures: Part12_Fig19_Plots.r.
#
# Usage (from the SEA-AD folder, or set SEAAD_ROOT), after Figures 15 to 18
#   source("Part12_rerun/Part12_Fig19_Analysis.r")
#
# Outputs (Part12_rerun/11_Fig19/<section>/): stats, readouts, null draws, logs.
# =============================================================================

source(file.path(Sys.getenv("SEAAD_ROOT", getwd()), "Part12_rerun", "Part12_03_Config.r"))

SCRIPT <- "Part12_Fig19_Analysis.r"
FIG_DIR <- "11_Fig19"   # every Figure 19 section folder sits inside this one
if (!isTRUE(as.integer(P12_CFG$config_version) >= 16L))
    stop(sprintf("config is version %s; Figure 19 needs version 16 or later",
                 P12_CFG$config_version), call. = FALSE)
B5_CFG_PATH <- file.path(P12_RERUN, "Part12_Fig19_config.json")
if (!file.exists(B5_CFG_PATH)) stop("Battery 5 config not found: ", B5_CFG_PATH, call. = FALSE)
B5 <- jsonlite::fromJSON(B5_CFG_PATH, simplifyVector = TRUE, simplifyDataFrame = FALSE,
                         simplifyMatrix = FALSE)
if (!identical(B5$figure_dir, FIG_DIR))
    stop("Part12_Fig19_config.json names another figure folder", call. = FALSE)

# Must be run from Part12_rerun.
local({
    src <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
    want <- file.path(P12_RERUN, SCRIPT)
    if (!is.null(src) && !identical(tolower(normalizePath(src, winslash = "/", mustWork = FALSE)),
                                    tolower(normalizePath(want, winslash = "/", mustWork = FALSE))))
        stop(sprintf("run this file from Part12_rerun: save it there, then source(\"%s\")", want),
             call. = FALSE)
})

# == Battery 5 bookkeeping ====================================================
b5_seed <- function(step) {
    code <- B5$seeds$step_codes[[step]]
    if (is.null(code)) p12_abort(sprintf("no Battery 5 step code for '%s'", step))
    as.integer(B5$seeds$master_seed) + as.integer(code)
}

# Logs the Battery 5 config and readings file; called after p12_start().
b5_header <- function() {
    rf <- p12_path(B5$readings_file)
    p12_say(sprintf("Battery 5 config: %s | version %s", p12_rel(B5_CFG_PATH),
                    B5$battery5_version))
    p12_say(sprintf("Battery 5 readings: %s%s", p12_rel(rf),
                    if (file.exists(rf)) "" else " | NOT FOUND"))
    p12_chk(file.exists(rf), "Battery 5 readings file present", p12_rel(rf))
    p12_say("")
}

# Names a Battery 5 pre-set reading (same rule as p12_reading(), other file).
b5_reading <- function(panel) {
    f <- p12_path(B5$readings_file)
    if (!file.exists(f)) p12_abort(paste0("Battery 5 readings file not found: ", f))
    x <- readLines(f, warn = FALSE, encoding = "UTF-8")
    head <- grep(sprintf("^### %s[. ]", panel), x)
    if (length(head) != 1L)
        p12_abort(sprintf("Battery 5 readings file has %d sections for %s", length(head), panel))
    p12_say(sprintf("Pre-set reading: %s | %s", sub("^### ", "", x[head]), p12_rel(f)))
}

# == Helpers copied from the Figure 16 and 17 analysis files ==================
p12_rs <- function(x, y) {
    if (length(x) < 3L) return(NA_real_)
    rx <- rank(x); ry <- rank(y)
    sx <- stats::sd(rx); sy <- stats::sd(ry)
    if (!is.finite(sx) || !is.finite(sy) || sx == 0 || sy == 0) return(NA_real_)
    stats::cor(rx, ry)
}

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

b5_write_table <- function(df, path) {
    df$run_id <- rep(P12$run_id, nrow(df))   # rep() so an empty table still writes
    if (!p12_ascii_ok(df[, vapply(df, is.character, TRUE), drop = FALSE]))
        p12_abort(sprintf("%s: a text field is not ASCII", basename(path)))
    p12_output(path)
    utils::write.csv(df, path, row.names = FALSE)
}

# Retested rows carry panel "19_<src>" and family "19|<family>".
b5_panel <- function(src) paste0("19_", src)
b5_family <- function(f) paste0("19|", f)

# == Inputs shared by the sections ============================================
b5_units <- function() {
    U <- utils::read.csv(p12_upstream_file("05_Units", "05_units.csv"), check.names = FALSE,
                         colClasses = c(donor_key = "character", region = "character"))
    p12_chk(all(U$labels_version == P12_CFG$labels$labels_version),
            "unit table carries this labels version", sprintf("%d units", nrow(U)))
    U
}
b5_donors <- function() {
    DN <- utils::read.csv(p12_upstream_file("05_Units", "05_donors.csv"), check.names = FALSE,
                          colClasses = c(donor_key = "character", donor_id = "character"))
    p12_chk(all(DN$labels_version == P12_CFG$labels$labels_version),
            "donor table carries this labels version", sprintf("%d donor rows", nrow(DN)))
    DN
}
b5_source_stats <- function(src) {
    e <- B5$exposed[[src]]
    if (is.null(e)) p12_abort(sprintf("no exposed entry for %s", src))
    S <- utils::read.csv(p12_upstream_file(e$folder, e$file), check.names = FALSE,
                         stringsAsFactors = FALSE)
    p12_chk(identical(names(S), P12_CFG$stats_columns) &&
                all(S$config_version == P12_CFG$config_version) &&
                all(S$labels_version == P12_CFG$labels$labels_version),
            sprintf("%s carries the configured columns, config and labels version", e$file),
            sprintf("%d rows", nrow(S)))
    S
}

# The section 26 outputs: the SA donors and the random-removal draws.
b5_sa_inputs <- function() {
    SA <- utils::read.csv(p12_upstream_file(file.path(FIG_DIR, "19_SADonors"), "SA_donors_19.csv"),
                          check.names = FALSE, colClasses = c(donor_key = "character",
                                                              donor_id = "character"))
    FR <- utils::read.csv(p12_upstream_file(file.path(FIG_DIR, "19_SADonors"), "Frames_19.csv"),
                          check.names = FALSE, stringsAsFactors = FALSE)
    DR <- utils::read.csv(p12_upstream_file(file.path(FIG_DIR, "19_SADonors"), "NullDraws_19.csv"),
                          check.names = FALSE, stringsAsFactors = FALSE,
                          colClasses = c(removed = "character"))
    sa <- SA$donor_key[SA$severely_affected]
    p12_chk(length(sa) == as.integer(B5$sa$expected_n), "SA donors read from section 26",
            sprintf("%d donors", length(sa)))
    # Draws as a list per frame: a character matrix, one row per draw.
    draws <- list()
    for (g in FR$frame) {
        d <- DR[DR$frame == g, ]
        d <- d[order(d$draw), ]
        k <- FR$k[FR$frame == g]
        draws[[g]] <- if (!nrow(d)) matrix(character(0), 0L, k) else
            do.call(rbind, lapply(strsplit(d$removed, ";", fixed = TRUE), function(z)
                if (length(z) == k) z else rep(NA_character_, k)))
        if (nrow(d) && anyNA(draws[[g]]))
            p12_abort(sprintf("frame %s: a draw does not hold k = %d donors", g, k))
    }
    list(SA = SA, sa = sa, frames = FR, draws = draws,
         high = SA$donor_key[SA$high_adnc])
}

# == The reading ==============================================================
# The change on removing a set: s x (e_removed - e), s = sign(e).
b5_change <- function(e_removed, e) {
    s <- sign(e); s[!is.finite(s) | s == 0] <- NA_real_
    s * (e_removed - e)
}

# Row readout: `R` one row per retested statistic, `N` rows x draws of
# random-removal estimates (NA where untestable).
b5_readout <- function(R, N) {
    tail_q <- as.numeric(B5$null$tail)
    alpha <- as.numeric(P12_CFG$inference$alpha)
    R$delta_sa <- b5_change(R$est_wo, R$est_all)
    DN <- if (nrow(R)) b5_change(N, matrix(R$est_all, nrow(N), ncol(N))) else N
    q <- function(v, pr) { v <- v[is.finite(v)]; if (length(v) < 20L) NA_real_ else
        stats::quantile(v, pr, names = FALSE, type = 7) }
    R$null_n <- if (nrow(R)) rowSums(is.finite(DN)) else integer(0)
    R$null_q_tail <- if (nrow(R)) apply(DN, 1L, q, tail_q) else numeric(0)
    R$null_q50 <- if (nrow(R)) apply(DN, 1L, q, 0.5) else numeric(0)
    R$null_q95 <- if (nrow(R)) apply(DN, 1L, q, 0.95) else numeric(0)
    R$null_est_lo <- R$est_all + sign(R$est_all) * R$null_q_tail
    R$null_est_hi <- R$est_all + sign(R$est_all) * R$null_q95
    R$p_null <- if (nrow(R)) vapply(seq_len(nrow(R)), function(i) {
        v <- DN[i, ]; v <- v[is.finite(v)]
        if (length(v) < 20L || !is.finite(R$delta_sa[i])) NA_real_ else
            (1 + sum(v <= R$delta_sa[i] + 1e-12)) / (1 + length(v))
    }, 0) else numeric(0)
    sig_a <- R$tested_all & is.finite(R$fdr_all) & R$fdr_all < alpha
    sig_w <- R$tested_wo & is.finite(R$fdr_wo) & R$fdr_wo < alpha
    same <- sign(R$est_all) == sign(R$est_wo)
    beyond <- is.finite(R$delta_sa) & is.finite(R$null_q_tail) & R$delta_sa < R$null_q_tail
    R$class <- ifelse(!R$tested_all, "Untestable",
               ifelse(sig_a & !R$tested_wo, "Untestable without the eleven",
               ifelse(sig_a & sig_w & same, "Holds",
               ifelse(sig_a & sig_w, "Reversed",
               ifelse(sig_a & !R$null_ok, "Lost, null untestable",
               ifelse(sig_a & beyond, "Lost, beyond random removal",
               ifelse(sig_a, "Lost, within random removal",
               ifelse(sig_w, "Gained", "Neither"))))))))
    R
}

# Class counts of a readout, for the log.
b5_log_classes <- function(R, by) {
    K <- B5$classes
    for (g in unique(R[[by]])) {
        r <- R[R[[by]] == g, ]
        tb <- table(factor(r$class, K))
        p12_say(sprintf("    %-28s %s", g, paste(sprintf("%s %d", names(tb)[tb > 0], tb[tb > 0]),
                                                    collapse = " | ")))
    }
}

# All-donor rows against the Battery 1 to 4 stats file, on the key columns.
b5_reproduce <- function(label, S, X, with_p = FALSE, tol = 1e-8) {
    k1 <- paste(S$arm, S$level, S$quantity, S$measure)
    k2 <- paste(X$arm, X$level, X$quantity, X$measure)
    j <- match(k2, k1)
    miss <- sum(is.na(j))
    d_e <- abs(X$estimate - S$estimate[j]); d_e[is.na(X$estimate) & is.na(S$estimate[j])] <- 0
    bad_e <- sum(!(is.finite(d_e) & d_e <= tol * pmax(1, abs(X$estimate))), na.rm = TRUE)
    nu <- if ("n_units" %in% names(X)) X$n_units else rep(NA_integer_, nrow(X))
    bad_n <- sum(!((is.na(nu) & is.na(S$n_units[j])) | (!is.na(nu) & !is.na(S$n_units[j]) &
                                                           nu == S$n_units[j])), na.rm = TRUE) +
        sum(X$n_donors != S$n_donors[j], na.rm = TRUE)
    bad_p <- 0L
    if (with_p) {
        d_p <- abs(X$p - S$p[j]); d_p[is.na(X$p) & is.na(S$p[j])] <- 0
        bad_p <- sum(!(is.finite(d_p) & d_p <= tol * pmax(1e-300, abs(X$p)) + 1e-15), na.rm = TRUE)
    }
    p12_chk(miss == 0L && bad_e == 0L && bad_n == 0L && bad_p == 0L,
            sprintf("%s: all-donor rows reproduce the Battery 1 to 4 stats file", label),
            sprintf("%d rows | %d unmatched | %d estimates, %d counts%s differ | largest estimate difference %.1e",
                    nrow(X), miss, bad_e, bad_n, if (with_p) sprintf(", %d p", bad_p) else "",
                    suppressWarnings(max(c(0, d_e), na.rm = TRUE))))
}

# =============================================================================
# SECTION 26: the SA donors, the frames and the random-removal draws
# Checks the Battery 1 to 4 logs, reads the SA flag
# (eleven donors, all High ADNC) and builds the frames and the draws.
# =============================================================================
p12_section_26 <- function() {
    STEP <- "26"
    D <- p12_step_dirs(STEP, "19_SADonors", parent = FIG_DIR)
    p12_clean(D$out)
    p12_start(STEP, SCRIPT, D)
    b5_header()
    F_SA <- file.path(D$out, "SA_donors_19.csv")
    F_FR <- file.path(D$out, "Frames_19.csv")
    F_DR <- file.path(D$out, "NullDraws_19.csv")

    main <- function() {
        R_NULL <- as.integer(B5$null$draws)
        b5_reading("Battery 5 definitions")

        p12_say("SECTION A -- the Battery 1 to 4 logs")
        sl <- B5$source_logs[names(B5$source_logs) != "_note"]
        for (k in names(sl)) {
            lg <- file.path(P12_RERUN, sl[[k]], sprintf("Part12_%s_log.txt", basename(sl[[k]])))
            if (!file.exists(lg)) { p12_chk(FALSE, sprintf("%s log present", k), p12_rel(lg)); next }
            x <- readLines(lg, warn = FALSE)
            scr <- sub("^Script: (\\S+).*$", "\\1", grep("^Script: ", x, value = TRUE)[1])
            ok_scr <- !is.na(scr) && file.exists(file.path(P12_RERUN, scr))
            p12_chk(any(x == "VERDICT: PASS") && ok_scr,
                    sprintf("%s log: PASS, and its analysis file is in Part12_rerun", k),
                    sprintf("%s %s", scr, if (ok_scr) "present" else "NOT FOUND"))
        }
        if (P12$n[["FAIL"]] > 0L) p12_abort("a Battery 1 to 4 log is missing or did not pass; nothing computed")
        p12_say("")

        p12_say("SECTION B -- the SA flag")
        TC <- P12_CFG$table_columns
        tv <- tolower(B5$sa$true_values); fv <- tolower(B5$sa$false_values)
        OB <- list()
        for (a in P12_ARMS) for (ds in P12_DATASETS) {
            f <- p12_log_input(p12_input("obs_dataset", arm = a, dataset = ds),
                               paste(a, ds, "obs_dataset"))
            o <- p12_read_csv_cols(f, c(donor = TC$obs_dataset$donor, sa = B5$sa$obs_column),
                                   paste(a, ds, "obs_dataset"))
            v <- tolower(trimws(ifelse(is.na(o$sa), "", o$sa)))
            OB[[length(OB) + 1L]] <- data.frame(donor_key = p12_donor_key(o$donor), v = v,
                                                stringsAsFactors = FALSE)
        }
        OB <- unique(do.call(rbind, OB))
        unk <- setdiff(unique(OB$v), c(tv, fv))
        p12_chk(!length(unk), "every SA value is a listed true or false value",
                if (length(unk)) paste("unlisted:", paste(unk, collapse = ", ")) else
                    sprintf("values: %s", paste(sort(unique(OB$v)), collapse = ", ")))
        OB$flag <- OB$v %in% tv
        nf <- tapply(OB$flag, OB$donor_key, function(z) length(unique(z)))
        p12_chk(all(nf == 1L), "the SA flag is constant within donor",
                sprintf("%d donors | %d with both values", length(nf), sum(nf > 1L)))
        flag <- tapply(OB$flag, OB$donor_key, any)

        DN <- b5_donors()
        HA <- B5$high_adnc
        hc <- p12_col(DN, HA$column, "donor table")
        dk <- sort(unique(DN$donor_key))
        adnc <- tapply(DN[[hc]], DN$donor_key, function(z) { z <- unique(z[is.finite(z)])
            if (length(z) == 1L) z else NA_real_ })[dk]
        p12_chk(all(dk %in% names(flag)), "every donor of the donor table has an SA value",
                sprintf("%d of %d donors", sum(dk %in% names(flag)), length(dk)))
        p12_chk(all(is.na(adnc) | adnc %in% as.numeric(HA$allowed_values)),
                sprintf("%s takes only the listed values", HA$column),
                sprintf("values: %s | missing %d", paste(sort(unique(adnc[is.finite(adnc)])),
                                                          collapse = ", "), sum(!is.finite(adnc))))
        sa <- names(flag)[flag]
        p12_chk(length(sa) == as.integer(B5$sa$expected_n),
                sprintf("exactly %s donors are flagged", B5$sa$expected_n),
                sprintf("%d flagged", length(sa)))
        p12_chk(all(sa %in% dk), "every SA donor is in the donor table",
                sprintf("%d of %d", sum(sa %in% dk), length(sa)))
        p12_chk(all(is.finite(adnc[sa]) & adnc[sa] == as.numeric(HA$value)),
                sprintf("every SA donor has %s %s (High)", HA$column, HA$value),
                sprintf("%s", paste(sprintf("%s:%s", sa, adnc[sa]), collapse = " ")))
        ids <- tapply(DN$donor_id, DN$donor_key, `[`, 1)
        p12_note("SA donors", paste(sort(ids[sa]), collapse = ", "))
        arms_of <- tapply(DN$arm, DN$donor_key, function(z) paste(P12_ARMS[P12_ARMS %in% z], collapse = ";"))
        SA <- data.frame(donor_key = dk, donor_id = unname(ids[dk]),
                         severely_affected = dk %in% sa,
                         adnc_rank = unname(adnc), high_adnc = is.finite(adnc) & adnc == HA$value,
                         arms = unname(arms_of[dk]), stringsAsFactors = FALSE)
        p12_note("donors", sprintf("%d | SA %d | High ADNC %d | High ADNC not SA %d", nrow(SA),
                                   sum(SA$severely_affected), sum(SA$high_adnc),
                                   sum(SA$high_adnc & !SA$severely_affected)))
        if (P12$n[["FAIL"]] > 0L) p12_abort("SA checks failed; nothing drawn")
        p12_say("")

        p12_say("SECTION C -- frames")
        U <- b5_units()
        J <- utils::read.csv(p12_upstream_file("09_Fig17/17AB_Histology", "Joined_17AB_units.csv"),
                             check.names = FALSE, colClasses = c(donor_key = "character",
                                                                 region = "character"))
        p12_chk(nrow(J) == nrow(U) && identical(J$donor_key, U$donor_key) && identical(J$arm, U$arm),
                "the joined histology table holds the unit table's units in its order",
                sprintf("%d units", nrow(J)))
        FR <- list()
        for (a in P12_ARMS) {
            FR[[paste("U", a, sep = "|")]] <- sort(unique(U$donor_key[U$arm == a]))
            FR[[paste("J19", a, sep = "|")]] <-
                sort(unique(J$donor_key[J$arm == a & (J$joined_histology | J$joined_cps)]))
            FR[[paste("JHU", a, sep = "|")]] <-
                sort(unique(J$donor_key[J$arm == a & J$joined_histology]))
        }
        for (a in P12_ARMS) for (ds in P12_DATASETS)
            FR[[paste("DN", a, ds, sep = "|")]] <- sort(unique(DN$donor_key[DN$arm == a & DN$dataset == ds]))
        high_non <- SA$donor_key[SA$high_adnc & !SA$severely_affected]
        fr_rows <- list(); dr_rows <- list()
        for (i in seq_along(FR)) {
            g <- names(FR)[i]; don <- FR[[g]]
            k <- sum(don %in% sa); pool <- don[don %in% high_non]
            status <- if (k == 0L) "no SA donors in the frame" else
                if (length(pool) < k) "untestable: pool smaller than k" else "testable"
            fr_rows[[i]] <- data.frame(frame = g, n_donors = length(don), k = k,
                                       n_pool = length(pool), null_status = status,
                                       stringsAsFactors = FALSE)
            if (status == "testable") {
                set.seed(b5_seed(STEP) + i)
                rm <- vapply(seq_len(R_NULL), function(r)
                    paste(sort(pool[sample.int(length(pool), k)]), collapse = ";"), "")
                dr_rows[[length(dr_rows) + 1L]] <- data.frame(frame = g, draw = seq_len(R_NULL),
                                                              removed = rm, stringsAsFactors = FALSE)
            }
            p12_say(sprintf("  %-16s donors %3d | SA %2d | pool %3d | %s%s", g, length(don), k,
                            length(pool), status,
                            if (status == "testable") sprintf(" | %d distinct draws of %d",
                                length(unique(dr_rows[[length(dr_rows)]]$removed)), R_NULL) else ""))
        }
        FRt <- do.call(rbind, fr_rows)
        DRt <- if (length(dr_rows)) do.call(rbind, dr_rows) else
            data.frame(frame = character(0), draw = integer(0), removed = character(0))
        p12_chk(all(FRt$null_status != "untestable: pool smaller than k"),
                "every frame with SA donors has a pool of at least k donors",
                paste(FRt$frame[FRt$null_status == "untestable: pool smaller than k"], collapse = ", "),
                warn_only = TRUE)
        p12_chk(!any(unlist(strsplit(DRt$removed, ";", fixed = TRUE)) %in% sa),
                "no draw removes an SA donor", sprintf("%d draws", nrow(DRt)))
        p12_say("")

        p12_say("SECTION D -- write")
        if (P12$n[["FAIL"]] > 0L) p12_abort("checks failed; nothing written")
        b5_write_table(SA, F_SA)
        b5_write_table(FRt, F_FR)
        b5_write_table(DRt, F_DR)
        p12_chk(TRUE, "tables written", sprintf("%d donors | %d frames | %d draws", nrow(SA),
                                                nrow(FRt), nrow(DRt)))
    }

    p12_run(main)
}

# =============================================================================
# SECTION 27: 19A, 16A, 16B and 16F between donors without the eleven
# Code of sections 13 and 17; random removal gives estimates only.
# =============================================================================

# Donor share of variance of y with donor codes dd (1..nd), as section 13.
b5_share <- function(y, dd) {
    nd <- max(c(0L, dd)); n_d <- tabulate(dd, nd)
    s <- rowsum(y, dd, reorder = TRUE)[, 1]
    ssb <- sum(s^2 / n_d[n_d > 0]) - length(y) * mean(y)^2
    ssb / sum((y - mean(y))^2)
}

# Median split-half rs over the donor columns `keep`, as section 13.
b5_splithalf <- function(A1, A2, keep) {
    if (length(keep) < 3L) return(NA_real_)
    stats::median(vapply(seq_len(nrow(A1)), function(s) p12_rs(A1[s, keep], A2[s, keep]), 0),
                  na.rm = TRUE)
}

# 16B bootstrap on one chunk of resamples. Copied from Part12_Fig16_Analysis.r.
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

# 16F bootstrap on one chunk of resamples. Copied from Part12_Fig16_Analysis.r.
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

# Long table of random-removal estimates for writing.
b5_null_long <- function(R, N) {
    if (!nrow(R) || !ncol(N)) return(data.frame(row_id = character(0), draw = integer(0),
                                                estimate = numeric(0)))
    data.frame(row_id = rep(R$row_id, times = ncol(N)), draw = rep(seq_len(ncol(N)), each = nrow(N)),
               estimate = as.vector(N), stringsAsFactors = FALSE)
}

# The 16B halves of the all-donor run, from section 13's split stream.
b5_16b_splits <- function(U, QB, MIND, N_SPLITS) {
    set.seed(p12_seed("13") + 1L)
    split_state <- get(".Random.seed", envir = globalenv())
    SPL <- list()
    for (a in P12_ARMS) for (q in QB) {
        u <- U[U$arm == a & is.finite(U[[q]]), ]
        kk <- table(u$donor_key)
        dons <- names(kk)[kk >= 2L]
        u <- u[u$donor_key %in% dons, ]
        A1 <- A2 <- NULL
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
        }
        SPL[[paste(a, q)]] <- list(dons = dons, A1 = A1, A2 = A2, n_units = nrow(u),
                                   units_of = table(u$donor_key))
    }
    SPL
}

# 16F between-donor tasks on the donors `dset`, as section 17.
b5_16f_tasks <- function(U, dset, QS, FIELD, MIND) {
    tasks <- list(); meta <- list()
    for (a in P12_ARMS) for (q in names(QS)) {
        u <- U[U$arm == a & U$donor_key %in% dset, ]
        x <- u[[QS[[q]]]]; y <- u[[paste0("lib_", FIELD)]]
        ok <- is.finite(x) & is.finite(y)
        if (!any(ok)) { x <- y <- numeric(0); d <- character(0) } else {
            mx <- tapply(x[ok], u$donor_key[ok], mean); my <- tapply(y[ok], u$donor_key[ok], mean)
            d <- names(mx); x <- as.numeric(mx); y <- as.numeric(my[d])
        }
        tasks[[length(tasks) + 1L]] <- list(x = x, y = y, d = match(d, dset), min_d = MIND, key = d)
        meta[[length(meta) + 1L]] <- data.frame(arm = a, q = q, n_units = length(x),
                                                n_donors = length(unique(d)))
    }
    list(tasks = tasks, meta = do.call(rbind, meta))
}

# rs of each 16F task less the donors in `drop`.
b5_16f_obs <- function(TK, MIND, drop = character(0)) {
    vapply(seq_along(TK$tasks), function(i) {
        t <- TK$tasks[[i]]
        dr <- if (is.list(drop)) drop[[TK$meta$arm[i]]] else drop
        r <- which(!(t$key %in% dr))
        if (length(unique(t$d[r])) >= MIND) p12_rs(t$x[r], t$y[r]) else NA_real_
    }, 0)
}

# 16F rows from task rs and, when given, the bootstrap matrix.
b5_16f_rows <- function(obs, BS, meta, panel, famfun, QS, LV, FIELD) {
    rows <- list()
    key <- function(a, q) which(meta$arm == a & meta$q == q)
    for (a in P12_ARMS) {
        for (q in names(QS)) {
            i <- key(a, q); ok <- is.finite(obs[i])
            ci <- if (ok && !is.null(BS)) p12_boot_ci(BS[, i]) else c(lo = NA_real_, hi = NA_real_)
            rows[[length(rows) + 1L]] <- p12_row(panel, a, LV, q, FIELD, meta$n_units[i],
                                                 meta$n_donors[i], obs[i], ci[["lo"]], ci[["hi"]],
                                                 if (ok && !is.null(BS)) p12_boot_p(BS[, i]) else NA,
                                                 famfun(paste("16F", a, LV, "rs", sep = "|")))
        }
        im <- key(a, "module")
        con <- list(module_minus_mito = abs(obs[im]) - abs(obs[key(a, "mito")]),
                    module_minus_depth = obs[im] - obs[key(a, "depth")])
        cbs <- list(module_minus_mito = if (!is.null(BS)) abs(BS[, im]) - abs(BS[, key(a, "mito")]),
                    module_minus_depth = if (!is.null(BS)) BS[, im] - BS[, key(a, "depth")])
        if (a != "HKG8") {
            ih <- key("HKG8", "module")
            con$arm_minus_hkg8 <- obs[im] - obs[ih]
            cbs$arm_minus_hkg8 <- if (!is.null(BS)) BS[, im] - BS[, ih]
        }
        for (cn in names(con)) {
            ok <- is.finite(con[[cn]])
            ci <- if (ok && !is.null(BS)) p12_boot_ci(cbs[[cn]]) else c(lo = NA_real_, hi = NA_real_)
            rows[[length(rows) + 1L]] <- p12_row(panel, a, LV, cn, FIELD, meta$n_units[im],
                                                 meta$n_donors[im], con[[cn]], ci[["lo"]], ci[["hi"]],
                                                 if (ok && !is.null(BS)) p12_boot_p(cbs[[cn]]) else NA,
                                                 famfun(paste("16F", a, LV, "contrast", sep = "|")))
        }
    }
    do.call(rbind, rows)
}

# 16F estimates as a plain vector, for the random-removal and influence loops.
b5_16f_est <- function(obs, meta, QS) {
    key <- function(a, q) which(meta$arm == a & meta$q == q)
    out <- numeric(0)
    for (a in P12_ARMS) {
        out <- c(out, obs[vapply(names(QS), function(q) key(a, q), 0L)])
        im <- key(a, "module")
        out <- c(out, abs(obs[im]) - abs(obs[key(a, "mito")]), obs[im] - obs[key(a, "depth")])
        if (a != "HKG8") out <- c(out, obs[im] - obs[key("HKG8", "module")])
    }
    unname(out)
}

p12_section_27 <- function() {
    STEP <- "27"
    D <- p12_step_dirs(STEP, "19_UnitPanels", parent = FIG_DIR)
    p12_clean(D$out)
    p12_start(STEP, SCRIPT, D, pkgs = "parallel")
    b5_header()
    F_A <- file.path(D$out, "Stats_19_16A.csv")
    F_B <- file.path(D$out, "Stats_19_16B.csv")
    F_F <- file.path(D$out, "Stats_19_16F.csv")
    F_R <- file.path(D$out, "Readout_19_Units.csv")
    F_N <- file.path(D$out, "Null_19_Units.csv")
    QA <- c("module_post_pct", "module_resid_tq4", "mito_log1p", "contamination", "ctrl_pct_null")
    QB <- c("module_post_pct", "module_resid_tq4", "mito_log1p")
    QS <- c(module = "module_post_pct", mito = "mito_log1p", ncells = "log_n_cells",
            depth = "log_median_umi")
    N_SPLITS <- 200L
    LV <- "between donors"; FIELD <- "intronic"

    main <- function() {
        R_PERM <- as.integer(P12_CFG$inference$permutation_R)
        B <- as.integer(P12_CFG$inference$bootstrap_B)
        MIND <- as.integer(P12_CFG$floors$donors_per_test)
        b5_reading("Battery 5 definitions")
        b5_reading("19A")

        p12_say("SECTION A -- inputs")
        U <- b5_units()
        miss <- setdiff(c(QA, QS, paste0("lib_", FIELD)), names(U))
        if (length(miss)) p12_abort(paste("unit table lacks:", paste(miss, collapse = ", ")))
        SA_ <- b5_sa_inputs(); sa <- SA_$sa
        S16A <- b5_source_stats("16A"); S16B <- b5_source_stats("16B"); S16F <- b5_source_stats("16F")
        R_NULL <- as.integer(B5$null$draws)
        drw <- function(a) SA_$draws[[paste("U", a, sep = "|")]]
        null_ok <- function(a) SA_$frames$null_status[SA_$frames$frame == paste("U", a, sep = "|")] != "untestable: pool smaller than k"
        cl <- p12_cluster(B)
        if (!is.null(cl)) on.exit(parallel::stopCluster(cl), add = TRUE)
        if (P12$n[["FAIL"]] > 0L) p12_abort("input checks failed; nothing computed")
        p12_say("")

        RD <- list(); NL <- list()
        push_readout <- function(src, a, lv, q, m, S, X, tested_all, tested_wo, nv, nk,
                                 nok = null_ok(a)) {
            k1 <- paste(S$arm, S$level, S$quantity, S$measure)
            j <- match(paste(a, lv, q, m), k1)
            RD[[length(RD) + 1L]] <<- data.frame(
                row_id = paste(src, a, lv, q, m, sep = "|"), panel = src, arm = a, level = lv,
                quantity = q, measure = m, frame = paste("U", a, sep = "|"), k_removed = nk,
                est_all = S$estimate[j], fdr_all = S$fdr[j], n_donors_all = S$n_donors[j],
                est_wo = X$estimate, lo_wo = X$lo, hi_wo = X$hi, fdr_wo = X$fdr,
                n_donors_wo = X$n_donors, tested_all = tested_all, tested_wo = tested_wo,
                null_ok = nok, stringsAsFactors = FALSE)
            NL[[length(NL) + 1L]] <<- nv
        }
        k_of <- function(a) SA_$frames$k[SA_$frames$frame == paste("U", a, sep = "|")]

        # == 16A ==============================================================
        p12_say(sprintf("SECTION B -- 16A, donor share against %d within-region shuffles", R_PERM))
        rep_a <- list(); rows_a <- list(); nulls_a <- list()
        set.seed(b5_seed(STEP))
        for (a in P12_ARMS) for (q in QA) {
            u <- U[U$arm == a & is.finite(U[[q]]), ]
            don <- as.integer(factor(u$donor_key)); nd <- max(c(0L, don))
            est0 <- if (nd >= MIND && isTRUE(stats::var(u[[q]]) > 0)) b5_share(u[[q]], don) else NA_real_
            rep_a[[length(rep_a) + 1L]] <- p12_row("16A", a, "units", q, "donor share of variance",
                                                   nrow(u), nd, est0, family = "x")
            w <- u[!(u$donor_key %in% sa), ]
            y <- w[[q]]; dw <- as.integer(factor(w$donor_key)); ndw <- max(c(0L, dw))
            est <- p <- lo <- hi <- NA_real_
            if (ndw >= MIND && isTRUE(stats::var(y) > 0)) {
                est <- b5_share(y, dw)
                reg <- split(seq_along(y), w$region)
                nul <- numeric(R_PERM)
                for (r in seq_len(R_PERM)) {
                    dd <- dw
                    for (g in reg) if (length(g) > 1L) dd[g] <- dw[g][sample.int(length(g))]
                    nul[r] <- b5_share(y, dd)
                }
                p <- (sum(nul >= est - 1e-12) + 1) / (R_PERM + 1)
                lo <- stats::median(nul); hi <- stats::quantile(nul, 0.95, names = FALSE)
            }
            rows_a[[length(rows_a) + 1L]] <- p12_row(b5_panel("16A"), a, "units", q,
                                                     "donor share of variance", nrow(w), ndw, est,
                                                     lo, hi, p, b5_family(paste("16A", a, sep = "|")))
            # random removal: estimate only
            dm <- drw(a)
            nv <- rep(NA_real_, R_NULL)
            if (nrow(dm)) for (r in seq_len(nrow(dm))) {
                z <- u[!(u$donor_key %in% dm[r, ]), ]
                dz <- as.integer(factor(z$donor_key))
                if (max(c(0L, dz)) >= MIND && isTRUE(stats::var(z[[q]]) > 0)) nv[r] <- b5_share(z[[q]], dz)
            }
            nulls_a[[length(nulls_a) + 1L]] <- nv
        }
        b5_reproduce("16A", S16A, do.call(rbind, rep_a))
        A <- p12_bh_by_family(do.call(rbind, rows_a))
        i <- 0L
        for (a in P12_ARMS) for (q in QA) {
            i <- i + 1L
            j <- which(S16A$arm == a & S16A$quantity == q)
            push_readout("16A", a, "units", q, "donor share of variance", S16A, A[i, ],
                         is.finite(S16A$p[j]), is.finite(A$p[i]), nulls_a[[i]], k_of(a))
            p12_say(sprintf("  %-6s %-17s donors %3d -> %3d | share %5.1f%% -> %5.1f%% | FDR %s -> %s",
                            a, q, S16A$n_donors[j], A$n_donors[i], 100 * S16A$estimate[j],
                            100 * A$estimate[i], p12_fmt(S16A$fdr[j], 1 / (R_PERM + 1)),
                            p12_fmt(A$fdr[i], 1 / (R_PERM + 1))))
        }
        p12_say("")

        # == 16B ==============================================================
        p12_say(sprintf("SECTION C -- 16B, split-half rs over %d splits (B = %d)", N_SPLITS, B))
        SPL <- b5_16b_splits(U, QB, MIND, N_SPLITS)
        rep_b <- list(); rows_b <- list(); nulls_b <- list()
        for (a in P12_ARMS) for (q in QB) {
            sp <- SPL[[paste(a, q)]]
            rep_b[[length(rep_b) + 1L]] <- p12_row("16B", a, "donors", q, "split-half rs, median over splits",
                                                   sp$n_units, length(sp$dons),
                                                   if (is.null(sp$A1)) NA_real_ else
                                                       b5_splithalf(sp$A1, sp$A2, seq_along(sp$dons)),
                                                   family = "x")
        }
        b5_reproduce("16B", S16B, do.call(rbind, rep_b))
        i <- 0L
        for (a in P12_ARMS) for (q in QB) {
            i <- i + 1L
            sp <- SPL[[paste(a, q)]]
            keep <- which(!(sp$dons %in% sa))
            est <- lo <- hi <- p <- NA_real_
            if (!is.null(sp$A1) && length(keep) >= MIND) {
                A1 <- sp$A1[, keep, drop = FALSE]; A2 <- sp$A2[, keep, drop = FALSE]
                est <- b5_splithalf(A1, A2, seq_along(keep))
                set.seed(b5_seed(STEP) + 1L)
                idx <- matrix(sample.int(length(keep), length(keep) * B, replace = TRUE), nrow = B)
                bs <- p12_boot_run(cl, idx, p12_16b_chunk, A = A1, Bm = A2)
                ci <- p12_boot_ci(bs); lo <- ci[["lo"]]; hi <- ci[["hi"]]; p <- p12_boot_p(bs)
            }
            nu <- if (length(keep)) sum(sp$units_of[sp$dons[keep]]) else 0L
            rows_b[[i]] <- p12_row(b5_panel("16B"), a, "donors", q, "split-half rs, median over splits",
                                   nu, length(keep), est, lo, hi, p, b5_family(paste("16B", a, sep = "|")))
            dm <- drw(a)
            nv <- rep(NA_real_, R_NULL)
            if (!is.null(sp$A1) && nrow(dm)) for (r in seq_len(nrow(dm))) {
                kp <- which(!(sp$dons %in% dm[r, ]))
                if (length(kp) >= MIND) nv[r] <- b5_splithalf(sp$A1, sp$A2, kp)
            }
            nulls_b[[i]] <- nv
        }
        Bt <- p12_bh_by_family(do.call(rbind, rows_b))
        i <- 0L
        for (a in P12_ARMS) for (q in QB) {
            i <- i + 1L
            j <- which(S16B$arm == a & S16B$quantity == q)
            push_readout("16B", a, "donors", q, "split-half rs, median over splits", S16B, Bt[i, ],
                         is.finite(S16B$p[j]), is.finite(Bt$p[i]), nulls_b[[i]], k_of(a))
            p12_say(sprintf("  %-6s %-17s donors %3d -> %3d | rs %+.3f -> %+.3f [%+.3f, %+.3f] | FDR %s -> %s",
                            a, q, S16B$n_donors[j], Bt$n_donors[i], S16B$estimate[j], Bt$estimate[i],
                            Bt$lo[i], Bt$hi[i], p12_fmt(S16B$fdr[j], p12_boot_floor(B)),
                            p12_fmt(Bt$fdr[i], p12_boot_floor(B))))
        }
        p12_say("")

        # == 16F between donors ===============================================
        p12_say(sprintf("SECTION D -- 16F between donors (field %s), donor bootstrap (B = %d)", FIELD, B))
        dons <- sort(unique(U$donor_key))
        T0 <- b5_16f_tasks(U, dons, QS, FIELD, MIND)
        b5_reproduce("16F", S16F[S16F$level == LV & S16F$measure == FIELD, ],
                     b5_16f_rows(b5_16f_obs(T0, MIND), NULL, T0$meta, "16F", identity, QS, LV, FIELD))
        dons_w <- dons[!(dons %in% sa)]
        TW <- b5_16f_tasks(U, dons_w, QS, FIELD, MIND)
        p12_chk(isTRUE(all.equal(b5_16f_obs(TW, MIND), b5_16f_obs(T0, MIND, sa), tolerance = 1e-12)),
                "16F: tasks rebuilt without the eleven equal the all-donor tasks less their rows",
                sprintf("%d tasks", length(TW$tasks)))
        set.seed(b5_seed(STEP) + 2L)
        idx <- matrix(sample.int(length(dons_w), length(dons_w) * B, replace = TRUE), nrow = B)
        BS <- p12_boot_run(cl, idx, p12_16fg_chunk, tasks = TW$tasks, nd = length(dons_w))
        Ft <- p12_bh_by_family(b5_16f_rows(b5_16f_obs(TW, MIND), BS, TW$meta, b5_panel("16F"), b5_family,
                                           QS, LV, FIELD))
        p12_chk(length(b5_16f_est(b5_16f_obs(TW, MIND), TW$meta, QS)) == nrow(Ft) &&
                    isTRUE(all.equal(b5_16f_est(b5_16f_obs(TW, MIND), TW$meta, QS), Ft$estimate)),
                "16F: the plain estimate vector follows the stats rows", sprintf("%d rows", nrow(Ft)))
        # Random removal: draw r of each arm's frame (both frames for a contrast).
        NF <- matrix(NA_real_, nrow(Ft), R_NULL)
        for (r in seq_len(R_NULL)) {
            drop_of <- lapply(stats::setNames(P12_ARMS, P12_ARMS), function(a) {
                dm <- drw(a); if (nrow(dm) >= r) dm[r, ] else character(0) })
            est <- b5_16f_est(b5_16f_obs(T0, MIND, drop_of), T0$meta, QS)
            # a frame without a testable null contributes NA
            for (a in P12_ARMS) if (!nrow(drw(a))) est[Ft$arm == a] <- NA_real_
            NF[, r] <- est
        }
        # HKG8 with SA donors but too small a pool leaves arm - HKG8 without a null
        hk_bad <- !null_ok("HKG8")
        for (i in seq_len(nrow(Ft))) {
            a <- Ft$arm[i]; q <- Ft$quantity[i]
            j <- which(S16F$arm == a & S16F$level == LV & S16F$measure == FIELD & S16F$quantity == q)
            nv <- NF[i, ]
            if (q == "arm_minus_hkg8" && hk_bad) nv[] <- NA_real_
            push_readout("16F", a, LV, q, FIELD, S16F[S16F$level == LV & S16F$measure == FIELD, ], Ft[i, ],
                         is.finite(S16F$p[j]), is.finite(Ft$p[i]), nv, k_of(a),
                         nok = null_ok(a) && !(q == "arm_minus_hkg8" && hk_bad))
        }
        for (a in P12_ARMS) {
            g <- Ft[Ft$arm == a, ]
            m <- g[g$quantity == "module", ]
            j <- which(S16F$arm == a & S16F$level == LV & S16F$measure == FIELD & S16F$quantity == "module")
            p12_say(sprintf("  %-6s Module rs %+.3f -> %+.3f [%+.3f, %+.3f] | FDR %s -> %s | donors %d -> %d",
                            a, S16F$estimate[j], m$estimate, m$lo, m$hi,
                            p12_fmt(S16F$fdr[j], p12_boot_floor(B)), p12_fmt(m$fdr, p12_boot_floor(B)),
                            S16F$n_donors[j], m$n_donors))
        }
        p12_say("")

        p12_say("SECTION E -- the reading")
        RDt <- do.call(rbind, RD)
        N <- do.call(rbind, NL)
        for (i in seq_len(nrow(RDt))) if (!RDt$null_ok[i] || RDt$k_removed[i] == 0L) N[i, ] <- NA_real_
        RDt <- b5_readout(RDt, N)
        b5_log_classes(transform(RDt, pa = paste(panel, arm)), "pa")
        p12_say("")

        p12_say("SECTION F -- write")
        if (P12$n[["FAIL"]] > 0L) p12_abort("checks failed; stats not written")
        p12_write_stats(p12_finalise(A), F_A)
        p12_write_stats(p12_finalise(Bt), F_B)
        p12_write_stats(p12_finalise(Ft), F_F)
        b5_write_table(RDt, F_R)
        b5_write_table(b5_null_long(RDt, N), F_N)
        p12_chk(TRUE, "stats and readout written",
                sprintf("16A %d | 16B %d | 16F %d | readout %d rows", nrow(A), nrow(Bt), nrow(Ft), nrow(RDt)))
    }

    p12_run(main)
}

# =============================================================================
# SECTION 28: 19B and 19C, 17A, 17B and 17C without the eleven
# Code of sections 19 and 20, on the joined histology table of section 19.
# =============================================================================

# == Copied from Part12_Fig17_Analysis.r ======================================
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

# Donor bootstrap per donor group, as p12_17_boot(), with the step's seed.
b5_17_boot <- function(tasks, group_of, dons_by_group, seed, B, cl, ncol_out) {
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
        set.seed(seed + 1000L * i)
        idx <- matrix(sample.int(nd, nd * B, replace = TRUE), nrow = B, ncol = nd)
        BS[, cols] <- p12_boot_run(cl, idx, p12_17_chunk, tasks = tk, nd = nd, ncol_out = loc - 1L,
                                   stat = p12_17_stat, min_d_boot = BOOT_MIN_D)
    }
    BS
}

# Random-removal estimates for one chunk of draws; returns draws x columns.
b5_17_null_chunk <- function(rs, tasks, keys, draws, frame_of, ncol_out, stat) {
    out <- matrix(NA_real_, length(rs), ncol_out)
    for (i in seq_along(rs)) {
        r <- rs[i]
        for (j in seq_along(tasks)) {
            t <- tasks[[j]]
            if (!t$testable) next
            dm <- draws[[frame_of[j]]]
            if (nrow(dm) < r) next
            keep <- which(!(keys[[j]] %in% dm[r, ]))
            v <- stat(t, keep)
            out[i, t$col:(t$col + length(v) - 1L)] <- v
        }
    }
    out
}

p12_section_28 <- function() {
    STEP <- "28"
    D <- p12_step_dirs(STEP, "19_Histology", parent = FIG_DIR)
    p12_clean(D$out)
    p12_start(STEP, SCRIPT, D, pkgs = "parallel")
    b5_header()
    F_A <- file.path(D$out, "Stats_19_17A.csv")
    F_B <- file.path(D$out, "Stats_19_17B.csv")
    F_C <- file.path(D$out, "Stats_19_17C.csv")
    F_V <- file.path(D$out, "Verdicts_19_17C.csv")
    F_R <- file.path(D$out, "Readout_19_Histology.csv")
    F_N <- file.path(D$out, "Null_19_Histology.csv")
    QS <- c("module_post_pct", "mito_log1p", "cyto_share", "contamination", "log_n_cells")
    LA <- "within donor"; LB <- "between donors"; LC <- "within donor"

    main <- function() {
        B <- as.integer(P12_CFG$inference$bootstrap_B)
        MIND <- as.integer(P12_CFG$floors$donors_per_test)
        MINU <- as.integer(P12_CFG$floors$finite_units_within_donor)
        R_NULL <- as.integer(B5$null$draws)
        HI <- P12_CFG$histology
        IBA <- vapply(HI$iba1_17A, function(z) z$measure, "")
        CPS <- vapply(HI$cps, function(z) z$column, "")
        MEAS <- unlist(lapply(HI$families, function(f) f$measures), use.names = FALSE)
        b5_reading("Battery 5 definitions")
        b5_reading("19B")
        b5_reading("19C")

        p12_say("SECTION A -- inputs")
        J <- utils::read.csv(p12_upstream_file("09_Fig17/17AB_Histology", "Joined_17AB_units.csv"),
                             check.names = FALSE, colClasses = c(donor_key = "character",
                                                                 region = "character"))
        miss <- setdiff(c(QS, MEAS, CPS, "joined_histology", "joined_cps"), names(J))
        if (length(miss)) p12_abort(paste("joined table lacks:", paste(miss, collapse = ", ")))
        SA_ <- b5_sa_inputs(); sa <- SA_$sa
        S17A <- b5_source_stats("17A"); S17B <- b5_source_stats("17B"); S17C <- b5_source_stats("17C")
        cl <- p12_cluster(B)
        if (!is.null(cl)) on.exit(parallel::stopCluster(cl), add = TRUE)
        if (P12$n[["FAIL"]] > 0L) p12_abort("input checks failed; nothing computed")
        fr <- SA_$frames
        fstat <- function(g) fr$null_status[fr$frame == g]
        fk <- function(g) fr$k[fr$frame == g]
        p12_say("")

        # Tasks of one run. `drop`: donors left out of the frames and the data.
        build <- function(drop) {
            Jx <- J[!(J$donor_key %in% drop), ]
            D19 <- lapply(stats::setNames(P12_ARMS, P12_ARMS), function(a)
                sort(unique(Jx$donor_key[Jx$arm == a & (Jx$joined_histology | Jx$joined_cps)])))
            DHU <- lapply(stats::setNames(P12_ARMS, P12_ARMS), function(a)
                sort(unique(Jx$donor_key[Jx$arm == a & Jx$joined_histology])))
            tasks <- list(); meta <- list(); col <- 1L
            push <- function(t, panel, a, lv, q, m, grp, dons) {
                t$col <- col
                tasks[[length(tasks) + 1L]] <<- t
                meta[[length(meta) + 1L]] <<- data.frame(panel = panel, arm = a, level = lv, q = q,
                                                         measure = m, group = grp,
                                                         n_units = t$n_units, n_donors = t$n_donors,
                                                         col = col, stringsAsFactors = FALSE)
                col <<- col + if (t$kind == "duel") 2L else 1L
            }
            for (a in P12_ARMS) {
                u <- Jx[Jx$arm == a, ]
                for (q in QS) for (m in IBA)
                    push(p12_17_task("rs", list(x = u[[q]], y = u[[m]]), u$donor_key, LA, D19[[a]], MIND, MINU),
                         "17A", a, LA, q, m, paste("J19", a, sep = "|"), D19[[a]])
                for (q in QS) for (m in CPS)
                    push(p12_17_task("rs", list(x = u[[q]], y = u[[m]]), u$donor_key, LB, D19[[a]], MIND, MINU),
                         "17B", a, LB, q, m, paste("J19", a, sep = "|"), D19[[a]])
                for (m in MEAS)
                    push(p12_17_task("duel", list(y = u[[m]], x = u$module_post_pct, z = u$mito_log1p),
                                     u$donor_key, LC, DHU[[a]], MIND, MINU),
                         "17C", a, LC, "duel", m, paste("JHU", a, sep = "|"), DHU[[a]])
            }
            MT <- do.call(rbind, meta)
            DBG <- c(stats::setNames(D19, paste("J19", names(D19), sep = "|")),
                     stats::setNames(DHU, paste("JHU", names(DHU), sep = "|")))
            OBS <- lapply(tasks, function(t) if (t$testable) p12_17_stat(t, seq_along(t$d)) else
                rep(NA_real_, if (t$kind == "duel") 2L else 1L))
            list(tasks = tasks, MT = MT, D = DBG, OBS = OBS, ncol = col - 1L)
        }

        # Stats rows from one run.
        rows_of <- function(R, BS, pfx) {
            famfun <- if (pfx) b5_family else identity
            pan <- function(p) if (pfx) b5_panel(p) else p
            A <- list(); Bt <- list(); C <- list(); V <- list()
            ci_of <- function(k, off = 0L) if (is.null(BS)) c(lo = NA_real_, hi = NA_real_) else
                p12_boot_ci(BS[, R$MT$col[k] + off])
            p_of <- function(k, off = 0L) if (is.null(BS)) NA_real_ else p12_boot_p(BS[, R$MT$col[k] + off])
            for (k in seq_len(nrow(R$MT))) {
                mt <- R$MT[k, ]; ob <- R$OBS[[k]]; ok <- all(is.finite(ob))
                if (mt$panel %in% c("17A", "17B")) {
                    ci <- if (ok) ci_of(k) else c(lo = NA_real_, hi = NA_real_)
                    row <- p12_row(pan(mt$panel), mt$arm, mt$level, mt$q, mt$measure, mt$n_units, mt$n_donors,
                                   ob, ci[["lo"]], ci[["hi"]], if (ok) p_of(k) else NA,
                                   famfun(paste(mt$panel, mt$arm, mt$level, mt$q, sep = "|")))
                    if (mt$panel == "17A") A[[length(A) + 1L]] <- row else Bt[[length(Bt) + 1L]] <- row
                } else {
                    t <- R$tasks[[k]]
                    cm <- if (ok) ci_of(k) else c(lo = NA, hi = NA)
                    ct <- if (ok) ci_of(k, 1L) else c(lo = NA, hi = NA)
                    C[[length(C) + 1L]] <- p12_row(pan("17C"), mt$arm, mt$level, "beta_module", mt$measure,
                                                   mt$n_units, mt$n_donors, ob[1], cm[["lo"]], cm[["hi"]],
                                                   if (ok) p_of(k) else NA,
                                                   famfun(paste("17C", mt$arm, mt$level, "module", sep = "|")))
                    C[[length(C) + 1L]] <- p12_row(pan("17C"), mt$arm, mt$level, "beta_mito", mt$measure,
                                                   mt$n_units, mt$n_donors, ob[2], ct[["lo"]], ct[["hi"]],
                                                   if (ok) p_of(k, 1L) else NA,
                                                   famfun(paste("17C", mt$arm, mt$level, "mito", sep = "|")))
                    r12 <- rsm <- rst <- NA_real_
                    if (t$testable) { r12 <- p12_rs(t$x, t$z); rsm <- p12_rs(t$x, t$y); rst <- p12_rs(t$z, t$y) }
                    vif <- if (is.finite(r12) && abs(r12) < 1) 1 / (1 - r12^2) else if (is.finite(r12)) Inf else NA
                    for (qq in c("vif", "rs_module_marginal", "rs_mito_marginal"))
                        C[[length(C) + 1L]] <- p12_row(pan("17C"), mt$arm, mt$level, qq, mt$measure,
                                                       mt$n_units, mt$n_donors,
                                                       c(vif = vif, rs_module_marginal = rsm,
                                                         rs_mito_marginal = rst)[[qq]],
                                                       family = famfun("17C|descriptive"))
                    V[[length(V) + 1L]] <- data.frame(arm = mt$arm, level = mt$level, measure = mt$measure,
                                                      n_donors = mt$n_donors, fitted = ok, vif = vif,
                                                      beta_module = ob[1], beta_mito = ob[2],
                                                      stringsAsFactors = FALSE)
                }
            }
            A <- p12_bh_by_family(do.call(rbind, A)); Bt <- p12_bh_by_family(do.call(rbind, Bt))
            C <- p12_bh_by_family(do.call(rbind, C)); V <- do.call(rbind, V)
            gf <- function(q) C$fdr[match(paste(V$arm, V$level, V$measure, q),
                                          paste(C$arm, C$level, C$measure, C$quantity))]
            V$fdr_module <- gf("beta_module"); V$fdr_mito <- gf("beta_mito")
            V$verdict <- p12_17_verdict(V$fdr_module, V$fdr_mito, V$vif, V$fitted)
            list(A = A, B = Bt, C = C, V = V)
        }

        p12_say("SECTION B -- all-donor reproduction")
        R0 <- build(character(0))
        X0 <- rows_of(R0, NULL, FALSE)
        b5_reproduce("17A", S17A[S17A$level == LA, ], X0$A)
        b5_reproduce("17B", S17B[S17B$level == LB, ], X0$B)
        b5_reproduce("17C", S17C[S17C$level == LC, ], X0$C)
        p12_say("")

        p12_say(sprintf("SECTION C -- without the eleven, donor bootstrap per frame (B = %d)", B))
        RW <- build(sa)
        # The same statistics from the all-donor tasks less the SA rows.
        sub_obs <- lapply(seq_along(R0$tasks), function(k) {
            t <- R0$tasks[[k]]; g <- R0$MT$group[k]
            keys <- R0$D[[g]][t$d]
            if (!t$testable) return(rep(NA_real_, if (t$kind == "duel") 2L else 1L))
            p12_17_stat(t, which(!(keys %in% sa)))
        })
        a1 <- unlist(sub_obs); a2 <- unlist(RW$OBS)
        p12_chk(length(a1) == length(a2) && isTRUE(all.equal(a1, a2, tolerance = 1e-12)),
                "tasks rebuilt without the eleven equal the all-donor tasks less their rows",
                sprintf("%d values", length(a1)))
        t0 <- Sys.time()
        BS <- b5_17_boot(RW$tasks, RW$MT$group, RW$D, b5_seed(STEP), B, cl, RW$ncol)
        p12_note("bootstrap time", sprintf("%.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
        XW <- rows_of(RW, BS, TRUE)
        p12_say("")

        p12_say(sprintf("SECTION D -- random removal (%d draws per frame)", R_NULL))
        keys <- lapply(seq_along(R0$tasks), function(k) R0$D[[R0$MT$group[k]]][R0$tasks[[k]]$d])
        t0 <- Sys.time()
        chunks <- p12_chunks(R_NULL, if (is.null(cl)) 1L else length(cl))
        res <- if (is.null(cl)) lapply(chunks, b5_17_null_chunk, tasks = R0$tasks, keys = keys,
                                       draws = SA_$draws, frame_of = R0$MT$group, ncol_out = R0$ncol,
                                       stat = p12_17_stat) else
            parallel::parLapply(cl, chunks, b5_17_null_chunk, tasks = R0$tasks, keys = keys,
                                draws = SA_$draws, frame_of = R0$MT$group, ncol_out = R0$ncol,
                                stat = p12_17_stat)
        NB <- do.call(rbind, res)   # draws x columns
        p12_note("random-removal time", sprintf("%.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
        p12_say("")

        p12_say("SECTION E -- the reading")
        RD <- list(); NL <- list()
        add <- function(src, S, X, k_col_off) {
            for (i in seq_len(nrow(X))) {
                x <- X[i, ]
                j <- match(paste(x$arm, x$level, x$quantity, x$measure),
                           paste(S$arm, S$level, S$quantity, S$measure))
                if (src == "17C") {
                    k <- which(R0$MT$panel == "17C" & R0$MT$arm == x$arm & R0$MT$measure == x$measure)
                    off <- if (x$quantity == "beta_mito") 1L else 0L
                } else {
                    k <- which(R0$MT$panel == src & R0$MT$arm == x$arm & R0$MT$q == x$quantity &
                                   R0$MT$measure == x$measure)
                    off <- 0L
                }
                g <- R0$MT$group[k]
                nv <- NB[, R0$MT$col[k] + off]
                if (fstat(g) != "testable") nv[] <- NA_real_
                RD[[length(RD) + 1L]] <<- data.frame(
                    row_id = paste(src, x$arm, x$level, x$quantity, x$measure, sep = "|"), panel = src,
                    arm = x$arm, level = x$level, quantity = x$quantity, measure = x$measure, frame = g,
                    k_removed = fk(g), est_all = S$estimate[j], fdr_all = S$fdr[j],
                    n_donors_all = S$n_donors[j], est_wo = x$estimate, lo_wo = x$lo, hi_wo = x$hi,
                    fdr_wo = x$fdr, n_donors_wo = x$n_donors, tested_all = is.finite(S$p[j]),
                    tested_wo = is.finite(x$p),
                    null_ok = fstat(g) != "untestable: pool smaller than k", stringsAsFactors = FALSE)
                NL[[length(NL) + 1L]] <<- nv
            }
        }
        add("17A", S17A[S17A$level == LA, ], XW$A)
        add("17B", S17B[S17B$level == LB, ], XW$B)
        add("17C", S17C[S17C$level == LC, ], XW$C[XW$C$quantity %in% c("beta_module", "beta_mito"), ])
        RDt <- b5_readout(do.call(rbind, RD), do.call(rbind, NL))
        N <- do.call(rbind, NL)
        b5_log_classes(transform(RDt, pa = paste(panel, arm)), "pa")

        # 17C verdicts, all donors against without the eleven
        V0 <- X0$V
        gs <- function(q) S17C$fdr[match(paste(V0$arm, V0$level, V0$measure, q),
                                         paste(S17C$arm, S17C$level, S17C$measure, S17C$quantity))]
        V0$fdr_module <- gs("beta_module"); V0$fdr_mito <- gs("beta_mito")
        V0$verdict <- p12_17_verdict(V0$fdr_module, V0$fdr_mito, V0$vif, V0$fitted)
        VW <- XW$V
        VV <- data.frame(arm = V0$arm, level = V0$level, measure = V0$measure,
                         n_donors_all = V0$n_donors, vif_all = V0$vif, verdict_all = V0$verdict,
                         n_donors_wo = VW$n_donors, vif_wo = VW$vif, verdict_wo = VW$verdict,
                         stringsAsFactors = FALSE)
        p12_chk(identical(paste(V0$arm, V0$measure), paste(VW$arm, VW$measure)),
                "17C verdict tables are in the same model order", sprintf("%d models", nrow(VV)))
        for (a in P12_ARMS) {
            g <- VV[VV$arm == a, ]
            p12_say(sprintf("  17C %-6s verdicts kept %d of %d | %s", a, sum(g$verdict_all == g$verdict_wo),
                            nrow(g), paste(sprintf("%s -> %s", g$verdict_all, g$verdict_wo)[
                                g$verdict_all != g$verdict_wo], collapse = "; ")))
        }
        p12_say("")

        p12_say("SECTION F -- write")
        if (P12$n[["FAIL"]] > 0L) p12_abort("checks failed; stats not written")
        p12_write_stats(p12_finalise(XW$A), F_A)
        p12_write_stats(p12_finalise(XW$B), F_B)
        p12_write_stats(p12_finalise(XW$C), F_C)
        b5_write_table(VV, F_V)
        b5_write_table(RDt, F_R)
        b5_write_table(b5_null_long(RDt, N), F_N)
        p12_chk(TRUE, "stats and readout written",
                sprintf("17A %d | 17B %d | 17C %d | readout %d rows", nrow(XW$A), nrow(XW$B), nrow(XW$C),
                        nrow(RDt)))
    }

    p12_run(main)
}

# =============================================================================
# SECTION 29: 19B and 19C, 15A, 17D and 18A without the eleven
# cov4 OLS as sections 07, 20 and 22; random removal refits by lm.fit,
# checked against lm() on the without-eleven run.
# =============================================================================

# == Model fits copied from sections 07, 20 and 22 ============================
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

p12_fit_duel <- function(d, y, covs, min_n) {
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
    list(n = nrow(d), fitted = TRUE,
         b_m = cf[".m", "Estimate"], lo_m = ci[".m", 1], hi_m = ci[".m", 2], p_m = cf[".m", "Pr(>|t|)"],
         b_t = cf[".t", "Estimate"], lo_t = ci[".t", 1], hi_t = ci[".t", 2], p_t = cf[".t", "Pr(>|t|)"],
         vif = unname(p12_vif(cbind(d$.m, d$.t))[1]))
}

p12_fit_18a <- function(d, y, covs, sec, min_n) {
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

# One model for the fast refit: complete-case rows, outcome, covariates and
# predictors in formula order (predictors standardised on each refit).
b5_spec <- function(d, y, covs, preds, check, min_n, out_cols = NULL) {
    cols <- unique(c(y, preds, covs))
    d <- d[stats::complete.cases(d[, cols, drop = FALSE]), , drop = FALSE]
    list(key = d$donor_key, y = as.numeric(d[[y]]),
         C = as.matrix(d[, covs, drop = FALSE]), P = as.matrix(d[, preds, drop = FALSE]),
         check = check, min_n = min_n, out = out_cols)
}

# Standardised predictor coefficients on the rows `keep`, by lm.fit.
b5_fast_fit <- function(sp, keep) {
    na <- rep(NA_real_, ncol(sp$P))
    if (sum(keep) < sp$min_n) return(na)
    y <- sp$y[keep]
    if (!isTRUE(stats::sd(y) > 0)) return(na)
    P <- sp$P[keep, , drop = FALSE]
    sdp <- apply(P, 2L, stats::sd)
    if (!all(is.finite(sdp[sp$check]) & sdp[sp$check] > 0)) return(na)
    P <- sweep(sweep(P, 2L, colMeans(P)), 2L, sdp, "/")
    X <- cbind(1, sp$C[keep, , drop = FALSE], P)
    f <- stats::lm.fit(X, (y - mean(y)) / stats::sd(y))
    cf <- f$coefficients[(ncol(X) - ncol(P) + 1L):ncol(X)]
    ifelse(is.na(cf), NA_real_, cf)
}

# Random-removal refits for one chunk of draws. Returns draws x outputs.
b5_clin_null_chunk <- function(rs, specs, draws, frame_of, n_out, fit) {
    out <- matrix(NA_real_, length(rs), n_out)
    for (i in seq_along(rs)) {
        r <- rs[i]
        for (j in seq_along(specs)) {
            dm <- draws[[frame_of[j]]]
            if (nrow(dm) < r) next
            sp <- specs[[j]]
            out[i, sp$out] <- fit(sp, !(sp$key %in% dm[r, ]))
        }
    }
    out
}

p12_section_29 <- function() {
    STEP <- "29"
    D <- p12_step_dirs(STEP, "19_Clinical", parent = FIG_DIR)
    p12_clean(D$out)
    p12_start(STEP, SCRIPT, D, pkgs = "parallel")
    b5_header()
    F_15 <- file.path(D$out, "Stats_19_15A.csv")
    F_17 <- file.path(D$out, "Stats_19_17D.csv")
    F_18 <- file.path(D$out, "Stats_19_18A.csv")
    F_V17 <- file.path(D$out, "Verdicts_19_17D.csv")
    F_V18 <- file.path(D$out, "Verdicts_19_18A.csv")
    F_R <- file.path(D$out, "Readout_19_Clinical.csv")
    F_N <- file.path(D$out, "Null_19_Clinical.csv")
    F18 <- P12_CFG$figure18

    main <- function() {
        MIND <- as.integer(P12_CFG$floors$donors_per_test)
        MINC <- as.integer(P12_CFG$floors$cells_per_group_per_donor)
        DET <- as.numeric(P12_CFG$nuclear_panels$detect_min)
        MING <- as.integer(F18$min_panel_genes)
        R_NULL <- as.integer(B5$null$draws)
        CL <- P12_CFG$clinical
        COVS <- unlist(CL$cov4)
        OUT <- vapply(CL$outcomes, function(o) o$column, "")
        OCLASS <- vapply(CL$outcomes, function(o) o$class, "")
        PAN <- F18$panels
        PK <- vapply(PAN, function(p) p$key, "")
        SECOND <- c(module = "module_post_pct", mito = "mito_log1p")
        b5_reading("Battery 5 definitions")
        b5_reading("19B")
        b5_reading("19C")

        p12_say("SECTION A -- inputs")
        DN <- b5_donors()
        miss <- setdiff(c(COVS, OUT, "module_post_pct", "mito_log1p", "contam_super_nonmodule"), names(DN))
        if (length(miss)) p12_abort(paste("donor table lacks:", paste(miss, collapse = ", ")))
        L <- utils::read.csv(p12_upstream_file("04_Labels", "04_cell_labels.csv.gz"), check.names = FALSE,
                             colClasses = c(barcode = "character", donor_key = "character"))
        X <- methods::as(Matrix::readMM(p12_upstream_file("02_Extract", "02_counts.mtx.gz")),
                         "CsparseMatrix")
        genes <- readLines(p12_upstream_file("02_Extract", "02_counts_genes.tsv"))
        cells <- readLines(p12_upstream_file("02_Extract", "02_counts_cells.tsv"))
        SETS <- utils::read.csv(p12_upstream_file("02_Extract", "02_matched_sets.csv"),
                                stringsAsFactors = FALSE)
        p12_chk(all(L$labels_version == P12_CFG$labels$labels_version),
                "labels carry this labels version", "checked")
        p12_chk(nrow(X) == length(cells) && ncol(X) == length(genes),
                "count matrix, gene list and cell list agree", sprintf("%d x %d", nrow(X), ncol(X)))
        j0 <- match(L$barcode, cells)
        p12_chk(!anyNA(j0), "every labelled cell is in the count matrix", sprintf("%d cells", nrow(L)))
        SA_ <- b5_sa_inputs(); sa <- SA_$sa
        S15 <- b5_source_stats("15A"); S17 <- b5_source_stats("17D"); S18 <- b5_source_stats("18A")
        fr <- SA_$frames
        fstat <- function(g) fr$null_status[fr$frame == g]
        fk <- function(g) fr$k[fr$frame == g]
        cl <- p12_cluster(R_NULL)
        if (!is.null(cl)) on.exit(parallel::stopCluster(cl), add = TRUE)
        if (P12$n[["FAIL"]] > 0L) p12_abort("input checks failed; nothing computed")
        p12_say("")

        p12_say("SECTION B -- donor panel scores (section 22, all donors, held fixed)")
        scores <- list()
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
                key <- paste(a, ds, PK[k], sep = "|")
                if (sum(det) < MING || !length(okd)) { scores[[key]] <- NULL; next }
                Md <- M[, det, drop = FALSE]
                Md <- Matrix::Diagonal(x = 1e4 / tot) %*% Md
                Md@x <- log1p(Md@x)
                v <- Matrix::rowMeans(Md)
                scores[[key]] <- tapply(v[don %in% okd], don[don %in% okd], stats::median)
            }
        }
        p12_note("panel scores", sprintf("%d of %d arm, region set and panel combinations scored",
                                         length(scores), length(P12_ARMS) * length(P12_DATASETS) * length(PAN)))
        rm(X); invisible(gc())
        p12_say("")

        # All three panels for one donor set. `drop`: donors removed.
        run <- function(drop, pfx) {
            pan <- function(p) if (pfx) b5_panel(p) else p
            famf <- if (pfx) b5_family else identity
            r15 <- list(); r17 <- list(); v17 <- list(); r18 <- list(); v18 <- list()
            for (a in P12_ARMS) for (ds in P12_DATASETS) {
                d <- DN[DN$arm == a & DN$dataset == ds & !(DN$donor_key %in% drop), ]
                for (k in names(OUT)) {
                    y <- OUT[[k]]
                    cv <- if (y %in% COVS) setdiff(COVS, y) else COVS
                    u <- p12_fit_clinical(d, y, "module_post_pct", cv, min_n = MIND)
                    us <- p12_fit_clinical(d[is.finite(d$contam_super_nonmodule), ], y, "module_post_pct", cv,
                                           min_n = MIND)
                    ad <- p12_fit_clinical(d, y, "module_post_pct", cv, extra = "contam_super_nonmodule",
                                           min_n = MIND)
                    f5 <- function(q, fit, est, lo, hi, p, fam)
                        r15[[length(r15) + 1L]] <<- p12_row(pan("15A"), a, ds, q, k, NA, fit[["n"]], est,
                                                            lo, hi, p, famf(fam))
                    f5("module_beta_unadjusted", u, u[["beta"]], u[["lo"]], u[["hi"]], u[["p"]],
                       paste("15A", a, ds, "unadjusted", sep = "|"))
                    f5("module_beta_contam_adjusted", ad, ad[["beta"]], ad[["lo"]], ad[["hi"]], ad[["p"]],
                       paste("15A", a, ds, "contamination-adjusted", sep = "|"))
                    fd <- paste("15A", a, ds, "dR2", sep = "|")
                    f5("dR2_unadjusted", u, u[["dr2"]], NA, NA, NA, fd)
                    f5("dR2_unadjusted_same_donors", us, us[["dr2"]], NA, NA, NA, fd)
                    f5("dR2_contam_adjusted", ad, ad[["dr2"]], NA, NA, NA, fd)
                    ret <- if (is.finite(us[["dr2"]]) && us[["dr2"]] > 0) ad[["dr2"]] / us[["dr2"]] else NA
                    f5("retention", ad, ret, NA, NA, NA, fd)

                    f <- p12_fit_duel(d, y, cv, MIND)
                    g <- function(z) if (isTRUE(f$fitted)) f[[z]] else NA_real_
                    r17[[length(r17) + 1L]] <- p12_row(pan("17D"), a, ds, "beta_module", k, NA, f$n, g("b_m"),
                                                       g("lo_m"), g("hi_m"), g("p_m"),
                                                       famf(paste("17D", a, ds, "module", sep = "|")))
                    r17[[length(r17) + 1L]] <- p12_row(pan("17D"), a, ds, "beta_mito", k, NA, f$n, g("b_t"),
                                                       g("lo_t"), g("hi_t"), g("p_t"),
                                                       famf(paste("17D", a, ds, "mito", sep = "|")))
                    r17[[length(r17) + 1L]] <- p12_row(pan("17D"), a, ds, "vif", k, NA, f$n, g("vif"),
                                                       family = famf("17D|descriptive"))
                    v17[[length(v17) + 1L]] <- data.frame(arm = a, dataset = ds, outcome = k, n_donors = f$n,
                                                          fitted = isTRUE(f$fitted), vif = g("vif"),
                                                          stringsAsFactors = FALSE)
                }
                for (kk in seq_along(PAN)) {
                    sc <- scores[[paste(a, ds, PK[kk], sep = "|")]]
                    d$.panel <- if (is.null(sc)) NA_real_ else as.numeric(sc[d$donor_key])
                    for (sp in names(SECOND)) for (o in names(OUT)) {
                        y <- OUT[[o]]
                        cv <- if (y %in% COVS) setdiff(COVS, y) else COVS
                        f <- if (is.null(sc)) list(n = 0L, fitted = FALSE) else p12_fit_18a(d, y, cv, SECOND[[sp]], MIND)
                        g <- function(z) if (isTRUE(f$fitted)) f[[z]] else NA_real_
                        fam <- paste("18A", a, ds, PK[kk], sp, sep = "|")
                        ms <- paste(sp, PK[kk], o, sep = "|")
                        r18[[length(r18) + 1L]] <- p12_row(pan("18A"), a, ds, "beta_second", ms, NA, f$n, g("b_s"),
                                                           p = g("p_s"), family = famf(paste(fam, "second", sep = "|")))
                        r18[[length(r18) + 1L]] <- p12_row(pan("18A"), a, ds, "beta_panel", ms, NA, f$n, g("b_n"),
                                                           p = g("p_n"), family = famf(paste(fam, "panel", sep = "|")))
                        r18[[length(r18) + 1L]] <- p12_row(pan("18A"), a, ds, "vif", ms, NA, f$n, g("vif"),
                                                           family = famf("18A|descriptive"))
                        v18[[length(v18) + 1L]] <- data.frame(arm = a, dataset = ds, predictor = sp, panel = PK[kk],
                                                              outcome = o, n_donors = f$n,
                                                              fitted = isTRUE(f$fitted), vif = g("vif"),
                                                              stringsAsFactors = FALSE)
                    }
                }
            }
            A <- p12_bh_by_family(do.call(rbind, r15)); B7 <- p12_bh_by_family(do.call(rbind, r17))
            C8 <- p12_bh_by_family(do.call(rbind, r18))
            V7 <- do.call(rbind, v17); V8 <- do.call(rbind, v18)
            list(A = A, B = B7, C = C8, V7 = V7, V8 = V8)
        }
        # Verdicts from a stats table (FDR) and the model table (fit, VIF).
        verdict17 <- function(S, V) {
            gd <- function(q) S$fdr[match(paste(V$arm, V$dataset, V$outcome, q),
                                          paste(S$arm, S$level, S$measure, S$quantity))]
            p12_17_verdict(gd("beta_module"), gd("beta_mito"), V$vif, V$fitted)
        }
        verdict18 <- function(S, V) {
            k <- function(q) match(paste(V$arm, V$dataset, paste(V$predictor, V$panel, V$outcome, sep = "|"), q),
                                   paste(S$arm, S$level, S$measure, S$quantity))
            p12_18_verdict(S$fdr[k("beta_second")], S$fdr[k("beta_panel")], V$vif, V$fitted)
        }

        p12_say("SECTION C -- all-donor reproduction")
        X0 <- run(character(0), FALSE)
        b5_reproduce("15A", S15, X0$A, with_p = TRUE)
        b5_reproduce("17D", S17, X0$B, with_p = TRUE)
        b5_reproduce("18A", S18, X0$C, with_p = TRUE)
        p12_say("")

        p12_say("SECTION D -- without the eleven")
        XW <- run(sa, TRUE)
        p12_say("")

        p12_say(sprintf("SECTION E -- random removal (%d draws per frame, fast refit)", R_NULL))
        specs <- list(); frame_of <- character(0); rid <- character(0); n_out <- 0L
        addspec <- function(sp, g, ids) {
            sp$out <- n_out + seq_along(ids)
            n_out <<- n_out + length(ids)
            specs[[length(specs) + 1L]] <<- sp
            frame_of <<- c(frame_of, g); rid <<- c(rid, ids)
        }
        for (a in P12_ARMS) for (ds in P12_DATASETS) {
            g <- paste("DN", a, ds, sep = "|")
            d <- DN[DN$arm == a & DN$dataset == ds, ]
            for (k in names(OUT)) {
                y <- OUT[[k]]
                cv <- if (y %in% COVS) setdiff(COVS, y) else COVS
                addspec(b5_spec(d, y, cv, "module_post_pct", 1L, MIND),
                        g, paste("15A", a, ds, "module_beta_unadjusted", k, sep = "|"))
                addspec(b5_spec(d, y, c(cv, "contam_super_nonmodule"), "module_post_pct", 1L, MIND),
                        g, paste("15A", a, ds, "module_beta_contam_adjusted", k, sep = "|"))
                addspec(b5_spec(d, y, cv, c("module_post_pct", "mito_log1p"), 1:2, MIND),
                        g, paste("17D", a, ds, c("beta_module", "beta_mito"), k, sep = "|"))
            }
            for (kk in seq_along(PAN)) {
                sc <- scores[[paste(a, ds, PK[kk], sep = "|")]]
                if (is.null(sc)) next
                d$.panel <- as.numeric(sc[d$donor_key])
                for (sp in names(SECOND)) for (o in names(OUT)) {
                    y <- OUT[[o]]
                    cv <- if (y %in% COVS) setdiff(COVS, y) else COVS
                    addspec(b5_spec(d, y, cv, c(SECOND[[sp]], ".panel"), 1:2, MIND), g,
                            paste("18A", a, ds, c("beta_second", "beta_panel"),
                                  paste(sp, PK[kk], o, sep = "|"), sep = "|"))
                }
            }
        }
        # The fast refit against lm() on the without-eleven run.
        fast_sa <- unlist(lapply(specs, function(sp) b5_fast_fit(sp, !(sp$key %in% sa))))
        XWall <- rbind(transform(XW$A, src = "15A"), transform(XW$B, src = "17D"), transform(XW$C, src = "18A"))
        lm_sa <- XWall$estimate[match(rid, paste(XWall$src, XWall$arm, XWall$level, XWall$quantity,
                                                 XWall$measure, sep = "|"))]
        dfz <- abs(fast_sa - lm_sa); dfz[is.na(fast_sa) & is.na(lm_sa)] <- 0
        p12_chk(all(is.finite(dfz) & dfz < 1e-8),
                "the fast refit equals lm() on the without-eleven run",
                sprintf("%d coefficients | %d differ | largest difference %.1e", length(dfz),
                        sum(!(is.finite(dfz) & dfz < 1e-8)), suppressWarnings(max(dfz, na.rm = TRUE))))
        t0 <- Sys.time()
        chunks <- p12_chunks(R_NULL, if (is.null(cl)) 1L else length(cl))
        res <- if (is.null(cl)) lapply(chunks, b5_clin_null_chunk, specs = specs, draws = SA_$draws,
                                       frame_of = frame_of, n_out = n_out, fit = b5_fast_fit) else
            parallel::parLapply(cl, chunks, b5_clin_null_chunk, specs = specs, draws = SA_$draws,
                                frame_of = frame_of, n_out = n_out, fit = b5_fast_fit)
        NB <- do.call(rbind, res)   # draws x coefficients
        p12_note("random-removal time", sprintf("%.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
        p12_say("")

        p12_say("SECTION F -- the reading")
        tested <- rbind(transform(XW$A[XW$A$quantity %in% c("module_beta_unadjusted", "module_beta_contam_adjusted"), ],
                                  src = "15A"),
                        transform(XW$B[XW$B$quantity %in% c("beta_module", "beta_mito"), ], src = "17D"),
                        transform(XW$C[XW$C$quantity %in% c("beta_second", "beta_panel"), ], src = "18A"))
        SS <- rbind(transform(S15, src = "15A"), transform(S17, src = "17D"), transform(S18, src = "18A"))
        kt <- paste(tested$src, tested$arm, tested$level, tested$quantity, tested$measure, sep = "|")
        j <- match(kt, paste(SS$src, SS$arm, SS$level, SS$quantity, SS$measure, sep = "|"))
        jn <- match(kt, rid)
        fg <- paste("DN", tested$arm, tested$level, sep = "|")
        RDt <- data.frame(row_id = kt, panel = tested$src, arm = tested$arm, level = tested$level,
                          quantity = tested$quantity, measure = tested$measure, frame = fg,
                          k_removed = vapply(fg, fk, 0L), est_all = SS$estimate[j], fdr_all = SS$fdr[j],
                          n_donors_all = SS$n_donors[j], est_wo = tested$estimate, lo_wo = tested$lo,
                          hi_wo = tested$hi, fdr_wo = tested$fdr, n_donors_wo = tested$n_donors,
                          tested_all = is.finite(SS$p[j]), tested_wo = is.finite(tested$p),
                          null_ok = vapply(fg, fstat, "") != "untestable: pool smaller than k",
                          stringsAsFactors = FALSE, row.names = NULL)
        N <- matrix(NA_real_, nrow(RDt), R_NULL)
        has <- !is.na(jn)
        N[has, ] <- t(NB[, jn[has], drop = FALSE])
        N[vapply(fg, fstat, "") != "testable", ] <- NA_real_
        RDt <- b5_readout(RDt, N)
        b5_log_classes(transform(RDt, pa = paste(panel, arm)), "pa")

        # Verdicts, all donors (Battery 3 and 4 FDR) against without the eleven
        V17 <- data.frame(X0$V7[, c("arm", "dataset", "outcome")], n_donors_all = X0$V7$n_donors,
                          vif_all = X0$V7$vif, verdict_all = verdict17(S17, X0$V7),
                          n_donors_wo = XW$V7$n_donors, vif_wo = XW$V7$vif,
                          verdict_wo = verdict17(XW$B, XW$V7), stringsAsFactors = FALSE)
        V18 <- data.frame(X0$V8[, c("arm", "dataset", "predictor", "panel", "outcome")],
                          n_donors_all = X0$V8$n_donors, vif_all = X0$V8$vif,
                          verdict_all = verdict18(S18, X0$V8), n_donors_wo = XW$V8$n_donors,
                          vif_wo = XW$V8$vif, verdict_wo = verdict18(XW$C, XW$V8), stringsAsFactors = FALSE)
        for (nm in c("17D", "18A")) {
            V <- if (nm == "17D") V17 else V18
            for (a in P12_ARMS) {
                g <- V[V$arm == a, ]
                ch <- table(paste(g$verdict_all, "->", g$verdict_wo)[g$verdict_all != g$verdict_wo])
                p12_say(sprintf("  %s %-6s verdicts kept %d of %d%s", nm, a, sum(g$verdict_all == g$verdict_wo),
                                nrow(g), if (length(ch)) paste0(" | ", paste(sprintf("%s: %d", names(ch), ch),
                                                                             collapse = "; ")) else ""))
            }
        }
        p12_say("")

        p12_say("SECTION G -- write")
        if (P12$n[["FAIL"]] > 0L) p12_abort("checks failed; stats not written")
        p12_write_stats(p12_finalise(XW$A), F_15)
        p12_write_stats(p12_finalise(XW$B), F_17)
        p12_write_stats(p12_finalise(XW$C), F_18)
        b5_write_table(V17, F_V17)
        b5_write_table(V18, F_V18)
        b5_write_table(RDt, F_R)
        b5_write_table(b5_null_long(RDt, N), F_N)
        p12_chk(TRUE, "stats and readout written",
                sprintf("15A %d | 17D %d | 18A %d | readout %d rows", nrow(XW$A), nrow(XW$B), nrow(XW$C),
                        nrow(RDt)))
    }

    p12_run(main)
}

# =============================================================================
# SECTION 30: leave-one-donor-out influence on 16A, 16B and 16F
# Each donor removed in turn; no test. Tables only.
# =============================================================================
p12_section_30 <- function() {
    STEP <- "30"
    D <- p12_step_dirs(STEP, "19_Influence", parent = FIG_DIR)
    p12_clean(D$out)
    p12_start(STEP, SCRIPT, D)
    b5_header()
    F_I <- file.path(D$out, "Influence_19.csv")
    F_S <- file.path(D$out, "Influence_19_summary.csv")
    QA <- c("module_post_pct", "module_resid_tq4", "mito_log1p", "contamination", "ctrl_pct_null")
    QB <- c("module_post_pct", "module_resid_tq4", "mito_log1p")
    QS <- c(module = "module_post_pct", mito = "mito_log1p", ncells = "log_n_cells",
            depth = "log_median_umi")
    N_SPLITS <- 200L
    LV <- "between donors"; FIELD <- "intronic"

    main <- function() {
        MIND <- as.integer(P12_CFG$floors$donors_per_test)
        b5_reading("Influence")
        p12_say("SECTION A -- inputs")
        U <- b5_units()
        SA_ <- b5_sa_inputs(); SAt <- SA_$SA
        RD <- utils::read.csv(p12_upstream_file(file.path(FIG_DIR, "19_UnitPanels"), "Readout_19_Units.csv"),
                              check.names = FALSE, stringsAsFactors = FALSE)
        dons <- sort(unique(U$donor_key))
        p12_say("")

        p12_say(sprintf("SECTION B -- %d donors removed in turn", length(dons)))
        # All-donor objects
        SPL <- b5_16b_splits(U, QB, MIND, N_SPLITS)
        T0 <- b5_16f_tasks(U, dons, QS, FIELD, MIND)
        f_rows <- b5_16f_rows(b5_16f_obs(T0, MIND), NULL, T0$meta, "16F", identity, QS, LV, FIELD)
        share_of <- function(a, q, drop) {
            u <- U[U$arm == a & is.finite(U[[q]]) & !(U$donor_key %in% drop), ]
            dd <- as.integer(factor(u$donor_key))
            if (max(c(0L, dd)) >= MIND && isTRUE(stats::var(u[[q]]) > 0)) b5_share(u[[q]], dd) else NA_real_
        }
        est_all <- function(drop) {
            e16a <- unlist(lapply(P12_ARMS, function(a) vapply(QA, function(q) share_of(a, q, drop), 0)))
            e16b <- unlist(lapply(P12_ARMS, function(a) vapply(QB, function(q) {
                sp <- SPL[[paste(a, q)]]
                if (is.null(sp$A1)) return(NA_real_)
                kp <- which(!(sp$dons %in% drop))
                if (length(kp) >= MIND) b5_splithalf(sp$A1, sp$A2, kp) else NA_real_
            }, 0)))
            e16f <- b5_16f_est(b5_16f_obs(T0, MIND, drop), T0$meta, QS)
            c(e16a, e16b, e16f)
        }
        KEY <- rbind(
            data.frame(panel = "16A", arm = rep(P12_ARMS, each = length(QA)), level = "units",
                       quantity = rep(QA, length(P12_ARMS)), measure = "donor share of variance"),
            data.frame(panel = "16B", arm = rep(P12_ARMS, each = length(QB)), level = "donors",
                       quantity = rep(QB, length(P12_ARMS)), measure = "split-half rs, median over splits"),
            data.frame(panel = "16F", arm = f_rows$arm, level = LV, quantity = f_rows$quantity, measure = FIELD))
        KEY$row_id <- paste(KEY$panel, KEY$arm, KEY$level, KEY$quantity, KEY$measure, sep = "|")
        E0 <- est_all(character(0))
        j <- match(KEY$row_id, RD$row_id)
        d0 <- abs(E0 - RD$est_all[j]); d0[is.na(E0) & is.na(RD$est_all[j])] <- 0
        p12_chk(!anyNA(j) && all(is.finite(d0) & d0 < 1e-8),
                "all-donor estimates equal section 27's readout",
                sprintf("%d rows | largest difference %.1e", nrow(KEY), suppressWarnings(max(d0, na.rm = TRUE))))
        # The donors each row is computed from.
        tk <- function(a, q) T0$tasks[[which(T0$meta$arm == a & T0$meta$q == q)]]$key
        FRM <- lapply(seq_len(nrow(KEY)), function(i) {
            a <- KEY$arm[i]; q <- KEY$quantity[i]
            switch(KEY$panel[i],
                   "16A" = unique(U$donor_key[U$arm == a & is.finite(U[[q]])]),
                   "16B" = SPL[[paste(a, q)]]$dons,
                   "16F" = switch(q,
                                  module_minus_mito = union(tk(a, "module"), tk(a, "mito")),
                                  module_minus_depth = union(tk(a, "module"), tk(a, "depth")),
                                  arm_minus_hkg8 = union(tk(a, "module"), tk("HKG8", "module")),
                                  tk(a, q)))
        })
        rows <- list()
        for (dk in dons) {
            E <- est_all(dk)
            inf <- vapply(FRM, function(f) dk %in% f, TRUE)
            rows[[length(rows) + 1L]] <- data.frame(KEY[inf, ], donor_key = dk,
                                                    est_all = E0[inf], est_minus_donor = E[inf],
                                                    change = b5_change(E[inf], E0[inf]),
                                                    stringsAsFactors = FALSE, row.names = NULL)
        }
        IN <- do.call(rbind, rows)
        IN$severely_affected <- IN$donor_key %in% SA_$sa
        IN$high_adnc <- IN$donor_key %in% SA_$high
        IN <- IN[order(IN$row_id, IN$change), ]
        IN$rank_by_change <- stats::ave(IN$change, IN$row_id, FUN = function(v) rank(v, ties.method = "min",
                                                                                   na.last = "keep"))
        p12_say("")

        p12_say("SECTION C -- summary per row")
        SM <- do.call(rbind, lapply(split(IN, IN$row_id), function(g) {
            s <- g$severely_affected; h <- g$high_adnc & !s
            k <- sum(s)
            data.frame(row_id = g$row_id[1], panel = g$panel[1], arm = g$arm[1], quantity = g$quantity[1],
                       n_donors = nrow(g), n_sa = k,
                       sa_ranks = paste(sort(g$rank_by_change[s]), collapse = ";"),
                       sa_in_top_k = if (k) sum(g$rank_by_change[s] <= k, na.rm = TRUE) else NA_integer_,
                       median_change_sa = if (k) stats::median(g$change[s], na.rm = TRUE) else NA_real_,
                       median_change_high_non_sa = if (any(h)) stats::median(g$change[h], na.rm = TRUE) else NA_real_,
                       median_change_other = if (any(!s & !h)) stats::median(g$change[!s & !h], na.rm = TRUE) else NA_real_,
                       sum_change_sa = if (k) sum(g$change[s], na.rm = TRUE) else NA_real_,
                       stringsAsFactors = FALSE)
        }))
        SM$change_sa_together <- RD$delta_sa[match(SM$row_id, RD$row_id)]
        rownames(SM) <- NULL
        for (i in which(SM$quantity == "module_post_pct" | SM$quantity == "module"))
            p12_say(sprintf("  %-45s SA %2d | SA ranks %s | median change SA %+.4f, other High %+.4f, rest %+.4f",
                            SM$row_id[i], SM$n_sa[i], SM$sa_ranks[i], SM$median_change_sa[i],
                            SM$median_change_high_non_sa[i], SM$median_change_other[i]))
        p12_say("")

        p12_say("SECTION D -- write")
        if (P12$n[["FAIL"]] > 0L) p12_abort("checks failed; tables not written")
        b5_write_table(IN, F_I)
        b5_write_table(SM, F_S)
        p12_chk(TRUE, "tables written", sprintf("influence %d rows | summary %d rows", nrow(IN), nrow(SM)))
    }

    p12_run(main)
}

# =============================================================================
# RUN, in order; stop at the first section that does not PASS
# =============================================================================
local({
    # Clear every Figure 19 section folder first, so no stale output remains.
    for (d in list.dirs(file.path(P12_RERUN, FIG_DIR), recursive = FALSE)) p12_clean(d)
    for (code in c("26", "27", "28", "29", "30")) {
        v <- get(paste0("p12_section_", code))()
        if (!identical(v, "PASS"))
            stop(sprintf("section %s ended with verdict %s; later sections not run. See its log.",
                         code, v), call. = FALSE)
    }
    message("Figure 19 analysis: sections 26 to 30 all PASS. Next: run Part12_Fig19_Plots.r")
})
