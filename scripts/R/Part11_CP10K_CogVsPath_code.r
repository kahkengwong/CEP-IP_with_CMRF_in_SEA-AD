# =============================================================================
# Part11_CP10K_CogVsPath_code.r
# SEA-AD CEP-IP/CMRF, Part 11, cognition against pathology (Figure 9D)
#   Compares Module and non-Module donor percentages on cognition (MMSE, CASI)
#   and pathology (ADNC, CERAD, Braak, Thal), from each arm's Part 10h
#   results: the contrast b(Module) - b(non-Module) in one model (section B),
#   the shared-component ratio (A), the cognition-to-pathology balance (C) and
#   the contrast on each outcome's own donors (D).
#
# Usage (from the SEA-AD folder, or set SEAAD_ROOT), after Part 10h
#   source("Part11_CP10K_CogVsPath_code.r")
#
# Outputs (Part11_Identity/Part11_CogVsPath/): figures (PDF, PNG), workbook,
#   results RDS and log.
# =============================================================================
suppressPackageStartupMessages({
    library(ggplot2)
})

# == CONFIGURATION ============================================================
# READ_ROOT holds the Part 10h outputs; OUT_ROOT is where this part writes.
if (!exists("READ_ROOT")) READ_ROOT <- Sys.getenv("SEAAD_ROOT", getwd())
if (!exists("OUT_ROOT"))  OUT_ROOT  <- READ_ROOT

OUT_DIR     <- file.path(OUT_ROOT, "Part11_Identity", "Part11_CogVsPath")
FILE_PREFIX <- "Part11_CP10K"

# Part 10h results file per arm; the untagged main-arm name is a fallback.
ARM_PATTERNS <- c(PVM3_all     = "_PVM_3_",
                  PVM1_all     = "_PVM_1_",
                  HKG8_matched = "_HKG8_")
MAIN_ARM       <- "Mic13_PVM3"
MAIN_BASENAME  <- "Part10h_CP10K_M4_Prognosis_results.rds"
MAIN_BASE_LEGACY <- "Part10h_CP10K_Prognosis_results.rds"

ARM_ORDER <- c("Mic13_PVM3", "PVM3_all", "PVM1_all", "HKG8_matched")
ARM_LABEL <- c(Mic13_PVM3   = "Mic13 PVM_3 (M4-top)",
               PVM3_all     = "PVM_3, all cells",
               PVM1_all     = "PVM_1, all cells",
               HKG8_matched = "HKG8 control")
CONTROL_ARM <- "HKG8_matched"

DS_ORDER <- c("NeoC", "NonNeoC")

COL_MOD <- "pct_post"        # module Post-IP
COL_NM  <- "pct_excl_post"   # non-module Post-IP
COV4    <- c("Age_num", "APOE_e4_dose", "Sex_num", "PMI_num")

# orient makes "more disease" positive on every outcome.
OUTCOMES <- data.frame(
    col    = c("MMSE_num", "CASI_num",
               "ADNC_rank", "CERAD_num", "Braak_num", "Thal_num"),
    label  = c("MMSE", "CASI", "ADNC", "CERAD", "Braak", "Thal"),
    domain = c("cognition", "cognition",
               "pathology", "pathology", "pathology", "pathology"),
    orient = c(-1, -1, 1, 1, 1, 1),
    stringsAsFactors = FALSE)

N_BOOT      <- 2000L
BOOT_SEED   <- 20260810L
MIN_DONORS  <- 30L

# Section A: no ratio when |denominator| is below this floor.
DENOM_FLOOR <- 0.10
# Section C: floor on the mean pathology beta.
PATH_FLOOR  <- 0.05

CI_LO <- 0.025
CI_HI <- 0.975
POWER_TARGET <- 0.80

# Figure colours.
COG_COL   <- "#A63A5A"
PATH_COL  <- "#1F3864"
MODULE_COL <- "#1F3864"
NONMOD_COL <- "#9AA5B1"
NULL_COL  <- "grey45"
# Section A intervals are clipped to this range on the figure.
SCR_CLIP  <- c(-250, 500)

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)


# == GUARDRAIL LEDGER =========================================================
.gr <- list()
gr11 <- function(id, status, detail) {
    if (length(detail) != 1L || is.na(detail))
        detail <- paste("[detail unavailable: the message was built from an",
                        "empty value]")
    .gr[[length(.gr) + 1L]] <<- list(id = id, status = status, detail = detail)
    cat(sprintf("  [%-4s]  %-12s %s\n", status, id, detail))
}

cat("\n", strrep("=", 78), "\n", sep = "")
cat("  Part 11 -- cognition versus pathology, module against non-module\n")
cat("  ", format(Sys.time()), "\n", sep = "")
cat("  Read : ", READ_ROOT, "\n", sep = "")
cat("  Write: ", OUT_DIR, "\n", sep = "")
cat(strrep("=", 78), "\n", sep = "")


# == HELPERS ==================================================================
zscore <- function(x) {
    s <- stats::sd(x)
    if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
    (x - mean(x)) / s
}

.cov_rhs <- paste(paste0("`", COV4, "`"), collapse = " + ")

# One marginal model: outcome on cov4 plus one standardised percentage.
.beta_marginal <- function(d, ycol, orient, pcol) {
    dd    <- d
    dd$.y <- orient * zscore(dd[[ycol]])
    dd$.p <- zscore(dd[[pcol]])
    if (all(is.na(dd$.y)) || all(is.na(dd$.p))) return(NA_real_)
    fit <- try(stats::lm(stats::as.formula(paste0(".y ~ ", .cov_rhs, " + .p")),
                         data = dd), silent = TRUE)
    if (inherits(fit, "try-error")) return(NA_real_)
    cf <- stats::coef(fit)
    if (!".p" %in% names(cf)) return(NA_real_)
    unname(cf[".p"])
}

# One joint model: both percentages together, plus the exact standard error of
# their difference from the fitted covariance matrix.
.joint <- function(d, ycol, orient) {
    dd      <- d
    dd$.y   <- orient * zscore(dd[[ycol]])
    dd$.mod <- zscore(dd[[COL_MOD]])
    dd$.nm  <- zscore(dd[[COL_NM]])
    if (all(is.na(dd$.y)) || all(is.na(dd$.mod)) || all(is.na(dd$.nm)))
        return(NULL)
    fit <- try(stats::lm(
        stats::as.formula(paste0(".y ~ ", .cov_rhs, " + .mod + .nm")),
        data = dd), silent = TRUE)
    if (inherits(fit, "try-error")) return(NULL)
    cf <- stats::coef(fit)
    if (!all(c(".mod", ".nm") %in% names(cf))) return(NULL)
    V  <- stats::vcov(fit)
    se <- sqrt(V[".mod", ".mod"] + V[".nm", ".nm"] - 2 * V[".mod", ".nm"])
    df <- stats::df.residual(fit)
    dl <- unname(cf[".mod"] - cf[".nm"])
    tv <- NA_real_
    if (is.finite(se) && se > 0) tv <- dl / se
    pv <- NA_real_
    if (is.finite(tv)) pv <- 2 * stats::pt(-abs(tv), df)
    list(b_mod = unname(cf[".mod"]), b_nm = unname(cf[".nm"]),
         delta = dl, se = se, df = df, p = pv)
}

# Section A on one donor frame.
.scr <- function(d, ycol, orient) {
    b_mod <- .beta_marginal(d, ycol, orient, COL_MOD)
    b_nm  <- .beta_marginal(d, ycol, orient, COL_NM)
    r     <- suppressWarnings(stats::cor(d[[COL_MOD]], d[[COL_NM]]))
    den   <- r * b_mod
    ratio <- NA_real_
    if (is.finite(den) && abs(den) >= DENOM_FLOOR && is.finite(b_nm))
        ratio <- b_nm / den
    list(b_mod = b_mod, b_nm = b_nm, r = r, denom = den, ratio = ratio)
}

# Section C on one donor frame.
.balance <- function(d) {
    one <- function(pcol) {
        b <- rep(NA_real_, nrow(OUTCOMES))
        for (k in seq_len(nrow(OUTCOMES)))
            b[k] <- .beta_marginal(d, OUTCOMES$col[k], OUTCOMES$orient[k],
                                   pcol)
        cog  <- mean(b[OUTCOMES$domain == "cognition"], na.rm = TRUE)
        path <- mean(b[OUTCOMES$domain == "pathology"], na.rm = TRUE)
        ratio <- NA_real_
        if (is.finite(path) && path >= PATH_FLOOR) ratio <- cog / path
        list(cog = cog, path = path, ratio = ratio)
    }
    list(module = one(COL_MOD), nonmodule = one(COL_NM))
}

.q <- function(v) {
    ok <- v[is.finite(v)]
    if (length(ok) < 200L) return(c(NA_real_, NA_real_))
    unname(stats::quantile(ok, c(CI_LO, CI_HI)))
}


# == STEP 1: locate the Part 10h results ======================================
cat("\nSTEP 1 -- Part 10h results objects\n")
all_rds <- list.files(READ_ROOT, pattern = "_Prognosis_results\\.rds$",
                      full.names = TRUE, recursive = TRUE)
arm_file <- stats::setNames(rep(NA_character_, length(ARM_ORDER)), ARM_ORDER)
claimed  <- character(0)
for (a in names(ARM_PATTERNS)) {
    hit <- all_rds[grepl(ARM_PATTERNS[[a]], basename(all_rds), fixed = TRUE)]
    if (length(hit)) {
        arm_file[[a]] <- hit[1]
        claimed <- c(claimed, hit)
    }
}
main_hit <- all_rds[basename(all_rds) == MAIN_BASENAME]
.main_src <- "tagged"
if (!length(main_hit)) {
    # Fall back to the untagged name.
    main_hit <- all_rds[basename(all_rds) == MAIN_BASE_LEGACY]
    if (length(main_hit)) .main_src <- "LEGACY UNTAGGED"
}
if (length(main_hit)) {
    arm_file[[MAIN_ARM]] <- main_hit[1]
    claimed <- c(claimed, main_hit)
}
.legacy_present <- any(basename(all_rds) == MAIN_BASE_LEGACY)
.main_txt <- sprintf(" | main arm from the %s name '%s'%s", .main_src,
                     basename(main_hit[1]),
                     if (.main_src == "tagged" && .legacy_present)
                         paste0("; '", MAIN_BASE_LEGACY,
                                "' also present, not used")
                     else if (.main_src != "tagged")
                         paste0("; '", MAIN_BASENAME, "' not found")
                     else "")
loose <- setdiff(all_rds, claimed)
.loose_txt <- ""
if (length(loose))
    .loose_txt <- sprintf(paste(" | %d unclaimed file(s) ignored (%s)"),
                          length(loose),
                          paste(basename(loose), collapse = ", "))

gr11("GR-11-01",
     if (any(is.na(arm_file)) || nzchar(.loose_txt) ||
         .main_src != "tagged" || .legacy_present) "WARN" else "PASS",
     sprintf(paste("%d of %d arm RDS resolved under %s (%s)%s%s"),
             sum(!is.na(arm_file)), length(ARM_ORDER), READ_ROOT,
             paste(sprintf("%s: %s", names(arm_file),
                           ifelse(is.na(arm_file), "ABSENT",
                                  basename(arm_file))), collapse = "; "),
             .main_txt, .loose_txt))
if (all(is.na(arm_file)))
    stop("No Part 10h results RDS found under ", READ_ROOT, call. = FALSE)


# == STEP 2: donor frames =====================================================
cat("\nSTEP 2 -- donor frames\n")
NEED <- c(COL_MOD, COL_NM, COV4, OUTCOMES$col)
frames <- list()
skips  <- character(0)

for (a in ARM_ORDER) {
    if (is.na(arm_file[[a]])) next
    obj <- try(readRDS(arm_file[[a]]), silent = TRUE)
    if (inherits(obj, "try-error")) {
        skips <- c(skips, sprintf("%s: RDS unreadable", a)); next
    }
    dss <- names(obj$datasets)
    if (is.null(dss)) {
        skips <- c(skips, sprintf("%s: no datasets slot", a)); next
    }
    for (ds in dss) {
        pd <- obj$datasets[[ds]]$prop_df
        if (is.null(pd)) {
            skips <- c(skips, sprintf("%s/%s: no prop_df", a, ds)); next
        }
        gone <- NEED[!NEED %in% names(pd)]
        if (length(gone)) {
            skips <- c(skips, sprintf("%s/%s: missing %s", a, ds,
                                      paste(gone, collapse = ", ")))
            next
        }
        d <- pd[stats::complete.cases(pd[, NEED]), , drop = FALSE]
        if (nrow(d) < MIN_DONORS) {
            skips <- c(skips, sprintf("%s/%s: %d donors, below %d",
                                      a, ds, nrow(d), MIN_DONORS))
            next
        }
        frames[[paste(a, ds, sep = "|")]] <-
            list(arm = a, dataset = ds, d = d, full = pd)
        cat(sprintf("  %-13s %-8s %3d donors complete on all six outcomes\n",
                    a, ds, nrow(d)))
    }
}

gr11("GR-11-02", if (length(skips)) "WARN" else "PASS",
     sprintf(paste("%d arm/dataset frame(s) on the common donor set; %s"),
             length(frames),
             if (length(skips)) paste("skipped:",
                                      paste(skips, collapse = "; "))
             else "nothing skipped"))
if (!length(frames)) stop("No usable donor frames.", call. = FALSE)

set.seed(BOOT_SEED)


# == STEP 3, SECTION A: shared-component ratio ================================
cat("\nSTEP 3 -- section A, shared-component ratio\n")
A_rows <- list()
for (key in names(frames)) {
    fr <- frames[[key]]; d <- fr$d; n <- nrow(d)
    for (k in seq_len(nrow(OUTCOMES))) {
        obs <- .scr(d, OUTCOMES$col[k], OUTCOMES$orient[k])
        bs <- rep(NA_real_, N_BOOT); nfl <- 0L
        for (b in seq_len(N_BOOT)) {
            st <- .scr(d[sample.int(n, n, replace = TRUE), , drop = FALSE],
                       OUTCOMES$col[k], OUTCOMES$orient[k])
            if (is.finite(st$denom) && abs(st$denom) < DENOM_FLOOR)
                nfl <- nfl + 1L
            bs[b] <- st$ratio
        }
        ci <- .q(bs)
        A_rows[[length(A_rows) + 1L]] <- data.frame(
            arm = fr$arm, dataset = fr$dataset, outcome = OUTCOMES$label[k],
            domain = OUTCOMES$domain[k], n_donors = n,
            r = obs$r, beta_module = obs$b_mod, beta_nonmodule = obs$b_nm,
            denom = obs$denom, scr_pct = 100 * obs$ratio,
            ci_lo_pct = 100 * ci[1], ci_hi_pct = 100 * ci[2],
            boot_below_floor = nfl,
            covers_100 = if (all(is.finite(ci)))
                (100 * ci[1] <= 100 && 100 * ci[2] >= 100) else NA,
            stringsAsFactors = FALSE)
    }
    cat(sprintf("  %-13s %-8s done\n", fr$arm, fr$dataset))
}
if (!length(A_rows)) stop("Section A produced no rows.", call. = FALSE)
A <- do.call(rbind, A_rows)

# A row is unstable when more than 10% of resamples fall below the floor.
UNSTABLE_FRAC <- 0.10
A$unstable <- A$boot_below_floor > UNSTABLE_FRAC * N_BOOT

.unstable <- sum(A$unstable, na.rm = TRUE)
.cov_real <- A[A$arm != CONTROL_ARM & A$domain == "cognition" &
               !is.na(A$covers_100), , drop = FALSE]
gr11("GR-11-03", if (.unstable > 0L) "WARN" else "PASS",
     sprintf(paste("%d of %d ratio rows unstable (denominator below %.2f in",
                   "more than %.0f%% of resamples); boot_below_floor non-zero",
                   "on %d rows; real-arm cognition rows covering 100%%: %d of",
                   "%d (%d without a finite interval)"),
             .unstable, nrow(A), DENOM_FLOOR, 100 * UNSTABLE_FRAC,
             sum(A$boot_below_floor > 0L, na.rm = TRUE),
             sum(.cov_real$covers_100), nrow(.cov_real),
             sum(A$arm != CONTROL_ARM & A$domain == "cognition" &
                 is.na(A$covers_100))))


# == STEP 4, SECTION B: the contrast ==========================================
cat("\nSTEP 4 -- section B, the contrast\n")
B_rows <- list()
for (key in names(frames)) {
    fr <- frames[[key]]; d <- fr$d; n <- nrow(d)
    for (k in seq_len(nrow(OUTCOMES))) {
        f0 <- .joint(d, OUTCOMES$col[k], OUTCOMES$orient[k])
        if (is.null(f0)) next
        bs <- rep(NA_real_, N_BOOT)
        for (b in seq_len(N_BOOT)) {
            f2 <- .joint(d[sample.int(n, n, replace = TRUE), , drop = FALSE],
                         OUTCOMES$col[k], OUTCOMES$orient[k])
            if (!is.null(f2)) bs[b] <- f2$delta
        }
        bci   <- .q(bs)
        tcrit <- stats::qt(CI_HI, f0$df)
        # Minimum effect this row could have detected at POWER_TARGET.
        mdes  <- (stats::qt(CI_HI, f0$df) +
                  stats::qt(POWER_TARGET, f0$df)) * f0$se
        B_rows[[length(B_rows) + 1L]] <- data.frame(
            arm = fr$arm, dataset = fr$dataset, outcome = OUTCOMES$label[k],
            domain = OUTCOMES$domain[k], n_donors = n,
            b_module = f0$b_mod, b_nonmodule = f0$b_nm,
            delta = f0$delta, se = f0$se,
            exact_lo = f0$delta - tcrit * f0$se,
            exact_hi = f0$delta + tcrit * f0$se,
            p_exact = f0$p, boot_lo = bci[1], boot_hi = bci[2],
            mdes80 = mdes,
            stringsAsFactors = FALSE)
    }
    cat(sprintf("  %-13s %-8s done\n", fr$arm, fr$dataset))
}
if (!length(B_rows)) stop("Section B produced no rows.", call. = FALSE)
B <- do.call(rbind, B_rows)

# BH within arm x dataset over the six outcomes.
B$p_bh <- NA_real_
for (a in unique(B$arm)) {
    for (ds in unique(B$dataset[B$arm == a])) {
        s <- B$arm == a & B$dataset == ds & is.finite(B$p_exact)
        if (any(s)) B$p_bh[s] <- stats::p.adjust(B$p_exact[s], "BH")
    }
}
B$verdict <- ifelse(!is.finite(B$p_bh), "not computed",
             ifelse(B$p_bh < 0.05, "module > non-module",
             ifelse(abs(B$delta) < B$mdes80, "no difference detectable",
                    "no difference")))

.sig <- B[is.finite(B$p_bh) & B$p_bh < 0.05, , drop = FALSE]
gr11("GR-11-04", if (nrow(.sig) == 0L) "WARN" else "PASS",
     sprintf(paste("%d of %d contrast rows significant after BH within",
                   "arm x dataset (%s)"),
             nrow(.sig), nrow(B),
             if (nrow(.sig))
                 paste(sprintf("%s/%s %s", .sig$arm, .sig$dataset,
                               .sig$outcome), collapse = "; ")
             else "none"))

.uninf <- sum(B$verdict == "no difference detectable", na.rm = TRUE)
gr11("GR-11-08", if (.uninf > 0L) "WARN" else "PASS",
     sprintf(paste("%d of %d rows non-significant with |delta| below the",
                   "minimum effect detectable at %.0f%% power"),
             .uninf, nrow(B), 100 * POWER_TARGET))

# The control must be null on the contrast.
.ctl <- B[B$arm == CONTROL_ARM & is.finite(B$p_bh), , drop = FALSE]
.ctl_sig <- .ctl[.ctl$p_bh < 0.05, , drop = FALSE]
gr11("GR-11-07",
     if (!nrow(.ctl)) "WARN" else if (nrow(.ctl_sig)) "WARN" else "PASS",
     sprintf(paste("control arm %s: %d of %d contrast rows significant after",
                   "BH%s"),
             CONTROL_ARM, nrow(.ctl_sig), nrow(.ctl),
             if (nrow(.ctl_sig))
                 sprintf(" (%s)",
                         paste(sprintf("%s %s", .ctl_sig$dataset,
                                       .ctl_sig$outcome), collapse = "; "))
             else ""))


# == STEP 5, SECTION C: the balance ===========================================
cat("\nSTEP 5 -- section C, cognition-to-pathology balance\n")
C_rows <- list()
for (key in names(frames)) {
    fr <- frames[[key]]; d <- fr$d; n <- nrow(d)
    obs <- .balance(d)
    dif <- NA_real_
    if (is.finite(obs$module$ratio) && is.finite(obs$nonmodule$ratio))
        dif <- obs$module$ratio - obs$nonmodule$ratio
    bm <- rep(NA_real_, N_BOOT); bn <- bm; bd <- bm; nfl <- 0L
    for (b in seq_len(N_BOOT)) {
        st <- .balance(d[sample.int(n, n, replace = TRUE), , drop = FALSE])
        if (is.finite(st$module$path) && st$module$path < PATH_FLOOR)
            nfl <- nfl + 1L
        bm[b] <- st$module$ratio
        bn[b] <- st$nonmodule$ratio
        if (is.finite(st$module$ratio) && is.finite(st$nonmodule$ratio))
            bd[b] <- st$module$ratio - st$nonmodule$ratio
    }
    cm <- .q(bm); cn <- .q(bn); cd <- .q(bd)
    okd <- bd[is.finite(bd)]
    pb <- NA_real_
    if (length(okd) >= 200L) pb <- 2 * min(mean(okd <= 0), mean(okd >= 0))
    C_rows[[length(C_rows) + 1L]] <- data.frame(
        arm = fr$arm, dataset = fr$dataset, n_donors = n,
        mod_cog = obs$module$cog, mod_path = obs$module$path,
        mod_ratio = obs$module$ratio, mod_lo = cm[1], mod_hi = cm[2],
        nm_cog = obs$nonmodule$cog, nm_path = obs$nonmodule$path,
        nm_ratio = obs$nonmodule$ratio, nm_lo = cn[1], nm_hi = cn[2],
        diff = dif, diff_lo = cd[1], diff_hi = cd[2], p_boot = pb,
        boot_path_below_floor = nfl,
        stringsAsFactors = FALSE)
    cat(sprintf("  %-13s %-8s module %.2f  non-module %.2f  diff %.2f\n",
                fr$arm, fr$dataset, obs$module$ratio, obs$nonmodule$ratio,
                dif))
}
if (!length(C_rows)) stop("Section C produced no rows.", call. = FALSE)
CBAL <- do.call(rbind, C_rows)
# Recorded, not applied: the family is one test per arm x dataset.
CBAL$p_bh_if_pooled <- stats::p.adjust(CBAL$p_boot, "BH")

.cs <- CBAL[is.finite(CBAL$p_boot) & CBAL$p_boot < 0.05, , drop = FALSE]
gr11("GR-11-05", if (nrow(.cs) == 0L) "WARN" else "PASS",
     sprintf(paste("balance difference significant on %d of %d arm/dataset",
                   "rows (%s); pooled BH over %d rows: %d significant"),
             nrow(.cs), nrow(CBAL),
             if (nrow(.cs)) paste(sprintf("%s/%s", .cs$arm, .cs$dataset),
                                  collapse = "; ") else "none",
             nrow(CBAL),
             sum(CBAL$p_bh_if_pooled < 0.05, na.rm = TRUE)))


# == STEP 6, SECTION D: donor-set sensitivity =================================
cat("\nSTEP 6 -- section D, per-outcome donor sets\n")
D_rows <- list()
for (key in names(frames)) {
    fr <- frames[[key]]; pd <- fr$full
    for (k in seq_len(nrow(OUTCOMES))) {
        cols <- c(OUTCOMES$col[k], COL_MOD, COL_NM, COV4)
        dk <- pd[stats::complete.cases(pd[, cols]), , drop = FALSE]
        if (nrow(dk) < MIN_DONORS) next
        f0 <- .joint(dk, OUTCOMES$col[k], OUTCOMES$orient[k])
        if (is.null(f0)) next
        tcrit <- stats::qt(CI_HI, f0$df)
        D_rows[[length(D_rows) + 1L]] <- data.frame(
            arm = fr$arm, dataset = fr$dataset, outcome = OUTCOMES$label[k],
            domain = OUTCOMES$domain[k], n_donors = nrow(dk),
            b_module = f0$b_mod, b_nonmodule = f0$b_nm,
            delta = f0$delta, se = f0$se,
            exact_lo = f0$delta - tcrit * f0$se,
            exact_hi = f0$delta + tcrit * f0$se,
            p_exact = f0$p,
            stringsAsFactors = FALSE)
    }
}
D <- NULL
if (length(D_rows)) D <- do.call(rbind, D_rows)

if (is.null(D)) {
    gr11("GR-11-06", "WARN",
         paste("section D produced no rows"))
} else {
    D$p_bh <- NA_real_
    for (a in unique(D$arm)) {
        for (ds in unique(D$dataset[D$arm == a])) {
            s_ <- D$arm == a & D$dataset == ds & is.finite(D$p_exact)
            if (any(s_)) D$p_bh[s_] <- stats::p.adjust(D$p_exact[s_], "BH")
        }
    }
    .mg <- merge(
        B[, c("arm", "dataset", "outcome", "delta", "p_bh", "n_donors")],
        D[, c("arm", "dataset", "outcome", "delta", "p_bh", "n_donors")],
        by = c("arm", "dataset", "outcome"), suffixes = c("_common", "_own"))
    if (!nrow(.mg)) {
        gr11("GR-11-06", "WARN",
             paste("no row shared by the common and per-outcome donor sets"))
    } else {
        .flip <- sum(sign(.mg$delta_common) != sign(.mg$delta_own),
                     na.rm = TRUE)
        .dis <- sum((.mg$p_bh_common < 0.05) != (.mg$p_bh_own < 0.05),
                    na.rm = TRUE)
        .gain <- .mg$n_donors_own - .mg$n_donors_common
        gr11("GR-11-06",
             if (.flip > 0L || .dis > 2L) "WARN" else "PASS",
             sprintf(paste("common against per-outcome donor sets: %d sign",
                           "flip(s), %d verdict disagreement(s) of %d rows;",
                           "donors gained %d to %d"),
                     .flip, .dis, nrow(.mg),
                     min(.gain, na.rm = TRUE), max(.gain, na.rm = TRUE)))
    }
}


# == STEP 7: figures ==========================================================
cat("\nSTEP 7 -- figures\n")

.arm_present <- ARM_ORDER[ARM_ORDER %in% unique(B$arm)]
.arm_f <- function(x) factor(unname(ARM_LABEL[x]),
                             levels = unname(ARM_LABEL[.arm_present]))
.ds_f  <- function(x) factor(x, levels = DS_ORDER[DS_ORDER %in% unique(x)])
.out_f <- function(x) factor(x, levels = rev(OUTCOMES$label))
.dom_f <- function(x) factor(x, levels = c("cognition", "pathology"))
DOM_COL <- c(cognition = COG_COL, pathology = PATH_COL)

.base_theme <- function() {
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 13, hjust = .5),
          plot.subtitle = element_text(size = 8, colour = "grey30",
                                       hjust = .5),
          strip.text = element_text(face = "bold", colour = "white"),
          strip.background = element_rect(fill = "#1F3864"),
          panel.grid.minor = element_blank(),
          legend.position = "bottom",
          plot.background = element_rect(fill = "white", colour = NA))
}

# The device opens inside a function so that on.exit closes it.
.render <- function(p, base, w, h) {
    one <- function(dev_open, path) {
        tryCatch({
            dev_open(path)
            on.exit(grDevices::dev.off(), add = TRUE)
            print(p)
            TRUE
        }, error = function(e) {
            cat("    ", basename(path), " failed: ", conditionMessage(e),
                "\n", sep = ""); FALSE })
    }
    o1 <- one(function(f) grDevices::pdf(f, width = w, height = h),
              paste0(base, ".pdf"))
    o2 <- one(function(f) grDevices::png(f, width = w * 300, height = h * 300,
                                         res = 300),
              paste0(base, ".png"))
    ok <- o1 && o2
    if (ok) cat(sprintf("  %s\n", basename(base)))
    ok
}

.n_fig <- 0L
.fig_h <- 2.0 + 1.55 * length(.arm_present)

# == F1: the contrast forest ==================================================
BF <- B
BF$arm_f <- .arm_f(BF$arm); BF$dataset_f <- .ds_f(BF$dataset)
BF$outcome_f <- .out_f(BF$outcome); BF$domain_f <- .dom_f(BF$domain)
if (anyNA(BF$arm_f) || anyNA(BF$outcome_f) || anyNA(BF$domain_f))
    stop("F1: an unmapped factor level would be dropped silently.",
         call. = FALSE)
BF$star <- ifelse(is.finite(BF$p_bh) & BF$p_bh < 0.05, "*", "")
# Filled point and * mark BH p < 0.05; hollow is not significant.
BF$sig_f <- factor(is.finite(BF$p_bh) & BF$p_bh < 0.05, levels = c(TRUE, FALSE))

# Interval and point in one layer, so the shape stays on its own row.
p1 <- ggplot(BF, aes(x = delta, y = outcome_f, colour = domain_f,
                     group = domain_f)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_pointrange(aes(xmin = exact_lo, xmax = exact_hi, shape = sig_f),
                    linewidth = .7, size = .55, fatten = 4.2,
                    orientation = "y") +
    geom_text(aes(label = star), nudge_y = .28, size = 5, show.legend = FALSE) +
    facet_grid(arm_f ~ dataset_f) +
    scale_shape_manual(values = c(`TRUE` = 19, `FALSE` = 1),
                       name = "FDR < 0.05",
                       labels = c(`TRUE` = "yes", `FALSE` = "no"),
                       drop = FALSE, na.translate = FALSE) +
    scale_colour_manual(values = DOM_COL, name = NULL) +
    labs(title = "Module against non-module, same model, same donors",
         subtitle = paste("Standardised coefficient difference",
                          "b(module) - b(non-module) with its exact 95%",
                          "interval. Outcomes oriented so positive means more",
                          "disease.\nFilled point and * mark BH p < 0.05",
                          "within",
                          "arm and dataset; hollow is not significant. A",
                          "non-significant row is not evidence of equality --",
                          "see the power figure."),
         x = "b(module) - b(non-module), standardised", y = NULL) +
    .base_theme()
if (.render(p1, file.path(OUT_DIR, paste0(FILE_PREFIX, "_F1_Contrast")),
            9, .fig_h + 1)) .n_fig <- .n_fig + 1L

# == F2: paired coefficients ==================================================
LB <- rbind(
    data.frame(BF[, c("arm_f", "dataset_f", "outcome_f", "domain_f")],
               population = "module", beta = BF$b_module,
               stringsAsFactors = FALSE),
    data.frame(BF[, c("arm_f", "dataset_f", "outcome_f", "domain_f")],
               population = "non-module", beta = BF$b_nonmodule,
               stringsAsFactors = FALSE))
LB$population <- factor(LB$population, levels = c("module", "non-module"))

p2 <- ggplot() +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_segment(data = BF,
                 aes(x = b_nonmodule, xend = b_module,
                     y = outcome_f, yend = outcome_f),
                 colour = "grey65", linewidth = .6) +
    geom_point(data = LB, aes(x = beta, y = outcome_f, fill = population),
               shape = 21, size = 2.6, colour = "white", stroke = .5) +
    facet_grid(arm_f ~ dataset_f) +
    scale_fill_manual(values = c(module = MODULE_COL,
                                 `non-module` = NONMOD_COL), name = NULL) +
    labs(title = "Where each population's association sits",
         subtitle = paste("Both standardised coefficients from the SAME joint",
                          "model. The gap between the two points is the",
                          "contrast in the previous figure."),
         x = "standardised coefficient (positive = more disease)", y = NULL) +
    .base_theme()
if (.render(p2, file.path(OUT_DIR, paste0(FILE_PREFIX, "_F2_Paired")),
            9, .fig_h + 1)) .n_fig <- .n_fig + 1L

# == F3: the balance ==========================================================
CB <- CBAL
CB$arm_f <- .arm_f(CB$arm); CB$dataset_f <- .ds_f(CB$dataset)
LC <- rbind(
    data.frame(CB[, c("arm_f", "dataset_f")], population = "module",
               ratio = CB$mod_ratio, lo = CB$mod_lo, hi = CB$mod_hi,
               stringsAsFactors = FALSE),
    data.frame(CB[, c("arm_f", "dataset_f")], population = "non-module",
               ratio = CB$nm_ratio, lo = CB$nm_lo, hi = CB$nm_hi,
               stringsAsFactors = FALSE))
LC$population <- factor(LC$population, levels = c("module", "non-module"))
CB$lab <- sprintf("difference %.2f [%.2f, %.2f]\nbootstrap p %.3f",
                  CB$diff, CB$diff_lo, CB$diff_hi, CB$p_boot)

p3 <- ggplot(LC, aes(x = population, y = ratio, colour = population)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_pointrange(aes(ymin = lo, ymax = hi), linewidth = .7, size = .5) +
    geom_text(data = CB, aes(x = 1.5, y = Inf, label = lab),
              inherit.aes = FALSE, vjust = 1.25, size = 2.6,
              colour = "grey20") +
    facet_grid(arm_f ~ dataset_f) +
    scale_colour_manual(values = c(module = MODULE_COL,
                                   `non-module` = NONMOD_COL), name = NULL) +
    labs(title = "Cognition effect relative to pathology effect",
         subtitle = paste("Mean oriented cognition beta divided by mean",
                          "pathology beta, with a donor bootstrap interval.",
                          "Above 1 means the population is cognition-weighted.",
                          "\nThe dissociation is the DIFFERENCE between the",
                          "two points, printed in each panel."),
         x = NULL, y = "cognition / pathology") +
    .base_theme()
if (.render(p3, file.path(OUT_DIR, paste0(FILE_PREFIX, "_F3_Balance")),
            8, .fig_h + 1.2)) .n_fig <- .n_fig + 1L

# == F4: the shared-component ratio ===========================================
AF <- A
AF$arm_f <- .arm_f(AF$arm); AF$dataset_f <- .ds_f(AF$dataset)
AF$outcome_f <- .out_f(AF$outcome); AF$domain_f <- .dom_f(AF$domain)
AF$lo_c <- pmax(AF$ci_lo_pct, SCR_CLIP[1])
AF$hi_c <- pmin(AF$ci_hi_pct, SCR_CLIP[2])

p4 <- ggplot(AF, aes(x = scr_pct, y = outcome_f, colour = domain_f)) +
    geom_vline(xintercept = 100, linetype = "dashed", colour = "grey40") +
    geom_linerange(aes(xmin = lo_c, xmax = hi_c,
                       alpha = ifelse(unstable, "unstable", "usable")),
                   linewidth = .7) +
    geom_point(size = 2.2) +
    facet_grid(arm_f ~ dataset_f) +
    scale_colour_manual(values = DOM_COL, name = NULL) +
    scale_alpha_manual(values = c(usable = 1, unstable = .3), name = NULL) +
    coord_cartesian(xlim = SCR_CLIP) +
    labs(title = "The shared-component ratio, reported with its interval",
         subtitle = paste("Non-module coefficient as a percentage of what its",
                          "correlation with the module predicts. 100% means",
                          "the association is entirely inherited.\nIntervals",
                          "are clipped to the panel; faded rows had an",
                          "unstable denominator. This is why the ratio is",
                          "reported and not led on."),
         x = "non-module beta as % of r x module beta", y = NULL) +
    .base_theme()
if (.render(p4, file.path(OUT_DIR, paste0(FILE_PREFIX, "_F4_SharedRatio")),
            9, .fig_h + 1)) .n_fig <- .n_fig + 1L

# == F5: what each row could have detected ====================================
PW <- BF
PW$abs_delta <- abs(PW$delta)
p5 <- ggplot(PW, aes(y = outcome_f)) +
    geom_col(aes(x = mdes80), fill = "grey85", width = .62) +
    geom_point(aes(x = abs_delta, colour = domain_f), size = 2.4) +
    facet_grid(arm_f ~ dataset_f) +
    scale_colour_manual(values = DOM_COL, name = NULL) +
    labs(title = "What each row could have detected",
         subtitle = paste("Grey bar is the smallest difference detectable at",
                          sprintf("%.0f%% power", 100 * POWER_TARGET),
                          "given that row's standard error. A point INSIDE",
                          "the bar\nmeans a non-significant result there is",
                          "uninformative, not a negative finding."),
         x = "absolute contrast, standardised", y = NULL) +
    .base_theme()
if (.render(p5, file.path(OUT_DIR, paste0(FILE_PREFIX, "_F5_Power")),
            9, .fig_h + 1)) .n_fig <- .n_fig + 1L


# == STEP 8: output ===========================================================
cat("\nSTEP 8 -- output\n")
wb <- list(Contrast = B,
           Balance = CBAL,
           Shared_Ratio = A,
           Power = B[, c("arm", "dataset", "outcome", "domain", "n_donors",
                         "delta", "se", "mdes80", "p_bh", "verdict")])
if (!is.null(D)) wb[["Donor_Sensitivity"]] <- D

xlsx_ok <- FALSE
if (requireNamespace("writexl", quietly = TRUE)) {
    xp <- file.path(OUT_DIR, paste0(FILE_PREFIX, "_CogVsPath.xlsx"))
    tryCatch({ writexl::write_xlsx(wb, xp); xlsx_ok <- TRUE
               cat(sprintf("  %s\n", basename(xp))) },
             error = function(e)
                 cat("  workbook failed:", conditionMessage(e), "\n"))
} else {
    cat("  writexl not installed; the workbook was not written\n")
}

saveRDS(list(contrast = B, balance = CBAL, shared_ratio = A,
             donor_sensitivity = D,
             meta = list(created = format(Sys.time()),
                         schema = "Part11 v1",
                         read_root = READ_ROOT, out_dir = OUT_DIR,
                         arms = names(arm_file)[!is.na(arm_file)],
                         n_boot = N_BOOT, boot_seed = BOOT_SEED,
                         denom_floor = DENOM_FLOOR, path_floor = PATH_FLOOR,
                         power_target = POWER_TARGET)),
        file.path(OUT_DIR, paste0(FILE_PREFIX, "_CogVsPath_results.rds")))

gr11("GR-11-09", if (.n_fig < 5L || !xlsx_ok) "WARN" else "PASS",
     sprintf(paste("%d of 5 figures written as PDF and PNG; workbook %s; RDS",
                   "written to %s"),
             .n_fig, if (xlsx_ok) "written" else "NOT written", OUT_DIR))


# == LEDGER ===================================================================
.lines <- c(strrep("=", 72),
            "Part 11 - Cognition vs Pathology - Guardrail Log [CP10K]",
            paste("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
            paste("Read:", READ_ROOT),
            paste("Write:", OUT_DIR),
            strrep("=", 72), "",
            vapply(.gr, function(g)
                sprintf("  [%-4s]  %-12s %s", g$status, g$id, g$detail),
                character(1)),
            "", strrep("=", 72),
            sprintf("  Total : %d", length(.gr)),
            sprintf("  PASS  : %d",
                    sum(vapply(.gr, function(g) g$status == "PASS",
                               logical(1)))),
            sprintf("  WARN  : %d",
                    sum(vapply(.gr, function(g) g$status == "WARN",
                               logical(1)))),
            sprintf("  FAIL  : %d",
                    sum(vapply(.gr, function(g) g$status == "FAIL",
                               logical(1)))),
            strrep("=", 72))
.log_path <- file.path(OUT_DIR,
                       paste0(FILE_PREFIX, "_CogVsPath_Guardrail_Log.txt"))
writeLines(.lines, .log_path)
cat("\n", paste(.lines[(length(.lines) - 5):length(.lines)], collapse = "\n"),
    "\n", sep = "")
cat("\nOutput directory:", OUT_DIR, "\n")
