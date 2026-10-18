# =============================================================================
# Part13b_CP10K_M4RandomSubset_code.r
# SEA-AD CEP-IP/CMRF, Part 13b, M4 cells against random PVM_3 cells
#   Compares each donor's M4 cells with random PVM_3 arm cells of the same
#   number, drawn per donor and region (B_NULL draws, one L'Ecuyer-CMRG stream
#   per draw), so the M4 share is identical in every draw. Every cell carries
#   its PVM_3 arm labels; the M4 arm's own labels are a reference only.
#
# Usage (from the SEA-AD folder, or set SEAAD_ROOT), after Part 10h, Part 13a
#   and Part 12 step 04
#   source("Part13b_CP10K_M4RandomSubset_code.r")
#
# Outputs (Part13b_M4RandomSubset/): pre-set readings, supplementary table,
#   null draws (CSV), supplementary figure (PDF, PNG), legend draft,
#   guardrail log and detail, console log.
# =============================================================================
# Settings take their default unless a list P13B_CFG exists before sourcing,
# e.g. P13B_CFG <- list(N_CORES = 4L). Other session variables are ignored.
cfg <- function(name, default) {
    o <- if (exists("P13B_CFG", envir = globalenv(), inherits = FALSE)) get("P13B_CFG", envir = globalenv()) else list()
    if (name %in% names(o)) o[[name]] else default
}

# == 0. CONFIGURATION =========================================================
CEP_HOME    <- cfg("CEP_HOME", Sys.getenv("SEAAD_ROOT", getwd()))
OUT_DIR     <- cfg("OUT_DIR", file.path(CEP_HOME, "Part13b_M4RandomSubset"))
if (!dir.exists(CEP_HOME)) stop("CEP_HOME does not exist: ", CEP_HOME, call. = FALSE)
if (!startsWith(normalizePath(dirname(OUT_DIR), winslash = "/", mustWork = FALSE),
                normalizePath(CEP_HOME, winslash = "/", mustWork = FALSE)))
    stop("OUT_DIR must sit inside CEP_HOME: ", OUT_DIR, call. = FALSE)
FILE_PREFIX <- "Part13b_CP10K_M4RandomSubset"
SEARCH_EXCLUDE <- cfg("SEARCH_EXCLUDE", c("Previous_version", "Archive", "archive"))

COMP_RDS         <- cfg("COMP_RDS", NA)      # NA = search CEP_HOME for the basename
COMP_BASENAME    <- "Part10h_CP10K_M4_Comparator_results.rds"
P10H_M4_RDS      <- cfg("P10H_M4_RDS", NA)
P10H_M4_BASENAME <- "Part10h_CP10K_M4_Prognosis_results.rds"
P10H_P3_RDS      <- cfg("P10H_P3_RDS", NA)
P10H_P3_BASENAME <- "Part10h_CP10K_PVM_3_Prognosis_results.rds"

# Cell label file of the Part 12 rerun (its mito_pct gives the Figure 14D %MT).
LAB_FILE       <- cfg("LAB_FILE", file.path(CEP_HOME, "Part12_rerun", "04_Labels", "04_cell_labels.csv.gz"))
LAB_BASENAME   <- "04_cell_labels.csv.gz"
# NA = found automatically; a partial list overrides, e.g. list(barcode = "cell_id").
LAB_COLS   <- utils::modifyList(list(donor = NA, dataset = NA, arm = NA, barcode = NA, mito = NA, region = NA),
                                cfg("LAB_COLS", list()))
ARM_VALUES <- utils::modifyList(list(M4 = c("M4", "Mic13_PVM3", "Mic13"),
                                     PVM3 = c("PVM3_all", "PVM_3", "PVM3", "PVM3_bulk")),
                                cfg("ARM_VALUES", list()))
# NULL = the label rule is found by reproduction; or set by hand, e.g.
#   list(single = list(col = "group", values = "Module post-IP"))
LAB_RULE <- cfg("LAB_RULE", NULL)

# Cell-level label check against the Validation 1 workbooks (Step4_Refined_Cells).
V1_CHECK  <- cfg("V1_CHECK", TRUE)
V1_FILES  <- cfg("V1_FILES", list())   # e.g. list(M4_NeoC = "<path to the Validation 1 workbook>")
V1_TAGS   <- c(M4 = "M4", PVM3 = "PVM_3")
V1_MODULE_CLASSES <- c("Original + Refined", "Recovered")   # Part 10h MODULE_CLASSES

MT_FILE       <- cfg("MT_FILE", file.path(CEP_HOME, "Part13_M4Share", "Part13_M4_pctMT_for14D.csv"))
MT_COLS       <- list(donor = "donor_id", dataset = "dataset", value = "pct_mt")

ARM_TAG_EXPECT    <- "Mic13_PVM3"
DATASETS          <- c("NeoC", "NonNeoC")
N_DONORS_EXPECT   <- cfg("N_DONORS_EXPECT", c(NeoC = 84L, NonNeoC = 83L))
N_M4_CELLS_EXPECT <- cfg("N_M4_CELLS_EXPECT", c(NeoC = 2130L, NonNeoC = 1987L))   # Part 9 GOI-positive, M4
N_P3_CELLS_EXPECT <- cfg("N_P3_CELLS_EXPECT", c(NeoC = 13331L, NonNeoC = 11132L)) # Part 9 GOI-positive, PVM_3
MIN_M4_CELLS_FRAC <- 5L
STRATA            <- cfg("STRATA", "donor_region")   # "donor_region" or "donor"
B_NULL            <- cfg("B_NULL", 10000L)
NULL_SEED         <- cfg("NULL_SEED", 20261011L)
NULL_MIN_VALID    <- 0.99   # share of draws in which every model was fitted on the observed donors
# Workers: all logical cores but one.
N_CORES           <- cfg("N_CORES", max(1L, parallel::detectCores(logical = TRUE) - 1L, na.rm = TRUE))
# Worker reproduction tolerance: a multithreaded BLAS can change the last digits.
TOL_REPRO         <- 1e-9

MIN_DONORS <- 20L; VIF_MAX <- 5; ALPHA <- 0.05
GUARD_STRICT         <- cfg("GUARD_STRICT", TRUE)
ALLOW_PRESET_REWRITE <- cfg("ALLOW_PRESET_REWRITE", FALSE)
TOL_EXACT <- 1e-8
TOL_PUB   <- 1.5e-3
# Published Figure 14D rows (M4, Module post-IP% + %MT): beta, CI low, CI high, n, FDR.
PUB_14D <- cfg("PUB_14D", list(c("NeoC", "MMSE score", -0.529, -0.808, -0.251, 74, 0.003),
                                c("NonNeoC", "Cognitive status", 0.415, 0.161, 0.670, 80, 0.015),
                                c("NonNeoC", "MMSE score", -0.337, -0.575, -0.099, 76, 0.028)))

# == 1. LIBRARIES, CONSOLE LOG, LEDGER ========================================
suppressPackageStartupMessages({
    library(openxlsx); library(ggplot2); library(patchwork)
})
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
STAMP    <- format(Sys.time(), "%Y%m%d_%H%M%S")
CONSOLE  <- file.path(OUT_DIR, sprintf("%s_console_log_%s.txt", FILE_PREFIX, STAMP))
LEDGER   <- file.path(OUT_DIR, paste0(FILE_PREFIX, "_Guardrail_Log.txt"))
GDETAIL  <- file.path(OUT_DIR, paste0(FILE_PREFIX, "_Guardrail_Detail.txt"))
PRESET   <- file.path(OUT_DIR, paste0(FILE_PREFIX, "_preset_readings.txt"))
OUT_XLSX <- file.path(OUT_DIR, paste0(FILE_PREFIX, "_SuppTable.xlsx"))
OUT_CSV  <- file.path(OUT_DIR, paste0(FILE_PREFIX, "_null_draws.csv"))
OUT_PDF  <- file.path(OUT_DIR, paste0(FILE_PREFIX, "_SuppFig.pdf"))
OUT_PNG  <- file.path(OUT_DIR, paste0(FILE_PREFIX, "_SuppFig.png"))
OUT_LEG  <- file.path(OUT_DIR, paste0(FILE_PREFIX, "_legend_draft.txt"))

.gr <- list(); .gd <- list()
gr13b <- function(id, status, detail) {
    if (length(detail) != 1L || is.na(detail)) detail <- "[detail unavailable]"
    .gr[[length(.gr) + 1L]] <<- data.frame(id = id, status = status, detail = detail,
                                           stringsAsFactors = FALSE)
    cat(sprintf("  [%-4s]  %-10s %s\n", status, id, detail))
    if (status == "FAIL" && GUARD_STRICT)
        stop(sprintf("%s FAILED: %s", id, detail), call. = FALSE)
}
gcmp <- function(id, item, mine, ref, tol) {
    mine <- as.numeric(mine); ref <- as.numeric(ref)
    if (!length(mine) || !length(ref)) { mine <- NA_real_; ref <- NA_real_ }
    ok <- !is.na(mine) & !is.na(ref) & abs(mine - ref) <= tol
    .gd[[length(.gd) + 1L]] <<- data.frame(id = id, item = item, mine = mine, ref = ref,
                                           abs_diff = abs(mine - ref), tol = tol, pass = ok,
                                           stringsAsFactors = FALSE)
    ok
}
write_ledger <- function(final_note) {
    G <- if (length(.gr)) do.call(rbind, .gr) else
        data.frame(id = character(0), status = character(0), detail = character(0))
    G <- G[order(G$id), , drop = FALSE]
    lines <- c(strrep("=", 72),
               "Part 13b - M4 cells against random PVM_3 cells - Guardrail Log [CP10K]",
               paste("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
               paste("Root:", CEP_HOME), paste("Write:", OUT_DIR),
               strrep("=", 72), "",
               if (nrow(G)) sprintf("  [%-4s]  %-10s %s", G$status, G$id, G$detail) else "  (none)",
               "", strrep("=", 72),
               sprintf("  Total : %d", nrow(G)),
               sprintf("  PASS  : %d", sum(G$status == "PASS")),
               sprintf("  WARN  : %d", sum(G$status == "WARN")),
               sprintf("  FAIL  : %d", sum(G$status == "FAIL")),
               sprintf("  SKIP  : %d", sum(G$status == "SKIP")),
               sprintf("  Run   : %s", final_note),
               strrep("=", 72))
    writeLines(lines, LEDGER)
    if (length(.gd)) utils::write.table(do.call(rbind, .gd), GDETAIL, sep = "\t", quote = FALSE,
                                        row.names = FALSE)
    invisible(G)
}

.stale <- c(OUT_XLSX, OUT_CSV, OUT_PDF, OUT_PNG, OUT_LEG, GDETAIL, LEDGER)
invisible(file.remove(.stale[file.exists(.stale)]))
.con <- file(CONSOLE, open = "wt")
sink(.con, split = TRUE)
.old_width <- options(width = 200)
cat(strrep("=", 78), "\n  Part 13b - M4 cells against random PVM_3 cells of the same number\n  ",
    format(Sys.time()), "\n  ", R.version.string, "\n  Root : ", CEP_HOME, "\n  Write: ", OUT_DIR, "\n",
    strrep("=", 78), "\n", sep = "")

# == 2. HELPERS ===============================================================
`%||%` <- function(a, b) if (is.null(a)) b else a
num <- function(x) suppressWarnings(as.numeric(x))
fmt <- function(x, d = 3) ifelse(is.na(x), "NA", formatC(x, format = "f", digits = d))
fmt_p <- function(x) ifelse(is.na(x), "NA", ifelse(x < 0.001, formatC(x, format = "e", digits = 2),
                                                    formatC(x, format = "f", digits = 3)))
bq <- function(v) sprintf("`%s`", v)
.z <- function(x) {                                   # comparator .z
    x <- num(x); s <- stats::sd(x, na.rm = TRUE)
    if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
    (x - mean(x, na.rm = TRUE)) / s
}
.safe_vif <- function(fit) {                          # comparator .safe_vif
    X <- stats::model.matrix(fit)[, -1, drop = FALSE]
    if (ncol(X) < 2L) return(stats::setNames(rep(1, ncol(X)), colnames(X)))
    v <- vapply(seq_len(ncol(X)), function(j) {
        r2 <- suppressWarnings(summary(stats::lm(X[, j] ~ X[, -j, drop = FALSE]))$r.squared)
        if (!is.finite(r2) || r2 >= 1) NA_real_ else 1 / (1 - r2)
    }, numeric(1))
    stats::setNames(v, colnames(X))
}
find_file <- function(path, base, required = TRUE) {
    if (!is.na(path)) {
        if (!file.exists(path)) stop("File not found: ", path, call. = FALSE)
        return(normalizePath(path, winslash = "/"))
    }
    hit <- list.files(CEP_HOME, pattern = paste0("^", gsub(".", "\\.", base, fixed = TRUE), "$"),
                      recursive = TRUE, full.names = TRUE)
    drop <- Reduce(`|`, lapply(SEARCH_EXCLUDE, function(s) grepl(s, hit, fixed = TRUE)),
                   rep(FALSE, length(hit)))
    if (any(drop)) cat(sprintf("  %s: %d copy(ies) under excluded folders ignored\n", base, sum(drop)))
    hit <- hit[!drop]
    if (length(hit) == 1L) return(normalizePath(hit, winslash = "/"))
    if (!length(hit) && !required) return(NA_character_)
    stop(sprintf("%s: %d match(es) under %s%s. Set the path in Section 0.", base, length(hit),
                 CEP_HOME, if (length(hit)) paste0(": ", paste(hit, collapse = "; ")) else ""),
         call. = FALSE)
}
sum_by <- function(v, g, n) {                         # sum of v per integer group 1..n
    out <- numeric(n)
    if (!length(v)) return(out)
    r <- rowsum(v, g, reorder = TRUE)
    out[as.integer(rownames(r))] <- r[, 1]
    out
}

# Reference fit: lm() with the comparator conventions; checks the fast fitter.
ref_fit <- function(D0, y, preds, covs) {
    covs <- setdiff(covs, y)
    vars <- unique(c(y, preds, covs))
    D <- as.data.frame(lapply(D0[, vars, drop = FALSE], num)); names(D) <- vars
    keep <- stats::complete.cases(D) & Reduce(`&`, lapply(D, is.finite))
    D <- D[keep, , drop = FALSE]; n <- nrow(D)
    if (n < MIN_DONORS) return(NULL)
    D[[y]] <- .z(D[[y]]); for (p in preds) D[[p]] <- .z(D[[p]])
    if (anyNA(D[, c(y, preds)])) return(NULL)
    fit <- stats::lm(stats::as.formula(paste(bq(y), "~", paste(bq(c(preds, covs)), collapse = " + "))), data = D)
    cf <- summary(fit)$coefficients; rownames(cf) <- gsub("`", "", rownames(cf))
    ci <- suppressMessages(stats::confint(fit)); rownames(ci) <- gsub("`", "", rownames(ci))
    vf <- .safe_vif(fit); names(vf) <- gsub("`", "", names(vf))
    list(n = n, beta = cf[preds, 1], se = cf[preds, 2], p = cf[preds, 4],
         ci_lo = ci[preds, 1], ci_hi = ci[preds, 2], vif = unname(vf[preds]))
}

# Fast fit for the null loop: the same model by QR. yz standardised outcome,
# P raw predictors (standardised here), C unscaled covariates.
ff_fit <- function(yz, P, C) {
    n <- length(yz); q <- ncol(P)
    Pz <- P
    for (j in seq_len(q)) {
        s <- stats::sd(P[, j]); if (!is.finite(s) || s == 0) return(NULL)
        Pz[, j] <- (P[, j] - mean(P[, j])) / s
    }
    X <- cbind(1, Pz, C)
    qx <- qr(X)
    if (qx$rank < ncol(X)) return(NULL)
    b <- qr.coef(qx, yz)
    res <- qr.resid(qx, yz)
    df <- n - ncol(X)
    s2 <- sum(res^2) / df
    V <- matrix(0, ncol(X), ncol(X))
    V[qx$pivot, qx$pivot] <- chol2inv(qr.R(qx))
    se <- sqrt(diag(V) * s2)
    idx <- 2:(q + 1)
    tt <- b[idx] / se[idx]
    qt975 <- stats::qt(0.975, df)
    R <- stats::cor(X[, -1, drop = FALSE])
    vif <- tryCatch(diag(solve(R))[seq_len(q)], error = function(e) rep(NA_real_, q))
    vif[!is.finite(vif)] <- NA_real_
    list(n = n, beta = unname(b[idx]), se = unname(se[idx]), p = unname(2 * stats::pt(-abs(tt), df)),
         ci_lo = unname(b[idx] - qt975 * se[idx]), ci_hi = unname(b[idx] + qt975 * se[idx]), vif = unname(vif))
}

VERDICTS <- c("Module only", "%MT only", "Both", "Neither", "Unreadable")
verdict_code <- function(vmax, q1, q2) {
    if (!is.finite(vmax) || vmax >= VIF_MAX || is.na(q1) || is.na(q2)) return(5)
    if (q1 < ALPHA && q2 < ALPHA) return(3)
    if (q1 < ALPHA) return(1)
    if (q2 < ALPHA) return(2)
    4
}

# == 2b. CELL SETS, FITS AND READOUTS (top level, for the workers) ============
CLIN_LAB <- c("Cognitive status", "MMSE score", "CASI score", "CERAD score", "Braak stage",
              "Thal phase", "ADNC rank", "CAA score", "APOE e4 dose")
COG_LAB <- CLIN_LAB[1:3]
SHORT <- c("Cognitive status" = "CogS", "MMSE score" = "MMSE", "CASI score" = "CASI")
EXP_SIGN <- c("Cognitive status" = 1, "MMSE score" = -1, "CASI score" = -1)

# Builds every cell set, observed or random, from one region set's frame.
agg_set <- function(P, sel) {
    g <- P$p$dix[sel]
    cnt <- tabulate(g, P$nD); a <- tabulate(g[P$p$.post[sel]], P$nD); ms <- sum_by(P$p$.mt[sel], g, P$nD)
    list(x = 100 * a / P$d$n_total, frac = ifelse(cnt >= MIN_M4_CELLS_FRAC, 100 * a / pmax(cnt, 1), NA_real_),
         mt = ifelse(cnt > 0, log1p(ms / pmax(cnt, 1)), NA_real_), share = 100 * cnt / P$d$n_total, cnt = cnt)
}

fit_one <- function(Y, lab, preds, V) {
    yy <- Y[[lab]]
    Pm <- do.call(cbind, V[preds])
    ok <- yy$base & (rowSums(!is.finite(Pm)) == 0)
    n <- sum(ok); if (n < MIN_DONORS) return(NULL)
    yz <- .z(yy$y[ok]); if (anyNA(yz)) return(NULL)        # outcome without variance: not fitted, counted
    r <- ff_fit(yz, Pm[ok, , drop = FALSE], yy$C[ok, , drop = FALSE])
    if (!is.null(r)) r$rows <- which(ok)
    r
}
# One readout vector per cell set: every quantity compared with the null.
readout <- function(Y, V) {
    out <- c(); ns <- c(); nfail <- 0L
    for (lab in COG_LAB) {
        s <- SHORT[[lab]]
        a <- fit_one(Y, lab, "x", V)
        c1 <- fit_one(Y, lab, c("share", "frac"), V)
        d1 <- fit_one(Y, lab, c("frac", "mt"), V)
        nfail <- nfail + is.null(a) + is.null(c1) + is.null(d1)
        out[paste0("bA_", s)] <- if (is.null(a)) NA else a$beta[1]
        out[paste0("bC_", s)] <- if (is.null(c1)) NA else c1$beta[2]
        out[paste0("vC_", s)] <- if (is.null(c1)) NA else if (anyNA(c1$vif)) Inf else max(c1$vif)
        out[paste0("bD_", s)] <- if (is.null(d1)) NA else d1$beta[1]
        out[paste0("bDm_", s)] <- if (is.null(d1)) NA else d1$beta[2]
        out[paste0("vD_", s)] <- if (is.null(d1)) NA else if (anyNA(d1$vif)) Inf else max(d1$vif)
        ns[paste0(c("MA_", "MC_", "MD_"), s)] <- c(a$n %||% NA, c1$n %||% NA, d1$n %||% NA)
    }
    bb <- lapply(CLIN_LAB, function(lab) fit_one(Y, lab, c("x", "mt"), V))
    nfail <- nfail + sum(vapply(bb, is.null, logical(1)))
    px <- vapply(bb, function(r) if (is.null(r)) NA_real_ else r$p[1], 0)
    pm <- vapply(bb, function(r) if (is.null(r)) NA_real_ else r$p[2], 0)
    qx <- px; qm <- pm
    qx[is.finite(px)] <- stats::p.adjust(px[is.finite(px)], "BH")
    qm[is.finite(pm)] <- stats::p.adjust(pm[is.finite(pm)], "BH")
    for (i in seq_along(CLIN_LAB)) if (CLIN_LAB[i] %in% COG_LAB) {
        s <- SHORT[[CLIN_LAB[i]]]; r <- bb[[i]]
        vm <- if (is.null(r)) NA_real_ else if (anyNA(r$vif)) Inf else max(r$vif)
        out[paste0("bB_", s)] <- if (is.null(r)) NA else r$beta[1]
        out[paste0("qB_", s)] <- qx[i]; out[paste0("qBm_", s)] <- qm[i]; out[paste0("vB_", s)] <- vm
        out[paste0("cB_", s)] <- verdict_code(vm, qx[i], qm[i])
    }
    for (i in seq_along(CLIN_LAB)) ns[paste0("MB_", CLIN_LAB[i])] <- bb[[i]]$n %||% NA
    kx <- is.finite(V$x) & is.finite(V$mt); kf <- is.finite(V$frac) & is.finite(V$mt)
    out["r_x_mt"] <- stats::cor(V$x[kx], V$mt[kx])
    out["r_frac_mt"] <- stats::cor(V$frac[kf], V$mt[kf])
    out["mean_frac"] <- mean(V$frac, na.rm = TRUE)
    attr(out, "ns") <- ns; attr(out, "nfail") <- nfail
    out
}

# One random subset: sample within strata, build the cell set, read it out.
# Each draw has its own L'Ecuyer-CMRG stream, so core count does not matter.
draw_one <- function(P) {
    u <- stats::runif(nrow(P$p))
    o <- order(P$p$six, u)                       # within stratum, random order
    rk <- sequence(P$size_s)                     # rank within stratum, strata in order
    o[rk <= rep(P$k_s, P$size_s)]
}
null_one <- function(seed, P, ns_obs) {
    assign(".Random.seed", seed, envir = globalenv())
    sel <- draw_one(P)
    V <- agg_set(P, sel)
    dev <- max(abs(V$cnt - P$d$k)) + max(abs(tabulate(P$p$six[sel], P$S) - P$k_s))
    V$share <- P$obs$share                       # identical by construction; checked by `dev`
    r <- readout(P$Y, V)
    nm <- sum(!mapply(identical, attr(r, "ns"), ns_obs[names(attr(r, "ns"))]))
    nf <- attr(r, "nfail"); attributes(r) <- list(names = names(r))
    list(r = r, dev = dev, dup = anyDuplicated(sel) > 0, foreign = sum(sel < 1L | sel > nrow(P$p)), nm = nm, nf = nf)
}
null_chunk <- function(seeds, P, ns_obs) lapply(seeds, null_one, P = P, ns_obs = ns_obs)
WORKER_OBJECTS <- c("num", ".z", "ff_fit", "sum_by", "verdict_code", "agg_set", "fit_one", "readout", "draw_one",
                    "null_one", "null_chunk", "%||%", "MIN_DONORS", "VIF_MAX", "ALPHA", "MIN_M4_CELLS_FRAC",
                    "CLIN_LAB", "COG_LAB", "SHORT")

# == 3. PRE-SET READINGS ======================================================
PRESET_TEXT <- c(
    "Part13b_CP10K_M4RandomSubset pre-set readings",
    "Fixed before the first run. A changed text stops the run unless ALLOW_PRESET_REWRITE = TRUE.",
    "",
    "Question: are the M4 results a property of the M4 cells, or would any PVM_3 cells of the same number per donor give the same results?",
    "",
    "Definitions",
    "- Pool: the PVM_3 arm cells (Micro-PVM_3 cells with MT-CO3 > 0). The M4 share counts only MT-CO3 > 0 cells, so the pool holds the same condition. The M4 cells are part of the pool.",
    "- Labels: every cell, M4 or random, carries its PVM_3 arm labels (PVM_3 Module and PVM_3 IP_y), so the M4 cells and the random cells are defined the same way. The M4 arm's own labels are reported as a reference without a reading.",
    "- Random subset: in each donor and region (each donor when the label file has no region column), as many pool cells as the donor has M4 cells there, drawn without replacement. The M4 share is the same in every draw.",
    sprintf("- Draws: B = %s per region set, seed %d.", format(B_NULL, big.mark = ","), NULL_SEED),
    "- Per cell set and donor: Module post-IP% = Module post-IP cells in the set / the donor's Micro-PVM_3 cells x 100 (Part 10h denominator); fraction = Module post-IP cells in the set / cells in the set x 100, donors with at least 5 M4 cells; %MT = log1p of the mean per-cell %MT over the set (as Figure 14D).",
    "- Models at cov4, outcome and predictors standardized after the complete-case filter: MA outcome ~ Module post-IP% (as Figure 9A); MB outcome ~ Module post-IP% + %MT (as Figure 14D, BH over the nine clinical outcomes); MC outcome ~ M4 share + fraction (as Part 13 S08); MD outcome ~ fraction + %MT (as Part 13 S15).",
    "- Empirical p: p_upper = (1 + draws >= observed) / (B + 1); p_lower = (1 + draws <= observed) / (B + 1); two-sided p = min(1, 2 x the smaller).",
    "- BH FDR within region set: each model coefficient over the three cognition outcomes; the three collinearity readouts (r of Module post-IP% with %MT, r of the fraction with %MT, mean fraction) as one family. Every reading uses the FDR.",
    "",
    "Primary test: NeoC MMSE, the MA coefficient of the M4 cells (PVM_3 labels) against the random subsets.",
    "I1 Identity adds: FDR < 0.05 and the M4 coefficient lies further in the expected direction (negative) than the median of the random subsets.",
    "I2 Identity weakens: FDR < 0.05 and the M4 coefficient lies on the other side of that median.",
    "I0 No evidence that identity matters: FDR >= 0.05.",
    "",
    "Second test: NeoC, the donor-level Pearson r between Module post-IP% and %MT (donors with M4 cells), M4 cells (PVM_3 labels) against the random subsets.",
    "S1 M4 cells less collinear with %MT than random subsets: FDR < 0.05 and r below the median of the random subsets.",
    "S2 M4 cells more collinear with %MT: FDR < 0.05 and r above that median.",
    "S0 No difference: FDR >= 0.05.",
    "",
    "Reported without a reading: NonNeoC, CASI and cognitive status, MB, MC, MD, the mean fraction, the r of the fraction with %MT, and every comparison that uses the M4 arm's own labels.",
    "Limit: the random subsets are matched on the number of cells per donor and region only. A difference can reflect any property that separates Mic.13-like cells from other PVM_3 cells, for example sequencing depth."
)

# ============================================================================
run_part13b <- function() {
    # == STEP 1: pre-set readings =============================================
    cat("\nSTEP 1 -- pre-set readings\n")
    if (file.exists(PRESET) && !ALLOW_PRESET_REWRITE) {
        same <- identical(readLines(PRESET, warn = FALSE), PRESET_TEXT)
        gr13b("GR-13b-01", if (same) "PASS" else "FAIL",
              sprintf("pre-set readings %s the file written by an earlier run (%s)",
                      if (same) "match" else "differ from", basename(PRESET)))
    } else {
        writeLines(PRESET_TEXT, PRESET)
        gr13b("GR-13b-01", "PASS", sprintf("pre-set readings written before any model: %s",
                                           basename(PRESET)))
    }

    # == STEP 2: inputs =======================================================
    cat("\nSTEP 2 -- inputs\n")
    comp_path <- find_file(COMP_RDS, COMP_BASENAME)
    pm4_path  <- find_file(P10H_M4_RDS, P10H_M4_BASENAME)
    pp3_path  <- find_file(P10H_P3_RDS, P10H_P3_BASENAME)
    lab_path  <- if (!is.na(LAB_FILE) && file.exists(LAB_FILE)) normalizePath(LAB_FILE, winslash = "/") else
        find_file(NA, LAB_BASENAME)
    if (is.na(MT_FILE) || !file.exists(MT_FILE)) stop("MT_FILE not found: ", MT_FILE, call. = FALSE)
    inputs <- c(comparator = comp_path, prognosis_M4 = pm4_path, prognosis_PVM3 = pp3_path,
                labels = lab_path, mt = normalizePath(MT_FILE, winslash = "/"))
    for (k in names(inputs)) cat(sprintf("  %-15s %s\n", k, inputs[[k]]))
    RC  <- readRDS(comp_path)
    PM4 <- readRDS(pm4_path)
    PP3 <- readRDS(pp3_path)
    if (is.null(RC$donor_frame) || is.null(RC$meta)) stop("The comparator RDS carries no donor_frame or meta.", call. = FALSE)
    # Exact names only: `$` partial-matches and errors on a non-list field.
    get_in <- function(o, path) {
        for (k in path) { if (!is.list(o) || !k %in% names(o)) return(NULL); o <- o[[k, exact = TRUE]] }
        o
    }
    arm_tag_of <- function(o) {
        for (pth in list("arm_tag", c("meta", "arm_tag"), c("config", "ARM_TAG"), c("run", "arm_tag"))) {
            f <- get_in(o, pth)
            if (is.character(f) && length(f) == 1L) return(f)
        }
        NA_character_
    }
    t_m4 <- arm_tag_of(PM4); t_p3 <- arm_tag_of(PP3)
    arm_ok <- identical(RC$meta$arm_tag, ARM_TAG_EXPECT) &&
        (is.na(t_m4) || identical(t_m4, ARM_TAG_EXPECT)) &&
        (is.na(t_p3) || t_p3 %in% c("PVM3_all", "PVM3_bulk"))
    p10_same <- is.null(RC$meta$part10h_rds) || identical(basename(RC$meta$part10h_rds), basename(pm4_path))
    gr13b("GR-13b-02", if (!(arm_ok && p10_same)) "FAIL" else "PASS",
          sprintf(paste("comparator arm_tag '%s' (expected '%s'); M4 Prognosis arm tag '%s'; PVM_3 Prognosis arm tag '%s'",
                        "(NA = not stored; the file is identified by its basename %s); comparator built from '%s', this run",
                        "reads '%s'"),
                  RC$meta$arm_tag %||% "unrecorded", ARM_TAG_EXPECT, t_m4, t_p3, basename(pp3_path),
                  basename(RC$meta$part10h_rds %||% "unrecorded"), basename(pm4_path)))

    # == STEP 3: donor frames and denominators ================================
    cat("\nSTEP 3 -- donor frames and denominators\n")
    DF <- as.data.frame(RC$donor_frame)
    need <- c(".donor", "dataset", "pct_post", "pct_pre", "pct_excl_post", "pct_excl_pre", "n_total")
    miss <- setdiff(need, names(DF))
    if (length(miss)) stop("donor_frame lacks: ", paste(miss, collapse = ", "), call. = FALSE)
    DF$.donor <- as.character(DF$.donor); DF$dataset <- as.character(DF$dataset)
    for (k in need[-(1:2)]) DF[[k]] <- num(DF[[k]])
    OUT_MAP <- unlist(RC$meta$outcomes)
    COV4 <- unlist(RC$meta$covariate_columns$cov4)
    lost <- setdiff(CLIN_LAB, names(OUT_MAP))
    if (length(lost)) stop("Outcomes absent from the comparator's meta$outcomes: ", paste(lost, collapse = ", "), call. = FALSE)
    if (length(COV4) != 4L || !all(COV4 %in% names(DF)))
        stop("cov4 does not resolve to four donor_frame columns.", call. = FALSE)
    DF$share <- DF$pct_pre + DF$pct_post + DF$pct_excl_pre + DF$pct_excl_post
    k_raw <- DF$share * DF$n_total / 100
    a_raw <- DF$pct_post * DF$n_total / 100
    DF$k <- round(k_raw); DF$a_nat <- round(a_raw)
    n_ds <- as.integer(table(factor(DF$dataset, DATASETS)))
    dup <- any(duplicated(DF[, c(".donor", "dataset")]))
    tot_k <- as.integer(tapply(DF$k, factor(DF$dataset, DATASETS), sum))
    whole <- max(abs(c(k_raw - DF$k, a_raw - DF$a_nat)))
    gr13b("GR-13b-03", if (!dup && all(n_ds == N_DONORS_EXPECT) && isTRUE(whole <= 1e-6) &&
                           all(tot_k == N_M4_CELLS_EXPECT) && all(DF$n_total >= 1)) "PASS" else "FAIL",
          sprintf(paste("M4 donor frame: donors %s (expected %s), duplicated rows %s; M4 cells (share x n_total / 100)",
                        "%s (expected %s); M4 cells and Module post-IP cells are whole numbers (max deviation %.1e)"),
                  paste(n_ds, collapse = "/"), paste(N_DONORS_EXPECT, collapse = "/"), dup,
                  paste(tot_k, collapse = "/"), paste(N_M4_CELLS_EXPECT, collapse = "/"), whole))

    prop_of <- function(R, what) {
        out <- do.call(rbind, lapply(DATASETS, function(ds) {
            p <- R$datasets[[ds]]$prop_df
            if (is.null(p)) stop(what, ": datasets$", ds, "$prop_df is absent.", call. = FALSE)
            p <- as.data.frame(p)
            nd <- c("donor_id", "n_total", "pct_post", "pct_pre", "pct_excl_post", "pct_excl_pre")
            if (length(setdiff(nd, names(p)))) stop(what, " prop_df lacks: ", paste(setdiff(nd, names(p)), collapse = ", "),
                                                   call. = FALSE)
            data.frame(.donor = as.character(p$donor_id), dataset = ds, n_total = num(p$n_total),
                       pct_post = num(p$pct_post), pct_pre = num(p$pct_pre), pct_excl_post = num(p$pct_excl_post),
                       pct_excl_pre = num(p$pct_excl_pre), stringsAsFactors = FALSE)
        }))
        out$share <- out$pct_pre + out$pct_post + out$pct_excl_pre + out$pct_excl_post
        out
    }
    PF4 <- prop_of(PM4, "M4 Prognosis"); PF3 <- prop_of(PP3, "PVM_3 Prognosis")
    key <- function(d) paste(d$dataset, d$.donor, sep = "|")
    i4 <- match(key(DF), key(PF4)); i3 <- match(key(DF), key(PF3))
    same_set4 <- setequal(key(DF), key(PF4)); same_set3 <- setequal(key(DF), key(PF3))
    d_m4 <- if (anyNA(i4)) Inf else max(abs(c(DF$pct_post - PF4$pct_post[i4], DF$n_total - PF4$n_total[i4],
                                              DF$share - PF4$share[i4])))
    d_n3 <- if (anyNA(i3)) Inf else max(abs(DF$n_total - PF3$n_total[i3]))
    gr13b("GR-13b-04", if (same_set4 && same_set3 && isTRUE(d_m4 <= TOL_EXACT) && isTRUE(d_n3 <= TOL_EXACT)) "PASS" else "FAIL",
          sprintf(paste("one donor set and one denominator: comparator donor_frame and the M4 Prognosis prop_df agree on",
                        "donors (%s), pct_post, n_total and share (max |diff| %.2e); the PVM_3 Prognosis prop_df has the",
                        "same donors (%s) and the same Micro-PVM_3 denominator n_total (max |diff| %.2e)"),
                  same_set4, d_m4, same_set3, d_n3))
    DF$k3 <- round(PF3$share[i3] * PF3$n_total[i3] / 100)
    DF$a3_ref <- round(PF3$pct_post[i3] * PF3$n_total[i3] / 100)
    whole3 <- max(abs(c(PF3$share[i3] * PF3$n_total[i3] / 100 - DF$k3, PF3$pct_post[i3] * PF3$n_total[i3] / 100 - DF$a3_ref)))
    tot_k3 <- as.integer(tapply(DF$k3, factor(DF$dataset, DATASETS), sum))
    gr13b("GR-13b-05", if (isTRUE(whole3 <= 1e-6) && all(tot_k3 == N_P3_CELLS_EXPECT) && all(DF$k <= DF$k3)) "PASS" else "FAIL",
          sprintf(paste("PVM_3 arm cells from the PVM_3 prop_df: %s (expected %s), whole numbers (max deviation %.1e);",
                        "donors with more M4 cells than PVM_3 arm cells: %d"),
                  paste(tot_k3, collapse = "/"), paste(N_P3_CELLS_EXPECT, collapse = "/"), whole3, sum(DF$k > DF$k3)))

    # == STEP 4: the cell label file ==========================================
    cat("\nSTEP 4 -- cell label file\n")
    L <- if (requireNamespace("data.table", quietly = TRUE))
        as.data.frame(data.table::fread(lab_path, colClasses = "character", showProgress = FALSE, na.strings = c("", "NA")))
        else utils::read.csv(lab_path, colClasses = "character", check.names = FALSE, na.strings = c("", "NA"))
    cat(sprintf("  %d row(s), %d column(s): %s\n", nrow(L), ncol(L), paste(names(L), collapse = ", ")))
    pick_col <- function(given, cands, test, what) {
        if (!is.na(given)) {
            if (!given %in% names(L)) stop(what, " column '", given, "' is not in the label file.", call. = FALSE)
            return(given)
        }
        hit <- cands[cands %in% names(L)]
        hit <- hit[vapply(hit, function(h) isTRUE(test(L[[h]])), logical(1))]
        if (!length(hit)) return(NA_character_)
        hit[1]
    }
    all_donors <- unique(DF$.donor)
    C_DS <- pick_col(LAB_COLS$dataset, c("dataset", "region_set", "set", "Dataset"),
                     function(v) mean(v %in% DATASETS) > 0.99, "dataset")
    C_DN <- pick_col(LAB_COLS$donor, c("donor_id", "donor", "donor_key", "Donor ID", "DonorID"),
                     function(v) mean(unique(v) %in% all_donors) > 0.9, "donor")
    C_ARM <- pick_col(LAB_COLS$arm, c("arm", "arm_tag", "Arm", "ARM_TAG"),
                      function(v) any(v %in% unlist(ARM_VALUES)), "arm")
    C_MT <- pick_col(LAB_COLS$mito, c("mito_pct", "pct_mt", "percent.mt", "pct_mito", "mito_frac", "percent_mito"),
                     function(v) mean(is.finite(num(v))) > 0.99, "mito")
    if (any(is.na(c(C_DS, C_DN, C_ARM, C_MT))))
        stop(sprintf("label file: could not resolve dataset '%s', donor '%s', arm '%s', mito '%s'. Set LAB_COLS.",
                     C_DS, C_DN, C_ARM, C_MT), call. = FALSE)
    arm_vals <- unique(L[[C_ARM]])
    V_M4 <- intersect(ARM_VALUES$M4, arm_vals)[1]; V_P3 <- intersect(ARM_VALUES$PVM3, arm_vals)[1]
    cat(sprintf("  arm values: %s\n", paste(sprintf("%s (%d)", names(table(L[[C_ARM]])), as.integer(table(L[[C_ARM]]))),
                                            collapse = ", ")))
    if (is.na(V_M4) || is.na(V_P3))
        stop("label file: no M4 or no PVM_3 arm value among: ", paste(arm_vals, collapse = ", "), ". Set ARM_VALUES.",
             call. = FALSE)
    L <- L[L[[C_ARM]] %in% c(V_M4, V_P3) & L[[C_DS]] %in% DATASETS, , drop = FALSE]
    is_m4 <- L[[C_ARM]] == V_M4
    # Barcode: a named column first, else the unique column linking most M4 rows to PVM_3.
    bc_score <- function(h, strict = TRUE) {
        v <- L[[h]]; if (anyNA(v)) return(-1)
        k1 <- paste(L[[C_DS]], L[[C_ARM]], v)
        if (strict && anyDuplicated(k1)) return(-1)
        m4k <- paste(L[[C_DS]][is_m4], v[is_m4]); p3k <- paste(L[[C_DS]][!is_m4], v[!is_m4])
        mean(m4k %in% p3k)
    }
    C_BC <- LAB_COLS$barcode
    if (is.na(C_BC)) {
        named <- intersect(c("barcode", "cell_id", "obs_name", "cell", "Barcode", "obs_names", "index", "cell_barcode"), names(L))
        sc <- vapply(named, bc_score, numeric(1), strict = FALSE)   # duplicates checked below
        C_BC <- if (length(sc) && max(sc) > 0.5) named[which.max(sc)] else {
            other <- setdiff(names(L), c(C_DS, C_DN, C_ARM, C_MT))
            other <- other[vapply(other, function(h) length(unique(L[[h]])) > 0.5 * nrow(L) / 2, logical(1))]
            sc2 <- vapply(other, bc_score, numeric(1))
            if (length(sc2) && max(sc2) > 0.5) other[which.max(sc2)] else NA_character_
        }
    }
    if (is.na(C_BC)) stop("label file: no barcode column links the M4 rows to the PVM_3 rows. Set LAB_COLS$barcode.", call. = FALSE)
    REGION_NAMES <- c("MTG", "DFC", "A1C", "V1C", "AnG", "STG", "ITG", "FI", "HIP", "LEC", "MEC", "PFC", "SSp", "S1C", "M1C", "IPL", "CgG")
    C_RG <- pick_col(LAB_COLS$region, c("region", "Region", "brain_region", "roi", "ROI", "Brain Region", "region_name"),
                     function(v) length(unique(v)) >= 2 && !anyNA(v), "region")
    if (is.na(C_RG)) {
        other <- setdiff(names(L), c(C_DS, C_DN, C_ARM, C_MT, C_BC))
        rg <- other[vapply(other, function(h) { u <- unique(L[[h]]); length(u) >= 2 && length(u) <= 40 &&
            sum(u %in% REGION_NAMES) >= 2 && !anyNA(L[[h]]) }, logical(1))]
        if (length(rg)) C_RG <- rg[1]
    }
    use_region <- STRATA == "donor_region" && !is.na(C_RG)
    L$.ds <- L[[C_DS]]; L$.dn <- trimws(L[[C_DN]]); L$.bc <- L[[C_BC]]; L$.mt <- num(L[[C_MT]])
    L$.rg <- if (!is.na(C_RG)) L[[C_RG]] else "all"
    cat(sprintf("  columns: dataset '%s', donor '%s', arm '%s' (M4 '%s', PVM_3 '%s'), barcode '%s', %%MT '%s', region '%s'\n",
                C_DS, C_DN, C_ARM, V_M4, V_P3, C_BC, C_MT, if (is.na(C_RG)) "none" else C_RG))
    if ("definition" %in% names(L)) cat(sprintf("  definition: %s\n", paste(unique(L$definition), collapse = ", ")))
    if ("labels_version" %in% names(L)) cat(sprintf("  labels_version: %s\n", paste(unique(L$labels_version), collapse = ", ")))
    dupk <- sum(duplicated(paste(L$.ds, L[[C_ARM]], L$.bc)))
    bad_dn <- sum(!paste(L$.ds, L$.dn, sep = "|") %in% key(DF))
    bad_mt <- sum(!is.finite(L$.mt) | L$.mt < 0 | L$.mt > 100)
    gr13b("GR-13b-06", if (dupk == 0 && bad_dn == 0 && bad_mt == 0 &&
                           (STRATA == "donor" || use_region)) "PASS" else
              if (dupk == 0 && bad_dn == 0 && bad_mt == 0) "WARN" else "FAIL",
          sprintf(paste("label file rows: duplicated (region set, arm, barcode) %d; rows whose donor and region set are not",
                        "in the donor frame %d; %%MT missing or outside 0-100 %d; strata %s%s"),
                  dupk, bad_dn, bad_mt, if (use_region) "donor x region" else "donor",
                  if (STRATA == "donor_region" && !use_region) " (no region column; donor used)" else ""))

    # Cell counts per donor reproduce the Part 10h arm cells.
    cnt_rows <- function(sel) {
        t <- table(factor(paste(L$.ds[sel], L$.dn[sel], sep = "|"), levels = key(DF)))
        as.integer(t)
    }
    n4 <- cnt_rows(is_m4); n3 <- cnt_rows(!is_m4)
    d4 <- max(abs(n4 - DF$k)); d3 <- max(abs(n3 - DF$k3))
    gr13b("GR-13b-07", if (d4 == 0 && d3 == 0) "PASS" else "FAIL",
          sprintf(paste("label rows per donor equal the Part 10h arm cells on every donor: M4 rows %s against M4 share",
                        "x n_total %s (max |diff| %d); PVM_3 rows %s against the PVM_3 prop_df %s (max |diff| %d)"),
                  paste(tapply(n4, DF$dataset, sum)[DATASETS], collapse = "/"), paste(tot_k, collapse = "/"), d4,
                  paste(tapply(n3, DF$dataset, sum)[DATASETS], collapse = "/"), paste(tot_k3, collapse = "/"), d3))

    # == STEP 5: the Module post-IP label, found by reproduction ==============
    cat("\nSTEP 5 -- Module post-IP label (the rule that reproduces Part 10h pct_post on both arms)\n")
    excl <- c(C_DS, C_DN, C_ARM, C_BC, C_MT, C_RG, ".ds", ".dn", ".bc", ".mt", ".rg")
    lowcard <- setdiff(names(L), excl)
    lowcard <- lowcard[vapply(lowcard, function(h) { u <- unique(stats::na.omit(L[[h]])); length(u) >= 2 && length(u) <= 30 },
                              logical(1))]
    for (h in lowcard) {
        t <- sort(table(L[[h]][is_m4], useNA = "ifany"), decreasing = TRUE)
        cat(sprintf("  M4 rows, %s: %s\n", h, paste(sprintf("%s (%d)", names(t), as.integer(t))[seq_len(min(8, length(t)))],
                                                    collapse = ", ")))
    }
    TRUE_TOK <- c("true", "t", "1", "yes", "y")
    is_mod_val <- function(v) tolower(trimws(v)) %in% c(TRUE_TOK, "module", "in_module") |
        grepl("^(original \\+ refined|recovered)$", trimws(v), ignore.case = TRUE) |
        (grepl("module", v, ignore.case = TRUE) & !grepl("non|post|pre|not|excl", v, ignore.case = TRUE))
    is_post_val <- function(v, colname) (grepl("post", v, ignore.case = TRUE) & !grepl("pre", v, ignore.case = TRUE)) |
        (grepl("post", colname, ignore.case = TRUE) & tolower(trimws(v)) %in% TRUE_TOK)
    cands <- list()
    add_c <- function(desc, f) cands[[length(cands) + 1L]] <<- list(desc = desc, f = f)
    if (!is.null(LAB_RULE)) {
        if (!is.null(LAB_RULE$single)) {
            s <- LAB_RULE$single
            add_c(sprintf("%s in {%s} (set by hand)", s$col, paste(s$values, collapse = ", ")),
                  local({ cc <- s$col; vv <- s$values; function(D) !is.na(D[[cc]]) & D[[cc]] %in% vv }))
        } else {
            m <- LAB_RULE$module; s <- LAB_RULE$side
            add_c(sprintf("%s in {%s} and %s in {%s} (set by hand)", m$col, paste(m$values, collapse = ", "), s$col,
                          paste(s$values, collapse = ", ")),
                  local({ mc <- m$col; mv <- m$values; sc <- s$col; sv <- s$values
                      function(D) !is.na(D[[mc]]) & D[[mc]] %in% mv & !is.na(D[[sc]]) & D[[sc]] %in% sv }))
        }
    } else {
        mods <- list(); sides <- list()
        for (h in lowcard) {
            u <- unique(stats::na.omit(L[[h]]))
            sp <- u[grepl("post", u, ignore.case = TRUE) & grepl("mod", u, ignore.case = TRUE) &
                    !grepl("non", u, ignore.case = TRUE)]
            if (length(sp)) add_c(sprintf("%s in {%s}", h, paste(sp, collapse = ", ")),
                                  local({ cc <- h; vv <- sp; function(D) !is.na(D[[cc]]) & D[[cc]] %in% vv }))
            if (grepl("post", h, ignore.case = TRUE) && grepl("mod", h, ignore.case = TRUE) && !grepl("non", h, ignore.case = TRUE)) {
                tv <- u[tolower(trimws(u)) %in% TRUE_TOK]
                if (length(tv)) add_c(sprintf("%s in {%s}", h, paste(tv, collapse = ", ")),
                                      local({ cc <- h; vv <- tv; function(D) !is.na(D[[cc]]) & D[[cc]] %in% vv }))
            }
            mv <- u[is_mod_val(u)]; if (length(mv) && length(mv) < length(u)) mods[[h]] <- mv
            pv <- u[is_post_val(u, h)]; if (length(pv) && length(pv) < length(u)) sides[[h]] <- pv
        }
        num_cols <- setdiff(names(L), excl)
        cd <- num_cols[grepl("cdfg", num_cols, ignore.case = TRUE)]
        th <- num_cols[grepl("^ip_y$|^threshold$|^IP_y$", num_cols, ignore.case = TRUE)]
        for (a in names(mods)) {
            for (b in setdiff(names(sides), a))
                add_c(sprintf("%s in {%s} and %s in {%s}", a, paste(mods[[a]], collapse = ", "), b,
                              paste(sides[[b]], collapse = ", ")),
                      local({ mc <- a; mv <- mods[[a]]; sc <- b; sv <- sides[[b]]
                          function(D) !is.na(D[[mc]]) & D[[mc]] %in% mv & !is.na(D[[sc]]) & D[[sc]] %in% sv }))
            for (c1 in cd) for (t1 in th)
                add_c(sprintf("%s in {%s} and %s >= %s", a, paste(mods[[a]], collapse = ", "), c1, t1),
                      local({ mc <- a; mv <- mods[[a]]; cc <- c1; tc <- t1
                          function(D) { x <- num(D[[cc]]); y <- num(D[[tc]])
                              !is.na(D[[mc]]) & D[[mc]] %in% mv & is.finite(x) & is.finite(y) & x >= y } }))
        }
    }
    if (!length(cands)) stop("No candidate Module post-IP rule in the label file. Set LAB_RULE.",
                             call. = FALSE)
    a_nat_ref <- DF$a_nat; a3_ref <- DF$a3_ref
    ev <- lapply(cands, function(cn) {
        lv <- cn$f(L); lv[is.na(lv)] <- FALSE
        c4 <- cnt_rows(is_m4 & lv); c3 <- cnt_rows(!is_m4 & lv)
        list(desc = cn$desc, lv = lv, dev4 = max(abs(c4 - a_nat_ref)), dev3 = max(abs(c3 - a3_ref)))
    })
    CAND <- data.frame(rule = vapply(ev, `[[`, "", "desc"), max_dev_M4 = vapply(ev, `[[`, 0, "dev4"),
                       max_dev_PVM3 = vapply(ev, `[[`, 0, "dev3"), stringsAsFactors = FALSE)
    CAND$reproduces <- CAND$max_dev_M4 == 0 & CAND$max_dev_PVM3 == 0
    print(CAND, row.names = FALSE)
    okc <- which(CAND$reproduces)
    distinct <- if (length(okc)) okc[!duplicated(lapply(ev[okc], function(e) which(e$lv)))] else integer(0)
    gr13b("GR-13b-08", if (length(distinct) == 1L) "PASS" else "FAIL",
          sprintf(paste("Module post-IP label: %d candidate rule(s), %d reproduce the Part 10h Module post-IP cells",
                        "(pct_post x n_total / 100) on every donor of BOTH arms, %d distinct cell set(s)%s"),
                  nrow(CAND), length(okc), length(distinct),
                  if (length(distinct) == 1L) paste0("; rule used: ", CAND$rule[distinct]) else
                      "; set LAB_RULE"))
    L$.post <- ev[[distinct]]$lv
    RULE_USED <- CAND$rule[distinct]

    # Cell-level check: Module class from cell_class, side from group_refined,
    # then IP_status.
    V1TAB <- NULL
    if (isTRUE(V1_CHECK)) {
        side_of <- function(z) { z <- as.character(z)
            ifelse(is.na(z), NA_character_, ifelse(grepl("Post-IP", z, fixed = TRUE), "Post-IP",
                   ifelse(grepl("Pre-IP", z, fixed = TRUE), "Pre-IP", NA_character_))) }
        for (a in names(V1_TAGS)) for (ds in DATASETS) {
            nm <- paste(a, ds, sep = "_")
            base <- sprintf("Part10e_CP10K_%s_V1_%s_SK_interaction_analysis.xlsx", V1_TAGS[[a]], ds)
            f <- V1_FILES[[nm]] %||% tryCatch(find_file(NA, base, required = FALSE), error = function(e) "AMBIGUOUS")
            if (is.na(f) || identical(f, "AMBIGUOUS")) {
                V1TAB <- rbind(V1TAB, data.frame(arm = a, dataset = ds, file = base,
                    found = FALSE, status = if (is.na(f)) "not found" else "several copies; set V1_FILES",
                    rows = NA, matched = NA, unlabelled = NA, mismatched = NA, stringsAsFactors = FALSE)); next }
            w <- openxlsx::read.xlsx(f, sheet = "Step4_Refined_Cells")
            if (!all(c("barcode", "cell_class") %in% names(w))) stop(base, ": Step4_Refined_Cells lacks barcode or cell_class.", call. = FALSE)
            sd1 <- if ("group_refined" %in% names(w)) side_of(w$group_refined) else rep(NA_character_, nrow(w))
            if ("IP_status" %in% names(w)) { s2 <- side_of(w$IP_status); sd1[is.na(sd1)] <- s2[is.na(sd1)] }
            mod <- w$cell_class %in% V1_MODULE_CLASSES
            known <- !mod | !is.na(sd1)                       # non-Module cells are never Module post-IP
            v1post <- mod & sd1 %in% "Post-IP"
            sel <- which(L$.ds == ds & L[[C_ARM]] == (if (a == "M4") V_M4 else V_P3))
            jx <- match(L$.bc[sel], as.character(w$barcode))
            ok <- !is.na(jx) & known[pmax(jx, 1L)]
            V1TAB <- rbind(V1TAB, data.frame(arm = a, dataset = ds, file = basename(f), found = TRUE, status = "read",
                rows = length(sel),
                matched = sum(!is.na(jx)), unlabelled = sum(!is.na(jx) & !ok),
                mismatched = sum(L$.post[sel][ok] != v1post[jx[ok]]), stringsAsFactors = FALSE))
        }
        print(V1TAB, row.names = FALSE)
    }
    if (is.null(V1TAB) || !any(V1TAB$found)) {
        gr13b("GR-13b-19", "WARN", paste("Validation 1 workbooks not read (not found, several copies, or V1_CHECK = FALSE);",
                                         "labels checked by per-donor counts only"))
    } else {
        f1 <- V1TAB[V1TAB$found, ]
        st <- if (any(f1$mismatched > 0)) "FAIL" else
            if (!all(V1TAB$found) || any(f1$matched < f1$rows) || any(f1$unlabelled > 0)) "WARN" else "PASS"
        gr13b("GR-13b-19", st,
              sprintf(paste("cell-level labels against Part 10h's Validation 1 cell table (Module class + inherited side):",
                            "%d of 4 workbooks found; label cells matched by barcode %d of %d; without an inherited side %d;",
                            "cells whose Module post-IP label differs %d"),
                      sum(V1TAB$found), sum(f1$matched), sum(f1$rows), sum(f1$unlabelled), sum(f1$mismatched)))
    }

    # == STEP 6: the M4 cells inside the PVM_3 pool ===========================
    cat("\nSTEP 6 -- the M4 cells inside the PVM_3 pool\n")
    M4 <- L[is_m4, c(".ds", ".dn", ".bc", ".mt", ".rg", ".post")]
    P3 <- L[!is_m4, c(".ds", ".dn", ".bc", ".mt", ".rg", ".post")]
    j <- match(paste(M4$.ds, M4$.bc), paste(P3$.ds, P3$.bc))
    n_miss <- sum(is.na(j)); jj <- j[!is.na(j)]; mm <- M4[!is.na(j), ]
    dn_bad <- sum(mm$.dn != P3$.dn[jj]); rg_bad <- sum(mm$.rg != P3$.rg[jj])
    mt_dev <- if (length(jj)) max(abs(mm$.mt - P3$.mt[jj])) else NA_real_
    gr13b("GR-13b-09", if (n_miss == 0 && dn_bad == 0 && rg_bad == 0 && isTRUE(mt_dev <= 1e-12)) "PASS" else "FAIL",
          sprintf(paste("every M4 cell is a PVM_3 pool cell of the same region set: %d of %d M4 cells absent from the",
                        "pool; of those found, %d carry another donor, %d another region, and the per-cell",
                        "%%MT differs by at most %.1e"),
                  n_miss, nrow(M4), dn_bad, rg_bad, mt_dev))
    M4$p3_row <- j
    M4$post3 <- ifelse(is.na(j), NA, P3$.post[pmax(j, 1L)])
    AGREE <- do.call(rbind, lapply(DATASETS, function(ds) {
        s <- M4$.ds == ds & !is.na(M4$post3)
        tb <- table(factor(M4$.post[s], c(FALSE, TRUE)), factor(M4$post3[s], c(FALSE, TRUE)))
        po <- sum(diag(tb)) / sum(tb); pe <- sum(rowSums(tb) * colSums(tb)) / sum(tb)^2
        data.frame(dataset = ds, m4_cells = sum(s), post_both = tb[2, 2], post_M4_only = tb[2, 1],
                   post_PVM3_only = tb[1, 2], post_neither = tb[1, 1], agreement = po,
                   kappa = (po - pe) / (1 - pe), stringsAsFactors = FALSE)
    }))
    cat("  M4 cells: Module post-IP under the M4 arm's own labels against the PVM_3 labels\n")
    print(transform(AGREE, agreement = fmt(agreement), kappa = fmt(kappa)), row.names = FALSE)
    # Agreement of the PVM_3 labels with the M4 arm's own labels on the M4 cells.
    gr13b("GR-13b-20", if (all(is.finite(AGREE$kappa)) && all(AGREE$kappa >= 0.4)) "PASS" else "WARN",
          sprintf(paste("M4 cells, Module post-IP under the M4 arm's own labels against the PVM_3 labels: agreement %s,",
                        "kappa %s (NeoC/NonNeoC; WARN below 0.4)"),
                  paste(fmt(AGREE$agreement), collapse = "/"), paste(fmt(AGREE$kappa), collapse = "/")))

    # == STEP 7: %MT of the M4 cells ==========================================
    cat("\nSTEP 7 -- %MT of the M4 cells\n")
    mt_raw <- utils::read.csv(MT_FILE, check.names = FALSE, stringsAsFactors = FALSE)
    if (length(setdiff(unlist(MT_COLS), names(mt_raw)))) stop("MT_FILE lacks: ",
        paste(setdiff(unlist(MT_COLS), names(mt_raw)), collapse = ", "), call. = FALSE)
    MTF <- data.frame(k = paste(mt_raw[[MT_COLS$dataset]], mt_raw[[MT_COLS$donor]], sep = "|"),
                      v = num(mt_raw[[MT_COLS$value]]), stringsAsFactors = FALSE)
    if (anyDuplicated(MTF$k)) stop("MT_FILE has duplicated donor rows.", call. = FALSE)
    mt_lab <- tapply(M4$.mt, factor(paste(M4$.ds, M4$.dn, sep = "|"), levels = key(DF)), mean)
    DF$mt_lab <- ifelse(DF$k > 0, log1p(as.numeric(mt_lab)), NA_real_)
    DF$mt_file <- MTF$v[match(key(DF), MTF$k)]
    both <- is.finite(DF$mt_lab) & is.finite(DF$mt_file)
    mt_d <- if (any(both)) max(abs(DF$mt_lab[both] - DF$mt_file[both])) else NA_real_
    only_lab <- sum(is.finite(DF$mt_lab) & !is.finite(DF$mt_file)); only_file <- sum(!is.finite(DF$mt_lab) & is.finite(DF$mt_file))
    extra <- sum(!MTF$k %in% key(DF))
    gr13b("GR-13b-10", if (isTRUE(mt_d <= 1e-9) && only_lab == 0 && only_file == 0 && extra == 0) "PASS" else "FAIL",
          sprintf(paste("log1p of the mean per-cell %%MT over the M4 label rows reproduces the Figure 14D %%MT file on %d",
                        "donors (max |diff| %.1e); donors with a value on one side only: %d / %d; file rows matching no donor %d"),
                  sum(both), mt_d, only_lab, only_file, extra))
    DF$mt <- DF$mt_lab

    # == STEP 8: pool and observed cell sets ==================================
    cat("\nSTEP 8 -- pool, strata and observed cell sets\n")
    PREP <- list(); obs_chk <- numeric(0); cap_bad <- 0L
    for (ds in DATASETS) {
        d <- DF[DF$dataset == ds, , drop = FALSE]
        D <- d$.donor; nD <- length(D)
        p <- P3[P3$.ds == ds, , drop = FALSE]
        p$dix <- match(p$.dn, D)
        skey <- if (use_region) paste(p$.dn, p$.rg, sep = "|") else p$.dn
        lev <- sort(unique(skey))
        p$six <- match(skey, lev)
        S <- length(lev)
        o <- order(p$six); p <- p[o, , drop = FALSE]          # pool sorted by stratum
        row_in_p <- match(seq_len(nrow(P3)), which(P3$.ds == ds)[o])  # P3 row -> sorted pool row
        m <- M4[M4$.ds == ds & !is.na(M4$p3_row), , drop = FALSE]
        m_idx <- row_in_p[m$p3_row]
        size_s <- tabulate(p$six, S); k_s <- tabulate(p$six[m_idx], S)
        cap_bad <- cap_bad + sum(k_s > size_s)
        PREP[[ds]] <- list(d = d, D = D, nD = nD, p = p, S = S, size_s = size_s, k_s = k_s, m_idx = m_idx)
        obs <- agg_set(PREP[[ds]], m_idx)
        # The same quantities by a direct count, to check agg_set on the observed set.
        mm4 <- M4[M4$.ds == ds, , drop = FALSE]
        a_dir <- as.integer(table(factor(mm4$.dn[mm4$post3 %in% TRUE], levels = D)))
        mt_dir <- tapply(mm4$.mt, factor(mm4$.dn, levels = D), mean)
        obs_chk <- c(obs_chk, max(abs(obs$cnt - d$k)), max(abs(obs$x - 100 * a_dir / d$n_total)),
                     max(abs(obs$share - d$share)), max(abs(obs$mt - log1p(as.numeric(mt_dir))), na.rm = TRUE),
                     max(abs(obs$mt - d$mt), na.rm = TRUE), sum(is.finite(obs$mt) != is.finite(d$mt)))
        nat <- list(x = d$pct_post, frac = ifelse(d$k >= MIN_M4_CELLS_FRAC, 100 * d$pct_post / d$share, NA_real_),
                    mt = d$mt, share = d$share)
        PREP[[ds]]$obs <- obs; PREP[[ds]]$nat <- nat
        cat(sprintf("  %-8s pool %d cells, %d strata, %d M4 cells (%d donors with M4 cells); pool cells per stratum median %d\n",
                    ds, nrow(p), S, length(m_idx), sum(d$k > 0), as.integer(stats::median(size_s))))
    }
    # Expected share of drawn cells that are M4 cells: sum(k^2 / pool) / sum(k).
    ovl <- vapply(DATASETS, function(ds) { P <- PREP[[ds]]; i <- P$k_s > 0
        sum(P$k_s[i]^2 / P$size_s[i]) / sum(P$k_s[i]) }, 0)
    full <- vapply(DATASETS, function(ds) { P <- PREP[[ds]]; sum(P$k_s > 0 & P$k_s == P$size_s) }, 0)
    gr13b("GR-13b-21", if (all(ovl <= 0.4)) "PASS" else "WARN",
          sprintf(paste("expected share of drawn cells that are M4 cells %s (NeoC/NonNeoC; WARN above 0.4); strata in",
                        "which every pool cell is an M4 cell %s"),
                  paste(fmt(ovl), collapse = "/"), paste(full, collapse = "/")))
    gr13b("GR-13b-11", if (cap_bad == 0 && isTRUE(max(obs_chk) <= 1e-9)) "PASS" else "FAIL",
          sprintf(paste("observed cell sets built by the draw function equal a direct count of the M4 rows: M4 cells,",
                        "Module post-IP%% under PVM_3 labels, M4 share, %%MT and the %%MT donor pattern (max |diff| %.1e);",
                        "strata with more M4 cells than pool cells: %d"), max(obs_chk), cap_bad))

    # == STEP 9: models on the observed data ==================================
    cat("\nSTEP 9 -- models on the observed cell sets\n")
    MODELS <- list(MA = list(preds = "x", outs = COG_LAB, label = "Module post-IP% (as Figure 9A)"),
                   MB = list(preds = c("x", "mt"), outs = CLIN_LAB, label = "Module post-IP% + %MT (as Figure 14D)"),
                   MC = list(preds = c("share", "frac"), outs = COG_LAB, label = "M4 share + fraction (as Part 13 S08)"),
                   MD = list(preds = c("frac", "mt"), outs = COG_LAB, label = "fraction + %MT (as Part 13 S15)"))
    prep_y <- function(d) lapply(CLIN_LAB, function(lab) {
        y <- num(d[[OUT_MAP[[lab]]]]); covs <- setdiff(COV4, OUT_MAP[[lab]])
        Cm <- as.matrix(as.data.frame(lapply(d[, covs, drop = FALSE], num)))
        list(y = y, C = Cm, base = is.finite(y) & (rowSums(!is.finite(Cm)) == 0))
    })
    for (ds in DATASETS) { PREP[[ds]]$Y <- prep_y(PREP[[ds]]$d); names(PREP[[ds]]$Y) <- CLIN_LAB }
    OBS <- list(); NAT <- list(); FITCHK <- list()
    for (ds in DATASETS) {
        P <- PREP[[ds]]
        OBS[[ds]] <- readout(P$Y, P$obs); NAT[[ds]] <- readout(P$Y, P$nat)
        # Fast fitter against lm() on every observed model.
        for (set in c("obs", "nat")) for (mn in names(MODELS)) for (lab in MODELS[[mn]]$outs) {
            V <- P[[set]]
            dd <- P$d; for (pp in MODELS[[mn]]$preds) dd[[paste0(".p_", pp)]] <- V[[pp]]
            rf <- ref_fit(dd, OUT_MAP[[lab]], paste0(".p_", MODELS[[mn]]$preds), COV4)
            fr <- fit_one(P$Y, lab, MODELS[[mn]]$preds, V)
            it <- paste(ds, set, mn, lab)
            if (is.null(rf) || is.null(fr)) { FITCHK[[length(FITCHK) + 1L]] <- is.null(rf) && is.null(fr); next }
            FITCHK[[length(FITCHK) + 1L]] <- c(
                gcmp("GR-13b-12", paste(it, "n"), fr$n, rf$n, 0),
                gcmp("GR-13b-12", paste(it, "beta"), fr$beta, rf$beta, 1e-10),
                gcmp("GR-13b-12", paste(it, "se"), fr$se, rf$se, 1e-10),
                gcmp("GR-13b-12", paste(it, "p"), fr$p, rf$p, 1e-10),
                gcmp("GR-13b-12", paste(it, "ci"), c(fr$ci_lo, fr$ci_hi), c(rf$ci_lo, rf$ci_hi), 1e-10),
                gcmp("GR-13b-12", paste(it, "vif"), fr$vif, rf$vif, 1e-6))
        }
    }
    fc <- unlist(FITCHK)
    gr13b("GR-13b-12", if (length(fc) && all(fc)) "PASS" else "FAIL",
          sprintf("fast fitter equals lm() with the comparator .z and .safe_vif on every observed model: %d of %d values",
                  sum(fc), length(fc)))

    # Reproduction of stored and published M4 values with the M4 arm's own labels.
    ok <- logical(0)
    for (ds in DATASETS) {
        mv <- PM4$datasets[[ds]]$mv_df
        if (is.null(mv)) next
        mv <- mv[mv$tier == "all" & mv$covariate_set == "cov4" & mv$predictor == "pct_post" & mv$outcome %in% COG_LAB, , drop = FALSE]
        for (jr in seq_len(nrow(mv))) {
            r <- fit_one(PREP[[ds]]$Y, mv$outcome[jr], "x", PREP[[ds]]$nat); it <- paste(ds, mv$outcome[jr], "Figure 9A")
            ok <- c(ok, gcmp("GR-13b-13", paste(it, "n"), r$n, mv$n_donors[jr], 0),
                    gcmp("GR-13b-13", paste(it, "beta"), r$beta[1], mv$beta_std[jr], 1e-6))
        }
    }
    n9a <- length(ok)
    if (!is.null(PUB_14D) && length(PUB_14D)) for (v in PUB_14D) {
        ds <- v[1]; lab <- v[2]; x <- as.numeric(v[3:7]); it <- paste(ds, lab, "Figure 14D")
        Y <- PREP[[ds]]$Y; V <- PREP[[ds]]$nat
        bb <- lapply(CLIN_LAB, function(l) fit_one(Y, l, c("x", "mt"), V))
        px <- vapply(bb, function(r) if (is.null(r)) NA_real_ else r$p[1], 0)
        q <- px; q[is.finite(px)] <- stats::p.adjust(px[is.finite(px)], "BH")
        r <- bb[[match(lab, CLIN_LAB)]]
        ok <- c(ok, gcmp("GR-13b-13", paste(it, "beta"), r$beta[1], x[1], TOL_PUB),
                gcmp("GR-13b-13", paste(it, "CI low"), r$ci_lo[1], x[2], TOL_PUB),
                gcmp("GR-13b-13", paste(it, "CI high"), r$ci_hi[1], x[3], TOL_PUB),
                gcmp("GR-13b-13", paste(it, "n"), r$n, x[4], 0),
                gcmp("GR-13b-13", paste(it, "FDR"), q[match(lab, CLIN_LAB)], x[5], TOL_PUB))
    }
    gr13b("GR-13b-13", if (n9a > 0 && length(ok) > n9a && all(ok)) "PASS" else "FAIL",
          sprintf(paste("with the M4 arm's own labels this script reproduces Part 10h Multivariate (cov4, pct_post,",
                        "cognition; %d values) and the published Figure 14D rows (%d values): %d of %d"),
                  n9a, length(ok) - n9a, sum(ok), length(ok)))

    # == STEP 10: random subsets ==============================================
    cat(sprintf("\nSTEP 10 -- random subsets: B = %d per region set, seed %d\n", B_NULL, NULL_SEED))
    # One L'Ecuyer-CMRG stream per draw, NeoC draws first, then NonNeoC.
    old_kind <- RNGkind()
    RNGkind("L'Ecuyer-CMRG"); set.seed(NULL_SEED)
    st <- .Random.seed; SEEDS <- vector("list", length(DATASETS) * B_NULL)
    for (i in seq_along(SEEDS)) { SEEDS[[i]] <- st; st <- parallel::nextRNGStream(st) }
    obs_ns <- lapply(OBS, attr, "ns")
    n_cores <- max(1L, min(as.integer(N_CORES), B_NULL), na.rm = TRUE)
    cl <- NULL
    if (n_cores > 1L) {
        cl <- parallel::makeCluster(n_cores, type = "PSOCK")
        on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
        parallel::clusterExport(cl, WORKER_OBJECTS, envir = globalenv())
    }
    cat(sprintf("  %d worker(s) (%d logical cores detected)\n", n_cores, parallel::detectCores()))
    NULLS <- list(); valid <- list(); RAW <- list()
    inv <- c(count = 0, dup = 0, foreign = 0, nmis = 0, nfail = 0)
    t0 <- Sys.time()
    for (k in seq_along(DATASETS)) {
        ds <- DATASETS[k]; P <- PREP[[ds]]
        sd_ds <- SEEDS[(k - 1L) * B_NULL + seq_len(B_NULL)]
        res <- if (is.null(cl)) null_chunk(sd_ds, P, obs_ns[[ds]]) else {
            ch <- split(sd_ds, cut(seq_len(B_NULL), min(B_NULL, n_cores * 4L), labels = FALSE))
            unlist(parallel::parLapply(cl, ch, null_chunk, P = P, ns_obs = obs_ns[[ds]]), recursive = FALSE)
        }
        if (length(res) != B_NULL) stop("STEP 10: ", length(res), " draws returned for ", ds, ", expected ", B_NULL, call. = FALSE)
        RAW[[ds]] <- res
        NULLS[[ds]] <- do.call(rbind, lapply(res, `[[`, "r"))
        valid[[ds]] <- vapply(res, function(z) z$nm == 0 && z$nf == 0, logical(1))
        inv["count"] <- max(inv["count"], vapply(res, `[[`, 0, "dev"))
        inv["dup"] <- inv["dup"] + sum(vapply(res, `[[`, FALSE, "dup"))
        inv["foreign"] <- inv["foreign"] + sum(vapply(res, `[[`, 0, "foreign"))
        inv["nmis"] <- inv["nmis"] + sum(vapply(res, `[[`, 0, "nm"))
        inv["nfail"] <- inv["nfail"] + sum(vapply(res, `[[`, 0, "nf"))
        cat(sprintf("    %s: %d draws (%.0f s)\n", ds, B_NULL, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    }
    # The first draws again in this session; the workers must give the same numbers.
    n_chk <- min(20L, B_NULL); rep_dev <- 0
    for (k in seq_along(DATASETS)) {
        ds <- DATASETS[k]
        for (b in seq_len(n_chk)) {
            a <- null_one(SEEDS[[(k - 1L) * B_NULL + b]], PREP[[ds]], obs_ns[[ds]])$r
            z <- RAW[[ds]][[b]]$r
            # Same names, NA and Inf patterns; finite values within TOL_REPRO.
            fa <- is.finite(a); fz <- is.finite(z)
            same <- identical(names(a), names(z)) && identical(is.na(a), is.na(z)) && identical(fa, fz) &&
                all(a[is.infinite(a)] == z[is.infinite(z)])
            rep_dev <- max(rep_dev, if (!same) Inf else if (any(fa)) max(abs(a[fa] - z[fa])) else 0)
        }
    }
    do.call(RNGkind, as.list(old_kind))
    gr13b("GR-13b-22", if (rep_dev <= TOL_REPRO) "PASS" else "FAIL",
          sprintf(paste("parallel draws reproduce exactly when rerun one by one in this session (first %d draws per region",
                        "set: same NA and Inf pattern, finite values max |diff| %s, tolerance %s); %d worker(s), one",
                        "L'Ecuyer-CMRG stream per draw from seed %d"),
                  n_chk, format(rep_dev), format(TOL_REPRO), n_cores, NULL_SEED))
    obs_fail <- sum(vapply(OBS, attr, 0L, "nfail"))
    vmin <- min(vapply(valid, mean, 0))
    gr13b("GR-13b-14", if (inv["count"] == 0 && inv["dup"] == 0 && inv["foreign"] == 0 && obs_fail == 0) "PASS" else "FAIL",
          sprintf(paste("every draw holds each donor's M4 cell number in every stratum, so the M4 share is identical",
                        "(max |deviation| %d over %d draws x 2 region sets); draws with duplicated cells %d; out-of-pool cells",
                        "%d; observed models not fitted %d"),
                  as.integer(inv["count"]), B_NULL, as.integer(inv["dup"]), as.integer(inv["foreign"]), obs_fail))
    gr13b("GR-13b-15", if (vmin >= NULL_MIN_VALID) "PASS" else "FAIL",
          sprintf(paste("draws in which every model was fitted on exactly the observed donors: %s%% at the lower region",
                        "set (minimum %s%%); model-n mismatches %d, failed fits %d; %.0f s"),
                  fmt(100 * vmin, 2), fmt(100 * NULL_MIN_VALID, 0), as.integer(inv["nmis"]), as.integer(inv["nfail"]),
                  as.numeric(difftime(Sys.time(), t0, units = "secs"))))

    # == STEP 11: summaries and readings ======================================
    cat("\nSTEP 11 -- observed against the random subsets\n")
    RO <- names(OBS[[1]])
    cont <- RO[!grepl("^cB_", RO)]
    COLL_RO <- c("r_x_mt", "r_frac_mt", "mean_frac")
    # Readouts that are tested against the null; VIFs and within-draw FDRs are descriptive.
    is_test <- function(ro) grepl("^(bA|bB|bC|bD|bDm)_", ro) | ro %in% COLL_RO
    emp <- function(nul, o) {
        if (!is.finite(o)) return(c(NA, NA, NA))
        nul <- nul[is.finite(nul)]; B <- length(nul)
        pu <- (1 + sum(nul >= o)) / (B + 1); pl <- (1 + sum(nul <= o)) / (B + 1)
        c(pu, pl, min(1, 2 * min(pu, pl)))
    }
    SUM <- do.call(rbind, lapply(DATASETS, function(ds) do.call(rbind, lapply(cont, function(ro) {
        nul <- NULLS[[ds]][valid[[ds]], ro]; nul_f <- nul[is.finite(nul)]
        o <- OBS[[ds]][[ro]]; nt <- NAT[[ds]][[ro]]
        e1 <- if (is_test(ro)) emp(nul, o) else c(NA, NA, NA)
        lab <- names(SHORT)[match(sub("^[^_]+_", "", ro), SHORT)]
        data.frame(dataset = ds, readout = ro, outcome = if (is.na(lab)) NA_character_ else lab,
                   observed_M4_PVM3labels = o, M4_own_labels = nt, null_mean = mean(nul_f), null_sd = stats::sd(nul_f),
                   null_median = stats::median(nul_f), null_q025 = unname(stats::quantile(nul_f, .025)),
                   null_q975 = unname(stats::quantile(nul_f, .975)),
                   p_upper = e1[1], p_lower = e1[2], p_two = e1[3],
                   draws_used = length(nul_f), stringsAsFactors = FALSE)
    }))))
    SUM$bh_family <- ifelse(SUM$readout %in% COLL_RO, paste(SUM$dataset, "collinearity"),
                     ifelse(is_test(SUM$readout), paste(SUM$dataset, sub("_[^_]+$", "", SUM$readout)), NA_character_))
    SUM$fdr <- NA_real_
    for (f in unique(stats::na.omit(SUM$bh_family))) {
        i <- which(SUM$bh_family %in% f & is.finite(SUM$p_two))
        if (length(i)) SUM$fdr[i] <- stats::p.adjust(SUM$p_two[i], "BH")
    }
    fam_n <- table(SUM$bh_family[is.finite(SUM$p_two)])
    if (any(fam_n != 3L)) stop("BH family sizes differ from 3: ", paste(names(fam_n), fam_n, collapse = "; "), call. = FALSE)
    VER <- do.call(rbind, lapply(DATASETS, function(ds) do.call(rbind, lapply(COG_LAB, function(lab) {
        ro <- paste0("cB_", SHORT[[lab]]); nul <- NULLS[[ds]][valid[[ds]], ro]
        fr <- as.numeric(table(factor(nul, 1:5))) / length(nul)
        data.frame(dataset = ds, outcome = lab, model = "MB (as Figure 14D)",
                   observed_M4_PVM3labels = VERDICTS[OBS[[ds]][[ro]]], M4_own_labels = VERDICTS[NAT[[ds]][[ro]]],
                   null_module_only = fr[1], null_mt_only = fr[2], null_both = fr[3], null_neither = fr[4],
                   null_unreadable = fr[5], stringsAsFactors = FALSE)
    }))))
    show <- SUM[!is.na(SUM$outcome) & grepl("^b", SUM$readout) | SUM$readout %in% c("r_x_mt", "r_frac_mt", "mean_frac"), ]
    print(transform(show[, c("dataset", "readout", "observed_M4_PVM3labels", "M4_own_labels", "null_median", "null_q025",
                             "null_q975", "p_two", "fdr")],
                    observed_M4_PVM3labels = fmt(observed_M4_PVM3labels), M4_own_labels = fmt(M4_own_labels),
                    null_median = fmt(null_median), null_q025 = fmt(null_q025), null_q975 = fmt(null_q975),
                    p_two = fmt_p(p_two), fdr = fmt_p(fdr)), row.names = FALSE)
    print(transform(VER, null_module_only = fmt(null_module_only, 3), null_mt_only = fmt(null_mt_only, 3),
                    null_both = fmt(null_both, 3), null_neither = fmt(null_neither, 3),
                    null_unreadable = fmt(null_unreadable, 3)), row.names = FALSE)

    row_of <- function(ds, ro) SUM[SUM$dataset == ds & SUM$readout == ro, , drop = FALSE]
    pr <- row_of("NeoC", "bA_MMSE")
    if (nrow(pr) != 1L || !is.finite(pr$observed_M4_PVM3labels)) {
        reading_I <- "not computable (the primary model was not fitted)"; i_ok <- FALSE
    } else {
        stronger <- EXP_SIGN[["MMSE score"]] * (pr$observed_M4_PVM3labels - pr$null_median) > 0
        reading_I <- if (pr$fdr < ALPHA && stronger) "I1 Identity adds" else
            if (pr$fdr < ALPHA) "I2 Identity weakens" else "I0 No evidence that identity matters"
        i_ok <- TRUE
    }
    sr <- row_of("NeoC", "r_x_mt")
    if (nrow(sr) != 1L || !is.finite(sr$observed_M4_PVM3labels)) {
        reading_S <- "not computable"; s_ok <- FALSE
    } else {
        lower <- sr$observed_M4_PVM3labels < sr$null_median
        reading_S <- if (sr$fdr < ALPHA && lower) "S1 M4 cells less collinear with %MT than random subsets" else
            if (sr$fdr < ALPHA) "S2 M4 cells more collinear with %MT" else "S0 No difference"
        s_ok <- TRUE
    }
    gr13b("GR-13b-16", if (i_ok) "PASS" else "FAIL",
          sprintf(paste("Primary reading (NeoC MMSE, MA): M4 cells beta %s against random subsets median %s",
                        "[95%% of draws %s to %s], two-sided p %s, FDR %s. %s. Reference without a test: the M4 cells",
                        "under the M4 arm's own labels give beta %s"),
                  fmt(pr$observed_M4_PVM3labels), fmt(pr$null_median), fmt(pr$null_q025), fmt(pr$null_q975),
                  fmt_p(pr$p_two), fmt_p(pr$fdr), reading_I, fmt(pr$M4_own_labels)))
    gr13b("GR-13b-17", if (s_ok) "PASS" else "FAIL",
          sprintf(paste("Second reading (NeoC, r of Module post-IP%% with %%MT): M4 cells r %s against random subsets",
                        "median %s [%s to %s], two-sided p %s, FDR %s. %s. Reference without a test: own labels r %s"),
                  fmt(sr$observed_M4_PVM3labels), fmt(sr$null_median), fmt(sr$null_q025), fmt(sr$null_q975),
                  fmt_p(sr$p_two), fmt_p(sr$fdr), reading_S, fmt(sr$M4_own_labels)))

    # == STEP 12: outputs =====================================================
    cat("\nSTEP 12 -- outputs\n")
    DR <- do.call(rbind, lapply(DATASETS, function(ds) data.frame(dataset = ds, draw = seq_len(B_NULL),
                                                                  valid = valid[[ds]], NULLS[[ds]], check.names = FALSE)))
    utils::write.csv(DR, OUT_CSV, row.names = FALSE)
    FRAME <- do.call(rbind, lapply(DATASETS, function(ds) {
        P <- PREP[[ds]]; d <- P$d
        data.frame(donor_id = d$.donor, dataset = ds, n_total = d$n_total, m4_cells = d$k, pvm3_pool_cells = d$k3,
                   m4_share = d$share, module_post_m4_own_labels = d$a_nat, module_post_pvm3_labels_on_m4 =
                       round(P$obs$x * d$n_total / 100), x_m4_own_labels = P$nat$x, x_m4_pvm3_labels = P$obs$x,
                   frac_m4_own_labels = P$nat$frac, frac_m4_pvm3_labels = P$obs$frac, mt_m4 = d$mt,
                   pvm3_module_post_pct = PF3$pct_post[match(paste(ds, d$.donor), paste(PF3$dataset, PF3$.donor))],
                   stringsAsFactors = FALSE)
    }))
    DESIGN <- data.frame(item = c("Pool", "Labels", "Strata", "Draws", "Seed", "Label rule", "Label file columns",
                                  "Models", "Empirical p"),
                         value = c("PVM_3 arm cells (MT-CO3 > 0), M4 cells included",
                                   "PVM_3 arm Module and IP_y for every cell; the M4 arm's own labels as a reference",
                                   if (use_region) "donor x region" else "donor",
                                   sprintf("%d per region set", B_NULL), sprintf("%d (L'Ecuyer-CMRG, one stream per draw; results do not depend on the number of cores)", NULL_SEED),
                                   RULE_USED,
                                   sprintf("dataset %s; donor %s; arm %s (M4 %s, PVM_3 %s); barcode %s; %%MT %s; region %s",
                                           C_DS, C_DN, C_ARM, V_M4, V_P3, C_BC, C_MT, if (is.na(C_RG)) "none" else C_RG),
                                   paste(sprintf("%s: %s", names(MODELS), vapply(MODELS, `[[`, "", "label")), collapse = "; "),
                                   "p_upper = (1 + draws >= observed)/(B + 1); p_lower = (1 + draws <= observed)/(B + 1); two-sided = min(1, 2 x smaller); BH FDR within region set per coefficient over the three cognition outcomes, and over the three collinearity readouts"),
                         stringsAsFactors = FALSE)
    STRATA_TAB <- do.call(rbind, lapply(DATASETS, function(ds) {
        P <- PREP[[ds]]; first <- match(seq_len(P$S), P$p$six)
        data.frame(dataset = ds, donor_id = P$p$.dn[first], region = if (use_region) P$p$.rg[first] else "all",
                   pool_cells = P$size_s, m4_cells = P$k_s, stringsAsFactors = FALSE)
    }))
    READ <- data.frame(item = c("Primary reading (NeoC MMSE, MA)", "Second reading (NeoC, r Module post-IP% with %MT)",
                                "Label rule", "Pre-set text"),
                       value = c(reading_I, reading_S, RULE_USED, paste(PRESET_TEXT, collapse = "\n")),
                       stringsAsFactors = FALSE)
    INPUTS <- data.frame(item = names(inputs), path = unname(inputs), stringsAsFactors = FALSE)
    README <- data.frame(sheet = c("A_Design", "B_Frame", "C_Strata", "D_Summary", "E_Verdicts_MB", "F_Label_Agreement",
                                   "G_Label_Rules", "G2_Label_V1_Check", "H_Readings", "I_Guardrails", "J_Guardrail_Detail", "K_Inputs"),
        content = c("Pool, labels, strata, draws, seed, the label rule and columns used, the models and the empirical p.",
                    "Per donor and region set: Micro-PVM_3 cells, M4 cells, PVM_3 pool cells, M4 share, Module post-IP cells and Module post-IP% of the M4 cells under the M4 arm's own labels and under PVM_3 labels, the fraction, %MT, and the PVM_3 arm's own Module post-IP%.",
                    "Every stratum: pool cells and M4 cells drawn.",
                    "Every readout: the M4 cells under PVM_3 labels, the M4 cells under their own labels (reference, no test), the random-subset distribution (mean, SD, median, 2.5 and 97.5 percentiles), empirical p (upper, lower, two-sided) and BH FDR (bh_family: within region set, each coefficient over the three cognition outcomes; the three collinearity readouts as one family). bA = beta of Module post-IP% in MA; bB = beta of Module post-IP% in MB, with qB, qBm (within-draw FDR of Module post-IP% and of %MT) and vB (largest VIF) descriptive; bC = beta of the fraction in MC, vC its largest VIF; bD, bDm = beta of the fraction and of %MT in MD, vD the largest VIF; r_x_mt and r_frac_mt = donor-level Pearson r with %MT; mean_frac = mean fraction over donors with at least 5 M4 cells. Descriptive readouts carry no p.",
                    "MB verdicts (as Figure 14D) of the M4 cells and the share of random subsets with each verdict.",
                    "M4 cells: Module post-IP under the M4 arm's own labels against the PVM_3 labels.",
                    "Every candidate label rule and its largest per-donor deviation from the Part 10h Module post-IP cells on each arm.",
                    "Cell-level check of the Module post-IP label against Part 10h's Validation 1 cell table, per arm and region set.",
                    "Pre-set readings and the readings of this run.",
                    "Guardrail ledger (also written as a .txt).",
                    "Every value compared by a guardrail (also written as a .txt).",
                    "Input files. The draws are in the _null_draws.csv file."),
        stringsAsFactors = FALSE)
    if (is.null(V1TAB)) V1TAB <- data.frame(note = "V1_CHECK = FALSE")
    G <- do.call(rbind, .gr); GD <- if (length(.gd)) do.call(rbind, .gd) else data.frame()
    wb <- list(README = README, A_Design = DESIGN, B_Frame = FRAME, C_Strata = STRATA_TAB, D_Summary = SUM,
               E_Verdicts_MB = VER, F_Label_Agreement = AGREE, G_Label_Rules = CAND, G2_Label_V1_Check = V1TAB, H_Readings = READ,
               I_Guardrails = G, J_Guardrail_Detail = GD, K_Inputs = INPUTS)
    wb <- lapply(wb, function(x) { if (is.data.frame(x)) for (k in names(x))
        if (is.numeric(x[[k]])) x[[k]][!is.finite(x[[k]]) & !is.na(x[[k]])] <- NA; x })
    openxlsx::write.xlsx(wb, OUT_XLSX, overwrite = TRUE, firstRow = TRUE)
    ok_fig <- tryCatch({ draw_figure(SUM, NULLS, valid, SHORT, COG_LAB); TRUE },
                       error = function(e) { cat("  Figure failed: ", conditionMessage(e), "\n", sep = ""); FALSE })
    .lc <- file(OUT_LEG, open = "w", encoding = "UTF-8"); writeLines(LEGEND_TXT, .lc); close(.lc)
    dr_ok <- file.exists(OUT_CSV) && nrow(utils::read.csv(OUT_CSV, nrows = 5)) > 0
    gr13b("GR-13b-18", if (ok_fig && dr_ok && all(file.exists(c(OUT_XLSX, OUT_PDF, OUT_PNG, OUT_LEG)))) "PASS" else "FAIL",
          sprintf("outputs written: %s, %s, %s, %s, %s", basename(OUT_XLSX), basename(OUT_CSV), basename(OUT_PDF),
                  basename(OUT_PNG), basename(OUT_LEG)))
    for (f in c(OUT_XLSX, OUT_CSV, OUT_PDF, OUT_PNG, OUT_LEG)) if (file.exists(f)) cat(sprintf("  %s\n", basename(f)))
    paste0(reading_I, " | ", reading_S)
}

# == FIGURE ===================================================================
# Part 12 arm palette: M4 green, PVM_3 amber; own labels neutral dark.
COL13B <- list(m4 = "#1B9E77", pvm3 = "#F2A950", own = "#3B3B3B", strip = "#262626", grid = "#EBEBEB",
               hist_fill = "#F9D8A8", hist_border = "#B9762A")
HIST_BINS <- 45L
# geom_label border: `linewidth` from ggplot2 3.5.0, `label.size` before.
label_noborder <- if (utils::packageVersion("ggplot2") >= "3.5.0") list(linewidth = 0) else list(label.size = 0)
p13b_theme <- function() {
    theme_bw(base_size = 7) +
        theme(panel.grid.minor = element_blank(),
              panel.grid.major = element_line(colour = COL13B$grid, linewidth = 0.25),
              strip.background = element_rect(fill = COL13B$strip, colour = COL13B$strip),
              strip.text = element_text(colour = "white", face = "bold", size = 6.5),
              plot.title = element_text(face = "bold", size = 7.5, hjust = 0),
              axis.title = element_text(size = 6.5), axis.text = element_text(size = 6, colour = "black"),
              legend.position = "bottom", legend.text = element_text(size = 6), legend.key.size = unit(3, "mm"),
              legend.margin = margin(0, 0, 0, 0), plot.tag = element_text(face = "bold", size = 9))
}
draw_figure <- function(SUM, NULLS, valid, SHORT, COG_LAB) {
    th <- p13b_theme()
    LINES <- c("M4 cells, PVM_3 labels" = COL13B$m4, "M4 cells, own labels" = COL13B$own)
    panel <- function(ros, labeller, title, xlab, flip_cogs = TRUE) {
        h <- do.call(rbind, lapply(names(NULLS), function(ds) do.call(rbind, lapply(ros, function(ro) {
            v <- NULLS[[ds]][valid[[ds]], ro]
            fl <- if (flip_cogs && grepl("_CogS$", ro)) -1 else 1
            data.frame(dataset = ds, readout = ro, value = fl * v[is.finite(v)], stringsAsFactors = FALSE)
        }))))
        s <- SUM[SUM$readout %in% ros, ]
        fl <- ifelse(flip_cogs & grepl("_CogS$", s$readout), -1, 1)
        ln <- rbind(data.frame(dataset = s$dataset, readout = s$readout, value = fl * s$observed_M4_PVM3labels,
                               series = "M4 cells, PVM_3 labels"),
                    data.frame(dataset = s$dataset, readout = s$readout, value = fl * s$M4_own_labels,
                               series = "M4 cells, own labels"))
        ln$series <- factor(ln$series, names(LINES))
        tx <- data.frame(dataset = s$dataset, readout = s$readout, txt = sprintf("FDR = %s", fmt_p(s$fdr)))
        lv <- unname(labeller[ros])
        flev <- paste0(rep(lv, each = length(NULLS)), ", ", rep(names(NULLS), length(lv)))
        addf <- function(z) { z$facet <- factor(paste0(labeller[z$readout], ", ", z$dataset), flev); z }
        h <- addf(h); ln <- addf(ln); tx <- addf(tx)
        ggplot(h, aes(x = value)) +
            geom_histogram(bins = HIST_BINS, fill = COL13B$hist_fill, colour = COL13B$hist_border, linewidth = 0.2) +
            geom_vline(data = ln, aes(xintercept = value, colour = series, linetype = series), linewidth = 0.55) +
            do.call(geom_label, c(list(data = tx, mapping = aes(label = txt), x = -Inf, y = Inf, hjust = -0.05, vjust = 1.15,
                                       size = 1.9, fill = "white", label.padding = unit(0.6, "mm"), inherit.aes = FALSE),
                                  label_noborder)) +
            scale_colour_manual(values = LINES, name = NULL, drop = FALSE) +
            scale_linetype_manual(values = c("M4 cells, PVM_3 labels" = "solid", "M4 cells, own labels" = "22"),
                                  name = NULL, drop = FALSE) +
            facet_wrap(~ facet, ncol = 2, scales = "free") +
            labs(x = xlab, y = "Random subsets (draws)", title = title) + th
    }
    cog <- function(prefix) stats::setNames(paste0(ifelse(COG_LAB == "Cognitive status", "CogS\u2020", SHORT[COG_LAB])),
                                            paste0(prefix, SHORT[COG_LAB]))
    pA <- panel(names(cog("bA_")), cog("bA_"), "Module post-IP% on cognition (MA, as Figure 9A)",
                "Standardized \u03b2 of Module post-IP%")
    pB <- panel(names(cog("bC_")), cog("bC_"), "Fraction at a fixed M4 share (MC, as Part 13 S08)",
                "Standardized \u03b2 of the fraction")
    labC <- c(r_x_mt = "r, Module post-IP% with %MT", r_frac_mt = "r, fraction with %MT", mean_frac = "Mean fraction (%)")
    pC <- panel(names(labC), labC, "Collinearity with %MT and the fraction level", "Value", flip_cogs = FALSE)
    fig <- pA / pB / pC + plot_layout(heights = c(1, 1, 1), guides = "collect") +
        plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")") &
        theme(legend.position = "bottom")
    dev_pdf <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
    ggsave(OUT_PDF, fig, width = 180, height = 240, units = "mm", device = dev_pdf)
    ggsave(OUT_PNG, fig, width = 180, height = 240, units = "mm", dpi = 300, bg = "white")
}

LEGEND_TXT <- c(
    paste0("Supplementary Figure [n+2]. Are the M4 results a property of the M4 cells? Each donor's M4 cells were compared with random PVM_3 cells of the same number, drawn without replacement within each donor and region from the PVM_3 arm cells (Micro-PVM_3 cells with MT-CO3 >0, the M4 cells included), so the M4 share is the same in every draw (", format(B_NULL, big.mark = ","), " draws per region set). Every cell, M4 or random, carries its PVM_3 arm labels (PVM_3 Module and PVM_3 IP_y), so the M4 cells and the random cells are defined the same way. Module post-IP%, Module post-IP cells in the set as a percentage of the donor's Micro-PVM_3 cells; fraction, Module post-IP cells as a percentage of the cells in the set (donors with at least 5 M4 cells); %MT, log1p of the mean per-cell %MT over the set, as in Figure 14D. Models at cov4 (age at death, sex, PMI and APOE \u03b54 dose), outcome and predictors standardized after the complete-case filter. Histograms (light amber), the random subsets; solid green line, the M4 cells under PVM_3 labels; dashed dark line, the M4 cells under the M4 arm's own labels (its own Module and IP_y), shown for reference only, because the random subsets carry PVM_3 labels. FDR, BH-corrected two-sided empirical p of the M4 cells under PVM_3 labels against the random subsets; the empirical p is min(1, 2 x the smaller of (1 + draws at or above) / (B + 1) and (1 + draws at or below) / (B + 1)); BH within region set over the three cognition outcomes for each coefficient (panels A and B) and over the three panel C readouts. The random subsets are matched on cell number per donor and region only, so a difference can reflect any property that separates Mic.13-like cells from other PVM_3 cells, including sequencing depth."),
    "(A) Standardized \u03b2 of Module post-IP% on cognition, as Figure 9A. \u2020Cognitive status is sign-reversed so that lower \u03b2 means worse cognition throughout.",
    "(B) Standardized \u03b2 of the fraction with the M4 share in the model, as Part 13 S08; the M4 share is identical in every draw, so only the fraction varies. Cognitive status sign-reversed as in panel A.",
    "(C) Donor-level Pearson r of Module post-IP% with %MT (donors with M4 cells) and of the fraction with %MT (donors with at least 5 M4 cells), and the mean fraction over donors with at least 5 M4 cells. Values, the Figure 14D model (Module post-IP% and %MT in one model) and its verdicts for the random subsets are in the supplementary table.")

# == RUN ======================================================================
.reading <- tryCatch(
    withCallingHandlers(run_part13b(), warning = function(w) {
        cat("  WARNING: ", conditionMessage(w), "\n", sep = ""); invokeRestart("muffleWarning")
    }),
    error = function(e) {
        cat("\nRun stopped: ", conditionMessage(e), "\n", sep = "")
        structure(conditionMessage(e), class = "p13b_error")
    })
.note <- if (inherits(.reading, "p13b_error")) paste("STOPPED -", .reading) else paste("completed; readings:", .reading)
write_ledger(.note)
cat("\nGuardrail ledger : ", LEDGER, "\nGuardrail detail : ", GDETAIL, "\nConsole log      : ", CONSOLE,
    "\nRun: ", .note, "\n\n", sep = "")
cat(utils::capture.output(utils::sessionInfo()), sep = "\n")
options(.old_width)
sink(); close(.con)
