# =============================================================================
# Part13_CP10K_M4Share_code.r
# SEA-AD CEP-IP/CMRF, Part 13, M4 share (Supplementary Figure S20, Table S18)
#   Tests whether the donor's M4 share (M4 cells / Micro-PVM_3 cells, the
#   denominator of Module post-IP%) explains the M4 Module post-IP% results.
#   Adds the share to the Figure 9, 10 and 14D models (M4 arm, NeoC and
#   NonNeoC), with a donor bootstrap of the change in the Module coefficient,
#   and a post-hoc model of the within-M4 fraction against %MT.
#
# Usage (from the SEA-AD folder, or set SEAAD_ROOT), after Part 10h, the
# Part 10h comparator and Part13a_CP10K_Find14DpctMT_code.r
#   source("Part13_CP10K_M4Share_code.r")
#
# Outputs (Part13_M4Share/): supplementary table (xlsx), two figures (PDF,
#   PNG), pre-set and post-hoc readings, legend draft, guardrail log and
#   console log.
# =============================================================================

# == 0. Configuration ========================================================
if (!exists("CEP_HOME")) CEP_HOME <- Sys.getenv("SEAAD_ROOT", getwd())
OUT_DIR      <- file.path(CEP_HOME, "Part13_M4Share")
FILE_PREFIX  <- "Part13_CP10K_M4Share"

COMP_RDS      <- NA   # NA = search CEP_HOME for COMP_BASENAME
COMP_BASENAME <- "Part10h_CP10K_M4_Comparator_results.rds"
P10H_RDS      <- NA   # NA = the path in the comparator's meta, else search
P10H_BASENAME <- "Part10h_CP10K_M4_Prognosis_results.rds"
CVP_RDS       <- NA   # NA = search; absent gives a WARN
CVP_BASENAME  <- "Part10h_CP10K_M4_CogVsPath_results.rds"
SEARCH_EXCLUDE <- c("Previous_version", "Archive", "archive")

# Figure 14D %MT from Part13a_CP10K_Find14DpctMT_code.r: log1p of the donor mean of
# per-cell %MT over the M4 cells. "MitoFrac_num" uses the donor_frame column instead.
MT_FILE      <- file.path(OUT_DIR, "Part13_M4_pctMT_for14D.csv")
MT_SHEET     <- NA    # for an xlsx; NA = first sheet
MT_COLS      <- list(donor = "donor_id", dataset = "dataset", value = "pct_mt")
MT_ARM_COL   <- NA    # set when the file holds several arms, e.g. "arm"
MT_ARM_VALUE <- c("M4", "Mic13_PVM3")
MT_TRANSFORM <- "none"  # the file already holds log1p values; do not transform again

ARM_TAG_EXPECT    <- "Mic13_PVM3"
DATASETS          <- c("NeoC", "NonNeoC")
N_DONORS_EXPECT   <- c(NeoC = 84L, NonNeoC = 83L)
N_M4_CELLS_EXPECT <- c(NeoC = 2130L, NonNeoC = 1987L)   # Methods S3, Table 3
MIN_M4_CELLS_FRAC <- 5L  # within-M4 fraction only for donors with at least this many M4 cells

# Optional sensitivity: M4 share from all M4-selected Micro-PVM_3 cells (MT-CO3 = 0
# included). CSV, one row per donor and region set; absent donors count as 0.
SEL_FILE <- NA
SEL_COLS <- list(donor = "donor_id", dataset = "dataset", n = "n_m4_selected")
N_SEL_CELLS_EXPECT <- c(NeoC = 2626L, NonNeoC = 2329L)  # Results: M4-selected PVM_3 cells

MIN_DONORS   <- 20L      # comparator MIN_DONORS
VIF_MAX      <- 5
ALPHA        <- 0.05
B_BOOT       <- 10000L
BOOT_SEED    <- 20261010L
BOOT_MIN_VALID <- 0.95   # share of replicates that must be full rank
GUARD_STRICT <- TRUE
ALLOW_PRESET_REWRITE <- FALSE
TOL_EXACT    <- 1e-8     # identities and unrounded statistics
TOL_ROUND4   <- 6e-5     # statistics stored at 4 decimals
TOL_ROUND2   <- 6e-3     # statistics stored at 2 decimals
TOL_PUB      <- 1.5e-3   # values published at 3 decimals

# == 1. Libraries, console log, ledger =======================================
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
OUT_PDF  <- file.path(OUT_DIR, paste0(FILE_PREFIX, "_SuppFig.pdf"))
OUT_PNG  <- file.path(OUT_DIR, paste0(FILE_PREFIX, "_SuppFig.png"))
OUT_LEG  <- file.path(OUT_DIR, paste0(FILE_PREFIX, "_legend_draft.txt"))
OUT_PDF2 <- file.path(OUT_DIR, paste0(FILE_PREFIX, "_SuppFig2.pdf"))
OUT_PNG2 <- file.path(OUT_DIR, paste0(FILE_PREFIX, "_SuppFig2.png"))
PH_PRESET <- file.path(OUT_DIR, paste0(FILE_PREFIX, "_posthoc_readings.txt"))

.gr <- list(); .gd <- list()
gr13 <- function(id, status, detail) {
    if (length(detail) != 1L || is.na(detail)) detail <- "[detail unavailable]"
    .gr[[length(.gr) + 1L]] <<- data.frame(id = id, status = status, detail = detail,
                                           stringsAsFactors = FALSE)
    cat(sprintf("  [%-4s]  %-9s %s\n", status, id, detail))
    if (status == "FAIL" && GUARD_STRICT)
        stop(sprintf("%s FAILED: %s", id, detail), call. = FALSE)
}
# Element-wise comparison; every value is kept for the detail file.
gcmp <- function(id, item, mine, ref, tol) {
    mine <- as.numeric(mine); ref <- as.numeric(ref)
    ok <- !is.na(mine) & !is.na(ref) & abs(mine - ref) <= tol
    .gd[[length(.gd) + 1L]] <<- data.frame(id = id, item = item, mine = mine, ref = ref,
                                           abs_diff = abs(mine - ref), tol = tol, pass = ok,
                                           stringsAsFactors = FALSE)
    ok
}
write_ledger <- function(final_note) {
    G <- if (length(.gr)) do.call(rbind, .gr) else
        data.frame(id = character(0), status = character(0), detail = character(0))
    G <- G[order(G$id), , drop = FALSE]   # run order is in the console log
    lines <- c(strrep("=", 72),
               "Part 13 - M4 share - Guardrail Log [CP10K]",
               paste("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
               paste("Root:", CEP_HOME), paste("Write:", OUT_DIR),
               strrep("=", 72), "",
               if (nrow(G)) sprintf("  [%-4s]  %-9s %s", G$status, G$id, G$detail) else "  (none)",
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

# Remove outputs of an earlier run first.
.stale <- c(OUT_XLSX, OUT_PDF, OUT_PNG, OUT_PDF2, OUT_PNG2, OUT_LEG, GDETAIL, LEDGER)
invisible(file.remove(.stale[file.exists(.stale)]))
.con <- file(CONSOLE, open = "wt")
sink(.con, split = TRUE)
cat(strrep("=", 78), "\n  Part 13 - M4 share\n  ", format(Sys.time()), "\n  ",
    R.version.string, "\n  Root : ", CEP_HOME, "\n  Write: ", OUT_DIR, "\n",
    strrep("=", 78), "\n", sep = "")

# == 2. Helpers ==============================================================
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
    if (any(drop)) cat(sprintf("  %s: %d copy(ies) under excluded folders ignored: %s\n",
                               base, sum(drop), paste(hit[drop], collapse = "; ")))
    hit <- hit[!drop]
    if (length(hit) == 1L) return(normalizePath(hit, winslash = "/"))
    if (!length(hit) && !required) return(NA_character_)
    stop(sprintf("%s: %d match(es) under %s%s. Set the path in Section 0.", base, length(hit),
                 CEP_HOME, if (length(hit)) paste0(": ", paste(hit, collapse = "; ")) else ""),
         call. = FALSE)
}

# One OLS fit: outcome and predictors standardised after the complete-case
# filter, covariates unscaled, the outcome removed from its own covariates.
fit_ols <- function(d, y, preds, covs) {
    covs <- setdiff(covs, y)
    vars <- unique(c(y, preds, covs))
    D <- as.data.frame(lapply(d[, vars, drop = FALSE], num))
    names(D) <- vars
    keep <- stats::complete.cases(D)
    D <- D[keep, , drop = FALSE]
    n <- nrow(D)
    if (n < MIN_DONORS) return(NULL)
    D[[y]] <- .z(D[[y]]); for (p in preds) D[[p]] <- .z(D[[p]])
    if (anyNA(D[, c(y, preds)])) return(NULL)
    rhs <- c(preds, covs)
    r2 <- function(terms, resp = y) {
        if (!length(terms)) return(0)
        summary(stats::lm(stats::as.formula(paste(bq(resp), "~", paste(bq(terms), collapse = " + "))),
                          data = D))$r.squared
    }
    fit <- stats::lm(stats::as.formula(paste(bq(y), "~", paste(bq(rhs), collapse = " + "))), data = D)
    cf <- summary(fit)$coefficients
    rn <- gsub("`", "", rownames(cf))
    ci <- suppressMessages(stats::confint(fit)); rownames(ci) <- gsub("`", "", rownames(ci))
    vf <- .safe_vif(fit); names(vf) <- gsub("`", "", names(vf))
    # A singular predictor correlation matrix gives NA (read as unreadable).
    vp <- if (length(preds) > 1L) tryCatch(diag(solve(stats::cor(D[, preds]))),
                                           error = function(e) stats::setNames(rep(NA_real_, length(preds)), preds))
          else stats::setNames(1, preds)
    R2 <- summary(fit)$r.squared; R2c <- r2(covs)
    n_alias <- sum(is.na(stats::coef(fit)))
    rows <- lapply(preds, function(p) {
        i <- match(p, rn)
        data.frame(predictor = p, n = n,
                   beta = if (is.na(i)) NA_real_ else cf[i, 1],
                   se = if (is.na(i)) NA_real_ else cf[i, 2],
                   p_value = if (is.na(i)) NA_real_ else cf[i, 4],
                   ci_lo = if (is.na(i)) NA_real_ else ci[p, 1],
                   ci_hi = if (is.na(i)) NA_real_ else ci[p, 2],
                   dR2_unique = R2 - r2(setdiff(rhs, p)), dR2_block = R2 - R2c,
                   vif_full = unname(vf[p]), vif_pred = unname(vp[p]),
                   n_aliased = n_alias, stringsAsFactors = FALSE)
    })
    list(rows = do.call(rbind, rows), donors = d$.donor[keep])
}

# == 3. Pre-set readings =====================================================
PRESET_TEXT <- c(
    "Part13_CP10K_M4Share pre-set readings",
    "Fixed before the first run. A changed text stops the run unless ALLOW_PRESET_REWRITE = TRUE.",
    "",
    "Definitions",
    "- M4 share: the donor's M4 arm cells (M4-selected Micro-PVM_3 cells with MT-CO3 > 0, every CMRF class) as a percentage of the donor's Micro-PVM_3 cells, per region set = pct_pre + pct_post + pct_excl_pre + pct_excl_post of Part 10h prop_df (the four Figure 9 populations). Donors without M4 cells score 0, as in Figure 9.",
    "- Within-M4 fraction: 100 x pct_post / M4 share, the Module post-IP cells as a percentage of the donor's M4 cells; computed only for donors with at least 5 M4 cells in the region set, otherwise excluded.",
    "- With Module post-IP% and the M4 share in one model, the Module coefficient is the association of a larger Module post-IP fraction of the M4 cells at a fixed M4 share (Module post-IP% = M4 share x within-M4 fraction).",
    "- Holds: the Module post-IP% coefficient reaches FDR < 0.05 with the expected sign (negative for MMSE and CASI, positive for cognitive status) and the full-model VIF of every predictor is < 5.",
    "- BH families: each model with the M4 share uses the family of the figure it extends (Fig 10B: the 15 outcomes; Fig 9B: the cognition outcomes per block; Fig 14D: the 9 clinical outcomes; Fig 10C/10D: MMSE and CASI per block and model). Bootstrap p-values: BH within region set over the nine cognition tests.",
    "",
    "Primary test: M4 NeoC MMSE, in three models",
    "  m1 = cov4 + Module post-IP% + M4 share (S03, Fig 10B family)",
    "  m2 = cov4 + Braak + Thal + CERAD + Module post-IP% + M4 share (S07, Fig 9B family)",
    "  m3 = cov4 + Module post-IP% + %MT + M4 share (S10, Fig 14D family)",
    "",
    "Readings (first that applies)",
    "R0 Not separable: a predictor VIF >= 5 in m1, m2 or m3. The question is then read on m1b = cov4 + M4 share + within-M4 fraction (S08): 'state beyond abundance' if the fraction holds (FDR < 0.05, negative, VIF < 5); 'abundance only' if only the share reaches FDR < 0.05; otherwise 'inconclusive'.",
    "R1 Not explained by the M4 share: Module post-IP% holds in m1, m2 and m3, and in m1 either the donor-bootstrap 95% CI of the change in its coefficient includes 0 or the attenuation is 25% or less.",
    "R2 Partly explained by the M4 share: Module post-IP% holds in at least one of m1 to m3 and R1 is not met.",
    "R3 Explained by the M4 share: Module post-IP% holds in none of m1 to m3 and the M4 share reaches FDR < 0.05 in at least one.",
    "R4 Inconclusive: neither Module post-IP% nor the M4 share reaches FDR < 0.05 in m1 to m3.",
    "",
    "Sensitivity without a reading: when SEL_FILE is supplied, m1 to m3 are refitted with the M4 share built from all M4-selected Micro-PVM_3 cells (MT-CO3 = 0 included) in place of the MT-CO3+ share (S03s, S07s, S10s).",
    "Reported without a reading: every other outcome, NonNeoC, the Figure 14D verdicts that change, the decomposition models (Fig 10C/10D) and m1b when R0 does not apply."
)

# Post-hoc readings for S15-S17, fixed before these models are first fitted.
PH_TEXT <- c(
    "Part13_CP10K_M4Share post-hoc readings: Module state against %MT within the same M4 cells",
    "Added 2026-10-10 after the first Part 13 run had been read (pre-set reading R0; m1b 'state beyond abundance'). This test is POST-HOC. The text is fixed before its models are first fitted; a changed text stops the run unless ALLOW_PRESET_REWRITE = TRUE.",
    "",
    "Question: within the same M4 cells, does the Module state carry cognition beyond the cells' %MT?",
    "- Both predictors come from the donor's M4 cells in the region set: the within-M4 fraction (Module post-IP cells as a percentage of M4 cells) and %MT (log1p of the donor mean of per-cell %MT over the M4 cells, as Figure 14D). Donors with at least 5 M4 cells.",
    "- p1 = cov4 + within-M4 fraction + %MT (S15; the 9 clinical outcomes; BH within region set and predictor, as Figure 14D).",
    "- p2 = cov4 + Braak + Thal + CERAD + within-M4 fraction + %MT (S16; cognition; BH within region set and predictor over the cognition outcomes, as Figure 9B).",
    "- p3 = p1 + M4 share (S17; sensitivity, no reading).",
    "- Holds: FDR < 0.05 with the expected sign (negative for MMSE and CASI, positive for cognitive status) and every full-model VIF < 5.",
    "",
    "Primary test: NeoC MMSE in p1. Readings (first that applies)",
    "P0 Not separable: a predictor VIF >= 5 in p1.",
    "P1 Module state beyond %MT: the fraction holds in p1 and %MT does not reach FDR < 0.05.",
    "P2 Both: the fraction holds in p1 and %MT reaches FDR < 0.05.",
    "P3 %MT beyond the Module state: %MT reaches FDR < 0.05 in p1 and the fraction does not hold.",
    "P4 Neither: neither reaches FDR < 0.05 in p1.",
    "",
    "Reported without a reading: p2, p3, NonNeoC and the other outcomes."
)

# ============================================================================
run_part13 <- function() {
    # == STEP 1: pre-set readings ============================================
    cat("\nSTEP 1 -- pre-set readings\n")
    if (file.exists(PRESET) && !ALLOW_PRESET_REWRITE) {
        same <- identical(readLines(PRESET, warn = FALSE), PRESET_TEXT)
        gr13("GR-13-01", if (same) "PASS" else "FAIL",
             sprintf("pre-set readings %s the file written by an earlier run (%s)",
                     if (same) "match" else "differ from", basename(PRESET)))
    } else {
        writeLines(PRESET_TEXT, PRESET)
        gr13("GR-13-01", "PASS", sprintf("pre-set readings written before any model: %s",
                                         basename(PRESET)))
    }

    # == STEP 2: inputs ======================================================
    cat("\nSTEP 2 -- inputs\n")
    comp_path <- find_file(COMP_RDS, COMP_BASENAME)
    RC <- readRDS(comp_path)
    if (is.null(RC$donor_frame) || is.null(RC$meta))
        stop("The comparator RDS carries no donor_frame or meta.", call. = FALSE)
    p10h_path <- P10H_RDS
    if (is.na(p10h_path) && !is.null(RC$meta$part10h_rds) && file.exists(RC$meta$part10h_rds))
        p10h_path <- RC$meta$part10h_rds
    p10h_path <- find_file(p10h_path, P10H_BASENAME)
    R10 <- readRDS(p10h_path)
    cvp_path <- find_file(CVP_RDS, CVP_BASENAME, required = FALSE)
    RV <- if (!is.na(cvp_path)) readRDS(cvp_path) else NULL
    MT_FROM_FRAME <- identical(MT_FILE, "MitoFrac_num")
    if (!MT_FROM_FRAME && (is.na(MT_FILE) || !file.exists(MT_FILE)))
        stop("MT_FILE is not set or not found. It must hold the Figure 14D donor %MT of the M4 cells.",
             call. = FALSE)
    inputs <- c(comparator = comp_path, part10h = p10h_path, cogvspath = cvp_path,
                mt = if (MT_FROM_FRAME) NA_character_ else MT_FILE)
    for (k in names(inputs)) if (!is.na(inputs[[k]]))
        cat(sprintf("  %-10s %s\n", k, inputs[[k]]))
    arm_ok <- identical(RC$meta$arm_tag, ARM_TAG_EXPECT)
    p10_same <- is.null(RC$meta$part10h_rds) ||
        identical(basename(RC$meta$part10h_rds), basename(p10h_path))
    gr13("GR-13-02", if (arm_ok && p10_same) "PASS" else "FAIL",
         sprintf(paste("comparator arm_tag '%s' (expected '%s'), schema '%s'; it was built from '%s' and",
                       "this run reads '%s'; CogVsPath RDS %s"),
                 RC$meta$arm_tag %||% "unrecorded", ARM_TAG_EXPECT, RC$meta$schema %||% "unrecorded",
                 basename(RC$meta$part10h_rds %||% "unrecorded"), basename(p10h_path),
                 if (is.null(RV)) "not found" else basename(cvp_path)))

    # == STEP 3: donor frame =================================================
    cat("\nSTEP 3 -- donor frame\n")
    DF <- as.data.frame(RC$donor_frame)
    need <- c(".donor", "dataset", "module", "pct_post", "pct_pre", "pct_excl_post", "pct_excl_pre",
              "pct_module_all", "pct_excl", "pct_any_post", "n_total", "B1", "B2", "B4", "B5")
    miss <- setdiff(need, names(DF))
    if (length(miss)) stop("donor_frame lacks: ", paste(miss, collapse = ", "), call. = FALSE)
    DF$.donor <- as.character(DF$.donor); DF$dataset <- as.character(DF$dataset)
    OUT_MAP <- unlist(RC$meta$outcomes)
    CONTROL <- RC$meta$control_outcomes
    COV4 <- unlist(RC$meta$covariate_columns$cov4)
    CLIN_LAB <- c("Cognitive status", "MMSE score", "CASI score", "CERAD score", "Braak stage",
                  "Thal phase", "ADNC rank", "CAA score", "APOE e4 dose")
    COG_LAB  <- CLIN_LAB[1:3]
    BTC_LAB  <- c("Braak stage", "Thal phase", "CERAD score")
    lost <- setdiff(CLIN_LAB, names(OUT_MAP))
    if (length(lost)) stop("Outcomes absent from the comparator's meta$outcomes: ",
                           paste(lost, collapse = ", "), call. = FALSE)
    if (length(COV4) != 4L || !all(COV4 %in% names(DF)))
        stop("cov4 does not resolve to four donor_frame columns: ", paste(COV4, collapse = ", "),
             call. = FALSE)
    BTC <- unname(OUT_MAP[BTC_LAB])
    OUT15 <- names(OUT_MAP)
    SHORT <- c("Cognitive status" = "CogS", "MMSE score" = "MMSE", "CASI score" = "CASI",
               "CERAD score" = "CERAD", "Braak stage" = "Braak", "Thal phase" = "Thal",
               "ADNC rank" = "ADNC", "CAA score" = "CAA", "APOE e4 dose" = "APOE")
    EXP_SIGN <- c("Cognitive status" = 1, "MMSE score" = -1, "CASI score" = -1)
    GROUP <- function(lab) ifelse(lab %in% COG_LAB, "Cognition",
                           ifelse(lab == "APOE e4 dose", "APOE",
                           ifelse(lab %in% CONTROL, "Control", "Pathology")))
    ds_ok <- identical(sort(unique(DF$dataset)), sort(DATASETS))
    n_ds <- as.integer(table(factor(DF$dataset, DATASETS)))
    dup <- any(duplicated(DF[, c(".donor", "dataset")]))
    gr13("GR-13-03", if (ds_ok && !dup && all(n_ds == N_DONORS_EXPECT)) "PASS" else "FAIL",
         sprintf("datasets %s; donors %s (expected %s); duplicated donor rows: %s; outcomes %d (%d controls); cov4 = %s",
                 paste(unique(DF$dataset), collapse = "/"), paste(n_ds, collapse = "/"),
                 paste(N_DONORS_EXPECT, collapse = "/"), dup, length(OUT15), length(CONTROL),
                 paste(COV4, collapse = " + ")))

    # == STEP 4: formula consistency =========================================
    cat("\nSTEP 4 -- formula consistency (numerators and denominators)\n")
    for (k in c("module", "pct_post", "pct_pre", "pct_excl_post", "pct_excl_pre", "pct_module_all",
                "pct_excl", "pct_any_post", "n_total", "B1")) DF[[k]] <- num(DF[[k]])
    d1 <- max(abs(DF$module - DF$pct_post))
    gr13("GR-13-04", if (isTRUE(d1 <= TOL_EXACT)) "PASS" else "FAIL",
         sprintf("comparator module == Part 10h pct_post on every donor: max |diff| %.2e", d1))
    d2 <- max(abs(DF$pct_any_post - (DF$pct_post + DF$pct_excl_post)))
    d3 <- max(abs(DF$pct_module_all - (DF$pct_pre + DF$pct_post)))
    gr13("GR-13-05", if (isTRUE(max(d2, d3) <= TOL_EXACT)) "PASS" else "FAIL",
         sprintf(paste("Part 10h identities on one denominator: pct_any_post = pct_post + pct_excl_post",
                       "(max |diff| %.2e); pct_module_all = pct_pre + pct_post (%.2e)"), d2, d3))
    DF$share <- DF$pct_pre + DF$pct_post + DF$pct_excl_pre + DF$pct_excl_post
    d4 <- max(abs(DF$share - (DF$pct_module_all + DF$pct_excl)))
    n_m4_raw <- DF$share * DF$n_total / 100
    d5 <- max(abs(n_m4_raw - round(n_m4_raw)))
    bad <- sum(!(DF$n_total >= 1) | DF$share > 100 + 1e-9 | DF$pct_post > DF$share + 1e-9)
    gr13("GR-13-06", if (isTRUE(d4 <= TOL_EXACT && d5 <= 1e-6 && bad == 0)) "PASS" else "FAIL",
         sprintf(paste("M4 share = the four Figure 9 populations = pct_module_all + pct_excl (max |diff| %.2e);",
                       "share x n_total / 100 is a whole cell count on every donor (max deviation %.2e);",
                       "donors with n_total < 1, share > 100 or Module > share: %d"), d4, d5, bad))
    DF$n_m4 <- round(n_m4_raw)
    tot <- as.integer(tapply(DF$n_m4, factor(DF$dataset, DATASETS), sum))
    gr13("GR-13-07", if (all(tot == N_M4_CELLS_EXPECT)) "PASS" else "FAIL",
         sprintf("M4 arm cells rebuilt from the shares: %s (expected %s, Methods S3 Table 3)",
                 paste(tot, collapse = "/"), paste(N_M4_CELLS_EXPECT, collapse = "/")))
    if ("n_micro" %in% names(DF)) {
        d6 <- max(abs(DF$B1 * num(DF$n_micro) / 100 - DF$n_total), na.rm = TRUE)
        gr13("GR-13-08", if (isTRUE(d6 <= 1e-6)) "PASS" else "FAIL",
             sprintf(paste("denominator n_total = the donor's Micro-PVM_3 cells: CP (B1) x microglial cells",
                           "/ 100 reproduces n_total to %.2e"), d6))
    } else {
        gr13("GR-13-08", "WARN", "n_micro absent from donor_frame; the denominator was not cross-checked against CP")
    }
    DF$frac <- ifelse(DF$n_m4 >= MIN_M4_CELLS_FRAC, 100 * DF$pct_post / DF$share, NA_real_)
    cat(sprintf("  within-M4 fraction: %d donor-region-set rows with 1 to %d M4 cells excluded, %d with 0 M4 cells\n",
                sum(DF$n_m4 > 0 & DF$n_m4 < MIN_M4_CELLS_FRAC), MIN_M4_CELLS_FRAC - 1L, sum(DF$n_m4 == 0)))
    # S02r: the M4 share on the S08 donors only (at least MIN_M4_CELLS_FRAC M4 cells).
    DF$share_r <- ifelse(is.finite(DF$frac), DF$share, NA_real_)
    HAS_SEL <- !is.na(SEL_FILE)
    if (HAS_SEL) {
        if (!file.exists(SEL_FILE)) stop("SEL_FILE not found: ", SEL_FILE, call. = FALSE)
        sel <- utils::read.csv(SEL_FILE, check.names = FALSE, stringsAsFactors = FALSE)
        miss <- setdiff(unlist(SEL_COLS), names(sel))
        if (length(miss)) stop("SEL_FILE lacks: ", paste(miss, collapse = ", "), call. = FALSE)
        sel <- data.frame(.donor = as.character(sel[[SEL_COLS$donor]]), dataset = as.character(sel[[SEL_COLS$dataset]]),
                          n_sel = num(sel[[SEL_COLS$n]]), stringsAsFactors = FALSE)
        if (any(duplicated(sel[, c(".donor", "dataset")]))) stop("SEL_FILE has duplicated donor rows.", call. = FALSE)
        unk <- sum(!paste(sel$.donor, sel$dataset) %in% paste(DF$.donor, DF$dataset))
        DF$n_sel <- sel$n_sel[match(paste(DF$.donor, DF$dataset), paste(sel$.donor, sel$dataset))]
        DF$n_sel[is.na(DF$n_sel)] <- 0
        DF$share_sel <- 100 * DF$n_sel / DF$n_total
        tot_sel <- as.integer(tapply(DF$n_sel, factor(DF$dataset, DATASETS), sum))
        bs1 <- sum(DF$n_sel < DF$n_m4); bs2 <- sum(DF$n_sel > DF$n_total)
        gr13("GR-13-21", if (bs1 == 0 && bs2 == 0 && unk == 0 && all(tot_sel == N_SEL_CELLS_EXPECT)) "PASS" else "FAIL",
             sprintf(paste("all-M4-selected share: totals %s (expected %s); donors with fewer selected cells than",
                           "MT-CO3+ M4 cells %d; with more selected cells than Micro-PVM_3 cells %d; file rows",
                           "matching no donor %d; r with the MT-CO3+ share %s"),
                     paste(tot_sel, collapse = "/"), paste(N_SEL_CELLS_EXPECT, collapse = "/"), bs1, bs2, unk,
                     fmt(stats::cor(DF$share_sel, DF$share))))
    } else {
        gr13("GR-13-21", "SKIP", "SEL_FILE not supplied; the all-M4-selected share sensitivity was not run")
    }

    if (MT_FROM_FRAME) {
        if (!"MitoFrac_num" %in% names(RC$donor_frame)) stop("donor_frame has no MitoFrac_num.", call. = FALSE)
        mt_raw <- data.frame(RC$donor_frame$.donor, RC$donor_frame$dataset, RC$donor_frame$MitoFrac_num)
        names(mt_raw) <- unlist(MT_COLS)
    } else mt_raw <- if (grepl("\\.xlsx?$", MT_FILE, ignore.case = TRUE))
        openxlsx::read.xlsx(MT_FILE, sheet = if (is.na(MT_SHEET)) 1 else MT_SHEET) else
        utils::read.csv(MT_FILE, check.names = FALSE, stringsAsFactors = FALSE)
    if (!is.na(MT_ARM_COL)) mt_raw <- mt_raw[mt_raw[[MT_ARM_COL]] %in% MT_ARM_VALUE, , drop = FALSE]
    miss <- setdiff(unlist(MT_COLS), names(mt_raw))
    if (length(miss)) stop("MT_FILE lacks: ", paste(miss, collapse = ", "), "; its columns are ",
                           paste(names(mt_raw), collapse = ", "), call. = FALSE)
    MT <- data.frame(.donor = as.character(mt_raw[[MT_COLS$donor]]),
                     dataset = as.character(mt_raw[[MT_COLS$dataset]]),
                     mt = num(mt_raw[[MT_COLS$value]]), stringsAsFactors = FALSE)
    if (MT_TRANSFORM == "log1p") MT$mt <- log1p(MT$mt)
    if (any(duplicated(MT[, c(".donor", "dataset")])))
        stop("MT_FILE has duplicated donor rows; set MT_ARM_COL.", call. = FALSE)
    DF <- merge(DF, MT, by = c(".donor", "dataset"), all.x = TRUE, sort = FALSE)
    DF <- DF[order(DF$dataset, DF$.donor), ]
    mt_zero <- sum(!is.na(DF$mt) & DF$n_m4 == 0)
    mt_miss <- sum(is.na(DF$mt) & DF$n_m4 > 0)
    r_mito <- if ("MitoFrac_num" %in% names(DF))
        suppressWarnings(stats::cor(DF$mt, num(DF$MitoFrac_num), use = "complete.obs")) else NA_real_
    gr13("GR-13-09", if (mt_zero == 0 && !isTRUE(r_mito > 0.999)) "PASS" else "WARN",
         sprintf(paste("%%MT of the M4 cells: %d donor(s) without M4 cells carry a value and %d with M4 cells",
                       "lack one; r with the Micro-PVM_3-wide MitoFrac_num %s"),
                 mt_zero, mt_miss, fmt(r_mito)))

    # == STEP 4b: post-hoc readings, fixed before the post-hoc models are fitted ====
    if (file.exists(PH_PRESET) && !ALLOW_PRESET_REWRITE) {
        same <- identical(readLines(PH_PRESET, warn = FALSE), PH_TEXT)
        gr13("GR-13-22", if (same) "PASS" else "FAIL",
             sprintf("post-hoc readings %s the file written by an earlier run (%s)",
                     if (same) "match" else "differ from", basename(PH_PRESET)))
    } else {
        writeLines(PH_TEXT, PH_PRESET)
        gr13("GR-13-22", "PASS", sprintf("post-hoc readings written before the post-hoc models: %s",
                                         basename(PH_PRESET)))
    }

    # == STEP 5: models ======================================================
    cat("\nSTEP 5 -- models\n")
    SPECS <- list(
        list(id = "S01", label = "Module alone",                fig = "9A/10A", preds = "module",                   block = "cov4", outs = OUT15),
        list(id = "S02", label = "M4 share alone",              fig = "new",    preds = "share",                    block = "cov4", outs = OUT15),
        list(id = "S03", label = "Module + M4 share",           fig = "10B+",   preds = c("module", "share"),       block = "cov4", outs = OUT15),
        list(id = "S04", label = "Module alone",                fig = "9B",     preds = "module",                   block = "cov4", outs = COG_LAB),
        list(id = "S05", label = "Module alone",                fig = "9B",     preds = "module",                   block = "+BTC", outs = COG_LAB),
        list(id = "S06", label = "Module + M4 share",           fig = "9B+",    preds = c("module", "share"),       block = "cov4", outs = COG_LAB),
        list(id = "S07", label = "Module + M4 share",           fig = "9B+",    preds = c("module", "share"),       block = "+BTC", outs = COG_LAB),
        list(id = "S08", label = "M4 share + within-M4 fraction", fig = "new",  preds = c("share", "frac"),         block = "cov4", outs = OUT15),
        list(id = "S02r", label = "M4 share alone, S08 donors", fig = "new",   preds = "share_r",                  block = "cov4", outs = OUT15),
        list(id = "S09", label = "Module + %MT",                fig = "14D",    preds = c("module", "mt"),          block = "cov4", outs = CLIN_LAB),
        list(id = "S10", label = "Module + %MT + M4 share",     fig = "14D+",   preds = c("module", "mt", "share"), block = "cov4", outs = CLIN_LAB),
        list(id = "S11", label = "Module + CP",                 fig = "10B",    preds = c("module", "B1"),          block = "cov4", outs = OUT15),
        list(id = "S12", label = "Module + AUCell",             fig = "10B",    preds = c("module", "B2"),          block = "cov4", outs = OUT15),
        list(id = "S13", label = "Module + UCell",              fig = "10B",    preds = c("module", "B4"),          block = "cov4", outs = OUT15),
        list(id = "S14", label = "Module + AMS",                fig = "10B",    preds = c("module", "B5"),          block = "cov4", outs = OUT15))
    DEC <- list(
        list(id = "D01", label = "CP + AUCell (Fig 10C)",       fig = "10C",  preds = c("module", "B1", "B2"),          score = "B2"),
        list(id = "D02", label = "M4 share + AUCell",           fig = "10C+", preds = c("module", "share", "B2"),       score = "B2"),
        list(id = "D03", label = "CP + M4 share + AUCell",      fig = "10C+", preds = c("module", "B1", "share", "B2"), score = "B2"),
        list(id = "D04", label = "CP + UCell (Fig 10D)",        fig = "10D",  preds = c("module", "B1", "B4"),          score = "B4"),
        list(id = "D05", label = "CP + AMS (Fig 10D)",          fig = "10D",  preds = c("module", "B1", "B5"),          score = "B5"),
        list(id = "D06", label = "M4 share + UCell",            fig = "10D+", preds = c("module", "share", "B4"),       score = "B4"),
        list(id = "D07", label = "M4 share + AMS",              fig = "10D+", preds = c("module", "share", "B5"),       score = "B5"),
        list(id = "D08", label = "AUCell only (Fig 10D label)", fig = "10D",  preds = c("module", "B2"),                score = "B2"),
        list(id = "D09", label = "UCell only (Fig 10D label)",  fig = "10D",  preds = c("module", "B4"),                score = "B4"),
        list(id = "D10", label = "AMS only (Fig 10D label)",    fig = "10D",  preds = c("module", "B5"),                score = "B5"))
    SPECS <- c(SPECS, list(   # POST-HOC, see PH_TEXT
        list(id = "S15", label = "Within-M4 fraction + %MT (post-hoc)", fig = "14D like-for-like", preds = c("frac", "mt"), block = "cov4", outs = CLIN_LAB),
        list(id = "S16", label = "Within-M4 fraction + %MT (post-hoc)", fig = "9B like-for-like",  preds = c("frac", "mt"), block = "+BTC", outs = COG_LAB),
        list(id = "S17", label = "Within-M4 fraction + %MT + M4 share (post-hoc)", fig = "14D like-for-like sens", preds = c("frac", "mt", "share"), block = "cov4", outs = CLIN_LAB)))
    if (HAS_SEL) SPECS <- c(SPECS, list(
        list(id = "S03s", label = "Module + all-selected M4 share", fig = "10B+ sens", preds = c("module", "share_sel"), block = "cov4", outs = OUT15),
        list(id = "S07s", label = "Module + all-selected M4 share", fig = "9B+ sens",  preds = c("module", "share_sel"), block = "+BTC", outs = COG_LAB),
        list(id = "S10s", label = "Module + %MT + all-selected M4 share", fig = "14D+ sens", preds = c("module", "mt", "share_sel"), block = "cov4", outs = CLIN_LAB)))
    for (b in c("cov4", "+BTC")) for (s in DEC) {
        s$id <- paste0(s$id, if (b == "cov4") "a" else "b"); s$block <- b
        s$outs <- c("MMSE score", "CASI score"); SPECS[[length(SPECS) + 1L]] <- s
    }
    covs_of <- function(block) if (block == "+BTC") c(COV4, BTC) else COV4

    RES <- list(); n_try <- 0L; n_fit <- 0L
    for (s in SPECS) for (ds in DATASETS) {
        d <- DF[DF$dataset == ds, , drop = FALSE]
        for (lab in s$outs) {
            n_try <- n_try + 1L
            r <- fit_ols(d, OUT_MAP[[lab]], s$preds, covs_of(s$block))
            if (is.null(r)) next
            n_fit <- n_fit + 1L
            RES[[length(RES) + 1L]] <- cbind(
                data.frame(spec = s$id, model = s$label, extends = s$fig, block = s$block,
                           dataset = ds, outcome = lab, group = GROUP(lab),
                           score = s$score %||% NA_character_, stringsAsFactors = FALSE), r$rows)
        }
    }
    RES <- do.call(rbind, RES)
    RES$fdr <- NA_real_
    key <- paste(RES$spec, RES$dataset, RES$predictor, sep = "|")
    for (k in unique(key)) {
        i <- which(key == k & is.finite(RES$p_value))
        if (length(i)) RES$fdr[i] <- stats::p.adjust(RES$p_value[i], method = "BH")
    }
    cat(sprintf("  %d of %d model(s) fitted (%d below %d donors)\n", n_fit, n_try, n_try - n_fit, MIN_DONORS))
    gv <- function(spec, ds, lab, pred, what = "beta") {
        i <- which(RES$spec == spec & RES$dataset == ds & RES$outcome == lab & RES$predictor == pred)
        if (length(i) != 1L) NA_real_ else RES[[what]][i]
    }
    vmax_of <- function(spec, ds, lab) {
        v <- RES$vif_full[RES$spec == spec & RES$dataset == ds & RES$outcome == lab]
        # No fitted model: NA. A predictor with NA VIF: Inf (read as unreadable).
        if (!length(v)) NA_real_ else if (anyNA(v)) Inf else max(v)
    }

    # == STEP 6: reproduction guardrails =====================================
    cat("\nSTEP 6 -- reproduction of the stored and published estimates\n")
    ok <- logical(0)
    for (ds in DATASETS) {
        mv <- R10$datasets[[ds]]$mv_df
        if (is.null(mv)) next
        mv <- mv[mv$tier == "all" & mv$covariate_set == "cov4" & mv$predictor == "pct_post", , drop = FALSE]
        for (j in seq_len(nrow(mv))) {
            lab <- mv$outcome[j]; if (!lab %in% OUT15 || !is.finite(mv$beta_std[j])) next
            it <- paste(ds, lab)
            ok <- c(ok, gcmp("GR-13-10", paste(it, "n"), gv("S01", ds, lab, "module", "n"), mv$n_donors[j], 0),
                    gcmp("GR-13-10", paste(it, "beta"), gv("S01", ds, lab, "module"), mv$beta_std[j], 1e-6),
                    gcmp("GR-13-10", paste(it, "dR2"), gv("S01", ds, lab, "module", "dR2_block"), mv$delta_R2[j], 1e-6),
                    gcmp("GR-13-10", paste(it, "FDR"), gv("S01", ds, lab, "module", "fdr"), mv$p_adj[j], 1e-6))
        }
    }
    gr13("GR-13-10", if (length(ok) && all(ok)) "PASS" else "FAIL",
         sprintf("Part 10h Multivariate (cov4, pct_post, 15 outcomes; Figure 9A) reproduced on %d of %d values",
                 sum(ok), length(ok)))
    ok <- logical(0); MG <- RC$marginal; DU <- RC$duels
    if (is.data.frame(MG)) {
        mg <- MG[MG$block == "cov4" & MG$key == "module" & MG$outcome %in% OUT15, , drop = FALSE]
        for (j in seq_len(nrow(mg))) {
            ds <- mg$dataset[j]; lab <- mg$outcome[j]; it <- paste(ds, lab, "marginal module")
            ok <- c(ok, gcmp("GR-13-11", paste(it, "n"), gv("S01", ds, lab, "module", "n"), mg$n[j], 0),
                    gcmp("GR-13-11", paste(it, "beta"), gv("S01", ds, lab, "module"), mg$beta[j], TOL_ROUND4),
                    gcmp("GR-13-11", paste(it, "dR2"), gv("S01", ds, lab, "module", "dR2_block"), mg$delta_r2[j], TOL_ROUND4),
                    gcmp("GR-13-11", paste(it, "FDR"), gv("S01", ds, lab, "module", "fdr"), mg$p_adj[j], TOL_EXACT))
        }
    }
    vmap <- c(module_vs_B1 = "S11", module_vs_B2 = "S12", module_vs_B4 = "S13", module_vs_B5 = "S14")
    if (is.data.frame(DU)) {
        du <- DU[DU$duel %in% names(vmap) & DU$outcome %in% OUT15, , drop = FALSE]
        for (j in seq_len(nrow(du))) {
            sp <- vmap[[du$duel[j]]]; rv <- sub("module_vs_", "", du$duel[j])
            ds <- du$dataset[j]; lab <- du$outcome[j]; it <- paste(ds, lab, du$duel[j])
            ok <- c(ok, gcmp("GR-13-11", paste(it, "beta module"), gv(sp, ds, lab, "module"), du$beta_module[j], TOL_ROUND4),
                    gcmp("GR-13-11", paste(it, "beta rival"), gv(sp, ds, lab, rv), du$beta_rival[j], TOL_ROUND4),
                    gcmp("GR-13-11", paste(it, "FDR module"), gv(sp, ds, lab, "module", "fdr"), du$p_module_adj[j], TOL_EXACT),
                    gcmp("GR-13-11", paste(it, "FDR rival"), gv(sp, ds, lab, rv, "fdr"), du$p_rival_adj[j], TOL_EXACT),
                    gcmp("GR-13-11", paste(it, "VIF"), vmax_of(sp, ds, lab), du$vif_joint[j], TOL_ROUND2))
        }
    }
    gr13("GR-13-11", if (length(ok) && all(ok)) "PASS" else "FAIL",
         sprintf(paste("comparator Marginal (module, cov4) and Duels vs CP, AUCell, UCell, AMS (Figures 10A-10B)",
                       "reproduced on %d of %d values"), sum(ok), length(ok)))
    if (!is.null(RV) && is.data.frame(RV$models)) {
        ok <- logical(0)
        mo <- RV$models[RV$models$key == "module" & RV$models$block %in% c("cov4", "+BTC") &
                        RV$models$outcome %in% COG_LAB, , drop = FALSE]
        fam <- unique(RV$models$outcome[RV$models$key == "module" & RV$models$block == "cov4"])
        for (j in seq_len(nrow(mo))) {
            sp <- if (mo$block[j] == "cov4") "S04" else "S05"
            ds <- mo$dataset[j]; lab <- mo$outcome[j]; it <- paste(ds, lab, mo$block[j])
            ok <- c(ok, gcmp("GR-13-12", paste(it, "n"), gv(sp, ds, lab, "module", "n"), mo$n[j], 0),
                    gcmp("GR-13-12", paste(it, "beta"), gv(sp, ds, lab, "module"), mo$beta[j], TOL_ROUND4),
                    gcmp("GR-13-12", paste(it, "FDR"), gv(sp, ds, lab, "module", "fdr"), mo$p_adj[j], TOL_EXACT))
        }
        gr13("GR-13-12", if (length(ok) && all(ok)) "PASS" else "FAIL",
             sprintf(paste("CogVsPath Models (module, cov4 and +BTC; Figure 9B) reproduced on %d of %d values;",
                           "its BH family holds: %s"), sum(ok), length(ok), paste(fam, collapse = ", ")))
    } else {
        gr13("GR-13-12", "WARN", "CogVsPath RDS not found; the Figure 9B estimates and family were not cross-checked")
    }
    pub10d <- list(
        NeoC = list(B2 = c(-.427, -.073, -.362, -.163, -.420, -.025, -.320, -.098),
                    B4 = c(-.423, -.077, -.350, -.179, -.418, -.027, -.312, -.109),
                    B5 = c(-.389, -.128, -.312, -.241, -.383, -.083, -.272, -.178)),
        NonNeoC = list(B2 = c(-.212, -.270, -.177, -.256, -.213, -.196, -.191, -.137),
                       B4 = c(-.201, -.286, -.163, -.280, -.204, -.211, -.181, -.157),
                       B5 = c(-.213, -.297, -.134, -.380, -.212, -.224, -.144, -.262)))
    pub10d_fdr <- list(NeoC = list(B2 = c(.006, .015, .0096, .039), B4 = c(.008, .022, .012, .048),
                                   B5 = c(.011, .028, .017, .070)),
                       NonNeoC = list(B2 = c(.052, .052, .086, .086), B4 = c(.065, .065, .101, .101),
                                      B5 = c(.081, .081, .150, .150)))
    pub10c_fdr <- list(NeoC = c(.011, .030, .014, .053), NonNeoC = c(.254, .284, .225, .225))
    withcp <- c(B2 = "D01", B4 = "D04", B5 = "D05"); nocp <- c(B2 = "D08", B4 = "D09", B5 = "D10")
    tests4 <- list(c("a", "MMSE score"), c("a", "CASI score"), c("b", "MMSE score"), c("b", "CASI score"))
    ok <- logical(0)
    for (ds in DATASETS) for (sc in names(withcp)) for (k in 1:4) {
        sp <- paste0(withcp[[sc]], tests4[[k]][1]); sq <- paste0(nocp[[sc]], tests4[[k]][1]); lab <- tests4[[k]][2]
        it <- paste(ds, sc, tests4[[k]][1], lab)
        ok <- c(ok, gcmp("GR-13-13", paste(it, "10D beta module"), gv(sp, ds, lab, "module"), pub10d[[ds]][[sc]][2 * k - 1], TOL_PUB),
                gcmp("GR-13-13", paste(it, "10D beta score"), gv(sp, ds, lab, sc), pub10d[[ds]][[sc]][2 * k], TOL_PUB),
                gcmp("GR-13-13", paste(it, "10D label FDR"), gv(sq, ds, lab, "module", "fdr"), pub10d_fdr[[ds]][[sc]][k], TOL_PUB))
    }
    for (ds in DATASETS) for (k in 1:4)
        ok <- c(ok, gcmp("GR-13-13", paste(ds, "10C", tests4[[k]][1], tests4[[k]][2], "FDR"),
                         gv(paste0("D01", tests4[[k]][1]), ds, tests4[[k]][2], "module", "fdr"), pub10c_fdr[[ds]][k], TOL_PUB))
    pub10c_vif <- list(NeoC = c(3.36, 3.16, 3.57, 3.36), NonNeoC = c(2.32, 2.38, 2.33, 2.39))
    okv <- logical(0)
    for (ds in DATASETS) for (k in 1:4)
        okv <- c(okv, gcmp("GR-13-20", paste(ds, "10C", tests4[[k]][1], tests4[[k]][2], "largest VIF"),
                           vmax_of(paste0("D01", tests4[[k]][1]), ds, tests4[[k]][2]), pub10c_vif[[ds]][k], TOL_ROUND2))
    gr13("GR-13-20", if (all(okv)) "PASS" else "WARN",
         sprintf(paste("published Figure 10C largest VIF reproduced on %d of %d boxes"), sum(okv), length(okv)))
    gr13("GR-13-13", if (all(ok)) "PASS" else "FAIL",
         sprintf(paste("published Figure 10C/10D coefficients and FDRs (BH over MMSE and CASI per block and model)",
                       "reproduced on %d of %d values"), sum(ok), length(ok)))
    pub14d <- list(c("NeoC", "MMSE score", -0.529, -0.808, -0.251, 74, 0.003),
                   c("NonNeoC", "Cognitive status", 0.415, 0.161, 0.670, 80, 0.015),
                   c("NonNeoC", "MMSE score", -0.337, -0.575, -0.099, 76, 0.028))
    ok <- logical(0)
    for (v in pub14d) {
        ds <- v[1]; lab <- v[2]; x <- as.numeric(v[3:7]); it <- paste(ds, lab, "14D")
        ok <- c(ok, gcmp("GR-13-14", paste(it, "beta"), gv("S09", ds, lab, "module"), x[1], TOL_PUB),
                gcmp("GR-13-14", paste(it, "CI low"), gv("S09", ds, lab, "module", "ci_lo"), x[2], TOL_PUB),
                gcmp("GR-13-14", paste(it, "CI high"), gv("S09", ds, lab, "module", "ci_hi"), x[3], TOL_PUB),
                gcmp("GR-13-14", paste(it, "n"), gv("S09", ds, lab, "module", "n"), x[4], 0),
                gcmp("GR-13-14", paste(it, "FDR"), gv("S09", ds, lab, "module", "fdr"), x[5], TOL_PUB))
    }
    gr13("GR-13-14", if (all(ok)) "PASS" else "FAIL",
         sprintf(paste("published Figure 14D estimates (beta, 95%% CI, n, FDR over the 9 clinical outcomes)",
                       "reproduced on %d of %d values"),
                 sum(ok), length(ok)))

    # == STEP 7: verdicts and hold status ====================================
    cat("\nSTEP 7 -- verdicts\n")
    pair_tab <- function(spec, second, main = "module") {
        a <- RES[RES$spec == spec & RES$predictor == main, ]
        b <- RES[RES$spec == spec & RES$predictor == second, c("dataset", "outcome", "beta", "fdr")]
        names(b)[3:4] <- c("beta_second", "fdr_second")
        w <- merge(a, b, by = c("dataset", "outcome"), sort = FALSE)
        w$vif_max <- mapply(vmax_of, spec, w$dataset, w$outcome)
        w$verdict <- ifelse(!is.finite(w$vif_max) | w$vif_max >= VIF_MAX | is.na(w$fdr) | is.na(w$fdr_second),
                            "Unreadable",
                     ifelse(w$fdr < ALPHA & w$fdr_second < ALPHA, "Both",
                     ifelse(w$fdr < ALPHA, "Module only",
                     ifelse(w$fdr_second < ALPHA, "Second only", "Neither"))))
        w$main <- main; w$second <- second
        w
    }
    V <- rbind(pair_tab("S03", "share"), pair_tab("S07", "share"), pair_tab("S11", "B1"),
               pair_tab("S09", "mt"), pair_tab("S10", "mt"))
    if (HAS_SEL) V <- rbind(V, pair_tab("S03s", "share_sel"), pair_tab("S07s", "share_sel"), pair_tab("S10s", "mt"))
    V <- rbind(V, pair_tab("S15", "mt", "frac"), pair_tab("S16", "mt", "frac"), pair_tab("S17", "mt", "frac"))
    shv <- RES[RES$spec == "S10" & RES$predictor == "share", c("spec", "dataset", "outcome", "beta", "fdr")]
    names(shv)[4:5] <- c("beta_share", "fdr_share")
    V <- merge(V, shv, by = c("spec", "dataset", "outcome"), all.x = TRUE, sort = FALSE)
    if (is.data.frame(DU)) {
        a <- V[V$spec == "S11", c("dataset", "outcome", "verdict", "vif_max")]
        b <- DU[DU$duel == "module_vs_B1", c("dataset", "outcome", "verdict")]
        chk <- merge(a, b, by = c("dataset", "outcome"), suffixes = c("_mine", "_comp"))
        comp <- c("module only" = "Module only", "rival only" = "Second only", both = "Both",
                  neither = "Neither")[chk$verdict_comp]
        agree <- ifelse(chk$vif_max < VIF_MAX, chk$verdict_mine == comp, chk$verdict_mine == "Unreadable")
        agree[is.na(agree)] <- FALSE
        gr13("GR-13-15", if (length(agree) && all(agree)) "PASS" else "FAIL",
             sprintf(paste("Module vs CP verdicts (Figure 10B) agree with the comparator's Duels on %d of %d",
                           "outcomes"),
                     sum(agree), length(agree)))
    } else gr13("GR-13-15", "WARN", "comparator Duels absent; the verdict rule was not cross-checked")

    DD <- RES[grepl("^D", RES$spec), ]
    H <- do.call(rbind, lapply(split(DD, paste(DD$spec, DD$dataset, DD$outcome, sep = "|")), function(t) {
        m <- t[t$predictor == "module", ]; sc <- t[t$predictor == t$score[1], ]
        vmax <- suppressWarnings(max(t$vif_full))
        st <- if (!is.finite(vmax) || vmax >= VIF_MAX) "unreadable" else
            if (isTRUE(m$fdr < ALPHA && sign(m$beta) == EXP_SIGN[[m$outcome]])) "HOLDS" else "doesn't hold"
        data.frame(spec = m$spec, model = m$model, extends = m$extends, block = m$block,
                   dataset = m$dataset, outcome = m$outcome, n = m$n, beta_module = m$beta,
                   beta_score = if (nrow(sc)) sc$beta else NA_real_, fdr_module = m$fdr,
                   vif_max = vmax, vif_pred_max = max(t$vif_pred), status = st, stringsAsFactors = FALSE)
    }))
    rownames(H) <- NULL

    # == STEP 8: donor bootstrap =============================================
    cat("\nSTEP 8 -- donor bootstrap\n")
    boot_job <- function(ds, lab, base, block, seed) {
        y <- OUT_MAP[[lab]]; covs <- setdiff(covs_of(block), y)
        vars <- unique(c(y, base, "share", covs))
        d <- DF[DF$dataset == ds, , drop = FALSE]
        D <- as.data.frame(lapply(d[, vars, drop = FALSE], num)); names(D) <- vars
        D <- D[stats::complete.cases(D), , drop = FALSE]; n <- nrow(D)
        X0 <- cbind(`(Intercept)` = 1, as.matrix(D[, c(base, covs), drop = FALSE]))
        X1 <- cbind(`(Intercept)` = 1, as.matrix(D[, c(base, "share", covs), drop = FALSE]))
        one <- function(i) {
            yv <- .z(D[[y]][i]); x0 <- X0[i, , drop = FALSE]; x1 <- X1[i, , drop = FALSE]
            for (p in base) x0[, p] <- .z(x0[, p])
            for (p in c(base, "share")) x1[, p] <- .z(x1[, p])
            if (anyNA(yv) || anyNA(x1[, c(base, "share")])) return(c(NA, NA, 1))
            f0 <- stats::lm.fit(x0, yv); f1 <- stats::lm.fit(x1, yv)
            if (f0$rank < ncol(x0) || f1$rank < ncol(x1)) return(c(NA, NA, 1))
            c(f0$coefficients[["module"]], f1$coefficients[["module"]], 0)
        }
        obs <- one(seq_len(n))
        set.seed(seed)
        bb <- t(vapply(seq_len(B_BOOT), function(b) one(sample.int(n, n, replace = TRUE)), numeric(3)))
        okb <- bb[, 3] == 0
        dl <- bb[okb, 2] - bb[okb, 1]; at <- 100 * (1 - bb[okb, 2] / bb[okb, 1])
        data.frame(dataset = ds, outcome = lab, block = block,
                   set = if ("mt" %in% base) "cov4 + %MT" else block, n = n,
                   beta_without = obs[1], beta_with = obs[2], delta = obs[2] - obs[1],
                   delta_lo = unname(stats::quantile(dl, .025)), delta_hi = unname(stats::quantile(dl, .975)),
                   atten_pct = 100 * (1 - obs[2] / obs[1]),
                   atten_lo = unname(stats::quantile(at, .025)), atten_hi = unname(stats::quantile(at, .975)),
                   p_boot = max(2 * min(mean(dl <= 0), mean(dl >= 0)), 2 / B_BOOT),
                   B_valid = sum(okb), B_rank_skipped = sum(!okb), stringsAsFactors = FALSE)
    }
    jobs <- list()
    for (ds in DATASETS) for (lab in COG_LAB) {
        jobs[[length(jobs) + 1L]] <- list(ds, lab, "module", "cov4", "S04", "S06")
        jobs[[length(jobs) + 1L]] <- list(ds, lab, "module", "+BTC", "S05", "S07")
        jobs[[length(jobs) + 1L]] <- list(ds, lab, c("module", "mt"), "cov4", "S09", "S10")
    }
    t0 <- Sys.time()
    BT <- do.call(rbind, lapply(seq_along(jobs), function(j) {
        a <- jobs[[j]]; r <- boot_job(a[[1]], a[[2]], a[[3]], a[[4]], BOOT_SEED + j)
        r$spec_without <- a[[5]]; r$spec_with <- a[[6]]; r
    }))
    BT$fdr_boot <- stats::ave(BT$p_boot, BT$dataset, FUN = function(p) stats::p.adjust(p, "BH"))
    cat(sprintf("  %d job(s) x B = %d in %.1f s\n", length(jobs), B_BOOT,
                as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    ref0 <- mapply(gv, BT$spec_without, BT$dataset, BT$outcome, "module")
    ref1 <- mapply(gv, BT$spec_with, BT$dataset, BT$outcome, "module")
    d7 <- max(abs(c(BT$beta_without - ref0, BT$beta_with - ref1)))
    vmin <- min(BT$B_valid) / B_BOOT
    gr13("GR-13-16", if (isTRUE(d7 <= 1e-8) && vmin >= BOOT_MIN_VALID) "PASS" else "FAIL",
         sprintf(paste("bootstrap point estimates equal the fitted models (max |diff| %.2e); full-rank replicates",
                       "%s%% of %d at the lowest job (rank-deficient replicates skipped and counted)"),
                 d7, fmt(100 * vmin, 1), B_BOOT))

    # == STEP 9: readings ====================================================
    cat("\nSTEP 9 -- readings\n")
    holds <- function(spec, ds, lab) {
        b <- gv(spec, ds, lab, "module"); f <- gv(spec, ds, lab, "module", "fdr"); v <- vmax_of(spec, ds, lab)
        c(holds = isTRUE(f < ALPHA && sign(b) == EXP_SIGN[[lab]] && v < VIF_MAX), vif_ok = isTRUE(v < VIF_MAX))
    }
    mm <- c(m1 = "S03", m2 = "S07", m3 = "S10")
    hv <- sapply(mm, function(s) holds(s, "NeoC", "MMSE score"))
    shs <- sapply(mm, function(s) isTRUE(gv(s, "NeoC", "MMSE score", "share", "fdr") < ALPHA))
    b1 <- BT[BT$dataset == "NeoC" & BT$outcome == "MMSE score" & BT$set == "cov4", ]
    small <- isTRUE(b1$delta_lo <= 0 && b1$delta_hi >= 0) || isTRUE(b1$atten_pct <= 25)
    reading <- if (!all(hv["vif_ok", ])) {
        fh <- isTRUE(gv("S08", "NeoC", "MMSE score", "frac", "fdr") < ALPHA &&
                     gv("S08", "NeoC", "MMSE score", "frac") < 0 &&
                     vmax_of("S08", "NeoC", "MMSE score") < VIF_MAX)
        sh <- isTRUE(gv("S08", "NeoC", "MMSE score", "share", "fdr") < ALPHA)
        paste("R0 Not separable; m1b reading:",
              if (fh) "state beyond abundance" else if (sh) "abundance only" else "inconclusive")
    } else if (all(hv["holds", ]) && small) "R1 Not explained by the M4 share" else
        if (any(hv["holds", ])) "R2 Partly explained by the M4 share" else
        if (any(shs)) "R3 Explained by the M4 share" else "R4 Inconclusive"
    PRIM <- do.call(rbind, lapply(names(mm), function(k) {
        s <- mm[[k]]; L <- "MMSE score"
        data.frame(model = k, spec = s, n = gv(s, "NeoC", L, "module", "n"), beta_module = gv(s, "NeoC", L, "module"),
                   ci_lo = gv(s, "NeoC", L, "module", "ci_lo"), ci_hi = gv(s, "NeoC", L, "module", "ci_hi"),
                   fdr_module = gv(s, "NeoC", L, "module", "fdr"), beta_share = gv(s, "NeoC", L, "share"),
                   fdr_share = gv(s, "NeoC", L, "share", "fdr"), vif_max = vmax_of(s, "NeoC", L),
                   module_holds = unname(hv["holds", k]), stringsAsFactors = FALSE)
    }))
    print(transform(PRIM, beta_module = fmt(beta_module), ci_lo = fmt(ci_lo), ci_hi = fmt(ci_hi),
                    fdr_module = fmt_p(fdr_module), beta_share = fmt(beta_share), fdr_share = fmt_p(fdr_share),
                    vif_max = fmt(vif_max, 2)), row.names = FALSE)
    cat(sprintf("  m1 change in beta %s (95%% CI %s to %s), attenuation %s%%, bootstrap p %s\n",
                fmt(b1$delta), fmt(b1$delta_lo), fmt(b1$delta_hi), fmt(b1$atten_pct, 1), fmt_p(b1$p_boot)))
    gr13("GR-13-17", if (all(hv["vif_ok", ])) "PASS" else "WARN",
         sprintf("primary models readable (full-model VIF < 5): m1 %s, m2 %s, m3 %s. Pre-set reading: %s",
                 hv["vif_ok", 1], hv["vif_ok", 2], hv["vif_ok", 3], reading))
    CHG <- merge(V[V$spec == "S09", c("dataset", "outcome", "verdict")],
                 V[V$spec == "S10", c("dataset", "outcome", "verdict")],
                 by = c("dataset", "outcome"), suffixes = c("_14D", "_with_share"))
    CHG <- CHG[CHG$verdict_14D != CHG$verdict_with_share, , drop = FALSE]
    cat(sprintf("  Figure 14D verdicts that change with the M4 share: %d\n", nrow(CHG)))
    if (nrow(CHG)) print(CHG, row.names = FALSE)

    # == STEP 9b: post-hoc, the Module state against %MT within the M4 cells ====
    cat("\nSTEP 9b -- post-hoc: within-M4 fraction against %MT of the same M4 cells\n")
    fr <- DF[is.finite(DF$frac), , drop = FALSE]
    # Check frac from whole cell counts: frac = 100 x Module post-IP cells / M4 cells.
    n_post_raw <- fr$pct_post * fr$n_total / 100
    d8 <- if (nrow(fr)) max(abs(fr$frac - 100 * round(n_post_raw) / fr$n_m4), abs(n_post_raw - round(n_post_raw)))
          else NA_real_
    rng_bad <- sum(fr$frac < -1e-9 | fr$frac > 100 + 1e-9)
    nchk <- do.call(rbind, lapply(DATASETS, function(ds) do.call(rbind, lapply(CLIN_LAB, function(lab) {
        d <- DF[DF$dataset == ds, , drop = FALSE]
        cv <- setdiff(COV4, OUT_MAP[[lab]])
        k <- is.finite(d$frac) & is.finite(d$mt) & is.finite(num(d[[OUT_MAP[[lab]]]])) &
            stats::complete.cases(as.data.frame(lapply(d[, cv, drop = FALSE], num)))
        data.frame(dataset = ds, outcome = lab, n_expected = sum(k), n_fitted = gv("S15", ds, lab, "frac", "n"),
                   min_m4_cells = if (any(k)) min(d$n_m4[k]) else NA_real_,
                   n_14D = gv("S09", ds, lab, "module", "n"),
                   stringsAsFactors = FALSE)
    }))))
    # A model below MIN_DONORS is not fitted (n_fitted NA).
    n_match <- ifelse(is.na(nchk$n_fitted), nchk$n_expected < MIN_DONORS, nchk$n_expected == nchk$n_fitted)
    n_sub <- is.na(nchk$n_fitted) | (!is.na(nchk$n_14D) & nchk$n_fitted <= nchk$n_14D)
    # S02r must use exactly the S08 donors.
    r_same <- unlist(lapply(DATASETS, function(ds) lapply(OUT15, function(lab) {
        a1 <- gv("S02r", ds, lab, "share_r", "n"); a2 <- gv("S08", ds, lab, "share", "n")
        (is.na(a1) && is.na(a2)) || isTRUE(a1 == a2) })))
    cat("  M4 share alone on the S08 donors (S02r) against the share term in S08, cognition:\n")
    for (ds in DATASETS) for (lab in COG_LAB)
        cat(sprintf("    %-8s %-17s S02r beta %s FDR %s | S08 share beta %s FDR %s | n %s\n", ds, lab,
                    fmt(gv("S02r", ds, lab, "share_r")), fmt_p(gv("S02r", ds, lab, "share_r", "fdr")),
                    fmt(gv("S08", ds, lab, "share")), fmt_p(gv("S08", ds, lab, "share", "fdr")),
                    gv("S08", ds, lab, "share", "n")))
    frame_ok <- all(r_same) && isTRUE(d8 <= 1e-6) && rng_bad == 0 && all(n_match) &&
        all(nchk$min_m4_cells >= MIN_M4_CELLS_FRAC, na.rm = TRUE) && all(n_sub)
    gr13("GR-13-23", if (frame_ok) "PASS" else "FAIL",
         sprintf(paste("post-hoc frame: within-M4 fraction = 100 x Module post-IP cells / M4 cells, rebuilt from whole cell counts (max |diff| %s); fraction",
                       "outside 0-100 on %d rows; fitted n equals donors with >= %d M4 cells, %%MT, outcome and cov4",
                       "on %d of %d models (MMSE n %s, Figure 14D frame %s); no model has more donors than its Figure 14D model; S02r on the S08 donors in %d of %d models"),
                 formatC(d8, format = "e", digits = 2), rng_bad, MIN_M4_CELLS_FRAC,
                 sum(n_match), nrow(nchk),
                 paste(nchk$n_fitted[nchk$outcome == "MMSE score"], collapse = "/"),
                 paste(nchk$n_14D[nchk$outcome == "MMSE score"], collapse = "/"), sum(r_same), length(r_same)))
    PHT <- do.call(rbind, lapply(c("S15", "S16", "S17"), function(sp) do.call(rbind, lapply(DATASETS, function(ds)
        do.call(rbind, lapply(COG_LAB, function(lab) {
            if (is.na(gv(sp, ds, lab, "frac", "n"))) return(NULL)
            data.frame(spec = sp, dataset = ds, outcome = lab, n = gv(sp, ds, lab, "frac", "n"),
                       beta_frac = gv(sp, ds, lab, "frac"), ci_lo_frac = gv(sp, ds, lab, "frac", "ci_lo"),
                       ci_hi_frac = gv(sp, ds, lab, "frac", "ci_hi"), fdr_frac = gv(sp, ds, lab, "frac", "fdr"),
                       beta_mt = gv(sp, ds, lab, "mt"), ci_lo_mt = gv(sp, ds, lab, "mt", "ci_lo"),
                       ci_hi_mt = gv(sp, ds, lab, "mt", "ci_hi"), fdr_mt = gv(sp, ds, lab, "mt", "fdr"),
                       vif_max = vmax_of(sp, ds, lab),
                       verdict = V$verdict[V$spec == sp & V$dataset == ds & V$outcome == lab][1],
                       stringsAsFactors = FALSE)
        }))))))
    PHT$verdict <- sub("^Module only$", "Fraction only", sub("^Second only$", "%MT only", PHT$verdict))
    print(transform(PHT, beta_frac = fmt(beta_frac), fdr_frac = fmt_p(fdr_frac), beta_mt = fmt(beta_mt),
                    fdr_mt = fmt_p(fdr_mt), vif_max = fmt(vif_max, 2), ci_lo_frac = NULL, ci_hi_frac = NULL,
                    ci_lo_mt = NULL, ci_hi_mt = NULL), row.names = FALSE)
    pv <- vmax_of("S15", "NeoC", "MMSE score")
    pbf <- gv("S15", "NeoC", "MMSE score", "frac"); pff <- gv("S15", "NeoC", "MMSE score", "frac", "fdr")
    pfm <- gv("S15", "NeoC", "MMSE score", "mt", "fdr")
    p_vif_ok <- isTRUE(pv < VIF_MAX)
    f_sig <- isTRUE(pff < ALPHA); f_holds <- f_sig && isTRUE(sign(pbf) == EXP_SIGN[["MMSE score"]]) && p_vif_ok
    m_sig <- isTRUE(pfm < ALPHA)
    ph_reading <- if (!p_vif_ok) "P0 Not separable" else
        if (f_holds && !m_sig) "P1 Module state beyond %MT" else
        if (f_holds && m_sig) "P2 Both" else
        if (m_sig) "P3 %MT beyond the Module state" else
        if (f_sig) "P4 Neither (the fraction reaches FDR < 0.05 with the wrong sign)" else "P4 Neither"
    gr13("GR-13-24", if (p_vif_ok) "PASS" else "WARN",
         sprintf(paste("post-hoc primary model readable (NeoC MMSE, p1, full-model VIF %s): fraction beta %s, FDR %s;",
                       "%%MT FDR %s. Post-hoc reading: %s"),
                 fmt(pv, 2), fmt(pbf), fmt_p(pff), fmt_p(pfm), ph_reading))

    # == STEP 10: collinearity ===============================================
    pairs <- list(c("share", "module", "M4 share ~ Module post-IP%"), c("share", "mt", "M4 share ~ %MT"),
                  c("share", "B1", "M4 share ~ CP"), c("share", "B2", "M4 share ~ AUCell"),
                  c("share", "frac", "M4 share ~ within-M4 fraction"),
                  c("module", "mt", "Module post-IP% ~ %MT (Fig 14E)"),
                  c("frac", "mt", "Within-M4 fraction ~ %MT (post-hoc)"))
    if (HAS_SEL) pairs <- c(pairs, list(c("share_sel", "share", "All-selected share ~ M4 share"),
                                        c("share_sel", "module", "All-selected share ~ Module post-IP%")))
    COLL <- do.call(rbind, lapply(DATASETS, function(ds) do.call(rbind, lapply(pairs, function(pp) {
        d <- DF[DF$dataset == ds, ]; x <- num(d[[pp[1]]]); y <- num(d[[pp[2]]]); k <- is.finite(x) & is.finite(y)
        r <- stats::cor(x[k], y[k])
        data.frame(dataset = ds, pair = pp[3], n = sum(k), r = r,
                   r_s = stats::cor(x[k], y[k], method = "spearman"), vif_pair = 1 / (1 - r^2),
                   stringsAsFactors = FALSE)
    }))))
    print(transform(COLL, r = fmt(r), r_s = fmt(r_s), vif_pair = fmt(vif_pair, 2)), row.names = FALSE)

    # == STEP 11: workbook ===================================================
    cat("\nSTEP 11 -- supplementary table\n")
    SPEC_TAB <- do.call(rbind, lapply(SPECS, function(s) data.frame(
        spec = s$id, model = s$label, extends = s$fig, block = s$block,
        predictors = paste(s$preds, collapse = " + "), covariates = paste(covs_of(s$block), collapse = " + "),
        outcomes_in_BH_family = paste(s$outs, collapse = "; "), stringsAsFactors = FALSE)))
    FRAME <- DF[, intersect(c(".donor", "dataset", "n_total", "n_m4", "share", "share_r", "n_sel", "share_sel", "pct_post", "pct_pre",
                              "pct_excl_post", "pct_excl_pre", "frac", "mt", "B1", "B2", "B4", "B5",
                              unname(OUT_MAP), COV4), names(DF))]
    names(FRAME)[1] <- "donor_id"
    README <- data.frame(sheet = c("A_Specs", "B_Frame", "C_Models", "D_Verdicts", "E_Decomposition",
                                   "F_Bootstrap", "G_Collinearity", "H_Readings", "I_Guardrails",
                                   "J_Guardrail_Detail", "K_Inputs", "L_PostHoc"),
        content = c(
            "Every model: predictors, covariates, the outcome set that forms its BH family, and the figure it extends (+ = with the M4 share).",
            "M4 donor frame per region set: Micro-PVM_3 cells (n_total), M4 cells (n_m4), M4 share, the four Figure 9 populations, within-M4 fraction, %MT, CP and gene scores, outcomes, cov4.",
            "Every fitted model, one row per predictor: n, standardised beta, SE, p, 95% CI (t), FDR, unique dR2, block dR2, full-model and predictor-only VIF.",
            "Verdicts at FDR < 0.05: Module vs M4 share (S03, S07), Module vs CP (S11, Fig 10B), Module vs %MT without (S09, Fig 14D) and with the M4 share (S10); S03s, S07s, S10s use the all-selected share when SEL_FILE is supplied. Unreadable = full-model VIF >= 5 or not computable.",
            "Figure 10C/10D models with and without the M4 share: HOLDS, doesn't hold, unreadable.",
            "Donor bootstrap of the Module coefficient without and with the M4 share: change, attenuation, percentile 95% CIs, two-sided p (floor 2/B), FDR within region set; rank-deficient replicates counted.",
            "Pearson and Spearman r between the M4 share and the other predictors; VIF = 1/(1 - r^2).",
            "Pre-set reading, the three primary models and the Figure 14D verdicts that change.",
            "Guardrail ledger (also written as a .txt).",
            "Every value compared by a guardrail (also written as a .txt).",
            "Input files.",
            "POST-HOC (added after the first run): within-M4 fraction and %MT of the same M4 cells in one model, cov4 (S15), +BTC (S16) and with the M4 share (S17); cognition rows with both coefficients, VIF and verdict. All outcomes are in C_Models and D_Verdicts."),
        stringsAsFactors = FALSE)
    INPUTS <- data.frame(item = names(inputs), path = unname(inputs),
                         stringsAsFactors = FALSE)
    READ <- rbind(data.frame(item = "Pre-set reading", value = reading),
                  data.frame(item = paste0(PRIM$model, " (", PRIM$spec, ")"),
                             value = sprintf("n %d, beta %s [%s, %s], FDR %s; share beta %s, FDR %s; VIF %s; holds %s",
                                             PRIM$n, fmt(PRIM$beta_module), fmt(PRIM$ci_lo), fmt(PRIM$ci_hi),
                                             fmt_p(PRIM$fdr_module), fmt(PRIM$beta_share), fmt_p(PRIM$fdr_share),
                                             fmt(PRIM$vif_max, 2), PRIM$module_holds)),
                  data.frame(item = "M4 share on the S08 donors (no reading)",
                             value = paste(sprintf("%s %s: share alone (S02r) beta %s, FDR %s; share with the fraction (S08) beta %s, FDR %s",
                                                   rep(DATASETS, each = 3), rep(SHORT[COG_LAB], 2),
                                                   fmt(mapply(gv, "S02r", rep(DATASETS, each = 3), rep(COG_LAB, 2), "share_r")),
                                                   fmt_p(mapply(gv, "S02r", rep(DATASETS, each = 3), rep(COG_LAB, 2), "share_r", "fdr")),
                                                   fmt(mapply(gv, "S08", rep(DATASETS, each = 3), rep(COG_LAB, 2), "share")),
                                                   fmt_p(mapply(gv, "S08", rep(DATASETS, each = 3), rep(COG_LAB, 2), "share", "fdr"))),
                                           collapse = "; ")),
                  data.frame(item = "m1 bootstrap", value = sprintf("change %s [%s, %s], attenuation %s%%",
                             fmt(b1$delta), fmt(b1$delta_lo), fmt(b1$delta_hi), fmt(b1$atten_pct, 1))),
                  if (HAS_SEL) data.frame(item = "Sensitivity, all-selected share (no reading)",
                      value = paste(sprintf("%s: Module beta %s, FDR %s, VIF %s", c("S03s", "S07s", "S10s"),
                                            fmt(sapply(c("S03s", "S07s", "S10s"), function(x) gv(x, "NeoC", "MMSE score", "module"))),
                                            fmt_p(sapply(c("S03s", "S07s", "S10s"), function(x) gv(x, "NeoC", "MMSE score", "module", "fdr"))),
                                            fmt(sapply(c("S03s", "S07s", "S10s"), function(x) vmax_of(x, "NeoC", "MMSE score")), 2)),
                                    collapse = "; ")) else NULL,
                  data.frame(item = "Figure 14D verdicts changed", value = if (nrow(CHG))
                      paste(sprintf("%s %s: %s -> %s", CHG$dataset, CHG$outcome, CHG$verdict_14D,
                                    CHG$verdict_with_share), collapse = "; ") else "none"),
                  data.frame(item = "Pre-set text", value = paste(PRESET_TEXT, collapse = "\n")),
                  data.frame(item = "POST-HOC reading (fraction vs %MT, NeoC MMSE, p1)",
                             value = sprintf("%s; n %s, fraction beta %s, FDR %s; %%MT beta %s, FDR %s; VIF %s", ph_reading,
                                             gv("S15", "NeoC", "MMSE score", "frac", "n"), fmt(pbf), fmt_p(pff),
                                             fmt(gv("S15", "NeoC", "MMSE score", "mt")), fmt_p(pfm), fmt(pv, 2))),
                  data.frame(item = "Post-hoc text", value = paste(PH_TEXT, collapse = "\n")))
    n_rows <- sum(!duplicated(RES[, c("spec", "dataset", "outcome")]))
    gr13("GR-13-18", if (n_rows == n_fit) "PASS" else "FAIL",
         sprintf("every fitted model is in the workbook: %d of %d", n_rows, n_fit))
    G <- do.call(rbind, .gr); GD <- if (length(.gd)) do.call(rbind, .gd) else data.frame()
    wb <- list(README = README, A_Specs = SPEC_TAB, B_Frame = FRAME, C_Models = RES, D_Verdicts = V,
               E_Decomposition = H, F_Bootstrap = BT, G_Collinearity = COLL, H_Readings = READ,
               I_Guardrails = G, J_Guardrail_Detail = GD, K_Inputs = INPUTS, L_PostHoc = PHT)
    # Inf (a non-computable VIF) is written as an empty cell.
    wb <- lapply(wb, function(x) { if (is.data.frame(x)) for (k in names(x))
        if (is.numeric(x[[k]])) x[[k]][!is.finite(x[[k]]) & !is.na(x[[k]])] <- NA; x })
    openxlsx::write.xlsx(wb, OUT_XLSX, overwrite = TRUE, firstRow = TRUE)

    # == STEP 12: figure =====================================================
    cat("\nSTEP 12 -- supplementary figure\n")
    ok_fig <- tryCatch({ draw_figure(RES, V, H, BT, COLL, SHORT, COG_LAB); TRUE },
                       error = function(e) { cat("  Figure failed: ", conditionMessage(e), "\n", sep = ""); FALSE })
    ok_fig2 <- tryCatch({ draw_figure2(RES, V, DF, SHORT, COG_LAB, OUT_MAP[["MMSE score"]]); TRUE },
                        error = function(e) { cat("  Figure 2 failed: ", conditionMessage(e), "\n", sep = ""); FALSE })
    writeLines(c(LEGEND_TXT, "", LEGEND2_TXT), OUT_LEG)
    gr13("GR-13-25", if (ok_fig2 && file.exists(OUT_PDF2) && file.exists(OUT_PNG2)) "PASS" else "FAIL",
         sprintf("post-hoc figure written: %s, %s", basename(OUT_PDF2), basename(OUT_PNG2)))
    gr13("GR-13-19", if (ok_fig && file.exists(OUT_PDF) && file.exists(OUT_PNG)) "PASS" else "FAIL",
         sprintf("outputs written: %s, %s, %s, %s", basename(OUT_XLSX), basename(OUT_PDF), basename(OUT_PNG),
                 basename(OUT_LEG)))
    for (f in c(OUT_XLSX, OUT_PDF, OUT_PNG, OUT_PDF2, OUT_PNG2, OUT_LEG)) if (file.exists(f))
        cat(sprintf("  %s\n", basename(f)))
    paste0(reading, " | post-hoc: ", ph_reading)
}

# == Figure ==================================================================
COL <- list(original = "#3B3B3B", variant = "#6A3D9A", strip = "#262626", cg = "#24456F",
            vifline = "#C0392B", grid = "#EBEBEB")
VERDICT_FILL <- c("Module only" = "#2D2D2D", "Second only" = "#1275D8", "Both" = "#7B7B7B",
                  "Neither" = "#D9D9D9", "Unreadable" = "#FBE7A2")
HOLD_FILL <- c("HOLDS" = "#2D2D2D", "doesn't hold" = "#D9D9D9", "unreadable" = "#FBE7A2")

p13_theme <- function() {
    theme_bw(base_size = 7) +
        theme(panel.grid.minor = element_blank(),
              panel.grid.major = element_line(colour = COL$grid, linewidth = 0.25),
              strip.background = element_rect(fill = COL$strip, colour = COL$strip),
              strip.text = element_text(colour = "white", face = "bold", size = 7),
              plot.title = element_text(face = "bold", size = 7.5, hjust = 0),
              axis.title = element_text(size = 6.5), axis.text = element_text(size = 6, colour = "black"),
              legend.position = "bottom", legend.title = element_text(size = 6.5, face = "bold"),
              legend.text = element_text(size = 6), legend.key.size = unit(2.6, "mm"),
              legend.margin = margin(0, 0, 0, 0), plot.tag = element_text(face = "bold", size = 9))
}

# Verdict bars shared by both figures.
verdict_bars <- function(v, ylab, grp, title, second_label, ylevels, glevels, th, first_label = "Module only") {
    v$ylab <- ylab; v$grp <- grp
    cnt <- stats::aggregate(list(k = rep(1, nrow(v))), by = list(grp = v$grp, ylab = v$ylab, verdict = v$verdict), sum)
    cnt$verdict <- factor(cnt$verdict, names(VERDICT_FILL))
    lv <- levels(cnt$verdict)
    lv[lv == "Second only"] <- second_label; lv[lv == "Module only"] <- first_label
    levels(cnt$verdict) <- lv
    fills <- VERDICT_FILL
    names(fills)[names(fills) == "Second only"] <- second_label
    names(fills)[names(fills) == "Module only"] <- first_label
    cnt$ylab <- factor(cnt$ylab, rev(ylevels)); cnt$grp <- factor(cnt$grp, glevels)
    cnt$txt <- ifelse(as.character(cnt$verdict) %in% c(first_label, second_label, "Both"), "white", "black")
    ggplot(cnt, aes(x = k, y = ylab, fill = verdict, group = verdict)) +
        geom_col(position = position_stack(reverse = TRUE), colour = "black", linewidth = 0.25, width = 0.7,
                 show.legend = TRUE) +
        geom_text(aes(label = k, colour = txt), position = position_stack(vjust = 0.5, reverse = TRUE), size = 2.1) +
        scale_colour_identity() +
        scale_fill_manual(values = fills, drop = FALSE, name = "Verdict (FDR<0.05; Unreadable, VIF \u22655)") +
        scale_x_continuous(breaks = function(l) seq(0, floor(l[2]), 1), expand = expansion(mult = c(0, 0.02))) +
        facet_grid(. ~ grp, scales = "free_x") +
        labs(x = "Number of outcomes", y = NULL, title = title) + th +
        theme(panel.grid.major.y = element_blank(), legend.direction = "horizontal",
              strip.text = element_text(colour = "white", face = "bold", size = 6),
              panel.spacing.x = unit(2.5, "mm")) +
        guides(fill = guide_legend(nrow = 3, title.position = "top"))
}

draw_figure <- function(RES, V, H, BT, COLL, SHORT, COG_LAB) {
    th <- p13_theme()
    sig <- function(x) !is.na(x) & x < 0.05
    short_cog <- function(lab) ifelse(lab == "Cognitive status", "CogS\u2020", SHORT[lab])

    # (A) Module beta without and with the share
    rowdef <- expand.grid(set = c("cov4", "+BTC", "cov4 + %MT"), outcome = COG_LAB, stringsAsFactors = FALSE)
    rowdef$y <- seq_len(nrow(rowdef)) + (match(rowdef$outcome, COG_LAB) - 1) * 0.6
    rowdef$label <- paste0(short_cog(rowdef$outcome), ", ", rowdef$set)
    pick <- function(spec, set, series) {
        x <- RES[RES$spec == spec & RES$predictor == "module" & RES$outcome %in% COG_LAB, ]
        x$set <- set; x$series <- series; x
    }
    fa <- rbind(pick("S04", "cov4", "Without M4 share"), pick("S06", "cov4", "With M4 share"),
                pick("S05", "+BTC", "Without M4 share"), pick("S07", "+BTC", "With M4 share"),
                pick("S09", "cov4 + %MT", "Without M4 share"), pick("S10", "cov4 + %MT", "With M4 share"))
    fa <- merge(fa, rowdef[, c("set", "outcome", "y")], by = c("set", "outcome"))
    fl <- ifelse(fa$outcome == "Cognitive status", -1, 1)
    fa$b <- fa$beta * fl
    fa$lo <- ifelse(fl < 0, -fa$ci_hi, fa$ci_lo); fa$hi <- ifelse(fl < 0, -fa$ci_lo, fa$ci_hi)
    fa$series <- factor(fa$series, c("Without M4 share", "With M4 share"))
    fa$yy <- fa$y + ifelse(fa$series == "Without M4 share", -0.2, 0.2)
    sc <- c("Without M4 share" = COL$original, "With M4 share" = COL$variant)
    # Readability as in the verdicts: the largest full-model VIF of the model.
    vm <- stats::aggregate(list(vif_max = RES$vif_full), by = list(spec = RES$spec, dataset = RES$dataset,
                                                                  outcome = RES$outcome),
                           FUN = function(v) if (anyNA(v)) Inf else max(v))
    fa <- merge(fa, vm, by = c("spec", "dataset", "outcome"), all.x = TRUE)
    fa$readable <- factor(ifelse(is.finite(fa$vif_max) & fa$vif_max < VIF_MAX, "VIF < 5", "VIF \u22655 (unreadable)"),
                          c("VIF < 5", "VIF \u22655 (unreadable)"))
    fa$fillc <- ifelse(sig(fa$fdr) & fa$readable == "VIF < 5", unname(sc[as.character(fa$series)]), "white")
    ysc <- scale_y_reverse(breaks = rowdef$y, labels = rowdef$label, expand = expansion(add = 0.6))
    pA <- ggplot(fa) +
        geom_vline(xintercept = 0, linetype = "dashed", colour = "grey45", linewidth = 0.3) +
        geom_segment(aes(x = lo, xend = hi, y = yy, yend = yy, colour = series), linewidth = 0.45) +
        geom_point(aes(x = b, y = yy, colour = series, fill = fillc, shape = readable), size = 1.7, stroke = 0.5) +
        scale_colour_manual(values = sc, name = NULL) + scale_fill_identity() + ysc +
        scale_shape_manual(values = c("VIF < 5" = 21, "VIF \u22655 (unreadable)" = 23), name = NULL, drop = FALSE) +
        facet_grid(. ~ dataset) +
        labs(x = "Standardized \u03b2 of Module post-IP% (95% CI)", y = NULL,
             title = "Module post-IP% \u03b2, without and with the M4 share") + th +
        guides(colour = guide_legend(order = 1, override.aes = list(fill = unname(sc), shape = 21)),
               shape = guide_legend(order = 2, override.aes = list(fill = "white", colour = "black"))) +
        theme(legend.box = "vertical", legend.spacing.y = unit(0.5, "mm"))

    # (B) change in beta, donor bootstrap
    fb <- merge(BT, rowdef[, c("set", "outcome", "y")], by = c("set", "outcome"))
    fl <- ifelse(fb$outcome == "Cognitive status", -1, 1)
    fb$d <- fb$delta * fl
    fb$dlo <- ifelse(fl < 0, -fb$delta_hi, fb$delta_lo); fb$dhi <- ifelse(fl < 0, -fb$delta_lo, fb$delta_hi)
    fb$fillc <- ifelse(sig(fb$fdr_boot), COL$variant, "white")
    pB <- ggplot(fb) +
        geom_vline(xintercept = 0, linetype = "dashed", colour = "grey45", linewidth = 0.3) +
        geom_segment(aes(x = dlo, xend = dhi, y = y, yend = y), colour = COL$variant, linewidth = 0.45) +
        geom_point(aes(x = d, y = y, fill = fillc), shape = 21, colour = COL$variant, size = 1.7, stroke = 0.5) +
        scale_fill_identity() + ysc + facet_grid(. ~ dataset) +
        scale_x_continuous(n.breaks = 4) +
        labs(x = "\u0394\u03b2 = \u03b2 with - \u03b2 without (bootstrap 95% CI)", y = NULL,
             title = "Change in \u03b2 with the M4 share") +
        th + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

    # (C), (D) verdict bars
    vC <- V[V$spec %in% c("S11", "S03", "S07") & V$group %in% c("Cognition", "Pathology"), ]
    gC <- ifelse(vC$group == "Cognition", paste0("Cognition\n", vC$block), "Pathology")
    yC <- paste0(vC$dataset, ifelse(vC$spec == "S11", ", vs CP (10B)", ", vs M4 share"))
    pC <- verdict_bars(vC, yC, gC, "Module post-IP% vs CP or M4 share", "CP or M4 share only",
               c("NeoC, vs CP (10B)", "NeoC, vs M4 share", "NonNeoC, vs CP (10B)", "NonNeoC, vs M4 share"),
               c("Cognition\ncov4", "Cognition\n+BTC", "Pathology"), th)
    vD <- V[V$spec %in% c("S09", "S10"), ]
    yD <- paste0(vD$dataset, ifelse(vD$spec == "S10", ", + M4 share", ", Fig 14D"))
    pD <- verdict_bars(vD, yD, vD$group, "Module vs %MT (Fig 14D) \u00b1 M4 share", "%MT only",
               c("NeoC, Fig 14D", "NeoC, + M4 share", "NonNeoC, Fig 14D", "NonNeoC, + M4 share"),
               c("Cognition", "Pathology", "APOE"), th)

    # (E) decomposition tiles
    cols_E <- c("D01", "D02", "D03", "D06", "D07")
    fe <- H[substr(H$spec, 1, 3) %in% cols_E, ]
    fe$col <- factor(substr(fe$spec, 1, 3), cols_E,
                     labels = c("CP + AUCell\n(Fig 10C)", "M4 share\n+ AUCell", "CP + M4 share\n+ AUCell",
                                "M4 share\n+ UCell", "M4 share\n+ AMS"))
    fe$row <- factor(paste0(fe$block, "; ", SHORT[fe$outcome]),
                     rev(c("cov4; MMSE", "cov4; CASI", "+BTC; MMSE", "+BTC; CASI")))
    fe$status <- factor(fe$status, names(HOLD_FILL))
    fe$lab <- sprintf("%s\n\u03b2 %s vs %s\nFDR %s, VIF %s", fe$status, fmt(fe$beta_module), fmt(fe$beta_score),
                      fmt_p(fe$fdr_module), fmt(fe$vif_max, 2))
    fe$txt <- ifelse(fe$status == "HOLDS", "white", "black")
    pE <- ggplot(fe, aes(x = col, y = row)) +
        geom_tile(aes(fill = status), colour = "black", linewidth = 0.3, width = 0.94, height = 0.9) +
        geom_text(aes(label = lab, colour = txt), size = 1.75, lineheight = 0.85) +
        scale_fill_manual(values = HOLD_FILL, drop = FALSE, name = NULL,
                          labels = c("HOLDS (FDR<0.05, expected sign, VIF<5)", "doesn't hold", "unreadable (VIF \u22655)")) +
        scale_colour_identity() + scale_x_discrete(position = "top") + facet_grid(dataset ~ .) +
        labs(x = NULL, y = NULL, title = "Decomposition (Figures 10C-D): Module post-IP% with abundance and a gene score") +
        th + theme(panel.grid = element_blank(), axis.ticks = element_blank(),
                   strip.text.y = element_text(angle = -90))

    # (F) collinearity
    COLL$pair <- factor(COLL$pair, rev(unique(COLL$pair)))
    pF <- ggplot(COLL, aes(y = pair)) +
        geom_vline(xintercept = sqrt(0.8), linetype = "dashed", colour = COL$vifline, linewidth = 0.35) +
        geom_segment(aes(x = 0, xend = abs(r), yend = pair), colour = COL$cg, linewidth = 0.5) +
        geom_point(aes(x = abs(r)), colour = COL$cg, size = 1.8) +
        geom_text(aes(x = abs(r), label = paste0("VIF=", fmt(vif_pair, 1))), hjust = -0.2, size = 1.9) +
        scale_x_continuous(limits = c(0, 1.18), breaks = seq(0, 1, 0.25)) + facet_grid(. ~ dataset) +
        labs(x = "|r| (Pearson, donor level)", y = NULL,
             title = "Collinearity of the M4 share (dashed line: VIF 5, |r| = 0.894)") + th

    fig <- (pA + pB + plot_layout(widths = c(0.63, 0.37))) / (pC + pD + plot_layout(widths = c(0.5, 0.5))) /
        pE / pF + plot_layout(heights = c(1.1, 0.62, 1.45, 0.6)) +
        plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")")
    dev_pdf <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
    ggsave(OUT_PDF, fig, width = 180, height = 270, units = "mm", device = dev_pdf)
    ggsave(OUT_PNG, fig, width = 180, height = 270, units = "mm", dpi = 300, bg = "white")
}

LEGEND_TXT <- c(
    "Supplementary Figure [n]. Does the donor's M4 share explain the M4 results? M4 only. M4 share: the donor's M4 cells (M4-selected Micro-PVM_3 cells with MT-CO3 >0) as a percentage of the donor's Micro-PVM_3 cells, per region set; Module post-IP% has the same denominator, whereas %MT is taken from the M4 cells only (log1p of the donor mean of per-cell %MT over the M4 cells, as in Figure 14D). Outcome and predictors were standardized after the complete-case filter; cov4, age at death, sex, PMI and APOE \u03b54 dose; +BTC, cov4 plus Braak, Thal and CERAD. Unreadable, a full-model VIF \u22655.",
    "(A) Standardized \u03b2 of Module post-IP% on cognition without (dark) and with (violet) the M4 share: cov4 and +BTC (as Figure 9B; BH within region set and block over the cognition outcomes) and cov4 with %MT (as Figure 14D; BH over the nine clinical outcomes). With the M4 share in the model, the coefficient is the association of a larger Module post-IP fraction of the donor's M4 cells at a fixed M4 share, because Module post-IP% = M4 share x within-M4 fraction. \u2020Cognitive status is sign-reversed so that lower \u03b2 means worse cognition throughout. Bars are 95% CIs from the t distribution; filled symbols, FDR <0.05 in a readable model; diamonds, full-model VIF \u22655 (unreadable, never filled).",
    "(B) Change in the Module post-IP% coefficient when the M4 share is added (\u03b2 with minus \u03b2 without; cognitive status sign-reversed as in panel A), with percentile 95% CIs from a donor bootstrap (B = 10,000) in which both models were refitted on each resample; filled circles, FDR <0.05 (BH within region set over the nine tests).",
    "(C) Module post-IP% with CP (as Figure 10B) or with the M4 share in one model: outcomes counted by verdict at FDR <0.05, at cov4 for cognition and pathology (BH within region set and predictor over the 15 outcomes of Figure 10) and at +BTC for cognition (BH over the cognition outcomes).",
    "(D) Module post-IP% and %MT in one model, as Figure 14D, without and with the M4 share as a further predictor; verdicts on the Module and %MT coefficients as Figure 14D (BH over the nine clinical outcomes).",
    "(E) Decomposition as Figures 10C and 10D: Module post-IP% with CP and AUCell (original), with the M4 share in place of CP, with both, and with the M4 share and UCell or AddModuleScore. HOLDS, FDR <0.05 with the expected sign and every VIF <5 (BH within region set, block and model over MMSE and CASI); each box gives the \u03b2 of Module post-IP% and of the gene score, the FDR and the largest VIF.",
    "(F) Pearson |r| between the M4 share and Module post-IP%, %MT, CP, AUCell and the within-M4 fraction (donors with at least 5 M4 cells), and between Module post-IP% and %MT (as Figure 14E), on donors with both values; VIF = 1/(1 - r\u00b2). Values are in the supplementary table.")

# == Figure 2 (post-hoc): the Module state against %MT within the M4 cells ====
draw_figure2 <- function(RES, V, DF, SHORT, COG_LAB, mmse_col) {
    th <- p13_theme()
    sig <- function(x) !is.na(x) & x < 0.05
    short_cog <- function(lab) ifelse(lab == "Cognitive status", "CogS\u2020", SHORT[lab])
    SET2 <- c(S09 = "Module + %MT (Fig 14D)", S15 = "fraction + %MT, cov4", S16 = "fraction + %MT, +BTC")
    STATE <- c(S09 = "module", S15 = "frac", S16 = "frac")
    SER <- c("Module state (Module post-IP% or within-M4 fraction)", "%MT of the M4 cells")
    sc <- stats::setNames(c(COL$variant, COL$original), SER)

    # (A) both coefficients in one model, cognition
    rowdef <- expand.grid(set = unname(SET2), outcome = COG_LAB, stringsAsFactors = FALSE)
    rowdef$y <- seq_len(nrow(rowdef)) + (match(rowdef$outcome, COG_LAB) - 1) * 0.6
    rowdef$label <- paste0(short_cog(rowdef$outcome), ", ", rowdef$set)
    vm <- stats::aggregate(list(vif_max = RES$vif_full), by = list(spec = RES$spec, dataset = RES$dataset,
                                                                  outcome = RES$outcome),
                           FUN = function(v) if (anyNA(v)) Inf else max(v))
    fa <- do.call(rbind, lapply(names(SET2), function(sp) {
        x <- RES[RES$spec == sp & RES$outcome %in% COG_LAB & RES$predictor %in% c(STATE[[sp]], "mt"), ]
        x$set <- SET2[[sp]]
        x$series <- ifelse(x$predictor == "mt", SER[2], SER[1])
        x
    }))
    fa <- merge(fa, rowdef[, c("set", "outcome", "y")], by = c("set", "outcome"))
    fa <- merge(fa, vm, by = c("spec", "dataset", "outcome"), all.x = TRUE)
    fl <- ifelse(fa$outcome == "Cognitive status", -1, 1)
    fa$b <- fa$beta * fl
    fa$lo <- ifelse(fl < 0, -fa$ci_hi, fa$ci_lo); fa$hi <- ifelse(fl < 0, -fa$ci_lo, fa$ci_hi)
    fa$series <- factor(fa$series, SER)
    fa$yy <- fa$y + ifelse(fa$series == SER[1], -0.2, 0.2)
    fa$readable <- factor(ifelse(is.finite(fa$vif_max) & fa$vif_max < VIF_MAX, "VIF < 5", "VIF \u22655 (unreadable)"),
                          c("VIF < 5", "VIF \u22655 (unreadable)"))
    fa$fillc <- ifelse(sig(fa$fdr) & fa$readable == "VIF < 5", unname(sc[as.character(fa$series)]), "white")
    pA <- ggplot(fa) +
        geom_vline(xintercept = 0, linetype = "dashed", colour = "grey45", linewidth = 0.3) +
        geom_segment(aes(x = lo, xend = hi, y = yy, yend = yy, colour = series), linewidth = 0.45) +
        geom_point(aes(x = b, y = yy, colour = series, fill = fillc, shape = readable), size = 1.7, stroke = 0.5) +
        scale_colour_manual(values = sc, name = NULL) + scale_fill_identity() +
        scale_shape_manual(values = c("VIF < 5" = 21, "VIF \u22655 (unreadable)" = 23), name = NULL, drop = FALSE) +
        scale_y_reverse(breaks = rowdef$y, labels = rowdef$label, expand = expansion(add = 0.6)) +
        facet_grid(. ~ dataset) +
        labs(x = "Standardized \u03b2 in one model (95% CI)", y = NULL,
             title = "Module state and %MT of the same M4 cells in one model (post-hoc)") + th +
        guides(colour = guide_legend(order = 1, override.aes = list(fill = unname(sc), shape = 21)),
               shape = guide_legend(order = 2, override.aes = list(fill = "white", colour = "black"))) +
        theme(legend.box = "vertical", legend.spacing.y = unit(0.5, "mm"))

    # (B) verdicts
    vB <- V[V$spec %in% c("S09", "S15", "S16") & V$group %in% c("Cognition", "Pathology", "APOE"), ]
    gB <- ifelse(vB$spec == "S16", "Cognition\n+BTC", vB$group)
    yB <- paste0(vB$dataset, ifelse(vB$spec == "S09", ", Module (14D)", ", fraction"))
    pB <- verdict_bars(vB, yB, gB, "Module state vs %MT, one model", "%MT only",
                       c("NeoC, Module (14D)", "NeoC, fraction", "NonNeoC, Module (14D)", "NonNeoC, fraction"),
                       c("Cognition", "Cognition\n+BTC", "Pathology", "APOE"), th,
                       first_label = "State only (Module or fraction)")

    # (C) the two quantities per donor
    dc <- DF[is.finite(DF$frac) & is.finite(DF$mt), c("dataset", "frac", "mt", mmse_col)]
    names(dc)[4] <- "mmse"; dc$mmse <- num(dc$mmse)
    lab <- do.call(rbind, lapply(split(dc, dc$dataset), function(d)
        data.frame(dataset = d$dataset[1], txt = sprintf("r = %.2f, n = %d", stats::cor(d$frac, d$mt), nrow(d)))))
    pC <- ggplot(dc, aes(x = frac, y = mt)) +
        geom_point(aes(fill = mmse), shape = 21, colour = "grey30", size = 1.6, stroke = 0.25) +
        geom_text(data = lab, aes(label = txt), x = -Inf, y = Inf, hjust = -0.08, vjust = 1.4, size = 2.1,
                  inherit.aes = FALSE) +
        scale_fill_gradient(low = COL$variant, high = "#F2F2F2", na.value = "white", name = "MMSE") +
        facet_grid(. ~ dataset) +
        labs(x = "Within-M4 fraction (%)", y = "%MT of the M4 cells (log1p)",
             title = "Per donor (\u22655 M4 cells)") + th +
        theme(legend.key.width = unit(5, "mm"))

    fig <- pA / (pB + pC + plot_layout(widths = c(0.55, 0.45))) +
        plot_layout(heights = c(1.0, 0.8)) +
        plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")")
    dev_pdf <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
    ggsave(OUT_PDF2, fig, width = 180, height = 190, units = "mm", device = dev_pdf)
    ggsave(OUT_PNG2, fig, width = 180, height = 190, units = "mm", dpi = 300, bg = "white")
}

LEGEND2_TXT <- c(
    "Supplementary Figure [n+1]. Post-hoc test: the Module state against %MT within the same M4 cells. Added after the M4 share analysis had been read, with its own readings file fixed before these models were first fitted. Within-M4 fraction: the donor's Module post-IP cells as a percentage of the donor's M4 cells, per region set, in donors with at least 5 M4 cells; %MT: log1p of the donor mean of per-cell %MT over the same M4 cells, as in Figure 14D. Both come from the same cells, so neither carries the abundance of M4 cells. A Module post-IP cell is one with cDFG, the mean of the eight mitochondrially encoded partner genes, at or above IP_y, so the fraction and %MT partly measure the same reads; a Both or %MT only verdict can reflect that overlap. Outcome and predictors were standardized after the complete-case filter; covariates as Supplementary Figure [n].",
    "(A) Standardized \u03b2 of the Module state (violet: Module post-IP% in the Figure 14D model, or the within-M4 fraction) and of %MT (dark) in one model, for cognition at cov4 (BH within region set and predictor over the nine clinical outcomes, as Figure 14D) and at +BTC (BH over the cognition outcomes, as Figure 9B). \u2020Cognitive status is sign-reversed so that lower \u03b2 means worse cognition throughout. Bars are 95% CIs from the t distribution; filled symbols, FDR <0.05 in a readable model; diamonds, full-model VIF \u22655 (unreadable, never filled).",
    "(B) Verdicts at FDR <0.05 for the Module state and %MT in one model: Module post-IP% (Figure 14D) and the within-M4 fraction, for cognition, pathology and APOE \u03b54 dose at cov4 and for cognition at +BTC. Unreadable, full-model VIF \u22655.",
    "(C) Within-M4 fraction against %MT of the M4 cells per donor, shaded by MMSE (white, missing MMSE); Pearson r and the number of donors in each panel. Values are in the supplementary table (sheets C_Models, D_Verdicts and L_PostHoc).")

# == Run =====================================================================
.reading <- tryCatch(
    withCallingHandlers(run_part13(), warning = function(w) {
        cat("  WARNING: ", conditionMessage(w), "\n", sep = ""); invokeRestart("muffleWarning")
    }),
    error = function(e) {
        cat("\nRun stopped: ", conditionMessage(e), "\n", sep = "")
        structure(conditionMessage(e), class = "p13_error")
    })
.note <- if (inherits(.reading, "p13_error")) paste("STOPPED -", .reading) else paste("completed; reading:", .reading)
write_ledger(.note)
cat("\nGuardrail ledger : ", LEDGER, "\nGuardrail detail : ", GDETAIL, "\nConsole log      : ", CONSOLE,
    "\nRun: ", .note, "\n\n", sep = "")
cat(utils::capture.output(utils::sessionInfo()), sep = "\n")
sink(); close(.con)
