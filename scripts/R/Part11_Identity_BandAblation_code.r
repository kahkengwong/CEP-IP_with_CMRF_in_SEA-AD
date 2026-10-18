# =============================================================================
# Part11_Identity_BandAblation_code.r
# SEA-AD CEP-IP/CMRF, Part 11, boundary ablation (Figure 11A and 11C, Table S12)
#   Compares four Module boundaries on the same cells: A0 OLS corridor, A1
#   penalised spline with constant width, A2 spline with locally estimated
#   width, and A3 CMRF (Validation 1 cell_class). Pre/Post uses one rule for
#   all four. Reports set overlap, donor associations with clinical outcomes
#   and negative controls, a scorecard, and repeated-measures rung tests.
#
# Usage (from the SEA-AD folder, or set SEAAD_ROOT), after Part 10e
#   Validation 1 and Part 10h
#   source("Part11_Identity_BandAblation_code.r")
#
# Outputs (Part11_Identity/Part11_Identity_BandAblation_P10hDenom/): figures
#   (PDF, PNG), workbook, results RDS and log.
# =============================================================================

library(readxl)
library(writexl)
library(ggplot2)
library(mgcv)
library(grid)          # unit(), for panel spacing in the association figure

CEP_HOME    <- Sys.getenv("SEAAD_ROOT", getwd())
V1_ROOT     <- file.path(CEP_HOME, "Part10e_Validation1")
EXPORT_ROOT <- file.path(CEP_HOME, "R_export_CP10K")
P10H_ROOT   <- file.path(CEP_HOME, "Part10h_Prognosis")
# Donor percentage denominator. "supertype": the donor's supertype cells
# (Part 10h prop_df n_total), donors with no arm cell score 0. "arm": the
# donor's GOI-positive cells on the arm. The modes write to different folders.
DENOM_MODE <- "supertype"
# TRUE: clinical outcomes and negative controls form separate BH families.
BH_SPLIT_NEG <- TRUE
OUT_DIR     <- file.path(CEP_HOME, "Part11_Identity",
                         if (identical(DENOM_MODE, "supertype"))
                             "Part11_Identity_BandAblation_P10hDenom"
                         else "Part11_Identity_BandAblation")

V1DIR  <- c(M4 = "Mic13_PVM3", PVM_3_all = "PVM_3",
            PVM_1_all = "PVM_1", HKG8 = "HKG8")
ARMDIR <- c(M4 = "Mic13_PVM3", PVM_3_all = "PVM3_all",
            PVM_1_all = "PVM1_all", HKG8 = "HKG8_matched")
STYPE  <- c(M4 = "Micro-PVM_3-SEAAD", PVM_3_all = "Micro-PVM_3-SEAAD",
            PVM_1_all = "Micro-PVM_1", HKG8 = "Micro-PVM_3-SEAAD")
P10HFILE <- c(M4 = "Part10h_CP10K_M4_Prognosis.xlsx",
              PVM_3_all = "Part10h_CP10K_PVM_3_Prognosis.xlsx",
              PVM_1_all = "Part10h_CP10K_PVM_1_Prognosis.xlsx",
              HKG8 = "Part10h_CP10K_HKG8_Prognosis.xlsx")
# Pinned Part 10h results RDS (the only frame with CASI); other groups are
# found by the fallbacks below.
P10HRDS <- list(
    M4 = file.path(P10H_ROOT, "Part10h_Prognosis_Mic13_PVM3",
                   "Part10h_CP10K_M4_Prognosis_results.rds"))
REGS <- list(NeoC = c("DFC", "MTG", "AnG", "STG", "V1C", "ITG"),
             NonNeoC = c("HIP", "MEC", "LEC", "FI"))

MODULE_CLASSES <- c("Original + Refined", "Recovered")
# A0 to A2 coverage is matched to A3's module fraction per group and dataset,
# so all four select the same number of cells. A number forces a fixed value.
COVERAGE_MODE <- "match_A3"
N_BINS_A2  <- 20L     # equal-count u bins for the local width
MIN_DONORS <- 25L     # below this a regression row is flagged
ALPHA      <- 0.05
RUNGS      <- c("A0_OLS", "A1_spline_fixed", "A2_spline_local", "A3_CMRF")
# Figure constants: group order, corridor colour, FDR < ALPHA colour.
GRP_ORDER    <- c("M4", "PVM_3", "PVM_1", "HKG8")
CORRIDOR_COL <- "#CC0000"
SIG_COL      <- "#6A3D9A"
# Also print each figure to the current device (FALSE for an unattended run).
SHOW_PLOTS <- TRUE
show_plot  <- function(p) if (isTRUE(SHOW_PLOTS)) print(p)

# Full outcome panel: clinical outcomes and negative controls, kept apart.
OUT_CLIN <- c("Cognitive status" = "CogStatus_num", "MMSE score" = "MMSE_num",
              "CASI score" = "CASI_num", "CERAD score" = "CERAD_num",
              "Braak stage" = "Braak_num", "Thal phase" = "Thal_num",
              "ADNC rank" = "ADNC_rank", "CAA score" = "CAA_num",
              "APOE e4 dose" = "APOE_e4_dose")
OUT_NEG  <- c("Age at death" = "Age_num", "Sex" = "Sex_num",
              "Brain pH" = "BrainPH_num", "PMI" = "PMI_num",
              "RIN" = "RIN_num", "Fresh brain weight" = "BrainWeight_num")
COG_LABS  <- c("Cognitive status", "MMSE score", "CASI score")
PATH_LABS <- c("CERAD score", "Braak stage", "Thal phase", "ADNC rank",
               "CAA score")
FIG_GROUP  <- "M4"    # which panel the band-shape figure draws
FIG_DS     <- "NeoC"

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
# Log written by append, so a partial run keeps its log.
LOG <- file.path(OUT_DIR, paste0("Part11_BandAblation_log_",
                                 format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
cat("", file = LOG)
say <- function(...) {
    m <- paste0(...)
    cat(m, "\n", sep = "")
    cat(m, "\n", sep = "", file = LOG, append = TRUE)
}

# Column picker: prefer a numeric-capable column, then a `_num` variant, then
# the shortest name. need_num = FALSE for columns that may be a factor.
pick <- function(df, pat, need_num = TRUE) {
    nms <- if (is.data.frame(df)) names(df) else df
    i <- grep(pat, nms, ignore.case = TRUE)
    if (length(i) == 0L) return(NA_character_)
    cand <- nms[i]
    if (is.data.frame(df) && need_num) {
        okn <- vapply(cand, function(cn) {
            x <- suppressWarnings(as.numeric(as.character(df[[cn]])))
            sum(is.finite(x)) >= 0.8 * length(x) &&
                stats::sd(x, na.rm = TRUE) > 0
        }, logical(1))
        if (any(okn)) cand <- cand[okn]
    }
    numish <- grepl("_num$", cand, ignore.case = TRUE)
    if (any(numish)) cand <- cand[numish]
    cand[order(nchar(cand))][1]
}
jac <- function(a, b) {
    u <- length(union(a, b))
    if (u == 0L) NA_real_ else length(intersect(a, b)) / u
}
num_or_factor <- function(x) {
    n <- suppressWarnings(as.numeric(as.character(x)))
    if (sum(is.finite(n)) >= 0.8 * length(x)) n else factor(as.character(x))
}

donor_map <- function(grp, ds) {
    sub <- if (identical(ds, "NeoC")) "neocortex" else "nonneocortex"
    st  <- gsub("[-/]", "_", STYPE[[grp]])
    bc <- character(0); dn <- character(0)
    for (rg in REGS[[ds]]) {
        f <- file.path(EXPORT_ROOT, ARMDIR[[grp]], sub, "gam_data",
                       paste0(rg, "_", st, "_selected.csv"))
        if (!file.exists(f)) next
        x <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
        if (!all(c("barcode", "donor") %in% names(x))) next
        bc <- c(bc, x$barcode); dn <- c(dn, as.character(x$donor))
    }
    if (!length(bc)) return(NULL)
    k <- !duplicated(bc)
    stats::setNames(dn[k], bc[k])
}

# Local quantile of |residual| over equal-count u bins; fit, then evaluate at
# any u.
local_width_fit <- function(u, r, cov, nb = N_BINS_A2) {
    flat <- as.numeric(stats::quantile(abs(r), cov, na.rm = TRUE))
    br <- unique(stats::quantile(u, probs = seq(0, 1, length.out = nb + 1),
                                 na.rm = TRUE))
    if (length(br) < 3) return(list(ctr = NULL, wid = NULL, flat = flat))
    b   <- cut(u, breaks = br, include.lowest = TRUE)
    ctr <- as.numeric(tapply(u, b, stats::median))
    wid <- as.numeric(tapply(abs(r), b, stats::quantile, probs = cov,
                             na.rm = TRUE))
    ok <- is.finite(ctr) & is.finite(wid)
    if (sum(ok) < 2) return(list(ctr = NULL, wid = NULL, flat = flat))
    list(ctr = ctr[ok], wid = wid[ok], flat = flat)
}
local_width_at <- function(fit, x) {
    if (is.null(fit$ctr)) return(rep(fit$flat, length(x)))
    stats::approx(fit$ctr, fit$wid, xout = x, rule = 2)$y
}

# A3 envelope for the figure: within equal-count u bins, the min and max v of
# A3 cells, with their median as the spine. No smoothing.
a3_envelope <- function(u, v, keep, nb = N_BINS_A2) {
    if (sum(keep) < 2L) return(NULL)
    br <- unique(stats::quantile(u, probs = seq(0, 1, length.out = nb + 1),
                                 na.rm = TRUE))
    if (length(br) < 3) return(NULL)
    b  <- cut(u, breaks = br, include.lowest = TRUE)
    uk <- u[keep]; vk <- v[keep]; bk <- b[keep]
    ctr <- as.numeric(tapply(uk, bk, stats::median))
    lo  <- as.numeric(tapply(vk, bk, min))
    hi  <- as.numeric(tapply(vk, bk, max))
    md  <- as.numeric(tapply(vk, bk, stats::median))
    ok  <- is.finite(ctr) & is.finite(lo) & is.finite(hi) & is.finite(md)
    if (sum(ok) < 2L) return(NULL)
    data.frame(u = ctr[ok], lo = lo[ok], hi = hi[ok], fit = md[ok])
}

# =============================================================================
# STEP 1: build the four cell sets
# =============================================================================
say(strrep("=", 76))
say("Part 11 Identity -- band ablation. Started ", format(Sys.time()))
say(strrep("=", 76))

cells <- list(); ovl <- list(); fitcurve <- list()
cover_tab <- list()
for (grp in names(V1DIR)) {
    for (ds in c("NeoC", "NonNeoC")) {
        dpath <- file.path(V1_ROOT, paste0("Part10e_Validation1_",
                                           V1DIR[[grp]]))
        f <- list.files(dpath, full.names = TRUE,
                        pattern = paste0("_V1_", ds,
                                         "_SK_interaction_analysis\\.xlsx$"))
        if (!length(f)) { say("  MISSING: ", grp, " ", ds); next }
        d <- as.data.frame(read_excel(f[1], "Step4_Refined_Cells"))
        need <- c("barcode", "cell_class", "u_all", "v_all", "cDFG_expr")
        miss <- setdiff(need, names(d))
        if (length(miss)) {
            say("  ", grp, " ", ds, " lacks: ", paste(miss, collapse = ", "))
            next
        }
        d <- d[stats::complete.cases(d[, c("u_all", "v_all", "cDFG_expr")]) &
               !is.na(d$cell_class), ]

        # == Pre/Post, one rule for all rungs ===============================
        lab <- ifelse(grepl("Post-IP", d$group_refined, fixed = TRUE), "Post",
               ifelse(grepl("Pre-IP",  d$group_refined, fixed = TRUE), "Pre",
               ifelse(grepl("Post-IP", d$IP_status, fixed = TRUE), "Post",
               ifelse(grepl("Pre-IP",  d$IP_status, fixed = TRUE), "Pre",
                      NA_character_))))
        hi <- suppressWarnings(max(d$cDFG_expr[!is.na(lab) & lab == "Pre"]))
        lo <- suppressWarnings(min(d$cDFG_expr[!is.na(lab) & lab == "Post"]))
        ip_y <- if (is.finite(hi) && is.finite(lo) && hi < lo) (hi + lo) / 2
                else stats::median(d$cDFG_expr)
        d$side <- ifelse(d$cDFG_expr < ip_y, "Pre", "Post")

        u <- d$u_all; v <- d$v_all
        fd <- data.frame(u = u, v = v)
        d$A3_CMRF <- d$cell_class %in% MODULE_CLASSES
        cover <- if (identical(COVERAGE_MODE, "match_A3")) mean(d$A3_CMRF)
                 else as.numeric(COVERAGE_MODE)
        if (!is.finite(cover) || cover <= 0 || cover >= 1) cover <- 0.95

        m0 <- stats::lm(v ~ u, data = fd)
        r0 <- as.numeric(stats::residuals(m0))
        w0 <- as.numeric(stats::quantile(abs(r0), cover, na.rm = TRUE))
        m1 <- mgcv::gam(v ~ s(u, bs = "tp"), data = fd, method = "REML")
        r1 <- as.numeric(stats::residuals(m1))
        w1 <- as.numeric(stats::quantile(abs(r1), cover, na.rm = TRUE))
        lw <- local_width_fit(u, r1, cover)
        w2 <- local_width_at(lw, u)

        d$A0_OLS          <- abs(r0) <= w0
        d$A1_spline_fixed <- abs(r1) <= w1
        d$A2_spline_local <- abs(r1) <= w2

        dm <- donor_map(grp, ds)
        d$donor <- if (is.null(dm)) NA_character_ else unname(dm[d$barcode])
        d$group <- grp; d$dataset <- ds
        cells[[paste(grp, ds)]] <- d

        say(sprintf(paste0("  %-10s %-8s n=%5d IP_y=%.4f cover=%.3f |",
                           " module n: A0 %5d  A1 %5d  A2 %5d  A3 %5d"),
                    grp, ds, nrow(d), ip_y, cover, sum(d$A0_OLS),
                    sum(d$A1_spline_fixed), sum(d$A2_spline_local),
                    sum(d$A3_CMRF)))
        cover_tab[[length(cover_tab) + 1L]] <- data.frame(
            group = grp, dataset = ds, n_cells = nrow(d), ip_y = ip_y,
            target_coverage = cover, width_A0 = w0, width_A1 = w1,
            width_A2_min = min(w2), width_A2_max = max(w2),
            n_A0 = sum(d$A0_OLS), n_A1 = sum(d$A1_spline_fixed),
            n_A2 = sum(d$A2_spline_local), n_A3 = sum(d$A3_CMRF),
            stringsAsFactors = FALSE)

        for (i in seq_along(RUNGS)) for (j in seq_along(RUNGS)) {
            if (j <= i) next
            ovl[[length(ovl) + 1L]] <- data.frame(
                group = grp, dataset = ds, a = RUNGS[i], b = RUNGS[j],
                n_a = sum(d[[RUNGS[i]]]), n_b = sum(d[[RUNGS[j]]]),
                jaccard = jac(d$barcode[d[[RUNGS[i]]]],
                              d$barcode[d[[RUNGS[j]]]]),
                stringsAsFactors = FALSE)
        }
        if (identical(grp, FIG_GROUP) && identical(ds, FIG_DS)) {
            # lo and hi per rung, since A3's envelope is not symmetric.
            ug  <- seq(min(u), max(u), length.out = 300)
            f0  <- as.numeric(stats::predict(m0, data.frame(u = ug)))
            f1  <- as.numeric(stats::predict(m1, data.frame(u = ug)))
            w2g <- local_width_at(lw, ug)
            fc <- rbind(
                data.frame(rung = "A0_OLS", u = ug,
                           fit = f0, lo = f0 - w0, hi = f0 + w0),
                data.frame(rung = "A1_spline_fixed", u = ug,
                           fit = f1, lo = f1 - w1, hi = f1 + w1),
                data.frame(rung = "A2_spline_local", u = ug,
                           fit = f1, lo = f1 - w2g, hi = f1 + w2g))
            e3 <- a3_envelope(u, v, d$A3_CMRF)
            if (!is.null(e3))
                fc <- rbind(fc, data.frame(rung = "A3_CMRF", u = e3$u,
                                           fit = e3$fit, lo = e3$lo,
                                           hi = e3$hi))
            else
                say("  A3 envelope not drawable; panel omitted")
            fitcurve[[1]] <- fc
        }
    }
}
if (!length(cells)) stop("Nothing was read; check V1_ROOT.", call. = FALSE)
OVL <- do.call(rbind, ovl)

say("")
# The three real groups and the HKG8 control are reported separately.
say("MEDIAN JACCARD AGAINST A3_CMRF, by rung")
for (tag in c("three real groups", "HKG8 control")) {
    for (rg in RUNGS[1:3]) {
        z <- OVL[OVL$b == "A3_CMRF" & OVL$a == rg &
                 ((tag == "HKG8 control") == (OVL$group == "HKG8")), ]
        if (!nrow(z)) next
        say(sprintf("  %-18s %-16s median %.3f  range %.3f - %.3f", tag, rg,
                    stats::median(z$jaccard), min(z$jaccard), max(z$jaccard)))
    }
}

# =============================================================================
# STEP 2: donor percentages and the clinical table
# =============================================================================
donor_pct <- list()
for (nm in names(cells)) {
    d <- cells[[nm]]
    d <- d[!is.na(d$donor), ]
    if (!nrow(d)) next
    tot <- table(d$donor)
    for (rg in RUNGS) {
        dn <- names(tot)
        cnt <- function(k) {
            tb <- table(d$donor[k])
            v <- as.integer(tb[match(dn, names(tb))])
            v[is.na(v)] <- 0L
            v
        }
        nt <- as.integer(tot[dn])
        donor_pct[[length(donor_pct) + 1L]] <- data.frame(
            group = d$group[1], dataset = d$dataset[1], rung = rg, donor = dn,
            n_arm = nt,
            n_mod_post = cnt(d[[rg]] & d$side == "Post"),
            n_nm_post = cnt(!d[[rg]] & d$side == "Post"),
            n_denominator = nt,
            pct_mod_post = 100 * cnt(d[[rg]] & d$side == "Post") / nt,
            pct_nm_post = 100 * cnt(!d[[rg]] & d$side == "Post") / nt,
            stringsAsFactors = FALSE)
    }
}
DP <- do.call(rbind, donor_pct)

# Find each Part 10h workbook by exact file name under CEP_HOME.
clin <- list(); clin_path <- list()
for (grp in names(P10HFILE)) {
    pat <- paste0("^", gsub(".", "\\.", P10HFILE[[grp]], fixed = TRUE), "$")
    f <- list.files(CEP_HOME, recursive = TRUE, full.names = TRUE,
                    pattern = pat)
    if (!length(f)) {
        say("  clinical MISSING for ", grp, " (", P10HFILE[[grp]], ")")
        next
    }
    say("  clinical ", grp, ": ", f[1])
    clin_path[[grp]] <- f[1]
    t <- as.data.frame(read_excel(f[1], "Donor_Table"))
    if (!length(clin)) say("  Donor_Table columns: ",
                           paste(names(t), collapse = " | "))
    clin[[grp]] <- t
}
if (!length(clin)) stop("No Donor_Table found under ", P10H_ROOT, call. = FALSE)
# Part 8's donor_clinical.csv, merged by donor ID on the Donor_Table route.
dc <- NULL
dcf <- list.files(EXPORT_ROOT, recursive = TRUE, full.names = TRUE,
                  pattern = "^donor_clinical\\.csv$")
if (length(dcf)) {
    dc <- read.csv(dcf[1], stringsAsFactors = FALSE, check.names = FALSE)
    say("  donor_clinical: ", dcf[1])
    say("  donor_clinical columns: ", paste(names(dc), collapse = " | "))
} else {
    say("  donor_clinical.csv not found under ", EXPORT_ROOT)
}

merge_clinical <- function(ct, cd_name) {
    if (is.null(dc)) return(ct)
    dcid <- pick(dc, "^donor", need_num = FALSE)
    if (is.na(dcid)) return(ct)
    a <- as.character(ct[[cd_name]])
    b <- as.character(dc[[dcid]])
    n_ov <- length(intersect(a, b))
    add <- setdiff(names(dc), c(dcid, names(ct)))
    if (!length(add) || n_ov == 0L) {
        say("  clinical merge: ", n_ov, " donor(s) overlap, ", length(add),
            " new column(s) -- nothing added")
        return(ct)
    }
    ct2 <- ct
    idx <- match(a, b)
    for (cn in add) ct2[[cn]] <- dc[[cn]][idx]
    say("  clinical merge: ", n_ov, " of ", length(a),
        " donor(s) matched, ", length(add), " column(s) added")
    ct2
}

# CASI is read from Part 10h's results RDS (datasets[[ds]]$prop_df), per group
# and dataset; Donor_Table is the fallback.
prop_frames <- list()
missing_pf  <- character(0)
for (grp in names(P10HFILE)) {
    if (is.null(clin_path[[grp]]) &&
        (is.null(P10HRDS[[grp]]) || !file.exists(P10HRDS[[grp]]))) next
    # Routes in order: pinned path, beside the workbook, the group's Part 10h
    # folder, then the implied name anywhere under CEP_HOME.
    route_rds <- NA_character_
    rf <- character(0)
    if (!is.null(P10HRDS[[grp]]) && file.exists(P10HRDS[[grp]])) {
        rf <- P10HRDS[[grp]]; route_rds <- "pinned"
    }
    if (!length(rf) && !is.null(clin_path[[grp]])) {
        rf <- list.files(dirname(clin_path[[grp]]), full.names = TRUE,
                         pattern = "_results\\.rds$")
        if (length(rf)) route_rds <- "beside the workbook"
    }
    if (!length(rf) && !is.null(V1DIR[[grp]])) {
        gdir <- file.path(P10H_ROOT, paste0("Part10h_Prognosis_",
                                            V1DIR[[grp]]))
        if (dir.exists(gdir)) {
            rf <- list.files(gdir, full.names = TRUE,
                             pattern = "_results\\.rds$")
            if (length(rf)) route_rds <- "group Part 10h folder"
        }
    }
    if (!length(rf) && !is.null(clin_path[[grp]])) {
        want <- sub("\\.xlsx$", "_results.rds", basename(clin_path[[grp]]))
        rf <- list.files(CEP_HOME, recursive = TRUE, full.names = TRUE,
                         pattern = paste0("^", gsub(".", "\\.", want,
                                                    fixed = TRUE), "$"))
        if (length(rf)) route_rds <- "name search under CEP_HOME"
    }
    if (length(rf) > 1L)
        say("  NOTE: ", length(rf), " results RDS matched for ", grp,
            "; using ", basename(rf[1]))
    if (!length(rf)) {
        say("  WARN: no Part 10h results RDS for ", grp,
            "; using Donor_Table, without CASI")
        missing_pf <- c(missing_pf, grp)
        next
    }
    r <- try(readRDS(rf[1]), silent = TRUE)
    if (inherits(r, "try-error") || is.null(r[["datasets"]])) {
        say("  Part 10h RDS unreadable or has no `datasets`: ", rf[1])
        next
    }
    for (ds in names(r[["datasets"]])) {
        pf <- r[["datasets"]][[ds]][["prop_df"]]
        if (is.null(pf) || !nrow(pf)) next
        prop_frames[[paste(grp, ds)]] <- as.data.frame(pf)
    }
    say("  Part 10h prop_df for ", grp, ": ",
        paste(names(r[["datasets"]]), collapse = ", "), " from ",
        rf[1], " (", route_rds, ")")
}
if (length(missing_pf))
    say("  WARN: CASI unavailable for ", paste(missing_pf, collapse = ", "),
        "; ", 2L * length(missing_pf), " row(s) fewer")
if (length(prop_frames))
    say("  prop_df columns: ",
        paste(names(prop_frames[[1]]), collapse = " | "))

# Columns are resolved on a prop_df when one exists (it carries CASI).
CT <- if (length(prop_frames)) prop_frames[[1]] else clin[[1]]
cd <- pick(CT, "^donor", need_num = FALSE)
if (!length(prop_frames) && !is.na(cd)) CT <- merge_clinical(CT, cd)
CO <- list(MMSE = pick(CT, "mmse"), CASI = pick(CT, "casi"))
CV4 <- c(age = pick(CT, "age"), apoe = pick(CT, "apoe"),
         sex = pick(CT, "^sex", need_num = FALSE), pmi = pick(CT, "^pmi"))
BTC <- c(braak = pick(CT, "braak"), thal = pick(CT, "thal"),
         cerad = pick(CT, "cerad"))
CV4 <- CV4[!is.na(CV4)]; BTC <- BTC[!is.na(BTC)]
say("  resolved: donor=", cd, " | outcomes=", paste(unlist(CO), collapse = ","),
    " | cov4=", paste(unlist(CV4), collapse = ","),
    " | BTC=", paste(unlist(BTC), collapse = ","))
if (is.na(cd) || all(is.na(unlist(CO))))
    stop("Could not resolve donor or outcome columns; see the list above.",
         call. = FALSE)

# =============================================================================
# STEP 2b: donor percentages on Part 10h's denominator
# =============================================================================
DENOM_CHECK <- NULL
if (identical(DENOM_MODE, "supertype")) {
    say("")
    say("DENOMINATOR: donor's supertype cells (Part 10h prop_df n_total)")
    new_dp <- list(); chk <- list()
    for (grp in unique(DP$group)) {
        for (ds in unique(DP$dataset[DP$group == grp])) {
            pf <- prop_frames[[paste(grp, ds)]]
            if (is.null(pf) || !all(c(cd, "n_total") %in% names(pf)))
                stop("DENOM_MODE = \"supertype\" needs Part 10h's prop_df with ",
                     "`n_total` for ", grp, " ", ds, ". Set DENOM_MODE <- ",
                     "\"arm\" to run without it.", call. = FALSE)
            ref <- data.frame(donor = as.character(pf[[cd]]),
                              n_denominator = suppressWarnings(
                                  as.numeric(pf[["n_total"]])),
                              stringsAsFactors = FALSE)
            ref <- ref[is.finite(ref$n_denominator) &
                       ref$n_denominator > 0, , drop = FALSE]
            for (rg in RUNGS) {
                dp <- DP[DP$group == grp & DP$dataset == ds & DP$rung == rg,
                         c("donor", "n_arm", "n_mod_post", "n_nm_post")]
                lost <- setdiff(dp$donor, ref$donor)
                m <- merge(ref, dp, by = "donor", all.x = TRUE)
                for (v in c("n_arm", "n_mod_post", "n_nm_post"))
                    m[[v]][is.na(m[[v]])] <- 0L
                m$pct_mod_post <- 100 * m$n_mod_post / m$n_denominator
                m$pct_nm_post  <- 100 * m$n_nm_post  / m$n_denominator
                m$group <- grp; m$dataset <- ds; m$rung <- rg
                new_dp[[length(new_dp) + 1L]] <- m[, c("group", "dataset",
                    "rung", "donor", "n_arm", "n_mod_post", "n_nm_post",
                    "n_denominator", "pct_mod_post", "pct_nm_post")]
                # At A3 the rebuilt percentage should reproduce Part 10h's
                # pct_post.
                if (identical(rg, "A3_CMRF") && "pct_post" %in% names(pf)) {
                    pp <- suppressWarnings(as.numeric(
                        pf[["pct_post"]][match(m$donor,
                                               as.character(pf[[cd]]))]))
                    ok <- is.finite(pp) & is.finite(m$pct_mod_post)
                    chk[[length(chk) + 1L]] <- data.frame(
                        group = grp, dataset = ds, n_donors = nrow(m),
                        n_zero_arm_cells = sum(m$n_arm == 0),
                        n_arm_donors_not_in_prop_df = length(lost),
                        rs_vs_pct_post = if (sum(ok) > 2)
                            suppressWarnings(stats::cor(m$pct_mod_post[ok],
                                pp[ok], method = "spearman")) else NA_real_,
                        max_abs_diff_pp = if (any(ok))
                            max(abs(m$pct_mod_post[ok] - pp[ok]))
                            else NA_real_,
                        stringsAsFactors = FALSE)
                }
            }
        }
    }
    DP <- do.call(rbind, new_dp); rownames(DP) <- NULL
    DENOM_CHECK <- do.call(rbind, chk)
    if (!is.null(DENOM_CHECK)) for (i in seq_len(nrow(DENOM_CHECK)))
        say(sprintf(paste0("  %-10s %-8s donors %3d (no arm cell: %2d) | A3 vs ",
                           "Part 10h pct_post: rs %.4f, max diff %.3f pp%s"),
                    DENOM_CHECK$group[i], DENOM_CHECK$dataset[i],
                    DENOM_CHECK$n_donors[i],
                    DENOM_CHECK$n_zero_arm_cells[i],
                    DENOM_CHECK$rs_vs_pct_post[i],
                    DENOM_CHECK$max_abs_diff_pp[i],
                    if (isTRUE(DENOM_CHECK$rs_vs_pct_post[i] < 0.99))
                        "  <-- CHECK" else ""))
} else {
    say("")
    say("DENOMINATOR: donor's cells on the arm")
}

# =============================================================================
# STEP 3: regressions
# =============================================================================
fit_row <- function(df, yv, preds, covs) {
    keep <- c(yv, preds, covs)
    # Two patterns can resolve to the same column; keep each name once.
    keep <- unique(keep[!is.na(keep) & keep %in% names(df)])
    z <- df[, keep, drop = FALSE]
    # Outcome made numeric before complete.cases.
    z[[yv]] <- suppressWarnings(as.numeric(as.character(z[[yv]])))
    z <- z[stats::complete.cases(z), , drop = FALSE]
    if (nrow(z) < 10) return(NULL)
    for (cn in setdiff(names(z), c(yv, preds)))
        z[[cn]] <- num_or_factor(z[[cn]])
    drop <- vapply(names(z), function(cn) {
        x <- z[[cn]]
        if (is.factor(x)) nlevels(droplevels(x)) < 2 else
            !is.finite(stats::sd(x)) || stats::sd(x) == 0
    }, logical(1))
    z <- z[, !drop, drop = FALSE]
    if (!all(c(yv, preds) %in% names(z))) return(NULL)
    z[[yv]] <- as.numeric(scale(z[[yv]]))
    for (p in preds) z[[p]] <- as.numeric(scale(z[[p]]))
    cv <- setdiff(names(z), c(yv, preds))
    # Covariates are backticked (they carry spaces); predictors are not.
    rhs <- paste(c(paste0("`", cv, "`"), preds), collapse = " + ")
    m <- try(stats::lm(stats::as.formula(paste0("`", yv, "` ~ ", rhs)), z),
             silent = TRUE)
    if (inherits(m, "try-error")) return(NULL)
    s <- summary(m)$coefficients
    # Delta R2 against the model without the predictors (their joint
    # contribution when there are two).
    rhs0 <- paste0("`", cv, "`", collapse = " + ")
    m0f <- if (nzchar(rhs0))
        try(stats::lm(stats::as.formula(paste0("`", yv, "` ~ ", rhs0)), z),
            silent = TRUE) else NULL
    r2f <- summary(m)$r.squared
    dr2 <- if (is.null(m0f) || inherits(m0f, "try-error")) NA_real_
           else r2f - summary(m0f)$r.squared
    out <- NULL
    for (p in preds) {
        if (!p %in% rownames(s)) next
        out <- rbind(out, data.frame(predictor = p, n = nrow(z),
                                     beta = s[p, 1], p_value = s[p, 4],
                                     delta_r2 = dr2,
                                     stringsAsFactors = FALSE))
    }
    out
}

# Drop a covariate that is the outcome.
drop_self <- function(covs, yv) covs[covs != yv]

assoc <- list(); enc <- list(); clin_route <- list()
for (grp in unique(DP$group)) {
    for (ds in unique(DP$dataset[DP$group == grp])) {
        pf <- prop_frames[[paste(grp, ds)]]
        route <- if (!is.null(pf)) "prop_df" else "Donor_Table"
        ct <- if (!is.null(pf)) pf else
              if (!is.null(clin[[grp]])) clin[[grp]] else CT
        if (is.null(pf) && cd %in% names(ct)) ct <- merge_clinical(ct, cd)
        if (!cd %in% names(ct)) {
            say("  SKIP ", grp, " ", ds, ": no `", cd, "` column (",
                route, ")")
            next
        }
        names(ct)[names(ct) == cd] <- "donor"
        ct$donor <- as.character(ct$donor)
        clin_route[[length(clin_route) + 1L]] <- data.frame(
            group = grp, dataset = ds, route = route, n_donors = nrow(ct),
            has_MMSE = !is.na(CO$MMSE) && CO$MMSE %in% names(ct),
            has_CASI = !is.na(CO$CASI) && CO$CASI %in% names(ct),
            stringsAsFactors = FALSE)
        for (rg in RUNGS) {
            dp <- DP[DP$group == grp & DP$dataset == ds & DP$rung == rg, ]
            mg <- merge(dp, ct, by = "donor")
            for (on in names(CO)) {
                yv <- CO[[on]]
                if (is.na(yv) || !yv %in% names(mg)) next
                for (blk in c("cov4", "BTC")) {
                    cv <- if (identical(blk, "cov4")) unlist(CV4)
                          else c(unlist(CV4), unlist(BTC))
                    r <- fit_row(mg, yv, "pct_mod_post", cv)
                    if (!is.null(r)) {
                        r$group <- grp; r$dataset <- ds; r$rung <- rg
                        r$outcome <- on; r$block <- blk
                        assoc[[length(assoc) + 1L]] <- r
                    }
                    e <- fit_row(mg, yv, c("pct_mod_post", "pct_nm_post"), cv)
                    if (!is.null(e)) {
                        e$group <- grp; e$dataset <- ds; e$rung <- rg
                        e$outcome <- on; e$block <- blk
                        enc[[length(enc) + 1L]] <- e
                    }
                }
            }
        }
    }
}
A <- do.call(rbind, assoc); E <- do.call(rbind, enc)
for (nm in c("A", "E")) {
    x <- get(nm)
    if (is.null(x)) next
    x$FDR <- stats::p.adjust(x$p_value, "BH")
    x$bh_family <- paste("Benjamini-Hochberg over all", nrow(x),
                         "rows of this sheet")
    x$underpowered <- x$n < MIN_DONORS
    assign(nm, x)
}

if (!is.null(A)) {
    say("")
    say("MODULE % -> COGNITION, +BTC block, standardised beta (BH)")
    z <- A[A$block == "BTC", ]
    for (i in seq_len(nrow(z)))
        say(sprintf("  %-10s %-8s %-16s %-5s beta %+.3f  BH %.4f  n=%d",
                    z$group[i], z$dataset[i], z$rung[i], z$outcome[i],
                    z$beta[i], z$FDR[i], z$n[i]))
}

# =============================================================================
# STEP 3a: the FULL outcome panel, clinical and negative control
# =============================================================================
# Every rung against every outcome, module only, non-module only and both in
# one model. BH within group x dataset x rung x block across outcomes.
OUT_ALL <- c(OUT_CLIN, OUT_NEG)
OUT_ALL <- OUT_ALL[!is.na(OUT_ALL)]
full <- list()
for (grp in unique(DP$group)) {
    for (ds in unique(DP$dataset[DP$group == grp])) {
        pf <- prop_frames[[paste(grp, ds)]]
        ct <- if (!is.null(pf)) pf else
              if (!is.null(clin[[grp]])) clin[[grp]] else CT
        if (!cd %in% names(ct)) next
        names(ct)[names(ct) == cd] <- "donor"
        ct$donor <- as.character(ct$donor)
        have <- OUT_ALL[OUT_ALL %in% names(ct)]
        if (!length(have)) next
        for (rgx in RUNGS) {
            dp <- DP[DP$group == grp & DP$dataset == ds & DP$rung == rgx, ]
            mg <- merge(dp, ct, by = "donor")
            for (lab in names(have)) {
                yv <- have[[lab]]
                blocks <- if (lab %in% COG_LABS) c("cov4", "BTC") else "cov4"
                for (blk in blocks) {
                    cvs <- if (identical(blk, "cov4")) unlist(CV4)
                           else c(unlist(CV4), unlist(BTC))
                    cvs <- drop_self(cvs, yv)
                    # A non-module-only model gives the non-module's own
                    # delta R2.
                    r <- fit_row(mg, yv, "pct_mod_post", cvs)
                    nn <- fit_row(mg, yv, "pct_nm_post", cvs)
                    e <- fit_row(mg, yv, c("pct_mod_post", "pct_nm_post"), cvs)
                    for (tag in c("module_only", "nonmodule_only",
                                  "encompassing")) {
                        z <- switch(tag, module_only = r,
                                    nonmodule_only = nn, encompassing = e)
                        if (is.null(z)) next
                        z$group <- grp; z$dataset <- ds; z$rung <- rgx
                        z$outcome <- lab; z$block <- blk; z$model <- tag
                        z$panel <- if (lab %in% names(OUT_NEG))
                            "negative control" else if (lab %in% COG_LABS)
                            "cognition" else "pathology"
                        full[[length(full) + 1L]] <- z
                    }
                }
            }
        }
    }
}
FA <- do.call(rbind, full)
if (!is.null(FA)) {
    FA$FDR <- NA_real_
    key <- paste(FA$group, FA$dataset, FA$rung, FA$block, FA$model,
                 FA$predictor)
    # Clinical outcomes and negative controls in separate families.
    if (isTRUE(BH_SPLIT_NEG))
        key <- paste(key, FA$panel == "negative control")
    for (k in unique(key)) {
        i <- key == k
        FA$FDR[i] <- stats::p.adjust(FA$p_value[i], "BH")
    }
    FA$significant <- !is.na(FA$FDR) & FA$FDR < ALPHA
    FA$bh_family <- if (isTRUE(BH_SPLIT_NEG))
        paste("BH within group x dataset x rung x block x model x term,",
              "clinical and negative-control outcomes separately")
        else "BH within group x dataset x rung x block x model x term"
}

# == Scorecard ================================================================
#   C1 clinical outcomes the module hits                       (higher better)
#   C2 negative controls the module hits                       (lower better)
#   C3 clinical outcomes the non-module also hits, encompassing (lower better)
#   C4 cognition rows significant at +BTC                      (higher better)
# HKG8 is scored separately.
score <- NULL
if (!is.null(FA)) {
    mk_score <- function(df, tag) {
        do.call(rbind, lapply(RUNGS, function(rgx) {
            a <- df[df$rung == rgx, ]
            mo <- a$model == "module_only" & a$predictor == "pct_mod_post"
            en <- a$model == "encompassing"
            data.frame(
                set = tag, rung = rgx,
                C1_clinical_hits = sum(mo & a$block == "cov4" &
                                       a$panel != "negative control" &
                                       a$significant),
                C1_of = sum(mo & a$block == "cov4" &
                            a$panel != "negative control"),
                C2_control_hits = sum(mo & a$panel == "negative control" &
                                      a$significant),
                C2_of = sum(mo & a$panel == "negative control"),
                C3_nonmodule_hits = sum(en & a$predictor == "pct_nm_post" &
                                        a$block == "cov4" &
                                        a$panel != "negative control" &
                                        a$significant),
                C3_of = sum(en & a$predictor == "pct_nm_post" &
                            a$block == "cov4" &
                            a$panel != "negative control"),
                C4_cog_BTC = sum(mo & a$block == "BTC" &
                                 a$panel == "cognition" & a$significant),
                C4_of = sum(mo & a$block == "BTC" & a$panel == "cognition"),
                median_abs_beta_nonmodule = stats::median(abs(
                    a$beta[en & a$predictor == "pct_nm_post" &
                           a$panel != "negative control"]), na.rm = TRUE),
                stringsAsFactors = FALSE)
        }))
    }
    score <- rbind(mk_score(FA[FA$group != "HKG8", ], "three real groups"),
                   mk_score(FA[FA$group == "HKG8", ], "HKG8 control"))
    say("")
    say("SCORECARD (C1 and C4 higher is better; C2 and C3 lower is better)")
    for (i in seq_len(nrow(score)))
        say(sprintf(paste0("  %-18s %-16s C1 %2d/%2d  C2 %2d/%2d  C3 %2d/%2d",
                           "  C4 %2d/%2d  |beta| nm %.3f"),
                    score$set[i], score$rung[i], score$C1_clinical_hits[i],
                    score$C1_of[i], score$C2_control_hits[i], score$C2_of[i],
                    score$C3_nonmodule_hits[i], score$C3_of[i],
                    score$C4_cog_BTC[i],
                    score$C4_of[i], score$median_abs_beta_nonmodule[i]))
}

# =============================================================================
# STEP 3b: do the four rungs give different donor percentages?
# =============================================================================
# Same donors across rungs: repeated-measures ANOVA if Shapiro-Wilk passes for
# every rung, else Friedman; paired t or Wilcoxon signed-rank post hoc (by
# Shapiro-Wilk on the differences), BH across all post hoc rows.
sw_ok <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) < 3 || length(x) > 5000) return(NA)
    p <- try(stats::shapiro.test(x)$p.value, silent = TRUE)
    if (inherits(p, "try-error")) NA else p >= 0.05
}
sw_p <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) < 3 || length(x) > 5000) return(NA_real_)
    p <- try(stats::shapiro.test(x)$p.value, silent = TRUE)
    if (inherits(p, "try-error")) NA_real_ else p
}

omni <- list(); post <- list()
for (grp in unique(DP$group)) for (ds in unique(DP$dataset[DP$group == grp])) {
    z <- DP[DP$group == grp & DP$dataset == ds, ]
    w <- stats::reshape(z[, c("donor", "rung", "pct_mod_post")],
                        idvar = "donor", timevar = "rung", direction = "wide")
    names(w) <- sub("^pct_mod_post\\.", "", names(w))
    rg <- intersect(RUNGS, names(w))
    w <- w[stats::complete.cases(w[, rg, drop = FALSE]), , drop = FALSE]
    if (nrow(w) < 5 || length(rg) < 2) next
    swp <- vapply(rg, function(r) sw_p(w[[r]]), numeric(1))
    normal <- all(!is.na(swp) & swp >= 0.05)
    M <- as.matrix(w[, rg, drop = FALSE])
    if (normal) {
        lg <- data.frame(y = as.numeric(M),
                         rung = factor(rep(rg, each = nrow(M))),
                         donor = factor(rep(w$donor, times = length(rg))))
        a <- try(stats::anova(stats::lm(y ~ rung + donor, lg)), silent = TRUE)
        op <- if (inherits(a, "try-error")) NA_real_ else a$`Pr(>F)`[1]
        otest <- "repeated-measures ANOVA"
    } else {
        fr <- try(stats::friedman.test(M), silent = TRUE)
        op <- if (inherits(fr, "try-error")) NA_real_ else fr$p.value
        otest <- "Friedman"
    }
    omni[[length(omni) + 1L]] <- data.frame(
        group = grp, dataset = ds, n_donors = nrow(w),
        shapiro_min_p = min(swp, na.rm = TRUE),
        all_normal = normal, omnibus_test = otest, omnibus_p = op,
        stringsAsFactors = FALSE)
    for (i in seq_along(rg)) for (j in seq_along(rg)) {
        if (j <= i) next
        dif <- M[, i] - M[, j]
        pn <- sw_ok(dif)
        use <- if (isTRUE(pn)) "paired t" else "Wilcoxon signed-rank"
        pv <- if (isTRUE(pn))
            try(stats::t.test(M[, i], M[, j], paired = TRUE)$p.value,
                silent = TRUE)
        else
            suppressWarnings(try(stats::wilcox.test(M[, i], M[, j],
                                                    paired = TRUE,
                                                    exact = FALSE)$p.value,
                                 silent = TRUE))
        post[[length(post) + 1L]] <- data.frame(
            group = grp, dataset = ds, a = rg[i], b = rg[j],
            n_donors = nrow(M), shapiro_p_diff = sw_p(dif),
            post_hoc_test = use,
            median_diff = stats::median(dif, na.rm = TRUE),
            p_value = if (inherits(pv, "try-error")) NA_real_ else
                as.numeric(pv),
            stringsAsFactors = FALSE)
    }
}
OM <- do.call(rbind, omni); PH <- do.call(rbind, post)
if (!is.null(OM)) OM$FDR <- stats::p.adjust(OM$omnibus_p, "BH")
if (!is.null(PH)) {
    PH$FDR <- stats::p.adjust(PH$p_value, "BH")
    PH$significant <- !is.na(PH$FDR) & PH$FDR < ALPHA
    say("")
    say("RUNG COMPARISON, donor module % (FDR < ", ALPHA, " significant)")
    for (i in seq_len(nrow(OM)))
        say(sprintf("  %-10s %-8s %-24s p %.4g  FDR %.4g  (Shapiro min %.3g)",
                    OM$group[i], OM$dataset[i], OM$omnibus_test[i],
                    OM$omnibus_p[i], OM$FDR[i], OM$shapiro_min_p[i]))
    say(sprintf("  post hoc: %d of %d pair(s) at FDR < %.2f",
                sum(PH$significant, na.rm = TRUE), nrow(PH), ALPHA))
}

# =============================================================================
# STEP 4: figures
# =============================================================================
pts <- cells[[paste(FIG_GROUP, FIG_DS)]]
if (length(fitcurve) && !is.null(pts)) {
    FC <- fitcurve[[1]]
    FC$rung <- factor(FC$rung, levels = RUNGS)
    p1 <- ggplot() +
        geom_point(data = pts, aes(x = u_all, y = v_all, colour = A3_CMRF),
                   size = 0.25, alpha = 0.5) +
        geom_ribbon(data = FC, aes(x = u, ymin = lo, ymax = hi),
                    fill = NA, colour = CORRIDOR_COL, linewidth = 0.3) +
        geom_line(data = FC, aes(x = u, y = fit), linewidth = 0.3) +
        scale_colour_manual(values = c(`TRUE` = "#05335C",
                                       `FALSE` = "#BEBEBE"),
                            labels = c(`FALSE` = "non-module",
                                       `TRUE`  = "module (A3)"),
                            name = NULL) +
        facet_wrap(~ rung, nrow = 1) +
        labs(x = "u", y = "v",
             title = paste0("Band shapes, ", FIG_GROUP, " ", FIG_DS),
             subtitle = paste("Red = each corridor. A0 to A2 are",
                              "fitted here; A3 is the observed envelope of the",
                              "CMRF classification read from Validation 1.")) +
        theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank())
    b <- file.path(OUT_DIR, "Part11_BandAblation_Shapes")
    pdf(paste0(b, ".pdf"), width = 14, height = 4); print(p1); dev.off()
    png(paste0(b, ".png"), width = 4200, height = 1200, res = 300)
    print(p1); dev.off()
    show_plot(p1)
}

# == Figure: delta R2 by outcome and rung =====================================
if (!is.null(FA)) {
    d4 <- FA[FA$model == "module_only" & FA$block == "cov4" &
             FA$group == FIG_GROUP & FA$dataset == FIG_DS, , drop = FALSE]
    d4n <- FA[FA$model == "nonmodule_only" & FA$block == "cov4" &
              FA$group == FIG_GROUP & FA$dataset == FIG_DS, , drop = FALSE]
    if (nrow(d4) > 0) {
        d4$cells <- "Module Post-IP"
        d4n$cells <- "Non-module Post-IP"
        d4b <- rbind(d4, d4n)
        d4b$outcome <- factor(d4b$outcome,
                              levels = c(names(OUT_CLIN), names(OUT_NEG)))
        d4b$panel <- factor(d4b$panel,
                            levels = c("cognition", "pathology",
                                       "negative control"))
        d4b$star <- ifelse(!d4b$significant, "",
                    ifelse(d4b$FDR < 0.001, "***",
                    ifelse(d4b$FDR < 0.01, "**", "*")))
        p4 <- ggplot(d4b, aes(x = outcome, y = delta_r2, fill = cells)) +
            geom_col(position = position_dodge(width = 0.8), width = 0.7) +
            geom_text(aes(label = star),
                      position = position_dodge(width = 0.8),
                      vjust = -0.2, size = 2.2) +
            scale_fill_manual(values = c(`Module Post-IP` = "#2D2D2D",
                                         `Non-module Post-IP` = "#C8A165"),
                              name = NULL) +
            facet_grid(rung ~ panel, scales = "free_x", space = "free_x") +
            labs(x = NULL, y = "delta R2 over cov4",
                 title = paste0("Every outcome, every corridor: ",
                                FIG_GROUP, " ", FIG_DS),
                 subtitle = paste("FDR < 0.05 *, < 0.01 **, < 0.001 ***.",
                                  "The negative controls are the panel that",
                                  "separates the corridors.")) +
            theme_bw(base_size = 8) +
            theme(panel.grid.minor = element_blank(),
                  strip.text.y = element_text(angle = 0),
                  axis.text.x = element_text(angle = 40, hjust = 1))
        bf <- file.path(OUT_DIR, "Part11_BandAblation_FullPanel")
        pdf(paste0(bf, ".pdf"), width = 12, height = 8); print(p4); dev.off()
        png(paste0(bf, ".png"), width = 3600, height = 2400, res = 300)
        print(p4); dev.off()
        show_plot(p4)
    }
}

if (!is.null(A)) {
    # Labels PVM_3 rather than PVM_3_all; facet order M4, PVM_3, PVM_1, HKG8.
    A$lab <- paste(sub("_all$", "", A$group), A$dataset)
    lab_order <- as.vector(t(outer(GRP_ORDER, c("NeoC", "NonNeoC"), paste)))
    A2f <- A[A$block == "BTC", ]
    # NA FDR counts as not significant.
    A2f$sig <- !is.na(A2f$FDR) & A2f$FDR < ALPHA
    A2f$lab <- factor(A2f$lab,
                      levels = intersect(lab_order, unique(A2f$lab)))
    A2f$rung <- factor(A2f$rung, levels = RUNGS,
                       labels = gsub("_", " ", RUNGS))
    p2 <- ggplot(A2f, aes(x = beta, y = rung, colour = sig)) +
        geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey35") +
        geom_segment(aes(x = 0, xend = beta, y = rung, yend = rung),
                     linewidth = 0.7, lineend = "round") +
        geom_point(size = 2.2) +
        scale_colour_manual(values = c(`TRUE` = SIG_COL,
                                       `FALSE` = "#C4C4C4"),
                            labels = c(`FALSE` = "ns",
                                       `TRUE`  = paste0("FDR < ", ALPHA)),
                            name = NULL) +
        scale_y_discrete(expand = ggplot2::expansion(add = 0.5)) +
        facet_grid(outcome ~ lab) +
        labs(x = "standardised beta, module % -> outcome (+BTC)", y = NULL,
             title = "Does the boundary definition change the association?",
             subtitle = paste("Four corridors, one downstream model,",
                              "covariates plus Braak, Thal and CERAD.")) +
        theme_bw(base_size = 8) +
        theme(panel.grid.minor = element_blank(),
              panel.grid.major.y = element_line(linewidth = 0.2,
                                                colour = "grey93"),
              panel.spacing.x = grid::unit(2.5, "pt"),
              panel.spacing.y = grid::unit(2.5, "pt"),
              strip.background = element_rect(fill = "grey96",
                                              colour = "grey70"),
              strip.text = element_text(size = 7.5),
              strip.text.y = element_text(angle = 0),
              axis.ticks.y = element_blank(),
              legend.position = "bottom",
              legend.margin = margin(t = -4),
              plot.title = element_text(face = "bold", size = 10))
    b <- file.path(OUT_DIR, "Part11_BandAblation_Assoc")
    pdf(paste0(b, ".pdf"), width = 12, height = 3.9); print(p2); dev.off()
    png(paste0(b, ".png"), width = 3600, height = 1170, res = 300)
    print(p2); dev.off()
    show_plot(p2)
}

sheets <- list(Overlap = OVL, Donor_Assoc = A, Encompassing = E,
               Full_Panel = FA, Scorecard = score,
               Rung_Omnibus = OM, Rung_PostHoc = PH,
               Clinical_Route = do.call(rbind, clin_route),
               Coverage = do.call(rbind, cover_tab), Donor_Pct = DP,
               Denominator_Check = DENOM_CHECK)
sheets <- sheets[!vapply(sheets, is.null, logical(1))]
write_xlsx(sheets, file.path(OUT_DIR, "Part11_Identity_BandAblation.xlsx"))
saveRDS(list(sheets = sheets, schema = "M11BA v1",
             generated = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")),
        file.path(OUT_DIR, "Part11_Identity_BandAblation_results.rds"))

say("")
say(strrep("=", 76))
say("Output: ", OUT_DIR)
for (f in sort(list.files(OUT_DIR)))
    say(sprintf("  %-52s %6.0f KB", f,
                file.size(file.path(OUT_DIR, f)) / 1e3))
say("COMPLETE ", format(Sys.time()))
say(strrep("=", 76))
