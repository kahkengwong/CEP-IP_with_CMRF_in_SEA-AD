# =============================================================================
# Part11_BandAblation_Fig11Export_code.r
# SEA-AD CEP-IP/CMRF, Part 11, boundary lines and cell membership (Figure 11B, Table S12)
#   Repeats step 1 of Part11_Identity_BandAblation_code.r with the same code,
#   inputs and settings, and writes the cell-level boundary membership and the
#   Figure 11B boundary curves. The recomputed counts, IP_y, coverage, widths
#   and Jaccard values are checked against the band-ablation workbook.
#
# Usage (from the SEA-AD folder, or set SEAAD_ROOT), after
#   Part11_Identity_BandAblation_code.r
#   source("Part11_BandAblation_Fig11Export_code.r")
#
# Outputs (Part11_Identity/Part11_BandAblation_Fig11Export/): cells, curves and
#   set-summary CSVs, log.
# =============================================================================

library(readxl)
library(mgcv)

CEP_HOME    <- Sys.getenv("SEAAD_ROOT", getwd())
V1_ROOT     <- file.path(CEP_HOME, "Part10e_Validation1")
EXPORT_ROOT <- file.path(CEP_HOME, "R_export_CP10K")
RUN_XLSX    <- file.path(CEP_HOME, "Part11_Identity",
                         "Part11_Identity_BandAblation_P10hDenom",
                         "Part11_Identity_BandAblation.xlsx")
OUT_DIR     <- file.path(CEP_HOME, "Part11_Identity",
                         "Part11_BandAblation_Fig11Export")

# Settings copied from Part11_Identity_BandAblation_code.r.
V1DIR  <- c(M4 = "Mic13_PVM3", PVM_3_all = "PVM_3",
            PVM_1_all = "PVM_1", HKG8 = "HKG8")
ARMDIR <- c(M4 = "Mic13_PVM3", PVM_3_all = "PVM3_all",
            PVM_1_all = "PVM1_all", HKG8 = "HKG8_matched")
STYPE  <- c(M4 = "Micro-PVM_3-SEAAD", PVM_3_all = "Micro-PVM_3-SEAAD",
            PVM_1_all = "Micro-PVM_1", HKG8 = "Micro-PVM_3-SEAAD")
REGS <- list(NeoC = c("DFC", "MTG", "AnG", "STG", "V1C", "ITG"),
             NonNeoC = c("HIP", "MEC", "LEC", "FI"))
MODULE_CLASSES <- c("Original + Refined", "Recovered")
COVERAGE_MODE  <- "match_A3"
N_BINS_A2      <- 20L
RUNGS          <- c("A0_OLS", "A1_spline_fixed", "A2_spline_local", "A3_CMRF")
FIG_GROUP      <- "M4"
FIG_DS         <- "NeoC"
TOL            <- 1e-9

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
LOG <- file.path(OUT_DIR, "Part11_BandAblation_Fig11Export_log.txt")
cat("", file = LOG)
say <- function(...) {
    m <- paste0(...)
    cat(m, "\n", sep = "")
    cat(m, "\n", sep = "", file = LOG, append = TRUE)
}
n_fail <- 0L
check <- function(id, ok, msg) {
    if (!isTRUE(ok)) n_fail <<- n_fail + 1L
    say(sprintf("  [%s] %s: %s", if (isTRUE(ok)) "PASS" else "FAIL", id, msg))
}

# == Helpers copied from the band-ablation script ============================
jac <- function(a, b) {
    u <- length(union(a, b))
    if (u == 0L) NA_real_ else length(intersect(a, b)) / u
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

say(strrep("=", 76))
say("Part11_BandAblation_Fig11Export. Started ", format(Sys.time()))
say("R ", R.version$major, ".", R.version$minor, "; mgcv ",
    as.character(utils::packageVersion("mgcv")))
say(strrep("=", 76))

# == Step 1, as in the band-ablation script ==================================
cells <- list(); ovl <- list(); cover_tab <- list(); fitcurve <- NULL
for (grp in names(V1DIR)) {
    for (ds in c("NeoC", "NonNeoC")) {
        dpath <- file.path(V1_ROOT, paste0("Part10e_Validation1_",
                                           V1DIR[[grp]]))
        f <- list.files(dpath, full.names = TRUE,
                        pattern = paste0("_V1_", ds,
                                         "_SK_interaction_analysis\\.xlsx$"))
        if (!length(f)) { say("  missing: ", grp, " ", ds); next }
        d <- as.data.frame(read_excel(f[1], "Step4_Refined_Cells"))
        need <- c("barcode", "cell_class", "u_all", "v_all", "cDFG_expr")
        miss <- setdiff(need, names(d))
        if (length(miss)) {
            say("  ", grp, " ", ds, " lacks: ", paste(miss, collapse = ", "))
            next
        }
        d <- d[stats::complete.cases(d[, c("u_all", "v_all", "cDFG_expr")]) &
               !is.na(d$cell_class), ]
        lab <- ifelse(grepl("Post-IP", d$group_refined, fixed = TRUE), "Post",
               ifelse(grepl("Pre-IP",  d$group_refined, fixed = TRUE), "Pre",
               ifelse(grepl("Post-IP", d$IP_status, fixed = TRUE), "Post",
               ifelse(grepl("Pre-IP",  d$IP_status, fixed = TRUE), "Pre",
                      NA_character_))))
        hi <- suppressWarnings(max(d$cDFG_expr[!is.na(lab) & lab == "Pre"]))
        lo <- suppressWarnings(min(d$cDFG_expr[!is.na(lab) & lab == "Post"]))
        ip_y <- if (is.finite(hi) && is.finite(lo) && hi < lo) (hi + lo) / 2
                else stats::median(d$cDFG_expr)
        ip_rule <- if (is.finite(hi) && is.finite(lo) && hi < lo)
                       "midpoint of highest Pre-IP and lowest Post-IP cDFG"
                   else "median cDFG (fallback)"
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
        cells[[paste(grp, ds)]] <- data.frame(
            group = grp, dataset = ds, barcode = d$barcode, donor = d$donor,
            cell_class = d$cell_class, cDFG = d$cDFG_expr, side = d$side,
            u = u, v = v,
            fit_OLS = as.numeric(stats::fitted(m0)),
            fit_spline = as.numeric(stats::fitted(m1)),
            width_local = w2,
            in_OLS = d$A0_OLS, in_fixed_spline = d$A1_spline_fixed,
            in_local_spline = d$A2_spline_local, in_CMRF = d$A3_CMRF,
            stringsAsFactors = FALSE)
        cover_tab[[length(cover_tab) + 1L]] <- data.frame(
            group = grp, dataset = ds, n_cells = nrow(d), ip_y = ip_y,
            ip_y_rule = ip_rule, target_coverage = cover, width_A0 = w0,
            width_A1 = w1, width_A2_min = min(w2), width_A2_max = max(w2),
            n_A0 = sum(d$A0_OLS), n_A1 = sum(d$A1_spline_fixed),
            n_A2 = sum(d$A2_spline_local), n_A3 = sum(d$A3_CMRF),
            n_no_donor = sum(is.na(d$donor)), stringsAsFactors = FALSE)
        for (i in seq_along(RUNGS)) for (j in seq_along(RUNGS)) {
            if (j <= i) next
            ovl[[length(ovl) + 1L]] <- data.frame(
                group = grp, dataset = ds, a = RUNGS[i], b = RUNGS[j],
                jaccard = jac(d$barcode[d[[RUNGS[i]]]],
                              d$barcode[d[[RUNGS[j]]]]),
                stringsAsFactors = FALSE)
        }
        if (identical(grp, FIG_GROUP) && identical(ds, FIG_DS)) {
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
            fitcurve <- data.frame(group = grp, dataset = ds, fc,
                                   stringsAsFactors = FALSE)
        }
        say(sprintf(paste0("  %-10s %-8s n=%5d IP_y=%.6f cover=%.4f |",
                           " A0 %5d A1 %5d A2 %5d A3 %5d"),
                    grp, ds, nrow(d), ip_y, cover, sum(d$A0_OLS),
                    sum(d$A1_spline_fixed), sum(d$A2_spline_local),
                    sum(d$A3_CMRF)))
    }
}
CELLS <- do.call(rbind, cells); rownames(CELLS) <- NULL
COV   <- do.call(rbind, cover_tab)
OVL   <- do.call(rbind, ovl)

# == Checks against the band-ablation workbook ===============================
say("")
say("Checks against ", RUN_XLSX)
check("X0", file.exists(RUN_XLSX), "run workbook found")
if (file.exists(RUN_XLSX)) {
    rc <- as.data.frame(read_excel(RUN_XLSX, "Coverage"))
    ro <- as.data.frame(read_excel(RUN_XLSX, "Overlap"))
    m <- merge(COV, rc, by = c("group", "dataset"), suffixes = c("", ".run"))
    check("X1", nrow(m) == 8L && nrow(COV) == 8L,
          sprintf("%d of 8 sets recomputed and matched to the run", nrow(m)))
    for (v in c("n_cells", "n_A0", "n_A1", "n_A2", "n_A3"))
        check("X2", all(m[[v]] == m[[paste0(v, ".run")]]),
              sprintf("%s equal in all sets", v))
    for (v in c("ip_y", "target_coverage", "width_A0", "width_A1",
                "width_A2_min", "width_A2_max")) {
        dv <- max(abs(m[[v]] - m[[paste0(v, ".run")]]))
        check("X3", dv < TOL, sprintf("%s max |diff| %.2e", v, dv))
    }
    mo <- merge(OVL, ro, by = c("group", "dataset", "a", "b"),
                suffixes = c("", ".run"))
    dj <- max(abs(mo$jaccard - mo$jaccard.run))
    check("X4", nrow(mo) == 48L && dj < TOL,
          sprintf("%d of 48 Jaccard values matched, max |diff| %.2e",
                  nrow(mo), dj))
}
check("X5", !is.null(fitcurve) && all(RUNGS %in% fitcurve$rung),
      sprintf("Figure 11B curves for %s %s: %d rows, rungs %s", FIG_GROUP,
              FIG_DS, if (is.null(fitcurve)) 0L else nrow(fitcurve),
              if (is.null(fitcurve)) "-" else
                  paste(unique(fitcurve$rung), collapse = ", ")))

# == Write ===================================================================
files <- c(cells = "BandAblation_cells.csv",
           curves = "BandAblation_Fig11B_curves.csv",
           summary = "BandAblation_set_summary.csv")
write.csv(CELLS, file.path(OUT_DIR, files[["cells"]]), row.names = FALSE)
if (!is.null(fitcurve))
    write.csv(fitcurve, file.path(OUT_DIR, files[["curves"]]),
              row.names = FALSE)
write.csv(COV, file.path(OUT_DIR, files[["summary"]]), row.names = FALSE)

say("")
say("Written to ", OUT_DIR)
for (f in files) {
    p <- file.path(OUT_DIR, f)
    if (file.exists(p))
        say(sprintf("  %-34s %8.0f KB", f, file.size(p) / 1e3))
}
say(sprintf("  cells: %d rows (%d without a donor)", nrow(CELLS),
            sum(is.na(CELLS$donor))))
say("")
say(sprintf("VERDICT: %s (%d FAIL)", if (n_fail == 0L) "PASS" else "FAIL",
            n_fail))
say("Finished ", format(Sys.time()))
