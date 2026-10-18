# =============================================================================
# Part10e_CP10K_Validation2_UMAP_code.r
# SEA-AD CEP-IP/CMRF, Part 10e, Validation 2: CMRF classes in transcriptomic
#   space (Figure 8B and 8C, Table S10)
#   Builds a Monocle3 UMAP from the Part 8 raw count export, colours it by
#   CMRF class and IP x TREP group, and compares Recovered and Original-only
#   cells with the four Original+Refined subgroups (TOST, PERMANOVA and
#   PERMDISP, iLISI, kNN proportions, kBET) at the cell and donor units. Also
#   computes the Figure 11 statistics: the Pre-IP to Post-IP TREP projection
#   in PCA and copula space and the bridge distance.
#
# Usage (from the SEA-AD folder, or set SEAAD_ROOT), after
#   Part10e_CP10K_Validation1_code.r
#   source("Part10e_CP10K_Validation2_UMAP_code.r")                         # all arms
#   CEP_ARM <- "PVM1_all"; source("Part10e_CP10K_Validation2_UMAP_code.r")  # one arm
#   Needs CEP_arms.R, R_config.R and R_config_HKG8.R in the same folder.
#
# Outputs (Part10e_Validation2/Part10e_Validation2_<arm>/): figures (PDF, PNG),
#   interactive 3D HTML, summary workbook, results RDS, guardrail ledger and
#   logs.
# =============================================================================


# =============================================================================
# SECTION 0: Libraries and shared arm definitions
# =============================================================================

library(Matrix)
library(monocle3)
library(ggplot2)
library(dplyr)
library(readxl)
library(writexl)
library(gridExtra)
library(grid)
library(scales)
library(FNN)

set.seed(42)

# vegan (PERMANOVA, PERMDISP), plotly and htmlwidgets (3D HTML) are optional;
# the matching output is skipped if one is missing.
.v2_vegan_ok <- tryCatch({ library(vegan); TRUE }, error = function(e) {
    cat("NOTE: vegan not installed -- PERMANOVA and PERMDISP will be skipped.\n")
    FALSE })
.v2_plotly_ok <- tryCatch({ library(plotly); TRUE }, error = function(e) {
    cat("NOTE: plotly not installed -- interactive 3D HTML will be skipped.\n")
    FALSE })
.v2_widgets_ok <- tryCatch({ library(htmlwidgets); TRUE }, error = function(e) {
    cat("NOTE: htmlwidgets not installed -- interactive 3D HTML will be skipped.\n")
    FALSE })

if (!exists("CEP_HOME")) CEP_HOME <- Sys.getenv("SEAAD_ROOT", getwd())
.p10v2_harness <- file.path(CEP_HOME, "CEP_arms.R")
if (!file.exists(.p10v2_harness))
    stop("CEP_arms.R not found at: ", .p10v2_harness,
         "\n  It holds the arm table shared by Parts 9 onwards.", call. = FALSE)
source(.p10v2_harness)

PART_TAG     <- "Part10e"              # file-name tag; "_V2" is appended
PART_DIR_TAG <- "Part10e_Validation2"  # results-folder tag

# == Embedding settings =======================================================
# The embedding source is the Part 8 raw count export.
V2_EMBED_SOURCE       <- "counts_raw_selected"
V2_NUM_DIM            <- 50L    # PCA dimensions fed into UMAP
V2_MIN_CELLS_PER_GENE <- 3L     # a gene must be detected in this many cells
V2_SEED               <- 42L

# == Analysis settings ========================================================
V2_KNN_K              <- 50L    # neighbourhood size for the proportion and kBET tests

# Interactive 3D scatter: one marker size for every group.
V2_3D_MARKER_SIZE <- 3.5
# Marker size of the Non-selected and Excluded background.
V2_3D_BG_SIZE     <- V2_3D_MARKER_SIZE
# Marker outline width, all traces.
V2_3D_OUTLINE_W   <- 0.20
# Draw the two 2D UMAP figures (the 2D embedding is always computed).
V2_MAKE_2D_FIGS   <- FALSE

# Point size of the static 3D projection panels, by cell count.
v2_proj_pt <- function(n)
    if (n > 20000) 0.60 else if (n > 10000) 0.80 else
    if (n >  5000) 1.00 else if (n >  1000) 1.40 else 2.00
# Static panel aesthetics, matching Validation1's CEP-IP scatter;
# 1.63 keeps Validation1's 2.2/1.35 ratio for the large markers.
V2_PROJ_ALPHA     <- 0.58    # Validation1's OrigRef circle alpha
V2_PROJ_RECOV_A   <- 0.93    # Validation1's Recovered diamond alpha
V2_PROJ_OO_A      <- 0.88    # Validation1's Original-only triangle alpha
V2_PROJ_STROKE    <- 0.15    # Validation1's white border width, both shapes
V2_PROJ_BIG_MULT  <- 1.63
# Groups drawn in the three-projection panel; the by-group figure shows all six.
V2_PROJ_GROUPS    <- c("Pre-IP TREP", "Post-IP TREP")
V2_PROJ_PDF_W     <- 21
V2_PROJ_PDF_H     <- 8
V2_PROJ_PNG_W     <- 4200
V2_PROJ_PNG_H     <- 1600
V2_PROJ_PNG_RES   <- 240
# Extra figure: one U1-U2 panel per group over a grey backdrop.
V2_PROJ_BY_GROUP  <- TRUE
V2_ILISI_PERPLEXITY   <- 30L    # iLISI Gaussian kernel perplexity; k = 3 x perplexity
V2_KBET_ALPHA         <- 0.05
V2_KBET_SAMPLE        <- 500L   # cells tested per pair; NULL tests all
V2_PERM_N             <- 999L   # PERMANOVA and PERMDISP permutations
V2_PERM_MAX_N         <- 500L   # cap per group in a PERMANOVA pair
V2_TOST_DELTA_FRAC    <- 0.5    # equivalence margin in pooled SD
V2_MIN_GROUP          <- 5L     # a group smaller than this is not tested
V2_ORIGREF_GROUP_COL  <- "group_refined"   # "group_orig" for the Part 10e groups
V2_WRITE_HTML         <- TRUE

# == Donor handling ===========================================================
# TOST and PERMANOVA are also run with each donor as one observation
# (the inverse-donor weighting).
V2_DONOR_UNIT_TESTS    <- TRUE
# == Figure 11 statistics =====================================================
# Copula space: spine residual and SK-plane projection (method checks, set by
# the selection rule). Transcriptome: position on the Pre-IP TREP -> Post-IP
# TREP axis in the 50-PC PCA (0 = Pre centroid, 1 = Post centroid).
# Tests: Shapiro-Wilk; two groups t-test or Mann-Whitney U; more than two
# ANOVA with Holm-Sidak or Kruskal-Wallis with Dunn. Donor unit (repeated
# measures on complete blocks): repeated-measures ANOVA or Friedman, Page's L
# for the trend, paired tests post hoc. BH: m = 3 pairwise within dataset x
# unit x statistic; m = 2 across datasets for omnibus and trend.
V2_STAT_ALPHA      <- 0.05  # normality and significance cut
V2_SW_MAX_N        <- 5000L # shapiro.test aborts above 5000
V2_STAT_MIN_GROUP  <- 5L    # cells or donors a group needs to be tested
V2_STAT_MIN_BLOCKS <- 5L    # complete donors the blocked tests need
# Group order, low to high, for the trend tests.
V2_PROJ_LEVELS     <- c("Pre-IP TREP", "Recovered", "Post-IP TREP")
V2_RESID_LEVELS    <- c("Original + Refined", "Recovered",
                        "Original only", "Excluded")
V2_RESID_COL       <- "resid_abs"   # Validation1's |v - g(u)|
# Run the Pre-IP -> Post-IP projection in the copula plane too (method check).
V2_PROJ_COPULA     <- TRUE
# The copula plane is the SKCS one, (u_SK, v_SK) = (F(Si_norm), F(Ki)), read
# from Validation1's RDS and checked against a recompute from the Step 4 sheet.
V2_COPULA_SK_FIELDS <- c("copula_u_sk", "copula_v_sk")
V2_COPULA_SK_RECALC <- c("spear_norm", "kend_contrib")
V2_COPULA_SK_TOL    <- 1e-8
# Spine residual in the summary table (it has no counterpart on the SK axes).
V2_RESID_IN_SUMMARY <- FALSE
# Bridge distance: per cell, the larger of the distances to the nearest Pre-IP
# TREP and nearest Post-IP TREP cell, on the SK axes. Recovered is compared
# with Excluded.
V2_BRIDGE_ENABLE  <- TRUE
V2_BRIDGE_LEVELS  <- c("Pre-IP TREP", "Recovered", "Post-IP TREP", "Excluded")
# SKCS cuts, to record which cut each Recovered cell fails (same as Part 10e).
V2_SI_NORM_CUT    <- 0.05
V2_KI_CUT         <- 0.5

V2_MIN_DONORS          <- 3L   # donors a group needs before it is tested
V2_MIN_CELLS_PER_DONOR <- 3L   # cells a donor needs to contribute a mean
V2_TVD_PERM            <- 200L # label permutations for the TVD null

# == The four CMRF classes ====================================================
CMRF_CLASS_LEVS <- c("Original + Refined", "Recovered",
                     "Original only",      "Excluded")

# Part 10e CP10K class colours.
VP_COL_ORIGREF <- "#009E73"
VP_COL_RECOV   <- "#7C3AED"
VP_COL_ORIONLY <- "#A0522D"
VP_COL_EXCL    <- "#8A9BAD"

CMRF_CLASS_COLS <- c(
    "Original + Refined" = VP_COL_ORIGREF,
    "Recovered"          = VP_COL_RECOV,
    "Original only"      = VP_COL_ORIONLY,
    "Excluded"           = VP_COL_EXCL
)

# The four IP x TREP groups, identical to Validation1's four_cols.
FOUR_COLS <- c(
    "Pre-IP TREP"      = "#E69F00",
    "Pre-IP non-TREP"  = "#1B3A8C",
    "Post-IP TREP"     = "#009E73",
    "Post-IP non-TREP" = "#D94F4F"
)
# Levels of the mapped colour scale.
GROUP5_LEVS    <- c("Non-selected", names(FOUR_COLS), "Unknown subgroup")
SCALE_COLS_VAL <- c("Non-selected" = "gray60", FOUR_COLS,
                    "Unknown subgroup" = "#7A7A7A")

RECOV_COL_04 <- VP_COL_RECOV     # Recovered diamond, shape 23
ORIGONLY_COL <- "#8B6914"        # Original-only inverted triangle, shape 25


# =============================================================================
# SECTION 0d: Pure helpers
# =============================================================================
# These depend only on their arguments.

# Coerce a logical, numeric or character flag column to logical.
as_flag <- function(x) {
    if (is.logical(x)) return(x)
    if (is.numeric(x)) return(x != 0)
    y <- toupper(trimws(as.character(x)))
    out <- rep(NA, length(y))
    out[y %in% c("TRUE", "T", "1", "YES")]  <- TRUE
    out[y %in% c("FALSE", "F", "0", "NO")]  <- FALSE
    out
}

theme_v2 <- function(base_size = 10) {
    theme_minimal(base_size = base_size) +
        theme(
            plot.title       = element_text(face = "bold", size = base_size + 1),
            plot.subtitle    = element_text(size = base_size - 1, colour = "gray30",
                                            lineheight = 1.25),
            axis.title       = element_text(face = "bold", size = base_size - 1),
            axis.text        = element_text(size = base_size - 2),
            panel.grid.major = element_line(colour = "gray90"),
            panel.grid.minor = element_blank(),
            panel.border     = element_rect(colour = "black", fill = NA,
                                            linewidth = 0.35),
            plot.background  = element_rect(fill = "white", colour = NA),
            panel.background = element_rect(fill = "white", colour = NA),
            legend.position  = "right"
        )
}

show_plot <- function(p) {
    tryCatch(print(p), error = function(e) invisible(NULL))
}

# Draw a grob in the Plots pane.
show_grob <- function(g) {
    tryCatch({
        grid::grid.newpage()
        grid::grid.draw(g)
    }, error = function(e) invisible(NULL))
}

save_pp <- function(p, base, pdf_w = 10, pdf_h = 9,
                    png_w = 2200, png_h = 2000, res = 200) {
    tryCatch({
        pdf(paste0(base, ".pdf"), width = pdf_w, height = pdf_h)
        print(p); dev.off()
    }, error = function(e) {
        tryCatch(dev.off(), error = function(e2) invisible(NULL))
        cat("  [save_pp] PDF error:", conditionMessage(e), "\n")
    })
    tryCatch({
        png(paste0(base, ".png"), width = png_w, height = png_h, res = res)
        print(p); dev.off()
    }, error = function(e) {
        tryCatch(dev.off(), error = function(e2) invisible(NULL))
        cat("  [save_pp] PNG error:", conditionMessage(e), "\n")
    })
}

save_grob <- function(g, base, pdf_w = 18, pdf_h = 7,
                      png_w = 3800, png_h = 1500, res = 200) {
    tryCatch({
        pdf(paste0(base, ".pdf"), width = pdf_w, height = pdf_h)
        grid.newpage(); grid.draw(g); dev.off()
    }, error = function(e) {
        tryCatch(dev.off(), error = function(e2) invisible(NULL))
        cat("  [save_grob] PDF error:", conditionMessage(e), "\n")
    })
    tryCatch({
        png(paste0(base, ".png"), width = png_w, height = png_h, res = res)
        grid.newpage(); grid.draw(g); dev.off()
    }, error = function(e) {
        tryCatch(dev.off(), error = function(e2) invisible(NULL))
        cat("  [save_grob] PNG error:", conditionMessage(e), "\n")
    })
}


# == v2_build_cds: Monocle3 CDS from a genes x cells matrix, PCA run once =====
v2_build_cds <- function(count_mat, class_df, num_dim = V2_NUM_DIM,
                         seed = V2_SEED) {

    cat(sprintf("    CDS: %d gene(s) x %d cell(s)\n",
                nrow(count_mat), ncol(count_mat)))

    cell_meta <- data.frame(
        cell_id    = colnames(count_mat),
        cell_class = as.character(class_df[colnames(count_mat), "cell_class"]),
        stringsAsFactors = FALSE)
    # Every cell must have a class.
    if (anyNA(cell_meta$cell_class))
        stop(sum(is.na(cell_meta$cell_class)),
             " cell(s) in the count matrix have no cell_class in the label ",
             "frame. The two are keyed on barcode and must cover the same ",
             "cells.", call. = FALSE)
    rownames(cell_meta) <- cell_meta$cell_id

    gene_meta <- data.frame(
        gene_id         = rownames(count_mat),
        gene_short_name = rownames(count_mat),
        stringsAsFactors = FALSE)
    rownames(gene_meta) <- gene_meta$gene_id

    cds <- new_cell_data_set(
        expression_data = count_mat,
        cell_metadata   = cell_meta,
        gene_metadata   = gene_meta)

    # Cap the PCA dimensions for irlba.
    n_pca <- max(2L, min(as.integer(num_dim),
                         floor(nrow(count_mat) / 2L) - 1L,
                         ncol(count_mat) - 1L))
    set.seed(seed)
    cds <- preprocess_cds(cds, num_dim = n_pca, method = "PCA", verbose = FALSE)
    cat(sprintf("    PCA: %d dimension(s)\n", n_pca))
    list(cds = cds, n_pca = n_pca)
}

# == v2_reduce: add a UMAP of the given dimensionality ========================
v2_reduce <- function(cds, umap_dims = 2L, seed = V2_SEED) {
    set.seed(seed)
    reduce_dimension(cds,
                     reduction_method  = "UMAP",
                     preprocess_method = "PCA",
                     max_components    = as.integer(umap_dims),
                     verbose           = FALSE)
}

# == v2_umap_coords: UMAP coordinate frame from a CDS =========================
v2_umap_coords <- function(cds_obj) {
    rd_names <- reducedDimNames(cds_obj)
    umap_key <- rd_names[grepl("umap", rd_names, ignore.case = TRUE)]
    if (length(umap_key) == 0L) {
        cat("    WARNING: no UMAP reducedDim. Available:",
            paste(rd_names, collapse = ", "), "\n")
        return(NULL)
    }
    umap_mat <- reducedDim(cds_obj, umap_key[1L])
    colnames(umap_mat) <- paste0("UMAP_", seq_len(ncol(umap_mat)))
    df_out <- as.data.frame(umap_mat)

    # Barcodes come from the CDS columns.
    bc <- colnames(cds_obj)
    if (is.null(bc) || length(bc) != nrow(df_out))
        stop(sprintf(paste0("UMAP has %d row(s) but the CDS carries %s cell",
                            " name(s); cannot align barcodes."),
                     nrow(df_out),
                     if (is.null(bc)) "no" else as.character(length(bc))),
             call. = FALSE)
    df_out$barcode    <- bc
    df_out$cell_class <- as.character(colData(cds_obj)$cell_class)
    rownames(df_out)  <- df_out$barcode
    df_out
}

# == v2_donor_means: donor-level means of the UMAP axes =======================
# Equal to the inverse-donor-weighted mean over cells.
v2_donor_means <- function(df, axes, donor_col = "donor",
                           min_cells = V2_MIN_CELLS_PER_DONOR) {
    if (is.null(df) || nrow(df) == 0L)  return(NULL)
    if (!donor_col %in% names(df))      return(NULL)
    d <- as.character(df[[donor_col]])
    keep <- !is.na(d) & nzchar(d)
    if (!any(keep)) return(NULL)
    d  <- d[keep]
    df <- df[keep, , drop = FALSE]

    n_by <- table(d)
    ok   <- names(n_by)[n_by >= as.integer(min_cells)]
    if (length(ok) == 0L) return(NULL)
    sel <- d %in% ok
    d   <- d[sel]
    df  <- df[sel, , drop = FALSE]

    out <- data.frame(donor = sort(unique(d)), stringsAsFactors = FALSE)
    for (ax in axes)
        out[[ax]] <- as.numeric(tapply(as.numeric(df[[ax]]), d, mean)[out$donor])
    out$n_cells <- as.integer(table(d)[out$donor])
    out
}

# == Figure 11 statistics helpers (base R) ====================================

# Console p formatting; NA for a missing or non-finite p.
v2_p_fmt <- function(p) {
    if (is.null(p) || length(p) != 1L || !is.finite(p)) return("NA")
    if (p < 1e-4) "<1e-04" else sprintf("%.4g", p)
}

# Shapiro-Wilk p, subsampled above V2_SW_MAX_N; NA when the test cannot run.
v2_sw <- function(x, max_n = V2_SW_MAX_N, seed = V2_SEED) {
    x <- as.numeric(x); x <- x[is.finite(x)]
    if (length(x) < 3L || length(unique(x)) < 3L) return(NA_real_)
    if (length(x) > max_n) {
        set.seed(seed); x <- x[sample.int(length(x), max_n)]
    }
    tryCatch(stats::shapiro.test(x)$p.value, error = function(e) NA_real_)
}

# Holm-Sidak step-down adjustment.
v2_sidak_holm <- function(p) {
    m <- length(p); if (m == 0L) return(numeric(0))
    o <- order(p); pa <- rep(NA_real_, m); run <- 0
    for (i in seq_len(m)) {
        v   <- 1 - (1 - p[o[i]])^(m - i + 1L)
        run <- max(run, v); pa[o[i]] <- min(1, run)
    }
    pa
}

# Mann-Whitney U count of (a < b), ties at half, from ranks.
.v2_u_lt <- function(a, b) {
    r <- rank(c(a, b))
    sum(r[(length(a) + 1L):(length(a) + length(b))]) -
        length(b) * (length(b) + 1) / 2
}

# Jonckheere-Terpstra test, normal approximation with tie correction; the
# factor level order of g is the hypothesised order. Two-sided.
v2_jt <- function(x, g) {
    ok <- is.finite(x) & !is.na(g)
    x  <- x[ok]; g <- droplevels(factor(g[ok], levels = levels(g)))
    lv <- levels(g); k <- length(lv)
    if (k < 2L) return(list(J = NA_real_, z = NA_real_, p = NA_real_))
    n <- as.numeric(table(g)[lv]); N <- sum(n)
    if (any(n < 2L)) return(list(J = NA_real_, z = NA_real_, p = NA_real_))
    J <- 0
    for (i in seq_len(k - 1L)) for (j in (i + 1L):k)
        J <- J + .v2_u_lt(x[g == lv[i]], x[g == lv[j]])
    mu <- (N^2 - sum(n^2)) / 4
    tb <- as.numeric(table(x))
    A  <- N * (N - 1) * (2 * N + 5) - sum(n * (n - 1) * (2 * n + 5)) -
        sum(tb * (tb - 1) * (2 * tb + 5))
    B  <- sum(n * (n - 1) * (n - 2)) * sum(tb * (tb - 1) * (tb - 2))
    C  <- sum(n * (n - 1)) * sum(tb * (tb - 1))
    v  <- A / 72 + B / (36 * N * (N - 1) * (N - 2)) + C / (8 * N * (N - 1))
    if (!is.finite(v) || v <= 0)
        return(list(J = J, z = NA_real_, p = NA_real_))
    z <- (J - mu) / sqrt(v)
    list(J = J, z = z, p = 2 * stats::pnorm(-abs(z)))
}

# Dunn's post hoc for Kruskal-Wallis, tie-corrected (independent groups).
v2_dunn <- function(x, g) {
    ok <- is.finite(x) & !is.na(g)
    x  <- x[ok]; g <- droplevels(factor(g[ok]))
    lv <- levels(g); k <- length(lv); N <- length(x)
    if (k < 2L) return(NULL)
    r    <- rank(x)
    tb   <- as.numeric(table(r)); ties <- sum(tb^3 - tb)
    sig  <- (N * (N + 1)) / 12 - ties / (12 * (N - 1))
    out  <- NULL
    for (i in seq_len(k - 1L)) for (j in (i + 1L):k) {
        ni <- sum(g == lv[i]); nj <- sum(g == lv[j])
        z  <- (mean(r[g == lv[i]]) - mean(r[g == lv[j]])) /
            sqrt(sig * (1 / ni + 1 / nj))
        out <- rbind(out, data.frame(
            group_a = lv[i], group_b = lv[j], post_hoc = "Dunn",
            stat = z, p_value = 2 * stats::pnorm(-abs(z)),
            n_a = ni, n_b = nj, stringsAsFactors = FALSE))
    }
    out
}

# Page's L on a blocks x treatments matrix, columns in the hypothesised order;
# incomplete blocks are dropped.
v2_page <- function(m) {
    m <- m[stats::complete.cases(m), , drop = FALSE]
    b <- nrow(m); k <- ncol(m)
    if (b < 2L || k < 3L)
        return(list(L = NA_real_, z = NA_real_, p = NA_real_, n_block = b))
    Rj <- colSums(t(apply(m, 1L, rank)))
    L  <- sum(seq_len(k) * Rj)
    mu <- b * k * (k + 1)^2 / 4
    sd <- sqrt(b * k^2 * (k + 1) * (k^2 - 1) / 144)
    if (!is.finite(sd) || sd <= 0)
        return(list(L = L, z = NA_real_, p = NA_real_, n_block = b))
    z <- (L - mu) / sd
    list(L = L, z = z, p = 2 * stats::pnorm(-abs(z)), n_block = b)
}

# One-way repeated-measures ANOVA on complete blocks.
v2_rm_anova <- function(m) {
    m <- m[stats::complete.cases(m), , drop = FALSE]
    b <- nrow(m); k <- ncol(m)
    if (b < 2L || k < 2L)
        return(list(F = NA_real_, p = NA_real_, n_block = b))
    gm   <- mean(m)
    ss_t <- b * sum((colMeans(m) - gm)^2)
    ss_b <- k * sum((rowMeans(m) - gm)^2)
    ss_e <- sum((m - gm)^2) - ss_t - ss_b
    df_t <- k - 1L; df_e <- (b - 1L) * (k - 1L)
    if (df_e < 1L || !is.finite(ss_e) || ss_e <= 0)
        return(list(F = NA_real_, p = NA_real_, n_block = b))
    Fv <- (ss_t / df_t) / (ss_e / df_e)
    list(F = Fv, p = stats::pf(Fv, df_t, df_e, lower.tail = FALSE),
         n_block = b)
}

# Per-group descriptives.
v2_group_desc <- function(x, g, levs) {
    do.call(rbind, lapply(levs, function(l) {
        v <- as.numeric(x[!is.na(g) & g == l]); v <- v[is.finite(v)]
        data.frame(group = l, n = length(v),
                   median = if (length(v)) stats::median(v) else NA_real_,
                   q25 = if (length(v)) unname(stats::quantile(v, .25)) else NA_real_,
                   q75 = if (length(v)) unname(stats::quantile(v, .75)) else NA_real_,
                   mean = if (length(v)) mean(v) else NA_real_,
                   sd   = if (length(v) > 1L) stats::sd(v) else NA_real_,
                   stringsAsFactors = FALSE)
    }))
}

# Test protocol. unit = "cell" (independent groups) or "donor" (repeated
# measures on complete blocks).
v2_protocol <- function(x, g, block = NULL, levs, unit,
                        alpha = V2_STAT_ALPHA) {
    g   <- factor(as.character(g), levels = levs)
    ok  <- is.finite(x) & !is.na(g)
    x   <- x[ok]; g <- g[ok]
    if (!is.null(block)) block <- as.character(block)[ok]
    desc <- v2_group_desc(x, g, levs)
    out  <- list(desc = desc, ph = NULL,
                 omnibus = data.frame(
                     unit = unit, n_group = length(levs),
                     normal = NA, shapiro_min_p = NA_real_,
                     omnibus_test = NA_character_, omnibus_stat = NA_real_,
                     omnibus_p = NA_real_, trend_test = NA_character_,
                     trend_stat = NA_real_, trend_p = NA_real_,
                     n_unit = NA_integer_, note = NA_character_,
                     stringsAsFactors = FALSE))

    if (identical(unit, "donor")) {
        if (is.null(block)) {
            out$omnibus$note <- "no donor column"; return(out)
        }
        # Donor means per group, complete donors only.
        agg <- tapply(x, list(block, g), mean)
        agg <- agg[, levs, drop = FALSE]
        m   <- agg[stats::complete.cases(agg), , drop = FALSE]
        out$omnibus$n_unit <- nrow(m)
        if (nrow(m) < V2_STAT_MIN_BLOCKS) {
            out$omnibus$note <- sprintf(
                "only %d complete donor block(s); minimum %d",
                nrow(m), V2_STAT_MIN_BLOCKS)
            return(out)
        }
        sw <- vapply(seq_len(ncol(m)), function(j) v2_sw(m[, j]), numeric(1))
        normal <- all(is.finite(sw)) && all(sw >= alpha)
        out$omnibus$shapiro_min_p <- suppressWarnings(min(sw, na.rm = TRUE))
        out$omnibus$normal <- normal
        if (normal) {
            a <- v2_rm_anova(m)
            out$omnibus$omnibus_test <- "repeated-measures ANOVA"
            out$omnibus$omnibus_stat <- a$F; out$omnibus$omnibus_p <- a$p
        } else {
            fr <- tryCatch(stats::friedman.test(m),
                           error = function(e) NULL)
            out$omnibus$omnibus_test <- "Friedman"
            if (!is.null(fr)) {
                out$omnibus$omnibus_stat <- unname(fr$statistic)
                out$omnibus$omnibus_p    <- fr$p.value
            }
        }
        pg <- v2_page(m)
        out$omnibus$trend_test <- "Page's L"
        out$omnibus$trend_stat <- pg$z; out$omnibus$trend_p <- pg$p
        # Post hoc: paired tests.
        ph <- NULL
        for (i in seq_len(ncol(m) - 1L)) for (j in (i + 1L):ncol(m)) {
            a <- m[, i]; b <- m[, j]
            if (normal) {
                tt <- tryCatch(stats::t.test(a, b, paired = TRUE),
                               error = function(e) NULL)
                nm <- "paired t-test"
            } else {
                tt <- tryCatch(suppressWarnings(
                    stats::wilcox.test(a, b, paired = TRUE, exact = FALSE)),
                    error = function(e) NULL)
                nm <- "paired Wilcoxon signed-rank"
            }
            ph <- rbind(ph, data.frame(
                unit = unit, group_a = levs[i], group_b = levs[j],
                post_hoc = nm,
                stat = if (is.null(tt)) NA_real_ else unname(tt$statistic),
                p_value = if (is.null(tt)) NA_real_ else tt$p.value,
                n_a = nrow(m), n_b = nrow(m), stringsAsFactors = FALSE))
        }
        if (!is.null(ph)) {
            ph$p_holm_sidak <- if (normal) v2_sidak_holm(ph$p_value) else NA_real_
            ph$p_BH_FDR     <- stats::p.adjust(ph$p_value, "BH")
            ph$bh_family_n  <- nrow(ph)
        }
        out$ph <- ph
        return(out)
    }

    # == cell unit, independent groups ========================================
    keep <- levs[table(g)[levs] >= V2_STAT_MIN_GROUP]
    out$omnibus$n_unit <- length(x)
    if (length(keep) < 2L) {
        out$omnibus$note <- "fewer than two groups reach the minimum size"
        return(out)
    }
    gk <- droplevels(factor(as.character(g), levels = keep))
    xk <- x[!is.na(gk)]; gk <- gk[!is.na(gk)]
    sw <- vapply(levels(gk), function(l) v2_sw(xk[gk == l]), numeric(1))
    normal <- all(is.finite(sw)) && all(sw >= alpha)
    out$omnibus$shapiro_min_p <- suppressWarnings(min(sw, na.rm = TRUE))
    out$omnibus$normal <- normal
    if (nlevels(gk) == 2L) {
        if (normal) {
            tt <- stats::t.test(xk ~ gk)
            out$omnibus$omnibus_test <- "t-test"
        } else {
            tt <- suppressWarnings(stats::wilcox.test(xk ~ gk, exact = FALSE))
            out$omnibus$omnibus_test <- "Mann-Whitney U"
        }
        out$omnibus$omnibus_stat <- unname(tt$statistic)
        out$omnibus$omnibus_p    <- tt$p.value
    } else if (normal) {
        av <- summary(stats::aov(xk ~ gk))[[1]]
        out$omnibus$omnibus_test <- "ANOVA"
        out$omnibus$omnibus_stat <- av[["F value"]][1]
        out$omnibus$omnibus_p    <- av[["Pr(>F)"]][1]
    } else {
        kw <- stats::kruskal.test(xk, gk)
        out$omnibus$omnibus_test <- "Kruskal-Wallis"
        out$omnibus$omnibus_stat <- unname(kw$statistic)
        out$omnibus$omnibus_p    <- kw$p.value
    }
    if (nlevels(gk) >= 3L) {
        jt <- v2_jt(xk, gk)
        out$omnibus$trend_test <- "Jonckheere-Terpstra"
        out$omnibus$trend_stat <- jt$z; out$omnibus$trend_p <- jt$p
    }
    ph <- NULL
    lk <- levels(gk)
    for (i in seq_len(length(lk) - 1L)) for (j in (i + 1L):length(lk)) {
        a <- xk[gk == lk[i]]; b <- xk[gk == lk[j]]
        if (normal) {
            tt <- tryCatch(stats::t.test(a, b), error = function(e) NULL)
            nm <- "t-test (Holm-Sidak)"
        } else {
            tt <- tryCatch(suppressWarnings(
                stats::wilcox.test(a, b, exact = FALSE)),
                error = function(e) NULL)
            nm <- "Dunn"
        }
        ph <- rbind(ph, data.frame(
            unit = unit, group_a = lk[i], group_b = lk[j], post_hoc = nm,
            stat = if (is.null(tt)) NA_real_ else unname(tt$statistic),
            p_value = if (is.null(tt)) NA_real_ else tt$p.value,
            n_a = length(a), n_b = length(b), stringsAsFactors = FALSE))
    }
    if (!normal && nlevels(gk) >= 3L) {
        dn <- v2_dunn(xk, gk)
        if (!is.null(dn)) {
            dn$unit <- unit
            ph <- dn[, c("unit", "group_a", "group_b", "post_hoc",
                         "stat", "p_value", "n_a", "n_b")]
        }
    }
    if (!is.null(ph)) {
        ph$p_holm_sidak <- if (normal) v2_sidak_holm(ph$p_value) else NA_real_
        ph$p_BH_FDR     <- stats::p.adjust(ph$p_value, "BH")
        ph$bh_family_n  <- nrow(ph)
    }
    out$ph <- ph
    out
}

# Distance from each row of X to its nearest row of Y, in blocks.
v2_nn_dist <- function(X, Y, block = 2000L) {
    if (is.null(X) || is.null(Y) || nrow(X) == 0L || nrow(Y) == 0L)
        return(rep(NA_real_, if (is.null(X)) 0L else nrow(X)))
    sy  <- rowSums(Y^2)
    out <- rep(NA_real_, nrow(X))
    for (i in seq(1L, nrow(X), by = block)) {
        j  <- i:min(i + block - 1L, nrow(X))
        Xb <- X[j, , drop = FALSE]
        d2 <- outer(rowSums(Xb^2), sy, "+") - 2 * (Xb %*% t(Y))
        d2[d2 < 0] <- 0                     # rounding can give a tiny negative
        out[j] <- sqrt(apply(d2, 1L, min))
    }
    out
}

# Bridge distance: distance to the further of the two arms. A reference cell
# is measured against the other members of its own arm.
v2_bridge_dist <- function(XY, grp, pre_lab, post_lab) {
    if (is.null(XY)) return(NULL)
    i_pre  <- which(!is.na(grp) & grp == pre_lab)
    i_post <- which(!is.na(grp) & grp == post_lab)
    if (length(i_pre) < 2L || length(i_post) < 2L) return(NULL)
    d_pre  <- rep(NA_real_, nrow(XY)); d_post <- d_pre
    own_pre  <- seq_len(nrow(XY)) %in% i_pre
    own_post <- seq_len(nrow(XY)) %in% i_post
    d_pre[!own_pre]   <- v2_nn_dist(XY[!own_pre, , drop = FALSE],
                                    XY[i_pre, , drop = FALSE])
    d_post[!own_post] <- v2_nn_dist(XY[!own_post, , drop = FALSE],
                                    XY[i_post, , drop = FALSE])
    # A reference cell's distance to its own arm is its nearest other member.
    for (idx in list(list(i_pre, TRUE), list(i_post, FALSE))) {
        ii <- idx[[1L]]
        if (length(ii) < 2L) next
        dd <- vapply(seq_along(ii), function(k)
            min(sqrt(colSums((t(XY[ii[-k], , drop = FALSE]) -
                                  as.numeric(XY[ii[k], ]))^2))), numeric(1))
        if (isTRUE(idx[[2L]])) d_pre[ii] <- dd else d_post[ii] <- dd
    }
    pmax(d_pre, d_post)
}

# Pre-IP TREP -> Post-IP TREP axis, rescaled so Pre is 0 and Post is 1
# (not clamped).
v2_prepost_proj <- function(pca_mat, barcodes, grp, pre_lab, post_lab) {
    if (is.null(pca_mat) || is.null(barcodes)) return(NULL)
    idx <- match(barcodes, rownames(pca_mat))
    if (anyNA(idx)) return(NULL)
    X    <- as.matrix(pca_mat)[idx, , drop = FALSE]
    i_pre  <- which(!is.na(grp) & grp == pre_lab)
    i_post <- which(!is.na(grp) & grp == post_lab)
    if (length(i_pre) < V2_STAT_MIN_GROUP ||
        length(i_post) < V2_STAT_MIN_GROUP) return(NULL)
    c_pre  <- colMeans(X[i_pre, , drop = FALSE])
    c_post <- colMeans(X[i_post, , drop = FALSE])
    ax     <- c_post - c_pre
    d2     <- sum(ax^2)
    if (!is.finite(d2) || d2 <= 0) return(NULL)
    tproj <- as.numeric((X %*% ax - sum(c_pre * ax)) / d2)
    # Off-axis distance: norm of the component perpendicular to the axis.
    resid <- sweep(X, 2L, c_pre, "-") - outer(tproj, ax)
    list(proj = tproj, off = sqrt(rowSums(resid^2)),
         sep = sqrt(d2), n_pre = length(i_pre), n_post = length(i_post))
}

# == v2_merge_groups: attach class and group columns, derive plot groups ======
# plot_group5 holds the levels of the colour scale.
v2_merge_groups <- function(df_umap, rc_src, group_col = V2_ORIGREF_GROUP_COL) {
    if (is.null(df_umap)) return(NULL)

    need <- c("barcode", "cell_class", "is_orig_skcs", "is_recovered", group_col)
    absent <- setdiff(need, names(rc_src))
    if (length(absent) > 0)
        stop("Validation1 sheet lacks: ", paste(absent, collapse = ", "),
             call. = FALSE)

    rc_merge <- data.frame(
        barcode      = as.character(rc_src$barcode),
        is_orig_skcs = as_flag(rc_src$is_orig_skcs),
        is_recovered = as_flag(rc_src$is_recovered),
        group_ref    = as.character(rc_src[[group_col]]),
        stringsAsFactors = FALSE)
    # Donor travels with the coordinates.
    rc_merge$donor <- if ("donor" %in% names(rc_src))
        as.character(rc_src$donor) else NA_character_

    merged <- merge(df_umap, rc_merge, by = "barcode", all.x = TRUE)

    # An Original+Refined cell outside the four subgroups is "Unknown subgroup".
    merged$plot_group <- dplyr::case_when(
        !is.na(merged$is_recovered) & merged$is_recovered      ~ "Recovered",
        merged$cell_class == "Original only"                   ~ "Original only",
        merged$cell_class == "Original + Refined" &
            merged$group_ref %in% names(FOUR_COLS)             ~ merged$group_ref,
        merged$cell_class == "Original + Refined"              ~ "Unknown subgroup",
        TRUE                                                   ~ "Non-selected")
    n_unknown_grp <- sum(merged$plot_group == "Unknown subgroup", na.rm = TRUE)
    if (n_unknown_grp > 0L)
        cat(sprintf(paste0("    %d Original+Refined cell(s) have a `%s` outside",
                           " the four IP x TREP subgroups (excluded from tests).\n"),
                    n_unknown_grp, group_col))
    attr(merged, "n_unknown_group") <- n_unknown_grp
    merged$plot_group <- factor(
        merged$plot_group,
        levels = c("Non-selected", names(FOUR_COLS),
                   "Original only", "Recovered", "Unknown subgroup"))

    merged$plot_group5 <- factor(
        ifelse(as.character(merged$plot_group) %in% GROUP5_LEVS,
               as.character(merged$plot_group), "Non-selected"),
        levels = GROUP5_LEVS)

    merged$cell_class_f <- factor(merged$cell_class, levels = CMRF_CLASS_LEVS)
    rownames(merged) <- merged$barcode
    merged
}

# == v2_plotly3d: interactive 3D scatter, one trace per group =================
v2_plotly3d <- function(df_3d, title_text) {

    trace_specs <- list(
        # Excluded cells are drawn by the next trace.
        list(label = "Non-selected",
             mask  = df_3d$plot_group == "Non-selected" &
                 df_3d$cell_class != "Excluded",
             col   = "#C8C8C8", sym = "circle",
             sz = V2_3D_BG_SIZE, op = 0.25),
        list(label = "Excluded",
             mask  = df_3d$cell_class == "Excluded",
             col   = VP_COL_EXCL, sym = "circle",
             sz = V2_3D_BG_SIZE, op = 0.40),
        list(label = "Original only",
             mask  = df_3d$cell_class == "Original only",
             col   = ORIGONLY_COL, sym = "square",
             sz = V2_3D_MARKER_SIZE, op = 0.88, keep_empty = TRUE),
        # Recovered is drawn before the four subgroups.
        list(label = "Recovered",
             mask  = df_3d$cell_class == "Recovered",
             col   = RECOV_COL_04, sym = "diamond",
             sz = V2_3D_MARKER_SIZE, op = 0.80, keep_empty = TRUE),
        list(label = "Pre-IP TREP",
             mask  = df_3d$plot_group == "Pre-IP TREP",
             col   = unname(FOUR_COLS["Pre-IP TREP"]),
             sym = "circle", sz = V2_3D_MARKER_SIZE, op = 0.80,
             keep_empty = TRUE),
        list(label = "Pre-IP non-TREP",
             mask  = df_3d$plot_group == "Pre-IP non-TREP",
             col   = unname(FOUR_COLS["Pre-IP non-TREP"]),
             sym = "circle", sz = V2_3D_MARKER_SIZE, op = 0.80,
             keep_empty = TRUE),
        list(label = "Post-IP TREP",
             mask  = df_3d$plot_group == "Post-IP TREP",
             col   = unname(FOUR_COLS["Post-IP TREP"]),
             sym = "circle", sz = V2_3D_MARKER_SIZE, op = 0.80,
             keep_empty = TRUE),
        list(label = "Post-IP non-TREP",
             mask  = df_3d$plot_group == "Post-IP non-TREP",
             col   = unname(FOUR_COLS["Post-IP non-TREP"]),
             sym = "circle", sz = V2_3D_MARKER_SIZE, op = 0.80,
             keep_empty = TRUE))

    fig <- plotly::plot_ly()
    .empty_lab <- character(0)
    for (tr in trace_specs) {
        m <- tr$mask
        m[is.na(m)] <- FALSE
        sub <- df_3d[m, , drop = FALSE]
        # An empty group keeps its legend entry: one size-0 point at the centre.
        n_sub <- nrow(sub)
        if (n_sub == 0L) {
            if (!isTRUE(tr$keep_empty)) next
            .empty_lab <- c(.empty_lab, tr$label)
            xs <- mean(df_3d$UMAP_1, na.rm = TRUE)
            ys <- mean(df_3d$UMAP_2, na.rm = TRUE)
            zs <- mean(df_3d$UMAP_3, na.rm = TRUE)
            if (!is.finite(xs)) xs <- 0
            if (!is.finite(ys)) ys <- 0
            if (!is.finite(zs)) zs <- 0
            mk_sz <- 0
            hov   <- "skip"
            hover_txt <- ""
        } else {
            xs <- sub$UMAP_1; ys <- sub$UMAP_2; zs <- sub$UMAP_3
            mk_sz <- tr$sz
            hov   <- "text"
            hover_txt <- paste0("Group: ",       tr$label,
                                "<br>Class: ",   sub$cell_class,
                                "<br>Barcode: ", sub$barcode)
        }
        fig <- plotly::add_trace(
            fig,
            x = xs, y = ys, z = zs,
            type = "scatter3d", mode = "markers",
            name = sprintf("%s (%d)", tr$label, n_sub),
            marker = list(color = tr$col, symbol = tr$sym, size = mk_sz,
                          opacity = tr$op,
                          line = list(color = "white",
                                      width = V2_3D_OUTLINE_W)),
            text = hover_txt, hoverinfo = hov, showlegend = TRUE)
    }
    if (length(.empty_lab) > 0L)
        cat(sprintf(paste0("    3D scatter: %d group(s) empty on this dataset",
                           " and drawn as a zero-count legend entry: %s\n"),
                    length(.empty_lab), paste(.empty_lab, collapse = ", ")))
    plotly::layout(
        fig,
        title  = list(text = title_text, x = 0.01),
        scene  = list(xaxis = list(title = "UMAP 1"),
                      yaxis = list(title = "UMAP 2"),
                      zaxis = list(title = "UMAP 3"),
                      aspectmode = "cube",
                      camera = list(eye = list(x = 1.5, y = 1.5, z = 1.1)),
                      bgcolor = "white"),
        # Uniform legend swatches, cube aspect and a fixed camera.
        legend = list(title = list(text = "<b>Cell group</b>"),
                      itemsizing = "constant"))
}


# =============================================================================
# SECTION 0f: Statistical helpers
# =============================================================================

# Shapiro-Wilk p, subsampled to 5000.
v2_sw_p <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) < 3L) return(NA_real_)
    if (length(x) > 5000L) { set.seed(V2_SEED); x <- sample(x, 5000L) }
    tryCatch(shapiro.test(x)$p.value, error = function(e) NA_real_)
}

# == v2_run_tost: two one-sided tests for equivalence on one axis =============
# Margin delta = delta_frac x pooled SD.
v2_run_tost <- function(x1, x2, delta_frac = V2_TOST_DELTA_FRAC) {
    x1 <- x1[!is.na(x1)]; x2 <- x2[!is.na(x2)]
    n1 <- length(x1); n2 <- length(x2)
    if (n1 < 3L || n2 < 3L)
        return(list(p_tost = NA_real_, p_lower = NA_real_, p_upper = NA_real_,
                    delta = NA_real_, delta_frac = delta_frac,
                    pooled_sd = NA_real_, mean_diff = NA_real_, equiv = NA,
                    n1 = n1, n2 = n2, sw_p1 = NA_real_, sw_p2 = NA_real_))

    pooled_sd <- sqrt(((n1 - 1) * var(x1) + (n2 - 1) * var(x2)) / (n1 + n2 - 2))
    delta     <- delta_frac * pooled_sd
    mean_diff <- mean(x1) - mean(x2)

    t_lower <- tryCatch(t.test(x1, x2, mu = -delta, alternative = "greater"),
                        error = function(e) NULL)
    t_upper <- tryCatch(t.test(x1, x2, mu =  delta, alternative = "less"),
                        error = function(e) NULL)
    p_lower <- if (!is.null(t_lower)) t_lower$p.value else NA_real_
    p_upper <- if (!is.null(t_upper)) t_upper$p.value else NA_real_
    p_tost  <- if (!is.na(p_lower) && !is.na(p_upper))
        max(p_lower, p_upper) else NA_real_

    list(p_tost = round(p_tost, 6), p_lower = round(p_lower, 6),
         p_upper = round(p_upper, 6), delta = round(delta, 4),
         delta_frac = delta_frac, pooled_sd = round(pooled_sd, 4),
         mean_diff = round(mean_diff, 4),
         equiv = !is.na(p_tost) && p_tost < 0.05,
         n1 = n1, n2 = n2,
         sw_p1 = round(v2_sw_p(x1), 6), sw_p2 = round(v2_sw_p(x2), 6))
}

# == v2_permanova_pair: PERMANOVA (adonis2) and PERMDISP (betadisper) =========
# Each group is subsampled to max_n; the n used is reported.
v2_permanova_pair <- function(focal_mat, ref_mat, n_perm = V2_PERM_N,
                              max_n = V2_PERM_MAX_N, seed = V2_SEED) {

    n_focal_all <- nrow(focal_mat); n_ref_all <- nrow(ref_mat)
    na_result <- list(R2 = NA_real_, F_stat = NA_real_, p_value = NA_real_,
                      disp_F = NA_real_, disp_p = NA_real_,
                      avg_dist_focal = NA_real_, avg_dist_ref = NA_real_,
                      n_focal = n_focal_all, n_ref = n_ref_all,
                      n_focal_used = 0L, n_ref_used = 0L,
                      subsampled = FALSE, n_perm = n_perm, note = "")

    if (!.v2_vegan_ok)   { na_result$note <- "vegan not installed"; return(na_result) }
    if (n_focal_all < 3L || n_ref_all < 3L) {
        na_result$note <- "n < 3"; return(na_result) }

    set.seed(seed)
    sub_rows <- function(m) {
        if (nrow(m) <= max_n) return(m)
        m[sort(sample.int(nrow(m), max_n)), , drop = FALSE]
    }
    fm <- sub_rows(as.matrix(focal_mat))
    rm_ <- sub_rows(as.matrix(ref_mat))
    subsampled <- (nrow(fm) < n_focal_all) || (nrow(rm_) < n_ref_all)

    comb_mat <- rbind(fm, rm_)
    grp_vec  <- factor(c(rep("focal", nrow(fm)), rep("ref", nrow(rm_))))
    # One distance matrix for both tests.
    d_comb <- dist(comb_mat, method = "euclidean")

    set.seed(seed)
    perm_res <- tryCatch(
        vegan::adonis2(d_comb ~ grp_vec, permutations = n_perm, by = "terms"),
        error = function(e) {
            cat(sprintf("      PERMANOVA error: %s\n", conditionMessage(e)))
            NULL })
    if (is.null(perm_res)) {
        na_result$note <- "adonis2 failed"
        na_result$n_focal_used <- nrow(fm); na_result$n_ref_used <- nrow(rm_)
        na_result$subsampled <- subsampled
        return(na_result)
    }

    set.seed(seed)
    disp_res <- tryCatch({
        bd <- vegan::betadisper(d_comb, grp_vec)
        pt <- vegan::permutest(bd, permutations = n_perm)
        list(F_stat  = round(pt$tab[["F"]][1], 4),
             p_value = round(pt$tab[["Pr(>F)"]][1], 6),
             avg_dist_focal = round(mean(bd$distances[grp_vec == "focal"]), 4),
             avg_dist_ref   = round(mean(bd$distances[grp_vec == "ref"]),   4))
    }, error = function(e) {
        cat(sprintf("      PERMDISP error: %s\n", conditionMessage(e)))
        list(F_stat = NA_real_, p_value = NA_real_,
             avg_dist_focal = NA_real_, avg_dist_ref = NA_real_) })

    list(R2      = round(perm_res$R2[1], 6),
         F_stat  = round(perm_res$F[1], 4),
         p_value = round(perm_res[["Pr(>F)"]][1], 6),
         disp_F  = disp_res$F_stat,
         disp_p  = disp_res$p_value,
         avg_dist_focal = disp_res$avg_dist_focal,
         avg_dist_ref   = disp_res$avg_dist_ref,
         n_focal = n_focal_all, n_ref = n_ref_all,
         n_focal_used = nrow(fm), n_ref_used = nrow(rm_),
         subsampled = subsampled, n_perm = n_perm, note = "")
}

# == v2_compute_ilisi: integration LISI over Gaussian-weighted kNN ============
# With two groups, 1 is no mixing and 2 is a 50/50 neighbourhood.
v2_compute_ilisi <- function(umap_mat, labels, perplexity = V2_ILISI_PERPLEXITY) {

    umap_mat <- as.matrix(umap_mat)
    labels   <- as.character(labels)
    n        <- nrow(umap_mat)
    k        <- min(as.integer(perplexity) * 3L, n - 1L)

    if (n < 10L || k < 3L)
        return(list(per_cell = rep(NA_real_, n), mean = NA_real_,
                    median = NA_real_, by_group = NULL,
                    n_groups = length(unique(labels)), k = k))

    knn_res <- FNN::get.knn(umap_mat, k = k)
    nn_idx  <- knn_res$nn.index
    nn_dist <- knn_res$nn.dist
    unique_labels <- sort(unique(labels))

    lisi_vals <- vapply(seq_len(n), function(i) {
        dists   <- nn_dist[i, ]
        sigma   <- max(median(dists), 1e-10)
        weights <- exp(-dists^2 / (2 * sigma^2))
        weights <- weights / sum(weights)
        nb_labels <- labels[nn_idx[i, ]]
        props <- vapply(unique_labels, function(lbl) sum(weights[nb_labels == lbl]),
                        numeric(1L))
        1.0 / sum(props^2)
    }, numeric(1L))

    by_grp <- data.frame(
        group = unique_labels,
        n = as.integer(vapply(unique_labels, function(l) sum(labels == l),
                              integer(1L))),
        mean_iLISI = round(vapply(unique_labels, function(l)
            mean(lisi_vals[labels == l], na.rm = TRUE), numeric(1L)), 4),
        median_iLISI = round(vapply(unique_labels, function(l)
            median(lisi_vals[labels == l], na.rm = TRUE), numeric(1L)), 4),
        pct_gt_1p5 = round(vapply(unique_labels, function(l)
            100 * mean(lisi_vals[labels == l] > 1.5, na.rm = TRUE), numeric(1L)), 2),
        stringsAsFactors = FALSE)

    list(per_cell = round(lisi_vals, 4),
         mean     = round(mean(lisi_vals, na.rm = TRUE), 4),
         median   = round(median(lisi_vals, na.rm = TRUE), 4),
         by_group = by_grp,
         n_groups = length(unique_labels),
         k        = k)
}

# == v2_neighbour_metrics: reference fraction among k neighbours ==============
# Enrichment = that fraction / the reference group's global share.
v2_neighbour_metrics <- function(focal_mat, ref_mat, focal_label, ref_label,
                                 k = V2_KNN_K) {
    focal_mat <- as.matrix(focal_mat)
    ref_mat   <- as.matrix(ref_mat)
    n_focal   <- nrow(focal_mat); n_ref <- nrow(ref_mat)
    n_total   <- n_focal + n_ref

    if (n_focal < 3L || n_ref < 3L)
        return(list(mean_ref_frac = NA_real_, pct_gt_40 = NA_real_,
                    pct_gt_50 = NA_real_, global_prop = NA_real_,
                    mean_enrichment = NA_real_, median_enrichment = NA_real_,
                    k = k))

    global_prop <- n_ref / n_total
    k_eff       <- min(as.integer(k), n_total - 1L)
    comb_mat    <- rbind(focal_mat, ref_mat)
    labels      <- c(rep(focal_label, n_focal), rep(ref_label, n_ref))
    knn_res     <- FNN::get.knn(comb_mat, k = k_eff)

    ref_fracs <- vapply(seq_len(n_focal), function(i)
        sum(labels[knn_res$nn.index[i, ]] == ref_label) / k_eff, numeric(1L))
    enrichment <- ref_fracs / global_prop

    list(mean_ref_frac     = round(mean(ref_fracs), 4),
         pct_gt_40         = round(100 * mean(ref_fracs > 0.40), 2),
         pct_gt_50         = round(100 * mean(ref_fracs > 0.50), 2),
         global_prop       = round(global_prop, 4),
         mean_enrichment   = round(mean(enrichment), 4),
         median_enrichment = round(median(enrichment), 4),
         k = k_eff)
}

# == v2_compute_kbet: simplified kBET (chi-square per cell) ===================
# Cells with an expected count below 1 return NA and are counted.
v2_compute_kbet <- function(focal_mat, ref_mat, focal_label, ref_label,
                            k = V2_KNN_K, alpha = V2_KBET_ALPHA,
                            n_sample = V2_KBET_SAMPLE) {
    focal_mat <- as.matrix(focal_mat)
    ref_mat   <- as.matrix(ref_mat)
    n_focal   <- nrow(focal_mat); n_ref <- nrow(ref_mat)
    n_total   <- n_focal + n_ref

    if (n_focal < 3L || n_ref < 3L)
        return(list(rejection_rate = NA_real_, mean_pval = NA_real_,
                    median_pval = NA_real_, n_tested = 0L, n_undefined = 0L,
                    k = k, alpha = alpha, note = "n < 3"))

    global_props <- c(n_focal / n_total, n_ref / n_total)
    k_eff        <- min(as.integer(k), n_total - 1L)
    comb_mat     <- rbind(focal_mat, ref_mat)
    labels       <- c(rep(focal_label, n_focal), rep(ref_label, n_ref))
    unique_labs  <- c(focal_label, ref_label)
    knn_res      <- FNN::get.knn(comb_mat, k = k_eff)

    test_idx <- seq_len(n_total)
    if (!is.null(n_sample) && n_total > n_sample) {
        set.seed(V2_SEED)
        test_idx <- sort(sample.int(n_total, n_sample))
    }

    expected <- global_props * k_eff
    if (any(expected < 1))
        return(list(rejection_rate = NA_real_, mean_pval = NA_real_,
                    median_pval = NA_real_, n_tested = 0L,
                    n_undefined = length(test_idx), k = k_eff, alpha = alpha,
                    note = sprintf(paste0("expected count %.2f below 1 at k = %d;",
                                          " chi-square undefined"),
                                   min(expected), k_eff)))

    pvals <- vapply(test_idx, function(i) {
        nn_labs  <- labels[knn_res$nn.index[i, ]]
        observed <- vapply(unique_labs, function(l) sum(nn_labs == l), integer(1L))
        tryCatch(suppressWarnings(
            chisq.test(observed, p = global_props)$p.value),
            error = function(e) NA_real_)
    }, numeric(1L))

    ok <- !is.na(pvals)
    list(rejection_rate = if (any(ok)) round(mean(pvals[ok] < alpha), 4) else NA_real_,
         mean_pval      = if (any(ok)) round(mean(pvals[ok]), 4) else NA_real_,
         median_pval    = if (any(ok)) round(median(pvals[ok]), 4) else NA_real_,
         n_tested       = sum(ok),
         n_undefined    = sum(!ok),
         k = k_eff, alpha = alpha, note = "")
}


# =============================================================================
# SECTION 0e: RDS size control
# =============================================================================
# Model environments are stripped and plot objects dropped before saving.
RDS_KEEP_PLOTS <- FALSE

rds_strip_gam_env <- function(m) {
    if (!inherits(m, c("gam", "glm", "lm"))) return(m)
    if (inherits(m$formula, "formula"))      environment(m$formula)      <- baseenv()
    if (inherits(m$pred.formula, "formula")) environment(m$pred.formula) <- baseenv()
    if (!is.null(m$terms))  attr(m$terms,  ".Environment") <- baseenv()
    if (!is.null(m$pterms)) attr(m$pterms, ".Environment") <- baseenv()
    if (!is.null(m$model) && !is.null(attr(m$model, "terms")))
        attr(attr(m$model, "terms"), ".Environment") <- baseenv()
    for (i in seq_along(m$smooth)) {
        if (inherits(m$smooth[[i]]$form, "formula"))
            environment(m$smooth[[i]]$form) <- baseenv()
        if (!is.null(m$smooth[[i]]$terms))
            attr(m$smooth[[i]]$terms, ".Environment") <- baseenv()
    }
    m
}

.rds_is_plot <- function(x)
    inherits(x, c("ggplot", "gg", "ggproto", "plotly", "htmlwidget",
                  "gtable", "grob", "gTree", "gList", "trellis",
                  "recordedplot"))

rds_slim <- function(x, keep_plots = RDS_KEEP_PLOTS) {
    if (inherits(x, c("gam", "glm", "lm"))) return(rds_strip_gam_env(x))
    if (!isTRUE(keep_plots) && .rds_is_plot(x)) return(NULL)
    if (is.list(x) && !is.data.frame(x)) {
        att <- attributes(x)
        out <- lapply(x, rds_slim, keep_plots = keep_plots)
        attributes(out) <- att
        return(out)
    }
    x
}

rds_save <- function(obj, path, label = "RDS") {
    saveRDS(rds_slim(obj), path)
    cat(sprintf("  %s saved: %s  (%.1f MB%s)\n", label, basename(path),
                file.size(path) / 1e6,
                if (isTRUE(RDS_KEEP_PLOTS)) ", plots kept" else ", plots dropped"))
    invisible(path)
}


# =============================================================================
# SECTION 0b: The pipeline for one arm
# =============================================================================

run_one_arm <- function(arm_tag, config_file, out_dir, file_tag) {

    config_path <- file.path(CEP_HOME, config_file)
    if (!file.exists(config_path))
        stop("Configuration not found: ", config_path, call. = FALSE)

    .p10v2_preamble <- character(0)
    assign("CEP_ARM", arm_tag, envir = globalenv())
    on.exit(suppressWarnings(rm("CEP_ARM", envir = globalenv())), add = TRUE)

    .p10v2_preamble <- tryCatch(
        utils::capture.output({
            cat("\n", strrep("=", 76), "\n", sep = "")
            cat("  ARM: ", arm_tag, "   config: ", config_file, "\n", sep = "")
            cat(strrep("=", 76), "\n", sep = "")
            source(config_path)
        }),
        error = function(e) {
            cat("\n  ARM: ", arm_tag, "   config: ", config_file, "\n", sep = "")
            stop(e)
        })

    .p10v2_needed <- c("ARM_TAG", "R_EXPORT_ROOT", "GOI", "DFG",
                       "DFG_COMPOSITE_LABEL", "PRIMARY_STYPE", "NORM_TAG",
                       "NORM_TARGET_SUM", "NORM_CEILING", "NEOC_DIR",
                       "NONNC_DIR", "NEOC_REGIONS", "NONNC_REGIONS")
    .p10v2_absent <- .p10v2_needed[!vapply(.p10v2_needed, exists, logical(1))]
    if (length(.p10v2_absent) > 0)
        stop("The configuration loaded from\n    ", config_path,
             "\n  does not define: ", paste(.p10v2_absent, collapse = ", "),
             "\n  This is the signature of a pre-migration config.", call. = FALSE)
    if (!identical(ARM_TAG, arm_tag))
        stop("The configuration reports ARM_TAG '", ARM_TAG,
             "' but this arm is '", arm_tag, "'.", call. = FALSE)

    .p10v2_preamble <- c(.p10v2_preamble, utils::capture.output(
        cat("  Config OK: ", ARM_TAG, " | ", GOI, " + ", length(DFG), " DFG | ",
            NORM_TAG, "\n", sep = "")))

    set.seed(V2_SEED)

    DATASETS <- c("NeoC", "NonNeoC")

    # == Input paths ==========================================================
    V1_DIR    <- cep_arm_dir_for("Part10e_Validation1", arm_tag)
    V1_PREFIX <- paste0(cep_out_prefix("Part10e", NORM_TAG,
                                       cep_file_tag_for(arm_tag)), "_V1")
    p10v2_v1_xlsx <- function(ds)
        file.path(V1_DIR, paste0(V1_PREFIX, "_", ds,
                                 "_SK_interaction_analysis.xlsx"))
    p10v2_v1_rds <- function(ds)
        file.path(V1_DIR, paste0(V1_PREFIX, "_", ds, "_results.rds"))

    # == Output ===============================================================
    OUT_DIR <- out_dir
    dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
    FILE_PREFIX <- paste0(cep_out_prefix(PART_TAG, NORM_TAG,
                                         cep_file_tag_for(arm_tag)), "_V2")

    # == Logging ==============================================================
    LOG_MESSAGES_V2 <- TRUE
    LOG_PATH_V2 <- file.path(OUT_DIR, paste0(FILE_PREFIX, "_console_log_",
                                             format(Sys.time(), "%Y%m%d_%H%M%S"),
                                             ".txt"))
    MSG_PATH_V2 <- file.path(OUT_DIR, paste0(FILE_PREFIX,
                                             "_console_log_messages.tmp"))
    while (sink.number() > 0) sink()
    if (sink.number(type = "message") != 2L) sink(type = "message")
    # Set before on.exit() registers.
    .p10v2_completed <- FALSE
    .p10v2_log_con <- file(LOG_PATH_V2, open = "wt")
    .p10v2_msg_con <- if (isTRUE(LOG_MESSAGES_V2))
        file(MSG_PATH_V2, open = "wt") else NULL
    sink(.p10v2_log_con, append = FALSE, split = TRUE)
    if (!is.null(.p10v2_msg_con)) sink(.p10v2_msg_con, type = "message")
    # Exit handler: four stages, each in its own tryCatch.
    on.exit({
        tryCatch({
            if (!isTRUE(.p10v2_completed))
                cat("\n", strrep("!", 72),
                    "\n  Arm did not complete; the error is printed to the console.\n",
                    strrep("!", 72), "\n", sep = "")
        }, error = function(e) invisible(NULL))
        tryCatch({
            if (sink.number(type = "message") != 2L) sink(type = "message")
            while (sink.number() > 0) sink()
        }, error = function(e) invisible(NULL))
        tryCatch(if (!is.null(.p10v2_msg_con)) close(.p10v2_msg_con),
                 error = function(e) invisible(NULL))
        tryCatch(close(.p10v2_log_con), error = function(e) invisible(NULL))
        tryCatch({
            if (!is.null(.p10v2_msg_con) && file.exists(MSG_PATH_V2)) {
                .msg_txt <- readLines(MSG_PATH_V2, warn = FALSE)
                if (length(.msg_txt) > 0)
                    cat(c("", strrep("=", 72),
                          "  MESSAGE STREAM (warnings, messages, errors)",
                          strrep("=", 72), .msg_txt, ""),
                        sep = "\n", file = LOG_PATH_V2, append = TRUE)
                unlink(MSG_PATH_V2)
            }
        }, error = function(e) invisible(NULL))
    }, add = TRUE)

    .p10v2_warn_opt  <- getOption("warn")
    options(warn = 1)
    on.exit(options(warn = .p10v2_warn_opt), add = TRUE)
    cat(sprintf("Part10e Validation2 (UMAP) - log started: %s\n\n", Sys.time()))
    if (length(.p10v2_preamble) > 0)
        cat(paste(.p10v2_preamble, collapse = "\n"), "\n\n", sep = "")

    cat("\n", strrep("=", 72), "\n")
    cat("  Part10e Validation2 - CMRF classes in transcriptomic space\n")
    cat("  Arm          :", ARM_TAG, "\n")
    cat("  GOI          :", GOI, "\n")
    cat("  Normalisation:", NORM_TAG, " log1p(UMIs/total x ",
        format(NORM_TARGET_SUM, scientific = FALSE), ")\n")
    cat("  Embedding    : Monocle3 UMAP on", V2_EMBED_SOURCE, "|",
        V2_NUM_DIM, "PCA dimensions\n")
    cat("  OrigRef group:", V2_ORIGREF_GROUP_COL, "\n")
    cat("  V1 input     :", V1_DIR, "\n")
    cat("  Output dir   :", OUT_DIR, "\n")
    cat(strrep("=", 72), "\n\n")

    INCOMPLETE_MARKER_V2 <- file.path(OUT_DIR, "_RUN_INCOMPLETE.txt")
    writeLines(c("This arm did not finish.", "",
                 paste("Arm     :", arm_tag),
                 paste("Part    :", PART_DIR_TAG),
                 paste("Started :", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), "",
                 "Any files in this folder are from a partial run.",
                 "This file is removed automatically when the arm completes."),
               INCOMPLETE_MARKER_V2)


    # =========================================================================
    # SECTION A0: Guardrail ledger
    # =========================================================================
    GR_LOG_V2  <- file.path(OUT_DIR, paste0(FILE_PREFIX, "_Guardrail_Log.txt"))
    .grv2_log  <- list()
    .grv2_pass <- 0L; .grv2_fail <- 0L; .grv2_warn <- 0L

    grv2 <- function(check_id, status, detail) {
        # A zero-length detail would vanish from console and ledger.
        if (length(detail) != 1L || is.na(detail))
            detail <- paste0("[detail unavailable]")
        # An unrecognised status is recorded as FAIL.
        if (!status %in% c("PASS", "FAIL", "WARN")) {
            detail <- paste0("[status '", as.character(status)[1],
                             "' is not PASS/FAIL/WARN] ", detail)
            status <- "FAIL"
        }
        tag <- switch(status, PASS = "[PASS]", FAIL = "[FAIL]", WARN = "[WARN]")
        cat(sprintf("  %s  %-18s  %s\n", tag, check_id, detail))
        .grv2_log[[length(.grv2_log) + 1L]] <<- list(
            id = check_id, status = status, detail = detail)
        if (status == "PASS") .grv2_pass <<- .grv2_pass + 1L
        if (status == "FAIL") .grv2_fail <<- .grv2_fail + 1L
        if (status == "WARN") .grv2_warn <<- .grv2_warn + 1L
    }

    grv2_write <- function() {
        con <- file(GR_LOG_V2, "wt")
        # Close the connection on any exit.
        on.exit(try(close(con), silent = TRUE), add = TRUE)
        writeLines(c(strrep("=", 72),
                     paste0("Part 10e Validation2 - Guardrail Log [", NORM_TAG, "]"),
                     paste("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
                     paste("Arm :", ARM_TAG),
                     paste("Embedding source:", V2_EMBED_SOURCE),
                     strrep("=", 72), ""), con)
        for (item in .grv2_log)
            writeLines(sprintf("  [%s]  %-18s  %s", item$status, item$id,
                               item$detail), con)
        writeLines(c("", strrep("=", 72),
                     sprintf("  Total : %d", .grv2_pass + .grv2_fail + .grv2_warn),
                     sprintf("  PASS  : %d", .grv2_pass),
                     sprintf("  WARN  : %d", .grv2_warn),
                     sprintf("  FAIL  : %d", .grv2_fail),
                     if (.grv2_fail == 0 && .grv2_warn == 0)
                         "\n  Result: no failures, no warnings"
                     else if (.grv2_fail == 0)
                         "\n  Result: no failures; review the warnings above"
                     else "\n  Result: failures present",
                     strrep("=", 72)), con)
        cat("  Guardrail log written:", GR_LOG_V2, "\n")
    }

    cat("  SECTION A0 - PRE-ANALYSIS GUARDRAILS\n\n")

    for (ds in DATASETS) {
        .p <- p10v2_v1_xlsx(ds)
        .r <- p10v2_v1_rds(ds)
        grv2(paste0("GR-V2-01-", ds),
             if (file.exists(.p) || file.exists(.r)) "PASS" else "FAIL",
             if (file.exists(.p))
                 sprintf("Validation1 workbook: %s",
                         file.path(basename(dirname(.p)), basename(.p)))
             else if (file.exists(.r))
                 sprintf("workbook absent; will use RDS: %s", basename(.r))
             else sprintf("neither found under %s", V1_DIR))
    }

    for (.d in list(list(tag = "NeoC", dir = NEOC_DIR),
                    list(tag = "NonNeoC", dir = NONNC_DIR))) {
        cdir <- file.path(.d$dir, "counts", "raw_selected")
        need <- file.path(cdir, c("matrix.mtx", "matrix_barcodes.tsv",
                                  "matrix_features.tsv"))
        grv2(paste0("GR-V2-02-", .d$tag),
             if (all(file.exists(need))) "PASS" else "FAIL",
             if (all(file.exists(need)))
                 sprintf("raw_selected present: %s", cdir)
             else sprintf("raw_selected incomplete: %s",
                          paste(basename(need[!file.exists(need)]),
                                collapse = ", ")))
    }

    cat(sprintf("\n  Pre-analysis guardrails: %d PASS / %d WARN / %d FAIL\n\n",
                .grv2_pass, .grv2_warn, .grv2_fail))


    # =========================================================================
    # SECTION A1: Donor lookup
    # =========================================================================
    # Donor by barcode from the per-region gam_data CSVs.
    .p10v2_donor_cache <- list()
    p10v2_donor_map <- function(ds) {
        if (!is.null(.p10v2_donor_cache[[ds]])) return(.p10v2_donor_cache[[ds]])
        gam_dir <- if (identical(ds, "NeoC")) NEOC_DIR else NONNC_DIR
        regs    <- if (identical(ds, "NeoC")) NEOC_REGIONS else NONNC_REGIONS
        stype   <- gsub("-", "_", gsub("/", "_", PRIMARY_STYPE))
        bc <- character(0); don <- character(0)
        for (reg in regs) {
            gp <- file.path(gam_dir, "gam_data",
                            paste0(reg, "_", stype, "_selected.csv"))
            if (!file.exists(gp)) next
            df <- read.csv(gp, stringsAsFactors = FALSE, check.names = FALSE)
            if (!all(c("barcode", "donor") %in% colnames(df)))
                stop(sprintf("Region %s: gam_data lacks barcode or donor.", reg),
                     call. = FALSE)
            bc  <- c(bc,  df$barcode)
            don <- c(don, as.character(df$donor))
        }
        if (length(bc) == 0L)
            stop("No gam_data CSVs found for ", ds, call. = FALSE)
        # Duplicated barcodes must agree on the donor.
        dup <- duplicated(bc)
        if (any(dup)) {
            .first <- setNames(don[!dup], bc[!dup])
            .clash <- unique(bc[dup][don[dup] != .first[bc[dup]]])
            .clash <- .clash[!is.na(.clash)]
            if (length(.clash) > 0)
                stop(sprintf(paste0("%s: %d barcode(s) appear in more than one",
                                    " region with different donors, e.g. %s."),
                             ds, length(.clash),
                             paste(utils::head(.clash, 3), collapse = ", ")),
                     call. = FALSE)
            cat(sprintf(paste0("    %d duplicated barcode(s) across regions,",
                               " all agreeing on the donor; keeping one each\n"),
                        sum(dup)))
            bc <- bc[!dup]; don <- don[!dup]
        }
        m <- setNames(don, bc)
        .p10v2_donor_cache[[ds]] <<- m
        m
    }


    # =========================================================================
    # SECTION A2: Label and count loaders
    # =========================================================================

    # Labels from the Validation1 workbook (Step4_Refined_Cells), RDS as fallback;
    # columns read by name.
    p10v2_load_labels <- function(ds) {
        xp <- p10v2_v1_xlsx(ds)
        rp <- p10v2_v1_rds(ds)
        rc <- NULL; src <- NA_character_

        if (file.exists(xp)) {
            rc <- tryCatch(
                as.data.frame(readxl::read_excel(xp,
                                                 sheet = "Step4_Refined_Cells"),
                              stringsAsFactors = FALSE),
                error = function(e) {
                    cat("    Validation1 workbook unreadable:",
                        conditionMessage(e), "\n"); NULL })
            if (!is.null(rc)) src <- basename(xp)
        }
        if (is.null(rc) && file.exists(rp)) {
            v1 <- tryCatch(readRDS(rp), error = function(e) NULL)
            if (!is.null(v1) && !is.null(v1$rc_v3)) {
                rc  <- as.data.frame(v1$rc_v3, stringsAsFactors = FALSE)
                src <- basename(rp)
            }
        }
        if (is.null(rc))
            stop("No Validation1 output for ", ds, " under ", V1_DIR,
                 call. = FALSE)

        need <- c("barcode", "cell_class", "is_orig_skcs", "is_recovered",
                  V2_ORIGREF_GROUP_COL)
        absent <- setdiff(need, names(rc))
        if (length(absent) > 0)
            stop("Step4_Refined_Cells in ", src, " lacks: ",
                 paste(absent, collapse = ", "), call. = FALSE)

        rc$barcode <- as.character(rc$barcode)
        if (anyNA(rc$barcode) || any(!nzchar(rc$barcode)))
            stop("Step4_Refined_Cells in ", src, " has empty barcodes.",
                 call. = FALSE)
        if (anyDuplicated(rc$barcode) > 0L)
            stop("Step4_Refined_Cells in ", src, " has duplicated barcodes.",
                 call. = FALSE)

        rc$cell_class <- as.character(rc$cell_class)
        # is_recovered and cell_class == "Recovered" must agree.
        .rc_flag  <- as_flag(rc$is_recovered)
        .rc_class <- rc$cell_class == "Recovered"
        .rc_diff  <- sum(!is.na(.rc_flag) & .rc_flag != .rc_class) +
                     sum(is.na(.rc_flag) & .rc_class)
        if (.rc_diff > 0)
            stop(.rc_diff, " cell(s) in ", src, " disagree between ",
                 "`is_recovered` and `cell_class == \"Recovered\"`.",
                 call. = FALSE)
        n_all  <- nrow(rc)
        rc_goi <- rc[rc$cell_class %in% CMRF_CLASS_LEVS, , drop = FALSE]
        rownames(rc_goi) <- rc_goi$barcode

        cat(sprintf("    Labels from %s: %d row(s), %d in the four classes\n",
                    src, n_all, nrow(rc_goi)))
        for (cl in CMRF_CLASS_LEVS)
            cat(sprintf("      %-22s %d\n", cl, sum(rc_goi$cell_class == cl)))

        list(rc = rc_goi, n_all = n_all,
             n_dropped = n_all - nrow(rc_goi), source = src)
    }

    # Raw counts: cells x genes on disk, genes x cells for Monocle3.
    p10v2_load_counts <- function(ds, keep_bc) {
        ds_dir <- if (identical(ds, "NeoC")) NEOC_DIR else NONNC_DIR
        cdir   <- file.path(ds_dir, "counts", "raw_selected")
        mtx_p  <- file.path(cdir, "matrix.mtx")
        bar_p  <- file.path(cdir, "matrix_barcodes.tsv")
        fea_p  <- file.path(cdir, "matrix_features.tsv")
        if (!all(file.exists(c(mtx_p, bar_p, fea_p))))
            stop("raw_selected counts incomplete under: ", cdir, call. = FALSE)

        X  <- as(readMM(mtx_p), "CsparseMatrix")
        bc <- readLines(bar_p, warn = FALSE)
        ft <- read.table(fea_p, sep = "\t", header = FALSE,
                         stringsAsFactors = FALSE)
        if (ncol(ft) < 2L)
            stop("matrix_features.tsv has ", ncol(ft),
                 " column(s); gene symbols are expected in column 2.",
                 call. = FALSE)
        if (nrow(X) != length(bc) || ncol(X) != nrow(ft))
            stop(sprintf(paste0("raw_selected dimensions disagree: matrix %d x %d,",
                                " barcodes %d, features %d"),
                         nrow(X), ncol(X), length(bc), nrow(ft)), call. = FALSE)
        rownames(X) <- bc
        colnames(X) <- make.unique(as.character(ft[, 2]))

        missing_bc <- setdiff(keep_bc, bc)
        if (length(missing_bc) > 0)
            stop(sprintf(paste0("%d labelled barcode(s) are absent from",
                                " raw_selected for %s (first: %s)."),
                         length(missing_bc), ds, missing_bc[1]), call. = FALSE)

        raw_int <- {
            v <- as.numeric(X@x[seq_len(min(10000L, length(X@x)))])
            all(abs(v - round(v)) < 0.01)
        }
        raw_max <- max(X)

        Xs <- Matrix::t(X[keep_bc, , drop = FALSE])   # genes x cells
        rm(X); gc(verbose = FALSE)

        # Per-gene detection count from the sparse row index.
        n_gene_all <- nrow(Xs)
        # Column-compressed format and no stored zeros before counting.
        Xs  <- Matrix::drop0(as(Xs, "CsparseMatrix"))
        det <- tabulate(Xs@i + 1L, nbins = nrow(Xs))
        Xs  <- Xs[det >= V2_MIN_CELLS_PER_GENE, , drop = FALSE]
        cat(sprintf(paste0("    Counts: %d cell(s); genes %d -> %d after the",
                           " >= %d-cell detection filter\n"),
                    ncol(Xs), n_gene_all, nrow(Xs), V2_MIN_CELLS_PER_GENE))

        list(mat = Xs, n_gene_all = n_gene_all, n_gene_kept = nrow(Xs),
             raw_integer_like = raw_int, raw_max = raw_max)
    }


    # =========================================================================
    # SECTION B: Main function (runs per dataset)
    # =========================================================================

    run_V2_dataset <- function(dataset_tag) {

        cat("\n", strrep("=", 72), "\n")
        cat("  Part10e Validation2:", dataset_tag, " [", ARM_TAG, "]\n")
        cat(strrep("=", 72), "\n\n")

        val_tabs <- list()

        # == STEP 1: labels ===================================================
        cat("  STEP 1 - CMRF class and group labels\n")
        lab    <- p10v2_load_labels(dataset_tag)
        rc     <- lab$rc
        n_cell <- nrow(rc)
        if (n_cell < 10L)
            stop("Only ", n_cell, " labelled cell(s) for ", dataset_tag,
                 call. = FALSE)

        # Donor per cell.
        donor_map <- p10v2_donor_map(dataset_tag)
        rc$donor  <- unname(donor_map[rc$barcode])
        n_no_don  <- sum(is.na(rc$donor))
        cat(sprintf("    Donors: %d distinct | %d cell(s) with no donor in gam_data\n",
                    length(unique(rc$donor[!is.na(rc$donor)])), n_no_don))

        # == STEP 2: counts ===================================================
        cat("\n  STEP 2 - Raw counts from the Part 8 export\n")
        cnt <- p10v2_load_counts(dataset_tag, rc$barcode)

        # == STEP 3: CDS and UMAP =============================================
        cat("\n  STEP 3 - Monocle3 CDS and UMAP\n")
        built  <- v2_build_cds(cnt$mat, rc, num_dim = V2_NUM_DIM, seed = V2_SEED)
        # PCA from the CDS, kept for the Pre->Post projection; rows named by barcode.
        pca_mat <- tryCatch({
            .rn <- reducedDimNames(built$cds)
            .pk <- .rn[toupper(.rn) == "PCA"]
            if (length(.pk) == 0L) NULL else {
                .m <- reducedDim(built$cds, .pk[1L])
                rownames(.m) <- colnames(built$cds)
                .m
            }
        }, error = function(e) NULL)
        if (is.null(pca_mat))
            cat("    WARNING: no PCA reducedDim; the Pre->Post projection",
                "will be skipped\n")
        cds_2d <- tryCatch(v2_reduce(built$cds, 2L, V2_SEED),
                           error = function(e) {
                               cat("    ERROR 2D UMAP:", conditionMessage(e), "\n")
                               NULL })
        cds_3d <- tryCatch(v2_reduce(built$cds, 3L, V2_SEED),
                           error = function(e) {
                               cat("    ERROR 3D UMAP:", conditionMessage(e), "\n")
                               NULL })
        if (is.null(cds_2d))
            stop("2D UMAP failed for ", dataset_tag, call. = FALSE)

        embed_label <- sprintf(
            "Monocle3 UMAP on raw counts (%s) | %d gene(s), %d PCA dimension(s)",
            NORM_TAG, cnt$n_gene_kept, built$n_pca)

        # == STEP 4: coordinates ==============================================
        cat("\n  STEP 4 - UMAP coordinates\n")
        df_umap2d <- v2_merge_groups(v2_umap_coords(cds_2d), rc)
        # Unknown-subgroup count, returned for the guardrail ledger.
        n_unknown_group <- {
            .a <- attr(df_umap2d, "n_unknown_group")
            if (is.null(.a)) NA_integer_ else as.integer(.a)[1]
        }
        df_umap3d <- if (!is.null(cds_3d))
            v2_merge_groups(v2_umap_coords(cds_3d), rc) else NULL
        cat(sprintf("    2D: %d cell(s)\n", nrow(df_umap2d)))
        if (!is.null(df_umap3d))
            cat(sprintf("    3D: %d cell(s)\n", nrow(df_umap3d)))

        has_3d <- !is.null(df_umap3d) &&
            all(c("UMAP_1", "UMAP_2", "UMAP_3") %in% names(df_umap3d))

        # == STEP 5: figure 1, four CMRF classes ==============================
        cat("\n  STEP 5 - Figures\n")
        # The two UMAP2D figures are drawn only if V2_MAKE_2D_FIGS is TRUE.
        if (!isTRUE(V2_MAKE_2D_FIGS))
            cat("    UMAP2D figures skipped (V2_MAKE_2D_FIGS = FALSE)\n")
        sub_txt <- sprintf("%s | %s\nOrigRef subgroups from %s",
                           ARM_TAG, embed_label, V2_ORIGREF_GROUP_COL)

        if (isTRUE(V2_MAKE_2D_FIGS)) {
        df_cc <- df_umap2d
        df_cc$cell_class_f <- factor(df_cc$cell_class, levels = CMRF_CLASS_LEVS)
        p_cc <- ggplot(df_cc, aes(x = UMAP_1, y = UMAP_2,
                                  colour = cell_class_f)) +
            geom_point(size = 0.9, alpha = 0.65, shape = 16) +
            scale_colour_manual(values = CMRF_CLASS_COLS, name = "CMRF class",
                                drop = FALSE, limits = CMRF_CLASS_LEVS,
                                guide = guide_legend(
                                    override.aes = list(
                                        size  = rep(2.4, length(CMRF_CLASS_LEVS)),
                                        alpha = rep(0.85, length(CMRF_CLASS_LEVS))))) +
            labs(title = sprintf("Validation2 - CMRF classes in UMAP space: %s",
                                 dataset_tag),
                 subtitle = sub_txt, x = "UMAP 1", y = "UMAP 2") +
            theme_v2() +
            theme(aspect.ratio = 1)
        save_pp(p_cc, file.path(OUT_DIR, paste0(FILE_PREFIX, "_", dataset_tag,
                                                "_UMAP2D_CellClass")))
        show_plot(p_cc)
        cat("    UMAP2D_CellClass saved\n")
        }

        # == STEP 5b: figure 2, Part 10e five-group colouring =================
        # Only plot_group5 is mapped to colour.
        if (isTRUE(V2_MAKE_2D_FIGS)) {
        df_nonsel   <- df_umap2d[df_umap2d$plot_group == "Non-selected", ]
        df_origonly <- df_umap2d[df_umap2d$cell_class == "Original only", ]
        df_origref  <- df_umap2d[df_umap2d$cell_class == "Original + Refined", ]
        df_recov    <- df_umap2d[df_umap2d$cell_class == "Recovered", ]

        p_grp <- ggplot() +
            geom_point(data = df_nonsel,
                       aes(x = UMAP_1, y = UMAP_2, colour = plot_group5),
                       size = 0.45, alpha = 0.45, shape = 16) +
            geom_point(data = df_origonly,
                       aes(x = UMAP_1, y = UMAP_2),
                       colour = "white", fill = ORIGONLY_COL,
                       shape = 25, size = 2.2, alpha = 0.88, stroke = 0.5) +
            geom_point(data = df_origref,
                       aes(x = UMAP_1, y = UMAP_2, colour = plot_group5),
                       size = 1.35, alpha = 0.58, shape = 16) +
            geom_point(data = df_recov,
                       aes(x = UMAP_1, y = UMAP_2),
                       colour = "white", fill = RECOV_COL_04,
                       shape = 23, size = 2.0, alpha = 0.93, stroke = 0.2) +
            scale_colour_manual(
                values = SCALE_COLS_VAL, name = "Cell group", drop = FALSE,
                limits = GROUP5_LEVS,
                # One override entry per level in limits.
                guide  = guide_legend(override.aes = list(
                    size  = c(1.8, rep(2.2, length(GROUP5_LEVS) - 1L)),
                    shape = rep(16L, length(GROUP5_LEVS)),
                    alpha = c(0.65, rep(0.75, length(GROUP5_LEVS) - 1L))))) +
            labs(title = sprintf("Validation2 - IP x TREP groups in UMAP space: %s",
                                 dataset_tag),
                 subtitle = paste0(sub_txt,
                                   "\nOriginal only = tan inverted triangle | ",
                                   "Recovered = violet diamond"),
                 x = "UMAP 1", y = "UMAP 2") +
            theme_v2() +
            theme(aspect.ratio = 1)
        save_pp(p_grp, file.path(OUT_DIR, paste0(FILE_PREFIX, "_", dataset_tag,
                                                 "_UMAP2D_GroupStyle")))
        show_plot(p_grp)
        cat("    UMAP2D_GroupStyle saved\n")
        }

        # == STEP 6: 3D projections ===========================================
        if (has_3d) {
            n3d <- nrow(df_umap3d)
            ps3 <- v2_proj_pt(n3d)
            cat(sprintf("    Projection point size %.2f for %d cell(s)\n",
                        ps3, n3d))

            # Original+Refined groups in V2_PROJ_GROUPS and Recovered only.
            make_proj <- function(df, xcol, ycol, panel_title) {
                d_or <- df[df$cell_class == "Original + Refined" &
                               as.character(df$plot_group) %in% V2_PROJ_GROUPS,
                           , drop = FALSE]
                d_rc <- df[df$cell_class == "Recovered", , drop = FALSE]
                ggplot() +
                    geom_point(data = d_or,
                               aes(x = .data[[xcol]], y = .data[[ycol]],
                                   colour = plot_group5),
                               size = ps3, alpha = V2_PROJ_ALPHA, shape = 16) +
                    geom_point(data = d_rc,
                               aes(x = .data[[xcol]], y = .data[[ycol]]),
                               colour = "white", fill = RECOV_COL_04,
                               shape = 23, size = ps3 * V2_PROJ_BIG_MULT,
                               alpha = V2_PROJ_RECOV_A,
                               stroke = V2_PROJ_STROKE) +
                    scale_colour_manual(values = SCALE_COLS_VAL,
                                        limits = GROUP5_LEVS,
                                        drop = FALSE, guide = "none") +
                    labs(title = panel_title,
                         x = gsub("UMAP_", "UMAP ", xcol),
                         y = gsub("UMAP_", "UMAP ", ycol)) +
                    theme_v2(base_size = 9) +
                    theme(aspect.ratio = 1, legend.position = "none")
            }

            g3 <- gridExtra::arrangeGrob(
                make_proj(df_umap3d, "UMAP_1", "UMAP_2", "3D UMAP: U1-U2"),
                make_proj(df_umap3d, "UMAP_1", "UMAP_3", "3D UMAP: U1-U3"),
                make_proj(df_umap3d, "UMAP_2", "UMAP_3", "3D UMAP: U2-U3"),
                nrow = 1,
                top = grid::textGrob(
                    sprintf(paste0("Validation2 3D UMAP projections: %s [%s]",
                                   "  --  %s and Recovered only"),
                            dataset_tag, ARM_TAG,
                            paste(V2_PROJ_GROUPS, collapse = ", ")),
                    gp = grid::gpar(fontface = "bold", fontsize = 12),
                    x = 0.01, hjust = 0))
            save_grob(g3, file.path(OUT_DIR, paste0(FILE_PREFIX, "_", dataset_tag,
                                                    "_UMAP3D_Projections")),
                      pdf_w = V2_PROJ_PDF_W, pdf_h = V2_PROJ_PDF_H,
                      png_w = V2_PROJ_PNG_W, png_h = V2_PROJ_PNG_H,
                      res   = V2_PROJ_PNG_RES)
            show_grob(g3)
            cat("    UMAP3D_Projections saved and drawn\n")

            # == One panel per group, common grey backdrop ====================
            if (isTRUE(V2_PROJ_BY_GROUP)) {
                # Panel order.
                .pg_lev <- c("Pre-IP TREP", "Recovered", "Post-IP TREP",
                             "Post-IP non-TREP", "Pre-IP non-TREP",
                             "Original only")
                d_mod <- df_umap3d[df_umap3d$cell_class %in%
                                       c("Original + Refined", "Recovered",
                                         "Original only"), , drop = FALSE]
                d_mod$panel_grp <- ifelse(
                    d_mod$cell_class == "Recovered", "Recovered",
                    ifelse(d_mod$cell_class == "Original only",
                           "Original only", as.character(d_mod$plot_group)))
                d_mod <- d_mod[d_mod$panel_grp %in% .pg_lev, , drop = FALSE]
                # Empty groups keep an empty panel.
                d_mod$panel_grp <- factor(d_mod$panel_grp, levels = .pg_lev)
                # The backdrop has no panel_grp, so it repeats in every facet.
                d_bg <- d_mod[, setdiff(names(d_mod), "panel_grp"), drop = FALSE]
                .pg_cols <- c(FOUR_COLS,
                              "Recovered"     = unname(RECOV_COL_04),
                              "Original only" = unname(ORIGONLY_COL))[.pg_lev]
                .n_lab <- table(d_mod$panel_grp)
                .pg_ttl <- sprintf("%s (%d)", .pg_lev, as.integer(.n_lab))
                levels(d_mod$panel_grp) <- .pg_ttl
                # Three layers so Recovered keeps its diamond and Original-only its triangle.
                .circ <- d_mod[!(d_mod$cell_class %in%
                                     c("Recovered", "Original only")), ,
                               drop = FALSE]
                .diam <- d_mod[d_mod$cell_class == "Recovered", , drop = FALSE]
                .tri  <- d_mod[d_mod$cell_class == "Original only", ,
                               drop = FALSE]
                p_bg <- ggplot() +
                    geom_point(data = d_bg,
                               aes(x = UMAP_1, y = UMAP_2),
                               colour = "gray86", size = ps3, alpha = 0.55,
                               shape = 16) +
                    geom_point(data = .circ,
                               aes(x = UMAP_1, y = UMAP_2, colour = panel_grp),
                               size = ps3 * 1.5, alpha = V2_PROJ_ALPHA,
                               shape = 16) +
                    geom_point(data = .tri,
                               aes(x = UMAP_1, y = UMAP_2),
                               colour = "white", fill = ORIGONLY_COL,
                               shape = 25, size = ps3 * V2_PROJ_BIG_MULT,
                               alpha = V2_PROJ_OO_A, stroke = V2_PROJ_STROKE) +
                    geom_point(data = .diam,
                               aes(x = UMAP_1, y = UMAP_2),
                               colour = "white", fill = RECOV_COL_04,
                               shape = 23, size = ps3 * V2_PROJ_BIG_MULT,
                               alpha = V2_PROJ_RECOV_A,
                               stroke = V2_PROJ_STROKE) +
                    scale_colour_manual(
                        values = setNames(unname(.pg_cols), .pg_ttl),
                        drop = FALSE, guide = "none") +
                    facet_wrap(~ panel_grp, nrow = 2, drop = FALSE) +
                    labs(title = sprintf(paste0("Validation2 3D UMAP by group",
                                                " (U1-U2): %s [%s]"),
                                         dataset_tag, ARM_TAG),
                         subtitle = paste0("Grey = every other module cell, ",
                                           "repeated in each panel | ",
                                           "Recovered = diamond, ",
                                           "Original-only = triangle"),
                         x = "UMAP 1", y = "UMAP 2") +
                    theme_v2(base_size = 9) +
                    theme(aspect.ratio = 1, legend.position = "none")
                save_pp(p_bg,
                        file.path(OUT_DIR, paste0(FILE_PREFIX, "_", dataset_tag,
                                                  "_UMAP3D_ByGroup")),
                        pdf_w = 13, pdf_h = 9,
                        png_w = 3120, png_h = 2160, res = 240)
                show_plot(p_bg)
                cat("    UMAP3D_ByGroup saved\n")
            }

            if (isTRUE(V2_WRITE_HTML) && .v2_plotly_ok && .v2_widgets_ok) {
                html_path <- file.path(
                    OUT_DIR, paste0(FILE_PREFIX, "_", dataset_tag,
                                    "_UMAP3D_Interactive.html"))
                ok <- tryCatch({
                    fig <- v2_plotly3d(
                        df_umap3d,
                        sprintf("Validation2 3D UMAP: %s [%s]<br><sup>%s</sup>",
                                dataset_tag, ARM_TAG, embed_label))
                    htmlwidgets::saveWidget(fig, file = html_path,
                                            selfcontained = TRUE)
                    TRUE
                }, error = function(e) {
                    cat("    selfcontained HTML failed (", conditionMessage(e),
                        "); retrying unpacked\n", sep = "")
                    tryCatch({
                        fig <- v2_plotly3d(
                            df_umap3d,
                            sprintf("Validation2 3D UMAP: %s [%s]",
                                    dataset_tag, ARM_TAG))
                        htmlwidgets::saveWidget(fig, file = html_path,
                                                selfcontained = FALSE)
                        TRUE
                    }, error = function(e2) {
                        cat("    interactive HTML skipped:",
                            conditionMessage(e2), "\n"); FALSE })
                })
                if (isTRUE(ok)) cat("    UMAP3D_Interactive.html saved\n")
            } else if (isTRUE(V2_WRITE_HTML)) {
                cat("    interactive HTML skipped (plotly or htmlwidgets absent)\n")
            }
        } else {
            cat("    3D UMAP unavailable; projections and HTML skipped\n")
        }

        # == STEP 7: donor composition ========================================
        cat("\n  STEP 7 - Donor composition of the four classes\n")
        don_tab <- as.data.frame(
            table(cell_class = rc$cell_class, donor = rc$donor),
            stringsAsFactors = FALSE)
        don_tab <- don_tab[don_tab$Freq > 0, , drop = FALSE]
        names(don_tab)[names(don_tab) == "Freq"] <- "n_cells"
        don_tab$dataset <- dataset_tag
        don_tab <- don_tab[, c("dataset", "cell_class", "donor", "n_cells")]
        val_tabs[["Donor_By_Class"]] <- don_tab

        donor_share <- function(cls) {
            d <- rc$donor[rc$cell_class == cls & !is.na(rc$donor)]
            if (length(d) == 0L) return(NULL)
            p <- table(d) / length(d)
            p
        }
        p_ref <- donor_share("Original + Refined")

        v2_tvd <- function(pa, pb) {
            lev <- union(names(pa), names(pb))
            a <- as.numeric(pa[lev]); a[is.na(a)] <- 0
            b <- as.numeric(pb[lev]); b[is.na(b)] <- 0
            0.5 * sum(abs(a - b))
        }

        # TVD against the 95th percentile of its own null (cell_class permuted,
        # donor fixed).
        .tvd_keep  <- !is.na(rc$donor)
        .tvd_don   <- rc$donor[.tvd_keep]
        .tvd_class <- rc$cell_class[.tvd_keep]
        tvd_null_p95 <- function(cls) {
            if (identical(cls, "Original + Refined")) return(NA_real_)
            if (sum(.tvd_class == cls) < 2L ||
                sum(.tvd_class == "Original + Refined") < 2L) return(NA_real_)
            set.seed(V2_SEED)
            v <- vapply(seq_len(V2_TVD_PERM), function(i) {
                perm <- sample(.tvd_class)
                a <- table(.tvd_don[perm == cls])
                b <- table(.tvd_don[perm == "Original + Refined"])
                v2_tvd(a / sum(a), b / sum(b))
            }, numeric(1L))
            round(as.numeric(quantile(v, 0.95, names = FALSE)), 4)
        }

        don_sum <- do.call(rbind, lapply(CMRF_CLASS_LEVS, function(cls) {
            p <- donor_share(cls)
            n <- sum(rc$cell_class == cls)
            if (is.null(p))
                return(data.frame(dataset = dataset_tag, cell_class = cls,
                                  n_cells = n, n_donors = 0L,
                                  top_donor_share = NA_real_,
                                  TVD_vs_OrigRef = NA_real_,
                                  TVD_null_p95 = NA_real_,
                                  stringsAsFactors = FALSE))
            tvd <- if (is.null(p_ref)) NA_real_ else round(v2_tvd(p, p_ref), 4)
            data.frame(dataset = dataset_tag, cell_class = cls, n_cells = n,
                       n_donors = length(p),
                       top_donor_share = round(max(p), 4),
                       TVD_vs_OrigRef = tvd,
                       TVD_null_p95 = tvd_null_p95(cls),
                       stringsAsFactors = FALSE)
        }))
        val_tabs[["Donor_Composition"]] <- don_sum
        for (i in seq_len(nrow(don_sum)))
            cat(sprintf("    %-22s n=%5d | donors=%3d | top share=%s | TVD=%s (null p95 %s)\n",
                        don_sum$cell_class[i], don_sum$n_cells[i],
                        don_sum$n_donors[i],
                        if (is.na(don_sum$top_donor_share[i])) "NA"
                        else sprintf("%.3f", don_sum$top_donor_share[i]),
                        if (is.na(don_sum$TVD_vs_OrigRef[i])) "NA"
                        else sprintf("%.3f", don_sum$TVD_vs_OrigRef[i]),
                        if (is.na(don_sum$TVD_null_p95[i])) "NA"
                        else sprintf("%.3f", don_sum$TVD_null_p95[i])))

        # == STEP 8: statistics ===============================================
        cat("\n  STEP 8 - Equivalence, dispersion and local mixing\n")

        df_cmp   <- if (has_3d) df_umap3d else df_umap2d
        umap_axes <- if (has_3d) c("UMAP_1", "UMAP_2", "UMAP_3")
                     else c("UMAP_1", "UMAP_2")
        cat(sprintf("    Using the %s embedding (%s)\n",
                    if (has_3d) "3D" else "2D",
                    paste(umap_axes, collapse = ", ")))

        # group_ref travels with the coordinates.
        origref    <- df_cmp[df_cmp$cell_class == "Original + Refined", ,
                             drop = FALSE]
        origref    <- origref[!is.na(origref$group_ref) &
                                  origref$group_ref %in% names(FOUR_COLS), ,
                              drop = FALSE]
        recov_umap <- df_cmp[df_cmp$cell_class == "Recovered", , drop = FALSE]
        oo_umap    <- df_cmp[df_cmp$cell_class == "Original only", , drop = FALSE]

        cat(sprintf("    OrigRef with a subgroup label: %d | Recovered: %d | Original only: %d\n",
                    nrow(origref), nrow(recov_umap), nrow(oo_umap)))
        for (g in names(FOUR_COLS))
            cat(sprintf("      %-22s %d\n", g, sum(origref$group_ref == g)))

        focal_sets <- list("Recovered" = recov_umap, "Original only" = oo_umap)
        ref_names  <- names(FOUR_COLS)

        # Donor-unit tests need a donor on at least some cells.
        donor_unit_ok <- isTRUE(V2_DONOR_UNIT_TESTS) &&
            "donor" %in% names(df_cmp) && any(!is.na(df_cmp$donor))
        cat(sprintf("    Donor-unit tests: %s\n",
                    if (donor_unit_ok)
                        sprintf(paste0("ON (>= %d donor(s) per group, >= %d",
                                       " cell(s) per donor)"),
                                V2_MIN_DONORS, V2_MIN_CELLS_PER_DONOR)
                    else "OFF"))

        tost_rows <- list(); perm_rows <- list(); ilisi_rows <- list()
        nbr_rows  <- list(); kbet_rows <- list()

        for (focal_name in names(focal_sets)) {
            focal_data <- focal_sets[[focal_name]]
            if (nrow(focal_data) < V2_MIN_GROUP) {
                cat(sprintf("    %s: %d cell(s), below V2_MIN_GROUP; skipped\n",
                            focal_name, nrow(focal_data)))
                next
            }
            focal_mat <- as.matrix(focal_data[, umap_axes, drop = FALSE])

            for (ref_grp in ref_names) {
                ref_cells <- origref[origref$group_ref == ref_grp, , drop = FALSE]
                if (nrow(ref_cells) < V2_MIN_GROUP) next
                ref_mat <- as.matrix(ref_cells[, umap_axes, drop = FALSE])

                # Donor-level means: one row per donor.
                focal_dm <- if (donor_unit_ok)
                    v2_donor_means(focal_data, umap_axes) else NULL
                ref_dm   <- if (donor_unit_ok)
                    v2_donor_means(ref_cells,  umap_axes) else NULL
                # NA when the donor unit is off.
                n_don_f  <- if (!donor_unit_ok) NA_integer_
                            else if (is.null(focal_dm)) 0L else nrow(focal_dm)
                n_don_r  <- if (!donor_unit_ok) NA_integer_
                            else if (is.null(ref_dm))   0L else nrow(ref_dm)
                donor_pair_ok <- !is.na(n_don_f) && !is.na(n_don_r) &&
                    n_don_f >= V2_MIN_DONORS && n_don_r >= V2_MIN_DONORS

                # TOST, one axis at a time, at each available unit
                for (ax in umap_axes) {
                    tost_units <- list(list(u = "cell",
                                            f = focal_data[[ax]],
                                            r = ref_cells[[ax]]))
                    if (donor_pair_ok)
                        tost_units[[2L]] <- list(u = "donor",
                                                 f = focal_dm[[ax]],
                                                 r = ref_dm[[ax]])
                    for (tu in tost_units) {
                        rt <- v2_run_tost(tu$f, tu$r)
                        tost_rows[[length(tost_rows) + 1L]] <- data.frame(
                            dataset = dataset_tag, unit = tu$u,
                            focal = focal_name, reference = ref_grp, axis = ax,
                            delta_frac = rt$delta_frac, delta = rt$delta,
                            pooled_sd = rt$pooled_sd, mean_diff = rt$mean_diff,
                            p_tost = rt$p_tost, equiv = rt$equiv,
                            p_lower = rt$p_lower, p_upper = rt$p_upper,
                            n_focal = rt$n1, n_ref = rt$n2,
                            n_donor_focal = n_don_f, n_donor_ref = n_don_r,
                            SW_p_focal = rt$sw_p1, SW_p_ref = rt$sw_p2,
                            stringsAsFactors = FALSE)
                        if (identical(focal_name, "Recovered"))
                            cat(sprintf("      TOST %-5s %-18s %s: p=%s delta=%s (%s)\n",
                                        tu$u, ref_grp, gsub("UMAP_", "U", ax),
                                        if (is.na(rt$p_tost)) "NA" else sprintf("%.4f", rt$p_tost),
                                        if (is.na(rt$delta)) "NA" else sprintf("%.3f", rt$delta),
                                        if (isTRUE(rt$equiv)) "equivalent" else "not equivalent"))
                    }
                }

                # PERMANOVA and PERMDISP at each unit; the donor unit is not subsampled.
                perm_units <- list(list(u = "cell", f = focal_mat, r = ref_mat,
                                        mx = V2_PERM_MAX_N))
                if (donor_pair_ok)
                    perm_units[[2L]] <- list(
                        u  = "donor",
                        f  = as.matrix(focal_dm[, umap_axes, drop = FALSE]),
                        r  = as.matrix(ref_dm[,   umap_axes, drop = FALSE]),
                        mx = .Machine$integer.max)
                for (pu in perm_units) {
                    rp <- v2_permanova_pair(pu$f, pu$r, max_n = pu$mx)
                    perm_rows[[length(perm_rows) + 1L]] <- data.frame(
                        dataset = dataset_tag, unit = pu$u,
                        focal = focal_name, reference = ref_grp,
                        n_axes = length(umap_axes),
                        n_focal = rp$n_focal, n_ref = rp$n_ref,
                        n_focal_used = rp$n_focal_used,
                        n_ref_used = rp$n_ref_used,
                        subsampled = rp$subsampled,
                        R2 = rp$R2, F_stat = rp$F_stat, p_value = rp$p_value,
                        disp_F = rp$disp_F, disp_p = rp$disp_p,
                        avg_dist_focal = rp$avg_dist_focal,
                        avg_dist_ref = rp$avg_dist_ref,
                        n_perm = rp$n_perm, note = rp$note,
                        stringsAsFactors = FALSE)
                    if (identical(focal_name, "Recovered"))
                        cat(sprintf("      PERMANOVA %-5s %-18s R2=%s p=%s | PERMDISP F=%s p=%s%s\n",
                                    pu$u, ref_grp,
                                    if (is.na(rp$R2)) "NA" else sprintf("%.4f", rp$R2),
                                    if (is.na(rp$p_value)) "NA" else sprintf("%.4f", rp$p_value),
                                    if (is.na(rp$disp_F)) "NA" else sprintf("%.2f", rp$disp_F),
                                    if (is.na(rp$disp_p)) "NA" else sprintf("%.4f", rp$disp_p),
                                    if (isTRUE(rp$subsampled)) "  [subsampled]" else ""))
                }

                # iLISI on the combined pair
                comb_mat <- rbind(focal_mat, ref_mat)
                comb_lab <- c(rep(focal_name, nrow(focal_mat)),
                              rep(ref_grp,    nrow(ref_mat)))
                li <- tryCatch(v2_compute_ilisi(comb_mat, comb_lab),
                               error = function(e) {
                                   cat("      iLISI error:", conditionMessage(e), "\n")
                                   NULL })
                if (!is.null(li) && !is.null(li$by_group)) {
                    for (g in seq_len(nrow(li$by_group))) {
                        bg <- li$by_group[g, ]
                        ilisi_rows[[length(ilisi_rows) + 1L]] <- data.frame(
                            dataset = dataset_tag, focal = focal_name,
                            reference = ref_grp, cell_group = bg$group,
                            n = bg$n, mean_iLISI = bg$mean_iLISI,
                            median_iLISI = bg$median_iLISI,
                            pct_gt_1p5 = bg$pct_gt_1p5,
                            overall_mean = li$mean, overall_median = li$median,
                            max_iLISI = li$n_groups, k = li$k,
                            stringsAsFactors = FALSE)
                    }
                    if (identical(focal_name, "Recovered"))
                        cat(sprintf("      iLISI     %-18s mean=%.3f median=%.3f (max %d)\n",
                                    ref_grp, li$mean, li$median, li$n_groups))
                }

                # Neighbour proportions and enrichment
                nm <- tryCatch(v2_neighbour_metrics(focal_mat, ref_mat,
                                                    focal_name, ref_grp),
                               error = function(e) NULL)
                if (!is.null(nm)) {
                    nbr_rows[[length(nbr_rows) + 1L]] <- data.frame(
                        dataset = dataset_tag, focal = focal_name,
                        reference = ref_grp, n_focal = nrow(focal_mat),
                        n_ref = nrow(ref_mat), k = nm$k,
                        global_prop_ref = nm$global_prop,
                        mean_ref_frac = nm$mean_ref_frac,
                        pct_gt_40 = nm$pct_gt_40, pct_gt_50 = nm$pct_gt_50,
                        mean_enrichment = nm$mean_enrichment,
                        median_enrichment = nm$median_enrichment,
                        stringsAsFactors = FALSE)
                    if (identical(focal_name, "Recovered"))
                        cat(sprintf("      kNN       %-18s ref_frac=%.3f global=%.3f enrich=%.2fx\n",
                                    ref_grp, nm$mean_ref_frac, nm$global_prop,
                                    nm$mean_enrichment))
                }

                # kBET
                kb <- tryCatch(v2_compute_kbet(focal_mat, ref_mat,
                                               focal_name, ref_grp),
                               error = function(e) NULL)
                if (!is.null(kb)) {
                    kbet_rows[[length(kbet_rows) + 1L]] <- data.frame(
                        dataset = dataset_tag, focal = focal_name,
                        reference = ref_grp, n_focal = nrow(focal_mat),
                        n_ref = nrow(ref_mat), k = kb$k, alpha = kb$alpha,
                        rejection_rate = kb$rejection_rate,
                        mean_pval = kb$mean_pval, median_pval = kb$median_pval,
                        n_tested = kb$n_tested, n_undefined = kb$n_undefined,
                        note = kb$note, stringsAsFactors = FALSE)
                    if (identical(focal_name, "Recovered"))
                        cat(sprintf("      kBET      %-18s reject=%s%s\n",
                                    ref_grp,
                                    if (is.na(kb$rejection_rate)) "NA"
                                    else sprintf("%.3f", kb$rejection_rate),
                                    if (nzchar(kb$note))
                                        paste0("  [", kb$note, "]") else ""))
                }
            }
        }

        # BH within unit.
        bh_by_unit <- function(p_vec, unit_vec) {
            out <- rep(NA_real_, length(p_vec))
            for (u in unique(unit_vec)) {
                i <- unit_vec == u
                out[i] <- round(p.adjust(p_vec[i], method = "BH"), 6)
            }
            out
        }

        if (length(tost_rows) > 0) {
            td <- do.call(rbind, tost_rows)
            td$p_tost_BH <- bh_by_unit(td$p_tost, td$unit)
            td$equiv_BH  <- !is.na(td$p_tost_BH) & td$p_tost_BH < 0.05
            val_tabs[["TOST_Equivalence"]] <- td
            for (u in unique(td$unit)) {
                i <- td$unit == u
                cat(sprintf("    TOST [%s]: %d test(s), %d equivalent unadjusted, %d after BH\n",
                            u, sum(i), sum(td$equiv[i], na.rm = TRUE),
                            sum(td$equiv_BH[i], na.rm = TRUE)))
            }
        }
        if (length(perm_rows) > 0) {
            pd <- do.call(rbind, perm_rows)
            pd$p_BH_FDR  <- bh_by_unit(pd$p_value, pd$unit)
            pd$disp_p_BH <- bh_by_unit(pd$disp_p,  pd$unit)
            pd$inseparable_BH <- !is.na(pd$p_BH_FDR) & pd$p_BH_FDR >= 0.05 &
                !is.na(pd$disp_p_BH) & pd$disp_p_BH >= 0.05
            val_tabs[["PERMANOVA_PERMDISP"]] <- pd
            for (u in unique(pd$unit)) {
                i <- pd$unit == u
                cat(sprintf(paste0("    PERMANOVA [%s]: %d test(s) | centroid ns %d",
                                   " | dispersion ns %d | both ns %d\n"),
                            u, sum(i),
                            sum(i & !is.na(pd$p_BH_FDR) & pd$p_BH_FDR >= 0.05),
                            sum(i & !is.na(pd$disp_p_BH) & pd$disp_p_BH >= 0.05),
                            sum(pd$inseparable_BH[i], na.rm = TRUE)))
            }
        }
        if (length(ilisi_rows) > 0)
            val_tabs[["iLISI"]] <- do.call(rbind, ilisi_rows)
        if (length(nbr_rows) > 0)
            val_tabs[["Neighbour_Proportions"]] <- do.call(rbind, nbr_rows)
        if (length(kbet_rows) > 0)
            val_tabs[["kBET"]] <- do.call(rbind, kbet_rows)

        # == STEP 9: overview =================================================
        n_by_class <- vapply(CMRF_CLASS_LEVS,
                             function(cl) sum(rc$cell_class == cl), integer(1L))
        # == STEP 9: Figure 11 statistics =====================================
        cat("\n  STEP 9 - Figure 11 statistics\n")
        fig11_desc <- NULL; fig11_omni <- NULL; fig11_ph <- NULL
        .proj_note <- "not computed"; .n_blocks_proj <- NA_integer_
        .n_blocks_res <- NA_integer_

        # Both statistics are built on df_cmp.
        .base <- df_cmp
        .grp_proj <- ifelse(!is.na(.base$cell_class) &
                                .base$cell_class == "Recovered", "Recovered",
                            as.character(.base$plot_group))
        .grp_res  <- as.character(.base$cell_class)

        .add <- function(stat, space, res) {
            if (is.null(res)) return(invisible(NULL))
            d <- res$desc; d$statistic <- stat; d$space <- space
            d$unit <- res$omnibus$unit
            fig11_desc <<- rbind(fig11_desc, d)
            o <- res$omnibus; o$statistic <- stat; o$space <- space
            fig11_omni <<- rbind(fig11_omni, o)
            if (!is.null(res$ph)) {
                p <- res$ph; p$statistic <- stat; p$space <- space
                fig11_ph <<- rbind(fig11_ph, p)
            }
        }

        # == Transcriptome space: Pre -> Post projection ======================
        pr <- v2_prepost_proj(pca_mat, .base$barcode, .grp_proj,
                              V2_PROJ_LEVELS[1L], V2_PROJ_LEVELS[3L])
        if (is.null(pr)) {
            cat("    Pre->Post projection SKIPPED (no PCA, or an arm below",
                V2_STAT_MIN_GROUP, "cells)\n")
        } else {
            .base$prepost     <- pr$proj
            .base$prepost_off <- pr$off
            .proj_note <- sprintf(
                "axis from %d Pre-IP TREP to %d Post-IP TREP cell(s), centroid separation %.3f",
                pr$n_pre, pr$n_post, pr$sep)
            cat("    ", .proj_note, "\n", sep = "")
            for (u in c("donor", "cell")) {
                r <- v2_protocol(.base$prepost, .grp_proj,
                                 block = if (u == "donor") .base$donor else NULL,
                                 levs = V2_PROJ_LEVELS, unit = u)
                if (u == "donor") .n_blocks_proj <- r$omnibus$n_unit
                .add("PrePost_projection", "transcriptome (PCA)", r)
                cat(sprintf("      %-5s omnibus %s p=%s | trend %s p=%s\n", u,
                            r$omnibus$omnibus_test, v2_p_fmt(r$omnibus$omnibus_p),
                            r$omnibus$trend_test, v2_p_fmt(r$omnibus$trend_p)))
                # Off-axis distance (omnibus only; no trend expected).
                ro <- v2_protocol(.base$prepost_off, .grp_proj,
                                  block = if (u == "donor") .base$donor else NULL,
                                  levs = V2_PROJ_LEVELS, unit = u)
                .add("PrePost_offaxis", "transcriptome (PCA)", ro)
                cat(sprintf("      %-5s off-axis omnibus %s p=%s\n", u,
                            ro$omnibus$omnibus_test,
                            v2_p_fmt(ro$omnibus$omnibus_p)))
            }
        }

        # == Copula space: the same projection on the SKCS axes ===============
        # (u_SK, v_SK) from Validation1's RDS, recomputed from the Step 4 sheet as a check.
        .sk_src <- "none"; .sk_dev <- NA_real_
        .uv_sk <- NULL
        if (isTRUE(V2_PROJ_COPULA)) {
            .rp_sk <- p10v2_v1_rds(dataset_tag)
            if (file.exists(.rp_sk)) {
                .v1 <- tryCatch(readRDS(.rp_sk), error = function(e) NULL)
                if (!is.null(.v1) &&
                    all(V2_COPULA_SK_FIELDS %in% names(.v1))) {
                    .us <- .v1[[V2_COPULA_SK_FIELDS[1L]]]
                    .vs <- .v1[[V2_COPULA_SK_FIELDS[2L]]]
                    if (!is.null(names(.us)) && !is.null(names(.vs))) {
                        .m <- cbind(u_sk = unname(.us[rc$barcode]),
                                    v_sk = unname(.vs[rc$barcode]))
                        rownames(.m) <- rc$barcode
                        .uv_sk <- .m[stats::complete.cases(.m), , drop = FALSE]
                        .sk_src <- "Validation1 RDS"
                    }
                }
            }
            # Recompute: average ranks over GOI-positive cells / (n + 1), as in Validation1.
            if (all(V2_COPULA_SK_RECALC %in% names(rc))) {
                .pos <- !is.na(rc$cell_class) & rc$cell_class != "GOI-negative"
                .n_c <- sum(.pos)
                .ur <- rep(NA_real_, nrow(rc)); .vr <- .ur
                .ur[.pos] <- rank(rc[[V2_COPULA_SK_RECALC[1L]]][.pos],
                                  ties.method = "average") / (.n_c + 1)
                .vr[.pos] <- rank(rc[[V2_COPULA_SK_RECALC[2L]]][.pos],
                                  ties.method = "average") / (.n_c + 1)
                .mr <- cbind(u_sk = .ur, v_sk = .vr)
                rownames(.mr) <- rc$barcode
                .mr <- .mr[stats::complete.cases(.mr), , drop = FALSE]
                if (is.null(.uv_sk)) {
                    .uv_sk  <- .mr
                    .sk_src <- "recomputed from the Step 4 sheet"
                } else {
                    .sh <- intersect(rownames(.uv_sk), rownames(.mr))
                    if (length(.sh) > 0L)
                        .sk_dev <- max(abs(.uv_sk[.sh, ] - .mr[.sh, ]))
                }
            }
        }
        if (!is.null(.uv_sk)) {
            prc <- v2_prepost_proj(.uv_sk, .base$barcode, .grp_proj,
                                   V2_PROJ_LEVELS[1L], V2_PROJ_LEVELS[3L])
            if (is.null(prc)) {
                cat("    Copula (u_SK, v_SK) projection SKIPPED\n")
            } else {
                .base$prepost_cop     <- prc$proj
                .base$prepost_cop_off <- prc$off
                .proj_note_cop <- sprintf(
                    "copula (u_SK, v_SK) axis from %s, centroid separation %.3f%s",
                    .sk_src, prc$sep,
                    if (is.finite(.sk_dev))
                        sprintf(", max |RDS - recomputed| = %.2e", .sk_dev)
                    else "")
                cat("    ", .proj_note_cop, "\n", sep = "")
                for (u in c("donor", "cell")) {
                    # Lateral displacement from the arm axis.
                    r <- v2_protocol(.base$prepost_cop, .grp_proj,
                                     block = if (u == "donor") .base$donor else NULL,
                                     levs = V2_PROJ_LEVELS, unit = u)
                    .add("PrePost_projection_copula",
                         "copula SK (axis-shape caveat)", r)
                    ro <- v2_protocol(.base$prepost_cop_off, .grp_proj,
                                      block = if (u == "donor") .base$donor else NULL,
                                      levs = V2_PROJ_LEVELS, unit = u)
                    .add("Lateral_displacement_copula",
                         "copula SK (method check)", ro)
                    cat(sprintf("      %-5s lateral omnibus %s p=%s\n", u,
                                ro$omnibus$omnibus_test,
                                v2_p_fmt(ro$omnibus$omnibus_p)))
                }

                # Bridge distance on the SK coordinates, and the SKCS cut each Recovered cell fails.
                if (isTRUE(V2_BRIDGE_ENABLE)) {
                    .bxy <- .uv_sk[match(.base$barcode, rownames(.uv_sk)), ,
                                   drop = FALSE]
                    .grp_br <- ifelse(!is.na(.base$cell_class) &
                                          .base$cell_class %in%
                                          c("Recovered", "Excluded"),
                                      .base$cell_class,
                                      as.character(.base$plot_group))
                    .bd <- v2_bridge_dist(.bxy, .grp_br,
                                          V2_PROJ_LEVELS[1L], V2_PROJ_LEVELS[3L])
                    if (is.null(.bd)) {
                        cat("    Bridge distance SKIPPED (an arm below 2 cells)\n")
                    } else {
                        .base$bridge <- .bd
                        for (u in c("donor", "cell")) {
                            rb <- v2_protocol(.base$bridge, .grp_br,
                                              block = if (u == "donor") .base$donor else NULL,
                                              levs = V2_BRIDGE_LEVELS, unit = u)
                            .add("Bridge_distance_copula",
                                 "copula SK (Recovered vs Excluded)", rb)
                            cat(sprintf("      %-5s bridge omnibus %s p=%s\n",
                                        u, rb$omnibus$omnibus_test,
                                        v2_p_fmt(rb$omnibus$omnibus_p)))
                        }
                    }
                    # Si only, Ki only or both.
                    if (all(c("spear_norm", "kend_contrib") %in% names(rc))) {
                        .si <- setNames(as.numeric(rc$spear_norm), rc$barcode)
                        .ki <- setNames(as.numeric(rc$kend_contrib), rc$barcode)
                        .rb <- .base$barcode[!is.na(.base$cell_class) &
                                                 .base$cell_class == "Recovered"]
                        .fs <- unname(.si[.rb]) <= V2_SI_NORM_CUT
                        .fk <- unname(.ki[.rb]) <= V2_KI_CUT
                        .ok <- !is.na(.fs) & !is.na(.fk)
                        .n  <- sum(.ok)
                        if (.n > 0L) {
                            .cut_tab <- data.frame(
                                dataset = dataset_tag,
                                n_recovered_scored = .n,
                                fail_Si_only = sum(.fs[.ok] & !.fk[.ok]),
                                fail_Ki_only = sum(!.fs[.ok] & .fk[.ok]),
                                fail_both    = sum(.fs[.ok] & .fk[.ok]),
                                fail_neither = sum(!.fs[.ok] & !.fk[.ok]),
                                Si_norm_cut  = V2_SI_NORM_CUT,
                                Ki_cut       = V2_KI_CUT,
                                stringsAsFactors = FALSE)
                            .cut_tab$pct_fail_Si_any <- round(
                                100 * (.cut_tab$fail_Si_only + .cut_tab$fail_both) / .n, 1)
                            .cut_tab$pct_fail_Ki_any <- round(
                                100 * (.cut_tab$fail_Ki_only + .cut_tab$fail_both) / .n, 1)
                            val_tabs[["Fig11_SKCS_Cuts"]] <- .cut_tab
                            cat(sprintf(paste0("    SKCS cut failed by the %d",
                                               " Recovered cell(s): Si only %d,",
                                               " Ki only %d, both %d, neither %d",
                                               " (Si-any %.1f%%, Ki-any %.1f%%)\n"),
                                        .n, .cut_tab$fail_Si_only,
                                        .cut_tab$fail_Ki_only, .cut_tab$fail_both,
                                        .cut_tab$fail_neither,
                                        .cut_tab$pct_fail_Si_any,
                                        .cut_tab$pct_fail_Ki_any))
                        }
                    }
                }
            }
        } else if (isTRUE(V2_PROJ_COPULA)) {
            cat("    Copula (u_SK, v_SK) projection SKIPPED: neither",
                paste(V2_COPULA_SK_FIELDS, collapse = "/"), "in the",
                "Validation1 RDS nor",
                paste(V2_COPULA_SK_RECALC, collapse = "/"),
                "in the Step 4 sheet\n")
        }

        # == Copula space: the spine residual =================================
        if (!V2_RESID_COL %in% names(rc)) {
            cat("    Spine residual SKIPPED (Validation1 wrote no",
                V2_RESID_COL, "column)\n")
        } else {
            .rmap <- setNames(as.numeric(rc[[V2_RESID_COL]]), rc$barcode)
            .base$resid <- unname(.rmap[.base$barcode])
            for (u in c("donor", "cell")) {
                r <- v2_protocol(.base$resid, .grp_res,
                                 block = if (u == "donor") .base$donor else NULL,
                                 levs = V2_RESID_LEVELS, unit = u)
                if (u == "donor") .n_blocks_res <- r$omnibus$n_unit
                .add("Spine_residual", "copula (method check)", r)
                cat(sprintf("      %-5s omnibus %s p=%s | trend %s p=%s\n", u,
                            r$omnibus$omnibus_test, v2_p_fmt(r$omnibus$omnibus_p),
                            r$omnibus$trend_test, v2_p_fmt(r$omnibus$trend_p)))
            }
        }

        # BH across the two datasets is applied in Section D.
        if (!is.null(fig11_desc)) val_tabs[["Fig11_Descriptives"]] <- fig11_desc
        if (!is.null(fig11_omni)) val_tabs[["Fig11_Omnibus"]]      <- fig11_omni
        if (!is.null(fig11_ph))   val_tabs[["Fig11_PostHoc"]]      <- fig11_ph

        val_tabs[["Fig11_Methods"]] <- data.frame(
            item = c("normality", "two groups", ">2 groups independent",
                     ">2 groups repeated measures", "ordered trend independent",
                     "ordered trend repeated measures",
                     "post hoc independent", "post hoc repeated measures",
                     "BH family pairwise", "BH family omnibus and trend",
                     "significance", "projection axis",
                     "complete donor blocks, projection",
                     "complete donor blocks, residual",
                     "copula projection axis"),
            value = c(
                sprintf("Shapiro-Wilk, alpha %.2f, subsampled above %d",
                        V2_STAT_ALPHA, V2_SW_MAX_N),
                "t-test or Mann-Whitney U",
                "ANOVA with Holm-Sidak, or Kruskal-Wallis with Dunn",
                "repeated-measures ANOVA, or Friedman (complete blocks)",
                "Jonckheere-Terpstra",
                "Page's L",
                "Holm-Sidak (normal) or Dunn (non-normal)",
                paste("paired t-test (normal) or paired Wilcoxon signed-rank",
                      "(non-normal)"),
                "m = 3, within one dataset x unit x statistic",
                "m = 2, across the two datasets (applied at the arm level)",
                "BH FDR < 0.05",
                .proj_note,
                as.character(.n_blocks_proj), as.character(.n_blocks_res),
                if (exists(".proj_note_cop", inherits = FALSE))
                    .proj_note_cop else "not computed"),
            stringsAsFactors = FALSE)

        val_tabs[["Overview"]] <- data.frame(
            dataset            = dataset_tag,
            arm                = ARM_TAG,
            label_source       = lab$source,
            embed_source       = V2_EMBED_SOURCE,
            embedding          = embed_label,
            n_rows_in_sheet    = lab$n_all,
            n_dropped_GOI_neg  = lab$n_dropped,
            n_cells_embedded   = n_cell,
            n_genes_input      = cnt$n_gene_all,
            n_genes_kept       = cnt$n_gene_kept,
            n_pca              = built$n_pca,
            n_orig_ref         = unname(n_by_class["Original + Refined"]),
            n_recovered        = unname(n_by_class["Recovered"]),
            n_orig_only        = unname(n_by_class["Original only"]),
            n_excluded         = unname(n_by_class["Excluded"]),
            origref_group_col  = V2_ORIGREF_GROUP_COL,
            donor_unit_tests   = donor_unit_ok,
            min_donors         = V2_MIN_DONORS,
            min_cells_per_donor = V2_MIN_CELLS_PER_DONOR,
            n_cells_no_donor   = n_no_don,
            # The 2D embedding is always built.
            umap2d_built       = TRUE,
            umap2d_figs        = isTRUE(V2_MAKE_2D_FIGS),
            umap3d_built       = has_3d,
            test_axes          = paste(umap_axes, collapse = "+"),
            stringsAsFactors   = FALSE)

        # The CDS objects are not returned.
        invisible(list(
            dataset_tag = dataset_tag,
            n_cells = n_cell, n_by_class = n_by_class,
            n_genes_kept = cnt$n_gene_kept, n_pca = built$n_pca,
            raw_integer_like = cnt$raw_integer_like, raw_max = cnt$raw_max,
            label_source = lab$source, embed_label = embed_label,
            has_3d = has_3d, umap_axes = umap_axes,
            n_unknown_group = n_unknown_group,
            n_origref_labelled = nrow(origref),
            n_recovered_tested = nrow(recov_umap),
            n_orig_only_tested = nrow(oo_umap),
            donor_summary = don_sum, n_no_donor = n_no_don,
            df_umap2d = df_umap2d, df_umap3d = df_umap3d,
            val_tabs = val_tabs))
    }


    # =========================================================================
    # SECTION C: Execute
    # =========================================================================

    all_results <- vector("list", length(DATASETS))
    names(all_results) <- DATASETS
    for (ds in DATASETS) {
        # Single-bracket assignment keeps a NULL result in the list.
        all_results[ds] <- list(tryCatch(
            run_V2_dataset(ds),
            error = function(e) {
                cat("  ERROR in dataset", ds, ":", conditionMessage(e), "\n")
                NULL
            }))
    }


    # =========================================================================
    # SECTION D: Save RDS and workbook
    # =========================================================================
    cat("\n", strrep("=", 72), "\n")
    cat("  SAVING RDS + EXCEL\n")
    cat(strrep("=", 72), "\n\n")

    # == The m = 2 BH family and the compact Figure 11 table ==================
    # BH across the two datasets for each omnibus and trend statistic.
    .f11 <- list()
    for (ds in DATASETS) {
        r <- all_results[[ds]]
        if (is.null(r) || is.null(r$val_tabs[["Fig11_Omnibus"]])) next
        o <- r$val_tabs[["Fig11_Omnibus"]]; o$dataset <- ds
        .f11[[ds]] <- o
    }
    if (length(.f11) > 0) {
        f11 <- do.call(rbind, .f11)
        key <- paste(f11$statistic, f11$unit, sep = "|")
        f11$omnibus_p_BH <- NA_real_; f11$trend_p_BH <- NA_real_
        for (k in unique(key)) {
            i <- which(key == k)
            f11$omnibus_p_BH[i] <- stats::p.adjust(f11$omnibus_p[i], "BH")
            f11$trend_p_BH[i]   <- stats::p.adjust(f11$trend_p[i], "BH")
        }
        f11$bh_family_n <- as.integer(table(key)[key])
        f11 <- f11[, c("dataset", "statistic", "space", "unit", "normal",
                       "shapiro_min_p", "omnibus_test", "omnibus_stat",
                       "omnibus_p", "omnibus_p_BH", "trend_test", "trend_stat",
                       "trend_p", "trend_p_BH", "n_unit", "bh_family_n",
                       "note")]
        arm_f11_omnibus <- f11

        # Compact table: one row per dataset; pick() returns NA for a missing entry.
        pick <- function(tb, cond, col) {
            if (is.null(tb) || !col %in% names(tb)) return(NA)
            i <- which(cond)
            if (length(i) == 0L) return(NA)
            tb[[col]][i[1L]]
        }
        rows <- list()
        for (ds in DATASETS) {
            r <- all_results[[ds]]
            if (is.null(r)) next
            de <- r$val_tabs[["Fig11_Descriptives"]]
            fo <- f11[f11$dataset == ds, , drop = FALSE]
            md <- function(stat, grp) {
                if (is.null(de)) return(NA_real_)
                i <- which(de$statistic == stat & de$unit == "cell" &
                               de$group == grp)
                if (length(i) == 0L) NA_real_ else de$median[i[1L]]
            }
            fp <- function(stat, unit, col) {
                i <- which(fo$statistic == stat & fo$unit == unit)
                if (length(i) == 0L) NA_real_ else fo[[col]][i[1L]]
            }
            pd <- r$val_tabs[["PERMANOVA_PERMDISP"]]
            nb <- r$val_tabs[["Neighbour_Proportions"]]
            kb <- r$val_tabs[["kBET"]]
            # Filter on unit only when the table has a unit column.
            .rp <- function(tb, col) {
                if (is.null(tb) || !all(c("focal", "reference") %in% names(tb)))
                    return(NA)
                cond <- !is.na(tb$focal) & tb$focal == "Recovered" &
                    !is.na(tb$reference) & tb$reference == "Post-IP TREP"
                if ("unit" %in% names(tb))
                    cond <- cond & !is.na(tb$unit) & tb$unit == "donor"
                pick(tb, cond, col)
            }
            rows[[ds]] <- data.frame(
                arm = ARM_TAG, dataset = ds,
                # Spine residual only if V2_RESID_IN_SUMMARY.
                resid_median_OrigRef   = if (isTRUE(V2_RESID_IN_SUMMARY))
                    md("Spine_residual", "Original + Refined") else NA_real_,
                resid_median_Recovered = if (isTRUE(V2_RESID_IN_SUMMARY))
                    md("Spine_residual", "Recovered") else NA_real_,
                resid_median_Excluded  = if (isTRUE(V2_RESID_IN_SUMMARY))
                    md("Spine_residual", "Excluded") else NA_real_,
                resid_trend_p_BH_donor = if (isTRUE(V2_RESID_IN_SUMMARY))
                    fp("Spine_residual", "donor", "trend_p_BH") else NA_real_,
                cop_proj_median_PreIPTREP  = md("PrePost_projection_copula", "Pre-IP TREP"),
                cop_proj_median_Recovered  = md("PrePost_projection_copula", "Recovered"),
                cop_proj_median_PostIPTREP = md("PrePost_projection_copula", "Post-IP TREP"),
                cop_proj_trend_p_donor     = fp("PrePost_projection_copula", "donor", "trend_p"),
                cop_proj_trend_p_BH_donor  = fp("PrePost_projection_copula", "donor", "trend_p_BH"),
                cop_proj_n_donor_blocks    = fp("PrePost_projection_copula", "donor", "n_unit"),
                cop_off_median_Recovered_lat = md("Lateral_displacement_copula", "Recovered"),
                cop_off_median_PreIPTREP_lat = md("Lateral_displacement_copula", "Pre-IP TREP"),
                cop_off_median_PostIPTREP_lat = md("Lateral_displacement_copula", "Post-IP TREP"),
                cop_lat_omnibus_p_BH_donor = fp("Lateral_displacement_copula", "donor", "omnibus_p_BH"),
                cop_bridge_median_PreIPTREP  = md("Bridge_distance_copula", "Pre-IP TREP"),
                cop_bridge_median_Recovered  = md("Bridge_distance_copula", "Recovered"),
                cop_bridge_median_PostIPTREP = md("Bridge_distance_copula", "Post-IP TREP"),
                cop_bridge_median_Excluded   = md("Bridge_distance_copula", "Excluded"),
                cop_bridge_omnibus_p_BH_donor = fp("Bridge_distance_copula", "donor", "omnibus_p_BH"),
                off_median_PreIPTREP       = md("PrePost_offaxis", "Pre-IP TREP"),
                off_median_Recovered       = md("PrePost_offaxis", "Recovered"),
                off_median_PostIPTREP      = md("PrePost_offaxis", "Post-IP TREP"),
                off_omnibus_p_BH_donor     = fp("PrePost_offaxis", "donor", "omnibus_p_BH"),
                proj_median_PreIPTREP  = md("PrePost_projection", "Pre-IP TREP"),
                proj_median_Recovered  = md("PrePost_projection", "Recovered"),
                proj_median_PostIPTREP = md("PrePost_projection", "Post-IP TREP"),
                proj_trend_p_donor     = fp("PrePost_projection", "donor", "trend_p"),
                proj_trend_p_BH_donor  = fp("PrePost_projection", "donor", "trend_p_BH"),
                proj_n_donor_blocks    = fp("PrePost_projection", "donor", "n_unit"),
                permanova_R2_RecVsPost   = suppressWarnings(as.numeric(.rp(pd, "R2"))),
                permanova_pBH_RecVsPost  = suppressWarnings(as.numeric(.rp(pd, "p_BH_FDR"))),
                nbr_enrich_RecVsPost     = suppressWarnings(as.numeric(.rp(nb, "mean_enrichment"))),
                kBET_reject_RecVsPost    = suppressWarnings(as.numeric(.rp(kb, "rejection_rate"))),
                stringsAsFactors = FALSE)
        }
        arm_f11_summary <- if (length(rows)) do.call(rbind, rows) else NULL
    } else {
        arm_f11_omnibus <- NULL; arm_f11_summary <- NULL
    }

    all_tabs <- list()
    for (ds in DATASETS) {
        r <- all_results[[ds]]
        if (is.null(r)) { cat("  [", ds, "] skipped (no results)\n"); next }

        rds_path <- file.path(OUT_DIR,
                              paste0(FILE_PREFIX, "_", ds, "_results.rds"))
        rds_save(r, rds_path, label = sprintf("[%s] RDS", ds))

        for (tnm in names(r$val_tabs)) {
            sheet <- paste0(ds, "_", tnm)
            # Sheet names: 31 characters at most, kept unique.
            if (nchar(sheet) > 31L) sheet <- substr(sheet, 1L, 31L)
            if (sheet %in% names(all_tabs)) {
                k <- 2L
                repeat {
                    cand <- paste0(substr(sheet, 1L, 31L - nchar(as.character(k)) - 1L),
                                   "_", k)
                    if (!cand %in% names(all_tabs)) break
                    k <- k + 1L
                }
                cat(sprintf("  sheet name '%s' collided; written as '%s'\n",
                            sheet, cand))
                sheet <- cand
            }
            all_tabs[[sheet]] <- r$val_tabs[[tnm]]
        }
    }

    # Arm-level sheets.
    if (!is.null(arm_f11_omnibus))
        all_tabs[["Fig11_Omnibus_BH"]] <- arm_f11_omnibus
    if (!is.null(arm_f11_summary))
        all_tabs[["Fig11_Summary"]] <- arm_f11_summary

    if (length(all_tabs) > 0) {
        xlsx_out <- file.path(OUT_DIR,
                              paste0(FILE_PREFIX, "_Validation2_Summary.xlsx"))
        tryCatch({
            write_xlsx(all_tabs, path = xlsx_out)
            cat(sprintf("  Excel saved: %s (%d sheet(s))\n",
                        basename(xlsx_out), length(all_tabs)))
        }, error = function(e)
            cat("  Excel error:", conditionMessage(e), "\n"))
    } else {
        cat("  WARNING: no tables to write\n")
    }


    # =========================================================================
    # SECTION E: Post-analysis guardrails and summary
    # =========================================================================
    cat("\n  POST-ANALYSIS GUARDRAILS\n\n")

    for (ds in DATASETS) {
        r <- all_results[[ds]]
        if (is.null(r)) {
            grv2(paste0("GR-V2-03-", ds), "FAIL", "Dataset returned NULL")
            next
        }

        grv2(paste0("GR-V2-03-", ds), "PASS",
             sprintf(paste0("embedded %d cell(s): OrigRef %d, Recovered %d, ",
                            "Orig-only %d, Excluded %d"),
                     r$n_cells,
                     unname(r$n_by_class["Original + Refined"]),
                     unname(r$n_by_class["Recovered"]),
                     unname(r$n_by_class["Original only"]),
                     unname(r$n_by_class["Excluded"])))

        # The count matrix must be integer-like.
        grv2(paste0("GR-V2-04-", ds),
             if (isTRUE(r$raw_integer_like) && r$raw_max > 20) "PASS" else "FAIL",
             sprintf("raw_selected max = %.0f, integer-like = %s | %d gene(s), %d PCA dim",
                     r$raw_max, if (isTRUE(r$raw_integer_like)) "TRUE" else "FALSE",
                     r$n_genes_kept, r$n_pca))

        # A 2D embedding is recorded.
        grv2(paste0("GR-V2-05-", ds),
             if (isTRUE(r$has_3d)) "PASS" else "WARN",
             sprintf("test axes: %s", paste(r$umap_axes, collapse = ", ")))

        # Donor composition of the Recovered class.
        ds_don  <- r$donor_summary
        rec_row <- ds_don[ds_don$cell_class == "Recovered", , drop = FALSE]
        top_share <- if (nrow(rec_row) == 1L) rec_row$top_donor_share[1] else NA_real_
        tvd       <- if (nrow(rec_row) == 1L) rec_row$TVD_vs_OrigRef[1] else NA_real_
        tvd_null  <- if (nrow(rec_row) == 1L) rec_row$TVD_null_p95[1]   else NA_real_
        # TVD against its own null.
        tvd_ok <- !is.na(tvd) && !is.na(tvd_null) && tvd <= tvd_null
        grv2(paste0("GR-V2-06-", ds),
             if (!is.na(top_share) && top_share < 0.5 && tvd_ok) "PASS" else "WARN",
             sprintf(paste0("Recovered donors: n = %s | largest donor share = %s",
                            " | TVD vs OrigRef = %s against a null p95 of %s%s"),
                     if (nrow(rec_row) == 1L) as.character(rec_row$n_donors[1]) else "NA",
                     if (is.na(top_share)) "NA" else sprintf("%.3f", top_share),
                     if (is.na(tvd)) "NA" else sprintf("%.3f", tvd),
                     if (is.na(tvd_null)) "NA" else sprintf("%.3f", tvd_null),
                     if (r$n_no_donor > 0)
                         sprintf(" | %d cell(s) without a donor", r$n_no_donor)
                     else ""))

        # Whether the pairwise tests had enough cells.
        grv2(paste0("GR-V2-07-", ds),
             if (r$n_recovered_tested >= V2_MIN_GROUP &&
                 r$n_origref_labelled >= V2_MIN_GROUP) "PASS" else "WARN",
             sprintf(paste0("tested: Recovered %d, Original-only %d, OrigRef with",
                            " a subgroup label %d (minimum %d)"),
                     r$n_recovered_tested, r$n_orig_only_tested,
                     r$n_origref_labelled, V2_MIN_GROUP))

        # Whether any cell-unit PERMANOVA pair was subsampled.
        pd <- r$val_tabs[["PERMANOVA_PERMDISP"]]
        pd_cell <- if (is.null(pd)) NULL else pd[pd$unit == "cell", , drop = FALSE]
        grv2(paste0("GR-V2-08-", ds),
             if (is.null(pd_cell) || nrow(pd_cell) == 0L) "WARN"
             else if (any(pd_cell$subsampled, na.rm = TRUE)) "WARN" else "PASS",
             if (is.null(pd_cell) || nrow(pd_cell) == 0L)
                 "PERMANOVA table absent (vegan missing or no pair qualified)"
             else sprintf("%d of %d cell-unit PERMANOVA pair(s) subsampled to %d per group",
                          sum(pd_cell$subsampled, na.rm = TRUE), nrow(pd_cell),
                          V2_PERM_MAX_N))

        # Whether the donor-unit tests ran.
        td     <- r$val_tabs[["TOST_Equivalence"]]
        n_td_d <- if (is.null(td)) 0L else sum(td$unit == "donor")
        n_pd_d <- if (is.null(pd)) 0L else sum(pd$unit == "donor")
        # Every Original+Refined cell must carry one of the four subgroups.
        .nu <- if (is.null(r$n_unknown_group)) NA_integer_
               else as.integer(r$n_unknown_group)[1]
        grv2(paste0("GR-V2-10-", ds),
             if (!is.finite(.nu)) "WARN" else if (.nu == 0L) "PASS" else "FAIL",
             if (!is.finite(.nu))
                 "the unknown-subgroup count was not recorded for this dataset"
             else if (.nu == 0L)
                 sprintf(paste0("every Original+Refined cell carries one of the",
                                " four IP x TREP subgroups from `%s`"),
                         V2_ORIGREF_GROUP_COL)
             else sprintf(paste0("%d Original+Refined cell(s) have a `%s`",
                                 " outside the four IP x TREP subgroups."),
                         as.integer(.nu), V2_ORIGREF_GROUP_COL))

        # == Figure 11 guardrails =============================================
        .fm <- r$val_tabs[["Fig11_Methods"]]
        .gv <- function(k) {
            if (is.null(.fm)) return(NA_character_)
            i <- which(.fm$item == k)
            if (length(i) == 0L) NA_character_ else .fm$value[i[1L]]
        }
        .axis <- .gv("projection axis")
        grv2(paste0("GR-V2-11-", ds),
             if (is.na(.axis) || identical(.axis, "not computed")) "WARN" else "PASS",
             sprintf("Pre->Post projection axis: %s",
                     if (is.na(.axis)) "not recorded" else .axis))

        .nb_pr <- suppressWarnings(as.integer(.gv("complete donor blocks, projection")))
        .nb_rs <- suppressWarnings(as.integer(.gv("complete donor blocks, residual")))
        # Complete donor blocks for the donor-unit tests (WARN below the minimum).
        grv2(paste0("GR-V2-12-", ds),
             if (!is.na(.nb_pr) && .nb_pr >= V2_STAT_MIN_BLOCKS) "PASS" else "WARN",
             sprintf(paste0("complete donor blocks: projection %s, spine",
                            " residual %s (minimum %d)"),
                     ifelse(is.na(.nb_pr), "n/a", as.character(.nb_pr)),
                     ifelse(is.na(.nb_rs), "n/a", as.character(.nb_rs)),
                     V2_STAT_MIN_BLOCKS))

        # SK copula axes: source and, when both routes ran, their agreement.
        .axis_cop <- .gv("copula projection axis")
        .cop_ok <- !is.na(.axis_cop) && !identical(.axis_cop, "not computed")
        grv2(paste0("GR-V2-13-", ds),
             if (.cop_ok) "PASS" else "WARN",
             sprintf(paste0("copula method check on (u_SK, v_SK): %s"),
                     if (.cop_ok) .axis_cop else "not computed"))

        grv2(paste0("GR-V2-09-", ds),
             if (!isTRUE(V2_DONOR_UNIT_TESTS)) "WARN"
             else if (n_td_d > 0L && n_pd_d > 0L) "PASS" else "WARN",
             sprintf(paste0("donor-unit (inverse-donor-weighted) tests: %d TOST,",
                            " %d PERMANOVA | thresholds >= %d donor(s) per group,",
                            " >= %d cell(s) per donor%s"),
                     n_td_d, n_pd_d, V2_MIN_DONORS, V2_MIN_CELLS_PER_DONOR,
                     if (!isTRUE(V2_DONOR_UNIT_TESTS))
                         " | V2_DONOR_UNIT_TESTS is FALSE" else ""))
    }

    grv2_write()
    cat(sprintf("\n  Final guardrail tally: %d PASS / %d WARN / %d FAIL\n",
                .grv2_pass, .grv2_warn, .grv2_fail))
    if (.grv2_fail > 0) {
        cat("  Failures present.\n")
    } else if (.grv2_warn > 0) {
        cat(sprintf("  %d passed, %d warning(s), no failures.\n",
                    .grv2_pass, .grv2_warn))
    } else {
        cat(sprintf("  All %d checks passed.\n", .grv2_pass))
    }

    cat("\n", strrep("=", 72), "\n")
    cat("  Part10e Validation2 COMPLETE [", ARM_TAG, "]\n")
    cat("  Output directory:", OUT_DIR, "\n")
    cat("  Files generated:\n")
    .transient <- c("_RUN_INCOMPLETE.txt", basename(MSG_PATH_V2))
    for (f in setdiff(sort(list.files(OUT_DIR)), .transient)) {
        sz <- round(file.size(file.path(OUT_DIR, f)) / 1e3, 0)
        cat(sprintf("    %-62s  %5d KB\n", f, sz))
    }
    cat(strrep("=", 72), "\n")

    cat("\n  Log written to:", LOG_PATH_V2, "\n")
    cat(sprintf("\nPart10e Validation2 - log ended: %s\n", Sys.time()))

    .p10v2_completed <- TRUE
    if (sink.number(type = "message") != 2L) sink(type = "message")
    while (sink.number() > 0) sink()
    suppressWarnings(file.remove(INCOMPLETE_MARKER_V2))

    invisible(list(arm = arm_tag, dir = OUT_DIR, pass = .grv2_pass,
                   warn = .grv2_warn, fail = .grv2_fail))
}


# =============================================================================
# SECTION 0c: Driver
# =============================================================================

cep_run_arms(PART_DIR_TAG, run_one_arm)
