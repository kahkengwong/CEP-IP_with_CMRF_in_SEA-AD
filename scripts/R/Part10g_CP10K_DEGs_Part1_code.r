# =============================================================================
# Part10g_CP10K_DEGs_Part1_code.r
# SEA-AD CEP-IP/CMRF, Part 10g, differential expression in the CMRF module
#   (Table S11)
#   Post-IP against Pre-IP module cells, each module half against the
#   non-module cells and the supertype background, and the non-module
#   contrasts. Every comparison is run at the cell unit (presto Wilcoxon) and
#   the donor unit (paired pseudobulk, one-sample Wilcoxon signed-rank on the
#   per-donor log2 fold change, BH across genes). Donor and library are nested,
#   so batch is absorbed by the donor unit.
#
# Usage (from the SEA-AD folder, or set SEAAD_ROOT), after
#   Part10e_CP10K_Validation1_code.r
#   source("Part10g_CP10K_DEGs_Part1_code.r")                         # all arms
#   CEP_ARM <- "PVM1_all"; source("Part10g_CP10K_DEGs_Part1_code.r")  # one arm
#
# Outputs (Part10g_DEGs/Part10g_DEGs_<arm>/): donor-unit workbook, optional
#   cell-unit workbook, figures (PDF, PNG), per-dataset results RDS, guardrail
#   ledger and logs.
# =============================================================================


# =============================================================================
# SECTION 0: Libraries and shared arm definitions
# =============================================================================

library(Matrix)
library(readxl)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(writexl)
library(grid)
library(tibble)

set.seed(42)

.v10g_presto_ok <- requireNamespace("presto", quietly = TRUE)
if (!.v10g_presto_ok)
    cat("NOTE: presto is not installed; the cell-unit Wilcoxon uses a base-R rank-sum.\n",
        sep = "")

if (!exists("CEP_HOME")) CEP_HOME <- Sys.getenv("SEAAD_ROOT", getwd())
.p10g_harness <- file.path(CEP_HOME, "CEP_arms.R")
if (!file.exists(.p10g_harness))
    stop("CEP_arms.R not found at: ", .p10g_harness, call. = FALSE)
source(.p10g_harness)

PART_TAG     <- "Part10g"
PART_DIR_TAG <- "Part10g_DEGs"

# == IP arm assignment ========================================================
# IP_BAND is kept at 0; no IP-zone exclusion exists, so a positive value stops.
IP_BAND <- 0.00
if (!is.finite(IP_BAND) || IP_BAND < 0)
    stop("IP_BAND must be a finite, non-negative number; got: ", IP_BAND,
         call. = FALSE)
if (IP_BAND > 0)
    stop("IP_BAND is set to ", IP_BAND, " but no IP-zone exclusion exists; set it to 0.",
         call. = FALSE)

# == DEG parameters ===========================================================
DEG_MIN_PCT      <- 0.05
DEG_LOGFC_THRESH <- 0
# Significance |log2FC| thresholds, calibrated for the CP10K scale.
DEG_SIG_LOGFC    <- 0.25   # head-to-head
DEG_SIG_LOGFC_BG <- 0.125  # arm vs other cells
DEG_SIG_PADJ     <- 0.05
DEG_MIN_CELLS_FG <- 10L
DEG_MIN_CELLS_BG <- 10L
# Below this size ratio the up/down guardrail keeps its down-only rule.
GR11_SIZE_RATIO  <- 1.5
# Below this many significant genes the up/down ratio is not assessed.
GR11_MIN_SIG     <- 10L
VOLCANO_TOP_N    <- 20L
DOT_TOP_N        <- 30L

# == Donor unit ===============================================================
V10G_DONOR_UNIT          <- TRUE
V10G_MIN_CELLS_PER_DONOR <- 3L    # per donor, in EACH arm, to be paired
V10G_MIN_DONORS          <- 5L    # paired donors needed to test at all
V10G_MIN_DONOR_PCT       <- 0.05  # gene detected in this share of paired donors
V10G_SEED                <- 42L
# The cell unit is always computed; this controls only its workbook and figures.
V10G_CELL_OUTPUTS        <- FALSE

# Keep the per-donor fold-change matrix (needed for per-gene intervals).
V10G_STORE_DONOR_FC <- TRUE

# Contrasts it is stored for.
V10G_FC_CONTRASTS   <- c("PostIP_vs_NMPost", "PreIP_vs_NMPre",
                         "PostIP_vs_PreIP",  "PostIP_vs_NMAll")

# Genes per block in the base-R rank-sum fallback (bounds memory).
V10G_RANKSUM_BLOCK       <- 2000L

# == Arm-vs-background ========================================================
# Background cells come from Part 8's dataset-wide raw matrix, normalised as
# ln(CP10K + 1) and checked against deg_support/norm.mtx on the selected cells.
V10G_BG_ENABLE    <- TRUE
V10G_BG_MAX_CELLS <- 150000L  # refuse to read a dataset matrix beyond this
V10G_NORM_TOL     <- 1e-6     # max |recomputed - exported| on the selected cells

# == CRS parameters ===========================================================
PCT_FILTER <- 0.10
TOP_N_CRS  <- 1000L

# == Colours ==================================================================
COL_PREIP_FILL    <- "#1B3A6B"
COL_PREIP_BORDER  <- "#0D1F3C"
COL_POSTIP_FILL   <- "#2D2D2D"
COL_POSTIP_BORDER <- "#000000"
COL_IPZONE_FILL   <- "#888888"
COL_NS   <- "#E0E0E0"
COL_UP   <- colorRampPalette(c("#FFFFAA", "#FFFF66", "#FFD700", "#FFA500",
                               "#FF6633", "#E63946", "#CC0000", "#990000",
                               "#660000"))(100)
COL_DOWN <- colorRampPalette(c("#CCFFE6", "#66FFB2", "#33CCCC", "#20B2AA",
                               "#1E90FF", "#2E86C1", "#4682B4", "#1A3D7C",
                               "#05335C"))(100)

MT_GENES <- c("MT-CO3", "MT-ATP6", "MT-ND3", "MT-CO2", "MT-ND4",
              "MT-ND2", "MT-ND1", "MT-CYB", "MT-CO1", "MT-ND5",
              "MT-ND4L", "MT-ATP8", "MTRNR2L12")

MODULE_CLASSES <- c("Original + Refined", "Recovered")

# == Non-module cells =========================================================
# Classes of Step4_Refined_Cells outside MODULE_CLASSES, derived per dataset;
# the expected names are used only to report a mismatch.
NONMODULE_CLASSES_EXPECTED <- c("Original only", "Excluded")

V10G_NONMOD_ENABLE <- TRUE

# Comparison tags; each must be unique in its first 17 characters (the
# 31-character sheet-name limit).
V10G_NONMOD_CMP <- c("NMPost_vs_NMPre", "PostIP_vs_NMPost", "PreIP_vs_NMPre")
# Pooled non-module group (both sides of the IP); both module halves use this pool.
V10G_NONMOD_CMP2 <- c("PostIP_vs_NMAll", "PreIP_vs_NMAll")
# Pooled non-module group against the shared background (BG significance gate).
V10G_NONMOD_BG  <- "NMAll_vs_BG"
V10G_CMP_ALL    <- c("PostIP_vs_PreIP", "PostIP_vs_BG", "PreIP_vs_BG",
                     V10G_NONMOD_CMP, V10G_NONMOD_CMP2, V10G_NONMOD_BG)
if (any(nchar(V10G_CMP_ALL) > 17L) ||
    anyDuplicated(substr(V10G_CMP_ALL, 1L, 17L)) > 0L)
    stop("Comparison tags must be at most 17 characters and unique in their ",
         "first 17, or mk_sheet() will collapse two sheets into one.\n  Got: ",
         paste(V10G_CMP_ALL, collapse = ", "), call. = FALSE)


# =============================================================================
# SECTION 0d: Pure helpers
# =============================================================================

as_flag <- function(x) {
    if (is.logical(x)) return(x)
    if (is.numeric(x)) return(x != 0)
    y <- toupper(trimws(as.character(x)))
    out <- rep(NA, length(y))
    out[y %in% c("TRUE", "T", "1", "YES")] <- TRUE
    out[y %in% c("FALSE", "F", "0", "NO")] <- FALSE
    out
}

# A subset condition must have no NA (an NA gives an all-NA phantom row).
require_no_na <- function(x, nm) {
    if (anyNA(x))
        stop(sum(is.na(x)), " missing value(s) in `", nm, "`.", call. = FALSE)
    x
}

show_plot <- function(p) {
    tryCatch({
        if (inherits(p, "ggplot")) print(p) else { grid.newpage(); grid.draw(p) }
    }, error = function(e) invisible(NULL))
}

# =============================================================================
# Seurat-convention fold change on a log1p CP10K matrix:
#   log2(mean(expm1(x_fg)) + 1) - log2(mean(expm1(x_bg)) + 1)
# =============================================================================
v10g_log2fc <- function(mat_fg, mat_bg) {
    f <- mat_fg; f@x <- expm1(f@x)
    b <- mat_bg; b@x <- expm1(b@x)
    log2(Matrix::rowMeans(f) + 1) - log2(Matrix::rowMeans(b) + 1)
}

# =============================================================================
# CP10K recipe as R_config.R documents it. X is cells x genes raw counts;
# the divisor is each cell's total over all genes.
# =============================================================================
p10g_cp10k <- function(X, target = NORM_TARGET_SUM) {
    X <- as(X, "CsparseMatrix")
    tot <- Matrix::rowSums(X)
    tot[tot == 0] <- 1
    # Per-row divisor applied through @i (row index of each stored value).
    Y <- X
    Y@x <- log1p(Y@x / tot[Y@i + 1L] * target)
    Y
}

# =============================================================================
# Check the CP10K recipe against deg_support/norm.mtx on the selected cells;
# the background comparisons are skipped if they do not agree.
# =============================================================================
p10g_verify_cp10k <- function(raw_dir, expr_norm, tol = V10G_NORM_TOL,
                              n_cells = 400L, n_genes = 2000L) {
    mp <- file.path(raw_dir, "matrix.mtx")
    bp <- file.path(raw_dir, "matrix_barcodes.tsv")
    fp <- file.path(raw_dir, "matrix_features.tsv")
    if (!all(file.exists(c(mp, bp, fp))))
        return(list(ok = FALSE, dev = NA_real_,
                    note = paste0("raw_selected not found under ", raw_dir)))
    X  <- as(readMM(mp), "CsparseMatrix")
    bc <- readLines(bp, warn = FALSE)
    ft <- read.table(fp, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
    if (nrow(X) != length(bc) || ncol(X) != nrow(ft))
        return(list(ok = FALSE, dev = NA_real_,
                    note = "raw_selected dimensions disagree with its barcodes or features"))
    rownames(X) <- bc
    colnames(X) <- make.unique(as.character(ft[, 2]))
    N <- Matrix::t(p10g_cp10k(X))          # genes x cells, CP10K
    rm(X); gc(verbose = FALSE)

    cc <- intersect(colnames(N), colnames(expr_norm))
    gg <- intersect(rownames(N), rownames(expr_norm))
    if (length(cc) < 10L || length(gg) < 10L)
        return(list(ok = FALSE, dev = NA_real_,
                    note = "raw_selected and norm.mtx share too few cells or genes"))
    set.seed(V10G_SEED)
    cc <- if (length(cc) > n_cells) sample(cc, n_cells) else cc
    gg <- if (length(gg) > n_genes) sample(gg, n_genes) else gg
    dev <- max(abs(as.matrix(N[gg, cc, drop = FALSE]) -
                   as.matrix(expr_norm[gg, cc, drop = FALSE])))
    list(ok = is.finite(dev) && dev <= tol, dev = dev,
         n_cells = length(cc), n_genes = length(gg), note = "")
}

# =============================================================================
# Non-M4-top PRIMARY_STYPE background on the CP10K scale: PRIMARY_STYPE cells
# of obs_dataset.csv outside the selected set, from Part 8's raw matrix.
# =============================================================================
p10g_load_background <- function(ds_dir_path, shared_ds_dir, sel_bc,
                                 stype_col, primary_stype, donor_col,
                                 max_cells = V10G_BG_MAX_CELLS) {
    obs_p <- file.path(ds_dir_path, "metadata", "obs_dataset.csv")
    if (!file.exists(obs_p))
        return(list(mat = NULL, note = paste0("obs_dataset.csv not found: ", obs_p)))
    obs <- read.csv(obs_p, stringsAsFactors = FALSE, check.names = FALSE)
    if (!("barcode" %in% names(obs))) obs$barcode <- obs[[1]]
    if (!(stype_col %in% names(obs)))
        return(list(mat = NULL,
                    note = sprintf("obs_dataset.csv has no `%s` column", stype_col)))
    # Background donors come from obs_dataset.csv.
    dcol <- c(donor_col, "donor", "Donor ID")
    dcol <- dcol[dcol %in% names(obs)][1]
    bg_donor <- if (is.na(dcol)) NULL
                else setNames(as.character(obs[[dcol]]), obs$barcode)
    st_all <- obs$barcode[!is.na(obs[[stype_col]]) &
                              obs[[stype_col]] == primary_stype]
    bg_bc  <- setdiff(st_all, sel_bc)
    # Report how many selected cells are PRIMARY_STYPE.
    n_sel_stype <- sum(sel_bc %in% st_all)
    cat(sprintf("    %s in the dataset: %d | selected: %d (%d of them %s) | background: %d\n",
                primary_stype, length(st_all), length(sel_bc), n_sel_stype,
                primary_stype, length(bg_bc)))
    if (n_sel_stype < length(sel_bc))
        cat(sprintf("    NOTE: %d selected cell(s) are not %s in obs_dataset.csv.\n",
                    length(sel_bc) - n_sel_stype, primary_stype))
    if (length(bg_bc) < DEG_MIN_CELLS_BG)
        return(list(mat = NULL, n_bg = length(bg_bc),
                    note = sprintf("only %d non-selected %s cell(s)",
                                   length(bg_bc), primary_stype)))

    rd <- file.path(shared_ds_dir, "counts", "raw_dataset")
    mp <- file.path(rd, "matrix.mtx")
    bp <- file.path(rd, "matrix_barcodes.tsv")
    fp <- file.path(rd, "matrix_features.tsv")
    if (!all(file.exists(c(mp, bp, fp))))
        return(list(mat = NULL, n_bg = length(bg_bc),
                    note = paste0("raw_dataset not found under ", rd)))
    bc <- readLines(bp, warn = FALSE)
    if (length(bc) > max_cells)
        return(list(mat = NULL, n_bg = length(bg_bc),
                    note = sprintf(paste0("raw_dataset holds %d cell(s), above the",
                                          " V10G_BG_MAX_CELLS cap of %d"), length(bc), max_cells)))
    have <- intersect(bg_bc, bc)
    if (length(have) < DEG_MIN_CELLS_BG)
        return(list(mat = NULL, n_bg = length(bg_bc),
                    note = sprintf("only %d of %d background cell(s) are in raw_dataset",
                                   length(have), length(bg_bc))))
    ft <- read.table(fp, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
    cat(sprintf("    raw_dataset: %d cell(s) x %d gene(s); reading\n",
                length(bc), nrow(ft)))
    X <- as(readMM(mp), "CsparseMatrix")
    if (nrow(X) != length(bc) || ncol(X) != nrow(ft))
        return(list(mat = NULL, n_bg = length(bg_bc),
                    note = "raw_dataset dimensions disagree with its barcodes or features"))
    rownames(X) <- bc
    colnames(X) <- make.unique(as.character(ft[, 2]))
    # Subset cells (not genes) before normalising; a gene subset would change the divisor.
    X <- X[have, , drop = FALSE]
    N <- Matrix::t(p10g_cp10k(X))
    rm(X); gc(verbose = FALSE)
    cat(sprintf("    background matrix: %d gene(s) x %d cell(s)\n",
                nrow(N), ncol(N)))
    dv <- if (is.null(bg_donor)) NULL else unname(bg_donor[colnames(N)])
    if (is.null(dv))
        cat(sprintf(paste0("    WARNING: obs_dataset.csv has no donor column",
                           " (looked for %s); BG comparisons are cell-unit only\n"),
                    paste(unique(c(donor_col, "donor", "Donor ID")),
                          collapse = ", ")))
    else
        cat(sprintf("    background donors: %d distinct | %d cell(s) without one\n",
                    length(unique(dv[!is.na(dv)])), sum(is.na(dv))))
    list(mat = N, n_bg = length(bg_bc), donor = dv, donor_col = dcol,
         n_stype = length(st_all), n_sel_stype = n_sel_stype, note = "")
}

# =============================================================================
# Cell-unit rank-sum: presto when available, otherwise a base-R rank-sum with
# tie correction (same p-values and AUC).
# =============================================================================
v10g_ranksum <- function(mat, is_fg) {
    if (.v10g_presto_ok) {
        grp <- ifelse(is_fg, "fg", "bg")
        wa  <- presto::wilcoxauc(mat, grp)
        wa  <- wa[wa$group == "fg", , drop = FALSE]
        return(data.frame(gene = as.character(wa$feature),
                          pval = as.numeric(wa$pval),
                          auc  = as.numeric(wa$auc),
                          stringsAsFactors = FALSE))
    }
    n1 <- sum(is_fg); n2 <- sum(!is_fg); n <- n1 + n2
    ng <- nrow(mat)
    mu  <- n1 * n2 / 2

    # One gene block at a time; ranks are within a gene, so values are unchanged.
    blk <- max(1L, as.integer(V10G_RANKSUM_BLOCK))
    auc <- numeric(ng)
    p   <- numeric(ng)
    for (s in seq.int(1L, ng, by = blk)) {
        e   <- min(s + blk - 1L, ng)
        idx <- s:e
        dm  <- as.matrix(mat[idx, , drop = FALSE])
        rk  <- t(apply(dm, 1L, rank))                 # ranks within each gene
        if (is.null(dim(rk))) rk <- matrix(rk, nrow = length(idx))
        R1  <- rowSums(rk[, is_fg, drop = FALSE])
        U1  <- R1 - n1 * (n1 + 1) / 2
        tie <- apply(dm, 1L, function(v) {
            tt <- table(v); sum(tt^3 - tt)
        })
        sdU <- sqrt((n1 * n2 / 12) * ((n + 1) - tie / (n * (n - 1))))
        auc[idx] <- U1 / (n1 * n2)
        p[idx]   <- 2 * pnorm(-abs((U1 - mu) / sdU))
        rm(dm, rk); gc(FALSE)
    }
    p[!is.finite(p)] <- 1
    data.frame(gene = rownames(mat), pval = p, auc = auc,
               stringsAsFactors = FALSE)
}

# == Smallest p the donor unit can return =====================================
# Normal approximation with continuity correction, as in v10g_pseudobulk.
v10g_sr_p_floor <- function(nd) {
    nd <- suppressWarnings(as.integer(nd))
    if (!isTRUE(is.finite(nd)) || nd < 3L) return(NA_real_)
    z <- (nd * (nd + 1) / 2 - nd * (nd + 1) / 4 - 0.5) /
         sqrt(nd * (nd + 1) * (2 * nd + 1) / 24)
    2 * stats::pnorm(-abs(z))
}

# =============================================================================
# Donor unit: paired pseudobulk. Donors with at least min_cells cells in both
# arms; per donor and gene the Seurat-convention log2 fold change; one-sample
# Wilcoxon signed-rank against zero, BH over genes. Effect size: matched-pairs
# rank-biserial correlation (W+ - W-) / (W+ + W-).
# =============================================================================
v10g_pseudobulk <- function(mat, is_fg, donor,
                            min_cells  = V10G_MIN_CELLS_PER_DONOR,
                            min_donors = V10G_MIN_DONORS,
                            min_pct    = V10G_MIN_DONOR_PCT) {

    d <- as.character(donor)
    ok <- !is.na(d) & nzchar(d)
    n_fg_by <- table(d[ok &  is_fg])
    n_bg_by <- table(d[ok & !is_fg])
    all_d   <- union(names(n_fg_by), names(n_bg_by))
    cnt_fg  <- setNames(as.integer(n_fg_by[all_d]), all_d); cnt_fg[is.na(cnt_fg)] <- 0L
    cnt_bg  <- setNames(as.integer(n_bg_by[all_d]), all_d); cnt_bg[is.na(cnt_bg)] <- 0L
    paired  <- all_d[cnt_fg >= min_cells & cnt_bg >= min_cells]

    info <- list(n_donors_total = length(all_d),
                 n_donors_paired = length(paired),
                 n_donors_fg_only = sum(cnt_fg >= min_cells & cnt_bg < min_cells),
                 n_donors_bg_only = sum(cnt_bg >= min_cells & cnt_fg < min_cells),
                 min_cells = min_cells)

    if (length(paired) < min_donors) {
        info$note <- sprintf(paste0("only %d donor(s) have >= %d cell(s) in",
                                    " both arms; %d needed"),
                             length(paired), min_cells, min_donors)
        return(list(res = NULL, info = info))
    }

    # per-donor fold change, one column per donor
    fc <- vapply(paired, function(dd) {
        i_f <- which(ok & is_fg  & d == dd)
        i_b <- which(ok & !is_fg & d == dd)
        v10g_log2fc(mat[, i_f, drop = FALSE], mat[, i_b, drop = FALSE])
    }, numeric(nrow(mat)))
    if (is.null(dim(fc))) fc <- matrix(fc, nrow = nrow(mat))
    rownames(fc) <- rownames(mat)

    # a gene must move in enough donors to be testable at all
    det <- rowMeans(fc != 0)
    keep <- det >= min_pct
    info$n_genes_in  <- nrow(fc)
    info$n_genes_out <- sum(keep)
    if (!any(keep)) {
        info$note <- "no gene reached the donor detection floor"
        return(list(res = NULL, info = info))
    }
    fc <- fc[keep, , drop = FALSE]

    nd <- ncol(fc)
    stat <- t(apply(fc, 1L, function(v) {
        v <- v[is.finite(v)]
        nz <- v[v != 0]
        if (length(nz) < 3L) return(c(NA_real_, NA_real_, NA_real_))
        wt <- suppressWarnings(tryCatch(
            stats::wilcox.test(nz, mu = 0, alternative = "two.sided",
                               exact = FALSE, correct = TRUE),
            error = function(e) NULL))
        p <- if (is.null(wt)) NA_real_ else as.numeric(wt$p.value)
        r  <- rank(abs(nz))
        Wp <- sum(r[nz > 0]); Wm <- sum(r[nz < 0])
        rrb <- if ((Wp + Wm) > 0) (Wp - Wm) / (Wp + Wm) else NA_real_
        c(p, rrb, mean(v))
    }))

    list(res = data.frame(gene = rownames(fc),
                          pval = stat[, 1], r_rb = stat[, 2],
                          avg_log2FC = stat[, 3],
                          n_donors = nd,
                          stringsAsFactors = FALSE),
         info = info, fc = fc)
}


# =============================================================================
# SECTION 0e: Plotting
# =============================================================================

theme_deg <- function(base_size = 10) {
    theme_minimal(base_size = base_size) +
        theme(
            plot.title       = element_text(face = "bold", size = base_size + 1),
            plot.subtitle    = element_text(size = base_size - 1, colour = "grey35",
                                            lineheight = 1.2),
            axis.title       = element_text(face = "bold", size = base_size - 1),
            axis.text        = element_text(size = base_size - 2),
            panel.grid.major = element_line(colour = "grey92"),
            panel.grid.minor = element_blank(),
            panel.border     = element_rect(colour = "black", fill = NA,
                                            linewidth = 0.3),
            plot.background  = element_rect(fill = "white", colour = NA),
            panel.background = element_rect(fill = "white", colour = NA),
            legend.position  = "right"
        )
}

# Save to PDF and PNG (no console print)
save_pdf_png <- function(p, base, w = 14, h = 10, dpi = 200) {
    is_grid <- inherits(p, "grob") || inherits(p, "gtable")
    draw_fn <- if (is_grid) function() { grid.newpage(); grid.draw(p) }
               else function() print(p)
    pdf(paste0(base, ".pdf"), width = w, height = h)
    draw_fn(); dev.off()
    png(paste0(base, ".png"), width = w * dpi, height = h * dpi, res = dpi)
    draw_fn(); dev.off()
}

# == Volcano theme ============================================================

theme_vol <- function(aspect = NULL) {
    theme_bw(base_size = 14) +
        theme(
            plot.title       = element_text(face = "bold", size = 13,
                                            hjust = 0.5),
            plot.subtitle    = element_text(size = 9, colour = "grey35",
                                            hjust = 0.5, lineheight = 1.2),
            axis.title       = element_text(face = "bold", size = 11),
            axis.text        = element_text(colour = "#303030", size = 9),
            panel.grid.major = element_line(colour = "grey88", linewidth = 0.3),
            panel.grid.minor = element_blank(),
            panel.border     = element_rect(colour = "#303030", fill = NA,
                                            linewidth = 0.5),
            plot.background  = element_rect(fill = "white", colour = NA),
            panel.background = element_rect(fill = "white", colour = NA),
            legend.position  = "none",
            aspect.ratio     = aspect
        )
}

# == CRS-gradient colours for volcano points ==================================
# Significant genes are binned by CRS into the 100-stop gradient; NS genes grey.
assign_volcano_colours <- function(df, fc_thresh = DEG_SIG_LOGFC) {
    df$point_colour <- COL_NS

    idx_up <- which(df$significant & df$avg_log2FC > fc_thresh)
    if (length(idx_up) > 0) {
        ss <- df$CRS[idx_up]
        lo <- min(ss); hi <- max(ss)
        if (isTRUE(lo == hi)) {
            df$point_colour[idx_up] <- COL_UP[50]
        } else {
            bins <- as.integer(cut(ss, breaks = seq(lo, hi, length.out = 101),
                                   include.lowest = TRUE))
            bins[is.na(bins)] <- 1L
            df$point_colour[idx_up] <- COL_UP[bins]
        }
    }

    idx_dn <- which(df$significant & df$avg_log2FC < -fc_thresh)
    if (length(idx_dn) > 0) {
        ss <- df$CRS[idx_dn]
        lo <- min(ss); hi <- max(ss)
        if (isTRUE(lo == hi)) {
            df$point_colour[idx_dn] <- COL_DOWN[50]
        } else {
            bins <- as.integer(cut(ss, breaks = seq(lo, hi, length.out = 101),
                                   include.lowest = TRUE))
            bins[is.na(bins)] <- 1L
            df$point_colour[idx_dn] <- COL_DOWN[bins]
        }
    }
    df
}

# == Plot 1: module cell counts by arm ========================================
plot_arm_counts <- function(rc, ds) {
    # Module cells (OriRef + Recovered) per arm.
    module <- rc[rc$cell_class %in% MODULE_CLASSES & !is.na(rc$arm), ]
    arms <- if (IP_BAND > 0) c("Pre-IP", "IP zone", "Post-IP")
            else c("Pre-IP", "Post-IP")
    cnt_df <- data.frame(
        arm     = arms,
        n_cells = vapply(arms, function(a)
            sum(module$arm == a, na.rm = TRUE), integer(1L)),
        stringsAsFactors = FALSE
    )
    cnt_df$arm_f <- factor(cnt_df$arm, levels = arms)

    ggplot(cnt_df, aes(x = arm_f, y = n_cells, fill = arm_f)) +
        geom_col(width = 0.65, colour = "grey20", linewidth = 0.3) +
        geom_text(aes(label = n_cells), vjust = -0.5, size = 3.5,
                  fontface = "bold") +
        scale_fill_manual(
            values = c("Pre-IP"  = COL_PREIP_FILL,
                       "IP zone" = COL_IPZONE_FILL,
                       "Post-IP" = COL_POSTIP_FILL),
            guide = "none") +
        scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
        labs(title    = sprintf("Module cell counts by IP arm - %s", ds),
             subtitle = sprintf(
                 "OriRef + Recovered | arms inherited: Pre-IP when cDFG < IP_y, never split on GOI"),
             x = "Arm", y = "Number of cells") +
        theme_deg()
}

# == Plot 2: DEG count bar ====================================================
plot_deg_bar <- function(deg_df, ds) {
    if (is.null(deg_df) || nrow(deg_df) == 0) return(NULL)

    sig <- deg_df[deg_df$significant, ]
    cnt_df <- data.frame(
        direction = c("Up in Post-IP", "Up in Pre-IP"),
        n_degs    = c(sum(sig$avg_log2FC > 0, na.rm = TRUE),
                       sum(sig$avg_log2FC < 0, na.rm = TRUE)),
        stringsAsFactors = FALSE
    )
    cnt_df$direction_f <- factor(cnt_df$direction,
                                  levels = c("Up in Post-IP", "Up in Pre-IP"))

    ggplot(cnt_df, aes(x = direction_f, y = n_degs, fill = direction_f)) +
        geom_col(width = 0.55, colour = "grey20", linewidth = 0.3) +
        geom_text(aes(label = n_degs), vjust = -0.5, size = 4,
                  fontface = "bold") +
        scale_fill_manual(
            values = c("Up in Post-IP" = COL_POSTIP_FILL,
                       "Up in Pre-IP"  = COL_PREIP_FILL),
            guide = "none") +
        scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
        labs(title    = sprintf("Significant DEGs - %s", ds),
             subtitle = sprintf(
                 "Post-IP vs Pre-IP | BH-FDR < %.2f & |log2FC| >= %.1f",
                 DEG_SIG_PADJ, DEG_SIG_LOGFC),
             x = NULL, y = "Number of DEGs") +
        theme_deg()
}

# == Plot 3: volcano ==========================================================
# CRS-gradient fill; top-5 up and down genes in the subtitle.
plot_volcano <- function(deg_df, ds, fc_thresh = DEG_SIG_LOGFC) {
    if (is.null(deg_df) || nrow(deg_df) == 0) return(NULL)

    df <- deg_df
    df$neg_log10_fdr <- -log10(pmax(df$p_val_adj_BH, 1e-300))

    # CRS-gradient colour assignment
    df <- assign_volcano_colours(df, fc_thresh = fc_thresh)

    n_up <- sum(df$significant & df$avg_log2FC >  fc_thresh, na.rm = TRUE)
    n_dn <- sum(df$significant & df$avg_log2FC < -fc_thresh, na.rm = TRUE)

    # Genes with no finite BH or log2FC are counted in the subtitle.
    n_unplot <- sum(!is.finite(df$neg_log10_fdr) | !is.finite(df$avg_log2FC))

    x_lo <- min(df$avg_log2FC,      na.rm = TRUE)
    x_hi <- max(df$avg_log2FC,      na.rm = TRUE)
    if (isTRUE(x_lo == x_hi)) { x_lo <- x_lo - 0.1; x_hi <- x_hi + 0.1 }

    y_hi  <- max(df$neg_log10_fdr, na.rm = TRUE)
    y_top <- y_hi * 1.06

    # Keep the significance line inside the panel.
    y_lim <- max(y_top * 1.04, -log10(DEG_SIG_PADJ) * 1.15)
    y_ann <- max(y_top, y_lim * 0.96)

    # Top-5 gene names by volcano_score for subtitle
    df$volcano_score <- df$neg_log10_fdr * abs(df$avg_log2FC)
    sig_up <- df[df$significant & df$avg_log2FC > fc_thresh, ]
    sig_dn <- df[df$significant & df$avg_log2FC < -fc_thresh, ]
    top5_up <- head(sig_up[order(-sig_up$volcano_score), ], 5)$gene
    top5_dn <- head(sig_dn[order(-sig_dn$volcano_score), ], 5)$gene
    up_str <- if (length(top5_up) > 0) paste(top5_up, collapse = ", ") else "none"
    dn_str <- if (length(top5_dn) > 0) paste(top5_dn, collapse = ", ") else "none"

    sub_text <- sprintf(
        "X = log2FC  |  Y = -log10(BH-FDR)  |  colour intensity = CRS  |  |log2FC| >= %.2f  &  BH-FDR < %.2f%s\nUp: %s\nDown: %s",
        fc_thresh, DEG_SIG_PADJ,
        if (n_unplot > 0)
            sprintf("  |  %d gene(s) not drawn (no finite BH or log2FC)",
                    n_unplot) else "",
        up_str, dn_str)

    p <- ggplot(df, aes(x = avg_log2FC, y = neg_log10_fdr,
                         fill = point_colour)) +
        geom_point(alpha = 0.9, size = 3.5, shape = 21,
                   stroke = 0.2, colour = "#303030") +
        scale_fill_identity() +
        geom_vline(xintercept =  fc_thresh, linetype = "dashed",
                   colour = "#303030", linewidth = 0.5) +
        geom_vline(xintercept = -fc_thresh, linetype = "dashed",
                   colour = "#303030", linewidth = 0.5) +
        geom_hline(yintercept = -log10(DEG_SIG_PADJ), linetype = "dashed",
                   colour = "#303030", linewidth = 0.5) +
        xlim(x_lo, x_hi) +
        ylim(0, y_lim) +
        # n= corner annotations
        annotate("text", x = x_hi, y = y_ann,
                 label = paste0("n=", n_up), size = 5, hjust = 1,
                 fontface = "bold") +
        annotate("text", x = x_lo, y = y_ann,
                 label = paste0("n=", n_dn), size = 5, hjust = 0,
                 fontface = "bold") +
        labs(
            title    = ds,
            subtitle = sub_text,
            x = expression(bold(paste("log"[2], " FC"))),
            y = expression(bold(-log[10](FDR)))
        ) +
        theme_vol(aspect = 1)
    p
}

# == Plot 4: top-genes dot plot (pct x mean expr) =============================
plot_dotplot <- function(deg_df, expr_mat, fg_barcodes, bg_barcodes,
                          ds, top_n = DOT_TOP_N) {
    if (is.null(deg_df) || nrow(deg_df) == 0) return(NULL)

    sig <- deg_df[deg_df$significant, ]
    if (nrow(sig) == 0) sig <- head(deg_df[order(deg_df$p_val_adj_BH), ], top_n)
    if (!"neg_log10p" %in% names(sig))
        sig$neg_log10p <- -log10(pmax(sig$p_val_adj_BH, 1e-300))
    sig$vs <- sig$neg_log10p * abs(sig$avg_log2FC)

    n_half   <- ceiling(top_n / 2)
    top_up   <- head(sig[sig$avg_log2FC > 0, ][
        order(-sig$vs[sig$avg_log2FC > 0]), ], n_half)
    top_down <- head(sig[sig$avg_log2FC < 0, ][
        order(-sig$vs[sig$avg_log2FC < 0]), ], n_half)
    top_genes <- c(top_up$gene, top_down$gene)
    if (length(top_genes) == 0) return(NULL)

    fg_bc <- intersect(fg_barcodes, colnames(expr_mat))
    bg_bc <- intersect(bg_barcodes, colnames(expr_mat))


    dot_rows <- lapply(top_genes, function(g) {
        if (!(g %in% rownames(expr_mat))) return(NULL)
        fg_expr <- as.numeric(expr_mat[g, fg_bc])
        bg_expr <- as.numeric(expr_mat[g, bg_bc])
        data.frame(
            gene      = g,
            group     = c("Post-IP", "Pre-IP"),
            pct_expr  = c(mean(fg_expr > 0) * 100,
                           mean(bg_expr > 0) * 100),
            mean_expr = c(mean(fg_expr), mean(bg_expr)),
            stringsAsFactors = FALSE)
    })
    dot_df <- do.call(rbind, Filter(Negate(is.null), dot_rows))
    if (is.null(dot_df) || nrow(dot_df) == 0) return(NULL)

    dot_df$gene_f  <- factor(dot_df$gene, levels = rev(unique(dot_df$gene)))
    dot_df$group_f <- factor(dot_df$group, levels = c("Post-IP", "Pre-IP"))

    ggplot(dot_df, aes(x = group_f, y = gene_f,
                        size = pct_expr, colour = mean_expr)) +
        geom_point() +
        scale_size_continuous(range = c(0.5, 7), name = "% expressed") +
        scale_colour_gradient2(
            low = COL_PREIP_FILL, mid = "grey90", high = COL_POSTIP_FILL,
            midpoint = median(dot_df$mean_expr, na.rm = TRUE),
            name = "Mean expr\n(log-norm)") +
        labs(title    = sprintf("Top-DEG dot plot - %s", ds),
             subtitle = sprintf(
                 "Post-IP vs Pre-IP | top %d genes by -log10(BH-p) x |log2FC|",
                 top_n),
             x = NULL, y = "Gene") +
        theme_deg() +
        theme(axis.text.y = element_text(size = 7))
}

# == Plot 5: CRS distribution among eligible genes ============================
plot_crs_distribution <- function(deg_df, ds) {
    if (is.null(deg_df) || nrow(deg_df) == 0) return(NULL)
    elig <- deg_df[deg_df$eligible, ]
    if (nrow(elig) < 5) return(NULL)

    elig$dir <- ifelse(elig$avg_log2FC > 0, "Up in Post-IP", "Up in Pre-IP")

    ggplot(elig, aes(x = CRS, fill = dir)) +
        geom_histogram(bins = 40, alpha = 0.7, colour = "grey30",
                        linewidth = 0.2, position = "identity") +
        scale_fill_manual(
            values = c("Up in Post-IP" = COL_POSTIP_FILL,
                       "Up in Pre-IP"  = COL_PREIP_FILL),
            name = NULL) +
        labs(title    = sprintf("CRS distribution among eligible DEGs - %s", ds),
             subtitle = sprintf(
                 "Eligible = significant & pct_filter >= %.0f%% | N = %d",
                 PCT_FILTER * 100, nrow(elig)),
             x = "Composite Rank Score (CRS)", y = "Count") +
        theme_deg()
}



# =============================================================================
# CRS (composite rank score)
# =============================================================================
add_crs_columns <- function(df) {
    if (is.null(df) || nrow(df) == 0) return(df)
    df$avg_log2FC   <- as.numeric(df$avg_log2FC)
    df$p_val_adj_BH <- as.numeric(df$p_val_adj_BH)
    df$r_rb         <- as.numeric(df$r_rb)
    df$pct_fg       <- as.numeric(df$pct_fg)
    df$pct_bg       <- as.numeric(df$pct_bg)
    df$significant  <- as.logical(df$significant)
    n <- nrow(df)

    df$pct_filter_pass <- (df$pct_fg >= PCT_FILTER | df$pct_bg >= PCT_FILTER)
    df$is_mt_gene      <- df$gene %in% MT_GENES

    padj_raw     <- rank(df$p_val_adj_BH, ties.method = "average", na.last = TRUE)
    df$padj_rank <- round((n + 1) - padj_raw, 6)
    df$auc_rank  <- round(rank(abs(df$r_rb), ties.method = "average",
                               na.last = TRUE), 6)
    df$fc_rank   <- round(rank(abs(df$avg_log2FC), ties.method = "average",
                               na.last = TRUE), 6)

    # With one testable gene CRS is degenerate; it is still computed.
    if (n == 1L)
        cat("    NOTE: 1 testable gene; CRS is degenerate for this comparison\n")
    s_auc  <- -log(1 - df$auc_rank / (n + 1))
    s_fc   <- -log(1 - df$fc_rank  / (n + 1))
    s_max  <- log(n + 1)
    df$CRS <- round((s_auc + s_fc) / (2 * s_max), 6)
    df$CRS_signed <- round(sign(df$avg_log2FC) * df$CRS, 6)

    df$eligible <- !is.na(df$significant) & df$significant &
                   !is.na(df$pct_filter_pass) & df$pct_filter_pass
    df$top1000_up <- FALSE; df$top1000_down <- FALSE
    df$rank_up <- NA_integer_; df$rank_down <- NA_integer_
    idx_up   <- which(df$eligible & !is.na(df$avg_log2FC) & df$avg_log2FC > 0)
    idx_down <- which(df$eligible & !is.na(df$avg_log2FC) & df$avg_log2FC < 0)
    if (length(idx_up) > 0) {
        top_up <- idx_up[order(-df$CRS[idx_up])[seq_len(min(TOP_N_CRS, length(idx_up)))]]
        df$top1000_up[top_up] <- TRUE
        df$rank_up[top_up]    <- rank(-df$CRS[top_up], ties.method = "min")
    }
    if (length(idx_down) > 0) {
        top_dn <- idx_down[order(-df$CRS[idx_down])[seq_len(min(TOP_N_CRS, length(idx_down)))]]
        df$top1000_down[top_dn] <- TRUE
        df$rank_down[top_dn]    <- rank(-df$CRS[top_dn], ties.method = "min")
    }
    df
}


# =============================================================================
# SECTION 0b: The pipeline for one arm
# =============================================================================

run_one_arm <- function(arm_tag, config_file, out_dir, file_tag) {

    config_path <- file.path(CEP_HOME, config_file)
    if (!file.exists(config_path))
        stop("Configuration not found: ", config_path, call. = FALSE)

    .p10g_preamble <- character(0)
    assign("CEP_ARM", arm_tag, envir = globalenv())
    on.exit(suppressWarnings(rm("CEP_ARM", envir = globalenv())), add = TRUE)

    .p10g_preamble <- tryCatch(
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

    .p10g_needed <- c("ARM_TAG", "R_EXPORT_ROOT", "GOI", "DFG", "NORM_TAG",
                      "PRIMARY_STYPE", "NEOC_DIR", "NONNC_DIR",
                      "NEOC_REGIONS", "NONNC_REGIONS",
                      "SHARED_DIR", "NORM_TARGET_SUM", "STYPE_COL", "DONOR_COL")
    .p10g_absent <- .p10g_needed[!vapply(.p10g_needed, exists, logical(1))]
    if (length(.p10g_absent) > 0)
        stop("The configuration loaded from\n    ", config_path,
             "\n  does not define: ", paste(.p10g_absent, collapse = ", "), call. = FALSE)
    if (!identical(ARM_TAG, arm_tag))
        stop("The configuration reports ARM_TAG '", ARM_TAG,
             "' but this arm is '", arm_tag, "'.", call. = FALSE)

    .p10g_preamble <- c(.p10g_preamble, utils::capture.output(
        cat("  Config OK: ", ARM_TAG, " | ", GOI, " + ", length(DFG), " DFG | ",
            NORM_TAG, "\n", sep = "")))

    set.seed(V10G_SEED)
    DATASETS <- c("NeoC", "NonNeoC")

    # == Inputs ===============================================================
    V1_DIR    <- cep_arm_dir_for("Part10e_Validation1", arm_tag)
    V1_PREFIX <- paste0(cep_out_prefix("Part10e", NORM_TAG,
                                       cep_file_tag_for(arm_tag)), "_V1")
    p10g_v1_xlsx <- function(ds)
        file.path(V1_DIR, paste0(V1_PREFIX, "_", ds,
                                 "_SK_interaction_analysis.xlsx"))
    p10g_v1_rds <- function(ds)
        file.path(V1_DIR, paste0(V1_PREFIX, "_", ds, "_results.rds"))
    ds_dir <- function(ds) if (identical(ds, "NeoC")) NEOC_DIR else NONNC_DIR

    # == IP_y =================================================================
    # A cell is Pre-IP when its cDFG lies below IP_y = f(IP_x); GOI is never split.
    p10g_iprs_rds <- function() {
        # Same four candidate locations as Parts 10b to 10e.
        cand <- c(
            file.path(cep_arm_dir_for("Part10b_SKCS", arm_tag),
                      paste0(cep_out_prefix("Part10b", NORM_TAG,
                                            cep_file_tag_for(arm_tag)),
                             "_IPRS_results.rds")),
            file.path(cep_arm_dir_for("Part10b_IPRS", arm_tag),
                      paste0(cep_out_prefix("Part10b_IPRS", NORM_TAG,
                                            cep_file_tag_for(arm_tag)),
                             "_results.rds")),
            file.path(cep_arm_dir_for("Part10b", arm_tag),
                      paste0(cep_out_prefix("Part10b", NORM_TAG,
                                            cep_file_tag_for(arm_tag)),
                             "_IPRS_results.rds")),
            file.path(R_EXPORT_ROOT, "Part10b_IPRS",
                      "Part10b_scran_IPRS_results.rds"))
        hit <- cand[file.exists(cand)]
        if (length(hit) == 0L) return(list(path = NA_character_, obj = NULL,
                                           note = paste(basename(cand),
                                                        collapse = " / ")))
        list(path = hit[1], obj = tryCatch(readRDS(hit[1]),
                                           error = function(e) NULL), note = "")
    }
    .p10g_iprs <- p10g_iprs_rds()
    p10g_IP_y <- function(ds) {
        o <- .p10g_iprs$obj
        if (is.null(o)) return(list(value = NA_real_, field = NA_character_))
        for (f in c("inflection_points_DFG_y", "inflection_points_DFG_y_full")) {
            v <- o[[f]]
            if (!is.null(v) && ds %in% names(v)) {
                x <- suppressWarnings(as.numeric(v[[ds]]))
                if (is.finite(x)) return(list(value = x, field = f))
            }
        }
        list(value = NA_real_, field = NA_character_)
    }
    p10g_IP_x <- function(ds) {
        o <- .p10g_iprs$obj
        if (is.null(o)) return(NA_real_)
        v <- o[["inflection_points_DFG"]]
        if (!is.null(v) && ds %in% names(v))
            return(suppressWarnings(as.numeric(v[[ds]])))
        NA_real_
    }

    OUT_DIR <- out_dir
    dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
    FILE_PREFIX <- cep_out_prefix(PART_TAG, NORM_TAG, cep_file_tag_for(arm_tag))

    # == Logging ==============================================================
    LOG_PATH_10G <- file.path(OUT_DIR, paste0(FILE_PREFIX, "_console_log_",
                                              format(Sys.time(), "%Y%m%d_%H%M%S"),
                                              ".txt"))
    MSG_PATH_10G <- file.path(OUT_DIR, paste0(FILE_PREFIX,
                                              "_console_log_messages.tmp"))
    while (sink.number() > 0) sink()
    if (sink.number(type = "message") != 2L) sink(type = "message")
    .p10g_log_con <- file(LOG_PATH_10G, open = "wt")
    .p10g_msg_con <- file(MSG_PATH_10G, open = "wt")
    sink(.p10g_log_con, append = FALSE, split = TRUE)
    sink(.p10g_msg_con, type = "message")
    on.exit({
        tryCatch({
            if (!isTRUE(.p10g_completed))
                cat("\n", strrep("!", 72),
                    "\n  Arm did not complete; the error is printed to the console.\n",
                    strrep("!", 72), "\n", sep = "")
            if (sink.number(type = "message") != 2L) sink(type = "message")
            while (sink.number() > 0) sink()
            close(.p10g_msg_con); close(.p10g_log_con)
            if (file.exists(MSG_PATH_10G)) {
                .msg <- readLines(MSG_PATH_10G, warn = FALSE)
                if (length(.msg) > 0)
                    cat(c("", strrep("=", 72),
                          "  MESSAGE STREAM (warnings, messages, errors)",
                          strrep("=", 72), .msg, ""),
                        sep = "\n", file = LOG_PATH_10G, append = TRUE)
                unlink(MSG_PATH_10G)
            }
        }, error = function(e) invisible(NULL))
    }, add = TRUE)

    .p10g_completed <- FALSE
    .p10g_warn_opt  <- getOption("warn")
    .p10g_pm_opt    <- getOption("warnPartialMatchDollar")
    options(warn = 1)
    # Warn on partial `$` matches.
    options(warnPartialMatchDollar = TRUE)
    on.exit({
        options(warn = .p10g_warn_opt)
        options(warnPartialMatchDollar = .p10g_pm_opt)
    }, add = TRUE)
    cat(sprintf("Part10g DEGs Part 1 - log started: %s\n\n", Sys.time()))
    if (length(.p10g_preamble) > 0)
        cat(paste(.p10g_preamble, collapse = "\n"), "\n\n", sep = "")

    cat("\n", strrep("=", 72), "\n")
    cat("  Part10g - Post-IP vs Pre-IP DEGs within the CMRF module\n")
    cat("  Arm          :", ARM_TAG, "\n")
    cat("  GOI          :", GOI, "|", NORM_TAG, "\n")
    cat("  IP arms      : inherited (group_refined -> IP_status -> cDFG < IP_y);",
        "never split on GOI\n")
    cat("  Units        : cell (comparison) and donor (PRIMARY, paired pseudobulk)\n")
    cat("  Donor pairing: >=", V10G_MIN_CELLS_PER_DONOR,
        "cell(s) per arm, >=", V10G_MIN_DONORS, "paired donors\n")
    cat("  V1 input     :", V1_DIR, "\n")
    cat("  Output dir   :", OUT_DIR, "\n")
    cat(strrep("=", 72), "\n\n")
    cat("  NOTE: library nests within donor; batch and donor are not separable.\n\n")

    INCOMPLETE_MARKER_10G <- file.path(OUT_DIR, "_RUN_INCOMPLETE.txt")
    writeLines(c("This arm did not finish.", "",
                 paste("Arm     :", arm_tag),
                 paste("Part    :", PART_DIR_TAG),
                 paste("Started :", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), "",
                 "Any files in this folder are from a partial run.",
                 "This file is removed automatically when the arm completes."),
               INCOMPLETE_MARKER_10G)

    # =========================================================================
    # Guardrail ledger
    # =========================================================================
    GR_LOG_10G <- file.path(OUT_DIR, paste0(FILE_PREFIX, "_Guardrail_Log.txt"))
    .gr_log <- list(); .gr_pass <- 0L; .gr_fail <- 0L; .gr_warn <- 0L
    gr10g <- function(id, status, detail) {
        # A detail that is not a single string is replaced, so the check still appears.
        if (length(detail) != 1L || is.na(detail))
            detail <- "[detail unavailable]"
        tag <- switch(status, PASS = "[PASS]", FAIL = "[FAIL]", WARN = "[WARN]")
        cat(sprintf("  %s  %-20s  %s\n", tag, id, detail))
        .gr_log[[length(.gr_log) + 1L]] <<- list(id = id, status = status,
                                                 detail = detail)
        if (status == "PASS") .gr_pass <<- .gr_pass + 1L
        if (status == "FAIL") .gr_fail <<- .gr_fail + 1L
        if (status == "WARN") .gr_warn <<- .gr_warn + 1L
    }
    gr10g_write <- function() {
        con <- file(GR_LOG_10G, "wt")
        writeLines(c(strrep("=", 72),
                     paste0("Part 10g DEGs Part 1 - Guardrail Log [", NORM_TAG, "]"),
                     paste("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
                     paste("Arm :", ARM_TAG),
                     "Primary inferential unit: donor (paired pseudobulk)",
                     strrep("=", 72), ""), con)
        for (it in .gr_log)
            writeLines(sprintf("  [%s]  %-20s  %s", it$status, it$id, it$detail), con)
        writeLines(c("", strrep("=", 72),
                     sprintf("  Total : %d", .gr_pass + .gr_fail + .gr_warn),
                     sprintf("  PASS  : %d", .gr_pass),
                     sprintf("  WARN  : %d", .gr_warn),
                     sprintf("  FAIL  : %d", .gr_fail),
                     if (.gr_fail == 0 && .gr_warn == 0)
                         "\n  Result: no failures, no warnings"
                     else if (.gr_fail == 0)
                         "\n  Result: no failures; warnings present"
                     else "\n  Result: failures present",
                     strrep("=", 72)), con)
        close(con)
        cat("  Guardrail log written:", GR_LOG_10G, "\n")
    }

    cat("  PRE-ANALYSIS GUARDRAILS\n\n")
    for (ds in DATASETS) {
        x <- p10g_v1_xlsx(ds); r <- p10g_v1_rds(ds)
        gr10g(paste0("GR-10G-01-", ds),
              if (file.exists(x) && file.exists(r)) "PASS" else "FAIL",
              if (file.exists(x) && file.exists(r))
                  sprintf("Validation1 workbook and RDS present in %s",
                          basename(V1_DIR))
              else sprintf("missing: %s",
                           paste(basename(c(x, r)[!file.exists(c(x, r))]),
                                 collapse = ", ")))
        dd <- file.path(ds_dir(ds), "deg_support")
        need <- file.path(dd, c("norm.mtx", "norm_barcodes.tsv",
                                "norm_features.tsv"))
        gr10g(paste0("GR-10G-02-", ds),
              if (all(file.exists(need))) "PASS" else "FAIL",
              if (all(file.exists(need)))
                  sprintf("deg_support present: %s", dd)
              else sprintf("deg_support incomplete: %s",
                           paste(basename(need[!file.exists(need)]),
                                 collapse = ", ")))
    }
    gr10g("GR-10G-03-presto",
          if (.v10g_presto_ok) "PASS" else "WARN",
          if (.v10g_presto_ok) "presto available for the cell unit"
          else "presto absent; the cell unit uses the base-R rank-sum fallback (same p-values, slower)")
    cat("\n")

    # =========================================================================
    # Donor lookup
    # =========================================================================
    .p10g_donor_cache <- list()
    p10g_donor_map <- function(ds) {
        if (!is.null(.p10g_donor_cache[[ds]])) return(.p10g_donor_cache[[ds]])
        gdir <- ds_dir(ds)
        regs <- if (identical(ds, "NeoC")) NEOC_REGIONS else NONNC_REGIONS
        st   <- gsub("-", "_", gsub("/", "_", PRIMARY_STYPE))
        bc <- character(0); dn <- character(0)
        for (rg in regs) {
            gp <- file.path(gdir, "gam_data", paste0(rg, "_", st, "_selected.csv"))
            if (!file.exists(gp)) next
            d <- read.csv(gp, stringsAsFactors = FALSE, check.names = FALSE)
            if (!all(c("barcode", "donor") %in% names(d)))
                stop(sprintf("Region %s: gam_data lacks barcode or donor.", rg),
                     call. = FALSE)
            bc <- c(bc, d$barcode); dn <- c(dn, as.character(d$donor))
        }
        if (length(bc) == 0L) stop("No gam_data CSVs for ", ds, call. = FALSE)
        k <- !duplicated(bc)
        m <- setNames(dn[k], bc[k])
        .p10g_donor_cache[[ds]] <<- m
        m
    }

    # =========================================================================
    # Normalised matrix from Part 8's deg_support export
    # =========================================================================
    p10g_load_norm <- function(ds) {
        dd  <- file.path(ds_dir(ds), "deg_support")
        mp  <- file.path(dd, "norm.mtx")
        bp  <- file.path(dd, "norm_barcodes.tsv")
        fp  <- file.path(dd, "norm_features.tsv")
        if (!all(file.exists(c(mp, bp, fp))))
            stop("deg_support incomplete under: ", dd, call. = FALSE)
        X  <- as(readMM(mp), "CsparseMatrix")      # cells x genes on disk
        bc <- readLines(bp, warn = FALSE)
        ft <- read.table(fp, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
        if (ncol(ft) < 2L)
            stop("norm_features.tsv has ", ncol(ft),
                 " column(s); gene symbols are expected in column 2.",
                 call. = FALSE)
        if (nrow(X) != length(bc) || ncol(X) != nrow(ft))
            stop(sprintf("norm.mtx is %d x %d but there are %d barcode(s) and %d feature(s)",
                         nrow(X), ncol(X), length(bc), nrow(ft)), call. = FALSE)
        rownames(X) <- bc
        colnames(X) <- make.unique(as.character(ft[, 2]))
        mx <- max(X)
        E  <- Matrix::t(X)                          # genes x cells
        rm(X); gc(verbose = FALSE)
        cat(sprintf("    norm.mtx: %d gene(s) x %d cell(s) | max = %.4f\n",
                    nrow(E), ncol(E), mx))
        list(mat = E, max = mx)
    }

    # =========================================================================
    # One comparison, at both units
    # =========================================================================
    # Fold-change gate: DEG_SIG_LOGFC_BG for tags ending in _vs_BG, else DEG_SIG_LOGFC.
    sig_thr <- function(tg)
        if (grepl("_vs_BG$", tg)) DEG_SIG_LOGFC_BG else DEG_SIG_LOGFC

    run_comparison <- function(expr, donor_vec, fg_bc, bg_bc, ds, tag,
                               fg_lab, bg_lab, sig_logfc) {

        fg <- intersect(fg_bc, colnames(expr))
        bg <- intersect(bg_bc, colnames(expr))
        fg <- setdiff(fg, bg)
        n_fg <- length(fg); n_bg <- length(bg)
        cat(sprintf("\n  [%s] %s: fg n=%d | bg n=%d\n", ds, tag, n_fg, n_bg))
        if (n_fg < DEG_MIN_CELLS_FG || n_bg < DEG_MIN_CELLS_BG) {
            cat(sprintf("    SKIP: below the %d/%d cell minimum\n",
                        DEG_MIN_CELLS_FG, DEG_MIN_CELLS_BG))
            return(NULL)
        }

        cells <- c(fg, bg)
        is_fg <- c(rep(TRUE, n_fg), rep(FALSE, n_bg))
        mat   <- expr[, cells, drop = FALSE]
        dvec  <- unname(donor_vec[cells])

        pct_fg <- Matrix::rowMeans(mat[, is_fg,  drop = FALSE] > 0)
        pct_bg <- Matrix::rowMeans(mat[, !is_fg, drop = FALSE] > 0)
        fc_all <- v10g_log2fc(mat[, is_fg, drop = FALSE],
                              mat[, !is_fg, drop = FALSE])
        keep <- (pct_fg >= DEG_MIN_PCT) | (pct_bg >= DEG_MIN_PCT)
        if (DEG_LOGFC_THRESH > 0) keep <- keep & (abs(fc_all) >= DEG_LOGFC_THRESH)
        cat(sprintf("    Detection filter: %d -> %d gene(s)\n",
                    length(keep), sum(keep)))
        if (!any(keep)) { cat("    SKIP: no gene passed the filter\n"); return(NULL) }
        mat <- mat[keep, , drop = FALSE]

        mean_fg <- Matrix::rowMeans(mat[, is_fg,  drop = FALSE])
        mean_bg <- Matrix::rowMeans(mat[, !is_fg, drop = FALSE])

        out <- list()

        # == Cell unit ========================================================
        t0 <- proc.time()[["elapsed"]]
        rs <- v10g_ranksum(mat, is_fg)
        cat(sprintf("    cell unit: %d gene(s) in %.1f s\n",
                    nrow(rs), proc.time()[["elapsed"]] - t0))
        g <- rs$gene
        padj <- p.adjust(rs$pval, method = "BH")
        cell_df <- data.frame(
            gene = g, dataset = ds, comparison = tag, unit = "cell",
            n_fg = n_fg, n_bg = n_bg, n_donors = NA_integer_,
            avg_log2FC = round(as.numeric(fc_all[g]), 5),
            pct_fg = round(as.numeric(pct_fg[g]), 4),
            pct_bg = round(as.numeric(pct_bg[g]), 4),
            mean_expr_fg = round(as.numeric(mean_fg[g]), 5),
            mean_expr_bg = round(as.numeric(mean_bg[g]), 5),
            p_val = signif(rs$pval, 5), p_val_adj_BH = signif(padj, 5),
            auc = round(rs$auc, 6), r_rb = round(2 * (rs$auc - 0.5), 6),
            stringsAsFactors = FALSE)
        cell_df$significant <- !is.na(cell_df$p_val_adj_BH) &
            cell_df$p_val_adj_BH < DEG_SIG_PADJ &
            abs(cell_df$avg_log2FC) >= sig_logfc
        cell_df$direction <- ifelse(cell_df$avg_log2FC > 0,
                                    paste0("up_in_", fg_lab),
                                    paste0("up_in_", bg_lab))
        cell_df <- add_crs_columns(cell_df)
        cell_df <- cell_df[order(cell_df$p_val_adj_BH,
                                 -abs(cell_df$avg_log2FC)), ]
        rownames(cell_df) <- NULL
        out$cell <- cell_df
        cat(sprintf("    cell unit: %d significant (BH < %.2f, |FC| >= %.2f)\n",
                    sum(cell_df$significant, na.rm = TRUE), DEG_SIG_PADJ,
                    sig_logfc))

        # == Donor unit (primary) =============================================
        if (isTRUE(V10G_DONOR_UNIT) && any(!is.na(dvec))) {
            t0 <- proc.time()[["elapsed"]]
            pb <- v10g_pseudobulk(mat, is_fg, dvec)
            cat(sprintf("    donor unit: %d donor(s) total, %d paired (>= %d cell(s) each arm), %d fg-only, %d bg-only\n",
                        pb$info$n_donors_total, pb$info$n_donors_paired,
                        pb$info$min_cells, pb$info$n_donors_fg_only,
                        pb$info$n_donors_bg_only))
            if (is.null(pb$res)) {
                cat("    donor unit SKIPPED:", pb$info$note, "\n")
                out$donor_info <- pb$info
            } else {
                gd <- pb$res$gene
                padj_d <- p.adjust(pb$res$pval, method = "BH")
                don_df <- data.frame(
                    gene = gd, dataset = ds, comparison = tag, unit = "donor",
                    n_fg = n_fg, n_bg = n_bg, n_donors = pb$res$n_donors,
                    avg_log2FC = round(pb$res$avg_log2FC, 5),
                    pct_fg = round(as.numeric(pct_fg[gd]), 4),
                    pct_bg = round(as.numeric(pct_bg[gd]), 4),
                    mean_expr_fg = round(as.numeric(mean_fg[gd]), 5),
                    mean_expr_bg = round(as.numeric(mean_bg[gd]), 5),
                    p_val = signif(pb$res$pval, 5),
                    p_val_adj_BH = signif(padj_d, 5),
                    auc = NA_real_, r_rb = round(pb$res$r_rb, 6),
                    stringsAsFactors = FALSE)
                don_df$significant <- !is.na(don_df$p_val_adj_BH) &
                    don_df$p_val_adj_BH < DEG_SIG_PADJ &
                    abs(don_df$avg_log2FC) >= sig_logfc
                don_df$direction <- ifelse(don_df$avg_log2FC > 0,
                                           paste0("up_in_", fg_lab),
                                           paste0("up_in_", bg_lab))
                don_df <- add_crs_columns(don_df)
                don_df <- don_df[order(don_df$p_val_adj_BH,
                                       -abs(don_df$avg_log2FC)), ]
                rownames(don_df) <- NULL
                out$donor <- don_df
                if (isTRUE(V10G_STORE_DONOR_FC) &&
                    tag %in% V10G_FC_CONTRASTS && !is.null(pb$fc)) {
                    # genes x paired donors, dimnames carried through so the
                    # matrix can be joined to don_df by gene downstream
                    out$donor_fc <- pb$fc
                    cat(sprintf("    donor_fc stored: %d gene(s) x %d donor(s)\n",
                                nrow(pb$fc), ncol(pb$fc)))
                }
                out$donor_info <- pb$info
                cat(sprintf("    donor unit: %d significant of %d testable, on %d donor(s)\n",
                            sum(don_df$significant, na.rm = TRUE), nrow(don_df),
                            pb$res$n_donors[1]))
            }
        }

        # == Testable universes, one per unit =================================
        out$universe <- do.call(rbind, lapply(names(out)[names(out) %in%
                                                         c("cell", "donor")],
            function(u) data.frame(dataset = ds, comparison = tag, unit = u,
                                   gene = out[[u]]$gene,
                                   stringsAsFactors = FALSE)))
        out
    }


    # =========================================================================
    # SECTION B: Main function (per dataset)
    # =========================================================================

    run_V10G_dataset <- function(ds) {

        cat("\n", strrep("=", 72), "\n")
        cat("  Part10g:", ds, " [", ARM_TAG, "]\n")
        cat(strrep("=", 72), "\n\n")

        # == Step 1: Validation1 labels =======================================
        cat("  STEP 1 - Validation1 cell annotations\n")
        xp <- p10g_v1_xlsx(ds)
        if (!file.exists(xp))
            stop("Validation1 workbook not found: ", xp, call. = FALSE)
        rc <- as.data.frame(
            readxl::read_excel(xp, sheet = "Step4_Refined_Cells"),
            stringsAsFactors = FALSE)
        # barcode and cell_class are required; u_all is optional (reference only).
        need <- c("barcode", "cell_class")
        miss <- setdiff(need, names(rc))
        if (length(miss) > 0)
            stop("Step4_Refined_Cells lacks: ", paste(miss, collapse = ", "),
                 call. = FALSE)
        rc$barcode    <- as.character(rc$barcode)
        rc$cell_class <- as.character(rc$cell_class)
        rc$u_all <- if ("u_all" %in% names(rc))
            suppressWarnings(as.numeric(rc$u_all)) else NA_real_
        if (anyDuplicated(rc$barcode) > 0L)
            stop("Step4_Refined_Cells has duplicated barcodes.", call. = FALSE)

        # Keep only rows that can be classified.
        n_class_na <- sum(is.na(rc$cell_class))
        if (n_class_na > 0)
            cat(sprintf("    %d row(s) with a missing cell_class dropped\n",
                        n_class_na))
        rc <- rc[!is.na(rc$cell_class), , drop = FALSE]
        for (cls in sort(unique(rc$cell_class)))
            cat(sprintf("    %-25s : %d\n", cls, sum(rc$cell_class == cls)))
        is_module <- rc$cell_class %in% MODULE_CLASSES
        cat(sprintf("    Module (OrigRef + Recovered) : %d\n", sum(is_module)))

        # == Step 2: IP_x, u_IP and IP_y ======================================
        cat("\n  STEP 2 - Inflection point (IP_x reported, IP_y classifies)\n")
        # The Validation 1 RDS supplies only u_IP (reference); a missing file is not fatal.
        rp <- p10g_v1_rds(ds)
        u_ip <- NA_real_
        if (!file.exists(rp)) {
            cat(sprintf("    NOTE: Validation1 RDS not found (%s); GOI-axis contrast skipped.\n",
                        basename(rp)))
        } else {
            v1 <- tryCatch(readRDS(rp), error = function(e) NULL)
            if (is.null(v1))
                cat(sprintf("    NOTE: Validation1 RDS could not be read (%s); GOI-axis contrast skipped.\n",
                            basename(rp)))
            else u_ip <- suppressWarnings(as.numeric(v1$u_IP)[1])
            rm(v1); gc(verbose = FALSE)
        }
        # u_IP is reference only; an absent or odd value is reported and skipped.
        if (length(u_ip) != 1L || !is.finite(u_ip) || u_ip <= 0 || u_ip >= 1) {
            cat(sprintf("    NOTE: Validation1 gave no usable u_IP (%s); the",
                        " GOI-axis contrast line will be skipped\n"),
                if (length(u_ip) == 1L) as.character(u_ip) else "absent")
            u_ip <- NA_real_
        }

        ipy  <- p10g_IP_y(ds)
        ip_y <- ipy$value
        ip_x <- p10g_IP_x(ds)
        cat(sprintf("    IP_x (reported, GOI axis)  = %s   [R(x)=0 crossing; not the classifier]\n",
                    if (is.finite(ip_x)) sprintf("%.6f", ip_x) else "unavailable"))
        cat(sprintf("    u_IP (copula rank of IP_x) = %s   [reference only]\n",
                    if (is.finite(u_ip)) sprintf("%.6f", u_ip) else "unavailable"))
        cat(sprintf("    IP_y = f(IP_x) (cDFG axis) = %s   [classifier: Pre-IP when cDFG < IP_y]\n",
                    if (is.finite(ip_y)) sprintf("%.6f", ip_y) else "unavailable"))
        cat(sprintf("    IP_y source: %s\n",
                    if (is.finite(ip_y))
                        sprintf("%s[[\"%s\"]] in %s", ipy$field, ds,
                                basename(.p10g_iprs$path))
                    else if (is.na(.p10g_iprs$path))
                        paste0("IPRS RDS not found; looked for ",
                               .p10g_iprs$note)
                    else if (is.null(.p10g_iprs$obj))
                        paste0("IPRS RDS found but unreadable: ",
                               basename(.p10g_iprs$path))
                    else sprintf(paste0("IPRS RDS %s carries no",
                                        " inflection_points_DFG_y[[\"%s\"]]"),
                                 basename(.p10g_iprs$path), ds)))

        # == Step 3: IP arms, inherited from Validation 1 =====================
        # Order: group_refined, then IP_status, then cDFG at IP_y (never the GOI axis).
        cat("\n  STEP 3 - IP arms (inherited from Part 10e / Validation 1)\n")

        pick_arm <- function(v) ifelse(is.na(v), NA_character_,
                        ifelse(grepl("Post-IP", v, fixed = TRUE), "Post-IP",
                        ifelse(grepl("Pre-IP",  v, fixed = TRUE), "Pre-IP",
                               NA_character_)))
        arm_gr <- if ("group_refined" %in% names(rc))
            pick_arm(as.character(rc$group_refined)) else rep(NA_character_, nrow(rc))
        arm_ip <- if ("IP_status" %in% names(rc))
            pick_arm(as.character(rc$IP_status)) else rep(NA_character_, nrow(rc))
        # Fallback and cross-check: cDFG_expr < IP_y -> Pre-IP (Validation 1's rule).
        has_cdfg <- "cDFG_expr" %in% names(rc) && is.finite(ip_y)
        arm_y <- if (has_cdfg) {
            cv <- suppressWarnings(as.numeric(rc$cDFG_expr))
            ifelse(is.na(cv), NA_character_,
                   ifelse(cv < ip_y, "Pre-IP", "Post-IP"))
        } else rep(NA_character_, nrow(rc))

        rc$arm        <- ifelse(!is.na(arm_gr), arm_gr,
                         ifelse(!is.na(arm_ip), arm_ip, arm_y))
        rc$arm_source <- ifelse(!is.na(arm_gr), "group_refined",
                         ifelse(!is.na(arm_ip), "IP_status",
                         ifelse(!is.na(arm_y), "cDFG_at_IP_y", NA_character_)))
        rc$arm_source[is.na(rc$arm)] <- NA_character_

        mod_i    <- rc$cell_class %in% MODULE_CLASSES
        arm_src  <- table(factor(rc$arm_source[mod_i],
                                 levels = c("group_refined", "IP_status",
                                            "cDFG_at_IP_y")))
        n_y_fall <- as.integer(arm_src[["cDFG_at_IP_y"]])
        n_unlab  <- sum(mod_i & is.na(rc$arm))

        # Agreement of the inherited label with cDFG_expr < IP_y.
        cmp <- mod_i & !is.na(rc$arm) & !is.na(arm_y) &
               rc$arm_source != "cDFG_at_IP_y"
        arm_agree <- if (any(cmp)) sum(rc$arm[cmp] == arm_y[cmp]) / sum(cmp)
                     else NA_real_
        arm_n_cmp <- sum(cmp)
        cat(sprintf("    Arm source: group_refined %d | IP_status %d | cDFG@IP_y fallback %d | unlabelled %d\n",
                    as.integer(arm_src[["group_refined"]]),
                    as.integer(arm_src[["IP_status"]]), n_y_fall, n_unlab))
        cat(sprintf("    Inherited label vs `cDFG_expr < IP_y`: %s on %d cell(s)\n",
                    if (is.na(arm_agree))
                        "not checked (no cDFG_expr column, or no IP_y)"
                    else sprintf("%.2f%% agreement", 100 * arm_agree), arm_n_cmp))
        if (!is.na(arm_agree) && arm_agree < 1)
            cat("    They differ; the inherited label is used.\n")

        # Stop rather than split on the GOI axis when no arm can be assigned.
        if (sum(mod_i) == 0L)
            stop("No module cell (Original + Refined or Recovered) in ",
                 "Step4_Refined_Cells for ", ds, ".", call. = FALSE)
        if (n_unlab == sum(mod_i))
            stop("No module cell could be assigned an IP arm for ", ds,
                 " (cDFG_expr column present: ",
                 "cDFG_expr" %in% names(rc), "; IP_y finite: ", is.finite(ip_y),
                 ").", call. = FALSE)

        # GOI-axis split at u_IP, reported only.
        arm_u <- if (is.finite(u_ip))
            ifelse(is.na(rc$u_all), NA_character_,
                   ifelse(rc$u_all < u_ip, "Pre-IP", "Post-IP"))
        else rep(NA_character_, nrow(rc))
        cmp_u <- mod_i & !is.na(rc$arm) & !is.na(arm_u)
        arm_agree_u <- if (any(cmp_u))
            sum(rc$arm[cmp_u] == arm_u[cmp_u]) / sum(cmp_u) else NA_real_
        cat(sprintf("    For reference only: the GOI-axis split at u_IP would disagree on %s of %d cell(s) (%s)\n",
                    if (is.na(arm_agree_u)) "NA"
                    else as.character(as.integer(round((1 - arm_agree_u) * sum(cmp_u)))),
                    sum(cmp_u),
                    if (is.na(arm_agree_u)) "not comparable"
                    else sprintf("%.2f%% agreement", 100 * arm_agree_u)))
        for (cls in MODULE_CLASSES) {
            i <- rc$cell_class == cls
            if (!any(i)) next
            cat(sprintf("    %-25s : Pre-IP=%d  Post-IP=%d  unlabelled=%d\n",
                        cls,
                        sum(rc$arm[i] == "Pre-IP",  na.rm = TRUE),
                        sum(rc$arm[i] == "Post-IP", na.rm = TRUE),
                        sum(is.na(rc$arm[i]))))
        }

        # Non-module cells: arm from the cDFG rule; provenance printed separately.
        nonmod_i <- !is.na(rc$cell_class) &
                    !(rc$cell_class %in% MODULE_CLASSES)
        nm_classes <- sort(unique(as.character(rc$cell_class[nonmod_i])))
        nm_src <- table(factor(rc$arm_source[nonmod_i],
                               levels = c("group_refined", "IP_status",
                                          "cDFG_at_IP_y")))
        nm_unlab <- sum(nonmod_i & is.na(rc$arm))
        cat(sprintf("    Non-module cells: %d in %d class(es): %s\n",
                    sum(nonmod_i), length(nm_classes),
                    if (length(nm_classes)) paste(nm_classes, collapse = ", ")
                    else "none"))
        .nm_unexpected <- setdiff(nm_classes, NONMODULE_CLASSES_EXPECTED)
        if (length(.nm_unexpected))
            cat(sprintf("    NOTE: %s not in NONMODULE_CLASSES_EXPECTED.\n",
                        paste(.nm_unexpected, collapse = ", ")))
        cat(sprintf(paste0("    Non-module arm source: group_refined %d |",
                           " IP_status %d | cDFG@IP_y %d | unlabelled %d\n"),
                    as.integer(nm_src[["group_refined"]]),
                    as.integer(nm_src[["IP_status"]]),
                    as.integer(nm_src[["cDFG_at_IP_y"]]), nm_unlab))
        for (cls in nm_classes) {
            i <- !is.na(rc$cell_class) & rc$cell_class == cls
            cat(sprintf("    %-25s : Pre-IP=%d  Post-IP=%d  unlabelled=%d\n",
                        cls,
                        sum(rc$arm[i] == "Pre-IP",  na.rm = TRUE),
                        sum(rc$arm[i] == "Post-IP", na.rm = TRUE),
                        sum(is.na(rc$arm[i]))))
        }

        # == Step 4: expression and donor =====================================
        cat("\n  STEP 4 - Normalised expression and donor\n")
        nm   <- p10g_load_norm(ds)
        expr <- nm$mat
        dmap <- p10g_donor_map(ds)
        donor_vec <- setNames(unname(dmap[colnames(expr)]), colnames(expr))
        n_no_donor <- sum(is.na(donor_vec))
        cat(sprintf("    Donors: %d distinct | %d cell(s) without one\n",
                    length(unique(donor_vec[!is.na(donor_vec)])), n_no_donor))

        mod_bc  <- rc$barcode[is_module & !is.na(rc$arm)]
        in_mat  <- intersect(mod_bc, colnames(expr))
        # Check the CP10K recipe before normalising any background cell.
        nv <- p10g_verify_cp10k(file.path(ds_dir(ds), "counts", "raw_selected"),
                                expr)
        cat(sprintf("    CP10K recipe vs deg_support/norm.mtx: %s\n",
                    if (is.finite(nv$dev))
                        sprintf("max |dev| = %.3e on %d gene(s) x %d cell(s)%s",
                                nv$dev, nv$n_genes, nv$n_cells,
                                if (nv$ok) "" else "  -- above tolerance")
                    else paste0("not verified: ", nv$note)))

        cat(sprintf("    Module cells with expression: %d of %d\n",
                    length(in_mat), length(mod_bc)))
        if (length(in_mat) < DEG_MIN_CELLS_FG + DEG_MIN_CELLS_BG)
            stop("Only ", length(in_mat),
                 " module cell(s) are present in norm.mtx.", call. = FALSE)

        post_bc <- rc$barcode[is_module & !is.na(rc$arm) & rc$arm == "Post-IP"]
        pre_bc  <- rc$barcode[is_module & !is.na(rc$arm) & rc$arm == "Pre-IP"]

        # Non-module halves (selected cells, so expected in norm.mtx); coverage is printed.
        nm_post_bc <- rc$barcode[nonmod_i & !is.na(rc$arm) &
                                 rc$arm == "Post-IP"]
        nm_pre_bc  <- rc$barcode[nonmod_i & !is.na(rc$arm) &
                                 rc$arm == "Pre-IP"]
        # Pooled non-module group.
        nm_all_bc  <- c(nm_pre_bc, nm_post_bc)
        nm_in_mat  <- length(intersect(c(nm_pre_bc, nm_post_bc),
                                       colnames(expr)))
        cat(sprintf(paste0("    Non-module cells with expression: %d of %d",
                           " (Pre-IP %d, Post-IP %d)\n"),
                    nm_in_mat, length(nm_pre_bc) + length(nm_post_bc),
                    length(nm_pre_bc), length(nm_post_bc)))
        # Donors clearing the per-donor floor, printed before Step 5.
        .nm_donor_count <- function(b) {
            b <- intersect(b, colnames(expr))
            if (!length(b)) return(0L)
            tb <- table(donor_vec[b])
            as.integer(sum(tb >= V10G_MIN_CELLS_PER_DONOR))
        }
        nm_donors_pre  <- .nm_donor_count(nm_pre_bc)
        nm_donors_post <- .nm_donor_count(nm_post_bc)
        nm_donors_all  <- .nm_donor_count(nm_all_bc)
        cat(sprintf(paste0("    Donors with >= %d non-module cell(s): Pre-IP",
                           " %d | Post-IP %d | pooled %d  (>= %d needed to",
                           " pair)\n"),
                    V10G_MIN_CELLS_PER_DONOR, nm_donors_pre, nm_donors_post,
                    nm_donors_all, V10G_MIN_DONORS))

        # == Step 5: comparison 1, the head-to-head ===========================
        cat("\n  STEP 5 - Comparisons\n")
        res <- list()
        res[["PostIP_vs_PreIP"]] <- run_comparison(
            expr, donor_vec, post_bc, pre_bc, ds, "PostIP_vs_PreIP",
            "PostIP", "PreIP", DEG_SIG_LOGFC)

        # == The five non-module contrasts ====================================
        # Each result is captured first, since res[[tag]] <- NULL deletes the element.
        nm_note    <- NA_character_
        nm_skipped <- character(0)
        if (!isTRUE(V10G_NONMOD_ENABLE)) {
            nm_note <- paste0("V10G_NONMOD_ENABLE is FALSE; the five",
                              " non-module contrasts were not attempted")
            cat("  ", nm_note, "\n", sep = "")
        } else if (sum(nonmod_i) == 0L) {
            nm_note <- paste0("Step4_Refined_Cells holds no cell outside",
                              " MODULE_CLASSES, so there is no non-module",
                              " population to contrast")
            cat("  ", nm_note, "\n", sep = "")
        } else {
            cat("\n  Non-module contrasts\n")
            nm_spec <- list(
                list(tag = V10G_NONMOD_CMP[1], fg = nm_post_bc,
                     bg = nm_pre_bc, fl = "NonModPostIP", bl = "NonModPreIP"),
                list(tag = V10G_NONMOD_CMP[2], fg = post_bc,
                     bg = nm_post_bc, fl = "PostIP",
                     bl = "NonModPostIP"),
                list(tag = V10G_NONMOD_CMP[3], fg = pre_bc,
                     bg = nm_pre_bc, fl = "PreIP",
                     bl = "NonModPreIP"),
                # Both module halves against the pooled non-module group.
                list(tag = V10G_NONMOD_CMP2[1], fg = post_bc,
                     bg = nm_all_bc, fl = "PostIP", bl = "NonModAll"),
                list(tag = V10G_NONMOD_CMP2[2], fg = pre_bc,
                     bg = nm_all_bc, fl = "PreIP", bl = "NonModAll"))
            for (sp in nm_spec) {
                .nm_res <- run_comparison(expr, donor_vec, sp$fg, sp$bg, ds,
                                          sp$tag, sp$fl, sp$bl, DEG_SIG_LOGFC)
                if (is.null(.nm_res)) nm_skipped <- c(nm_skipped, sp$tag)
                else res[[sp$tag]] <- .nm_res
            }
            nm_note <- sprintf(paste0("non-module = %d selected cell(s) that",
                                      " failed the CMRF band (%s); Pre-IP %d,",
                                      " Post-IP %d, %d with expression; %d of",
                                      " 5 contrast(s) ran"),
                               sum(nonmod_i),
                               if (length(nm_classes))
                                   paste(nm_classes, collapse = " + ")
                               else "no class",
                               length(nm_pre_bc), length(nm_post_bc),
                               nm_in_mat,
                               length(c(V10G_NONMOD_CMP, V10G_NONMOD_CMP2)) -
                                   length(nm_skipped))
        }

        # == Arm vs the non-M4-top background =================================
        # Built only if the CP10K recipe check passed.
        bg_note <- NA_character_; n_bg_cells <- 0L
        if (!isTRUE(V10G_BG_ENABLE)) {
            bg_note <- "V10G_BG_ENABLE is FALSE; the three background comparisons were not attempted"
            cat("  ", bg_note, "\n", sep = "")
        } else if (!isTRUE(nv$ok)) {
            bg_note <- sprintf(paste0("the CP10K recipe did not reproduce",
                                      " deg_support/norm.mtx (max |dev| = %s); background",
                                      " comparisons skipped"),
                               if (is.finite(nv$dev)) sprintf("%.3e", nv$dev)
                               else paste0("not computed: ", nv$note))
            cat("  ", bg_note, "\n", sep = "")
        } else {
            cat("\n  Background (non-M4-top ", PRIMARY_STYPE, ")\n", sep = "")
            bgl <- p10g_load_background(ds_dir(ds),
                                        file.path(SHARED_DIR, basename(ds_dir(ds))),
                                        colnames(expr), STYPE_COL, PRIMARY_STYPE,
                                        DONOR_COL)
            if (is.null(bgl$mat)) {
                bg_note <- sprintf("background not built: %s", bgl$note)
                cat("  ", bg_note, "\n", sep = "")
            } else {
                gg <- intersect(rownames(expr), rownames(bgl$mat))
                if (length(gg) < 100L) {
                    bg_note <- sprintf(paste0("only %d gene(s) shared between",
                                              " norm.mtx and raw_dataset; background",
                                              " comparisons skipped"), length(gg))
                    cat("  ", bg_note, "\n", sep = "")
                } else {
                    expr_bg <- cbind(expr[gg, , drop = FALSE],
                                     bgl$mat[gg, , drop = FALSE])
                    bg_bc   <- colnames(bgl$mat)
                    n_bg_cells <- length(bg_bc)
                    # Donors: selected cells from gam_data, background from obs_dataset.csv.
                    donor_bg <- setNames(unname(dmap[colnames(expr_bg)]),
                                         colnames(expr_bg))
                    if (!is.null(bgl$donor))
                        donor_bg[bg_bc] <- bgl$donor[match(bg_bc,
                                                           colnames(bgl$mat))]
                    n_bg_no_donor <- sum(is.na(donor_bg[bg_bc]))
                    cat(sprintf("    combined donors: %d distinct | %d background cell(s) without one\n",
                                length(unique(donor_bg[!is.na(donor_bg)])),
                                n_bg_no_donor))
                    cat(sprintf("    combined matrix: %d gene(s) x %d cell(s) (%d background)\n",
                                nrow(expr_bg), ncol(expr_bg), n_bg_cells))
                    # Capture each result first (res[[tag]] <- NULL deletes the element).
                    .bg_skipped <- character(0)
                    .bg_po <- run_comparison(
                        expr_bg, donor_bg, post_bc, bg_bc, ds, "PostIP_vs_BG",
                        "PostIP", "BG", DEG_SIG_LOGFC_BG)
                    if (is.null(.bg_po))
                        .bg_skipped <- c(.bg_skipped, "PostIP_vs_BG")
                    else res[["PostIP_vs_BG"]] <- .bg_po
                    .bg_pr <- run_comparison(
                        expr_bg, donor_bg, pre_bc, bg_bc, ds, "PreIP_vs_BG",
                        "PreIP", "BG", DEG_SIG_LOGFC_BG)
                    if (is.null(.bg_pr))
                        .bg_skipped <- c(.bg_skipped, "PreIP_vs_BG")
                    else res[["PreIP_vs_BG"]] <- .bg_pr
                    # Pooled non-module group against the same background, if it exists.
                    .nmbg_note <- NA_character_
                    if (!isTRUE(V10G_NONMOD_ENABLE) ||
                        length(nm_all_bc) == 0L) {
                        .nmbg_note <- if (!isTRUE(V10G_NONMOD_ENABLE))
                            paste("V10G_NONMOD_ENABLE is FALSE, so",
                                  "NMAll_vs_BG was not attempted either")
                        else paste("no non-module cell, so",
                                   "NMAll_vs_BG was not attempted")
                        cat("    ", .nmbg_note, "\n", sep = "")
                    } else {
                        .nmbg <- run_comparison(
                            expr_bg, donor_bg, nm_all_bc, bg_bc, ds,
                            V10G_NONMOD_BG, "NonModAll", "BG",
                            DEG_SIG_LOGFC_BG)
                        if (is.null(.nmbg))
                            .bg_skipped <- c(.bg_skipped, V10G_NONMOD_BG)
                        else res[[V10G_NONMOD_BG]] <- .nmbg
                    }
                    rm(expr_bg); gc(verbose = FALSE)
                    bg_note <- sprintf(paste0("background = %d non-selected %s",
                                              " cell(s) of %d in the dataset,",
                                              " from raw_dataset, normalised on",
                                              " the verified CP10K recipe; %d",
                                              " gene(s) shared with norm.mtx;",
                                              " %d of %d selected cell(s) are",
                                              " %s; %d of 3 background",
                                              " contrast(s) ran%s"),
                                       n_bg_cells, PRIMARY_STYPE, bgl$n_stype,
                                       length(gg), bgl$n_sel_stype,
                                       ncol(expr), PRIMARY_STYPE,
                                       sum(c("PostIP_vs_BG", "PreIP_vs_BG",
                                             V10G_NONMOD_BG) %in% names(res)),
                                       paste0(
                                           if (is.na(.nmbg_note)) ""
                                           else paste0(" -- ", .nmbg_note),
                                           if (length(.bg_skipped))
                                               paste0(" -- returned nothing: ",
                                                      paste(.bg_skipped,
                                                            collapse = ", "))
                                           else ""))
                }
            }
        }

        # == Step 6: figures ==================================================
        cat("\n  STEP 6 - Figures\n")
        p_arm <- plot_arm_counts(rc[is_module, , drop = FALSE], ds)
        save_pdf_png(p_arm, file.path(OUT_DIR, paste0(FILE_PREFIX, "_", ds,
                                                      "_arm_counts")),
                     w = 8, h = 6)
        show_plot(p_arm)
        # Figures for every comparison, each at its own fold-change gate.
        n_fig <- 0L
        for (tg in names(res)) {
            if (is.null(res[[tg]])) next
            thr <- sig_thr(tg)
            for (u in c("cell", "donor")) {
                # Cell-unit figures only when V10G_CELL_OUTPUTS is TRUE.
                if (identical(u, "cell") && !isTRUE(V10G_CELL_OUTPUTS)) next
                d <- res[[tg]][[u]]
                if (is.null(d) || nrow(d) == 0) next
                lab  <- sprintf("%s %s [%s unit]", ds, tg, u)
                base <- file.path(OUT_DIR,
                                  paste0(FILE_PREFIX, "_", ds, "_", tg, "_", u))
                for (fn in list(
                    list(p = plot_volcano(d, lab, thr), s = "_volcano",
                         w = 11, h = 9),
                    list(p = plot_deg_bar(d, lab), s = "_counts", w = 8, h = 6),
                    list(p = plot_crs_distribution(d, lab), s = "_CRS",
                         w = 9, h = 6))) {
                    if (is.null(fn$p)) next
                    save_pdf_png(fn$p, paste0(base, fn$s), w = fn$w, h = fn$h)
                    show_plot(fn$p); n_fig <- n_fig + 1L
                }
                # The dot plot needs both cell sets in `expr` (not the background contrasts).
                .dp_pair <- if (identical(tg, "PostIP_vs_PreIP"))
                                list(post_bc, pre_bc)
                            else if (identical(tg, V10G_NONMOD_CMP[1]))
                                list(nm_post_bc, nm_pre_bc)
                            else if (identical(tg, V10G_NONMOD_CMP[2]))
                                list(post_bc, nm_post_bc)
                            else if (identical(tg, V10G_NONMOD_CMP[3]))
                                list(pre_bc, nm_pre_bc)
                            else if (identical(tg, V10G_NONMOD_CMP2[1]))
                                list(post_bc, nm_all_bc)
                            else if (identical(tg, V10G_NONMOD_CMP2[2]))
                                list(pre_bc, nm_all_bc)
                            else NULL
                if (!is.null(.dp_pair)) {
                    dp <- tryCatch(plot_dotplot(d, expr, .dp_pair[[1]],
                                                .dp_pair[[2]], lab),
                                   error = function(e) {
                                       cat("    dot plot skipped:",
                                           conditionMessage(e), "\n"); NULL })
                    if (!is.null(dp)) {
                        save_pdf_png(dp, paste0(base, "_dotplot"), w = 10, h = 10)
                        show_plot(dp); n_fig <- n_fig + 1L
                    }
                }
            }
        }
        cat(sprintf("    %d figure(s) written across %d comparison(s)\n",
                    n_fig, sum(!vapply(res, is.null, logical(1)))))

        rm(expr); gc(verbose = FALSE)

        # == Summary ==========================================================
        summ <- do.call(rbind, lapply(names(res), function(tg) {
            r <- res[[tg]]
            if (is.null(r)) return(NULL)
            do.call(rbind, lapply(c("cell", "donor"), function(u) {
                d <- r[[u]]
                if (is.null(d)) return(NULL)
                data.frame(
                    dataset = ds, comparison = tg, unit = u,
                    n_fg = d$n_fg[1], n_bg = d$n_bg[1],
                    n_donors = d$n_donors[1],
                    n_tested = nrow(d),
                    n_significant = sum(d$significant, na.rm = TRUE),
                    n_up = sum(d$significant & d$avg_log2FC > 0, na.rm = TRUE),
                    n_down = sum(d$significant & d$avg_log2FC < 0, na.rm = TRUE),
                    sig_logfc = sig_thr(tg),
                    sig_padj = DEG_SIG_PADJ,
                    # Axis and rule that split the arms.
                    class_axis = "cDFG at IP_y",
                    class_rule = "Pre-IP when cDFG_expr < IP_y",
                    IP_y = if (is.finite(ip_y)) round(ip_y, 6) else NA_real_,
                    IP_x_reported = if (is.finite(ip_x)) round(ip_x, 6)
                                    else NA_real_,
                    IP_y_field = ipy$field,
                    IP_y_source = if (is.na(.p10g_iprs$path)) NA_character_
                                  else basename(.p10g_iprs$path),
                    arm_from_group_refined = as.integer(arm_src[["group_refined"]]),
                    arm_from_IP_status     = as.integer(arm_src[["IP_status"]]),
                    arm_from_cDFG_at_IP_y  = n_y_fall,
                    arm_unlabelled         = n_unlab,
                    arm_vs_cDFG_at_IP_y    = if (is.na(arm_agree)) NA_real_
                                             else round(arm_agree, 6),
                    stringsAsFactors = FALSE)
            }))
        }))

        # Calibration and drop-out diagnostics for every comparison at its own gate.
        fc_diag <- list()
        for (tg in names(res)) {
            if (is.null(res[[tg]])) next
            thr <- sig_thr(tg)
            for (u in c("cell", "donor")) {
                d <- res[[tg]][[u]]
                if (is.null(d) || nrow(d) == 0) next
                fcv <- abs(d$avg_log2FC[is.finite(d$avg_log2FC)])
                sg  <- d[!is.na(d$significant) & d$significant, , drop = FALSE]
                dn  <- sg[sg$avg_log2FC < 0, , drop = FALSE]
                # Donor-unit resolution: theoretical floor beside the observed minimum p.
                pfl <- NULL
                if (identical(u, "donor")) {
                    .nd <- suppressWarnings(as.integer(d$n_donors[1]))
                    .pv <- d$p_val[is.finite(d$p_val)]
                    .pa <- d$p_val_adj_BH[is.finite(d$p_val_adj_BH)]
                    .pmin <- if (length(.pv)) min(.pv) else NA_real_
                    pfl <- list(
                        n_donors = if (is.finite(.nd)) .nd else NA_integer_,
                        p_floor  = v10g_sr_p_floor(.nd),
                        p_min    = .pmin,
                        n_at_min = if (is.finite(.pmin))
                            sum(.pv <= .pmin * (1 + 1e-9)) else NA_integer_,
                        padj_min = if (length(.pa)) min(.pa) else NA_real_)
                }
                # Up side and group sizes, for the up/down guardrail.
                up  <- sg[sg$avg_log2FC > 0, , drop = FALSE]
                fc_diag[[paste0(tg, "|", u)]] <- list(
                    comparison = tg, unit = u, thr = thr,
                    n_tested   = nrow(d),
                    p_floor    = pfl,
                    frac_over  = if (length(fcv)) mean(fcv >= thr) else NA_real_,
                    med_absfc  = if (length(fcv)) median(fcv) else NA_real_,
                    n_sig      = nrow(sg),
                    n_up       = nrow(up), n_down = nrow(dn),
                    n_fg       = suppressWarnings(as.numeric(d$n_fg[1])),
                    n_bg       = suppressWarnings(as.numeric(d$n_bg[1])),
                    up_pct     = if (nrow(up))
                                     median(pmax(up$pct_fg, up$pct_bg))
                                 else NA_real_,
                    down_pct   = if (nrow(dn)) median(pmax(dn$pct_fg, dn$pct_bg))
                                 else NA_real_,
                    med_expr   = median(pmax(d$mean_expr_fg, d$mean_expr_bg),
                                        na.rm = TRUE))
            }
        }
        for (k in names(fc_diag)) {
            f <- fc_diag[[k]]
            cat(sprintf("    [%s %s] median |log2FC| %.4f | %.1f%% of genes over the %.3f gate | median max expr %.3f\n",
                        f$comparison, f$unit, f$med_absfc, 100 * f$frac_over,
                        f$thr, f$med_expr))
        }

        invisible(list(dataset = ds, u_IP = u_ip, n_no_donor = n_no_donor,
                       n_module = sum(is_module), bg_note = bg_note,
                       norm_verify = nv, n_bg_cells = n_bg_cells,
                       arm_src = arm_src, arm_agree = arm_agree,
                       arm_n_cmp = arm_n_cmp, n_y_fall = n_y_fall,
                       n_unlab = n_unlab, ip_y = ip_y, ip_x = ip_x,
                       ip_y_field = ipy$field, ip_y_path = .p10g_iprs$path,
                       arm_agree_u = arm_agree_u, has_cdfg = has_cdfg,
                       fc_diag = fc_diag,
                       nm_note = nm_note, nm_classes = nm_classes,
                       n_nonmod = sum(nonmod_i), nm_src = nm_src,
                       nm_unlab = nm_unlab, nm_in_mat = nm_in_mat,
                       n_nm_pre = length(nm_pre_bc),
                       n_nm_post = length(nm_post_bc),
                       nm_donors_pre = nm_donors_pre,
                       nm_donors_post = nm_donors_post,
                       nm_donors_all = nm_donors_all,
                       nm_skipped = nm_skipped,
                       norm_max = nm$max, res = res, summary = summ,
                       donor_info = lapply(res, function(r)
                           if (is.null(r)) NULL else r$donor_info)))
    }


    # =========================================================================
    # SECTION C: Execute
    # =========================================================================
    all_results <- vector("list", length(DATASETS))
    names(all_results) <- DATASETS
    for (ds in DATASETS) {
        all_results[ds] <- list(tryCatch(
            run_V10G_dataset(ds),
            error = function(e) {
                cat("  ERROR in dataset", ds, ":", conditionMessage(e), "\n")
                NULL
            }))
    }


    # =========================================================================
    # SECTION D: Workbook
    # =========================================================================
    ok_ds <- DATASETS[!vapply(all_results, is.null, logical(1L))]

    # =========================================================================
    # Two workbooks split by unit: <prefix>_DEGs_Part1_donor.xlsx (primary) and
    # <prefix>_DEGs_Part1_cell_diagnostic.xlsx (cell unit, diagnostic). Every
    # sheet carries a `unit` column.
    # =========================================================================
    cat("\n", strrep("=", 72), "\n")
    cat("  SAVING RESULTS\n")
    cat(strrep("=", 72), "\n\n")

    # One sheet accumulator per workbook; a name collision gets a numeric suffix.
    new_tabs <- function() {
        e <- new.env(parent = emptyenv())
        e$tabs <- list()
        list(
            add = function(nm, df) {
                if (is.null(df) || nrow(df) == 0) return(invisible(NULL))
                sh <- if (nchar(nm) > 31L) substr(nm, 1L, 31L) else nm
                if (sh %in% names(e$tabs)) {
                    k <- 2L
                    repeat {
                        cand <- paste0(substr(sh, 1L,
                                              31L - nchar(as.character(k)) - 1L),
                                       "_", k)
                        if (!cand %in% names(e$tabs)) break
                        k <- k + 1L
                    }
                    cat(sprintf("  sheet '%s' collided; written as '%s'\n",
                                sh, cand))
                    sh <- cand
                }
                e$tabs[[sh]] <- df
                invisible(NULL)
            },
            get = function() e$tabs)
    }

    # Sheet names fit 31 characters; the comparison is shortened so the
    # unit suffix survives.
    mk_sheet <- function(d, tg, u) {
        suf  <- paste0("_", u)
        room <- 31L - nchar(d) - 1L - nchar(suf)
        if (room < 1L) return(paste0(substr(d, 1L, 31L - nchar(suf)), suf))
        paste0(d, "_", substr(tg, 1L, room), suf)
    }

    wb_paths <- c(donor = file.path(OUT_DIR, paste0(FILE_PREFIX,
                                                    "_DEGs_Part1_donor.xlsx")),
                  cell  = file.path(OUT_DIR, paste0(FILE_PREFIX,
                                                    "_DEGs_Part1_cell_diagnostic.xlsx")))
    # Workbooks written; the cell unit stays in the RDS either way.
    wb_units <- if (isTRUE(V10G_CELL_OUTPUTS)) c("donor", "cell") else "donor"
    wb_chk <- NULL

    if (length(ok_ds) > 0) {
        summ_all <- do.call(rbind, lapply(ok_ds,
                                          function(d) all_results[[d]]$summary))
        # Testable universes (enrichment uses the donor workbook's).
        univ_all <- do.call(rbind, lapply(ok_ds, function(d)
            do.call(rbind, lapply(all_results[[d]]$res, function(r)
                if (is.null(r)) NULL else r$universe))))

        wb <- list()
        for (u in wb_units) wb[[u]] <- new_tabs()
        for (u in wb_units) {
            # The notes sheet carries a `unit` column like every other sheet.
            wb[[u]]$add("Read_Me", data.frame(
                unit = u, arm = ARM_TAG,
                inferential_unit = if (identical(u, "donor"))
                    "donor (paired pseudobulk, one difference per donor)"
                else "cell (rank-sum over cells)",
                downstream_use = if (identical(u, "donor"))
                    "primary result"
                else "diagnostic only (pseudoreplicated)",
                note = if (identical(u, "donor"))
                    "n is the number of paired donors"
                else "n is the number of cells",
                stringsAsFactors = FALSE))
            wb[[u]]$add("Summary",
                        summ_all[!is.na(summ_all$unit) & summ_all$unit == u, ,
                                 drop = FALSE])
            for (d in ok_ds) for (tg in names(all_results[[d]]$res)) {
                r <- all_results[[d]]$res[[tg]]
                if (is.null(r)) next
                wb[[u]]$add(mk_sheet(d, tg, u), r[[u]])
            }
            wb[[u]]$add("Testable_Universe",
                        univ_all[!is.na(univ_all$unit) & univ_all$unit == u, ,
                                 drop = FALSE])
        }

        # Leakage check: every sheet must carry a `unit` column equal to its workbook's unit.
        wb_chk <- lapply(wb_units, function(u) {
            tt  <- wb[[u]]$get()
            bad <- character(0)
            for (nm in names(tt)) {
                dd <- tt[[nm]]
                if (is.null(dd[["unit"]])) {
                    bad <- c(bad, paste0(nm, "(no unit column)"))
                    next
                }
                uu <- unique(as.character(dd[["unit"]]))
                uu <- uu[!is.na(uu)]
                if (!all(uu == u))
                    bad <- c(bad, sprintf("%s(%s)", nm,
                                          paste(setdiff(uu, u), collapse = "/")))
            }
            list(unit = u, n_sheets = length(tt),
                 n_rows = sum(vapply(tt, nrow, integer(1L))),
                 bad = bad, file = basename(wb_paths[[u]]))
        })
        names(wb_chk) <- wb_units

        # == Stored donor_fc matrices (comparisons sit under all_results[[d]]$res) ===
        {
            fcs <- unlist(lapply(ok_ds, function(d) {
                rr <- all_results[[d]]$res
                if (!is.list(rr) || length(rr) == 0L) return(character(0))
                out <- vapply(names(rr), function(tg) {
                    r <- rr[[tg]]
                    if (!is.list(r) || is.null(r[["donor_fc"]]))
                        return(NA_character_)
                    sprintf("%s/%s:%dx%d", d, tg, nrow(r[["donor_fc"]]),
                            ncol(r[["donor_fc"]]))
                }, character(1L))
                out[!is.na(out)]
            }))
            want <- length(ok_ds) * length(V10G_FC_CONTRASTS)
            gr10g("GR-10G-45",
                  if (length(fcs) == want) "PASS" else "WARN",
                  sprintf(paste0("donor_fc: %d matrix(es) of %d expected",
                                 " (%d dataset(s) x %d contrast(s)) | %s"),
                          length(fcs), want, length(ok_ds),
                          length(V10G_FC_CONTRASTS),
                          if (length(fcs)) paste(fcs, collapse = ", ") else
                          "none stored"))
        }

        for (d in ok_ds) {
            rp <- file.path(OUT_DIR, paste0(FILE_PREFIX, "_", d, "_results.rds"))
            saveRDS(all_results[[d]], rp)
            cat(sprintf("  [%s] RDS saved: %s (%.1f MB)\n", d, basename(rp),
                        file.size(rp) / 1e6))
        }
        for (u in wb_units) {
            tt <- wb[[u]]$get()
            if (length(tt) == 0L) {
                cat(sprintf("  %s workbook: nothing to write\n", u))
                next
            }
            tryCatch({
                write_xlsx(tt, path = wb_paths[[u]])
                cat(sprintf("  Excel saved: %-46s (%d sheet(s), %s)\n",
                            basename(wb_paths[[u]]), length(tt),
                            if (identical(u, "donor")) "primary"
                            else "diagnostic only"))
            }, error = function(e)
                cat("  Excel error:", conditionMessage(e), "\n"))
        }
    } else {
        cat("  No dataset completed; nothing to write.\n")
    }


    # =========================================================================
    # SECTION E: Post-analysis guardrails
    # =========================================================================
    cat("\n  POST-ANALYSIS GUARDRAILS\n\n")
    for (ds in DATASETS) {
        r <- all_results[[ds]]
        if (is.null(r)) {
            gr10g(paste0("GR-10G-04-", ds), "FAIL", "Dataset returned NULL")
            next
        }
        gr10g(paste0("GR-10G-04-", ds), "PASS",
              sprintf("u_IP = %.6f | %d module cell(s) | %d cell(s) without a donor",
                      r$u_IP, r$n_module, r$n_no_donor))

        # the normalised matrix must be on the CP10K scale, not raw counts
        gr10g(paste0("GR-10G-05-", ds),
              if (is.finite(r$norm_max) && r$norm_max > 0 &&
                  r$norm_max <= log1p(1e4) + 1e-4) "PASS" else "FAIL",
              sprintf("norm.mtx max = %.4f (CP10K ceiling log1p(10000) = %.4f)",
                      r$norm_max, log1p(1e4)))

        h2h <- r$res[["PostIP_vs_PreIP"]]
        # The donor unit is primary and must have run; lookups use exact [[ ]].
        di      <- if (is.null(h2h)) NULL else h2h[["donor_info"]]
        h2h_don <- if (is.null(h2h)) NULL else h2h[["donor"]]
        h2h_cel <- if (is.null(h2h)) NULL else h2h[["cell"]]
        gr10g(paste0("GR-10G-06-", ds),
              if (!is.null(h2h_don)) "PASS" else "WARN",
              if (is.null(h2h)) "head-to-head comparison did not run"
              else if (!is.null(h2h_don))
                  sprintf("donor unit ran on %d paired donor(s) of %d (>= %d cell(s) per arm); %d fg-only, %d bg-only",
                          h2h_don$n_donors[1], di$n_donors_total,
                          di$min_cells, di$n_donors_fg_only, di$n_donors_bg_only)
              else sprintf("donor unit did not run: %s; cell unit only",
                           if (is.null(di) || is.null(di$note))
                               "reason not recorded" else di$note))

        # Significant counts at both units.
        if (!is.null(h2h_cel) && !is.null(h2h_don)) {
            nc <- sum(h2h_cel$significant, na.rm = TRUE)
            nd <- sum(h2h_don$significant, na.rm = TRUE)
            gr10g(paste0("GR-10G-07-", ds), "PASS",
                  sprintf(paste0("significant genes: %d at the cell unit, %d at",
                                 " the donor unit%s"),
                          nc, nd,
                          if (nd > 0) sprintf(" (%.1fx)", nc / nd)
                          else " (none survive the donor unit)"))
        }

        # Donor unit present on every comparison.
        cmp_names <- names(r$res)
        du <- vapply(cmp_names, function(tg) {
            x <- r$res[[tg]]
            !is.null(x) && !is.null(x[["donor"]])
        }, logical(1L))
        miss <- cmp_names[!du]
        gr10g(paste0("GR-10G-15-", ds),
              if (length(cmp_names) > 0L && all(du)) "PASS" else "WARN",
              if (length(cmp_names) == 0L) "no comparison ran"
              else sprintf(paste0("donor unit present on %d of %d comparison(s):",
                                  " %s%s"),
                           sum(du), length(du), paste(cmp_names[du], collapse = ", "),
                           if (length(miss))
                               sprintf(" | missing on %s (cell unit only)",
                                       paste(miss, collapse = ", "))
                           else ""))

        # The five non-module contrasts: FAIL if none ran or coverage is short; WARN if
        # the donor unit is missing on any.
        nm_ran <- intersect(c(V10G_NONMOD_CMP, V10G_NONMOD_CMP2), names(r$res))
        nm_dn  <- vapply(nm_ran, function(tg)
            !is.null(r$res[[tg]]) && !is.null(r$res[[tg]][["donor"]]),
            logical(1L))
        nm_cov_ok <- isTRUE(r$nm_in_mat >= r$n_nm_pre + r$n_nm_post)
        nm_status <- if (!isTRUE(V10G_NONMOD_ENABLE)) "WARN"
                     else if (isTRUE(r$n_nonmod == 0L)) "WARN"
                     else if (length(nm_ran) == 0L || !nm_cov_ok) "FAIL"
                     else if (!all(nm_dn)) "WARN"
                     else "PASS"
        gr10g(paste0("GR-10G-17-", ds), nm_status,
              sprintf(paste0("non-module %d cell(s) in %s | Pre-IP %d /",
                             " Post-IP %d | %d of %d in norm.mtx | %d of 5",
                             " contrast(s) ran (%s)%s | donors with >= %d",
                             " cell(s): Pre-IP %d, Post-IP %d, pooled %s",
                             " against the %d needed to pair%s"),
                      r$n_nonmod,
                      if (length(r$nm_classes))
                          paste(r$nm_classes, collapse = " + ")
                      else "no class",
                      r$n_nm_pre, r$n_nm_post, r$nm_in_mat,
                      r$n_nm_pre + r$n_nm_post,
                      length(nm_ran),
                      if (length(nm_ran)) paste(nm_ran, collapse = ", ")
                      else "none",
                      if (length(r$nm_skipped))
                          sprintf(" | skipped below the cell minimum: %s",
                                  paste(r$nm_skipped, collapse = ", "))
                      else "",
                      V10G_MIN_CELLS_PER_DONOR, r$nm_donors_pre,
                      r$nm_donors_post,
                      if (is.null(r$nm_donors_all)) "not recorded"
                      else format(r$nm_donors_all),
                      V10G_MIN_DONORS,
                      if (length(nm_ran) && !all(nm_dn))
                          sprintf(paste0(" | donor unit missing on %s",
                                         " (cell unit only)"),
                                  paste(nm_ran[!nm_dn], collapse = ", "))
                      else ""))

        # IP_y must be present; IP_x is reported only.
        gr10g(paste0("GR-10G-14-", ds),
              if (isTRUE(is.finite(r$ip_y))) "PASS" else "FAIL",
              if (isTRUE(is.finite(r$ip_y)))
                  sprintf(paste0("IP_y = %.6f from `%s` in %s | IP_x = %s",
                                 " (reported, GOI axis, not the classifier)"),
                          r$ip_y, r$ip_y_field, basename(r$ip_y_path),
                          if (is.finite(r$ip_x)) sprintf("%.6f", r$ip_x)
                          else "unavailable")
              else sprintf("IP_y not available (%s)",
                           if (is.na(r$ip_y_path)) "no IPRS RDS found"
                           else paste0("no inflection_points_DFG_y in ",
                                       basename(r$ip_y_path))))

        # Arm provenance and agreement of the inherited label with cDFG < IP_y:
        # FAIL below 0.999, WARN if unchecked or not all from group_refined.
        asrc <- r$arm_src
        .agree_bad  <- isTRUE(!is.na(r$arm_agree) && r$arm_agree < 0.999)
        .no_check   <- !isTRUE(r$has_cdfg) || is.na(r$arm_agree)
        .degraded   <- isTRUE(r$n_y_fall > 0L) || isTRUE(r$n_unlab > 0L) ||
                       !isTRUE(as.integer(asrc[["group_refined"]]) > 0)
        gr10g(paste0("GR-10G-09-", ds),
              if (.agree_bad) "FAIL"
              else if (.no_check || .degraded) "WARN" else "PASS",
              sprintf(paste0("arm inherited: group_refined %d, IP_status %d,",
                             " cDFG@IP_y fallback %d, unlabelled %d | inherited",
                             " vs `cDFG < IP_y`: %s on %d cell(s) | GOI-axis",
                             " split at u_IP would agree %s%s%s"),
                      as.integer(asrc[["group_refined"]]),
                      as.integer(asrc[["IP_status"]]), r$n_y_fall, r$n_unlab,
                      if (is.na(r$arm_agree)) "not checked"
                      else sprintf("%.2f%%", 100 * r$arm_agree),
                      r$arm_n_cmp,
                      if (is.na(r$arm_agree_u)) "NA"
                      else sprintf("%.2f%%", 100 * r$arm_agree_u),
                      if (.agree_bad)
                          " -- inherited labels disagree with cDFG at this IP_y"
                      else "",
                      if (!isTRUE(r$has_cdfg))
                          " -- cDFG_expr or IP_y unavailable; cross-check not run"
                      else if (.degraded)
                          " -- not every module cell took its arm from group_refined"
                      else ""))

        # Fold-change gate calibration: 2-40% of tested genes should clear it.
        for (k in names(r$fc_diag)) {
            f <- r$fc_diag[[k]]
            gr10g(paste0("GR-10G-10-", ds, "-", f$comparison, "-", f$unit),
                  if (is.finite(f$frac_over) &&
                      f$frac_over >= 0.02 && f$frac_over <= 0.40) "PASS" else "WARN",
                  sprintf(paste0("[%s %s] %.2f%% of %d tested gene(s) clear",
                                 " |log2FC| >= %.3f (median |log2FC| %.4f,",
                                 " median max expr %.3f; expected 2-40%%)"),
                          f$comparison, f$unit, 100 * f$frac_over, f$n_tested,
                          f$thr, f$med_absfc, f$med_expr))
        }

        # Drop-out artefact risk: skew toward the larger group (the down-only rule when
        # the groups are within GR11_SIZE_RATIO), with the detection rate.
        for (k in names(r$fc_diag)) {
            f <- r$fc_diag[[k]]
            # Below GR11_MIN_SIG significant genes the ratio is not assessed.
            if (f$n_sig < GR11_MIN_SIG) {
                gr10g(paste0("GR-10G-11-", ds, "-", f$comparison, "-", f$unit),
                      "PASS",
                      sprintf(paste0("[%s %s] %d significant: %d up / %d down.",
                                     " %s"),
                              f$comparison, f$unit, f$n_sig, f$n_up, f$n_down,
                              if (f$n_sig == 0L)
                                  "No significant gene, so there is no skew"
                              else sprintf(paste("Fewer than %d significant",
                                                 "gene(s), so the skew ratio is",
                                                 "not assessable"),
                                           GR11_MIN_SIG)))
                next
            }
            # Group sizes as length-1 values.
            .nf <- if (length(f$n_fg)) f$n_fg[1] else NA_real_
            .nb <- if (length(f$n_bg)) f$n_bg[1] else NA_real_
            .sized <- isTRUE(is.finite(.nf)) && isTRUE(is.finite(.nb)) &&
                      .nf > 0 && .nb > 0 &&
                      max(.nf, .nb) / min(.nf, .nb) >= GR11_SIZE_RATIO
            .toward_fg <- .sized && .nf > .nb
            if (!.sized) {
                .maj <- f$n_down; .min <- f$n_up
                .maj_pct <- f$down_pct; .lab <- "down"
                .why <- sprintf(paste("groups within %.1fx (%s against %s),",
                                      "so the original down-only rule",
                                      "applies"),
                                GR11_SIZE_RATIO,
                                if (isTRUE(is.finite(.nf))) format(.nf)
                                else "NA",
                                if (isTRUE(is.finite(.nb))) format(.nb)
                                else "NA")
            } else if (.toward_fg) {
                .maj <- f$n_up; .min <- f$n_down
                .maj_pct <- f$up_pct; .lab <- "up"
                .why <- sprintf(paste("the foreground is the larger group",
                                      "(%s against %s), so the artefact",
                                      "direction is up"),
                                format(.nf), format(.nb))
            } else {
                .maj <- f$n_down; .min <- f$n_up
                .maj_pct <- f$down_pct; .lab <- "down"
                .why <- sprintf(paste("the background is the larger group",
                                      "(%s against %s), so the artefact",
                                      "direction is down"),
                                format(.nb), format(.nf))
            }
            skew <- if (.min > 0) .maj / .min else Inf
            # An empty artefact side has no median detection.
            .det <- if (.maj == 0L)
                        sprintf(paste("no gene(s) on the %s side, so no",
                                      "detection rate to report"), .lab)
                    else if (is.na(.maj_pct))
                        sprintf(paste("those %s genes carry no detection",
                                      "rate"), .lab)
                    else
                        sprintf(paste("those %s genes are detected in %.1f%%",
                                      "of cells at best (median)"),
                                .lab, 100 * .maj_pct)
            gr10g(paste0("GR-10G-11-", ds, "-", f$comparison, "-", f$unit),
                  if (skew <= 3 && (is.na(.maj_pct) || .maj_pct >= 0.10))
                      "PASS" else "WARN",
                  sprintf(paste0("[%s %s] %d significant: %d up / %d down;",
                                 " %s. Skew toward %s = %.2f; %s"),
                          f$comparison, f$unit, f$n_sig, f$n_up, f$n_down,
                          .why, .lab, skew, .det))
        }

        # Donor-unit resolution floor: smallest attainable signed-rank p; WARN when the
        # observed minimum BH is above alpha (resolution-limited).
        for (k in names(r$fc_diag)) {
            f  <- r$fc_diag[[k]]
            if (!identical(f$unit, "donor")) next
            pf <- f$p_floor
            # A floor that cannot be computed still gets a row.
            if (is.null(pf) || !isTRUE(is.finite(pf$p_floor))) {
                gr10g(paste0("GR-10G-12-", ds, "-", f$comparison), "WARN",
                      sprintf(paste0("[%s donor] the signed-rank resolution",
                                     " floor could not be computed"),
                              f$comparison))
                next
            }
            padj_one    <- pf$p_floor * f$n_tested
            g_need      <- ceiling(padj_one / DEG_SIG_PADJ)
            floor_tight <- isTRUE(padj_one > DEG_SIG_PADJ)
            limited     <- if (isTRUE(is.finite(pf$padj_min)))
                               isTRUE(pf$padj_min > DEG_SIG_PADJ) else floor_tight
            tail_txt <-
                if (limited && isTRUE(is.finite(pf$padj_min)))
                    sprintf(paste0(" Resolution-limited: no gene reached",
                                   " BH < %.2f (observed minimum %.2e); one gene",
                                   " at the floor reaches only BH = %.2e",
                                   " (%d tied gene(s) needed)."),
                            DEG_SIG_PADJ, pf$padj_min, padj_one, g_need)
                else if (limited)
                    sprintf(paste0(" Resolution-limited: observed minimum BH",
                                   " unavailable; one gene at the floor reaches",
                                   " only BH = %.2e."),
                            padj_one)
                else if (floor_tight)
                    sprintf(paste0(" Single-gene BH at the floor %.2e is above",
                                   " alpha = %.2f; observed minimum %.2e; not",
                                   " resolution-limited."),
                            padj_one, DEG_SIG_PADJ, pf$padj_min)
                else ""
            gr10g(paste0("GR-10G-12-", ds, "-", f$comparison),
                  if (limited) "WARN" else "PASS",
                  sprintf(paste0("[%s donor] %d paired donor(s): smallest",
                                 " attainable signed-rank p = %.2e; one gene",
                                 " alone reaches BH = %.2e over %d gene(s).",
                                 " Observed: %s gene(s) tie at p = %s, reaching",
                                 " BH = %s.%s"),
                          f$comparison, pf$n_donors, pf$p_floor, padj_one,
                          f$n_tested,
                          if (is.finite(pf$n_at_min)) format(pf$n_at_min) else "NA",
                          if (is.finite(pf$p_min)) sprintf("%.2e", pf$p_min) else "NA",
                          if (is.finite(pf$padj_min)) sprintf("%.2e", pf$padj_min)
                          else "NA",
                          tail_txt))
        }

        # CP10K recipe check: WARN if raw_selected is absent, FAIL if outside tolerance.
        nv <- r$norm_verify
        gr10g(paste0("GR-10G-13-", ds),
              if (isTRUE(nv$ok)) "PASS"
              else if (!is.finite(nv$dev)) "WARN" else "FAIL",
              if (is.finite(nv$dev))
                  sprintf(paste0("CP10K recomputed from raw_selected vs",
                                 " deg_support/norm.mtx: max |dev| = %.3e on %d",
                                 " gene(s) x %d cell(s), tolerance %.1e"),
                          nv$dev, nv$n_genes, nv$n_cells, V10G_NORM_TOL)
              else paste0("recipe not verified: ", nv$note))

        # Arm-vs-background comparisons: did they run, and on what.
        bg_ran <- !is.null(r$res[["PostIP_vs_BG"]]) ||
                  !is.null(r$res[["PreIP_vs_BG"]])
        gr10g(paste0("GR-10G-08-", ds),
              if (bg_ran) "PASS" else "WARN",
              if (is.na(r$bg_note)) "no background note recorded" else r$bg_note)
    }

    # Leakage check, arm level: a wrong-unit row or a sheet with no `unit` column FAILs.
    if (is.null(wb_chk)) {
        gr10g("GR-10G-16", "WARN",
              "no workbook was written, so the unit split could not be checked")
    } else {
        for (u in c("donor", "cell")) {
            # A suppressed cell workbook reports rather than warns.
            if (!u %in% wb_units) {
                # Cell-unit files left from an earlier run are flagged.
                .stale <- list.files(OUT_DIR, pattern = "_cell_")
                gr10g(paste0("GR-10G-16-", u),
                      if (length(.stale) == 0L) "PASS" else "WARN",
                      sprintf(paste0("%s: not written (V10G_CELL_OUTPUTS = FALSE);",
                                     " the %s unit is in %s_<ds>_results.rds.%s"),
                              basename(wb_paths[[u]]), u, FILE_PREFIX,
                              if (length(.stale) == 0L) ""
                              else sprintf(paste0(" Stale: %d file(s) matching",
                                                  " \"_cell_\" from an earlier run: %s"),
                                           length(.stale),
                                           paste(head(sort(.stale), 6L),
                                                 collapse = ", "))))
                next
            }
            w <- wb_chk[[u]]
            gr10g(paste0("GR-10G-16-", u),
                  if (length(w$bad) == 0L) "PASS" else "FAIL",
                  if (length(w$bad) == 0L)
                      sprintf(paste0("%s: %d sheet(s), %d row(s), every sheet",
                                     " carries unit = \"%s\" and nothing else%s"),
                              w$file, w$n_sheets, w$n_rows, u,
                              if (identical(u, "donor"))
                                  " -- primary"
                              else " -- diagnostic only")
                  else
                      sprintf(paste0("%s: %d sheet(s) carry a unit other than",
                                     " \"%s\" or no unit column at all: %s"),
                              w$file, length(w$bad), u,
                              paste(w$bad, collapse = ", ")))
        }
    }

    gr10g_write()
    cat(sprintf("\n  Final guardrail tally: %d PASS / %d WARN / %d FAIL\n",
                .gr_pass, .gr_warn, .gr_fail))

    cat("\n", strrep("=", 72), "\n")
    cat("  Part10g DEGs Part 1 complete [", ARM_TAG, "]\n")
    cat("  Output directory:", OUT_DIR, "\n")
    cat("  Files generated:\n")
    .tr <- c("_RUN_INCOMPLETE.txt", basename(MSG_PATH_10G))
    for (f in setdiff(sort(list.files(OUT_DIR)), .tr)) {
        sz <- round(file.size(file.path(OUT_DIR, f)) / 1e3, 0)
        cat(sprintf("    %-62s  %5d KB\n", f, sz))
    }
    cat(strrep("=", 72), "\n")
    cat("\n  Log written to:", LOG_PATH_10G, "\n")
    cat(sprintf("\nPart10g DEGs Part 1 - log ended: %s\n", Sys.time()))

    .p10g_completed <- TRUE
    if (sink.number(type = "message") != 2L) sink(type = "message")
    while (sink.number() > 0) sink()
    suppressWarnings(file.remove(INCOMPLETE_MARKER_10G))

    invisible(list(arm = arm_tag, dir = OUT_DIR, pass = .gr_pass,
                   warn = .gr_warn, fail = .gr_fail))
}


# =============================================================================
# SECTION 0c: Driver
# =============================================================================

cep_run_arms(PART_DIR_TAG, run_one_arm)
