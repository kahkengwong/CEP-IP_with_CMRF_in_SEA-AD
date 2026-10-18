# =============================================================================
# Part10g_CP10K_DEGs_Part1_code_S11_followup_analyses.r
# SEA-AD CEP-IP/CMRF, Part 10g, follow-up DE analyses (Table S11)
#   F1 control-test survivors in both region sets; F2 phagocytic, lysosomal
#   and homeostatic genes; F3 short neuron-specific genes; F4 UMIs per nucleus
#   Post-IP against Pre-IP within donor; F5 panel UMIs per nucleus (amount and
#   share) on the F4 cells.
#
# Usage (from the SEA-AD folder, or set SEAAD_ROOT), after Part 10g Part 1,
#   Part12_04_Labels.r, with Table S11 in Manuscript/Supp_Tables/ (F1 reads
#   its sheet B)
#   source("Part10g_CP10K_DEGs_Part1_code_S11_followup_analyses.r")
#   Optional: FU_F4_IPY, a named numeric IP_y per arm and region set
#   ("M4.NeoC", "PVM_3.NonNeoC", ...), checked against the Module gap.
#
# Outputs (SuppTable_S11_build/): S11_F1 to S11_F5 CSV tables and log.
# =============================================================================

if (!exists("CEP_HOME")) CEP_HOME <- Sys.getenv("SEAAD_ROOT", getwd())
stopifnot(dir.exists(CEP_HOME))
if (!requireNamespace("readxl", quietly = TRUE))
    stop("Package readxl is required.", call. = FALSE)

FU_OUT  <- file.path(CEP_HOME, "SuppTable_S11_build")
dir.create(FU_OUT, showWarnings = FALSE)
FU_LOG  <- file.path(FU_OUT, paste0("S11_followup_log_",
                                    format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
S11_XLSX <- file.path(CEP_HOME, "Manuscript", "Supp_Tables",
                      "Supplementary_Table_S11_DEGs_ModvsNonMod.xlsx")

# == Shared settings and helpers =============================================
FU_ARMS     <- c(M4 = "Mic13_PVM3", PVM_3 = "PVM_3", PVM_1 = "PVM_1",
                 HKG8 = "HKG8")
FU_MAIN     <- c("M4", "PVM_3", "PVM_1")
FU_DATASETS <- c("NeoC", "NonNeoC")

fu_hdr <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n",
                          strrep("=", 78), "\n", sep = "")
fu_lgl <- function(x) toupper(trimws(as.character(x))) %in% "TRUE"

# Read the eight Part 1 RDS once.
fu_read_rds <- function() {
    R <- list()
    for (a in names(FU_ARMS)) for (ds in FU_DATASETS) {
        d1 <- file.path(CEP_HOME, "Part10g_DEGs",
                        paste0("Part10g_DEGs_", FU_ARMS[[a]]))
        f  <- list.files(d1, pattern = paste0("_", ds, "_results\\.rds$"),
                         full.names = TRUE)
        if (length(f) != 1L) stop(sprintf("%s %s: %d RDS found in %s",
                                          a, ds, length(f), d1), call. = FALSE)
        R[[paste(a, ds)]] <- readRDS(f)
        cat(sprintf("  %-13s %s\n", paste(a, ds), basename(f)))
    }
    R
}

# Long table: genes x arms x region sets x contrasts, from the donor frames.
fu_lookup <- function(R, genes, contrasts) {
    out <- list()
    for (a in names(FU_ARMS)) for (ds in FU_DATASETS) {
        o <- R[[paste(a, ds)]]
        for (tg in contrasts) {
            d <- o$res[[tg]]$donor
            if (!is.data.frame(d))
                stop(sprintf("%s %s %s: no donor frame", a, ds, tg), call. = FALSE)
            k <- match(genes, d$gene)
            out[[length(out) + 1L]] <- data.frame(
                gene = genes, arm = a, region_set = ds, contrast = tg,
                n_donors = as.numeric(d$n_donors[1]),
                log2FC = round(as.numeric(d$avg_log2FC[k]), 3),
                BH = signif(as.numeric(d$p_val_adj_BH[k]), 3),
                significant = fu_lgl(d$significant[k]),
                stringsAsFactors = FALSE)
        }
    }
    do.call(rbind, out)
}

# Wide print: one row per gene, one column per arm; * = significant.
fu_print_wide <- function(L, genes, tg, ds) {
    s <- L[L$contrast == tg & L$region_set == ds, ]
    w <- data.frame(gene = genes, stringsAsFactors = FALSE)
    for (a in names(FU_ARMS)) {
        sa <- s[s$arm == a, ]; k <- match(genes, sa$gene)
        w[[a]] <- ifelse(is.na(sa$log2FC[k]), "NA",
                         paste0(formatC(sa$log2FC[k], format = "f", digits = 3),
                                ifelse(sa$significant[k], "*", "")))
    }
    nd <- vapply(names(FU_ARMS), function(a) s$n_donors[s$arm == a][1], numeric(1))
    cat(sprintf("\n-- %s, %s (paired donors: %s) --\n", tg, ds,
                paste(sprintf("%s %d", names(FU_ARMS), nd), collapse = ", ")))
    print(w, row.names = FALSE)
}

# =============================================================================
# F1. Control-test survivors in both region sets
# =============================================================================
fu_F1 <- function(R) {
    fu_hdr("F1. NON-MITOCHONDRIAL POST-IP vs PRE-IP SURVIVORS, BOTH REGION SETS")
    if (!file.exists(S11_XLSX))
        stop("Table S11 workbook not found: ", S11_XLSX, call. = FALSE)
    B <- as.data.frame(readxl::read_excel(S11_XLSX, sheet = "B_Control_test"))
    sv <- B[B$comparison == "PostIP_vs_PreIP" & fu_lgl(B$passes_control_test) &
            !fu_lgl(B$mitochondrial),
            c("gene", "region_set", "direction_on_PVM_3")]
    names(sv) <- c("gene", "passed_in", "direction")
    if (anyDuplicated(sv$gene))
        stop("a survivor passed in both region sets; F1 assumes one",
             call. = FALSE)
    cat(sprintf("  %d survivors read from sheet B (%s)\n", nrow(sv),
                paste(sprintf("%s %d", names(table(sv$passed_in)),
                              table(sv$passed_in)), collapse = ", ")))
    genes <- sv$gene

    L <- fu_lookup(R, genes, c("PostIP_vs_PreIP", "NMPost_vs_NMPre"))
    for (tg in c("PostIP_vs_PreIP", "NMPost_vs_NMPre"))
        for (ds in FU_DATASETS) fu_print_wide(L, genes, tg, ds)

    # Same direction in the other region set (Module Post-IP vs Pre-IP)
    M <- L[L$contrast == "PostIP_vs_PreIP", ]
    S <- sv
    S$other_region <- ifelse(S$passed_in == "NeoC", "NonNeoC", "NeoC")
    sg <- ifelse(S$direction == "up", 1, -1)
    for (a in names(FU_ARMS)) for (side in c("passed", "other")) {
        rs <- if (side == "passed") S$passed_in else S$other_region
        v  <- M$log2FC[match(paste(genes, a, rs),
                             paste(M$gene, M$arm, M$region_set))]
        S[[paste0("log2FC_", a, "_", side)]] <- v
    }
    same <- sapply(FU_MAIN, function(a)
        sign(S[[paste0("log2FC_", a, "_other")]]) == sg)
    same[is.na(same)] <- FALSE
    S$main_arms_same_direction_other <- rowSums(same)
    # Smallest |log2FC| over the main arms in the other region set, so a
    # same-sign value near zero is visible.
    S$min_abs_log2FC_main_other <- apply(abs(as.matrix(
        S[, paste0("log2FC_", FU_MAIN, "_other")])), 1, min)
    S$all_main_arms_same_direction_other <- S$main_arms_same_direction_other == 3L
    S$HKG8_opposite_other <- sign(S$log2FC_HKG8_other) == -sg
    S$HKG8_opposite_other[is.na(S$HKG8_opposite_other)] <- FALSE
    S$sig_other_any_main_arm <- sapply(seq_len(nrow(S)), function(i)
        any(M$significant[M$gene == genes[i] & M$region_set == S$other_region[i] &
                          M$arm %in% FU_MAIN]))

    cat("\n  Summary (Module Post-IP vs Pre-IP, region set where the gene did not pass):\n")
    print(S[, c("gene", "passed_in", "direction",
                paste0("log2FC_", FU_MAIN, "_other"), "log2FC_HKG8_other",
                "main_arms_same_direction_other", "min_abs_log2FC_main_other",
                "HKG8_opposite_other",
                "sig_other_any_main_arm")], row.names = FALSE)
    cat(sprintf(paste0("\n  Same direction on all three main arms in the other",
                       " region set: %d of %d (%s)\n"),
                sum(S$all_main_arms_same_direction_other), nrow(S),
                paste(S$gene[S$all_main_arms_same_direction_other], collapse = ", ")))
    cat(sprintf("  Not on all three: %s\n",
                paste(S$gene[!S$all_main_arms_same_direction_other], collapse = ", ")))
    cat(sprintf("  Gene x main-arm values in the same direction: %d of %d\n",
                sum(same), length(same)))
    cat(sprintf("  HKG8 opposite in the other region set: %d of %d\n",
                sum(S$HKG8_opposite_other), nrow(S)))
    cat(sprintf("  Significant on any main arm in the other region set: %d of %d\n",
                sum(S$sig_other_any_main_arm), nrow(S)))

    write.csv(L, file.path(FU_OUT, "S11_F1_survivors_both_regions.csv"),
              row.names = FALSE)
    write.csv(S, file.path(FU_OUT, "S11_F1_survivors_summary.csv"),
              row.names = FALSE)
    invisible(list(long = L, summary = S))
}

# =============================================================================
# F2. Phagocytic / lysosomal and homeostatic genes
# =============================================================================
FU_PHAGO <- c("CD68", "LAMP1", "LAMP2", "CTSB", "CTSD", "CTSS", "CTSZ", "TREM2",
              "TYROBP", "AXL", "MERTK", "CD36", "LPL", "GPNMB", "ITGAX",
              "CLEC7A", "MSR1")
FU_HOMEO <- c("P2RY12", "CX3CR1", "TMEM119")

fu_F2 <- function(R) {
    fu_hdr("F2. PHAGOCYTIC / LYSOSOMAL AND HOMEOSTATIC GENES")
    genes <- c(FU_PHAGO, FU_HOMEO)
    cat("  Phagocytic / lysosomal:", paste(FU_PHAGO, collapse = ", "), "\n")
    cat("  Homeostatic           :", paste(FU_HOMEO, collapse = ", "), "\n")
    L <- fu_lookup(R, genes, c("PostIP_vs_PreIP", "PostIP_vs_NMPost"))
    for (tg in c("PostIP_vs_PreIP", "PostIP_vs_NMPost"))
        for (ds in FU_DATASETS) fu_print_wide(L, genes, tg, ds)
    n_na <- sum(is.na(L$log2FC))
    cat(sprintf("\n  Gene x arm x region x contrast values not in a frame (NA): %d of %d\n",
                n_na, nrow(L)))
    write.csv(L, file.path(FU_OUT, "S11_F2_phago_lookup.csv"), row.names = FALSE)
    invisible(L)
}

# =============================================================================
# F3. Short neuron-specific genes against the long neuronal survivors
# =============================================================================
FU_SHORT_NEURO <- c("SNAP25", "NRGN", "SLC17A7", "GAD1", "GAD2", "STMN2",
                    "SYP", "VSNL1", "NEFL", "RBFOX3", "ENC1")

fu_F3 <- function(R, long_genes) {
    fu_hdr("F3. SHORT NEURON-SPECIFIC GENES vs THE LONG NEURONAL SURVIVORS")
    cat("  Short neuron-specific:", paste(FU_SHORT_NEURO, collapse = ", "), "\n")
    cat("  Long survivors (F1)  :", paste(long_genes, collapse = ", "), "\n")
    tg2 <- c("PostIP_vs_PreIP", "NMPost_vs_NMPre")
    L <- fu_lookup(R, c(FU_SHORT_NEURO, long_genes), tg2)
    L$panel <- ifelse(L$gene %in% FU_SHORT_NEURO, "short_neuronal",
                      "long_survivor")
    for (tg in tg2) for (ds in FU_DATASETS)
        fu_print_wide(L, FU_SHORT_NEURO, tg, ds)

    # Per arm, region set and contrast: how each panel moves.
    S <- do.call(rbind, lapply(split(L, list(L$panel, L$arm, L$region_set,
                                             L$contrast), drop = TRUE),
        function(s) data.frame(
            panel = s$panel[1], arm = s$arm[1], region_set = s$region_set[1],
            contrast = s$contrast[1], n_genes = nrow(s),
            n_detected = sum(!is.na(s$log2FC)),
            median_log2FC = if (all(is.na(s$log2FC))) NA_real_ else
                round(stats::median(s$log2FC, na.rm = TRUE), 3),
            n_up = sum(s$log2FC > 0, na.rm = TRUE),
            n_sig_up = sum(s$significant & s$log2FC > 0, na.rm = TRUE),
            n_sig_down = sum(s$significant & s$log2FC < 0, na.rm = TRUE),
            stringsAsFactors = FALSE)))
    S$arm <- factor(S$arm, levels = names(FU_ARMS))
    S <- S[order(match(S$contrast, tg2), S$region_set, S$arm, S$panel), ]
    S$arm <- as.character(S$arm); rownames(S) <- NULL
    cat("\n  Panel summary (median_log2FC over detected genes):\n")
    print(S, row.names = FALSE)
    write.csv(L, file.path(FU_OUT, "S11_F3_short_neuronal_lookup.csv"),
              row.names = FALSE)
    write.csv(S, file.path(FU_OUT, "S11_F3_short_neuronal_summary.csv"),
              row.names = FALSE)
    invisible(list(long = L, summary = S))
}

# =============================================================================
# F4. Total UMIs per nucleus, Post-IP against Pre-IP, within donor
# =============================================================================
# UMIs per nucleus from the Part 8 raw counts (total and non-MT), joined to
# the labels by barcode. Column candidates; FU_F4_COLS_OVERRIDE adds names.
FU_F4_COLS <- list(
    barcode    = c("barcode", "cell_id", "cell", "obs_names", "index"),
    arm        = c("arm", "arm_label", "ARM_TAG", "arm_tag"),
    region_set = c("dataset", "region_set", "ds", "region_group"),
    donor      = c("donor", "donor_id", "Donor ID", "donor_name"),
    module     = c("is_module", "module", "in_module", "module_member"),
    ip_side    = c("ip_status", "ip_side", "IP_side", "IP_status", "side",
                   "post_ip", "is_post_ip"),
    composite  = c("composite", "composite_score", "arm_composite", "score"))
if (exists("FU_F4_COLS_OVERRIDE") && is.list(FU_F4_COLS_OVERRIDE))
    for (k in names(FU_F4_COLS_OVERRIDE))
        FU_F4_COLS[[k]] <- unique(c(FU_F4_COLS_OVERRIDE[[k]], FU_F4_COLS[[k]]))
FU_F4_MIN_CELLS <- 3L     # per donor and group, as Part 1's pairing rule
FU_F4_ARM_TAGS  <- c(M4 = "Mic13_PVM3", PVM_3 = "PVM3_all",
                     PVM_1 = "PVM1_all", HKG8 = "HKG8_matched")
# F5 calibration panel: long microglial/PVM genes, in no arm's composite.
FU_MG_OWN <- c("CSF1R", "P2RY12", "CX3CR1", "DOCK8", "ARHGAP15", "INPP5D")

fu_pick <- function(nm, cands) {
    hit <- nm[match(tolower(cands), tolower(nm))]; hit <- hit[!is.na(hit)]
    if (length(hit)) return(hit[1])
    # Fallback: a name that ends with the candidate (e.g. after a BOM).
    for (cc in tolower(cands)) {
        k <- which(endsWith(tolower(nm), cc))
        if (length(k)) return(nm[k[1]])
    }
    NA_character_
}
fu_arm_std <- function(x) {
    u <- toupper(gsub("[^A-Za-z0-9]", "", as.character(x)))
    out <- rep(NA_character_, length(u))
    out[u %in% c("M4", "MIC13PVM3", "MIC13")] <- "M4"
    out[u %in% c("PVM3", "PVM3ALL", "PVM3BULK")] <- "PVM_3"
    out[u %in% c("PVM1", "PVM1ALL", "PVM1BULK")] <- "PVM_1"
    out[u %in% c("HKG8", "HKG8MATCHED")] <- "HKG8"
    out
}
fu_rs_std <- function(x) {
    u <- toupper(gsub("[^A-Za-z]", "", as.character(x)))
    ifelse(u %in% c("NEOC", "NEOCORTEX"), "NeoC",
    ifelse(u %in% c("NONNEOC", "NONNEOCORTEX"), "NonNeoC", NA_character_))
}
fu_lgl_map <- function(v, true_re, false_re) {
    if (is.logical(v)) return(v)
    u <- tolower(gsub("[^A-Za-z0-9]", "", as.character(v)))
    ifelse(u %in% c("true", "1") | grepl(true_re, u), TRUE,
    ifelse(u %in% c("false", "0") | grepl(false_re, u), FALSE, NA))
}

# UMIs per barcode from one raw_selected folder: total, non-mitochondrial
# and, for each gene panel in `panels` (named list), the panel's UMIs.
fu_umi_from_raw <- function(d, panels = list()) {
    mtx <- list.files(d, pattern = "^matrix\\.mtx(\\.gz)?$", full.names = TRUE)
    bcf <- list.files(d, pattern = "barcode", full.names = TRUE)
    ftf <- list.files(d, pattern = "feature|gene", full.names = TRUE)
    if (length(mtx) != 1L || length(bcf) != 1L || length(ftf) != 1L)
        stop(sprintf("F4: %s needs one matrix.mtx, one barcodes file and one features file (found %d, %d, %d)",
                     d, length(mtx), length(bcf), length(ftf)), call. = FALSE)
    M  <- Matrix::readMM(mtx)
    bc <- utils::read.delim(bcf, header = FALSE, stringsAsFactors = FALSE)[[1]]
    ft <- utils::read.delim(ftf, header = FALSE, stringsAsFactors = FALSE)
    # Symbol column: the one with the most MT- genes (10x has id, name, type).
    nmt <- vapply(ft, function(z) sum(grepl("^MT-", z)), numeric(1))
    sym <- ft[[if (any(nmt > 0)) which.max(nmt) else ncol(ft)]]
    cells_by_row <- nrow(M) == length(bc) && ncol(M) == length(sym)
    cells_by_col <- ncol(M) == length(bc) && nrow(M) == length(sym)
    if (!cells_by_row && !cells_by_col)
        stop(sprintf("F4: matrix %d x %d does not match %d barcodes and %d features in %s",
                     nrow(M), ncol(M), length(bc), length(sym), d), call. = FALSE)
    if (cells_by_col) M <- Matrix::t(M)
    M  <- methods::as(M, "CsparseMatrix")
    mt <- grepl("^MT-", sym)
    tot <- Matrix::rowSums(M)
    nm  <- if (any(mt)) Matrix::rowSums(M[, !mt, drop = FALSE]) else rep(NA_real_, length(bc))
    out <- data.frame(barcode = bc, total = as.numeric(tot), non_mito = as.numeric(nm),
                      n_mt_genes = sum(mt), stringsAsFactors = FALSE)
    found <- character(0)
    for (p in names(panels)) {
        k <- sym %in% panels[[p]]
        out[[paste0("p_", p)]] <- if (any(k))
            as.numeric(Matrix::rowSums(M[, k, drop = FALSE])) else NA_real_
        found[p] <- sprintf("%s %d/%d", p, length(unique(sym[k])),
                            length(panels[[p]]))
    }
    attr(out, "panels_found") <- found
    out
}

# Non-Module cells: Post-IP when composite >= IP_y. IP_y comes from FU_F4_IPY
# or is bracketed by the Module cells, whose labels follow the same rule.
fu_split_nonmodule <- function(L) {
    rows <- list()
    for (a in names(FU_ARMS)) for (ds in FU_DATASETS) {
        k  <- which(L$arm == a & L$rs == ds)
        if (!length(k)) next
        km <- k[L$mod[k] %in% TRUE & !is.na(L$post[k]) & is.finite(L$comp[k])]
        if (!any(L$post[km]) || all(L$post[km])) {
            cat(sprintf("  %-5s %-7s no Module cells on one side of the IP; non-Module cells not split\n",
                        a, ds)); next
        }
        lo <- max(L$comp[km][!L$post[km]]); hi <- min(L$comp[km][L$post[km]])
        if (!(lo < hi))
            stop(sprintf(paste0("F4: %s %s Module Pre-IP and Post-IP composites",
                                " overlap (max Pre %.6f, min Post %.6f), so the",
                                " labels do not follow composite >= IP_y"),
                         a, ds, lo, hi), call. = FALSE)
        key <- paste(a, ds, sep = ".")
        ipy <- if (exists("FU_F4_IPY") && key %in% names(FU_F4_IPY))
            as.numeric(FU_F4_IPY[[key]]) else NA_real_
        if (is.finite(ipy) && !(ipy > lo && ipy <= hi))
            stop(sprintf(paste0("F4: FU_F4_IPY[\"%s\"] = %.6f is outside the",
                                " Module gap (%.6f, %.6f]"), key, ipy, lo, hi),
                 call. = FALSE)
        kn <- k[L$mod[k] %in% FALSE & is.na(L$post[k])]
        cmp <- L$comp[kn]
        side <- if (is.finite(ipy)) ifelse(cmp >= ipy, TRUE, FALSE) else
            ifelse(cmp >= hi, TRUE, ifelse(cmp <= lo, FALSE, NA))
        side[!is.finite(cmp)] <- NA
        L$post[kn] <- side
        rows[[key]] <- data.frame(
            arm = a, region_set = ds,
            ip_y_source = if (is.finite(ipy)) "FU_F4_IPY" else "Module gap",
            ip_y = if (is.finite(ipy)) ipy else NA_real_,
            module_pre_max = lo, module_post_min = hi,
            nonmodule_cells = length(kn), nonmodule_post = sum(side %in% TRUE),
            nonmodule_pre = sum(side %in% FALSE),
            nonmodule_in_gap = sum(is.finite(cmp) & is.na(side)),
            nonmodule_no_composite = sum(!is.finite(cmp)))
    }
    S <- do.call(rbind, rows); rownames(S) <- NULL
    list(L = L, split = S)
}

fu_F4 <- function(panels = list()) {
    fu_hdr("F4. UMIs PER NUCLEUS, POST-IP vs PRE-IP, WITHIN DONOR")
    # Labels from Part12_04_Labels.r
    f <- if (exists("FU_F4_FILE")) FU_F4_FILE else
        list.files(file.path(CEP_HOME, "Part12_rerun", "04_Labels"),
                   pattern = "^04_cell_labels\\.csv(\\.gz)?$",
                   full.names = TRUE)
    if (length(f) != 1L)
        stop(sprintf(paste0("F4: expected one 04_cell_labels.csv(.gz) under",
                            " Part12_rerun/04_Labels, found %d%s. Set FU_F4_FILE."), length(f),
                     if (length(f)) paste0(": ", paste(f, collapse = "; ")) else ""),
             call. = FALSE)
    cat("  Labels:", f, "\n")
    # fileEncoding drops a byte-order mark from the first column name.
    X <- utils::read.csv(f, check.names = FALSE, stringsAsFactors = FALSE,
                         fileEncoding = "UTF-8-BOM")
    names(X) <- trimws(sub("^(\ufeff|\u00ef\u00bb\u00bf)", "", names(X)))
    nm <- names(X)
    col <- vapply(FU_F4_COLS, fu_pick, character(1), nm = nm)
    need <- c("barcode", "arm", "region_set", "donor", "module", "ip_side")
    miss <- need[is.na(col[need])]
    if (length(miss)) {
        cat("  Columns in the file (escaped, so hidden characters show):\n   ",
            paste(encodeString(nm, quote = "'"), collapse = ", "), "\n")
        stop("F4: no column found for: ", paste(miss, collapse = ", "),
             ". Set FU_F4_COLS_OVERRIDE to name them.", call. = FALSE)
    }
    cat("  Columns used:", paste(sprintf("%s = '%s'", names(col)[!is.na(col)],
                                         col[!is.na(col)]), collapse = "; "), "\n")
    cat("  Raw values seen (check the mapping):\n")
    for (k in c("arm", "region_set", "module", "ip_side"))
        cat(sprintf("    %-10s: %s\n", k, paste(utils::head(sort(unique(as.character(
            X[[col[[k]]]]))), 12), collapse = ", ")))
    L <- data.frame(
        barcode = if (!is.na(col[["barcode"]])) as.character(X[[col[["barcode"]]]]) else NA,
        arm = fu_arm_std(X[[col[["arm"]]]]), rs = fu_rs_std(X[[col[["region_set"]]]]),
        donor = as.character(X[[col[["donor"]]]]),
        mod = fu_lgl_map(X[[col[["module"]]]], "^module|^yes|refined|recovered",
                         "^nonmodule|^no|originalonly|excluded"),
        post = fu_lgl_map(X[[col[["ip_side"]]]], "post", "pre"),
        comp = if (!is.na(col[["composite"]]))
            suppressWarnings(as.numeric(X[[col[["composite"]]]])) else NA_real_,
        stringsAsFactors = FALSE)
    cat(sprintf("  Labels mapped (NA counts): arm %d, region set %d, module %d, IP side %d of %d\n",
                sum(is.na(L$arm)), sum(is.na(L$rs)), sum(is.na(L$mod)),
                sum(is.na(L$post)), nrow(L)))
    cat(sprintf("  IP side NA: %d Module, %d non-Module\n",
                sum(is.na(L$post) & L$mod %in% TRUE),
                sum(is.na(L$post) & L$mod %in% FALSE)))

    # Non-Module IP side from the composite
    SPLIT <- NULL
    if (any(is.na(L$post) & L$mod %in% FALSE)) {
        if (is.na(col[["composite"]]) || all(!is.finite(L$comp)))
            cat("  No composite column, so non-Module cells cannot be split at IP_y\n")
        else {
            r <- fu_split_nonmodule(L); L <- r$L; SPLIT <- r$split
            cat("\n  Non-Module cells split at IP_y on the composite:\n")
            print(SPLIT, row.names = FALSE)
            cat("\n")
            write.csv(SPLIT, file.path(FU_OUT, "S11_F4_nonmodule_split.csv"),
                      row.names = FALSE)
        }
    }

    # UMIs per nucleus, total and non-MT, from the Part 8 raw counts
    measures <- c("total", "non_mito")
    if (!requireNamespace("Matrix", quietly = TRUE))
        stop("F4: package Matrix is required", call. = FALSE)
    L$total <- NA_real_; L$non_mito <- NA_real_
    for (p in names(panels)) L[[paste0("p_", p)]] <- NA_real_
    root <- file.path(CEP_HOME, "R_export_CP10K")
    for (a in names(FU_F4_ARM_TAGS)) {
        adir <- file.path(root, FU_F4_ARM_TAGS[[a]])
        raws <- list.dirs(adir, recursive = TRUE, full.names = TRUE)
        raws <- raws[basename(raws) == "raw_selected" &
                     basename(dirname(raws)) == "counts"]
        for (ds in FU_DATASETS) {
            d <- raws[fu_rs_std(basename(dirname(dirname(raws)))) %in% ds]
            k <- which(L$arm == a & L$rs == ds)
            if (length(d) != 1L) {
                cat(sprintf("  %-5s %-7s raw_selected folders found: %d under %s - skipped\n",
                            a, ds, length(d), adir)); next
            }
            U <- fu_umi_from_raw(d, panels)
            m <- match(L$barcode[k], U$barcode)
            L$total[k] <- U$total[m]; L$non_mito[k] <- U$non_mito[m]
            for (p in names(panels))
                L[[paste0("p_", p)]][k] <- U[[paste0("p_", p)]][m]
            cat(sprintf("  %-5s %-7s %s: %d barcodes, %d MT- genes, labels matched %d of %d%s\n",
                        a, ds, basename(dirname(dirname(d))), nrow(U),
                        U$n_mt_genes[1], sum(!is.na(m)), length(k),
                        if (length(panels)) paste0("; panel genes found: ",
                            paste(attr(U, "panels_found"), collapse = ", ")) else ""))
            if (length(k) && mean(!is.na(m)) < 0.9)
                stop(sprintf(paste0("F4: only %.1f%% of %s %s labels matched a",
                                    " raw_selected barcode; check the barcode",
                                    " format"), 100 * mean(!is.na(m)), a, ds),
                     call. = FALSE)
        }
    }
    if (all(is.na(L$non_mito))) measures <- "total"
    ok <- !is.na(L$arm) & !is.na(L$rs) & !is.na(L$mod) & !is.na(L$post) &
          is.finite(L$total) & L$total > 0
    cat(sprintf("  Cells usable: %d of %d\n", sum(ok), nrow(L)))
    if (!any(ok)) stop("F4: no usable cells after mapping", call. = FALSE)
    CELLS <- L[ok, ]
    D <- CELLS
    D$total <- log2(D$total)
    D$non_mito <- ifelse(is.finite(D$non_mito) & D$non_mito > 0, log2(D$non_mito), NA)

    by_d <- list(); summ <- list()
    for (ms in measures)
    for (a in names(FU_ARMS)) for (ds in FU_DATASETS) for (m in c(TRUE, FALSE)) {
        s <- D[D$arm == a & D$rs == ds & D$mod == m & is.finite(D[[ms]]), ]
        if (!nrow(s)) next
        s$lu <- s[[ms]]
        rows <- lapply(split(s, s$donor), function(x) {
            n1 <- sum(x$post); n0 <- sum(!x$post)
            if (n1 < FU_F4_MIN_CELLS || n0 < FU_F4_MIN_CELLS) return(NULL)
            data.frame(measure = ms, arm = a, region_set = ds,
                       cells = if (m) "Module" else "non-Module",
                       donor = x$donor[1], n_post = n1, n_pre = n0,
                       med_log2UMI_post = stats::median(x$lu[x$post]),
                       med_log2UMI_pre  = stats::median(x$lu[!x$post]))
        })
        dd <- do.call(rbind, rows)
        if (is.null(dd) || !nrow(dd)) next
        dd$log2_ratio <- dd$med_log2UMI_post - dd$med_log2UMI_pre
        by_d[[length(by_d) + 1L]] <- dd
        p <- if (nrow(dd) >= 5L) suppressWarnings(stats::wilcox.test(
            dd$log2_ratio, exact = FALSE, correct = TRUE)$p.value) else NA_real_
        summ[[length(summ) + 1L]] <- data.frame(
            measure = ms, arm = a, region_set = ds, cells = dd$cells[1],
            n_donors = nrow(dd),
            median_log2_ratio_post_vs_pre = round(stats::median(dd$log2_ratio), 3),
            fold_post_vs_pre = round(2^stats::median(dd$log2_ratio), 3),
            n_donors_post_lower = sum(dd$log2_ratio < 0),
            p_signed_rank = signif(p, 3))
    }
    if (!length(summ)) stop("F4: no arm and region set had paired donors",
                            call. = FALSE)
    BD <- do.call(rbind, by_d); S <- do.call(rbind, summ)
    S$BH <- NA_real_
    for (ms in unique(S$measure)) for (cc in unique(S$cells)) {
        k <- S$measure == ms & S$cells == cc & is.finite(S$p_signed_rank)
        S$BH[k] <- signif(stats::p.adjust(S$p_signed_rank[k], "BH"), 3)
    }
    S <- S[order(match(S$measure, measures), S$cells != "Module",
                 match(S$arm, names(FU_ARMS)), S$region_set), ]
    rownames(S) <- NULL
    cat(sprintf("\n  Within donor (>= %d cells per group), Post-IP minus Pre-IP median log2 UMIs:\n",
                FU_F4_MIN_CELLS))
    print(S, row.names = FALSE)
    write.csv(BD, file.path(FU_OUT, "S11_F4_umi_by_donor.csv"), row.names = FALSE)
    write.csv(S, file.path(FU_OUT, "S11_F4_umi_summary.csv"), row.names = FALSE)
    invisible(list(by_donor = BD, summary = S, split = SPLIT, cells = CELLS))
}

# =============================================================================
# F5. Panel UMIs per nucleus: amount and share, Post-IP against Pre-IP
# =============================================================================
# Per donor: amount = mean panel UMIs per nucleus; share = panel UMIs over
# non-MT UMIs. Post/Pre log2 ratio, median over donors, signed-rank test, BH.
fu_F5 <- function(F4, panels) {
    fu_hdr("F5. PANEL UMIs PER NUCLEUS, POST-IP vs PRE-IP, WITHIN DONOR")
    if (is.null(F4) || is.null(F4$cells))
        stop("F5 needs the F4 cells; F4 did not run", call. = FALSE)
    C <- F4$cells
    for (p in names(panels))
        cat(sprintf("  %-15s: %s\n", p, paste(panels[[p]], collapse = ", ")))
    by_d <- list(); summ <- list()
    for (p in names(panels)) for (m in c(TRUE, FALSE))
    for (a in names(FU_ARMS)) for (ds in FU_DATASETS) {
        pc <- paste0("p_", p)
        s <- C[C$arm == a & C$rs == ds & C$mod == m & is.finite(C[[pc]]) &
               is.finite(C$non_mito) & C$non_mito > 0, ]
        if (!nrow(s)) next
        n_zero <- 0L
        rows <- lapply(split(s, s$donor), function(x) {
            n1 <- sum(x$post); n0 <- sum(!x$post)
            if (n1 < FU_F4_MIN_CELLS || n0 < FU_F4_MIN_CELLS) return(NULL)
            a1 <- mean(x[[pc]][x$post]); a0 <- mean(x[[pc]][!x$post])
            if (a1 == 0 || a0 == 0) { n_zero <<- n_zero + 1L; return(NULL) }
            o1 <- mean(x$non_mito[x$post]); o0 <- mean(x$non_mito[!x$post])
            data.frame(panel = p, cells = if (m) "Module" else "non-Module",
                       arm = a, region_set = ds, donor = x$donor[1],
                       n_post = n1, n_pre = n0,
                       panel_umi_per_nucleus_post = a1,
                       panel_umi_per_nucleus_pre = a0,
                       share_post = sum(x[[pc]][x$post]) / sum(x$non_mito[x$post]),
                       share_pre  = sum(x[[pc]][!x$post]) / sum(x$non_mito[!x$post]),
                       log2_amount = log2(a1 / a0),
                       log2_own_rna = log2(o1 / o0))
        })
        dd <- do.call(rbind, rows)
        if (is.null(dd) || !nrow(dd)) next
        dd$log2_share <- log2(dd$share_post / dd$share_pre)
        by_d[[length(by_d) + 1L]] <- dd
        sr <- function(v) if (length(v) >= 5L) suppressWarnings(stats::wilcox.test(
            v, exact = FALSE, correct = TRUE)$p.value) else NA_real_
        summ[[length(summ) + 1L]] <- data.frame(
            panel = p, cells = dd$cells[1], arm = a, region_set = ds,
            n_donors = nrow(dd), n_donors_zero = n_zero,
            umi_per_nucleus_pre = round(stats::median(dd$panel_umi_per_nucleus_pre), 2),
            fold_own_rna = round(2^stats::median(dd$log2_own_rna), 3),
            fold_amount = round(2^stats::median(dd$log2_amount), 3),
            p_amount = signif(sr(dd$log2_amount), 3),
            fold_share = round(2^stats::median(dd$log2_share), 3),
            p_share = signif(sr(dd$log2_share), 3),
            stringsAsFactors = FALSE)
    }
    if (!length(summ)) stop("F5: no arm and region set had paired donors",
                            call. = FALSE)
    BD <- do.call(rbind, by_d); S <- do.call(rbind, summ)
    S$BH_amount <- NA_real_; S$BH_share <- NA_real_
    for (p in unique(S$panel)) for (cc in unique(S$cells)) {
        k <- S$panel == p & S$cells == cc
        S$BH_amount[k] <- signif(stats::p.adjust(S$p_amount[k], "BH"), 3)
        S$BH_share[k]  <- signif(stats::p.adjust(S$p_share[k], "BH"), 3)
    }
    S <- S[, c("panel", "cells", "arm", "region_set", "n_donors", "n_donors_zero",
               "umi_per_nucleus_pre", "fold_own_rna", "fold_amount", "p_amount",
               "BH_amount", "fold_share", "p_share", "BH_share")]
    rownames(S) <- NULL
    cat(sprintf("\n  Within donor (>= %d cells per group), medians over donors:\n",
                FU_F4_MIN_CELLS))
    for (cc in c("Module", "non-Module")) {
        cat(sprintf("\n  %s cells:\n", cc))
        print(S[S$cells == cc, setdiff(names(S), "cells")], row.names = FALSE)
    }
    write.csv(BD, file.path(FU_OUT, "S11_F5_panels_by_donor.csv"), row.names = FALSE)
    write.csv(S, file.path(FU_OUT, "S11_F5_panels_summary.csv"), row.names = FALSE)
    invisible(list(by_donor = BD, summary = S))
}

# =============================================================================
# Run
# =============================================================================
.fu_con <- file(FU_LOG, open = "wt")
.fu_opt <- options(width = 200)
sink(.fu_con, split = TRUE)
tryCatch({
    fu_hdr("S11 FOLLOW-UP ANALYSES")
    cat("Date     :", format(Sys.time()), "\n")
    cat("CEP_HOME :", CEP_HOME, "\n")
    cat("R        :", R.version.string, "\n\n")
    cat("Inputs:\n")
    FU_RDS <- fu_read_rds()
    FU_F1  <- fu_F1(FU_RDS)
    FU_F2  <- fu_F2(FU_RDS)
    FU_F3  <- fu_F3(FU_RDS, FU_F1$summary$gene)
    FU_PANELS <- list(
        long_neuronal  = FU_F1$summary$gene[FU_F1$summary$direction == "up"],
        short_neuronal = FU_SHORT_NEURO,
        microglial_own = FU_MG_OWN)
    FU_F4  <- tryCatch(fu_F4(FU_PANELS), error = function(e) {
        cat("\n  F4 not run:", conditionMessage(e), "\n"); NULL })
    FU_F5  <- tryCatch(fu_F5(FU_F4, FU_PANELS), error = function(e) {
        cat("\n  F5 not run:", conditionMessage(e), "\n"); NULL })
    fu_hdr("DONE")
    cat("Outputs in", FU_OUT, "\n")
}, error = function(e) cat("\nSTOPPED:", conditionMessage(e), "\n"),
   finally = {
       while (sink.number() > 0) sink()
       try(close(.fu_con), silent = TRUE)
       options(.fu_opt)
       cat("Log written to:", FU_LOG, "\n")
   })
