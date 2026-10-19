################################################################################
# build_SuppTable_S11.r
#
# Builds SuppTable_S11_DEGs_ModvsNonMod.xlsx, the source workbook of
# Supplementary Table S11: differential expression between Module and
# non-Module cells, and between Post-IP and Pre-IP cells, on the four arms
# (M4, PVM_3, PVM_1, HKG8) and the two region sets (NeoC, NonNeoC). Sheet A
# counts the significant genes of the Part 10g Part 1 donor-unit contrasts.
# Sheets B to D hold the test of the PVM_3 DEGs against HKG8, the axis check
# on non-Module cells and the Module agreement between arms (Part 10g Steps
# 7b, 7c, 8 and 9). Sheets E and F hold the survivors in both region sets and
# the neuronal gene panels. Sheets G and H hold UMIs per nucleus and panel
# UMIs per nucleus from the Part 8 raw counts. Every count in sheet A is
# recomputed from the Part 1 donor frames and compared with Part 1's own
# Summary, and sheets B, C, E and F are checked against the RDS gene by gene.
# Any FAIL stops the build before the workbook is written.
# finalize_supp_tables.py turns this workbook into Supplementary Table S11;
# the Contents sheet written here is replaced there.
#
# Usage (from the SEA-AD root folder, or with SEAAD_ROOT set):
#   Rscript build_SuppTable_S11.r
#   or in R: CEP_HOME <- "/data"; source("build_SuppTable_S11.r")
# CEP_HOME defaults to SEAAD_ROOT, else the current folder.
# Packages: readxl, writexl, Matrix.
#
# Inputs (under CEP_HOME; nothing is modified):
#   Part10g_DEGs/Part10g_DEGs_<arm>/*_<NeoC|NonNeoC>_results.rds
#       (Part10g_CP10K_DEGs_Part1_code.r)
#   Part10e_Validation1/Part10e_Validation1_<arm>/*_V1_<NeoC|NonNeoC>_SK_interaction_analysis.xlsx
#       (Part10e_CP10K_Validation1_code.r)
#   Part10g_DEGs/Step7to9_rebuild/Part10g_Step7b_targeted_tests_verdict.csv
#   Part10g_DEGs/Step7to9_rebuild/Part10g_Step7c_PvP_targeted_tests.csv
#   Part10g_DEGs/Step7to9_rebuild/Part10g_Step8_control_overlap.csv
#   Part10g_DEGs/Step7to9_rebuild/Part10g_Step9_ambient_check.csv
#       (Part10g_Step7to9_rebuild_code.r)
#   Part12_rerun/04_Labels/04_cell_labels.csv.gz   (Part12_04_Labels.r)
#   R_export_CP10K/<ARM_TAG>/<region set>/counts/raw_selected/
#       (Part 8 data export: matrix.mtx, barcodes and features)
#   <arm> is the folder name Mic13_PVM3, PVM_3, PVM_1 or HKG8; <ARM_TAG> is
#   Mic13_PVM3, PVM3_all, PVM1_all or HKG8_matched.
#
# Output (CEP_HOME/Manuscript/Supp_Tables/build/):
#   SuppTable_S11_DEGs_ModvsNonMod.xlsx, with the sheets Contents, A_Contrasts,
#   B_Control_test, B_Summary, C_Axis_check, D_Control_validity,
#   E_Survivors_regions, F_Neuronal_summary, F_Neuronal_genes, G_UMI_summary,
#   G_UMI_donors, G_IP_split, H_Panel_summary and H_Panel_donors. Bold header
#   only; no fills, borders or colours.
#   SuppTable_S11_build_log_<timestamp>.txt, ending with VERDICT: PASS or FAIL.
#
# Two steps restore the layout of the published table from the rebuilt Step
# 7 to 9 files. The p-values and FDR of the Step 7b verdict file are rounded
# as in the file the published table was built from: values below 1e-4 to
# three significant digits, other values to nine decimals. Rows of sheet C
# that passed the control test are ordered by FDR_M4 within region set.
################################################################################

if (!exists("CEP_HOME")) CEP_HOME <- Sys.getenv("SEAAD_ROOT", getwd())
stopifnot(dir.exists(CEP_HOME))
for (pk in c("readxl", "writexl", "Matrix"))
    if (!requireNamespace(pk, quietly = TRUE))
        stop("Package ", pk, " is required.", call. = FALSE)

S11_OUT  <- file.path(CEP_HOME, "Manuscript", "Supp_Tables", "build")
dir.create(S11_OUT, showWarnings = FALSE, recursive = TRUE)
.ts      <- format(Sys.time(), "%Y%m%d_%H%M%S")
.log_txt <- file.path(S11_OUT, paste0("SuppTable_S11_build_log_", .ts, ".txt"))
S11_XLSX <- file.path(S11_OUT, "SuppTable_S11_DEGs_ModvsNonMod.xlsx")

.old_opt <- options(width = 200, warn = 1, stringsAsFactors = FALSE)
.log_con <- file(.log_txt, open = "wt")
sink(.log_con, split = TRUE)

s11_build <- function() {

    # ---- Settings ----------------------------------------------------------
    DEG_DIR  <- file.path(CEP_HOME, "Part10g_DEGs")
    V1_DIR   <- file.path(CEP_HOME, "Part10e_Validation1")
    ARMS     <- c(M4 = "Mic13_PVM3", PVM_3 = "PVM_3", PVM_1 = "PVM_1",
                  HKG8 = "HKG8")             # label = folder subdir
    DATASETS <- c("NeoC", "NonNeoC")
    CONTRASTS <- c("PostIP_vs_NMPost", "PreIP_vs_NMPre", "PostIP_vs_PreIP",
                   "NMPost_vs_NMPre")
    GROUPS <- list(
        PostIP_vs_NMPost = c("Module Post-IP", "non-Module Post-IP"),
        PreIP_vs_NMPre   = c("Module Pre-IP", "non-Module Pre-IP"),
        PostIP_vs_PreIP  = c("Module Post-IP", "Module Pre-IP"),
        NMPost_vs_NMPre  = c("non-Module Post-IP", "non-Module Pre-IP"))
    MITO_RE  <- "^MT-|^MTRNR2L|^TSPO$"   # circular for this module
    ALPHA    <- 0.05
    TOL      <- 1e-3                     # CSV values are rounded
    MODULE_CLASSES    <- c("Original + Refined", "Recovered")
    NONMODULE_CLASSES <- c("Original only", "Excluded")
    # Arms tested against HKG8. The control test uses the two-sided
    # FDR2_<arm> columns of Step 7b/7c (BH within contrast and region set);
    # checks S11-11 and S11-12 assume this.
    CTRL_ARMS <- c("M4", "PVM_1")
    B_CONTRASTS <- c("PostIP_vs_NMPost", "PreIP_vs_NMPre", "PostIP_vs_PreIP")
    # Steps 7b, 7c, 8 and 9 as written by Part10g_Step7to9_rebuild_code.r
    S79_DIR  <- file.path(DEG_DIR, "Step7to9_rebuild")
    # Cell labels as written by Part12_04_Labels.r
    F_LABELS <- file.path(CEP_HOME, "Part12_rerun", "04_Labels",
                          "04_cell_labels.csv.gz")
    # Rounding of the p-values and FDR in the Step 7b verdict file the
    # published table was built from: below 1e-4 three significant digits,
    # otherwise nine decimals.
    round_7b <- function(x) {
        x <- suppressWarnings(as.numeric(x))
        ifelse(is.finite(x) & abs(x) < 1e-4, signif(x, 3), round(x, 9))
    }

    # ---- Check ledger ------------------------------------------------------
    LEDGER <- list()
    chk <- function(id, status, detail) {
        if (length(detail) != 1L || is.na(detail)) detail <- "[no detail]"
        cat(sprintf("  [%s]  %-7s %s\n", status, id, detail))
        LEDGER[[length(LEDGER) + 1L]] <<- data.frame(id = id, status = status,
                                                     detail = detail)
    }
    hdr <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n",
                           strrep("=", 78), "\n", sep = "")
    as_lgl <- function(x) toupper(trimws(as.character(x))) %in% "TRUE"
    num <- function(x) suppressWarnings(as.numeric(x))

    # Smallest two-sided p the signed-rank test can give with n paired
    # donors, all differences of one sign. Normal approximation with
    # continuity correction, as used by Part 1 (exact = FALSE, correct = TRUE).
    sr_p_floor <- function(n) {
        if (!is.finite(n) || n < 1) return(NA_real_)
        mu <- n * (n + 1) / 4
        sg <- sqrt(n * (n + 1) * (2 * n + 1) / 24)
        2 * stats::pnorm(-(n * (n + 1) / 2 - mu - 0.5) / sg)
    }

    hdr("SuppTable_S11 BUILD")
    cat("Date     :", format(Sys.time()), "\n")
    cat("CEP_HOME :", CEP_HOME, "\n")
    cat("R        :", R.version.string, "\n")

    # =========================================================================
    # 1. Read Part 1 RDS and Validation 1 cell classes
    # =========================================================================
    hdr("1. INPUTS")
    RDS <- list(); V1 <- list()
    for (lab in names(ARMS)) {
        sub <- ARMS[[lab]]
        d1  <- file.path(DEG_DIR, paste0("Part10g_DEGs_", sub))
        dv  <- file.path(V1_DIR, paste0("Part10e_Validation1_", sub))
        for (ds in DATASETS) {
            key <- paste(lab, ds)
            f <- list.files(d1, pattern = paste0("_", ds, "_results\\.rds$"),
                            full.names = TRUE)
            if (length(f) != 1L) {
                chk("S11-01", "FAIL", sprintf("%s: %d Part 1 RDS found in %s",
                                              key, length(f), d1))
                next
            }
            RDS[[key]] <- readRDS(f)
            cat(sprintf("  %-13s %s (%s)\n", key, basename(f),
                        format(file.info(f)$mtime, "%Y-%m-%d %H:%M")))
            fv <- list.files(dv, pattern = paste0("_V1_", ds,
                                                  "_SK_interaction_analysis\\.xlsx$"),
                             full.names = TRUE)
            if (length(fv) == 1L) {
                rc <- readxl::read_excel(fv, sheet = "Step4_Refined_Cells")
                V1[[key]] <- table(factor(rc$cell_class,
                                          levels = c(MODULE_CLASSES,
                                                     NONMODULE_CLASSES)))
                rm(rc)
            } else {
                chk("S11-01", "WARN", sprintf("%s: Validation 1 workbook not found",
                                              key))
            }
        }
    }
    chk("S11-01", if (length(RDS) == 8L) "PASS" else "FAIL",
        sprintf("%d of 8 Part 1 RDS read; %d of 8 Validation 1 workbooks read",
                length(RDS), length(V1)))

    # =========================================================================
    # 2. Sheet A: one row per arm, region set and contrast
    # =========================================================================
    hdr("2. SHEET A (contrasts)")
    A <- list(); n_cnt_bad <- 0L; n_cell_bad <- 0L
    for (key in names(RDS)) {
        o <- RDS[[key]]; lab <- sub(" .*", "", key); ds <- sub(".* ", "", key)
        sm <- o[["summary"]]
        for (tg in CONTRASTS) {
            d <- o[["res"]][[tg]][["donor"]]
            if (!is.data.frame(d) || nrow(d) == 0L) {
                chk("S11-02", "FAIL", sprintf("%s %s: no donor frame", key, tg))
                next
            }
            sig  <- as_lgl(d$significant)
            mito <- grepl(MITO_RE, d$gene)
            nd   <- num(d$n_donors[1])
            pfl  <- sr_p_floor(nd)
            p_adj <- num(d$p_val_adj_BH)
            row <- data.frame(
                arm = lab, region_set = ds, contrast = tg,
                group_1 = GROUPS[[tg]][1], group_2 = GROUPS[[tg]][2],
                n_cells_group_1 = num(d$n_fg[1]),
                n_cells_group_2 = num(d$n_bg[1]),
                n_paired_donors = nd, n_genes_tested = nrow(d),
                log2FC_gate = NA_real_,
                n_significant = sum(sig),
                n_significant_up = sum(sig & num(d$avg_log2FC) > 0),
                n_significant_down = sum(sig & num(d$avg_log2FC) < 0),
                n_significant_mito = sum(sig & mito),
                n_significant_nonmito = sum(sig & !mito),
                min_BH = if (all(is.na(p_adj))) NA_real_ else
                             signif(min(p_adj, na.rm = TRUE), 4),
                p_floor_signed_rank = signif(pfl, 4),
                single_gene_BH_floor = signif(min(1, pfl * nrow(d)), 4))
            row$resolution_limited <- is.finite(row$single_gene_BH_floor) &&
                row$single_gene_BH_floor > ALPHA && row$n_significant == 0L
            # Compare with Part 1's own Summary row
            if (is.data.frame(sm)) {
                k <- which(sm$comparison == tg &
                           (if ("unit" %in% names(sm)) sm$unit == "donor" else TRUE))
                if (length(k) == 1L) {
                    row$log2FC_gate <- num(sm$sig_logfc[k])
                    ok <- num(sm$n_significant[k]) == row$n_significant &&
                          num(sm$n_up[k]) == row$n_significant_up &&
                          num(sm$n_down[k]) == row$n_significant_down &&
                          num(sm$n_donors[k]) == row$n_paired_donors &&
                          num(sm$n_tested[k]) == row$n_genes_tested
                    if (!isTRUE(ok)) {
                        n_cnt_bad <- n_cnt_bad + 1L
                        cat(sprintf("    mismatch with Part 1 Summary: %s %s\n",
                                    key, tg))
                    }
                } else n_cnt_bad <- n_cnt_bad + 1L
            } else n_cnt_bad <- n_cnt_bad + 1L
            A[[length(A) + 1L]] <- row
        }
        # Cell counts: module and non-module totals against Validation 1
        r_pp <- o[["res"]][["PostIP_vs_PreIP"]][["donor"]]
        r_nm <- o[["res"]][["NMPost_vs_NMPre"]][["donor"]]
        n_mod <- num(r_pp$n_fg[1]) + num(r_pp$n_bg[1])
        n_non <- num(r_nm$n_fg[1]) + num(r_nm$n_bg[1])
        ok_rds <- n_mod == num(o[["n_module"]]) && n_non == num(o[["n_nonmod"]])
        ok_v1 <- TRUE
        if (!is.null(V1[[key]])) {
            ok_v1 <- n_mod == sum(V1[[key]][MODULE_CLASSES]) &&
                     n_non == sum(V1[[key]][NONMODULE_CLASSES])
        }
        if (!isTRUE(ok_rds && ok_v1)) n_cell_bad <- n_cell_bad + 1L
        cat(sprintf("  %-13s Module %d (Pre %s / Post %s) | non-Module %d (Pre %s / Post %s)%s\n",
                    key, n_mod, r_pp$n_bg[1], r_pp$n_fg[1], n_non,
                    r_nm$n_bg[1], r_nm$n_fg[1],
                    if (isTRUE(ok_rds && ok_v1)) "" else "  <- MISMATCH"))
    }
    A <- do.call(rbind, A)
    A$arm <- factor(A$arm, levels = names(ARMS))
    A$contrast <- factor(A$contrast, levels = CONTRASTS)
    A <- A[order(A$contrast, A$arm, A$region_set), ]
    A$arm <- as.character(A$arm); A$contrast <- as.character(A$contrast)
    rownames(A) <- NULL
    chk("S11-02", if (nrow(A) == 32L && n_cnt_bad == 0L) "PASS" else "FAIL",
        sprintf("%d of 32 rows; %d disagree with Part 1's Summary", nrow(A),
                n_cnt_bad))
    chk("S11-03", if (n_cell_bad == 0L) "PASS" else "FAIL",
        sprintf(paste0("Module and non-Module cell totals agree with the RDS",
                       " and Validation 1 on %d of %d arm/region sets"),
                length(RDS) - n_cell_bad, length(RDS)))
    print(A[, c("arm", "region_set", "contrast", "n_paired_donors",
                "n_significant", "n_significant_mito", "n_significant_nonmito",
                "min_BH", "resolution_limited")], row.names = FALSE)

    # =========================================================================
    # 3. Sheet B: genes significant on PVM_3, tested against HKG8
    # =========================================================================
    hdr("3. SHEET B (test against the HKG8 control)")
    f7b <- file.path(S79_DIR, "Part10g_Step7b_targeted_tests_verdict.csv")
    f7c <- file.path(S79_DIR, "Part10g_Step7c_PvP_targeted_tests.csv")
    for (f in c(f7b, f7c)) if (!file.exists(f)) stop("not found: ", f,
                                                     call. = FALSE)
    s7b <- utils::read.csv(f7b, check.names = FALSE)
    s7c <- utils::read.csv(f7c, check.names = FALSE)
    cat(sprintf("  Step 7b: %s, %d rows\n  Step 7c: %s, %d rows\n", f7b,
                nrow(s7b), f7c, nrow(s7c)))
    # Round the Step 7b p-values and FDR as in the published table.
    n_rnd <- 0L
    for (cn in grep("^(p[12]_|FDR)", names(s7b), value = TRUE)) {
        v0 <- num(s7b[[cn]]); v1 <- round_7b(v0)
        n_rnd <- n_rnd + sum(is.finite(v0) & v0 != v1)
        s7b[[cn]] <- v1
    }
    cat(sprintf(paste0("  Step 7b p-values and FDR rounded (below 1e-4: three",
                       " significant digits; otherwise nine decimals): %d",
                       " values changed\n"), n_rnd))
    # Row order of the published table within each contrast and region set:
    # verdict (wrong direction, pass, one arm only, not significant), then
    # non-mitochondrial genes before mitochondrial ones, then file order.
    V7B <- c("wrong direction", "pass", "one arm only", "not significant")
    blk <- paste(s7b$comparison, s7b$dataset)
    s7b <- s7b[order(match(blk, unique(blk)), match(s7b$verdict, V7B),
                     as.logical(s7b$circular), seq_len(nrow(s7b))), ]
    cat("  Step 7b rows ordered by verdict, then non-mitochondrial before mitochondrial\n")
    s7c$verdict <- NA_character_
    keepB <- c("gene", "comparison", "dataset", "direction_on_PVM_3",
               "circular", "also_called_by_HKG8",
               "n_PVM_3", "med_PVM_3", "n_M4", "med_M4", "n_PVM_1", "med_PVM_1",
               "n_HKG8", "med_HKG8", "test_M4", "test_PVM_1",
               "concordant_M4", "concordant_PVM_1",
               "FDR1_M4", "FDR1_PVM_1", "FDR2_M4", "FDR2_PVM_1",
               "pass_one_sided", "pass_two_sided", "verdict")
    miss <- setdiff(keepB, union(names(s7b), names(s7c)))
    if (length(miss)) stop("Step 7b/7c CSV lacks: ", paste(miss, collapse = ", "),
                           call. = FALSE)
    B <- rbind(s7b[, keepB], s7c[, keepB])
    names(B)[names(B) == "dataset"] <- "region_set"
    B$mitochondrial <- grepl(MITO_RE, B$gene)
    n_flag_diff <- sum(as_lgl(B$circular) != B$mitochondrial)
    chk("S11-04", if (n_flag_diff == 0L) "PASS" else "WARN",
        sprintf(paste0("Step 7b/7c 'circular' flag agrees with the mitochondrial",
                       " pattern on %d of %d genes"), nrow(B) - n_flag_diff,
                nrow(B)))

    # ---- One control-test rule for every row -------------------------------
    for (a in CTRL_ARMS) {
        raw <- toupper(trimws(as.character(B[[paste0("concordant_", a)]])))
        con <- raw %in% "TRUE"
        fdr <- num(B[[paste0("FDR2_", a)]])
        sig <- is.finite(fdr) & fdr < ALPHA
        B[[paste0("concordant_", a)]] <- con
        B[[paste0("FDR_", a)]]  <- fdr
        B[[paste0("pass_", a)]] <- sig & con
        # 'wrong direction' only when concordance is known to be FALSE; a
        # missing concordance value is not evidence of the opposite sign.
        B[[paste0(".wrong_", a)]] <- sig & raw %in% "FALSE"
    }
    B$passes_control_test <- B$pass_M4 & B$pass_PVM_1
    v_7b <- B$verdict
    B$verdict <- ifelse(B$passes_control_test, "pass",
                 ifelse(B$.wrong_M4 | B$.wrong_PVM_1, "wrong direction",
                 ifelse(B$pass_M4 | B$pass_PVM_1, "one arm only",
                        "not significant")))
    has7b <- !is.na(v_7b) & nzchar(v_7b)
    n_v_bad <- sum(B$verdict[has7b] != v_7b[has7b])
    chk("S11-11", if (n_v_bad == 0L) "PASS" else "FAIL",
        sprintf(paste0("verdict from the control-test rule equals the Step 7b verdict on",
                       " %d of %d Module vs non-Module rows"),
                sum(has7b) - n_v_bad, sum(has7b)))
    n_2s_bad <- sum(B$passes_control_test != as_lgl(B$pass_two_sided))
    chk("S11-12", if (n_2s_bad == 0L) "PASS" else "FAIL",
        sprintf(paste0("passes_control_test equals Step 7b/7c pass_two_sided on",
                       " %d of %d rows"), nrow(B) - n_2s_bad, nrow(B)))
    # FDR1 is not written. Report how often it is larger than FDR2 for a
    # gene in the PVM_3 direction, which a one-sided version of the same
    # test with the same BH family cannot give.
    n_odd <- 0L; n_con <- 0L
    for (a in CTRL_ARMS) {
        f1 <- num(B[[paste0("FDR1_", a)]]); f2 <- num(B[[paste0("FDR2_", a)]])
        cc <- B[[paste0("concordant_", a)]] & is.finite(f1) & is.finite(f2)
        n_con <- n_con + sum(cc); n_odd <- n_odd + sum(f1[cc] > f2[cc] + 1e-12)
    }
    chk("S11-13", "NOTE",
        sprintf(paste0("FDR1 > FDR2 on %d of %d concordant gene x arm values;",
                       " FDR1 is not used"), n_odd, n_con))
    n_1s <- sum(as_lgl(B$pass_one_sided) != B$passes_control_test)
    if (n_1s > 0L) {
        cat("  Rows where the one-sided Step 7b/7c pass differs from the two-sided rule:\n")
        print(B[as_lgl(B$pass_one_sided) != B$passes_control_test,
                c("gene", "comparison", "region_set", "pass_one_sided",
                  "passes_control_test")], row.names = FALSE)
    }

    # ---- Every PVM_3 DEG is in sheet B --------------------------------------
    B$tested_against_HKG8 <- TRUE
    B$not_tested_reason   <- NA_character_
    add <- list(); n_set_bad <- 0L
    for (tg in B_CONTRASTS) for (ds in DATASETS) {
        o <- RDS[[paste("PVM_3", ds)]]; if (is.null(o)) next
        d <- o[["res"]][[tg]][["donor"]]
        rds_sig <- d$gene[as_lgl(d$significant)]
        csv_g   <- B$gene[B$comparison == tg & B$region_set == ds]
        a_only <- setdiff(csv_g, rds_sig); r_only <- setdiff(rds_sig, csv_g)
        cat(sprintf("  %-17s %-8s CSV %4d genes | PVM_3 significant %4d | CSV only %d | RDS only %d\n",
                    tg, ds, length(csv_g), length(rds_sig), length(a_only),
                    length(r_only)))
        if (length(a_only)) n_set_bad <- n_set_bad + 1L
        if (!length(r_only)) next
        fc3 <- o[["res"]][[tg]][["donor_fc"]]
        hk  <- RDS[[paste("HKG8", ds)]][["res"]][[tg]][["donor"]]
        for (g in r_only) {
            absent <- character(0)
            for (arm in c(CTRL_ARMS, "HKG8")) {
                da <- RDS[[paste(arm, ds)]][["res"]][[tg]][["donor"]]
                if (is.null(da) || !g %in% da$gene) absent <- c(absent, arm)
            }
            rsn <- if (length(absent))
                       paste("not tested on", paste(absent, collapse = ", "),
                             "for this contrast") else
                       "tested on all arms but absent from the Step 7b/7c output"
            r <- B[1, , drop = FALSE]; r[1, ] <- NA
            k <- match(g, d$gene)
            r$gene <- g; r$comparison <- tg; r$region_set <- ds
            lfc_k <- num(d$avg_log2FC[k])
            r$direction_on_PVM_3 <- if (!is.finite(lfc_k)) NA_character_ else
                                    if (lfc_k > 0) "up" else "down"
            r$mitochondrial <- grepl(MITO_RE, g)
            r$also_called_by_HKG8 <- !is.null(hk) && g %in% hk$gene &&
                as_lgl(hk$significant[match(g, hk$gene)])
            r$n_PVM_3 <- num(d$n_donors[k])
            r$med_PVM_3 <- if (is.matrix(fc3) && g %in% rownames(fc3))
                round(stats::median(fc3[g, ], na.rm = TRUE), 4) else NA_real_
            r$concordant_M4 <- FALSE; r$concordant_PVM_1 <- FALSE
            r$pass_M4 <- FALSE; r$pass_PVM_1 <- FALSE
            r$passes_control_test <- FALSE; r$verdict <- "not tested"
            r$tested_against_HKG8 <- FALSE; r$not_tested_reason <- rsn
            add[[length(add) + 1L]] <- r
            cat(sprintf("    added %-12s %-17s %-8s %s\n", g, tg, ds, rsn))
        }
    }
    if (length(add)) B <- rbind(B, do.call(rbind, add))
    chk("S11-05", if (n_set_bad == 0L) "PASS" else "FAIL",
        "every gene in sheet B is significant on PVM_3 for its contrast in the RDS")
    n_cov_bad <- 0L
    for (tg in B_CONTRASTS) for (ds in DATASETS) {
        nA <- A$n_significant[A$arm == "PVM_3" & A$contrast == tg &
                              A$region_set == ds]
        nB <- sum(B$comparison == tg & B$region_set == ds)
        if (length(nA) != 1L || nA != nB) n_cov_bad <- n_cov_bad + 1L
    }
    chk("S11-14", if (n_cov_bad == 0L) "PASS" else "FAIL",
        sprintf(paste0("sheet B rows equal the PVM_3 significant counts of",
                       " sheet A on %d of 6 contrast x region sets; %d untested",
                       " DEGs added"), 6L - n_cov_bad, length(add)))

    # Medians of the per-donor fold changes against the stored donor_fc.
    n_med <- 0L; n_med_ok <- 0L
    for (arm in c("PVM_3", "M4", "PVM_1", "HKG8")) for (ds in DATASETS) {
        o <- RDS[[paste(arm, ds)]]; if (is.null(o)) next
        for (tg in unique(B$comparison)) {
            fc <- o[["res"]][[tg]][["donor_fc"]]
            if (!is.matrix(fc)) next
            ii <- which(B$comparison == tg & B$region_set == ds)
            g  <- B$gene[ii]; m <- num(B[[paste0("med_", arm)]][ii])
            have <- g %in% rownames(fc) & is.finite(m)
            if (!any(have)) next
            mr <- apply(fc[g[have], , drop = FALSE], 1, stats::median, na.rm = TRUE)
            n_med <- n_med + sum(have)
            n_med_ok <- n_med_ok + sum(abs(mr - m[have]) <= TOL, na.rm = TRUE)
        }
    }
    chk("S11-06", if (n_med > 0L && n_med_ok == n_med) "PASS" else "WARN",
        sprintf(paste0("med_<arm> reproduced from the stored donor_fc within %g",
                       " on %d of %d gene x arm values"), TOL, n_med_ok, n_med))

    B$comparison <- factor(B$comparison, levels = B_CONTRASTS)
    B <- B[order(B$comparison, B$region_set, !B$tested_against_HKG8), ]
    B$comparison <- as.character(B$comparison)
    B <- B[, c("gene", "comparison", "region_set", "direction_on_PVM_3",
               "mitochondrial", "also_called_by_HKG8", "tested_against_HKG8",
               "not_tested_reason",
               "n_PVM_3", "med_PVM_3", "n_M4", "med_M4", "n_PVM_1", "med_PVM_1",
               "n_HKG8", "med_HKG8", "test_M4", "test_PVM_1",
               "concordant_M4", "concordant_PVM_1", "FDR_M4", "FDR_PVM_1",
               "pass_M4", "pass_PVM_1", "passes_control_test", "verdict")]
    rownames(B) <- NULL

    # ---- B_Summary -----------------------------------------------------------
    BS <- do.call(rbind, lapply(B_CONTRASTS, function(tg)
        do.call(rbind, lapply(DATASETS, function(ds) {
            s <- B[B$comparison == tg & B$region_set == ds, ]
            mt <- s$mitochondrial
            data.frame(
                contrast = tg, region_set = ds,
                PVM_3_DEGs = nrow(s), PVM_3_DEGs_mito = sum(mt),
                PVM_3_DEGs_nonmito = sum(!mt),
                tested_against_HKG8 = sum(s$tested_against_HKG8),
                pass_M4_mito = sum(s$pass_M4 & mt),
                pass_M4_nonmito = sum(s$pass_M4 & !mt),
                pass_PVM_1_mito = sum(s$pass_PVM_1 & mt),
                pass_PVM_1_nonmito = sum(s$pass_PVM_1 & !mt),
                pass_both_mito = sum(s$passes_control_test & mt),
                pass_both_nonmito = sum(s$passes_control_test & !mt),
                nonmito_genes_pass_both = paste(s$gene[s$passes_control_test & !mt],
                                                collapse = ", "))
        }))))
    rownames(BS) <- NULL
    cat("\n  Genes that pass against HKG8 (two-sided FDR < 0.05, PVM_3 direction):\n")
    print(BS, row.names = FALSE)

    # =========================================================================
    # 4. Sheet C: the axis check (Step 9)
    # =========================================================================
    hdr("4. SHEET C (non-Module cells crossing the inflection point)")
    f9 <- file.path(S79_DIR, "Part10g_Step9_ambient_check.csv")
    if (!file.exists(f9)) stop("not found: ", f9, call. = FALSE)
    C <- utils::read.csv(f9, check.names = FALSE)
    names(C)[names(C) == "dataset"] <- "region_set"
    n_c <- 0L; n_c_ok <- 0L
    for (arm in c("PVM_3", "M4", "PVM_1", "HKG8")) for (ds in DATASETS) {
        o <- RDS[[paste(arm, ds)]]; if (is.null(o)) next
        ii <- which(C$region_set == ds)
        for (pr in list(c("mod_", "PostIP_vs_PreIP"),
                        c("non_", "NMPost_vs_NMPre"))) {
            col <- paste0(pr[1], arm)
            if (!col %in% names(C)) next
            d <- o[["res"]][[pr[2]]][["donor"]]
            v <- num(d$avg_log2FC)[match(C$gene[ii], d$gene)]
            m <- num(C[[col]][ii]); have <- is.finite(v) & is.finite(m)
            n_c <- n_c + sum(have)
            n_c_ok <- n_c_ok + sum(abs(v[have] - m[have]) <= TOL)
        }
    }
    chk("S11-07", if (n_c > 0L && n_c_ok == n_c) "PASS" else "FAIL",
        sprintf(paste0("sheet C log2FC values equal the Part 1 donor frames",
                       " (PostIP_vs_PreIP and NMPost_vs_NMPre) on %d of %d"),
                n_c_ok, n_c))
    names(C) <- sub("^mod_", "log2FC_ModulePostVsPre_", names(C))
    names(C) <- sub("^non_sig_", "significant_nonModulePostVsPre_", names(C))
    names(C) <- sub("^non_", "log2FC_nonModulePostVsPre_", names(C))
    # Plain set labels, and the non-defining genes must be the
    # non-mitochondrial Post-IP vs Pre-IP survivors of sheet B.
    is_def <- grepl("defining", C$set, ignore.case = TRUE)
    C$set  <- ifelse(is_def, "defining gene", "passed control test")
    n_c_set_bad <- 0L
    for (ds in DATASETS) {
        gC <- C$gene[!is_def & C$region_set == ds]
        gB <- B$gene[B$comparison == "PostIP_vs_PreIP" & B$region_set == ds &
                     B$passes_control_test & !B$mitochondrial]
        if (!setequal(gC, gB)) {
            n_c_set_bad <- n_c_set_bad + 1L
            cat(sprintf("    %s: sheet C only: %s | sheet B only: %s\n", ds,
                        paste(setdiff(gC, gB), collapse = ", "),
                        paste(setdiff(gB, gC), collapse = ", ")))
        }
    }
    chk("S11-15", if (n_c_set_bad == 0L) "PASS" else "FAIL",
        sprintf(paste0("sheet C genes equal the non-mitochondrial Post-IP vs",
                       " Pre-IP survivors of sheet B in %d of 2 region sets"),
                2L - n_c_set_bad))
    # Row order of the published table: per region set, the genes that
    # passed the control test by FDR_M4 of sheet B (smallest first), then
    # the defining genes in the order of the Step 9 file.
    kB <- which(B$comparison == "PostIP_vs_PreIP")
    fdr_c <- B$FDR_M4[kB][match(paste(C$gene, C$region_set),
                                paste(B$gene[kB], B$region_set[kB]))]
    fdr_c[is_def] <- 0
    if (anyNA(fdr_c))
        stop("sheet C: FDR_M4 missing in sheet B for a gene that passed",
             call. = FALSE)
    C <- C[order(match(C$region_set, DATASETS), is_def, fdr_c,
                 seq_len(nrow(C))), ]
    rownames(C) <- NULL
    cat("  Sheet C rows ordered by region set, then FDR_M4 within the passed genes\n")
    print(utils::head(C, 40), row.names = FALSE)

    # =========================================================================
    # 4b. Sheets E and F: survivors in both region sets; neuronal panels
    # =========================================================================
    hdr("4b. SHEETS E AND F (survivors in both region sets; neuronal panels)")
    MAIN_ARMS <- c("M4", "PVM_3", "PVM_1")
    SHORT_NEURO <- c("SNAP25", "NRGN", "SLC17A7", "GAD1", "GAD2", "STMN2",
                     "SYP", "VSNL1", "NEFL", "RBFOX3", "ENC1")
    lookup <- function(genes, tgs) {
        out <- list()
        for (lab in names(ARMS)) for (ds in DATASETS) {
            o <- RDS[[paste(lab, ds)]]; if (is.null(o)) next
            for (tg in tgs) {
                d <- o[["res"]][[tg]][["donor"]]; k <- match(genes, d$gene)
                out[[length(out) + 1L]] <- data.frame(
                    gene = genes, arm = lab, region_set = ds, contrast = tg,
                    n_paired_donors = num(d$n_donors[1]),
                    log2FC = round(num(d$avg_log2FC[k]), 3),
                    BH = signif(num(d$p_val_adj_BH[k]), 3),
                    significant = as_lgl(d$significant[k]))
            }
        }
        do.call(rbind, out)
    }

    # ---- Sheet E -------------------------------------------------------------
    sv <- B[B$comparison == "PostIP_vs_PreIP" & B$passes_control_test &
            !B$mitochondrial, c("gene", "region_set", "direction_on_PVM_3")]
    names(sv) <- c("gene", "passed_in", "direction_on_PVM_3")
    if (anyDuplicated(sv$gene))
        stop("a survivor passed in both region sets; sheet E assumes one",
             call. = FALSE)
    sv$other_region <- ifelse(sv$passed_in == "NeoC", "NonNeoC", "NeoC")
    LE <- lookup(sv$gene, "PostIP_vs_PreIP")
    E <- sv
    for (lab in names(ARMS)) for (side in c("passed", "other")) {
        rs <- if (side == "passed") E$passed_in else E$other_region
        k  <- match(paste(E$gene, lab, rs), paste(LE$gene, LE$arm, LE$region_set))
        E[[paste0("log2FC_", lab, "_", side)]] <- LE$log2FC[k]
        if (side == "other") {
            E[[paste0("BH_", lab, "_other")]] <- LE$BH[k]
            E[[paste0("significant_", lab, "_other")]] <- LE$significant[k]
        }
    }
    sg <- ifelse(E$direction_on_PVM_3 == "up", 1, -1)
    same <- sapply(MAIN_ARMS, function(a)
        sign(E[[paste0("log2FC_", a, "_other")]]) == sg)
    same[is.na(same)] <- FALSE
    E$main_arms_same_direction_other <- as.integer(rowSums(same))
    E$min_abs_log2FC_main_other <- apply(abs(as.matrix(
        E[, paste0("log2FC_", MAIN_ARMS, "_other")])), 1, min)
    E$HKG8_opposite_other <- !is.na(E$log2FC_HKG8_other) &
        sign(E$log2FC_HKG8_other) == -sg
    E$PVM_3_BH_below_0.05_other <- !is.na(E$BH_PVM_3_other) &
        E$BH_PVM_3_other < ALPHA
    rownames(E) <- NULL
    # Cross-check: the passed-region Module log2FC equals sheet C's value.
    n_e <- 0L; n_e_ok <- 0L
    for (lab in names(ARMS)) {
        cc <- paste0("log2FC_ModulePostVsPre_", lab)
        if (!cc %in% names(C)) next
        k <- match(paste(E$gene, E$passed_in), paste(C$gene, C$region_set))
        v <- num(C[[cc]][k]); e <- E[[paste0("log2FC_", lab, "_passed")]]
        have <- is.finite(v) & is.finite(e)
        n_e <- n_e + sum(have); n_e_ok <- n_e_ok + sum(abs(v[have] - e[have]) <= TOL)
    }
    chk("S11-16", if (n_e > 0L && n_e_ok == n_e) "PASS" else "FAIL",
        sprintf(paste0("sheet E passed-region log2FC equals sheet C on %d of %d",
                       " gene x arm values"), n_e_ok, n_e))
    cat(sprintf(paste0("  Same direction on M4, PVM_3 and PVM_1 in the other",
                       " region set: %d of %d (not: %s); gene x main-arm values:",
                       " %d of %d; HKG8 opposite: %d of %d; PVM_3 BH < 0.05 in",
                       " the other region set: %d of %d\n"),
                sum(E$main_arms_same_direction_other == 3L), nrow(E),
                paste(E$gene[E$main_arms_same_direction_other < 3L], collapse = ", "),
                sum(same), length(same), sum(E$HKG8_opposite_other), nrow(E),
                sum(E$PVM_3_BH_below_0.05_other), nrow(E)))

    # ---- Sheets F ------------------------------------------------------------
    tgF <- c("PostIP_vs_PreIP", "NMPost_vs_NMPre")
    FG <- lookup(c(SHORT_NEURO, sv$gene), tgF)
    FG$panel <- ifelse(FG$gene %in% SHORT_NEURO, "short neuron-specific",
                       "control-test survivor")
    FG <- FG[, c("panel", names(FG)[names(FG) != "panel"])]
    FS <- do.call(rbind, lapply(split(FG, list(FG$panel, FG$arm, FG$region_set,
                                               FG$contrast), drop = TRUE),
        function(x) data.frame(
            contrast = x$contrast[1], region_set = x$region_set[1],
            arm = x$arm[1], panel = x$panel[1], n_genes = nrow(x),
            n_detected = sum(!is.na(x$log2FC)),
            median_log2FC = if (all(is.na(x$log2FC))) NA_real_ else
                round(stats::median(x$log2FC, na.rm = TRUE), 3),
            n_up = sum(x$log2FC > 0, na.rm = TRUE),
            n_significant_up = sum(x$significant & x$log2FC > 0, na.rm = TRUE),
            n_significant_down = sum(x$significant & x$log2FC < 0, na.rm = TRUE),
            n_BH_below_0.05_up = sum(x$BH < ALPHA & x$log2FC > 0, na.rm = TRUE))))
    FS <- FS[order(match(FS$contrast, tgF), FS$region_set,
                   match(FS$arm, names(ARMS)), FS$panel), ]
    rownames(FS) <- NULL
    FG <- FG[order(match(FG$contrast, tgF), FG$region_set,
                   match(FG$arm, names(ARMS)), FG$panel), ]
    rownames(FG) <- NULL
    # Cross-check: survivor rows of F equal sheet E (Module Post-IP vs Pre-IP).
    fk <- FG$panel == "control-test survivor" & FG$contrast == "PostIP_vs_PreIP"
    ek <- match(paste(FG$gene[fk], FG$arm[fk], FG$region_set[fk]),
                paste(LE$gene, LE$arm, LE$region_set))
    n_f_bad <- sum(abs(FG$log2FC[fk] - LE$log2FC[ek]) > TOL, na.rm = TRUE)
    chk("S11-17", if (n_f_bad == 0L) "PASS" else "FAIL",
        sprintf("sheet F survivor rows equal sheet E on %d of %d values",
                sum(fk) - n_f_bad, sum(fk)))
    chk("S11-18", "NOTE",
        sprintf(paste0("short neuron-specific genes not in a donor frame (NA):",
                       " %d of %d gene x arm x region x contrast values"),
                sum(is.na(FG$log2FC[FG$panel == "short neuron-specific"])),
                sum(FG$panel == "short neuron-specific")))
    print(FS[FS$contrast == "PostIP_vs_PreIP", ], row.names = FALSE)

    # =========================================================================
    # 4c. Sheets G and H: UMIs per nucleus and panel UMIs
    # =========================================================================
    hdr("4c. SHEETS G AND H (UMIs per nucleus; panel UMIs per nucleus)")
    # Cell groups: Part12_rerun/04_Labels/04_cell_labels.csv.gz. Module
    # cells carry ip_status; non-Module cells are split at IP_y on the
    # labels' composite (composite >= IP_y is Post-IP). IP_y is not in the labels file, so it is bracketed
    # by the highest Module Pre-IP and the lowest Module Post-IP composite;
    # non-Module cells inside that gap are assigned so that the non-Module
    # Post-IP and Pre-IP counts equal Part 1's (highest composites to Post-IP).
    # UMIs: Part 8 raw counts, R_export_CP10K/<ARM_TAG>/<region set>/counts/
    # raw_selected (matrix.mtx, barcodes, features).
    ARM_TAGS <- c(M4 = "Mic13_PVM3", PVM_3 = "PVM3_all", PVM_1 = "PVM1_all",
                  HKG8 = "HKG8_matched")
    MIN_CELLS <- 3L     # per donor and side, as Part 1's pairing rule
    MG_OWN <- c("CSF1R", "P2RY12", "CX3CR1", "DOCK8", "ARHGAP15", "INPP5D")
    PANELS <- list(
        long_neuronal  = sv$gene[sv$direction_on_PVM_3 == "up"],
        short_neuronal = SHORT_NEURO,
        microglial_own = MG_OWN)
    arm_std <- function(x) {
        u <- toupper(gsub("[^A-Za-z0-9]", "", as.character(x)))
        out <- rep(NA_character_, length(u))
        out[u %in% c("M4", "MIC13PVM3", "MIC13")] <- "M4"
        out[u %in% c("PVM3", "PVM3ALL")] <- "PVM_3"
        out[u %in% c("PVM1", "PVM1ALL")] <- "PVM_1"
        out[u %in% c("HKG8", "HKG8MATCHED")] <- "HKG8"
        out
    }
    rs_std <- function(x) {
        u <- toupper(gsub("[^A-Za-z]", "", as.character(x)))
        ifelse(u %in% c("NEOC", "NEOCORTEX"), "NeoC",
        ifelse(u %in% c("NONNEOC", "NONNEOCORTEX"), "NonNeoC", NA_character_))
    }

    # ---- Labels ---------------------------------------------------------------
    # 04_cell_labels.csv.gz (Part12_04_Labels.r): one row per cell and arm,
    # with barcode, arm, dataset, donor_id, is_module, ip_status (Post or Pre
    # on Module cells, empty otherwise) and composite.
    f_lab <- F_LABELS
    if (!file.exists(f_lab)) stop("not found: ", f_lab, call. = FALSE)
    X <- utils::read.csv(f_lab, check.names = FALSE,
                         colClasses = c(barcode = "character",
                                        donor_id = "character"))
    need <- c("barcode", "arm", "dataset", "donor_id", "is_module",
              "ip_status", "composite")
    if (!all(need %in% names(X)))
        stop("labels file lacks columns: ",
             paste(setdiff(need, names(X)), collapse = ", "), call. = FALSE)
    cat(sprintf("  Labels: %s (%s), %d cells\n", f_lab,
                format(file.info(f_lab)$mtime, "%Y-%m-%d %H:%M"), nrow(X)))
    L <- data.frame(barcode = as.character(X$barcode), arm = arm_std(X$arm),
                    rs = rs_std(X$dataset), donor = trimws(as.character(X$donor_id)),
                    mod = as_lgl(X$is_module),
                    post = ifelse(tolower(X$ip_status) == "post", TRUE,
                           ifelse(tolower(X$ip_status) == "pre", FALSE, NA)),
                    comp = num(X$composite))
    rm(X)
    if (anyNA(L$arm) || anyNA(L$rs))
        stop("labels: arm or dataset values not recognised", call. = FALSE)

    # ---- Non-Module IP side; Module and non-Module counts against Part 1 ----
    SPL <- list(); n_cnt19_bad <- 0L
    for (lab in names(ARMS)) for (ds in DATASETS) {
        key <- paste(lab, ds); o <- RDS[[key]]; if (is.null(o)) next
        pp <- o[["res"]][["PostIP_vs_PreIP"]][["donor"]]
        nm <- o[["res"]][["NMPost_vs_NMPre"]][["donor"]]
        k  <- which(L$arm == lab & L$rs == ds)
        km <- k[L$mod[k]]; kn <- k[!L$mod[k]]
        lo <- max(L$comp[km][L$post[km] %in% FALSE], na.rm = TRUE)
        hi <- min(L$comp[km][L$post[km] %in% TRUE], na.rm = TRUE)
        cn <- L$comp[kn]
        side <- ifelse(cn >= hi, TRUE, ifelse(cn <= lo, FALSE, NA))
        n_gap <- sum(is.na(side))
        need_post <- num(nm$n_fg[1]) - sum(side %in% TRUE)
        ok_gap <- is.finite(need_post) && need_post >= 0 && need_post <= n_gap &&
                  lo < hi && !anyNA(cn)
        if (ok_gap && n_gap > 0L) {
            g <- which(is.na(side))
            g <- g[order(cn[g], decreasing = TRUE)]
            side[g] <- seq_along(g) <= need_post
        }
        L$post[kn] <- side
        r <- data.frame(
            arm = lab, region_set = ds,
            module_pre_max_composite = round(lo, 6),
            module_post_min_composite = round(hi, 6),
            module_post = sum(L$post[km] %in% TRUE),
            module_pre = sum(L$post[km] %in% FALSE),
            part1_module_post = num(pp$n_fg[1]), part1_module_pre = num(pp$n_bg[1]),
            nonmodule_cells = length(kn),
            nonmodule_in_gap = n_gap,
            nonmodule_gap_to_post = if (ok_gap) need_post else NA_real_,
            nonmodule_post = sum(side %in% TRUE),
            nonmodule_pre = sum(side %in% FALSE),
            part1_nonmodule_post = num(nm$n_fg[1]),
            part1_nonmodule_pre = num(nm$n_bg[1]))
        r$agrees_with_part1 <- ok_gap &&
            r$module_post == r$part1_module_post && r$module_pre == r$part1_module_pre &&
            r$nonmodule_post == r$part1_nonmodule_post &&
            r$nonmodule_pre == r$part1_nonmodule_pre
        if (!r$agrees_with_part1) n_cnt19_bad <- n_cnt19_bad + 1L
        SPL[[key]] <- r
    }
    GS <- do.call(rbind, SPL); rownames(GS) <- NULL
    chk("S11-19", if (n_cnt19_bad == 0L && nrow(GS) == 8L) "PASS" else "FAIL",
        sprintf(paste0("labels give Part 1's Module and non-Module Post-IP and",
                       " Pre-IP counts on %d of %d arm/region sets (%d",
                       " non-Module cells in the IP_y gap)"),
                nrow(GS) - n_cnt19_bad, nrow(GS), sum(GS$nonmodule_in_gap)))
    print(GS, row.names = FALSE)

    # ---- UMIs from the raw counts ----------------------------------------------
    if (!requireNamespace("Matrix", quietly = TRUE))
        stop("package Matrix is required for sheets G and H", call. = FALSE)
    umi_from_raw <- function(d) {
        mtx <- list.files(d, pattern = "^matrix\\.mtx(\\.gz)?$", full.names = TRUE)
        bcf <- list.files(d, pattern = "barcode", full.names = TRUE)
        ftf <- list.files(d, pattern = "feature|gene", full.names = TRUE)
        if (length(mtx) != 1L || length(bcf) != 1L || length(ftf) != 1L)
            stop(sprintf("%s needs one matrix.mtx, barcodes and features file", d),
                 call. = FALSE)
        M  <- Matrix::readMM(mtx)
        bc <- utils::read.delim(bcf, header = FALSE)[[1]]
        ft <- utils::read.delim(ftf, header = FALSE)
        nmt <- vapply(ft, function(z) sum(grepl("^MT-", z)), numeric(1))
        sym <- ft[[if (any(nmt > 0)) which.max(nmt) else ncol(ft)]]
        if (ncol(M) == length(bc) && nrow(M) == length(sym)) M <- Matrix::t(M)
        if (nrow(M) != length(bc) || ncol(M) != length(sym))
            stop(sprintf("matrix in %s does not match its barcodes and features", d),
                 call. = FALSE)
        M  <- methods::as(M, "CsparseMatrix")
        mt <- grepl("^MT-", sym)
        out <- data.frame(barcode = bc, total = as.numeric(Matrix::rowSums(M)),
                          non_mito = as.numeric(Matrix::rowSums(M[, !mt, drop = FALSE])))
        miss <- character(0)
        for (p in names(PANELS)) {
            kk <- sym %in% PANELS[[p]]
            miss <- c(miss, setdiff(PANELS[[p]], sym))
            out[[paste0("p_", p)]] <- as.numeric(Matrix::rowSums(M[, kk, drop = FALSE]))
        }
        attr(out, "n_mt") <- sum(mt); attr(out, "missing") <- unique(miss)
        out
    }
    L$total <- NA_real_; L$non_mito <- NA_real_
    for (p in names(PANELS)) L[[paste0("p_", p)]] <- NA_real_
    n_match_bad <- 0L; panel_missing <- character(0)
    root <- file.path(CEP_HOME, "R_export_CP10K")
    for (lab in names(ARMS)) {
        adir <- file.path(root, ARM_TAGS[[lab]])
        raws <- list.dirs(adir, recursive = TRUE, full.names = TRUE)
        raws <- raws[basename(raws) == "raw_selected" &
                     basename(dirname(raws)) == "counts"]
        for (ds in DATASETS) {
            d <- raws[rs_std(basename(dirname(dirname(raws)))) %in% ds]
            k <- which(L$arm == lab & L$rs == ds)
            if (length(d) != 1L) {
                n_match_bad <- n_match_bad + 1L
                cat(sprintf("  %-5s %-7s %d raw_selected folders under %s\n",
                            lab, ds, length(d), adir)); next
            }
            U <- umi_from_raw(d)
            m <- match(L$barcode[k], U$barcode)
            for (cc in setdiff(names(U), "barcode")) L[[cc]][k] <- U[[cc]][m]
            if (anyNA(m)) n_match_bad <- n_match_bad + 1L
            panel_missing <- unique(c(panel_missing, attr(U, "missing")))
            cat(sprintf("  %-5s %-7s %s: %d barcodes, %d MT- genes, labels matched %d of %d\n",
                        lab, ds, basename(dirname(dirname(d))), nrow(U),
                        attr(U, "n_mt"), sum(!is.na(m)), length(k)))
        }
    }
    chk("S11-20", if (n_match_bad == 0L) "PASS" else "FAIL",
        sprintf(paste0("every labelled cell found in the Part 8 raw counts on",
                       " %d of 8 arm/region sets"), 8L - n_match_bad))
    chk("S11-21", if (!length(panel_missing)) "PASS" else "WARN",
        if (!length(panel_missing)) "all panel genes present in the raw counts"
        else paste("panel genes absent from the raw counts:",
                   paste(panel_missing, collapse = ", ")))

    CE <- L[!is.na(L$post) & is.finite(L$total) & L$total > 0 &
            is.finite(L$non_mito) & L$non_mito > 0, ]
    cat(sprintf("  Cells used: %d of %d\n", nrow(CE), nrow(L)))
    sr <- function(v) if (length(v) >= 5L) suppressWarnings(stats::wilcox.test(
        v, exact = FALSE, correct = TRUE)$p.value) else NA_real_
    cells_lab <- function(m) if (m) "Module" else "non-Module"

    # ---- Sheet G: UMIs per nucleus -----------------------------------------
    gd <- list(); gs <- list()
    for (ms in c("total", "non_mito")) for (m in c(TRUE, FALSE))
    for (lab in names(ARMS)) for (ds in DATASETS) {
        s <- CE[CE$arm == lab & CE$rs == ds & CE$mod == m, ]
        if (!nrow(s)) next
        lu <- log2(s[[ms]])
        rows <- lapply(split(seq_len(nrow(s)), s$donor), function(i) {
            p1 <- s$post[i]
            if (sum(p1) < MIN_CELLS || sum(!p1) < MIN_CELLS) return(NULL)
            data.frame(measure = ms, cells = cells_lab(m), arm = lab,
                       region_set = ds, donor = s$donor[i[1]],
                       n_post = sum(p1), n_pre = sum(!p1),
                       median_log2_umi_post = stats::median(lu[i][p1]),
                       median_log2_umi_pre  = stats::median(lu[i][!p1]))
        })
        dd <- do.call(rbind, rows); if (is.null(dd)) next
        dd$log2_ratio <- dd$median_log2_umi_post - dd$median_log2_umi_pre
        gd[[length(gd) + 1L]] <- dd
        gs[[length(gs) + 1L]] <- data.frame(
            measure = ms, cells = cells_lab(m), arm = lab, region_set = ds,
            n_paired_donors = nrow(dd),
            median_log2_ratio = round(stats::median(dd$log2_ratio), 3),
            fold_post_vs_pre = round(2^stats::median(dd$log2_ratio), 3),
            n_donors_post_lower = sum(dd$log2_ratio < 0),
            p_signed_rank = signif(sr(dd$log2_ratio), 3))
    }
    GD <- do.call(rbind, gd); G <- do.call(rbind, gs)
    G$BH <- NA_real_
    for (ms in unique(G$measure)) for (cc in unique(G$cells)) {
        k <- G$measure == ms & G$cells == cc
        G$BH[k] <- signif(stats::p.adjust(G$p_signed_rank[k], "BH"), 3)
    }
    rownames(G) <- NULL; rownames(GD) <- NULL
    GD[, c("median_log2_umi_post", "median_log2_umi_pre", "log2_ratio")] <-
        round(GD[, c("median_log2_umi_post", "median_log2_umi_pre", "log2_ratio")], 4)
    print(G[G$measure == "non_mito", ], row.names = FALSE)

    # ---- Sheet H: panel UMIs per nucleus, amount and share -------------------
    hd <- list(); hs <- list()
    for (p in names(PANELS)) for (m in c(TRUE, FALSE))
    for (lab in names(ARMS)) for (ds in DATASETS) {
        pc <- paste0("p_", p)
        s <- CE[CE$arm == lab & CE$rs == ds & CE$mod == m & is.finite(CE[[pc]]), ]
        if (!nrow(s)) next
        n_zero <- 0L
        rows <- lapply(split(seq_len(nrow(s)), s$donor), function(i) {
            p1 <- s$post[i]
            if (sum(p1) < MIN_CELLS || sum(!p1) < MIN_CELLS) return(NULL)
            x <- s[[pc]][i]; o <- s$non_mito[i]
            a1 <- mean(x[p1]); a0 <- mean(x[!p1])
            if (a1 == 0 || a0 == 0) { n_zero <<- n_zero + 1L; return(NULL) }
            data.frame(panel = p, cells = cells_lab(m), arm = lab,
                       region_set = ds, donor = s$donor[i[1]],
                       n_post = sum(p1), n_pre = sum(!p1),
                       panel_umi_per_nucleus_post = a1,
                       panel_umi_per_nucleus_pre = a0,
                       share_post = sum(x[p1]) / sum(o[p1]),
                       share_pre = sum(x[!p1]) / sum(o[!p1]),
                       log2_own_rna = log2(mean(o[p1]) / mean(o[!p1])),
                       log2_amount = log2(a1 / a0))
        })
        dd <- do.call(rbind, rows); if (is.null(dd)) next
        dd$log2_share <- log2(dd$share_post / dd$share_pre)
        hd[[length(hd) + 1L]] <- dd
        hs[[length(hs) + 1L]] <- data.frame(
            panel = p, cells = cells_lab(m), arm = lab, region_set = ds,
            n_paired_donors = nrow(dd), n_donors_zero = n_zero,
            panel_umi_per_nucleus_pre = round(stats::median(dd$panel_umi_per_nucleus_pre), 2),
            fold_own_rna = round(2^stats::median(dd$log2_own_rna), 3),
            fold_amount = round(2^stats::median(dd$log2_amount), 3),
            p_amount = signif(sr(dd$log2_amount), 3),
            fold_share = round(2^stats::median(dd$log2_share), 3),
            p_share = signif(sr(dd$log2_share), 3))
    }
    HD <- do.call(rbind, hd); H <- do.call(rbind, hs)
    H$BH_amount <- NA_real_; H$BH_share <- NA_real_
    for (p in unique(H$panel)) for (cc in unique(H$cells)) {
        k <- H$panel == p & H$cells == cc
        H$BH_amount[k] <- signif(stats::p.adjust(H$p_amount[k], "BH"), 3)
        H$BH_share[k]  <- signif(stats::p.adjust(H$p_share[k], "BH"), 3)
    }
    H <- H[, c("panel", "cells", "arm", "region_set", "n_paired_donors",
               "n_donors_zero", "panel_umi_per_nucleus_pre", "fold_own_rna",
               "fold_amount", "p_amount", "BH_amount", "fold_share", "p_share",
               "BH_share")]
    rownames(H) <- NULL; rownames(HD) <- NULL
    hcols <- c("panel_umi_per_nucleus_post", "panel_umi_per_nucleus_pre",
               "share_post", "share_pre", "log2_own_rna", "log2_amount",
               "log2_share")
    HD[, hcols] <- signif(HD[, hcols], 5)
    # Calibration: the microglial panel should follow the cell's own RNA.
    hc <- H[H$panel == "microglial_own" & H$cells == "Module" &
            H$arm %in% MAIN_ARMS, ]
    chk("S11-22", "NOTE",
        sprintf(paste0("microglial panel, Module cells, main arms: fold_amount",
                       " %.2f to %.2f against fold_own_rna %.2f to %.2f"),
                min(hc$fold_amount), max(hc$fold_amount),
                min(hc$fold_own_rna), max(hc$fold_own_rna)))
    print(H[H$cells == "Module", ], row.names = FALSE)

    # =========================================================================
    # 5. Sheet D: control validity (Step 8)
    # =========================================================================
    hdr("5. SHEET D (agreement of Module assignment between arms)")
    f8 <- file.path(S79_DIR, "Part10g_Step8_control_overlap.csv")
    if (!file.exists(f8)) stop("not found: ", f8, call. = FALSE)
    D <- utils::read.csv(f8, check.names = FALSE)
    names(D)[names(D) == "dataset"] <- "region_set"
    ok_d <- all(num(D$both_module) + num(D$a_only_module) +
                num(D$b_only_module) + num(D$neither_module) == num(D$n_shared))
    chk("S11-08", if (ok_d) "PASS" else "FAIL",
        "both + a only + b only + neither = shared cells on every row")
    print(D, row.names = FALSE)

    # =========================================================================
    # 6. Contents
    # =========================================================================
    # Short Contents sheet in the layout of the other supplementary tables:
    # title, one note, a blank row, then one row per sheet (sheet, rows,
    # description). The title is written as the sheet's header row, so it is
    # the bold row; writexl cannot bold a later row or set column widths.
    DATA <- list(A_Contrasts = A, B_Control_test = B, B_Summary = BS,
                 C_Axis_check = C, D_Control_validity = D,
                 E_Survivors_regions = E, F_Neuronal_summary = FS,
                 F_Neuronal_genes = FG, G_UMI_summary = G, G_UMI_donors = GD,
                 G_IP_split = GS, H_Panel_summary = H, H_Panel_donors = HD)
    DESC <- c(
        A_Contrasts = "Four comparisons per arm and region set: cells per group, paired donors, genes tested, significant genes by direction and mitochondrial flag, the smallest BH, and the smallest p and BH attainable with that donor count.",
        B_Control_test = "Genes significant on PVM_3, with the test of their per-donor log2FC on M4 and on PVM_1 against HKG8; passes_control_test is TRUE when both arms differ from HKG8 at two-sided FDR <0.05 in the direction seen on PVM_3.",
        B_Summary = "PVM_3 DEGs per comparison and region set, and those that pass against HKG8 on M4, on PVM_1 and on both, mitochondrial or not.",
        C_Axis_check = "The 14 non-mitochondrial genes that passed on post-IP against pre-IP, with their log2FC in Module cells and in non-Module cells crossing the same IP; the nine defining genes are listed for scale.",
        D_Control_validity = "Cohen's kappa for Module membership between arms, on shared cells.",
        E_Survivors_regions = "The 14 genes in both region sets on all four arms, and whether each changed in the same direction in the region set where it did not pass.",
        F_Neuronal_summary = "Eleven short neuron-specific genes beside the 14 genes, in Module and in non-Module cells crossing the IP. The pseudocount compresses small log2FC values, so only signs are comparable between panels.",
        F_Neuronal_genes = "Gene-level values behind F_Neuronal_summary.",
        G_UMI_summary = "Total and non-mitochondrial UMIs per nucleus, post-IP against pre-IP within donor, in Module and in non-Module cells.",
        G_UMI_donors = "Donor values behind G_UMI_summary.",
        G_IP_split = "Assignment of cells to post-IP and pre-IP for sheets G and H, checked against the Part 1 cell counts.",
        H_Panel_summary = "UMIs per nucleus (amount) and share of non-mitochondrial UMIs for the long neuronal genes, the short neuron-specific genes and six microglial genes, post-IP against pre-IP within donor.",
        H_Panel_donors = "Donor values behind H_Panel_summary.")
    stopifnot(identical(names(DESC), names(DATA)))
    TITLE <- paste("Supplementary Table S11. Differential expression between",
                   "Module and non-Module cells, and between post-IP and",
                   "pre-IP cells.")
    NOTE <- paste("Paired pseudobulk per donor (Wilcoxon signed-rank test;",
                  "BH across genes within each comparison); significant:",
                  "BH <0.05 and |log2FC| >= 0.25. M4, PVM_3 and PVM_1 carry",
                  "the MT-CO3 Module; HKG8 is the housekeeping control arm.",
                  "NeoC and NonNeoC were analysed separately. Mitochondrial",
                  "genes (MT-, MTRNR2L, TSPO) are flagged because the Module",
                  "is defined on them.")
    CT <- data.frame(
        c1 = c(NOTE, "", "Sheet", names(DATA)),
        c2 = c("", "", "Rows", vapply(DATA, function(x) as.character(nrow(x)),
                                      character(1))),
        c3 = c("", "", "Description", unname(DESC)))
    names(CT) <- c(TITLE, "", "")
    nonascii <- grepl("[^ -~]", c(TITLE, unlist(CT)))
    chk("S11-09", if (!any(nonascii)) "PASS" else "WARN",
        sprintf("Contents text is plain ASCII (%d non-ASCII entries)", sum(nonascii)))

    # =========================================================================
    # 7. Write
    # =========================================================================
    hdr("7. WRITE")
    # Remove the workbook of an earlier run first, so a blocked or failed
    # build cannot leave a stale file at the output path.
    if (file.exists(S11_XLSX)) {
        unlink(S11_XLSX)
        if (file.exists(S11_XLSX))
            stop("cannot remove the earlier ", basename(S11_XLSX),
                 "; close it in Excel and run again", call. = FALSE)
        cat("  Removed the workbook of an earlier run.\n")
    }
    LG <- do.call(rbind, LEDGER)
    n_fail <- sum(LG$status == "FAIL")
    if (n_fail > 0L) {
        cat("  ", n_fail, " FAIL(s): workbook NOT written.\n", sep = "")
        return(invisible(LG))
    }
    # writexl gives bold headers and no fills or borders. openxlsx is not
    # used because its files carry drawing relationships with no drawing
    # part, which openpyxl cannot open.
    sheets <- c(list(Contents = CT), DATA)
    writexl::write_xlsx(sheets, S11_XLSX, format_headers = TRUE)

    # Read back and compare row counts
    back <- vapply(names(sheets), function(nm)
        nrow(readxl::read_excel(S11_XLSX, sheet = nm)), integer(1))
    want <- vapply(sheets, nrow, integer(1))
    ok_back <- identical(unname(back), unname(want))
    chk("S11-10", if (ok_back) "PASS" else "FAIL",
        paste(sprintf("%s %d", names(back), back), collapse = " | "))
    if (!ok_back) {
        unlink(S11_XLSX)
        cat("  Read-back row counts differ: workbook removed.\n")
    } else cat("  Written:", S11_XLSX, "\n")

    LG <- do.call(rbind, LEDGER)
    hdr("SUMMARY")
    cat(sprintf("  PASS %d | WARN %d | FAIL %d | NOTE %d\n",
                sum(LG$status == "PASS"), sum(LG$status == "WARN"),
                sum(LG$status == "FAIL"), sum(LG$status == "NOTE")))
    invisible(LG)
}

# PASS only when the build ran to the end with no FAIL and the workbook
# was written and read back.
.lg <- tryCatch(s11_build(),
                error = function(e) {
                    cat("\nSTOPPED:", conditionMessage(e), "\n")
                    NULL
                })
.pass <- is.data.frame(.lg) && !any(.lg$status == "FAIL") &&
         any(.lg$id == "S11-10" & .lg$status == "PASS") && file.exists(S11_XLSX)
cat(sprintf("\nVERDICT: %s\n", if (.pass) "PASS" else "FAIL"))
while (sink.number() > 0) sink()
try(close(.log_con), silent = TRUE)
options(.old_opt)
cat("Log written to:", .log_txt, "\n")
