##############################################################################
# build_SuppTable_S10_step0_Fig8_cells.r
#
# Writes the eight per-arm Figure 8 cell files that build_SuppTable_S10.py
# reads for sheet A of Supplementary Table S10: the 2D UMAP coordinates, CMRF
# class and refit group of every GOI-positive cell (from the Validation 2
# results RDS) and its copula score coordinates u_SK and v_SK (from the
# Validation 1 results RDS). No value is computed; the columns are copied and
# joined on the barcode. Run it once before build_SuppTable_S10.py.
#
# Usage
#   Rscript build_SuppTable_S10_step0_Fig8_cells.r            (root: $SEAAD_ROOT, else the current folder)
#   SEAAD_ROOT=/data Rscript build_SuppTable_S10_step0_Fig8_cells.r
#   In R: CEP_HOME <- "/data"; source("build_SuppTable_S10_step0_Fig8_cells.r")
#
# Inputs (read only), written by the deposited pipeline:
#   Part10e_Validation2/Part10e_Validation2_<folder>/Part10e_CP10K_<tag>_V2_<ds>_results.rds
#       written by Part10e_CP10K_Validation2_UMAP_code.r; element df_umap2d
#       (barcode, UMAP_1, UMAP_2, cell_class, group_ref, ...)
#   Part10e_Validation1/Part10e_Validation1_<folder>/Part10e_CP10K_<tag>_V1_<ds>_results.rds
#       written by Part10e_CP10K_Validation1_code.r; elements copula_u_sk and
#       copula_v_sk (named vectors, barcode -> u_SK, v_SK), the same fields
#       Validation 2 reads for its copula projection
#   <tag>/<folder>: M4/Mic13_PVM3, PVM_3/PVM_3, PVM_1/PVM_1, HKG8/HKG8;
#   <ds>: NeoC, NonNeoC
#
# Output, in <root>/Manuscript/Supp_Tables/build/ (created if missing):
#   Part10e_CP10K_<tag>_V2_<ds>_Fig8_cells.csv, columns barcode, UMAP_1,
#   UMAP_2, group_ref, cell_class, u_SK, v_SK in the row order of df_umap2d
#   (existing files are replaced)
#   build_SuppTable_S10_step0_Fig8_cells_log.txt, ending with VERDICT: PASS or
#   VERDICT: FAIL; build_SuppTable_S10.py requires PASS
#
# Checks per arm and region set:
#   F8-01  the V2 and V1 results RDS exist and hold df_umap2d, copula_u_sk
#          and copula_v_sk
#   F8-02  every output column can be produced
#   F8-03  barcodes are unique and every cell has u_SK and v_SK
##############################################################################

if (!exists("CEP_HOME")) CEP_HOME <- Sys.getenv("SEAAD_ROOT", getwd())
stopifnot(dir.exists(CEP_HOME))
OUT <- file.path(CEP_HOME, "Manuscript", "Supp_Tables", "build")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
.log <- file.path(OUT, "build_SuppTable_S10_step0_Fig8_cells_log.txt")
.con <- file(.log, open = "wt")
sink(.con, split = TRUE)

ARMS <- c(M4 = "Mic13_PVM3", PVM_3 = "PVM_3", PVM_1 = "PVM_1", HKG8 = "HKG8")
DATASETS <- c("NeoC", "NonNeoC")
SK_FIELDS <- c("copula_u_sk", "copula_v_sk")
FIXED_COLS <- c("barcode", "UMAP_1", "UMAP_2", "group_ref", "cell_class", "u_SK", "v_SK")
LEDGER <- character(0)
chk <- function(id, ok, what) {
    s <- sprintf("  [%s] %s %s", if (isTRUE(ok)) "PASS" else "FAIL", id, what)
    cat(s, "\n")
    LEDGER <<- c(LEDGER, if (isTRUE(ok)) "PASS" else "FAIL")
    invisible(isTRUE(ok))
}

cat("build_SuppTable_S10_step0_Fig8_cells.r\n")
cat("Date    :", format(Sys.time()), "\n")
cat("CEP_HOME:", CEP_HOME, "\n")
cat("Output  :", OUT, "\n")
cat("R       :", R.version.string, "\n")

for (tag in names(ARMS)) {
    folder <- ARMS[[tag]]
    v2_dir <- file.path(CEP_HOME, "Part10e_Validation2", paste0("Part10e_Validation2_", folder))
    v1_dir <- file.path(CEP_HOME, "Part10e_Validation1", paste0("Part10e_Validation1_", folder))
    for (ds in DATASETS) {
        cat(sprintf("\n%s %s\n", tag, ds))
        f_v2 <- file.path(v2_dir, sprintf("Part10e_CP10K_%s_V2_%s_results.rds", tag, ds))
        f_v1 <- file.path(v1_dir, sprintf("Part10e_CP10K_%s_V1_%s_results.rds", tag, ds))
        f_out <- file.path(OUT, sprintf("Part10e_CP10K_%s_V2_%s_Fig8_cells.csv", tag, ds))
        cat("    V2 RDS:", f_v2, "\n    V1 RDS:", f_v1, "\n")
        if (!chk("F8-01", file.exists(f_v2) && file.exists(f_v1),
                 sprintf("V2 RDS %s, V1 RDS %s",
                         if (file.exists(f_v2)) "found" else "MISSING",
                         if (file.exists(f_v1)) "found" else "MISSING"))) next

        r2 <- readRDS(f_v2)
        d <- r2$df_umap2d
        r1 <- readRDS(f_v1)
        has_sk <- all(SK_FIELDS %in% names(r1))
        if (!chk("F8-01", !is.null(d) && has_sk,
                 sprintf("df_umap2d in the V2 RDS: %s; %s in the V1 RDS: %s",
                         if (is.null(d)) "no" else "yes", paste(SK_FIELDS, collapse = ", "),
                         if (has_sk) "yes" else "no"))) next
        d <- as.data.frame(d, stringsAsFactors = FALSE)
        d$barcode <- as.character(d$barcode)
        for (k in names(d)) if (is.factor(d[[k]])) d[[k]] <- as.character(d[[k]])
        d$u_SK <- unname(r1[[SK_FIELDS[1]]][d$barcode])
        d$v_SK <- unname(r1[[SK_FIELDS[2]]][d$barcode])
        rm(r1, r2); invisible(gc())

        miss <- setdiff(FIXED_COLS, names(d))
        chk("F8-02", length(miss) == 0L,
            sprintf("columns %s%s", paste(FIXED_COLS, collapse = ", "),
                    if (length(miss)) paste0(" | cannot produce: ", paste(miss, collapse = ", ")) else ""))
        if (length(miss)) next
        new <- d[, FIXED_COLS, drop = FALSE]; rownames(new) <- NULL
        ok <- !anyDuplicated(new$barcode) && !anyNA(new$u_SK) && !anyNA(new$v_SK)
        chk("F8-03", ok, sprintf("%d cells, barcodes unique, u_SK and v_SK present for every cell",
                                 nrow(new)))
        if (!ok) next
        utils::write.csv(new, f_out, row.names = FALSE)
        cat("    written:", f_out, "\n")
    }
}

n_fail <- sum(LEDGER == "FAIL")
n_expected <- 4L * length(ARMS) * length(DATASETS)   # four checks per arm and region set
cat(sprintf("\nChecks: %d PASS, %d FAIL\n", sum(LEDGER == "PASS"), n_fail))
cat("Log:", .log, "\n")
cat(if (n_fail > 0L || length(LEDGER) < n_expected) "VERDICT: FAIL\n" else "VERDICT: PASS\n")
sink(); close(.con)
