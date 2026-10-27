# R11_Pineda_nuclei_table.R
# Collects the R02 cell-calling results of external validation 1 (Pineda et al. 2024, BA9 batch
# 210526) into one nuclei table and one donor table for the Python steps (R12 and R13). This is the
# Pineda counterpart of the nuclei table that R12 writes for Rexach. No correlation or other
# statistic is computed here.
#
# Steps:
#   1. Reads R02_cell_calling/R02_<LIB>.rds for the 63 libraries of the manifest.
#   2. Checks each library: the kept counts match the kept rows of the nuclei table (names and total
#      UMIs), and the gene table is the same in every library.
#   3. Joins the libraries to the S01 donor table (sheets sample_mapping and donor_information):
#      one BA9 library per donor, and the diagnosis in the library name agrees with the S01 Group
#      (C9ALS counts as ALS and C9FTLD as FTLD).
#   4. Writes the nuclei table (all called nuclei, with the gene-floor, doublet, and keep flags),
#      the gene table, the donor table, and a summary workbook.
#
# The counts themselves are not written again. R12 and R13 read the counts of the kept nuclei from
# the Cell Ranger raw .h5 files in counts/, and check them against the total UMIs written here.
#
# Usage (R on Windows):
#   source("C:/Users/Wong/Desktop/SEA-AD/External_Dataset_1/R11_Pineda_nuclei_table.R")
# Inputs, in ROOT: Pineda_BA9_210526_manifest.tsv, R02_cell_calling/R02_<LIB>.rds, and the S01
#   workbook (S01_sample_information.xlsx, from Synapse syn51105515), or its path set in S01.
# Outputs, in ROOT/R11_nuclei/:
#   R11_Pineda_nuclei.csv.gz   one row per called nucleus
#   R11_Pineda_genes.csv.gz    gene_id and symbol, in the order of the Cell Ranger matrices
#   R11_Pineda_donors.csv      one row per library (donor), with the S01 covariates
#   R11_Pineda_summary.xlsx    sheets Libraries, Donors, By_disease, Settings, and Session
#   R11_Pineda_log.txt         console output

# The data folder: PINEDA_ROOT when set (the container sets it to /data), otherwise the Windows folder.
ROOT <- Sys.getenv("PINEDA_ROOT", "C:/Users/Wong/Desktop/SEA-AD/External_Dataset_1")
# S01 workbook. NULL searches ROOT for a .xlsx whose name starts with "S01"; otherwise give the path.
S01 <- NULL

SET <- list(n_libraries = 63L, s01_region = "BA9",
            group_to_disease = c(ALS = "ALS", C9ALS = "ALS", FTLD = "FTLD", C9FTLD = "FTLD", PN = "PN"))

MANIFEST <- file.path(ROOT, "Pineda_BA9_210526_manifest.tsv")
IN02 <- file.path(ROOT, "R02_cell_calling")
OUT <- file.path(ROOT, "R11_nuclei")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
LOG <- file.path(OUT, "R11_Pineda_log.txt")
F_S01 <- if (is.null(S01)) list.files(ROOT, pattern = "^S01.*\\.xlsx$", full.names = TRUE, ignore.case = TRUE) else S01

need <- c("Matrix", "openxlsx")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Install these packages first: ", paste(miss, collapse = ", "))

say <- function(...) cat(sprintf(...), "\n", sep = "")

main <- function() {
t0 <- Sys.time()
say("R11_Pineda_nuclei_table.R started %s", format(t0, "%Y-%m-%d %H:%M:%S"))
say("ROOT %s", ROOT)

# ---- 1. Inputs -----------------------------------------------------------------------------------
if (!file.exists(MANIFEST)) stop("Manifest not found: ", MANIFEST, call. = FALSE)
if (length(F_S01) != 1L || !file.exists(F_S01))
    stop("S01 workbook not found in ROOT (or more than one): ", paste(basename(F_S01), collapse = ", "),
         ". Set S01 at the top of the script.", call. = FALSE)
man <- utils::read.delim(MANIFEST, stringsAsFactors = FALSE, check.names = FALSE)
libs <- man[[1]]
if (length(libs) != SET$n_libraries || anyDuplicated(libs))
    stop(sprintf("The manifest holds %d libraries (%d unique); expected %d.", length(libs),
                 length(unique(libs)), SET$n_libraries), call. = FALSE)
rds <- file.path(IN02, paste0("R02_", libs, ".rds"))
if (!all(file.exists(rds)))
    stop(sprintf("%d of %d R02 files are missing, for example %s. Run R02 on all libraries first.",
                 sum(!file.exists(rds)), length(rds), basename(rds[!file.exists(rds)][1])), call. = FALSE)
say("Manifest: %d libraries | R02 files: all present", length(libs))

# ---- 2. Libraries ----------------------------------------------------------------------------------
nuc_all <- vector("list", length(libs)); lib_rows <- vector("list", length(libs))
genes0 <- NULL
last <- Sys.time()
for (i in seq_along(libs)) {
    r <- readRDS(rds[i])
    nuc <- r$nuclei; K <- r$counts; g <- r$genes
    if (is.null(genes0)) genes0 <- g
    same_genes <- identical(g$id, genes0$id) && identical(g$symbol, genes0$symbol)
    kept <- nuc[nuc$keep, ]
    names_ok <- identical(colnames(K), paste0(libs[i], "_", kept$barcode))
    umis_ok <- isTRUE(all(Matrix::colSums(K) == kept$total_umis))
    lib_ok <- all(nuc$library == libs[i])
    if (!(same_genes && names_ok && umis_ok && lib_ok))
        stop(sprintf("%s: gene table same as the first library %s | kept names match %s | kept total UMIs match %s | library column %s",
                     libs[i], same_genes, names_ok, umis_ok, lib_ok), call. = FALSE)
    nuc$cell_id <- paste0(libs[i], "_", nuc$barcode)
    nuc_all[[i]] <- nuc
    lr <- r$library
    lr$checks <- "kept counts and nuclei table agree; gene table as in the first library"
    lib_rows[[i]] <- lr
    rm(r, K); gc(verbose = FALSE)
    if (as.numeric(difftime(Sys.time(), last, units = "secs")) >= 30 || i == length(libs)) {
        say("  %d of %d libraries read | %.1f min", i, length(libs), as.numeric(difftime(Sys.time(), t0, units = "mins")))
        last <- Sys.time()
    }
}
nuc <- do.call(rbind, nuc_all); rm(nuc_all)
libtab <- do.call(rbind, lib_rows)
if (anyDuplicated(nuc$cell_id)) stop("Duplicated cell IDs across libraries.", call. = FALSE)
i_mt <- grep("^MT-", genes0$symbol)
say("Genes: %d (MT- genes %d) | called nuclei: %d | kept: %d", nrow(genes0), length(i_mt), nrow(nuc), sum(nuc$keep))
if (length(i_mt) != 13L) stop("Expected 13 MT- genes.", call. = FALSE)

# ---- 3. Donors from S01 -------------------------------------------------------------------------------
sm <- openxlsx::read.xlsx(F_S01, sheet = "sample_mapping")
di <- openxlsx::read.xlsx(F_S01, sheet = "donor_information", check.names = FALSE, sep.names = " ")
sm <- sm[sm$Sample_ID %in% libs, ]
dn <- data.frame(library_id = libs, stringsAsFactors = FALSE)
dn$donor_number <- as.character(as.integer(sm$Donor[match(libs, sm$Sample_ID)]))
dn$region_S01 <- sm$Region[match(libs, sm$Sample_ID)]
dn$disease_from_name <- sub("^[0-9]+_([A-Za-z]+)_[0-9]+_.*$", "\\1", libs)
dn$donor_id <- sub("^[0-9]+_([A-Za-z]+_[0-9]+)_.*$", "\\1", libs)
di$donor_number <- as.character(as.integer(di$Donor))
m <- match(dn$donor_number, di$donor_number)
col_of <- function(pattern) {
    k <- grep(pattern, names(di), ignore.case = TRUE)
    if (length(k) != 1L) stop("S01 donor_information: column matching '", pattern, "' found ", length(k), " times.",
                              call. = FALSE)
    names(di)[k]
}
num <- function(v) suppressWarnings(as.numeric(as.character(v)))
dn$group_S01 <- di$Group[m]
dn$disease <- unname(SET$group_to_disease[dn$group_S01])
dn$sex <- di$Sex[m]
dn$race <- di$Race[m]
dn$age_at_death <- num(di[[col_of("^Age of death")]][m])
dn$age_at_death_S01_text <- as.character(di[[col_of("^Age of death")]][m])
dn$pmi <- num(di[[col_of("^PMI")]][m])
dn$pmi_S01_text <- as.character(di[[col_of("^PMI")]][m])
dn$braak_stage <- num(di[[col_of("^Braak")]][m])
dn$thal_phase <- num(di[[col_of("^Thal")]][m])
dn$age_onset_cognitive <- num(di[[col_of("onset cognitive")]][m])
dn$disease_duration <- num(di[[col_of("^Disease duration")]][m])

checks <- c(
    all_mapped = all(!is.na(dn$donor_number)) && all(!is.na(m)),
    one_per_donor = !anyDuplicated(dn$donor_number) && !anyDuplicated(dn$donor_id),
    region_BA9 = all(dn$region_S01 == SET$s01_region),
    known_group = all(dn$group_S01 %in% names(SET$group_to_disease)),
    name_agrees = isTRUE(all(dn$disease == dn$disease_from_name)),
    sex_FM = all(dn$sex %in% c("F", "M")))
say("S01 (%s): %s", basename(F_S01), paste(sprintf("%s %s", names(checks), checks), collapse = " | "))
if (!all(checks)) {
    bad <- dn[is.na(m) | dn$disease != dn$disease_from_name | !(dn$sex %in% c("F", "M")),
              c("library_id", "donor_number", "group_S01", "disease", "disease_from_name", "sex")]
    print(bad, row.names = FALSE)
    stop("The S01 join failed one of the checks above.", call. = FALSE)
}
say("Donors by disease: %s", paste(sprintf("%s %d", names(table(dn$disease)), as.integer(table(dn$disease))),
                                    collapse = ", "))
say("Age at death missing: %s | PMI missing: %s",
    paste(dn$donor_id[is.na(dn$age_at_death)], collapse = ", "), paste(dn$donor_id[is.na(dn$pmi)], collapse = ", "))

nuc$library_id <- nuc$library
nuc$donor_id <- dn$donor_id[match(nuc$library_id, dn$library_id)]
nuc$disease <- dn$disease[match(nuc$library_id, dn$library_id)]
nuc$sex <- dn$sex[match(nuc$library_id, dn$library_id)]
keep_cols <- c("cell_id", "library_id", "donor_id", "disease", "sex", "barcode", "total_umis", "n_genes",
               "mt_umis", "pct_mt", "emptydrops_FDR", "pass_gene_floor", "scDblFinder_score",
               "scDblFinder_class", "keep")
nuc <- nuc[, keep_cols]

# ---- 4. Outputs ----------------------------------------------------------------------------------------
write_gz <- function(d, f) {
    tmp <- paste0(f, ".tmp")
    con <- gzfile(tmp, "w"); utils::write.csv(d, con, row.names = FALSE); close(con)
    if (!file.rename(tmp, f)) stop("could not rename ", tmp, call. = FALSE)
}
write_gz(nuc, file.path(OUT, "R11_Pineda_nuclei.csv.gz"))
write_gz(data.frame(gene_id = genes0$id, symbol = genes0$symbol), file.path(OUT, "R11_Pineda_genes.csv.gz"))
utils::write.csv(dn, file.path(OUT, "R11_Pineda_donors.csv"), row.names = FALSE)

k <- nuc[nuc$keep, ]
by_dis <- do.call(rbind, lapply(split(k, k$disease), function(d) data.frame(
    disease = d$disease[1], donors = length(unique(d$donor_id)), nuclei_kept = nrow(d),
    median_per_donor = stats::median(table(d$donor_id)), median_umis = stats::median(d$total_umis),
    median_genes = stats::median(d$n_genes), median_pct_mt = stats::median(d$pct_mt),
    max_pct_mt = max(d$pct_mt))))
print(by_dis, row.names = FALSE, digits = 4)
libtab$donor_id <- dn$donor_id[match(libtab$library, dn$library_id)]
libtab$disease <- dn$disease[match(libtab$library, dn$library_id)]
settings <- data.frame(setting = c(names(SET), "S01 file"),
                       value = c(vapply(SET, function(v) paste(names(v), v, sep = if (is.null(names(v))) "" else "=",
                                                               collapse = ", "), character(1)), basename(F_S01)),
                       stringsAsFactors = FALSE)
session <- data.frame(line = utils::capture.output(utils::sessionInfo()))
wb <- openxlsx::createWorkbook()
bold <- openxlsx::createStyle(textDecoration = "bold")
for (nm in c("Libraries", "Donors", "By_disease", "Settings", "Session")) {
    d <- switch(nm, Libraries = libtab, Donors = dn, By_disease = by_dis, Settings = settings, Session = session)
    openxlsx::addWorksheet(wb, nm)
    openxlsx::writeData(wb, nm, d, headerStyle = bold)
    openxlsx::freezePane(wb, nm, firstRow = TRUE)
}
openxlsx::saveWorkbook(wb, file.path(OUT, "R11_Pineda_summary.xlsx"), overwrite = TRUE)
say("Written to %s: R11_Pineda_nuclei.csv.gz, R11_Pineda_genes.csv.gz, R11_Pineda_donors.csv, R11_Pineda_summary.xlsx", OUT)
say("Done in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))
invisible(TRUE)
}

zz <- file(LOG, open = "wt")
sink(zz, split = TRUE)
tryCatch(main(), error = function(e) say("STOPPED: %s", conditionMessage(e)),
         finally = { sink(); close(zz) })
