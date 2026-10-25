# R15_Rexach_diagnosis_contrasts.R
# Diagnosis contrasts of external validation 2 (Rexach et al. 2024, insula), Section 6 of the
# readings (external_validation/External2_Rexach_Preset_Readings.md, committed in 195439d, with the
# dated notes on the Section 4 review split and on the Section 6 covariates). Reported as observed,
# without a reading.
#
#   Predictor   Module post-IP% per donor (R14 donor table: Module post-IP cells / all Micro-PVM_3
#               nuclei of the donor x 100), the 33 donors at or above the floor of 5 nuclei.
#   Model       cepip_associate(): OLS on z-scored variables, with age at death, sex (F = 0, M = 1),
#               and PMI as covariates. This is cov3: the SEA-AD cov4 covariates without APOE e4
#               dose, which Table S1 does not give.
#   Contrasts   All three diagnoses combined against controls (min_donors = 20), and AD, bvFTD, and
#               PSP each against controls (min_donors = 15). BH across the four contrasts is
#               reported beside the p values.
#   Arms        Only arms whose Section 4 gate in R14 allows the contrasts. HKG8 is reported beside
#               PVM_3 when its gate allows it, and otherwise as not available, with the reason.
# Also reported: the SEA-AD share of Module post-IP cells above 8% MT (readings Section 1, beside
# Section 3), from the SEA-AD file Part12c_cell_metrics.csv, PVM_3 arm.
#
# Age at death and PMI come from Table S1 of the paper, sheet 5.Metadata_perRNAlibrary, insula
# libraries. CELLxGENE gives a range of ages for 27 of the 37 donors and no PMI.
#
# Guardrails (ledger R15_readings/R15_Rexach_Guardrail_Log.txt, PASS / WARN / FAIL):
#   GR-01  input files present               GR-07  arms allowed by the Section 4 gate
#   GR-02  R14 finished without a FAIL        GR-08  donors per contrast and min_donors
#   GR-03  Table S1 insula libraries cover    GR-09  models fitted with finite estimates
#          the CELLxGENE donors               GR-10  SEA-AD share above 8% MT
#   GR-04  sex agrees with CELLxGENE           GR-11  R warnings during the run
#   GR-05  ages agree with CELLxGENE           GR-12  output files written
#   GR-06  diagnoses agree with CELLxGENE      GR-09b PMI known for every donor
# A FAIL in GR-01 to GR-06 stops the run before any model is fitted.
#
# Inputs, in External_Dataset_2:
#   R14_readings/R14_Rexach_results.rds and R14_readings/R14_Rexach_Guardrail_Log.txt
#   R12_mapmycells_input/R12_Rexach_nuclei.csv.gz (donor metadata of all insula nuclei)
#   Table S1 (mmc1.xlsx, or a name starting with "Table S1" or "Table_S1"), or its path set in TABLE_S1
#   Part12c_cell_metrics.csv (copied from Fedora: SEA-AD/Part12_archive/Part12c_DoubletsSoupAndCargo/)
# Outputs, in R15_readings/:
#   R15_Rexach_contrasts.xlsx, R15_Rexach_Guardrail_Log.txt, and R15_Rexach_log.txt
#
# Run in R on Windows:
#   source("C:/Users/Wong/Desktop/SEA-AD/External_Dataset_2/R15_Rexach_diagnosis_contrasts.R")

# The data folder: REXACH_ROOT when set (the container sets it to /data), otherwise the Windows folder.
ROOT <- Sys.getenv("REXACH_ROOT", "C:/Users/Wong/Desktop/SEA-AD/External_Dataset_2")
# Table S1 of the paper. NULL searches ROOT for mmc1.xlsx (the name Cell gives Table S1) or a .xlsx
# whose name starts with "Table S1" or "Table_S1"; otherwise give the full path.
TABLE_S1 <- NULL

SET <- list(model = "cov3", covariates = c("age_at_death", "sex", "pmi"), min_donors_combined = 20L,
            min_donors_each = 15L, mt_cut_seaad = 8, seaad_arm = "PVM_3",
            seaad_post_group = "module Post")
# CELLxGENE disease label -> Table S1 clinical_dx label and contrast name.
DX <- data.frame(cxg = c("normal", "Alzheimer disease", "Pick disease", "progressive supranuclear palsy"),
                 s1 = c("Control", "AD", "bvFTD", "PSP-S"),
                 short = c("Control", "AD", "bvFTD", "PSP"), stringsAsFactors = FALSE)

OUT <- file.path(ROOT, "R15_readings")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
LOG <- file.path(OUT, "R15_Rexach_log.txt")
GR_LOG <- file.path(OUT, "R15_Rexach_Guardrail_Log.txt")
INCOMPLETE <- file.path(OUT, "_RUN_INCOMPLETE.txt")
F_RDS <- file.path(ROOT, "R14_readings", "R14_Rexach_results.rds")
F_GR14 <- file.path(ROOT, "R14_readings", "R14_Rexach_Guardrail_Log.txt")
# All insula nuclei (R12), so every donor is present, also donors without Micro-PVM_3 nuclei.
F_NUC <- file.path(ROOT, "R12_mapmycells_input", "R12_Rexach_nuclei.csv.gz")
F_S1 <- if (is.null(TABLE_S1)) list.files(ROOT, pattern = "^(Table[ _.]*S1.*|mmc1)\\.xlsx$", full.names = TRUE,
                                           ignore.case = TRUE) else TABLE_S1
F_P12C <- file.path(ROOT, "Part12c_cell_metrics.csv")

need <- c("cepip", "openxlsx")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Install these packages first: ", paste(miss, collapse = ", "))
suppressPackageStartupMessages(library(cepip))

say <- function(...) cat(sprintf(...), "\n", sep = "")
GR <- new.env()
GR$log <- list(); GR$warnings <- character(0)
gr <- function(id, status, detail) {
    if (length(detail) != 1L || is.na(detail)) detail <- "[detail unavailable]"
    say("  [%s]  %-16s  %s", status, id, detail)
    GR$log[[length(GR$log) + 1L]] <- list(id = id, status = status, detail = detail)
    invisible(status)
}
gr_count <- function(s) sum(vapply(GR$log, function(x) x$status == s, logical(1)))
gr_write <- function(finished) {
    n_p <- gr_count("PASS"); n_w <- gr_count("WARN"); n_f <- gr_count("FAIL")
    writeLines(c(strrep("=", 72), "R15 Rexach diagnosis contrasts - Guardrail Log",
                 paste("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
                 "Readings: external_validation/External2_Rexach_Preset_Readings.md, Section 6",
                 strrep("=", 72), "",
                 vapply(GR$log, function(x) sprintf("  [%s]  %-16s  %s", x$status, x$id, x$detail),
                        character(1)),
                 "", strrep("=", 72),
                 sprintf("  Total : %d", n_p + n_w + n_f), sprintf("  PASS  : %d", n_p),
                 sprintf("  WARN  : %d", n_w), sprintf("  FAIL  : %d", n_f),
                 if (!finished) "\n  Result: the run did not finish; see R15_Rexach_log.txt"
                 else if (n_f == 0 && n_w == 0) "\n  Result: no failures, no warnings"
                 else if (n_f == 0) "\n  Result: no failures; review the warnings above"
                 else "\n  Result: failures present",
                 strrep("=", 72)), GR_LOG)
}
stop_if_fail <- function() if (gr_count("FAIL") > 0)
    stop("Pre-analysis guardrail failure; no model was fitted. See the guardrail log.", call. = FALSE)

main <- function() {
t0 <- Sys.time()
si <- sessionInfo()
say("R15_Rexach_diagnosis_contrasts.R started %s", format(t0))
say("R %s | cepip %s | BLAS: %s | LAPACK: %s", getRversion(), packageVersion("cepip"),
    if (is.null(si$BLAS)) "not reported" else si$BLAS, if (is.null(si$LAPACK)) "not reported" else si$LAPACK)

# ---- Pre-analysis guardrails ----------------------------------------------------------------
say("\n== Pre-analysis guardrails")
files <- c(R14_results = F_RDS, R14_guardrails = F_GR14, R12_nuclei = F_NUC,
           Table_S1 = if (length(F_S1) == 1L) F_S1 else NA_character_)
ok <- !is.na(files) & file.exists(files)
gr("GR-01-inputs", if (all(ok)) "PASS" else "FAIL",
   if (all(ok)) "R14 results and guardrail log, R12 nuclei table, and Table S1 are present"
   else paste("missing:", paste(names(files)[!ok], collapse = ", "),
              if (length(F_S1) > 1L) "(more than one Table S1 .xlsx in the folder)" else "",
              "| .xlsx files in ROOT:", paste(basename(list.files(ROOT, pattern = "\\.xlsx$")), collapse = ", ")))
stop_if_fail()

r14 <- readRDS(F_RDS)
g14 <- readLines(F_GR14, warn = FALSE)
r14_fail <- suppressWarnings(as.integer(sub(".*FAIL\\s*:\\s*", "", grep("^\\s*FAIL\\s*:", g14, value = TRUE))))
r14_done <- !any(grepl("did not finish", g14))
gr("GR-02-R14", if (r14_done && length(r14_fail) == 1L && r14_fail == 0L) "PASS" else "FAIL",
   sprintf("R14 guardrail log: finished %s | FAIL count %s", r14_done,
           if (length(r14_fail)) r14_fail else "not found"))
stop_if_fail()

nuc <- utils::read.csv(gzfile(F_NUC), stringsAsFactors = FALSE)
cxg <- unique(nuc[, c("donor_id", "disease", "sex", "development_stage")])
one_row <- !anyDuplicated(cxg$donor_id)
cxg$autopsy_id <- sub("^R24_", "", cxg$donor_id)
S1_SHEET <- "5.Metadata_perRNAlibrary"
sheets_s1 <- openxlsx::getSheetNames(F_S1)
if (!(S1_SHEET %in% sheets_s1)) {
    gr("GR-03-TableS1", "FAIL", sprintf("%s has no sheet %s (sheets: %s); this is not Table S1",
                                        basename(F_S1), S1_SHEET, paste(sheets_s1, collapse = ", ")))
    stop_if_fail()
}
s1 <- openxlsx::read.xlsx(F_S1, sheet = S1_SHEET)
ins <- s1[tolower(s1$region) %in% "insula", c("library_id", "autopsy_id", "clinical_dx", "age", "sex", "pmi")]
ins$age <- suppressWarnings(as.numeric(ins$age))
cover_ok <- one_row && all(cxg$autopsy_id %in% ins$autopsy_id) && !anyDuplicated(ins$autopsy_id)
gr("GR-03-TableS1", if (cover_ok) "PASS" else "FAIL",
   sprintf("%s: insula libraries %d | CELLxGENE donors %d (one metadata row each: %s), all found %s | in Table S1 only: %s",
           basename(F_S1), nrow(ins), length(unique(cxg$donor_id)), one_row, all(cxg$autopsy_id %in% ins$autopsy_id),
           paste(setdiff(ins$autopsy_id, cxg$autopsy_id), collapse = ", ")))
stop_if_fail()
m <- merge(cxg, ins, by = "autopsy_id", suffixes = c("_cxg", "_s1"))

sex_ok <- all(toupper(substr(m$sex_cxg, 1, 1)) == toupper(m$sex_s1))
gr("GR-04-sex", if (sex_ok) "PASS" else "FAIL",
   sprintf("sex agrees for %d of %d donors", sum(toupper(substr(m$sex_cxg, 1, 1)) == toupper(m$sex_s1)), nrow(m)))

exact <- suppressWarnings(as.numeric(sub("^(\\d+)-year-old.*$", "\\1", m$development_stage)))
is_exact <- grepl("^\\d+-year-old", m$development_stage)
in_range <- ifelse(is_exact, exact == m$age,
            ifelse(grepl("^60-79", m$development_stage), m$age >= 60 & m$age <= 79,
            ifelse(grepl("^80 year-old and over", m$development_stage), m$age >= 80, NA)))
age_ok <- all(is.finite(m$age)) && all(in_range %in% TRUE)
gr("GR-05-age", if (age_ok) "PASS" else "FAIL",
   sprintf("Table S1 ages: %d of %d donors numeric | exact CELLxGENE ages equal Table S1 for %d of %d | ranges consistent for %d of %d",
           sum(is.finite(m$age)), nrow(m), sum(in_range[is_exact] %in% TRUE), sum(is_exact),
           sum(in_range[!is_exact] %in% TRUE), sum(!is_exact)))

m$dx_expected <- DX$s1[match(m$disease, DX$cxg)]
dx_ok <- all(m$clinical_dx == m$dx_expected)
gr("GR-06-diagnosis", if (dx_ok) "PASS" else "FAIL",
   sprintf("clinical_dx agrees with the CELLxGENE disease for %d of %d donors", sum(m$clinical_dx == m$dx_expected),
           nrow(m)))
stop_if_fail()

# ---- Arms allowed by the Section 4 gate ----------------------------------------------------------
g4 <- r14$gate4
allowed <- g4$arm[startsWith(g4$diagnosis_contrasts, "run")]
conditional <- g4$arm[grepl("conditional", g4$diagnosis_contrasts)]
gr("GR-07-gate", if ("PVM_3" %in% allowed) "PASS" else "WARN",
   paste(sprintf("%s: %s", g4$arm, g4$diagnosis_contrasts), collapse = " | "))

# ---- Contrasts ------------------------------------------------------------------------------------------
say("\n== Diagnosis contrasts (Section 6, reported as observed)")
meta <- data.frame(donor = m$donor_id, dx = DX$short[match(m$disease, DX$cxg)],
                   age_at_death = m$age, sex = ifelse(toupper(substr(m$sex_cxg, 1, 1)) == "M", 1, 0),
                   pmi = suppressWarnings(as.numeric(m$pmi)), stringsAsFactors = FALSE)
contrasts <- list(All_vs_Control = c("AD", "bvFTD", "PSP"), AD_vs_Control = "AD",
                  bvFTD_vs_Control = "bvFTD", PSP_vs_Control = "PSP")
rows <- list(); donor_tabs <- list()
for (arm in c("PVM_3", "HKG8")) {
    if (!(arm %in% allowed)) {
        for (ct in names(contrasts)) rows[[length(rows) + 1L]] <- data.frame(
            arm = arm, contrast = ct, n = NA, n_disease = NA, n_control = NA, beta_std = NA, se = NA,
            ci_low = NA, ci_high = NA, p_value = NA, delta_R2 = NA,
            note = g4$diagnosis_contrasts[g4$arm == arm], stringsAsFactors = FALSE)
        next
    }
    dt <- merge(r14$donors[[arm]][, c("donor", "n_total", "n_PostIP", "pct_PostIP")], meta, by = "donor")
    n14 <- nrow(r14$donors[[arm]])
    gr(sprintf("GR-08-%s-donors", arm), if (nrow(dt) == n14 && n14 == length(r14$floor_donors)) "PASS" else "FAIL",
       sprintf("R14 donor table %d rows | floor donors %d | with Table S1 metadata %d", n14,
               length(r14$floor_donors), nrow(dt)))
    if (nrow(dt) != n14) stop("Donors of the R14 table are missing from the metadata.", call. = FALSE)
    donor_tabs[[arm]] <- dt[order(dt$dx, dt$donor), ]
    for (ct in names(contrasts)) {
        d <- dt[dt$dx %in% c("Control", contrasts[[ct]]), ]
        d$diagnosis <- as.numeric(d$dx != "Control")
        md <- if (ct == "All_vs_Control") SET$min_donors_combined else SET$min_donors_each
        a <- cepip_associate(d, "pct_PostIP", outcomes = "diagnosis", covariates = SET$covariates,
                             min_donors = md)
        n_cc <- a$n
        gr(sprintf("GR-08-%s-%s", arm, ct), if (n_cc >= md) "PASS" else "WARN",
           sprintf("%d complete donors (%d disease, %d control) | min_donors %d", n_cc,
                   sum(d$diagnosis == 1), sum(d$diagnosis == 0), md))
        rows[[length(rows) + 1L]] <- data.frame(
            arm = arm, contrast = ct, n = n_cc, n_disease = sum(d$diagnosis == 1),
            n_control = sum(d$diagnosis == 0), beta_std = a$beta_std, se = a$se, ci_low = a$ci_low,
            ci_high = a$ci_high, p_value = a$p_value, delta_R2 = a$delta_R2,
            note = if (arm %in% conditional) "conditional on the Section 4 review" else "",
            stringsAsFactors = FALSE)
    }
}
res <- do.call(rbind, rows)
res$BH_over_4_contrasts <- NA_real_
for (arm in unique(res$arm)) {
    k <- res$arm == arm & is.finite(res$p_value)
    if (any(k)) res$BH_over_4_contrasts[k] <- stats::p.adjust(res$p_value[k], method = "BH")
}
fit_ok <- all(is.finite(res$beta_std[res$arm %in% allowed]))
gr("GR-09-models", if (fit_ok) "PASS" else "WARN",
   sprintf("%d of %d models of the allowed arms have finite estimates",
           sum(is.finite(res$beta_std[res$arm %in% allowed])), sum(res$arm %in% allowed)))
print(res[, c("arm", "contrast", "n", "beta_std", "ci_low", "ci_high", "p_value", "BH_over_4_contrasts", "note")],
      row.names = FALSE, digits = 4)
if (length(intersect(c("PVM_3", "HKG8"), allowed)) == 2L) {
    dif <- merge(res[res$arm == "PVM_3", c("contrast", "beta_std")], res[res$arm == "HKG8", c("contrast", "beta_std")],
                 by = "contrast", suffixes = c("_PVM_3", "_HKG8"))
    dif$difference <- dif$beta_std_PVM_3 - dif$beta_std_HKG8
} else dif <- data.frame(item = "PVM_3 minus HKG8", value = "not available: HKG8 did not pass the Section 4 gate")

summ <- do.call(rbind, lapply(names(donor_tabs), function(a) {
    d <- donor_tabs[[a]]
    do.call(rbind, lapply(split(d, d$dx), function(z) data.frame(
        arm = a, diagnosis = z$dx[1], donors = nrow(z), median_pct_PostIP = stats::median(z$pct_PostIP),
        min = min(z$pct_PostIP), max = max(z$pct_PostIP), median_age = stats::median(z$age_at_death),
        male = sum(z$sex == 1), median_pmi = stats::median(z$pmi), min_pmi = min(z$pmi),
        max_pmi = max(z$pmi), stringsAsFactors = FALSE)))
}))
if (!is.null(summ)) print(summ, row.names = FALSE, digits = 4)

pmi_ok <- all(is.finite(meta$pmi))
gr("GR-09b-PMI", if (pmi_ok) "PASS" else "WARN",
   sprintf("PMI known for %d of %d donors (Table S1)", sum(is.finite(meta$pmi)), nrow(meta)))

# ---- SEA-AD share of Module post-IP cells above 8% MT ----------------------------------------------
say("\n== SEA-AD: Module post-IP cells above %g%% MT (expected loss to the authors' cut)", SET$mt_cut_seaad)
seaad <- data.frame()
if (file.exists(F_P12C)) {
    p12 <- utils::read.csv(F_P12C, stringsAsFactors = FALSE)
    p <- p12[p12$arm == SET$seaad_arm, ]
    post <- p$group == SET$seaad_post_group
    for (cut in c(1, 5, SET$mt_cut_seaad, 10)) seaad <- rbind(seaad, data.frame(
        cutoff_pct_MT = cut, module_postIP_cells = sum(post), module_postIP_above = sum(p$mito_pct[post] > cut),
        share_module_postIP_above = mean(p$mito_pct[post] > cut), other_arm_cells = sum(!post),
        other_above = sum(p$mito_pct[!post] > cut), share_other_above = mean(p$mito_pct[!post] > cut)))
    lv <- unique(p$labels_version)
    s8 <- seaad[seaad$cutoff_pct_MT == SET$mt_cut_seaad, ]
    gr("GR-10-SEAAD_MT8", if (length(lv) == 1L && lv == "2026-09-06_module-pool") "PASS" else "WARN",
       sprintf("SEA-AD PVM_3: %d of %d Module post-IP cells above %g%% MT (%.2f%%) | labels %s",
               s8$module_postIP_above, s8$module_postIP_cells, SET$mt_cut_seaad,
               100 * s8$share_module_postIP_above, paste(lv, collapse = ", ")))
    print(seaad, row.names = FALSE, digits = 4)
} else {
    gr("GR-10-SEAAD_MT8", "WARN", paste("not computed:", basename(F_P12C), "not found in", ROOT))
}

# ---- Save -------------------------------------------------------------------------------------------------
say("\n== Saving")
age_check <- data.frame(donor = m$donor_id, disease = m$disease, sex_cxg = m$sex_cxg, sex_TableS1 = m$sex_s1,
                        development_stage = m$development_stage, age_TableS1 = m$age, pmi_TableS1 = m$pmi,
                        library_TableS1 = m$library_id, consistent = in_range, stringsAsFactors = FALSE)
settings <- data.frame(setting = c(names(SET), "Table S1 file", "SEA-AD file"),
                       value = c(vapply(SET, function(v) paste(v, collapse = ", "), character(1)),
                                 basename(F_S1), basename(F_P12C)), stringsAsFactors = FALSE)
session <- data.frame(item = c("started", "finished", "R", "cepip", "cepip RemoteSha", "openxlsx", "BLAS", "LAPACK"),
                      value = c(format(t0), format(Sys.time()), as.character(getRversion()),
                                as.character(packageVersion("cepip")),
                                if (is.null(packageDescription("cepip")$RemoteSha)) NA else packageDescription("cepip")$RemoteSha,
                                as.character(packageVersion("openxlsx")),
                                if (is.null(si$BLAS)) NA else si$BLAS, if (is.null(si$LAPACK)) NA else si$LAPACK),
                      stringsAsFactors = FALSE)
sheets <- list("01_Contrasts" = res,
               "02_PVM3_minus_HKG8" = dif, "03_By_diagnosis" = summ,
               "04_Donors_PVM_3" = donor_tabs$PVM_3, "05_Donors_HKG8" = donor_tabs$HKG8,
               "06_SEAAD_MT_cut" = seaad, "07_Age_sex_check" = age_check, "08_Section4_gate" = g4,
               "09_Settings" = settings, "10_Session" = session)
sheets <- sheets[!vapply(sheets, function(d) is.null(d) || !nrow(d), logical(1))]
wb <- openxlsx::createWorkbook()
bold <- openxlsx::createStyle(textDecoration = "bold")
for (nm in names(sheets)) {
    openxlsx::addWorksheet(wb, nm)
    openxlsx::writeData(wb, nm, sheets[[nm]], headerStyle = bold)
    openxlsx::freezePane(wb, nm, firstRow = TRUE)
}
xlsx_path <- file.path(OUT, "R15_Rexach_contrasts.xlsx")
openxlsx::saveWorkbook(wb, xlsx_path, overwrite = TRUE)

gr("GR-11-R_warnings", if (length(GR$warnings)) "WARN" else "PASS",
   if (length(GR$warnings)) sprintf("%d R warning(s): %s", length(GR$warnings),
                                    paste(unique(GR$warnings), collapse = " || "))
   else "no R warnings during the run")
gr("GR-12-outputs", if (file.exists(xlsx_path)) "PASS" else "FAIL", basename(xlsx_path))
say("\nGuardrail tally: %d PASS / %d WARN / %d FAIL", gr_count("PASS"), gr_count("WARN"), gr_count("FAIL"))
say("Done in %.1f s", as.numeric(difftime(Sys.time(), t0, units = "secs")))
invisible(TRUE)
}

writeLines(paste("R15 started", format(Sys.time()), "- this file is removed when the run finishes."), INCOMPLETE)
old_out <- file.path(OUT, "R15_Rexach_contrasts.xlsx")
invisible(file.remove(old_out[file.exists(old_out)]))
zz <- file(LOG, open = "wt")
sink(zz, split = TRUE)
finished <- FALSE
tryCatch({
    withCallingHandlers(main(),
        warning = function(w) {
            GR$warnings <- c(GR$warnings, conditionMessage(w))
            say("R warning: %s", conditionMessage(w))
            invokeRestart("muffleWarning")
        },
        error = function(e) say("ERROR: %s", conditionMessage(e)))
    finished <- TRUE
}, finally = {
    gr_write(finished)
    say("Guardrail log written: %s", GR_LOG)
    sink(); close(zz)
    if (finished) unlink(INCOMPLETE)
})
