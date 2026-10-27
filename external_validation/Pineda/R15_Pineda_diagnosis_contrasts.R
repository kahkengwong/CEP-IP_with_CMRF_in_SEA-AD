# R15_Pineda_diagnosis_contrasts.R
# Diagnosis contrasts of external validation 1 (Pineda et al. 2024, BA9 batch 210526). This is the
# Pineda counterpart of R15_Rexach_diagnosis_contrasts.R. The model code is the same; the differences
# are listed under "Pineda changes" below.
#
#   Predictor   Module post-IP% per donor (R14 donor table: Module post-IP cells / all Micro-PVM_3
#               nuclei of the donor x 100), the donors at or above the floor of 5 nuclei.
#   Model       cepip_associate(): OLS on z-scored variables, diagnosis (PN = 0, disease = 1), with
#               age at death and sex (F = 0, M = 1) as covariates (cov2), and with PMI added (cov3,
#               the donors with PMI recorded; PMI is missing for 7 of the 16 PN donors). S01 has no
#               APOE, RIN, brain pH, or brain weight. Complete cases only (FTLD donor 226 has no age
#               at death).
#   Contrasts   FTLD against PN and ALS against PN, each with cov2 and with cov3, with the
#               Benjamini-Hochberg FDR over the four contrasts of each arm. Run only for arms that
#               the IP gate in R14 allows.
#   Arms        HKG8 is fitted with the same models beside PVM_3 when its gate allows it, with the
#               PVM_3 minus HKG8 difference, and is otherwise reported as not available.
#
# Pineda changes from R15_Rexach_diagnosis_contrasts.R:
#   - Donor metadata from S01 (sample_mapping and donor_information) in place of Table S1 and
#     CELLxGENE. The checks compare the S01 group with the diagnosis in the library name, and the
#     sex of S01 with the sex carried through R11 to R14.
#   - The four Pineda contrasts above, in place of the four Rexach contrasts.
#   - The SEA-AD share above 8% MT is not computed. It applies to the 8% MT cut of the Rexach
#     authors' QC.
#
# Guardrails (ledger R15_contrasts/R15_Pineda_Guardrail_Log.txt, PASS / WARN / FAIL):
#   GR-01  input files present               GR-07  arms allowed by the IP gate
#   GR-02  R14 finished without a FAIL        GR-08  donors per contrast and min_donors
#   GR-03  S01 covers the 63 libraries,       GR-09  models fitted with finite estimates
#          one library per donor              GR-11  R warnings during the run
#   GR-04  sex agrees with the R14 donors      GR-12  output files written
#   GR-05  ages: missing only for the
#          expected donors
#   GR-06  S01 group agrees with the library name
# A FAIL in GR-01 to GR-06 stops the run before any model is fitted.
#
# Inputs, in External_Dataset_1:
#   R14_core/R14_Pineda_results.rds and R14_core/R14_Pineda_Guardrail_Log.txt
#   The S01 workbook (S01_sample_information.xlsx, Synapse syn51105515), or its path set in S01
#   Pineda_BA9_210526_manifest.tsv
# Outputs, in R15_contrasts/:
#   R15_Pineda_contrasts.xlsx, R15_Pineda_Guardrail_Log.txt, and R15_Pineda_log.txt
#
# Run in R on Windows:
#   source("C:/Users/Wong/Desktop/SEA-AD/External_Dataset_1/R15_Pineda_diagnosis_contrasts.R")

# The data folder: PINEDA_ROOT when set (the container sets it to /data), otherwise the Windows folder.
ROOT <- Sys.getenv("PINEDA_ROOT", "C:/Users/Wong/Desktop/SEA-AD/External_Dataset_1")
# S01 workbook. NULL searches ROOT for a .xlsx whose name starts with "S01"; otherwise give the path.
S01 <- NULL

SET <- list(model = "cov2", covariates = c("age_at_death", "sex"),
            covariates_pmi = c("age_at_death", "sex", "pmi"), min_donors = 20L,
            # Donors without age at death in S01 ("NR").
            age_missing_expected = "FTLD_226")
GROUP_TO_DX <- c(ALS = "ALS", C9ALS = "ALS", FTLD = "FTLD", C9FTLD = "FTLD", PN = "PN")
# Contrasts: the disease group against PN, and the covariates.
CONTRASTS <- list(
    FTLD_vs_PN = list(dx = "FTLD", covariates = SET$covariates),
    FTLD_vs_PN_PMI = list(dx = "FTLD", covariates = SET$covariates_pmi),
    ALS_vs_PN = list(dx = "ALS", covariates = SET$covariates),
    ALS_vs_PN_PMI = list(dx = "ALS", covariates = SET$covariates_pmi))

OUT <- file.path(ROOT, "R15_contrasts")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
LOG <- file.path(OUT, "R15_Pineda_log.txt")
GR_LOG <- file.path(OUT, "R15_Pineda_Guardrail_Log.txt")
INCOMPLETE <- file.path(OUT, "_RUN_INCOMPLETE.txt")
F_RDS <- file.path(ROOT, "R14_core", "R14_Pineda_results.rds")
F_GR14 <- file.path(ROOT, "R14_core", "R14_Pineda_Guardrail_Log.txt")
F_MAN <- file.path(ROOT, "Pineda_BA9_210526_manifest.tsv")
F_S01 <- if (is.null(S01)) list.files(ROOT, pattern = "^S01.*\\.xlsx$", full.names = TRUE, ignore.case = TRUE) else S01

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
    writeLines(c(strrep("=", 72), "R15 Pineda diagnosis contrasts - Guardrail Log",
                 paste("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
                 strrep("=", 72), "",
                 vapply(GR$log, function(x) sprintf("  [%s]  %-16s  %s", x$status, x$id, x$detail),
                        character(1)),
                 "", strrep("=", 72),
                 sprintf("  Total : %d", n_p + n_w + n_f), sprintf("  PASS  : %d", n_p),
                 sprintf("  WARN  : %d", n_w), sprintf("  FAIL  : %d", n_f),
                 if (!finished) "\n  Result: the run did not finish; see R15_Pineda_log.txt"
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
say("R15_Pineda_diagnosis_contrasts.R started %s", format(t0))
say("R %s | cepip %s | BLAS: %s | LAPACK: %s", getRversion(), packageVersion("cepip"),
    if (is.null(si$BLAS)) "not reported" else si$BLAS, if (is.null(si$LAPACK)) "not reported" else si$LAPACK)

# ---- Pre-analysis guardrails ----------------------------------------------------------------
say("\n== Pre-analysis guardrails")
files <- c(R14_results = F_RDS, R14_guardrails = F_GR14, manifest = F_MAN,
           S01 = if (length(F_S01) == 1L) F_S01 else NA_character_)
ok <- !is.na(files) & file.exists(files)
gr("GR-01-inputs", if (all(ok)) "PASS" else "FAIL",
   if (all(ok)) "R14 results and guardrail log, the manifest, and S01 are present"
   else paste("missing:", paste(names(files)[!ok], collapse = ", "),
              if (length(F_S01) > 1L) "(more than one S01 .xlsx in the folder)" else ""))
stop_if_fail()

r14 <- readRDS(F_RDS)
g14 <- readLines(F_GR14, warn = FALSE)
r14_fail <- suppressWarnings(as.integer(sub(".*FAIL\\s*:\\s*", "", grep("^\\s*FAIL\\s*:", g14, value = TRUE))))
r14_done <- !any(grepl("did not finish", g14))
gr("GR-02-R14", if (r14_done && length(r14_fail) == 1L && r14_fail == 0L) "PASS" else "FAIL",
   sprintf("R14 guardrail log: finished %s | FAIL count %s", r14_done,
           if (length(r14_fail)) r14_fail else "not found"))
stop_if_fail()

libs <- utils::read.delim(F_MAN, stringsAsFactors = FALSE, check.names = FALSE)[[1]]
sm <- openxlsx::read.xlsx(F_S01, sheet = "sample_mapping")
di <- openxlsx::read.xlsx(F_S01, sheet = "donor_information", check.names = FALSE, sep.names = " ")
col_of <- function(pattern) {
    k <- grep(pattern, names(di), ignore.case = TRUE)
    if (length(k) != 1L) stop("S01 donor_information: column matching '", pattern, "' found ", length(k),
                              " times.", call. = FALSE)
    names(di)[k]
}
num <- function(v) suppressWarnings(as.numeric(as.character(v)))
meta <- data.frame(library_id = libs, stringsAsFactors = FALSE)
meta$donor_number <- as.character(as.integer(sm$Donor[match(libs, sm$Sample_ID)]))
meta$donor <- sub("^[0-9]+_([A-Za-z]+_[0-9]+)_.*$", "\\1", libs)
meta$dx_name <- sub("^[0-9]+_([A-Za-z]+)_[0-9]+_.*$", "\\1", libs)
m <- match(meta$donor_number, as.character(as.integer(di$Donor)))
meta$group_S01 <- di$Group[m]
meta$dx <- unname(GROUP_TO_DX[meta$group_S01])
meta$sex_S01 <- di$Sex[m]
meta$age_at_death <- num(di[[col_of("^Age of death")]][m])
meta$age_S01_text <- as.character(di[[col_of("^Age of death")]][m])
meta$pmi <- num(di[[col_of("^PMI")]][m])
meta$pmi_S01_text <- as.character(di[[col_of("^PMI")]][m])
cover_ok <- length(libs) == 63L && all(!is.na(meta$donor_number)) && all(!is.na(m)) &&
    !anyDuplicated(meta$donor_number) && all(r14$floor_donors %in% meta$donor)
gr("GR-03-S01", if (cover_ok) "PASS" else "FAIL",
   sprintf("%s: %d manifest libraries | mapped to an S01 donor %d | one library per donor %s | R14 floor donors found %d of %d",
           basename(F_S01), length(libs), sum(!is.na(m)), !anyDuplicated(meta$donor_number),
           sum(r14$floor_donors %in% meta$donor), length(r14$floor_donors)))
stop_if_fail()

d14 <- r14$donors$PVM_3
sx <- meta$sex_S01[match(d14$donor, meta$donor)]
sex_ok <- !is.null(d14) && all(sx %in% c("F", "M")) && all(sx == d14$sex)
gr("GR-04-sex", if (sex_ok) "PASS" else "FAIL",
   sprintf("S01 sex agrees with the R14 donor table for %d of %d donors", sum(sx == d14$sex, na.rm = TRUE),
           if (is.null(d14)) 0L else nrow(d14)))

age_na <- sort(meta$donor[!is.finite(meta$age_at_death)])
age_rng <- range(meta$age_at_death, na.rm = TRUE)
age_ok <- identical(age_na, sort(SET$age_missing_expected)) && age_rng[1] >= 18 && age_rng[2] <= 110
gr("GR-05-age", if (age_ok) "PASS" else "FAIL",
   sprintf("age at death numeric for %d of %d donors (range %.1f to %.1f) | missing: %s (expected: %s)",
           sum(is.finite(meta$age_at_death)), nrow(meta), age_rng[1], age_rng[2],
           paste(age_na, collapse = ", "), paste(SET$age_missing_expected, collapse = ", ")))

dx_ok <- all(meta$dx == meta$dx_name, na.rm = FALSE) && !anyNA(meta$dx)
gr("GR-06-diagnosis", if (isTRUE(dx_ok)) "PASS" else "FAIL",
   sprintf("S01 group (C9 groups merged) agrees with the library name for %d of %d libraries | %s",
           sum(meta$dx == meta$dx_name, na.rm = TRUE), nrow(meta),
           paste(sprintf("%s %d", names(table(meta$group_S01)), as.integer(table(meta$group_S01))),
                 collapse = ", ")))
stop_if_fail()

# ---- Arms allowed by the IP gate ----------------------------------------------------------------
g4 <- r14$gate4
allowed <- g4$arm[startsWith(g4$diagnosis_contrasts, "run")]
conditional <- g4$arm[grepl("conditional", g4$diagnosis_contrasts)]
gr("GR-07-gate", if ("PVM_3" %in% allowed) "PASS" else "WARN",
   paste(sprintf("%s: %s", g4$arm, g4$diagnosis_contrasts), collapse = " | "))

# ---- Contrasts ------------------------------------------------------------------------------------------
say("\n== Diagnosis contrasts")
mt <- data.frame(donor = meta$donor, dx = meta$dx, age_at_death = meta$age_at_death,
                 sex = ifelse(meta$sex_S01 == "M", 1, 0), pmi = meta$pmi, stringsAsFactors = FALSE)
rows <- list(); donor_tabs <- list()
for (arm in c("PVM_3", "HKG8")) {
    if (!(arm %in% allowed)) {
        for (ct in names(CONTRASTS)) rows[[length(rows) + 1L]] <- data.frame(
            arm = arm, contrast = ct,
            covariates = paste(CONTRASTS[[ct]]$covariates, collapse = ", "), n = NA, n_disease = NA,
            n_control = NA, beta_std = NA, se = NA, ci_low = NA, ci_high = NA, p_value = NA, delta_R2 = NA,
            note = g4$diagnosis_contrasts[g4$arm == arm], stringsAsFactors = FALSE)
        next
    }
    dt <- merge(r14$donors[[arm]][, c("donor", "n_total", "n_PostIP", "pct_PostIP")], mt, by = "donor")
    n14 <- nrow(r14$donors[[arm]])
    gr(sprintf("GR-08-%s-donors", arm), if (nrow(dt) == n14 && n14 == length(r14$floor_donors)) "PASS" else "FAIL",
       sprintf("R14 donor table %d rows | floor donors %d | with S01 metadata %d", n14,
               length(r14$floor_donors), nrow(dt)))
    if (nrow(dt) != n14) stop("Donors of the R14 table are missing from the metadata.", call. = FALSE)
    donor_tabs[[arm]] <- dt[order(dt$dx, dt$donor), ]
    for (ct in names(CONTRASTS)) {
        cs <- CONTRASTS[[ct]]
        d <- dt[dt$dx %in% c("PN", cs$dx), ]
        d$diagnosis <- as.numeric(d$dx != "PN")
        cc <- stats::complete.cases(d[, c("pct_PostIP", "diagnosis", cs$covariates)])
        nd <- sum(d$diagnosis[cc] == 1); nc <- sum(d$diagnosis[cc] == 0)
        a <- cepip_associate(d, "pct_PostIP", outcomes = "diagnosis", covariates = cs$covariates,
                             min_donors = SET$min_donors)
        gr(sprintf("GR-08-%s-%s", arm, ct), if (isTRUE(a$n >= SET$min_donors) && a$n == sum(cc)) "PASS" else "WARN",
           sprintf("%d complete donors (%d %s, %d PN) | left out for missing covariates: %s | min_donors %d",
                   sum(cc), nd, cs$dx, nc, if (any(!cc)) paste(d$donor[!cc], collapse = ", ") else "none",
                   SET$min_donors))
        rows[[length(rows) + 1L]] <- data.frame(
            arm = arm, contrast = ct, covariates = paste(cs$covariates, collapse = ", "),
            n = a$n, n_disease = nd, n_control = nc, beta_std = a$beta_std, se = a$se, ci_low = a$ci_low,
            ci_high = a$ci_high, p_value = a$p_value, delta_R2 = a$delta_R2,
            note = if (arm %in% conditional) "conditional on the IP gate review" else "",
            stringsAsFactors = FALSE)
    }
}
res <- do.call(rbind, rows)
# Benjamini-Hochberg FDR over the four contrasts of each arm (FTLD and ALS against PN, each with age at
# death and sex, and with PMI added), as for the four Rexach contrasts.
res$FDR <- NA_real_
for (a_ in unique(res$arm)) {
    k <- res$arm == a_ & is.finite(res$p_value)
    if (any(k)) res$FDR[k] <- stats::p.adjust(res$p_value[k], method = "BH")
}
fit_ok <- all(is.finite(res$beta_std[res$arm %in% allowed]))
gr("GR-09-models", if (fit_ok) "PASS" else "WARN",
   sprintf("%d of %d models of the allowed arms have finite estimates",
           sum(is.finite(res$beta_std[res$arm %in% allowed])), sum(res$arm %in% allowed)))
print(res[, c("arm", "contrast", "n", "beta_std", "ci_low", "ci_high", "p_value", "FDR", "note")],
      row.names = FALSE, digits = 4)
if (length(intersect(c("PVM_3", "HKG8"), allowed)) == 2L) {
    dif <- merge(res[res$arm == "PVM_3", c("contrast", "beta_std")], res[res$arm == "HKG8", c("contrast", "beta_std")],
                 by = "contrast", suffixes = c("_PVM_3", "_HKG8"))
    dif$difference <- dif$beta_std_PVM_3 - dif$beta_std_HKG8
} else dif <- data.frame(item = "PVM_3 minus HKG8", value = "not available: an arm did not pass the IP gate")

summ <- do.call(rbind, lapply(names(donor_tabs), function(a) {
    d <- donor_tabs[[a]]
    do.call(rbind, lapply(split(d, d$dx), function(z) data.frame(
        arm = a, diagnosis = z$dx[1], donors = nrow(z), median_pct_PostIP = stats::median(z$pct_PostIP),
        min = min(z$pct_PostIP), max = max(z$pct_PostIP), median_age = stats::median(z$age_at_death, na.rm = TRUE),
        age_missing = sum(!is.finite(z$age_at_death)), male = sum(z$sex == 1),
        median_pmi = stats::median(z$pmi, na.rm = TRUE), pmi_missing = sum(!is.finite(z$pmi)),
        stringsAsFactors = FALSE)))
}))
if (!is.null(summ)) print(summ, row.names = FALSE, digits = 4)

# ---- Save -------------------------------------------------------------------------------------------------
say("\n== Saving")
settings <- data.frame(setting = c(names(SET), "contrasts", "S01 file"),
                       value = c(vapply(SET, function(v) paste(v, collapse = ", "), character(1)),
                                 paste(sprintf("%s: %s vs PN, %s", names(CONTRASTS),
                                               vapply(CONTRASTS, `[[`, "", "dx"),
                                               vapply(CONTRASTS, function(x) paste(x$covariates, collapse = " + "), "")),
                                       collapse = " | "),
                                 basename(F_S01)), stringsAsFactors = FALSE)
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
               "06_Metadata_check" = meta, "07_IP_gate" = g4,
               "08_Settings" = settings, "09_Session" = session)
sheets <- sheets[!vapply(sheets, function(d) is.null(d) || !nrow(d), logical(1))]
wb <- openxlsx::createWorkbook()
bold <- openxlsx::createStyle(textDecoration = "bold")
for (nm in names(sheets)) {
    openxlsx::addWorksheet(wb, nm)
    openxlsx::writeData(wb, nm, sheets[[nm]], headerStyle = bold)
    openxlsx::freezePane(wb, nm, firstRow = TRUE)
}
xlsx_path <- file.path(OUT, "R15_Pineda_contrasts.xlsx")
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
old_out <- file.path(OUT, "R15_Pineda_contrasts.xlsx")
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
