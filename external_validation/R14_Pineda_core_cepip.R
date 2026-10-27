# R14_Pineda_core_cepip.R
# Core cepip run on external validation 1 (Pineda et al. 2024, BA9 batch 210526), following the
# committed readings (external_validation/External1_Pineda_Preset_Readings.md, commit in
# READINGS_COMMIT below). This is R14_Rexach_core_cepip.R with the Pineda inputs; the analysis
# code is the same. The differences are listed under "Pineda changes" below.
#
#   Section 3  Do the eight DFGs replicate? MT-CO3 screen on all Micro-PVM_3 nuclei (donors with
#              >= 5 nuclei), tier 1 (rs >= 0.6, tau_b >= 0.5) and tier 2 (rs >= 0.5, tau_b >= 0.4).
#   Section 4  IP gate of each arm: IPRS tier, the C1 to C5 flag, and, for Moderate, Weak, or a
#              component at 0, the CMPB 2026 Section 3.12 mosaic comparison (Fisher's exact test,
#              pre-IP and post-IP split on the GOI axis as in CMPB 2026).
#   Section 6  Observed items that need no ages: within-positive rs and tau_b, the agnostic screen
#              at both tiers, pipeline counts, the same-cell control (cDFG against cHKG9), the
#              known-filter test (1%, 5%, and 10% MT, and the authors' published QC), the share of
#              Module post-IP cells that are Micro-PVM_3 in both mappings, per-donor Module
#              post-IP%, and DEGs.
# The association with diagnosis (Section 5) needs the S01 ages and is run by R15, only for arms
# that pass the Section 4 gate.
#
# Pineda changes from R14_Rexach_core_cepip.R:
#   - Inputs from R13_Pineda_count_checks.py; EXPECT holds the Pineda count checks (fill from the
#     R13 log before the run; GR-07 fails while a value is NA).
#   - Known-filter test: cutoffs 1%, 5%, and 10% MT (no %MT cut was applied to the Pineda reads),
#     and the share removed by the authors' published QC (nuclei absent from the authors' metadata).
#   - GR-18 reports the %MT range of the PVM_3 arm. There is no authors' cut to check against.
#   - GR-13: a cepip error in the HKG8 control arm is a WARN, not a FAIL (readings dated note of
#     2026-10-26). HKG8 is then reported as not available, and R15 runs Section 5 on PVM_3. An error
#     in PVM_3 or HKG_control is still a FAIL (panel B of R16 needs HKG_control).
#
# Arms (readings Section 1):
#   PVM_3         MT-CO3 against cDFG (mean of the eight DFGs), Micro-PVM_3 (Deep Generative
#                 Mapping) nuclei with MT-CO3 > 0, all donors
#   HKG8          ACTB against cHKG8, PVM_3 arm nuclei with ACTB > 0
#   HKG_control   MT-CO3 against cHKG9 on the PVM_3 arm nuclei (same-cell control; CEP-IP only, as in
#                 the SEA-AD reproduction script)
# Expression is ln(CP10K + 1), with CP10K from each nucleus's total UMIs over all genes.
#
# Guardrails (ledger R14_readings/R14_Pineda_Guardrail_Log.txt, PASS / WARN / FAIL):
#   GR-01  input files present              GR-11  Section 3 screen cells = floor donors' nuclei
#   GR-02  matrix and tables consistent      GR-12  Section 3 pass flags match rs and tau_b
#   GR-03  raw counts integer, non-negative  GR-13  cepip run completed (each arm; WARN for HKG8)
#   GR-04  total UMIs equal R13              GR-14  IP found and IPRS tier (each arm)
#   GR-05  ln(CP10K + 1) pins 10,000         GR-15  CMRF classes cover the arm; Module =
#   GR-06  arm genes present and unique             OriRef + Recovered (each arm)
#   GR-07  counts equal the committed        GR-16  Validation1 present (each arm)
#          count checks (R13)                GR-17  donor table consistent (each arm)
#   GR-08  arms nested and GOI-positive      GR-18  PVM_3 arm %MT range (no cut applied)
#   GR-09  cepip version and commit          GR-19  R warnings during the run
#   GR-10  settings equal the readings       GR-20  output files written
# A FAIL in GR-01 to GR-10 stops the run before any statistic. R14_readings/_RUN_INCOMPLETE.txt
# exists while the run is going and is removed when it finishes.
#
# Input (from R13, Python in WSL): R13_count_checks/R13_Pineda_MicroPVM3_counts.mtx.gz,
#   R13_Pineda_MicroPVM3_genes.csv.gz, and R13_Pineda_MicroPVM3_nuclei.csv.gz.
# Outputs, one flat folder per arm plus one for the readings:
#   R14_PVM_3/, R14_HKG8/, R14_HKG_control/   cepip_save() tables, plots, and settings
#   R14_readings/R14_Pineda_readings.xlsx     the Section 3, 4, and 6 results
#   R14_readings/R14_Pineda_Guardrail_Log.txt guardrail ledger
#   R14_readings/R14_Pineda_log.txt           console output
#   R14_readings/R14_Pineda_results.rds       the cepip objects and donor tables, for R15
#
# Run in R on Windows:
#   source("C:/Users/Wong/Desktop/SEA-AD/External_Dataset_1/R14_Pineda_core_cepip.R")

ROOT    <- "C:/Users/Wong/Desktop/SEA-AD/External_Dataset_1"
WORKERS <- 6L
# Commit of the committed Pineda readings (filled at the commit).
READINGS_COMMIT <- "[fill at commit]"

# The method settings of the committed readings (Section 1). GR-10 checks them against READINGS.
SET <- list(
    goi = "MT-CO3",
    dfg8 = c("MT-ATP6", "MT-CO1", "MT-CO2", "MT-CYB", "MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4"),
    hkg_goi = "ACTB",
    hk8 = c("EEF1A1", "RPL13A", "B2M", "YWHAZ", "FTL", "PSAP", "RAB7A", "HNRNPA2B1"),
    donor_floor = 5L,
    tier1 = c(rs = 0.6, tau = 0.5), tier2 = c(rs = 0.5, tau = 0.4), n_pass = 5L,
    ip_k = 10L, ip_axis = "y", kappa = 0.5, bw = 0.1, repeats = 40L,
    min_arm_cells = 50L, mt_cutoffs = c(1, 5, 10), review_alpha = 0.05,
    rng_kind = c("Mersenne-Twister", "Inversion", "Rejection"))
# The same settings as written in the committed readings, kept separately so that GR-10 detects an
# edit of SET.
READINGS <- list(goi = "MT-CO3", dfg8 = c("MT-ATP6", "MT-CO1", "MT-CO2", "MT-CYB", "MT-ND1", "MT-ND2",
                                          "MT-ND3", "MT-ND4"),
                 hkg_goi = "ACTB", hk8 = c("EEF1A1", "RPL13A", "B2M", "YWHAZ", "FTL", "PSAP", "RAB7A",
                                           "HNRNPA2B1"),
                 donor_floor = 5L, tier1 = c(rs = 0.6, tau = 0.5), tier2 = c(rs = 0.5, tau = 0.4),
                 n_pass = 5L, ip_k = 10L, ip_axis = "y", kappa = 0.5, bw = 0.1, repeats = 40L,
                 min_arm_cells = 50L, mt_cutoffs = c(1, 5, 10), review_alpha = 0.05)
# Commits of cepip 0.1.0 with the code of the readings (4c4f70c is the CRAN submission).
COMMITS <- c("fe4b269", "4c4f70c")
# Counts from the committed count checks (readings Section 2, R13 log sections D and E, run of
# 2026-10-27). GR-07 fails if the inputs give other counts.
EXPECT <- list(micro_pvm3 = 3442L, pvm3_arm = 2795L, hkg8_arm = 1632L,
               floor_donors = 59L,
               floor_by_disease = c("ALS" = 22L, "FTLD" = 23L, "PN" = 14L))

IN  <- file.path(ROOT, "R13_count_checks")
OUT <- file.path(ROOT, "R14_readings")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
LOG <- file.path(OUT, "R14_Pineda_log.txt")
GR_LOG <- file.path(OUT, "R14_Pineda_Guardrail_Log.txt")
INCOMPLETE <- file.path(OUT, "_RUN_INCOMPLETE.txt")
FILES_IN <- file.path(IN, c("R13_Pineda_MicroPVM3_counts.mtx.gz", "R13_Pineda_MicroPVM3_genes.csv.gz",
                            "R13_Pineda_MicroPVM3_nuclei.csv.gz"))

need <- c("cepip", "Matrix", "openxlsx")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Install these packages first: ", paste(miss, collapse = ", "))
suppressPackageStartupMessages({ library(cepip); library(Matrix) })
do.call(RNGkind, as.list(SET$rng_kind))

say <- function(...) cat(sprintf(...), "\n", sep = "")
# cepip progress messages go to the console and the log.
quiet_msgs <- function(expr) withCallingHandlers(expr, message = function(m) {
    cat(conditionMessage(m)); invokeRestart("muffleMessage") })

# ---- Guardrail ledger ------------------------------------------------------------------------------
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
    writeLines(c(strrep("=", 72),
                 "R14 Pineda core cepip - Guardrail Log",
                 paste("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
                 sprintf("Readings: external_validation/External1_Pineda_Preset_Readings.md (commit %s)",
                         READINGS_COMMIT),
                 "Normalisation: ln(CP10K + 1) from the raw counts of R13",
                 strrep("=", 72), "",
                 vapply(GR$log, function(x) sprintf("  [%s]  %-16s  %s", x$status, x$id, x$detail),
                        character(1)),
                 "", strrep("=", 72),
                 sprintf("  Total : %d", n_p + n_w + n_f), sprintf("  PASS  : %d", n_p),
                 sprintf("  WARN  : %d", n_w), sprintf("  FAIL  : %d", n_f),
                 if (!finished) "\n  Result: the run did not finish; see R14_Pineda_log.txt"
                 else if (n_f == 0 && n_w == 0) "\n  Result: no failures, no warnings"
                 else if (n_f == 0) "\n  Result: no failures; review the warnings above"
                 else "\n  Result: failures present",
                 strrep("=", 72)), GR_LOG)
}
stop_if_fail <- function() if (gr_count("FAIL") > 0)
    stop("Pre-analysis guardrail failure; no statistic was computed. See the guardrail log.",
         call. = FALSE)

main <- function() {
t0 <- Sys.time()
say("R14_Pineda_core_cepip.R started %s | workers %d", format(t0), WORKERS)
say("BLAS: %s | LAPACK: %s", extSoftVersion()[["BLAS"]], La_library())

# ---- Pre-analysis guardrails ----------------------------------------------------------------------------
say("\n== Pre-analysis guardrails")
ok_in <- file.exists(FILES_IN)
gr("GR-01-inputs", if (all(ok_in)) "PASS" else "FAIL",
   if (all(ok_in)) "the three R13 export files are present"
   else paste("missing:", paste(basename(FILES_IN[!ok_in]), collapse = ", ")))
stop_if_fail()

M <- Matrix::readMM(gzfile(FILES_IN[1]))
M <- methods::as(M, "CsparseMatrix")
genes <- utils::read.csv(gzfile(FILES_IN[2]), stringsAsFactors = FALSE)
nuc <- utils::read.csv(gzfile(FILES_IN[3]), stringsAsFactors = FALSE)
as_flag <- function(v) toupper(as.character(v)) %in% c("TRUE", "T", "1")
for (f in c("is_pvm3_dgm", "is_pvm3_hm", "in_pvm3_arm", "in_hkg8_arm")) nuc[[f]] <- as_flag(nuc[[f]])
# In the authors' published QC: TRUE or FALSE, or NA for libraries without the authors' metadata.
nuc$in_author_qc <- ifelse(is.na(nuc$in_author_qc) | nuc$in_author_qc %in% "", NA, as_flag(nuc$in_author_qc))
nuc$donor_id <- as.character(nuc$donor_id)

dims_ok <- nrow(M) == nrow(genes) && ncol(M) == nrow(nuc) && !anyDuplicated(nuc$cell_id)
gr("GR-02-tables", if (dims_ok) "PASS" else "FAIL",
   sprintf("matrix %d genes x %d nuclei | gene table %d rows | nuclei table %d rows | duplicated cell IDs %d",
           nrow(M), ncol(M), nrow(genes), nrow(nuc), sum(duplicated(nuc$cell_id))))
stop_if_fail()
rownames(M) <- make.unique(genes$symbol)
colnames(M) <- nuc$cell_id

int_ok <- all(M@x >= 0) && all(M@x == round(M@x))
gr("GR-03-counts", if (int_ok) "PASS" else "FAIL",
   sprintf("raw counts: min %s | max %s | all whole numbers %s", min(M@x), max(M@x),
           all(M@x == round(M@x))))
lib <- Matrix::colSums(M)
umi_dev <- max(abs(lib - nuc$total_umis))
gr("GR-04-total_UMIs", if (umi_dev == 0) "PASS" else "FAIL",
   sprintf("largest difference from the R13 total UMIs: %s", umi_dev))
stop_if_fail()
N <- M
N@x <- log1p(N@x / rep(lib, diff(N@p)) * 1e4)   # ln(CP10K + 1)
E1 <- N; E1@x <- expm1(E1@x)
cs_dev <- max(abs(Matrix::colSums(E1) - 1e4)); rm(E1)
gr("GR-05-normalisation", if (cs_dev < 1e-6 * 1e4 && max(N@x) <= log1p(1e4) + 1e-9) "PASS" else "FAIL",
   sprintf("expm1 column sums pin 10,000 within %.2e | max %.4f (ceiling %.4f)", cs_dev, max(N@x),
           log1p(1e4)))

needed <- c(SET$goi, SET$dfg8, SET$hkg_goi, SET$hk8)
dup <- intersect(needed, genes$symbol[duplicated(genes$symbol)])
absent <- setdiff(needed, rownames(M))
gr("GR-06-genes", if (!length(dup) && !length(absent)) "PASS" else "FAIL",
   sprintf("%d arm genes | missing: %s | duplicated symbols: %s | duplicated symbols elsewhere: %d",
           length(needed), if (length(absent)) paste(absent, collapse = ", ") else "none",
           if (length(dup)) paste(dup, collapse = ", ") else "none", sum(duplicated(genes$symbol))))

pvm3_all <- which(nuc$is_pvm3_dgm)
n_pvm3 <- table(nuc$donor_id[pvm3_all])
floor_donors <- names(n_pvm3)[n_pvm3 >= SET$donor_floor]
donor_info <- unique(data.frame(donor = nuc$donor_id[pvm3_all], disease = nuc$disease[pvm3_all],
                                sex = nuc$sex[pvm3_all], stringsAsFactors = FALSE))
by_dis <- table(donor_info$disease[donor_info$donor %in% floor_donors])
got <- c(micro_pvm3 = length(pvm3_all), pvm3_arm = sum(nuc$in_pvm3_arm),
         hkg8_arm = sum(nuc$in_hkg8_arm), floor_donors = length(floor_donors))
exp_v <- unlist(EXPECT[c("micro_pvm3", "pvm3_arm", "hkg8_arm", "floor_donors")])
dis_ok <- isTRUE(all(names(EXPECT$floor_by_disease) %in% names(by_dis)) &&
    all(as.integer(by_dis[names(EXPECT$floor_by_disease)]) == EXPECT$floor_by_disease))
filled <- !anyNA(exp_v) && !anyNA(EXPECT$floor_by_disease)
gr("GR-07-count_checks", if (filled && isTRUE(all(got == exp_v)) && dis_ok) "PASS" else "FAIL",
   sprintf("%sMicro-PVM_3 %d (%s) | PVM_3 arm %d (%s) | HKG8 arm %d (%s) | floor donors %d (%s) | %s; committed values in brackets",
           if (filled) "" else "EXPECT is not filled from the R13 log | ",
           got[1], exp_v[1], got[2], exp_v[2], got[3], exp_v[3], got[4], exp_v[4],
           paste(sprintf("%s %d", names(by_dis), as.integer(by_dis)), collapse = ", ")))

arm_cells <- list(PVM_3 = which(nuc$in_pvm3_arm), HKG8 = which(nuc$in_hkg8_arm))
goi_raw <- as.numeric(M[SET$goi, ]); hkg_raw <- as.numeric(M[SET$hkg_goi, ])
nest_ok <- all(arm_cells$PVM_3 %in% pvm3_all) && all(arm_cells$HKG8 %in% arm_cells$PVM_3) &&
    all(goi_raw[arm_cells$PVM_3] > 0) && all(hkg_raw[arm_cells$HKG8] > 0) &&
    !any(goi_raw[setdiff(pvm3_all, arm_cells$PVM_3)] > 0) &&
    !any(hkg_raw[setdiff(arm_cells$PVM_3, arm_cells$HKG8)] > 0)
gr("GR-08-arms", if (nest_ok) "PASS" else "FAIL",
   "HKG8 arm within the PVM_3 arm within Micro-PVM_3; MT-CO3 > 0 exactly in the PVM_3 arm; ACTB > 0 exactly in the HKG8 arm")

sha <- packageDescription("cepip")$RemoteSha
ver_ok <- as.character(packageVersion("cepip")) == "0.1.0"
gr("GR-09-cepip", if (!ver_ok) "FAIL" else if (is.null(sha)) "WARN"
   else if (any(startsWith(sha, COMMITS))) "PASS" else "FAIL",
   sprintf("cepip %s | installed from commit %s | readings commit fe4b269 (code identical: 4c4f70c)",
           packageVersion("cepip"), if (is.null(sha)) "not recorded by the installer" else sha))

same_set <- identical(SET[names(READINGS)], READINGS)
gr("GR-10-settings", if (same_set) "PASS" else "FAIL",
   if (same_set) sprintf("ip_k %d | ip_axis %s | kappa %.2f | bw %.2f | repeats %d | tiers %s and %s | floor %d",
                         SET$ip_k, SET$ip_axis, SET$kappa, SET$bw, SET$repeats,
                         paste(SET$tier1, collapse = "/"), paste(SET$tier2, collapse = "/"),
                         SET$donor_floor)
   else paste("SET differs from the committed readings in:", paste(names(READINGS)[!mapply(identical,
       SET[names(READINGS)], READINGS)], collapse = ", ")))
stop_if_fail()
say("Pre-analysis guardrails passed: %d PASS, %d WARN", gr_count("PASS"), gr_count("WARN"))

dense <- function(g, cells) t(as.matrix(N[g, cells, drop = FALSE]))
has_pcapp <- requireNamespace("pcaPP", quietly = TRUE)
ktau <- if (has_pcapp) pcaPP::cor.fk else function(a, b) stats::cor(a, b, method = "kendall")
say("Kendall tau-b by %s", if (has_pcapp) "pcaPP::cor.fk" else
    "stats::cor (install pcaPP for speed; same values)")

# ---- Section 3: DFG replication ------------------------------------------------------------------
say("\n== Section 3. Do the eight DFGs replicate?")
X3 <- N[, pvm3_all]
d3 <- nuc$donor_id[pvm3_all]
screens <- list()
for (tr in c("tier1", "tier2")) {
    th <- SET[[tr]]
    screens[[tr]] <- quiet_msgs(cepip_screen(X3, goi = SET$goi, donor = d3, min_pct = 0,
                                             min_mean = 0, min_donor_cells = SET$donor_floor,
                                             rs = th[["rs"]], tau = th[["tau"]], verbose = TRUE))
}
cells3 <- screens$tier1$cells
exp3 <- nuc$cell_id[pvm3_all][d3 %in% floor_donors]
gr("GR-11-screen_cells", if (setequal(cells3, exp3) && setequal(screens$tier2$cells, exp3)) "PASS" else "FAIL",
   sprintf("screen cells %d | Micro-PVM_3 nuclei of the %d floor donors %d", length(cells3),
           length(floor_donors), length(exp3)))
xg <- as.numeric(N[SET$goi, cells3])
dfg <- data.frame(DFG = SET$dfg8, stringsAsFactors = FALSE)
dfg$rs <- vapply(SET$dfg8, function(g) stats::cor(xg, as.numeric(N[g, cells3]), method = "spearman"),
                 numeric(1))
dfg$tau_b <- vapply(SET$dfg8, function(g) ktau(xg, as.numeric(N[g, cells3])), numeric(1))
for (tr in c("tier1", "tier2")) {
    p <- screens[[tr]]$partners
    dfg[[paste0("passes_", tr)]] <- dfg$DFG %in% p$partner[p$GOI == SET$goi & p$direction == "positive"]
}
flag_ok <- all(dfg$passes_tier1 == (dfg$rs >= SET$tier1[["rs"]] - 1e-12 & dfg$tau_b >= SET$tier1[["tau"]] - 1e-12)) &&
    all(dfg$passes_tier2 == (dfg$rs >= SET$tier2[["rs"]] - 1e-12 & dfg$tau_b >= SET$tier2[["tau"]] - 1e-12))
gr("GR-12-screen_flags", if (flag_ok) "PASS" else "FAIL",
   "the screen's pass flags equal rs and tau_b computed directly on the same cells, at both tiers")
n1 <- sum(dfg$passes_tier1); n2 <- sum(dfg$passes_tier2)
reading3 <- if (n1 >= SET$n_pass) "The module replicates at the SEA-AD thresholds (tier 1)." else
    if (n2 >= SET$n_pass) "The module replicates at the CEP-IP lower thresholds (tier 2)." else
        "The module does not replicate at either tier in this dataset."
print(dfg, row.names = FALSE, digits = 4)
say("Pass tier 1: %d of 8 | pass tier 2: %d of 8 | needed: %d", n1, n2, SET$n_pass)
say("READING (Section 3): %s", reading3)
gate3 <- screens$tier1$gate
say("MT-CO3 gate values: share expressing %.4f | mean over expressing nuclei %.4f",
    gate3$pct_expr[1], gate3$mean_expressing[1])

# ---- Section 4: arms and IP gate --------------------------------------------------------------------
# A cepip error in one arm is recorded (GR-13) and the other arms still run.
run_arm <- function(arm, cells, x_gene, y_genes, cmrf) {
    say("\n== Arm %s: %d nuclei, %d donors", arm, length(cells), length(unique(nuc$donor_id[cells])))
    x <- as.numeric(N[x_gene, cells])
    y <- cepip_composite(dense(y_genes, cells), y_genes)
    r <- tryCatch(quiet_msgs(
        if (cmrf) cepip(x, y, donor = nuc$donor_id[cells], id = nuc$cell_id[cells], cmrf = TRUE,
                        ip_k = SET$ip_k, ip_axis = SET$ip_axis, kappa = SET$kappa, bw = SET$bw,
                        repeats = SET$repeats, workers = WORKERS)
        else cepip(x, y, donor = nuc$donor_id[cells], id = nuc$cell_id[cells], workers = WORKERS)),
        error = function(e) e)
    if (inherits(r, "error")) {
        # HKG8 is the control arm: its failure is reported, and the PVM_3 analysis goes ahead.
        gr(paste0("GR-13-", arm), if (arm == "HKG8") "WARN" else "FAIL",
           paste("cepip stopped:", conditionMessage(r),
                 if (arm == "HKG8") "| HKG8 is reported as not available" else ""))
        return(NULL)
    }
    gr(paste0("GR-13-", arm), "PASS", sprintf("cepip completed on %d cells | status: %s | ip_axis %s",
                                               r$fit$n_cells, r$status, r$settings$ip_axis))
    ip <- r$ip
    # The same-cell control is compared at the GAM step only, so it has no IP guardrail.
    if (cmrf) gr(paste0("GR-14-", arm),
       if (!isTRUE(ip$has_ip) || ip$iprs$tier == "Failed") "WARN"
       else if (ip$iprs$tier %in% c("Very Strong", "Strong") && !grepl("C=0", ip$iprs$flag)) "PASS" else "WARN",
       if (!isTRUE(ip$has_ip)) paste("no IP:", ip$status)
       else sprintf("IP_x %.4f | IP_y %.4f | IPRS %.3f (%s) | flag: %s", ip$ip_x, ip$ip_y,
                    ip$iprs$IPRS, ip$iprs$tier, ip$iprs$flag))
    if (cmrf && !is.null(r$cmrf)) {
        cc <- r$cmrf$counts
        mod_ok <- sum(cc) == length(cells) &&
            sum(r$cmrf$cells$in_module) == cc[["Original + Refined"]] + cc[["Recovered"]]
        gr(paste0("GR-15-", arm), if (mod_ok) "PASS" else "FAIL",
           sprintf("%s | Module %d | arm cells %d",
                   paste(sprintf("%s %d", names(cc), as.integer(cc)), collapse = " | "),
                   sum(r$cmrf$cells$in_module), length(cells)))
        gr(paste0("GR-16-", arm), if (!is.null(r$validate)) "PASS" else "WARN",
           if (!is.null(r$validate)) sprintf("Validation1: k %d | DE %.4f | label agreement %.1f%%",
                                             r$validate$fit$k, r$validate$de, 100 * r$validate$agreement)
           else "Validation1 not available")
    } else if (cmrf) {
        gr(paste0("GR-15-", arm), "WARN", "CMRF not run (the run stopped after the IP step)")
    }
    r
}

# The CMPB 2026 Section 3.12 mosaic comparison on the curve that defines the IP: the share of TREP
# cells above the main GAM curve, post-IP against pre-IP (two-sided Fisher's exact test). CMPB 2026
# splits pre-IP from post-IP on the GOI axis (x < IP_x), so the review uses that split. The split on
# the cDFG axis (y < IP_y, the analysis labels) is reported beside it.
mosaic <- function(res, axis = c("x", "y")) {
    axis <- match.arg(axis)
    na <- list(pre_n = NA, pre_above = NA, post_n = NA, post_above = NA, share_pre = NA,
               share_post = NA, p = NA)
    if (is.null(res) || !isTRUE(res$ip$has_ip)) return(na)
    ip <- res$ip
    e <- ip$ep
    trep <- ip$is_top
    lab <- if (axis == "x") ifelse(e$x < ip$ip_x, "Pre-IP", "Post-IP") else ip$labels
    side <- lab[trep]
    above <- (e$y > e$fitted)[trep]
    pre_n <- sum(side == "Pre-IP"); post_n <- sum(side == "Post-IP")
    pre_a <- sum(above & side == "Pre-IP"); post_a <- sum(above & side == "Post-IP")
    p <- if (pre_n > 0 && post_n > 0)
        stats::fisher.test(matrix(c(post_a, post_n - post_a, pre_a, pre_n - pre_a), 2))$p.value else NA
    list(pre_n = pre_n, pre_above = pre_a, post_n = post_n, post_above = post_a,
         share_pre = if (pre_n) pre_a / pre_n else NA, share_post = if (post_n) post_a / post_n else NA,
         p = p)
}

gate_row <- function(arm, res, n_cells) {
    mo <- mosaic(res, "x"); my <- mosaic(res, "y")
    if (is.null(res)) {
        i <- list(C1 = NA, C2 = NA, C3 = NA, C4 = NA, C5 = NA, IPRS = NA, tier = "not available")
        has_ip <- FALSE; ip_x <- NA; ip_y <- NA
    } else {
        i <- res$ip$iprs; has_ip <- isTRUE(res$ip$has_ip); ip_x <- res$ip$ip_x; ip_y <- res$ip$ip_y
    }
    any_zero <- has_ip && any(unlist(i[c("C1", "C2", "C3", "C4", "C5")]) == 0)
    decision <- if (is.null(res)) "not run: the cepip run failed" else
        if (n_cells < SET$min_arm_cells) "not run: fewer than 50 GOI-positive cells" else
        if (!has_ip) "not run: no IP found" else
            if (i$tier == "Failed") "not run: IPRS Failed" else
                if (i$tier %in% c("Very Strong", "Strong") && !any_zero) "run" else {
                    ok <- isTRUE(mo$share_post > mo$share_pre) && isTRUE(mo$p < SET$review_alpha)
                    if (ok) "run, labelled as conditional on the review" else
                        "not run: review not confirmed"
                }
    data.frame(arm = arm, n_cells = n_cells, has_IP = has_ip, IP_x = ip_x, IP_y = ip_y,
               C1 = i$C1, C2 = i$C2, C3 = i$C3, C4 = i$C4, C5 = i$C5, IPRS = i$IPRS, tier = i$tier,
               any_component_zero = any_zero,
               review_needed = !(i$tier %in% c("Very Strong", "Strong")) || any_zero,
               TREP_pre_IP = mo$pre_n, TREP_pre_IP_above_curve = mo$pre_above,
               TREP_post_IP = mo$post_n, TREP_post_IP_above_curve = mo$post_above,
               share_above_pre = mo$share_pre, share_above_post = mo$share_post,
               fisher_p = mo$p, review_split = "GOI axis (x < IP_x), as CMPB 2026",
               share_above_pre_ysplit = my$share_pre, share_above_post_ysplit = my$share_post,
               fisher_p_ysplit = my$p, diagnosis_contrasts = decision, stringsAsFactors = FALSE)
}

res <- list()
res$PVM_3 <- run_arm("PVM_3", arm_cells$PVM_3, SET$goi, SET$dfg8, cmrf = TRUE)
res$HKG8 <- run_arm("HKG8", arm_cells$HKG8, SET$hkg_goi, SET$hk8, cmrf = TRUE)
res$HKG_control <- run_arm("HKG_control", arm_cells$PVM_3, SET$goi, c(SET$hkg_goi, SET$hk8),
                           cmrf = FALSE)
gate4 <- rbind(gate_row("PVM_3", res$PVM_3, length(arm_cells$PVM_3)),
               gate_row("HKG8", res$HKG8, length(arm_cells$HKG8)))
say("\n== Section 4. IP gate")
print(gate4[, c("arm", "n_cells", "IPRS", "tier", "any_component_zero", "share_above_pre",
                "share_above_post", "fisher_p", "diagnosis_contrasts")], row.names = FALSE, digits = 4)
for (k in seq_len(nrow(gate4)))
    say("READING (Section 4, %s): IPRS %s (%s); diagnosis contrasts %s.", gate4$arm[k],
        format(round(gate4$IPRS[k], 3), nsmall = 3), gate4$tier[k], gate4$diagnosis_contrasts[k])

# ---- Section 6: observed items ------------------------------------------------------------------------
say("\n== Section 6. Reported as observed")
ac <- arm_cells$PVM_3
xa <- as.numeric(N[SET$goi, ac]); ya <- cepip_composite(dense(SET$dfg8, ac), SET$dfg8)
within <- data.frame(item = c("rs MT-CO3 vs cDFG, PVM_3 arm (MT-CO3 > 0)",
                              "tau_b MT-CO3 vs cDFG, PVM_3 arm (MT-CO3 > 0)"),
                     value = c(stats::cor(xa, ya, method = "spearman"), ktau(xa, ya)))
print(within, row.names = FALSE, digits = 4)

agn <- list(); agn_rank <- data.frame()
for (tr in c("tier1", "tier2")) {
    th <- SET[[tr]]
    say("Agnostic screen, %s", tr)
    agn[[tr]] <- quiet_msgs(cepip_screen(X3, goi = NULL, donor = d3, min_donor_cells = SET$donor_floor,
                                         rs = th[["rs"]], tau = th[["tau"]], verbose = TRUE))
    rk <- agn[[tr]]$ranking
    pos <- if (nrow(rk)) match(SET$goi, rk$GOI) else NA
    agn_rank <- rbind(agn_rank, data.frame(
        tier = tr, GOIs_screened = sum(agn[[tr]]$gate$passes), GOIs_with_DFGs = nrow(rk),
        MT_CO3_rank = pos, MT_CO3_n_DFG = if (is.na(pos)) 0L else rk$n_DFG[pos],
        top_GOI = if (nrow(rk)) rk$GOI[1] else NA, top_n_DFG = if (nrow(rk)) rk$n_DFG[1] else NA,
        stringsAsFactors = FALSE))
}
print(agn_rank, row.names = FALSE)

done <- res[!vapply(res, is.null, logical(1))]
pipe <- do.call(rbind, lapply(names(done), function(a) {
    s <- summary(done[[a]]); data.frame(arm = a, item = s$item, value = s$value, stringsAsFactors = FALSE)
}))
de_of <- function(r) if (is.null(r)) NA_real_ else r$fit$de
same_cell <- data.frame(item = c("DE, MT-CO3 vs cDFG (PVM_3 arm)", "DE, MT-CO3 vs cHKG9 (same cells)",
                                 "DE difference (cDFG minus cHKG9)"),
                        value = c(de_of(res$PVM_3), de_of(res$HKG_control),
                                  de_of(res$PVM_3) - de_of(res$HKG_control)))
print(same_cell, row.names = FALSE, digits = 4)

known <- data.frame(); both_map <- data.frame()
if (!is.null(res$PVM_3$cmrf)) {
    g <- cepip_groups(res$PVM_3)
    m <- match(g$id, nuc$cell_id)
    post <- g$group == "PostIP"
    for (cut in SET$mt_cutoffs) known <- rbind(known, data.frame(
        filter = sprintf("%%MT above %g%%", cut),
        cutoff_pct_MT = cut, module_postIP_cells = sum(post),
        module_postIP_removed = sum(nuc$pct_mt[m][post] > cut),
        share_module_postIP_removed = mean(nuc$pct_mt[m][post] > cut),
        other_arm_cells = sum(!post), other_arm_removed = sum(nuc$pct_mt[m][!post] > cut),
        share_other_removed = mean(nuc$pct_mt[m][!post] > cut)))
    # The authors' published QC: nuclei absent from the authors' metadata count as removed. Only
    # nuclei of libraries with the authors' metadata enter this row.
    qa <- nuc$in_author_qc[m]
    ok <- !is.na(qa)
    known <- rbind(known, data.frame(
        filter = sprintf("authors' published QC (%d of %d arm cells with the authors' metadata)",
                         sum(ok), length(ok)),
        cutoff_pct_MT = NA_real_, module_postIP_cells = sum(post & ok),
        module_postIP_removed = sum(post & ok & !qa, na.rm = TRUE),
        share_module_postIP_removed = if (any(post & ok)) mean(!qa[post & ok]) else NA_real_,
        other_arm_cells = sum(!post & ok), other_arm_removed = sum(!post & ok & !qa, na.rm = TRUE),
        share_other_removed = if (any(!post & ok)) mean(!qa[!post & ok]) else NA_real_))
    both_map <- data.frame(item = c("Module post-IP cells", "also Micro-PVM_3 by Hierarchical Mapping",
                                    "share"),
                           value = c(sum(post), sum(nuc$is_pvm3_hm[m][post]),
                                     mean(nuc$is_pvm3_hm[m][post])))
    print(known, row.names = FALSE, digits = 4)
    print(both_map, row.names = FALSE, digits = 4)
}

# Per-donor Module post-IP%: denominator all Micro-PVM_3 nuclei of the donor (MT-CO3 = 0 included);
# donors at or above the floor only. Ages are added by R15.
donor_tab <- function(arm, r) {
    if (is.null(r$cmrf)) return(NULL)
    dp <- cepip_donor(r, n_total = n_pvm3)
    dp <- dp[dp$donor %in% floor_donors, ]
    out <- merge(donor_info, dp, by = "donor")
    tot_ok <- nrow(out) == length(floor_donors) &&
        all(out$n_total == as.integer(n_pvm3[out$donor])) &&
        all(out$pct_PostIP >= 0 & out$pct_PostIP <= 100) &&
        all(abs(out$n_PostIP + out$n_PreIP + out$n_NMPost + out$n_NMPre -
                as.integer(table(factor(nuc$donor_id[arm_cells[[arm]]], levels = out$donor)))) == 0)
    gr(paste0("GR-17-", arm), if (tot_ok) "PASS" else "FAIL",
       sprintf("%d donors | denominator = all Micro-PVM_3 nuclei of the donor | group counts sum to the donor's arm cells | Module post-IP%% from %.2f to %.2f",
               nrow(out), min(out$pct_PostIP), max(out$pct_PostIP)))
    out
}
donors <- list(PVM_3 = donor_tab("PVM_3", res$PVM_3), HKG8 = donor_tab("HKG8", res$HKG8))

mt_max <- max(nuc$pct_mt[ac])
gr("GR-18-pct_MT", if (is.finite(mt_max) && mt_max >= 0 && mt_max <= 100) "PASS" else "FAIL",
   sprintf("PVM_3 arm (no %%MT cut applied): median %%MT %.3f | max %.3f | nuclei above 1%%, 5%%, and 10%%: %d, %d, %d",
           stats::median(nuc$pct_mt[ac]), mt_max, sum(nuc$pct_mt[ac] > 1), sum(nuc$pct_mt[ac] > 5),
           sum(nuc$pct_mt[ac] > 10)))


deg <- NULL
if (!is.null(res$PVM_3$cmrf)) {
    say("DEGs, PVM_3 arm (cepip_deg defaults, the S8 thresholds)")
    E <- N[, ac]
    E <- E[Matrix::rowSums(E > 0) >= 10, ]
    deg <- quiet_msgs(cepip_deg(res$PVM_3, E, verbose = TRUE))
}

# ---- Save -----------------------------------------------------------------------------------------------------
say("\n== Saving")
arm_dirs <- c(PVM_3 = "R14_PVM_3", HKG8 = "R14_HKG8", HKG_control = "R14_HKG_control")
quiet_msgs({
    if (!is.null(res$PVM_3))
        cepip_save(res$PVM_3, file.path(ROOT, arm_dirs[["PVM_3"]]), arm = arm_dirs[["PVM_3"]],
                   extras = list(all = list(screen = screens$tier1, donor = donors$PVM_3, deg = deg)))
    if (!is.null(res$HKG8))
        cepip_save(res$HKG8, file.path(ROOT, arm_dirs[["HKG8"]]), arm = arm_dirs[["HKG8"]],
                   extras = list(all = list(donor = donors$HKG8)))
    if (!is.null(res$HKG_control))
        cepip_save(res$HKG_control, file.path(ROOT, arm_dirs[["HKG_control"]]),
                   arm = arm_dirs[["HKG_control"]])
})

reading_tab <- data.frame(section = c("3", "4", "4"),
                          reading = c(reading3,
                                      sprintf("%s: IPRS %s (%s); diagnosis contrasts %s.", gate4$arm,
                                              format(round(gate4$IPRS, 3), nsmall = 3), gate4$tier,
                                              gate4$diagnosis_contrasts)),
                          stringsAsFactors = FALSE)
settings <- data.frame(setting = names(SET), value = vapply(SET, function(v) paste(v, collapse = ", "),
                                                            character(1)), stringsAsFactors = FALSE)
session <- data.frame(item = c("started", "finished", "R", "cepip", "cepip RemoteSha", "mgcv", "Matrix",
                               "openxlsx", "pcaPP", "BLAS", "LAPACK", "workers", "readings commit"),
                      value = c(format(t0), format(Sys.time()), as.character(getRversion()),
                                as.character(packageVersion("cepip")), if (is.null(sha)) NA else sha,
                                as.character(packageVersion("mgcv")), as.character(packageVersion("Matrix")),
                                as.character(packageVersion("openxlsx")),
                                if (has_pcapp) as.character(packageVersion("pcaPP")) else "not installed",
                                extSoftVersion()[["BLAS"]], La_library(), WORKERS, READINGS_COMMIT),
                      stringsAsFactors = FALSE)
sheets <- list("01_Readings" = reading_tab, "02_S3_DFGs" = dfg, "03_S3_MTCO3_gate" = gate3,
               "04_S4_IP_gate" = gate4, "05_S6_within_positive" = within,
               "06_S6_agnostic_rank" = agn_rank,
               "07_S6_agnostic_tier1" = agn$tier1$ranking, "08_S6_agnostic_tier2" = agn$tier2$ranking,
               "09_S6_pipeline" = pipe, "10_S6_same_cell" = same_cell,
               "11_S6_known_filter" = known, "12_S6_both_mappings" = both_map,
               "13_Donors_PVM_3" = donors$PVM_3, "14_Donors_HKG8" = donors$HKG8,
               "15_Settings" = settings, "16_Session" = session)
sheets <- sheets[!vapply(sheets, function(d) is.null(d) || !nrow(d), logical(1))]
wb <- openxlsx::createWorkbook()
bold <- openxlsx::createStyle(textDecoration = "bold")
for (nm in names(sheets)) {
    openxlsx::addWorksheet(wb, nm)
    openxlsx::writeData(wb, nm, sheets[[nm]], headerStyle = bold)
    openxlsx::freezePane(wb, nm, firstRow = TRUE)
}
xlsx_path <- file.path(OUT, "R14_Pineda_readings.xlsx")
rds_path <- file.path(OUT, "R14_Pineda_results.rds")
openxlsx::saveWorkbook(wb, xlsx_path, overwrite = TRUE)
saveRDS(list(res = res, donors = donors, gate4 = gate4, reading3 = reading3, dfg = dfg,
             n_pvm3 = n_pvm3, floor_donors = floor_donors, settings = SET), rds_path)

gr("GR-19-R_warnings", if (length(GR$warnings)) "WARN" else "PASS",
   if (length(GR$warnings)) sprintf("%d R warning(s): %s", length(GR$warnings),
                                    paste(unique(GR$warnings), collapse = " || "))
   else "no R warnings during the run")
expected_out <- c(xlsx_path, rds_path,
                  file.path(ROOT, arm_dirs[names(res)[!vapply(res, is.null, logical(1))]], "00_settings.txt"))
gr("GR-20-outputs", if (all(file.exists(expected_out))) "PASS" else "FAIL",
   sprintf("%d of %d expected files written", sum(file.exists(expected_out)), length(expected_out)))
say("\nGuardrail tally: %d PASS / %d WARN / %d FAIL", gr_count("PASS"), gr_count("WARN"), gr_count("FAIL"))
say("Done in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))
invisible(TRUE)
}

# The console output also goes to the log. R warnings are logged as they occur and counted in GR-19.
# The guardrail ledger is written whether or not the run finishes.
writeLines(paste("R14 started", format(Sys.time()), "- this file is removed when the run finishes."),
           INCOMPLETE)
# Results of an earlier run are removed, so a failed run cannot leave old readings in place.
old_out <- file.path(OUT, c("R14_Pineda_readings.xlsx", "R14_Pineda_results.rds"))
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
