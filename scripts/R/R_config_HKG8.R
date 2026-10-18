# =============================================================================
# R_config_HKG8.R
# SEA-AD CEP-IP/CMRF, settings for the HKG8 control arm (Parts 9 onwards)
#   Housekeeping control: GOI ACTB and eight HKGs as the composite (EEF1A1,
#   RPL13A, B2M, YWHAZ, FTL, PSAP, RAB7A, HNRNPA2B1), on the primary arm's cells
#   with ACTB > 0 (arm HKG8_matched). Used in place of R_config.R for this arm;
#   checked against the arm's Part8_export_record.csv.
#
# Usage: sourced by the analysis scripts, not run directly.
#
# Outputs: none.
# =============================================================================

R_EXPORT_ROOT <- paste0(Sys.getenv("SEAAD_ROOT", getwd()), "/R_export_CP10K")

# This file defines one arm; a different CEP_ARM from the caller stops the run.
ARM_TAG <- "HKG8_matched"

if (exists("CEP_ARM") && !is.null(CEP_ARM) && !identical(CEP_ARM, ARM_TAG)) {
  stop(sprintf(paste0("R_config_HKG8.R defines arm '%s' but the caller asked ",
                      "for '%s'.\n  Module 1 arms are defined by R_config.R."),
               ARM_TAG, CEP_ARM))
}

# == Dataset paths (same cells as the primary arm, HKG panel) =================
ARM_DIR    <- file.path(R_EXPORT_ROOT, ARM_TAG)
NEOC_DIR   <- file.path(ARM_DIR, "neocortex")
NONNC_DIR  <- file.path(ARM_DIR, "nonneocortex")
SHARED_DIR <- file.path(R_EXPORT_ROOT, "shared")

# == GOI and DFG definitions ==================================================
# ACTB as GOI; the other eight HKGs form the composite.
GOI <- "ACTB"

DFG <- c("EEF1A1", "RPL13A", "B2M", "YWHAZ",
         "FTL",    "PSAP",   "RAB7A", "HNRNPA2B1")
DFG_NEOC  <- DFG   # same set
DFG_NONNC <- DFG   # same set

# Full gene set (GOI + DFGs = 9 genes, the same as HKG_GENES)
MOD1_GENES <- c(GOI, DFG)

# Genes required to be > 0 for a cell to be selected (Part 8 SELECT_GENES).
SELECT_GENES <- c("MT-CO3", "ACTB")

# == Labels for output file naming ============================================
DFG_COMPOSITE_LABEL <- "cHKG8"
RUN_SUFFIX          <- "_HKG8"

# == Region definitions =======================================================
NEOC_REGIONS  <- c("DFC", "MTG", "AnG", "STG", "V1C", "ITG")
NONNC_REGIONS <- c("HIP", "MEC", "LEC", "FI")

# == Supertype of primary interest ============================================
PRIMARY_STYPE <- "Micro-PVM_3-SEAAD"

# == Normalisation constants ==================================================
NORM_TARGET_SUM <- 10000
NORM_CEILING    <- log1p(NORM_TARGET_SUM)   # 9.2104
NORM_TAG        <- "CP10K"

# == Housekeeping panel =======================================================
# GOI + DFG equals HKG_GENES in this arm.
HKG_GENES <- c("ACTB", "EEF1A1", "RPL13A", "B2M", "YWHAZ",
               "FTL",  "PSAP", "RAB7A", "HNRNPA2B1")

# == Biological validation gene panels ========================================
GLYC_GENES   <- c("HK2", "PKM")
DAM_GENES    <- c("TREM2", "APOE", "CLEC7A", "CST7", "TSPO",
                   "TYROBP", "SPP1",  "LPL",  "CD9", "GPNMB")
INFLAM_GENES <- c("NLRP3", "CASP1", "IL1B", "PYCARD")
MIC13UP_GENES <- c("PTPRG", "GLDN",  "MYO1E", "DPYD",    "PPARG",
                    "ALCAM", "CTSD",  "MSR1",  "ITGAX",   "OLR1",
                    "LIPA",  "ASAH1", "STARD13")
MIC13DN_GENES <- c("CX3CR1", "FOXP2", "PRKCA", "ITPR2", "TMEM163")
BIO_GENES    <- unique(c(GLYC_GENES, DAM_GENES, INFLAM_GENES,
                          MIC13UP_GENES, MIC13DN_GENES))


# == Ordinal maps =============================================================
# Observed levels are strings, so they are mapped explicitly.
BRAAK_MAP <- c("Braak 0" = 0, "Braak I" = 1, "Braak II" = 2, "Braak III" = 3,
               "Braak IV" = 4, "Braak V" = 5, "Braak VI" = 6)
CERAD_MAP <- c("Absent" = 0, "Sparse" = 1, "Moderate" = 2, "Frequent" = 3)
THAL_MAP  <- c("Thal 0" = 0, "Thal 1" = 1, "Thal 2" = 2,
               "Thal 3" = 3, "Thal 4" = 4, "Thal 5" = 5)
CAA_MAP   <- c("Not identified" = 0, "Mild" = 1, "Moderate" = 2, "Severe" = 3)

# Apply a map to a character vector. Unmapped values stop the run rather than
# becoming NA.
map_ordinal <- function(x, map, label = "ordinal") {
  x   <- as.character(x)
  out <- unname(map[x])
  bad <- unique(x[is.na(out) & !is.na(x)])
  if (length(bad)) {
    stop(sprintf("%s: %d value(s) not in the map: %s",
                 label, length(bad), paste(bad, collapse = ", ")))
  }
  out
}

# == Validate against Part 8's export record ==================================
# Stops the run if the settings above disagree with Part8_export_record.csv.
validate_export_record <- function(arm_dir = ARM_DIR) {
  rec_path <- file.path(arm_dir, "Part8_export_record.csv")
  if (!file.exists(rec_path)) {
    stop(paste("Export record not found:", rec_path,
               "- run Part 8 for this arm first."))
  }
  rec <- read.csv(rec_path, stringsAsFactors = FALSE, check.names = FALSE)
  get_field <- function(f) {
    v <- rec$value[rec$field == f]
    if (!length(v)) stop(paste("Export record has no field:", f))
    as.character(v[1])
  }
  chk <- function(field, expected, label) {
    actual <- get_field(field)
    if (!identical(actual, expected)) {
      stop(sprintf(paste0("Config/export mismatch on %s.\n",
                          "  config: %s\n  export: %s\n",
                          "  Rerun Part 8 for arm '%s', or correct the config."),
                   label, expected, actual, ARM_TAG))
    }
  }
  chk("arm_tag",     ARM_TAG,                     "ARM_TAG")
  chk("goi",         GOI,                         "GOI")
  chk("dfg",         paste(DFG, collapse = ";"),  "DFG")
  chk("panel_label", DFG_COMPOSITE_LABEL,         "DFG_COMPOSITE_LABEL")
  chk("run_suffix",  RUN_SUFFIX,                  "RUN_SUFFIX")

  fmt_map <- function(m) paste(paste0(names(m), "=", m), collapse = ";")
  for (nm in c("braak", "cerad", "thal", "caa")) {
    chk(paste0(nm, "_map"),
        fmt_map(get(paste0(toupper(nm), "_MAP"))),
        paste0(toupper(nm), "_MAP"))
  }

  cat("  Export record validated:", basename(rec_path), "\n")
  cat("    h5ad    :", get_field("h5ad"), "\n")
  cat("    release :", get_field("release"), "\n")
  cat("    selected:", get_field("n_selected_neoc"), "NeoC /",
      get_field("n_selected_nonneoc"), "NonNeoC\n")
  invisible(rec)
}

# == Function to load GAM data ================================================
# check.names = FALSE keeps hyphenated gene columns intact.
load_gam_data <- function(dataset_dir, region, supertype = PRIMARY_STYPE) {
  fname <- paste0(region, "_",
                  gsub("-", "_", gsub("/", "_", supertype)),
                  "_selected.csv")
  fpath <- file.path(dataset_dir, "gam_data", fname)
  if (!file.exists(fpath)) {
    stop(paste("GAM data file not found:", fpath))
  }
  df <- read.csv(fpath, stringsAsFactors = FALSE, check.names = FALSE)
  rownames(df) <- df$barcode
  return(df)
}

# == Function to load raw count matrix ========================================
load_raw_matrix <- function(dataset_dir, selected_only = TRUE) {
  count_dir <- if (selected_only) {
    file.path(dataset_dir, "counts", "raw_selected")
  } else {
    file.path(SHARED_DIR, basename(dataset_dir), "counts", "raw_dataset")
  }
  mtx_path  <- file.path(count_dir, "matrix.mtx")
  bar_path  <- file.path(count_dir, "matrix_barcodes.tsv")
  feat_path <- file.path(count_dir, "matrix_features.tsv")
  library(Matrix)
  X        <- readMM(mtx_path)
  barcodes <- readLines(bar_path)
  features <- read.table(feat_path, sep = "\t", header = FALSE,
                          stringsAsFactors = FALSE)
  rownames(X) <- barcodes
  colnames(X) <- features[, 2]   # gene symbols in column 2
  return(list(X = X, barcodes = barcodes, features = features))
}

# == Function to load cell metadata ===========================================
load_metadata <- function(dataset_dir, selected_only = TRUE) {
  fname <- if (selected_only) "obs_selected.csv" else "obs_dataset.csv"
  df    <- read.csv(file.path(dataset_dir, "metadata", fname),
                    stringsAsFactors = FALSE, check.names = FALSE)
  rownames(df) <- df$barcode
  return(df)
}

# == Column name constants ====================================================
REGION_COL  <- "Brain Region"
STYPE_COL   <- "Supertype"
DONOR_COL   <- "Donor ID"
ADNC_COL    <- "Overall AD neuropathological Change"
BRAAK_COL   <- "Braak"

# == GAM fitting parameters ===================================================
GAM_GAMMA    <- 1.5
GAM_METHOD   <- "REML"
GAM_SELECT   <- TRUE
GAM_BS       <- "tp"
GAM_MAX_K    <- 10
PRSS_N_ITER  <- 100
PRSS_EARLY_STOP <- 20

cat("R_config_HKG8.R loaded (HKG specificity control: ACTB + 8 cHKGs).\n")
cat("  ARM_TAG   :", ARM_TAG, "\n")
cat("  ARM_DIR   :", ARM_DIR, "\n")
cat("  NEOC_DIR  :", NEOC_DIR, "\n")
cat("  NONNC_DIR :", NONNC_DIR, "\n")
cat("  SHARED_DIR:", SHARED_DIR, "\n")
cat("  GOI       :", GOI, "\n")
cat("  DFG       :", paste(DFG, collapse=", "), "\n")
cat("  DFG_LABEL :", DFG_COMPOSITE_LABEL, "\n")
cat("  RUN_SUFFIX:", RUN_SUFFIX, "\n")
cat("  Selected  :", paste(paste0(SELECT_GENES, ">0"), collapse=" AND "), "\n")
cat("  Mod1 genes:", paste(MOD1_GENES, collapse=", "),
    "(", length(MOD1_GENES), ")\n")
validate_export_record()
