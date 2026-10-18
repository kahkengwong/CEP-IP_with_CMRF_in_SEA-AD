# =============================================================================
# R_config.R
# SEA-AD CEP-IP/CMRF, settings for the three Module 1 arms (Parts 9 onwards)
#   Module 1: GOI MT-CO3 and eight DFGs (MT-ATP6, MT-CO2, MT-ND3, MT-ND4,
#   MT-CYB, MT-ND1, MT-ND2, MT-CO1). Sets the export paths for the arm named in
#   CEP_ARM (Mic13_PVM3, PVM3_all or PVM1_all; Mic13_PVM3 when unset), the
#   regions, constants, ordinal maps and loaders, and checks them against the
#   arm's Part8_export_record.csv.
#
# Usage: sourced by the analysis scripts, not run directly.
#
# Outputs: none.
# =============================================================================

R_EXPORT_ROOT <- paste0(Sys.getenv("SEAAD_ROOT", getwd()), "/R_export_CP10K")

# == Arm selection ============================================================
# Module 1 arms. The caller sets CEP_ARM before sourcing; unset gives the
# primary arm. The HKG8 control arm uses R_config_HKG8.R.
ARM_TAG_DEFAULT <- "Mic13_PVM3"
ARM_TAG_ALLOWED <- c("Mic13_PVM3", "PVM3_all", "PVM1_all")

# Braces keep the top-level if/else parseable in a sourced file.
ARM_TAG <- if (exists("CEP_ARM") && !is.null(CEP_ARM)) {
  CEP_ARM
} else {
  ARM_TAG_DEFAULT
}

if (!ARM_TAG %in% ARM_TAG_ALLOWED) {
  stop(sprintf(paste0("R_config.R: arm '%s' is not one of the Module 1 arms ",
                      "(%s).\n  The HKG8 control is defined by ",
                      "R_config_HKG8.R."),
               ARM_TAG, paste(ARM_TAG_ALLOWED, collapse = ", ")))
}

ARM_DIR    <- file.path(R_EXPORT_ROOT, ARM_TAG)
NEOC_DIR   <- file.path(ARM_DIR, "neocortex")
NONNC_DIR  <- file.path(ARM_DIR, "nonneocortex")
SHARED_DIR <- file.path(R_EXPORT_ROOT, "shared")

# == GOI and DFG definitions ==================================================
# One gene set for both NeoC and NonNeoC.
GOI <- "MT-CO3"

DFG <- c("MT-ATP6", "MT-CO2", "MT-ND3", "MT-ND4",
         "MT-CYB",  "MT-ND1", "MT-ND2", "MT-CO1")
DFG_NEOC  <- DFG   # same set
DFG_NONNC <- DFG   # same set

# Full Module 1 gene set (GOI + DFGs = 9 genes)
MOD1_GENES <- c(GOI, DFG)

# Genes required to be > 0 for a cell to be selected (Part 8 SELECT_GENES).
SELECT_GENES <- c("MT-CO3")

# == Labels for output file naming ============================================
# RUN_SUFFIX: none for the primary arm, its own suffix for each all-cells arm.
DFG_COMPOSITE_LABEL <- "cDFG"
RUN_SUFFIX <- switch(ARM_TAG,
                     "Mic13_PVM3" = "",
                     "PVM3_all"   = "_PVM3all",
                     "PVM1_all"   = "_PVM1all",
                     paste0("_", ARM_TAG))

# == Region definitions =======================================================
NEOC_REGIONS  <- c("DFC", "MTG", "AnG", "STG", "V1C", "ITG")
NONNC_REGIONS <- c("HIP", "MEC", "LEC", "FI")

# == Supertype of primary interest ============================================
# PVM1_all selects Micro-PVM_1; the other two arms Micro-PVM_3-SEAAD.
PRIMARY_STYPE <- if (ARM_TAG == "PVM1_all") "Micro-PVM_1" else "Micro-PVM_3-SEAAD"

# == Normalisation constants ==================================================
# adata.X is log1p(UMIs / total_counts * 1e4); Part 8 exports it unchanged.
NORM_TARGET_SUM <- 10000
NORM_CEILING    <- log1p(NORM_TARGET_SUM)   # 9.2104
NORM_TAG        <- "CP10K"

# == Housekeeping genes (negative-control panel) ==============================
# In the HKG8 control, ACTB is the GOI and the other eight form cHKG8.
HKG_GENES <- c("ACTB", "EEF1A1", "RPL13A", "B2M", "YWHAZ",
               "FTL",  "PSAP", "RAB7A", "HNRNPA2B1")

# == Biological validation gene panels ========================================
# Part 10 appends these genes to the pooled SKCS scores when BIO_GENES exists.
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
# check.names = FALSE keeps hyphenated gene columns intact (MT-CO3, not MT.CO3).
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
# Selected cells are in the arm tree; the dataset-wide matrix is under SHARED_DIR.
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

# == Column name constants (as written by Part 8) =============================
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

cat("R_config.R loaded (Module 1: MT-CO3 + 8 DFGs, \u22658/10 regions).\n")
cat("  ARM_TAG   :", ARM_TAG, "\n")
cat("  ARM_DIR   :", ARM_DIR, "\n")
cat("  NEOC_DIR  :", NEOC_DIR, "\n")
cat("  NONNC_DIR :", NONNC_DIR, "\n")
cat("  SHARED_DIR:", SHARED_DIR, "\n")
cat("  GOI       :", GOI, "\n")
cat("  DFG_NEOC  :", paste(DFG_NEOC, collapse=", "), "\n")
cat("  DFG_NONNC :", paste(DFG_NONNC, collapse=", "), "\n")
cat("  Composite :", DFG_COMPOSITE_LABEL, "\n")
cat("  RUN_SUFFIX:", if (nchar(RUN_SUFFIX)) RUN_SUFFIX else "(none)", "\n")
cat("  Supertype :", PRIMARY_STYPE, "\n")
cat("  Mod1 genes:", paste(MOD1_GENES, collapse=", "),
    "(", length(MOD1_GENES), ")\n")
cat("  HKG_GENES :", paste(HKG_GENES, collapse=", "),
    "(", length(HKG_GENES), ")\n")
validate_export_record()
