# =============================================================================
# Part10h_CP10K_Prognosis_code.r
# SEA-AD CEP-IP/CMRF, Part 10h, donor-level clinical correlation of the CMRF
#   groups (Figure 9A, Supplementary Figure S14, Table S12)
#   For each arm and dataset: takes the Pre-IP/Post-IP labels from Validation 1
#   (or the cDFG rule at IP_y), computes donor percentages of the module and
#   non-module groups, and tests them against the clinical variables
#   (Spearman, Steiger, Jonckheere-Terpstra, nested F) and in multivariate
#   models with the covariate sets cov4, cov7, covM, covT, covP and covE.
#
# Usage (from the SEA-AD folder, or set SEAAD_ROOT), after
#   Part10e_CP10K_Validation1_code.r
#   source("Part10h_CP10K_Prognosis_code.r")                         # all arms
#   CEP_ARM <- "PVM1_all"; source("Part10h_CP10K_Prognosis_code.r")  # one arm
#   Needs CEP_arms.R, R_config.R and R_config_HKG8.R in the same folder.
#
# Outputs (Part10h_Prognosis/Part10h_Prognosis_<arm>/): figures (PDF, PNG),
#   workbook, results RDS (every table and figure), guardrail ledger and logs.
# =============================================================================

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(writexl)
library(patchwork)

# The seed is set inside run_one_arm, after the arm's config is sourced.

if (!exists("CEP_HOME")) CEP_HOME <- Sys.getenv("SEAAD_ROOT", getwd())
.p10h_harness <- file.path(CEP_HOME, "CEP_arms.R")
if (!file.exists(.p10h_harness))
    stop("CEP_arms.R not found at: ", .p10h_harness,
         "\n  It holds the arm table shared by Parts 9 onwards.", call. = FALSE)
source(.p10h_harness)

PART_TAG     <- "Part10h"
PART_DIR_TAG <- "Part10h_Prognosis"

# == IP_BAND (retired; must stay 0) ===========================================
IP_BAND <- 0.00
if (!is.finite(IP_BAND) || IP_BAND < 0)
    stop("IP_BAND must be a finite, non-negative number; got: ", IP_BAND,
         call. = FALSE)
if (IP_BAND > 0)
    stop("IP_BAND is set to ", IP_BAND, " but the exclusion mechanism it used ",
         "to control no longer exists.\n  The u_IP split was replaced with ",
         "the inherited cDFG label at IP_y, so no cell is assigned\n  to an ",
         "\"IP zone\" and none is excluded. To reinstate a band, define it on ",
         "the cDFG axis\n  around IP_y.", call. = FALSE)

# == CMRF classes =============================================================
CMRF_CLASS_LEVS <- c("Original + Refined", "Recovered", "Original only",
                     "Excluded")
MODULE_CLASSES  <- c("Original + Refined", "Recovered")

# Non-module comparison pool: cells of this supertype and donor that are not in
# the final module (Excluded and Original only).
NONMODULE_CLASSES <- c("Excluded", "Original only")

# The CMRF band has 95% coverage of the SKCS-core, so about 5% of the core
# (Original only) is expected outside it; checked per arm.
CMRF_MISS_NOMINAL <- 0.05
CMRF_MISS_TOL     <- 0.05   # generous: this flags a different residual, not noise

# == Arm groups ===============================================================
ARM_GROUPS <- list(
    CMRF_core_PreIP  = list(classes = MODULE_CLASSES, arms = "Pre-IP",
                            label = "Module Pre-IP"),
    CMRF_core_PostIP = list(classes = MODULE_CLASSES, arms = "Post-IP",
                            label = "Module Post-IP"),
    # The unsplit module; the restricted model in the nested F test
    # (pct_module_all = pct_pre + pct_post).
    CMRF_core_all    = list(classes = MODULE_CLASSES, arms = NULL,
                            label = "Module (all, no IP split)"),
    # All classified cells, split at IP_y only (no CMRF);
    # pct_any_post = pct_post + pct_excl_post exactly.
    AllCells_PostIP  = list(classes = CMRF_CLASS_LEVS, arms = "Post-IP",
                            label = "All cells Post-IP (no CMRF)"),
    OR_PreIP         = list(classes = "Original + Refined", arms = "Pre-IP",
                            label = "Orig+Ref Pre-IP"),
    OR_PostIP        = list(classes = "Original + Refined", arms = "Post-IP",
                            label = "Orig+Ref Post-IP"),
    Recov            = list(classes = "Recovered", arms = NULL,
                            label = "Recovered (all)"),
    # Original only, split by the cDFG rule at IP_y.
    Ori_only_PreIP   = list(classes = "Original only", arms = "Pre-IP",
                            label = "Original only Pre-IP"),
    Ori_only_PostIP  = list(classes = "Original only", arms = "Post-IP",
                            label = "Original only Post-IP"),
    Ori_only         = list(classes = "Original only", arms = NULL,
                            label = "Original only (all)"),
    # Sensitivity arm: the module with Original-only cells put back.
    ModPlusOri_PreIP = list(classes = c(MODULE_CLASSES, "Original only"),
                            arms = "Pre-IP",
                            label = "Module + Ori-only Pre-IP"),
    ModPlusOri_PostIP = list(classes = c(MODULE_CLASSES, "Original only"),
                            arms = "Post-IP",
                            label = "Module + Ori-only Post-IP"),
    # Non-module cells; their arm comes from the cDFG rule at IP_y.
    Excl_PreIP       = list(classes = NONMODULE_CLASSES, arms = "Pre-IP",
                            label = "Non-module Pre-IP"),
    Excl_PostIP      = list(classes = NONMODULE_CLASSES, arms = "Post-IP",
                            label = "Non-module Post-IP"),
    Excl_all         = list(classes = NONMODULE_CLASSES, arms = NULL,
                            label = "Non-module (all)")
)
# The predictors carried into the multivariate models and the Steiger
# contrasts, and the prop_df column each lands in.
PRED_GROUPS <- c(pct_post      = "CMRF_core_PostIP",
                 pct_pre       = "CMRF_core_PreIP",
                 # The unsplit module.
                 pct_module_all = "CMRF_core_all",
                 # The no-CMRF baseline.
                 pct_any_post   = "AllCells_PostIP",
                 pct_excl_post = "Excl_PostIP",
                 pct_excl_pre  = "Excl_PreIP",
                 pct_excl      = "Excl_all",
                 # Sensitivity arms.
                 pct_mpo_post  = "ModPlusOri_PostIP",
                 pct_mpo_pre   = "ModPlusOri_PreIP")
PRED_LABEL <- c(pct_post = "Module Post-IP", pct_pre = "Module Pre-IP",
                pct_module_all = "Module (all, no IP split)",
                pct_any_post = "All cells Post-IP (no CMRF)",
                pct_excl_post = "Non-module Post-IP",
                pct_excl_pre = "Non-module Pre-IP",
                pct_excl = "Non-module (all)",
                pct_mpo_post = "Module + Ori-only Post-IP",
                pct_mpo_pre  = "Module + Ori-only Pre-IP")
PRED_COL <- c(pct_post = "#B2182B", pct_pre = "#2166AC",
              pct_module_all = "#762A83",
              pct_any_post = "#1B7837",
              pct_excl_post = "#E08214", pct_excl_pre = "#8073AC",
              pct_excl = "#4D4D4D",
              pct_mpo_post = "#D6604D", pct_mpo_pre = "#4393C3")

# == Delta-R2 figure ==========================================================
# Bar order per outcome; must match PRED_LABEL (checked by a guardrail).
MV_ARM_ORDER <- c("Non-module Pre-IP", "Module Pre-IP",
                  "Non-module Post-IP", "Module Post-IP")
MV_ARM_FILL <- c("Non-module Pre-IP"  = "#A8C8E0",   # light blue
                 "Module Pre-IP"      = "#5A7FA6",   # light navy
                 "Non-module Post-IP" = "#B08968",   # dark tan
                 "Module Post-IP"     = "#2D2D2D")   # darkest gray
MV_ARM_BORDER <- c("Non-module Pre-IP"  = "#6E9BBF",
                   "Module Pre-IP"      = "#2F4F70",
                   "Non-module Post-IP" = "#7F5F44",
                   "Module Post-IP"     = "#000000")

# Fixed clinical order: cognition, neuropathology, genetic risk, then the
# demographic and tissue-quality variables.
MV_OUTCOME_ORDER <- c(
    "Cognitive status", "MMSE score", "CASI score",
    "CERAD score", "Braak stage", "Thal phase", "ADNC rank", "CAA score",
    "APOE e4 dose",
    "LATE stage", "Lewy body stage", "Atherosclerosis", "Arteriolosclerosis",
    "Microinfarcts",
    "Age at death", "Sex", "Brain pH", "PMI", "RIN", "Fresh brain weight",
    "Years of education", "Months from last MMSE",
    "Mito fraction (donor mean)", "UMIs per cell (donor mean)",
    "Genes detected (donor mean)", "Doublet score (donor mean)",
    "Intronic:exonic ratio (donor mean)",
    "Sequencing saturation (donor mean)")

# Outcomes plotted as negative controls (demographic and tissue-quality
# variables the arm proportion is not expected to predict).
MV_CONTROL_OUTCOMES <- c("Age at death", "Sex", "Brain pH", "PMI", "RIN",
                         "Fresh brain weight")
MV_PANEL_LABEL <- c(clinical = "Clinical outcomes",
                    control  = "Tissue-quality and demographic variables (negative controls)")

MV_FIGURE_SETS <- c("cov4", "covM", "covT")

# Arm-term VIF above which a covariate set gets no figure.
MV_VIF_FIGURE_MAX <- 10

MV_COV_TITLE <- list(
    cov4 = "Age + APOE e4 + Sex + PMI",
    cov7 = "Age + APOE e4 + Sex + PMI + RIN + Brain pH + Brain weight",
    covM = "Age + APOE e4 + Sex + PMI + mito fraction + UMIs per cell",
    covT = "the core + the technical block (mito, genes, intronic:exonic, doublet)",
    covP = "the core + co-pathology (LATE, Lewy, atherosclerosis, arteriolosclerosis, microinfarcts)",
    covE = "the core + cognitive reserve (years of education, months since last MMSE)")
MV_COV_NOTE <- list(
    cov4 = paste0("THE PRIMARY MODEL. Quote effect sizes from here: 16.8 ",
                  "donors per predictor, max VIF 1.19-1.33.\nThe bounds are ",
                  "this set (under-adjusted) against covT (over-adjusted); ",
                  "cov7 sits between them and is in the workbook."),
    cov7 = paste0("Partly OVER-adjusted: RIN tracks ADNC (-0.2295) and fresh ",
                  "brain weight more so (-0.3261), so both sit partly on the ",
                  "causal path. Workbook only -- it moves the arm ",
                  "coefficient very little, and covT is the informative ",
                  "upper bound."),
    covM = paste0("Adjusts for the donor-mean mitochondrial fraction ",
                  "and sequencing depth over the same cells that form the ",
                  "proportion denominator.\nAn arm effect surviving here is ",
                  "not a restatement of mitochondrial content. The mito ",
                  "fraction correlates with the arm exposure at rho +0.83 ",
                  "(NeoC) and +0.72 (NonNeoC), so a large attenuation here ",
                  "is expected and is not by itself a failure."),
    covT = paste0("THE OVER-ADJUSTED BOUND. The technical block at once, ",
                  "including the intronic:exonic ambient index. This set ",
                  "dropped UMIs and saturation: with genes detected they are ",
                  "three measures of library depth, and all three together ",
                  "drove the VIF to 14.\nCheck this set's VIF ",
                  "before presenting the panel. An effect surviving from cov4 ",
                  "to here is not a technical artefact."),
    covP = paste0("Co-pathology. Measured |rho| with the arm ",
                  "exposure is <= 0.15 for Lewy and the vascular scores, so ",
                  "this is NOT confounding control -- it answers the separate ",
                  "question of whether the module tracks cognition ",
                  "independently of co-pathology."),
    covE = paste0("Cognitive reserve. Years of education correlates ",
                  "with the arm exposure at only -0.08 (NeoC) and -0.02 ",
                  "(NonNeoC), so it is NOT a confounder in this cohort ",
                  "despite being the expected one. Included because it is ",
                  "expected and cheap, and reported as measured."))

# == JT panel figure ==========================================================
# Tertile fills for continuous panels; severity ramp for stacked ordinal panels.
JT_TERTILE_FILL <- c("T1" = "#E8923A", "T2" = "#C4A35A", "T3" = "#00796B")
JT_ORDINAL_VARS <- c("Thal phase", "CERAD score", "CAA score",
                     "ADNC rank", "Braak stage", "Cognitive status")
.jt_sev_palette <- function(n) {
    all_cols <- c("#F5EDD6", "#EDD9A8", "#EDCB8C", "#E0A070",
                  "#D4886B", "#B8605F", "#8B4049")
    if (!is.finite(n) || n < 1L) return(character(0))
    if (n == 1L) return(all_cols[1])
    if (n <= length(all_cols))
        all_cols[round(seq(1, length(all_cols), length.out = n))]
    else
        grDevices::colorRampPalette(c("#F5EDD6", "#EDCB8C", "#D4886B",
                                      "#8B4049"))(n)
}

# * / ** / *** from a BH-adjusted p, blank when absent or not significant.
.sig_star <- function(p) {
    ifelse(!is.finite(p), "",
    ifelse(p < 0.001, "***",
    ifelse(p < 0.01,  "**",
    ifelse(p < 0.05,  "*", ""))))
}
# Steiger contrasts (dependent correlations sharing the clinical variable);
# a positive Z means the first member has the larger rho.
STEIGER_PAIRS <- list(
    c("pct_post",      "pct_pre"),        # within the module, the original ask
    c("pct_post",      "pct_excl_post"),  # module vs excluded, same arm
    c("pct_pre",       "pct_excl_pre"),
    c("pct_excl_post", "pct_excl_pre"),   # does the split matter off-module
    # Sensitivity arm.
    c("pct_post",      "pct_mpo_post"),
    # Split against unsplit module.
    c("pct_post",      "pct_module_all"),
    # Module against all cells split at the same IP.
    c("pct_post",      "pct_any_post"),
    c("pct_pre",       "pct_mpo_pre")
)

# == Clinical parameters ======================================================
# `expected`: +1 higher group % -> higher value, -1 -> lower value. Recorded
# and printed only, never used to select a result.
CLINICAL_PARAMS <- list(
    list(label = "Thal phase",      col = "Thal_num",      expected = +1L,
         role = "AD pathology", note = "Abeta staging, Arabic-coded 0-5"),
    list(label = "CERAD score",     col = "CERAD_num",     expected = +1L,
         role = "AD pathology",
         note = "Abeta plaque burden; WORD-valued, digit parser fails 84/84"),
    list(label = "CAA score",       col = "CAA_num",       expected = +1L,
         role = "AD pathology",
         note = "Cerebrovascular amyloid; map adopted here, see CAA_MAP"),
    list(label = "ADNC rank",       col = "ADNC_rank",     expected = +1L,
         role = "AD pathology",
         note = "Composite; WORD-valued, digit parser fails 84/84"),
    list(label = "Braak stage",     col = "Braak_num",     expected = +1L,
         role = "AD pathology",
         note = "Tau spread; ROMAN-coded, digit parser fails 82/84"),
    list(label = "APOE e4 dose",    col = "APOE_e4_dose",  expected = +1L,
         role = "Genetic risk", note = "e4 alleles 0/1/2 from APOE Genotype"),
    list(label = "Cognitive status", col = "CogStatus_num", expected = +1L,
         role = "Cognition", note = "No dementia 0, MCI 1, Dementia 2"),
    list(label = "MMSE score",      col = "MMSE_num",      expected = -1L,
         role = "Cognition, inverse AD axis",
         note = "NOT a negative control - MMSE falls with AD"),
    # CASI is a further cognitive outcome. MoCA (15 of 84 donors) is loaded
    # for the completeness report but is not an outcome.
    list(label = "CASI score",      col = "CASI_num",      expected = -1L,
         role = "Cognition, inverse AD axis",
         note = "100-point battery; falls with AD, like MMSE"),
    list(label = "Age at death",    col = "Age_num",       expected = -1L,
         role = "Demographic, inverse AD axis",
         note = "NOT a negative control - age is an AD risk factor")
)
VAR_ORDER_DISPLAY <- vapply(CLINICAL_PARAMS, `[[`, character(1), "label")

# =============================================================================
# Extra donor variables from obs_dataset.csv
# =============================================================================
# The maps below use the observed levels. The literal string "nan" is a
# missing value and becomes NA (NA_TOKENS applies to every ordinal).
NA_TOKENS <- c("nan", "na", "n/a", "none given", "not applicable",
               "unknown", "unclassifiable", "")

# LATE: 57 of 84 donors carry a stage, so this is not a rare co-pathology in
# this cohort. "Unclassifiable" (1 donor) is an NA token, not a stage.
LATE_MAP <- c("not identified" = 0L, "none" = 0L,
              "late stage 1" = 1L, "stage 1" = 1L, "1" = 1L,
              "late stage 2" = 2L, "stage 2" = 2L, "2" = 2L,
              "late stage 3" = 3L, "stage 3" = 3L, "3" = 3L)

# Lewy: both "not identified" spellings are stage 0.
LEWY_MAP <- c("not identified (olfactory bulb not assessed)" = 0L,
              "not identified (olfactory bulb assessed)" = 0L,
              "not identified" = 0L, "none" = 0L,
              "olfactory bulb only" = 1L,
              "amygdala-predominant" = 2L,
              "brainstem-predominant" = 3L,
              "limbic (transitional)" = 4L,
              "neocortical (diffuse)" = 5L)

# Vascular: no donor has a "none" level in this cohort.
VASC_MAP <- c("none" = 0L, "not identified" = 0L,
              "mild" = 1L, "moderate" = 2L, "severe" = 3L)

# == Donor-level variables from obs_dataset.csv ===============================
# agg "const": one value per donor (conflicts are an error).
# agg "mean":  per-cell quantity averaged over the arm's PRIMARY_STYPE cells.
DONOR_EXTRA <- list(
    list(out = "CASI_num",     agg = "const", type = "numeric",
         patterns = c("Last CASI Score")),
    list(out = "MOCA_num",     agg = "const", type = "numeric",
         patterns = c("Last MOCA Score")),
    list(out = "MMSEint_num",  agg = "const", type = "numeric",
         patterns = c("Interval from last MMSE in months")),
    list(out = "CASIint_num",  agg = "const", type = "numeric",
         patterns = c("Interval from last CASI in months")),
    list(out = "EduYears_num", agg = "const", type = "numeric",
         patterns = c("Years of education")),
    list(out = "LATE_num",     agg = "const", type = "ordinal",
         patterns = c("LATE"), map = LATE_MAP),
    list(out = "LewyFull_num", agg = "const", type = "ordinal",
         patterns = c("Highest Lewy Body Disease"), map = LEWY_MAP),
    list(out = "Athero_num",   agg = "const", type = "ordinal",
         patterns = c("Atherosclerosis"), map = VASC_MAP),
    list(out = "Arteriolo_num", agg = "const", type = "ordinal",
         patterns = c("Arteriolosclerosis"), map = VASC_MAP),
    list(out = "Microinf_num", agg = "const", type = "numeric",
         patterns = c("Total microinfarcts in screening sections")),
    list(out = "MitoFrac_num", agg = "mean", type = "numeric",
         patterns = c("Fraction mitochondrial UMIs", "fraction_mito",
                      "pct_mito", "percent.mt", "mito_fraction")),
    list(out = "nUMI_num",     agg = "mean", type = "numeric",
         patterns = c("Number of UMIs", "n_umis", "nCount_RNA",
                      "total_counts", "n_counts")),
    list(out = "nGene_num",    agg = "mean", type = "numeric",
         patterns = c("Genes detected", "n_genes", "nFeature_RNA")),
    list(out = "Doublet_num",  agg = "mean", type = "numeric",
         patterns = c("Doublet score", "doublet_score")),
    list(out = "Satur_num",    agg = "mean", type = "numeric",
         patterns = c("GEX_sequencing_saturation")),
    # Cohort membership; scanned, in no block. The intronic:exonic ratio below
    # is the ambient index (a low ratio is a high-ambient library).
    list(out = "Study_num",    agg = "const", type = "ordinal",
         patterns = c("Primary Study Name"),
         map = c("act" = 0L, "adrc clinical core" = 1L)),
    list(out = "IntronExon_num", agg = "mean", type = "ratio",
         num = "GEX_Reads_mapped_confidently_to_intronic_regions",
         den = "GEX_Reads_mapped_confidently_to_exonic_regions")
)

# == Covariate blocks =========================================================
# A small core plus one block at a time, with the arm coefficient read across
# them.
COV_SET_4 <- c("age", "apoe", "sex", "pmi")
COV_SET_7 <- c("age", "apoe", "sex", "pmi", "rin", "ph", "weight")
COV_SET_M <- c("age", "apoe", "sex", "pmi", "mito", "numi")
# covT: UMIs per cell and saturation left out (with genes detected, three
# measures of library depth).
COV_SET_T <- c("age", "apoe", "sex", "pmi", "mito", "ngene",
               "intronexon", "doublet")
COV_SET_P <- c("age", "apoe", "sex", "pmi", "late", "lewyfull", "athero",
               "arteriolo", "microinf")
COV_SET_E <- c("age", "apoe", "sex", "pmi", "eduyears", "mmseint")

COV_MAP <- list(
    age    = "Age_num",        apoe   = "APOE_e4_dose",
    sex    = "Sex_num",        pmi    = "PMI_num",
    rin    = "RIN_num",        ph     = "BrainPH_num",
    weight = "BrainWeight_num",
    mito   = "MitoFrac_num",   numi   = "nUMI_num",
    ngene  = "nGene_num",      doublet = "Doublet_num",
    intronexon = "IntronExon_num", satur = "Satur_num",
    late   = "LATE_num",       lewyfull = "LewyFull_num",
    athero = "Athero_num",     arteriolo = "Arteriolo_num",
    microinf = "Microinf_num",
    eduyears = "EduYears_num", mmseint = "MMSEint_num",
    casiint  = "CASIint_num", study = "Study_num"
)
COV_LABEL <- c(age = "Age at death", apoe = "APOE e4 dose", sex = "Sex",
               pmi = "PMI", rin = "RIN", ph = "Brain pH",
               weight = "Fresh brain weight",
               mito = "Mito fraction (donor mean)",
               numi = "UMIs per cell (donor mean)",
               ngene = "Genes detected (donor mean)",
               doublet = "Doublet score (donor mean)",
               intronexon = "Intronic:exonic ratio (donor mean)",
               satur = "Sequencing saturation (donor mean)",
               late = "LATE stage", lewyfull = "Lewy body stage",
               athero = "Atherosclerosis", arteriolo = "Arteriolosclerosis",
               microinf = "Microinfarcts",
               eduyears = "Years of education",
               mmseint = "Months from last MMSE",
               casiint = "Months from last CASI",
               study = "Primary study")

# Smallest core coefficient for which a percentage change is reported;
# below it the scan reports the absolute shift only.
MV_MIN_BETA_FOR_PCT <- 0.10

# Leave-one-in scan: each covariate is added alone to the core.
COV_SCAN <- c("rin", "ph", "weight", "mito", "numi", "ngene", "doublet",
              "intronexon", "satur", "late", "lewyfull", "athero",
              "arteriolo", "microinf", "eduyears", "mmseint", "casiint",
              "study")
# casiint is scanned only (69 of 84 donors complete).

# == Predictor head-to-head ===================================================
# Q1: does the IP split earn its place? Nested 1-df F test of
#     y ~ cov + pct_module_all against y ~ cov + pct_pre + pct_post.
# Q2: does the predictor beat a rival (e.g. mito fraction)? Non-nested duel:
#     the verdict comes from the encompassing model y ~ cov + A + B; the
#     Vuong test and a donor-bootstrap CI on the delta-R2 difference are
#     reported. Collinear rivals (fitted-value correlation above
#     PHH_DUEL_FIT_COR) are reported as indistinguishable.
PHH_ENABLE   <- TRUE
PHH_BOOT     <- 1000L      # donor resamples for the delta-R2 difference CI
PHH_DUEL_FIT_COR <- 0.9999 # above this the two models are one model
PHH_SEED     <- 20260805L  # the duel is reproducible run to run
PHH_MIN_N    <- 20L        # same floor as the multivariate
# The duel is bootstrapped on the clinical outcomes only.
PHH_OUTCOMES <- c("Cognitive status", "MMSE score", "CASI score",
                  "CERAD score", "Braak stage", "Thal phase", "ADNC rank",
                  "CAA score")
# Duels: `a` against a rival given by `b_pred` (a PRED_GROUPS entry) or `b_cov`
# (a COV_MAP key); exactly one is set.
# == Restriction tests ========================================================
# The unsplit predictor is the exact sum of the two split ones, so each test is
# a nested 1-df F test; the identity is checked at run time.
PHH_RESTRICTIONS <- list(
    list(id = "ip_split", unsplit = "pct_module_all",
         split = c("pct_pre", "pct_post"),
         question = "does splitting the module at IP_y beat not splitting it"),
    list(id = "cmrf_split", unsplit = "pct_any_post",
         split = c("pct_post", "pct_excl_post"),
         question = paste("does the CMRF module boundary beat splitting every",
                          "cell of the supertype at the same IP_y")))

# == Placebo-threshold profile ================================================
# The split is recomputed at a grid of cDFG thresholds. p_all counts every grid
# threshold that matches or beats IP_y; p_far excludes thresholds within
# PHH_THRESH_EXCL quantile points of IP_y.
PHH_THRESH_ENABLE <- TRUE
PHH_THRESH_GRID   <- seq(0.05, 0.95, by = 0.01)  # quantiles of module cDFG
PHH_THRESH_EXCL   <- 0.10   # quantile half-width excluded from `p_far`

PHH_DUELS <- list(
    list(id = "module_vs_mito", a = "pct_post", b_cov = "mito",
         question = paste("does the module beat the donor's raw",
                          "mitochondrial fraction")),
    list(id = "module_vs_unsplit", a = "pct_post", b_pred = "pct_module_all",
         question = paste("does splitting the module at the IP beat not",
                          "splitting it")),
    list(id = "module_vs_nonmodule", a = "pct_post", b_pred = "pct_excl_post",
         question = "does the module boundary beat the cells it excluded"))
# Duel declarations are checked at load.
local({
    # `%||%` is not used (base only from R 4.4.0).
    .ids <- vapply(PHH_DUELS, function(d)
        if (is.null(d$id)) NA_character_ else as.character(d$id), character(1))
    if (anyNA(.ids) || anyDuplicated(.ids))
        stop("PHH_DUELS: every duel needs a unique `id`", call. = FALSE)
    for (d in PHH_DUELS) {
        if (is.null(d$a) || !(d$a %in% names(PRED_GROUPS)))
            stop("PHH_DUELS[", d$id, "]: `a` must name a PRED_GROUPS entry; ",
                 "got: ", if (is.null(d$a)) "NULL" else d$a, call. = FALSE)
        if ((!is.null(d$b_pred)) + (!is.null(d$b_cov)) != 1L)
            stop("PHH_DUELS[", d$id, "]: set exactly one of `b_pred` and ",
                 "`b_cov`", call. = FALSE)
        if (!is.null(d$b_pred) && !(d$b_pred %in% names(PRED_GROUPS)))
            stop("PHH_DUELS[", d$id, "]: `b_pred` must name a PRED_GROUPS ",
                 "entry; got: ", d$b_pred, call. = FALSE)
        if (!is.null(d$b_cov) && is.null(COV_MAP[[d$b_cov]]))
            stop("PHH_DUELS[", d$id, "]: `b_cov` must name a COV_MAP key; ",
                 "got: ", d$b_cov, call. = FALSE)
        if (!is.null(d$b_pred) && identical(d$b_pred, d$a))
            stop("PHH_DUELS[", d$id, "]: a predictor cannot duel itself",
                 call. = FALSE)
    }
    if (!is.finite(PHH_BOOT) || PHH_BOOT < 100L)
        stop("PHH_BOOT must be at least 100; the CI is a percentile interval ",
             "and below that it is noise. Got: ", PHH_BOOT, call. = FALSE)
})

# == Declared ordinal maps ====================================================
# Keys are matched on tolower(trimws(x)); `.map_ordinal()` reports unmatched
# levels.
ADNC_MAP <- c("not ad" = 0L, "no ad" = 0L, "reference" = 0L,
              "low" = 1L, "intermediate" = 2L, "high" = 3L)
CERAD_MAP <- c("absent" = 0L, "sparse" = 1L, "moderate" = 2L, "frequent" = 3L,
               # numeric spellings, in case a re-export writes codes
               "0" = 0L, "1" = 1L, "2" = 2L, "3" = 3L)
BRAAK_MAP <- c("braak 0" = 0L, "0" = 0L, "none" = 0L,
               "braak i" = 1L, "i" = 1L, "braak ii" = 2L, "ii" = 2L,
               "braak iii" = 3L, "iii" = 3L, "braak iv" = 4L, "iv" = 4L,
               "braak v" = 5L, "v" = 5L, "braak vi" = 6L, "vi" = 6L)
THAL_MAP <- c("thal 0" = 0L, "0" = 0L, "thal 1" = 1L, "1" = 1L,
              "thal 2" = 2L, "2" = 2L, "thal 3" = 3L, "3" = 3L,
              "thal 4" = 4L, "4" = 4L, "thal 5" = 5L, "5" = 5L)
# CAA: evenly spaced 4-level ordinal.
CAA_MAP <- c("not identified" = 0L, "none" = 0L, "absent" = 0L,
             "mild" = 1L, "moderate" = 2L, "severe" = 3L,
             "0" = 0L, "1" = 1L, "2" = 2L, "3" = 3L)
COGSTATUS_MAP <- c("no dementia" = 0L, "reference" = 0L,
                   "mild cognitive impairment" = 1L, "mci" = 1L,
                   "dementia" = 2L)

# == Minimum sample sizes =====================================================
MIN_DONORS      <- 5L    # donors with valid clinical + group data
MIN_CELLS_GRP   <- 10L   # cells in an arm group, summed over donors
MIN_MV_DONORS   <- 20L   # donors needed to fit a multivariate model at all
JT_MIN_DONORS   <- 20L   # donors needed to cut tertiles and run the trend test
# Donors with fewer module cells than this are dropped from the JT ordering.
JT_MIN_CELLS_PER_DONOR <- 5L

# == Results RDS ==============================================================
# One RDS per arm: every table, the per-dataset objects, the declarations, the
# guardrail ledger and the ggplot objects.
SAVE_RESULTS_RDS    <- TRUE
SAVE_FIGURE_OBJECTS <- TRUE


# =============================================================================
# SECTION 0: SHARED HELPERS (arm-independent)
# =============================================================================

# Apply a declared map and report unmatched levels. Values are forced to
# character first (as.numeric on a factor returns the level index). Two passes:
# an exact match after folding Unicode whitespace, then a loose match on
# alphanumerics only for levels the exact pass missed (ambiguous loose keys are
# refused).
.norm_key <- function(z) {
    z <- as.character(z)
    # Strip non-breaking space, zero-width space and BOM bytes first, on the
    # untouched string (enc2utf8 would rewrite them).
    z <- gsub("\xc2\xa0|\xe2\x80\x8b|\xef\xbb\xbf", " ", z, useBytes = TRUE)
    # A bare 0xA0 (latin1 non-breaking space) only in non-UTF-8 strings.
    .lat <- !is.na(z) & !validUTF8(z)
    if (any(.lat))
        z[.lat] <- gsub("\xa0", " ", z[.lat], useBytes = TRUE)
    # Mark remaining valid UTF-8 strings so the Unicode pass can run.
    ok <- !is.na(z) & validUTF8(z)
    if (any(ok)) Encoding(z[ok]) <- "UTF-8"
    z <- tryCatch(gsub("[\\p{Zs}\\p{Cf}]+", " ", z, perl = TRUE),
                  error = function(e) z)
    z <- gsub("[[:space:]]+", " ", z)
    trimws(tolower(z))
}

.loose_key <- function(z) gsub("[^a-z0-9]", "", .norm_key(z))

.map_ordinal <- function(x, map, var_label) {
    raw  <- .norm_key(x)                 # what the exact pass compares
    key  <- raw
    out  <- unname(map[key])
    out  <- suppressWarnings(as.numeric(out))      # never trust a named vector

    # Loose pass, for the levels the exact pass could not place.
    .lm_key <- .loose_key(names(map))
    .lm_val <- as.integer(map)
    .amb    <- unique(.lm_key[duplicated(.lm_key) &
                              !duplicated(paste(.lm_key, .lm_val))])
    .keepl  <- !(.lm_key %in% .amb) & !duplicated(.lm_key)
    loose_map <- setNames(.lm_val[.keepl], .lm_key[.keepl])
    n_loose <- 0L
    if (any(is.na(out) & nzchar(raw))) {
        i <- is.na(out) & nzchar(raw)
        v <- suppressWarnings(as.numeric(unname(loose_map[.loose_key(x)[i]])))
        n_loose <- sum(!is.na(v))
        out[i] <- v
    }

    # `raw` (normalised) is the level the audit reports.
    raw <- .norm_key(x)
    # What a digit parser would have produced.
    naive <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", raw)))
    lv <- sort(unique(raw[!is.na(raw) & nzchar(raw)]))
    # A column that exists but is empty returns one explicit audit row.
    if (length(lv) == 0L)
        return(list(value = rep(NA_real_, length(raw)),
                    audit = data.frame(
                        variable = var_label, level = "(no non-empty value)",
                        n_donors = 0L, mapped_to = "COLUMN PRESENT BUT EMPTY",
                        naive_digit_parse = NA_character_,
                        raw_escaped = NA_character_,
                        stringsAsFactors = FALSE, row.names = NULL),
                    n_present = 0L, n_unmapped = 0L, n_naive_fail = 0L,
                    n_loose = 0L, n_ambiguous_loose = 0L))
    tab <- data.frame(
        variable   = var_label,
        level      = lv,
        n_donors   = vapply(lv, function(l) sum(raw == l, na.rm = TRUE),
                            integer(1)),
        mapped_to  = vapply(lv, function(l) {
            v <- out[match(l, raw)]
            if (is.na(v)) "UNMAPPED" else as.character(v)
        }, character(1)),
        naive_digit_parse = vapply(lv, function(l) {
            v <- naive[match(l, raw)]
            if (is.na(v)) "NA" else as.character(v)
        }, character(1)),
        # Escaped bytes of the level as it arrived, before normalisation.
        raw_escaped = vapply(lv, function(l) {
            z <- as.character(x)[match(l, raw)]
            if (is.na(z)) "NA" else encodeString(z)
        }, character(1)),
        stringsAsFactors = FALSE, row.names = NULL)
    n_present  <- sum(!is.na(raw) & nzchar(raw))
    n_unmapped <- sum(!is.na(raw) & nzchar(raw) & is.na(out))
    n_naive_na <- sum(!is.na(raw) & nzchar(raw) & is.na(naive))
    list(value = out, audit = tab, n_present = n_present,
         n_unmapped = n_unmapped, n_naive_fail = n_naive_na,
         n_loose = n_loose, n_ambiguous_loose = length(.amb))
}

# APOE genotype -> e4 dose ("3/4", "E3/E4"). Unparseable genotypes return NA.
apoe_to_e4_dose <- function(apoe_str) {
    s <- toupper(trimws(as.character(apoe_str)))
    s <- gsub("(^|[/|-])\\s*E", "\\1", s)     # "E3/E4" -> "3/4"
    s <- gsub("\\s", "", s)
    out <- rep(NA_integer_, length(s))
    ok  <- grepl("^[2-4][/|-][2-4]$", s)
    if (any(ok)) {
        alleles <- strsplit(s[ok], "[/|-]")
        out[ok] <- vapply(alleles, function(a) sum(a == "4"), integer(1))
    }
    out
}

.num_col <- function(df, nm) {
    if (!(nm %in% names(df))) return(rep(NA_real_, nrow(df)))
    suppressWarnings(as.numeric(trimws(as.character(df[[nm]]))))
}

# Case- and punctuation-insensitive column finder.
.find_col <- function(df, patterns) {
    nm <- names(df)
    key <- tolower(gsub("[^a-z0-9]", "", tolower(nm)))
    for (p in patterns) {
        pk <- tolower(gsub("[^a-z0-9]", "", tolower(p)))
        hit <- which(key == pk)
        if (length(hit) > 0) return(nm[hit[1]])
    }
    # Fallback: one whole token of the column name must equal the pattern.
    toks <- lapply(nm, function(z)
        tolower(Filter(nzchar, strsplit(z, "[^A-Za-z0-9]+")[[1]])))
    for (p in patterns) {
        pk <- tolower(gsub("[^A-Za-z0-9]", "", p))
        if (!nzchar(pk)) next
        hit <- which(vapply(toks, function(tk) pk %in% tk, logical(1)))
        if (length(hit) > 0) return(nm[hit[1]])
    }
    NA_character_
}

# Steiger test for two dependent, overlapping correlations (Meng, Rosenthal
# and Rubin 1992). Inputs are Spearman rho, so the Fisher-z variance is
# inflated by 1.06 (Fieller, Hartley and Pearson 1957).
SPEARMAN_Z_VAR_FACTOR <- 1.06

steiger_z <- function(r_1, r_2, r_12, n, spearman = TRUE) {
    out <- list(Z = NA_real_, p = NA_real_, delta_rho = NA_real_, n = n,
                note = "")
    if (any(!is.finite(c(r_1, r_2)))) {
        out$note <- "a correlation with the clinical variable is missing"
        return(out)
    }
    # The gap is reported even when the test cannot run.
    out$delta_rho <- r_1 - r_2
    if (length(n) != 1L || !is.finite(n) || n < 10L) {
        out$note <- sprintf("n = %s < 10",
                            if (length(n) == 1L) as.character(n) else "NA")
        return(out)
    }
    if (!is.finite(r_12)) {
        out$note <- "the correlation between the two predictors is missing"
        return(out)
    }
    # |r_12| at 1: the difference is not identifiable.
    if (abs(r_12) >= 1 - 1e-6) {
        out$note <- paste0("the two predictors are perfectly correlated ",
                           "(r = ", sprintf("%.6f", r_12), "); the difference ",
                           "in their correlations is not identifiable")
        return(out)
    }
    # Clamp just inside +-1 so atanh stays finite; noted on the row.
    cl <- function(r) max(min(r, 1 - 1e-8), -1 + 1e-8)
    if (max(abs(c(r_1, r_2))) >= 1)
        out$note <- "a correlation reached +-1 and was clamped for the transform"
    z1 <- atanh(cl(r_1)); z2 <- atanh(cl(r_2))
    r_bar  <- (cl(r_1) + cl(r_2)) / 2
    f_star <- (1 - cl(r_12)) / (2 * (1 - r_bar^2))
    f_star <- min(f_star, 1)
    h <- (1 - f_star * r_bar^2) / (1 - r_bar^2)
    v <- if (isTRUE(spearman)) SPEARMAN_Z_VAR_FACTOR else 1
    denom <- sqrt(2 * (1 - cl(r_12)) * h * v / (n - 3))
    if (!is.finite(denom) || denom < 1e-12) {
        out$note <- paste(c(out$note,
                            "the two predictors are collinear; Z is undefined"),
                          collapse = "; ")
        return(out)
    }
    out$Z <- (z1 - z2) / denom
    out$p <- 2 * stats::pnorm(-abs(out$Z))
    out
}

# Jonckheere-Terpstra trend test, donor level, with the Lehmann (1975) tie
# correction. Donors are cut into tertiles of their module Post-IP fraction.
jt_test <- function(y, g, k) {
    out <- list(T = NA_real_, E = NA_real_, V = NA_real_, Z = NA_real_,
                p = NA_real_, n = 0L, nk = integer(k), n_tie_groups = NA_integer_,
                note = "")
    ok <- is.finite(y) & !is.na(g)
    y <- y[ok]; g <- as.integer(g[ok])
    N <- length(y)
    if (N == 0L) { out$note <- "no donor has both a group and a value"; return(out) }
    # Group labels must be a contiguous 1..k for tabulate().
    nk <- tabulate(g, nbins = k)
    out$n <- N; out$nk <- nk
    if (sum(nk > 0L) < 2L) {
        out$note <- "fewer than two non-empty ordered groups"
        return(out)
    }
    TT <- 0
    for (i in seq_len(k - 1L)) {
        yi <- y[g == i]
        if (length(yi) == 0L) next
        for (j in seq(i + 1L, k)) {
            yj <- y[g == j]
            if (length(yj) == 0L) next
            d <- outer(yj, yi, "-")
            TT <- TT + sum(d > 0) + 0.5 * sum(d == 0)
        }
    }
    E_T <- (N^2 - sum(nk^2)) / 4
    tj <- as.integer(table(y))          # tie-group sizes in the POOLED sample
    out$n_tie_groups <- length(tj)
    if (N < 3L) { out$note <- "N < 3; the variance is undefined"; return(out) }
    A <- N * (N - 1) * (2 * N + 5) - sum(nk * (nk - 1) * (2 * nk + 5)) -
         sum(tj * (tj - 1) * (2 * tj + 5))
    B <- sum(nk * (nk - 1) * (nk - 2)) * sum(tj * (tj - 1) * (tj - 2))
    C <- sum(nk * (nk - 1)) * sum(tj * (tj - 1))
    V_T <- A / 72 +
           B / (36 * N * (N - 1) * (N - 2)) +
           C / (8 * N * (N - 1))
    out$T <- TT; out$E <- E_T; out$V <- V_T
    if (!is.finite(V_T) || V_T <= 0) {
        out$note <- "variance is zero or undefined; every value is tied"
        return(out)
    }
    out$Z <- (TT - E_T) / sqrt(V_T)
    out$p <- 2 * stats::pnorm(-abs(out$Z))
    out
}

# Tertile labels sized to the groups that exist.
.jt_x_labels <- function(k) {
    if (!is.finite(k) || k < 1L) return(character(0))
    if (k == 1L) return("T1")
    out <- paste0("T", seq_len(k))
    out[1]  <- paste0(out[1],  "\nPre-IP")
    out[k]  <- paste0(out[k],  "\nPost-IP")
    if (k > 2L) out[2:(k - 1L)] <- paste0(out[2:(k - 1L)], "\nmid")
    out
}

theme_h10 <- function(base_size = 10) {
    theme_bw(base_size = base_size) +
        theme(panel.grid.minor = element_blank(),
              plot.title    = element_text(face = "bold", size = base_size + 2),
              plot.subtitle = element_text(size = base_size - 1,
                                           colour = "gray30"),
              plot.caption  = element_text(size = base_size - 2.5,
                                           colour = "gray35", hjust = 0),
              strip.background = element_rect(fill = "gray95"),
              legend.position  = "right")
}

# Spearman with an explicit n; the tie warning is suppressed and the tie state
# recorded.
.spearman <- function(x, y, min_n) {
    keep <- is.finite(x) & is.finite(y)
    n <- sum(keep)
    if (n < min_n)
        return(list(rho = NA_real_, p = NA_real_, n = n, ties = NA,
                    note = sprintf("n = %d < %d", n, min_n)))
    xs <- x[keep]; ys <- y[keep]
    if (sd(xs) < 1e-12 || sd(ys) < 1e-12)
        return(list(rho = NA_real_, p = NA_real_, n = n, ties = NA,
                    note = "zero variance in one vector"))
    ct <- suppressWarnings(cor.test(xs, ys, method = "spearman", exact = FALSE))
    list(rho = unname(ct$estimate), p = ct$p.value, n = n,
         ties = anyDuplicated(xs) > 0L || anyDuplicated(ys) > 0L, note = "")
}


# =============================================================================
# SECTION 1: THE ARM RUNNER
# =============================================================================

# Worker called by cep_run_arms() with four arguments; `file_tag` is not read
# (the prefix comes from cep_file_tag_for(arm_tag)).
run_one_arm <- function(arm_tag, config_file, out_dir, file_tag) {

    # == Arm configuration ====================================================
    config_path <- file.path(CEP_HOME, config_file)
    if (!file.exists(config_path))
        stop("Configuration not found: ", config_path, call. = FALSE)

    assign("CEP_ARM", arm_tag, envir = globalenv())
    on.exit(suppressWarnings(rm("CEP_ARM", envir = globalenv())), add = TRUE)

    # The config prints as it loads, and the sinks are not open yet, so its
    # output is captured and replayed into the log once they are.
    .p10h_preamble <- tryCatch(
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

    .p10h_needed <- c("ARM_TAG", "R_EXPORT_ROOT", "GOI", "NORM_TAG",
                      "PRIMARY_STYPE", "NEOC_DIR", "NONNC_DIR", "SHARED_DIR",
                      "STYPE_COL", "DONOR_COL")
    .p10h_missing <- .p10h_needed[!vapply(.p10h_needed, exists, logical(1))]
    if (length(.p10h_missing) > 0)
        stop("The configuration loaded from\n    ", config_path,
             "\n  does not define: ", paste(.p10h_missing, collapse = ", "),
             "\n  This is the signature of a pre-migration config.",
             call. = FALSE)
    if (!identical(ARM_TAG, arm_tag))
        stop("The configuration reports ARM_TAG '", ARM_TAG,
             "' but this arm is '", arm_tag, "'.", call. = FALSE)
    .p10h_preamble <- c(.p10h_preamble, utils::capture.output(
        cat("  Config OK: ", ARM_TAG, " | ", GOI, " | ", NORM_TAG, "\n",
            sep = "")))

    DATASETS <- c("NeoC", "NonNeoC")
    set.seed(42)

    OUT_DIR <- out_dir
    dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
    FILE_PREFIX <- cep_out_prefix(PART_TAG, NORM_TAG, cep_file_tag_for(arm_tag))

    # == Logging ==============================================================
    LOG_PATH_10H <- file.path(OUT_DIR, paste0(FILE_PREFIX, "_console_log_",
                                              format(Sys.time(), "%Y%m%d_%H%M%S"),
                                              ".txt"))
    MSG_PATH_10H <- file.path(OUT_DIR, paste0(FILE_PREFIX,
                                              "_console_log_messages.tmp"))
    while (sink.number() > 0) sink()
    if (sink.number(type = "message") != 2L) sink(type = "message")
    # If the second connection fails, the first is closed.
    .p10h_log_con <- file(LOG_PATH_10H, open = "wt")
    .p10h_msg_con <- tryCatch(file(MSG_PATH_10H, open = "wt"),
                              error = function(e) {
                                  try(close(.p10h_log_con), silent = TRUE)
                                  stop("Could not open the message log at ",
                                       MSG_PATH_10H, ": ", conditionMessage(e),
                                       call. = FALSE)
                              })
    sink(.p10h_log_con, append = FALSE, split = TRUE)
    sink(.p10h_msg_con, type = "message")
    # Replay the config load into the log.
    if (length(.p10h_preamble) > 0)
        cat(paste(.p10h_preamble, collapse = "\n"), "\n\n", sep = "")
    .p10h_completed <- FALSE
    on.exit({
        tryCatch({
            if (!isTRUE(.p10h_completed))
                cat("\n", strrep("!", 72),
                    "\n  Arm did not complete; the error is printed to the console.\n",
                    strrep("!", 72), "\n", sep = "")
            if (sink.number(type = "message") != 2L) sink(type = "message")
            while (sink.number() > 0) sink()
            close(.p10h_msg_con); close(.p10h_log_con)
            if (file.exists(MSG_PATH_10H)) {
                .msg <- readLines(MSG_PATH_10H, warn = FALSE)
                if (length(.msg) > 0)
                    cat(c("", strrep("=", 72),
                          "  MESSAGE STREAM (warnings, messages, errors)",
                          strrep("=", 72), .msg, ""),
                        sep = "\n", file = LOG_PATH_10H, append = TRUE)
                unlink(MSG_PATH_10H)
            }
        }, error = function(e) invisible(NULL))
    }, add = TRUE)

    .p10h_warn_opt <- getOption("warn")
    .p10h_pm_opt   <- getOption("warnPartialMatchDollar")
    options(warn = 1)
    options(warnPartialMatchDollar = TRUE)
    on.exit({
        options(warn = .p10h_warn_opt)
        options(warnPartialMatchDollar = .p10h_pm_opt)
    }, add = TRUE)

    INCOMPLETE_MARKER <- file.path(OUT_DIR, "_RUN_INCOMPLETE.txt")
    writeLines(c("This arm did not finish.", "",
                 paste("Arm     :", arm_tag),
                 paste("Part    :", PART_DIR_TAG),
                 paste("Started :", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), "",
                 "Any files in this folder are from a partial run.",
                 "This file is removed automatically when the arm completes."),
               INCOMPLETE_MARKER)

    # == Guardrail ledger =====================================================
    GR_LOG_10H <- file.path(OUT_DIR, paste0(FILE_PREFIX, "_Guardrail_Log.txt"))
    .gr_log <- character(0); .gr_pass <- 0L; .gr_warn <- 0L; .gr_fail <- 0L
    gr10h <- function(check_id, status, detail) {
        # Zero-length or NA detail would print nothing; replaced.
        if (length(detail) == 0L || (length(detail) == 1L && is.na(detail)))
            detail <- "[detail unavailable: the value passed was NA or zero-length]"
        if (length(detail) > 1L) detail <- paste(detail, collapse = " ")
        line <- sprintf("  [%-4s]  %-22s %s", status, check_id, detail)
        cat(line, "\n", sep = "")
        .gr_log  <<- c(.gr_log, line)
        if (status == "PASS") .gr_pass <<- .gr_pass + 1L
        if (status == "FAIL") .gr_fail <<- .gr_fail + 1L
        if (status == "WARN") .gr_warn <<- .gr_warn + 1L
    }
    gr10h_write <- function() {
        writeLines(c(sprintf("Part10h Prognosis guardrail log  [%s]", arm_tag),
                     sprintf("Written: %s", Sys.time()), strrep("=", 72),
                     .gr_log, strrep("=", 72),
                     sprintf("PASS %d | WARN %d | FAIL %d",
                             .gr_pass, .gr_warn, .gr_fail)),
                   GR_LOG_10H)
    }

    # == Paths ================================================================
    V1_DIR    <- cep_arm_dir_for("Part10e_Validation1", arm_tag)
    V1_PREFIX <- paste0(cep_out_prefix("Part10e", NORM_TAG,
                                       cep_file_tag_for(arm_tag)), "_V1")
    p10h_v1_xlsx <- function(ds)
        file.path(V1_DIR, paste0(V1_PREFIX, "_", ds,
                                 "_SK_interaction_analysis.xlsx"))
    p10h_v1_rds <- function(ds)
        file.path(V1_DIR, paste0(V1_PREFIX, "_", ds, "_results.rds"))
    ds_dir <- function(ds) if (identical(ds, "NeoC")) NEOC_DIR else NONNC_DIR

    # == IP_y, the classifying quantity =======================================
    # Candidate IPRS RDS files, searched in this order.
    p10h_iprs_rds <- function() {
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
        if (length(hit) == 0L)
            return(list(path = NA_character_, obj = NULL,
                        note = paste(basename(cand), collapse = " / ")))
        list(path = hit[1],
             obj = tryCatch(readRDS(hit[1]), error = function(e) NULL),
             note = "")
    }
    .p10h_iprs <- p10h_iprs_rds()
    p10h_IP_y <- function(ds) {
        o <- .p10h_iprs$obj
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
    p10h_IP_x <- function(ds) {
        o <- .p10h_iprs$obj
        if (is.null(o)) return(NA_real_)
        v <- o[["inflection_points_DFG"]]
        if (!is.null(v) && ds %in% names(v))
            return(suppressWarnings(as.numeric(v[[ds]])))
        NA_real_
    }

    cat(sprintf("Part10h Prognosis - log started: %s\n\n", Sys.time()))
    cat("\n", strrep("=", 72), "\n")
    cat("  Part10h - donor-level clinical correlation for the CMRF arm groups\n")
    cat("  Arm          :", ARM_TAG, "\n")
    cat("  GOI          :", GOI, "|", NORM_TAG, "\n")
    cat("  Supertype    :", PRIMARY_STYPE, "\n")
    cat("  IP arms      : inherited (group_refined -> IP_status -> cDFG < IP_y);",
        "never split on GOI\n")
    cat("  Covariates   : 4-covariate (Age, APOE, Sex, PMI) AND 7-covariate\n")
    cat("                 (+ RIN, Brain pH, Fresh brain weight), side by side\n")
    cat("  V1 input     :", V1_DIR, "\n")
    cat("  Output dir   :", OUT_DIR, "\n")
    cat(strrep("=", 72), "\n\n")
    cat("  NOTE: the 4- and 7-covariate models are bounds; the ambient index\n")
    cat("        and batch are not adjusted for (both nest in donor).\n\n")
    # =========================================================================
    # STEP 1: Donor clinical table
    # =========================================================================
    cat("\nSTEP 1 -- donor clinical table\n")

    # S5b_donor_clinical.csv first; the legacy donor_clinical.csv is the
    # fallback (tissue-quality columns then come from obs_all_cells.csv).
    clin_cand <- c(file.path(SHARED_DIR, "S5b_donor_clinical.csv"),
                   file.path(SHARED_DIR, "donor_clinical.csv"),
                   file.path(R_EXPORT_ROOT, "shared", "S5b_donor_clinical.csv"),
                   file.path(R_EXPORT_ROOT, "shared", "donor_clinical.csv"))
    clin_hit <- clin_cand[file.exists(clin_cand)]
    if (length(clin_hit) == 0L)
        stop("No donor clinical table found. Looked for:\n  ",
             paste(clin_cand, collapse = "\n  "), call. = FALSE)
    donor_clinical <- read.csv(clin_hit[1], stringsAsFactors = FALSE,
                               check.names = FALSE)
    cat(sprintf("  Source : %s\n", basename(clin_hit[1])))
    cat(sprintf("  Donors : %d rows, %d columns\n",
                nrow(donor_clinical), ncol(donor_clinical)))

    dclin_donor_col <- .find_col(donor_clinical,
                                 c("Donor ID", "donor_id", "DonorID", "donor"))
    if (is.na(dclin_donor_col))
        stop("The donor clinical table has no donor identifier column.",
             call. = FALSE)
    donor_clinical$donor_id <- trimws(as.character(
        donor_clinical[[dclin_donor_col]]))
    if (anyDuplicated(donor_clinical$donor_id) > 0L)
        stop("The donor clinical table has duplicated donor IDs; one row per ",
             "donor is required.", call. = FALSE)
    # A table with no usable donor ID stops here.
    n_usable_donor <- sum(!is.na(donor_clinical$donor_id) &
                          nzchar(donor_clinical$donor_id))
    if (nrow(donor_clinical) == 0L || n_usable_donor == 0L)
        stop("The donor clinical table has ", nrow(donor_clinical),
             " row(s) and ", n_usable_donor, " usable donor ID(s).\n  ",
             basename(clin_hit[1]), " is present but carries no donor this ",
             "script can join on.", call. = FALSE)

    # =========================================================================
    # STEP 1a: Ordinal parsing
    # =========================================================================
    cat("\nSTEP 1a -- ordinal parsing (declared maps, not digit extraction)\n")

    ord_specs <- list(
        list(out = "ADNC_rank",    map = ADNC_MAP,      label = "ADNC",
             pats = c("Overall AD neuropathological Change", "ADNC",
                      "AD neuropathological")),
        list(out = "Braak_num",    map = BRAAK_MAP,     label = "Braak",
             pats = c("Braak", "Braak stage")),
        list(out = "Thal_num",     map = THAL_MAP,      label = "Thal",
             pats = c("Thal", "Thal phase")),
        list(out = "CERAD_num",    map = CERAD_MAP,     label = "CERAD",
             pats = c("CERAD score", "CERAD")),
        list(out = "CAA_num",      map = CAA_MAP,       label = "CAA",
             pats = c("Overall CAA Score", "CAA")),
        list(out = "CogStatus_num", map = COGSTATUS_MAP, label = "Cognitive status",
             pats = c("Cognitive Status", "CogStatus", "Cognitive"))
    )

    ordinal_audit <- list()
    for (sp in ord_specs) {
        src <- .find_col(donor_clinical, sp$pats)
        if (is.na(src)) {
            donor_clinical[[sp$out]] <- rep(NA_real_, nrow(donor_clinical))
            cat(sprintf("  %-16s source column not found (looked for: %s)\n",
                        sp$label, paste(sp$pats, collapse = ", ")))
            ordinal_audit[[sp$label]] <- data.frame(
                variable = sp$label, level = NA_character_, n_donors = NA_integer_,
                mapped_to = "SOURCE COLUMN NOT FOUND",
                naive_digit_parse = NA_character_, raw_escaped = NA_character_,
                stringsAsFactors = FALSE)
            next
        }
        m <- .map_ordinal(donor_clinical[[src]], sp$map, sp$label)
        donor_clinical[[sp$out]] <- m$value
        ordinal_audit[[sp$label]] <- m$audit
        cat(sprintf(paste0("  %-16s <- \"%s\"  |  %d level(s), %d/%d mapped%s",
                           "  |  a digit parser would fail %d/%d\n"),
                    sp$label, src, nrow(m$audit),
                    m$n_present - m$n_unmapped, m$n_present,
                    if (m$n_loose > 0)
                        sprintf(" (%d rescued by the loose pass)",
                                m$n_loose) else "",
                    m$n_naive_fail, m$n_present))
        if (m$n_unmapped > 0) {
            .u <- m$audit$mapped_to == "UNMAPPED"
            cat(sprintf("    Unmapped levels: %s\n",
                        paste(m$audit$level[.u], collapse = " | ")))
            # Escaped bytes of the unmapped levels.
            cat(sprintf("    As received    : %s\n",
                        paste(m$audit$raw_escaped[.u], collapse = " | ")))
        }
    }
    ordinal_audit_df <- do.call(rbind, ordinal_audit)
    rownames(ordinal_audit_df) <- NULL

    # == Continuous / derived clinical columns ================================
    apoe_src <- .find_col(donor_clinical, c("APOE Genotype", "APOE", "APOE4"))
    donor_clinical$APOE_e4_dose <- if (is.na(apoe_src)) rep(NA_real_,
                                                            nrow(donor_clinical))
                                   else apoe_to_e4_dose(donor_clinical[[apoe_src]])
    mmse_src <- .find_col(donor_clinical, c("Last MMSE Score", "MMSE"))
    donor_clinical$MMSE_num <- if (is.na(mmse_src)) rep(NA_real_,
                                                        nrow(donor_clinical))
                               else .num_col(donor_clinical, mmse_src)
    age_src <- .find_col(donor_clinical, c("Age at Death", "Age_at_Death", "Age"))
    # Age "90+" is read as 90; the count is printed.
    age_raw <- if (is.na(age_src)) rep(NA_character_, nrow(donor_clinical))
               else trimws(as.character(donor_clinical[[age_src]]))
    n_censored <- sum(grepl("\\+", age_raw), na.rm = TRUE)
    donor_clinical$Age_num <- suppressWarnings(
        as.numeric(gsub("[^0-9.]", "", age_raw)))
    donor_clinical$Age_censored <- grepl("\\+", age_raw)
    cat(sprintf("  %-16s <- \"%s\"  |  %d donor(s) censored as \"90+\" and read at the floor\n",
                "Age at death", if (is.na(age_src)) "not found" else age_src,
                n_censored))

    # =========================================================================
    # STEP 1b: Tissue-quality and demographic covariates
    # =========================================================================
    cat("\nSTEP 1b -- tissue-quality and demographic covariates\n")

    cov_specs <- list(
        list(out = "Sex_raw",         pats = c("Sex", "sex"),          num = FALSE),
        list(out = "PMI_num",         pats = c("PMI", "Post-Mortem Interval",
                                               "PostMortem Interval"), num = TRUE),
        list(out = "RIN_num",         pats = c("RIN", "RNA Integrity Number",
                                               "RNA Integrity"),       num = TRUE),
        list(out = "BrainPH_num",     pats = c("Brain pH", "Brain_pH", "pH"),
             num = TRUE),
        list(out = "BrainWeight_num", pats = c("Fresh Brain Weight",
                                               "Brain Weight"),        num = TRUE),
        list(out = "Lewy_raw",        pats = c("Highest Lewy Body Disease",
                                               "Lewy Body Disease", "Lewy"),
             num = FALSE),
        list(out = "Severe_raw",      pats = c("Severely Affected Donor",
                                               "Severely Affected"),   num = FALSE)
    )

    # Covariates missing from the clinical table are taken from
    # obs_all_cells.csv, one value per donor (conflicts are reported).
    obs_all_path <- c(file.path(SHARED_DIR, "obs_all_cells.csv"),
                      file.path(R_EXPORT_ROOT, "shared", "obs_all_cells.csv"))
    obs_all_path <- obs_all_path[file.exists(obs_all_path)]
    obs_extra <- NULL
    cov_conflicts <- character(0)
    if (length(obs_all_path) > 0) {
        # Read the header first, then only the columns needed.
        .obs_head <- tryCatch(
            read.csv(obs_all_path[1], nrows = 1L, stringsAsFactors = FALSE,
                     check.names = FALSE),
            error = function(e) NULL)
        .keep <- NULL
        if (!is.null(.obs_head)) {
            .d <- .find_col(.obs_head, c("Donor ID", "donor_id", "DonorID"))
            if (!is.na(.d)) {
                .keep <- .d
                for (cs in cov_specs) {
                    h <- .find_col(.obs_head, cs$pats)
                    if (!is.na(h)) .keep <- c(.keep, h)
                }
                .keep <- unique(.keep)
            }
        }
        obs_all <- NULL
        if (!is.null(.keep)) {
            obs_all <- tryCatch({
                if (requireNamespace("data.table", quietly = TRUE))
                    as.data.frame(
                        data.table::fread(obs_all_path[1], select = .keep,
                                          showProgress = FALSE),
                        stringsAsFactors = FALSE)
                else {
                    # Without data.table, skip the other columns.
                    cc <- ifelse(names(.obs_head) %in% .keep, NA, "NULL")
                    read.csv(obs_all_path[1], colClasses = cc,
                             stringsAsFactors = FALSE, check.names = FALSE)
                }
            }, error = function(e) {
                cat(sprintf("  obs_all_cells.csv could not be read: %s\n",
                            conditionMessage(e)))
                NULL
            })
        }
        rm(.obs_head)
        if (!is.null(obs_all)) {
            oa_donor <- .find_col(obs_all, c("Donor ID", "donor_id", "DonorID"))
            if (!is.na(oa_donor)) {
                want <- character(0)
                for (cs in cov_specs) {
                    hit <- .find_col(obs_all, cs$pats)
                    if (!is.na(hit)) want <- c(want, setNames(hit, cs$out))
                }
                if (length(want) > 0) {
                    sub <- obs_all[, c(oa_donor, unname(want)), drop = FALSE]
                    names(sub) <- c("donor_id", names(want))
                    sub$donor_id <- trimws(as.character(sub$donor_id))
                    sub <- unique(sub)
                    dup <- sub$donor_id[duplicated(sub$donor_id)]
                    if (length(dup) > 0) {
                        cov_conflicts <- unique(dup)
                        cat(sprintf(paste0("  WARNING: %d donor(s) carry more ",
                                           "than one value for a donor-level ",
                                           "covariate; first row kept. Donors: %s\n"),
                                    length(dup),
                                    paste(utils::head(dup, 5), collapse = ", ")))
                        sub <- sub[!duplicated(sub$donor_id), , drop = FALSE]
                    }
                    obs_extra <- sub
                }
            }
        }
        rm(obs_all); gc(verbose = FALSE)
        cat(sprintf("  obs_all_cells.csv: %s\n",
                    if (is.null(obs_extra)) "present but no usable donor key"
                    else sprintf("%d donor(s), %d covariate column(s)",
                                 nrow(obs_extra), ncol(obs_extra) - 1L)))
    } else {
        cat("  obs_all_cells.csv not found; covariates come from the clinical",
            " table only\n", sep = "")
    }

    cov_source <- character(0)
    for (cs in cov_specs) {
        src <- .find_col(donor_clinical, cs$pats)
        if (!is.na(src)) {
            donor_clinical[[cs$out]] <- if (cs$num) .num_col(donor_clinical, src)
                                        else trimws(as.character(
                                            donor_clinical[[src]]))
            cov_source[cs$out] <- paste0("clinical table: ", src)
        } else if (!is.null(obs_extra) && cs$out %in% names(obs_extra)) {
            v <- obs_extra[[cs$out]][match(donor_clinical$donor_id,
                                           obs_extra$donor_id)]
            donor_clinical[[cs$out]] <- if (cs$num)
                suppressWarnings(as.numeric(trimws(as.character(v))))
                else trimws(as.character(v))
            cov_source[cs$out] <- "obs_all_cells.csv"
        } else {
            donor_clinical[[cs$out]] <- if (cs$num) rep(NA_real_,
                                                        nrow(donor_clinical))
                                        else rep(NA_character_,
                                                 nrow(donor_clinical))
            cov_source[cs$out] <- "NOT AVAILABLE"
        }
    }

    # Sex: Female 0, Male 1, anything else NA.
    sx <- tolower(trimws(as.character(donor_clinical$Sex_raw)))
    donor_clinical$Sex_num <- ifelse(sx %in% c("male", "m"), 1L,
                              ifelse(sx %in% c("female", "f"), 0L, NA_integer_))
    # Lewy: "Not Identified" 0, named pathology 1, blank NA.
    lw <- trimws(as.character(donor_clinical$Lewy_raw))
    donor_clinical$Lewy_binary <- ifelse(
        is.na(lw) | !nzchar(lw), NA_integer_,
        ifelse(grepl("not identified|none|absent", lw, ignore.case = TRUE),
               0L, 1L))
    # Severely Affected Donor flag.
    sv <- tolower(trimws(as.character(donor_clinical$Severe_raw)))
    donor_clinical$Severe_flag <- ifelse(
        is.na(sv) | !nzchar(sv), NA_integer_,
        ifelse(sv %in% c("y", "yes", "true", "1", "severely affected"), 1L, 0L))

    for (cs in cov_specs)
        cat(sprintf("  %-16s : %s\n", cs$out, cov_source[cs$out]))
    .rng <- function(v) {
        v <- v[is.finite(v)]
        if (length(v) == 0) return("no finite values")
        sprintf("n=%d, range [%.2f, %.2f], median %.2f",
                length(v), min(v), max(v), median(v))
    }
    cat(sprintf("    PMI    : %s\n", .rng(donor_clinical$PMI_num)))
    cat(sprintf("    RIN    : %s\n", .rng(donor_clinical$RIN_num)))
    cat(sprintf("    Brain pH: %s\n", .rng(donor_clinical$BrainPH_num)))
    cat(sprintf("    Weight : %s\n", .rng(donor_clinical$BrainWeight_num)))
    cat(sprintf("    Sex    : %d female, %d male, %d unknown\n",
                sum(donor_clinical$Sex_num == 0, na.rm = TRUE),
                sum(donor_clinical$Sex_num == 1, na.rm = TRUE),
                sum(is.na(donor_clinical$Sex_num))))
    cat(sprintf("    Severely Affected Donor: %d flagged, %d not, %d unknown\n",
                sum(donor_clinical$Severe_flag == 1, na.rm = TRUE),
                sum(donor_clinical$Severe_flag == 0, na.rm = TRUE),
                sum(is.na(donor_clinical$Severe_flag))))

    # == Covariate completeness ===============================================
    .complete_n <- function(set) {
        cols <- unlist(COV_MAP[set], use.names = FALSE)
        have <- cols[cols %in% names(donor_clinical)]
        if (length(have) < length(cols)) return(NA_integer_)
        sum(complete.cases(donor_clinical[, have, drop = FALSE]))
    }
    n_cov4 <- .complete_n(COV_SET_4)
    n_cov7 <- .complete_n(COV_SET_7)
    cat(sprintf(paste0("  Complete cases: %s on the 4-covariate set, %s on the",
                       " 7-covariate set\n"),
                if (is.na(n_cov4)) "NA (a covariate column is missing)"
                else as.character(n_cov4),
                if (is.na(n_cov7)) "NA (a covariate column is missing)"
                else as.character(n_cov7)))
    cat("  Quote the realised n per model.\n")

    # =========================================================================
    # STEP 2: Per-dataset analysis
    # =========================================================================
    all_res <- list()

    # Every figure passes through save_pp(), which stores it here by basename.
    .fig_store <- list()

    for (ds in DATASETS) {
        cat(sprintf("\n%s\n  %s\n%s\n", strrep("=", 70), ds, strrep("=", 70)))
        res <- tryCatch({

        # == 2a: Validation 1 refined cells ===================================
        xp <- p10h_v1_xlsx(ds)
        if (!file.exists(xp))
            stop("Validation1 workbook not found: ", basename(xp), call. = FALSE)
        rc <- as.data.frame(
            readxl::read_excel(xp, sheet = "Step4_Refined_Cells"),
            stringsAsFactors = FALSE)
        need <- c("barcode", "cell_class")
        miss <- setdiff(need, names(rc))
        if (length(miss) > 0)
            stop("Step4_Refined_Cells lacks: ", paste(miss, collapse = ", "),
                 call. = FALSE)
        rc$barcode    <- as.character(rc$barcode)
        rc$cell_class <- as.character(rc$cell_class)
        if (anyDuplicated(rc$barcode) > 0L)
            stop("Step4_Refined_Cells has duplicated barcodes.", call. = FALSE)
        # u_all is optional (printed contrast only); cDFG_expr classifies.
        rc$u_all <- if ("u_all" %in% names(rc))
            suppressWarnings(as.numeric(rc$u_all)) else NA_real_
        rc$cDFG_expr <- if ("cDFG_expr" %in% names(rc))
            suppressWarnings(as.numeric(rc$cDFG_expr)) else NA_real_
        n_class_na <- sum(is.na(rc$cell_class))
        if (n_class_na > 0)
            cat(sprintf("  %d row(s) with a missing cell_class dropped\n",
                        n_class_na))
        rc <- rc[!is.na(rc$cell_class), , drop = FALSE]
        cat(sprintf("  Step4_Refined_Cells: %d cell(s)\n", nrow(rc)))
        for (cls in CMRF_CLASS_LEVS)
            cat(sprintf("    %-25s : %d\n", cls,
                        sum(rc$cell_class == cls, na.rm = TRUE)))

        # == 2b: the IP, and the arms =========================================
        ipy  <- p10h_IP_y(ds)
        ip_y <- ipy$value
        ip_x <- p10h_IP_x(ds)
        # u_IP is reference only; a missing Validation 1 RDS is a note.
        rp <- p10h_v1_rds(ds)
        u_ip <- NA_real_
        if (!file.exists(rp)) {
            cat(sprintf("  NOTE: Validation1 RDS not found (%s); the GOI-axis\n",
                        basename(rp)),
                "        contrast line is skipped. u_IP is reference only.\n",
                sep = "")
        } else {
            v1 <- tryCatch(readRDS(rp), error = function(e) NULL)
            if (is.null(v1))
                cat("  NOTE: Validation1 RDS could not be read; the GOI-axis",
                    " contrast line is skipped.\n", sep = "")
            else u_ip <- suppressWarnings(as.numeric(v1$u_IP)[1])
            rm(v1); gc(verbose = FALSE)
        }
        cat(sprintf("  IP_y = %s  (%s)   |   IP_x = %s (reported IP, not the split)\n",
                    if (is.finite(ip_y)) sprintf("%.6f", ip_y) else "NA",
                    if (is.na(ipy$field)) "no IPRS RDS field" else ipy$field,
                    if (is.finite(ip_x)) sprintf("%.6f", ip_x) else "NA"))

        # Inherit the arm in order; never fall back to the GOI axis.
        pick_arm <- function(z) {
            z <- as.character(z)
            ifelse(is.na(z), NA_character_,
            ifelse(grepl("Post-IP", z, fixed = TRUE), "Post-IP",
            ifelse(grepl("Pre-IP",  z, fixed = TRUE), "Pre-IP", NA_character_)))
        }
        arm <- rep(NA_character_, nrow(rc))
        arm_src <- rep(NA_character_, nrow(rc))
        if ("group_refined" %in% names(rc)) {
            a <- pick_arm(rc$group_refined)
            arm_src[is.na(arm) & !is.na(a)] <- "group_refined"
            arm[is.na(arm) & !is.na(a)] <- a[is.na(arm) & !is.na(a)]
        }
        if ("IP_status" %in% names(rc)) {
            a <- pick_arm(rc$IP_status)
            arm_src[is.na(arm) & !is.na(a)] <- "IP_status"
            arm[is.na(arm) & !is.na(a)] <- a[is.na(arm) & !is.na(a)]
        }
        if (any(is.na(arm)) && is.finite(ip_y) && any(is.finite(rc$cDFG_expr))) {
            i <- is.na(arm) & is.finite(rc$cDFG_expr)
            arm[i] <- ifelse(rc$cDFG_expr[i] < ip_y, "Pre-IP", "Post-IP")
            arm_src[i] <- "cDFG_at_IP_y"
        }
        rc$arm <- arm
        n_unlab <- sum(is.na(arm) & rc$cell_class %in% MODULE_CLASSES)
        # Only unclassifiable module cells stop the run.
        if (n_unlab > 0 && !is.finite(ip_y))
            stop(sprintf(paste0("%d module cell(s) cannot be classified: ",
                                "no usable group_refined or IP_status and no ",
                                "IP_y."), n_unlab),
                 call. = FALSE)
        .src_tab <- function(sel) table(factor(arm_src[sel],
                                              levels = c("group_refined",
                                                         "IP_status",
                                                         "cDFG_at_IP_y")))
        asrc      <- .src_tab(rc$cell_class %in% MODULE_CLASSES)
        asrc_excl <- .src_tab(rc$cell_class %in% NONMODULE_CLASSES)
        # Module and non-module sources are counted separately; non-module
        # cells normally take the cDFG route.
        n_y_fall_excl <- as.integer(asrc_excl[["cDFG_at_IP_y"]])
        n_y_fall <- as.integer(asrc[["cDFG_at_IP_y"]])
        n_unlab_excl <- sum(is.na(arm) & rc$cell_class %in% NONMODULE_CLASSES)
        cat(sprintf(paste0("  Arm source, module cells  : group_refined %d, ",
                           "IP_status %d, cDFG@IP_y %d, unlabelled %d\n"),
                    as.integer(asrc[["group_refined"]]),
                    as.integer(asrc[["IP_status"]]), n_y_fall, n_unlab))
        cat(sprintf(paste0("  Arm source, NON-MODULE cells: group_refined %d, ",
                           "IP_status %d, cDFG@IP_y %d, unlabelled %d\n"),
                    as.integer(asrc_excl[["group_refined"]]),
                    as.integer(asrc_excl[["IP_status"]]), n_y_fall_excl,
                    n_unlab_excl))
        if (n_unlab_excl > 0)
            cat(sprintf(paste0("    %d excluded cell(s) have no cDFG_expr and ",
                               "carry no arm; they enter Excluded (all) only.\n"),
                        n_unlab_excl))

        # Agreement of inherited labels with the cDFG rule (fallback cells left out).
        arm_agree <- NA_real_; arm_n_cmp <- 0L; arm_agree_u <- NA_real_
        if (is.finite(ip_y) && any(is.finite(rc$cDFG_expr))) {
            i <- !is.na(arm) & arm_src != "cDFG_at_IP_y" &
                 is.finite(rc$cDFG_expr)
            if (any(i)) {
                rule <- ifelse(rc$cDFG_expr[i] < ip_y, "Pre-IP", "Post-IP")
                arm_agree <- mean(arm[i] == rule)
                arm_n_cmp <- sum(i)
            }
        }
        if (is.finite(ip_x) && any(is.finite(rc$u_all)) && is.finite(u_ip)) {
            i <- !is.na(arm) & arm_src != "cDFG_at_IP_y" & is.finite(rc$u_all)
            if (any(i)) {
                rule_u <- ifelse(rc$u_all[i] < u_ip, "Pre-IP", "Post-IP")
                arm_agree_u <- mean(arm[i] == rule_u)
            }
        }
        cat(sprintf(paste0("  Inherited vs `cDFG < IP_y`: %s on %d cell(s)  |",
                           "  a GOI split at u_IP would agree %s\n"),
                    if (is.na(arm_agree)) "not checked"
                    else sprintf("%.2f%%", 100 * arm_agree), arm_n_cmp,
                    if (is.na(arm_agree_u)) "NA"
                    else sprintf("%.2f%%", 100 * arm_agree_u)))

        # == 2c: donor denominators ===========================================
        obs_ds_path <- file.path(ds_dir(ds), "metadata", "obs_dataset.csv")
        if (!file.exists(obs_ds_path))
            stop("obs_dataset.csv not found: ", obs_ds_path, call. = FALSE)
        obs_dataset <- read.csv(obs_ds_path, stringsAsFactors = FALSE,
                                check.names = FALSE)
        # Column names come from the arm's config.
        d_col <- .find_col(obs_dataset, c(DONOR_COL, "Donor ID", "donor_id",
                                          "DonorID"))
        b_col <- .find_col(obs_dataset, c("barcode", "Barcode", "cell_id"))
        s_col <- .find_col(obs_dataset, c(STYPE_COL, "Supertype", "supertype",
                                          "Supertype (non-expanded)"))
        if (is.na(d_col) || is.na(b_col) || is.na(s_col))
            stop("obs_dataset.csv lacks a donor, barcode or supertype column.",
                 call. = FALSE)
        obs_dataset[[b_col]] <- as.character(obs_dataset[[b_col]])
        obs_dataset[[d_col]] <- trimws(as.character(obs_dataset[[d_col]]))

        # Cell counts are reported; only internal consistency is checked.
        n_obs_total <- nrow(obs_dataset)
        obs_pri <- obs_dataset[obs_dataset[[s_col]] == PRIMARY_STYPE, ,
                               drop = FALSE]
        n_obs_pri <- nrow(obs_pri)
        cat(sprintf("  obs_dataset: %d cell(s) total, %d of supertype %s\n",
                    n_obs_total, n_obs_pri, PRIMARY_STYPE))
        if (n_obs_pri == 0L)
            stop("obs_dataset.csv holds no cell of supertype ", PRIMARY_STYPE,
                 ".", call. = FALSE)

        # Cells with no donor ID are counted and dropped.
        .d_bad <- is.na(obs_pri[[d_col]]) | !nzchar(obs_pri[[d_col]])
        if (any(.d_bad)) {
            cat(sprintf(paste0("  %d of %d %s cell(s) carry no donor ID and ",
                               "are left out.\n"),
                        sum(.d_bad), nrow(obs_pri), PRIMARY_STYPE))
            obs_pri <- obs_pri[!.d_bad, , drop = FALSE]
        }
        n_no_donor_obs <- sum(.d_bad)
        if (nrow(obs_pri) == 0L)
            stop("Every ", PRIMARY_STYPE, " cell in obs_dataset.csv lacks a ",
                 "donor ID; there is no denominator to compute against.",
                 call. = FALSE)
        donor_totals <- as.data.frame(
            table(donor_id = obs_pri[[d_col]]), stringsAsFactors = FALSE)
        names(donor_totals) <- c("donor_id", "n_total")
        donor_totals$n_total <- as.integer(donor_totals$n_total)
        cat(sprintf("  Donors with %s cells: %d\n", PRIMARY_STYPE,
                    nrow(donor_totals)))

        bc_donor <- obs_dataset[, c(b_col, d_col)]
        names(bc_donor) <- c("barcode", "donor_id")
        bc_donor <- bc_donor[!duplicated(bc_donor$barcode), , drop = FALSE]

        # Donor-level extra variables (DONOR_EXTRA). A variable that does not
        # resolve or has no variance stays NA and disables only the sets that
        # use it.
        donor_extra <- data.frame(donor_id = donor_totals$donor_id,
                                  stringsAsFactors = FALSE)
        extra_note <- character(0)
        extra_ok   <- character(0)
        for (.e in DONOR_EXTRA) {
            .src <- NULL; .from <- NA_character_
            if (identical(.e$type, "ratio")) {
                .cn <- .find_col(obs_pri, .e$num); .cx <- .find_col(obs_pri, .e$den)
                if (!is.na(.cn) && !is.na(.cx)) {
                    .den <- suppressWarnings(as.numeric(obs_pri[[.cx]]))
                    # A zero denominator is Inf, not a very large ratio.
                    .den[!is.finite(.den) | .den == 0] <- NA_real_
                    .src <- suppressWarnings(as.numeric(obs_pri[[.cn]])) / .den
                    .from <- paste(.cn, "/", .cx)
                }
            } else {
                .fc <- .find_col(obs_pri, .e$patterns)
                if (!is.na(.fc)) { .src <- obs_pri[[.fc]]; .from <- .fc }
            }
            if (is.null(.src)) {
                donor_extra[[.e$out]] <- NA_real_
                extra_note <- c(extra_note, paste0(.e$out, ": column absent"))
                next
            }
        # Ordinals: NA tokens removed, then mapped on the normalised key.
            if (identical(.e$type, "ordinal")) {
                .k <- .norm_key(.src)
                .k[.k %in% NA_TOKENS] <- NA_character_
                .v <- unname(.e$map[.k])
                .n_un <- sum(is.na(.v) & !is.na(.k))
                if (.n_un > 0L) {
                    .bad <- utils::head(unique(.k[is.na(.v) & !is.na(.k)]), 5)
                    extra_note <- c(extra_note,
                        sprintf("%s: %d cell(s) UNMAPPED (%s)", .e$out, .n_un,
                                paste(.bad, collapse = "; ")))
                }
            } else {
                .v <- suppressWarnings(as.numeric(as.character(.src)))
            }
            .agg <- if (identical(.e$agg, "mean")) function(z)
                        mean(z, na.rm = TRUE) else function(z) {
                            z <- z[!is.na(z)]
                            if (!length(z)) NA_real_ else z[1]
                        }
            .keep <- !is.na(obs_pri[[d_col]]) & nzchar(obs_pri[[d_col]])
            .ag <- stats::aggregate(list(v = .v[.keep]),
                                    by = list(donor_id = obs_pri[[d_col]][.keep]),
                                    FUN = .agg)
            donor_extra[[.e$out]] <- .ag$v[match(donor_extra$donor_id,
                                                 .ag$donor_id)]
            # A donor-constant variable must be constant within donor.
            if (identical(.e$agg, "const")) {
                .nc <- sum(vapply(split(.v[.keep], obs_pri[[d_col]][.keep]),
                                  function(z) {
                                      z <- z[!is.na(z)]
                                      length(unique(z)) > 1L
                                  }, logical(1)))
                if (.nc > 0L)
                    extra_note <- c(extra_note,
                        sprintf("%s: NOT donor-constant on %d donor(s)",
                                .e$out, .nc))
            }
            .n_fin <- sum(is.finite(donor_extra[[.e$out]]))
            .sd <- stats::sd(donor_extra[[.e$out]], na.rm = TRUE)
            if (.n_fin >= MIN_MV_DONORS && isTRUE(.sd > 0)) {
                extra_ok <- c(extra_ok, .e$out)
            } else {
                extra_note <- c(extra_note,
                    sprintf("%s: unusable (%d donor(s), sd %s)", .e$out, .n_fin,
                            if (is.finite(.sd)) sprintf("%.3g", .sd) else "NA"))
            }
            extra_note <- c(extra_note, sprintf("%s <- %s [%d donor(s)]",
                                                .e$out, .from, .n_fin))
        }
        # mito_cov_ok is reported by a guardrail.
        mito_cov_ok <- all(c("MitoFrac_num", "nUMI_num") %in% extra_ok)
        depth_note  <- paste(extra_note, collapse = " | ")
        cat(sprintf("  Donor-level extra variables: %d of %d usable\n",
                    length(extra_ok), length(DONOR_EXTRA)))
        for (.z in extra_note) cat(sprintf("    %s\n", .z))

        # One donor-level clinical frame read by every consumer.
        donor_clin_all <- merge(donor_clinical, donor_extra, by = "donor_id",
                                all = TRUE)
        stopifnot(all(vapply(CLINICAL_PARAMS, function(z) z$col,
                             character(1)) %in% names(donor_clin_all)))

        rc$donor_id <- bc_donor$donor_id[match(rc$barcode, bc_donor$barcode)]
        # A blank donor ID is set to NA.
        rc$donor_id[!is.na(rc$donor_id) & !nzchar(rc$donor_id)] <- NA_character_
        rm(obs_dataset, obs_pri, bc_donor); gc(verbose = FALSE)
        n_no_donor <- sum(is.na(rc$donor_id))
        if (n_no_donor > 0)
            cat(sprintf(paste0("  WARNING: %d of %d refined cell(s) have no ",
                               "donor in obs_dataset.csv and are left out.\n"),
                        n_no_donor, nrow(rc)))

        # == 2d: per-donor proportions per arm group ==========================
        # Denominator: all PRIMARY_STYPE cells of the donor; a donor with none in
        # a group is a real zero.
        grp_pct <- list(); grp_meta <- list()
        for (gn in names(ARM_GROUPS)) {
            g <- ARM_GROUPS[[gn]]
            sel <- rc$cell_class %in% g$classes & !is.na(rc$donor_id)
            if (!is.null(g$arms)) sel <- sel & !is.na(rc$arm) & rc$arm %in% g$arms
            n_cells <- sum(sel)
            if (n_cells < MIN_CELLS_GRP) {
                cat(sprintf("    SKIP %-17s : %d cell(s) < %d\n", gn, n_cells,
                            MIN_CELLS_GRP))
                grp_meta[[gn]] <- list(n_cells = n_cells, skipped = TRUE)
                next
            }
            cnt <- as.data.frame(table(donor_id = rc$donor_id[sel]),
                                 stringsAsFactors = FALSE)
            names(cnt) <- c("donor_id", "n_grp")
            m <- merge(donor_totals, cnt, by = "donor_id", all.x = TRUE)
            m$n_grp <- ifelse(is.na(m$n_grp), 0L, as.integer(m$n_grp))
            m$pct <- 100 * m$n_grp / m$n_total
            grp_pct[[gn]] <- m[, c("donor_id", "n_grp", "n_total", "pct")]
            grp_meta[[gn]] <- list(n_cells = n_cells, skipped = FALSE,
                                   n_donors = nrow(m))
            cat(sprintf("    %-17s : %6d cell(s) over %d donor(s), median %.3f%%\n",
                        gn, n_cells, nrow(m), median(m$pct)))
        }
        if (length(grp_pct) == 0L)
            stop("No arm group cleared the ", MIN_CELLS_GRP,
                 "-cell floor for ", ds, ".", call. = FALSE)

        # == 2e: Spearman, BH within group ====================================
        cat("\n  Spearman rho, BH-FDR within each group across the ",
            length(CLINICAL_PARAMS), " clinical variables\n", sep = "")
        rows <- list()
        for (gn in names(grp_pct)) {
            m <- merge(grp_pct[[gn]], donor_clin_all, by = "donor_id",
                       all.x = TRUE)
            gr_rows <- list()
            for (cp in CLINICAL_PARAMS) {
                y <- if (cp$col %in% names(m)) suppressWarnings(
                    as.numeric(m[[cp$col]])) else rep(NA_real_, nrow(m))
                st <- .spearman(m$pct, y, MIN_DONORS)
                gr_rows[[cp$label]] <- data.frame(
                    dataset = ds, group = gn,
                    group_label = ARM_GROUPS[[gn]]$label,
                    variable = cp$label, role = cp$role,
                    expected_direction = cp$expected,
                    rho = st$rho, p_value = st$p, n_donors = st$n,
                    ties_present = st$ties, note = st$note,
                    stringsAsFactors = FALSE)
            }
            gdf <- do.call(rbind, gr_rows)
            gdf$p_adj <- NA_real_
            ok <- is.finite(gdf$p_value)
            if (any(ok)) gdf$p_adj[ok] <- p.adjust(gdf$p_value[ok],
                                                   method = "BH")
            gdf$significant <- is.finite(gdf$p_adj) & gdf$p_adj < 0.05
            gdf$matches_expected <- is.finite(gdf$rho) &
                sign(gdf$rho) == sign(gdf$expected_direction)
            rows[[gn]] <- gdf
            cat(sprintf("    %-17s : %d/%d significant after BH\n",
                        gn, sum(gdf$significant, na.rm = TRUE), nrow(gdf)))
        }
        res_df <- do.call(rbind, rows)
        rownames(res_df) <- NULL

        # == 2f: donor table for the multivariate models ======================
        # A donor with no cell in an arm carries 0.
        prop_df <- donor_totals
        # A group under the cell floor leaves its column all-NA.
        pred_avail <- character(0)
        for (nm in names(PRED_GROUPS)) {
            gn <- PRED_GROUPS[[nm]]
            if (is.null(grp_pct[[gn]])) { prop_df[[nm]] <- NA_real_; next }
            prop_df[[nm]] <- grp_pct[[gn]]$pct[
                match(prop_df$donor_id, grp_pct[[gn]]$donor_id)]
            pred_avail <- c(pred_avail, nm)
        }
        cat(sprintf("  Predictors available: %s\n",
                    if (length(pred_avail)) paste(unname(PRED_LABEL[pred_avail]),
                                                  collapse = ", ")
                    else "NONE"))
        prop_df <- merge(prop_df, donor_clin_all, by = "donor_id", all.x = TRUE)
        cat(sprintf("\n  Donor table for the models: %d donor(s)\n",
                    nrow(prop_df)))

        # == 2g: multivariate models ==========================================
        # outcome ~ %arm + covariates (OLS); delta R2 is the increment the arm
        # proportion adds over the covariates. OUTCOMES = CLINICAL_PARAMS plus
        # the cov4/cov7 covariates as negative-control outcomes; an outcome is
        # never its own covariate.
        .cov_keys_all <- unique(c(COV_SET_4, COV_SET_7))
        .excl_for <- function(cl) {
            hit <- names(COV_MAP)[vapply(COV_MAP, identical, logical(1), cl)]
            if (length(hit)) hit[1] else NULL
        }
        OUTCOMES <- lapply(CLINICAL_PARAMS, function(cp)
            list(label = cp$label, col = cp$col, excl = .excl_for(cp$col)))
        for (k in .cov_keys_all) {
            cl <- COV_MAP[[k]]
            if (is.null(cl)) next
            if (cl %in% vapply(OUTCOMES, function(z) z$col, character(1))) next
            OUTCOMES[[length(OUTCOMES) + 1L]] <-
                list(label = unname(COV_LABEL[k]), col = cl, excl = k)
        }
        OUTCOMES <- Filter(function(z) z$col %in% names(prop_df), OUTCOMES)
        scale_safe <- function(x) {
            s <- sd(x, na.rm = TRUE)
            if (!is.finite(s) || s < 1e-12) return(x - mean(x, na.rm = TRUE))
            (x - mean(x, na.rm = TRUE)) / s
        }
        fit_one <- function(covs, dat) {
            keep_cols <- c("y", "xpct", covs)
            d <- dat[, keep_cols, drop = FALSE]
            d <- d[complete.cases(d), , drop = FALSE]
            n <- nrow(d)
            if (n < MIN_MV_DONORS)
                return(list(n = n, beta = NA_real_, se = NA_real_,
                            t = NA_real_, p = NA_real_, r2_full = NA_real_,
                            r2_cov = NA_real_, delta_r2 = NA_real_,
                            vif_max = NA_real_,
                            vif_arm = NA_real_,
                            note = sprintf("n = %d < %d", n, MIN_MV_DONORS)))
            for (cn in c("y", "xpct", covs)) d[[cn]] <- scale_safe(d[[cn]])
            if (sd(d$xpct) < 1e-12)
                return(list(n = n, beta = NA_real_, se = NA_real_, t = NA_real_,
                            p = NA_real_, r2_full = NA_real_, r2_cov = NA_real_,
                            delta_r2 = NA_real_, vif_max = NA_real_,
                            vif_arm = NA_real_,
                            note = "zero variance in the arm proportion"))
            f_full <- stats::as.formula(
                paste("y ~ xpct", if (length(covs)) paste("+", paste(covs,
                      collapse = " + ")) else ""))
            f_cov  <- if (length(covs))
                stats::as.formula(paste("y ~", paste(covs, collapse = " + ")))
                else stats::as.formula("y ~ 1")
            m_full <- tryCatch(stats::lm(f_full, data = d),
                               error = function(e) NULL)
            m_cov  <- tryCatch(stats::lm(f_cov, data = d),
                               error = function(e) NULL)
            if (is.null(m_full) || is.null(m_cov))
                return(list(n = n, beta = NA_real_, se = NA_real_, t = NA_real_,
                            p = NA_real_, r2_full = NA_real_, r2_cov = NA_real_,
                            delta_r2 = NA_real_, vif_max = NA_real_,
                            vif_arm = NA_real_,
                            note = "model did not fit"))
            co <- summary(m_full)$coefficients
            if (!("xpct" %in% rownames(co)))
                return(list(n = n, beta = NA_real_, se = NA_real_, t = NA_real_,
                            p = NA_real_,
                            r2_full = summary(m_full)$r.squared,
                            r2_cov = summary(m_cov)$r.squared,
                            delta_r2 = NA_real_, vif_max = NA_real_,
                            vif_arm = NA_real_,
                            note = "arm proportion aliased by the covariates"))
            # Largest VIF and the arm term's own VIF are both reported.
            vif_max <- NA_real_
            vif_arm <- NA_real_
            if (length(covs) > 0) {
                vifs <- vapply(c("xpct", covs), function(v) {
                    others <- setdiff(c("xpct", covs), v)
                    if (length(others) == 0) return(1)
                    r2 <- suppressWarnings(summary(stats::lm(
                        stats::as.formula(paste(v, "~",
                                                paste(others, collapse = " + "))),
                        data = d))$r.squared)
                    if (!is.finite(r2) || r2 >= 1) return(Inf)
                    1 / (1 - r2)
                }, numeric(1))
                vif_max <- max(vifs, na.rm = TRUE)
                vif_arm <- unname(vifs[["xpct"]])
            }
            # Terms dropped by lm() as aliased are recorded on the row.
            n_dropped <- sum(is.na(stats::coef(m_full)))
            note <- character(0)
            if (n_dropped > 0)
                note <- c(note, sprintf("%d term(s) dropped by lm() as aliased",
                                        n_dropped))
            if (!is.finite(vif_max))
                note <- c(note, "perfect collinearity among the predictors")
            else if (vif_max >= 10)
                note <- c(note, sprintf("max VIF %.1f: the predictors are collinear and this coefficient is unstable",
                                        vif_max))
            list(n = n, beta = unname(co["xpct", 1]), se = unname(co["xpct", 2]),
                 t = unname(co["xpct", 3]), p = unname(co["xpct", 4]),
                 r2_full = summary(m_full)$r.squared,
                 r2_cov  = summary(m_cov)$r.squared,
                 delta_r2 = summary(m_full)$r.squared -
                            summary(m_cov)$r.squared,
                 vif_max = vif_max, vif_arm = vif_arm,
                 note = paste(note, collapse = "; "))
        }

        mv_rows <- list()
        # Tiers: "all" is primary; drop_severe (Severely Affected Donors left out)
        # and plus_lewy (Lewy body binary added) are pre-specified sensitivity.
        tiers <- list(
            list(id = "all",         label = "All donors"),
            list(id = "drop_severe", label = "Severely Affected Donors dropped"),
            list(id = "plus_lewy",   label = "All donors, + Lewy body binary"))
        # Reset per dataset.
        .cov_sets_seen <- NULL
        for (tr in tiers) {
            base <- prop_df
            if (tr$id == "drop_severe") {
                if (all(is.na(base$Severe_flag))) {
                    cat("  Tier drop_severe: no Severely Affected Donor flag;",
                        " skipped\n", sep = "")
                    next
                }
                base <- base[!is.na(base$Severe_flag) & base$Severe_flag == 0L, ,
                             drop = FALSE]
            }
            # A block is added only when every covariate it needs resolved with
            # variance.
            .set_ok <- function(keys) all(vapply(keys, function(k) {
                cl <- COV_MAP[[k]]
                is.null(cl) || cl %in% c("Age_num", "APOE_e4_dose", "Sex_num",
                                         "PMI_num", "RIN_num", "BrainPH_num",
                                         "BrainWeight_num") || cl %in% extra_ok
            }, logical(1)))
            .cov_sets <- list(list(id = "cov4", set = COV_SET_4),
                              list(id = "cov7", set = COV_SET_7))
            for (.b in list(list(id = "covM", set = COV_SET_M),
                            list(id = "covT", set = COV_SET_T),
                            list(id = "covP", set = COV_SET_P),
                            list(id = "covE", set = COV_SET_E))) {
                if (.set_ok(.b$set))
                    .cov_sets[[length(.cov_sets) + 1L]] <- .b
                else
                    cat(sprintf("    block %s skipped: a covariate is absent or constant\n",
                                .b$id))
            }
            # The head-to-head block reuses these blocks.
            if (identical(tr$id, "all")) .cov_sets_seen <- .cov_sets
            for (cs in .cov_sets) {
                covset <- cs$set
                if (tr$id == "plus_lewy") covset <- c(covset, "lewy")
                cov_cols <- vapply(covset, function(k)
                    if (identical(k, "lewy")) "Lewy_binary" else COV_MAP[[k]],
                    character(1))
                have <- covset[cov_cols %in% names(base)]
                if (length(have) < length(covset)) {
                    cat(sprintf(paste0("  Tier %s / %s: missing covariate ",
                                       "column(s) %s; skipped\n"),
                                tr$id, cs$id,
                                paste(setdiff(covset, have), collapse = ", ")))
                    next
                }
                for (oc in OUTCOMES) {
                    if (!(oc$col %in% names(base))) next
                    active <- setdiff(covset, oc$excl)
                    dat <- data.frame(y = suppressWarnings(
                        as.numeric(base[[oc$col]])))
                    for (k in active) {
                        cc <- if (identical(k, "lewy")) "Lewy_binary"
                              else COV_MAP[[k]]
                        dat[[k]] <- suppressWarnings(as.numeric(base[[cc]]))
                    }
                    for (a in pred_avail) {
                        dat$xpct <- suppressWarnings(as.numeric(base[[a]]))
                        f <- fit_one(active, dat)
                        mv_rows[[length(mv_rows) + 1L]] <- data.frame(
                            dataset = ds, tier = tr$id, tier_label = tr$label,
                            covariate_set = cs$id,
                            n_covariates = length(active),
                            covariates = paste(unname(COV_LABEL[active]),
                                               collapse = " + "),
                            outcome = oc$label,
                            predictor = a,
                            arm = unname(PRED_LABEL[a]),
                            n_donors = f$n, beta_std = f$beta, se = f$se,
                            t_value = f$t, p_value = f$p,
                            R2_full = f$r2_full, R2_covariates_only = f$r2_cov,
                            delta_R2 = f$delta_r2, max_VIF = f$vif_max,
                            arm_VIF = f$vif_arm,
                            note = f$note, stringsAsFactors = FALSE)
                    }
                }
            }
        }
        mv_df <- if (length(mv_rows)) do.call(rbind, mv_rows) else NULL
        if (is.null(.cov_sets_seen))
            .cov_sets_seen <- list(list(id = "cov4", set = COV_SET_4))
        scan_df <- NULL
        if (!is.null(mv_df)) {
            # BH across outcomes within each tier x covariate set x predictor.
            mv_df$p_adj <- NA_real_
            key <- paste(mv_df$tier, mv_df$covariate_set, mv_df$predictor,
                         sep = "|")
            for (k in unique(key)) {
                i <- key == k & is.finite(mv_df$p_value)
                if (any(i)) mv_df$p_adj[i] <- p.adjust(mv_df$p_value[i],
                                                       method = "BH")
            }
            mv_df$significant <- is.finite(mv_df$p_adj) & mv_df$p_adj < 0.05

            # == Leave-one-in scan ============================================
            # Each candidate is added alone to the 4-covariate core; the arm
            # coefficient is reported beside its core value.
            scan_rows <- list()
            .scan_keys <- COV_SCAN[vapply(COV_SCAN, function(k) {
                cl <- COV_MAP[[k]]
                !is.null(cl) && cl %in% names(prop_df) &&
                    sum(is.finite(suppressWarnings(
                        as.numeric(prop_df[[cl]])))) >= MIN_MV_DONORS &&
                    isTRUE(stats::sd(suppressWarnings(
                        as.numeric(prop_df[[cl]])), na.rm = TRUE) > 0)
            }, logical(1))]
            for (oc in CLINICAL_PARAMS) {
                if (!(oc$col %in% names(prop_df))) next
                # An outcome is never its own covariate (core included).
                .oc_excl <- .excl_for(oc$col)
                for (pr in pred_avail) {
                    dat0 <- data.frame(
                        y = suppressWarnings(as.numeric(prop_df[[oc$col]])),
                        xpct = prop_df[[pr]], stringsAsFactors = FALSE)
                    for (k in unique(c(COV_SET_4, .scan_keys))) {
                        cl <- COV_MAP[[k]]
                        if (is.null(cl) || !(cl %in% names(prop_df))) next
                        dat0[[k]] <- suppressWarnings(as.numeric(prop_df[[cl]]))
                    }
                    base_keys <- COV_SET_4[COV_SET_4 %in% names(dat0)]
                    if (!is.null(.oc_excl))
                        base_keys <- setdiff(base_keys, .oc_excl)
                    b0 <- fit_one(base_keys, dat0)
                    if (!is.finite(b0$beta)) next
                    for (k in .scan_keys) {
                        if (!(k %in% names(dat0))) next
                        # An outcome cannot be its own covariate.
                        if (identical(COV_MAP[[k]], oc$col)) next
                        f1 <- fit_one(unique(c(base_keys, k)), dat0)
                        scan_rows[[length(scan_rows) + 1L]] <- data.frame(
                            dataset = ds, outcome = oc$label,
                            arm = unname(PRED_LABEL[pr]),
                            added_covariate = unname(COV_LABEL[k]),
                            covariate_key = k,
                            n_core = b0$n, beta_core = round(b0$beta, 4),
                            p_core = signif(b0$p, 4),
                            n_added = f1$n, beta_added = round(f1$beta, 4),
                            p_added = signif(f1$p, 4),
                            # Percentage change only when |core beta| >=
                            # MV_MIN_BETA_FOR_PCT.
                            pct_change_beta =
                                if (isTRUE(abs(b0$beta) >= MV_MIN_BETA_FOR_PCT))
                                    round(100 * (f1$beta - b0$beta) /
                                          abs(b0$beta), 1)
                                else NA_real_,
                            abs_change_beta = round(f1$beta - b0$beta, 4),
                            pct_suppressed =
                                isTRUE(abs(b0$beta) < MV_MIN_BETA_FOR_PCT),
                            delta_r2_added = round(f1$delta_r2, 5),
                            max_VIF = round(f1$vif_max, 2),
                            # A sign flip counts only when |core beta| >=
                            # MV_MIN_BETA_FOR_PCT; the raw flip is kept beside it.
                            sign_flipped = isTRUE(is.finite(f1$beta) &&
                                abs(b0$beta) >= MV_MIN_BETA_FOR_PCT &&
                                sign(f1$beta) != sign(b0$beta)),
                            sign_flipped_any = isTRUE(is.finite(f1$beta) &&
                                sign(f1$beta) != sign(b0$beta)),
                            lost_significance = isTRUE(is.finite(b0$p) &&
                                is.finite(f1$p) && b0$p < 0.05 && f1$p >= 0.05),
                            stringsAsFactors = FALSE)
                    }
                }
            }
            scan_df <- if (length(scan_rows)) do.call(rbind, scan_rows) else NULL
            if (!is.null(scan_df)) {
                .mp <- scan_df[scan_df$arm == "Module Post-IP", , drop = FALSE]
                cat(sprintf(paste0("  Leave-one-in scan: %d model(s) over %d",
                                   " candidate(s); Module Post-IP loses",
                                   " significance in %d of %d\n"),
                            nrow(scan_df), length(.scan_keys),
                            sum(.mp$lost_significance, na.rm = TRUE), nrow(.mp)))
            }

            for (cid in c("cov4", "cov7", "covM", "covT", "covP", "covE")) {
                i <- mv_df$tier == "all" & mv_df$covariate_set == cid
                if (!any(i)) next
                cat(sprintf(paste0("  Multivariate %s (all donors): %d/%d ",
                                   "significant after BH, n = %s\n"),
                            cid, sum(mv_df$significant[i], na.rm = TRUE), sum(i),
                            paste(range(mv_df$n_donors[i], na.rm = TRUE),
                                  collapse = "-")))
            }
        }

        # == Predictor head-to-head ===========================================
        # Q1 (nested F) and Q2 (non-nested duel).
        ipsplit_df <- NULL
        duel_df    <- NULL
        thresh_df  <- NULL
        curve_df   <- NULL
        ident_ok   <- list()
        if (isTRUE(PHH_ENABLE) && !is.null(mv_df)) {
            .phh_base <- prop_df
            # One helper for the F test, the duel and every bootstrap replicate
            # (same scaling and complete-case set).
            .phh_r2 <- function(d, xs, covs) {
                cols <- c("y", xs, covs)
                d <- d[, cols, drop = FALSE]
                d <- d[complete.cases(d), , drop = FALSE]
                if (nrow(d) < PHH_MIN_N) return(NULL)
                for (cn in cols) d[[cn]] <- scale_safe(d[[cn]])
                if (any(vapply(xs, function(v) sd(d[[v]]) < 1e-12, logical(1))))
                    return(NULL)
                f <- stats::as.formula(paste("y ~",
                    paste(c(xs, covs), collapse = " + ")))
                f0 <- if (length(covs))
                    stats::as.formula(paste("y ~", paste(covs, collapse = " + ")))
                    else stats::as.formula("y ~ 1")
                m  <- tryCatch(stats::lm(f,  data = d), error = function(e) NULL)
                m0 <- tryCatch(stats::lm(f0, data = d), error = function(e) NULL)
                if (is.null(m) || is.null(m0)) return(NULL)
                list(m = m, m0 = m0, n = nrow(d), d = d,
                     r2 = summary(m)$r.squared,
                     r2_cov = summary(m0)$r.squared)
            }
            # Returns NULL when a covariate column is missing (the block is skipped).
            .phh_dat <- function(oc, covset, extra) {
                if (!(oc$col %in% names(.phh_base))) return(NULL)
                active <- setdiff(covset, oc$excl)
                cols <- vapply(active, function(k) {
                    cc <- COV_MAP[[k]]
                    if (is.null(cc)) NA_character_ else cc
                }, character(1))
                if (anyNA(cols) || !all(cols %in% names(.phh_base)))
                    return(NULL)
                if (!all(unlist(extra) %in% names(.phh_base))) return(NULL)
                d <- data.frame(y = suppressWarnings(
                    as.numeric(.phh_base[[oc$col]])))
                for (k in active) d[[k]] <- suppressWarnings(
                    as.numeric(.phh_base[[COV_MAP[[k]]]]))
                for (nm in names(extra)) d[[nm]] <- suppressWarnings(
                    as.numeric(.phh_base[[extra[[nm]]]]))
                list(d = d, active = active)
            }

            # == Q1: is a boundary worth its degree of freedom? ===============
            ips_rows <- list()
            for (rs in PHH_RESTRICTIONS) {
              if (!all(c(rs$unsplit, rs$split) %in% pred_avail)) {
                  cat(sprintf(paste0("  Restriction %s: skipped, %s not",
                                     " available\n"), rs$id,
                              paste(setdiff(c(rs$unsplit, rs$split),
                                            pred_avail), collapse = " and ")))
                  next
              }
              # The sum identity is verified on the donor columns first; if it
              # fails the models are not nested and the restriction is skipped.
              .sum_dev <- max(abs(.phh_base[[rs$unsplit]] -
                                  rowSums(.phh_base[, rs$split, drop = FALSE])),
                              na.rm = TRUE)
              ident_ok[[rs$id]] <- .sum_dev
              if (!is.finite(.sum_dev) || .sum_dev > 1e-8) {
                  cat(sprintf(paste0("  Restriction %s: SKIPPED -- %s is not",
                                     " the sum of %s (max |dev| %.3e), so the",
                                     " models are not nested\n"),
                              rs$id, rs$unsplit,
                              paste(rs$split, collapse = " + "), .sum_dev))
                  next
              }
              for (cs in .cov_sets_seen) {
                for (oc in OUTCOMES) {
                    if (!(oc$col %in% names(.phh_base))) next
                    z <- .phh_dat(oc, cs$set,
                                  stats::setNames(
                                      as.list(c(rs$unsplit, rs$split)),
                                      c("x_all", "x_pre", "x_post")))
                    if (is.null(z)) next
                    # One complete-case set and one standardisation for both
                    # models; x_all is the raw total standardised.
                    full <- .phh_r2(z$d, c("x_all", "x_pre", "x_post"),
                                    z$active)
                    if (is.null(full)) next
                    dd3 <- full$d
                    m_full <- tryCatch(stats::lm(stats::as.formula(paste(
                        "y ~", paste(c("x_pre", "x_post", z$active),
                                     collapse = " + "))), data = dd3),
                        error = function(e) NULL)
                    rest2 <- tryCatch(stats::lm(stats::as.formula(paste(
                        "y ~", paste(c("x_all", z$active), collapse = " + "))),
                        data = dd3), error = function(e) NULL)
                    m_cov <- tryCatch(stats::lm(if (length(z$active))
                        stats::as.formula(paste("y ~",
                            paste(z$active, collapse = " + ")))
                        else stats::as.formula("y ~ 1"), data = dd3),
                        error = function(e) NULL)
                    if (is.null(m_full) || is.null(rest2) || is.null(m_cov))
                        next
                    av <- tryCatch(stats::anova(rest2, m_full),
                                   error = function(e) NULL)
                    if (is.null(av) || nrow(av) < 2L) next
                    co <- summary(m_full)$coefficients
                    .cf <- function(v, j) if (v %in% rownames(co)) co[v, j]
                                          else NA_real_
                    ips_rows[[length(ips_rows) + 1L]] <- data.frame(
                        dataset = ds, test = rs$id, question = rs$question,
                        unsplit_predictor = unname(PRED_LABEL[rs$unsplit]),
                        split_predictors = paste(unname(PRED_LABEL[rs$split]),
                                                 collapse = " + "),
                        covariate_set = cs$id,
                        covariates = paste(unname(COV_LABEL[z$active]),
                                           collapse = " + "),
                        outcome = oc$label, n_donors = full$n,
                        R2_covariates_only = summary(m_cov)$r.squared,
                        R2_unsplit = summary(rest2)$r.squared,
                        R2_split   = summary(m_full)$r.squared,
                        delta_R2_unsplit = summary(rest2)$r.squared -
                                           summary(m_cov)$r.squared,
                        delta_R2_split   = summary(m_full)$r.squared -
                                           summary(m_cov)$r.squared,
                        gain_from_split  = summary(m_full)$r.squared -
                                           summary(rest2)$r.squared,
                        beta_split1 = .cf("x_pre", 1),
                        p_split1 = .cf("x_pre", 4),
                        beta_split2 = .cf("x_post", 1),
                        p_split2 = .cf("x_post", 4),
                        F_value = av[2, "F"], df = av[2, "Df"],
                        p_value = av[2, "Pr(>F)"],
                        stringsAsFactors = FALSE)
                }
              }
            }
            ipsplit_df <- if (length(ips_rows)) do.call(rbind, ips_rows)
                          else NULL
            if (!is.null(ipsplit_df)) {
                # BH within test x covariate set across outcomes.
                ipsplit_df$p_adj <- NA_real_
                .fam <- paste(ipsplit_df$test, ipsplit_df$covariate_set,
                              sep = "|")
                for (k in unique(.fam)) {
                    i <- .fam == k & is.finite(ipsplit_df$p_value)
                    if (any(i)) ipsplit_df$p_adj[i] <-
                        p.adjust(ipsplit_df$p_value[i], method = "BH")
                }
                ipsplit_df$significant <- is.finite(ipsplit_df$p_adj) &
                                          ipsplit_df$p_adj < 0.05
                for (tid in unique(ipsplit_df$test)) {
                    .i4 <- ipsplit_df$covariate_set == "cov4" &
                           ipsplit_df$test == tid
                    cat(sprintf(paste0("  Restriction %s (cov4): the split",
                                       " beats the unsplit predictor on %d of",
                                       " %d outcome(s) after BH\n"), tid,
                                sum(ipsplit_df$significant[.i4], na.rm = TRUE),
                                sum(.i4)))
                }
            }

            # == Q2: non-nested duel, with a donor bootstrap ==================
            duel_rows <- list()
            # The RNG stream is saved before set.seed(PHH_SEED) and restored
            # after the block.
            .rng_saved <- if (exists(".Random.seed", envir = globalenv()))
                get(".Random.seed", envir = globalenv()) else NULL
            set.seed(PHH_SEED)
            for (du in PHH_DUELS) {
                if (!(du$a %in% pred_avail) ||
                    !(du$a %in% names(.phh_base))) next
                b_col <- if (!is.null(du$b_pred)) {
                             if (!(du$b_pred %in% pred_avail)) next
                             du$b_pred
                         } else COV_MAP[[du$b_cov]]
                if (is.null(b_col) || !(b_col %in% names(.phh_base))) next
                b_lab <- if (!is.null(du$b_pred)) unname(PRED_LABEL[du$b_pred])
                         else unname(COV_LABEL[du$b_cov])
                for (oc in OUTCOMES) {
                    if (!(oc$label %in% PHH_OUTCOMES)) next
                    if (!(oc$col %in% names(.phh_base))) next
                    # A rival covariate in the core is removed from the core for
                    # this duel (recorded in self_excluded).
                    core <- setdiff(COV_SET_4, oc$excl)
                    self_excl <- character(0)
                    if (!is.null(du$b_cov) && du$b_cov %in% core) {
                        self_excl <- du$b_cov
                        core <- setdiff(core, du$b_cov)
                    }
                    z <- .phh_dat(oc, core, c(xa = du$a, xb = b_col))
                    if (is.null(z)) next
                    both <- .phh_r2(z$d, c("xa", "xb"), z$active)
                    if (is.null(both)) next
                    dj <- both$d                    # the shared donor set
                    fa <- .phh_r2(dj, "xa", z$active)
                    fb <- .phh_r2(dj, "xb", z$active)
                    if (is.null(fa) || is.null(fb)) next
                    da <- fa$r2 - fa$r2_cov
                    db <- fb$r2 - fb$r2_cov
                    # Donor bootstrap on the difference.
                    bs <- rep(NA_real_, PHH_BOOT)
                    for (b in seq_len(PHH_BOOT)) {
                        ii <- sample.int(nrow(dj), nrow(dj), replace = TRUE)
                        ra <- .phh_r2(dj[ii, , drop = FALSE], "xa", z$active)
                        rb <- .phh_r2(dj[ii, , drop = FALSE], "xb", z$active)
                        if (!is.null(ra) && !is.null(rb))
                            bs[b] <- (ra$r2 - ra$r2_cov) - (rb$r2 - rb$r2_cov)
                    }
                    bs <- bs[is.finite(bs)]
                    ci <- if (length(bs) >= 100L)
                        stats::quantile(bs, c(0.025, 0.975), names = FALSE)
                        else c(NA_real_, NA_real_)
                    # Vuong statistic on the two models fitted on `dj` (equal
                    # parameter counts, no penalty).
                    rA <- stats::residuals(fa$m); rB <- stats::residuals(fb$m)
                    nn <- length(rA)
                    s2A <- sum(rA^2) / nn; s2B <- sum(rB^2) / nn
                    fitcor <- suppressWarnings(stats::cor(
                        stats::fitted(fa$m), stats::fitted(fb$m)))
                    mvec <- if (s2A > 0 && s2B > 0)
                        (-0.5 * log(2 * pi * s2A) - rA^2 / (2 * s2A)) -
                        (-0.5 * log(2 * pi * s2B) - rB^2 / (2 * s2B))
                        else rep(NA_real_, nn)
                    msd <- if (all(is.finite(mvec)))
                        sqrt(mean((mvec - mean(mvec))^2)) else NA_real_
                    degenerate <- (is.finite(fitcor) &&
                                   abs(fitcor) > PHH_DUEL_FIT_COR) ||
                                  !is.finite(msd) || msd < 1e-12
                    vz <- if (degenerate) NA_real_
                          else sqrt(nn) * mean(mvec) / msd
                    vp <- if (is.finite(vz)) 2 * stats::pnorm(-abs(vz))
                          else NA_real_
                    co <- summary(both$m)$coefficients
                    .cf <- function(v, j) if (v %in% rownames(co)) co[v, j]
                                          else NA_real_
                    .vf <- function(v) {
                        others <- setdiff(c("xa", "xb", z$active), v)
                        if (!length(others)) return(1)
                        r2 <- suppressWarnings(summary(stats::lm(
                            stats::as.formula(paste(v, "~",
                                paste(others, collapse = " + "))),
                            data = dj))$r.squared)
                        if (!is.finite(r2) || r2 >= 1) Inf else 1 / (1 - r2)
                    }
                    duel_rows[[length(duel_rows) + 1L]] <- data.frame(
                        dataset = ds, duel = du$id, question = du$question,
                        outcome = oc$label, n_donors = both$n,
                        predictor_A = unname(PRED_LABEL[du$a]),
                        predictor_B = b_lab,
                        core_covariates = paste(unname(COV_LABEL[z$active]),
                                                collapse = " + "),
                        self_excluded = if (length(self_excl))
                            unname(COV_LABEL[self_excl]) else "",
                        delta_R2_A = da, delta_R2_B = db,
                        difference = da - db,
                        boot_lo = ci[1], boot_hi = ci[2],
                        boot_n = length(bs),
                        # Descriptive only; the verdict is assigned below.
                        favours_bootstrap = if (!is.finite(ci[1]))
                                                "undetermined"
                                            else if (ci[1] > 0) "A"
                                            else if (ci[2] < 0) "B"
                                            else "neither",
                        vuong_z = vz, vuong_p = vp,
                        models_indistinguishable = degenerate,
                        fitted_correlation = fitcor,
                        joint_beta_A = .cf("xa", 1), joint_p_A = .cf("xa", 4),
                        joint_beta_B = .cf("xb", 1), joint_p_B = .cf("xb", 4),
                        joint_VIF_A = .vf("xa"), joint_VIF_B = .vf("xb"),
                        stringsAsFactors = FALSE)
                }
            }
            if (!is.null(.rng_saved))
                assign(".Random.seed", .rng_saved, envir = globalenv())
            duel_df <- if (length(duel_rows)) do.call(rbind, duel_rows)
                       else NULL
            if (!is.null(duel_df)) {
                # BH within duel across outcomes.
                for (cn in c("vuong_p", "joint_p_A", "joint_p_B")) {
                    duel_df[[paste0(cn, "_adj")]] <- NA_real_
                    for (k in unique(duel_df$duel)) {
                        i <- duel_df$duel == k & is.finite(duel_df[[cn]])
                        if (any(i)) duel_df[[paste0(cn, "_adj")]][i] <-
                            p.adjust(duel_df[[cn]][i], method = "BH")
                    }
                }
                # Verdict from the encompassing model: "A", "B", "both" or
                # "neither".
                duel_df$favours <- ifelse(
                    duel_df$models_indistinguishable, "indistinguishable",
                    ifelse(!is.finite(duel_df$joint_p_A_adj) |
                           !is.finite(duel_df$joint_p_B_adj), "undetermined",
                    ifelse(duel_df$joint_p_A_adj < 0.05 &
                           duel_df$joint_p_B_adj < 0.05, "both",
                    ifelse(duel_df$joint_p_A_adj < 0.05, "A",
                    ifelse(duel_df$joint_p_B_adj < 0.05, "B", "neither")))))
                # Vuong verdict, reported only.
                duel_df$favours_vuong <- ifelse(
                    duel_df$models_indistinguishable, "indistinguishable",
                    ifelse(!is.finite(duel_df$vuong_p_adj), "undetermined",
                    ifelse(duel_df$vuong_p_adj >= 0.05, "neither",
                    ifelse(duel_df$vuong_z > 0, "A", "B"))))
            }
            if (!is.null(duel_df)) for (dn in unique(duel_df$duel)) {
                i <- duel_df$duel == dn
                cat(sprintf(paste0("  Duel %s (encompassing): module only %d,",
                                   " rival only %d, both %d, neither %d of",
                                   " %d\n"),
                            dn, sum(duel_df$favours[i] == "A"),
                            sum(duel_df$favours[i] == "B"),
                            sum(duel_df$favours[i] == "both"),
                            sum(!(duel_df$favours[i] %in%
                                  c("A", "B", "both"))), sum(i)))
            }
        }

        # == Placebo-threshold profile ========================================
        # The split is recomputed from cell-level cDFG at a grid of thresholds,
        # refitted at each, and the gain in R2 read off as a curve.
        if (isTRUE(PHH_ENABLE) && isTRUE(PHH_THRESH_ENABLE) &&
            !is.null(mv_df) && is.finite(ip_y) &&
            all(c("pct_module_all", "pct_pre", "pct_post") %in% pred_avail)) {
            .mod <- rc$cell_class %in% MODULE_CLASSES & !is.na(rc$donor_id) &
                    is.finite(rc$cDFG_expr)
            if (sum(.mod) >= MIN_CELLS_GRP) {
                .cd  <- rc$cDFG_expr[.mod]
                .dn  <- as.character(rc$donor_id[.mod])
                # Percentages at threshold t on the same denominator; a donor
                # with no module cell carries 0.
                .pct_at <- function(t) {
                    ip <- .cd >= t
                    a <- table(factor(.dn[ip],  levels = donor_totals$donor_id))
                    b <- table(factor(.dn[!ip], levels = donor_totals$donor_id))
                    data.frame(
                        donor_id = donor_totals$donor_id,
                        x_post = 100 * as.integer(a) / donor_totals$n_total,
                        x_pre  = 100 * as.integer(b) / donor_totals$n_total,
                        stringsAsFactors = FALSE)
                }
                # Grid in quantiles of the module cells' cDFG; IP_y placed by its
                # empirical quantile.
                .grid  <- unname(stats::quantile(.cd, PHH_THRESH_GRID,
                                                 na.rm = TRUE, names = FALSE))
                .ip_q  <- mean(.cd < ip_y, na.rm = TRUE)
                .keep  <- abs(PHH_THRESH_GRID - .ip_q) > PHH_THRESH_EXCL
                .base2 <- merge(prop_df, .pct_at(ip_y), by = "donor_id",
                                all.x = TRUE)
                # Recomputing at IP_y must reproduce pct_post.
                .recon_dev <- max(abs(.base2$x_post - .base2$pct_post),
                                  na.rm = TRUE)
                th_rows <- list()
                # The full curve, long format.
                curve_rows <- list()
                .cov4 <- COV_SET_4
                for (oc in OUTCOMES) {
                    if (!(oc$label %in% PHH_OUTCOMES)) next
                    z0 <- .phh_dat(oc, .cov4, c(x_all = "pct_module_all"))
                    if (is.null(z0)) next
                    r0 <- .phh_r2(z0$d, "x_all", z0$active)
                    if (is.null(r0)) next
                    keep_id <- prop_df$donor_id[
                        stats::complete.cases(z0$d[, c("y", "x_all",
                                                       z0$active),
                                                   drop = FALSE])]
                    .gain <- function(t) {
                        pc <- .pct_at(t)
                        d  <- z0$d
                        d$x_pre  <- pc$x_pre[match(prop_df$donor_id,
                                                   pc$donor_id)]
                        d$x_post <- pc$x_post[match(prop_df$donor_id,
                                                    pc$donor_id)]
                        d <- d[prop_df$donor_id %in% keep_id, , drop = FALSE]
                        if (sd(d$x_post, na.rm = TRUE) < 1e-12) return(NA_real_)
                        f <- tryCatch(stats::lm(stats::as.formula(paste(
                            "y ~", paste(c("x_pre", "x_post", z0$active),
                                         collapse = " + "))), data = d),
                            error = function(e) NULL)
                        if (is.null(f)) return(NA_real_)
                        summary(f)$r.squared - r0$r2
                    }
                    g_grid <- vapply(.grid, .gain, numeric(1))
                    g_ip   <- .gain(ip_y)
                    ok <- is.finite(g_grid)
                    if (!is.finite(g_ip) || !any(ok)) next
                    # +1 in numerator and denominator (permutation p).
                    p_all <- (1 + sum(g_grid[ok] >= g_ip)) / (1 + sum(ok))
                    okf   <- ok & .keep
                    p_far <- if (any(okf))
                        (1 + sum(g_grid[okf] >= g_ip)) / (1 + sum(okf))
                        else NA_real_
                    th_rows[[length(th_rows) + 1L]] <- data.frame(
                        dataset = ds, outcome = oc$label, n_donors = r0$n,
                        IP_y = ip_y, IP_quantile = .ip_q,
                        gain_at_IP = g_ip,
                        best_gain = max(g_grid[ok]),
                        best_quantile = PHH_THRESH_GRID[ok][
                            which.max(g_grid[ok])],
                        best_threshold = .grid[ok][which.max(g_grid[ok])],
                        n_grid = sum(ok), n_grid_far = sum(okf),
                        p_all = p_all, p_far = p_far,
                        recon_dev_at_IP = .recon_dev,
                        stringsAsFactors = FALSE)
                    # One row per grid point plus one row for IP_y (is_IP);
                    # in_far_null marks the points p_far counted.
                    curve_rows[[length(curve_rows) + 1L]] <- data.frame(
                        dataset = ds, outcome = oc$label, n_donors = r0$n,
                        quantile = c(PHH_THRESH_GRID[ok], .ip_q),
                        threshold = c(.grid[ok], ip_y),
                        gain = c(g_grid[ok], g_ip),
                        is_IP = c(rep(FALSE, sum(ok)), TRUE),
                        in_far_null = c(.keep[ok], FALSE),
                        R2_covariates_only = r0$r2_cov,
                        R2_unsplit = r0$r2,
                        stringsAsFactors = FALSE)
                }
                thresh_df <- if (length(th_rows)) do.call(rbind, th_rows)
                             else NULL
                curve_df <- if (length(curve_rows))
                                do.call(rbind, curve_rows) else NULL
                if (!is.null(thresh_df)) {
                    thresh_df$p_far_adj <- NA_real_
                    i <- is.finite(thresh_df$p_far)
                    if (any(i)) thresh_df$p_far_adj[i] <-
                        p.adjust(thresh_df$p_far[i], method = "BH")
                    cat(sprintf(paste0("  Threshold profile: IP_y at quantile",
                                       " %.3f | median p_far %.3f over %d",
                                       " outcome(s) | recomputed pct_post at",
                                       " IP_y vs stored: max |dev| %.2e\n"),
                                .ip_q,
                                stats::median(thresh_df$p_far, na.rm = TRUE),
                                nrow(thresh_df), .recon_dev))
                }
            }
        }

        # == 2g2: Steiger, dependent correlations =============================
        # Module against non-module pairs share donors, supertype and
        # denominator.
        st_rows <- list()
        for (pr in STEIGER_PAIRS) {
            if (!all(pr %in% pred_avail)) {
                cat(sprintf("  Steiger %s vs %s: skipped, %s not available\n",
                            PRED_LABEL[pr[1]], PRED_LABEL[pr[2]],
                            paste(unname(PRED_LABEL[setdiff(pr, pred_avail)]),
                                  collapse = " and ")))
                next
            }
            for (cp in CLINICAL_PARAMS) {
                if (!(cp$col %in% names(prop_df))) next
                y  <- suppressWarnings(as.numeric(prop_df[[cp$col]]))
                x1 <- suppressWarnings(as.numeric(prop_df[[pr[1]]]))
                x2 <- suppressWarnings(as.numeric(prop_df[[pr[2]]]))
                # All three correlations on the same donors.
                keep <- is.finite(y) & is.finite(x1) & is.finite(x2)
                n <- sum(keep)
                s1  <- .spearman(x1[keep], y[keep], MIN_DONORS)
                s2  <- .spearman(x2[keep], y[keep], MIN_DONORS)
                s12 <- .spearman(x1[keep], x2[keep], MIN_DONORS)
                st  <- steiger_z(s1$rho, s2$rho, s12$rho, n)
                st_rows[[length(st_rows) + 1L]] <- data.frame(
                    dataset = ds, variable = cp$label, role = cp$role,
                    predictor_1 = unname(PRED_LABEL[pr[1]]),
                    predictor_2 = unname(PRED_LABEL[pr[2]]),
                    rho_1 = s1$rho, rho_2 = s2$rho,
                    rho_between_predictors = s12$rho,
                    delta_rho = st$delta_rho, Z = st$Z, p_value = st$p,
                    n_donors = n,
                    test = "Meng-Rosenthal-Rubin z, dependent overlapping",
                    spearman_variance_factor = SPEARMAN_Z_VAR_FACTOR,
                    note = st$note, stringsAsFactors = FALSE)
            }
        }
        st_df <- if (length(st_rows)) do.call(rbind, st_rows) else NULL
        if (!is.null(st_df)) {
            st_df$p_adj <- NA_real_
            kk <- paste(st_df$predictor_1, st_df$predictor_2, sep = " vs ")
            for (k1 in unique(kk)) {
                i <- kk == k1 & is.finite(st_df$p_value)
                if (any(i)) st_df$p_adj[i] <- p.adjust(st_df$p_value[i],
                                                       method = "BH")
            }
            st_df$significant <- is.finite(st_df$p_adj) & st_df$p_adj < 0.05
            for (k1 in unique(kk))
                cat(sprintf("  Steiger %-42s : %d/%d significant after BH\n",
                            k1, sum(st_df$significant[kk == k1], na.rm = TRUE),
                            sum(kk == k1)))
        }

        # == 2g3: Jonckheere-Terpstra trend test ==============================
        # Donors cut into tertiles of their Post-IP fraction; the non-module
        # ordering is the control.
        jt_rows <- list(); jt_donor <- list()
        jt_orderings <- list(
            list(id = "module",   classes = MODULE_CLASSES,
                 label = "Module post-IP fraction"),
            list(id = "excluded", classes = NONMODULE_CLASSES,
                 label = "Non-module post-IP fraction"))
        for (jo in jt_orderings) {
            sel  <- rc$cell_class %in% jo$classes & !is.na(rc$donor_id) &
                    !is.na(rc$arm) & rc$arm %in% c("Pre-IP", "Post-IP")
            if (!any(sel)) {
                cat(sprintf("  JT %s: no cell carries both a donor and an arm; skipped\n",
                            jo$id)); next
            }
            dc <- data.frame(
                donor_id = names(table(rc$donor_id[sel])),
                n_pre  = as.integer(table(factor(rc$donor_id[sel & rc$arm == "Pre-IP"],
                                                 levels = names(table(rc$donor_id[sel]))))),
                n_post = as.integer(table(factor(rc$donor_id[sel & rc$arm == "Post-IP"],
                                                 levels = names(table(rc$donor_id[sel]))))),
                stringsAsFactors = FALSE)
            dc$n_cells  <- dc$n_pre + dc$n_post
            # Per-donor cell floor, reported before and after.
            n_donor_all  <- nrow(dc)
            n_donor_thin <- sum(dc$n_cells < JT_MIN_CELLS_PER_DONOR)
            dc <- dc[dc$n_cells >= JT_MIN_CELLS_PER_DONOR, , drop = FALSE]
            if (n_donor_thin > 0)
                cat(sprintf(paste0("  JT %s: %d of %d donor(s) dropped for ",
                                   "fewer than %d cell(s)\n"),
                            jo$id, n_donor_thin, n_donor_all,
                            JT_MIN_CELLS_PER_DONOR))
            dc$ip_ratio <- dc$n_post / dc$n_cells
            if (nrow(dc) < JT_MIN_DONORS) {
                cat(sprintf(paste0("  JT %s: %d donor(s) clear the %d-cell ",
                                   "floor, below the %d-donor floor; skipped\n"),
                            jo$id, nrow(dc), JT_MIN_CELLS_PER_DONOR,
                            JT_MIN_DONORS)); next
            }
            br <- unique(stats::quantile(dc$ip_ratio, probs = c(0, 1/3, 2/3, 1),
                                         na.rm = TRUE))
            if (length(br) < 3L) {
                cat(sprintf(paste0("  JT %s: ip_ratio has %d distinct tertile ",
                                   "boundary(ies); not enough variation to ",
                                   "order donors. Skipped\n"),
                            jo$id, length(br))); next
            }
            grp <- cut(dc$ip_ratio, breaks = br, include.lowest = TRUE,
                       labels = FALSE)
            # Tertile labels relabelled to a contiguous 1..k in the same order.
            keep <- !is.na(grp)
            dc <- dc[keep, , drop = FALSE]; grp <- grp[keep]
            lv  <- sort(unique(grp))
            dc$jt_group <- match(grp, lv)
            k_g <- length(lv)
            if (k_g < 2L) {
                cat(sprintf("  JT %s: donors fall in one tertile only; skipped\n",
                            jo$id)); next
            }
            cat(sprintf("  JT %s -- tertiles of %s, %d donor(s):\n",
                        jo$id, jo$label, nrow(dc)))
            for (q in seq_len(k_g)) {
                dq <- dc[dc$jt_group == q, ]
                cat(sprintf("    T%d: n=%2d, ip_ratio [%.3f, %.3f]\n", q,
                            nrow(dq), min(dq$ip_ratio), max(dq$ip_ratio)))
            }
            jd <- merge(dc, donor_clin_all, by = "donor_id", all.x = TRUE)
            jd$ordering <- jo$id
            jt_donor[[jo$id]] <- jd
            for (cp in CLINICAL_PARAMS) {
                if (!(cp$col %in% names(jd))) next
                yv <- suppressWarnings(as.numeric(jd[[cp$col]]))
                o  <- jt_test(yv, jd$jt_group, k_g)
                meds <- vapply(seq_len(k_g), function(q) {
                    v <- yv[jd$jt_group == q & is.finite(yv)]
                    if (length(v)) stats::median(v) else NA_real_
                }, numeric(1))
                jt_rows[[length(jt_rows) + 1L]] <- data.frame(
                    dataset = ds, ordering = jo$id, ordering_label = jo$label,
                    variable = cp$label, role = cp$role,
                    expected_direction = cp$expected,
                    n_donors = o$n, n_groups = k_g,
                    n_donors_before_cell_floor = n_donor_all,
                    n_donors_dropped_cell_floor = n_donor_thin,
                    min_cells_per_donor = JT_MIN_CELLS_PER_DONOR,
                    group_sizes = paste(o$nk, collapse = "/"),
                    group_medians = paste(round(meds, 3), collapse = "/"),
                    n_tie_groups = o$n_tie_groups,
                    JT_statistic = o$T, JT_expected = o$E,
                    JT_variance_tie_corrected = o$V,
                    JT_Z = o$Z, JT_p = o$p, note = o$note,
                    stringsAsFactors = FALSE)
            }
        }
        jt_df <- if (length(jt_rows)) do.call(rbind, jt_rows) else NULL
        if (!is.null(jt_df)) {
            jt_df$JT_p_adj <- NA_real_
            for (od in unique(jt_df$ordering)) {
                i <- jt_df$ordering == od & is.finite(jt_df$JT_p)
                if (any(i)) jt_df$JT_p_adj[i] <- p.adjust(jt_df$JT_p[i],
                                                          method = "BH")
            }
            jt_df$significant <- is.finite(jt_df$JT_p_adj) & jt_df$JT_p_adj < 0.05
            for (od in unique(jt_df$ordering))
                cat(sprintf("  JT %-9s : %d/%d significant after BH\n", od,
                            sum(jt_df$significant[jt_df$ordering == od],
                                na.rm = TRUE),
                            sum(jt_df$ordering == od)))
        }

        # == 2h: figures ======================================================
        save_pp <- function(p, base, w = 12, h = 8, dpi = 200) {
            # Keep the plot object as well as the file (results RDS).
            if (isTRUE(SAVE_FIGURE_OBJECTS))
                .fig_store[[base]] <<- list(plot = p, width = w, height = h,
                                            dpi = dpi, dataset = ds)
            # Each render is wrapped; the device is closed by on.exit registered
            # after it opens.
            render <- function(open_fn, path) {
                tryCatch({
                    open_fn()
                    on.exit(grDevices::dev.off(), add = TRUE)
                    print(p)
                    TRUE
                }, error = function(e) {
                    cat(sprintf("    figure not written (%s): %s\n",
                                basename(path), conditionMessage(e)))
                    FALSE
                })
            }
            pdf_path <- file.path(OUT_DIR, paste0(base, ".pdf"))
            png_path <- file.path(OUT_DIR, paste0(base, ".png"))
            render(function() grDevices::pdf(pdf_path, width = w, height = h),
                   pdf_path)
            render(function() grDevices::png(png_path, width = w * dpi,
                                             height = h * dpi, res = dpi),
                   png_path)
            invisible(NULL)
        }

        # Provenance caption carried by every figure.
        fig_note <- sprintf(paste0(
            "Arms inherited: Pre-IP when cDFG < IP_y = %s; ",
            "never split on GOI. Reported IP x_IP = %s.\n",
            "Denominator: all %s cells per donor in %s. Arm: %s. ",
            "Spearman rho, BH-FDR within group across %d variables."),
            if (is.finite(ip_y)) sprintf("%.4f", ip_y) else "NA",
            if (is.finite(ip_x)) sprintf("%.4f", ip_x) else "NA",
            PRIMARY_STYPE, ds, ARM_TAG, length(CLINICAL_PARAMS))

        # Heatmap: groups x clinical variables. Fill capped at +-0.8 (the
        # printed number is the true rho); significant cells bold.
        hm <- res_df
        hm$group_label <- factor(hm$group_label,
                                 levels = vapply(ARM_GROUPS, `[[`,
                                                 character(1), "label"))
        hm$variable <- factor(hm$variable, levels = rev(VAR_ORDER_DISPLAY))
        hm$star <- .sig_star(hm$p_adj)
        hm$rho_cap <- pmax(pmin(hm$rho, 0.8), -0.8)
        hm$is_sig <- nzchar(hm$star)
        hm$lab <- ifelse(is.finite(hm$rho),
                         ifelse(hm$is_sig,
                                sprintf("%+.2f\n%s", hm$rho, hm$star),
                                sprintf("%+.2f", hm$rho)), "")
        .hm_breaks <- c(-0.8, -0.4, 0, 0.4, 0.8)
        .hm_cols   <- c("#053061", "#4393C3", "#F7F7F7", "#D6604D", "#67001F")
        p_hm <- ggplot(hm[!is.na(hm$group_label), ],
                       aes(x = group_label, y = variable, fill = rho_cap)) +
            geom_tile(colour = "white", linewidth = 0.5) +
            geom_text(aes(label = lab,
                          fontface = ifelse(is_sig, "bold", "plain")),
                      size = 2.8, lineheight = 0.85) +
            scale_fill_gradientn(
                colours  = .hm_cols,
                values   = scales::rescale(.hm_breaks, from = c(-0.8, 0.8)),
                limits   = c(-0.8, 0.8),
                na.value = "gray92",
                name     = "Spearman rho\n(fill capped \u00b10.8)",
                guide    = guide_colorbar(barwidth = 1.2, barheight = 8,
                                          ticks.colour = "black")) +
            labs(title = sprintf("Donor-level Spearman rho by CMRF arm group - %s",
                                 ds),
                 subtitle = paste0("Cell % of each group per donor vs clinical ",
                                   "variable. * BH q<0.05, ** q<0.01, *** q<0.001.\n",
                                   "MMSE and age are NOT negative controls: ",
                                   "MMSE falls with AD and age is a risk factor."),
                 caption = fig_note, x = NULL, y = NULL) +
            theme_h10() +
            theme(axis.text.x = element_text(angle = 30, hjust = 1),
                  panel.border = element_rect(colour = "grey40", fill = NA,
                                              linewidth = 0.5))
        save_pp(p_hm, paste0(FILE_PREFIX, "_", ds, "_Heatmap"), w = 11, h = 7)

        # Scatter grid: the module core arms against every clinical variable.
        sc <- list()
        for (cp in CLINICAL_PARAMS) {
            if (!(cp$col %in% names(prop_df))) next
            for (a in pred_avail) {
                d <- data.frame(x = prop_df[[a]],
                                y = suppressWarnings(as.numeric(prop_df[[cp$col]])))
                d <- d[is.finite(d$x) & is.finite(d$y), , drop = FALSE]
                if (nrow(d) < MIN_DONORS) next
                st <- .spearman(d$x, d$y, MIN_DONORS)
                sc[[length(sc) + 1L]] <- ggplot(d, aes(x, y)) +
                    geom_point(alpha = 0.65, size = 1.4,
                               colour = PRED_COL[a]) +
                    geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                                linewidth = 0.5, colour = "gray35") +
                    labs(title = sprintf("%s  |  %s", cp$label,
                                         unname(PRED_LABEL[a])),
                         subtitle = sprintf("rho = %.3f, p = %.3g, n = %d",
                                            st$rho, st$p, st$n),
                         x = "% of donor's cells", y = cp$label) +
                    theme_h10(base_size = 8)
            }
        }
        if (length(sc) > 0) {
            p_sc <- patchwork::wrap_plots(sc, ncol = 4) +
                patchwork::plot_annotation(
                    title = sprintf("Donor scatter - module and excluded arms vs clinical variables (%s)",
                                    ds),
                    caption = fig_note,
                    theme = theme(plot.caption = element_text(size = 7,
                                                              colour = "gray35",
                                                              hjust = 0)))
            save_pp(p_sc, paste0(FILE_PREFIX, "_", ds, "_Scatter"),
                    w = 16, h = 3.2 * ceiling(length(sc) / 4))
        }

        # Dumbbell per Steiger pair; q is the dependent-correlation test.
        if (!is.null(st_df)) {
            for (pr in STEIGER_PAIRS) {
                l1 <- unname(PRED_LABEL[pr[1]]); l2 <- unname(PRED_LABEL[pr[2]])
                dsub <- st_df[st_df$predictor_1 == l1 &
                              st_df$predictor_2 == l2, , drop = FALSE]
                dsub <- dsub[is.finite(dsub$rho_1) & is.finite(dsub$rho_2), ,
                             drop = FALSE]
                if (nrow(dsub) == 0L) next
                dsub$variable <- factor(dsub$variable,
                                        levels = rev(VAR_ORDER_DISPLAY))
                dsub$lab <- ifelse(is.finite(dsub$p_adj),
                                   sprintf("q=%.3g", dsub$p_adj), "")
                lg <- rbind(
                    data.frame(variable = dsub$variable, rho = dsub$rho_1,
                               which = l1, stringsAsFactors = FALSE),
                    data.frame(variable = dsub$variable, rho = dsub$rho_2,
                               which = l2, stringsAsFactors = FALSE))
                p_db <- ggplot() +
                    geom_segment(data = dsub,
                                 aes(y = variable, yend = variable,
                                     x = rho_2, xend = rho_1),
                                 colour = "gray65", linewidth = 0.9) +
                    geom_point(data = lg,
                               aes(x = rho, y = variable, colour = which),
                               size = 3) +
                    geom_text(data = dsub,
                              aes(x = pmax(rho_1, rho_2) + 0.04, y = variable,
                                  label = lab),
                              hjust = 0, size = 2.7, colour = "gray25") +
                    geom_vline(xintercept = 0, colour = "gray40",
                               linewidth = 0.4) +
                    scale_colour_manual(values = setNames(
                        c(PRED_COL[[pr[1]]], PRED_COL[[pr[2]]]), c(l1, l2)),
                        name = NULL) +
                    scale_x_continuous(limits = c(-1, 1.25),
                                       breaks = seq(-1, 1, 0.5)) +
                    labs(title = sprintf("%s vs %s - %s", l1, l2, ds),
                         subtitle = paste0("Dependent correlations sharing the ",
                                           "clinical variable; q is the ",
                                           "Meng-Rosenthal-Rubin z after BH ",
                                           "within this pair.\nSpearman rho is ",
                                           "not Pearson r, so the Fisher ",
                                           "variance carries the 1.06 Spearman ",
                                           "factor. An approximation either way."),
                         caption = fig_note, x = expression(rho), y = NULL) +
                    theme_h10()
                save_pp(p_db, sprintf("%s_%s_Dumbbell_%s_vs_%s", FILE_PREFIX,
                                      ds, pr[1], pr[2]), w = 10, h = 6)
            }
        }

        # JT panel: one panel per clinical variable, tertiles on x; ordinals as
        # stacked percentage bars, continuous variables as boxplots. The
        # annotation is the trend p.
        if (!is.null(jt_df) && length(jt_donor) > 0) {
            for (od in names(jt_donor)) {
                jd  <- jt_donor[[od]]
                sub <- jt_df[jt_df$ordering == od, , drop = FALSE]
                if (is.null(jd) || nrow(sub) == 0L) next
                .g <- suppressWarnings(as.integer(jd$jt_group))
                n_grp <- length(unique(.g[is.finite(.g)]))
                if (!is.finite(n_grp) || n_grp < 2L) next
                q_lab  <- paste0("T", seq_len(n_grp))
                t_disp <- .jt_x_labels(n_grp)
                panels <- list()
                for (cp in CLINICAL_PARAMS) {
                    if (!(cp$col %in% names(jd))) next
                    fr <- sub[sub$variable == cp$label, , drop = FALSE]
                    fq <- if (nrow(fr)) fr$JT_p_adj[1] else NA_real_
                    fz <- if (nrow(fr)) fr$JT_Z[1] else NA_real_
                    flab <- if (is.finite(fq) && is.finite(fz))
                        sprintf("JT Z=%+.2f  FDR=%s %s", fz,
                                if (fq < 0.001) sprintf("%.3e", fq)
                                else sprintf("%.4f", fq), .sig_star(fq))
                    else "not testable"
                    bdf <- data.frame(
                        y = suppressWarnings(as.numeric(jd[[cp$col]])),
                        group = factor(paste0("T", .g), levels = q_lab),
                        stringsAsFactors = FALSE)
                    bdf <- bdf[is.finite(bdf$y) & !is.na(bdf$group), ,
                               drop = FALSE]
                    if (nrow(bdf) == 0L) next
                    if (cp$label %in% JT_ORDINAL_VARS) {
                        bdf$y <- as.integer(round(bdf$y))
                        sev <- sort(unique(bdf$y))
                        tb <- as.data.frame(table(group = bdf$group, y = bdf$y),
                                            stringsAsFactors = FALSE)
                        tot <- stats::aggregate(list(tot = tb$Freq),
                                                by = list(group = tb$group),
                                                FUN = sum)
                        tb$pct <- 100 * tb$Freq /
                            tot$tot[match(tb$group, tot$group)]
                        tb <- tb[is.finite(tb$pct), , drop = FALSE]
                        tb$y     <- factor(as.integer(tb$y), levels = sev)
                        tb$group <- factor(tb$group, levels = q_lab)
                        sc <- .jt_sev_palette(length(sev))
                        names(sc) <- as.character(sev)
                        p <- ggplot(tb, aes(x = group, y = pct, fill = y)) +
                            geom_col(width = 0.7, colour = "white",
                                     linewidth = 0.3) +
                            scale_fill_manual(
                                values = sc, name = cp$label,
                                guide = guide_legend(reverse = TRUE, nrow = 1,
                                                     keywidth = 0.8,
                                                     keyheight = 0.8)) +
                            scale_y_continuous(
                                labels = function(x) paste0(x, "%"),
                                expand = expansion(mult = c(0, 0.14))) +
                            annotate("text", x = n_grp / 2 + 0.5, y = 105,
                                     label = flab, size = 2.8,
                                     colour = "grey20", fontface = "bold") +
                            labs(title = cp$label, x = NULL, y = "% of donors")
                    } else {
                        p <- ggplot(bdf, aes(x = group, y = y, fill = group)) +
                            geom_boxplot(width = 0.6, outlier.size = 1.2,
                                         alpha = 0.85, linewidth = 0.4,
                                         colour = "grey30") +
                            geom_jitter(width = 0.15, size = 0.9, alpha = 0.4,
                                        colour = "grey30") +
                            scale_fill_manual(values = JT_TERTILE_FILL,
                                              guide = "none") +
                            annotate("text", x = n_grp / 2 + 0.5, y = Inf,
                                     label = flab, vjust = 1.5, size = 2.8,
                                     colour = "grey20", fontface = "bold") +
                            labs(title = cp$label, x = NULL, y = cp$label)
                    }
                    p <- p +
                        scale_x_discrete(labels = t_disp) +
                        theme_bw(base_size = 9) +
                        theme(plot.title = element_text(face = "bold",
                                                        size = 10, hjust = 0.5),
                              legend.position = "bottom",
                              legend.title = element_text(size = 7),
                              legend.text  = element_text(size = 7),
                              legend.margin = margin(0, 0, 0, 0),
                              panel.grid.minor   = element_blank(),
                              panel.grid.major.x = element_blank(),
                              plot.background  = element_rect(fill = "white",
                                                              colour = NA),
                              panel.background = element_rect(fill = "white",
                                                              colour = NA))
                    panels[[length(panels) + 1L]] <- p
                }
                if (length(panels) < 2L) next
                n_rows_fig <- ceiling(length(panels) / 4)
                fig <- patchwork::wrap_plots(panels, ncol = 4) +
                    patchwork::plot_annotation(
                        title = sprintf(paste0("JT trend test: clinical scores",
                                               " by post-IP tertile \u2014 %s (%s)"),
                                        ds, od),
                        subtitle = paste0(
                            "Donors ordered by the ",
                            if (od == "module") "MODULE" else "EXCLUDED",
                            " post-IP fraction; positive Z = scores rise from T1",
                            " to T", n_grp, ".  FDR is BH within this ordering.\n",
                            "The variance carries the Lehmann tie correction, ",
                            "which these short ordinals need. Donors below the ",
                            "cell floor are dropped."),
                        caption = fig_note,
                        theme = theme(
                            plot.title    = element_text(face = "bold",
                                                         size = 13, hjust = 0.5),
                            plot.subtitle = element_text(size = 8.5,
                                                         colour = "grey40",
                                                         hjust = 0.5),
                            plot.caption  = element_text(size = 7,
                                                         colour = "grey55"),
                            plot.background = element_rect(fill = "white",
                                                           colour = NA)))
                save_pp(fig, sprintf("%s_%s_JT_%s", FILE_PREFIX, ds, od),
                        w = 14, h = max(7, n_rows_fig * 5 + 2))
            }
        }

        # == Threshold-profile figure =========================================
        # One panel per clinical outcome: gain in R2 from splitting the module
        # against the split point, with IP_y marked.
        if (!is.null(curve_df) && nrow(curve_df)) {
            cd <- curve_df[curve_df$outcome %in% PHH_OUTCOMES, , drop = FALSE]
            cd <- cd[is.finite(cd$gain) & is.finite(cd$quantile), ,
                     drop = FALSE]
            if (nrow(cd)) {
                cd$outcome <- factor(cd$outcome,
                                     levels = intersect(MV_OUTCOME_ORDER,
                                                        unique(cd$outcome)))
                grid_pts <- cd[!cd$is_IP, , drop = FALSE]
                ip_pts   <- cd[cd$is_IP, , drop = FALSE]
                # Shaded band from the same flag p_far counted.
                excl <- grid_pts[!grid_pts$in_far_null, , drop = FALSE]
                shade <- if (nrow(excl))
                    data.frame(outcome = factor(levels(cd$outcome),
                                                levels = levels(cd$outcome)),
                               xmin = min(excl$quantile),
                               xmax = max(excl$quantile),
                               stringsAsFactors = FALSE) else NULL
                best <- do.call(rbind, lapply(split(grid_pts, grid_pts$outcome),
                    function(g) if (!nrow(g)) NULL
                                else g[which.max(g$gain), , drop = FALSE]))
                p_tp <- ggplot(grid_pts, aes(x = quantile, y = gain))
                if (!is.null(shade))
                    p_tp <- p_tp + geom_rect(data = shade, inherit.aes = FALSE,
                        aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
                        fill = "grey85", alpha = 0.55)
                p_tp <- p_tp +
                    geom_hline(yintercept = 0, linewidth = 0.3,
                               colour = "grey60") +
                    geom_line(linewidth = 0.6, colour = "#4D4D4D")
                # Marker layers are added only when they have rows (data = NULL
                # would inherit the parent data).
                if (!is.null(best) && nrow(best))
                    p_tp <- p_tp + geom_point(data = best, size = 2.6,
                                              shape = 21, fill = "white",
                                              colour = "#4D4D4D", stroke = 0.7)
                if (nrow(ip_pts))
                    p_tp <- p_tp + geom_point(data = ip_pts, size = 3.4,
                                              shape = 23, fill = "#B2182B",
                                              colour = "black", stroke = 0.6)
                p_tp <- p_tp +
                    facet_wrap(~ outcome, scales = "free_y") +
                    labs(
                        title = sprintf(paste0("Where the module is split:",
                                               " gain in R2 across candidate",
                                               " cDFG thresholds - %s"), ds),
                        subtitle = paste0(
                            "X = split point as a quantile of the module",
                            " cells' cDFG. Y = R2 of `outcome ~ cov4 +",
                            " pre + post` minus `outcome ~ cov4 + module",
                            " total`.\nRed diamond = IP_y. Open circle = best",
                            " grid point. Shaded band = the neighbourhood of",
                            " IP_y excluded from p_far, since a split a hair",
                            " away is the same split.\nA broad plateau",
                            " containing IP_y and a flat line with edge spikes",
                            " are both 'not significant' and mean opposite",
                            " things."),
                        x = "split point (quantile of module cDFG)",
                        y = expression("gain in " * R^2 * " from splitting"),
                        caption = fig_note) +
                    theme_bw(base_size = 12) +
                    theme(
                        plot.title    = element_text(face = "bold", size = 13,
                                                     hjust = 0.5),
                        plot.subtitle = element_text(size = 8.5,
                                                     colour = "grey40",
                                                     hjust = 0.5),
                        plot.caption  = element_text(size = 7,
                                                     colour = "grey55"),
                        strip.background = element_rect(fill = "grey92",
                                                        colour = NA),
                        panel.grid.minor = element_blank(),
                        plot.background = element_rect(fill = "white",
                                                       colour = NA))
                save_pp(p_tp, sprintf("%s_%s_ThresholdProfile", FILE_PREFIX,
                                      ds), w = 13, h = 8)
            }
        }

        # Delta-R2 figure: one figure per covariate set, four bars per outcome
        # in a fixed arm order and a fixed clinical order; stars are the BH q of
        # the arm term.
        if (!is.null(mv_df)) {
            dd <- mv_df[mv_df$tier == "all" & is.finite(mv_df$delta_R2), ,
                        drop = FALSE]
            dd <- dd[dd$arm %in% MV_ARM_ORDER, , drop = FALSE]
            if (nrow(dd) > 0) {
                dd$arm_f <- factor(dd$arm, levels = MV_ARM_ORDER)
                # Outcomes not in the declared order go at the end.
                .lv <- c(MV_OUTCOME_ORDER,
                         setdiff(sort(unique(dd$outcome)), MV_OUTCOME_ORDER))
                .lv <- .lv[.lv %in% dd$outcome]
                dd$outcome_f <- factor(dd$outcome, levels = .lv)
                dd$star <- .sig_star(dd$p_adj)
                # Figures for MV_FIGURE_SETS only (all sets stay in the workbook),
                # one y scale across them, and no figure for a set whose arm VIF
                # is at or above MV_VIF_FIGURE_MAX.
                .arm_vif_of <- function(cid) {
                    v <- dd$arm_VIF[dd$covariate_set == cid]
                    v <- v[is.finite(v)]
                    if (!length(v)) NA_real_ else max(v)
                }
                .fig_sets <- intersect(MV_FIGURE_SETS, unique(dd$covariate_set))
                .fig_drop <- .fig_sets[vapply(.fig_sets, function(cid) {
                    va <- .arm_vif_of(cid)
                    isTRUE(is.finite(va) && va >= MV_VIF_FIGURE_MAX)
                }, logical(1))]
                if (length(.fig_drop)) {
                    for (cid in .fig_drop)
                        cat(sprintf(paste0("    DeltaR2 figure for %s suppressed",
                                           ": arm VIF %.2f >= %.1f.",
                                           " Its rows remain in the workbook.\n"),
                                    cid, .arm_vif_of(cid), MV_VIF_FIGURE_MAX))
                    .fig_sets <- setdiff(.fig_sets, .fig_drop)
                }
                y_max <- suppressWarnings(
                    max(dd$delta_R2[dd$covariate_set %in% .fig_sets],
                        na.rm = TRUE)) * 1.45
                if (!is.finite(y_max) || y_max <= 0) y_max <- 0.1
                for (cid in .fig_sets) {
                    di <- dd[dd$covariate_set == cid, , drop = FALSE]
                    if (nrow(di) == 0L) next
                    # Clinical and control panels; shared y scale.
                    di$panel_f <- factor(
                        ifelse(as.character(di$outcome) %in% MV_CONTROL_OUTCOMES,
                               "control", "clinical"),
                        levels = c("clinical", "control"),
                        labels = unname(MV_PANEL_LABEL[c("clinical", "control")]))
                    p_dr <- ggplot(di, aes(x = outcome_f, y = delta_R2,
                                           fill = arm_f, colour = arm_f)) +
                        geom_col(position = position_dodge(0.7), width = 0.6,
                                 linewidth = 0.5) +
                        scale_fill_manual(values = MV_ARM_FILL, name = NULL,
                                          drop = FALSE) +
                        scale_colour_manual(values = MV_ARM_BORDER, name = NULL,
                                            guide = "none", drop = FALSE) +
                        geom_text(aes(label = star,
                                      y = delta_R2 + y_max * 0.025),
                                  position = position_dodge(0.7),
                                  size = 3.2, colour = "grey15", vjust = 0,
                                  fontface = "bold") +
                        # Shared limit across the figures.
                        scale_y_continuous(
                            limits = c(0, y_max), expand = c(0, 0),
                            labels = function(x) sprintf("%.2f", x)) +
                        geom_hline(yintercept = 0, linewidth = 0.3,
                                   colour = "grey40") +
                        facet_grid(~ panel_f, scales = "free_x",
                                   space = "free_x") +
                        labs(
                            title = sprintf(
                                "Unique variance the arm proportion adds over %s \u2014 %s",
                                MV_COV_TITLE[[cid]], ds),
                            # The Delta is written as plotmath on the y axis (U+0394
                            # is dropped by the Windows graphics device).
                            subtitle = paste0(
                                "outcome ~ %arm + ", MV_COV_TITLE[[cid]],
                                " vs covariate-only  |  bars are the R\u00b2",
                                " increment;",
                                " * BH q<0.05  ** q<0.01  *** q<0.001\n",
                                sprintf(paste0("y axis is COMMON to the %s ",
                                               "figures, so heights compare ",
                                               "directly across them."),
                                        paste(.fig_sets, collapse = " / ")),
                                "\n", MV_COV_NOTE[[cid]]),
                            caption = fig_note, x = NULL,
                            y = expression(Delta * R^2 ~
                                "(full model " * minus * " covariate-only model)")) +
                        theme_bw(base_size = 11) +
                        theme(
                            plot.title    = element_text(face = "bold",
                                                         size = 12, hjust = 0.5),
                            plot.subtitle = element_text(size = 8.5,
                                                         colour = "grey40",
                                                         hjust = 0.5),
                            plot.caption  = element_text(size = 7,
                                                         colour = "grey55"),
                            axis.text.x   = element_text(angle = 35, hjust = 1,
                                                         size = 9),
                            panel.grid.minor   = element_blank(),
                            panel.grid.major.x = element_blank(),
                            legend.position = "top")
                    save_pp(p_dr, sprintf("%s_%s_DeltaR2_%s", FILE_PREFIX,
                                          ds, cid),
                            w = 13, h = 7)
                }
            }
        }

        list(dataset = ds, res_df = res_df, mv_df = mv_df, st_df = st_df,
             ipsplit_df = ipsplit_df, duel_df = duel_df,
             thresh_df = thresh_df, curve_df = curve_df,
             ident_ok = ident_ok,
             scan_df = scan_df,
             jt_df = jt_df, jt_donor = jt_donor,
             mito_cov_ok = mito_cov_ok, depth_note = depth_note,
             extra_ok = extra_ok,
             pred_avail = pred_avail, asrc_excl = asrc_excl,
             n_y_fall_excl = n_y_fall_excl, n_unlab_excl = n_unlab_excl,
             prop_df = prop_df,
             grp_meta = grp_meta, ip_y = ip_y, ip_x = ip_x, u_ip = u_ip,
             ipy_field = ipy$field, arm_src = asrc, n_y_fall = n_y_fall,
             n_unlab = n_unlab, arm_agree = arm_agree, arm_n_cmp = arm_n_cmp,
             arm_agree_u = arm_agree_u,
             has_cdfg = any(is.finite(rc$cDFG_expr)) && is.finite(ip_y),
             n_obs_total = n_obs_total, n_obs_pri = n_obs_pri,
             n_no_donor_obs = n_no_donor_obs,
             n_no_donor = n_no_donor, n_donor_totals = nrow(donor_totals),
             n_cells_rc = nrow(rc))

        }, error = function(e) {
            cat(sprintf("  Dataset failed (%s): %s\n", ds, conditionMessage(e)))
            NULL
        })
        all_res[ds] <- list(res)
    }

    # =========================================================================
    # STEP 3: Workbook
    # =========================================================================
    cat("\nSTEP 3 -- workbook\n")
    ok_ds <- names(all_res)[!vapply(all_res, is.null, logical(1))]
    if (length(ok_ds) == 0L)
        stop("No dataset completed for arm ", arm_tag, ".", call. = FALSE)

    sheets <- list()
    add_tab <- function(nm, df) {
        if (is.null(df) || nrow(df) == 0L) return(invisible(NULL))
        # Sheet names are capped at 31 characters; collisions are numbered.
        key <- substr(nm, 1, 31)
        if (!is.null(sheets[[key]])) {
            base_key <- substr(key, 1, 28)
            i <- 2L
            while (!is.null(sheets[[sprintf("%s_%d", base_key, i)]]) && i < 99L)
                i <- i + 1L
            key <- sprintf("%s_%d", base_key, i)
            cat(sprintf(paste0("  NOTE: sheet name \"%s\" collided after the ",
                               "31-character Excel cap; written as \"%s\"\n"),
                        nm, key))
        }
        sheets[[key]] <<- as.data.frame(df)
    }

    # Provenance sheet: classification axis, IP_y, source file and realised n.
    prov <- do.call(rbind, lapply(ok_ds, function(ds) {
        r <- all_res[[ds]]
        data.frame(
            dataset = ds,
            class_axis = "cDFG at IP_y",
            class_rule = "Pre-IP when cDFG_expr < IP_y",
            IP_y = if (is.finite(r$ip_y)) round(r$ip_y, 6) else NA_real_,
            IP_y_field = r$ipy_field,
            IP_y_source = if (is.na(.p10h_iprs$path)) NA_character_
                          else basename(.p10h_iprs$path),
            IP_x_reported = if (is.finite(r$ip_x)) round(r$ip_x, 6) else NA_real_,
            u_IP_reference_only = if (is.finite(r$u_ip)) round(r$u_ip, 6)
                                  else NA_real_,
            arm_from_group_refined = as.integer(r$arm_src[["group_refined"]]),
            arm_from_IP_status     = as.integer(r$arm_src[["IP_status"]]),
            arm_from_cDFG_at_IP_y  = r$n_y_fall,
            arm_unlabelled         = r$n_unlab,
            # Non-module provenance, kept apart from the module's.
            excl_from_group_refined = as.integer(r$asrc_excl[["group_refined"]]),
            excl_from_IP_status     = as.integer(r$asrc_excl[["IP_status"]]),
            excl_from_cDFG_at_IP_y  = r$n_y_fall_excl,
            excl_unlabelled         = r$n_unlab_excl,
            predictors_available    = paste(r$pred_avail, collapse = ", "),
            inherited_vs_cDFG_at_IP_y = if (is.na(r$arm_agree)) NA_real_
                                        else round(r$arm_agree, 6),
            GOI_split_at_u_IP_would_agree = if (is.na(r$arm_agree_u)) NA_real_
                                            else round(r$arm_agree_u, 6),
            n_cells_refined = r$n_cells_rc,
            n_cells_no_donor = r$n_no_donor,
            n_donors_with_supertype = r$n_donor_totals,
            supertype = PRIMARY_STYPE,
            obs_dataset_cells_total = r$n_obs_total,
            obs_dataset_cells_supertype = r$n_obs_pri,
            arm_tag = ARM_TAG, norm_tag = NORM_TAG,
            stringsAsFactors = FALSE)
    }))
    add_tab("Provenance", prov)

    add_tab("Spearman", do.call(rbind, lapply(ok_ds,
                                              function(d) all_res[[d]]$res_df)))
    add_tab("Steiger", do.call(rbind, lapply(ok_ds,
                                             function(d) all_res[[d]]$st_df)))
    # Restriction tests (IP split and CMRF boundary), told apart by `test`.
    add_tab("Restriction_Tests", do.call(rbind, lapply(ok_ds,
            function(d) all_res[[d]]$ipsplit_df)))
    add_tab("IP_Threshold_Profile", do.call(rbind, lapply(ok_ds,
            function(d) all_res[[d]]$thresh_df)))
    # The curve behind the threshold profile.
    add_tab("IP_Threshold_Curve", do.call(rbind, lapply(ok_ds,
            function(d) all_res[[d]]$curve_df)))
    add_tab("Predictor_Duel", do.call(rbind, lapply(ok_ds,
            function(d) all_res[[d]]$duel_df)))
    add_tab("Jonckheere_Terpstra",
            do.call(rbind, lapply(ok_ds, function(d) all_res[[d]]$jt_df)))
    add_tab("JT_donor_tertiles", do.call(rbind, lapply(ok_ds, function(d) {
        jd <- all_res[[d]]$jt_donor
        if (is.null(jd) || length(jd) == 0L) return(NULL)
        do.call(rbind, lapply(names(jd), function(od)
            jd[[od]][, intersect(c("donor_id", "ordering", "n_pre", "n_post",
                                   "n_cells", "ip_ratio", "jt_group"),
                                 names(jd[[od]])), drop = FALSE]))
    })))
    mv_all <- do.call(rbind, lapply(ok_ds, function(d) all_res[[d]]$mv_df))
    add_tab("Multivariate", mv_all)
    # One row per outcome x arm x candidate covariate.
    add_tab("LeaveOneIn_Scan",
            do.call(rbind, lapply(ok_ds, function(d) all_res[[d]]$scan_df)))
    # Which extra donor variables resolved on this arm, and why any did not.
    add_tab("Donor_Extra_Audit", do.call(rbind, lapply(ok_ds, function(d) {
        r <- all_res[[d]]
        if (is.null(r$depth_note)) return(NULL)
        data.frame(dataset = d,
                   n_usable = length(r$extra_ok),
                   usable = paste(r$extra_ok, collapse = ", "),
                   detail = r$depth_note, stringsAsFactors = FALSE)
    })))

    # cov4 and cov7 side by side.
    if (!is.null(mv_all)) {
        base_mv <- mv_all[mv_all$tier == "all", , drop = FALSE]
        if (nrow(base_mv) > 0) {
            # Keyed on dataset, outcome and predictor.
            k <- paste(base_mv$dataset, base_mv$outcome, base_mv$predictor,
                       sep = "|")
            cmp <- do.call(rbind, lapply(unique(k), function(kk) {
                i4 <- k == kk & base_mv$covariate_set == "cov4"
                i7 <- k == kk & base_mv$covariate_set == "cov7"
                if (!any(i4) && !any(i7)) return(NULL)
                g <- function(i, cn) if (any(i)) base_mv[[cn]][which(i)[1]]
                                     else NA
                data.frame(
                    dataset   = g(i4 | i7, "dataset"),
                    outcome   = g(i4 | i7, "outcome"),
                    predictor = g(i4 | i7, "predictor"),
                    arm       = g(i4 | i7, "arm"),
                    n_cov4 = g(i4, "n_donors"), beta_cov4 = g(i4, "beta_std"),
                    p_cov4 = g(i4, "p_value"), padj_cov4 = g(i4, "p_adj"),
                    dR2_cov4 = g(i4, "delta_R2"), maxVIF_cov4 = g(i4, "max_VIF"),
                    n_cov7 = g(i7, "n_donors"), beta_cov7 = g(i7, "beta_std"),
                    p_cov7 = g(i7, "p_value"), padj_cov7 = g(i7, "p_adj"),
                    dR2_cov7 = g(i7, "delta_R2"), maxVIF_cov7 = g(i7, "max_VIF"),
                    stringsAsFactors = FALSE)
            }))
            if (!is.null(cmp)) {
                cmp$beta_attenuation <- ifelse(
                    is.finite(cmp$beta_cov4) & is.finite(cmp$beta_cov7) &
                        abs(cmp$beta_cov4) > 1e-12,
                    round(1 - cmp$beta_cov7 / cmp$beta_cov4, 4), NA_real_)
                add_tab("Cov4_vs_Cov7", cmp)
            }
        }
    }

    # Ordinal audit: every level, its mapped value and the naive digit parse.
    add_tab("Ordinal_Audit", ordinal_audit_df)

    cov_audit <- data.frame(
        covariate = names(cov_source), source = unname(cov_source),
        stringsAsFactors = FALSE)
    cov_audit$n_complete <- vapply(cov_audit$covariate, function(cn)
        if (cn %in% names(donor_clinical))
            sum(!is.na(donor_clinical[[cn]])) else NA_integer_, integer(1))
    add_tab("Covariate_Audit", cov_audit)

    donor_out <- donor_clinical[, intersect(
        c("donor_id", "Thal_num", "CERAD_num", "CAA_num", "ADNC_rank",
          "Braak_num", "APOE_e4_dose", "CogStatus_num", "MMSE_num", "Age_num",
          "Age_censored", "Sex_num", "PMI_num", "RIN_num", "BrainPH_num",
          "BrainWeight_num", "Lewy_binary", "Severe_flag"),
        names(donor_clinical)), drop = FALSE]
    add_tab("Donor_Table", donor_out)

    for (ds in ok_ds) add_tab(paste0("Donor_pct_", ds),
                              all_res[[ds]]$prop_df[, intersect(
                                  c("donor_id", "n_total", names(PRED_GROUPS)),
                                  names(all_res[[ds]]$prop_df)), drop = FALSE])

    # Heatmap pivot, one sheet per dataset: rows = clinical variable,
    # columns = arm group, values = rho.
    for (ds in ok_ds) {
        rd <- all_res[[ds]]$res_df
        pv <- tryCatch(tidyr::pivot_wider(
            rd[, c("variable", "group_label", "rho")],
            names_from = "group_label", values_from = "rho"),
            error = function(e) NULL)
        if (!is.null(pv)) add_tab(paste0("Heatmap_rho_", ds), as.data.frame(pv))
    }

    xlsx_path <- file.path(OUT_DIR, paste0(FILE_PREFIX, "_Prognosis.xlsx"))
    tryCatch({
        writexl::write_xlsx(sheets, xlsx_path)
        cat(sprintf("  %d sheet(s) written to %s\n", length(sheets),
                    basename(xlsx_path)))
    }, error = function(e)
        cat(sprintf("  Workbook not written: %s\n", conditionMessage(e))))

    # =========================================================================
    # STEP 4: Guardrails
    # =========================================================================
    cat("\nSTEP 4 -- guardrails\n")

    .num1 <- function(z, d = 6L) if (length(z) != 1L || !is.finite(z)) "NA"
                                 else sprintf(paste0("%.", d, "f"), z)
    .ok1  <- function(z) length(z) == 1L && is.finite(z)

    # GR-10H-01: the IPRS RDS and both inflection-point fields.
    .have_iprs <- !is.na(.p10h_iprs$path)
    gr10h("GR-10H-01", if (.have_iprs) "PASS" else "FAIL",
          if (.have_iprs) sprintf("IPRS RDS: %s", basename(.p10h_iprs$path))
          else sprintf(paste0("no IPRS RDS found; looked for %s."),
                       .p10h_iprs$note))

    for (ds in ok_ds) {
        r <- all_res[[ds]]

        # GR-10H-02: IP_y present, and IP_x reported beside it.
        gr10h(paste0("GR-10H-02-", ds), if (.ok1(r$ip_y)) "PASS" else "FAIL",
              sprintf(paste0("IP_y = %s (field: %s) | IP_x reported = %s | ",
                             "u_IP reference only = %s"),
                      .num1(r$ip_y),
                      if (length(r$ipy_field) == 1L && !is.na(r$ipy_field))
                          r$ipy_field else "NONE",
                      .num1(r$ip_x), .num1(r$u_ip)))

        # GR-10H-05: classification axis; FAIL if the inherited labels agree
        # with the cDFG rule below 0.999.
        .agree_bad <- isTRUE(!is.na(r$arm_agree) && r$arm_agree < 0.999)
        .no_check  <- !isTRUE(r$has_cdfg) || is.na(r$arm_agree)
        .degraded  <- isTRUE(r$n_y_fall > 0L) || isTRUE(r$n_unlab > 0L) ||
                      !isTRUE(as.integer(r$arm_src[["group_refined"]]) > 0)
        gr10h(paste0("GR-10H-05-", ds),
              if (.agree_bad) "FAIL"
              else if (.no_check || .degraded) "WARN" else "PASS",
              sprintf(paste0("axis = cDFG at IP_y | arm from: group_refined %d,",
                             " IP_status %d, cDFG@IP_y %d, unlabelled %d | ",
                             "inherited vs `cDFG < IP_y`: %s on %d cell(s) | a",
                             " GOI split at u_IP would agree %s%s%s"),
                      as.integer(r$arm_src[["group_refined"]]),
                      as.integer(r$arm_src[["IP_status"]]),
                      r$n_y_fall, r$n_unlab,
                      if (is.na(r$arm_agree)) "not checked"
                      else sprintf("%.2f%%", 100 * r$arm_agree),
                      r$arm_n_cmp,
                      if (is.na(r$arm_agree_u)) "NA"
                      else sprintf("%.2f%%", 100 * r$arm_agree_u),
                      if (.agree_bad)
                          " -- the inherited labels do not match cDFG at this IP_y"
                      else "",
                      if (!isTRUE(r$has_cdfg))
                          " -- cDFG_expr or IP_y unavailable, so the cross-check could not run"
                      else if (.degraded)
                          " -- not every module cell took its arm from group_refined"
                      else ""))

        # GR-10H-04: barcode coverage (counts reported).
        .cov <- if (r$n_cells_rc > 0) 1 - r$n_no_donor / r$n_cells_rc else NA_real_
        gr10h(paste0("GR-10H-04-", ds),
              if (!.ok1(.cov)) "WARN" else if (.cov >= 0.999) "PASS" else "FAIL",
              sprintf(paste0("%d of %d refined cell(s) matched a donor in ",
                             "obs_dataset.csv (%s) | obs_dataset holds %d cell",
                             "(s), %d of supertype %s (%d with no donor ID), ",
                             "over %d donor(s)"),
                      r$n_cells_rc - r$n_no_donor, r$n_cells_rc,
                      if (.ok1(.cov)) sprintf("%.2f%%", 100 * .cov) else "NA",
                      r$n_obs_total, r$n_obs_pri, PRIMARY_STYPE,
                      r$n_no_donor_obs, r$n_donor_totals))

        # GR-10H-06: did any group fall under the cell floor.
        sk <- names(r$grp_meta)[vapply(r$grp_meta,
                                       function(z) isTRUE(z$skipped),
                                       logical(1))]
        gr10h(paste0("GR-10H-06-", ds),
              if (length(sk) == 0L) "PASS" else "WARN",
              if (length(sk) == 0L)
                  sprintf("all %d arm group(s) cleared the %d-cell floor",
                          length(ARM_GROUPS), MIN_CELLS_GRP)
              else sprintf("%d group(s) below the %d-cell floor and not tested: %s",
                           length(sk), MIN_CELLS_GRP, paste(sk, collapse = ", ")))

        # GR-10H-08: realised n and collinearity for every covariate set.
        # Verdict on the arm term's VIF: WARN at 5, FAIL at 10.
        if (is.null(r$mv_df)) {
            gr10h(paste0("GR-10H-08-", ds), "WARN",
                  "no multivariate model was fitted")
        } else {
            .set_i <- function(cid) r$mv_df$tier == "all" &
                r$mv_df$covariate_set == cid & is.finite(r$mv_df$n_donors)
            i4 <- .set_i("cov4")
            i7 <- .set_i("cov7")
            # n range over models that cleared the floor; the rest are counted.
            .fit <- function(i) i & is.finite(r$mv_df$n_donors) &
                                r$mv_df$n_donors >= MIN_MV_DONORS
            .n_rng <- function(i) {
                f <- .fit(i)
                if (!any(f)) return("none fitted")
                rg <- range(r$mv_df$n_donors[f])
                paste(if (rg[1] == rg[2]) rg[1] else paste(rg, collapse = "-"),
                      sprintf("(%d of %d model(s) fitted)", sum(f), sum(i)))
            }
            .vmax <- function(i) {
                v <- r$mv_df$max_VIF[.fit(i)]
                v <- v[is.finite(v)]
                if (!length(v)) "NA" else sprintf("%.2f", max(v))
            }
            # The arm term's own VIF.
            .varm <- function(i) {
                v <- r$mv_df$arm_VIF[.fit(i)]
                v <- v[is.finite(v)]
                if (!length(v)) NA_real_ else max(v)
            }
            # Declared sets absent from mv_df are named.
            .declared <- c("cov4", "cov7", "covM", "covT", "covP", "covE")
            .sets <- intersect(.declared, unique(r$mv_df$covariate_set))
            .absent <- setdiff(.declared, .sets)
            .rows <- character(0); .worst <- 0; .worst_set <- NA_character_
            .none <- character(0)
            for (cid in .sets) {
                ii <- .set_i(cid)
                if (!any(.fit(ii))) { .none <- c(.none, cid); next }
                v  <- suppressWarnings(max(r$mv_df$max_VIF[.fit(ii)],
                                           na.rm = TRUE))
                va <- .varm(ii)
                # Verdict follows the arm VIF; the maximum is printed too.
                if (is.finite(va) && va > .worst) { .worst <- va; .worst_set <- cid }
                .rows <- c(.rows, sprintf("%s n=%s armVIF=%s maxVIF=%s%s", cid,
                                          .n_rng(ii),
                                          if (is.finite(va)) sprintf("%.2f", va)
                                          else "NA", .vmax(ii),
                                          if (is.finite(va) && va >= 10) " FAIL"
                                          else if (is.finite(va) && va >= 5) " WARN"
                                          else if (is.finite(v) && v >= 10)
                                              " (max high, arm clean)" else ""))
            }
            gr10h(paste0("GR-10H-08-", ds),
                  if (length(.none) || length(.absent) || !length(.rows)) "WARN"
                  else if (.worst >= 10) "FAIL"
                  else if (.worst >= 5) "WARN" else "PASS",
                  sprintf(paste0("realised n and collinearity per covariate ",
                                 "set: %s | worst ARM VIF %s in %s%s%s%s"),
                          paste(.rows, collapse = " ; "),
                          if (.worst > 0) sprintf("%.2f", .worst) else "NA",
                          if (is.na(.worst_set)) "none" else .worst_set,
                          if (length(.none))
                              sprintf(" | no model fitted for: %s",
                                      paste(.none, collapse = ", ")) else "",
                          if (length(.absent))
                              sprintf(paste0(" | not run (a covariate was",
                                             " absent or constant): %s"),
                                      paste(.absent, collapse = ", ")) else "",
                          if (.worst >= 10)
                              " -- arm VIF >= 10: that set's figure is suppressed"
                          else if (.worst >= 5)
                              " -- arm VIF >= 5"
                          else ""))
        }
        # GR-10H-10: the non-module comparison set.
        .have_excl <- c("pct_excl_post", "pct_excl_pre", "pct_excl") %in%
                      r$pred_avail
        gr10h(paste0("GR-10H-10-", ds),
              if (all(.have_excl)) "PASS"
              else if (any(.have_excl)) "WARN" else "FAIL",
              sprintf(paste0("excluded comparison set: %d of 3 predictor(s) ",
                             "available (%s) | arms from cDFG@IP_y %d, ",
                             "unlabelled %d | %s"),
                      sum(.have_excl),
                      paste(c("Non-module Post-IP", "Non-module Pre-IP",
                              "Non-module (all)")[.have_excl], collapse = ", "),
                      r$n_y_fall_excl, r$n_unlab_excl,
                      if (all(.have_excl))
                          "excluded cells take their arm from the cDFG rule at IP_y"
                      else "module-vs-excluded contrast incomplete; Steiger pairs needing a missing predictor skipped"))

        # GR-10H-12: trend test, orderings and donors.
        .jt <- r$jt_df
        .njt_ord <- if (is.null(.jt)) 0L else length(unique(.jt$ordering))
        .jt_bad  <- if (is.null(.jt)) 0L else sum(!is.finite(.jt$JT_Z))
        gr10h(paste0("GR-10H-12-", ds),
              if (.njt_ord >= 2L && .jt_bad == 0L) "PASS"
              else if (.njt_ord >= 1L) "WARN" else "FAIL",
              sprintf(paste0("%d of 2 donor ordering(s) tested | %s row(s), %d ",
                             "not testable | tertile sizes %s | %d donor(s) ",
                             "dropped below the %d-cell evidence floor | ",
                             "variance tie-corrected (Lehmann 1975)%s"),
                      .njt_ord,
                      if (is.null(.jt)) "0" else as.character(nrow(.jt)),
                      .jt_bad,
                      # Tertile sizes from rows that produced a Z.
                      if (is.null(.jt) || !any(is.finite(.jt$JT_Z))) "none"
                      else paste(unique(.jt$group_sizes[is.finite(.jt$JT_Z)]),
                                 collapse = " ; "),
                      if (is.null(.jt) || nrow(.jt) == 0L) 0L
                      else max(.jt$n_donors_dropped_cell_floor, na.rm = TRUE),
                      if (is.null(.jt) || nrow(.jt) == 0L) JT_MIN_CELLS_PER_DONOR
                      else max(.jt$min_cells_per_donor, na.rm = TRUE),
                      if (.njt_ord < 2L)
                          " -- the non-module (control) ordering was not tested"
                      else ""))

        # GR-10H-11: Steiger contrasts. Rows that could not be computed (noted
        # rows) give WARN.
        .nst <- if (is.null(r$st_df)) 0L else nrow(r$st_df)
        .npair <- if (is.null(r$st_df)) 0L
                  else length(unique(paste(r$st_df$predictor_1,
                                           r$st_df$predictor_2)))
        .noted <- if (is.null(r$st_df)) logical(0)
                  else !is.na(r$st_df$note) & nzchar(r$st_df$note)
        .n_noted <- sum(.noted)
        .note_txt <- if (!.n_noted) "" else {
            .nd <- r$st_df[.noted, , drop = FALSE]
            .by <- tapply(seq_len(nrow(.nd)), .nd$note, function(ii)
                sprintf("%s (%d row(s) on %s)", .nd$note[ii[1]], length(ii),
                        paste(unique(paste(.nd$predictor_1[ii], "vs",
                                           .nd$predictor_2[ii])),
                              collapse = "; ")))
            sprintf(paste0(" -- not computed (not nulls): %s"),
                    paste(unlist(.by), collapse = " | "))
        }
        gr10h(paste0("GR-10H-11-", ds),
              if (.npair < 1L) "FAIL"
              else if (.npair < length(STEIGER_PAIRS) || .n_noted > 0L) "WARN"
              else "PASS",
              sprintf(paste0("%d of %d Steiger pair(s) tested over %d row(s) ",
                             "| test: Meng-Rosenthal-Rubin z for dependent ",
                             "overlapping correlations, Spearman variance ",
                             "factor %.2f | %d row(s) carry a note%s"),
                      .npair, length(STEIGER_PAIRS), .nst,
                      SPEARMAN_Z_VAR_FACTOR, .n_noted, .note_txt))
    }

    # GR-10H-03: the clinical table and the ordinal maps.
    n_unmapped_tot <- sum(ordinal_audit_df$mapped_to == "UNMAPPED", na.rm = TRUE)
    n_missing_src  <- sum(ordinal_audit_df$mapped_to %in%
                          c("SOURCE COLUMN NOT FOUND",
                            "COLUMN PRESENT BUT EMPTY"), na.rm = TRUE)
    gr10h("GR-10H-03",
          if (n_unmapped_tot > 0 || n_missing_src > 0) "FAIL" else "PASS",
          sprintf(paste0("%d donor(s) in the clinical table | %d ordinal ",
                         "level(s) audited | %d UNMAPPED | %d variable(s) with ",
                         "no usable source column | %d donor(s) carry conflicting ",
                         "donor-level covariates%s"),
                  nrow(donor_clinical), nrow(ordinal_audit_df),
                  n_unmapped_tot, n_missing_src, length(cov_conflicts),
                  if (n_unmapped_tot > 0)
                      " -- unmapped levels become NA"
                  else ""))

    # GR-10H-07: covariate completeness (cov7 should cost about one donor).
    gr10h("GR-10H-07",
          if (is.na(n_cov4) || is.na(n_cov7)) "WARN"
          else if (n_cov7 >= n_cov4 - 3L) "PASS" else "WARN",
          sprintf(paste0("complete cases: %s on the 4-covariate set, %s on the ",
                         "7-covariate set%s"),
                  if (is.na(n_cov4)) "NA" else as.character(n_cov4),
                  if (is.na(n_cov7)) "NA" else as.character(n_cov7),
                  if (!is.na(n_cov4) && !is.na(n_cov7) && n_cov7 < n_cov4 - 3L)
                      " -- larger than expected"
                  else ""))

    # GR-10H-13: mito/depth-adjusted set (covM); Module Post-IP effect on
    # cognitive outcomes against cov4.
    for (ds in DATASETS) {
        r <- all_res[[ds]]
        if (is.null(r) || is.null(r$mv_df)) next
        mm <- r$mv_df[r$mv_df$tier == "all" &
                      r$mv_df$covariate_set == "covM" &
                      r$mv_df$arm == "Module Post-IP", , drop = FALSE]
        m4 <- r$mv_df[r$mv_df$tier == "all" &
                      r$mv_df$covariate_set == "cov4" &
                      r$mv_df$arm == "Module Post-IP", , drop = FALSE]
        if (nrow(mm) == 0L) {
            gr10h(paste0("GR-10H-13-", ds), "WARN",
                  paste0("mito/depth-adjusted models did not run",
                         if (!is.null(r$depth_note) && nzchar(r$depth_note))
                             paste0(": ", r$depth_note) else ""))
            next
        }
        .cmp <- character(0)
        for (oc in c("Cognitive status", "MMSE score", "ADNC rank")) {
            a <- mm[mm$outcome == oc, , drop = FALSE]
            b <- m4[m4$outcome == oc, , drop = FALSE]
            if (nrow(a) == 0L || nrow(b) == 0L) next
            .cmp <- c(.cmp, sprintf("%s beta %+.2f%s (4-cov %+.2f%s)", oc,
                                    a$beta_std[1], .sig_star(a$p_adj[1]),
                                    b$beta_std[1], .sig_star(b$p_adj[1])))
        }
        .nsig <- sum(mm$significant, na.rm = TRUE)
        gr10h(paste0("GR-10H-13-", ds),
              if (.nsig > 0L) "PASS" else "WARN",
              sprintf(paste0("mito/depth-adjusted set ran on %d model(s), %d",
                             " significant after BH | %s | %s"),
                      nrow(mm), .nsig,
                      if (length(.cmp)) paste(.cmp, collapse = "; ")
                      else "no cognitive outcome fitted",
                      if (!is.null(r$depth_note)) r$depth_note else ""))
    }

    # GR-10H-16: every declared clinical variable reaches the donor frame, the
    # Spearman and the multivariate.
    for (ds in DATASETS) {
        r <- all_res[[ds]]
        if (is.null(r) || is.null(r$prop_df)) next
        want <- vapply(CLINICAL_PARAMS, function(z) z$col, character(1))
        # Not `labs`, which would shadow ggplot2::labs().
        .lab <- vapply(CLINICAL_PARAMS, function(z) z$label, character(1))
        no_col <- .lab[!want %in% names(r$prop_df)]
        no_mv  <- if (is.null(r$mv_df)) .lab
                  else .lab[!.lab %in% unique(r$mv_df$outcome)]
        no_sp  <- if (is.null(r$res_df)) .lab
                  else .lab[!.lab %in% unique(r$res_df$variable)]
        # Outcomes under the donor floor are left to GR-10H-14.
        n_low <- vapply(CLINICAL_PARAMS, function(z)
            if (z$col %in% names(r$prop_df))
                sum(is.finite(suppressWarnings(as.numeric(r$prop_df[[z$col]]))))
            else 0L, numeric(1))
        no_mv <- setdiff(no_mv, .lab[n_low < MIN_MV_DONORS])
        bad <- c(if (length(no_col)) paste0("absent from the donor frame: ",
                                            paste(no_col, collapse = ", ")),
                 if (length(no_sp))  paste0("absent from the Spearman: ",
                                            paste(no_sp, collapse = ", ")),
                 if (length(no_mv))  paste0("absent from the multivariate: ",
                                            paste(no_mv, collapse = ", ")))
        gr10h(paste0("GR-10H-16-", ds),
              if (length(bad) == 0L) "PASS" else "FAIL",
              if (length(bad) == 0L)
                  sprintf(paste0("all %d clinical variable(s) reach the donor",
                                 " frame, the Spearman and the multivariate;",
                                 " OUTCOMES derived from CLINICAL_PARAMS, %d",
                                 " outcome(s) modelled"),
                          length(CLINICAL_PARAMS),
                          if (is.null(r$mv_df)) 0L
                          else length(unique(r$mv_df$outcome)))
              else paste(bad, collapse = " | "))
    }

    # GR-10H-14: outcomes too sparse to model are declared.
    for (ds in DATASETS) {
        r <- all_res[[ds]]
        if (is.null(r) || is.null(r$prop_df)) next
        .sp <- character(0)
        # MoCA completeness is reported (not an outcome).
        if ("MOCA_num" %in% names(r$prop_df)) {
            .nm <- sum(is.finite(suppressWarnings(
                as.numeric(r$prop_df[["MOCA_num"]]))))
            .sp <- c(.sp, sprintf("MOCA score (%d donor(s), not an outcome)",
                                  .nm))
        }
        for (cp in CLINICAL_PARAMS) {
            if (!(cp$col %in% names(r$prop_df))) {
                .sp <- c(.sp, paste0(cp$label, " (column absent)")); next
            }
            .n <- sum(is.finite(suppressWarnings(
                as.numeric(r$prop_df[[cp$col]]))))
            if (.n < MIN_MV_DONORS)
                .sp <- c(.sp, sprintf("%s (%d donor(s))", cp$label, .n))
        }
        gr10h(paste0("GR-10H-14-", ds),
              if (length(.sp) == 0L) "PASS" else "WARN",
              if (length(.sp) == 0L)
                  sprintf("all %d clinical outcome(s) clear the %d-donor floor",
                          length(CLINICAL_PARAMS), MIN_MV_DONORS)
              else sprintf(paste0("below the %d-donor modelling floor and",
                                  " not fitted: %s"),
                           MIN_MV_DONORS, paste(.sp, collapse = "; ")))
    }

    # GR-10H-15: leave-one-in scan verdict.
    for (ds in DATASETS) {
        r <- all_res[[ds]]
        sc <- if (is.null(r)) NULL else r$scan_df
        if (is.null(sc) || nrow(sc) == 0L) {
            gr10h(paste0("GR-10H-15-", ds), "WARN",
                  "leave-one-in scan produced no rows")
            next
        }
        key <- sc[sc$arm == "Module Post-IP" &
                  sc$outcome %in% c("Cognitive status", "MMSE score",
                                    "ADNC rank", "CERAD score"), , drop = FALSE]
        lost <- key[!is.na(key$lost_significance) &
                    key$lost_significance, , drop = FALSE]
        flip     <- sum(key$sign_flipped, na.rm = TRUE)
        flip_any <- sum(key$sign_flipped_any, na.rm = TRUE)
        # Rank on reportable percentages only.
        .rk <- key[is.finite(key$pct_change_beta), , drop = FALSE]
        .worst <- if (nrow(.rk)) .rk[which.max(abs(.rk$pct_change_beta)), ,
                                     drop = FALSE] else NULL
        .n_supp <- sum(key$pct_suppressed, na.rm = TRUE)
        gr10h(paste0("GR-10H-15-", ds),
              if (nrow(lost) == 0L && flip == 0L) "PASS" else "WARN",
              sprintf(paste0("leave-one-in over %d candidate(s) on the four",
                             " headline outcomes: %d of %d model(s) lose",
                             " significance, %d flip sign from a core |beta|",
                             " >= %.2f (%d at any magnitude)%s"),
                      length(unique(sc$covariate_key)), nrow(lost), nrow(key),
                      flip, MV_MIN_BETA_FOR_PCT, flip_any,
                      paste0(
                          if (!is.null(.worst) && nrow(.worst) == 1L &&
                              is.finite(.worst$pct_change_beta))
                              sprintf(paste0(" | largest shift: %s on %s moves",
                                             " beta %+.1f%% (%+.3f absolute)"),
                                      .worst$added_covariate, .worst$outcome,
                                      .worst$pct_change_beta,
                                      .worst$abs_change_beta)
                          else " | no reportable percentage shift",
                          if (.n_supp > 0L)
                              sprintf(paste0(" | %d model(s) had a core beta",
                                             " below %.2f, so their percentage",
                                             " change is not reported"),
                                      .n_supp, MV_MIN_BETA_FOR_PCT)
                          else "")))
    }

    # GR-10H-17: figure arm labels must match PRED_LABEL, and each predictor
    # label must equal its group label.
    .grp_lab <- vapply(ARM_GROUPS, `[[`, character(1), "label")
    .lab_split <- names(PRED_GROUPS)[
        unname(PRED_LABEL) != unname(.grp_lab[PRED_GROUPS])]
    .fig_lab_missing <- setdiff(MV_ARM_ORDER, unname(PRED_LABEL))
    .fig_lab_unstyled <- setdiff(MV_ARM_ORDER,
                                 intersect(names(MV_ARM_FILL),
                                           names(MV_ARM_BORDER)))
    gr10h("GR-10H-17",
          if (length(.fig_lab_missing) || length(.fig_lab_unstyled) ||
              length(.lab_split)) "FAIL"
          else "PASS",
          if (length(.lab_split))
              sprintf(paste0("the same cell set carries two names: %s"),
                      paste(sprintf("%s is '%s' as a predictor but '%s' as a group",
                                    .lab_split, PRED_LABEL[.lab_split],
                                    .grp_lab[PRED_GROUPS[.lab_split]]),
                            collapse = " | "))
          else if (length(.fig_lab_missing) || length(.fig_lab_unstyled))
              sprintf(paste0("figure arm labels do not match the data:",
                             " not in PRED_LABEL: %s | missing a fill or",
                             " border: %s"),
                      paste(c(.fig_lab_missing, "none")[
                          seq_len(max(1L, length(.fig_lab_missing)))],
                          collapse = ", "),
                      paste(c(.fig_lab_unstyled, "none")[
                          seq_len(max(1L, length(.fig_lab_unstyled)))],
                          collapse = ", "))
          else sprintf(paste0("all %d plotted arm label(s) resolve in",
                              " PRED_LABEL and carry a fill and a border;",
                              " every predictor label equals its group label;",
                              " %d predictor(s) defined, %d plotted, %d group(s)"),
                       length(MV_ARM_ORDER), length(PRED_LABEL),
                       length(MV_ARM_ORDER), length(ARM_GROUPS)))

    # GR-10H-18: Original only as a fraction of the SKCS-core, against the
    # nominal 5% miss rate of the CMRF band.
    for (ds in DATASETS) {
        r <- all_res[[ds]]
        gm <- if (is.null(r)) NULL else r$grp_meta
        if (is.null(gm) || is.null(gm$Ori_only) || is.null(gm$OR_PreIP)) {
            gr10h(paste0("GR-10H-18-", ds), "WARN",
                  "Ori-only or Orig+Refined group counts unavailable")
            next
        }
        .n_ori  <- gm$Ori_only$n_cells
        .n_core <- gm$OR_PreIP$n_cells + gm$OR_PostIP$n_cells + .n_ori
        .frac   <- if (.n_core > 0) .n_ori / .n_core else NA_real_
        gr10h(paste0("GR-10H-18-", ds),
              if (!is.finite(.frac)) "WARN"
              else if (abs(.frac - CMRF_MISS_NOMINAL) <= CMRF_MISS_TOL) "PASS"
              else "WARN",
              sprintf(paste0("Ori-only is %d of %d SKCS-core cell(s) = %.2f%%;",
                             " the CMRF band carries %.0f%% coverage so the",
                             " expected residual is %.2f%% (tolerance +/-%.2f",
                             " points)%s"),
                      .n_ori, .n_core, 100 * .frac,
                      100 * (1 - CMRF_MISS_NOMINAL), 100 * CMRF_MISS_NOMINAL,
                      100 * CMRF_MISS_TOL,
                      if (is.finite(.frac) &&
                          abs(.frac - CMRF_MISS_NOMINAL) > CMRF_MISS_TOL)
                          " -- outside tolerance"
                      else ""))
    }

    # == GR-10H-19 and GR-10H-20 ==============================================
    # These report a result; a null result is WARN, not FAIL.
    for (ds in ok_ds) {
        r <- all_res[[ds]]
        # Filter to this test before the emptiness check.
        ips <- r$ipsplit_df
        if (!is.null(ips)) ips <- ips[ips$test == "ip_split", , drop = FALSE]
        if (is.null(ips) || !nrow(ips)) {
            gr10h(paste0("GR-10H-19-", ds), "WARN",
                  paste0("the IP-split test did not run (needs pct_module_all,",
                         " pct_pre and pct_post)"))
        } else {
            i4 <- ips$covariate_set == "cov4"
            n4 <- sum(i4); s4 <- sum(ips$significant[i4], na.rm = TRUE)
            .top <- ips[i4 & is.finite(ips$p_adj), , drop = FALSE]
            .top <- .top[order(.top$p_adj), , drop = FALSE]
            gr10h(paste0("GR-10H-19-", ds), if (s4 > 0L) "PASS" else "WARN",
                  sprintf(paste0("IP-split test (nested F, 1 df, cov4): the",
                                 " split beats the unsplit module on %d of %d",
                                 " outcome(s) after BH. Median gain in R2 from",
                                 " splitting %.4f.%s%s"),
                          s4, n4,
                          stats::median(ips$gain_from_split[i4], na.rm = TRUE),
                          if (nrow(.top))
                              sprintf(paste0(" Strongest: %s F=%.2f BH=%.3g",
                                             " (split dR2 %.3f vs unsplit %.3f)."),
                                      .top$outcome[1], .top$F_value[1],
                                      .top$p_adj[1], .top$delta_R2_split[1],
                                      .top$delta_R2_unsplit[1]) else "",
                          if (s4 > 0L) ""
                          else " No outcome needs the split on this arm."))
        }
        dl <- r$duel_df
        if (is.null(dl) || !nrow(dl)) {
            gr10h(paste0("GR-10H-20-", ds), "WARN",
                  paste0("the predictor duels did not run; no head-to-head",
                         " is available"))
        } else {
            .lines <- character(0); .lost <- character(0)
            for (dn in unique(dl$duel)) {
                i <- dl$duel == dn
                # "both" is counted separately.
                nA <- sum(dl$favours[i] == "A"); nB <- sum(dl$favours[i] == "B")
                nBo <- sum(dl$favours[i] == "both")
                .lines <- c(.lines,
                            sprintf(paste0("%s: module only %d, rival only %d,",
                                           " both %d, neither %d of %d"),
                                    dn, nA, nB, nBo,
                                    sum(!(dl$favours[i] %in%
                                          c("A", "B", "both"))), sum(i)))
                if (nB > nA) .lost <- c(.lost, dn)
            }
            # Joint VIF >= 5 on the module term (Inf counts, NA does not).
            .hiv <- sum(!is.na(dl$joint_VIF_A) & dl$joint_VIF_A >= 5)
            gr10h(paste0("GR-10H-20-", ds),
                  if (length(.lost) || .hiv > 0L) "WARN" else "PASS",
                  sprintf(paste0("predictor duels over the cov4 core",
                                 " (verdict: encompassing test, BH within a duel;",
                                 " Vuong and a %d-resample donor bootstrap",
                                 " reported). %s%s%s"),
                          PHH_BOOT, paste(.lines, collapse = " | "),
                          if (length(.lost))
                              sprintf(paste0(" | rival wins more outcomes",
                                             " than the module on: %s"),
                                      paste(.lost, collapse = ", ")) else "",
                          if (.hiv > 0L)
                              sprintf(paste0(" | %d row(s) carry a joint VIF >= 5",
                                             " on the module term"), .hiv)
                          else ""))
        }
    }

    # == GR-10H-21 and GR-10H-22 ==============================================
    for (ds in ok_ds) {
        r <- all_res[[ds]]
        # GR-10H-21: the CMRF boundary; the sum identity is reported first.
        dev <- if (is.null(r$ident_ok[["cmrf_split"]])) NA_real_
               else r$ident_ok[["cmrf_split"]]
        cm <- if (is.null(r$ipsplit_df)) NULL
              else r$ipsplit_df[r$ipsplit_df$test == "cmrf_split", ,
                                drop = FALSE]
        if (is.null(cm) || !nrow(cm)) {
            gr10h(paste0("GR-10H-21-", ds), "WARN",
                  sprintf(paste0("the CMRF-boundary test did not run.",
                                 " Sum identity max |dev|: %s"),
                          if (is.finite(dev)) sprintf("%.3e", dev)
                          else "not checked"))
        } else {
            i4 <- cm$covariate_set == "cov4"
            n4 <- sum(i4); s4 <- sum(cm$significant[i4], na.rm = TRUE)
            top <- cm[i4 & is.finite(cm$p_adj), , drop = FALSE]
            top <- top[order(top$p_adj), , drop = FALSE]
            gr10h(paste0("GR-10H-21-", ds), if (s4 > 0L) "PASS" else "WARN",
                  sprintf(paste0("CMRF-boundary test (nested F, 1 df, cov4):",
                                 " the module/non-module split beats one",
                                 " undivided Post-IP percentage on %d of %d",
                                 " outcome(s) after BH. Median gain in R2",
                                 " %.4f. Sum identity `pct_any_post ==",
                                 " pct_post + pct_excl_post` max |dev|",
                                 " %.3e.%s%s"),
                          s4, n4,
                          stats::median(cm$gain_from_split[i4], na.rm = TRUE),
                          dev,
                          if (nrow(top))
                              sprintf(paste0(" Strongest: %s F=%.2f BH=%.3g",
                                             " (split dR2 %.3f vs undivided",
                                             " %.3f)."),
                                      top$outcome[1], top$F_value[1],
                                      top$p_adj[1], top$delta_R2_split[1],
                                      top$delta_R2_unsplit[1]) else "",
                          if (s4 > 0L) ""
                          else " No outcome needs the CMRF boundary on this arm."))
        }
        # GR-10H-22: is IP_y a special threshold, or would any split do?
        th <- r$thresh_df
        if (is.null(th) || !nrow(th)) {
            gr10h(paste0("GR-10H-22-", ds), "WARN",
                  "the placebo-threshold profile did not run")
        } else {
            nsig <- sum(is.finite(th$p_far_adj) & th$p_far_adj < 0.05)
            recon <- max(th$recon_dev_at_IP, na.rm = TRUE)
            near <- mean(abs(th$best_quantile - th$IP_quantile), na.rm = TRUE)
            # A recomputation that does not reproduce the split FAILs.
            gr10h(paste0("GR-10H-22-", ds),
                  if (!is.finite(recon) || recon > 1e-6) "FAIL"
                  else if (nsig > 0L) "PASS" else "WARN",
                  sprintf(paste0("threshold profile over %d grid point(s) in",
                                 " quantiles of module cDFG, IP_y at quantile",
                                 " %.3f: IP_y beats a genuinely different",
                                 " threshold on %d of %d outcome(s) (p_far,",
                                 " BH, excluding +/-%.2f quantile around",
                                 " IP_y). Median p_far %.3f; the best grid",
                                 " threshold sits %.3f quantile point(s) from",
                                 " IP_y on average. Recomputed pct_post at",
                                 " IP_y vs the stored column: max |dev| %.2e",
                                 " (tolerance 1.0e-06).%s"),
                          max(th$n_grid, na.rm = TRUE), th$IP_quantile[1],
                          nsig, nrow(th), PHH_THRESH_EXCL,
                          stats::median(th$p_far, na.rm = TRUE), near, recon,
                          if (is.finite(recon) && recon > 1e-6)
                              " The sweep did not reproduce the stored split."
                          else if (nsig == 0L)
                              " No outcome singles IP_y out."
                          else ""))
        }
    }

    # GR-10H-09: declared limitations (always PASS).
    gr10h("GR-10H-09", "PASS",
          paste0("declared limitations: batch is not independently ",
                 "adjusted (library nests within donor); the intronic:exonic ",
                 "ratio in covT is a per-library ambient proxy, not a per-cell ",
                 "correction; the 4- and 7-covariate models are bounds"))

    gr10h_write()
    cat(sprintf("\n  Final guardrail tally: %d PASS / %d WARN / %d FAIL\n",
                .gr_pass, .gr_warn, .gr_fail))

    # =========================================================================
    # STEP 5: Results RDS (written last, with the guardrail ledger)
    # =========================================================================
    if (isTRUE(SAVE_RESULTS_RDS)) {
        cat("\nSTEP 5 -- reloadable results object\n")

        # Ledger stored as lines and as a frame.
        .gr_df <- tryCatch({
            m <- regmatches(.gr_log, regexec(
                "^\\s*\\[\\s*([A-Z]+)\\s*\\]\\s+(\\S+)\\s+(.*)$", .gr_log))
            keep <- vapply(m, length, integer(1)) == 4L
            if (!any(keep)) NULL else data.frame(
                status = vapply(m[keep], `[`, character(1), 2L),
                check_id = vapply(m[keep], `[`, character(1), 3L),
                detail = vapply(m[keep], `[`, character(1), 4L),
                stringsAsFactors = FALSE)
        }, error = function(e) NULL)

        p10h_results <- list(
            meta = list(
                part = PART_TAG, part_dir = PART_DIR_TAG,
                arm_tag = ARM_TAG, arm_requested = arm_tag,
                config_file = config_file, norm_tag = NORM_TAG,
                GOI = GOI, DFG = DFG, supertype = PRIMARY_STYPE,
                datasets = DATASETS, datasets_completed = ok_ds,
                file_prefix = FILE_PREFIX, out_dir = OUT_DIR,
                iprs_source = if (is.na(.p10h_iprs$path)) NA_character_
                              else .p10h_iprs$path,
                written = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                R_version = R.version.string,
                figures_stored = isTRUE(SAVE_FIGURE_OBJECTS),
                schema = "M10H-33 v1"),

            guardrails = list(lines = .gr_log, table = .gr_df,
                              pass = .gr_pass, warn = .gr_warn,
                              fail = .gr_fail),

            # Every declaration a figure or table reads.
            declarations = list(
                ARM_GROUPS = ARM_GROUPS, PRED_GROUPS = PRED_GROUPS,
                PRED_LABEL = PRED_LABEL, PRED_COL = PRED_COL,
                CMRF_CLASS_LEVS = CMRF_CLASS_LEVS,
                MODULE_CLASSES = MODULE_CLASSES,
                NONMODULE_CLASSES = NONMODULE_CLASSES,
                CLINICAL_PARAMS = CLINICAL_PARAMS,
                VAR_ORDER_DISPLAY = VAR_ORDER_DISPLAY,
                STEIGER_PAIRS = STEIGER_PAIRS,
                MV_ARM_ORDER = MV_ARM_ORDER, MV_ARM_FILL = MV_ARM_FILL,
                MV_ARM_BORDER = MV_ARM_BORDER,
                MV_OUTCOME_ORDER = MV_OUTCOME_ORDER,
                MV_CONTROL_OUTCOMES = MV_CONTROL_OUTCOMES,
                MV_PANEL_LABEL = MV_PANEL_LABEL,
                MV_FIGURE_SETS = MV_FIGURE_SETS,
                MV_VIF_FIGURE_MAX = MV_VIF_FIGURE_MAX,
                MV_COV_TITLE = MV_COV_TITLE, MV_COV_NOTE = MV_COV_NOTE,
                MV_MIN_BETA_FOR_PCT = MV_MIN_BETA_FOR_PCT,
                COV_MAP = COV_MAP, COV_LABEL = COV_LABEL,
                COV_SCAN = COV_SCAN,
                COV_SETS = list(cov4 = COV_SET_4, cov7 = COV_SET_7,
                                covM = COV_SET_M, covT = COV_SET_T,
                                covP = COV_SET_P, covE = COV_SET_E),
                DONOR_EXTRA = DONOR_EXTRA,
                JT_TERTILE_FILL = JT_TERTILE_FILL,
                JT_ORDINAL_VARS = JT_ORDINAL_VARS,
                PHH_RESTRICTIONS = PHH_RESTRICTIONS, PHH_DUELS = PHH_DUELS,
                PHH_OUTCOMES = PHH_OUTCOMES, PHH_BOOT = PHH_BOOT,
                PHH_SEED = PHH_SEED, PHH_MIN_N = PHH_MIN_N,
                PHH_THRESH_GRID = PHH_THRESH_GRID,
                PHH_THRESH_EXCL = PHH_THRESH_EXCL,
                thresholds = list(MIN_DONORS = MIN_DONORS,
                                  MIN_CELLS_GRP = MIN_CELLS_GRP,
                                  MIN_MV_DONORS = MIN_MV_DONORS,
                                  JT_MIN_DONORS = JT_MIN_DONORS,
                                  JT_MIN_CELLS_PER_DONOR =
                                      JT_MIN_CELLS_PER_DONOR,
                                  SPEARMAN_Z_VAR_FACTOR =
                                      SPEARMAN_Z_VAR_FACTOR,
                                  CMRF_MISS_NOMINAL = CMRF_MISS_NOMINAL,
                                  CMRF_MISS_TOL = CMRF_MISS_TOL)),

            donor_clinical = donor_clinical,
            ordinal_audit  = ordinal_audit_df,
            covariate_audit = cov_audit,
            covariate_source = cov_source,

            # The workbook, sheet for sheet.
            sheets = sheets,

            # Everything per dataset.
            datasets = all_res,

            figures = if (isTRUE(SAVE_FIGURE_OBJECTS)) .fig_store else NULL)

        rds_path <- file.path(OUT_DIR,
                              paste0(FILE_PREFIX, "_Prognosis_results.rds"))
        tryCatch({
            saveRDS(p10h_results, rds_path)
            cat(sprintf(paste0("  %s written: %d sheet(s), %d dataset(s),",
                               " %d figure(s), %.2f MB\n"),
                        basename(rds_path), length(sheets), length(all_res),
                        length(p10h_results$figures),
                        file.size(rds_path) / 1e6))
        }, error = function(e)
            cat(sprintf("  Results RDS not written: %s\n",
                        conditionMessage(e))))
    }

    cat("\n", strrep("=", 72), "\n")
    cat("  Part10h Prognosis COMPLETE [", ARM_TAG, "]\n")
    cat("  Output directory:", OUT_DIR, "\n")
    cat("  Files generated:\n")
    .tr <- c("_RUN_INCOMPLETE.txt", basename(MSG_PATH_10H))
    for (f in setdiff(sort(list.files(OUT_DIR)), .tr)) {
        sz <- round(file.size(file.path(OUT_DIR, f)) / 1e3, 0)
        cat(sprintf("    %-62s  %5d KB\n", f, sz))
    }
    cat(strrep("=", 72), "\n")
    cat("\n  Log written to:", LOG_PATH_10H, "\n")
    cat(sprintf("\nPart10h Prognosis - log ended: %s\n", Sys.time()))

    .p10h_completed <- TRUE
    if (sink.number(type = "message") != 2L) sink(type = "message")
    while (sink.number() > 0) sink()
    suppressWarnings(file.remove(INCOMPLETE_MARKER))

    invisible(list(arm = arm_tag, dir = OUT_DIR, pass = .gr_pass,
                   warn = .gr_warn, fail = .gr_fail))
}


# =============================================================================
# SECTION 2: DRIVER
# =============================================================================

cep_run_arms(PART_DIR_TAG, run_one_arm)
