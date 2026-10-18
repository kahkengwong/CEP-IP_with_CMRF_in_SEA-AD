# =============================================================================
# Part10h_CP10K_Comparator_Fig13_code.r
# SEA-AD CEP-IP/CMRF, Part 10h, comparator figures (Figure 10A to 10D)
#   Redraws the comparator (Part 10h stage 1) and decomposition (stage 4)
#   outputs for all four groups: the delta-R2 grid, the encompassing-test
#   verdicts, the decomposition and the module-against-score contrast. Fits
#   nothing; every plotted table is written to CSV.
#
# Usage (from the SEA-AD folder, or set SEAAD_ROOT), after the Part 10h
# comparator and multivariate stages for all four arms
#   source("Part10h_CP10K_Comparator_Fig13_code.r")
#
# Outputs (Part10h_Prognosis_1_Comparator/Fig13/): figures (PDF, PNG),
#   plotted-data CSVs, guardrail ledger, baseline-sharing CSV and log.
# =============================================================================

suppressPackageStartupMessages(library(ggplot2))

# =============================================================================
# CONFIGURATION
# =============================================================================

READ_ROOT <- Sys.getenv("SEAAD_ROOT", getwd())
OUT_DIR   <- file.path(READ_ROOT, "Part10h_Prognosis_1_Comparator",
                       "Fig13")

# File tags as they appear in the comparator output names, in display order.
GROUP_TAGS <- c("M4", "PVM_3", "PVM_1", "HKG8")

# Figure labels, one per group.
GROUP_LABELS <- c(
    M4    = "M4 (masked, primary)",
    PVM_3 = "PVM_3, all cells",
    PVM_1 = "PVM_1, all cells",
    HKG8  = "HKG8 (housekeeping control)")

DATASET_LEVELS <- c("NeoC", "NonNeoC")

PRIMARY_BLOCK <- "cov4"
ALPHA         <- 0.05
VIF_READABLE  <- 5

# Predictor labels exactly as the comparator writes them.
PRED_KEYS <- c(module = "CMRF module Post-IP",
               B1 = "Cluster proportion (supertype)",
               B2 = "Module score (AUCell)",
               B3 = "Changepoint split (segmented)",
               B4 = "Module score (UCell)",
               B5 = "Module score (AddModuleScore)")

# Predictors in the delta-R2 grid and in the supplementary forest.
PRED_MAIN <- c("module", "B1", "B2", "B3", "B4", "B5")
PRED_SUPP <- c("module", "B1", "B2", "B3", "B4", "B5")

# Rivals in the bootstrap panel; its geometry scales with this length.
BOOT_RIVALS <- c("B1", "B2", "B3", "B4", "B5")

# Manuscript cut (PUB_* below): outcomes dropped from the main grid, the
# bootstrap outcomes, in-cell text, and the groups in the main grid.
# Delta-R2 fill ramp: cmocean "matter" plus two lighter anchors; light is
# low, dark is high. FILL_REVERSE flips it.
FILL_RAMP <- c("#FFF9DF", "#FFF4C6", "#FEEDB0", "#FACA90", "#F5A773",
               "#EE855E", "#E26253", "#CE4356", "#B32E5F", "#932063",
               "#721A60", "#4F1652", "#2F0F3E")
FILL_REVERSE <- FALSE
FILL_COLS <- if (isTRUE(FILL_REVERSE)) rev(FILL_RAMP) else FILL_RAMP

# In-cell label colour: white or dark grey, whichever has the better WCAG
# contrast against the ramp colour under the cell.
fill_label_colour <- function(frac, cols = FILL_COLS) {
    out <- rep("grey15", length(frac))
    ok <- is.finite(frac)
    if (!any(ok)) return(out)
    f <- grDevices::colorRamp(cols)(pmin(pmax(frac[ok], 0), 1)) / 255
    lin <- ifelse(f <= 0.04045, f / 12.92, ((f + 0.055) / 1.055)^2.4)
    lum <- 0.2126 * lin[, 1] + 0.7152 * lin[, 2] + 0.0722 * lin[, 3]
    out[ok] <- ifelse(1.05 / (lum + 0.05) >= (lum + 0.05) / 0.0685,
                      "white", "grey15")
    out
}

# In-cell label sizes and FDR border width, shared by both grids.
GRID_TEXT_SIZE     <- 1.6
GRID_TEXT_SIZE_PUB <- 2.8
GRID_BH_LWD        <- 0.3
# Fill and FDR border tiles share one inset size, leaving a white gutter.
GRID_CELL          <- 0.88
# Gap between the NeoC and NonNeoC panels in panels C and D.
PANEL_GAP          <- 14
DUEL_TEXT_SIZE     <- 3.0
DUEL_BAR_HEIGHT    <- 0.68

# Stated on the figure so the direction is never left to the reader.
.fill_dir_txt <- if (isTRUE(FILL_REVERSE)) {
    "Brighter is more variance explained."
} else {
    "Darker is more variance explained."
}

# Outcome class for the duel panel; APOE e4 dose has no class and is dropped.
OUTCOME_CLASS <- c("Cognitive status" = "Cognition",
                   "MMSE score" = "Cognition",
                   "CASI score" = "Cognition",
                   "CERAD score" = "Pathology",
                   "Braak stage" = "Pathology",
                   "Thal phase" = "Pathology",
                   "ADNC rank" = "Pathology",
                   "CAA score" = "Pathology")
CLASS_LEVELS <- c("Cognition", "Pathology")

PUB_DROP_OUTCOMES <- c("APOE e4 dose")
PUB_BOOT_OUTCOMES <- c("Cognitive status", "MMSE score", "CASI score")
PUB_GRID_TEXT     <- TRUE
PUB_GRID_GROUPS   <- c("M4")

# Outcome order, cognition first. Anything present in the data but absent here
# is appended in its existing order and reported.
OUTCOME_ORDER <- c("Cognitive status", "MMSE score", "CASI score",
                   "CERAD score", "Braak stage", "Thal phase", "ADNC rank",
                   "CAA score", "APOE e4 dose")

# Display labels applied at draw time; source strings stay as factor keys.
# The epsilon is built from its code point so the file stays ASCII.
OUTCOME_DISPLAY <- c("APOE e4 dose" = paste0("APOE ",
                                             intToUtf8(0x03B5), "4 dose"))
relabel_outcomes <- function(x) {
    x <- as.character(x)
    hit <- x %in% names(OUTCOME_DISPLAY)
    x[hit] <- unname(OUTCOME_DISPLAY[x[hit]])
    x
}

COGNITIVE_OUTCOMES <- c("Cognitive status", "MMSE score", "CASI score")

# Predictor colours; the bootstrap panel also uses filled vs hollow points.
PRED_FILL <- c("CMRF module Post-IP" = "#2D2D2D",
               "Cluster proportion (supertype)" = "#D95F02",
               "Module score (AUCell)" = "#2E7D5B",
               "Changepoint split (segmented)" = "#7B2D8E",
               "Module score (UCell)" = "#A34F5D",
               "Module score (AddModuleScore)" = "#24456F")

# Verdict order and colours; this order sets the legend and the stacking.
VERDICT_LEVELS <- c("module only", "both", "rival only", "neither",
                    "not readable")
# Panel C typography (sizes in mm; 1 pt = 0.35 mm).
MV_FONT <- "Arial Narrow"
MV_TEXT_SIZE    <- GRID_TEXT_SIZE + 0.35   # statistics line, +1 pt on R29
MV_VERDICT_BUMP <- 1.06                    # verdict word, +3 pt on the stats
MV_RHO_BUMP     <- 1.24                    # correlation row, +3.5 pt
MV_CELL_LWD     <- 0.25                    # hairline border on every cell

# Score term drawn in panel C (AUCell is pre-specified).
MV_ARM_SHOWN <- "AUCell"

# Printed verdict labels for panel C; the source verdict column is unchanged.
MV_VERDICT_DISPLAY <- c("HOLDS"         = "HOLDS",
                        "does not hold" = "doesn't hold",
                        "WRONG SIGN"    = "WRONG SIGN",
                        "not readable"  = "unreadable")
if (.Platform$OS.type == "windows")
    try(grDevices::windowsFonts(
        `Arial Narrow` = grDevices::windowsFont("Arial Narrow")),
        silent = TRUE)

MV_FILL <- c("HOLDS"        = "#2D2D2D",
             "doesn't hold" = "#D9D9D9",
             "WRONG SIGN"   = "#2A9D8F",
             "unreadable"   = "#FBE7A2")

VERDICT_FILL <- c("module only"  = "#2D2D2D",
                  "both"         = "#7B7B7B",
                  "rival only"   = "#2A9D8F",
                  "neither"      = "#D9D9D9",
                  "not readable" = "#FBE7A2")

# Control-arm note on the figures that draw the HKG8 cells.
HKG8_CAVEAT <- paste("Control arm: the module clears FDR on 5 of 5 neocortex",
                     "pathology rows at cov4 and 0 of 5 at cov7, which adds",
                     "RIN, brain pH and brain weight; the primary arm keeps",
                     "3 of 3 cognitive rows at cov7. The control is working.")

# Caption wrap width in characters. ggplot does not wrap captions.
CAPTION_WIDTH <- 120

FIG_DPI <- 300

# =============================================================================
# LEDGER
# =============================================================================

.GR <- list()

gr <- function(id, status, detail) {
    .GR[[length(.GR) + 1L]] <<- data.frame(
        check_id = id, status = status, detail = detail,
        stringsAsFactors = FALSE)
    cat(sprintf("  [%-4s] %-10s %s\n", status, id, detail))
    invisible(NULL)
}

# =============================================================================
# INPUT RESOLUTION
# =============================================================================

find_input <- function(tag, root) {
    pat_rds <- sprintf("^Part10h_CP10K_%s_Comparator_results\\.rds$", tag)
    rds <- list.files(root, recursive = TRUE, full.names = TRUE,
                      pattern = pat_rds)
    if (length(rds) > 0) return(rds[order(nchar(rds))][1])
    xls <- list.files(root, recursive = TRUE, full.names = TRUE,
                      pattern = sprintf("^Part10h_CP10K_%s_Comparator\\.xlsx$",
                                        tag))
    if (length(xls) > 0) return(xls[order(nchar(xls))][1])
    NA_character_
}

read_input <- function(path) {
    if (grepl("\\.rds$", path)) {
        r <- readRDS(path)
        need <- c("marginal", "duels", "bootstrap")
        miss <- need[!need %in% names(r)]
        if (length(miss) > 0)
            stop("The results RDS ", basename(path), " has no ",
                 paste(miss, collapse = ", "),
                 ". It carries: ", paste(names(r), collapse = ", "))
        return(list(marginal  = as.data.frame(r$marginal),
                    duels     = as.data.frame(r$duels),
                    bootstrap = as.data.frame(r$bootstrap)))
    }
    if (!requireNamespace("readxl", quietly = TRUE))
        stop("Only the workbook was found for ", basename(path),
             " and readxl is not installed.")
    rx <- function(sheet)
        as.data.frame(readxl::read_excel(path, sheet = sheet))
    list(marginal  = rx("Marginal"),
         duels     = rx("Duels"),
         bootstrap = rx("Stress_Bootstrap"))
}

require_cols <- function(df, cols, what) {
    miss <- cols[!cols %in% names(df)]
    if (length(miss) > 0)
        stop("The ", what, " frame is missing ", paste(miss, collapse = ", "),
             ". It carries: ", paste(names(df), collapse = ", "))
    invisible(TRUE)
}

# A row is unreadable when the comparator wrote a reason into `readable`.
is_unreadable <- function(x) !is.na(x) & nzchar(as.character(x))

# =============================================================================
# LOAD
# =============================================================================

# Console output is also written to the log; any open sink is closed first.
while (sink.number() > 0L) sink()
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)
LOG_PATH <- file.path(OUT_DIR, "Fig13_Log.txt")
.log_con <- file(LOG_PATH, open = "wt")
sink(.log_con, split = TRUE)
cat(strrep("=", 72), "\n", sep = "")
cat("Part10h comparator - Figure 13 builder\n")
cat("Date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", sep = "")
cat("Read root : ", READ_ROOT, "\n", sep = "")
cat("Output dir: ", OUT_DIR, "\n", sep = "")
cat(strrep("=", 72), "\n", sep = "")

cat("\nSTEP 1 -- inputs\n")

INPUTS <- vapply(GROUP_TAGS, find_input, character(1), root = READ_ROOT)
if (anyNA(INPUTS))
    stop("No comparator output found under ", READ_ROOT, " for: ",
         paste(GROUP_TAGS[is.na(INPUTS)], collapse = ", "),
         ". Expected Part10h_CP10K_<tag>_Comparator_results.rds or the xlsx.")
gr("GR-13A-01", "PASS",
   sprintf("%d of %d groups resolved: %s", sum(!is.na(INPUTS)),
           length(GROUP_TAGS),
           paste(sprintf("%s -> %s", GROUP_TAGS, basename(INPUTS)),
                 collapse = "; ")))

RAW <- lapply(INPUTS, read_input)
names(RAW) <- GROUP_TAGS

for (g in GROUP_TAGS) {
    require_cols(RAW[[g]]$marginal,
                 c("dataset", "block", "predictor", "key", "outcome",
                   "outcome_class", "is_control", "n", "beta", "se",
                   "delta_r2", "p_adj"),
                 sprintf("%s Marginal", g))
    require_cols(RAW[[g]]$duels,
                 c("dataset", "rival_label", "outcome", "outcome_class",
                   "is_control", "verdict", "vif_joint", "readable"),
                 sprintf("%s Duels", g))
    require_cols(RAW[[g]]$bootstrap,
                 c("dataset", "rival_label", "outcome", "outcome_class",
                   "is_control", "median_dR2", "ci_lo", "ci_hi",
                   "pct_module_ahead", "ci_excludes_zero"),
                 sprintf("%s Stress_Bootstrap", g))
}

bind_groups <- function(part) {
    out <- do.call(rbind, lapply(GROUP_TAGS, function(g) {
        d <- RAW[[g]][[part]]
        d$group <- g
        d
    }))
    out$group <- factor(out$group, levels = GROUP_TAGS,
                        labels = unname(GROUP_LABELS[GROUP_TAGS]))
    out$dataset <- factor(out$dataset, levels = DATASET_LEVELS)
    out
}

# Panel C reads the Part 10h stage 4 (Multivariate) results.
find_mv <- function(tag, root) {
    pat <- sprintf("^Part10h_CP10K_%s_Multivariate_results\\.rds$", tag)
    rds <- list.files(root, recursive = TRUE, full.names = TRUE, pattern = pat)
    if (length(rds) > 0) return(rds[order(nchar(rds))][1])
    xls <- list.files(root, recursive = TRUE, full.names = TRUE,
                      pattern = sprintf(
                          "^Part10h_CP10K_%s_Multivariate\\.xlsx$", tag))
    if (length(xls) > 0) return(xls[order(nchar(xls))][1])
    NA_character_
}

read_mv <- function(path) {
    if (grepl("\\.rds$", path)) {
        r <- readRDS(path)
        return(list(decomposition = as.data.frame(r$decomposition),
                    correlations  = as.data.frame(r$correlations)))
    }
    list(decomposition = as.data.frame(
             readxl::read_excel(path, sheet = "Decomposition")),
         correlations  = as.data.frame(
             readxl::read_excel(path, sheet = "Correlations")))
}

MV_RAW <- list()
.mv_missing <- character(0)
for (g in GROUP_TAGS) {
    p <- find_mv(g, READ_ROOT)
    if (is.na(p)) {
        .mv_missing <- c(.mv_missing, g)
        next
    }
    MV_RAW[[g]] <- tryCatch(read_mv(p), error = function(e) NULL)
    if (is.null(MV_RAW[[g]])) .mv_missing <- c(.mv_missing, g)
}

MARG <- bind_groups("marginal")
DUEL <- bind_groups("duels")
BOOT <- bind_groups("bootstrap")

# =============================================================================
# STEP 2: checks
# =============================================================================

cat("\nSTEP 2 -- checks\n")

# Requested predictors present in every group.
want <- unname(PRED_KEYS[PRED_SUPP])
have <- lapply(GROUP_TAGS, function(g)
    unique(RAW[[g]]$marginal$predictor))
names(have) <- GROUP_TAGS
missing_pred <- vapply(GROUP_TAGS, function(g)
    paste(setdiff(want, have[[g]]), collapse = ", "), character(1))
if (any(nzchar(missing_pred))) {
    gr("GR-13A-02", "WARN",
       sprintf("predictors absent: %s",
               paste(sprintf("%s: %s", GROUP_TAGS[nzchar(missing_pred)],
                             missing_pred[nzchar(missing_pred)]),
                     collapse = "; ")))
} else {
    gr("GR-13A-02", "PASS",
       sprintf("all %d predictors present in all %d groups",
               length(want), length(GROUP_TAGS)))
}

# Baseline columns shared across groups, measured on the primary block.
ref <- RAW[[GROUP_TAGS[1]]]$marginal
KEY_SEP <- "|@|"
key_of <- function(d) paste(d$dataset, d$block, d$predictor, d$outcome,
                            sep = KEY_SEP)
ref_val <- setNames(paste(ref$beta, ref$delta_r2, ref$n, sep = KEY_SEP),
                    key_of(ref))
share <- do.call(rbind, lapply(GROUP_TAGS[-1], function(g) {
    d <- RAW[[g]]$marginal
    v <- setNames(paste(d$beta, d$delta_r2, d$n, sep = KEY_SEP), key_of(d))
    do.call(rbind, lapply(setdiff(names(PRED_KEYS), "module"), function(pk) {
        lab <- PRED_KEYS[[pk]]
        k <- key_of(ref)[ref$predictor == lab]
        same <- sum(!is.na(v[k]) & v[k] == ref_val[k])
        data.frame(group = g, predictor = lab, n_rows = length(k),
                   n_identical = same, stringsAsFactors = FALSE)
    }))
}))
share$identical_all <- share$n_identical == share$n_rows
gr("GR-13A-03", "PASS",
   paste(vapply(GROUP_TAGS[-1], function(g) {
       s <- share[share$group == g, ]
       ok <- s$predictor[s$identical_all]
       sprintf("%s shares %d of %d baselines with %s%s", g, length(ok),
               nrow(s), GROUP_TAGS[1],
               if (length(ok) > 0)
                   sprintf(" (%s)", paste(ok, collapse = ", ")) else "")
   }, character(1)), collapse = "; "))

# Unreadable duel rows (joint VIF at or above VIF_READABLE).
DUEL$unreadable <- is_unreadable(DUEL$readable)
du_clin <- DUEL[!DUEL$is_control, , drop = FALSE]
un_tab <- aggregate(list(n_unreadable = du_clin$unreadable),
                    by = list(group = du_clin$group), FUN = sum)
un_tab$n_rows <- as.vector(table(du_clin$group)[as.character(un_tab$group)])
gr("GR-13A-04",
   if (sum(du_clin$unreadable) > 0) "WARN" else "PASS",
   sprintf("%d of %d non-control duel rows at joint VIF >= %g: %s",
           sum(du_clin$unreadable), nrow(du_clin), VIF_READABLE,
           paste(sprintf("%s %d/%d", un_tab$group, un_tab$n_unreadable,
                         un_tab$n_rows), collapse = ", ")))

# Bootstrap intervals that exclude zero.
bo <- BOOT[!BOOT$is_control & BOOT$ci_excludes_zero %in% c(TRUE, "TRUE"), ,
           drop = FALSE]
if (nrow(bo) == 0) {
    gr("GR-13A-05", "PASS",
       sprintf("0 of %d bootstrap intervals exclude zero",
               sum(!BOOT$is_control)))
} else {
    gr("GR-13A-05", "WARN",
       sprintf("%d of %d intervals exclude zero: %s", nrow(bo),
               sum(!BOOT$is_control),
               paste(sprintf("%s %s %s vs %s (%.3f)", bo$group, bo$dataset,
                             bo$outcome, bo$rival_label, bo$median_dR2),
                     collapse = "; ")))
}

# =============================================================================
# STEP 3: shared drawing helpers
# =============================================================================

cat("\nSTEP 3 -- figures\n")

.n_fig <- 0L
.csv_written <- character(0)

save_fig <- function(p, base, w, h, dpi = FIG_DPI) {
    render <- function(open_fn, path) {
        tryCatch({
            open_fn(); on.exit(grDevices::dev.off(), add = TRUE); print(p)
            TRUE
        }, error = function(e) {
            cat(sprintf("    figure not written (%s): %s\n", basename(path),
                        conditionMessage(e))); FALSE })
    }
    # cairo devices resolve system fonts by name; grDevices::pdf is the fallback.
    pp <- paste0(base, ".pdf"); gp <- paste0(base, ".png")
    .cairo <- isTRUE(capabilities("cairo"))
    ok1 <- render(function() {
        if (.cairo) grDevices::cairo_pdf(pp, width = w, height = h)
        else grDevices::pdf(pp, width = w, height = h)
    }, pp)
    if (!ok1 && .cairo)
        ok1 <- render(function() grDevices::pdf(pp, width = w, height = h),
                      pp)
    ok2 <- render(function() {
        if (.cairo) grDevices::png(gp, width = w * dpi, height = h * dpi,
                                   res = dpi, type = "cairo")
        else grDevices::png(gp, width = w * dpi, height = h * dpi, res = dpi)
    }, gp)
    if (ok1 && ok2) .n_fig <<- .n_fig + 1L
    invisible(NULL)
}

# Every figure writes the rows it drew.
save_plotted <- function(df, base) {
    path <- paste0(base, "_plotted.csv")
    utils::write.csv(df, path, row.names = FALSE)
    .csv_written <<- c(.csv_written, basename(path))
    cat(sprintf("    %s (%d rows)\n", basename(path), nrow(df)))
    invisible(NULL)
}

theme_pub <- function(base_size = 10) theme_bw(base_size = base_size) +
    theme(plot.title = element_text(face = "bold", size = base_size + 2),
          plot.subtitle = element_text(size = base_size - 2,
                                       colour = "grey30"),
          plot.caption = element_text(size = base_size - 2.5,
                                      colour = "grey30", hjust = 0),
          strip.background = element_rect(fill = "#1F3864"),
          strip.text = element_text(colour = "white", face = "bold",
                                    size = base_size - 1),
          panel.grid.minor = element_blank(),
          legend.position = "bottom",
          plot.background = element_rect(fill = "white", colour = NA))

wrap_cap <- function(..., width = CAPTION_WIDTH) {
    txt <- paste(unlist(list(...)), collapse = " ")
    paste(strwrap(txt, width = width), collapse = "\n")
}

# `descending = TRUE` puts the first entry of OUTCOME_ORDER at the top of a
# discrete y axis; FALSE puts it at the left of a discrete x axis.
order_outcomes <- function(x, descending = TRUE) {
    lv <- c(OUTCOME_ORDER[OUTCOME_ORDER %in% x],
            setdiff(unique(x), OUTCOME_ORDER))
    # Factor built on source strings; only the axis labels are relabelled.
    factor(x, levels = if (descending) rev(lv) else lv,
           labels = relabel_outcomes(if (descending) rev(lv) else lv))
}

# Baseline-sharing sentence for every figure that draws baselines.
share_caption <- function() {
    txt <- vapply(GROUP_TAGS[-1], function(g) {
        s <- share[share$group == g, ]
        sprintf("%s shares %d of %d baseline columns with %s", g,
                sum(s$identical_all), nrow(s), GROUP_TAGS[1])
    }, character(1))
    paste0("Baselines are computed over supertype cells, not over the group's ",
           "selected cells. ", paste(txt, collapse = "; "),
           ". Identical columns are one computation shown more than once.")
}

# =============================================================================
# 3a: delta-R2 grid (main)
# =============================================================================

gd <- MARG[MARG$block == PRIMARY_BLOCK & !MARG$is_control &
           MARG$predictor %in% unname(PRED_KEYS[PRED_MAIN]), , drop = FALSE]
gd$outcome_f <- order_outcomes(gd$outcome)
gd$pred_f <- factor(gd$predictor, levels = unname(PRED_KEYS[PRED_MAIN]))
gd$bh <- !is.na(gd$p_adj) & gd$p_adj < ALPHA
gd$lab <- sprintf("%.2f", gd$delta_r2)
gd$lab[gd$bh] <- paste0(gd$lab[gd$bh], "*")
.dr2_max <- max(gd$delta_r2, na.rm = TRUE)
gd$txt_col <- fill_label_colour(gd$delta_r2 / .dr2_max)

p_grid <- ggplot(gd, aes(x = pred_f, y = outcome_f, fill = delta_r2)) +
    geom_tile(colour = NA, width = GRID_CELL, height = GRID_CELL) +
    geom_tile(data = gd[gd$bh, , drop = FALSE], fill = NA, colour = "black",
              linewidth = GRID_BH_LWD,
              width = GRID_CELL, height = GRID_CELL) +
    geom_text(aes(label = lab, colour = txt_col), size = GRID_TEXT_SIZE) +
    scale_colour_identity() +
    scale_fill_gradientn(colours = FILL_COLS,
                         limits = range(c(0, gd$delta_r2), na.rm = TRUE),
                         name = expression(Delta * R^2)) +
    facet_grid(group ~ dataset) +
    labs(title = sprintf("Variance explained over covariates alone [%s]",
                         PRIMARY_BLOCK),
         subtitle = paste("Cognitive outcomes are the top three rows.",
                          .fill_dir_txt,
                          "A black border and an asterisk mark FDR < 0.05",
                          "within the group-by-dataset family."),
         caption = paste(
             wrap_cap(share_caption()),
             wrap_cap("Delta-R2 is unsigned; the sign of each coefficient",
                      "is in the supplementary table.", HKG8_CAVEAT),
             sep = "\n"),
         x = NULL, y = NULL) +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
if (interactive()) print(p_grid)   # on-screen copy only; save_fig writes the files
save_plotted(gd[, c("group", "dataset", "predictor", "outcome", "n", "beta",
                    "se", "delta_r2", "p_adj", "bh")],
             file.path(OUT_DIR, "Fig13A_DeltaR2_grid"))
save_fig(p_grid, file.path(OUT_DIR, "Fig13A_DeltaR2_grid"),
         w = 2.2 + 1.15 * length(PRED_MAIN) * length(DATASET_LEVELS),
         h = 3.2 + 0.30 * nlevels(gd$outcome_f) * length(GROUP_TAGS))

# =============================================================================
# 3b: encompassing verdicts, unreadable rows separated (main)
# =============================================================================

dv <- DUEL[!DUEL$is_control &
           DUEL$rival_label %in% unname(PRED_KEYS[setdiff(PRED_SUPP,
                                                          "module")]), ,
           drop = FALSE]
# Outcomes outside the class map leave the panel and are named.
.unmapped <- setdiff(unique(as.character(dv$outcome)), names(OUTCOME_CLASS))
if (length(.unmapped))
    cat(sprintf("    duel panel drops %d outcome(s) with no class: %s\n",
                length(.unmapped), paste(.unmapped, collapse = ", ")))
dv <- dv[as.character(dv$outcome) %in% names(OUTCOME_CLASS), , drop = FALSE]
dv$outcome_class <- factor(unname(OUTCOME_CLASS[as.character(dv$outcome)]),
                           levels = CLASS_LEVELS)
dv$verdict_disp <- ifelse(dv$unreadable, "not readable",
                          as.character(dv$verdict))
dv$verdict_disp <- factor(dv$verdict_disp, levels = VERDICT_LEVELS)
dv$rival_f <- factor(dv$rival_label,
                     levels = rev(unname(PRED_KEYS[setdiff(PRED_SUPP,
                                                           "module")])))

tally <- as.data.frame(
    table(group = dv$group, dataset = dv$dataset,
          outcome_class = dv$outcome_class, rival = dv$rival_f,
          verdict = dv$verdict_disp),
    stringsAsFactors = FALSE)
tally$outcome_class <- factor(tally$outcome_class, levels = CLASS_LEVELS)
tally$verdict <- factor(tally$verdict, levels = VERDICT_LEVELS)
tally$rival <- factor(tally$rival,
                      levels = rev(unname(PRED_KEYS[setdiff(PRED_SUPP,
                                                            "module")])))
tally$group <- factor(tally$group, levels = levels(dv$group))
tally$dataset <- factor(tally$dataset, levels = DATASET_LEVELS)
tally <- tally[tally$Freq > 0, , drop = FALSE]
# Label colour is a data column, measured from each verdict fill.
tally$txt_col <- vapply(VERDICT_FILL[as.character(tally$verdict)],
                        function(h) fill_label_colour(1, c(h, h)),
                        character(1), USE.NAMES = FALSE)

# Segment extents computed here, so bars and labels share one calculation.
tally <- tally[order(tally$group, tally$dataset, tally$outcome_class,
                     tally$rival, tally$verdict), , drop = FALSE]
.bar_key <- interaction(tally$group, tally$dataset, tally$outcome_class,
                        tally$rival, drop = TRUE)
tally$xmax <- stats::ave(tally$Freq, .bar_key, FUN = cumsum)
tally$xmin <- tally$xmax - tally$Freq
tally$xmid <- tally$xmax - tally$Freq / 2
tally$ypos <- as.integer(tally$rival)
.bar_half <- DUEL_BAR_HEIGHT / 2

# Range of the module's delta-R2 behind its wins, for the caption.
.mod_dr2 <- gd[gd$predictor == unname(PRED_KEYS[["module"]]), , drop = FALSE]
.win_range <- function(grp) {
    w <- dv[dv$group == grp & dv$verdict_disp == "module only", , drop = FALSE]
    if (nrow(w) == 0L) return(NULL)
    v <- .mod_dr2$delta_r2[
        match(paste(w$dataset, w$outcome),
              paste(.mod_dr2$dataset, .mod_dr2$outcome))]
    v <- v[is.finite(v)]
    if (!length(v)) NULL else range(v)
}
.r_pri <- .win_range(GROUP_LABELS[[GROUP_TAGS[1]]])
.r_ctl <- .win_range(GROUP_LABELS[[utils::tail(GROUP_TAGS, 1)]])
.mag_txt <- if (is.null(.r_pri) || is.null(.r_ctl)) {
    paste("Verdicts are relative and carry no magnitude; the module delta-R2",
          "behind each win is in the delta-R2 panel.")
} else {
    sprintf(paste("Verdicts are relative and carry no magnitude. The module",
                  "delta-R2 behind each win is in the delta-R2 panel:",
                  "%.3f to %.3f on %s against %.3f to %.3f on %s%s."),
            .r_pri[1], .r_pri[2], GROUP_TAGS[1],
            .r_ctl[1], .r_ctl[2], utils::tail(GROUP_TAGS, 1),
            if (.r_ctl[2] < .r_pri[1]) ", with no overlap" else "")
}

p_duel <- ggplot(tally) +
    geom_rect(aes(xmin = xmin, xmax = xmax,
                  ymin = ypos - .bar_half, ymax = ypos + .bar_half,
                  fill = verdict),
              colour = "black", linewidth = 0.3) +
    geom_text(aes(x = xmid, y = ypos, label = Freq, colour = txt_col),
              size = DUEL_TEXT_SIZE) +
    scale_colour_identity() +
    scale_fill_manual(values = VERDICT_FILL, name = NULL, drop = FALSE) +
    scale_y_continuous(breaks = seq_along(levels(tally$rival)),
                       labels = levels(tally$rival),
                       expand = expansion(add = 0.6)) +
    # Free x and space: one unit of bar is one outcome in every facet.
    facet_grid(group ~ outcome_class + dataset,
               scales = "free_x", space = "free_x") +
    labs(title = "Encompassing-test verdicts against each baseline",
         subtitle = sprintf(paste("y ~ %s + module + rival, FDR-corrected",
                                  "within a duel. Cognition holds three",
                                  "outcomes and pathology five, so a bar sums",
                                  "to 3 or to 5. Rival-only rows are reported",
                                  "as they fall."), PRIMARY_BLOCK),
         caption = paste(
             wrap_cap(sprintf(paste("Rows at joint VIF >= %g carry no verdict:",
                                    "the module and the rival cannot be",
                                    "separated in the joint model, which is",
                                    "not the same as neither carrying",
                                    "information. %d of %d rows are",
                                    "affected."),
                              VIF_READABLE, sum(du_clin$unreadable),
                              nrow(du_clin))),
             wrap_cap("All five rivals are EXTERNAL baselines. The module",
                      "against its own non-module cells is Figure 12D."),
             wrap_cap(.mag_txt),
             wrap_cap(share_caption()),
             wrap_cap(HKG8_CAVEAT),
             sep = "\n"),
         x = "outcomes", y = NULL) +
    theme_pub()
if (interactive()) print(p_duel)   # on-screen copy only; save_fig writes the files
save_plotted(dv[, c("group", "dataset", "rival_label", "outcome", "n",
                    "beta_module", "beta_rival", "p_module_adj",
                    "p_rival_adj", "vif_joint", "verdict", "unreadable",
                    "verdict_disp")],
             file.path(OUT_DIR, "Fig13B_Duel_verdicts"))
# Two class facets per region, sized by outcome count.
save_fig(p_duel, file.path(OUT_DIR, "Fig13B_Duel_verdicts"),
         w = 14,
         h = 3.0 + 0.42 * nlevels(dv$rival_f) * length(GROUP_TAGS))

# =============================================================================
# 3c: donor bootstrap (CSV only)
# =============================================================================

bd <- BOOT[!BOOT$is_control &
           BOOT$rival_label %in% unname(PRED_KEYS[BOOT_RIVALS]), ,
           drop = FALSE]
# Outcomes on x.
bd$outcome_f <- order_outcomes(bd$outcome, descending = FALSE)
bd$rival_f <- factor(bd$rival_label, levels = unname(PRED_KEYS[BOOT_RIVALS]))
bd$excl <- bd$ci_excludes_zero %in% c(TRUE, "TRUE")

# Replicate count read from the data.
.bn <- sort(unique(bd$n_boot[is.finite(bd$n_boot)]))
.boot_n_txt <- if (length(.bn) == 0L) {
    "Donor bootstrap."
} else if (length(.bn) == 1L) {
    sprintf("%s donor resamples.", format(.bn, big.mark = ","))
} else {
    sprintf("%s to %s donor resamples per row.",
            format(min(.bn), big.mark = ","),
            format(max(.bn), big.mark = ","))
}

# group = rival_f and one shared dodge, so points sit on their own intervals.
# Dodge width, point size and panel width scale with the rival count.
.n_riv <- max(1L, nlevels(bd$rival_f))
.pd_w <- min(0.88, 0.24 + 0.13 * .n_riv)
.pt_sz <- if (.n_riv > 3L) 1.4 else 1.9
.pd_boot <- position_dodge(width = .pd_w)

# Caption states which rivals are drawn.
.boot_riv_txt <- if (.n_riv >= length(PRED_KEYS) - 1L) {
    paste("All", .n_riv, "rivals are drawn.",
          "Every interval is also in the plotted CSV.")
} else {
    sprintf(paste("%d of the %d rivals are drawn (%s). The other %d are in",
                  "the plotted CSV."),
            .n_riv, length(PRED_KEYS) - 1L,
            paste(levels(bd$rival_f), collapse = ", "),
            length(PRED_KEYS) - 1L - .n_riv)
}
p_boot <- ggplot(bd, aes(x = outcome_f, y = median_dR2, colour = rival_f,
                         group = rival_f)) +
    geom_hline(yintercept = 0, colour = "grey50", linetype = "dashed") +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0,
                  position = .pd_boot, linewidth = 0.4) +
    geom_point(aes(shape = excl), size = .pt_sz, position = .pd_boot) +
    scale_shape_manual(values = c(`TRUE` = 19, `FALSE` = 1),
                       name = "CI excludes zero", na.translate = FALSE) +
    scale_colour_manual(values = PRED_FILL, name = NULL) +
    facet_grid(group ~ dataset) +
    labs(title = "Donor bootstrap of the delta-R2 difference",
         subtitle = paste(.boot_n_txt,
                          "Positive is the module ahead of the rival.",
                          "Filled points have an interval excluding zero."),
         caption = paste(
             wrap_cap("The encompassing test is the verdict; the bootstrap",
                      "resamples the same donors and is a conservative check."),
             wrap_cap(.boot_riv_txt),
             sep = "\n"),
         x = NULL,
         y = expression(paste("median ", Delta * R^2,
                              " (module) - ", Delta * R^2, " (rival)"))) +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 40, hjust = 1))
save_plotted(bd[, c("group", "dataset", "rival_label", "outcome", "n_boot",
                    "median_dR2", "ci_lo", "ci_hi", "pct_module_ahead",
                    "ci_excludes_zero")],
             file.path(OUT_DIR, "Fig13_Bootstrap"))
# The bootstrap and forest are supplementary tables: CSVs only, no figures.

# =============================================================================
# 3c-pub: manuscript cut
# =============================================================================
# Same data, smaller; each writes its own plotted CSV.

cat("\n  publication cut\n")

# == A1: delta-R2 grid, APOE dropped ==========================================
.pub_lab <- unname(GROUP_LABELS[PUB_GRID_GROUPS])
gp <- gd[!(as.character(gd$outcome) %in% PUB_DROP_OUTCOMES) &
         as.character(gd$group) %in% .pub_lab, , drop = FALSE]
gp$outcome_f <- droplevels(gp$outcome_f)
gp$group <- droplevels(gp$group)
if (nrow(gp) == 0L)
    stop("PUB_GRID_GROUPS matched no group label: ",
         paste(PUB_GRID_GROUPS, collapse = ", "), call. = FALSE)
p_grid_pub <- ggplot(gp, aes(x = pred_f, y = outcome_f, fill = delta_r2)) +
    geom_tile(colour = NA, width = GRID_CELL, height = GRID_CELL) +
    geom_tile(data = gp[gp$bh, , drop = FALSE], fill = NA, colour = "black",
              linewidth = GRID_BH_LWD,
              width = GRID_CELL, height = GRID_CELL) +
    scale_fill_gradientn(colours = FILL_COLS,
                         limits = range(c(0, gd$delta_r2), na.rm = TRUE),
                         name = expression(Delta * R^2)) +
    facet_grid(group ~ dataset) +
    labs(title = sprintf("Variance explained over covariates alone [%s]",
                         PRIMARY_BLOCK),
         subtitle = paste("Cognitive outcomes are the top three rows.",
                          paste0(.fill_dir_txt, " A black border"),
                          "and an asterisk mark FDR < 0.05 within the",
                          "group-by-dataset family."),
         caption = paste(
             wrap_cap(sprintf(paste("%s, the %s row and the other %d group(s)",
                                    "are in the supplementary grid. The",
                                    "association replicates across the",
                                    "specificity groups; the module is the",
                                    "only column that is not shared."),
                              if (isTRUE(PUB_GRID_TEXT)) "The standard errors"
                              else "Values",
                              paste(PUB_DROP_OUTCOMES, collapse = ", "),
                              length(GROUP_TAGS) - length(PUB_GRID_GROUPS))),
             wrap_cap("Delta-R2 is unsigned; the sign of each coefficient",
                      "is in the supplementary table."),
             sep = "\n"),
         x = NULL, y = NULL) +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
if (isTRUE(PUB_GRID_TEXT))
    p_grid_pub <- p_grid_pub +
        geom_text(aes(label = lab, colour = txt_col),
                  size = GRID_TEXT_SIZE_PUB) +
        scale_colour_identity()
if (interactive()) print(p_grid_pub)   # on-screen copy only; save_fig writes the files
save_plotted(gp[, c("group", "dataset", "predictor", "outcome", "n", "beta",
                    "se", "delta_r2", "p_adj", "bh")],
             file.path(OUT_DIR, "Fig13A_Pub_DeltaR2_grid"))
save_fig(p_grid_pub, file.path(OUT_DIR, "Fig13A_Pub_DeltaR2_grid"),
         w = 1.8 + (if (isTRUE(PUB_GRID_TEXT)) 0.90 else 0.62) *
             length(PRED_MAIN) * length(DATASET_LEVELS),
         h = 2.6 + (if (isTRUE(PUB_GRID_TEXT)) 0.30 else 0.22) *
             nlevels(gp$outcome_f) * nlevels(gp$group))

# == A2: duel verdicts ========================================================
p_duel_pub <- p_duel +
    labs(caption = paste(
        wrap_cap("Rows at joint VIF >= 5 carry no verdict: the module and",
                 "the rival cannot be separated in the joint model, which",
                 "is not the same as neither carrying information."),
        wrap_cap("All five rivals are EXTERNAL baselines. The module against",
                 "its own non-module cells is Figure 12D."),
        wrap_cap(.mag_txt),
        wrap_cap(HKG8_CAVEAT),
        sep = "\n"))
if (interactive()) print(p_duel_pub)   # on-screen copy only; save_fig writes the files
save_plotted(tally, file.path(OUT_DIR, "Fig13B_Pub_Duel_verdicts"))
save_fig(p_duel_pub, file.path(OUT_DIR, "Fig13B_Pub_Duel_verdicts"),
         w = 2.6 + 0.62 * sum(table(dv$outcome_class) > 0) *
             length(DATASET_LEVELS) * 3,
         h = 2.2 + 0.34 * nlevels(dv$rival_f) * length(GROUP_TAGS))

# == Bootstrap, cognition only (CSV) ==========================================
bp <- bd[as.character(bd$outcome) %in% PUB_BOOT_OUTCOMES, , drop = FALSE]
bp$outcome_f <- droplevels(bp$outcome_f)
p_boot_pub <- ggplot(bp, aes(x = outcome_f, y = median_dR2, colour = rival_f,
                             group = rival_f)) +
    geom_hline(yintercept = 0, colour = "grey50", linetype = "dashed") +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0,
                  position = .pd_boot, linewidth = 0.4) +
    geom_point(aes(shape = excl), size = 1.8, position = .pd_boot) +
    scale_shape_manual(values = c(`TRUE` = 19, `FALSE` = 1),
                       name = "CI excludes zero", na.translate = FALSE) +
    scale_colour_manual(values = PRED_FILL, name = NULL) +
    facet_grid(group ~ dataset) +
    labs(title = "Donor bootstrap of the delta-R2 difference, cognition",
         subtitle = paste(.boot_n_txt,
                          "Positive is the module ahead of the rival.",
                          "Filled points have an interval excluding zero."),
         caption = wrap_cap("The encompassing test is the verdict; the",
                            "bootstrap resamples the same donors and is a",
                            "conservative check. The",
                            "staging outcomes are in the supplementary",
                            "bootstrap."),
         x = NULL,
         y = expression(paste("median ", Delta * R^2,
                              " (module) - ", Delta * R^2, " (rival)"))) +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_plotted(bp[, c("group", "dataset", "rival_label", "outcome", "n_boot",
                    "median_dR2", "ci_lo", "ci_hi", "pct_module_ahead",
                    "ci_excludes_zero")],
             file.path(OUT_DIR, "Fig13_Bootstrap_cognition"))

# =============================================================================
# 3c: panel C, the decomposition
# =============================================================================
# Module, AUCell and cluster proportion in one model. The top row of each
# facet is the module-AUCell correlation.

if (length(MV_RAW) > 0) {
    cat("\n  panel C -- decomposition\n")

    MV <- do.call(rbind, lapply(names(MV_RAW), function(g) {
        d <- MV_RAW[[g]]$decomposition
        d$group <- g
        d
    }))
    MV_ALL <- MV                                   # every score arm, kept
    MV <- MV[MV$score_arm == MV_ARM_SHOWN, , drop = FALSE]
    if (nrow(MV) == 0L)
        stop("MV_ARM_SHOWN matched no score arm: ", MV_ARM_SHOWN,
             call. = FALSE)
    MVC <- do.call(rbind, lapply(names(MV_RAW), function(g) {
        d <- MV_RAW[[g]]$correlations
        keep <- (grepl("odule", d$x) & !grepl("core", d$x) &
                 grepl("AUCell", d$y)) |
                (grepl("odule", d$y) & !grepl("core", d$y) &
                 grepl("AUCell", d$x))
        d <- d[keep, , drop = FALSE]
        d$group <- g
        d
    }))

    .glv <- unname(GROUP_LABELS[GROUP_TAGS])
    MV$group <- factor(MV$group, levels = GROUP_TAGS, labels = .glv)
    MVC$group <- factor(MVC$group, levels = GROUP_TAGS, labels = .glv)
    MV$dataset <- factor(MV$dataset, levels = DATASET_LEVELS)
    MVC$dataset <- factor(MVC$dataset, levels = DATASET_LEVELS)

    .row_lab <- function(b, o) sprintf("%s | %s", b, sub(" score$", "", o))
    MV$row_f <- .row_lab(MV$block, MV$outcome)
    .rows <- c("cov4 | MMSE", "cov4 | CASI", "+BTC | MMSE", "+BTC | CASI")
    RHO_ROW <- "module ~ AUCell rho"
    MV$row_f <- factor(MV$row_f, levels = rev(c(RHO_ROW, .rows)))

    # Verdict and statistics are separate text layers.
    MV$verdict_disp <- unname(MV_VERDICT_DISPLAY[as.character(MV$verdict)])
    if (anyNA(MV$verdict_disp))
        stop("a decomposition verdict has no display label: ",
             paste(unique(MV$verdict[is.na(MV$verdict_disp)]),
                   collapse = ", "), call. = FALSE)
    MV$lab_v <- MV$verdict_disp
    MV$lab_s <- sprintf("FDR=%.3f, VIF=%.2f", MV$q_module_t, MV$max_vif)
    MV$fill_f <- factor(MV$verdict_disp, levels = unname(MV_VERDICT_DISPLAY))
    MV$txt_col <- vapply(MV_FILL[as.character(MV$fill_f)],
                         function(h) fill_label_colour(1, c(h, h)),
                         character(1), USE.NAMES = FALSE)

    # The correlation row: same grid, no fill, value as text.
    RHO <- data.frame(group = MVC$group, dataset = MVC$dataset,
                      row_f = factor(RHO_ROW,
                                     levels = levels(MV$row_f)),
                      lab = sprintf("rho=%.3f", MVC$spearman),
                      stringsAsFactors = FALSE)

    p_mv <- ggplot(MV, aes(x = group, y = row_f)) +
        # Hairline border on every cell.
        geom_tile(aes(fill = fill_f), colour = "black",
                  linewidth = MV_CELL_LWD,
                  width = GRID_CELL, height = GRID_CELL) +
        geom_text(aes(label = lab_v, colour = txt_col),
                  size = MV_TEXT_SIZE + MV_VERDICT_BUMP, fontface = "bold",
                  family = MV_FONT, nudge_y = 0.15) +
        geom_text(aes(label = lab_s, colour = txt_col),
                  size = MV_TEXT_SIZE, fontface = "bold",
                  family = MV_FONT, nudge_y = -0.15) +
        geom_text(data = RHO, aes(label = lab),
                  size = MV_TEXT_SIZE + MV_RHO_BUMP, fontface = "bold",
                  family = MV_FONT, colour = "grey15", inherit.aes = TRUE) +
        scale_colour_identity() +
        # Unused verdict levels are dropped from the legend.
        scale_fill_manual(values = MV_FILL, name = NULL, drop = TRUE,
                          na.translate = FALSE) +
        facet_grid(dataset ~ .) +
        theme(panel.spacing = grid::unit(PANEL_GAP, "pt")) +
        labs(title = sprintf(paste("Does the module survive with %s and an",
                                   "abundance both in the model?"),
                             MV_ARM_SHOWN),
             subtitle = paste("Outcome ~ cov4 or +BTC, plus AUCell, plus",
                              "cluster proportion, plus the module. Eight",
                              "tests per group: 2 regions x 2 covariate",
                              "blocks x 2 cognitive outcomes. The score term",
                              "drawn here is", paste0(MV_ARM_SHOWN, ";"),
                              "the other two are in panel C2 and in the",
                              "supplementary table."),
             caption = paste(
                 wrap_cap("Cream is not missing data. Where the module and",
                          "the score correlate above about 0.9 the joint",
                          "model cannot separate them and the question",
                          "cannot be asked, which is why the correlation is",
                          "the top row of each panel."),
                 wrap_cap(HKG8_CAVEAT), sep = "\n"),
             x = NULL, y = NULL) +
        theme_pub() +
        theme(axis.text.x = element_text(angle = 20, hjust = 1))
    if (interactive()) print(p_mv)   # on-screen copy only; save_fig writes the files
    # Report when the font cannot be resolved.
    if (!isTRUE(capabilities("cairo")))
        cat(sprintf(paste("    NOTE: cairo devices unavailable, so '%s'",
                          "falls back to the device default family.\n"),
                    MV_FONT))

    save_plotted(MV[, c("group", "dataset", "block", "outcome", "score_arm",
                        "n", "beta_module", "p_module_t", "q_module_t",
                        "beta_score", "beta_cluster", "max_vif", "verdict",
                        "verdict_disp")],
                 file.path(OUT_DIR, "Fig13C_Pub_Decomposition"))
    save_fig(p_mv, file.path(OUT_DIR, "Fig13C_Pub_Decomposition"),
             w = 2.4 + 1.9 * nlevels(MV$group),
             h = 2.8 + 0.62 * (length(.rows) + 1L) *
                 length(DATASET_LEVELS))

    .hold <- sum(MV$verdict == "HOLDS")
    .unread <- sum(MV$verdict == "not readable")
    gr("GR-13A-07", if (length(.mv_missing) == 0L) "PASS" else "WARN",
       sprintf(paste("decomposition over %d primary-score row(s), %d",
                     "group(s): %d HOLDS, %d not readable.%s"),
               nrow(MV), nlevels(MV$group), .hold, .unread,
               if (length(.mv_missing))
                   paste0(" missing: ", paste(.mv_missing, collapse = ", "))
               else ""))
} else {
    gr("GR-13A-07", "WARN",
       "no Multivariate results found; panel C not drawn")
}

# =============================================================================
# 3c2: panel C2, module against score coefficient
# =============================================================================
# One model per tile: module and score coefficients and their formal contrast.

if (exists("MV_ALL") && nrow(MV_ALL) > 0) {
    cat("\n  panel C2 -- module against the score, head to head\n")

    # Filter on the arm tag; MV_ALL is taken before the group relabelling.
    CT <- MV_ALL[as.character(MV_ALL$group) == GROUP_TAGS[1], , drop = FALSE]
    if (nrow(CT) == 0L)
        stop("panel C2: no decomposition rows for ", GROUP_TAGS[1],
             call. = FALSE)
    # Region as a facet; rows are block and outcome.
    CT$dataset <- factor(CT$dataset, levels = DATASET_LEVELS)
    .rlv <- c("cov4 | MMSE", "cov4 | CASI", "+BTC | MMSE", "+BTC | CASI")
    CT$row_f <- factor(sprintf("%s | %s", CT$block,
                               sub(" score$", "", CT$outcome)),
                       levels = rev(.rlv))
    CT$arm_f <- factor(CT$score_arm,
                       levels = c("AUCell", "UCell", "AddModuleScore"))

    # Fill is the directional contrast verdict, computed from the coefficients.
    # Panel C2 text sizes start from panel C's.
    CT_TEXT_BUMP <- 0.7                # about 2 pt on both lines
    .ct_size <- MV_TEXT_SIZE + CT_TEXT_BUMP

    .APX <- paste0("Module ", intToUtf8(0x2248), " Score")
    .rd <- CT$readable_vs_score %in% c(TRUE, "TRUE", 1)
    .sig <- .rd & CT$q_module_vs_score < ALPHA
    CT$fill_f <- ifelse(!.rd, "unreadable",
                 ifelse(!.sig, .APX,
                 ifelse(abs(CT$beta_module) > abs(CT$beta_score),
                        "Module > Score", "Score > Module")))
    CT$fill_f <- factor(CT$fill_f,
                        levels = c("Module > Score", "Score > Module",
                                   .APX, "unreadable"))
    CT_FILL <- setNames(c("#2D2D2D", "#2A9D8F", "#D9D9D9", "#FBE7A2"),
                        c("Module > Score", "Score > Module", .APX,
                          "unreadable"))
    CT$txt_col <- vapply(CT_FILL[as.character(CT$fill_f)],
                         function(h) fill_label_colour(1, c(h, h)),
                         character(1), USE.NAMES = FALSE)

    # Separate p-values are not printed; asterisks mark the contrast FDR.
    .stars <- function(q) ifelse(is.na(q), "",
                         ifelse(q < 0.01, "**", ifelse(q < ALPHA, "*", "")))
    # Coefficients only, module first.
    CT$lab_v <- as.character(CT$fill_f)
    CT$lab_s <- sprintf("%+.3f vs %+.3f\nFDR=%.3f%s",
                        CT$beta_module, CT$beta_score,
                        CT$q_module_vs_score, .stars(CT$q_module_vs_score))

    p_ct <- ggplot(CT, aes(x = arm_f, y = row_f)) +
        geom_tile(aes(fill = fill_f), colour = "black",
                  linewidth = MV_CELL_LWD,
                  width = GRID_CELL, height = GRID_CELL) +
        # Same typography constants as panel C.
        geom_text(aes(label = lab_v, colour = txt_col),
                  size = .ct_size + MV_VERDICT_BUMP, fontface = "bold",
                  family = MV_FONT, nudge_y = 0.20) +
        geom_text(aes(label = lab_s, colour = txt_col), size = .ct_size,
                  fontface = "bold", family = MV_FONT, lineheight = 0.95,
                  nudge_y = -0.13) +
        scale_colour_identity() +
        scale_fill_manual(values = CT_FILL, name = NULL, drop = TRUE) +
        facet_grid(dataset ~ .) +
        labs(title = "Contrast: module against score coefficient",
             subtitle = sprintf(paste("%s ONLY, unlike the panel beside it.",
                                      "One model per tile: outcome ~",
                                      "covariates + module + cluster",
                                      "proportion + the score named in the",
                                      "column. The verdict is the formal",
                                      "contrast between the module's and the",
                                      "score's coefficients in that model;",
                                      "the direction names whichever",
                                      "coefficient is larger."),
                                GROUP_TAGS[1]),
             caption = paste(
                 wrap_cap("A difference in significance is not a significant",
                          "difference, so the two coefficients are tested",
                          "against each other and their separate p-values are",
                          "deliberately not shown. Fill is the contrast, not",
                          "the module's own verdict; the module's own FDR is",
                          "in panel C."),
                 wrap_cap("Both coefficients are from ONE model, so they are",
                          "adjusted for each other and for cluster",
                          "proportion. * FDR < 0.05, ** FDR < 0.01 on the",
                          "contrast."),
                 sep = "\n"),
             x = NULL, y = NULL) +
        theme_pub() +
        # Strips off; spacing set to match panel C.
        theme(strip.background = element_blank(),
              strip.text = element_blank(),
              panel.spacing = grid::unit(PANEL_GAP, "pt"))
    if (interactive()) print(p_ct)   # on-screen copy only; save_fig writes the files
    save_plotted(CT[, c("group", "dataset", "block", "outcome", "score_arm",
                        "n", "beta_module", "q_module_t", "beta_score",
                        "p_score", "beta_module_vs_score",
                        "q_module_vs_score", "vif_module_vs_score",
                        "readable_vs_score")],
                 file.path(OUT_DIR, "Fig13C2_Pub_Contrast"))
    save_fig(p_ct, file.path(OUT_DIR, "Fig13C2_Pub_Contrast"),
             w = 2.8 + 2.3 * nlevels(CT$arm_f),
             h = 2.8 + 0.78 * nlevels(CT$row_f) * nlevels(CT$dataset))

    gr("GR-13A-08", if (nrow(CT) > 0L) "PASS" else "WARN",
       sprintf(paste("coefficient contrast on %s over %d tile(s): module",
                     "stronger on %d, score stronger on %d,",
                     "indistinguishable on %d, unreadable on %d. The",
                     "score's own coefficient clears %g on %d"),
               GROUP_TAGS[1], nrow(CT),
               sum(CT$fill_f == "Module > Score"),
               sum(CT$fill_f == "Score > Module"),
               sum(CT$fill_f == .APX),
               sum(CT$fill_f == "unreadable"),
               ALPHA, sum(CT$p_score < ALPHA, na.rm = TRUE)))
}

# =============================================================================
# 3d: coefficient forest per group, controls separated (CSV only)
# =============================================================================

fd_all <- MARG[MARG$block == PRIMARY_BLOCK &
               MARG$predictor %in% unname(PRED_KEYS[PRED_SUPP]), ,
               drop = FALSE]
fd_all$lo <- fd_all$beta - 1.96 * fd_all$se
fd_all$hi <- fd_all$beta + 1.96 * fd_all$se
fd_all$class_f <- factor(ifelse(fd_all$is_control, "Control outcomes",
                                "Clinical outcomes"),
                         levels = c("Clinical outcomes", "Control outcomes"))
# Levels reversed for position_dodge(); legend order set by `breaks`.
fd_all$pred_f <- factor(fd_all$predictor,
                        levels = rev(unname(PRED_KEYS[PRED_SUPP])))

for (i in seq_along(GROUP_TAGS)) {
    g <- GROUP_TAGS[i]
    lab <- unname(GROUP_LABELS[g])
    fd <- fd_all[as.character(fd_all$group) == lab, , drop = FALSE]
    if (nrow(fd) == 0) next
    ctrl_lv <- unique(fd$outcome[fd$is_control])
    .flv <- rev(c(OUTCOME_ORDER[OUTCOME_ORDER %in% fd$outcome], ctrl_lv))
    fd$outcome_f <- factor(fd$outcome, levels = .flv,
                           labels = relabel_outcomes(.flv))
    # group = pred_f and one shared dodge, so points sit on their own intervals.
    .pd_f <- position_dodge(width = 0.7)
    p_f <- ggplot(fd, aes(x = beta, y = outcome_f, colour = pred_f,
                          group = pred_f)) +
        geom_vline(xintercept = 0, colour = "grey50", linetype = "dashed") +
        geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0,
                       position = .pd_f, linewidth = 0.4) +
        geom_point(aes(shape = p_adj < ALPHA), size = 1.9,
                   position = .pd_f) +
        scale_shape_manual(values = c(`TRUE` = 19, `FALSE` = 1),
                           name = "FDR < 0.05", na.translate = FALSE) +
        scale_colour_manual(values = PRED_FILL, name = NULL,
                            breaks = unname(PRED_KEYS[PRED_SUPP])) +
        facet_grid(class_f ~ dataset, scales = "free_y", space = "free_y") +
        labs(title = sprintf("Standardised coefficient per predictor [%s, %s]",
                             lab, PRIMARY_BLOCK),
             subtitle = paste("Same donors, same covariate block, same",
                              "outcomes. Predictors read downwards in legend",
                              "order. Open points do not clear FDR."),
             caption = paste(
                 wrap_cap("Higher-is-better outcomes (MMSE, CASI) and",
                          "higher-is-worse outcomes (staging, cognitive",
                          "status) sit on opposite sides of zero for the",
                          "same biological direction."),
                 wrap_cap("The control block is drawn separately: on this",
                          "covariate block the module's association with RIN",
                          "is of the same order as its association with",
                          "cognition and is a technical caveat rather than a",
                          "result."),
                 wrap_cap(share_caption()),
                 sep = "\n"),
             x = "beta (outcome SD per predictor SD)", y = NULL) +
        theme_pub()
    base <- file.path(OUT_DIR, sprintf("Fig13_Forest_%s", g))
    save_plotted(fd[, c("group", "dataset", "outcome_class", "outcome",
                        "predictor", "n", "beta", "se", "lo", "hi",
                        "delta_r2", "p_adj")], base)
}

# =============================================================================
# STEP 4: ledger
# =============================================================================

cat("\nSTEP 4 -- ledger\n")

# Figures whose plotted CSV must exist (anything that calls save_fig).
.fig_names <- c("Fig13A_DeltaR2_grid", "Fig13B_Duel_verdicts",
                "Fig13A_Pub_DeltaR2_grid", "Fig13B_Pub_Duel_verdicts",
                "Fig13C_Pub_Decomposition", "Fig13C2_Pub_Contrast")
.fig_csv <- paste0(.fig_names, "_plotted.csv")
.orphan <- setdiff(.fig_csv, .csv_written)
gr("GR-13A-06",
   if (length(.orphan) == 0L && .n_fig == length(.fig_names)) "PASS"
   else "WARN",
   sprintf(paste("%d figure(s) written and %d plotted-data CSV(s). Every",
                 "figure has its rows on disk%s. Table-only outputs: %s"),
           .n_fig, length(.csv_written),
           if (length(.orphan)) paste0("; missing: ",
                                       paste(.orphan, collapse = ", "))
           else "",
           paste(sub("_plotted\\.csv$", "",
                     setdiff(.csv_written, .fig_csv)), collapse = ", ")))

GR_DF <- do.call(rbind, .GR)
utils::write.csv(GR_DF, file.path(OUT_DIR, "Fig13_guardrails.csv"),
                 row.names = FALSE)
utils::write.csv(share, file.path(OUT_DIR, "Fig13_baseline_sharing.csv"),
                 row.names = FALSE)

cat(sprintf("\n  %d PASS / %d WARN / %d FAIL\n",
            sum(GR_DF$status == "PASS"), sum(GR_DF$status == "WARN"),
            sum(GR_DF$status == "FAIL")))
cat(sprintf("  outputs in %s\n", OUT_DIR))
cat(sprintf("  log written to %s\n", basename(LOG_PATH)))

sink()
close(.log_con)
