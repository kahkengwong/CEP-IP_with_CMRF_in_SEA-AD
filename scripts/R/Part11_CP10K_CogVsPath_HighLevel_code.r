# =============================================================================
# Part11_CP10K_CogVsPath_HighLevel_code.r
# SEA-AD CEP-IP/CMRF, Part 11, high-level summary panel (Figure 9E)
#   Draws the four cell populations (module or non-module, Pre-IP or Post-IP)
#   in copula space for M4 NeoC, from the Validation 1 refined-cell table, with
#   the IP line, a display spine and reference blocks for the panel questions
#   and the covariate blocks. Nothing is refitted.
#
# Usage (from the SEA-AD folder, or set SEAAD_ROOT), after Part 10e Validation1
#   source("Part11_CP10K_CogVsPath_HighLevel_code.r")
#
# Outputs (Part11_Identity/Part11_CogVsPath_HighLevel/): figure (PDF, PNG),
#   population counts and guardrails (CSV), log.
# =============================================================================
suppressPackageStartupMessages({
    library(readxl)
    library(ggplot2)
})
.HAS_MGCV <- requireNamespace("mgcv", quietly = TRUE)

# == CONFIGURATION ============================================================
CEP_HOME <- Sys.getenv("SEAAD_ROOT", getwd())
V1_ROOT  <- file.path(CEP_HOME, "Part10e_Validation1")
OUT_DIR  <- file.path(CEP_HOME, "Part11_Identity",
                      "Part11_CogVsPath_HighLevel")

# One group and region, matching the band-shape figure of Part 11 BandAblation.
FIG_GROUP <- "M4"
FIG_DS    <- "NeoC"
V1DIR     <- c(M4 = "Mic13_PVM3", PVM_3_all = "PVM_3",
               PVM_1_all = "PVM_1", HKG8 = "HKG8")
MODULE_CLASSES <- c("Original + Refined", "Recovered")

# Four population colours, matching the Part 10h DeltaR2 figure legend.
COL_NM_PRE   <- "#A8C8E0"   # Non-module Pre-IP
COL_MOD_PRE  <- "#5A7FA6"   # Module Pre-IP
COL_NM_POST  <- "#B08968"   # Non-module Post-IP
COL_MOD_POST <- "#2D2D2D"   # Module Post-IP, the predictor

IP_LINE_COL <- "#4A4A4A"
# Display spine: a GAM refit of v ~ s(u), not the CMRF spine.
SHOW_SPINE  <- TRUE
SPINE_COL   <- "#1A1A1A"
# Reference-block colours: rule lines, headings and body text.
RULE_COL    <- "#5A5A5A"
HEAD_COL    <- "#05335C"
BODY_COL    <- "#1A1A1A"

PT_SIZE  <- 0.30           # BandAblation uses 0.25 over three facets
PT_ALPHA <- 0.55
PANEL_SQUARE    <- TRUE    # force a square panel, as in the band-shape figure
SHOW_IN_CONSOLE <- TRUE    # also draw to the active device (the Plots pane)
FIG_W    <- 7.5
FIG_H    <- 8.5
TOP_FRAC <- 0.66           # share of the height given to the scatter
DPI      <- 300

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
LOG <- file.path(OUT_DIR, paste0("Part11_CogVsPath_HighLevel_log_",
                                 format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
say <- function(...) {
    m <- paste0(...)
    cat(m, "\n", sep = "")
    cat(m, "\n", sep = "", file = LOG, append = TRUE)
}
GR <- list()
gr <- function(id, status, msg) {
    GR[[length(GR) + 1L]] <<- data.frame(id = id, status = status,
                                         message = msg,
                                         stringsAsFactors = FALSE)
    say(sprintf("  [%-4s]  %-12s %s", status, id, msg))
}

say(strrep("=", 76))
say("Part 11 high-level summary panel. Started ", format(Sys.time()))
say("  Read : ", V1_ROOT)
say("  Write: ", OUT_DIR)
say(strrep("=", 76))

# =============================================================================
# STEP 1: Refined cells (same file, sheet and filter as Part 11 BandAblation)
# =============================================================================
say("")
say("STEP 1: refined cells")
dpath <- file.path(V1_ROOT, paste0("Part10e_Validation1_", V1DIR[[FIG_GROUP]]))
f <- list.files(dpath, full.names = TRUE,
                pattern = paste0("_V1_", FIG_DS,
                                 "_SK_interaction_analysis\\.xlsx$"))
if (!length(f))
    stop("No Validation 1 workbook for ", FIG_GROUP, " ", FIG_DS, " under ",
         dpath, call. = FALSE)
d <- as.data.frame(read_excel(f[1], "Step4_Refined_Cells"))
need <- c("barcode", "cell_class", "u_all", "v_all", "cDFG_expr")
miss <- setdiff(need, names(d))
if (length(miss))
    stop("Step4_Refined_Cells lacks: ", paste(miss, collapse = ", "),
         call. = FALSE)
d <- d[stats::complete.cases(d[, c("u_all", "v_all", "cDFG_expr")]) &
       !is.na(d$cell_class), , drop = FALSE]
gr("GR-11HL-01", "PASS",
   sprintf(paste("%s %s from %s: %d refined cell(s), %d column(s)"),
           FIG_GROUP, FIG_DS, basename(f[1]), nrow(d), ncol(d)))

# =============================================================================
# STEP 2: The IP and the module boundary
# =============================================================================
say("")
say("STEP 2: IP and module boundary")

# == IP on the cDFG axis, derived as in Part 11 BandAblation ================
lab <- ifelse(grepl("Post-IP", d$group_refined, fixed = TRUE), "Post",
       ifelse(grepl("Pre-IP",  d$group_refined, fixed = TRUE), "Pre",
       ifelse(grepl("Post-IP", d$IP_status, fixed = TRUE), "Post",
       ifelse(grepl("Pre-IP",  d$IP_status, fixed = TRUE), "Pre",
              NA_character_))))
hi <- suppressWarnings(max(d$cDFG_expr[!is.na(lab) & lab == "Pre"]))
lo <- suppressWarnings(min(d$cDFG_expr[!is.na(lab) & lab == "Post"]))
ip_y <- if (is.finite(hi) && is.finite(lo) && hi < lo) (hi + lo) / 2 else
        stats::median(d$cDFG_expr)
d$side <- ifelse(d$cDFG_expr < ip_y, "Pre", "Post")

# The horizontal IP line is valid only if v is monotone in cDFG; checked here.
.rho <- suppressWarnings(stats::cor(d$cDFG_expr, d$v_all, method = "spearman"))
.o <- order(d$cDFG_expr)
v_ip <- as.numeric(stats::approx(d$cDFG_expr[.o], d$v_all[.o], xout = ip_y,
                                 rule = 2, ties = "ordered")$y)
.agree <- mean((d$v_all < v_ip) == (d$cDFG_expr < ip_y))
gr("GR-11HL-02", if (is.finite(.rho) && .rho > 0.999 && .agree > 0.999)
                     "PASS" else "WARN",
   sprintf(paste("IP_y = %.4f maps to v = %.4f; Spearman(cDFG, v) = %.6f;",
                 "v split matches the cDFG split on %.2f%% of cells"),
           ip_y, v_ip, .rho, 100 * .agree))

# == Module boundary from the CMRF classification =========================
d$in_module <- d$cell_class %in% MODULE_CLASSES
d$pop <- factor(
    ifelse( d$in_module & d$side == "Post", "Module Post-IP",
    ifelse( d$in_module & d$side == "Pre",  "Module Pre-IP",
    ifelse(!d$in_module & d$side == "Post", "Non-module Post-IP",
                                            "Non-module Pre-IP"))),
    levels = c("Non-module Pre-IP", "Module Pre-IP",
               "Non-module Post-IP", "Module Post-IP"))
.tab <- table(d$pop)
gr("GR-11HL-03", if (sum(.tab) == nrow(d) && !anyNA(d$pop)) "PASS" else "FAIL",
   sprintf(paste("four populations over %d cell(s): %s"),
           nrow(d), paste(sprintf("%s %d", names(.tab), as.integer(.tab)),
                          collapse = "; ")))

POP_COL <- c(`Non-module Pre-IP`  = COL_NM_PRE,
             `Module Pre-IP`      = COL_MOD_PRE,
             `Non-module Post-IP` = COL_NM_POST,
             `Module Post-IP`     = COL_MOD_POST)
gr("GR-11HL-04",
   if (setequal(names(POP_COL), levels(d$pop))) "PASS" else "FAIL",
   sprintf(paste("colours drawn: %s"),
           paste(sprintf("%s %s", names(POP_COL), unname(POP_COL)),
                 collapse = "; ")))

# Interleaving: per u-slice, non-module cells whose v lies within the module's
# v range in that slice. Recorded only; nothing in the panel uses it.
interleave <- function(um, vm, un, vn, nb = 40L) {
    if (length(um) < 30L || !length(un)) return(NULL)
    br <- unique(stats::quantile(um, probs = seq(0, 1, length.out = nb + 1L),
                                 na.rm = TRUE))
    if (length(br) < 4L) return(NULL)
    bm <- cut(um, breaks = br, include.lowest = TRUE)
    bn <- cut(un, breaks = br, include.lowest = TRUE)
    lo <- tapply(vm, bm, min); hi <- tapply(vm, bm, max)
    inr <- rep(NA, length(un))
    k <- !is.na(bn)
    inr[k] <- vn[k] >= lo[bn[k]] & vn[k] <= hi[bn[k]]
    list(n_nonmodule = length(un),
         n_in_slice = sum(!is.na(bn)),
         n_interleaved = sum(inr, na.rm = TRUE),
         frac = mean(inr, na.rm = TRUE),
         median_width = stats::median(hi - lo, na.rm = TRUE),
         bins = sum(!is.na(lo)))
}
kp <- d$in_module & d$side == "Post"   # module Post-IP cells
kr <- d$in_module & d$side == "Pre"    # module Pre-IP cells
.ip <- interleave(d$u_all[kp], d$v_all[kp],
                  d$u_all[!d$in_module & d$side == "Post"],
                  d$v_all[!d$in_module & d$side == "Post"])
.ir <- interleave(d$u_all[kr], d$v_all[kr],
                  d$u_all[!d$in_module & d$side == "Pre"],
                  d$v_all[!d$in_module & d$side == "Pre"])
gr("GR-11HL-10", if (isTRUE(.ip$frac < 0.01) && isTRUE(.ir$frac < 0.01))
                     "PASS" else "WARN",
   sprintf(paste("interleaving (recorded only). Post-IP: %d of %d",
                 "non-module cells (%.1f%%) within the module v range, %d",
                 "slices, median height %.3f. Pre-IP: %d of %d (%.1f%%), %d",
                 "slices, median height %.3f"),
           .ip$n_interleaved, .ip$n_in_slice, 100 * .ip$frac, .ip$bins,
           .ip$median_width,
           .ir$n_interleaved, .ir$n_in_slice, 100 * .ir$frac, .ir$bins,
           .ir$median_width))

# == Fitted spine, a visual guide only ====================================
SP <- NULL
if (SHOW_SPINE && .HAS_MGCV) {
    fd <- data.frame(u = d$u_all, v = d$v_all)
    m1 <- mgcv::gam(v ~ s(u, bs = "tp"), data = fd, method = "REML")
    ug <- seq(min(fd$u), max(fd$u), length.out = 300)
    SP <- data.frame(u = ug,
                     fit = as.numeric(stats::predict(m1, data.frame(u = ug))))
}

# =============================================================================
# STEP 3: Figure
# =============================================================================
say("")
say("STEP 3: figure")

p_top <- ggplot() +
    geom_point(data = d, aes(x = u_all, y = v_all, colour = pop),
               size = PT_SIZE, alpha = PT_ALPHA) +
    geom_hline(yintercept = v_ip, linetype = "dashed",
               colour = IP_LINE_COL, linewidth = 0.35)
if (!is.null(SP))
    p_top <- p_top + geom_line(data = SP, aes(x = u, y = fit),
                               colour = SPINE_COL, linewidth = 0.35)
p_top <- p_top +
    annotate("text", x = 0.015, y = v_ip, label = "inflection point (IP)",
             hjust = 0, vjust = -0.6, size = 2.7, colour = IP_LINE_COL) +
    scale_colour_manual(values = POP_COL, name = NULL,
                        guide = guide_legend(nrow = 2, byrow = TRUE,
                            override.aes = list(size = 2, alpha = 1))) +
    scale_x_continuous(limits = c(0, 1), expand = expansion(mult = 0.02)) +
    scale_y_continuous(expand = expansion(mult = 0.03)) +
    labs(x = "u  (target gene, rank scale)",
         y = "v  (module composite, rank scale)",
         title = "(E) The four cell populations",
         subtitle = paste("One cell per point,", FIG_GROUP, FIG_DS,
                          "\nColour is the CMRF band crossed with the IP,",
                          "the same four colours as (A); the dashed line is",
                          "the IP. Each donor's value in (A) to (D) is the",
                          "percentage of their cells in one of these four",
                          "populations.")) +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", size = 11),
          plot.subtitle = element_text(size = 7.5, colour = "grey30"),
          # Legend at the bottom; the right margin stays clear for braces.
          legend.position = "bottom",
          legend.margin = margin(0, 0, 0, 0),
          legend.text = element_text(size = 7.5),
          plot.margin = margin(4, 6, 2, 4))
if (PANEL_SQUARE) p_top <- p_top + theme(aspect.ratio = 1)

# == Reference blocks: panel questions and covariate blocks ===============
BLK <- rbind(
  data.frame(row = 0, key = "What each panel asks", txt = "", head = TRUE,
             stringsAsFactors = FALSE),
  data.frame(row = 1, key = "(A)",
             txt = paste("Does the fraction predict cognition?  Yes, and",
                         "most of it sits past the IP."),
             head = FALSE),
  data.frame(row = 2, key = "(B)",
             txt = paste("Is it only plaque and tangle burden?  No. It",
                         "survives; the control does not."),
             head = FALSE),
  data.frame(row = 3, key = "(C)",
             txt = paste("Does the link run both ways?  No. Pathology drops",
                         "out once cognition is in."),
             head = FALSE),
  data.frame(row = 4, key = "(D)",
             txt = paste("Is it the module cells specifically?  On cognition",
                         "yes, on pathology no."),
             head = FALSE),
  data.frame(row = 5.6, key = "The covariate blocks", txt = "", head = TRUE),
  data.frame(row = 6.6, key = "cov4",
             txt = "Age at death, sex, APOE e4 dose, PMI. The primary model.",
             head = FALSE),
  data.frame(row = 7.6, key = "+ADNC",
             txt = "cov4 plus the ADNC composite.",
             head = FALSE),
  data.frame(row = 8.6, key = "+BTC",
             txt = "cov4 plus Braak, Thal and CERAD. The headline block.",
             head = FALSE),
  data.frame(row = 9.6, key = "+ALL",
             txt = "+BTC plus ADNC. Collinear; shown as sensitivity only.",
             head = FALSE))
# One column at full width, so a question and its verdict fit on one line.
BLK$x  <- 0.01
BLK$xt <- ifelse(BLK$head, 0.01, 0.075)
BLK$y  <- 10.2 - BLK$row

p_bot <- ggplot(BLK) +
    geom_segment(data = data.frame(x = 0.01, xe = 0.99,
                                   y = BLK$y[BLK$head] - 0.34),
                 aes(x = x, xend = xe, y = y, yend = y),
                 colour = RULE_COL, linewidth = 0.4, inherit.aes = FALSE) +
    geom_text(data = BLK[BLK$head, ],
              aes(x = x, y = y, label = key), hjust = 0, size = 3.0,
              fontface = "bold", colour = HEAD_COL) +
    geom_text(data = BLK[!BLK$head, ],
              aes(x = x, y = y, label = key), hjust = 0, size = 2.8,
              fontface = "bold", colour = BODY_COL) +
    geom_text(data = BLK[!BLK$head, ],
              aes(x = xt, y = y, label = txt), hjust = 0, size = 2.7,
              colour = BODY_COL) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0.2, 10.6), expand = FALSE) +
    theme_void() +
    theme(plot.margin = margin(2, 8, 4, 8))

draw <- function() {
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(
        2, 1, heights = grid::unit(c(TOP_FRAC, 1 - TOP_FRAC), "npc"))))
    print(p_top, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
    print(p_bot, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
    grid::popViewport()
}

# Draw once to the open device as well as to file.
if (SHOW_IN_CONSOLE) try(draw(), silent = TRUE)

base <- file.path(OUT_DIR, "Part11_CogVsPath_HighLevel")
.ok <- 0L
.try <- function(open_fn, path) tryCatch({
    open_fn(); on.exit(grDevices::dev.off(), add = TRUE); draw(); TRUE
}, error = function(e) {
    say("    not written (", basename(path), "): ", conditionMessage(e))
    FALSE })
if (.try(function() grDevices::pdf(paste0(base, ".pdf"), width = FIG_W,
                                   height = FIG_H), paste0(base, ".pdf")))
    .ok <- .ok + 1L
if (.try(function() grDevices::png(paste0(base, ".png"),
                                   width = FIG_W * DPI, height = FIG_H * DPI,
                                   res = DPI), paste0(base, ".png")))
    .ok <- .ok + 1L
gr("GR-11HL-05", if (.ok == 2L) "PASS" else "WARN",
   sprintf("%d of 2 device(s) written (PDF and PNG at %d dpi, %.1f x %.1f in)",
           .ok, DPI, FIG_W, FIG_H))

# =============================================================================
# STEP 4: Record (four-way population counts and guardrails)
# =============================================================================
POPS <- data.frame(population = names(.tab), n_cells = as.integer(.tab),
                   pct_of_cells = round(100 * as.numeric(.tab) / nrow(d), 2),
                   band = ifelse(grepl("^Module", names(.tab)),
                                 "module", "non-module"),
                   side = ifelse(grepl("Post-IP", names(.tab), fixed = TRUE),
                                 "Post-IP", "Pre-IP"),
                   stringsAsFactors = FALSE)
utils::write.csv(POPS, file.path(OUT_DIR,
                 "Part11_CogVsPath_HighLevel_populations.csv"),
                 row.names = FALSE)
utils::write.csv(do.call(rbind, GR), file.path(OUT_DIR,
                 "Part11_CogVsPath_HighLevel_guardrails.csv"),
                 row.names = FALSE)
gr("GR-11HL-06", "PASS",
   sprintf(paste("population counts and guardrail table written")))

say("")
say(strrep("=", 76))
say(sprintf("  Total: %d | PASS %d | WARN %d | FAIL %d",
            length(GR),
            sum(vapply(GR, function(x) x$status == "PASS", TRUE)),
            sum(vapply(GR, function(x) x$status == "WARN", TRUE)),
            sum(vapply(GR, function(x) x$status == "FAIL", TRUE))))
say("  Output: ", OUT_DIR)
for (fn in sort(list.files(OUT_DIR)))
    say(sprintf("    %-52s %8.1f KB", fn,
                file.size(file.path(OUT_DIR, fn)) / 1e3))
say(strrep("=", 76))
