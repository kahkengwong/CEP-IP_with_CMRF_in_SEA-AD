# =============================================================================
# Part10g_Step7to9_rebuild_code.r
# SEA-AD CEP-IP/CMRF, Part 10g, Steps 7b, 7c, 8 and 9 (targeted tests, control
#   overlap and ambient check)
#   Rebuilds Part10g_Step7b_targeted_tests.csv (and its _verdict version),
#   Part10g_Step7c_PvP_targeted_tests.csv, Part10g_Step8_control_overlap.csv
#   and Part10g_Step9_ambient_check.csv from the Part 10g Part 1 RDS files and
#   the Validation 1 workbooks, and compares each with the existing CSV if one
#   is present in Part10g_DEGs/.
#   Test settings: Welch t-test; Mann-Whitney U with exact = FALSE; Step 9 uses
#   avg_log2FC. Set CHECK_ALL_SETTINGS <- TRUE to rerun the search over the
#   candidate settings.
#
# Usage (from the SEA-AD folder, or set SEAAD_ROOT), after Part 10g Part 1
#   and Validation 1
#   source("Part10g_Step7to9_rebuild_code.r")
#
# Outputs (Part10g_DEGs/Step7to9_rebuild/): the five CSVs and log.
# =============================================================================

if (!exists("CEP_HOME"))
    CEP_HOME <- Sys.getenv("SEAAD_ROOT", getwd())
stopifnot(dir.exists(CEP_HOME))
if (!exists("CHECK_ALL_SETTINGS")) CHECK_ALL_SETTINGS <- FALSE
if (!requireNamespace("readxl", quietly = TRUE))
    stop("Package readxl is required.", call. = FALSE)

DEG_DIR <- file.path(CEP_HOME, "Part10g_DEGs")
V1_DIR  <- file.path(CEP_HOME, "Part10e_Validation1")
OUT_DIR <- file.path(DEG_DIR, "Step7to9_rebuild")
dir.create(OUT_DIR, showWarnings = FALSE)

.ts      <- format(Sys.time(), "%Y%m%d_%H%M%S")
.log_txt <- file.path(OUT_DIR, paste0("Part10g_Step7to9_rebuild_log_", .ts, ".txt"))
.old_opt <- options(width = 200, warn = 1, stringsAsFactors = FALSE)
.log_con <- file(.log_txt, open = "wt")
sink(.log_con, split = TRUE)

rebuild_7to9 <- function() {

    # == Settings =============================================================
    ARMS     <- c(M4 = "Mic13_PVM3", PVM_3 = "PVM_3", PVM_1 = "PVM_1",
                  HKG8 = "HKG8")              # label = folder subdir
    TEST_ARMS <- c("PVM_3", "M4", "PVM_1")    # each tested against HKG8
    CTRL_ARMS <- c("M4", "PVM_1")             # used for pass and verdict
    DATASETS <- c("NeoC", "NonNeoC")
    C7B      <- c("PostIP_vs_NMPost", "PreIP_vs_NMPre")
    C7C      <- "PostIP_vs_PreIP"
    CIRC_RE  <- "^MT-|^MTRNR2L"
    DEFINING <- c("MT-CO3", "MT-ATP6", "MT-CO2", "MT-ND3", "MT-ND4", "MT-CYB",
                  "MT-ND1", "MT-ND2", "MT-CO1")
    ALPHA    <- 0.05
    SW_ALPHA <- 0.05
    MODULE_CLASSES <- c("Original + Refined", "Recovered")
    TOL_RND  <- 1.5e-4   # for values stored rounded to 4 decimals
    TOL_P    <- 1e-6     # relative, for p-values and FDR
    # The _verdict CSV stores rounded p-values and FDR, hence a looser tolerance.
    TOL_P_VERDICT <- 5e-3

    F7B  <- file.path(DEG_DIR, "Part10g_Step7b_targeted_tests.csv")
    F7BV <- file.path(DEG_DIR, "Part10g_Step7b_targeted_tests_verdict.csv")
    F7C  <- file.path(DEG_DIR, "Part10g_Step7c_PvP_targeted_tests.csv")
    F8   <- file.path(DEG_DIR, "Part10g_Step8_control_overlap.csv")
    F9   <- file.path(DEG_DIR, "Part10g_Step9_ambient_check.csv")

    hdr <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n",
                           strrep("=", 78), "\n", sep = "")
    as_lgl <- function(x) toupper(trimws(as.character(x))) %in% "TRUE"
    num <- function(x) suppressWarnings(as.numeric(x))

    hdr("Part 10g Steps 7b, 7c, 8 and 9 rebuild")
    cat("Date     :", format(Sys.time()), "\n")
    cat("CEP_HOME :", CEP_HOME, "\n")
    cat("R        :", R.version.string, "\n")
    cat("Output   :", OUT_DIR, "\n")

    # =========================================================================
    # 1. Inputs
    # =========================================================================
    hdr("1. INPUTS")
    RDS <- list()
    for (lab in names(ARMS)) for (ds in DATASETS) {
        d1 <- file.path(DEG_DIR, paste0("Part10g_DEGs_", ARMS[[lab]]))
        f  <- list.files(d1, pattern = paste0("_", ds, "_results\\.rds$"),
                         full.names = TRUE)
        if (length(f) != 1L)
            stop(sprintf("%s %s: %d Part 1 RDS found in %s", lab, ds,
                         length(f), d1), call. = FALSE)
        RDS[[paste(lab, ds)]] <- readRDS(f)
        cat(sprintf("  %-13s %s (%s)\n", paste(lab, ds), basename(f),
                    format(file.info(f)$mtime, "%Y-%m-%d %H:%M")))
    }
    o1 <- RDS[[1]][["res"]][[C7C]]
    cat("\n  Structure check on", names(RDS)[1], C7C, "\n")
    cat("    donor columns :", paste(names(o1[["donor"]]), collapse = ", "), "\n")
    cat("    donor_fc      :", class(o1[["donor_fc"]])[1],
        paste(dim(o1[["donor_fc"]]), collapse = " x "), "\n")

    get_res <- function(arm, ds, tg) RDS[[paste(arm, ds)]][["res"]][[tg]]

    # =========================================================================
    # 2. Steps 7b and 7c
    # =========================================================================
    sw_ok <- function(x) {
        x <- x[is.finite(x)]
        if (length(x) < 3L || length(x) > 5000L || length(unique(x)) < 3L)
            return(FALSE)
        p <- tryCatch(stats::shapiro.test(x)$p.value, error = function(e) NA)
        is.finite(p) && p > SW_ALPHA
    }
    # One arm against HKG8 under one setting. Returns test name, p2, p1.
    test_one <- function(x, y, up, welch, mw_exact) {
        x <- x[is.finite(x)]; y <- y[is.finite(y)]
        if (length(x) < 2L || length(y) < 2L)
            return(list(test = NA_character_, p2 = NA_real_, p1 = NA_real_))
        alt1 <- if (up) "greater" else "less"
        if (sw_ok(x) && sw_ok(y)) {
            p2 <- tryCatch(stats::t.test(x, y, var.equal = !welch)$p.value,
                           error = function(e) NA_real_)
            p1 <- tryCatch(stats::t.test(x, y, var.equal = !welch,
                                         alternative = alt1)$p.value,
                           error = function(e) NA_real_)
            return(list(test = "t-test", p2 = p2, p1 = p1))
        }
        ex <- if (identical(mw_exact, "default")) NULL else FALSE
        p2 <- suppressWarnings(stats::wilcox.test(x, y, exact = ex)$p.value)
        p1 <- suppressWarnings(stats::wilcox.test(x, y, exact = ex,
                                                  alternative = alt1)$p.value)
        list(test = "Mann-Whitney U", p2 = p2, p1 = p1)
    }

    build_7 <- function(contrasts, welch, mw_exact) {
        rows <- list()
        for (tg in contrasts) for (ds in DATASETS) {
            r3 <- get_res("PVM_3", ds, tg); d3 <- r3[["donor"]]
            sig <- d3$gene[as_lgl(d3$significant)]
            fcs <- lapply(c(TEST_ARMS, "HKG8"),
                          function(a) get_res(a, ds, tg)[["donor_fc"]])
            names(fcs) <- c(TEST_ARMS, "HKG8")
            dfr <- lapply(c(TEST_ARMS, "HKG8"),
                          function(a) get_res(a, ds, tg)[["donor"]])
            names(dfr) <- c(TEST_ARMS, "HKG8")
            keep <- sig[vapply(sig, function(g) all(vapply(names(fcs), function(a)
                g %in% rownames(fcs[[a]]) && g %in% dfr[[a]]$gene, logical(1))),
                logical(1))]
            hk <- dfr[["HKG8"]]
            for (g in keep) {
                lfc <- num(d3$avg_log2FC[match(g, d3$gene)])
                up  <- lfc > 0
                y   <- fcs[["HKG8"]][g, ]
                r <- list(gene = g, comparison = tg, dataset = ds,
                          direction_on_PVM_3 = if (up) "up" else "down",
                          circular = grepl(CIRC_RE, g),
                          also_called_by_HKG8 = as_lgl(hk$significant[match(g, hk$gene)]),
                          n_HKG8 = sum(is.finite(y)),
                          med_HKG8 = round(stats::median(y, na.rm = TRUE), 4))
                for (a in TEST_ARMS) {
                    x  <- fcs[[a]][g, ]
                    tt <- test_one(x, y, up, welch, mw_exact)
                    mx <- stats::median(x, na.rm = TRUE)
                    my <- stats::median(y, na.rm = TRUE)
                    r[[paste0("n_", a)]]          <- sum(is.finite(x))
                    r[[paste0("med_", a)]]        <- round(mx, 4)
                    r[[paste0("test_", a)]]       <- tt$test
                    r[[paste0("p2_", a)]]         <- tt$p2
                    r[[paste0("p1_", a)]]         <- tt$p1
                    r[[paste0("concordant_", a)]] <- if (up) mx > my else mx < my
                }
                rows[[length(rows) + 1L]] <- as.data.frame(r, check.names = FALSE)
            }
        }
        D <- do.call(rbind, rows)
        grp <- paste(D$comparison, D$dataset)
        nc  <- !D$circular
        for (a in TEST_ARMS) for (s in c("2", "1")) {
            p <- D[[paste0("p", s, "_", a)]]
            D[[paste0("FDR", s, "_", a)]] <- ave(p, grp,
                FUN = function(v) stats::p.adjust(v, "BH"))
            v <- rep(NA_real_, nrow(D))
            v[nc] <- ave(p[nc], grp[nc], FUN = function(v) stats::p.adjust(v, "BH"))
            D[[paste0("FDR", s, "_nc_", a)]] <- v
        }
        pass <- function(s) Reduce(`&`, lapply(CTRL_ARMS, function(a)
            D[[paste0("FDR", s, "_", a)]] < ALPHA & D[[paste0("concordant_", a)]]))
        D$pass_two_sided <- pass("2") %in% TRUE
        D$pass_one_sided <- pass("1") %in% TRUE
        ord <- c("gene", "comparison", "dataset", "direction_on_PVM_3",
                 "circular", "also_called_by_HKG8", "n_HKG8", "med_HKG8",
                 unlist(lapply(TEST_ARMS, function(a) paste0(c("n_", "med_",
                     "test_", "p2_", "p1_", "concordant_"), a))),
                 unlist(lapply(TEST_ARMS, function(a) paste0(c("FDR2_", "FDR1_",
                     "FDR2_nc_", "FDR1_nc_"), a))),
                 "pass_two_sided", "pass_one_sided")
        D[, ord]
    }

    add_verdict <- function(D) {
        lab <- function(a) {
            s <- D[[paste0("FDR2_", a)]] < ALPHA; s[is.na(s)] <- FALSE
            c_ <- D[[paste0("concordant_", a)]]
            ifelse(!s, paste(a, "ns"),
                   ifelse(c_, paste(a, "sig, right way"), paste(a, "sig, WRONG way")))
        }
        sg <- sapply(CTRL_ARMS, function(a) {
            s <- D[[paste0("FDR2_", a)]] < ALPHA; s[is.na(s)] <- FALSE; s })
        rw <- sapply(CTRL_ARMS, function(a) D[[paste0("concordant_", a)]])
        n_right <- rowSums(sg & rw); n_wrong <- rowSums(sg & !rw)
        D$verdict <- ifelse(n_wrong > 0, "wrong direction",
                     ifelse(n_right == length(CTRL_ARMS), "pass",
                     ifelse(n_right > 0, "one arm only", "not significant")))
        D$verdict_detail <- paste(lab("M4"), "|", lab("PVM_1"))
        D
    }

    # == Comparison of a rebuilt table with the existing CSV ===================
    # An existing CSV is read only if present; without it the comparison is skipped.
    read_old <- function(f) if (file.exists(f)) read.csv(f, stringsAsFactors = FALSE) else NULL
    compare <- function(new, old, keys, label, quiet = FALSE, tol_p = TOL_P) {
        if (is.null(old)) {
            if (!quiet) cat(sprintf("\n  %s: no existing CSV; comparison skipped\n", label))
            return(invisible(0))
        }
        kn <- do.call(paste, c(new[keys], sep = "|"))
        ko <- do.call(paste, c(old[keys], sep = "|"))
        only_new <- setdiff(kn, ko); only_old <- setdiff(ko, kn)
        m <- match(ko, kn); have <- !is.na(m)
        n_bad_total <- length(only_new) + length(only_old)
        if (!quiet) {
            cat(sprintf("\n  %s: existing %d rows | rebuilt %d rows | in existing only %d | in rebuilt only %d\n",
                        label, nrow(old), nrow(new), length(only_old), length(only_new)))
            if (length(only_old)) cat("    existing only:", head(only_old, 8), "\n")
            if (length(only_new)) cat("    rebuilt only :", head(only_new, 8), "\n")
        }
        for (cn in setdiff(intersect(names(old), names(new)), keys)) {
            a <- old[[cn]][have]; b <- new[[cn]][m[have]]
            if (is.numeric(a) || is.numeric(b)) {
                a <- num(a); b <- num(b)
                if (grepl("^p[12]_|^FDR", cn)) {
                    bad <- !((is.na(a) & is.na(b)) |
                             (abs(a - b) <= tol_p * pmax(abs(a), 1e-300)))
                } else {
                    bad <- !((is.na(a) & is.na(b)) | (abs(a - b) <= TOL_RND))
                }
            } else {
                bad <- !(as.character(a) == as.character(b) |
                         (is.na(a) & is.na(b)))
                if (is.logical(old[[cn]]) || is.logical(new[[cn]]))
                    bad <- as_lgl(a) != as_lgl(b)
            }
            bad[is.na(bad)] <- TRUE
            n_bad_total <- n_bad_total + sum(bad)
            if (!quiet && any(bad)) {
                cat(sprintf("    %-22s %4d of %d differ", cn, sum(bad), length(bad)))
                k <- which(bad)[1:min(3, sum(bad))]
                cat("  e.g.", paste(sprintf("%s: %s vs %s", ko[have][k],
                    format(a[k], digits = 6), format(b[k], digits = 6)),
                    collapse = "; "), "\n")
            }
        }
        if (!quiet) cat(sprintf("    VERDICT: %s\n", if (n_bad_total == 0L)
            "MATCH" else paste(n_bad_total, "differences")))
        invisible(n_bad_total)
    }

    hdr("2. STEPS 7b AND 7c: choose the test settings")
    old7b <- read_old(F7B)
    old7c <- read_old(F7C)
    K7 <- c("gene", "comparison", "dataset")
    grid <- if (isTRUE(CHECK_ALL_SETTINGS))
                expand.grid(welch = c(TRUE, FALSE), mw_exact = c("default", "FALSE"),
                            stringsAsFactors = FALSE) else
                data.frame(welch = TRUE, mw_exact = "FALSE", stringsAsFactors = FALSE)
    score <- numeric(nrow(grid)); built <- list()
    for (i in seq_len(nrow(grid))) {
        b7 <- build_7(c(C7B, C7C), grid$welch[i], grid$mw_exact[i])
        built[[i]] <- b7
        score[i] <- compare(b7[b7$comparison %in% C7B, ], old7b, K7, "", TRUE) +
                    compare(b7[b7$comparison %in% C7C, ], old7c, K7, "", TRUE)
        cat(sprintf("  t-test %-7s | Mann-Whitney exact = %-7s | differences %d\n",
                    if (grid$welch[i]) "Welch" else "Student",
                    grid$mw_exact[i], score[i]))
    }
    best <- which.min(score)
    cat(sprintf("\n  Using: t-test %s, Mann-Whitney exact = %s\n",
                if (grid$welch[best]) "Welch" else "Student", grid$mw_exact[best]))
    b7  <- built[[best]]
    n7b <- b7[b7$comparison %in% C7B, ]; rownames(n7b) <- NULL
    n7c <- b7[b7$comparison %in% C7C, ]; rownames(n7c) <- NULL
    n7v <- add_verdict(n7b)

    hdr("3. STEPS 7b AND 7c: comparison with the existing CSVs")
    compare(n7b, old7b, K7, "Step 7b")
    compare(n7v, read_old(F7BV), K7, "Step 7b verdict",
            tol_p = TOL_P_VERDICT)
    compare(n7c, old7c, K7, "Step 7c")
    write.csv(n7b, file.path(OUT_DIR, basename(F7B)), row.names = FALSE)
    write.csv(n7v, file.path(OUT_DIR, basename(F7BV)), row.names = FALSE)
    write.csv(n7c, file.path(OUT_DIR, basename(F7C)), row.names = FALSE)

    # =========================================================================
    # 4. Step 8
    # =========================================================================
    hdr("4. STEP 8: Module membership overlap between arms")
    cells <- list()
    for (lab in names(ARMS)) for (ds in DATASETS) {
        dv <- file.path(V1_DIR, paste0("Part10e_Validation1_", ARMS[[lab]]))
        fv <- list.files(dv, pattern = paste0("_V1_", ds,
                                              "_SK_interaction_analysis\\.xlsx$"),
                         full.names = TRUE)
        if (length(fv) != 1L)
            stop(sprintf("%s %s: %d Validation 1 workbooks in %s", lab, ds,
                         length(fv), dv), call. = FALSE)
        rc <- as.data.frame(readxl::read_excel(fv, sheet = "Step4_Refined_Cells"))
        # Cell ID column: a named one if present, otherwise the first column
        # whose values are unique text.
        pref <- c("barcode", "cell_id", "exact_id", "cell", "cell_name",
                  "sample_id", "obs_names")
        idc <- intersect(pref, names(rc))
        if (!length(idc))
            idc <- names(rc)[vapply(rc, function(v) is.character(v) &&
                                    !anyDuplicated(v), logical(1))]
        if (!length(idc))
            stop(sprintf("%s %s: no cell ID column. Columns: %s", lab, ds,
                         paste(names(rc), collapse = ", ")), call. = FALSE)
        cells[[paste(lab, ds)]] <- data.frame(
            id = as.character(rc[[idc[1]]]),
            module = rc$cell_class %in% MODULE_CLASSES)
        cat(sprintf("  %-13s %6d cells | ID column %-12s | Module %6d\n",
                    paste(lab, ds), nrow(rc), idc[1],
                    sum(rc$cell_class %in% MODULE_CLASSES)))
    }
    pairs <- list(c("HKG8", "PVM_3"), c("HKG8", "M4"), c("M4", "PVM_3"))
    rows <- list()
    for (ds in DATASETS) for (pr in pairs) {
        A <- cells[[paste(pr[1], ds)]]; B <- cells[[paste(pr[2], ds)]]
        sh <- intersect(A$id, B$id)
        a <- A$module[match(sh, A$id)]; b <- B$module[match(sh, B$id)]
        both <- sum(a & b); ao <- sum(a & !b); bo <- sum(!a & b)
        nei <- sum(!a & !b); N <- length(sh)
        po <- (both + nei) / N
        pe <- ((both + ao) * (both + bo) + (bo + nei) * (ao + nei)) / N^2
        rows[[length(rows) + 1L]] <- data.frame(
            dataset = ds, arm_a = pr[1], arm_b = pr[2],
            n_a = nrow(A), n_b = nrow(B), n_shared = N,
            both_module = both, a_only_module = ao, b_only_module = bo,
            neither_module = nei,
            jaccard_module = round(both / (both + ao + bo), 4),
            kappa = round((po - pe) / (1 - pe), 4),
            pct_a_module_that_is_b_module = round(100 * both / (both + ao), 1))
    }
    n8 <- do.call(rbind, rows)
    print(n8, row.names = FALSE)
    compare(n8, read_old(F8),
            c("dataset", "arm_a", "arm_b"), "Step 8")
    write.csv(n8, file.path(OUT_DIR, basename(F8)), row.names = FALSE)

    # =========================================================================
    # 5. Step 9
    # =========================================================================
    hdr("5. STEP 9: axis genes on non-Module cells")
    old9 <- read_old(F9)
    hits <- n7c[!n7c$circular & n7c$pass_two_sided, c("gene", "dataset")]
    cat("  Non-circular Step 7c passes (two-sided):",
        paste(paste(hits$gene, hits$dataset), collapse = ", "), "\n")
    sel <- rbind(data.frame(hits, set = "neuronal hit from Step 7c"),
                 do.call(rbind, lapply(DATASETS, function(ds)
                     data.frame(gene = DEFINING, dataset = ds,
                                set = "defining (circular)"))))
    sel <- sel[order(match(sel$dataset, DATASETS)), ]
    # Candidate sources for mod_ and non_: each numeric column of the donor
    # frame, and the median of donor_fc.
    value_of <- function(arm, ds, tg, g, src) {
        r <- get_res(arm, ds, tg)
        if (src == "median_donor_fc") {
            fc <- r[["donor_fc"]]
            if (!g %in% rownames(fc)) return(NA_real_)
            return(round(stats::median(fc[g, ], na.rm = TRUE), 4))
        }
        d <- r[["donor"]]; k <- match(g, d$gene)
        if (is.na(k)) NA_real_ else round(num(d[[src]][k]), 4)
    }
    d0 <- get_res("PVM_3", "NeoC", C7C)[["donor"]]
    srcs <- if (isTRUE(CHECK_ALL_SETTINGS))
                c(names(d0)[vapply(d0, is.numeric, logical(1))], "median_donor_fc") else
                "avg_log2FC"
    build_9 <- function(src) {
        out <- sel
        for (a in c("PVM_3", "M4", "PVM_1", "HKG8")) {
            out[[paste0("mod_", a)]] <- mapply(function(g, ds)
                value_of(a, ds, C7C, g, src), out$gene, out$dataset)
            out[[paste0("non_", a)]] <- mapply(function(g, ds)
                value_of(a, ds, "NMPost_vs_NMPre", g, src), out$gene, out$dataset)
            out[[paste0("non_sig_", a)]] <- mapply(function(g, ds) {
                d <- get_res(a, ds, "NMPost_vs_NMPre")[["donor"]]
                as_lgl(d$significant[match(g, d$gene)]) }, out$gene, out$dataset)
        }
        rownames(out) <- NULL
        out
    }
    sc9 <- vapply(srcs, function(s) compare(build_9(s), old9,
                  c("gene", "dataset"), "", TRUE), numeric(1))
    for (s in srcs) cat(sprintf("  source %-20s differences %d\n", s, sc9[[s]]))
    best9 <- srcs[which.min(sc9)]
    cat("\n  Using source:", best9, "\n")
    n9 <- build_9(best9)
    compare(n9, old9, c("gene", "dataset"), "Step 9")
    write.csv(n9, file.path(OUT_DIR, basename(F9)), row.names = FALSE)

    hdr("DONE")
    cat("Rebuilt files and this log are in:\n  ", OUT_DIR, "\n", sep = "")
    invisible(NULL)
}

tryCatch(rebuild_7to9(),
         error = function(e) cat("\nERROR:", conditionMessage(e), "\n"),
         finally = {
             sink(); close(.log_con); options(.old_opt)
             cat("Log saved:", .log_txt, "\n")
         })
