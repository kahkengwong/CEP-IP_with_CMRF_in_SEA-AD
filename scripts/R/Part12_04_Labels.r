# =============================================================================
# Part12_04_Labels.r
# SEA-AD CEP-IP/CMRF, Part 12, step 04: cell labels
#   Checks the step 02 tables against the step 02 log, rebuilds the Module
#   labels per arm and dataset (Module = Original + Refined or Recovered,
#   Post-IP = composite at or above the Part 10b IP_y), assigns each cell to one
#   of four groups, joins donor, region, library and DecontX contamination, and
#   recomputes each arm's composite from the step 02 counts as a check.
#
# Usage (from the SEA-AD folder, or set SEAAD_ROOT), after step 02
#   source("Part12_rerun/Part12_04_Labels.r")
#
# Outputs (Part12_rerun/04_Labels/): 04_cell_labels.csv.gz, 04_module_cuts.csv, log.
# =============================================================================

source(file.path(Sys.getenv("SEAAD_ROOT", getwd()), "Part12_rerun", "Part12_03_Config.r"))

STEP <- "04"
D <- p12_step_dirs(STEP, "Labels")
p12_clean(c(D$tables, D$figures))
p12_start(STEP, "Part12_04_Labels.r", D)

F_LABELS <- file.path(D$tables, "04_cell_labels.csv.gz")
F_CUTS <- file.path(D$tables, "04_module_cuts.csv")

main <- function() {
    TC <- P12_CFG$table_columns
    LB <- P12_CFG$labels
    POST <- TC$ip_status_post
    PRE <- TC$ip_status_pre

    # == A. step 02 tables ====================================================
    p12_say("SECTION A -- step 02 tables copied from Fedora")
    d02 <- file.path(P12_RERUN, "02_Extract")
    log02 <- file.path(d02, "Part12_02_Extract_log.txt")
    if (!p12_chk(file.exists(log02), "step 02 log present", p12_rel(log02)))
        p12_abort("copy the whole 02_Extract folder (tables and log) from Fedora")
    want <- p12_log_outputs(log02)
    n_fail0 <- P12$n[["FAIL"]]
    p12_chk(length(want) == 6L, "step 02 log lists six outputs",
            sprintf("%d listed", length(want)))
    for (rel in names(want)) {
        f <- file.path(d02, basename(rel))
        p12_chk(file.exists(f), sprintf("%s listed in the step 02 log and present", basename(f)),
                p12_rel(f))
    }
    if (P12$n[["FAIL"]] > n_fail0) p12_abort("a step 02 table listed in the log is missing")

    C02 <- utils::read.csv(file.path(d02, "02_cells.csv.gz"),
                           check.names = FALSE,
                           colClasses = c(barcode = "character", donor_id = "character",
                                          region = "character", supertype = "character",
                                          library = "character"))
    genes02 <- readLines(file.path(d02, "02_counts_genes.tsv"))
    cells02 <- readLines(file.path(d02, "02_counts_cells.tsv"))
    X02 <- Matrix::readMM(file.path(d02, "02_counts.mtx.gz"))
    X02 <- methods::as(X02, "CsparseMatrix")
    p12_chk(nrow(X02) == length(cells02) && ncol(X02) == length(genes02) &&
                identical(cells02, as.character(C02$barcode)),
            "step 02 count matrix, cell list and cell table agree",
            sprintf("%d cells x %d genes", nrow(X02), ncol(X02)))
    p12_say("")

    # == B. inputs ============================================================
    p12_say("SECTION B -- inputs")
    f_dx <- p12_log_input(p12_input("decontx"), "DecontX (Part 7)")
    f_iprs <- f_prog <- list()
    for (a in P12_ARMS) {
        f_iprs[[a]] <- p12_log_input(p12_input("iprs_rds", arm = a), paste(a, "IPRS RDS"))
        f_prog[[a]] <- p12_log_input(p12_input("prognosis", arm = a), paste(a, "Prognosis"))
    }
    f_v1 <- f_obs <- list()
    for (a in P12_ARMS) for (ds in P12_DATASETS) {
        k <- paste(a, ds)
        f_v1[[k]] <- p12_log_input(p12_input("validation1", arm = a, dataset = ds),
                                   paste(k, "Validation 1"))
        f_obs[[k]] <- p12_log_input(p12_input("obs_dataset", arm = a, dataset = ds),
                                    paste(k, "obs_dataset"))
    }
    p12_say("")

    # DecontX: region, library and contamination per barcode.
    DX <- utils::read.csv(f_dx, check.names = FALSE, colClasses = "character")
    cx <- TC$decontx
    DX <- data.frame(barcode = as.character(DX[[p12_col(DX, cx$barcode, "DecontX")]]),
                     region = as.character(DX[[p12_col(DX, cx$region, "DecontX")]]),
                     library = as.character(DX[[p12_col(DX, cx$library, "DecontX")]]),
                     contamination = p12_num(DX[[p12_col(DX, cx$contamination, "DecontX")]]))
    p12_chk(!anyDuplicated(DX$barcode), "DecontX barcodes unique",
            sprintf("%d rows", nrow(DX)))

    # IPRS: IP_x and IP_y per arm and dataset.
    IPX <- IPY <- list()
    n_fail0 <- P12$n[["FAIL"]]
    for (a in P12_ARMS) {
        R <- readRDS(f_iprs[[a]])
        IPX[[a]] <- R[[TC$iprs$ip_x]]
        IPY[[a]] <- R[[TC$iprs$ip_y]]
        p12_chk(!is.null(IPY[[a]]) && all(P12_DATASETS %in% names(IPY[[a]])),
                sprintf("%s IPRS RDS carries IP_y for both datasets", a),
                if (is.null(IPY[[a]])) "no IP_y element" else
                    paste(names(IPY[[a]]), collapse = ", "))
    }
    if (P12$n[["FAIL"]] > n_fail0) p12_abort("an IPRS RDS lacks IP_y")

    # == C. the cut, per arm and dataset ======================================
    p12_say("SECTION C -- the Module cut (IP_y), per arm and dataset")
    cuts <- list(); labs <- list()
    for (a in P12_ARMS) for (ds in P12_DATASETS) {
        k <- paste(a, ds)
        cv <- TC$validation1
        V1 <- p12_read_sheet(f_v1[[k]], cv$sheet)
        what <- paste(k, "Validation 1")
        bc <- as.character(V1[[p12_col(V1, cv$barcode, what)]])
        cd <- p12_num(V1[[p12_col(V1, cv$composite, what)]])
        ipv <- tolower(trimws(as.character(V1[[p12_col(V1, cv$ip_status, what)]])))
        cls <- trimws(as.character(V1[[p12_col(V1, cv$cell_class, what)]]))
        is_mod <- cls %in% LB$module_classes
        p12_chk(!anyDuplicated(bc), sprintf("%s barcodes unique", k),
                sprintf("%d cells | %d Module", length(bc), sum(is_mod)))
        p12_chk(all(is.finite(cd)), sprintf("%s composite finite on every cell", k),
                sprintf("%d missing", sum(!is.finite(cd))))
        lev <- sort(unique(ipv[is_mod]))
        p12_chk(length(lev) == 2L && any(lev %in% POST) && any(lev %in% PRE),
                sprintf("%s IP_status on Module cells is Post or Pre", k),
                paste(lev, collapse = ", "))

        ipx <- as.numeric(IPX[[a]][[ds]]); ipy <- as.numeric(IPY[[a]][[ds]])
        if (!is.finite(ipy)) p12_abort(sprintf("%s has no finite IP_y", k))
        thr <- ipy

        # Prognosis Post-IP target: sum over donors of pct_post x n_total.
        cp <- TC$prognosis_donor_pct
        DP <- p12_read_sheet(f_prog[[a]], .p12_fill(cp$sheet, dataset = ds))
        wp <- paste(a, "Prognosis", ds)
        dp_don <- trimws(as.character(DP[[p12_col(DP, cp$donor, wp)]]))
        dp_exp <- round(p12_num(DP[[p12_col(DP, cp$pct_post, wp)]]) / 100 *
                        p12_num(DP[[p12_col(DP, cp$n_total, wp)]]))
        p12_chk(!anyDuplicated(dp_don), sprintf("%s Prognosis donors unique", k),
                sprintf("%d donors", length(dp_don)))
        target <- as.integer(sum(dp_exp, na.rm = TRUE))

        mod_ok <- is_mod & is.finite(cd)
        drv_post <- mod_ok & cd >= thr
        got <- sum(drv_post)
        p12_chk(abs(got - target) <= max(5, 0.01 * target),
                sprintf("%s cut reproduces the Prognosis Post-IP count", k),
                sprintf("IP_y %.6f -> %d Post-IP (target %d)", thr, got, target))

        # Top-N rule, as a check only: the N-th largest Module composite.
        cdm <- sort(cd[mod_ok], decreasing = TRUE)
        thr_topn <- if (target > 0 && target <= length(cdm)) cdm[target] else NA_real_
        n_between <- if (is.finite(thr_topn))
            sum(mod_ok & cd >= min(ipy, thr_topn) & cd < max(ipy, thr_topn)) else NA_integer_
        p12_chk(isTRUE(n_between == 0L),
                sprintf("%s IP_y and the top-N rule split the Module cells identically", k),
                sprintf("IP_x %.6f | IP_y %.6f | top-N %.6f | %s Module cell(s) between",
                        ipx, ipy, thr_topn, n_between))

        wb_post <- is_mod & ipv %in% POST
        n_mis <- sum(xor(drv_post, wb_post))
        p12_chk(n_mis == 0L, sprintf("%s labels match the workbook IP_status", k),
                sprintf("%d of %d Module cells disagree", n_mis, sum(is_mod)))
        p_lo <- suppressWarnings(min(cd[wb_post]))
        q_hi <- suppressWarnings(max(cd[is_mod & !wb_post]))
        p12_chk(is.finite(p_lo) && is.finite(q_hi) && p_lo > q_hi,
                sprintf("%s IP_status is a clean threshold on the composite", k),
                sprintf("Post min %.6f | Pre max %.6f", p_lo, q_hi))

        # Donor from obs_dataset; region, library, contamination from DecontX.
        co <- TC$obs_dataset
        OB <- utils::read.csv(f_obs[[k]], check.names = FALSE, colClasses = "character")
        wo <- paste(k, "obs_dataset")
        ob_bc <- as.character(OB[[p12_col(OB, co$barcode, wo)]])
        p12_chk(!anyDuplicated(ob_bc), sprintf("%s obs_dataset barcodes unique", k),
                sprintf("%d rows", length(ob_bc)))
        j <- match(bc, ob_bc)
        don_raw <- as.character(OB[[p12_col(OB, co$donor, wo)]])[j]
        don <- trimws(don_raw)
        ob_reg <- as.character(OB[[p12_col(OB, co$region, wo)]])[j]
        p12_chk(all(!is.na(don) & nzchar(don)), sprintf("%s donor joined for every cell", k),
                sprintf("%d of %d", sum(!is.na(don)), length(don)))
        mi <- match(bc, DX$barcode)
        p12_chk(all(!is.na(mi)), sprintf("%s every cell found in DecontX", k),
                sprintf("%d of %d", sum(!is.na(mi)), length(mi)))
        reg <- DX$region[mi]
        n_dis <- sum(reg != ob_reg, na.rm = TRUE)
        p12_chk(n_dis == 0L, sprintf("%s region agrees between DecontX and obs_dataset", k),
                sprintf("%d cell(s) disagree", n_dis))
        # Checked on the raw values, so stray spaces and non-ASCII are reported.
        dk <- unique(don_raw[!is.na(don_raw)])
        bad <- dk[grepl("[^ -~]", dk) | dk != trimws(dk)]
        p12_chk(!length(bad), sprintf("%s donor IDs are clean ASCII", k),
                sprintf("%d donors | %d bad", length(dk), length(bad)))

        # Per-donor Post-IP counts against the Prognosis percentages.
        g <- table(don[drv_post])
        kk <- union(dp_don, names(g))
        ev <- stats::setNames(dp_exp, dp_don)[kk]; ev[is.na(ev)] <- 0
        gv <- as.numeric(g[kk]); gv[is.na(gv)] <- 0
        dv <- abs(ev - gv)
        p12_note(sprintf("%s per-donor Post-IP counts against Prognosis", k),
                 sprintf("%d donors | max deviation %d | off by more than 1: %d",
                         length(kk), as.integer(max(c(0, dv))), sum(dv > 1)))

        group <- ifelse(is_mod, ifelse(drv_post, "Module post-IP", "Module pre-IP"),
                        ifelse(is.finite(cd) & cd >= thr, "non-Module post-IP",
                               "non-Module pre-IP"))
        cuts[[k]] <- data.frame(
            arm = a, dataset = ds, cut_token = LB$cut_token,
            labels_version = LB$labels_version, threshold = thr, ip_x = ipx,
            ip_y = ipy, threshold_topn = thr_topn, n_module_between_rules = n_between,
            target_post = target, reproduced_post = got, n_v1 = length(bc),
            n_module = sum(is_mod), n_nonmodule_post = sum(group == "non-Module post-IP"),
            n_nonmodule_pre = sum(group == "non-Module pre-IP"),
            mismatch_vs_workbook = n_mis, v1_file = basename(f_v1[[k]]))
        labs[[k]] <- data.frame(
            barcode = bc, arm = a, dataset = ds, donor_id = don,
            donor_key = p12_donor_key(don), region = reg, library = DX$library[mi],
            contamination = DX$contamination[mi], composite = cd,
            cell_class = cls, is_module = is_mod,
            ip_status = ifelse(is_mod, ifelse(drv_post, "Post", "Pre"), ""),
            group = group)
        p12_say("")
    }
    CUTS <- do.call(rbind, cuts); rownames(CUTS) <- NULL
    LABS <- do.call(rbind, labs); rownames(LABS) <- NULL

    # == D. join the step 02 table and recompute the composites ===============
    p12_say("SECTION D -- join to the step 02 tables")
    j <- match(LABS$barcode, C02$barcode)
    p12_chk(all(!is.na(j)), "every labelled cell is in the step 02 table",
            sprintf("%d of %d", sum(!is.na(j)), length(j)))
    if (anyNA(j)) p12_abort("labelled cells missing from step 02")
    want_st <- vapply(LABS$arm, function(a) P12_CFG$arms[[a]]$supertype, "")
    st_ok <- C02$supertype[j] == want_st
    p12_chk(all(st_ok), "every labelled cell has its arm's supertype in the h5ad",
            sprintf("%d mismatches", sum(!st_ok)))
    dup <- anyDuplicated(paste(LABS$arm, LABS$barcode))
    p12_chk(dup == 0L, "each cell appears once per arm", sprintf("%d rows", nrow(LABS)))
    same_reg <- LABS$region == C02$region[j]
    p12_chk(all(same_reg), "region agrees between DecontX and the h5ad",
            sprintf("%d cell(s) disagree", sum(!same_reg)))
    same_don <- LABS$donor_id == trimws(C02$donor_id[j])
    p12_chk(all(same_don), "donor agrees between obs_dataset and the h5ad",
            sprintf("%d cell(s) disagree", sum(!same_don)))

    # Composite recomputed from raw counts: mean over the arm's partner genes
    # of log1p(count / total UMIs x 10,000).
    tot <- C02$total_umis[j]
    LABS$composite_recomputed <- NA_real_
    for (a in P12_ARMS) {
        part <- P12_CFG$genes[[P12_CFG$arms[[a]]$partners]]
        gi <- match(part, genes02)
        if (anyNA(gi)) p12_abort(sprintf("%s partner genes missing from step 02", a))
        r <- which(LABS$arm == a)
        M <- as.matrix(X02[j[r], gi, drop = FALSE])
        LABS$composite_recomputed[r] <- rowMeans(log1p(M / tot[r] * 1e4))
        for (ds in P12_DATASETS) {
            rr <- which(LABS$arm == a & LABS$dataset == ds)
            rs <- suppressWarnings(stats::cor(LABS$composite[rr],
                                              LABS$composite_recomputed[rr],
                                              method = "spearman"))
            md <- max(abs(LABS$composite[rr] - LABS$composite_recomputed[rr]))
            p12_chk(isTRUE(rs >= 0.99),
                    sprintf("%s %s recomputed composite agrees with the workbook", a, ds),
                    sprintf("rs %.4f | max |difference| %.4g", rs, md), warn_only = TRUE)
        }
    }

    keep02 <- setdiff(names(C02), c("barcode", "donor_id", "region", "library"))
    LABS <- cbind(LABS, C02[j, keep02, drop = FALSE])
    LABS$labels_version <- LB$labels_version
    LABS$cut_token <- LB$cut_token
    LABS$config_version <- P12_CFG$config_version
    LABS$run_id <- P12$run_id
    rownames(LABS) <- NULL
    p12_say("")

    # == E. totals and write ==================================================
    p12_say("SECTION E -- totals")
    p12_chk(nrow(CUTS) == length(P12_ARMS) * length(P12_DATASETS),
            "every arm and dataset produced a cut", sprintf("%d units", nrow(CUTS)))
    p12_chk(all(nzchar(LABS$ip_status) == LABS$is_module),
            "ip_status is set on Module cells only",
            sprintf("%d Module cells", sum(LABS$is_module)))
    tb <- table(paste(LABS$arm, LABS$dataset), LABS$group)
    for (i in seq_len(nrow(tb)))
        p12_say(sprintf("  %-14s %s", rownames(tb)[i],
                        paste(sprintf("%s %d", colnames(tb), tb[i, ]), collapse = " | ")))
    p12_say(sprintf("  donors per arm: %s", paste(sprintf("%s %d", P12_ARMS,
            vapply(P12_ARMS, function(a) length(unique(LABS$donor_key[LABS$arm == a])), 0L)),
            collapse = " | ")))
    p12_chk(p12_ascii_ok(LABS[, c("arm", "dataset", "donor_id", "region", "library",
                                  "cell_class", "ip_status", "group")]),
            "label text fields are ASCII", "checked")
    if (P12$n[["FAIL"]] > 0L) p12_abort("checks failed; labels not written")

    # Registered before writing, so a partial write is renamed _REJECTED.
    p12_output(F_LABELS)
    p12_output(F_CUTS)
    con <- gzfile(F_LABELS, "wb")
    utils::write.csv(LABS, con, row.names = FALSE)
    close(con)
    utils::write.csv(CUTS, F_CUTS, row.names = FALSE)
    back <- utils::read.csv(F_LABELS, check.names = FALSE)
    p12_chk(nrow(back) == nrow(LABS) && ncol(back) == ncol(LABS),
            "labels written and read back",
            sprintf("%d rows x %d columns", nrow(back), ncol(back)))
}

p12_run(main)
