# =============================================================================
# Part12_05_Units.r
# SEA-AD CEP-IP/CMRF, Part 12 rerun, step 05, donor-region units and donors
#   Checks the step 04 labels and the library, donor and region nesting, then
#   builds the donor-region unit table and the donor table (Module post-IP%,
#   %MT, non-Module supertype contamination, cov4, outcomes and tissue-quality
#   fields) and adds the tissue-quality residuals on the unit table.
#
# Usage (from the SEA-AD folder, or set SEAAD_ROOT), after step 04
#   source("Part12_rerun/Part12_05_Units.r")
#
# Outputs (Part12_rerun/05_Units/): 05_units.csv, 05_donors.csv and log.
# =============================================================================

source(file.path(Sys.getenv("SEAAD_ROOT", getwd()), "Part12_rerun", "Part12_03_Config.r"))

STEP <- "05"
D <- p12_step_dirs(STEP, "Units")
p12_clean(c(D$tables, D$figures))
p12_start(STEP, "Part12_05_Units.r", D)

F_UNITS <- file.path(D$tables, "05_units.csv")
F_DONORS <- file.path(D$tables, "05_donors.csv")

main <- function() {
    TC <- P12_CFG$table_columns
    CL <- P12_CFG$clinical
    MIN_CELLS <- as.integer(P12_CFG$floors$cells_per_unit)

    # == A. Step 04 labels ===================================================
    p12_say("SECTION A -- step 04 labels")
    f_lab <- p12_upstream_file("04_Labels", "04_cell_labels.csv.gz")
    L <- utils::read.csv(f_lab, check.names = FALSE,
                         colClasses = c(barcode = "character", donor_id = "character",
                                        donor_key = "character", region = "character",
                                        library = "character", arm = "character",
                                        dataset = "character", group = "character"))
    p12_chk(all(L$labels_version == P12_CFG$labels$labels_version) &&
                all(L$cut_token == P12_CFG$labels$cut_token),
            "labels carry this config's labels version and cut token",
            paste(unique(L$labels_version), unique(L$cut_token)))
    # Labels depend only on the labels section of the config, checked above.
    p12_note("config version of the labels",
             sprintf("labels %s | this run %s", paste(unique(L$config_version), collapse = ","),
                     P12_CFG$config_version))
    L$post <- L$group == "Module post-IP"
    p12_say(sprintf("  %d labelled cells | %d arms | %d donors", nrow(L),
                    length(unique(L$arm)), length(unique(L$donor_key))))
    p12_say("")

    # == B. Nesting checks ===================================================
    p12_say("SECTION B -- libraries, donors, regions and datasets")
    cells <- L[!duplicated(L$barcode), ]
    dpl <- tapply(cells$donor_key, cells$library, function(x) length(unique(x)))
    p12_chk(all(dpl == 1L), "each library belongs to one donor",
            sprintf("%d libraries | %d span more than one donor", length(dpl), sum(dpl > 1L)))
    libcols <- grep("^lib_", names(L), value = TRUE)
    for (lc in libcols) {
        v <- tapply(cells[[lc]], cells$library, function(x) length(unique(x[!is.na(x)])))
        p12_chk(all(v <= 1L), sprintf("%s is constant within each library", lc),
                sprintf("%d libraries with more than one value", sum(v > 1L)))
    }
    rds <- tapply(L$dataset, L$region, function(x) length(unique(x)))
    p12_chk(all(rds == 1L), "each region belongs to one dataset",
            paste(sprintf("%s:%s", names(rds),
                          tapply(L$dataset, L$region, function(x) unique(x)[1])),
                  collapse = " "))
    p12_say("")

    # Every cell of the two supertypes, from step 02. Module post-IP% uses the
    # arm's supertype cells or its labelled cells as the denominator.
    f02 <- p12_upstream_file("02_Extract", "02_cells.csv.gz")
    C02 <- utils::read.csv(f02,
                           check.names = FALSE,
                           colClasses = c(barcode = "character", donor_id = "character",
                                          region = "character", supertype = "character"))
    reg_ds <- tapply(L$dataset, L$region, `[`, 1)
    C02$dataset <- reg_ds[C02$region]
    C02$donor_key <- p12_donor_key(C02$donor_id)
    DENOM <- P12_CFG$quantities$module_post_pct_denominator
    if (!(DENOM %in% c("supertype", "labelled")))
        p12_abort("quantities$module_post_pct_denominator must be supertype or labelled")
    p12_note("Module post-IP% denominator", DENOM)
    n_super_unit <- table(paste(C02$supertype, C02$dataset, C02$donor_key, C02$region,
                                sep = "|"))
    n_super_donor <- table(paste(C02$supertype, C02$dataset, C02$donor_key, sep = "|"))
    arm_st <- vapply(P12_ARMS, function(a) P12_CFG$arms[[a]]$supertype, "")
    p12_say("")

    # == C. Unit table =======================================================
    p12_say(sprintf("SECTION C -- donor-region units (at least %d cells of the arm)", MIN_CELLS))
    key <- paste(L$arm, L$dataset, L$donor_key, L$region, sep = "|")
    grp <- split(seq_len(nrow(L)), key)
    n_all <- lengths(grp)
    grp <- grp[n_all >= MIN_CELLS]
    med <- function(x) if (all(is.na(x))) NA_real_ else stats::median(x, na.rm = TRUE)
    mn <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
    # A library field missing for every cell is kept numeric (unit value NA).
    libnum <- libcols[vapply(libcols, function(lc)
        is.numeric(L[[lc]]) || (is.logical(L[[lc]]) && all(is.na(L[[lc]]))), TRUE)]
    libchr <- setdiff(libcols, libnum)
    U <- do.call(rbind, lapply(grp, function(i) {
        r <- L[i, ]
        base <- data.frame(
            arm = r$arm[1], dataset = r$dataset[1], donor_id = r$donor_id[1],
            donor_key = r$donor_key[1], region = r$region[1],
            n_cells = length(i), n_libs = length(unique(r$library)),
            libs = paste(sort(unique(r$library)), collapse = ";"),
            n_post = sum(r$post),
            mito_pct = mn(r$mito_pct),
            contamination = mn(r$contamination),
            cyto_share = mn(r$cyto_share),
            ctrl_pct = mn(r$ctrl_pct),
            ctrl_pct_null = mn(r$ctrl_pct_null_mean),
            median_umi = stats::median(r$total_umis))
        for (lc in libnum) base[[lc]] <- med(r[[lc]])
        for (lc in libchr) base[[lc]] <- paste(sort(unique(r[[lc]][!is.na(r[[lc]])])),
                                              collapse = ";")
        base
    }))
    rownames(U) <- NULL
    U$n_super <- as.integer(n_super_unit[paste(arm_st[U$arm], U$dataset, U$donor_key,
                                               U$region, sep = "|")])
    U$module_post_pct_labelled <- 100 * U$n_post / U$n_cells
    U$module_post_pct_supertype <- 100 * U$n_post / U$n_super
    U$module_post_pct <- if (DENOM == "supertype") U$module_post_pct_supertype else
        U$module_post_pct_labelled
    p12_chk(all(is.finite(U$n_super) & U$n_super >= U$n_cells),
            "every unit has at least as many supertype cells as labelled cells",
            sprintf("%d units", nrow(U)))
    U$mito_log1p <- log1p(U$mito_pct)
    U$log_n_cells <- log(U$n_cells)
    U$log_median_umi <- log(U$median_umi)
    U$n_units_donor <- stats::ave(U$n_cells, U$arm, U$dataset, U$donor_key, FUN = length)
    for (a in P12_ARMS) for (ds in P12_DATASETS) {
        i <- U$arm == a & U$dataset == ds
        k <- sprintf("%s|%s", a, ds)
        tot_units <- sum(startsWith(names(n_all), paste0(k, "|")))
        p12_say(sprintf("  %-6s %-8s %4d units of %4d donor-regions | %3d donors | %d cells",
                        a, ds, sum(i), tot_units, length(unique(U$donor_key[i])),
                        sum(U$n_cells[i])))
    }
    p12_chk(all(U$n_cells >= MIN_CELLS), "every unit meets the cell floor",
            sprintf("%d units", nrow(U)))
    p12_chk(all(is.finite(U$module_post_pct) & is.finite(U$mito_pct)),
            "Module post-IP% and %MT finite on every unit", "checked")
    p12_chk(mean(is.finite(U$contamination)) == 1, "contamination finite on every unit",
            sprintf("%d missing", sum(!is.finite(U$contamination))))
    p12_say("")

    # == D. Donor table ======================================================
    p12_say("SECTION D -- donor table")
    # Donor-level fields from the source arm's Prognosis Donor_Table.
    ct <- TC$donor_table
    f_dt <- p12_log_input(p12_input("prognosis", arm = ct$source_arm),
                          paste(ct$source_arm, "Prognosis"))
    DT <- p12_read_sheet(f_dt, ct$sheet)
    w <- paste(ct$source_arm, ct$sheet)
    dt_key <- p12_donor_key(DT[[p12_col(DT, ct$donor, w)]])
    p12_chk(!anyDuplicated(dt_key), "Donor_Table has one row per donor",
            sprintf("%d donors", length(dt_key)))
    need <- unique(c(unlist(CL$cov4),
                     vapply(CL$outcomes, function(o) o$column, "")[
                         vapply(CL$outcomes, function(o) is.null(o$source), TRUE)],
                     CL$tissue_quality$donor_table_fields))
    DV <- data.frame(donor_key = dt_key)
    for (v in need) DV[[v]] <- p12_num(DT[[p12_col(DT, v, w)]])

    # Fresh brain weight and CASI from obs_dataset; each constant within donor.
    of <- TC$obs_donor_fields
    OB <- do.call(rbind, lapply(P12_DATASETS, function(ds) {
        f <- p12_log_input(p12_input("obs_dataset", arm = ct$source_arm, dataset = ds),
                           paste(ct$source_arm, ds, "obs_dataset"))
        o <- p12_read_csv_cols(f, c(donor = TC$obs_dataset$donor,
                                    fbw = of$fresh_brain_weight, casi = of$casi),
                               paste(ct$source_arm, ds, "obs_dataset"))
        data.frame(donor_key = p12_donor_key(o$donor), fbw = p12_num(o$fbw),
                   casi = p12_num(o$casi))
    }))
    for (v in c("fbw", "casi")) {
        nv <- tapply(OB[[v]], OB$donor_key, function(x) length(unique(x[is.finite(x)])))
        p12_chk(all(nv <= 1L), sprintf("obs %s is constant within donor", v),
                sprintf("%d donors with more than one value", sum(nv > 1L)))
    }
    first <- function(x) { x <- x[is.finite(x)]; if (length(x)) x[1] else NA_real_ }
    obs_fbw <- tapply(OB$fbw, OB$donor_key, first)[DV$donor_key]
    # Donor_Table brain weight against the obs_dataset value, per donor.
    bw <- DV$BrainWeight_num
    both <- is.finite(bw) & is.finite(obs_fbw)
    p12_chk(all(abs(bw[both] - obs_fbw[both]) < 1e-6),
            "Donor_Table brain weight equals the obs_dataset Fresh Brain Weight",
            sprintf("%d donors compared | %d differ | Donor_Table %d, obs %d non-missing",
                    sum(both), sum(abs(bw[both] - obs_fbw[both]) >= 1e-6),
                    sum(is.finite(bw)), sum(is.finite(obs_fbw))), warn_only = TRUE)
    casi_col <- CL$outcomes[["CASI score"]]$column
    DV[[casi_col]] <- tapply(OB$casi, OB$donor_key, first)[DV$donor_key]

    # Per arm, dataset and donor: Module post-IP% and %MT.
    dk <- paste(L$arm, L$dataset, L$donor_key, sep = "|")
    DN <- data.frame(
        arm = tapply(L$arm, dk, `[`, 1), dataset = tapply(L$dataset, dk, `[`, 1),
        donor_key = tapply(L$donor_key, dk, `[`, 1),
        donor_id = tapply(L$donor_id, dk, `[`, 1),
        n_cells = as.integer(tapply(L$post, dk, length)),
        n_post = as.integer(tapply(L$post, dk, sum)),
        mito_pct = tapply(L$mito_pct, dk, mean, na.rm = TRUE))
    rownames(DN) <- NULL
    DN$n_super <- as.integer(n_super_donor[paste(arm_st[DN$arm], DN$dataset, DN$donor_key,
                                                 sep = "|")])
    DN$module_post_pct_labelled <- 100 * DN$n_post / DN$n_cells
    DN$module_post_pct_supertype <- 100 * DN$n_post / DN$n_super
    DN$module_post_pct <- if (DENOM == "supertype") DN$module_post_pct_supertype else
        DN$module_post_pct_labelled
    DN$mito_log1p <- log1p(DN$mito_pct)

    # Against Prognosis: Post-IP count (pct_post x n_total, compared as counts
    # because pct_post is rounded) and n_total against the denominator.
    cp <- TC$prognosis_donor_pct
    for (a in P12_ARMS) {
        f_pr <- p12_input("prognosis", arm = a)
        for (ds in P12_DATASETS) {
            DP <- p12_read_sheet(f_pr, .p12_fill(cp$sheet, dataset = ds))
            wp <- paste(a, "Prognosis", ds)
            pk <- p12_donor_key(DP[[p12_col(DP, cp$donor, wp)]])
            pp <- p12_num(DP[[p12_col(DP, cp$pct_post, wp)]])
            pn <- p12_num(DP[[p12_col(DP, cp$n_total, wp)]])
            i <- which(DN$arm == a & DN$dataset == ds)
            j <- match(DN$donor_key[i], pk)
            post_prog <- round(pp[j] * pn[j] / 100)
            d_post <- DN$n_post[i] != post_prog
            d_lab <- DN$n_cells[i] != pn[j]
            d_sup <- DN$n_super[i] != pn[j]
            d_use <- if (DENOM == "supertype") d_sup else d_lab
            dpp <- abs(DN$module_post_pct[i] - pp[j])
            p12_chk(all(!is.na(j)) && !any(d_post, na.rm = TRUE),
                    sprintf("%s %s donor Post-IP counts equal Prognosis", a, ds),
                    sprintf("%d donors | %d differ", length(i), sum(d_post, na.rm = TRUE)))
            p12_chk(all(!is.na(j)) && !any(d_use, na.rm = TRUE),
                    sprintf("%s %s Prognosis n_total equals the %s denominator", a, ds, DENOM),
                    sprintf(paste0("n_total differs from supertype cells for %d, from labelled ",
                                   "cells for %d | max |pct_post difference| %.2g pp"),
                            sum(d_sup, na.rm = TRUE), sum(d_lab, na.rm = TRUE),
                            max(c(0, dpp), na.rm = TRUE)))
        }
    }

    # Non-Module supertype contamination: donor mean DecontX over the arm's
    # supertype cells, excluding its Module cells; at least 5 cells.
    f_dx <- p12_input("decontx")
    DX <- p12_read_csv_cols(f_dx, c(barcode = TC$decontx$barcode,
                                    contamination = TC$decontx$contamination), "DecontX")
    DXc <- p12_num(DX$contamination)[match(C02$barcode, DX$barcode)]
    DN$contam_super_nonmodule <- NA_real_
    for (a in P12_ARMS) {
        st <- P12_CFG$arms[[a]]$supertype
        modbc <- L$barcode[L$arm == a & L$is_module]
        s <- C02$supertype == st & !(C02$barcode %in% modbc) & is.finite(DXc) &
            !is.na(C02$dataset)
        kk <- paste(a, C02$dataset[s], C02$donor_key[s], sep = "|")
        m <- tapply(DXc[s], kk, mean)
        n <- tapply(DXc[s], kk, length)
        m[n < 5] <- NA_real_
        i <- DN$arm == a
        DN$contam_super_nonmodule[i] <- m[paste(DN$arm[i], DN$dataset[i], DN$donor_key[i],
                                                sep = "|")]
    }
    p12_note("non-Module supertype contamination per donor",
             sprintf("%d of %d donor rows have 5 or more cells",
                     sum(is.finite(DN$contam_super_nonmodule)), nrow(DN)))

    jd <- match(DN$donor_key, DV$donor_key)
    p12_chk(all(!is.na(jd)), "every labelled donor is in the Donor_Table",
            sprintf("%d of %d donors", length(unique(DN$donor_key[!is.na(jd)])),
                    length(unique(DN$donor_key))))
    DN <- cbind(DN, DV[jd, setdiff(names(DV), "donor_key"), drop = FALSE])
    rownames(DN) <- NULL
    for (v in unique(c(unlist(CL$cov4), vapply(CL$outcomes, function(o) o$column, ""),
                       CL$tissue_quality$fields)))
        p12_say(sprintf("  %-22s %3d of %3d donors non-missing", v,
                        sum(is.finite(DV[[v]])), nrow(DV)))
    p12_say("")

    # == E. Tissue-quality residuals =========================================
    p12_say("SECTION E -- tissue-quality residuals (TQR), per arm over units")
    tq <- CL$tissue_quality$fields
    ju <- match(U$donor_key, DV$donor_key)
    TQ <- DV[ju, tq, drop = FALSE]
    names(TQ) <- paste0("tq_", make.names(tq))
    U <- cbind(U, TQ)
    tqc <- names(TQ)
    U$module_resid_tq4 <- NA_real_
    U$mito_resid_tq4 <- NA_real_
    for (a in P12_ARMS) {
        i <- which(U$arm == a)
        kp <- i[stats::complete.cases(U[i, c(tqc, "module_post_pct", "mito_log1p"),
                                         drop = FALSE])]
        if (length(kp) < 10L) {
            p12_note(sprintf("%s tissue-quality residuals not fitted", a),
                     sprintf("%d complete units", length(kp)))
            next
        }
        for (y in c("module_post_pct", "mito_log1p")) {
            f <- stats::as.formula(paste(y, "~", paste(tqc, collapse = " + ")))
            fit <- stats::lm(f, data = U[kp, ])
            out <- if (y == "module_post_pct") "module_resid_tq4" else "mito_resid_tq4"
            U[[out]][kp] <- as.numeric(stats::resid(fit))
            p12_say(sprintf("  %-6s %-16s %4d units | %3d donors | R2 %.3f", a, y, length(kp),
                            length(unique(U$donor_key[kp])), summary(fit)$r.squared))
        }
    }
    p12_note("units without all four tissue fields",
             sprintf("%d of %d", sum(!is.finite(U$module_resid_tq4)), nrow(U)))
    p12_say("")

    # == F. Write ============================================================
    p12_say("SECTION F -- write")
    U$labels_version <- DN$labels_version <- P12_CFG$labels$labels_version
    U$config_version <- DN$config_version <- P12_CFG$config_version
    U$run_id <- DN$run_id <- P12$run_id
    p12_chk(p12_ascii_ok(U[, c("arm", "dataset", "donor_id", "region", "libs")]) &&
                p12_ascii_ok(DN[, c("arm", "dataset", "donor_id")]),
            "text fields are ASCII", "checked")
    if (P12$n[["FAIL"]] > 0L) p12_abort("checks failed; tables not written")
    p12_output(F_UNITS)
    p12_output(F_DONORS)
    utils::write.csv(U, F_UNITS, row.names = FALSE)
    utils::write.csv(DN, F_DONORS, row.names = FALSE)
    bu <- utils::read.csv(F_UNITS, check.names = FALSE)
    bd <- utils::read.csv(F_DONORS, check.names = FALSE)
    p12_chk(nrow(bu) == nrow(U) && nrow(bd) == nrow(DN), "tables written and read back",
            sprintf("units %d rows | donors %d rows", nrow(bu), nrow(bd)))
}

p12_run(main)
