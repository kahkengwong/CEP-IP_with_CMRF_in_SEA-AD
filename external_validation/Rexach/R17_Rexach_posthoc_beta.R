# R17_Rexach_posthoc_beta.R
# Post hoc exploratory analyses of the Rexach diagnosis contrasts (R15, cov3 model). These are not in
# the pre-specified readings. They ask two questions of the four PVM_3 betas.
#
#   1. Is the beta stable?
#      a. Freedman-Lane permutation test of the Module post-IP% coefficient, with the cov3
#         covariates kept in the model (100,000 permutations). With covariates, the labels cannot
#         simply be relabelled, so the residuals of the covariate-only model are permuted.
#      b. Leave-one-donor-out (LODO) refits, with DFBETAS and Cook's distance from the full fit.
#      c. Donor bootstrap stratified by diagnosis (10,000 resamples), percentile and BCa 95% CI.
#      d. Firth logistic regression (bias-reduced, brglm2 AS_mean), odds ratio per SD of Module
#         post-IP% with Wald 95% CI.
#   2. Is it consistent with SEA-AD?
#      a. z test of the difference between each Rexach beta and the SEA-AD PVM_3 cognitive-status
#         beta (NeoC 0.465, NonNeoC 0.499, cov4). The SEA-AD SE is derived from its beta and p value.
#      b. Whether the SEA-AD beta lies inside the Rexach 95% CI.
#      c. Donors needed for 80% power at the SEA-AD beta, scaling the Rexach SE by the residual df
#         with the same group ratio and covariates. The SEA-AD beta is used as the target, because
#         power computed at the observed beta is a function of the observed p value.
#   The SEA-AD outcome is cognitive status and the Rexach outcome is diagnosis, so the comparison is
#   of standardized betas from models of the same form (cepip_associate), not of the same outcome.
#
# Model, as cepip_associate(): OLS on z-scored variables, diagnosis (control = 0, disease = 1) on
# Module post-IP% with age at death, sex, and PMI (cov3).
#
# Guardrails (ledger R17_posthoc/R17_Rexach_Guardrail_Log.txt, PASS / WARN / FAIL):
#   GR-01  input files present                  GR-07  LODO refits finite
#   GR-02  R15 finished without a FAIL           GR-08  bootstrap failures below 1%
#   GR-03  donor table complete                  GR-09  Firth fits converged
#   GR-04  refit reproduces the R15 contrasts    GR-10  donors needed found
#   GR-05  SEA-AD reference values read          GR-11  R warnings during the run
#   GR-06  permutation statistic reproduces lm   GR-12  output files written
# A FAIL in GR-01 to GR-05 stops the run before any analysis.
#
# Inputs, in External_Dataset_2:
#   R15_readings/R15_Rexach_contrasts.xlsx and R15_readings/R15_Rexach_Guardrail_Log.txt
#   The SEA-AD reference values shipped with cepip (reference/seaad_reference.csv)
# Outputs, in R17_posthoc/:
#   R17_Rexach_posthoc_beta.xlsx, R17_Rexach_Guardrail_Log.txt, and R17_Rexach_log.txt
#
# Packages: cepip, openxlsx, and brglm2 (install.packages("brglm2")).
# Run in R on Windows:
#   source("C:/Users/Wong/Desktop/SEA-AD/External_Dataset_2/R17_Rexach_posthoc_beta.R")

# The data folder: REXACH_ROOT when set (the container sets it to /data), otherwise the Windows folder.
ROOT <- Sys.getenv("REXACH_ROOT", "C:/Users/Wong/Desktop/SEA-AD/External_Dataset_2")

SET <- list(model = "cov3", covariates = c("age_at_death", "sex", "pmi"), seed = 20261025L,
            n_perm = 100000L, n_boot = 10000L, power = 0.80, alpha = 0.05, max_n = 5000L,
            # SEA-AD donors with a Module post-IP% in each region set (cepip reference, 07_donor).
            # The SE of the SEA-AD beta is derived with df = n - 6. Replace with the complete-case
            # n of the SEA-AD cov4 model if it differs.
            seaad_n = c(NeoC = 84L, NonNeoC = 83L))
CONTRASTS <- list(All_vs_Control = c("AD", "bvFTD", "PSP"), AD_vs_Control = "AD",
                  bvFTD_vs_Control = "bvFTD", PSP_vs_Control = "PSP")
EXPECT_DONORS <- c(AD = 10L, bvFTD = 8L, Control = 7L, PSP = 8L)

OUT <- file.path(ROOT, "R17_posthoc")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
LOG <- file.path(OUT, "R17_Rexach_log.txt")
GR_LOG <- file.path(OUT, "R17_Rexach_Guardrail_Log.txt")
INCOMPLETE <- file.path(OUT, "_RUN_INCOMPLETE.txt")
F_R15 <- file.path(ROOT, "R15_readings", "R15_Rexach_contrasts.xlsx")
F_GR15 <- file.path(ROOT, "R15_readings", "R15_Rexach_Guardrail_Log.txt")
XLSX <- file.path(OUT, "R17_Rexach_posthoc_beta.xlsx")

need <- c("cepip", "openxlsx", "brglm2")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Install these packages first: ", paste(miss, collapse = ", "))

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
    writeLines(c(strrep("=", 72), "R17 Rexach post hoc beta analyses - Guardrail Log",
                 paste("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
                 "Post hoc and exploratory; not part of the pre-specified readings",
                 strrep("=", 72), "",
                 vapply(GR$log, function(x) sprintf("  [%s]  %-16s  %s", x$status, x$id, x$detail),
                        character(1)),
                 "", strrep("=", 72),
                 sprintf("  Total : %d", n_p + n_w + n_f), sprintf("  PASS  : %d", n_p),
                 sprintf("  WARN  : %d", n_w), sprintf("  FAIL  : %d", n_f),
                 if (!finished) "\n  Result: the run did not finish; see R17_Rexach_log.txt"
                 else if (n_f == 0 && n_w == 0) "\n  Result: no failures, no warnings"
                 else if (n_f == 0) "\n  Result: no failures; review the warnings above"
                 else "\n  Result: failures present",
                 strrep("=", 72)), GR_LOG)
}
stop_if_fail <- function() if (gr_count("FAIL") > 0)
    stop("Guardrail failure before the analyses. See the guardrail log.", call. = FALSE)

f3 <- function(v) ifelse(is.finite(v), sprintf("%.3f", v), "NA")
fp <- function(p) ifelse(!is.finite(p), "NA", ifelse(p < 0.001, sprintf("%.1e", p), sprintf("%.3f", p)))

# z-score as cepip_associate does.
zs <- function(v) {
    v <- as.numeric(v); s <- stats::sd(v)
    if (!is.finite(s) || s < 1e-12) v - mean(v) else (v - mean(v)) / s
}
# Design for one contrast: z-scored outcome, predictor, and covariates.
design <- function(d) {
    list(y = zs(d$diagnosis), x = zs(d$pct_PostIP),
         C = sapply(SET$covariates, function(v) zs(d[[v]])))
}
# Fast OLS for the predictor coefficient (column 2). NA if the predictor is not estimable.
beta_x <- function(y, X) {
    q <- qr(X)
    if (q$rank < ncol(X) && !(2L %in% q$pivot[seq_len(q$rank)])) return(NA_real_)
    unname(qr.coef(q, y)[2])
}
fit_lm <- function(d) {
    z <- design(d)
    df <- data.frame(y = z$y, x = z$x, z$C)
    stats::lm(stats::as.formula(paste("y ~ x +", paste(SET$covariates, collapse = " + "))), data = df)
}

main <- function() {
t0 <- Sys.time()
si <- sessionInfo()
say("R17_Rexach_posthoc_beta.R started %s", format(t0))
say("Post hoc and exploratory; not part of the pre-specified readings.")
say("R %s | cepip %s | brglm2 %s", getRversion(), packageVersion("cepip"), packageVersion("brglm2"))

# ---- Guardrails on the inputs --------------------------------------------------------------------
say("\n== Guardrails on the inputs")
ok <- file.exists(c(F_R15, F_GR15))
gr("GR-01-inputs", if (all(ok)) "PASS" else "FAIL",
   if (all(ok)) "R15 contrasts workbook and guardrail log are present"
   else paste("missing:", paste(basename(c(F_R15, F_GR15))[!ok], collapse = ", ")))
stop_if_fail()
g15 <- readLines(F_GR15, warn = FALSE)
n_f15 <- suppressWarnings(as.integer(sub(".*FAIL\\s*:\\s*", "", grep("^\\s*FAIL\\s*:", g15, value = TRUE))))
done15 <- !any(grepl("did not finish", g15))
gr("GR-02-R15", if (done15 && length(n_f15) == 1L && n_f15 == 0L) "PASS" else "FAIL",
   sprintf("R15 guardrail log: finished %s | FAIL count %s", done15, if (length(n_f15)) n_f15 else "not found"))
stop_if_fail()

con15 <- openxlsx::read.xlsx(F_R15, "01_Contrasts")
con15 <- con15[con15$arm == "PVM_3", ]
don <- openxlsx::read.xlsx(F_R15, "04_Donors_PVM_3")
vars <- c("donor", "dx", "pct_PostIP", SET$covariates)
have <- all(vars %in% names(don))
cnt <- if (have) table(factor(don$dx, names(EXPECT_DONORS))) else integer(0)
don_ok <- have && all(stats::complete.cases(don[, vars])) && all(cnt == EXPECT_DONORS) && nrow(con15) == 4L
gr("GR-03-donors", if (don_ok) "PASS" else "FAIL",
   sprintf("%d donors | complete cov3 data %s | per group %s | PVM_3 contrasts %d", nrow(don),
           have && all(stats::complete.cases(don[, vars])),
           paste(sprintf("%s %d", names(cnt), cnt), collapse = ", "), nrow(con15)))
stop_if_fail()

dat <- lapply(names(CONTRASTS), function(ct) {
    d <- don[don$dx %in% c("Control", CONTRASTS[[ct]]), ]
    d$diagnosis <- as.numeric(d$dx != "Control")
    d[order(d$dx, d$donor), ]
})
names(dat) <- names(CONTRASTS)

fits <- lapply(dat, fit_lm)
obs <- do.call(rbind, lapply(names(CONTRASTS), function(ct) {
    f <- fits[[ct]]; co <- summary(f)$coefficients; ci <- stats::confint(f)["x", ]
    data.frame(contrast = ct, n = nrow(dat[[ct]]), n_disease = sum(dat[[ct]]$diagnosis == 1),
               n_control = sum(dat[[ct]]$diagnosis == 0), beta_std = co["x", 1], se = co["x", 2],
               t = co["x", 3], ci_low = ci[1], ci_high = ci[2], p_value = co["x", 4],
               df_resid = f$df.residual, stringsAsFactors = FALSE)
}))
rownames(obs) <- NULL
m15 <- match(obs$contrast, con15$contrast)
dev <- max(abs(c(obs$beta_std - con15$beta_std[m15], obs$ci_low - con15$ci_low[m15],
                 obs$p_value - con15$p_value[m15])))
gr("GR-04-reproduce", if (is.finite(dev) && dev < 1e-8) "PASS" else "FAIL",
   sprintf("4 cov3 refits reproduce the R15 beta, CI, and p | largest difference %.1e", dev))

ref_f <- system.file("reference", "seaad_reference.csv", package = "cepip")
ref <- if (nzchar(ref_f)) utils::read.csv(ref_f, stringsAsFactors = FALSE) else data.frame()
get_ref <- function(set, item) {
    v <- ref$value[ref$arm == "PVM_3" & ref$region_set == set & ref$step == "09_associate" & ref$item == item]
    if (length(v) == 1L) v else NA_real_
}
sea <- data.frame(set = c("NeoC", "NonNeoC"), stringsAsFactors = FALSE)
sea$beta <- vapply(sea$set, get_ref, numeric(1), item = "beta_std_Cognitive status")
sea$p <- vapply(sea$set, get_ref, numeric(1), item = "p_value_Cognitive status")
sea$n <- unname(SET$seaad_n[sea$set])
sea$df <- sea$n - 6L
sea$se <- abs(sea$beta) / stats::qt(1 - sea$p / 2, sea$df)
sea$ci_low <- sea$beta - stats::qt(0.975, sea$df) * sea$se
sea$ci_high <- sea$beta + stats::qt(0.975, sea$df) * sea$se
sea_ok <- all(is.finite(sea$se)) && all(abs(sea$beta - c(0.4645, 0.4986)) < 5e-4)
gr("GR-05-SEAAD", if (sea_ok) "PASS" else "FAIL",
   sprintf("SEA-AD cognitive status (cov4): NeoC beta %s (SE %s, n %d) | NonNeoC beta %s (SE %s, n %d)",
           f3(sea$beta[1]), f3(sea$se[1]), sea$n[1], f3(sea$beta[2]), f3(sea$se[2]), sea$n[2]))
stop_if_fail()

say("\n== Observed (cov3, as R15)")
print(data.frame(contrast = obs$contrast, n = obs$n, beta = f3(obs$beta_std), se = f3(obs$se),
                 ci = paste(f3(obs$ci_low), "to", f3(obs$ci_high)), p = fp(obs$p_value)), row.names = FALSE)

# ---- 1a. Freedman-Lane permutation ---------------------------------------------------------------
say("\n== 1a. Freedman-Lane permutation test (%s permutations)", format(SET$n_perm, big.mark = ","))
perm <- list(); t_dev <- 0
for (k in seq_along(dat)) {
    ct <- names(dat)[k]; z <- design(dat[[ct]]); n <- length(z$y)
    Xr <- cbind(1, z$C); Xf <- cbind(1, z$x, z$C)
    fit_r <- drop(Xr %*% solve(crossprod(Xr), crossprod(Xr, z$y))); e_r <- z$y - fit_r
    XtXi <- solve(crossprod(Xf)); H <- XtXi %*% t(Xf); pdf <- n - ncol(Xf)
    tstat <- function(Y) {
        B <- H %*% Y; R <- Y - Xf %*% B
        B[2, ] / sqrt(colSums(R^2) / pdf * XtXi[2, 2])
    }
    t_obs <- tstat(matrix(z$y))
    t_dev <- max(t_dev, abs(t_obs - obs$t[k]))
    set.seed(SET$seed + k)
    tp <- numeric(0); left <- SET$n_perm
    while (left > 0) {
        b <- min(left, 20000L)
        P <- replicate(b, sample.int(n))
        tp <- c(tp, tstat(matrix(fit_r, n, b) + matrix(e_r[P], n, b)))
        left <- left - b
    }
    p_perm <- (sum(abs(tp) >= abs(t_obs) - 1e-12) + 1) / (length(tp) + 1)
    perm[[ct]] <- data.frame(contrast = ct, n = n, t_obs = t_obs, n_perm = length(tp),
                             n_as_extreme = sum(abs(tp) >= abs(t_obs) - 1e-12), p_perm = p_perm,
                             p_ols = obs$p_value[k], stringsAsFactors = FALSE)
}
perm <- do.call(rbind, perm); rownames(perm) <- NULL
gr("GR-06-perm_stat", if (t_dev < 1e-8) "PASS" else "FAIL",
   sprintf("permutation t statistic equals the lm t for the observed data | largest difference %.1e", t_dev))
print(data.frame(contrast = perm$contrast, t = f3(perm$t_obs), p_perm = fp(perm$p_perm), p_ols = fp(perm$p_ols)),
      row.names = FALSE)

# ---- 1b. Leave-one-donor-out -----------------------------------------------------------------------
say("\n== 1b. Leave-one-donor-out (LODO)")
lodo <- list(); lodo_sum <- list(); jack <- list()
for (ct in names(dat)) {
    d <- dat[[ct]]; n <- nrow(d); f <- fits[[ct]]
    b_i <- vapply(seq_len(n), function(i) {
        z <- design(d[-i, ]); beta_x(z$y, cbind(1, z$x, z$C))
    }, numeric(1))
    jack[[ct]] <- b_i
    dfb <- stats::dfbetas(f)[, "x"]; ck <- stats::cooks.distance(f)
    b0 <- obs$beta_std[obs$contrast == ct]
    lodo[[ct]] <- data.frame(contrast = ct, donor = d$donor, dx = d$dx, pct_PostIP = d$pct_PostIP,
                            beta_without = b_i, change = b_i - b0, dfbetas = unname(dfb),
                            cooks_d = unname(ck), flag_dfbetas = abs(unname(dfb)) > 2 / sqrt(n),
                            flag_cooks = unname(ck) > 4 / n, stringsAsFactors = FALSE)
    w <- which.max(abs(b_i - b0))
    lodo_sum[[ct]] <- data.frame(contrast = ct, n = n, beta_std = b0, lodo_min = min(b_i), lodo_max = max(b_i),
                                lodo_positive = sum(b_i > 0), most_influential = d$donor[w],
                                most_influential_dx = d$dx[w], beta_without_it = b_i[w],
                                n_flag_dfbetas = sum(abs(dfb) > 2 / sqrt(n)), n_flag_cooks = sum(ck > 4 / n),
                                stringsAsFactors = FALSE)
}
lodo <- do.call(rbind, lodo); lodo_sum <- do.call(rbind, lodo_sum); rownames(lodo) <- rownames(lodo_sum) <- NULL
gr("GR-07-LODO", if (all(is.finite(lodo$beta_without))) "PASS" else "WARN",
   sprintf("%d of %d leave-one-donor-out (LODO) refits finite", sum(is.finite(lodo$beta_without)), nrow(lodo)))
print(data.frame(contrast = lodo_sum$contrast, beta = f3(lodo_sum$beta_std),
                 lodo_range = paste(f3(lodo_sum$lodo_min), "to", f3(lodo_sum$lodo_max)),
                 positive = sprintf("%d of %d", lodo_sum$lodo_positive, lodo_sum$n),
                 most_influential = sprintf("%s (%s): %s", lodo_sum$most_influential, lodo_sum$most_influential_dx,
                                            f3(lodo_sum$beta_without_it))), row.names = FALSE)

# ---- 1c. Stratified donor bootstrap -----------------------------------------------------------------
say("\n== 1c. Donor bootstrap stratified by diagnosis (%s resamples)", format(SET$n_boot, big.mark = ","))
boot <- list(); n_fail <- 0L; n_all <- 0L
for (k in seq_along(dat)) {
    ct <- names(dat)[k]; d <- dat[[ct]]; strata <- split(seq_len(nrow(d)), d$dx)
    set.seed(SET$seed + 100L + k)
    tb <- vapply(seq_len(SET$n_boot), function(b) {
        idx <- unlist(lapply(strata, function(ix) ix[sample.int(length(ix), length(ix), replace = TRUE)]),
                      use.names = FALSE)
        z <- design(d[idx, ]); beta_x(z$y, cbind(1, z$x, z$C))
    }, numeric(1))
    ok_b <- is.finite(tb); n_fail <- n_fail + sum(!ok_b); n_all <- n_all + length(tb); tb <- tb[ok_b]
    b0 <- obs$beta_std[k]; J <- jack[[ct]]; L <- mean(J) - J
    a <- sum(L^3) / (6 * sum(L^2)^1.5)
    z0 <- stats::qnorm(mean(tb < b0))
    zq <- stats::qnorm(c(0.025, 0.975))
    adj <- stats::pnorm(z0 + (z0 + zq) / (1 - a * (z0 + zq)))
    bca <- if (is.finite(z0)) unname(stats::quantile(tb, adj, type = 6)) else c(NA_real_, NA_real_)
    pct <- unname(stats::quantile(tb, c(0.025, 0.975), type = 6))
    boot[[ct]] <- data.frame(contrast = ct, n_boot = length(tb), n_failed = sum(!ok_b), beta_std = b0,
                             boot_mean = mean(tb), boot_sd = sd(tb), pct_low = pct[1], pct_high = pct[2],
                             bca_low = bca[1], bca_high = bca[2], bias_z0 = z0, acceleration = a,
                             share_above_0 = mean(tb > 0), ols_low = obs$ci_low[k], ols_high = obs$ci_high[k],
                             stringsAsFactors = FALSE)
}
boot <- do.call(rbind, boot); rownames(boot) <- NULL
gr("GR-08-bootstrap", if (n_fail / n_all < 0.01 && all(is.finite(boot$bca_low))) "PASS" else "WARN",
   sprintf("%d of %d resamples failed (%.2f%%) | BCa CI finite for %d of 4", n_fail, n_all,
           100 * n_fail / n_all, sum(is.finite(boot$bca_low))))
print(data.frame(contrast = boot$contrast, OLS = paste(f3(boot$ols_low), "to", f3(boot$ols_high)),
                 percentile = paste(f3(boot$pct_low), "to", f3(boot$pct_high)),
                 BCa = paste(f3(boot$bca_low), "to", f3(boot$bca_high)), above_0 = f3(boot$share_above_0)),
      row.names = FALSE)

# ---- 1d. Firth logistic regression --------------------------------------------------------------------
say("\n== 1d. Firth logistic regression (brglm2, AS_mean)")
# brglmFit works on adjusted responses, and glm's binomial step then warns "non-integer #successes".
# The outcome itself is 0/1, so only that warning is muffled here; it is counted in GR-09.
n_nonint <- 0L
firth <- do.call(rbind, lapply(names(dat), function(ct) {
    d <- dat[[ct]]
    df <- data.frame(y = d$diagnosis, x = zs(d$pct_PostIP), sapply(SET$covariates, function(v) zs(d[[v]])))
    f <- withCallingHandlers(
        stats::glm(stats::as.formula(paste("y ~ x +", paste(SET$covariates, collapse = " + "))),
                   family = stats::binomial(), data = df, method = brglm2::brglmFit, type = "AS_mean"),
        warning = function(w) if (grepl("non-integer #successes", conditionMessage(w), fixed = TRUE)) {
            n_nonint <<- n_nonint + 1L; invokeRestart("muffleWarning") })
    co <- summary(f)$coefficients
    est <- co["x", 1]; se <- co["x", 2]
    data.frame(contrast = ct, n = nrow(d), log_or = est, se = se, or_per_sd = exp(est),
               or_low = exp(est - 1.959964 * se), or_high = exp(est + 1.959964 * se),
               p_wald = co["x", 4], converged = isTRUE(f$converged), stringsAsFactors = FALSE)
}))
y01 <- all(vapply(dat, function(d) all(d$diagnosis %in% c(0, 1)), logical(1)))
gr("GR-09-Firth", if (all(firth$converged) && all(is.finite(firth$log_or)) && y01) "PASS" else "WARN",
   sprintf("%d of 4 Firth fits converged with finite estimates | outcome 0/1: %s | brglm2 non-integer warnings muffled: %d",
           sum(firth$converged & is.finite(firth$log_or)), y01, n_nonint))
print(data.frame(contrast = firth$contrast, OR_per_SD = f3(firth$or_per_sd),
                 CI = paste(f3(firth$or_low), "to", f3(firth$or_high)), p = fp(firth$p_wald)), row.names = FALSE)

# ---- 2a and 2b. Consistency with SEA-AD ----------------------------------------------------------------
say("\n== 2a, 2b. Difference from the SEA-AD beta (cognitive status, cov4)")
cmp <- do.call(rbind, lapply(seq_len(nrow(obs)), function(k) do.call(rbind, lapply(seq_len(nrow(sea)), function(j) {
    dz <- (obs$beta_std[k] - sea$beta[j]) / sqrt(obs$se[k]^2 + sea$se[j]^2)
    data.frame(contrast = obs$contrast[k], seaad_set = sea$set[j], beta_rexach = obs$beta_std[k],
               se_rexach = obs$se[k], beta_seaad = sea$beta[j], se_seaad = sea$se[j],
               difference = obs$beta_std[k] - sea$beta[j], z = dz, p_difference = 2 * stats::pnorm(-abs(dz)),
               rexach_ci_low = obs$ci_low[k], rexach_ci_high = obs$ci_high[k],
               seaad_inside_rexach_ci = sea$beta[j] >= obs$ci_low[k] & sea$beta[j] <= obs$ci_high[k],
               stringsAsFactors = FALSE)
}))))
print(data.frame(contrast = cmp$contrast, SEAAD = cmp$seaad_set, difference = f3(cmp$difference),
                 p = fp(cmp$p_difference), inside_CI = cmp$seaad_inside_rexach_ci), row.names = FALSE)

# ---- 2c. Donors needed at the SEA-AD beta ------------------------------------------------------------
say("\n== 2c. Donors needed for %.0f%% power at the SEA-AD beta", 100 * SET$power)
k_par <- 2L + length(SET$covariates)   # intercept, predictor, and covariates
pow_at <- function(n, se0, df0, beta) {
    dfn <- n - k_par; if (dfn < 1) return(0)
    se <- se0 * sqrt(df0 / dfn); tc <- stats::qt(1 - SET$alpha / 2, dfn); ncp <- beta / se
    1 - stats::pt(tc, dfn, ncp) + stats::pt(-tc, dfn, ncp)
}
need_n <- do.call(rbind, lapply(seq_len(nrow(obs)), function(k) do.call(rbind, lapply(seq_len(nrow(sea)), function(j) {
    ns <- (k_par + 2L):SET$max_n
    pw <- vapply(ns, pow_at, numeric(1), se0 = obs$se[k], df0 = obs$df_resid[k], beta = sea$beta[j])
    hit <- which(pw >= SET$power)[1]
    nn <- if (is.na(hit)) NA_integer_ else ns[hit]
    share_c <- obs$n_control[k] / obs$n[k]
    data.frame(contrast = obs$contrast[k], seaad_set = sea$set[j], target_beta = sea$beta[j],
               n_now = obs$n[k], power_now = pow_at(obs$n[k], obs$se[k], obs$df_resid[k], sea$beta[j]),
               donors_needed = nn, controls_needed = if (is.na(nn)) NA else ceiling(nn * share_c),
               disease_needed = if (is.na(nn)) NA else nn - ceiling(nn * share_c), stringsAsFactors = FALSE)
}))))
gr("GR-10-donors_needed", if (all(is.finite(need_n$donors_needed))) "PASS" else "WARN",
   sprintf("found for %d of %d contrast and target pairs (search to %d donors)",
           sum(is.finite(need_n$donors_needed)), nrow(need_n), SET$max_n))
print(data.frame(contrast = need_n$contrast, target = sprintf("%s %s", need_n$seaad_set, f3(need_n$target_beta)),
                 power_now = f3(need_n$power_now), donors_needed = need_n$donors_needed,
                 split = sprintf("%s control, %s disease", need_n$controls_needed, need_n$disease_needed)),
      row.names = FALSE)

# ---- Summary and save -----------------------------------------------------------------------------------
summ <- data.frame(contrast = obs$contrast, n = obs$n, beta_std = obs$beta_std, ols_ci_low = obs$ci_low,
                   ols_ci_high = obs$ci_high, p_ols = obs$p_value, p_freedman_lane = perm$p_perm,
                   bca_ci_low = boot$bca_low, bca_ci_high = boot$bca_high, lodo_min = lodo_sum$lodo_min,
                   lodo_max = lodo_sum$lodo_max, lodo_positive = lodo_sum$lodo_positive,
                   firth_or_per_sd = firth$or_per_sd, firth_or_low = firth$or_low, firth_or_high = firth$or_high,
                   firth_p = firth$p_wald,
                   p_diff_seaad_NeoC = cmp$p_difference[cmp$seaad_set == "NeoC"],
                   p_diff_seaad_NonNeoC = cmp$p_difference[cmp$seaad_set == "NonNeoC"],
                   donors_needed_NeoC_beta = need_n$donors_needed[need_n$seaad_set == "NeoC"],
                   donors_needed_NonNeoC_beta = need_n$donors_needed[need_n$seaad_set == "NonNeoC"],
                   stringsAsFactors = FALSE)
say("\n== Saving")
settings <- data.frame(setting = c(names(SET), "SEA-AD SE", "Permutation", "Bootstrap", "Firth", "Power"),
                       value = c(vapply(SET, function(v) paste(names(v), v, sep = if (is.null(names(v))) "" else "=",
                                                               collapse = ", "), character(1)),
                                 "|beta| / qt(1 - p/2, n - 6), from the cepip SEA-AD reference beta and p",
                                 "Freedman-Lane: residuals of the covariate-only model permuted; statistic t of Module post-IP%; p = (extreme + 1) / (permutations + 1)",
                                 "donors resampled within diagnosis; BCa with jackknife acceleration from the LODO refits",
                                 "glm with brglm2::brglmFit, type AS_mean; Wald 95% CI",
                                 "SE scaled by sqrt(df_now / df_new); noncentral t; same group ratio and covariates"),
                       stringsAsFactors = FALSE)
session <- data.frame(item = c("started", "finished", "R", "cepip", "cepip RemoteSha", "brglm2", "openxlsx"),
                      value = c(format(t0), format(Sys.time()), as.character(getRversion()),
                                as.character(packageVersion("cepip")),
                                if (is.null(packageDescription("cepip")$RemoteSha)) NA else packageDescription("cepip")$RemoteSha,
                                as.character(packageVersion("brglm2")), as.character(packageVersion("openxlsx"))),
                      stringsAsFactors = FALSE)
sheets <- list("01_Summary" = summ, "02_Observed_cov3" = obs, "03_Permutation" = perm, "04_LODO_summary" = lodo_sum,
               "05_LODO_donors" = lodo, "06_Bootstrap" = boot, "07_Firth" = firth, "08_SEAAD_difference" = cmp,
               "09_Donors_needed" = need_n, "10_SEAAD_reference" = sea, "11_Settings" = settings,
               "12_Session" = session)
wb <- openxlsx::createWorkbook()
bold <- openxlsx::createStyle(textDecoration = "bold")
for (nm in names(sheets)) {
    openxlsx::addWorksheet(wb, nm)
    openxlsx::writeData(wb, nm, sheets[[nm]], headerStyle = bold)
    openxlsx::freezePane(wb, nm, firstRow = TRUE)
}
openxlsx::saveWorkbook(wb, XLSX, overwrite = TRUE)

gr("GR-11-R_warnings", if (length(GR$warnings)) "WARN" else "PASS",
   if (length(GR$warnings)) sprintf("%d R warning(s): %s", length(GR$warnings),
                                    paste(unique(GR$warnings), collapse = " || "))
   else "no R warnings during the run")
gr("GR-12-outputs", if (file.exists(XLSX)) "PASS" else "FAIL", basename(XLSX))
say("\nGuardrail tally: %d PASS / %d WARN / %d FAIL", gr_count("PASS"), gr_count("WARN"), gr_count("FAIL"))
say("Done in %.1f s", as.numeric(difftime(Sys.time(), t0, units = "secs")))
invisible(TRUE)
}

writeLines(paste("R17 started", format(Sys.time()), "- this file is removed when the run finishes."), INCOMPLETE)
invisible(file.remove(XLSX[file.exists(XLSX)]))
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
