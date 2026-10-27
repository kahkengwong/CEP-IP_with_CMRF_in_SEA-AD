# External validation 1 (Pineda et al. 2024): pre-specified readings, core analyses

Status: committed 2026-10-27, after the count checks of Section 2 and before any statistic of Sections 3 to 5. After the commit, a reading is changed only with a dated note that states the reason.

Three readings are fixed: the replication of the eight DFGs (Section 3), the IP gate (Section 4), and the association with diagnosis (Section 5). All other results are reported as observed, as in External validation 2. Post hoc analyses, if any, are chosen after the core results and labelled post hoc.

## 1. Design and fixed settings

- Data. Pineda et al. 2024 (Cell), prefrontal cortex (BA9), batch 210526. Reads from SRA PRJNA1073234. Donor information from Synapse syn51105515 (S01).
- Donors. There are 63: 23 ALS, 24 FTLD, and 16 PN (pathologically normal). All donors enter the CEP-IP fit. The association with diagnosis compares FTLD with PN.
- Processing. Cell Ranger 10.1.0 with GRCh38-2020-A, counting intronic reads. Cells are called with emptyDrops (FDR <= 0.01). Nuclei with fewer than 500 genes and doublets (scDblFinder) are removed. No %MT cutoff is applied. Expression is ln(CP10K + 1).
- Cell types. MapMyCells Deep Generative Mapping to the SEA-AD taxonomy (CCN20230505). Hierarchical Mapping is a second mapping.
- Arms. The PVM_3 arm (primary) is the Micro-PVM_3 nuclei with MT-CO3 > 0. The HKG8 arm (control) is the PVM_3 arm nuclei with ACTB > 0, with ACTB plotted against the mean of eight housekeeping genes (EEF1A1, RPL13A, B2M, YWHAZ, FTL, PSAP, RAB7A, and HNRNPA2B1), as in SEA-AD.
- Method. `cepip(x, y, donor, id, cmrf = TRUE, ip_k = 10, ip_axis = "y", kappa = 0.5, bw = 0.1, repeats = 40)`, cepip 0.1.0 at commit fe4b26934d2f4ff5aa1331a56f0e40ab70d8234f. The y axis is the mean of the eight SEA-AD DFGs (MT-ATP6, MT-CO1, MT-CO2, MT-CYB, MT-ND1, MT-ND2, MT-ND3, and MT-ND4).
- Significance. FDR < 0.05 (Benjamini-Hochberg).
- Order in the paper. Pineda is external validation 1 and Rexach is external validation 2. The validations then end with the dataset that has AD donors, as SEA-AD does.

## 2. Count checks (no statistic of Sections 3 to 5)

Checks 1 to 4 were done on 2026-10-27 and check 5 on 2026-10-26.

1. Micro-PVM_3 nuclei. There are 3,442: 1,097 ALS, 1,699 FTLD, and 646 PN. Four donors have fewer than 5 (ALS_118, FTLD_209, PN_302, and PN_303), so 59 donors are analysed: 22 ALS, 23 FTLD, and 14 PN. Section 5 has 22 FTLD and 14 PN.
2. Depth. The PVM_3 arm has 2,795 nuclei and the HKG8 arm has 1,632. The median arm nucleus has 1,918 UMIs and 26 MT UMIs. At this depth, the expected r_s from SEA-AD (Stats_15H) is 0.534. It sets no threshold.
3. Doublets and cell labels. The median doublet rate per library is 8.7%. All but one of the authors' microglia are Microglia-PVM in our mapping. Hierarchical Mapping agrees with Deep Generative Mapping for 36% of the Micro-PVM_3 nuclei.
4. Taxonomy names. Micro-PVM_3-SEAAD and Micro-PVM_1 are present in both mappings.
5. Covariates. S01 has age at death and sex for every donor except FTLD_226, who has no age at death. PMI is missing for 4 FTLD and 7 PN donors. S01 has no APOE, RIN, brain pH, or brain weight.

## 3. Replication of the eight DFGs

- Estimate. `cepip_screen(goi = "MT-CO3", min_pct = 0, min_mean = 0)` on all Micro-PVM_3 nuclei of donors with >= 5 nuclei (zeros included, no subsample). The GOI gate is 0 because the GOI is fixed. r_s and tau_b of MT-CO3 with each DFG, and MT-CO3's gate values.
- Tiers. Tier 1: r_s >= 0.6 and tau_b >= 0.5 (SEA-AD). Tier 2: r_s >= 0.5 and tau_b >= 0.4 (CEP-IP, CMPB 2026). The screen is run once at each tier and not lowered below tier 2.
- Reading.
  - Five or more pass tier 1: replicates at the SEA-AD thresholds.
  - Otherwise, five or more pass tier 2: replicates at the CEP-IP lower thresholds.
  - Otherwise: does not replicate.
- Five follows SEA-AD Supplementary Table S3, where MT-CO3 had 7 of the eight as tier 1 partners in NeoC and 5 in NonNeoC (2,500-nucleus subsamples).
- Sections 4 and 5 run on the fixed composite whatever this reading is.

## 4. IP gate

The IPRS tiers of CEP-IP (CMPB 2026, Section 2.8) apply to each arm.
- Very Strong or Strong (IPRS >= 0.60). Section 5 is run.
- Moderate or Weak (0.20 to < 0.60), or any of C1 to C5 equal to 0. Review by the mosaic test of CMPB 2026 Section 3.12: the share of TREP cells above the main GAM curve in pre-IP against post-IP (split at x < IP_x; the cDFG-axis split reported beside), two-sided Fisher's exact test. Confirmed if the post-IP share is higher with p < 0.05; Section 5 is then run and labelled conditional on the review. Otherwise it is not run.
- Failed (< 0.20), no IP, or < 50 GOI-positive cells. Section 5 is not run.

## 5. Association with diagnosis (primary test)

- Donors. FTLD and PN donors with >= 5 Micro-PVM_3 nuclei, as in External validation 2. At least 10 donors per group are needed. If either group has fewer after the floor, Section 5 is reported as observed, as in External validation 2.
- Predictor. Module post-IP% = PVM_3 Module post-IP cells / all Micro-PVM_3 nuclei of the donor x 100 (`cepip_donor()`, `n_total` = all Micro-PVM_3 nuclei of the 63 donors), so donors with no Module post-IP cells enter with 0.
- Model. `cepip_associate(data, "pct_PostIP", outcomes = "diagnosis", covariates = c("age_at_death", "sex"), min_donors = 20)`, with PN = 0 and FTLD = 1, and sex F = 0 and M = 1. OLS with standardized variables, as for cognitive status in SEA-AD. Complete cases: FTLD donor 226 is left out, giving 23 FTLD and 16 PN before the nuclei floor.
- Covariates (cov2). Age at death and sex. S01 has no APOE, RIN, brain pH, or brain weight. PMI is left out because it is missing for 7 of 16 PN donors: only 9 PN would remain, below 10 per group, and controls would be removed more often than FTLD donors (7 of 16 against 4 of 24).
- Family. One test.
- Reading.
  - Beta > 0 with FDR < 0.05: a higher Module post-IP% is associated with FTLD, the SEA-AD direction.
  - Beta < 0 with FDR < 0.05: reported as its own outcome.
  - FDR >= 0.05: not detected in this sample. Stated with the SEA-AD donor requirement for PVM_3 MMSE (30 donors, Figure 15D, estimated for a continuous score).
- HKG8 is fitted beside, with the PVM_3 minus HKG8 difference. If the HKG8 cepip run stops with an error, HKG8 is reported as not available and Section 5 is run on PVM_3.

## 6. Reported as observed

- Screen. Within-positive r_s and tau_b of MT-CO3 with cDFG; the agnostic screen at both tiers; MT-CO3's rank.
- Pipeline. SKCS-core, OriRef, Recovered, and gate counts; Validation1; the same-cell control (cDFG against cHKG9); the share of Module post-IP cells that are Micro-PVM_3 in both mappings.
- Known-filter test. Share of Module post-IP cells and other arm cells removed by 1%, 5%, and 10% MT cutoffs, and by the authors' published QC.
- Section 5 with PMI added (18 FTLD and 9 PN after the nuclei floor).
- ALS against PN, with the Section 5 model and HKG8 beside.
- DEGs by `cepip_deg()` with the S8 thresholds.
- Every dataset examined for external validation, with its result, in a supplementary table.

## Dated notes

- 2026-10-26, before any Pineda statistic. Covariates set to age at death and sex (count check 5). Section 6 aligned with External validation 2, without a pooled estimate.
- 2026-10-26, before any Pineda statistic. A cepip error in the HKG8 control arm no longer stops the PVM_3 analysis: HKG8 is reported as not available (Section 5).
