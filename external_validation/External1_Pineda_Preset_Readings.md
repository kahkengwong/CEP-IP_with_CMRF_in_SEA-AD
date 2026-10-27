# External validation 1 (Pineda et al. 2024): pre-specified readings, core analyses

Status: committed 2026-10-27, after the count checks of Section 2 and before any statistic of Sections 3 to 5. After the commit, a reading is changed only with a dated note that states the reason.

Three readings are fixed: the replication of the eight DFGs (Section 3), the IP gate (Section 4), and the association with diagnosis (Section 5). All other results are reported as observed, as in External validation 2. Post hoc analyses, if any, are chosen after the core results and labelled post hoc.

## 1. Design and fixed settings

- Data. Pineda et al. 2024 (Cell), BA9, batch 210526. Raw reads from SRA PRJNA1073234 as SRA Lite files, converted with fastq-dump (the pilot library 210526_FTLD_228_snRNA-C8 from the ENA FASTQ). SRA Lite carries simplified base qualities, so Q30 reads 100% in every library. Donor metadata from Synapse syn51105515 (S01).
- Donors. 63: 23 ALS (6 C9ALS), 24 FTLD (11 C9FTLD), and 16 PN, as grouped in S01. All enter the CEP-IP fit. The association uses FTLD and PN.
- Processing. Cell Ranger 10.1.0, GRCh38-2020-A, intronic reads counted. Cells called by emptyDrops (lower = 100, FDR <= 0.01, 10,000 iterations), as in Part12_06. Nuclei with < 500 genes removed. No %MT cutoff. Doublets removed by scDblFinder per library, default settings. Expression is ln(CP10K + 1).
- Cell types. MapMyCells Deep Generative Mapping to the SEA-AD taxonomy (CCN20230505), no confidence cutoff. Hierarchical Mapping as a second mapping.
- Arms. PVM_3 (primary): Micro-PVM_3 nuclei with MT-CO3 > 0. HKG8 (control): PVM_3 arm nuclei with ACTB > 0, ACTB against cHKG8 (the mean of EEF1A1, RPL13A, B2M, YWHAZ, FTL, PSAP, RAB7A, and HNRNPA2B1), nested in the primary arm as in SEA-AD.
- Method. `cepip(x, y, donor, id, cmrf = TRUE, ip_k = 10, ip_axis = "y", kappa = 0.5, bw = 0.1, repeats = 40)`, cepip 0.1.0 at commit fe4b26934d2f4ff5aa1331a56f0e40ab70d8234f. Composite: mean of the eight SEA-AD DFGs (MT-ATP6, MT-CO1, MT-CO2, MT-CYB, MT-ND1, MT-ND2, MT-ND3, and MT-ND4).
- Significance. FDR < 0.05, Benjamini-Hochberg within each family.
- Order in the paper. Pineda (BA9; ALS and FTLD) is external validation 1, and Rexach (insula; AD, bvFTD, and PSP) is external validation 2, so that the validations end with the dataset that includes AD donors, as SEA-AD does.

## 2. Count checks (no statistic of Sections 3 to 5)

1. Micro-PVM_3 nuclei per donor, by diagnosis (done 2026-10-27). Deep Generative Mapping gives 3,442 Micro-PVM_3 nuclei: 1,097 ALS (median 27 per donor, range 1 to 311), 1,699 FTLD (median 56.5, range 3 to 436), and 646 PN (median 22.5, range 3 to 140). Four donors have < 5: ALS_118, FTLD_209, PN_302, and PN_303. So 59 donors enter Section 3 (22 ALS, 23 FTLD, and 14 PN). Section 5 has 22 FTLD and 14 PN after FTLD_226 is left out, Section 5 with PMI added has 18 FTLD and 9 PN, and ALS against PN has 22 ALS and 14 PN.
2. Median UMIs and median MT UMIs per PVM_3 arm nucleus, and the expected r_s at this depth from SEA-AD Stats_15H. It sets no threshold (done 2026-10-27). The PVM_3 arm has 2,795 nuclei (81.2% of the Micro-PVM_3 nuclei have MT-CO3 > 0), and the HKG8 arm has 1,632. Per PVM_3 arm nucleus, the median is 1,918 UMIs, 26 MT UMIs, and 4 MT-CO3 UMIs. The expected r_s is 0.534 (95% CI 0.529 to 0.538), from the SEA-AD PVM_3 arm thinned to 1,504 UMIs, the largest SEA-AD target depth at or below 1,918 (SEA-AD median 5,650 UMIs).
3. Doublets per library, and the agreement of the two mappings and of the authors' microglia labels (done 2026-10-27). scDblFinder doublets per library: median 8.72% (range 4.20% to 25.65%); ALS 10.21%, FTLD 7.94%, and PN 7.44%. Deep Generative Mapping and Hierarchical Mapping agree for 89.3% of kept nuclei at subclass and 51.5% at supertype. Of the 3,442 Micro-PVM_3 nuclei by Deep Generative Mapping, 1,237 (35.9%) are Micro-PVM_3 by Hierarchical Mapping. Of the authors' Micro nuclei, 12,212 of 12,213 are Microglia-PVM by Deep Generative Mapping. Of the Micro-PVM_3 nuclei, 2,189 (63.6%) are in the authors' published QC set, and 2,188 of these are Micro in the authors' labels.
4. Taxonomy names Micro-PVM_3-SEAAD and Micro-PVM_1 match the SEA-AD Multiregion release (done 2026-10-27). Both names are present in both mappings. Deep Generative Mapping version 2.0.2 (vae_cell_type_mapper_scvi) and Hierarchical Mapping version 1.7.4 (cell_type_mapper).
5. Covariates in S01 (done 2026-10-26). S01 has no APOE, RIN, brain pH, or brain weight. Age at death and sex: 23 FTLD (226 has no age at death) and 16 PN. PMI: 20 of 24 FTLD and 9 of 16 PN, so 19 FTLD and 9 PN have all three.

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
