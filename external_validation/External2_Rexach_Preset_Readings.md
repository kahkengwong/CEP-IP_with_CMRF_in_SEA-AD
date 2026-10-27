# External validation 2 (Rexach et al. 2024): pre-specified readings, core analyses

Status: COMMITTED, 2026-10-25, after the count checks of Section 2 and before any statistic of Sections 3 and 4. After the commit, a reading is changed only with a dated note that states the reason.

Two readings are fixed: the replication of the eight DFGs (Section 3) and the IP gate (Section 4). All other results are reported as observed. Design and settings are those of External validation 1 (Pineda), except where stated.

## 1. Design and fixed settings

- Data. Rexach et al. 2024 (Cell 187:5753), insular cortex, 10x 3' v3. Raw counts (raw/X) from CELLxGENE dataset aa87f914-07e0-48aa-b915-8de906c95baf, insula cells only. Age at death and PMI from Table S1 (sheet 5.Metadata_perRNAlibrary), which gives both for all 37 donors and no APOE genotype.
- Donors. 37 (10 AD, 9 bvFTD, 10 PSP, and 8 controls), one insula library each. All enter the CEP-IP fit.
- Authors' processing, kept as published. Cell Ranger 3.0 (Ensembl 93); top barcodes by UMI; nuclei with < 200 genes or > 8% MT removed; DoubletFinder; low-quality subclusters removed.
- Our processing. Nuclei with < 500 genes removed. Expression is ln(CP10K + 1). No further cell calling, %MT cutoff, or doublet removal.
- Cell types. MapMyCells Deep Generative Mapping (vae_cell_type_mapper_scvi 2.0.2) to the SEA-AD taxonomy (CCN20230505), no confidence cutoff. Hierarchical Mapping (cell_type_mapper 1.7.4) as a second mapping.
- Arms. PVM_3 (primary): Micro-PVM_3 nuclei with MT-CO3 > 0. HKG8 (control): PVM_3 arm nuclei with ACTB > 0, ACTB against cHKG8.
- Method. `cepip(x, y, donor, id, cmrf = TRUE, ip_k = 10, ip_axis = "y", kappa = 0.5, bw = 0.1, repeats = 40)`, cepip 0.1.0 at commit fe4b26934d2f4ff5aa1331a56f0e40ab70d8234f. Composite: mean of the eight SEA-AD DFGs (MT-ATP6, MT-CO1, MT-CO2, MT-CYB, MT-ND1, MT-ND2, MT-ND3, and MT-ND4).
- Significance. FDR < 0.05, Benjamini-Hochberg within each family.
- Known limitation. The authors' 8% MT cut and percent_mito subcluster filter may remove module nuclei and lower r_s. A failed reading is reported as failed, with these filters named as a possible reason.

## 2. Count checks (R11 to R13)

1. 3,710 Micro-PVM_3 nuclei. Four donors have < 5 (R24_P2730, R24_P2603, R24_P2423, and R24_P2737). Donor floor: 5 nuclei, leaving 33 donors (10 AD, 8 bvFTD, 8 PSP, and 7 controls).
2. PVM_3 arm 3,417 nuclei; HKG8 arm 1,844. Median 3,047 UMIs and 1.41% MT per arm nucleus. Expected r_s at this depth (SEA-AD Stats_15H, 3,000 UMIs): 0.630 (0.626 to 0.632). It sets no threshold.
3. The 500-gene floor removed 2.4% of insula nuclei. Of the Micro-PVM_3 nuclei, 40.8% are Micro-PVM_3 by Hierarchical Mapping. All are microglia in the authors' labels.
4. Maximum %MT in the PVM_3 arm is 8.03%, which confirms the authors' cut.
5. Taxonomy names match the SEA-AD Multiregion release.

## 3. Do the eight DFGs replicate?

- Estimate. `cepip_screen(goi = "MT-CO3", min_pct = 0, min_mean = 0)` on all Micro-PVM_3 nuclei of the 33 floor donors. r_s and tau_b of MT-CO3 with each DFG.
- Tiers. Tier 1: r_s >= 0.6 and tau_b >= 0.5 (SEA-AD). Tier 2: r_s >= 0.5 and tau_b >= 0.4 (CEP-IP, CMPB 2026).
- Reading.
  - Five or more pass tier 1: replicates at the SEA-AD thresholds.
  - Otherwise, five or more pass tier 2: replicates at the CEP-IP lower thresholds.
  - Otherwise: does not replicate.
- Section 4 runs on the fixed composite whatever this reading is.

## 4. IP gate

The IPRS tiers of CEP-IP (CMPB 2026, Section 2.8) apply to each arm.
- Very Strong or Strong (IPRS >= 0.60). The diagnosis contrasts of Section 6 are run.
- Moderate or Weak (0.20 to < 0.60), or any of C1 to C5 equal to 0. Review by the mosaic test of CMPB 2026 Section 3.12: the share of TREP cells above the main GAM curve in pre-IP against post-IP (split at x < IP_x), two-sided Fisher's exact test. Confirmed if the post-IP share is higher with p < 0.05; the contrasts are then run and labelled conditional on the review. Otherwise they are not run.
- Failed (< 0.20), no IP, or < 50 GOI-positive cells. The contrasts are not run.

## 5. Association with diagnosis

No test is pre-specified, because the insula has fewer than 10 control donors (7 after the floor). The contrasts are reported in Section 6.

## 6. Reported as observed

- Screen. Within-positive r_s and tau_b of MT-CO3 with cDFG; the agnostic screen at both tiers; MT-CO3's rank.
- Pipeline. SKCS-core, OriRef, Recovered, and gate counts; Validation1; the same-cell control (cDFG against cHKG9); the share of Module post-IP cells that are Micro-PVM_3 in both mappings.
- Known-filter test. Share of Module post-IP cells and other arm cells removed by 1% and 5% MT cutoffs, beside the SEA-AD share of Module post-IP cells above 8% MT.
- Diagnosis contrasts. Predictor: Module post-IP% per donor (`cepip_donor()`, `n_total` = all Micro-PVM_3 nuclei), 33 floor donors. Model: `cepip_associate()` with age at death, sex (F = 0, M = 1), and PMI (cov3: the SEA-AD cov4 covariates without APOE e4 dose, which Table S1 does not give). Contrasts: all three diagnoses against controls (`min_donors = 20`), and AD, bvFTD, and PSP each against controls (`min_donors = 15`). HKG8 is fitted beside, with the PVM_3 minus HKG8 difference.
- DEGs by `cepip_deg()` with the S8 thresholds.

## Dated notes

- 2026-10-25, before any statistic of Sections 3 and 4. Section 4 review: pre-IP and post-IP are split on the GOI axis (x < IP_x), as in CMPB 2026 Section 3.12; the cDFG-axis split is reported beside it. Two-sided; confirmed only when the post-IP share is higher.
- 2026-10-25. Section 6 diagnosis contrasts: PMI added to the covariates (cov3), because Table S1 gives PMI for all donors and the SEA-AD model adjusted for PMI.
