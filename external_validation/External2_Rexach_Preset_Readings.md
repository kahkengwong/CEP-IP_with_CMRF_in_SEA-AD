# External validation 2 (Rexach et al. 2024): pre-specified readings, core analyses

Status: COMMITTED, 2026-10-25. Committed after the count checks of Section 2 (R11, R12, and R13) and before any statistic of Sections 3 and 4 was computed. After this commit, a reading is changed only with a dated note in this file that states the reason. The one [CHECK] item left (Table S1 ages) is data for the Section 6 diagnosis contrasts and is filled before that step.

This file fixes two readings: the replication of the eight DFGs (Section 3) and the IP gate (Section 4). The association with diagnosis is not a reading in this dataset, because the insula has fewer than 10 control donors (Section 5). Every other result of the core pipeline is reported as observed, without a reading. The design, arms, and method settings are those of External validation 1 (Pineda), except where Section 1 states otherwise.

## 1. Design and fixed settings

- Data. Rexach et al. 2024 (Cell 187:5753). Insular cortex, 10x Chromium Single Cell 3' v3. Counts from the CELLxGENE dataset aa87f914-07e0-48aa-b915-8de906c95baf (fixed version; collection c53573b2-eff4-4c5e-9ad0-b24d422dfd9b), insula cells only. Raw counts are read from raw/X (integer counts; X holds normalized values). Sex is taken from the obs column sex. The obs column development_stage gives an exact age for 10 donors and a range (60-79, or 80 and over) for 27 (R12), so age at death is taken from Table S1 [CHECK: Table S1 ages, needed only for the diagnosis contrasts of Section 6]. The library is the obs column library_id.
- Region and donors. 37 insula donors (R11): 10 AD, 9 bvFTD (CELLxGENE label "Pick disease"), 10 PSP, and 8 controls ("normal"); 20 male and 17 female. Each donor has one insula library, and each library one donor (R12). All donors enter the CEP-IP fit. The insula holds 228,671 cells, all 10x 3' v3, of which 16,976 are labelled microglial cell by the authors.
- Processing by the authors, kept as published (STAR Methods).
  - Cell Ranger 3.0 with Ensembl release 93. The paper does not state whether intronic reads were counted.
  - Cells are the top barcodes by UMI, with the number kept set from the nuclei loaded.
  - Nuclei with fewer than 200 genes or more than 8% MT counts are removed. Genes detected in fewer than 3 cells are removed.
  - Doublets are removed with DoubletFinder, run by region and diagnosis (expected rate 7.6%, adjusted for homotypic doublets).
  - Outlier samples, and subclusters judged non-representative, low-quality (including association with percent_mito), or ambiguous, are removed.
- Our processing. Nuclei with fewer than 500 genes are removed. Expression is ln(CP10K + 1) from the raw counts. No further cell calling, %MT cutoff, or doublet removal is applied.
- Differences from External validation 1. Cell calling (top barcodes against emptyDrops), Cell Ranger version and intron counting, the authors' 8% MT cut (no cut in Pineda), and the doublet tool (DoubletFinder against scDblFinder). These steps cannot be redone from the CELLxGENE copy.
- Cell types. MapMyCells Deep Generative Mapping (vae_cell_type_mapper_scvi 2.0.2) to the 10x Human MTG SEA-AD taxonomy (CCN20230505), with no confidence cutoff. SEA-AD Hierarchical Mapping (cell_type_mapper 1.7.4) is run as a second mapping. The arm is defined by Deep Generative Mapping.
- Arms.
  - PVM_3 (primary). Micro-PVM_3 nuclei with MT-CO3 > 0.
  - HKG8 (control). PVM_3 arm nuclei with ACTB > 0, with ACTB as the GOI and cHKG8 as the composite (the mean of EEF1A1, RPL13A, B2M, YWHAZ, FTL, PSAP, RAB7A, and HNRNPA2B1).
- Method settings. All are the SEA-AD values, passed explicitly: `cepip(x, y, donor, id, cmrf = TRUE, ip_k = 10, ip_axis = "y", kappa = 0.5, bw = 0.1, repeats = 40)`, with cepip 0.1.0 at commit fe4b26934d2f4ff5aa1331a56f0e40ab70d8234f. The composite is the mean of the eight SEA-AD DFGs (MT-ATP6, MT-CO1, MT-CO2, MT-CYB, MT-ND1, MT-ND2, MT-ND3, and MT-ND4).
- Significance. FDR < 0.05, Benjamini-Hochberg within the family named in each section.
- Known limitation, stated in advance. The authors' 8% MT cut and their percent_mito subcluster filter may remove module nuclei, and the 8% cut narrows the %MT range, which lowers r_s. A failed reading in Section 3 or 4 is reported as failed. These filters are named as a possible reason, and the reading is not changed. The share of SEA-AD Module post-IP cells above 8% MT, from the deposited SEA-AD outputs, is reported beside Section 3 as the expected loss to the authors' cut.

## 2. Count checks before commit (no statistic of Sections 3 and 4)

1. Micro-PVM_3 nuclei per donor, by diagnosis (R13). There are 3,710 Micro-PVM_3 nuclei. Per donor, the median is 126 in AD, 134 in bvFTD, 54.5 in PSP, and 33.5 in controls. Four donors have fewer than 5: R24_P2730 (bvFTD, 0), R24_P2603 (PSP, 0), R24_P2423 (PSP, 2), and R24_P2737 (control, 4). No donor has 5 to 9, so a floor of 5 or 10 keeps the same donors. The donor floor is 5 Micro-PVM_3 nuclei for Section 3 and the diagnosis contrasts: 33 donors (10 AD, 8 bvFTD, 8 PSP, and 7 controls).
2. Depth of the PVM_3 arm (R13). The arm has 3,417 nuclei (92.1% of Micro-PVM_3), and the HKG8 arm 1,844. Per arm nucleus, the median is 3,047 UMIs, 1,742 genes, 42 MT UMIs, and 6 MT-CO3 UMIs. Of the arm nuclei, 30.0% have 1 to 3 MT-CO3 UMIs. The expected r_s at this depth is read from SEA-AD Stats_15H: with SEA-AD PVM_3 nuclei downsampled to 3,000 UMIs, the median r_s of MT-CO3 over the eight partners is 0.630 (0.626 to 0.632). The depth row applies, because the Rexach arm has a higher MT fraction than SEA-AD (median 1.41% against 0.57%), not the lower fraction of the depth plus MT row. It sets no threshold.
3. The 500-gene floor and the mappings (R12 and R13). The floor removed 5,418 of 228,671 insula nuclei (2.4%), and 392 of 16,976 of the authors' microglia (2.3%). The two mappings agree on 94.9% of nuclei at subclass and 44.7% at supertype. Of the 3,710 Micro-PVM_3 nuclei by Deep Generative Mapping, 1,514 (40.8%) are Micro-PVM_3 by Hierarchical Mapping, and most of the rest are Micro-PVM_2_3-SEAAD (1,162) or Micro-PVM_2 (814). All 3,710 are microglial cell in the authors' labels. Deep Generative Mapping assigns 7,793 nuclei to Sst Chodl. The authors call 7,791 of them inhibitory interneurons, their median gene count is 4,565, and Hierarchical Mapping spreads them over other neuron subclasses. This does not touch the arms.
4. The maximum %MT per nucleus (R11 and R13). Among the authors' insula microglia, the maximum is 8.05% and 9 of 16,976 nuclei are above 8%, which confirms the authors' cut. In the PVM_3 arm, the median %MT is 1.41%, the maximum 8.03%, and 3 nuclei are above 8%.
5. The taxonomy names Micro-PVM_3-SEAAD and Micro-PVM_1 appear exactly as in the SEA-AD Multiregion release, in both mappings (R13).

## 3. Do the eight DFGs replicate?

- Estimate. `cepip_screen(goi = "MT-CO3", min_pct = 0, min_mean = 0)` on all Micro-PVM_3 nuclei of the 33 donors with >= 5 Micro-PVM_3 nuclei (zeros included; no subsample). The GOI gate is set to 0 because the GOI is fixed in advance. MT-CO3's own gate values (share expressing and mean over expressing nuclei) are reported. For each of the eight DFGs: r_s and tau_b with MT-CO3.
- Tiers. Tier 1 is r_s >= 0.6 and tau_b >= 0.5 (SEA-AD). Tier 2 is r_s >= 0.5 and tau_b >= 0.4 (CEP-IP, CMPB 2026). The screen is run once at each tier, because tau_b is computed only for pairs that reach the r_s threshold. The thresholds are not lowered below tier 2.
- Why five. In the SEA-AD pooled PVM_3 screens (Supplementary Table S3), MT-CO3 had 7 of the eight as partners in NeoC and 5 in NonNeoC at tier 1. The SEA-AD screens used 2,500-nucleus subsamples. Here all nuclei are used.
- Reading.
  - Five or more of the eight pass tier 1: the module replicates at the SEA-AD thresholds.
  - Otherwise, five or more pass tier 2: the module replicates at the CEP-IP lower thresholds.
  - Otherwise: the module does not replicate at either tier in this dataset.
- What a pass shows. A pass is expected in most datasets, because the 13 MT genes are transcribed together from mtDNA. A pass shows that the composite behaves as in SEA-AD. The evidence for the module itself comes from Section 4, the diagnosis contrasts, and the HKG8 control.
- Section 4 runs on the fixed composite whatever this reading is.

## 4. IP gate

- The IPRS tiers and rules of CEP-IP (CMPB 2026, Section 2.8) apply to each arm, as in External validation 1.
  - Very Strong or Strong (IPRS >= 0.60): the diagnosis contrasts of Section 6 are run.
  - Moderate or Weak (0.20 to below 0.60), or any of C1 to C5 equal to 0: the arm is flagged for review, as in CMPB 2026 Section 2.8. That section gives no numeric rule. It asks the review to confirm post-IP TREP enrichment above the main GAM curve. Here the review is the mosaic comparison of CMPB 2026 Section 3.12 (Fig. 9A), written as a count.
    - The shares of TREP cells above the main GAM curve in pre-IP and in post-IP, compared by Fisher's exact test.
    - Confirmed if the post-IP share is higher, with p < 0.05. Then the diagnosis contrasts are run, and their results are labelled as conditional on the review.
    - Not confirmed: the diagnosis contrasts are not run for that arm.
    - The R(x) plot, the main GAM plot, both shares, and the test are reported.
    - The test is lenient, because the IP is placed where R(x) of the TREP residuals turns positive. It guards against an IP that rests on few cells or on noise in R(x).
  - Failed (< 0.20), no IP found, or fewer than 50 GOI-positive cells: the diagnosis contrasts are not run for that arm.

## 5. Association with diagnosis

- The rule of External validation 1 needs at least 10 donors in each group. The insula has 8 controls (R11), and 7 after the donor floor (R13), so no diagnosis test is pre-specified. The contrasts are reported as observed in Section 6.

## 6. Reported as observed, without a reading

- Screen. The within-positive r_s and tau_b of MT-CO3 with cDFG; the agnostic screen at both tiers; and the rank of MT-CO3 among all GOIs.
- Pipeline. SKCS-core, OriRef, Recovered, and gate counts; Validation1 DE and label agreement; the same-cell control (DE of cDFG against cHKG9 on MT-CO3); and the share of Module post-IP cells that are Micro-PVM_3 in both mappings.
- Known-filter test. The share of Module post-IP cells and of other PVM_3 arm cells removed by %MT cutoffs of 1% and 5%. The 10% cutoff of External validation 1 is not used, because no nucleus is above the authors' 8% cut.
- Diagnosis contrasts. The predictor is Module post-IP% = PVM_3 Module post-IP cells / all Micro-PVM_3 nuclei of the donor x 100, from `cepip_donor()` with `n_total` set to all Micro-PVM_3 nuclei of all donors. The donors are the 33 above the floor of count check 1. The model is `cepip_associate()` with age at death and sex as covariates (sex F = 0 and M = 1), as in External validation 1.
  - All three diagnoses combined against controls (diagnosis control = 0 and disease = 1), with `min_donors = 20`.
  - AD, bvFTD, and PSP each against controls, with `min_donors = 15`, because the contrasts have 17 (AD) and 15 (bvFTD and PSP) donors after the donor floor.
  - HKG8 is fitted with the same models and reported beside, with the PVM_3 minus HKG8 coefficient difference.
- DEGs by `cepip_deg()` with the S8 thresholds.

## Dated notes

- 2026-10-25, before any statistic of Sections 3 and 4 was computed. Clarification of Section 4: the mosaic comparison splits pre-IP from post-IP on the GOI axis (x < IP_x), as CMPB 2026 Section 3.12 does. The split on the cDFG axis (y < IP_y) is reported beside it. The test is two-sided, and the review is confirmed only when the post-IP share is the higher one.
