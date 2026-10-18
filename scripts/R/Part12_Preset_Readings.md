# Part 12 rerun: pre-specified readings

Status: Figure 15 agreed 2026-09-21. Figure 16 agreed 2026-09-22. Figure 17 agreed 2026-09-22. Figure 18 agreed 2026-09-22. Figures 16 to 18 are added one battery at a time, before that battery runs. Each panel script names its reading at the top of its log. A reading is changed only with a dated note in this file stating the reason, and the Results report the change.

## Rules common to every panel

1. **Arms.** M4 is the primary arm and carries the reading. PVM_3 is the reference and PVM_1 is reported. HKG8 is the control. If M4 fails a floor on a panel, PVM_3 carries the reading and the Results say so.
2. **Significance.** A result is significant at FDR < 0.05, BH within the family named for the panel. The figure marks * at FDR < 0.05 and ** at FDR < 0.01. The text gives the exact FDR as `%.2e`, and a value at the resampling floor is written with "<=".
3. **Donor bootstrap.** B = 10,000, 95% percentile interval, p = min(1, 2(k+1)/(B+1)). The Gerrits panels resample donors and keep their samples together.
4. **Floors.** 20 cells per group per donor, 10 donors per test, 3 units with a finite value for the within-donor level. A test that fails a floor is reported as untestable, never dropped.
5. **Control reading.** A shift that appears on HKG8 with the same sign and FDR < 0.05 is a property of the construction or the preparation, not of the Module. A shift on the Module arms with none on HKG8 is attributed to the Module.
6. **Direction.** Every reading states the sign that supports each account. A result of the opposite sign is reported as its own outcome.
7. **Old results.** A result that differs from the archived run is not investigated or reconciled; the rulebook changed the methods on purpose. Old and new values meet only in the comparison table (D20). The one exception is a structural check (a cell, unit or donor count that depends only on the frozen labels) that fails: that is treated as a possible code defect and fixed before the panel's result is read.
   - **Note, 2026-09-22, after every figure had run.** The comparison table (D20) is dropped; old and new values are not compared anywhere and results are reported as they stand. No reading changes.

## Figure 15. Battery 1: Is it contamination?

### 15A. Does contamination explain the clinical association?

- **Model.** The Figure 12A clinical model: cov4 OLS on the donor Module post-IP%, region sets split, 9 outcomes, all Figure 12A patients. Refit with donor-mean supertype non-Module DecontX contamination added as a covariate.
- **Estimate.** For each arm, region set and outcome: the Module coefficient's FDR before and after adjustment, and retention = dR2 after / dR2 before.
- **Family.** Arm x region set x adjustment state, over the 9 outcomes.
- **Rows read.** The rows significant before adjustment. The primary rows are the M4 cognition rows (cognitive status, MMSE, CASI in each region set).
- **Reading.**
  - Every primary row significant before adjustment is still significant after: contamination does not account for the cognition association.
  - One or more primary rows lose significance, with retention below 0.5: contamination accounts for part of the association on those rows.
  - Rows that lose significance with retention of 0.5 or more are reported as borderline, because the loss is a shift across the threshold without a large drop in dR2.
- **Control.** HKG8 rows significant before adjustment are read the same way and reported beside M4.

### 15B. Do doublets drive the clinical association?

- **Model.** As 15A, refitted after Module post-IP nuclei with a doublet score at or above 0.15, 0.20 and 0.30 are removed and Module post-IP% is recomputed. The share of Module post-IP nuclei removed at each threshold is reported per arm.
- **Family.** Arm x region set x threshold, over the 9 outcomes.
- **Reading.**
  - Every primary row significant before is still significant at all three thresholds: doublets do not drive the association.
  - A primary row loses significance at 0.30, where only high-scoring nuclei are removed: doublets contribute to that row.
  - A primary row loses significance only at 0.15, with at least a quarter of Module post-IP nuclei removed: reported as a loss of cells and power, not attributed to doublets.

### 15C. Are Module post-IP cells more contaminated?

- **Estimate.** Per donor, the median DecontX contamination of each cell group; per contrast, the median over donors of the donor-level difference. Three contrasts: Module post-IP minus Module pre-IP (primary), Module post-IP minus non-Module, Module pre-IP minus non-Module. NeoC and NonNeoC pooled within donor. All four arms.
- **Family.** Within arm, over the 3 contrasts.
- **Reading.**
  - The primary contrast is at or below zero, or FDR >= 0.05: no evidence that Module post-IP cells are more contaminated.
  - The primary contrast is positive with FDR < 0.05 on M4, and HKG8 shows no such shift: Module post-IP cells carry more contamination, and 15D and 15E to 15G decide whether it drives the co-expression.
  - The same positive shift on HKG8: the difference comes from the post-IP construction and is not specific to the Module.

### 15D. Does co-expression rise with contamination?

- **Estimate.** Module cells of each arm, split into contamination tertiles within donor (D25) and pooled per tertile across donors. Spearman rs of MT-CO3 with the cDFG composite in each tertile (ACTB with cHKG8 on HKG8). Primary estimate: rs(top tertile) minus rs(bottom tertile).
- **Family.** The 4 arms, one test each.
- **Reading.**
  - The difference is at or below zero, or FDR >= 0.05: the co-expression does not rise with contamination, which an ambient source would predict.
  - The difference is positive with FDR < 0.05 on M4 and not on HKG8: contamination contributes to the co-expression.
  - Positive on both M4 and HKG8: the rise is a general effect of contamination on any composite and is not specific to the mitochondrial axis.

### 15E. Is the ambient pool more mitochondrial than the nuclei? (Gerrits)

- **Estimate.** Per sample, log of the ratio of ambient %MT (barcodes of 1 to 100 UMIs) to nuclear %MT (D16 cells). Median over samples, donor bootstrap.
- **Family.** One test.
- **Reading.**
  - Median log ratio > 0 with FDR < 0.05: the ambient pool is enriched for mitochondrial reads, so an ambient source of MT reads exists; 15F and 15G test whether it produces the co-expression.
  - At or below zero, or FDR >= 0.05: no mitochondrial enrichment in the ambient pool.

### 15F. Does co-expression survive removing ambient RNA and depth? (Gerrits)

- **Estimate.** On one cell set (D16) and one population (D27): per sample, the median over the 8 DFG partners of rs(MT-CO3, partner), raw and partial (soup index and log depth partialled out). Primary estimate: median over samples of the partial value. Secondary: the same per partner, and partial divided by raw, reported without a test.
- **Family.** Primary: raw and partial, 2 tests. Secondary: 8 partners, within raw and within partial.
- **Reading.**
  - Partial median > 0 with FDR < 0.05: the co-expression is present after ambient load and depth are removed.
  - Partial median at or below zero, or FDR >= 0.05, while raw is positive: the raw co-expression is accounted for by ambient load and depth.
  - Partial divided by raw is reported as the retained share, with no threshold.

### 15G. Does co-expression change with ambient load? (Gerrits)

- **Estimate.** Soup-index tertiles within sample. Per sample, the median over the 8 partners of rs(top tertile) minus rs(bottom tertile). Median over samples, donor bootstrap. Secondary: per partner.
- **Family.** Primary: one test. Secondary: 8 partners.
- **Reading.**
  - Difference < 0 with FDR < 0.05: co-expression is weaker in nuclei with more ambient RNA, so the ambient RNA dilutes it and does not produce it.
  - Difference > 0 with FDR < 0.05: co-expression is stronger with more ambient RNA, so ambient RNA contributes to it.
  - FDR >= 0.05: co-expression does not change with ambient load.

### 15H. Do depth and mitochondrial content explain the gap to Gerrits? (simulation)

- **Population.** The PVM_3 arm population (Micro-PVM_3 cells with MT-CO3 > 0), as the Methods state and as 16D uses.
- **Simulation.** Binomial thinning of each cell's counts to the nine target depths, with and without matching the mitochondrial fraction to Gerrits', 20 repeats. The median rs over repeats and its range are drawn. No test.
- **Reference.** The interquartile range of the Gerrits per-sample raw rs, read from the 15F stats file, so 15H runs after 15EG. Gerrits' depth and mitochondrial fraction come from the 06 output.
- **Reading.**
  - At Gerrits' depth with matched mitochondrial content, the SEA-AD rs falls inside the Gerrits interquartile range: depth and mitochondrial content account for the gap.
  - It stays above the range: depth and mitochondrial content do not account for the whole gap, and the remainder is attributed to protocol differences, not tested here.
  - Matching depth alone brings it inside the range: depth alone accounts for the gap.

### Supplementary Figure 10. Is the extra material the cell's own cytoplasm or another lineage's?

- **Estimate.** Paired within-donor differences, Module post-IP minus Module pre-IP, in four metrics: own-ER share, ribosomal share, NEAT1 share and the combined other-lineage index (12c per-dataset guard, D21). NeoC and NonNeoC pooled within donor. All four arms.
- **Family.** Within arm, over the 4 metrics.
- **Reading.**
  - Own-ER and ribosomal shares higher with FDR < 0.05, and the other-lineage index not higher: the extra material is the cell's own cytoplasm.
  - Other-lineage index higher with FDR < 0.05, and own-ER and ribosomal shares not higher: the extra material comes from other cells.
  - Both higher: both sources contribute; the sizes are reported side by side.
  - NEAT1 share lower with FDR < 0.05 is read as loss of nuclear RNA and reported beside the cytoplasmic result.
  - The same shifts on HKG8 are read by common rule 5.

## Figure 16. Battery 2: Is it the measurement?

Status: agreed 2026-09-22, with choices A1, A2, B1, C1, D1, E1 and F1 as recommended.

### 16A. Is Module post-IP% donor-structured?

- **Estimate.** Across units, the share of variance explained by donor (donor entered first, then region), for five quantities per arm: Module post-IP%, its TQR residual, log1p %MT, contamination, and the matched non-MT share (`ctrl_pct_null`). Null: donor labels shuffled within region, 10,000 times; p = (k+1)/(R+1). This is the section 4 exception: the question is whether donor structure exceeds what region imbalance produces, which a bootstrap interval cannot answer.
- **Family.** Within arm, over the 5 quantities (4 families of 5; the archived run used one family of 20).
- **Reading.**
  - Module post-IP% donor share above the null with FDR < 0.05 on M4: Module post-IP% is donor-structured.
  - The TQR residual also above the null: the donor structure is not explained by the four tissue-quality fields.
  - The same result on HKG8 is read by common rule 5: donor structure comes from the construction or the preparation, not from the Module.
  - %MT, contamination and `ctrl_pct_null` are reported beside the Module as comparators, with no reading of their own.
- **No equivalence test.** The archived 55% equivalence bound on HKG8 is dropped; rule 5 reads HKG8.

### 16B. Does a donor keep its value across regions?

- **Estimate.** For donors with at least two units, each split puts a random half of the donor's units in each half; the donor's mean in each half is taken; rs across donors between the halves. 200 random splits, fixed from the step seed; the estimate is the median rs over the splits. Donor bootstrap on that median (the same splits in every resample), 10-donor floor. All four arms.
- **Quantities.** Module post-IP%, its TQR residual and log1p %MT.
- **Family.** Within arm, over the quantities.
- **Reading.**
  - Module post-IP% split-half rs > 0 with FDR < 0.05 on M4: a donor keeps its Module post-IP% across regions.
  - HKG8 is read by rule 5. An arm below the donor floor is untestable, which may apply to M4 because of few units per donor.

### 16C. Do two arms label the same cells?

- **Estimate.** Cohen's kappa for Module membership on cells shared with the PVM_3 reference: M4 against PVM_3 (positive control) and HKG8 against PVM_3 (negative control). Pooled over regions with at least 20 shared cells; donor bootstrap interval. PVM_3 is the reference and PVM_1 shares no cells with it (D24).
- **Family.** None: two estimates read against fixed thresholds.
- **Reading.** Thresholds kept from the archived run.
  - M4 kappa above 0.60 and HKG8 |kappa| below 0.20: the Module is reproduced on a shared axis and not by the construction.
  - Otherwise the failing arm is named.
- The thresholds are read on the point estimate, as in the archived run; the interval is reported beside it. Per-region-set rows (NeoC, NonNeoC) are reported and drawn as in the archived figure.

### 16D. Does the MT-CO3 to cDFG relationship hold inside single donors?

- **Estimate.** On the PVM_3 arm population (MT-CO3-positive Micro-PVM_3 cells), per donor with at least 20 cells, rs of MT-CO3 with the cDFG composite (log1p CP10K), raw and with log depth partialled out. The median over donors, with a donor-bootstrap interval. Beside it, the pooled rs over all cells and the between-donor rs of donor means. Control pair: ACTB with the cHKG8 composite on the same cells.
- **Family.** No test; the median is read against fixed marks.
- **Reading.** Marks kept from the archived run, on the raw within-donor median.
  - At or above 0.75: the relationship holds inside donors and is not produced by pooling.
  - 0.40 to 0.75: it holds partly inside donors.
  - Below 0.40: the pooled relationship is mostly between donors.

### 16E. Is the offset share equal between the branches?

- **Estimate.** Per donor and branch (Module post-IP, Module pre-IP), OLS of MT UMIs on total UMIs; the offset share is the median intercept over donors divided by the fitted MT UMIs at the median depth. The contrast is post-IP minus pre-IP in percentage points, with a paired donor bootstrap (a donor's two fits stay together). 20 cells per donor and branch, 10 donors. No depth filter (D8). All four arms.
- **Family.** The 4 arms.
- **Reading.** Margin 10 pp, read on the 90% interval (the section 4 exception); the 95% interval gives the non-zero test.
  - The whole 90% interval inside +/-10 pp and the 95% interval covering zero: the offset share is the same in both branches, within the margin.
  - Inside the margin but non-zero: a difference exists and is smaller than 10 pp.
  - Non-zero and not inside the margin: the post-IP branch carries a different offset share, which would mean the Module selects the part of the MT signal that does not scale with depth.
  - Neither non-zero nor inside the margin: underpowered.
- **Note, 2026-09-22.** The reading said both that the family is the 4 arms and that the 95% interval gives the non-zero test. These disagree when the interval excludes zero but the FDR over the arms does not pass. Resolved in favour of common rule 2: non-zero means FDR < 0.05, BH over the arms; the margin still reads the 90% interval. This was decided after the first run, in which M4 had p 0.0412 (95% interval +1.0 to +37.2) and FDR 0.12; the Results report the change.

### 16F. Does Module post-IP% track the library intronic rate?

- **Estimate.** rs of Module post-IP% with the library intronic rate at two levels: between donors (donor means over units) and within donor (units centred on their donor mean; donors with at least 3 units). Comparators, as in the archived run: log1p %MT, log cell count and log median depth. Contrasts on the same donor draws: |Module rs| minus |%MT rs|, Module minus depth, and the arm minus HKG8. Donor bootstrap; donors drawn once from all arms together, so the arm-minus-HKG8 contrast keeps the donor overlap. All four arms.
- **Family.** Within arm and level: the rs rows, then the contrasts as their own family.
- **Reading.**
  - Module rs non-zero with FDR < 0.05 at a level: Module post-IP% tracks the intronic rate at that level.
  - The arm minus HKG8 contrast non-zero: the association is specific to the Module; otherwise specificity is not established.
  - The Module minus %MT contrast is reported beside it, with no reading of its own.

### 16G. Do sort purity or library construction order it?

- **Estimate.** As 16F, for each library field with no more than 25% missing: antisense rate, NeuN-positive fraction on FANS, share of cDNA above 400 bp, PCR cycles and sequencing saturation. log median depth is carried as a comparator, with the Module minus depth contrast.
- **Family.** As 16F, within arm and level over the fields.
- **Reading.** Four conditions, kept from the archived run. A field orders Module post-IP% only if all four hold within donor:
  1. the Module rs is non-zero after BH;
  2. the Module minus depth contrast excludes zero;
  3. the arm minus HKG8 contrast excludes zero;
  4. the Module rs is not matched by %MT: the 95% interval of |Module rs| minus |%MT rs| is not wholly below zero.
  - The between-donor rows are reported beside, with no condition of their own.

## Figure 17. Battery 3: Is Module post-IP% just %MT, for histology and clinical outcomes?

Status: agreed 2026-09-22, with choices M1, A1, B1, C1, D1, E1 and S1 as recommended.

### Definitions shared by 17A to 17E and Supplementary Figure 12

- **Units.** The step 05 unit table (arm, dataset, donor, region; at least 10 cells of the arm). NeoC and NonNeoC units are pooled within donor for 17A, 17B, 17C and Supplementary Figure 12, because a donor's regions span both region sets and neither set alone gives 10 donors with 3 regions on M4 or HKG8.
- **Quantities.** Module post-IP% (`module_post_pct`, supertype denominator, D28); %MT as `mito_log1p` (D5); cytoplasmic share (`cyto_share`); contamination (unit mean DecontX); log cell count (`log_n_cells`). All five are computed for every measure; each panel draws its own set.
- **Neuropathology measures (choice M1).** `sea-ad-quantitative-neuropathology-063026.xlsx`, sheet `Travaglini 2026 (multi-region)`. The 28 stain measures of that sheet: its 32 measure columns less the four hematoxylin nucleus measures, which are a counterstain and belong to no stain family. The 28 names are listed in the config (version 14) and read by exact name. Stain families, in drawing order: Amyloid 6E10 (9, including the two 6E10-IBA1 colocalisation measures), Tau AT8 (3), TDP43 (3), Synuclein (3), Astro GFAP (1), Microglia IBA1 (6), Neuron NeuN (3).
- **CPS scores.** `sea-ad_global_and_local_cps_050126.xlsx`, sheet `Travaglini 2026`: CPS_Local, CPS_Local_ABeta, CPS_Local_pTau (region-resolved) and CPS_Global, CPS_Global_ABeta, CPS_Global_pTau (one value per donor; checked constant within donor).
- **Joins.** Histology and CPS region names are mapped to snRNA-seq region names by the config region map only. The neuropathology sheet has one analysis region (Grey matter) and one row per donor and region (checked 2026-09-22), so no rows are averaged; the join stops if either no longer holds. ESOC and CaH have no sequencing and are left out; LEC units have no histology or CPS row. The log lists every unmapped region name and the share of units joined per arm.
- **Levels.**
  - Within donor: donors with at least 3 units finite on both variables; each variable centred on its donor mean (`p12_centre_within`); rs over the pooled centred values.
  - Between donors: donor means over units of both variables; for a global CPS score, the donor's value.
- **Inference.** Common rule 3: donors resampled, B = 10,000, indices drawn once per arm in the main process. Floor: 10 donors per test (common rule 4). A test below the floor is drawn as untestable.
- **Note, 2026-09-22, before any run.** "Once per arm" is made exact: indices are drawn once per arm and frame, so each resample draws as many donors as the frame has. The frames are the arm's histology-joined units (17A to 17C, S12 and the 17E histology-unit rows) and the arm's donor table per region set (the 17E clinical-grain rows). The 10-donor floor applies to the observed data; a resample needs 3 distinct donors.
- **Readability.** A two-predictor model is unreadable when the VIF of Module post-IP% and %MT, 1/(1 - r^2) from `p12_vif` on the two predictor columns of the fitted frame, is 5 or more. The VIF is taken on the full sample, not per resample.

### 17A. Does microglial shape track Module post-IP%, or just %MT?

- **Estimate.** Within-donor rs of each of the 6 Microglia IBA1 measures with Module post-IP%, %MT and the cytoplasmic share. All four arms. The between-donor rows are computed and written to the stats file but not drawn.
- **Family (choice A1).** Within arm, level and quantity, over the 6 IBA1 measures (the archived run used the 26 measures).
- **Reading.** The activation account predicts a positive rs with activated IBA1+ cells per area, and a negative rs with inactivated IBA1+ cells per area and with process area and process length per cell. IBA1+ cells per area and IBA1+ area (%) are abundance measures with no sign set in advance.
  - Module rs non-zero with FDR < 0.05 in the predicted direction on M4 for one or more shape measures: Module post-IP% tracks microglial activation within donor.
  - %MT rs of the same sign on the same measures: the marginal association is shared with %MT; 17C tests whether the two can be separated.
  - Cytoplasmic share rs of the same sign on the same measures: consistent with more cytoplasm carried over from activated cells; reported beside, with no reading of its own.
  - A result of the opposite sign is reported as its own outcome (common rule 6).
  - A null on an abundance measure is not evidence against the Module, which is a fraction.
- **Control.** HKG8 read by common rule 5.

### 17B. Does nearby pathology track Module post-IP%, or just %MT?

- **Estimate.** rs of each CPS score with each quantity. Within donor: the three local scores (the global scores are constant within donor). Between donors: all six. All four arms.
- **Quantities drawn (choice B1).** Module post-IP%, %MT, cytoplasmic share and log cell count, as the archived panel plotted them; the key names the cytoplasmic share correctly (the archived key said contamination). Contamination is in the stats file, not drawn.
- **Family.** Within arm, level and quantity, over the CPS scores of that level (3 within donor, 6 between donors).
- **Reading.** A positive rs supports the account that Module post-IP% follows pathology.
  - Module rs > 0 with FDR < 0.05 on M4 with CPS_Local within donor: the Module's regional variation follows regional severity.
  - Within donor not significant and between donors > 0 with FDR < 0.05: Module post-IP% is a donor-level marker of disease stage.
  - A sign reversal between the two levels is reported as its own outcome.
  - Log cell count or cytoplasmic share tracking the same score with the same sign: the association may reflect tissue recovery; reported beside the Module.
  - %MT is reported beside; the separation of the two is not tested on CPS.
- **Control.** HKG8 read by common rule 5.
- **Change.** The archived within-donor permutation null is dropped; the donor bootstrap gives p (common rule 3).

### 17C. Can Module post-IP% and %MT be separated on histology?

- **Model.** Per arm, level and measure: rank(measure) ~ rank(Module post-IP%) + rank(%MT), standardised, no covariates. Within donor, every variable is centred on its donor before ranking. Between donors, donor means. 28 measures, 2 levels, all four arms. %MT in its log form only; the archived raw and rank variants are dropped (D5).
- **Inference.** The two coefficients on the same donor resamples; p by common rule 3. Sign flip (the coefficient's sign differs from the marginal rs) is a column.
- **Family (choice C1).** Within arm, level and predictor, over the 28 measures. The four stain blocks (Amyloid 6E10; Tau and TDP43; Microglia IBA1 and astrocyte GFAP; NeuN and synuclein) are for drawing only.
- **Verdict per row.** Unreadable (VIF >= 5); otherwise Module only, %MT only, Both or Neither, from the two FDRs at 0.05.
- **Reading.** On M4 within donor:
  - One or more Module-only rows: Module post-IP% carries histological information beyond %MT on those measures.
  - No Module-only rows and one or more %MT-only rows: the histology association is carried by %MT. A %MT-only row with a Module sign flip is read as shared variance the model could not allocate, not as absence of a Module association.
  - All rows Neither or Unreadable: the data cannot separate the two on histology. This is not evidence for either.
  - The between-donor rows are read the same way and reported beside.
- **Control.** HKG8, where the two predictors are expected to be weakly correlated, is read by common rule 5.
- **Title.** The archived panel title said cov4; the histology model carries no covariates, and the title says so.

### 17D. Can Module post-IP% and %MT be separated on the clinical outcomes?

- **Model.** The 15A clinical model with both predictors: outcome ~ cov4 + Module post-IP% + %MT, donor table, region sets fitted separately, the 9 outcomes of Figure 12, APOE e4 dose leaving the covariates when it is the outcome. Same donors as 15A. OLS p (section 4 exception 1, as 15A).
- **Family (choice D1).** Within arm, region set and predictor, over the 9 outcomes, as 15A. Outcome classes (cognition 3, pathology 5, APOE 1) are for drawing only.
- **Verdict per row.** As 17C, with the VIF from the two predictor columns of the donor frame.
- **Reading.** The primary rows are the M4 cognition rows in each region set, as 15A.
  - A primary row Module only: Module post-IP% carries clinical information beyond %MT on that outcome.
  - A primary row %MT only: that clinical association is carried by %MT, with the sign-flip reading of 17C.
  - Neither or Unreadable: the two cannot be separated on that outcome.
- **Control.** HKG8 read by common rule 5. HKG8 now has its own %MT column (D5), so its verdicts change from the archived run.

### 17E. Why are so many models unreadable?

- **Estimate (choice E1).** |r| between Module post-IP% and %MT and the VIF, in exactly the frame each model uses:
  - within donor (histology units): Spearman, NeoC and NonNeoC pooled, as 17C;
  - between donors (histology units): Spearman on donor means, pooled, as 17C;
  - between donors (clinical grain): Pearson on the donor table, region sets split, as 17D.
  Donor-bootstrap intervals on |r|. No test (section 4 exception 5).
- **Reading.** Descriptive, read against the VIF 5 line (|r| 0.894).
  - Rows at or above the line on the real arms and below it on HKG8: collinearity is the cause of the Unreadable verdicts in 17C and 17D.
  - Rows below the line with Neither verdicts: the models are readable and the two predictors are not separable at this sample size.

### Supplementary Figure 12. How do the stain families track each quantity?

- **Estimate.** Per arm, level, quantity and stain family: the median over the family's measures of the rs of 17A's design (panel A). Panel B: the family median of |Module rs| minus |%MT rs|, on the same resamples. Measures below the floor leave the family median; the drawn count is (tested/total).
- **Drawn (choice S1).** Within donor, quantities Module post-IP%, %MT, contamination and log cell count, as the archived panel. Between-donor rows and the cytoplasmic share are in the stats file.
- **Family.** Within arm, level and quantity over the 7 stain families (panel A); within arm and level over the 7 families (panel B).
- **Reading.**
  - Panel B < 0 with FDR < 0.05: %MT tracks that stain family more strongly than Module post-IP%.
  - Panel B > 0 with FDR < 0.05: Module post-IP% tracks it more strongly.
  - A Module rs smaller than the %MT rs by about 0.80-fold is what a cut version of the same axis predicts; this is reported with no threshold.

## Figure 18. Battery 4: Can %MT replace the Module?

Status: agreed 2026-09-22, with choices A1, A2, B1, C1, D1, D2, E1 and S1 as recommended.

### Figure 18 definitions, shared by 18A to 18E

- **Arms and region sets.** All four arms. 18A, 18B and 18C keep NeoC and NonNeoC apart (18A is a clinical model, D11; 18B and 18C work on cells within a region set). 18D and 18E use the step 05 unit table with the region sets pooled within donor, as 16A and 16B do.
- **%MT.** The 13-gene %MT of the arm's own cells (D5): per cell for 18B and 18C, the donor mean in log1p form (`mito_log1p`) for 18A.
- **Module cells.** The arm's Module cells (step 04 `is_module`) against its non-Module labelled cells.
- **Inference.** Common rule 3 for the 18C test. 18A keeps OLS p (section 4 exception 1). 18B, 18D and 18E carry no test (exceptions 4 and 5). Simulation draws are made from the step seed per simulation cell, so the result does not depend on the worker count.

### 18A. Can a nuclear-encoded mitochondrial panel replace the Module?

- **Panels.** The 15 nuclear-encoded panels of the config (all nuclear-encoded mito, OXPHOS core, Complex I to V, assembly factors, biogenesis, mitoribosome, protein import, fission and fusion, mitophagy, TCA cycle, glycolysis) and the 300-gene matched reference (the null), from `02_matched_sets.csv`. Each panel keeps the genes detected in at least 10% of the arm's cells in that region set, and needs at least 3 such genes.
- **Donor score.** Per cell, the mean over the panel's genes of log1p(count / total UMIs x 10,000). Per donor, the median over the donor's labelled cells of the arm in that region set (choice A1: donors with at least 20 cells, common rule 4; the archived run used 5).
- **Model.** Per arm, region set, panel and second predictor (Module post-IP% or %MT): outcome ~ cov4 + second predictor + panel score, all standardised, the 9 outcomes of Figure 12, APOE e4 dose leaving the covariates when it is the outcome. 10-donor floor. OLS p.
- **Readability (choice A2).** VIF of the two predictors from the two predictor columns, as 17D; 5 or more is unreadable.
- **Family.** Within arm, region set, panel and predictor, over the 9 outcomes, each coefficient separately, as 17D.
- **Verdict.** Module only (or %MT only), NE panel only, Both, Neither, Unreadable, Untestable, from the two FDRs at 0.05.
- **Reading.** The primary rows are the M4 cognition rows (3 outcomes in each region set) on the 15 real panels.
  - Module only or Both on the primary rows, and no NE panel only row: no nuclear-encoded panel replaces the Module; the Module carries the association beyond each panel.
  - NE panel only on a primary row: that panel carries the association on that outcome and the Module stands in for it there.
  - The matched reference is the null. A real panel that wins no more rows than the reference has shown nothing. Reference wins are read as expression level, not mitochondrial biology.
  - The %MT tally is read the same way and reported beside the Module tally.
- **Control.** HKG8 read by common rule 5; its %MT is its own column (D5).

### 18B. Does any %MT threshold reproduce the Module cells? (descriptive)

- **Estimate.** Per arm and region set, cells pooled over donors: the densities of log10(%MT + 0.01) for Module and non-Module cells, the AUC of %MT for Module membership, the best single threshold by Youden's J, and Cohen's kappa at that threshold. No test (D19); 18C is the test.
- **Drawn (choice B1).** Both region sets, NeoC above NonNeoC. This takes in the NonNeoC densities that sat in Supplementary Figure 15.
- **Reading.** Bands kept from the archived run, on the real arms.
  - AUC below 0.60: %MT does not recover the Module cells.
  - AUC 0.60 to 0.80: %MT carries part of the partition, and cells at the same %MT fall on both sides.
  - AUC 0.80 or more: the Module is close to a %MT threshold.
  - Kappa at the best threshold is reported beside the AUC, with no band of its own.

### 18C. Does %MT distinguish Module cells within a donor?

- **Estimate.** Per arm, region set and donor with at least 20 Module and 20 non-Module cells (choice C1, common rule 4; the archived run used 10): the AUC of %MT for Module membership. The median over donors, with a donor-bootstrap interval; p for the median against 0.5 by common rule 3.
- **Family.** Within region set, over the 4 arms.
- **Reading.** M4 carries the reading; if M4 is below the donor floor, PVM_3 carries it (common rule 1).
  - FDR >= 0.05: %MT does not distinguish Module from non-Module cells within a donor.
  - Median above 0.5 with FDR < 0.05: Module cells carry more %MT within a donor; the size is read against the 18B bands.
  - Median below 0.5 with FDR < 0.05 is reported as its own outcome.
- **Control.** HKG8 read by common rule 5.

### 18D. How many donors does each effect need? (simulation)

- **Design.** Per arm, donors are subsampled without replacement at 10, 15, 20, 25, 30, 40, 50, 60, 70 and 83 donors (up to the arm's count), 300 draws per point, with the unit Module post-IP% rescaled so the donor share is 1.0, 0.75 or 0.5 times the observed share (the archived rescaling of the donor deviation).
- **Rules.**
  - R1, the 16A rule: the donor share of the subsample exceeds the 95th percentile of its null (donor labels shuffled within region; 20 subsets x 25 shuffles per donor count, median of the 95th percentiles).
  - R2, the 16B rule: split-half rs over 200 random splits of each donor's units (donors with 2 or more units, 10-donor floor) is 0.50 or more.
  - R3, D22 (choice D1): MMSE score ~ cov4 + donor Module post-IP%, standardised, Module coefficient OLS p < 0.05, where the donor value is the mean of the donor's rescaled units in the subsample (region sets pooled) and cov4 and MMSE come from the donor table.
  - The share of draws in which R2 could be computed is drawn beside R2.
- **Estimate.** The detection rate per rule and donor count, and the smallest donor count reaching 80%.
- **Drawn (choice D2).** One row: M4 at 1.0x and 0.75x, PVM_3 at 1.0x and 0.75x, as the archived figure laid them out (the archived order was PVM_3 first).
- **Reading.** No test.
  - The donor requirement is the largest of the three crossings at 0.75x on M4 and on PVM_3; a rule that does not reach 80% needs more donors than the cohort has.
  - HKG8 R1 should not reach 80% at any donor count. If it does, the curves measure sample size and 18D is not read.
- **Closure check.** At 1.0x with every donor, the donor share equals the 16A estimate for Module post-IP% on the same arm.

### 18E. What does raising the region floor cost? (descriptive)

- **Estimate.** Per arm and minimum regions per donor, from 3 to the largest count any donor has: the share of the arm's donors kept, the donor share on the kept donors (3-donor floor, as archived), and the split-half rs by the 16B rule on the kept donors (choice E1: 10-donor floor, common rule 4; the archived run used 15).
- **Reading.** The region requirement is the smallest minimum beyond which neither the donor share nor the split-half rs rises by more than 0.05 while retention falls. Kept donors at a high minimum are the donors sampled in many regions, not a random subset, and the Results say so.

### Supplementary Figure 17. Donor-count curves per arm

- **Estimate.** 18D for each arm at 1.0x, 0.75x and 0.5x, one file per arm (choice S1). Supplementary Figures 14, 15 and 16 are not drawn (section 7: 14 and 16 dropped with their 12j and 12r content, 15 with 12m; its NonNeoC densities are in 18B).
- **Reading.** As 18D, per arm.
