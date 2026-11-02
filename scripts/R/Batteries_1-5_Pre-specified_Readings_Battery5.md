# Part 12, Battery 5: pre-specified readings

Question: do the eleven severely affected (SA) donors drive the between-donor and clinical results of Batteries 1 to 4?

Written 2026-09-25, before any Battery 1 to 4 statistic was computed without these donors. The all-donor results of Batteries 1 to 4 were known when this file was written; the without-eleven results were not. Every section of Part12_Fig19_Analysis.r names the readings it uses from this file before it computes anything.

## Scope and status

### Battery 5 definitions.

- Why a separate battery. All 84 donors are primary in Batteries 1 to 4, and the eleven SA donors were pre-specified as a sensitivity analysis of the clinical association only (Supplementary Figure 10). The Methods state that the between-donor results may depend on them. Battery 5 tests that statement on the panels exposed to a donor-level shift, instead of leaving it open.
- Change of decision, recorded. On 2026-09-23 no further reruns of Batteries 1 to 4 were planned. Battery 5 does not rerun them: it reads their stats files and inputs, reproduces the all-donor estimates it retests, and computes new statistics only without the eleven and under random removal. It is recorded in the migration summary as a dated entry.
- Files left unchanged. Part12_01_config.json (version 16), Part12_Preset_Readings.md, Part12_03_Config.r and the Figure 15 to 18 analysis files are not edited, so they stay as they were run. Battery 5 has its own configuration (Part12_Fig19_config.json) and this readings file.
- What is held fixed. Cell labels, the inflection point, the units, the Module post-IP% of every unit and donor, the %MT, the nuclear panel genes of 18A and every donor panel score are those of the all-donor run. Only donor rows are removed from each statistic. Battery 5 therefore asks whether the eleven donors carry the statistics, not whether they shaped the Module definition.
- The SA flag. The SEA-AD column "Severely Affected Donor" of obs_dataset.csv, every arm and region set. Exactly eleven donors must be flagged, the flag must be constant within donor, and every flagged donor must have ADNC_rank 3 (High). Otherwise section 26 fails.
- Frames. A frame is the donor set a statistic is drawn from: the arm's donors in the unit table (16A, 16B, 16F); the arm's donors joined to histology or CPS (17A, 17B); the arm's donors joined to histology (17C); the arm's donor-table rows in one region set (15A, 17D, 18A). Removing the eleven removes those of them present in the frame; k is that number.
- Without the eleven. Each retested statistic is recomputed on the frame less its SA donors, with the method, floors and multiple-testing families of the all-donor run, and a new resampling seed (Part12_Fig19_config.json). FDR is BH within the same family as in the all-donor run.
- Random removal (the null). In each frame, k donors are drawn without replacement from the frame's High-ADNC donors that are not SA, 1,000 times; the draws are made once in section 26 and written, so every later section uses the same draws. A contrast between two arms uses draw r of each arm's frame. Each draw recomputes the estimate only.
- Change. For a row with all-donor estimate e and sign s = sign(e), the change on removing a set is s x (e_removed - e). A negative change is attenuation toward zero. The SA change is compared with the 1,000 random-removal changes of the same row; p_null = (1 + number of random changes at or below the SA change) / (1 + number of finite random changes), one-sided, descriptive.
- Reproduction check. Before anything is removed, each section recomputes the all-donor estimates it retests and compares them with the Battery 1 to 4 stats files (estimate, units and donors; and the OLS p for 15A, 17D and 18A). A difference above 1e-8 fails the section.
- Row classes. Significant means FDR < 0.05; the all-donor FDR is the one in the Battery 1 to 4 stats file.
  - Holds: significant with and without the eleven, same sign.
  - Reversed: significant with and without, opposite sign.
  - Lost, beyond random removal: significant with, not without, and the SA change is below the 5th percentile of the random-removal changes.
  - Lost, within random removal: significant with, not without, and the SA change is at or above that percentile.
  - Lost, null untestable: significant with, not without, and the frame has fewer non-SA High-ADNC donors than k.
  - Untestable without the eleven: significant with, below the donor floor without.
  - Gained: not significant with, significant without.
  - Neither: significant in neither.
  - Untestable: untestable in the all-donor run.
- Panel reading. A panel holds without the eleven when none of its rows is Reversed or Lost, beyond random removal. Otherwise the panel is reported as depending on the eleven for that number of rows, named. Rows lost within random removal are reported as lost with power, not as dependence. The counts of every class are reported for every panel, whatever they are.
- Placement. Figure 19 of the Part 12 rerun files (the manuscript number is set at assembly), in the main text, whatever the outcome. Every row of every retested panel is in the Figure 19 stats and readout files.

### 19A. Unit-level donor structure: 16A, 16B and 16F between donors (section 27)

- 16A: donor share of variance per arm and quantity (five quantities), against 10,000 within-region shuffles of the donor labels; BH within arm.
- 16B: split-half rs per arm and quantity (three quantities). The halves are those of the all-donor run (reproduced from its seed), so without the eleven each remaining donor keeps its halves; donor bootstrap on the remaining donors; BH within arm.
- 16F between donors: rs of donor means with the library intronic rate for Module post-IP%, log1p %MT, log cell count and log median depth, and the contrasts |Module rs| - |%MT rs|, Module - depth and arm - HKG8; donors drawn once over all arms; BH within arm, level and row type (rs, contrast).
- Reading. The rows the Battery 2 claim rests on are the Module post-IP% rows of 16A and 16B on M4, PVM_3 and PVM_1, and the Module rs and arm - HKG8 contrast of 16F. If these hold, the donor structure and its intronic-rate association do not rest on the eleven. If they are lost beyond random removal, the Methods sentence "may depend on them" becomes "depends on them" for the named rows. HKG8 is read as the control: its rows are expected not to hold in 16A and 16B, as in the all-donor run.

### 19B. Histology and clinical rows by class: 17A within donor, 17B between donors, 17C within donor, 15A, 17D, 18A (sections 28 and 29)

- 17A within donor: rs of each of the five unit quantities with the six IBA1 measures; BH within arm and quantity over the six.
- 17B between donors: rs of the five quantities with the six CPS scores; BH within arm and quantity over the six.
- 17C within donor: the two rank coefficients (Module post-IP% and %MT) for each of the 28 stain measures; BH within arm and predictor over the 28; VIF recomputed.
- 15A: the cov4 Module coefficient, unadjusted and adjusted for donor-mean contamination, per arm, region set and outcome; BH as in 15A. dR2 and retention are carried without a test.
- 17D: the two coefficients of the cov4 model with both predictors; BH as in 17D; VIF recomputed.
- 18A: the two coefficients of each cov4 model with the second predictor and a nuclear-encoded panel; panel genes and donor scores fixed at the all-donor values; BH as in 18A; VIF recomputed.
- Reading. Class counts per panel and arm, drawn in Figure 19B. 17B between donors carries the exposure named in the Methods; 17A, 17C, 15A, 17D and 18A carry moderate exposure. A panel is read by the panel reading in the definitions.

### 19C. Verdicts with and without the eleven: 17C within donor, 17D, 18A (sections 28 and 29)

- The verdict of each model (Module only, %MT only, Both, Neither, Unreadable, Untestable; for 18A second only, NE panel only, Both, Neither, Unreadable, Untestable) is recomputed without the eleven with the Battery 3 and 4 rules, and cross-tabulated against the all-donor verdict.
- Reading. The Battery 3 and 4 conclusions rest on the absence of "Module only" rows where %MT is readable, and on the unreadable PVM_3 and PVM_1 rows. If the cross-tabulation moves rows into "Module only" or out of "Unreadable", the conclusion is re-stated for those rows; if the verdicts stay on the diagonal or move only between "Neither" and "%MT only", the Battery 3 and 4 conclusions stand.

### Influence. Leave-one-donor-out: 16A, 16B, 16F between donors (section 30)

- Every donor of the cohort is removed in turn and the 16A, 16B and 16F between-donor estimates recomputed (no test). Reported per donor with its SA and High-ADNC status: the change s x (e_minus_donor - e), and per row the rank of each SA donor among all donors by change. Tables only; no figure.
- Reading. Descriptive. It shows whether the SA-removal change in 19A comes from a few SA donors or from the group, and whether non-SA donors are as influential.
