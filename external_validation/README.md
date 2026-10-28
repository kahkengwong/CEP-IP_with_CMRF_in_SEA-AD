# External Validation

The pipeline was run on two independent snRNA-seq datasets with the R package [cepip](https://github.com/kahkengwong/cepip) (0.1.0, commit fe4b269) at the SEA-AD settings. Each dataset was run on the PVM_3 arm (Micro-PVM_3 nuclei with *MT-CO3*>0), with the HKG8 arm as the control.

## 1. Contents of this folder
| Folder | Contents |
|---|---|
| [`Pineda/`](Pineda/) | External dataset 1: the download and Cell Ranger scripts (01 to 08), the scripts R02 to R17, the Table S20 builder, the container recipe, the run script `pineda-run`, the WSL step scripts, and [`README_container_Pineda.md`](Pineda/README_container_Pineda.md) |
| `Pineda/verification/` | The comparison of the author's container rerun with the original results, and `independent_test/` for the rerun on another computer |
| [`Rexach/`](Rexach/) | External dataset 2: the scripts R11 to R17, the Table S21 builder, the container recipe, the run script `rexach-run`, the WSL step scripts, and [`README_container_Rexach.md`](Rexach/README_container_Rexach.md) |
| `Rexach/verification/` | The comparison of the author's container rerun with the original results, and `independent_test/` for the rerun on another computer |

## 2. Datasets and results
| Item | External dataset 1: Pineda et al. (2024) | External dataset 2: Rexach et al. (2024) |
|---|---|---|
| **Region and diagnoses** | Prefrontal cortex (BA9); ALS, FTLD, and pathologically normal controls (PN) | Insular cortex; AD, bvFTD, PSP, and controls |
| **Donors** | 63 (49 with at least 5 Micro-PVM_3 nuclei: 13 ALS, 23 FTLD, and 13 PN) | 37 (33 with at least 5 Micro-PVM_3 nuclei: 10 AD, 8 bvFTD, 8 PSP, and 7 controls) |
| **Input** | Raw reads (SRA PRJNA1073234, SRA Lite) counted with Cell Ranger 10.1.0 (GRCh38-2020-A, intronic reads included); cells called by emptyDrops (FDR ≤0.01); nuclei with <500 genes and scDblFinder doublets removed; only nuclei kept by the authors' quality control (≥700 genes, ≤20% MT, and removal by cluster) | Raw counts from CELLxGENE, with the authors' quality control (nuclei with <200 genes or >8% MT removed); nuclei with <500 genes removed |
| **Cell types** | MapMyCells Deep Generative Mapping to the SEA-AD taxonomy (CCN20230505), after the authors' quality control | MapMyCells Deep Generative Mapping to the SEA-AD taxonomy (CCN20230505) |
| **DFG replication** | Replicated at the SEA-AD thresholds: 5 of 8 DFGs at tier 1 (r<sub>s</sub> ≥0.6 and τ<sub>b</sub> ≥0.5), and 2 more at tier 2 | Replicated at the SEA-AD thresholds: 5 of 8 DFGs at tier 1, and the other 3 at tier 2 |
| **IP gate** | PVM_3: IPRS 0.754 (strong tier), passed. HKG8: no IP (failed). | PVM_3: IPRS 0.758 (strong tier), passed. HKG8: no IP (IPRS 0.200, failed). |
| **Module** | 1,225 of 1,817 PVM_3 arm cells, of which 219 Recovered by CMRF | 2,578 of 3,417 PVM_3 arm cells, of which 430 Recovered by CMRF |
| **Module post-IP% against diagnosis** | Standardized β against PN, adjusted for age at death and sex: FTLD -0.122 and ALS -0.187. With PMI added: FTLD -0.019 and ALS -0.074. FDR=0.923 for all four contrasts. | Standardized β against controls, adjusted for age at death, sex, and PMI: all three diagnoses 0.277 (FDR=0.261), AD 0.478 (FDR=0.254), bvFTD 0.641 (FDR=0.254), and PSP 0.064 (FDR=0.882) |
| **Post hoc** | A permutation test, leave-one-donor-out refits, a donor bootstrap, Firth logistic regression, the difference from the SEA-AD β, and the power to detect the SEA-AD NeoC β (0.465) | The same as for external dataset 1 |
| **Preprint** | Figure 19 and Supplementary Table S20 | Figure 20 and Supplementary Table S21 |

## 3. Rerunning an external validation
Each external validation has its own container, built on the SEA-AD container. The steps, inputs, and outputs are in the README of its folder.

| Item | External dataset 1: Pineda et al. (2024) | External dataset 2: Rexach et al. (2024) |
|---|---|---|
| **Instructions** | [`README_container_Pineda.md`](Pineda/README_container_Pineda.md) | [`README_container_Rexach.md`](Rexach/README_container_Rexach.md) |
| **Container image** | `pineda-cepip_image.tar.gz` on Zenodo (https://doi.org/10.5281/zenodo.23653156), or built from `Pineda/Containerfile` | `rexach-cepip_image.tar.gz` on Zenodo (https://doi.org/10.5281/zenodo.23653156), or built from `Rexach/Containerfile` |
| **Inputs** | The 63 Cell Ranger raw count matrices, the MapMyCells outputs, and the other input files from [Hugging Face](https://huggingface.co/datasets/kahkengwong/CEP-IP_with_CMRF_in_SEA-AD), folder `external_validation/Pineda/`. The S01 workbook and the authors' per-library metadata from Synapse (syn51105515). | The h5ad from CELLxGENE. The MapMyCells outputs and the SEA-AD %MT table from [Hugging Face](https://huggingface.co/datasets/kahkengwong/CEP-IP_with_CMRF_in_SEA-AD), folder `external_validation/Rexach/`. Table S1 from the paper. |
| **Command** | `podman run --rm -v /path/to/Pineda_rerun:/data:Z pineda-cepip pineda-run all` | `podman run --rm -v /path/to/Rexach_rerun:/data:Z rexach-cepip rexach-run all` |
| **Run time** | R02 (cell calling): about 1.5 to 2 hours (108 minutes with 6 workers on Windows). R11 to R17 and Table S20: about 15 minutes in the container. `pineda-run r` (R14 to R17 and Table S20): about 5 minutes on the Windows 11 laptop of the independent test. | R11 to R17 and Table S21: about 8 minutes on the development machine, and about 10 minutes on the Windows 11 laptop of the independent test |

## 4. Reproducibility check
The external validations were checked as the SEA-AD pipeline ([Section 7 of the main README](../README.md#7-reproducibility-check)), with the same rules: numbers agree when their relative difference is within 1e-6, and text must match exactly. The tables that differ only in their run records (versions, paths, times, the BLAS library, and the number of workers) are counted as agreed.

| Check | External dataset 1: Pineda et al. (2024) | External dataset 2: Rexach et al. (2024) |
|---|---|---|
| **Container rerun by the author** (Podman in WSL, new folder with the inputs only) | From the Cell Ranger count matrices: the cell calls were identical in all 63 libraries, and 151 of 152 tables agreed with the original results. The other table is the R11 nuclei table, whose scDblFinder scores differ by at most 3.7e-4; no doublet class changed. | 113 of 113 tables agreed with the original results |
| **Independent test** (second person, Windows 11, Docker Desktop) | `pineda-run r` (R14 to R17 and Table S20, from the outputs of R13 to R13c): 93 of 93 tables agreed with the original results, and 93 of 93 with the author's container rerun (28 October 2026, image built at commit 4268f76) | 113 of 113 tables agreed with the original results, and 113 of 113 with the author's container rerun (26 October 2026, image built at commit 8b98869) |
| **Files** | `Pineda/verification/` and `Pineda/verification/independent_test/` | `Rexach/verification/` and `Rexach/verification/independent_test/` |

## 5. References
- S.S. Pineda, H. Lee, M.J. Ulloa-Navas, R.M. Linville, F.J. Garcia, et al., Single-cell dissection of the human motor and prefrontal cortices in ALS and FTLD, Cell, 187 (2024) 1971-1989. https://doi.org/10.1016/j.cell.2024.02.031
- J.E. Rexach, Y. Cheng, L. Chen, D. Polioudakis, L.-C. Lin, et al., Cross-disorder and disease-specific pathways in dementia revealed by single-cell genomics, Cell, 187 (2024) 5753-5774. https://doi.org/10.1016/j.cell.2024.08.019
