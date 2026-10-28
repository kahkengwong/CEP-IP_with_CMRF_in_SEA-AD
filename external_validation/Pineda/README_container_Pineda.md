# Container for external validation 1 (Pineda et al. 2024, BA9)

This container reruns external validation 1, from the Cell Ranger raw count matrices to Supplementary Table S20. It is built on the SEA-AD container (`seaad-cepip`, README_container.md), so R 4.6.1, the 280 R packages of `renv.lock`, Python 3.13.14, and the packages of `requirements-python.txt` are those used for the paper. It adds four things:

1) cepip 0.1.0 at commit fe4b26934d2f4ff5aa1331a56f0e40ab70d8234f (R code identical to 4c4f70c, the CRAN submission).
2) brglm2 1.1.1, for the Firth logistic regression of R17 (Posit Package Manager snapshot of 2026-10-25), as in the Rexach container.
3) scDblFinder 1.26.9, for the doublets of R02, and the 16 other packages of `pineda_r_packages.tsv`, at the versions of the R02 run on Windows (Bioconductor 3.23 and the same snapshot).
4) The scripts R02 to R17, `build_SuppTable_S20.py`, `compare_pineda.py`, and the run script `pineda-run`.

The build stops if any package of `renv.lock` or `pineda_r_packages.tsv` has another version (`check_lock_pineda.R`). The container holds no data. The Pineda folder is mounted at `/data` when the container starts, and every script reads and writes inside that folder.

## 1. Steps before the container

The raw reads (SRA PRJNA1073234, SRA Lite files) were downloaded and counted with Cell Ranger 10.1.0 on a Fedora workstation, with the scripts `01_setup.sh` to `08_batch_ncbi.sh` in this folder. Cell Ranger is not rerun in the container, because the reads total about 1.5 TB. The container starts from the 63 raw count matrices that Cell Ranger wrote.

| Script | Step |
|---|---|
| `01_setup.sh` | Checks Cell Ranger 10.1.0 and downloads the GRCh38-2020-A reference |
| `03_pilot_download.sh`, `02_pilot_count.sh` | The pilot library (210526_FTLD_228_snRNA-C8), from the ENA FASTQ |
| `07_ncbi_fallback.sh` | One library (210526_PN_301_snRNA-D7) from NCBI SRA Lite |
| `08_batch_ncbi.sh` | The other 61 libraries from NCBI SRA Lite, streamed to paired FASTQ |
| `05_batch_count.sh` | `cellranger count` of each library (`--chemistry=auto --include-introns=true --create-bam=false`), and `05_cellranger_metrics.csv` |

Tools: Cell Ranger 10.1.0, sra-tools 3.4.1 (fastq-dump and fasterq-dump), aria2 1.37.0, pigz 2.8, and mawk 1.3.4. `Pineda_BA9_210526_manifest.tsv` lists the 63 libraries and their SRA runs.

## 2. Build

In the repository, the build folder is `external_validation/Pineda/`. The SEA-AD image must be present first, built from the repository or loaded from Zenodo (`podman load -i seaad-cepip_image.tar.gz`). `podman images` should list `localhost/seaad-cepip`.

From the build folder:

```bash
podman build -t pineda-cepip -f Containerfile .
```

The build takes about 20 to 40 minutes, because the Bioconductor packages are compiled, and needs an internet connection (GitHub for cepip, Posit Package Manager for the CRAN packages, Bioconductor for scDblFinder and its packages). With Docker, if the SEA-AD image is listed under another name (for example `seaad-cepip:latest`):

```bash
docker build --build-arg BASE=seaad-cepip:latest -t pineda-cepip -f Containerfile .
```

The build log ends with `check_lock_pineda.R: PASS`, and lists cepip with its commit and the packages added beyond `renv.lock` and the R02 list.

## 3. Inputs

Place these in the data folder before running:

1) `counts/`, with the 63 Cell Ranger raw count matrices `<LIB>_raw_feature_bc_matrix.h5`. Read by R02, R12, R13, and R13c.
2) `05_cellranger_metrics.csv` (the Cell Ranger metrics of the 63 libraries, read by R02) and `Pineda_BA9_210526_manifest.tsv`.
3) The S01 workbook `S01_sample_information.xlsx` of Pineda et al. 2024 (Synapse syn51105515), read by R11 and R15.
4) `Pineda2024_metadata/`, with the authors' `<LIB>_col_metadata.tsv` file of each library (Synapse syn51105515). R12 and R13 keep the nuclei of the authors' published QC: a nucleus is kept when its barcode is in the authors' file.
5) `Stats_15H.csv`, the SEA-AD source data of Figure 15H, read by R13 for the expected r_s.
6) The MapMyCells results of the original run. MapMyCells is a web service of the Allen Institute and is not run in the container. The same file mapped twice gave another Deep Generative Mapping supertype for 0.8% of nuclei, so the rerun uses these files:
   - `R13_mapmycells_output/`: 12 csv files (six parts, each with Deep Generative Mapping and Hierarchical Mapping), the nuclei of the authors' QC. R13 checks that each file covers exactly the nuclei of its part in the rerun of R12, so the rerun stops if R12 does not reproduce the parts.
   - `R13_mapmycells_output_allKept/` and `R13_mapmycells_output_allKept_repeat/`: 26 csv files each, all R02 kept nuclei mapped before the authors' QC, and the same files mapped again. Read by R13b and R13c for sheets J and K of Table S20.
7) For `pineda-run r` only: the outputs of R13, R13b, and R13c (`pineda-run` lists them).

`pineda-run` lists the inputs it finds before the first step and stops if one is missing.

## 4. Run

```bash
podman run --rm -v /path/to/Pineda_rerun:/data:Z pineda-cepip pineda-run all
```

`:Z` lets the container use the folder on systems with SELinux, such as Fedora. With Docker Desktop on Windows, give the Windows path and leave out `:Z`. To compare the rerun with the original results, mount the original `External_Dataset_1` folder read-only at `/ref`, and the original `Manuscript/Supp_Tables` folder at `/ref_supp`, and add `--compare`:

```bash
podman run --rm -v /path/to/Pineda_rerun:/data:Z \
    -v /path/to/SEA-AD/External_Dataset_1:/ref:ro,Z \
    -v /path/to/SEA-AD/Manuscript/Supp_Tables:/ref_supp:ro,Z \
    pineda-cepip pineda-run all --compare
```

| Command | Runs | Needs |
|---|---|---|
| `pineda-run all` | R02 and R11 (R), R12 to R13c (Python), R14 to R17 (R), and Table S20 | section 3, items 1 to 6 |
| `pineda-run r` | R14 to R17 and Table S20 | section 3, items 2, 3, and 7 |
| `pineda-run list all` | prints the steps of a mode | - |

Settings, given with `-e` on the run line:

1) `PINEDA_FIGURE`. The main figure number in the Table S20 title and file name (default 19).
2) `PINEDA_WORKERS`. The libraries processed at the same time in R02 and the workers of `cepip()` in R14 (default 6). The results do not depend on the number of workers. Each R02 worker holds one raw matrix in memory and used up to 9 GB in the author's rerun, so 4 workers need about 36 GB. On a computer with less memory, use `-e PINEDA_WORKERS=2`. WSL gives Linux half of the Windows memory by default; a `.wslconfig` file in the Windows user folder with the lines `[wsl2]`, `memory=40GB`, and `swap=16GB` raises it after `wsl --shutdown`.

Each step's console output is saved to `pineda-run_logs/<step>.txt`, and the start, end, and result of every step to `pineda-run_logs/pineda-run_summary.txt`. The run stops at the first step that fails, because later steps read earlier outputs, and prints the command to resume from that step. A step fails when it exits with an error, or when its log shows a Python traceback, a failed check (`[FAIL]`, `VERDICT: FAIL`) or a stop (`STOPPED`, `ERROR:`, `Execution halted`). R02 also fails unless all 63 libraries finished.

### Step scripts for WSL

The three scripts below run the whole check in WSL with Podman, each with a log (`<script>_log.txt`) in the build folder:

1) `21_build_pineda_image.sh`. Builds `pineda-cepip` on the SEA-AD image. If the SEA-AD image is not loaded yet, give the Zenodo archive: `SEAAD_TAR=/path/to/seaad-cepip_image.tar.gz bash 21_build_pineda_image.sh`.
2) `22_prepare_pineda_rerun.sh`. Makes `~/Pineda_rerun` with the inputs of section 3, copied from an existing External_Dataset_1 folder (`ORIG=`), and checks every count and size.
3) `23_run_pineda_container.sh`. Runs `pineda-run all --compare`, with the reference results at `ORIG=` and the folder holding the original Table S20 at `SUPP=`. The run summary, the comparison, and the guardrail logs are copied to `OUTC=` (default `/mnt/c/Users/Wong/Downloads/pineda_results`). It uses 4 workers by default, because WSL has half the memory of Windows unless set otherwise.

## 5. Run order and outputs

1) R02_Pineda_cell_calling_QC.R. emptyDrops, the 500-gene floor, and scDblFinder, per library. `R02_cell_calling/`
2) R11_Pineda_nuclei_table.R. The nuclei and donor tables, with the S01 covariates. `R11_nuclei/`
3) R12_Pineda_mapping_input.py. The authors' QC and the MapMyCells input files. `R12_mapmycells_input/`
4) R13_Pineda_count_checks.py. The count checks and the counts of the Micro-PVM_3 nuclei. `R13_count_checks/`
5) R13b_Pineda_mapping_comparison.py, run twice: the repeat mapping, and the QC before mapping. `R13b_mapping_comparison_repeat/` and `R13b_mapping_comparison_QC/`
6) R13c_Pineda_noAuthorQC_arm.py. The PVM_3 arm of the run without the authors' QC. `R13c_noAuthorQC_arm/`
7) R14_Pineda_core_cepip.R. The DFG screen, the IP gate, and the cepip fits of the three arms. `R14_core/`, `R14_PVM_3/`, `R14_HKG8/`, and `R14_HKG_control/`
8) R15_Pineda_diagnosis_contrasts.R. The four diagnosis contrasts (cov2 and cov3). `R15_contrasts/`
9) R17_Pineda_posthoc_beta.R. The post hoc analyses of the betas. `R17_posthoc/`
10) R16_Pineda_main_figure.R. The panels of the main figure and their source data. `R16_Figure_Pineda/`
11) build_SuppTable_S20.py. Supplementary Table S20. `Supp_Tables/`

Each script writes a log next to its outputs, and R14 to R17 also write a guardrail log (PASS, WARN, and FAIL).

## 6. Comparison

`compare_pineda.py` uses the comparison of `compare_rexach.py` and of `verification/compare_rerun.py` of the SEA-AD analysis, with the same tolerances. Every CSV and Excel table (every sheet) is compared cell by cell. Numbers agree when their relative difference is within 1e-6, and text must match exactly. Cells that hold a path, a date, a time, or a BLAS library are counted separately, and the column `minutes` of the R02 summary is not compared. The R12 MapMyCells input files and the R13 count matrix are compared exactly (nuclei, genes, and every count). A Session sheet records the run (versions, paths, times, the BLAS library, and the number of workers), so a Session sheet that differs is reported as PROVENANCE ONLY.

The R02 cell calls are also compared nucleus by nucleus for each library, from the R11 nuclei table: the called nuclei, the emptyDrops FDR, the scDblFinder score and class, and the kept nuclei. The result is written to `compare/compare_pineda_all_R02_libraries.csv`.

The comparison is written to `compare/compare_pineda_<mode>_log.txt` (summary) and `compare/compare_pineda_<mode>_detail.csv` (one row per file and sheet). With `-e COMPARE_CELLS=1`, each differing cell is also listed.

`external_validation/Pineda/verification/` holds the comparison of the author's container rerun with the original results. The rerun ran in WSL on Windows 11 with Podman, from a new folder with the inputs of section 3, with 4 workers. The folder holds the run summary (`pineda-run_summary.txt`, which also records the earlier attempts that stopped when R02 ran out of memory; R02 resumed each time from the libraries already finished), `compare_pineda_all_log.txt`, `compare_pineda_all_detail.csv`, and `compare_pineda_all_R02_libraries.csv`. The cell calls were identical in all 63 libraries. All result tables agreed within a relative difference of 1e-6, except the scDblFinder scores of the R11 and R12 nuclei tables (largest difference 3.7e-4, no doublet class changed).

`external_validation/Pineda/verification/independent_test/` holds the same check for an independent rerun on another computer. It ran `pineda-run r` (R14 to R17 and Table S20, from the outputs of R13 to R13c) on Windows 11 with Docker Desktop, with the image built from this folder at commit 4268f76. The folder holds the run summary and two comparisons, each with its detail file. The first is with the original results (`compare_tester_vs_original_log.txt`, 93 of 93 tables agreed), and the second with the author's container rerun (`compare_tester_vs_container_log.txt`, 93 of 93 tables agreed).

## 7. Notes

1) Systems. The original R steps ran on Windows 11 (R 4.6.1, reference BLAS) and the Python steps in WSL. The container runs Linux, so some numbers can differ in the last digits across operating systems and linear algebra libraries. scDblFinder runs a PCA in each library, so a doublet call near its threshold can differ between the two systems; the R02 comparison shows any such nucleus.
2) Random numbers. Every script that draws random numbers uses fixed seeds (emptyDrops and scDblFinder in R02, the cepip CMRF repeats in R14, the permutations and the bootstrap in R17).
3) Python. The original Python steps (R12 to R13c) ran in WSL with Python 3.14.6. The container has Python 3.13.14 and the packages of the SEA-AD container. The R12 MapMyCells input files and the R13 count matrix of the rerun were identical to the original files.
4) Kendall tau-b. R14 computes tau-b with `pcaPP::cor.fk` when pcaPP is installed, as on Windows, and otherwise with `stats::cor`, which gives the same values more slowly. The container does not install pcaPP; the DFG sheets of R14 agreed within 1e-6.
5) Fonts. R16 draws the panels with the default sans font of the system, so the letter shapes in the container can differ slightly from the Windows panels. The figure values are in the source data and in Table S20.
