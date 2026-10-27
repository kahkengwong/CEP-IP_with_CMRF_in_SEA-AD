# Container for external validation 2 (Rexach et al. 2024, insula)

This container reruns external validation 2, from the CELLxGENE h5ad to Supplementary Table S21. It is built on the SEA-AD container (`seaad-cepip`, README_container.md), so R 4.6.1, the 280 R packages of `renv.lock`, Python 3.13.14, and the packages of `requirements-python.txt` are those used for the paper. It adds three things:

1) cepip 0.1.0 at commit fe4b26934d2f4ff5aa1331a56f0e40ab70d8234f (R code identical to 4c4f70c, the CRAN submission).
2) brglm2 1.1.1, for the Firth logistic regression of R17, with the packages it needs that are not already installed (Posit Package Manager snapshot of 2026-10-25).
3) The scripts R11 to R17, `build_SuppTable_S21.py`, `compare_rexach.py`, and the run script `rexach-run`.

The build stops if any package of `renv.lock` changes version (`check_lock.R`). The container holds no data. The Rexach folder is mounted at `/data` when the container starts, and every script reads and writes inside that folder.

## 1. Build

In the repository, the build folder is `external_validation/Rexach/`. The SEA-AD image must be present first, built from the repository or loaded from Zenodo (`podman load -i seaad-cepip_image.tar.gz`). `podman images` should list `localhost/seaad-cepip`.

The build folder holds:

```
Containerfile
check_lock.R
rexach-run.sh
compare_rexach.py
build_SuppTable_S21.py
R11_Rexach_cellxgene_check.py
R12_Rexach_mapping_input.py
R13_Rexach_count_checks.py
R14_Rexach_core_cepip.R
R15_Rexach_diagnosis_contrasts.R
R16_Rexach_main_figure.R
R17_Rexach_posthoc_beta.R
```

From that folder:

```bash
podman build -t rexach-cepip -f Containerfile .
```

The build takes a few minutes and needs an internet connection (GitHub for cepip, Posit Package Manager for brglm2). With Docker, if the SEA-AD image is listed under another name (for example `seaad-cepip:latest`):

```bash
docker build --build-arg BASE=seaad-cepip:latest -t rexach-cepip -f Containerfile .
```

The build log ends with `check_lock.R: PASS`, and lists cepip with its commit and the packages added beyond `renv.lock`.

The built image is also on Zenodo, as `rexach-cepip_image.tar.gz` in the same record as the SEA-AD image (https://doi.org/10.5281/zenodo.23653156). `podman load -i rexach-cepip_image.tar.gz` (or `docker load -i`) gives `localhost/rexach-cepip`, so the build can be skipped. The image holds the scripts of this folder at commit 8b98869.

## 2. Run

Make a new folder that holds only the inputs of section 3, for example `Rexach_rerun`. Then:

```bash
podman run --rm -v /path/to/Rexach_rerun:/data:Z rexach-cepip rexach-run all
```

`:Z` lets the container use the folder on systems with SELinux, such as Fedora. With Docker Desktop on Windows, give the Windows path and leave out `:Z`:

```bash
docker run --rm -v "C:\Users\<name>\Rexach_rerun:/data" rexach-cepip rexach-run all
```

To compare the rerun with the original results, mount the original `External_Dataset_2` folder read-only at `/ref`, and the original `Manuscript/Supp_Tables` folder at `/ref_supp`, and add `--compare`:

```bash
podman run --rm -v /path/to/Rexach_rerun:/data:Z \
    -v /path/to/SEA-AD/External_Dataset_2:/ref:ro,Z \
    -v /path/to/SEA-AD/Manuscript/Supp_Tables:/ref_supp:ro,Z \
    rexach-cepip rexach-run all --compare
```

| Command | Runs | Needs |
|---|---|---|
| `rexach-run all` | R11 to R13 (Python), R14 to R17 (R), and Table S21 | section 3, items 1 to 4 |
| `rexach-run r` | R14 to R17 and Table S21 | section 3, items 3 to 5 |
| `rexach-run list r` | prints the steps of a mode | - |

Settings, given with `-e` on the run line:

1) `REXACH_FIGURE`. The main figure number in the Table S21 title and file name (default 19).
2) `REXACH_WORKERS`. The workers of `cepip()` in R14 (default 6). The results do not depend on the number of workers. On a laptop with less memory, use `-e REXACH_WORKERS=2`.

Each step's console output is saved to `rexach-run_logs/<step>.txt`, and the start, end, and result of every step to `rexach-run_logs/rexach-run_summary.txt`. The run stops at the first step that fails, because later steps read earlier outputs, and prints the command to resume from that step, for example:

```bash
podman run --rm -v /path/to/Rexach_rerun:/data:Z rexach-cepip rexach-run all --from R15_Rexach_diagnosis_contrasts
```

A step fails when it exits with an error, or when its log shows a Python traceback, a failed check (`[FAIL]`, `VERDICT: FAIL`) or a stop (`STOPPED`, `ERROR:`, `Execution halted`).

## 3. Inputs

Place these in the data folder before running:

1) The CELLxGENE h5ad `aa87f914-07e0-48aa-b915-8de906c95baf.h5ad` (dataset version aa87f914-07e0-48aa-b915-8de906c95baf, collection c53573b2-eff4-4c5e-9ad0-b24d422dfd9b; 4,993,733,632 bytes), as the only `.h5ad` in the folder. Read by R11, R12, and R13. Download: https://datasets.cellxgene.cziscience.com/aa87f914-07e0-48aa-b915-8de906c95baf.h5ad
2) `R13_mapmycells_output/`, with the six MapMyCells csv files of the original run (three parts, each with Deep Generative Mapping and Hierarchical Mapping). MapMyCells is a web service of the Allen Institute and is not run in the container. R13 checks that each file covers exactly the nuclei of its part in the rerun of R12, so the rerun stops if R12 does not reproduce the parts. Download: the Hugging Face dataset https://huggingface.co/datasets/kahkengwong/CEP-IP_with_CMRF_in_SEA-AD, folder `external_validation/Rexach/R13_mapmycells_output/`.
3) `mmc1.xlsx`, Table S1 of Rexach et al. 2024 (Cell 187:5753), read by R15. It is the publisher's file and is not redistributed here. Download it from the Supplemental Information of the paper.
4) `Part12c_cell_metrics.csv`, the %MT per cell of the SEA-AD arms, read by R15. Download: the same Hugging Face dataset, folder `external_validation/Rexach/`.
5) For `rexach-run r` only: `R12_mapmycells_input/R12_Rexach_nuclei.csv.gz` and the three files of `R13_count_checks/` (`R13_Rexach_MicroPVM3_counts.mtx.gz`, `_genes.csv.gz`, and `_nuclei.csv.gz`).

`rexach-run` lists the inputs it finds before the first step and stops if one is missing.

In WSL, `22b_download_rexach_inputs.sh` downloads items 1, 2, and 4 into a new folder `~/Rexach_rerun`, with progress shown, resumes a stopped download when it is run again, and checks every file against the size of the original. Give the path of `mmc1.xlsx` with `MMC1=`:

```bash
MMC1=/mnt/c/Users/<name>/Downloads/mmc1.xlsx bash 22b_download_rexach_inputs.sh
```

## 3b. Step scripts for WSL

The three scripts below run the whole check in WSL with Podman, each with a log (`<script>_log.txt`) in the build folder:

1) `21_build_rexach_image.sh`. Builds `rexach-cepip` on the SEA-AD image. If the SEA-AD image is not loaded yet, give the Zenodo archive: `SEAAD_TAR=/path/to/seaad-cepip_image.tar.gz bash 21_build_rexach_image.sh`.
2) `22b_download_rexach_inputs.sh` (or `22_prepare_rexach_rerun.sh`, which copies the inputs from an existing External_Dataset_2 folder). Makes `~/Rexach_rerun` with the inputs of section 3.
3) `23_run_rexach_container.sh`. Runs `rexach-run all --compare`. It needs the reference results at `ORIG=` (an External_Dataset_2 folder with the original outputs) and `SUPP=` (the folder holding the original Table S21). The run summary, the comparison, and the guardrail logs are copied to `OUTC=` (default `/mnt/c/Users/Wong/Downloads/rexach_results`).

On Ubuntu 22.04 in WSL, Podman is installed with `sudo apt install --no-install-recommends podman uidmap slirp4netns fuse-overlayfs catatonit crun`.

## 4. Run order and outputs

1) R11_Rexach_cellxgene_check.py. Donors, diagnoses, and microglia of the insula. Outputs in the data folder.
2) R12_Rexach_mapping_input.py. The 500-gene floor and the MapMyCells input files. `R12_mapmycells_input/`
3) R13_Rexach_count_checks.py. The count checks and the counts of the Micro-PVM_3 nuclei. `R13_count_checks/`
4) R14_Rexach_core_cepip.R. The DFG screen, the IP gate, and the cepip fits of the three arms. `R14_readings/`, `R14_PVM_3/`, `R14_HKG8/`, and `R14_HKG_control/`
5) R15_Rexach_diagnosis_contrasts.R. The diagnosis contrasts (cov3). `R15_readings/`
6) R17_Rexach_posthoc_beta.R. The post hoc analyses of the betas. `R17_posthoc/`
7) R16_Rexach_main_figure.R. The panels of the main figure and their source data. `R16_Figure_Rexach/`
8) build_SuppTable_S21.py. Supplementary Table S21. `Supp_Tables/`

Each script writes a log and a guardrail log (PASS, WARN, and FAIL) next to its outputs.

## 5. Comparison

`compare_rexach.py` uses the comparison of `verification/compare_rerun.py` of the SEA-AD analysis, with the same tolerances. Every CSV and Excel table (every sheet) is compared cell by cell. Numbers agree when their relative difference is within 1e-6, and text must match exactly. Cells that hold a path, a date, a time, or a BLAS library are counted separately. The R12 MapMyCells input files and the R13 count matrix are compared exactly (nuclei, genes, and every count). Each table is given one result, from SAME to OUTSIDE 0.1%, as in README_container.md. A Session sheet whose only differences are in its text (versions, paths, times, and the BLAS library) is reported as PROVENANCE ONLY.

The comparison is written to `compare/compare_rexach_<mode>_log.txt` (summary) and `compare/compare_rexach_<mode>_detail.csv` (one row per file and sheet). With `-e COMPARE_CELLS=1`, each differing cell is also listed.

`external_validation/Rexach/verification/` holds the comparison of the author's container rerun, made from a new folder with the inputs of section 3, with the original results: `rexach-run_summary.txt` (the steps and their times), `compare_rexach_all_log.txt`, and `compare_rexach_all_detail.csv`.

`external_validation/Rexach/verification/independent_test/` holds the same check for an independent rerun on another computer. It ran on Windows 11 with Docker Desktop, with the image built from this folder at commit 8b98869. The folder holds the run summary and two comparisons, each with its detail file. The first is with the original results (`compare_tester_vs_original_log.txt`), and the second with the author's container rerun (`compare_tester_vs_container_log.txt`).

## 6. Notes

1) Systems. The original R steps ran on Windows 11 (R 4.6.1) and the Python steps in WSL. The container runs Linux, so some numbers can differ in the last digits across operating systems and linear algebra libraries.
2) Random numbers. Every script that draws random numbers uses fixed seeds (cepip CMRF repeats in R14, the permutations and the bootstrap in R17).
3) Fonts. R16 draws the panels with the default sans font of the system, so the letter shapes in the container can differ slightly from the Windows panels. The figure values are in the source data and in Table S21.
