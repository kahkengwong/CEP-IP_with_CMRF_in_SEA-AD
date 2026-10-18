# Container for the SEA-AD CEP-IP/CMRF analysis

The container holds the software used for the paper: R 4.6.1 with the packages in `renv.lock`, Python 3.13.14 with the packages in `requirements-python.txt`, and the deposited scripts. It does not hold any data. The SEA-AD folder is mounted at `/data` when the container starts, and every script reads and writes inside that folder.

The scripts find the project folder through the environment variable `SEAAD_ROOT`, which the container sets to `/data`. Outside the container, set `SEAAD_ROOT` to the SEA-AD folder, or start R or Python from that folder.

## 1. Build

Podman is used here (it is installed on Fedora). Docker accepts the same commands with `docker` in place of `podman`.

The build folder holds:

```
Containerfile
seaad-setup.sh
renv.lock
requirements-python.txt
scripts/R/       the deposited R scripts, the Part 12 JSON settings and the Part 12
                 pre-specified readings (Part12_Preset_Readings*.md)
scripts/Python/  the deposited Python scripts
```

From that folder:

```bash
podman build --no-cache -t seaad-cepip -f Containerfile .
```

The first build takes about one to two hours, mostly to compile the Bioconductor packages and monocle3. It needs an internet connection. The image is about 6 to 8 GB.

Instead of building it, the image can be loaded from the archived copy on Zenodo (`seaad-cepip_image.tar.gz`, 1.7 GB, https://doi.org/10.5281/zenodo.23446228), which holds the same software and scripts:
```bash
podman load -i seaad-cepip_image.tar.gz
```
This gives the image `localhost/seaad-cepip`; Docker accepts `docker load -i` in the same way.

For an HPC cluster that uses Apptainer instead of Podman:

```bash
podman save -o seaad-cepip.tar seaad-cepip
apptainer build seaad-cepip.sif docker-archive://seaad-cepip.tar
```

## 2. Start the container

```bash
podman run --rm -it --shm-size=16g -v /path/to/SEA-AD:/data:Z seaad-cepip
```

`:Z` lets the container use the folder on systems with SELinux, such as Fedora. `--shm-size=16g` is needed for Part 7 (DecontX), whose worker processes share the count matrix through shared memory; with the default size, Part 7 stops with a bus error. Inside the container, run once:

```bash
seaad-setup
```

This copies the shared settings files (`CEP_arms.R`, `R_config.R`, `R_config_HKG8.R`) into the data folder, and the Part 12 scripts, settings and pre-specified readings into `Part12_rerun`, if they are missing. The Part 12 scripts are run from `Part12_rerun`, because each one checks that its own copy and `Part12_03_Config.r` are there. It reports which inputs it finds and does not change any file that is already there.

## 3. Inputs

Place these in the data folder before running:

1) The SEA-AD immune h5ad, `SEAAD_Immune_HIP_MEC_LEC_ITG_MTG_FI_STG_DFC_AnG_V1C_RNAseq_final-nuclei.2026-06-22.h5ad`, read by Parts 1 to 8 and Part12_02.
2) The SEA-AD quantitative neuropathology workbook (`sea-ad-quantitative-neuropathology-063026.xlsx`) and the CPS workbook (`sea-ad_global_and_local_cps_050126.xlsx`), used by the Batteries.
3) The Gerrits et al. deposit (GSE148822) in `Gerrits/`, with the series matrix in `Gerrits/Metadata/`, used by Part12_06 and Battery 1.

To run from Part 9 without Parts 1 to 8, also place the Part 8 export (`R_export_CP10K/`), the Part 2 donor tables (`*_S5b_donor_clinical.csv`) and the Part 7 contamination table (`Part7_DecontX/decontX_contamination_all.csv`) in the data folder.

Two scripts start from the supplementary tables (Zenodo) instead of computing them. Place these copies before running them:
- `Part10g_CP10K_DEGs_Part1_code_S11_followup_analyses.r` reads `Manuscript/Supp_Tables/Supplementary_Table_S11_DEGs_ModvsNonMod.xlsx`.
- `Fig17_Summary_code.r` reads `Manuscript/Supp_Tables/Supplementary_Table_S12_Prognosis_&_Ablation.xlsx` and `Part12_rerun/07_Fig15/SuppTable_Fig15.xlsx`, `08_Fig16/SuppTable_Fig16.xlsx`, `09_Fig17/SuppTable_Fig17.xlsx`, `10_Fig18/SuppTable_Fig18.xlsx` and `11_Fig19/SuppTable_Fig19.xlsx` (all under `Part12_rerun`).

Each script lists its own inputs in its header.

## 4. Run order

Run every script from the data folder (`/data`). R scripts are run with `source()`, as in the original runs. Each script writes a log next to its outputs. Part 8 writes one arm per run; the arm is selected at the top of the script.

```bash
# Python (Parts 1 to 8, then the Part 12 extraction)
python /opt/seaad/Python/Part1a_NormForensics_code.py      # Parts 1a to 1g, in order
python /opt/seaad/Python/Part2_NeoC_CP10K_code.py
python /opt/seaad/Python/Part2_NonNeoC_CP10K_code.py
#   Parts 3 to 6b (Figure 1, Supplementary Figure S1, Tables S1 to S5)
python /opt/seaad/Python/Part7_DecontX.py
python /opt/seaad/Python/Part8_Data-export-CP10K_MT-CO3_8DFG.py

# R
Rscript -e 'source("/opt/seaad/R/Part9_CP10K_QC_GAM_code.r")'

# Part 12 (after seaad-setup has copied the files into Part12_rerun)
python Part12_rerun/Part12_02_Extract.py
Rscript -e 'source("Part12_rerun/Part12_04_Labels.r")'
```

The full order of the R scripts:

1) Part9_CP10K_QC_GAM_code.r
2) Part9b_CP10K_SiKi_code.r
3) Part10_CP10K_SKCS_code.r
4) Part10b_CP10K_IPRS_code.r
5) Part10b_CP10K_SKCS_TREP_code.r
6) Part10c_CP10K_Copula1_code.r
7) Part10d_CP10K_Copula2_code.r
8) Part10e_CP10K_CMRF_code.r
9) Part10e_CP10K_Validation1_code.r
10) Part10e_CP10K_Validation2_UMAP_code.r
11) Part10f_CP10K_Band_Formation_code.r
12) Part10f_CP10K_CMRF_RSkfCV_Mechanism_code.r
13) Part10g_CP10K_DEGs_Part1_code.r
14) Part10g_Step7to9_rebuild_code.r
15) Part10h_CP10K_Prognosis_code.r
16) Part10h_CP10K_Prognosis_1_Comparator_code.r
17) Part10h_CP10K_Prognosis_2_CogVsPath_code.r
18) Part10h_CP10K_Prognosis_4_Multivariate_code.r
19) Part10h_CP10K_Comparator_Fig13_code.r
20) Part11_CP10K_CogVsPath_code.r
21) Part11_CP10K_CogVsPath_HighLevel_code.r
22) Part11_Identity_BandAblation_code.r
23) Part11_BandAblation_Fig11Export_code.r
24) Part12_02_Extract.py (Python)
25) Part12_04_Labels.r, Part12_05_Units.r, Part12_06_GerritsCells.r
26) Part12_Fig15 to Part12_Fig19, each Analysis then Plots
27) Part12_ModulePostIP_Overview.r
28) Part10g_CP10K_DEGs_Part1_code_S11_followup_analyses.r (reads the Part 12 labels)
29) Part13a_CP10K_Find14DpctMT_code.r, Part13_CP10K_M4Share_code.r, Part13b_CP10K_M4RandomSubset_code.r
30) Fig17_Summary_code.r

`verification/run_rest.sh` in the repository runs steps 15 to 30 in this order (or from any named step), saves each console log in `console_logs/` and ends with the comparison against the original results mounted at `/win`. Each step is marked failed if it stops or its log ends with `VERDICT: FAIL`.

`verification/run_python.sh` runs Python Parts 1 to 8 in order (or from any named step), with Part 8 once per arm, saves each console log in `console_logs_py/`, and ends with `compare_python.py`, which compares the outputs with the original results mounted at `/ref`. It stops at the first step that fails.

## 5. Notes

1) Fonts. Two figure scripts use Arial and Arial Narrow, which are Windows fonts. The container uses Liberation Sans and Liberation Sans Narrow instead, which have the same letter widths, so the layout is unchanged and the letter shapes differ slightly.
2) Optional packages. `pals` and `svglite` were not installed for the original runs and are not in the container. The scripts then use the built-in ocean.curl colours and R's own SVG device, as in the original runs.
3) Memory in Part12_06. `Part12_06_GerritsCells.r` calls cells in parallel, and each worker holds one raw Gerrits matrix in memory. By default it uses all logical cores but one. On a machine with many cores and limited memory, set fewer workers, for example `Rscript -e 'P12_WORKERS <- 6L; source("Part12_rerun/Part12_06_GerritsCells.r")'`. Each sample has its own seed, so the number of workers does not change the results.
4) Package versions. The container provides the software environment used at the end of the project. Some earlier analysis steps were run on older versions of some R packages, so minor numerical differences are possible. The R analyses for the paper also ran on Windows and the Python parts on Fedora, and a different operating system or linear algebra library can change the last digits of some results.
