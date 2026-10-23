# Container for the SEA-AD CEP-IP/CMRF analysis

The container holds the software used for the paper: R 4.6.1 with the packages in `renv.lock`, Python 3.13.14 with the packages in `requirements-python.txt`, and the deposited scripts. It does not hold any data. The SEA-AD folder is mounted at `/data` when the container starts, and every script reads and writes inside that folder.

The scripts find the project folder through the environment variable `SEAAD_ROOT`, which the container sets to `/data`. Outside the container, set `SEAAD_ROOT` to the SEA-AD folder, or start R or Python from that folder.

## 1. Build

Podman is used here (it is installed on Fedora). Docker accepts the same commands with `docker` in place of `podman`.

The build folder holds:

```
Containerfile
seaad-setup.sh
seaad-run.sh
renv.lock
requirements-python.txt
scripts/R/          the deposited R scripts, the Part 12 JSON settings and the Part 12
                    pre-specified readings (Part12_Preset_Readings*.md)
scripts/Python/     the deposited Python scripts and the Mic.13 gene table (Mic13_genes.txt)
scripts/SuppTables/ the scripts that build Supplementary Tables S1 to S19
verification/       compare_rerun.py, compare_python.py and compare_gamdata.py
```

From that folder:

```bash
podman build --no-cache -t seaad-cepip -f Containerfile .
```

The first build takes about one to two hours, mostly to compile the Bioconductor packages and monocle3. It needs an internet connection. The image is about 6 to 8 GB. After a change to the scripts only, the same command without `--no-cache` rebuilds only the last layers.

Instead of building it, the image can be loaded from the archived copy on Zenodo (`seaad-cepip_image.tar.gz`, https://doi.org/10.5281/zenodo.23446228), which holds the same software and scripts:
```bash
podman load -i seaad-cepip_image.tar.gz
```
This gives the image `localhost/seaad-cepip`; Docker accepts `docker load -i` in the same way.

For an HPC cluster that uses Apptainer instead of Podman:

```bash
podman save -o seaad-cepip.tar seaad-cepip
apptainer build seaad-cepip.sif docker-archive://seaad-cepip.tar
```

## 2. Run the whole analysis with one command

```bash
podman run --rm --shm-size=16g -v /path/to/SEA-AD:/data:Z seaad-cepip seaad-run all
```

`:Z` lets the container use the folder on systems with SELinux, such as Fedora. `--shm-size=16g` is needed for Part 7 (DecontX), whose worker processes share the count matrix through shared memory; with the default size, Part 7 stops with a bus error.

`seaad-run all` runs `seaad-setup`, then every script in the order of section 4, from the SEA-AD h5ad to Figure 17 and Supplementary Tables S1 to S19. It takes about six hours on the development machine. Other modes run part of the analysis:

| Command | Runs | Needs |
|---|---|---|
| `seaad-run all` | Python Parts 1 to 8, R Parts 9 to 13, Figure 17, Tables S1 to S19 | the inputs of section 3, items 1 to 3 |
| `seaad-run python` | Python Parts 1 to 8, Tables S1 to S6 | the SEA-AD h5ad |
| `seaad-run r` | R Parts 9 to 13, Figure 17, Tables S7 to S19 | items 1 to 4 of section 3 |
| `seaad-run hf` | R Parts 9 to 11, Tables S7 to S10 and S12 | item 4 of section 3 only (the Hugging Face files) |
| `seaad-run tables` | Tables S1 to S19 from outputs already in the data folder | the outputs of the earlier parts |
| `seaad-run list r` | prints the steps of a mode | - |

Each step's console output is saved to `seaad-run_logs/<step>.txt`, and the start, end and duration of every step to `seaad-run_logs/seaad-run_summary.txt`. The run stops at the first step that fails, because later steps read earlier outputs, and prints the command to resume from that step, for example:

```bash
podman run --rm --shm-size=16g -v /path/to/SEA-AD:/data:Z seaad-cepip seaad-run all --from Part10e_CP10K_CMRF_code
```

A step fails when it exits with an error, or when its log shows a Python traceback, a failed check (`[FAIL]`, `GUARDRAIL FAILED`, `VERDICT: FAIL`) or a stop (`Run stopped`, `STOPPED`, `ERROR:`, `Execution halted`).

With `--compare`, the results are then compared with the original results mounted read-only at `/ref` (a folder holding the original Part folders and `Manuscript/Supp_Tables`, the published Supplementary Tables):

```bash
podman run --rm --shm-size=16g -v /path/to/SEA-AD:/data:Z -v /path/to/original/SEA-AD:/ref:ro,Z \
    seaad-cepip seaad-run all --compare
```

The comparisons are written to `compare/` (pipeline tables) and `Manuscript/Supp_Tables/compare_supp_tables_log.txt` (Supplementary Tables).

To work inside the container instead, start it with a shell and run `seaad-setup` once:

```bash
podman run --rm -it --shm-size=16g -v /path/to/SEA-AD:/data:Z seaad-cepip
seaad-setup
```

`seaad-setup` copies the shared settings files (`CEP_arms.R`, `R_config.R`, `R_config_HKG8.R`) and the Mic.13 gene table (`Mic13_genes.txt`) into the data folder, and the Part 12 scripts, settings and pre-specified readings into `Part12_rerun`, if they are missing. The Part 12 scripts are run from `Part12_rerun`, because each one checks that its own copy and `Part12_03_Config.r` are there. It reports which inputs it finds and does not change any file that is already there.

## 3. Inputs

Place these in the data folder before running:

1) The SEA-AD immune h5ad, `SEAAD_Immune_HIP_MEC_LEC_ITG_MTG_FI_STG_DFC_AnG_V1C_RNAseq_final-nuclei.2026-06-22.h5ad`, read by Parts 1 to 8 and Part12_02, and the SEA-AD donor metadata workbook, `MetaData/sea-ad_cohort_donor_metadata_072524.xlsx`, read by Parts 1d and 1e.
2) The SEA-AD quantitative neuropathology workbook (`sea-ad-quantitative-neuropathology-063026.xlsx`) and the CPS workbook (`sea-ad_global_and_local_cps_050126.xlsx`), used by the Batteries.
3) The Gerrits et al. deposit (GSE148822) in `Gerrits/`, with the series matrix in `Gerrits/Metadata/`, used by Part12_06 and Battery 1.
4) To run from Part 9 without Parts 1 to 8: the Part 8 export (`R_export_CP10K/`), the Part 2 donor tables (`*_S5b_donor_clinical.csv`) and the Part 7 contamination table (`Part7_DecontX/decontX_contamination_all.csv`), all on Hugging Face.

The Mic.13 gene table (`Mic13_genes.txt`, from Green et al. 2024, Extended Data Fig. 2c), read by Part 5b, is deposited with the Python scripts, and `seaad-setup` copies it into the data folder.

The Supplementary Tables are not inputs. The scripts in `scripts/SuppTables/` build all 19 of them from the pipeline outputs, and the two scripts that read tables (the Part 10g follow-up reads Table S11, and `Fig17_Summary_code.r` reads Table S12 and the five Part 12 workbooks of Tables S13 to S17) run after those tables are built.

Each script lists its own inputs in its header.

## 4. Run order

`seaad-run list all` prints the full order. Run every script from the data folder (`/data`). R scripts are run with `source()`, as in the original runs. Each script writes a log next to its outputs. Part 8 writes one arm per run; the arm is selected at the top of the script, and `seaad-run` runs it once per arm (PVM3_all, PVM1_all, Mic13_PVM3, HKG8_matched).

1) Python Parts 1a to 8 (`/opt/seaad/Python`), then Tables S1 to S6
2) Part9_CP10K_QC_GAM_code.r
3) Part9b_CP10K_SiKi_code.r
4) Part10_CP10K_SKCS_code.r
5) Part10b_CP10K_IPRS_code.r
6) Part10b_CP10K_SKCS_TREP_code.r
7) Part10c_CP10K_Copula1_code.r
8) Part10d_CP10K_Copula2_code.r
9) Part10e_CP10K_CMRF_code.r
10) Part10e_CP10K_Validation1_code.r
11) Part10e_CP10K_Validation2_UMAP_code.r
12) Part10f_CP10K_Band_Formation_code.r
13) Part10f_CP10K_CMRF_RSkfCV_Mechanism_code.r
14) Part10g_CP10K_DEGs_Part1_code.r
15) Part10g_Step7to9_rebuild_code.r
16) Part10h_CP10K_Prognosis_code.r
17) Part10h_CP10K_Prognosis_1_Comparator_code.r
18) Part10h_CP10K_Prognosis_2_CogVsPath_code.r
19) Part10h_CP10K_Prognosis_4_Multivariate_code.r
20) Part10h_CP10K_Comparator_Fig13_code.r
21) Part11_CP10K_CogVsPath_code.r
22) Part11_CP10K_CogVsPath_HighLevel_code.r
23) Part11_Identity_BandAblation_code.r
24) Part11_BandAblation_Fig11Export_code.r
25) Tables S7 to S10 and S12
26) Part12_02_Extract.py (Python)
27) Part12_04_Labels.r, Part12_05_Units.r, Part12_06_GerritsCells.r
28) Part12_Fig15 to Part12_Fig19, each Analysis then Plots
29) Part12_ModulePostIP_Overview.r
30) Table S11, then Part10g_CP10K_DEGs_Part1_code_S11_followup_analyses.r (reads the Part 12 labels and Table S11)
31) Part13a_CP10K_Find14DpctMT_code.r, Part13_CP10K_M4Share_code.r, Part13b_CP10K_M4RandomSubset_code.r
32) Table S18 (runs Part13c_CP10K_M4_SuppConsolidate_code.r, which also draws Supplementary Figures S20 and S21), then Tables S13 to S17
33) Fig17_Summary_code.r
34) Table S19

To run steps by hand, for example:

```bash
python /opt/seaad/Python/Part1a_NormForensics_code.py
Rscript -e 'source("/opt/seaad/R/Part9_CP10K_QC_GAM_code.r")'
python Part12_rerun/Part12_02_Extract.py
Rscript -e 'source("Part12_rerun/Part12_04_Labels.r")'
bash /opt/seaad/SuppTables/run_supp_tables.sh 7-10,12
```

`run_supp_tables.sh` builds the listed tables. Each table has a builder (`build_SuppTable_S<n>`), whose workbook and log go to `Manuscript/Supp_Tables/build/`. `finalize_supp_tables.py` then writes the published workbook (`Manuscript/Supp_Tables/Supplementary_Table_S<n>_*.xlsx`) with the published sheet names, the Contents sheet (from `supp_tables_contents.json`) and the published wording. With `PUBLISHED=/ref/Manuscript/Supp_Tables`, it also runs `compare_supp_tables.py` against the published tables.

`verification/run_rest.sh` and `verification/run_python.sh` in the repository are the run scripts used for the reproducibility check of the paper, before `seaad-run` was written.

## 5. Notes

1) Fonts. Two figure scripts use Arial and Arial Narrow, which are Windows fonts. The container uses Liberation Sans and Liberation Sans Narrow instead, which have the same letter widths, so the layout is unchanged and the letter shapes differ slightly.
2) Optional packages. `pals` and `svglite` were not installed for the original runs and are not in the container. The scripts then use the built-in ocean.curl colours and R's own SVG device, as in the original runs.
3) Memory in Part12_06. `Part12_06_GerritsCells.r` calls cells in parallel, and each worker holds one raw Gerrits matrix in memory. `seaad-run` uses 6 workers; set another number with `-e P12_WORKERS_06=4` on the `podman run` line. By hand: `Rscript -e 'P12_WORKERS <- 6L; source("Part12_rerun/Part12_06_GerritsCells.r")'`. Each sample has its own seed, so the number of workers does not change the results.
4) Package versions. The container provides the software environment used at the end of the project. Some earlier analysis steps were run on older versions of some R packages, so minor numerical differences are possible. The R analyses for the paper also ran on Windows and the Python parts on Fedora, and a different operating system or linear algebra library can change the last digits of some results.
5) Supplementary Tables from a rerun. Compared with the published tables, a container rerun differs in the values that already differ in the pipeline outputs: the Validation2 UMAP coordinates, kBET and neighbour enrichment in Table S10, the depth-thinning curve in Table S13 (sheet P), and the order of a few near-tied genes in Table S3 (sheets C and D). The Table S10 builder accepts these UMAP-based differences and records them in its log. Table S19 records the date each source workbook was last modified.
