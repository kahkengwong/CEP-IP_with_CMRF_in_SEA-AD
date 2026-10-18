# Run logs

## container_rerun/

Console logs of the container rerun on Fedora (Podman), one per script, from Part 9 to `Fig17_Summary_code.r`. The logs for Part 9 to Part 10g come from the batch runs made before `verification/run_rest.sh`. The logs from Part 10h onward come from `run_rest.sh`, and `run_rest_summary.txt` lists the status of each of its steps. Inside the container, the data folder is `/data` and the original results are `/win`. The run ID in the Part 12 logs ends with the host name.

## original_R/

Logs written by the original runs of the scripts from Part 9 to `Fig17_Summary_code.r`, in the same folder layout as the analysis folder. For each log that the container rerun wrote, the matching original log is included. Where a script wrote the same log at several times, the latest one is kept. Paths to the project folder are replaced with `<SEA-AD>`, and other user folder paths with `C:/Users/<user>` or `C:\Users\<user>`. No other text was changed.

## original_python/

Logs of the original runs of the Python Parts 1 to 8. The project folder path is replaced with `<SEA-AD>`, and the shell prompt user and host with `user@host`. No other text was changed.

Some scripts were renamed after they were run, so some log names differ from the script names.

| Log | Script | Note |
|---|---|---|
| `Part1a-1c_logs.txt` | `Part1a_NormForensics_code.py`, `Part1b_MT-audit_code.py`, `Part1c_WithinDonor_code.py` | |
| `Part1d_logs.txt` | `Part1d_Errata_code.py` | |
| `Part1e_logs.txt` | `Part1e_TechQC_code.py` | |
| `Part1f_logs.txt` | `Part1f_Ordinals_code.py` | |
| - | `Part1g_TrendReconcile_code.py` | No log was kept. |
| `Part2_NeoC_Logs.txt` | `Part2_NeoC_CP10K_code.py` | |
| `Part2_NonNeoC_Logs.txt` | `Part2_NonNeoC_CP10K_code.py` | |
| `Part2_M4Composition_log.txt` | `Part3_M4Composition_code.py` | Log named after the earlier script name. |
| `Part6c_NeoC_ScreenScatter_log.txt` | `Part4_NeoC_ScreenScatter.py` | Log named after the earlier script name. |
| `Part6c_NonNeoC_ScreenScatter_log.txt` | `Part4_NonNeoC_ScreenScatter.py` | Log named after the earlier script name. |
| - | `Part5_NeoC_DEGs.py` | No log was kept. |
| `Part5_logs.txt` | `Part5_NonNeoC_DEGs.py` | |
| `Part5b_logs.txt` | `Part5b_NeoC_Validation.py`, `Part5b_NonNeoC_Validation.py` | |
| `Part6_logs.txt` | `Part6_NeoC_GOI-DFG_code.py`, `Part6_NonNeoC_GOI-DFG_code.py` | |
| `Part6b_logs.txt` | `Part6b_NeoC_Pooled_GOI-DFG.py`, `Part6b_NonNeoC_Pooled_GOI-DFG.py` | Console capture of both runs. |
| `Part6b_NeoC_Pooled_GOI-DFG_log.txt` | `Part6b_NeoC_Pooled_GOI-DFG.py` | Log file written by the script. It also captured later steps run in the same session, including steps that are not part of this deposit. |
| `Part6b_NonNeoC_Pooled_GOI-DFG_log.txt` | `Part6b_NonNeoC_Pooled_GOI-DFG.py` | As above. |
| `Part7_logs.txt` | `Part7_DecontX.py` | |
| `Part8_Export_logs.txt` | `Part8_Data-export-CP10K_MT-CO3_8DFG.py` | One run per arm. At the time of these runs, the `PVM3_all` and `PVM1_all` folders were named `PVM3_bulk` and `PVM1_bulk`. |
| `Part8_Guardrail_Log_<arm>.txt` | `Part8_Data-export-CP10K_MT-CO3_8DFG.py` | Guardrail log of each arm, written as `Part8_Guardrail_Log.txt` in the arm folder of `R_export_CP10K`. |
