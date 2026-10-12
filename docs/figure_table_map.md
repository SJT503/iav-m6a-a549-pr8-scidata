# Figure and table map

Every figure and table in the Data Descriptor is produced by the scripts below, from the
files deposited in GEO (GSE336465) plus the provider quality-control reports in
`data/qc_evidence/`. Stages run in the order given; `run_all.sh` runs all of them.

| Output | Script(s) | Input | Notes |
|---|---|---|---|
| **Figure 1** — study design and data workflow | `41_figure1_schematic.py` | `figdata/fig2d_ip_matrix.csv` (probe count only) | schematic drawn entirely in code (matplotlib primitives and text); contains no raster or generated artwork |
| **Figure 2** — array and replicate quality | `30_export_figdata.R` → `40_figures.py` | `data/raw/Matrix.xlsx`, FE files, `data/qc_evidence/RNA_QC_gel_image_original.jpeg` | panel a is the provider's RNA gel; PCA is a standard centred PCA (`prcomp`) on the 40,955 probes |
| **Figure 3** — m6A quantity and expression | `10_limma.R`, `30_export_figdata.R` → `40_figures.py` | `data/raw/Matrix.xlsx` | limma + Benjamini–Hochberg, deposited normalisation |
| **Figure 4** — methylation fraction (%modified) | `12_detection_flags.py`, `13_pctmod_detected.R`, `30_export_figdata.R` → `40_figures.py` | FE files, `data/raw/Matrix.xlsx` | restricted to probes with a Sup signal above background on all six arrays |
| **Table 1** — data records, sizes and MD5 | `61_build_table_docx.py` | `data/raw_md5.txt`, `data/processed_md5.txt`, the deposited files themselves | checksums recomputed and compared against the manifest |
| **Table 2** — per-sample QC | `60_build_tables.py`, `61_build_table_docx.py` | `data/qc_evidence/*`, FE file headers | NanoDrop, spike-in recovery, labelling |
| **Table 3** — array layout and calibration | `60_build_tables.py`, `61_build_table_docx.py` | FE file headers, `data/raw/Matrix.xlsx` | per-array position, scan image, calibration constants |
| **Table 4** — interferon and ISG transcripts | `50_table4_data.py`, `61_build_table_docx.py` | `figdata/fig3_limma_ip_sup.csv`, FE files | one probe per gene: the one with the largest expression change |
| **Supplementary Table 1** — normalisation-robust %modified probes (n = 210) | `11_normalisation_sensitivity.R`, `12_detection_flags.py`, `13_pctmod_detected.R`, `60_build_tables.py` | `data/raw/Matrix.xlsx`, FE files | supplied as an xlsx file; definitions on the Legend sheet |

## Stage list

| Stage | Script | What it does | Needs |
|---|---|---|---|
| 00 | `00_prepare_inputs.py` | caches the provider result tables as CSV | `data/processed/` |
| 01 | `01_prepare_inputs_large.py` | streaming cache of the two largest tables | `data/processed/` |
| 02 | `02_prepare_inputs_wide.py` | wide cache of the mRNA quantity and expression tables | `data/processed/` |
| 03 | `03_probe_map_mrna.py` | probe → transcript / gene / m6A-region map, and the probe-filter rule | FE files, cache |
| 04 | `04_probe_map_lncrna_others.py` | probe maps for lncRNA and other non-coding RNA, and the matrix calibration factors | FE files, cache |
| 10 | `10_limma.R` | limma + BH on the m6A-quantity, expression and %modified layers | cache |
| 11 | `11_normalisation_sensitivity.R` | deposited / centred / joint-QN / separate-QN normalisations, decomposition of the offsets | cache |
| 12 | `12_detection_flags.py` | per-probe Sup and IP detection flags from the FE files | FE files, and the stage-11 normalisation results (`limma_pct_SI.csv`, `limma_pct_QNj.csv`, `pctmod_robust_set.csv`) |
| 13 | `13_pctmod_detected.R` | %modified restricted to probes with a measurable Sup signal; counterpart and sensitivity analyses | cache, flags |
| 20 | `20_pca_standard.R` | standard centred PCA (IP and Sup) | cache |
| 21 | `21_isg_detection.py` | detection status of the interferon transcripts in both channels | FE files, map |
| 30 | `30_export_figdata.R` | one small CSV per figure panel group | results, cache |
| 40 | `40_figures.py` | renders Figures 2–4 (PDF + PNG, 300 dpi) | figdata |
| 41 | `41_figure1_schematic.py` | draws Figure 1, the study-design schematic (PDF + PNG, 300 dpi) | figdata (probe count) |
| 50 | `50_table4_data.py` | Table 4 values, one probe per gene | figdata, FE files |
| 60 | `60_build_tables.py` | Supplementary Table 1 (xlsx) and the data for Tables 2–3 | results, cache, QC reports |
| 61 | `61_build_table_docx.py` | Table 1–4 as Word files | results, data root |
| 70 | `70_check_numbers.py` | re-derives every number quoted in the manuscript text and compares | all of the above |

`tools/make_manifest.py` regenerates `MANIFEST.md5`, which records checksums of the
distributed files. It is run by hand when the repository changes, not by `run_all.sh`.

## What the code does not do

The array-service provider (Aksomics, Shanghai, China) performed RNA quality control,
immunoprecipitation, labelling, hybridisation, scanning, Agilent Feature Extraction,
calibration-spike-in normalisation, probe filtering and the differential tables of the
original submission. Those steps are **not** reimplemented here: their outputs are the
deposited files that this repository reads. Everything reported in the revised manuscript
— limma differential analysis, %modified calculation and testing, standard PCA, replicate
correlations, the decomposition of the normalisation constants and the sensitivity
analyses, the interferon panel, and all figures and tables — is computed by the code in
this repository from those deposited files.
