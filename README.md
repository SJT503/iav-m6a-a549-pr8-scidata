# IAV m6A A549 PR8 — analysis code

Reproduces every analysis, figure and table in the Data Descriptor
**"An m6A epitranscriptomic microarray dataset of human A549 cells infected with influenza A virus PR8"**
(*Scientific Data*).

The dataset profiles human A549 cells under mock infection or influenza A virus PR8 infection
for 24 h at MOI 0.5, with three biological replicates per condition, on GPL25759
(Agilent-085753 Arraystar Human Epitranscriptomic microarray). Raw and processed data:
**GEO [GSE336465](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE336465)**.

Code released under the MIT License. If you use it, please cite the Data Descriptor and this
archive (see `CITATION.cff`).

---

## 1. What you need

| | |
|---|---|
| Python | 3.14.4 — numpy, pandas, scipy, openpyxl, matplotlib, Pillow, python-docx, PyMuPDF |
| R | 4.5.3 — limma (3.66.0), jsonlite (2.0.0) |
| Data | the GEO deposit (see step 2) |

```bash
conda env create -f environment.yml     # or install the packages listed in environment.yml
```

## 2. Get the data

Download the Series from GEO and arrange it as shown. File names must keep their
GEO names, because the configuration refers to them.

```bash
mkdir -p data/raw data/processed
# from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE336465
#   supplementary files  ->  data/processed/
#   the six Feature Extraction files (MOCK_1.txt ... PR8_3.txt) and Matrix.xlsx  ->  data/raw/
```

```
data/
├── raw/                 MOCK_1.txt … PR8_3.txt, Matrix.xlsx
├── processed/           mRNA_Quantity(flag_pass).xlsx, … (8 provider tables)
├── qc_evidence/         provider QC reports — already in this repository
├── raw_md5.txt          checksums of the uploaded raw files
└── processed_md5.txt    checksums of the uploaded processed files
```

Sample metadata — the mapping from sample name to GEO accession, condition, replicate,
array position and scan image — is in [`metadata/sample_metadata.csv`](metadata/sample_metadata.csv).

Check the download before running anything:

```bash
cd data && md5sum -c raw_md5.txt && md5sum -c processed_md5.txt
```

`data/qc_evidence/` holds four provider quality-control reports (NanoDrop quantification,
spike-in qPCR, the Agilent labelling report and the RNA gel image). They are **not** part of
the GEO deposit; they are shipped here because Figure 1d and Table 2 are built from them.

## 3. Configure

All paths and thresholds live in **`config/config.ini`**; no script contains a
machine-specific path. The defaults assume the layout above. If your data lives elsewhere,
change one line:

```ini
[paths]
data_root = /path/to/your/data      # relative paths resolve against the repository root
```

Sanity-check the configuration and the inputs:

```bash
python scripts/_config.py
```

It prints the resolved paths and lists any expected input file that is missing.

The R stages read the same `config/config.ini` through `scripts/_config.R`, so one file configures both languages.

## 4. Run

```bash
bash run_all.sh                 # everything, in order
bash run_all.sh 10 11           # only the stages whose number matches a prefix
```

Outputs: `results/` (statistical tables and summaries), `figdata/` (one small CSV per
figure panel group), `figures/` (Figures 1–4, PDF + PNG at 300 dpi), `tables/`
(Tables 1–4 as Word files, Supplementary Table 1 as xlsx).

[`docs/figure_table_map.md`](docs/figure_table_map.md) maps every figure and table to the
scripts that produce it, and lists what each stage does.

## 5. What the code does and does not do

The array-service provider (Aksomics, Shanghai, China) performed RNA quality control,
immunoprecipitation, labelling, hybridisation, scanning, Agilent Feature Extraction,
calibration-spike-in normalisation, probe filtering and the differential tables of the
original submission. Those steps are **not** reimplemented here; their outputs are the
deposited files this repository reads.

Everything reported in the manuscript is computed here from those deposited files:
limma differential analysis with Benjamini–Hochberg adjustment, the %modified calculation
and its testing on the logit scale, standard centred PCA, replicate correlations, the
decomposition of the array calibration constants, the normalisation sensitivity analyses,
the interferon panel, and all figures and tables.

## 6. Conventions that matter when reusing the data

**The provider's differential tables use the MOCK-versus-PR8 convention.** Fold change is
mock/PR8, so a transcript induced by infection has a fold change **below 1** and appears in
the `hypo` or `down` sheet. This is the opposite of the usual infected-versus-control
direction and is the single easiest way to misread the deposit.

| File | Sheet | Meaning |
|---|---|---|
| `mRNA_Quantity(diffmethylated).xlsx` | `hyper_MOCK_vs_PR8` | m6A quantity higher in mock, i.e. **lower after infection** |
| | `hypo_MOCK_vs_PR8` | m6A quantity lower in mock, i.e. **higher after infection** |
| `mRNA_ExpressionLevel(diffexpressed).xlsx` | `up_MOCK_vs_PR8` | expression higher in mock, i.e. **lower after infection** |
| | `down_MOCK_vs_PR8` | expression lower in mock, i.e. **higher after infection** |
| `lncRNA_Quantity(diffmethylated).xlsx`, `others_Quantity(diffmethylated).xlsx` | `hyper_` / `hypo_` | as for mRNA quantity |

Checked against the deposited values: IFNL1 (`ENST00000333625`) has a log2 fold change of
−6.64 in `hypo_MOCK_vs_PR8` and −5.21 in `down_MOCK_vs_PR8`, consistent with strong induction
by PR8 (about 2,700-fold in the expression layer; Table 4).

**%modified needs a measurable denominator.** %modified = 100 × 2^IP/(2^IP + 2^Sup) is
undefined when the Sup (Cy3) signal is at background. Analyses of %modified in the
manuscript are therefore restricted to the 25,279 probes whose Sup signal is above
background (`gIsWellAboveBG = 1`) on all six arrays; `results/probe_detection_flags.csv`
gives the per-array flags so that you can apply your own rule.

**Array-level offsets.** The provider's normalisation adds one constant per channel per
array. Because that constant is the same for every probe on an array, it shifts %modified
for all probes on that array, and the IP-minus-Sup constant was on average 0.26 log2 units
higher on the PR8 arrays. A global change in m6A stoichiometry therefore cannot be
separated from a difference in immunoprecipitation efficiency in this dataset, and global
shifts in %modified should not be read biologically. `results/norm_offsets_table.csv` lists
the constants; stage 11 quantifies the effect and reproduces the sensitivity analyses.

## 7. Repository layout

```
config/config.ini     every path and threshold
scripts/              the workflow, in run order (see docs/figure_table_map.md)
run_all.sh            runs all stages
metadata/             sample metadata (name, accession, condition, array, scan image)
data/                 inputs (GEO data, provider QC reports, checksum manifests)
docs/                 figure/table map, and the manuscript text used by the number check
                      (figure_table_map.md, manuscript_text.md)
results/  figdata/  figures/  tables/     outputs
```

## 8. Reference

Ritchie ME, et al. limma powers differential expression analyses for RNA-sequencing and
microarray studies. *Nucleic Acids Res.* 2015;43(7):e47.
Benjamini Y, Hochberg Y. Controlling the false discovery rate. *J R Stat Soc Series B*
1995;57(1):289–300.
