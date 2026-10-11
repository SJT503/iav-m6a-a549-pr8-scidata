# Data directory

The inputs to this repository are the files deposited in GEO **GSE336465**, plus four
provider quality-control reports that are shipped in `qc_evidence/`.

```
data/
├── raw/                 MOCK_1.txt … PR8_3.txt  (Feature Extraction exports)
│                        Matrix.xlsx            (normalised 40,955-probe matrix)
├── processed/           the eight provider result tables:
│                          mRNA_Quantity(flag_pass).xlsx
│                          mRNA_ExpressionLevel(flag_pass).xlsx
│                          lncRNA_Quantity(flag_pass).xlsx
│                          others_Quantity(flag_pass).xlsx
│                          mRNA_Quantity(diffmethylated).xlsx
│                          lncRNA_Quantity(diffmethylated).xlsx
│                          others_Quantity(diffmethylated).xlsx
│                          mRNA_ExpressionLevel(diffexpressed).xlsx
├── qc_evidence/         provider QC reports (shipped here, not in GEO):
│                          RNA_QC_nanoDrop_quantification.csv
│                          MeRIP_QC_spikein_pct_input.csv
│                          Agilent Labeling Report.pdf
│                          RNA_QC_gel_image_original.jpeg
├── raw_md5.txt          MD5 checksums of the uploaded raw files
└── processed_md5.txt    MD5 checksums of the uploaded processed files
```

`interim/` is created at run time by `scripts/00_prepare_inputs.py`; the CSV caches it
writes are regenerated from the xlsx files and are not distributed.

## Downloading

Download the Series supplementary files from
<https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE336465> and place them as shown
above, then verify:

```bash
cd data && md5sum -c raw_md5.txt && md5sum -c processed_md5.txt
```

GEO stores the six Feature Extraction files gzip-compressed and renames them to
`GSM<accession>_<sample>.txt.gz`; decompress them and restore the sample names
(`MOCK_1.txt` … `PR8_3.txt`) that the configuration and the metadata table use.

## Source

Array service and data generation: Aksomics (Shanghai, China), on the Arraystar Human
mRNA&lncRNA Epitranscriptomic Microarray (8×60K), GEO platform GPL25759. The provider
performed RNA quality control, immunoprecipitation, labelling, hybridisation, scanning,
Feature Extraction, calibration-spike-in normalisation, probe filtering and the
differential tables; those outputs are the deposited files above.
