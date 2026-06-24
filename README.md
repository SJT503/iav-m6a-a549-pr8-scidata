# IAV m6A A549 PR8 Scientific Data code release

This repository contains the analysis and figure-generation code used for the manuscript **"An m6A epitranscriptomic microarray dataset of human A549 cells infected with influenza A virus PR8"**.

## Dataset

The dataset profiles human A549 cells under mock infection or influenza A virus PR8 infection for 24 h at MOI 0.5, with three biological replicates per condition, using the Arraystar Human mRNA&lncRNA Epitranscriptomic Microarray platform.

The raw and processed data are deposited in GEO. The GEO Series accession and reviewer-access token will be added after GEO curation.

## Repository contents

- `scripts/` — Python scripts used for data anchoring, quality assessment, figure generation, manuscript/table generation, and reference checks.
- `docs/` — project governance and descriptor framing documents.
- `results/anchored_facts.md` — manuscript-facing factual anchors extracted from the source report files.
- `data/README.md` — expected input files and rerun instructions.
- `requirements.txt` — Python package requirements.

## Reproducibility notes

The scripts preserve the exact analysis provenance used during manuscript preparation. Several scripts contain local path constants at the top of the file. To rerun the workflow on another machine:

1. Download the GEO raw and processed files after the GEO Series is available.
2. Place them under a local data directory.
3. Edit the path constants near the top of each script to point to the downloaded files and desired output directory.
4. Install dependencies:

```bash
pip install -r requirements.txt
```

5. Run scripts in numerical order where applicable.

## Data and code separation

This repository intentionally excludes raw microarray files, processed spreadsheets, GEO upload archives, and generated binary outputs. These files are provided through GEO or regenerated from the scripts.

## License

Code is released under the MIT License.
