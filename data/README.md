# Data directory

Raw and processed data are not distributed in this code repository. They should be obtained from the GEO Series accession associated with the manuscript after curation/release.

Expected inputs include:

- `Matrix.xlsx` — normalized probe intensity matrix.
- Six Agilent Feature Extraction raw text files: `MOCK_1.txt`, `MOCK_2.txt`, `MOCK_3.txt`, `PR8_1.txt`, `PR8_2.txt`, `PR8_3.txt`.
- Processed Arraystar/Aksomics result tables used for descriptive analyses and figure generation.

The original scripts in `scripts/` preserve the analysis provenance and contain the local paths used during manuscript preparation. To rerun on a new machine, edit the path constants near the top of each script to point to the downloaded GEO files and a local output directory.
