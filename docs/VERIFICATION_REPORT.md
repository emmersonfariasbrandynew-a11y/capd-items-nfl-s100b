# Reproducibility verification

**Date:** 2026-11-08. **Environment:** R 4.6.1 on Windows 11 (see `sessionInfo.txt`).

**What was run:** the complete pipeline in `code/` exactly as deposited (`Rscript EXECUTAR_ANALISE.R`, scripts 00–09), on the restricted anonymized dataset `capd_items_biomarkers_anonymized.csv` (SHA-256 `763e2673267b58128d0de888f3be893170b158847b4d9fac376b1fe7ceb6e830`). Before running, the runner checked the SHA-256 of every script and of the dataset against `code/MANIFESTO_EXECUCAO.csv`. Seeds, resample numbers, test families, and thresholds were those of the analysis reported in the article. The run started at 21:43 and finished at 23:34 without errors; all bootstrap contrasts were computed with status "ok".

**Comparison:** every aggregated output table (31 CSV files) was compared byte by byte with the outputs of the analysis run reported in the article (`results/`):
- 30 of 31 tables were **identical**, including point estimates, bootstrap confidence intervals, bootstrap p values, and adjusted q values;
- Table T1 (sample characteristics) differs only in the age rows. The article reports median age in months, whereas the repository dataset provides age bands (<1, 1–5, ≥6 years). All other rows of T1 are identical.

**Code version verified:** relative to the analysis run, the deposited scripts differ only in (i) data input and age description adapted to the repository dataset, (ii) removal of code comments and addition of a short English header to each script, and (iii) edited console messages. The removal of comments and the edits to console messages were checked with the R parser: outside comments and printed text, the executable code is identical.

**Notes:**
- Bootstrap resampling depends on row order. The anonymized dataset therefore keeps the row order of the analysis dataset, while participant identifiers were replaced by random codes.
- Individual-level outputs (trajectory tables and figures) are produced locally by the pipeline but are not part of this repository.
