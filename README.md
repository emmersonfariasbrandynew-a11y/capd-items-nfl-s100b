# Plasma neurofilament light chain, S100 calcium-binding protein B, and individual delirium screening items in mechanically ventilated children

Analysis code, data dictionary, and aggregated results for the exploratory secondary analysis of a prospective pediatric intensive care cohort (Fundação Santa Casa de Misericórdia do Pará, Belém, Brazil; June 2022 to December 2025), submitted for publication.

## What this analysis does
For 34 mechanically ventilated children with an assessable Cornell Assessment of Pediatric Delirium (CAPD) on study days 1 and 5, the analysis:
- compares the Spearman correlations of plasma neurofilament light chain (NfL) and S100 calcium-binding protein B (S100B) with each of the eight CAPD items, as signed differences (Δρ = ρ(item, NfL) − ρ(item, S100B));
- uses studentized bootstrap-t confidence intervals with jackknife standard errors (5,000 resamples by participant; seed 20260818);
- adjusts within each family of tests (Holm for the total score; Benjamini–Hochberg, q <0.10, for the item families);
- includes secondary analyses (total score, temporal differences, changes, contrasts with pain and withdrawal scores) and three sensitivity analyses.

## Repository contents
| Folder or file | Content | Access |
|---|---|---|
| `code/` | R scripts (00–09), the runner `EXECUTAR_ANALISE.R`, and the integrity manifest `MANIFESTO_EXECUCAO.csv` (SHA-256) | open (MIT) |
| `docs/data_dictionary.csv` | Definition, coding, unit, and missing-data meaning of each deposited variable | open (CC BY 4.0) |
| `docs/VERIFICATION_REPORT.md` | Reproducibility verification of this repository version | open (CC BY 4.0) |
| `docs/sessionInfo.txt` | R version and package versions of the verification run | open |
| `docs/sessionInfo_analysis_run.txt` | R session of the analysis run reported in the article | open |
| `docs/SHA256SUMS_open_files.txt`, `docs/SHA256_restricted_dataset.txt` | SHA-256 checksums of the open files and of the restricted dataset | open |
| `results/` | Aggregated outputs (CSV tables of estimates) of the analysis run reported in the article | open (CC BY 4.0) |
| `data_restricted/` | Anonymized dataset (45 participants, 50 variables) | **restricted**: deposited as a separate restricted-access Zenodo record and never in the GitHub repository; access on request (see below) |

## Anonymization of the restricted dataset
- Only the variables required to reproduce this analysis are included.
- Participant identifiers were replaced by random codes; the linkage key is not deposited.
- Age is provided in bands (<1, 1–5, ≥6 years), together with the developmental age band required to apply the CAPD anchor points (0–6, >6–12, >12–24, >24 months); exact age, dates, and free text are not included.
- Variables from other analyses of the cohort are not included.

## Reproducing the analysis
1. Obtain the restricted dataset and place it at `data_restricted/capd_items_biomarkers_anonymized.csv`.
2. Install R 4.6.1 with the packages `boot`, `digest`, `ggplot2`, and `patchwork` (see `docs/sessionInfo.txt`).
3. From the `code/` folder, run `Rscript EXECUTAR_ANALISE.R`. The runner verifies the SHA-256 of every script and of the dataset before running, and stops on any mismatch. The verification run took about 2 hours on a desktop computer (single R process).

With the repository dataset, all inferential results reproduce those of the article. The only expected difference is in the descriptive age row of Table 1: the article reports the median age in months, while the repository dataset provides age bands. The files in `results/` are the outputs of the analysis run; rerunning the code yields age-band rows in Table 1. Relative to the analysis run, the analytic code is unchanged: data input and age description were adapted to the repository dataset, code comments were removed, and console messages were edited (see `docs/VERIFICATION_REPORT.md`).

Each script begins with a short English description. Object names and console messages are in Portuguese, the working language of the study team; this README and the data dictionary are in English.

## Data access
Individual data remain under the custody of Fundação Santa Casa de Misericórdia do Pará. The ethics approval (no. 5,276,259; CAAE 45317221.1.0000.5171) and informed consent did not authorize unrestricted public release. Requests for the restricted dataset are assessed by the institution under the applicable approvals, without prior guarantee of access.

## License
- Code: MIT License (`LICENSE`).
- Documentation and aggregated results: Creative Commons Attribution 4.0 International (`LICENSE-CC-BY-4.0.txt`).
- Restricted dataset: access and reuse under the terms set by the institution.

## Citation
See `CITATION.cff`. Please also cite the article.
