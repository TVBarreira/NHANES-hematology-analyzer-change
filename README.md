# Analyzer change distorts blood-count aging and inflammation measures in NHANES

Code and results for the manuscript of that title (Barreira, Harmon, Costa, Armstrong, Heffernan). Three scripts reproduce every number and the figure; `results/` holds the output files from the run that produced the submitted manuscript.

| Step | Script | What it does |
|---|---|---|
| 1 | `M01_build_cbc_frame.R` | Downloads the NHANES demographic, complete blood count and biochemistry files for the 11 cycles 1999-2000 to 2021-2023 (`nhanesA`), applies the exclusions, defines the adolescent referent (12-19 years, nonpregnant), computes the composite indices, and saves `data/m_cbc_frame.rds`. Downloads are cached per cycle. 20-40 min the first time. |
| 2 | `M02_measurement_analyses.R` | All estimation: the change in each blood-count measure at every adjacent-cycle boundary in adolescents; the percentile test (additive vs proportional); the correction factors; phenotypic age (BioAge `phenoage_calc`, trained in NHANES III, CRP omitted) and its change at the analyzer boundary; the counterfactual correction and chemistry placebo tests; composite measures observed and corrected in adults; pooled-vs-within-cycle quartile reassignment. About 1 min. |
| 3 | `F01_figure.py` | Draws the figure from `results/m_all_boundary_steps.csv` and `results/m_figure_panels.csv`. |

## Setup

- Create a project folder with empty `data/` and `results/` subfolders, and either run the scripts from it or point the environment variable `NHANES_BIOAGE_DIR` to it.
- R 4.5.1 with `nhanesA`, `dplyr`, `tidyr`, `tibble`, `survey` (4.5), `BioAge` (0.1.0; `remotes::install_github("dayoonkwon/BioAge")`), and optionally `NoSleepR`.
- Python 3 with `numpy` and `matplotlib` for the figure.

## Where each reported number comes from (all in `results/`)

The files in `results/` are from the runs of 18-19 October 2026 (R 4.5.1). Re-running the scripts reproduces them; the public version of `M02` additionally writes the corrected acceleration series against the observed age fit (`RDW + WBC (same reference)`), which `F01_figure.py` uses when present.

| In the letter | File |
|---|---|
| Table 1 and Supplementary Table 1: change at each boundary, largest other change | `m_all_boundary_steps.csv` |
| Table 1: correction factors | `m_correction_factors.csv` |
| Supplementary Table 2: change by percentile | `m_correction_form_quantiles.csv` |
| Table 2: composite measures, observed and after correction (adults) | `m_derived_boundary_steps.csv` |
| Table 2: adults changing quartile | `m_quartile_shift.csv` (rows with `population == "adults"`) |
| Phenotypic age step (2.63 y), share removed (92.0%), residual (0.21 y) | `m_pheno_steps.csv`, `m_counterfactual.csv` |
| Fig. 1a / Fig. 1b | `m_all_boundary_steps.csv` / `m_figure_panels.csv` |
| Sample flow | `m_sample_ledger.csv`, `m_build_log.csv` |
| Every number in the text, by token | `m_placeholders.csv` |

## Notes

- Survey weights: `WTMEC2YR` for every cycle, divided by the number of cycles pooled; masked strata and PSUs for variance.
- The analyzer boundary is 2011-2012 to 2013-2014 (Beckman Coulter MAXM to DxH 800). NCHS indicates the instruments were replaced during 2012, so the estimated step is conservative.
- Phenotypic age is computed with `phenoage_calc` trained in NHANES III and projected with that fit. Published coefficients with the CRP term deleted give a different, uncalibrated measure about 0.8 years higher in level.
- In `m_figure_panels.csv` the corrected acceleration series appears twice: as its own pooled residual (`RDW + WBC`), which re-centres on the pooled mean, and against the observed age fit (`RDW + WBC (same reference)`), which is what the figure shows.

## License

MIT (see `LICENSE`).

## Data

No participant-level data are included. All inputs are public NHANES files, downloaded by `M01` from the National Center for Health Statistics (https://wwwn.cdc.gov/nchs/nhanes/).
