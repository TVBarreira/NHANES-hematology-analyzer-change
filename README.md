# An uncorrected hematology analyzer change distorts blood-count measures of aging and inflammation in NHANES

Code and results for the manuscript of that title (Barreira, Harmon, Costa, Armstrong, Heffernan), submitted to *PLOS Medicine*. The scripts reproduce every number and Fig 1. The `m_*.csv` files in `results/` are the output files from the runs behind the manuscript.

| Step | Script | What it does |
|---|---|---|
| 1 | `M01_build_cbc_frame.R` | Downloads the NHANES demographic, complete blood count and biochemistry files for the 11 cycles 1999-2000 to 2021-2023 (`nhanesA`), applies the exclusions, defines the adolescent referent (12-19 years, nonpregnant), computes the composite indices, and saves `data/m_cbc_frame.rds`. Downloads are cached per cycle. 20-40 min the first time. |
| 2 | `M02_measurement_analyses.R` | All main estimation: the change in each blood-count measure at every adjacent-cycle boundary in adolescents; the percentile test (additive vs proportional); the correction factors; phenotypic age (BioAge `phenoage_calc`, trained in NHANES III, CRP omitted) and its change at the analyzer boundary; the counterfactual correction and chemistry placebo tests; composite measures observed and corrected in adults; pooled-vs-within-cycle quartile reassignment as measured. About 1 min. |
| 3 | `M03_verification.R` | The age-group-by-cycle interaction test; the headline estimates with the subgroup analysed as a survey domain; all 2021-2023 results repeated with the phlebotomy weight (`WTPH2YR`); quartile reassignment and top-quartile shares after correction. Writes `results/m_verification.csv`. 5-10 min. |
| 4 | `F01_figure_PLOS.py` | Draws Fig 1 (600 dpi TIFF, 7.5 in wide) from `m_all_boundary_steps.csv` and `m_figure_panels.csv` in the local `results/` folder. |

`F01_figure.py` is the earlier version of the same figure drawn for a research-letter format; it uses the same inputs.

## Setup

- Create a project folder with empty `data/` and `results/` subfolders, and either run the scripts from it or point the environment variable `NHANES_BIOAGE_DIR` to it.
- R 4.5.1 with `nhanesA`, `dplyr`, `tidyr`, `tibble`, `survey` (4.5), `BioAge` (0.1.0; `remotes::install_github("dayoonkwon/BioAge")`), and optionally `NoSleepR`.
- Python 3 with `numpy`, `matplotlib` and `Pillow` for the figure.

## Where each reported number comes from

| In the manuscript | File |
|---|---|
| Table 1 and S1 Table: change at each boundary, largest other change | `m_all_boundary_steps.csv` |
| Table 1: correction factors | `m_correction_factors.csv` |
| S2 Table: change by percentile | `m_correction_form_quantiles.csv` |
| Table 2: composite measures, observed and after correction (adults) | `m_derived_boundary_steps.csv` |
| Table 2: adults changing quartile, as measured | `m_quartile_shift.csv` (rows with `population == "adults"`) |
| Adults changing quartile after correction; share in the top pooled quartile before and from 2013, measured and corrected | `m_verification.csv` (section `D quartiles`, version `MEC weight (as in M02)`) |
| Age-group interaction for the phenotypic age step (p = 0.47) | `m_verification.csv` (section `A interaction`, first row) |
| Results with the 2021-2023 phlebotomy weight (e.g., Fig 1B 2021-2023: 1.05 vs 1.02 years) | `m_verification.csv` (section `C phlebotomy weight`) |
| Phenotypic age step (2.63 y), share removed (92.0%), residual (0.21 y) | `m_pheno_steps.csv`, `m_counterfactual.csv` |
| Fig 1A / Fig 1B | `m_all_boundary_steps.csv` / `m_figure_panels.csv` |
| Sample flow | `m_sample_ledger.csv`, `m_build_log.csv` |

## Notes

- Survey weights: `WTMEC2YR` for every cycle, divided by the number of cycles pooled; masked strata and PSUs for variance. `M03` repeats the analyses involving 2021-2023 with the phlebotomy weight that NCHS recommends for blood analytes in that cycle; results are essentially unchanged.
- The interaction term that `M02` writes (age as a continuous variable by cycle, across all 11 cycles) is a different test from the one the manuscript describes and reports; the reported test is in `M03`.
- The analyzer boundary is 2011-2012 to 2013-2014 (Beckman Coulter MAXM to DxH 800). NCHS indicates the instruments were replaced during 2012, so the estimated step is conservative.
- Phenotypic age is computed with `phenoage_calc` trained in NHANES III and projected with that fit. Published coefficients with the CRP term deleted give a different, uncalibrated measure about 0.8 years higher in level.
- In `m_figure_panels.csv` the corrected acceleration series appears as its own pooled residual (`RDW + WBC`), which re-centres on the pooled mean. The figure scripts express it against the observed age fit, so the two series coincide before 2013-2014.

## License

MIT (see `LICENSE`).

## Data

No participant-level data are included. All inputs are public NHANES files, downloaded by `M01` from the National Center for Health Statistics (https://wwwn.cdc.gov/nchs/nhanes/).
