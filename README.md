# pain-ipd-uhc-7cohort

Minimal reproducibility package for the seven-cohort harmonised individual-participant-data (IPD) pooled longitudinal analysis of baseline pain and incident ADL disability in older adults, with country-level Universal Health Coverage (UHC) as a pre-specified moderator.

> Companion code archive for the submitted manuscript. During double-anonymised peer review, author and institutional information are intentionally omitted.

---

## Repository contents

| Path | Content |
|---|---|
| `protocol/pre_analysis_protocol_v0.1.md` | Pre-specified analysis protocol (v0.1, 14 April 2026). |
| `scripts/_paths.R` | Portable path / helper stub (env variables). |
| `scripts/10_main_cox_rcs.R` | Primary analysis: cohort-stratified one-stage Cox. |
| `scripts/13b_meta_regression.R` | Random-effects REML meta-regression on UHC; leave-one-out. |
| `scripts/15_figures_main.R` | Figures 1–3 redraw code. |
| `aggregates/cohort_logHR_uhc.csv` | Seven-row country-level aggregate (log-HR, SE, UHC, log-GDP, baseline prev). |
| `LICENSE` | MIT (code). |
| `LICENSE-data.md` | CC-BY 4.0 (aggregate data). |

---

## Data availability (layered declaration)

This repository distributes **no individual-level data**. Data access responsibilities remain with the original cohort custodians.

| Tier | Content | Status | How to obtain |
|---|---|---|---|
| 1. Country-level aggregate | `aggregates/cohort_logHR_uhc.csv` (7 rows) | **Open (CC-BY 4.0)** | This repo |
| 2. Cohort-level intermediates | Per-cohort event counts, sensitivity tables | Not in repo | Reproducible by running scripts on source cohorts |
| 3. Individual participant data | Harmonised person-wave tables from CHARLS / ELSA / HRS / KLoSA / LASI / MHAS / SHARE | **Not distributable** under respective DUAs | Apply directly to each cohort's portal; use Gateway to Global Aging Data (g2aging.org) |

Country-level moderators (UHC Service Coverage Index, GDP per capita) are publicly sourced from the WHO 2023 UHC monitoring report and the World Bank.

---

## How to reproduce

1. Apply for IPD access from each of the seven cohorts (see `protocol/`).
2. Obtain Gateway-harmonised RAND-style files and place them under a local directory.
3. Set environment variables:
   ```bash
   export STUDY24_DATA_DIR=/your/local/path/to/harmonised_rds
   export STUDY24_TAB_DIR=/your/local/path/to/intermediate_csv
   ```
4. Run scripts in order: `10_main_cox_rcs.R` → `13b_meta_regression.R` → `15_figures_main.R`.
5. The meta-regression (β, UHC 95.6% R², leave-one-out) can be reproduced directly from `aggregates/cohort_logHR_uhc.csv` + `13b_meta_regression.R` **without any IPD access**.

**R version**: 4.4.0. Packages: `metafor`, `survival`, `rms`, `dplyr`, `tidyr`, `ggplot2`, `patchwork`.

---

## Caveats

- Minimal reproducibility package, not a full pipeline. Harmonisation, sensitivity, mediation, cognitive LMM, IPCW/MICE/IPTW, Fine-Gray and negative-control scripts are available from the corresponding author upon reasonable request.
- Meta-regression with k = 7 is statistically fragile; leave-one-out β 0.031–0.038 is the basis for inference stability. R² = 95.6% is not causal attribution.
- UHC direction (higher-UHC → stronger observed HR) is interpreted as ascertainment heterogeneity, not true risk elevation.
- External validation in TILDA, JSTAR, ELSI-Brazil, TLSA is warranted once Gateway modules become available.

---

## Citation

> The Study 24 Authors. Universal Health Coverage statistically accounts for 95.6% of the between-country heterogeneity in the pain–disability association among older adults: a seven-cohort harmonised individual participant data pooled analysis. Submitted, 2026.

Full citation (authors, DOI, journal) will be added upon acceptance.

---

## License

- Code (`scripts/*.R`): MIT — see `LICENSE`.
- Aggregate data (`aggregates/*.csv`): CC-BY 4.0 — see `LICENSE-data.md`.

---

## Contact

During double-anonymised peer review, please contact the journal editorial office. Post-acceptance, corresponding-author details will be added here.
