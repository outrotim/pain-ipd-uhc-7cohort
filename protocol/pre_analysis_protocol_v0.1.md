# Study 24 — PROTOCOL v0.1
# Pain, disability, cognitive decline and mortality in 7 ageing cohorts: a harmonised IPD pooled analysis of >250 000 participants

> **Drafted**: 2026-04-14
> **Status**: PROTOCOL v0.1 — pre-specified
> **Reporting**: PRISMA-IPD (Stewart 2015) + STROBE + GATHER

---

## 1. Background

Chronic pain prevalence in adults ≥50 years is 30–50% and is a leading contributor to geriatric disability and care burden (GBD 2021 ranked low back pain first among YLD causes in older adults). Existing evidence has three systematic gaps:

1. **Geographic bias**: 95% of longitudinal evidence comes from high-income white cohorts (HRS/ELSA); East Asian, South Asian, and Latin American older populations are underrepresented.
2. **Methodological uniformity**: Most studies use linear Cox/logistic models, missing biologically plausible non-linear pain–disability dose–response (thresholds, saturation, J-shape).
3. **Mediation pathway unquantified**: The relative contribution of depression, physical function and social participation to the pain → disability pathway has never been decomposed in a cross-cultural framework.

This study uses the Gateway to Global Aging Data (g2aging.org) harmonisation system to pool seven nationally representative ageing cohorts (CHARLS, HRS, ELSA, KLoSA, LASI, MHAS, SHARE) under a strict IPD pooled framework.

---

## 2. Objectives & hypotheses

### 2.1 Primary
Non-linear dose–response of pain intensity and multisite pain count to:
- (O1) incident ADL disability
- (O2) incident IADL disability
- (O3) cognitive decline rate + MCI/dementia incidence (subsample)
- (O4) all-cause mortality

### 2.2 Secondary
- (S1) Mediation proportion of depression, physical function (grip + gait), social participation
- (S2) Cross-country heterogeneity + meta-regression on GDP per capita and UHC index
- (S3) Sex, age, education subgroup interactions

### 2.3 Hypotheses
- H1: Baseline pain predicts incident ADL disability (HR > 1.5)
- H2: Pain severity shows dose–response
- H3: Country-level UHC moderates the pain–disability association (direction not pre-specified)

---

## 3. PICO-T

| Dimension | Definition |
|---|---|
| **P** | Community-dwelling adults ≥50 y, non-demented, with valid pain measurement |
| **E** | Pain intensity (NRS 0-10 or 4-level ordinal) + multisite pain count (0, 1, 2, 3, ≥4) |
| **C** | No pain reference / mildest category |
| **O** | ADL, IADL disability; cognitive decline; all-cause mortality |
| **T** | Baseline to last available wave (CHARLS 2020, HRS 2022, ELSA 2023, etc.) |

---

## 4. Data sources & authorisation

| Cohort | Portal |
|---|---|
| CHARLS | charls.pku.edu.cn |
| HRS | hrs.isr.umich.edu (RAND HRS Fat File) |
| ELSA | UK Data Service |
| KLoSA | survey.keis.or.kr |
| LASI | iipsindia.ac.in/lasi |
| MHAS | mhasweb.org |
| SHARE | share-eric.eu |

All seven datasets were accessed under their respective Data Use Agreements. Each constituent cohort obtained original ethics approval from its originating institution; this secondary analysis was determined exempt from additional institutional review.

---

## 5. Harmonisation plan

Core principle: Gateway to Global Aging Data naming scheme is primary; missing variables use operational harmonisation (equivalent local-variable mapping).

### 5.1 Exposure harmonisation

| Concept | Gateway name | CHARLS | HRS | ELSA | KLoSA | LASI | MHAS | SHARE |
|---|---|---|---|---|---|---|---|---|
| Pain present | `r*painfr` | da041 | R*PAIN | painfr | C062 | ht801 | AGE22_1 | ph010 |
| Pain intensity | `r*painlv` | da042 | R*PAINLV | painlv | C063 | ht801a | AGE22_2 | ph011 |
| Multisite | `r*painpt*` | da041_1..8 | local map | pain_sites | C062_1..8 | ht801b | AGE22_3 | ph089_* |

### 5.2 Outcomes
- ADL (6-item Katz): `r*adlfive` / `r*adl6a`
- IADL (5-6 item Lawton): `r*iadlfour` / `r*iadlza`
- Cognition: `r*cogtot` + subdomains (orientation/memory/executive)
- Death: `radyear` / `radmonth`

### 5.3 Mediators & covariates
- Depression (CES-D): `r*cesd` (8 or 10 item, Z-score)
- Grip strength: `r*mxgrip`
- Gait speed: `r*wspeed1`
- Social participation: cohort-specific scales → PCA first component
- Education: `raedyrs` + ISCED coding
- Income: within-country quintiles

### 5.4 QC
- Each harmonised variable reports missingness, distribution, extreme-value checks
- `data/harmonization_qc_report.html` generated

---

## 6. Statistical analysis plan

### 6.1 Primary
- Cohort-stratified one-stage Cox proportional hazards model, adjusting for baseline age and sex
- Two-stage REML random-effects meta-analysis of cohort-specific log-HRs as validation
- RCS (4 knots at 5/35/65/95 percentile) for non-linear dose–response

### 6.2 Moderator (H3, pre-specified)
- Random-effects meta-regression on UHC Service Coverage Index (WHO 2021, SDG 3.8.1), log-GDP per capita, baseline ADL prevalence
- Leave-one-out sensitivity for UHC
- Bonferroni correction for three moderators

### 6.3 Sensitivity block (8 pre-specified)
| # | Analysis | Method |
|---|---|---|
| S1 | Full-covariate | Cox + HTN/DM/HD/stroke/arthritis/smoking/education |
| S2 | Competing risk | Fine-Gray subdistribution HR (two directions) |
| S3 | IPTW | `WeightIt`, trimmed 1/99 |
| S4 | Multiple imputation | MICE (m=10 per cohort) |
| S5 | IPCW | Attrition-adjusted (3 cohorts ≥7 waves) |
| S6 | Time-varying exposure | Counting-process Cox |
| S7 | Washout ≥2 waves (~4 y) | Reverse-causation sensitivity |
| S8 | Washout ≥4 waves (~8 y) | Reverse-causation sensitivity |

### 6.4 Supporting analyses
- Regression-based mediation (product method, 500 bootstraps) for depression and grip strength
- Negative-control outcome: hearing worsening
- E-value (VanderWeele & Ding 2017)
- Cognitive LMM: fixed effects for pain × time
- Pooled PAF via inverse-variance-weighted Levin's formula
- Counterfactual UHC simulation (illustrative, not causally identified)

### 6.5 Software
R 4.4.0+. Packages: `rms`, `survival`, `metafor`, `mice`, `WeightIt`, `haven`, `gtsummary`, `ggplot2`. Code versioned in Git.

---

## 7. Reporting plan

### 7.1 Main tables
- Table 1: Baseline characteristics across 7 cohorts
- Table 2: Primary HR + 8 sensitivity analyses
- Table 3: Meta-regression moderators

### 7.2 Main figures
- Figure 1: Forest plots (ADL/IADL/mortality)
- Figure 2: Meta-regression bubble plot (UHC × log-HR)
- Figure 3: Global HR map
- Appendix figures: dose–response, KM, leave-one-out, counterfactual

---

## 8. Ethics & data governance

- All constituent studies obtained original ethics approval
- Secondary-analysis framework determined exempt from additional institutional review
- No individual-level data are redistributed
- Code will be archived under MIT License; aggregate products under CC-BY 4.0

---

## 9. Version history
- v0.1 (2026-04-14): Initial protocol draft, pre-specified prior to primary analysis execution

---

**Note**: Local-institution and author-identifying information has been redacted from this public copy during double-anonymised peer review.
