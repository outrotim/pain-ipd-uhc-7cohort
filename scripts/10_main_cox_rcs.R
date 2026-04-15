# =============================================================================
# Study 24: 10_main_cox_rcs.R
# Primary analysis -- Cox PH + RCS non-linear dose-response + cohort strata
# Outcomes: incident ADL / IADL / all-cause mortality (mortality in 13_mortality.R)
# Exposure: pain_severity_std (0-3) and pain_any_derived
# Strategy:
#   (A) Two-stage IPD: per-cohort Cox -> random-effects meta-analysis
#   (B) One-stage IPD: cohort-stratified baseline, common exposure effect
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(survival); library(rms)
  library(metafor); library(readr); library(broom); library(purrr)
})
source("_paths.R")  # defines DATA_DIR, TAB_DIR, FIG_DIR, log_step(); see _paths.R

d <- readRDS(file.path(DATA_DIR, "ipd_pooled.rds"))
log_step(sprintf("Loaded IPD: %d rows", nrow(d)), "10_cox")

# ---- Build survival dataset (person-level: first baseline -> first event or last wave) ----
surv_data <- d %>%
  filter(!is.na(age_years), !is.na(sex), !is.na(pain_any_derived)) %>%
  arrange(cohort, ID, wave) %>%
  group_by(cohort, ID) %>%
  mutate(
    baseline_age = first(age_years),
    baseline_pain_any = first(pain_any_derived),
    baseline_pain_sev = first(pain_severity_std),
    baseline_sex = first(sex),
    baseline_edu = first(edu_isced),
    adl_base = first(adl_disability),
    adl_event_wave = suppressWarnings(min(wave[adl_disability == 1], na.rm = TRUE)),
    adl_event = ifelse(is.finite(adl_event_wave) & adl_base == 0, 1L, 0L),
    iadl_base = first(iadl_disability),
    iadl_event_wave = suppressWarnings(min(wave[iadl_disability == 1], na.rm = TRUE)),
    iadl_event = ifelse(is.finite(iadl_event_wave) & iadl_base == 0, 1L, 0L),
    first_wave = first(wave),
    last_wave = last(wave),
    fu_waves = last_wave - first_wave + 1
  ) %>%
  ungroup()

log_step(sprintf("Survival data: %d person-waves", nrow(surv_data)), "10_cox")

# Person-level (one row per person)
person <- surv_data %>%
  distinct(cohort, ID, .keep_all = TRUE) %>%
  mutate(
    time_adl = ifelse(adl_event == 1,
                      adl_event_wave - first_wave,
                      last_wave - first_wave),
    time_iadl = ifelse(iadl_event == 1,
                       iadl_event_wave - first_wave,
                       last_wave - first_wave),
    .keep = "all"
  ) %>%
  filter(fu_waves > 0, !is.na(baseline_age), !is.na(baseline_pain_any))

log_step(sprintf("Person-level: %d", nrow(person)), "10_cox")

# ADL / IADL: condition on baseline-free
person_adl <- person %>% filter(adl_base == 0, time_adl > 0)
person_iadl <- person %>% filter(iadl_base == 0, time_iadl > 0)

log_step(sprintf("ADL analytic set: %d (events %d)", nrow(person_adl), sum(person_adl$adl_event)), "10_cox")
log_step(sprintf("IADL analytic set: %d (events %d)", nrow(person_iadl), sum(person_iadl$iadl_event)), "10_cox")

# ---- One-stage cohort-stratified Cox + RCS ---------------------------------
log_step("=== One-stage Cox + RCS ===", "10_cox")

run_stratified_cox <- function(dat, time_var, event_var, exposure_var, label) {
  f <- as.formula(sprintf("Surv(%s, %s) ~ %s + baseline_age + baseline_sex + strata(cohort)",
                          time_var, event_var, exposure_var))
  fit <- tryCatch(coxph(f, data = dat), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  s <- summary(fit)
  tibble(
    outcome = label,
    exposure = exposure_var,
    hr = round(s$coefficients[exposure_var, "exp(coef)"], 3),
    lci = round(s$conf.int[exposure_var, "lower .95"], 3),
    uci = round(s$conf.int[exposure_var, "upper .95"], 3),
    p = signif(s$coefficients[exposure_var, "Pr(>|z|)"], 3),
    n_obs = fit$n,
    n_event = fit$nevent
  )
}

res_main <- bind_rows(
  run_stratified_cox(person_adl %>% filter(!is.na(baseline_pain_any)),
                     "time_adl", "adl_event", "baseline_pain_any", "ADL disability"),
  run_stratified_cox(person_adl %>% filter(!is.na(baseline_pain_sev)),
                     "time_adl", "adl_event", "baseline_pain_sev", "ADL disability (severity)"),
  run_stratified_cox(person_iadl %>% filter(!is.na(baseline_pain_any)),
                     "time_iadl", "iadl_event", "baseline_pain_any", "IADL disability"),
  run_stratified_cox(person_iadl %>% filter(!is.na(baseline_pain_sev)),
                     "time_iadl", "iadl_event", "baseline_pain_sev", "IADL disability (severity)")
)
print(res_main)
write_csv(res_main, file.path(TAB_DIR, "cox_main_results.csv"))

# ---- RCS dose-response (severity 0-3) --------------------------------------
log_step("=== RCS dose-response (ADL) ===", "10_cox")

run_rcs <- function(dat, time_var, event_var, label) {
  dd <- dat %>%
    filter(!is.na(baseline_pain_sev), !is.na(baseline_age), !is.na(baseline_sex)) %>%
    as.data.frame()
  if (nrow(dd) < 100) return(NULL)
  dd_r <- dd
  dd_r$cohort_factor <- factor(dd_r$cohort)
  # Cox with RCS on pain severity; 3 knots since variable has only 4 levels (0-3)
  f <- as.formula(sprintf("Surv(%s, %s) ~ rcs(baseline_pain_sev, 3) + baseline_age + baseline_sex + strata(cohort_factor)",
                          time_var, event_var))
  fit <- tryCatch(coxph(f, data = dd_r), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  pred_grid <- expand.grid(
    baseline_pain_sev = 0:3,
    baseline_age = mean(dd_r$baseline_age),
    baseline_sex = "Female",
    cohort_factor = levels(dd_r$cohort_factor)[1]
  )
  lp <- predict(fit, pred_grid, type = "lp")
  tibble(
    outcome = label,
    pain_severity = 0:3,
    log_HR = round(lp - lp[1], 4),
    HR = round(exp(lp - lp[1]), 3)
  )
}

rcs_adl <- run_rcs(person_adl, "time_adl", "adl_event", "ADL disability")
rcs_iadl <- run_rcs(person_iadl, "time_iadl", "iadl_event", "IADL disability")
rcs_all <- bind_rows(rcs_adl, rcs_iadl)
print(rcs_all)
write_csv(rcs_all, file.path(TAB_DIR, "rcs_dose_response.csv"))

# ---- Two-stage: per-cohort Cox -> random-effects meta-analysis -------------
log_step("=== Two-stage meta-analysis ===", "10_cox")

run_cohort_cox <- function(dat, time_var, event_var, exposure_var) {
  dat %>%
    filter(!is.na(!!sym(exposure_var))) %>%
    group_by(cohort) %>%
    group_split() %>%
    map_dfr(function(dd) {
      if (nrow(dd) < 100 || sum(dd[[event_var]]) < 10) return(NULL)
      f <- as.formula(sprintf("Surv(%s, %s) ~ %s + baseline_age + baseline_sex",
                              time_var, event_var, exposure_var))
      fit <- tryCatch(coxph(f, data = dd), error = function(e) NULL)
      if (is.null(fit)) return(NULL)
      s <- summary(fit)
      coef_row <- tryCatch(s$coefficients[exposure_var, ], error = function(e) NULL)
      if (is.null(coef_row)) return(NULL)
      tibble(
        cohort = unique(dd$cohort),
        logHR = coef_row["coef"],
        SE = coef_row["se(coef)"],
        HR = round(exp(coef_row["coef"]), 3),
        lci = round(exp(coef_row["coef"] - 1.96 * coef_row["se(coef)"]), 3),
        uci = round(exp(coef_row["coef"] + 1.96 * coef_row["se(coef)"]), 3),
        n = fit$n,
        n_event = fit$nevent
      )
    })
}

meta_adl_any <- run_cohort_cox(person_adl, "time_adl", "adl_event", "baseline_pain_any")
meta_adl_sev <- run_cohort_cox(person_adl, "time_adl", "adl_event", "baseline_pain_sev")
meta_iadl_sev <- run_cohort_cox(person_iadl, "time_iadl", "iadl_event", "baseline_pain_sev")

summarize_meta <- function(tbl, label) {
  if (nrow(tbl) < 2) return(NULL)
  fit <- rma(yi = logHR, sei = SE, data = tbl, method = "REML")
  tibble(
    outcome = label,
    n_cohorts = fit$k,
    pooled_HR = round(exp(fit$beta[1, 1]), 3),
    lci = round(exp(fit$ci.lb), 3),
    uci = round(exp(fit$ci.ub), 3),
    I2 = round(fit$I2, 1),
    tau2 = round(fit$tau2, 4),
    p = signif(fit$pval, 3)
  )
}

meta_summary <- bind_rows(
  summarize_meta(meta_adl_any, "ADL disab ~ any pain"),
  summarize_meta(meta_adl_sev, "ADL disab ~ pain severity"),
  summarize_meta(meta_iadl_sev, "IADL disab ~ pain severity")
)
print(meta_summary)
write_csv(meta_summary, file.path(TAB_DIR, "meta_analysis_pooled.csv"))

write_csv(meta_adl_any, file.path(TAB_DIR, "cohort_cox_adl_any.csv"))
write_csv(meta_adl_sev, file.path(TAB_DIR, "cohort_cox_adl_sev.csv"))
write_csv(meta_iadl_sev, file.path(TAB_DIR, "cohort_cox_iadl_sev.csv"))

log_step("=== Primary analysis done ===", "10_cox")
cat("\nOutputs:\n")
cat("- cox_main_results.csv\n")
cat("- rcs_dose_response.csv\n")
cat("- meta_analysis_pooled.csv\n")
cat("- cohort_cox_adl_*.csv\n")
