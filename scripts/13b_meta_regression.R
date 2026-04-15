# =============================================================================
# Study 24: 13b_meta_regression.R
# Cohort-level meta-regression to decompose I^2 ~ 97% between-cohort heterogeneity.
# Moderators:
#   - GDP per capita (World Bank, mid-observation year)
#   - UHC service coverage index (WHO SDG 3.8.1, 2021)
#   - baseline ADL disability prevalence (within-data)
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(metafor); library(readr); library(ggplot2); library(stringr)
})
source("_paths.R")  # defines DATA_DIR, TAB_DIR, FIG_DIR, log_step(); see _paths.R

# ---- Cohort macro covariates (World Bank + WHO, all public) -----------------
cohort_macro <- tribble(
  ~cohort, ~country_region,     ~gdp_pc_2015, ~uhc_index_2021, ~income_group,
  "CHARLS", "China",              8050,        81,            "Upper-middle",
  "ELSA",   "UK",                 44450,       87,            "High",
  "HRS",    "USA",                56790,       85,            "High",
  "KLoSA",  "South Korea",        27230,       88,            "High",
  "LASI",   "India",              1590,        63,            "Lower-middle",
  "MHAS",   "Mexico",             9520,        75,            "Upper-middle",
  "SHARE",  "Europe (pooled)",    36000,       82,            "High"
)
write_csv(cohort_macro, file.path(TAB_DIR, "cohort_macro_covariates.csv"))

# ---- Load per-cohort Cox results (produced by 10_main_cox_rcs) ---------------
cox_adl <- read_csv(file.path(TAB_DIR, "cohort_cox_adl_any.csv"), show_col_types=FALSE) %>%
  mutate(outcome = "ADL disability")
cox_adl_sev <- read_csv(file.path(TAB_DIR, "cohort_cox_adl_sev.csv"), show_col_types=FALSE) %>%
  mutate(outcome = "ADL disability (severity)")
cox_mort <- read_csv(file.path(TAB_DIR, "mortality_per_cohort.csv"), show_col_types=FALSE) %>%
  mutate(outcome = "All-cause mortality")

# Baseline ADL prev (from pooled_cohort_summary)
coh_base <- read_csv(file.path(TAB_DIR, "pooled_cohort_summary.csv"), show_col_types=FALSE) %>%
  select(cohort, baseline_adl_disab = pct_adl_disab, baseline_pain = pct_pain, mean_age_coh = mean_age)

all_cohort_effects <- bind_rows(cox_adl, cox_adl_sev, cox_mort) %>%
  left_join(cohort_macro, by = "cohort") %>%
  left_join(coh_base, by = "cohort")

write_csv(all_cohort_effects, file.path(TAB_DIR, "meta_regression_input.csv"))

# ---- Meta-regression for each outcome ---------------------------------------
run_meta_regression <- function(dat, outcome_label) {
  dat <- dat %>% filter(!is.na(logHR), !is.na(SE))
  if (nrow(dat) < 3) return(NULL)

  out <- list()

  # (1) total heterogeneity
  m0 <- rma(yi = logHR, sei = SE, data = dat, method = "REML")
  out$total <- tibble(
    outcome = outcome_label,
    model = "Total (no moderator)",
    k = m0$k,
    pooled_HR = round(exp(m0$beta[1,1]), 3),
    lci = round(exp(m0$ci.lb), 3),
    uci = round(exp(m0$ci.ub), 3),
    I2 = round(m0$I2, 1),
    tau2 = round(m0$tau2, 4),
    Q_p = signif(m0$QEp, 3)
  )

  # (2) UHC index
  if (all(!is.na(dat$uhc_index_2021))) {
    m_uhc <- rma(yi = logHR, sei = SE, mods = ~ uhc_index_2021, data = dat, method = "REML")
    out$uhc <- tibble(
      outcome = outcome_label,
      model = "~ UHC index",
      k = m_uhc$k,
      I2_residual = round(m_uhc$I2, 1),
      R2 = round(m_uhc$R2, 1),
      beta_mod = round(m_uhc$beta[2,1], 4),
      mod_SE = round(m_uhc$se[2], 4),
      mod_p = signif(m_uhc$pval[2], 3)
    )
  }

  # (3) log-GDP per capita
  if (all(!is.na(dat$gdp_pc_2015))) {
    dat$log_gdp <- log(dat$gdp_pc_2015)
    m_gdp <- rma(yi = logHR, sei = SE, mods = ~ log_gdp, data = dat, method = "REML")
    out$gdp <- tibble(
      outcome = outcome_label,
      model = "~ log(GDP per capita)",
      k = m_gdp$k,
      I2_residual = round(m_gdp$I2, 1),
      R2 = round(m_gdp$R2, 1),
      beta_mod = round(m_gdp$beta[2,1], 4),
      mod_SE = round(m_gdp$se[2], 4),
      mod_p = signif(m_gdp$pval[2], 3)
    )
  }

  # (4) Baseline ADL prev
  if (all(!is.na(dat$baseline_adl_disab))) {
    m_base <- rma(yi = logHR, sei = SE, mods = ~ baseline_adl_disab, data = dat, method = "REML")
    out$baseline <- tibble(
      outcome = outcome_label,
      model = "~ Baseline ADL prev",
      k = m_base$k,
      I2_residual = round(m_base$I2, 1),
      R2 = round(m_base$R2, 1),
      beta_mod = round(m_base$beta[2,1], 4),
      mod_SE = round(m_base$se[2], 4),
      mod_p = signif(m_base$pval[2], 3)
    )
  }
  bind_rows(out)
}

mr_adl <- run_meta_regression(
  all_cohort_effects %>% filter(outcome == "ADL disability"),
  "ADL disability ~ any pain"
)
mr_mort <- run_meta_regression(
  all_cohort_effects %>% filter(outcome == "All-cause mortality"),
  "All-cause mortality ~ any pain"
)
mr_all <- bind_rows(mr_adl, mr_mort)
print(mr_all)
write_csv(mr_all, file.path(TAB_DIR, "meta_regression_results.csv"))

# ---- Leave-one-out for UHC (pre-specified sensitivity) ----------------------
run_loo <- function(dat, moderator) {
  dat <- dat %>% filter(!is.na(logHR), !is.na(.data[[moderator]]))
  map_dfr(seq_len(nrow(dat)), function(i) {
    excl <- dat$cohort[i]
    dd <- dat[-i, ]
    f <- as.formula(sprintf("~ %s", moderator))
    m <- tryCatch(rma(yi = logHR, sei = SE, mods = f, data = dd, method = "REML"),
                  error = function(e) NULL)
    if (is.null(m)) return(NULL)
    tibble(excluded = excl,
           beta = round(m$beta[2, 1], 4),
           SE   = round(m$se[2], 4),
           p    = signif(m$pval[2], 3),
           R2   = round(m$R2, 1))
  })
}
adl_eff <- all_cohort_effects %>% filter(outcome == "ADL disability")
loo_uhc <- run_loo(adl_eff, "uhc_index_2021")
print(loo_uhc)
write_csv(loo_uhc, file.path(TAB_DIR, "A3_loo_meta_regression.csv"))

# ---- Bubble plots -----------------------------------------------------------
make_bubble <- function(dat, moderator, mod_label, outcome_label, fname) {
  dat <- dat %>% filter(!is.na(logHR), !is.na(.data[[moderator]]))
  if (nrow(dat) < 3) return(NULL)
  p <- ggplot(dat, aes(x = .data[[moderator]], y = HR, size = 1/SE^2, color = cohort)) +
    geom_point(alpha = 0.7) +
    geom_text(aes(label = cohort), vjust = -1.2, size = 3.5, color = "black") +
    geom_hline(yintercept = 1, linetype = "dashed") +
    geom_smooth(method = "lm", se = TRUE, color = "grey30", linewidth = 0.7,
                show.legend = FALSE, inherit.aes = FALSE,
                aes(x = .data[[moderator]], y = HR, weight = 1/SE^2)) +
    scale_y_log10() +
    labs(title = sprintf("Meta-regression: %s ~ %s", outcome_label, mod_label),
         x = mod_label, y = "Cohort HR (log scale)", size = "Precision (1/SE^2)") +
    theme_minimal(base_size = 12) + theme(legend.position = "bottom")
  ggsave(file.path(FIG_DIR, fname), p, width = 9, height = 6, dpi = 300)
}

make_bubble(adl_eff, "uhc_index_2021", "UHC Index (2021)", "ADL disability", "fig_meta_reg_uhc.png")
make_bubble(adl_eff, "gdp_pc_2015", "GDP per capita 2015 (USD)", "ADL disability", "fig_meta_reg_gdp.png")
make_bubble(adl_eff, "baseline_adl_disab", "Baseline ADL prevalence (%)", "ADL disability", "fig_meta_reg_base.png")

log_step("=== meta-regression done ===", "13b_mr")
cat("\nOutputs:\n")
cat("- meta_regression_results.csv\n")
cat("- A3_loo_meta_regression.csv\n")
cat("- cohort_macro_covariates.csv\n")
cat("- figures/fig_meta_reg_*.png\n")
