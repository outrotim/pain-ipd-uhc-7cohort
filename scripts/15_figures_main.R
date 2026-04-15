# =============================================================================
# Study 24: 15_figures_main.R
# Main figures: Fig 1 forest plots / Fig 2 categorical dose-response /
#               Fig 3 cohort-level prevalence / Fig 4 Kaplan-Meier by severity
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
  library(readr); library(survival); library(scales); library(stringr)
})
source("_paths.R")  # defines DATA_DIR, TAB_DIR, FIG_DIR, log_step(); see _paths.R

d <- readRDS(file.path(DATA_DIR, "ipd_pooled.rds"))
cohort_cox_adl <- read_csv(file.path(TAB_DIR, "cohort_cox_adl_any.csv"), show_col_types = FALSE)
cohort_cox_adl_sev <- read_csv(file.path(TAB_DIR, "cohort_cox_adl_sev.csv"), show_col_types = FALSE)
cohort_cox_iadl_sev <- read_csv(file.path(TAB_DIR, "cohort_cox_iadl_sev.csv"), show_col_types = FALSE)
meta_pooled <- read_csv(file.path(TAB_DIR, "meta_analysis_pooled.csv"), show_col_types = FALSE)

theme_pub <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 13)
  )
cohort_colors <- c(CHARLS="#C3423F", ELSA="#2E86AB", HRS="#A23B72",
                   KLoSA="#F18F01", LASI="#048A81", MHAS="#9B5094", SHARE="#3E5641")

# ---- Fig 1: Forest plots (ADL~any / ADL~severity / IADL~severity) ----------
make_forest <- function(cohort_df, pool_row, title_text) {
  fx <- cohort_df %>%
    mutate(label = sprintf("%s (n=%s, ev=%s)",
                           cohort, format(n, big.mark=","), format(n_event, big.mark=","))) %>%
    arrange(HR)
  pool <- tibble(label = "Pooled (RE)",
                 cohort = "Pooled",
                 HR = pool_row$pooled_HR,
                 lci = pool_row$lci, uci = pool_row$uci,
                 is_pool = TRUE)
  fx$is_pool <- FALSE
  all <- bind_rows(fx %>% select(label, cohort, HR, lci, uci, is_pool), pool)
  all$label <- factor(all$label, levels = c(pool$label, rev(fx$label)))

  ggplot(all, aes(x = HR, y = label, color = cohort)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
    geom_point(aes(size = is_pool, shape = is_pool)) +
    geom_errorbarh(aes(xmin = lci, xmax = uci), height = 0.2) +
    geom_text(aes(label = sprintf("%.2f (%.2f-%.2f)", HR, lci, uci)),
              hjust = -0.2, vjust = -0.8, size = 3.2, color = "grey30") +
    scale_x_log10(breaks = c(0.5, 1, 1.5, 2, 3, 5), limits = c(0.5, 6)) +
    scale_size_manual(values = c(`FALSE` = 3, `TRUE` = 5), guide = "none") +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 18), guide = "none") +
    scale_color_manual(values = c(cohort_colors, Pooled = "black")) +
    labs(title = title_text, x = "Hazard Ratio (log scale)", y = NULL) +
    theme_pub + theme(legend.position = "none")
}

p1a <- make_forest(cohort_cox_adl,
                   meta_pooled %>% filter(outcome == "ADL disab ~ any pain"),
                   "A. ADL disability ~ baseline any pain")
p1b <- make_forest(cohort_cox_adl_sev,
                   meta_pooled %>% filter(outcome == "ADL disab ~ pain severity"),
                   "B. ADL disability ~ pain severity (per 1-unit)")
p1c <- make_forest(cohort_cox_iadl_sev,
                   meta_pooled %>% filter(outcome == "IADL disab ~ pain severity"),
                   "C. IADL disability ~ pain severity (per 1-unit)")

fig1 <- p1a / p1b / p1c + plot_annotation(
  title = "Figure 1. Cohort-specific and pooled Hazard Ratios",
  subtitle = "Random-effects IPD meta-analysis across 7 international cohorts",
  caption = "Age- and sex-adjusted Cox models; random-effects pooled HR with 95% CI"
)
ggsave(file.path(FIG_DIR, "fig1_forest_plots.png"), fig1, width = 10, height = 11, dpi = 300)
ggsave(file.path(FIG_DIR, "fig1_forest_plots.pdf"), fig1, width = 10, height = 11)
log_step("Fig 1 forest done", "15_fig")

# ---- Fig 2: Categorical dose-response (severity 0-3 -> HR) -----------------
baseline <- d %>% group_by(cohort, ID) %>% slice(1) %>% ungroup() %>%
  mutate(pain_sev_cat = factor(pain_severity_std, levels = 0:3,
                                labels = c("None","Mild","Moderate","Severe")))

surv <- d %>% arrange(cohort, ID, wave) %>%
  group_by(cohort, ID) %>%
  summarise(
    adl_base = first(adl_disability),
    adl_event = as.integer(any(adl_disability == 1, na.rm = TRUE) & first(adl_disability) == 0),
    iadl_base = first(iadl_disability),
    iadl_event = as.integer(any(iadl_disability == 1, na.rm = TRUE) & first(iadl_disability) == 0),
    time = n_distinct(wave),
    baseline_age = first(age_years),
    baseline_sex = first(sex),
    pain_sev = first(pain_severity_std),
    .groups = "drop"
  ) %>%
  mutate(pain_sev_cat = factor(pain_sev, levels = 0:3,
                                labels = c("None","Mild","Moderate","Severe")))

cat_cox <- function(outcome_var, base_var) {
  dd <- surv %>% filter(!is.na(pain_sev_cat), !is.na(baseline_age), !is.na(baseline_sex),
                        .data[[base_var]] == 0) %>%
    as.data.frame()
  if (nrow(dd) < 100) return(NULL)
  f <- as.formula(sprintf("Surv(time, %s) ~ pain_sev_cat + baseline_age + baseline_sex + strata(cohort)",
                          outcome_var))
  fit <- coxph(f, data = dd)
  s <- summary(fit)
  tibble(
    category = c("None (ref)", "Mild", "Moderate", "Severe"),
    HR = c(1, s$coefficients[grep("pain_sev_cat", rownames(s$coefficients)), "exp(coef)"]),
    lci = c(1, s$conf.int[grep("pain_sev_cat", rownames(s$conf.int)), "lower .95"]),
    uci = c(1, s$conf.int[grep("pain_sev_cat", rownames(s$conf.int)), "upper .95"])
  )
}

dr_adl <- cat_cox("adl_event", "adl_base") %>% mutate(outcome = "ADL disability")
dr_iadl <- cat_cox("iadl_event", "iadl_base") %>% mutate(outcome = "IADL disability")
dr <- bind_rows(dr_adl, dr_iadl) %>%
  mutate(category = factor(category, levels = c("None (ref)","Mild","Moderate","Severe")))

write_csv(dr, file.path(TAB_DIR, "dose_response_categorical.csv"))

p2 <- ggplot(dr, aes(x = category, y = HR, color = outcome, group = outcome)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  geom_point(size = 4) +
  geom_line(linewidth = 1) +
  geom_errorbar(aes(ymin = lci, ymax = uci), width = 0.15, linewidth = 0.9) +
  geom_text(aes(label = sprintf("%.2f\n(%.2f-%.2f)", HR, lci, uci)),
            vjust = -0.4, size = 3, color = "grey20") +
  scale_y_log10(breaks = c(1, 1.5, 2, 3, 5), limits = c(0.8, 5)) +
  scale_color_manual(values = c("ADL disability" = "#C3423F", "IADL disability" = "#2E86AB")) +
  labs(
    title = "Figure 2. Dose-response of pain severity on incident disability",
    subtitle = "One-stage stratified Cox model, age- and sex-adjusted",
    x = "Baseline pain severity", y = "Hazard Ratio (log scale)",
    color = "Outcome"
  ) + theme_pub
ggsave(file.path(FIG_DIR, "fig2_dose_response.png"), p2, width = 8, height = 6, dpi = 300)
log_step("Fig 2 dose-response done", "15_fig")

# ---- Fig 3: Cohort-level baseline prevalence -------------------------------
cohort_prev <- baseline %>%
  group_by(cohort) %>%
  summarise(
    n = n(),
    pct_pain = 100 * mean(pain_any_derived, na.rm = TRUE),
    pct_sev_mod_plus = 100 * mean(pain_severity_std >= 2, na.rm = TRUE),
    pct_adl = 100 * mean(adl_disability, na.rm = TRUE),
    pct_iadl = 100 * mean(iadl_disability, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(c(pct_pain, pct_sev_mod_plus, pct_adl, pct_iadl),
               names_to = "indicator", values_to = "pct") %>%
  mutate(indicator = recode(indicator,
                            pct_pain = "Any pain",
                            pct_sev_mod_plus = "Moderate/Severe",
                            pct_adl = "ADL disability",
                            pct_iadl = "IADL disability"))

p3 <- ggplot(cohort_prev, aes(x = cohort, y = pct, fill = indicator)) +
  geom_col(position = position_dodge(0.8), width = 0.75) +
  geom_text(aes(label = sprintf("%.1f", pct)),
            position = position_dodge(0.8), vjust = -0.4, size = 2.8) +
  scale_fill_brewer(palette = "Set1") +
  labs(
    title = "Figure 3. Baseline prevalence of pain and disability across cohorts",
    x = NULL, y = "Prevalence (%)", fill = NULL
  ) + theme_pub + theme(axis.text.x = element_text(angle = 0))
ggsave(file.path(FIG_DIR, "fig3_prevalence.png"), p3, width = 10, height = 5.5, dpi = 300)
log_step("Fig 3 prevalence done", "15_fig")

# ---- Fig 4 (appendix): Kaplan-Meier cumulative ADL by severity -------------
km_data <- surv %>%
  filter(!is.na(pain_sev_cat), adl_base == 0, time > 0) %>%
  as.data.frame()
fit_km <- survfit(Surv(time, adl_event) ~ pain_sev_cat, data = km_data)

km_df <- broom::tidy(fit_km) %>%
  mutate(strata = str_remove(strata, "pain_sev_cat="),
         cum_event = 1 - estimate,
         cum_lci = 1 - conf.high,
         cum_uci = 1 - conf.low)

p4 <- ggplot(km_df, aes(x = time, y = cum_event, color = strata, fill = strata)) +
  geom_step(linewidth = 1) +
  geom_ribbon(aes(ymin = cum_lci, ymax = cum_uci), alpha = 0.15, color = NA) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, NA)) +
  scale_color_manual(values = c("None"="#048A81", "Mild"="#F5B841",
                                 "Moderate"="#E67E22", "Severe"="#C3423F")) +
  scale_fill_manual(values = c("None"="#048A81", "Mild"="#F5B841",
                                "Moderate"="#E67E22", "Severe"="#C3423F")) +
  labs(
    title = "Appendix Fig. Cumulative incidence of ADL disability by pain severity",
    x = "Follow-up (survey waves)", y = "Cumulative incidence",
    color = "Baseline pain", fill = "Baseline pain"
  ) + theme_pub
ggsave(file.path(FIG_DIR, "fig_appendix_km_adl.png"), p4, width = 8, height = 6, dpi = 300)
log_step("KM done", "15_fig")

cat("\n===== Main figures done =====\n")
cat(sprintf("Output dir: %s\n", FIG_DIR))
list.files(FIG_DIR, pattern = "^fig")
