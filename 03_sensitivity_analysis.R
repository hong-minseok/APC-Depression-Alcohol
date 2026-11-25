source("02_main_analysis.R")

cohort_curve_summary <- function(fit, period_log_trend = 0) {
  decades <- (birth_years - 1979) / 10
  flat <- regTermTest(fit, ~phq5:period10 + phq5:cohort_curve,
                      null = c(period_log_trend, 0), method = "Wald")
  straight <- regTermTest(fit, ~phq5:cohort_curve, method = "Wald")
  list(n = nrow(fit$model),
       age = ratio_estimates(fit, effect_contrasts["age"]) * exp(-period_log_trend),
       cohort = ratio_estimates(fit, cohort_contrasts_at_age_40) * exp(-period_log_trend * decades),
       p_cohort_flat = flat$p, p_cohort_straight = straight$p)
}

arguments <- commandArgs(trailingOnly = TRUE)
data <- readRDS(arguments[1])
design <- analysis_design(data)
results <- list()
primary <- fit_model(design, primary_terms)
results$primary <- cohort_curve_summary(primary)

high_risk_design <- analysis_design(data, "high_risk_drinking")
results$outcome_definition <- cohort_curve_summary(fit_model(high_risk_design, primary_terms))

period_only <- fit_model(design, "phq5:period10")
period_log_trend <- coef(period_only)[["phq5:period10"]]
results$period_trend_assumption <- cohort_curve_summary(primary, period_log_trend)
results$period_trend_assumption$period <- ratio_estimates(
  period_only, list(period = c("phq5:period10" = 1)))

waves <- sort(unique(data$year))
results$one_survey_omitted <- setNames(lapply(waves, function(wave) {
  cohort_curve_summary(fit_model(subset(design, year != wave), primary_terms))
}), waves)

socioeconomic <- c("education", "income_quartile", "marital_status", "employment")
adjustments <- list(complete_cases = character(), socioeconomic = socioeconomic,
                    bmi_activity = c(socioeconomic, "bmi", "physical_activity"))
complete_design <- subset(design, complete.cases(design$variables[adjustments$bmi_activity]))
results$covariate_adjustment <- lapply(adjustments, function(covariates) {
  cohort_curve_summary(fit_model(complete_design, primary_terms, covariates))
})

results$high_risk_model_selection <- select_models(high_risk_design)$comparison
dir.create(arguments[2], recursive = TRUE, showWarnings = FALSE)
saveRDS(results, file.path(arguments[2], "sensitivity_analysis.rds"))
