library(survey)
library(splines)
options(survey.lonely.psu = "adjust", survey.adjust.domain.lonely = TRUE)

curvature_basis <- function(x, knots, boundaries) {
  spline <- function(x) cbind(1, ns(x, knots = knots, Boundary.knots = boundaries))
  grid <- seq(boundaries[1], boundaries[2])
  linear_trend <- cbind(1, grid - mean(grid))
  spline(x) %*% MASS::Null(crossprod(spline(grid), linear_trend))
}

primary_cohort_basis <- function(cohort) curvature_basis(cohort, 1976, c(1955, 2004))

analysis_design <- function(data, outcome = "weekly_hed") {
  data <- data[data$design_valid, ]
  data$outcome <- data[[outcome]]
  design <- svydesign(ids = ~psu, strata = ~stratum, weights = ~weight, data = data)
  design <- subset(design, age >= 20 & age <= 59 &
                   !is.na(phq9) & !is.na(outcome) & !is.na(sex))
  update(design,
    phq5 = phq9 / 5,
    sex_centered = sex - 1.5,
    cohort = year - age,
    age10 = (age - 40) / 10,
    period10 = (year - 2019) / 10,
    cohort10 = (cohort - 1979) / 10,
    survey_year = factor(year),
    baseline_age = ns(age, knots = c(30, 40, 50), Boundary.knots = c(20, 60)),
    sex_age = ns(age, knots = c(33, 46), Boundary.knots = c(20, 60)),
    cohort_curve = primary_cohort_basis(cohort))
}

baseline_formula <- outcome ~ baseline_age * survey_year + sex_centered +
  sex_centered:(sex_age + survey_year) + phq5 + phq5:sex_centered
primary_terms <- c("phq5:age10", "phq5:period10", "phq5:cohort_curve")
effect_contrasts <- list(age = c("phq5:age10" = 1, "phq5:period10" = 1),
                         cohort = c("phq5:period10" = 1))
birth_years <- c(1960, 1980, 2000)
cohort_contrasts_at_age_40 <- lapply(setNames(birth_years, birth_years), function(year) {
  c("phq5" = 1, "phq5:period10" = (year - 1979) / 10,
    "phq5:cohort_curve" = drop(primary_cohort_basis(year)))
})

fit_model <- function(design, time_terms, covariates = character()) {
  formula <- update(baseline_formula, reformulate(c(".", time_terms, covariates), response = "."))
  svyglm(formula, design = design, family = quasibinomial(),
         control = glm.control(epsilon = 1e-10, maxit = 60))
}

model_estimates <- function(fit) {
  list(coefficients = coef(fit), covariance = vcov(fit),
       residual_df = fit$df.residual, n = nrow(fit$model))
}

ratio_estimates <- function(fit, contrasts) {
  estimates <- svycontrast(fit, contrasts)
  exp(cbind(ratio = coef(estimates), confint(estimates, df = fit$df.residual)))
}

wald_p_value <- function(fit, contrasts) {
  estimates <- svycontrast(fit, contrasts)
  beta <- coef(estimates)
  covariance <- vcov(estimates)
  df <- length(beta)
  statistic <- drop(crossprod(beta, solve(covariance, beta))) / df
  pf(statistic, df, fit$df.residual, lower.tail = FALSE)
}

select_models <- function(design) {
  axes <- c(age = "age", period = "year", cohort = "cohort")
  boundaries <- list(age = c(20, 60), period = c(2014, 2024), cohort = c(1955, 2004))
  for (axis in names(axes)) {
    values <- design$variables[[axes[[axis]]]]
    for (number in 1:3) {
      knots <- round(quantile(values, seq_len(number) / (number + 1)))
      design$variables[[paste0(axis, number)]] <- curvature_basis(values, knots, boundaries[[axis]])
    }
  }
  grid <- expand.grid(age = 0:3, period = 0:3, cohort = 0:3)
  fits <- apply(grid, 1, function(knot_counts) {
    curves <- paste0("phq5:", names(knot_counts), knot_counts)[knot_counts > 0]
    fit_model(design, c("phq5:age10", "phq5:period10", curves))
  })
  grid$time_terms <- 2 + rowSums(grid)
  grid$AIC <- sapply(fits, AIC)["AIC", ]
  grid$delta_AIC <- grid$AIC - min(grid$AIC)
  eligible <- which(grid$delta_AIC <= 2)
  selected <- eligible[order(grid$time_terms[eligible], grid$AIC[eligible])][1]
  lowest <- which.min(grid$AIC)
  list(comparison = grid, selected = grid[selected, ], lowest_AIC = grid[lowest, ],
       selected_model = model_estimates(fits[[selected]]),
       lowest_AIC_model = model_estimates(fits[[lowest]]))
}

compare_time_models <- function(design, covariates = character()) {
  model_terms <- list(None = character(), Age = "phq5:age10", Period = "phq5:period10",
                      Cohort = c("phq5:cohort10", "phq5:cohort_curve"), APC = primary_terms)
  fits <- lapply(model_terms, fit_model, design = design, covariates = covariates)
  comparisons <- list(None_Age = c("None", "Age"), None_Period = c("None", "Period"),
                      None_Cohort = c("None", "Cohort"), None_APC = c("None", "APC"),
                      Age_APC = c("Age", "APC"), Period_APC = c("Period", "APC"),
                      Cohort_APC = c("Cohort", "APC"))
  tests <- lapply(comparisons, function(pair) {
    test <- anova(fits[[pair[1]]], fits[[pair[2]]], method = "Wald", test = "F")
    c(F = test$Ftest, df1 = test$df, df2 = test$ddf, p = test$p)
  })
  aic <- sapply(fits, AIC)["AIC", ]
  list(models = lapply(fits, model_estimates),
       effects = ratio_estimates(fits$APC, effect_contrasts),
       AIC = aic, delta_AIC = aic - min(aic), tests = do.call(rbind, tests))
}

period_test <- function(design) {
  waves <- sort(unique(design$variables$year))
  polynomial <- poly(waves, length(waves) - 1)
  design <- update(design, period_deviations = predict(polynomial, year)[, -1, drop = FALSE])
  fit <- fit_model(design, c(primary_terms, "phq5:period_deviations"))
  regTermTest(fit, ~phq5:period_deviations, method = "Wald")
}

cohort_curvature <- function(design) {
  fit <- fit_model(design, primary_terms)
  decade_slope <- function(years) {
    basis <- primary_cohort_basis(years)[, 1]
    10 * coef(lm(basis ~ years))[["years"]]
  }
  slope_change <- decade_slope(1980:2000) - decade_slope(1960:1979)
  contrast <- list(slope_change = c("phq5:cohort_curve" = slope_change))
  list(test = regTermTest(fit, ~phq5:cohort_curve, method = "Wald"),
       slope_change = ratio_estimates(fit, contrast))
}

sex_comparison <- function(design) {
  fit <- fit_model(design, c(primary_terms, paste0("sex_centered:", primary_terms)))
  contrasts <- c(effect_contrasts, cohort_contrasts_at_age_40,
                 list(curvature = c("phq5:cohort_curve" = 1)))
  at_sex <- function(contrast, sex_centered) {
    interaction <- setNames(sex_centered * contrast, paste0("sex_centered:", names(contrast)))
    c(contrast, interaction)
  }
  men <- lapply(contrasts, at_sex, sex_centered = -0.5)
  women <- lapply(contrasts, at_sex, sex_centered = 0.5)
  men_to_women <- Map(`-`, men, women)

  summarize <- function(row, n = NA_integer_) {
    list(n = n,
         age = ratio_estimates(fit, row["age"]),
         cohort = ratio_estimates(fit, row[names(cohort_contrasts_at_age_40)]),
         p_cohort_flat = wald_p_value(fit, row[c("cohort", "curvature")]),
         p_cohort_straight = wald_p_value(fit, row["curvature"]))
  }
  list(men = summarize(men, sum(fit$model$sex_centered == -0.5)),
       women = summarize(women, sum(fit$model$sex_centered == 0.5)),
       men_to_women = summarize(men_to_women))
}

if (sys.nframe() == 0L) {
  arguments <- commandArgs(trailingOnly = TRUE)
  data <- readRDS(arguments[1])
  design <- analysis_design(data)
  results <- list(
    model_selection = select_models(design),
    time_models = compare_time_models(design),
    period = period_test(design),
    cohort_curvature = cohort_curvature(design),
    sex_comparison = sex_comparison(design))
  dir.create(arguments[2], recursive = TRUE, showWarnings = FALSE)
  saveRDS(results, file.path(arguments[2], "main_analysis.rds"))
}
