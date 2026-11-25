library(haven)

valid_response <- function(x, values) {
  ifelse(x %in% values, as.numeric(x), NA_real_)
}

read_wave <- function(year, raw_directory) {
  file <- file.path(raw_directory, year, sprintf("hn%02d_all.sas7bdat", year %% 100))
  fields <- c("age", "sex", "wt_itvex", "kstrata", "psu",
              paste0("BP_PHQ_", 1:9), "BD1", "BD1_11", "BD2_1", "BD2_31", "BD2_32",
              "edu", "ho_incm", "marri_1", "marri_2", "EC1_1", "HE_BMI", "pa_aerobic")
  raw <- read_sas(file, col_select = tidyselect::all_of(fields), encoding = "CP949")
  sex <- valid_response(raw$sex, 1:2)
  phq_items <- sapply(raw[paste0("BP_PHQ_", 1:9)], valid_response, values = 0:3)
  lifetime_drinking <- valid_response(raw$BD1, 1:2)
  drinking_frequency <- valid_response(raw$BD1_11, 1:6)
  drinking_quantity <- valid_response(raw$BD2_1, 1:5)
  heavy_drinking_frequency <- valid_response(ifelse(sex == 1, raw$BD2_31, raw$BD2_32), 1:5)
  abstainer <- lifetime_drinking %in% 1 | drinking_frequency %in% 1
  high_quantity <- ifelse(sex == 1, drinking_quantity >= 4, drinking_quantity >= 3)
  high_risk_drinking <- as.integer(high_quantity & drinking_frequency >= 5)
  weekly_hed <- as.integer(heavy_drinking_frequency >= 4)
  inconsistent_skip <- (lifetime_drinking %in% 1 & drinking_frequency %in% 2:6) |
    (abstainer & (!is.na(drinking_quantity) | heavy_drinking_frequency %in% 2:5))
  high_risk_drinking[abstainer] <- weekly_hed[abstainer] <- 0L
  high_risk_drinking[inconsistent_skip] <- weekly_hed[inconsistent_skip] <- NA_integer_

  data.frame(
    year, age = as.numeric(raw$age), sex,
    phq9 = rowSums(phq_items), weekly_hed, high_risk_drinking,
    weight = as.numeric(raw$wt_itvex),
    stratum = paste(year, raw$kstrata, sep = ":"),
    psu = paste(year, raw$psu, sep = ":"),
    design_valid = is.finite(raw$wt_itvex) & raw$wt_itvex > 0 &
      !is.na(raw$kstrata) & nzchar(as.character(raw$kstrata)) &
      !is.na(raw$psu) & nzchar(as.character(raw$psu)),
    education = factor(raw$edu, levels = 4:1,
      labels = c("college", "high_school", "middle_or_less", "middle_or_less")),
    income_quartile = factor(raw$ho_incm, levels = 4:1),
    marital_status = factor(ifelse(raw$marri_1 %in% 2, "never_married",
      ifelse(raw$marri_2 %in% 1, "married",
        ifelse(raw$marri_2 %in% 2:4, "separated_widowed_divorced", NA))),
      levels = c("married", "never_married", "separated_widowed_divorced")),
    employment = factor(raw$EC1_1, levels = 1:2, labels = c("yes", "no")),
    bmi = relevel(cut(raw$HE_BMI, c(0, 18.5, 23, 25, Inf), right = FALSE,
      labels = c("underweight", "normal", "overweight", "obese")), ref = "normal"),
    physical_activity = factor(raw$pa_aerobic, levels = 0:1, labels = c("no", "yes"))
  )
}

arguments <- commandArgs(trailingOnly = TRUE)
data <- do.call(rbind, lapply(seq(2014, 2024, 2), read_wave, raw_directory = arguments[1]))
dir.create(dirname(arguments[2]), recursive = TRUE, showWarnings = FALSE)
saveRDS(data, arguments[2])
