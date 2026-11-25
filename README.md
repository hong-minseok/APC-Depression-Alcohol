# Generational weakening of the association between depressive symptoms and heavy drinking among adults in Korea: an age–period–cohort analysis, 2014–2024

This repository provides the analysis code for a submitted manuscript.

```text
01_preprocessing.R          Variable definitions and data preparation
02_main_analysis.R          APC model fitting, model comparison, and age and cohort effects by sex
03_sensitivity_analysis.R   Outcome definition, period trend assumption, one survey omitted, and covariate adjustment
```

Requires R with `haven`, `survey`, `tidyselect`, and `MASS` (an R recommended package).
Verified with R 4.6.1, haven 2.5.5, survey 4.5, MASS 7.3-65, and tidyselect 1.2.1.
From this directory, run:

```sh
Rscript 01_preprocessing.R /path/to/KNHANES results/data.rds
Rscript 02_main_analysis.R results/data.rds results
Rscript 03_sensitivity_analysis.R results/data.rds results
```

Place the KNHANES files for even years from 2014 to 2024 in year-specific folders (e.g., `2014/hn14_all.sas7bdat`). Survey data are distributed separately.

The scripts save two results files in `results/`:

- `main_analysis.rds` contains model comparisons and estimates of age, period, and cohort effects. Its `sex_comparison` element gives the estimates for men and women and their men-to-women ratios from one interaction model (Table S2).
- `sensitivity_analysis.rds` contains the high-risk drinking model comparisons (Table S1) and the four sensitivity analyses (Table S3).

MIT License. Copyright © 2026 Hong Minseok.
