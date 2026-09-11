## contents

| File | Description |
|---|---|
| `first_phase_recoding.Rmd` | Data preparation: imports the 2023 BRFSS public-use file, applies the inclusion and exclusion criteria, and recodes study variables, with the reasoning for each coding decision documented in the text. |
| `Analytic_code_Revised.R` | Main analysis: descriptive statistics, Welch's ANOVA and Tukey comparisons, the nested linear regression models (Models 0–5) with HC3 robust standard errors, logistic regression, and sensitivity analyses (negative binomial, survey-weighted, psychosocial subsets, Benjamini–Hochberg correction). Produces the manuscript tables, figures, and supplementary appendices. |

## How to reproduce

1. Download the 2023 BRFSS data (SAS Transport format) from
   https://www.cdc.gov/brfss/annual_data/annual_2023.html
   (analyses used the file downloaded in 2024; CDC may has since modified the posted file).
2. Run `first_phase_recoding.Rmd` to create the recoded analytic dataset.
3. Run `Analytic_code_Revised.R`.

Analyses were run in R 4.2.2. BRFSS data are not included in this repository.
