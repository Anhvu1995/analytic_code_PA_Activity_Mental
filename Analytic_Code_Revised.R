# =============================================================
# Final analytic (2023 BRFSS only)
# Changes implemented per reviewer and self-review:
#   1. "Intensity" language removed; framing is PA type, frequency, duration
#   2. Past tense throughout comments and narrative guidance
#   3. p-values formatted as p < 0.0001 / 0.001 / 0.01 / 0.05
#   4. Effect sizes: eta^2, epsilon^2, Cohen's d, standardised beta, OR, Nagelkerke R^2
#   5. Table 3 -> TukeyHSD pairwise results only (ANOVA omnibus in text)
#   6. Calisthenics removed from conclusions (p = 0.078, non-significant)
#   7. LRT Model 3 vs 2 explicitly flagged: PA modality adds no fit (p = 0.341)
#   8. Logistic regression (>=14 days) elevated to main Results section
#   9. Appendix A: Accuracy/Precision/Sensitivity columns removed;
#      RMSE, R^2, Adj-R^2, MAE, AIC, BIC retained
#  10. HC3 robust SEs reported as primary inference (BP confirmed heteroscedasticity)
#  11. Food stamps finding flagged as non-robust to weighting
#  12. Weighted model divergences (Hispanic, job loss) noted for Discussion
#  13. Figure 3 legend fixed with explicit fill labels
#  14. EXERANY2 yes/no preliminary analysis removed (no variation in analytic sample)
#  15. Gardening/yard work retained in model but flagged: do not interpret directionally
# -------------------------------------------------------------
#  SECTION H additions (this revision) - reviewer round 2:
#  H1. Attrition table: included vs. excluded, reproduced as a step-by-step
#      WATERFALL matching the actual Section 3 filter order            -> Supp. Table S5
#  H2. MENTHLTH histogram/density plot                                  -> Supp. Figure S3
#  H3. Negative binomial regression (sensitivity check)                 -> Supp. Table S6
#  H4. Psychosocial-predictor subset sensitivity analysis                -> Supp. Table S7
#  H5. Combined unweighted-vs-weighted side-by-side table                -> Supp. Table S8

# ---- 0. Packages ----
packages <- c(
  "haven", "dplyr", "tidyr", "ggplot2", "car", "lmtest",
  "sandwich", "nortest", "ggfortify", "broom", "stringr",
  "lm.beta", "effectsize", "rstatix", "FSA", "flextable",
  "tableone", "survey", "labelled", "gridExtra", "scales",
  "DescTools", "tibble", "forcats", "patchwork"
)
lapply(packages, library, character.only = TRUE)


# ---- 1. Load raw data ----
BRFSS_2023 <- read_xpt("LLCP2023.XPT")


# ---- 2. Select & rename ----
BRFSS_2023.small <- BRFSS_2023 %>%
  select(
    "_AGE80", "_IMPRACE", "_INCOMG1", "_LLCPWT", "_STATE",
    "LANDSEX2", "CELLSEX2",
    "EMPLOY1", "MENTHLTH", "EXERANY2", "ADDEPEV3",
    "EXRACT12", "PAFREQ1_", "PADUR1_", "_MINAC12", "STRFREQ_",
    "LSATISFY", "EMTSUPRT", "SDLONELY", "SDHSTRE1",
    "SDHEMPLY", "FOODSTMP", "SDHBILLS", "SDHUTILS", "SDHTRNSP"
  ) %>%
  mutate(SEX = coalesce(LANDSEX2, CELLSEX2)) %>%
  rename(
    AGE80   = `_AGE80`,
    INCOMG1 = `_INCOMG1`,
    IMPRACE = `_IMPRACE`,
    MINAC12 = `_MINAC12`,
    LLCPWT  = `_LLCPWT`,
    STATE   = `_STATE`
  ) %>%
  filter(EMPLOY1 == 6) %>%               # student status only
  mutate_at(
    vars(-MENTHLTH, -AGE80, -PAFREQ1_, -PADUR1_, -MINAC12, -STRFREQ_, -LLCPWT),
    as.factor
  ) %>%
  select(-LANDSEX2, -CELLSEX2, -EMPLOY1) %>%
  select(SEX, AGE80, IMPRACE, LLCPWT, everything())

cat("Step 2 check - Student participants (EMPLOY1 == 6):", nrow(BRFSS_2023.small), "\n")

# ---- 3. Recode ----
BRFSS_2023.recoded <- BRFSS_2023.small %>%
  
  # Sex
  filter(!SEX %in% c("7", "9")) %>%
  mutate(SEX = recode_factor(SEX,
                             `1` = "Male",
                             `2` = "Female",
                             `3` = "Unspecified or another gender identity")) %>%
  
  # Age: restrict to 18-24
  filter(AGE80 >= 18, AGE80 <= 24) %>%
  
  # Ethnicity
  mutate(IMPRACE = recode_factor(IMPRACE,
                                 `1` = "White, Non-Hispanic",
                                 `2` = "Black, Non-Hispanic",
                                 `3` = "Asian, Non-Hispanic",
                                 `4` = "American Indian/Alaskan Native, Non-Hispanic",
                                 `5` = "Hispanic",
                                 `6` = "Other race, Non-Hispanic")) %>%
  
  # Income
  filter(INCOMG1 != "9") %>%
  mutate(INCOMG1 = recode_factor(INCOMG1,
                                 `1` = "Less than $15,000",
                                 `2` = "$15,000 - $25,000",
                                 `3` = "$25,000 - $35,000",
                                 `4` = "$35,000 - $50,000",
                                 `5` = "$50,000 - $100,000",
                                 `6` = "$100,000 - $200,000",
                                 `7` = "$200,000 or more")) %>%
  
  # Outcome: continuous (recode 88 -> 0, drop 77/99)
  mutate(MENTHLTH = if_else(MENTHLTH == 88, 0, MENTHLTH)) %>%
  filter(!MENTHLTH %in% c(77, 99)) %>%
  
  # Outcome: binary - CDC "frequent mental distress" threshold (>=14 days)
  # This binary outcome is used in the supplementary logistic model to
  # produce odds ratios and address reviewer concern about modest R^2.
  mutate(MENTHLTH_BIN = as.integer(MENTHLTH >= 14)) %>%
  
  # Exercise (any): used only to confirm all analytic-sample respondents
  # exercised in the past month. Variable dropped after filtering.
  filter(EXERANY2 == "1") %>%
  select(-EXERANY2) %>%
  
  # Exercise type
  filter(!EXRACT12 %in% c("77", "99")) %>%
  mutate(EXRACT12 = recode_factor(EXRACT12,
                                  `1`  = "Walking",
                                  `2`  = "Running or jogging",
                                  `3`  = "Gardening or yard work",
                                  `4`  = "Bicycling",
                                  `5`  = "Aerobics",
                                  `6`  = "Calisthenics",
                                  `7`  = "Elliptical/EFX",
                                  `8`  = "Household activities",
                                  `9`  = "Weight lifting",
                                  `10` = "Yoga/Pilates/Tai Chi",
                                  `11` = "Other")) %>%
  
  # PA frequency and strength (stored as integer x 1000)
  filter(!PAFREQ1_ %in% 99000) %>%
  mutate(PAFREQ1_ = PAFREQ1_ / 1000) %>%
  filter(!STRFREQ_  %in% 99000) %>%
  mutate(STRFREQ_  = STRFREQ_  / 1000) %>%
  
  # Depression diagnosis
  filter(!ADDEPEV3 %in% c("7", "9")) %>%
  mutate(ADDEPEV3 = recode_factor(ADDEPEV3, `1` = "Yes", `2` = "No")) %>%
  
  # Life satisfaction
  filter(!LSATISFY %in% c("7", "9")) %>%
  mutate(LSATISFY = recode_factor(LSATISFY,
                                  `1` = "Very satisfied", `2` = "Satisfied",
                                  `3` = "Dissatisfied",   `4` = "Very dissatisfied")) %>%
  
  # Emotional support
  filter(!EMTSUPRT %in% c("7", "9")) %>%
  mutate(EMTSUPRT = recode_factor(EMTSUPRT,
                                  `1` = "Always", `2` = "Usually", `3` = "Sometimes",
                                  `4` = "Rarely",  `5` = "Never")) %>%
  
  # Loneliness
  filter(!SDLONELY %in% c("7", "9")) %>%
  mutate(SDLONELY = recode_factor(SDLONELY,
                                  `1` = "Always", `2` = "Usually", `3` = "Sometimes",
                                  `4` = "Rarely",  `5` = "Never")) %>%
  
  # Stress
  filter(!SDHSTRE1 %in% c("7", "9")) %>%
  mutate(SDHSTRE1 = recode_factor(SDHSTRE1,
                                  `1` = "Always", `2` = "Usually", `3` = "Sometimes",
                                  `4` = "Rarely",  `5` = "Never")) %>%
  
  # SES indicators
  filter(!SDHEMPLY %in% c("7", "9")) %>%
  mutate(SDHEMPLY = recode_factor(SDHEMPLY, `1` = "Yes", `2` = "No")) %>%
  filter(!FOODSTMP %in% c("7", "9")) %>%
  mutate(FOODSTMP = recode_factor(FOODSTMP, `1` = "Yes", `2` = "No")) %>%
  filter(!SDHBILLS %in% c("7", "9")) %>%
  mutate(SDHBILLS = recode_factor(SDHBILLS, `1` = "Yes", `2` = "No")) %>%
  filter(!SDHUTILS %in% c("7", "9")) %>%
  mutate(SDHUTILS = recode_factor(SDHUTILS, `1` = "Yes", `2` = "No")) %>%
  filter(!SDHTRNSP %in% c("7", "9")) %>%
  mutate(SDHTRNSP = recode_factor(SDHTRNSP, `1` = "Yes", `2` = "No")) %>%
  
  # Drop any remaining NAs (listwise deletion)
  na.omit() %>%
  
  # Exclude PA type groups with insufficient sample size for reliable inference:
  #   Household activities (n = 3), Aerobics (n = 13), Elliptical/EFX (n = 14).
  # The household activities coefficient (beta = 9.34, p = 0.01) in earlier --> runs was driven by only 3 observations, too highly sensitive to outliers.
  filter(!EXRACT12 %in% c("Household activities", "Aerobics", "Elliptical/EFX")) %>%
  mutate(EXRACT12 = droplevels(EXRACT12))

cat("Final analytic sample (after small-group exclusion):", nrow(BRFSS_2023.recoded), "\n")
# Should print 1,862 to match the manuscript and flow diagram.

# ---- 3b. Set reference categories (must run before any model is fit) ----
BRFSS_2023.recoded <- BRFSS_2023.recoded %>%
  mutate(
    SEX      = relevel(factor(SEX),      ref = "Male"),
    IMPRACE  = relevel(factor(IMPRACE),  ref = "White, Non-Hispanic"),
    INCOMG1  = relevel(factor(INCOMG1),  ref = "Less than $15,000"),
    EXRACT12 = relevel(factor(EXRACT12), ref = "Walking"),
    ADDEPEV3 = relevel(factor(ADDEPEV3), ref = "Yes"),
    LSATISFY = relevel(factor(LSATISFY), ref = "Very satisfied"),
    EMTSUPRT = relevel(factor(EMTSUPRT), ref = "Always"),
    SDLONELY = relevel(factor(SDLONELY), ref = "Always"),
    SDHSTRE1 = relevel(factor(SDHSTRE1), ref = "Always"),
    SDHEMPLY = relevel(factor(SDHEMPLY), ref = "Yes"),
    FOODSTMP = relevel(factor(FOODSTMP), ref = "Yes"),
    SDHBILLS = relevel(factor(SDHBILLS), ref = "Yes"),
    SDHUTILS = relevel(factor(SDHUTILS), ref = "Yes"),
    SDHTRNSP = relevel(factor(SDHTRNSP), ref = "Yes")
  )

# Save cleaned/recoded file
write.csv(BRFSS_2023.recoded, "BRFSS_2023.recoded.csv", row.names = FALSE)

# ---- 4. Variable labels ----
var_label(BRFSS_2023.recoded) <- list(
  SEX      = "Sex Category",
  AGE80    = "Age category",
  IMPRACE  = "Ethnicity",
  INCOMG1  = "Income Categories",
  MENTHLTH = "Days Mental Health Not Good (Past 30 Days)",
  MENTHLTH_BIN = "Frequent Mental Distress (>=14 days)",
  ADDEPEV3 = "(Ever told) to Have a Depressive Disorder",
  EXRACT12 = "Type of Physical Activity",
  PAFREQ1_ = "Physical Activity Frequency per Week",
  PADUR1_  = "Duration of One Session in Minutes",
  MINAC12  = "Total Minutes of Physical Activity per Week",
  STRFREQ_ = "Strength Activity Frequency per Week",
  LSATISFY = "Life Satisfaction",
  EMTSUPRT = "Adequate Social and Emotional Support",
  SDLONELY = "Feeling Loneliness",
  SDHSTRE1 = "Feeling Stress Within the Last 30 Days",
  SDHEMPLY = "Job Loss or Reduced Hours",
  FOODSTMP = "Received Food Stamps (Past 12 Months)",
  SDHBILLS = "Unable to Pay Bills",
  SDHUTILS = "Unable to Pay Utilities or Service Threatened",
  SDHTRNSP = "Lack of Transportation Affected Daily Needs"
)


# ---- 5. P-value formatter ----
fmt_p <- function(p) {
  case_when(
    p < 0.0001 ~ "p < 0.0001",
    p < 0.001  ~ "p < 0.001",
    p < 0.01   ~ "p < 0.01",
    p < 0.05   ~ "p < 0.05",
    TRUE        ~ paste0("p = ", round(p, 3))
  )
}


# ========================================================
# TABLE 1 & 2: Descriptive statistics
# ========================================================
cat_vars  <- c("SEX", "IMPRACE", "INCOMG1", "ADDEPEV3", "EXRACT12", "STATE",
               "LSATISFY", "EMTSUPRT", "SDLONELY", "SDHSTRE1",
               "SDHEMPLY", "FOODSTMP", "SDHBILLS", "SDHUTILS", "SDHTRNSP")
cont_vars <- c("AGE80", "MENTHLTH", "PAFREQ1_", "PADUR1_",
               "MINAC12", "STRFREQ_")

tbl1 <- CreateTableOne(vars = cat_vars,  data = BRFSS_2023.recoded, factorVars = cat_vars)
tbl2 <- CreateTableOne(vars = cont_vars, data = BRFSS_2023.recoded)

print(tbl1, showAllLevels = TRUE, varLabels = TRUE)
print(tbl2, varLabels = TRUE)


# ========================================================
# SECTION A: ANOVA - PA type vs. MENTHLTH
# ========================================================

fig1 <- ggplot(
  BRFSS_2023.recoded,
  aes(x    = fct_reorder(EXRACT12, MENTHLTH, .fun = median),
      y    = MENTHLTH,
      fill = EXRACT12)
) +
  geom_jitter(
    aes(color = EXRACT12),
    width  = 0.25,
    alpha  = 0.40,
    size   = 0.7,
    shape  = 16
  ) +
  geom_boxplot(
    width         = 0.45,
    outlier.shape = NA,
    alpha         = 0.55,
    color         = "grey25",
    linewidth     = 0.45
  ) +
  scale_fill_manual(
    guide  = "none",
    values = c(
      "Walking"              = "#4E79A7",
      "Running or jogging"   = "#F28E2B",
      "Gardening or yard work" = "#59A14F",
      "Bicycling"            = "#E15759",
      "Calisthenics"         = "#B07AA1",
      "Weight lifting"       = "#76B7B2",
      "Yoga/Pilates/Tai Chi" = "#FF9DA7",
      "Other"                = "#9C755F"
    )
  ) +
  scale_color_manual(
    guide  = "none",
    values = c(
      "Walking"              = "#4E79A7",
      "Running or jogging"   = "#F28E2B",
      "Gardening or yard work" = "#59A14F",
      "Bicycling"            = "#E15759",
      "Calisthenics"         = "#B07AA1",
      "Weight lifting"       = "#76B7B2",
      "Yoga/Pilates/Tai Chi" = "#FF9DA7",
      "Other"                = "#9C755F"
    )
  ) +
  scale_y_continuous(breaks = seq(0, 30, by = 5), limits = c(-1, 32)) +
  coord_flip() +
  labs(
    x = "Physical Activity Type",
    y = "Days of Poor Mental Health (past 30 days)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.y        = element_text(size = 9),
    axis.title.x       = element_text(size = 13),
    axis.title.y       = element_text(size = 13),
    plot.subtitle      = element_text(size = 7.5, color = "grey40"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank()
  )

ggsave("Figure1_boxplot_jitter_PA_type.png", fig1,
       width = 11, height = 6, dpi = 300)
print(fig1)


# --- One-way ANOVA ---
anova_result <- aov(MENTHLTH ~ EXRACT12, data = BRFSS_2023.recoded)
anova_summary <- summary(anova_result)
print(anova_summary)

eta2_anova <- eta_squared(anova_result)
print(eta2_anova)


# --- Assumption 1: Homoscedasticity (Levene's test) ---
levene_res <- leveneTest(MENTHLTH ~ EXRACT12, data = BRFSS_2023.recoded)
print(levene_res)

# --- Assumption 2: Independence (Durbin-Watson) ---
dw_res <- durbinWatsonTest(anova_result)
print(dw_res)

# --- Welch's ANOVA (primary, robust to heteroscedasticity) ---
welch_result <- oneway.test(MENTHLTH ~ EXRACT12,
                            data      = BRFSS_2023.recoded,
                            var.equal = FALSE)
print(welch_result)

# --- Kruskal-Wallis (robust to non-normality) ---
kw_result <- kruskal.test(MENTHLTH ~ EXRACT12, data = BRFSS_2023.recoded)
print(kw_result)

epsilon_sq <- as.numeric(kw_result$statistic) / (nrow(BRFSS_2023.recoded) - 1)
cat("Kruskal-Wallis epsilon-squared:", round(epsilon_sq, 3), "\n")


# ---- TABLE 3 (revised): TukeyHSD pairwise results ----
tukey_res <- TukeyHSD(anova_result, conf.level = 0.95)

tukey_df <- as.data.frame(tukey_res$EXRACT12) %>%
  rownames_to_column("Comparison") %>%
  rename(
    Mean_Difference = diff,
    Lower_95CI      = lwr,
    Upper_95CI      = upr,
    p_adj           = `p adj`
  ) %>%
  mutate(
    p_fmt = fmt_p(p_adj),
    sig   = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE          ~ ""
    )
  ) %>%
  filter(p_adj < 0.05) %>%
  arrange(p_adj)

print(tukey_df)

png("Figure2_TukeyHSD.png", width = 1400, height = 1000, res = 150)
par(mar = c(4, 16, 3, 2))
plot(tukey_res, las = 1, cex.axis = 0.6)
dev.off()

sig_pairs <- tukey_df$Comparison

cohens_d_manual <- lapply(sig_pairs, function(pair) {
  lvls  <- levels(BRFSS_2023.recoded$EXRACT12)
  grp_b <- NA_character_
  grp_a <- NA_character_
  for (l in lvls) {
    if (startsWith(pair, l)) {
      grp_b <- l
      grp_a <- substr(pair, nchar(l) + 2, nchar(pair))
      break
    }
  }
  if (is.na(grp_b) || !grp_a %in% lvls) return(NULL)
  
  sub_df <- BRFSS_2023.recoded %>%
    filter(EXRACT12 %in% c(grp_a, grp_b)) %>%
    mutate(EXRACT12 = droplevels(EXRACT12))
  
  counts <- table(sub_df$EXRACT12)
  if (any(counts < 2)) return(NULL)
  
  tryCatch({
    d <- cohens_d(sub_df, MENTHLTH ~ EXRACT12, var.equal = FALSE)
    d$comparison <- pair
    d
  }, error = function(e) NULL)
})
cohens_d_df <- bind_rows(Filter(Negate(is.null), cohens_d_manual))
print(cohens_d_df)


# ========================================================
# SECTION C: Linear regression - nested models
# ========================================================

model_0 <- lm(MENTHLTH ~ SEX + AGE80 + IMPRACE,
              data = BRFSS_2023.recoded)
summary(model_0)

model_1 <- lm(MENTHLTH ~ EXRACT12,
              data = BRFSS_2023.recoded)
summary(model_1)

model_2 <- lm(MENTHLTH ~ SEX + AGE80 + IMPRACE + EXRACT12,
              data = BRFSS_2023.recoded)
summary(model_2)

model_3 <- lm(MENTHLTH ~ SEX + AGE80 + IMPRACE + EXRACT12 +
                PAFREQ1_ + PADUR1_ + MINAC12 + STRFREQ_,
              data = BRFSS_2023.recoded)
summary(model_3)

model_4 <- lm(MENTHLTH ~ SEX + AGE80 + IMPRACE +
                EXRACT12 + PAFREQ1_ + PADUR1_ + MINAC12 + STRFREQ_ +
                ADDEPEV3 + LSATISFY + EMTSUPRT + SDLONELY + SDHSTRE1,
              data = BRFSS_2023.recoded)
summary(model_4)

model_5 <- lm(MENTHLTH ~ SEX + AGE80 + IMPRACE + INCOMG1 +
                EXRACT12 + PAFREQ1_ + PADUR1_ + MINAC12 + STRFREQ_ +
                ADDEPEV3 + LSATISFY + EMTSUPRT + SDLONELY + SDHSTRE1 +
                SDHEMPLY + FOODSTMP + SDHBILLS + SDHUTILS + SDHTRNSP,
              data = BRFSS_2023.recoded)
summary(model_5)


model5_tidy <- tidy(model_5, conf.int = TRUE) %>%
  mutate(p_fmt = fmt_p(p.value)) %>%
  select(term, estimate, std.error, statistic, conf.low, conf.high, p_fmt)

print(model5_tidy)


eta2_m5 <- eta_squared(model_5)
print(eta2_m5)


model_5_std  <- lm.beta(model_5)
std_coefs_m5 <- data.frame(
  term     = names(coef(model_5_std)),
  std_beta = coef(model_5_std)
) %>%
  filter(!is.na(std_beta), term != "(Intercept)") %>%
  arrange(desc(abs(std_beta)))

print(std_coefs_m5)


calculate_metrics <- function(model, data) {
  actual      <- data$MENTHLTH
  predictions <- predict(model, newdata = data)
  rmse        <- sqrt(mean((actual - predictions)^2, na.rm = TRUE))
  smry        <- summary(model)
  mae         <- mean(abs(actual - predictions), na.rm = TRUE)
  list(
    RMSE   = round(rmse, 3),
    R2     = round(smry$r.squared, 3),
    Adj_R2 = round(smry$adj.r.squared, 3),
    MAE    = round(mae, 3),
    AIC    = round(AIC(model), 1),
    BIC    = round(BIC(model), 1)
  )
}

model_list <- list(model_0, model_1, model_2, model_3, model_4, model_5)
metrics_df <- do.call(rbind, lapply(seq_along(model_list), function(i) {
  data.frame(Model = paste0("Model_", i - 1),
             as.data.frame(calculate_metrics(model_list[[i]],
                                             BRFSS_2023.recoded)))
}))
print(metrics_df)


lrt_1v0 <- anova(model_0, model_1, test = "F")
lrt_2v1 <- anova(model_1, model_2, test = "F")
lrt_3v2 <- anova(model_2, model_3, test = "F")
lrt_4v3 <- anova(model_3, model_4, test = "F")
lrt_5v4 <- anova(model_4, model_5, test = "F")

lrt_list <- list("1 vs 0" = lrt_1v0, "2 vs 1" = lrt_2v1,
                 "3 vs 2" = lrt_3v2, "4 vs 3" = lrt_4v3,
                 "5 vs 4" = lrt_5v4)
lapply(names(lrt_list), function(nm) {
  cat("\n--- LRT:", nm, "---\n")
  print(lrt_list[[nm]])
})


# ========================================================
# SECTION D: Logistic regression - binary "frequent distress"
# ========================================================

logit_full <- glm(
  MENTHLTH_BIN ~ SEX + AGE80 + IMPRACE + INCOMG1 +
    EXRACT12 + PAFREQ1_ + PADUR1_ + MINAC12 + STRFREQ_ +
    ADDEPEV3 + LSATISFY + EMTSUPRT + SDLONELY + SDHSTRE1 +
    SDHEMPLY + FOODSTMP + SDHBILLS + SDHUTILS + SDHTRNSP,
  data   = BRFSS_2023.recoded,
  family = binomial(link = "logit")
)
summary(logit_full)

or_table <- exp(cbind(OR = coef(logit_full), confint(logit_full)))
print(round(or_table, 3))

cat("Nagelkerke pseudo-R^2:", round(PseudoR2(logit_full, which = "Nagelkerke"), 3), "\n")

logit_tidy <- tidy(logit_full, conf.int = TRUE, exponentiate = TRUE) %>%
  mutate(p_fmt = fmt_p(p.value)) %>%
  select(term, estimate, conf.low, conf.high, p_fmt) %>%
  rename(OR = estimate, Lower_CI = conf.low, Upper_CI = conf.high)

print(logit_tidy)


# ========================================================
# SECTION E: Figures
# ========================================================

# --- Figure 4: Standardised effect sizes (Model 5) ---
n_m5 <- nobs(model_5)

std_plot_data <- std_coefs_m5 %>%
  filter(grepl(
    "EXRACT12|PAFREQ1_|PADUR1_|MINAC12|STRFREQ_|ADDEPEV3|LSATISFY|EMTSUPRT|SDLONELY|SDHSTRE1|SDHEMPLY|FOODSTMP|SDHBILLS|SDHUTILS|SDHTRNSP|SEX|IMPRACE",term)) %>%
  mutate(
    term_clean = str_replace_all(term, c(
      "ADDEPEV3No"                    = "Depression History - No",
      "LSATISFY"                      = "Life Satisfaction - ",
      "EMTSUPRT"                      = "Emotional Support - ",
      "SDLONELY"                      = "Loneliness - ",
      "SDHSTRE1"                      = "Stress - ",
      "EXRACT12"                      = "PA Type - ",
      "PAFREQ1_"                      = "PA Frequency (per wk)",
      "PADUR1_"                       = "Session Duration (min)",
      "MINAC12"                       = "Primary PA min (per wk)",
      "STRFREQ_"                      = "Strength Freq. (per wk)",
      "SDHEMPLY"                      = "Job Loss - ",
      "FOODSTMP"                      = "Food Stamps - ",
      "SDHBILLS"                      = "Unpaid Bills - ",
      "SDHUTILS"                      = "Utility Issues - ",
      "SDHTRNSP"                      = "Transport Barrier - ",
      "SEXFemale"                     = "Sex - Female",
      "SEXUnspecified or another gender identity" = "Sex - Nonbinary",
      "IMPRACEBlack, Non-Hispanic"    = "Race - Black",
      "IMPRACEAsian, Non-Hispanic"    = "Race - Asian",
      "IMPRACEHispanic"               = "Race - Hispanic",
      "IMPRACEOther race, Non-Hispanic" = "Race - Others",
      "IMPRACEAmerican Indian/Alaskan Native, Non-Hispanic" = "Race - American Indian/Alaskan Native"
    )),
    direction = ifelse(std_beta > 0, "Higher distress", "Lower distress")
  )

# FIX 1: max_abs now computed AFTER std_plot_data exists (was referenced
# before creation in the original script, which throws "object not found").
max_abs <- max(abs(std_plot_data$std_beta), na.rm = TRUE)

fig3 <- ggplot(std_plot_data,
               aes(x = reorder(term_clean, std_beta),
                   y = std_beta, fill = direction)) +
  geom_col(alpha = 0.85) +
  scale_fill_manual(
    values = c("Higher distress" = "#d73027", "Lower distress" = "#4575b4"),
    name   = "Association direction",
    labels = c("Higher distress", "Lower distress")
  ) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed", color = "grey30") +
  scale_y_continuous(limits = c(-max_abs, max_abs) * 1.05) +
  coord_flip() +
  labs(
    x = NULL,
    y = "Standardised beta coefficient"
  ) +
  theme_minimal(base_size = 10) +
  theme(axis.text.y     = element_text(size = 7.5),
        legend.position = "bottom",
        legend.title    = element_text(size = 8),
        legend.text     = element_text(size = 8))

print(fig3)


# --- Figure 5: Odds ratios - logistic model ---
term_labels <- c(
  "AGE80"                         = "Age",
  "INCOMG1"                       = "",
  "ADDEPEV3No"                    = "No",
  "LSATISFY"                      = "",
  "EMTSUPRT"                      = "",
  "SDLONELY"                      = "",
  "SDHSTRE1"                      = "",
  "EXRACT12"                      = "",
  "PAFREQ1_"                      = "PA Frequency (per wk)",
  "PADUR1_"                       = "Session Duration (min)",
  "MINAC12"                       = "Primary PA min (per wk)",
  "STRFREQ_"                      = "Strength Freq. (per wk)",
  "SDHEMPLY"                      = "Job Loss - ",
  "FOODSTMP"                      = "Food Stamps - ",
  "SDHBILLS"                      = "Unpaid Bills - ",
  "SDHUTILS"                      = "Utility Issues - ",
  "SDHTRNSP"                      = "Transport Barrier - ",
  "SEXFemale"                     = "Female",
  "SEXUnspecified or another gender identity" = "Nonbinary",
  "IMPRACEBlack, Non-Hispanic"    = "Black",
  "IMPRACEAsian, Non-Hispanic"    = "Asian",
  "IMPRACEHispanic"               = "Hispanic",
  "IMPRACEOther race, Non-Hispanic" = "Others",
  "IMPRACEAmerican Indian/Alaskan Native, Non-Hispanic" = "American Indian/Alaskan Native"
)

or_plot_data <- logit_tidy %>%
  filter(term != "(Intercept)") %>%
  mutate(
    term_clean = str_replace_all(term, term_labels),
    
    # Assign each predictor to a group (based on the raw BRFSS variable name)
    group = case_when(
      str_detect(term, "^(SEX)")                     ~ "Sex",
      str_detect(term, "^(AGE80)")                     ~ "Age",
      str_detect(term, "^(IMPRACE)")                     ~ "Race",
      str_detect(term, "^(INCOMG1)")                     ~ "Income",
      str_detect(term, "^(EXRACT12)")    ~ "Physical Activity\nType",
      str_detect(term, "^(PAFREQ1_|PADUR1_|MINAC12|STRFREQ_)")    ~ "Physical\nModality",
      str_detect(term, "^(ADDEPEV3)")  ~ "Depression Diagnosis",
      str_detect(term, "^(LSATISFY)")  ~ "Life Satisfaction",
      str_detect(term, "^(EMTSUPRT)")  ~ "Emotional Support",
      str_detect(term, "^(SDLONELY)")  ~ "Loneliness",
      str_detect(term, "^(SDHSTRE1)")  ~ "Stress",
      str_detect(term, "^(SDHEMPLY|FOODSTMP|SDHBILLS|SDHUTILS|SDHTRNSP)")  ~ "Social\nDeterminants",
      TRUE                                                                 ~ "Other"
    ),
    group = factor(group, levels = c("Sex",
                                     "Age",
                                     "Race",
                                     "Income",
                                     "Physical Activity\nType",
                                     "Physical\nModality",
                                     "Mental Health &\nSocial Support",
                                     "Depression Diagnosis",
                                     "Life Satisfaction",
                                     "Emotional Support",
                                     "Loneliness",
                                     "Stress",
                                     "Social\nDeterminants",
                                     "Other")),
    
    # CI text: 3 decimals for very narrow CIs, otherwise 2
    ci_text = ifelse(Upper_CI - Lower_CI < 0.05,
                     sprintf("%.3f – %.3f", Lower_CI, Upper_CI),
                     sprintf("%.2f – %.2f", Lower_CI, Upper_CI)),
    
    # Order rows by OR (within each group once faceted)
    term_clean = fct_reorder(term_clean, OR)
  ) %>%
  droplevels()

# Left panel: variable names + CI + p-value
table_data <- or_plot_data %>%
  mutate(variable = as.character(term_clean),
         p_fmt    = as.character(p_fmt)) %>%
  select(group, term_clean, variable, ci_text, p_fmt) %>%
  pivot_longer(c(variable, ci_text, p_fmt),
               names_to = "column", values_to = "value") %>%
  mutate(column = factor(column,
                         levels = c("variable", "ci_text", "p_fmt"),
                         labels = c("Variable", "95% CI", "p-value")))

# Right panel: table plot with group strips on the right
table_txt <- 3.2   
header_txt <- 10   

p_table <- ggplot(table_data, aes(x = column, y = term_clean, label = value)) +
  geom_text(data = filter(table_data, column == "Variable"),
            hjust = 0, nudge_x = -0.45, size = table_txt, colour = "grey20") +
  geom_text(data = filter(table_data, column != "Variable"),
            size = table_txt, colour = "grey20") +
  scale_x_discrete(position = "top") +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid      = element_blank(),
    axis.text.y     = element_blank(),
    axis.text.x.top = element_text(face = "bold", size = header_txt),
    strip.text      = element_blank(),
    panel.spacing   = unit(6, "pt")
  )

fig5 <- p_table + p_forest + plot_layout(widths = c(1.5, 1.8))
print(fig5)

# ========================================================
# SECTION F: Model diagnostics (Model 5)
# ========================================================

vif_vals <- vif(model_5)
print(vif_vals)

autoplot(model_5)

bp_test <- bptest(model_5)
print(bp_test)

cat("\n--- HC3 Robust Standard Errors (Model 5) ---\n")
print(coeftest(model_5, vcov = vcovHC(model_5, type = "HC3")))


# ========================================================
# SECTION G: Weighted sensitivity analysis (Appendix C)
# ========================================================
svy_design <- svydesign(ids = ~1, weights = ~LLCPWT, data = BRFSS_2023.recoded)

svy_model <- svyglm(
  MENTHLTH ~ SEX + AGE80 + IMPRACE + INCOMG1 +
    EXRACT12 + PAFREQ1_ + PADUR1_ + MINAC12 + STRFREQ_ +
    ADDEPEV3 + LSATISFY + EMTSUPRT + SDLONELY + SDHSTRE1 +
    SDHEMPLY + FOODSTMP + SDHBILLS + SDHUTILS + SDHTRNSP,
  design = svy_design
)
summary(svy_model)

svy_tidy <- tidy(svy_model, conf.int = TRUE) %>%
  mutate(p_fmt = fmt_p(p.value)) %>%
  select(term, estimate, std.error, conf.low, conf.high, p_fmt)
print(svy_tidy)

comparison <- left_join(
  tidy(model_5)   %>% select(term, unweighted_est = estimate),
  tidy(svy_model) %>% select(term, weighted_est   = estimate),
  by = "term"
) %>%
  mutate(
    direction_match = sign(unweighted_est) == sign(weighted_est),
    abs_diff        = round(abs(unweighted_est - weighted_est), 3)
  )
print(comparison)

# FIX 2: The original script had a second block here titled "Relevel to
# match Appendix B reference categories" that built
# `BRFSS_2023.recoded_svy` and then referenced an undefined
# `svy_model_releveled` two lines later -> guaranteed "object not found"
# error. BRFSS_2023.recoded was ALREADY releveled in Step 3b, before
# svy_model was ever fit, so svy_model/svy_tidy above already use the
# correct reference categories. That redundant block has been removed
# entirely rather than fixed, to avoid computing the same thing twice.


# ========================================================
# SECTION H: Reviewer round-2 sensitivity analyses
# Requires Sections 0-3 (data load/recode) and Section E (model_5)
# to have already been run in this session.
# ========================================================

# ---- H1. Attrition table: WATERFALL matching Section 3's actual filter
#      order, so the printed counts can be pasted directly into the flow
#      diagram / Appendix C1 without hand reconciliation. ----
# Reviewer Concerns #1 & #3: show, step by step, how many respondents were
# dropped at each stage and why, then compare included vs. excluded
# respondents on demographics available for both groups.

waterfall <- tibble(step = character(), n_before = integer(),
                    n_after = integer(), n_dropped = integer())

add_step <- function(df_before, df_after, label) {
  waterfall <<- bind_rows(waterfall, tibble(
    step       = label,
    n_before   = nrow(df_before),
    n_after    = nrow(df_after),
    n_dropped  = nrow(df_before) - nrow(df_after)
  ))
  df_after
}

pool <- BRFSS_2023.small %>% mutate(.row_id = row_number())

step_sex      <- pool %>% filter(!SEX %in% c("7", "9"))
pool          <- add_step(pool, step_sex, "SEX = 'Not sure'/'Refused'")

step_age      <- pool %>% filter(AGE80 >= 18, AGE80 <= 24)
pool          <- add_step(pool, step_age, "Age outside 18-24")

step_income   <- pool %>% filter(INCOMG1 != "9")
pool          <- add_step(pool, step_income, "INCOMG1 = 'Not sure'/'Refused'")

step_menthlth <- pool %>%
  mutate(MENTHLTH = if_else(MENTHLTH == 88, 0, MENTHLTH)) %>%
  filter(!MENTHLTH %in% c(77, 99))
pool          <- add_step(pool, step_menthlth, "MENTHLTH = 'Not sure'/'Refused'")

step_exerany  <- pool %>% filter(EXERANY2 == "1")
pool          <- add_step(pool, step_exerany, "EXERANY2 = No/Not sure/Refused")

step_exract   <- pool %>% filter(!EXRACT12 %in% c("77", "99")) %>%
  mutate(EXRACT12 = recode_factor(EXRACT12,
                                  `1`  = "Walking",
                                  `2`  = "Running or jogging",
                                  `3`  = "Gardening or yard work",
                                  `4`  = "Bicycling",
                                  `5`  = "Aerobics",
                                  `6`  = "Calisthenics",
                                  `7`  = "Elliptical/EFX",
                                  `8`  = "Household activities",
                                  `9`  = "Weight lifting",
                                  `10` = "Yoga/Pilates/Tai Chi",
                                  `11` = "Other"))
pool          <- add_step(pool, step_exract, "EXRACT12 = 'Not sure'/'Refused'")

step_pafreq   <- pool %>% filter(!PAFREQ1_ %in% 99000)
pool          <- add_step(pool, step_pafreq, "PAFREQ1_ = 'Not sure'/'Refused'")

step_strfreq  <- pool %>% filter(!STRFREQ_ %in% 99000)
pool          <- add_step(pool, step_strfreq, "STRFREQ_ = 'Not sure'/'Refused'")

step_addepev3 <- pool %>% filter(!ADDEPEV3 %in% c("7", "9"))
pool          <- add_step(pool, step_addepev3, "ADDEPEV3 = 'Not sure'/'Refused'")

step_lsatisfy <- pool %>% filter(!LSATISFY %in% c("7", "9"))
pool          <- add_step(pool, step_lsatisfy, "LSATISFY = 'Not sure'/'Refused'")

step_emtsuprt <- pool %>% filter(!EMTSUPRT %in% c("7", "9"))
pool          <- add_step(pool, step_emtsuprt, "EMTSUPRT = 'Not sure'/'Refused'")

step_sdlonely <- pool %>% filter(!SDLONELY %in% c("7", "9"))
pool          <- add_step(pool, step_sdlonely, "SDLONELY = 'Not sure'/'Refused'")

step_sdhstre1 <- pool %>% filter(!SDHSTRE1 %in% c("7", "9"))
pool          <- add_step(pool, step_sdhstre1, "SDHSTRE1 = 'Not sure'/'Refused'")

step_sdhemply <- pool %>% filter(!SDHEMPLY %in% c("7", "9"))
pool          <- add_step(pool, step_sdhemply, "SDHEMPLY = 'Not sure'/'Refused'")

step_foodstmp <- pool %>% filter(!FOODSTMP %in% c("7", "9"))
pool          <- add_step(pool, step_foodstmp, "FOODSTMP = 'Not sure'/'Refused'")

step_sdhbills <- pool %>% filter(!SDHBILLS %in% c("7", "9"))
pool          <- add_step(pool, step_sdhbills, "SDHBILLS = 'Not sure'/'Refused'")

step_sdhutils <- pool %>% filter(!SDHUTILS %in% c("7", "9"))
pool          <- add_step(pool, step_sdhutils, "SDHUTILS = 'Not sure'/'Refused'")

step_sdhtrnsp <- pool %>% filter(!SDHTRNSP %in% c("7", "9"))
pool          <- add_step(pool, step_sdhtrnsp, "SDHTRNSP = 'Not sure'/'Refused'")

step_naomit   <- pool %>% na.omit()
pool          <- add_step(pool, step_naomit, "NA omitted (listwise deletion, remaining true NAs)")

step_smallpa  <- pool %>% filter(!EXRACT12 %in% c("Household activities", "Aerobics", "Elliptical/EFX"))
pool          <- add_step(pool, step_smallpa, "PA groups w/ small sample size (n < 15)")

cat("\n--- Waterfall: paste this directly into Figure 1 / Appendix C1 ---\n")
print(waterfall, n = Inf)
cat("\nStudent participants (start):", nrow(BRFSS_2023.small), "\n")
cat("Final analytic cohort (end):  ", nrow(pool), "\n")
cat("Sanity check - should equal 1,862 and match nrow(BRFSS_2023.recoded):",
    nrow(pool) == nrow(BRFSS_2023.recoded), "\n")

# Age-restricted pool (matches "Adult-aged college students" box in Figure 1)
cat("\nAge-restricted pool (post age-filter, pre everything else):", nrow(step_age), "\n")

# Included vs. excluded comparison, on demographics available pre-listwise-deletion.
# Uses .row_id (added above) rather than matching on full row content, because
# content-based anti_join() cannot distinguish different respondents who share
# identical answers across all ~10 categorical variables - which is common
# enough in this data to meaningfully distort the excluded/included counts.
excluded_df <- step_age %>% filter(!.row_id %in% pool$.row_id) %>%
  mutate(attrition_group = "Excluded")
included_df <- pool %>% mutate(attrition_group = "Included")

attrition_compare <- bind_rows(
  included_df %>% select(SEX, AGE80, IMPRACE, INCOMG1, attrition_group),
  excluded_df %>% select(SEX, AGE80, IMPRACE, INCOMG1, attrition_group)
)

cat("\n--- Supplementary Table S5: Included vs. Excluded ---\n")
cat("Sanity check - excluded + included should equal age-restricted pool (",
    nrow(step_age), "):", nrow(excluded_df) + nrow(included_df) == nrow(step_age), "\n")
tbl_attrition <- CreateTableOne(
  vars = c("SEX", "AGE80", "IMPRACE", "INCOMG1"), strata = "attrition_group",
  data = attrition_compare, factorVars = c("SEX", "IMPRACE", "INCOMG1")
)
print(tbl_attrition, showAllLevels = TRUE, varLabels = TRUE)
# NOTE: this table compares on SEX/AGE80/IMPRACE/INCOMG1 only, since those
# are the only variables observed for BOTH included and excluded
# respondents at this stage (MENTHLTH is also available - add it if desired,
# matching the original Appendix C1 layout).


# ---- H2. MENTHLTH distribution: histogram + density overlay ----
menthlth_hist <- ggplot(BRFSS_2023.recoded, aes(x = MENTHLTH)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 1,
                 fill = "#4575b4", alpha = 0.75, boundary = 0) +
  geom_density(color = "#d73027", linewidth = 0.8) +
  geom_vline(xintercept = mean(BRFSS_2023.recoded$MENTHLTH), linetype = "dashed") +
  labs(x = "Poor mental health days (past 30 days)", y = "Density") +
  theme_minimal(base_size = 11)
print(menthlth_hist)

.skew <- function(x) {
  n <- length(x); m <- mean(x); s <- sd(x)
  (sum((x - m)^3) / n) / s^3
}
cat("\nMENTHLTH skewness:", round(.skew(BRFSS_2023.recoded$MENTHLTH), 3), "\n")
cat("MENTHLTH mean:", round(mean(BRFSS_2023.recoded$MENTHLTH), 2),
    " SD:", round(sd(BRFSS_2023.recoded$MENTHLTH), 2), "\n")


# ---- H3. Negative binomial regression (sensitivity check) ----
model_5_formula <- formula(model_5)
nb_model <- MASS::glm.nb(model_5_formula, data = BRFSS_2023.recoded)

cat("\n--- Supplementary Table S6: Negative Binomial Regression (Model 5 spec.) ---\n")
nb_tidy <- tidy(nb_model, conf.int = TRUE, exponentiate = FALSE) %>%
  mutate(p_fmt = fmt_p(p.value)) %>%
  select(term, estimate, std.error, conf.low, conf.high, p_fmt)
print(nb_tidy, n = Inf)

ols_vs_nb <- left_join(
  tidy(model_5) %>% select(term, ols_estimate = estimate, ols_p = p.value),
  tidy(nb_model) %>% select(term, nb_estimate = estimate, nb_p = p.value),
  by = "term"
) %>%
  mutate(
    direction_match = sign(ols_estimate) == sign(nb_estimate),
    sig_match       = (ols_p < 0.05) == (nb_p < 0.05)
  )
print(ols_vs_nb, n = Inf)


# ---- H4. Psychosocial-predictor subset sensitivity analysis ----
base_rhs <- "SEX + AGE80 + IMPRACE + INCOMG1 + EXRACT12 + PAFREQ1_ + PADUR1_ + MINAC12 + STRFREQ_ + ADDEPEV3 + EMTSUPRT + SDHEMPLY + FOODSTMP + SDHBILLS + SDHUTILS + SDHTRNSP"

subset_formulas <- list(
  "Stress + Loneliness only"       = paste("MENTHLTH ~", base_rhs, "+ SDHSTRE1 + SDLONELY"),
  "Stress + Life satisfaction only" = paste("MENTHLTH ~", base_rhs, "+ SDHSTRE1 + LSATISFY"),
  "All three (stress+loneliness+life satisfaction; = Model 5)" =
    paste("MENTHLTH ~", base_rhs, "+ SDHSTRE1 + SDLONELY + LSATISFY")
)

# FIX 4: broom::tidy() does not reliably support `coeftest` objects across
# broom versions. Replaced with a manual data.frame conversion of the
# coeftest matrix (which always works, regardless of broom version).
tidy_coeftest <- function(ct) {
  as.data.frame(unclass(ct)) %>%
    rownames_to_column("term") %>%
    rename(estimate = Estimate, std.error = `Std. Error`,
           statistic = `t value`, p.value = `Pr(>|t|)`)
}

cat("\n--- Supplementary Table S7: Psychosocial-Predictor Subset Sensitivity ---\n")
subset_results <- lapply(names(subset_formulas), function(nm) {
  fit <- lm(as.formula(subset_formulas[[nm]]), data = BRFSS_2023.recoded)
  ct  <- coeftest(fit, vcov = vcovHC(fit, type = "HC3"))
  tidy_coeftest(ct) %>%
    filter(term == "EXRACT12Running or jogging") %>%
    mutate(subset = nm, adj_r2 = summary(fit)$adj.r.squared)
})
subset_results_df <- bind_rows(subset_results) %>%
  mutate(p_fmt = fmt_p(p.value)) %>%
  select(subset, term, estimate, std.error, p_fmt, adj_r2)
print(subset_results_df)


# ---- H5. Combined unweighted vs. weighted side-by-side table ----
cat("\n--- Supplementary Table S8: Unweighted vs. Weighted Model 5 Estimates ---\n")
combined_uw_w <- left_join(
  tidy(model_5, conf.int = TRUE) %>%
    select(term, uw_estimate = estimate, uw_p = p.value,
           uw_low = conf.low, uw_high = conf.high),
  tidy(svy_model, conf.int = TRUE) %>%
    select(term, w_estimate = estimate, w_p = p.value,
           w_low = conf.low, w_high = conf.high),
  by = "term"
) %>%
  mutate(
    uw_p_fmt = fmt_p(uw_p), w_p_fmt = fmt_p(w_p),
    diverges = (uw_p < 0.05) != (w_p < 0.05)
  )
print(combined_uw_w, n = Inf)


### ADDITIONAL ANALYSIS ###
attrition_compare %>%
  filter(INCOMG1 != "9") %>% droplevels() %>%
  CreateTableOne(vars = "INCOMG1", strata = "attrition_group",
                 factorVars = "INCOMG1", data = .)

nested_tests <- list(
  "Model 2 vs 0 (adds PA type)"             = anova(model_0, model_2),
  "Model 2 vs 1 (adds demographics)"        = anova(model_1, model_2),
  "Model 3 vs 2 (adds PA modality)"         = anova(model_2, model_3),
  "Model 4 vs 3 (adds psychosocial)"        = anova(model_3, model_4),
  "Model 5 vs 4 (adds income and hardship)" = anova(model_4, model_5))

f_table <- purrr::imap_dfr(nested_tests, function(a, nm) tibble(
  Comparison = nm, df1 = a$Df[2], df2 = a$Res.Df[2],
  F = round(a$F[2], 2), p = fmt_p(a$`Pr(>F)`[2])))
print(f_table)

lmtest::waldtest(model_3, model_4, test = "F",
                 vcov = function(m) sandwich::vcovHC(m, type = "HC3"))

table(step_menthlth$EXERANY2, useNA = "ifany")

