# Figure 4 stats

# SET UP ENVIRONMENT----
# Load necessary libraries
library(dplyr); library(emmeans); library(multcompView); library(ggplot2)

# Shannon diversity - observed=============================
# Read in data
alpha_div_data <- read.csv("input_files/alpha_diversity_rare_data.csv", header = TRUE)

# create a "day within phase" column that restarts with each dilution rate
alpha_div_data <- alpha_div_data %>%
  arrange(model, days_after)

alpha_div_data <- alpha_div_data %>%
  group_by(model, dilution_rate) %>%
  mutate(
    day_within_phase = days_after - min(days_after)
  ) %>%
  ungroup()

alpha_div_data %>%
  select(SampleID, model, dilution_rate, days_after, day_within_phase) %>%
  arrange(model, days_after)

alpha_obs_model_GLMM1 <- lmer(Shannon ~ dilution_rate + (1 | model), data = alpha_div_data, REML = FALSE)
alpha_obs_model_GLMM2 <- lmer(Shannon ~ dilution_rate * day_within_phase + (1 | model), data = alpha_div_data, REML = FALSE)
alpha_obs_model_GLMM3 <- lmer(Shannon ~ dilution_rate + day_within_phase + (1 | model), data = alpha_div_data, REML = FALSE)
alpha_obs_model_null1 <- lmer(Shannon ~ 1 + (1 | model), data = alpha_div_data, REML = FALSE)

AICc(alpha_obs_model_GLMM1, alpha_obs_model_GLMM2, alpha_obs_model_GLMM3, alpha_obs_model_null1)
model.sel(alpha_obs_model_GLMM1, alpha_obs_model_GLMM2, alpha_obs_model_GLMM3, alpha_obs_model_null1)
summary(alpha_obs_model_GLMM1)

# Now to get posthoc information
emm_dilution <- emmeans(alpha_obs_model_GLMM1, ~ dilution_rate)
emm_dilution
pairs(emm_dilution, adjust = "tukey")
confint(pairs(emm_dilution, adjust = "tukey"))

# Community stability - observed==============================
# Read in data
stability_data <- read.csv("input_files/BRAY_DISS_updated.csv", header = TRUE)

# create a "day within phase" column that restarts with each dilution rate
stability_data <- stability_data %>%
  arrange(model, days_after)

stability_data <- stability_data %>%
  group_by(model, DILUTION_RATE) %>%
  mutate(
    day_within_phase = days_after - min(days_after)
  ) %>%
  ungroup()

stability_data %>%
  select(SAMPLE_ID, model, DILUTION_RATE, days_after, day_within_phase) %>%
  arrange(model, days_after)

beta_obs_model_GLMM1 <- lmer(SIM ~ DILUTION_RATE + (1 | model), data = stability_data, REML = FALSE)
beta_obs_model_GLMM2 <- lmer(SIM ~ DILUTION_RATE * day_within_phase + (1 | model), data = stability_data, REML = FALSE)
beta_obs_model_GLMM3 <- lmer(SIM ~ DILUTION_RATE + day_within_phase + (1 | model), data = stability_data, REML = FALSE)
beta_obs_model_null1 <- lmer(SIM ~ 1 + (1 | model), data = stability_data, REML = FALSE)

AICc(beta_obs_model_GLMM1, beta_obs_model_GLMM2, beta_obs_model_GLMM3, beta_obs_model_null1)
model.sel(beta_obs_model_GLMM1, beta_obs_model_GLMM2, beta_obs_model_GLMM3, beta_obs_model_null1)
summary(beta_obs_model_GLMM1)

# Now to get posthoc information
emm_dilution <- emmeans(beta_obs_model_GLMM1, ~ DILUTION_RATE)
emm_dilution
pairs(emm_dilution, adjust = "tukey")
confint(pairs(emm_dilution, adjust = "tukey"))


# Shannon diversity - modeled================================
# Read in data
alpha_modeled <- read.csv("input_files/alpha_diversity_modeled.csv", header = T)
head(alpha_modeled)

# Subset by distribution (uniform versus beta)
alpha_uniform <- alpha_modeled %>%
  filter(distribution == "uniform")

alpha_beta <- alpha_modeled %>%
  filter(distribution == "beta")

## Uniform==========
# create a "day within phase" column that restarts with each dilution rate
alpha_uniform <- alpha_uniform %>%
  arrange(model, days_after)

alpha_uniform <- alpha_uniform %>%
  group_by(model, dilution_rate) %>%
  mutate(
    day_within_phase = days_after - min(days_after)
  ) %>%
  ungroup()

alpha_uniform %>%
  select(SampleID, model, dilution_rate, days_after, day_within_phase) %>%
  arrange(model, days_after)

head(alpha_uniform)

alpha_unimod_GLMM1 <- lmer(shannon ~ dilution_rate + (1 | model), data = alpha_uniform, REML = FALSE)
alpha_unimod_GLMM2 <- lmer(shannon ~ dilution_rate * day_within_phase + (1 | model), data = alpha_uniform, REML = FALSE)
alpha_unimod_GLMM3 <- lmer(shannon ~ dilution_rate + day_within_phase + (1 | model), data = alpha_uniform, REML = FALSE)
alpha_unimod_null1 <- lmer(shannon ~ 1 + (1 | model), data = alpha_uniform, REML = FALSE)

AICc(alpha_unimod_GLMM1, alpha_unimod_GLMM2, alpha_unimod_GLMM3, alpha_unimod_null1)
model.sel(alpha_unimod_GLMM1, alpha_unimod_GLMM2, alpha_unimod_GLMM3, alpha_unimod_null1)
summary(alpha_unimod_GLMM2)

# Now to get posthoc information
emm_dilution <- emmeans(alpha_unimod_GLMM2, ~dilution_rate) 
pairs(emm_dilution, adjust = "tukey")

trends <- emtrends(alpha_unimod_GLMM2, ~ dilution_rate, var = "day_within_phase")
pairs(trends, adjust = "tukey")

emm_dilution
pairs(emm_dilution, adjust = "tukey")
confint(pairs(emm_dilution, adjust = "tukey"))

## Beta===========
# create a "day within phase" column that restarts with each dilution rate
alpha_beta <- alpha_beta %>%
  arrange(model, days_after)

alpha_beta <- alpha_beta %>%
  group_by(model, dilution_rate) %>%
  mutate(
    day_within_phase = days_after - min(days_after)
  ) %>%
  ungroup()

alpha_beta %>%
  select(SampleID, model, dilution_rate, days_after, day_within_phase) %>%
  arrange(model, days_after)

head(alpha_beta)

alpha_betamod_GLMM1 <- lmer(shannon ~ dilution_rate + (1 | model), data = alpha_beta, REML = FALSE)
alpha_betamod_GLMM2 <- lmer(shannon ~ dilution_rate * day_within_phase + (1 | model), data = alpha_beta, REML = FALSE)
alpha_betamod_GLMM3 <- lmer(shannon ~ dilution_rate + day_within_phase + (1 | model), data = alpha_beta, REML = FALSE)
alpha_betamod_null1 <- lmer(shannon ~ 1 + (1 | model), data = alpha_beta, REML = FALSE)

AICc(alpha_betamod_GLMM1, alpha_betamod_GLMM2, alpha_betamod_GLMM3, alpha_betamod_null1)
model.sel(alpha_betamod_GLMM1, alpha_betamod_GLMM2, alpha_betamod_GLMM3, alpha_betamod_null1)
summary(alpha_betamod_GLMM3)

# Now to get posthoc information
emm_dilution <- emmeans(alpha_betamod_GLMM3, ~dilution_rate) 
pairs(emm_dilution, adjust = "tukey")
confint(pairs(emm_dilution, adjust = "tukey"))

# Beta diversity - modeled ============================
# Read in data
stability_modeled <- read.csv("input_files/model_dissimilarity_data_within_dilutions.csv", header = T)

# Only sequential data
stability_modeled <- stability_modeled %>%
  filter(sequential == "YES")

# Subset by distribution (uniform versus beta)
stability_modeled_uniform <- stability_modeled %>%
  filter(distribution == "uniform")
head(stability_modeled_uniform)

stability_modeled_beta <- stability_modeled %>%
  filter(distribution == "beta")

## Uniform==============
# create a "day within phase" column that restarts with each dilution rate
stability_modeled_uniform <- stability_modeled_uniform %>%
  arrange(model, days_after)

stability_modeled_uniform <- stability_modeled_uniform %>%
  group_by(model, dilution_rate) %>%
  mutate(
    day_within_phase = days_after - min(days_after)
  ) %>%
  ungroup()

stability_modeled_uniform %>%
  select(sampleID, model, dilution_rate, days_after, day_within_phase) %>%
  arrange(model, days_after)

head(stability_modeled_uniform)

# Models
beta_unimod_GLMM1 <- lmer(sim ~ dilution_rate + (1 | model), data = stability_modeled_uniform, REML = FALSE)
beta_unimod_GLMM2 <- lmer(sim ~ dilution_rate * day_within_phase + (1 | model), data = stability_modeled_uniform, REML = FALSE)
beta_unimod_GLMM3 <- lmer(sim ~ dilution_rate + day_within_phase + (1 | model), data = stability_modeled_uniform, REML = FALSE)
beta_unimod_null1 <- lmer(sim ~ 1 + (1 | model), data = stability_modeled_uniform, REML = FALSE)

AICc(beta_unimod_GLMM1, beta_unimod_GLMM2, beta_unimod_GLMM3, beta_unimod_null1)
model.sel(beta_unimod_GLMM1, beta_unimod_GLMM2, beta_unimod_GLMM3, beta_unimod_null1)
summary(beta_unimod_GLMM3) # null is best supported model but if we go with additive model

# Now to get posthoc information
emm_dilution <- emmeans(beta_unimod_GLMM3, ~dilution_rate) 
pairs(emm_dilution, adjust = "tukey")
# suggests a difference between 0.01 and 1, but not between 0.1 and 10, and 1 and 10

## Beta==============
# create a "day within phase" column that restarts with each dilution rate
# create a "day within phase" column that restarts with each dilution rate
stability_modeled_beta <- stability_modeled_beta %>%
  arrange(model, days_after)

stability_modeled_beta <- stability_modeled_beta %>%
  group_by(model, dilution_rate) %>%
  mutate(
    day_within_phase = days_after - min(days_after)
  ) %>%
  ungroup()

stability_modeled_beta %>%
  select(sampleID, model, dilution_rate, days_after, day_within_phase) %>%
  arrange(model, days_after)

head(stability_modeled_beta)

# Models
beta_betamod_GLMM1 <- lmer(sim ~ dilution_rate + (1 | model), data = stability_modeled_beta, REML = FALSE)
beta_betamod_GLMM2 <- lmer(sim ~ dilution_rate * day_within_phase + (1 | model), data = stability_modeled_beta, REML = FALSE)
beta_betamod_GLMM3 <- lmer(sim ~ dilution_rate + day_within_phase + (1 | model), data = stability_modeled_beta, REML = FALSE)
beta_betamod_null1 <- lmer(sim ~ 1 + (1 | model), data = stability_modeled_beta, REML = FALSE)

AICc(beta_betamod_GLMM1, beta_betamod_GLMM2, beta_betamod_GLMM3, beta_betamod_null1)
model.sel(beta_betamod_GLMM1, beta_betamod_GLMM2, beta_betamod_GLMM3, beta_betamod_null1)
summary(beta_betamod_GLMM3)

# Now to get posthoc information
emm_dilution <- emmeans(beta_unimod_GLMM3, ~dilution_rate) 
pairs(emm_dilution, adjust = "tukey")
# suggests a difference between 0.01 and 1, but not between 0.1 and 10, and 1 and 10
