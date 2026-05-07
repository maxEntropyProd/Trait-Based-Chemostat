# Figure 2 stats

# SET UP ENVIRONMENT----
## Load necessary libraries----
library(dplyr); library(emmeans); library(multcompView); library(ggplot2)

## Read in data----
# Read in the data (dissolved O2 had a different output)
monitor_var <- read.csv("input_files/monitor_param_LONG.csv", header = T)
dissolved_oxy <- read.csv("input_files/dissolved_oxygen.csv", header = T)
nuts_data <- read.csv("input_files/new_nutrients.csv", header = TRUE)

# Models====================
## CO2 (%)----
# Create "time within phase" that restarts for each dilution rate within each model
monitor_var <- monitor_var %>%
  arrange(chemostat_ID, time) %>%
  group_by(chemostat_ID, dilution_rate) %>%
  mutate(time_within_phase = time - min(time, na.rm = TRUE)) %>%
  ungroup()

# Inspect to confirm it resets correctly
monitor_var %>%
  select(chemostat_ID, dilution_rate, time, time_within_phase) %>%
  arrange(chemostat_ID, time)

CO2_GLMM1 <- lmer(carbon_dio_perc ~ dilution_rate + (1 | chemostat_ID), data = monitor_var, REML = FALSE)
CO2_GLMM2 <- lmer(carbon_dio_perc ~ dilution_rate * time_within_phase + (1 | chemostat_ID), data = monitor_var, REML = FALSE)
CO2_GLMM3 <- lmer(carbon_dio_perc ~ dilution_rate + time_within_phase + (1 | chemostat_ID), data = monitor_var, REML = FALSE)
CO2_null <- lmer(carbon_dio_perc ~ 1 + (1 | chemostat_ID), data = monitor_var, REML = FALSE)

AICc(CO2_GLMM1, CO2_GLMM2, CO2_GLMM3, CO2_null)
model.sel(CO2_GLMM1, CO2_GLMM2, CO2_GLMM3, CO2_null)
summary(CO2_GLMM2)

# Now to get posthoc information
emm_dilution <- emmeans(CO2_GLMM2, ~ dilution_rate)
emm_dilution
pairs(emm_dilution, adjust = "tukey")
confint(pairs(emm_dilution, adjust = "tukey"))

# report emtrends due to significant interaction
co2_trends <- emtrends(CO2_GLMM2,~ dilution_rate, var = "time_within_phase")
co2_trends
pairs(co2_trends)

## O2 (%)----
# Can use the same dataset that already has a "time within phase" column added

O2_GLMM1 <- lmer(oxygen_perc ~ dilution_rate + (1 | chemostat_ID), data = monitor_var, REML = FALSE)
O2_GLMM2 <- lmer(oxygen_perc ~ dilution_rate * time_within_phase + (1 | chemostat_ID), data = monitor_var, REML = FALSE)
O2_GLMM3 <- lmer(oxygen_perc ~ dilution_rate + time_within_phase + (1 | chemostat_ID), data = monitor_var, REML = FALSE)
O2_null <- lmer(oxygen_perc ~ 1 + (1 | chemostat_ID), data = monitor_var, REML = FALSE)

AICc(O2_GLMM1, O2_GLMM2, O2_GLMM3, O2_null)
model.sel(O2_GLMM1, O2_GLMM2, O2_GLMM3, O2_null)
summary(O2_GLMM2)

# Now to get posthoc information
emm_dilution <- emmeans(O2_GLMM2, ~ dilution_rate)
emm_dilution
pairs(emm_dilution, adjust = "tukey")
confint(pairs(emm_dilution, adjust = "tukey"))

# report emtrends due to significant interaction
o2_trends <- emtrends(O2_GLMM2,~ dilution_rate, var = "time_within_phase")
o2_trends
pairs(o2_trends)

## pH ----
# Can use the same dataset that already has a "time within phase" column added

ph_GLMM1 <- lmer(ph ~ dilution_rate + (1 | chemostat_ID), data = monitor_var, REML = FALSE)
ph_GLMM2 <- lmer(ph ~ dilution_rate * time_within_phase + (1 | chemostat_ID), data = monitor_var, REML = FALSE)
ph_GLMM3 <- lmer(ph ~ dilution_rate + time_within_phase + (1 | chemostat_ID), data = monitor_var, REML = FALSE)
ph_null <- lmer(ph ~ 1 + (1 | chemostat_ID), data = monitor_var, REML = FALSE)

AICc(ph_GLMM1, ph_GLMM2, ph_GLMM3, ph_null)
model.sel(ph_GLMM1, ph_GLMM2, ph_GLMM3, ph_null)
summary(ph_GLMM2)

# Now to get posthoc information
emm_dilution <- emmeans(ph_GLMM2, ~ dilution_rate)
emm_dilution
pairs(emm_dilution, adjust = "tukey")
confint(pairs(emm_dilution, adjust = "tukey"))

# report emtrends due to significant interaction
ph_trends <- emtrends(ph_GLMM2,~ dilution_rate, var = "time_within_phase")
ph_trends
pairs(ph_trends)

## Dissolved oxy ----
# Create "time within phase" that restarts for each dilution rate within each model
dissolved_oxy <- dissolved_oxy %>%
  arrange(chemostat_ID, time) %>%
  group_by(chemostat_ID, dilution_rate) %>%
  mutate(time_within_phase = time - min(time, na.rm = TRUE)) %>%
  ungroup()

# Inspect to confirm it resets correctly
dissolved_oxy %>%
  select(chemostat_ID, dilution_rate, time, time_within_phase) %>%
  arrange(chemostat_ID, time)

do_GLMM1 <- lmer(diss_o2 ~ dilution_rate + (1 | chemostat_ID), data = dissolved_oxy, REML = FALSE)
do_GLMM2 <- lmer(diss_o2 ~ dilution_rate * time_within_phase + (1 | chemostat_ID), data = dissolved_oxy, REML = FALSE)
do_GLMM3 <- lmer(diss_o2 ~ dilution_rate + time_within_phase + (1 | chemostat_ID), data = dissolved_oxy, REML = FALSE)
do_null <- lmer(diss_o2 ~ 1 + (1 | chemostat_ID), data = dissolved_oxy, REML = FALSE)

AICc(do_GLMM1, do_GLMM2, do_GLMM3, do_null)
model.sel(do_GLMM1, do_GLMM2, do_GLMM3, do_null)
summary(do_GLMM2)

# Now to get posthoc information
emm_dilution <- emmeans(do_GLMM2, ~ dilution_rate)
emm_dilution
pairs(emm_dilution, adjust = "tukey")
confint(pairs(emm_dilution, adjust = "tukey"))

# report emtrends due to significant interaction
do_trends <- emtrends(do_GLMM2,~ dilution_rate, var = "time_within_phase")
do_trends
pairs(do_trends)

## Nitrate ----
# create a "day within phase" column that restarts with each dilution rate
nuts_data <- nuts_data %>%
  arrange(chemostat_ID, days_after)

nuts_data <- nuts_data %>%
  group_by(chemostat_ID, dilution_rate) %>%
  mutate(
    day_within_phase = days_after - min(days_after)
  ) %>%
  ungroup()

nuts_data %>%
  select(sample_ID, chemostat_ID, dilution_rate, days_after, day_within_phase) %>%
  arrange(chemostat_ID, days_after)

nitrate_GLMM1 <- lmer(nitrate_conc ~ dilution_rate + (1 | chemostat_ID), data = nuts_data, REML = FALSE)
nitrate_GLMM2 <- lmer(nitrate_conc ~ dilution_rate * day_within_phase + (1 | chemostat_ID), data = nuts_data, REML = FALSE)
nitrate_GLMM3 <- lmer(nitrate_conc ~ dilution_rate + day_within_phase + (1 | chemostat_ID), data = nuts_data, REML = FALSE)
nitrate_null <- lmer(nitrate_conc ~ 1 + (1 | chemostat_ID), data = nuts_data, REML = FALSE)

AICc(nitrate_GLMM1, nitrate_GLMM2, nitrate_GLMM3, nitrate_null)
model.sel(nitrate_GLMM1, nitrate_GLMM2, nitrate_GLMM3, nitrate_null)
summary(nitrate_GLMM3)

# Now to get posthoc information
emm_dilution <- emmeans(nitrate_GLMM3, ~ dilution_rate)
emm_dilution
pairs(emm_dilution, adjust = "tukey")
confint(pairs(emm_dilution, adjust = "tukey"))

## Phosphate ----
# Can use the same day within phase created in "nuts_data" on Ln 128
phosphate_GLMM1 <- lmer(phos_conc ~ dilution_rate + (1 | chemostat_ID), data = nuts_data, REML = FALSE)
phosphate_GLMM2 <- lmer(phos_conc ~ dilution_rate * day_within_phase + (1 | chemostat_ID), data = nuts_data, REML = FALSE)
phosphate_GLMM3 <- lmer(phos_conc ~ dilution_rate + day_within_phase + (1 | chemostat_ID), data = nuts_data, REML = FALSE)
phosphate_null <- lmer(phos_conc ~ 1 + (1 | chemostat_ID), data = nuts_data, REML = FALSE)

AICc(phosphate_GLMM1, phosphate_GLMM2, phosphate_GLMM3, phosphate_null)
model.sel(phosphate_GLMM1, phosphate_GLMM2, phosphate_GLMM3, phosphate_null)
summary(phosphate_GLMM3)

# Now to get posthoc information
emm_dilution <- emmeans(phosphate_GLMM3, ~ dilution_rate)
emm_dilution
pairs(emm_dilution, adjust = "tukey")
confint(pairs(emm_dilution, adjust = "tukey"))

## Ammonium ----
# Can use the same day within phase created in "nuts_data" on Ln 128
ammonium_GLMM1 <- lmer(ammonium_con ~ dilution_rate + (1 | chemostat_ID), data = nuts_data, REML = FALSE)
ammonium_GLMM2 <- lmer(ammonium_con ~ dilution_rate * day_within_phase + (1 | chemostat_ID), data = nuts_data, REML = FALSE)
ammonium_GLMM3 <- lmer(ammonium_con ~ dilution_rate + day_within_phase + (1 | chemostat_ID), data = nuts_data, REML = FALSE)
ammonium_null <- lmer(ammonium_con ~ 1 + (1 | chemostat_ID), data = nuts_data, REML = FALSE)

AICc(ammonium_GLMM1, ammonium_GLMM2, ammonium_GLMM3, ammonium_null)
model.sel(ammonium_GLMM1, ammonium_GLMM2, ammonium_GLMM3, ammonium_null)
summary(ammonium_GLMM3)

# Now to get posthoc information
emm_dilution <- emmeans(ammonium_GLMM3, ~ dilution_rate)
emm_dilution
pairs(emm_dilution, adjust = "tukey")
confint(pairs(emm_dilution, adjust = "tukey"))

## DOC ----
# Can use the same day within phase created in "nuts_data" on Ln 128
doc_GLMM1 <- lmer(doc_conc_mmol ~ dilution_rate + (1 | chemostat_ID), data = nuts_data, REML = FALSE)
doc_GLMM2 <- lmer(doc_conc_mmol ~ dilution_rate * day_within_phase + (1 | chemostat_ID), data = nuts_data, REML = FALSE)
doc_GLMM3 <- lmer(doc_conc_mmol ~ dilution_rate + day_within_phase + (1 | chemostat_ID), data = nuts_data, REML = FALSE)
doc_null <- lmer(doc_conc_mmol ~ 1 + (1 | chemostat_ID), data = nuts_data, REML = FALSE)

AICc(doc_GLMM1, doc_GLMM2, doc_GLMM3, doc_null)
model.sel(doc_GLMM1, doc_GLMM2, doc_GLMM3, doc_null)
summary(doc_GLMM1)
summary(doc_GLMM2) # weaker support for interactive model, but include for consistency?

# Now to get posthoc information
emm_dilution <- emmeans(doc_GLMM1, ~ dilution_rate)
emm_dilution
pairs(emm_dilution, adjust = "tukey")
confint(pairs(emm_dilution, adjust = "tukey"))

# report emtrends due to significant interaction (if we go with interactive model)
doc_trends <- emtrends(doc_GLMM2,~ dilution_rate, var = "day_within_phase")
doc_trends
pairs(doc_trends)

