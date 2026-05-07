# Figure 2 
# Environmental parameters from chemostats (pH, CO2, O2, Dissolved O2)
# Nutrients (nitrate, ammonium, phosphate, DOC)

# SET UP ENVIRONMENT----
# Load necessary libraries
library("ggplot2"); library(patchwork); library(lme4); library(emmeans)
library(sjstats); library(lmerTest); library(MuMIn); library(scales)

# Set figure themes
# This is a slightly different theme function from others 
pretty.theme <- function(){
  theme_bw()+
    theme(axis.text.x=element_text(size=14, color = "black", angle = 0),
          axis.text.y=element_text(size=14, color = "black"),
          axis.title.x=element_text(size=14, color = "black"),             
          axis.title.y=element_text(size=14, color = "black"),             
          panel.grid.major.x=element_blank(),                                          
          panel.grid.minor.x=element_blank(),
          panel.grid.minor.y=element_blank(),
          panel.grid.major.y=element_blank(),  
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), units = , "cm"),
          plot.title = element_text(size=14))
}
# Read in the data (dissolved O2 had a different output)
monitor_var <- read.csv("input_files/monitor_param_LONG.csv", header = T)
dissolved_oxy <- read.csv("input_files/dissolved_oxygen.csv", header = T)
nuts_data <- read.csv("input_files/new_nutrients.csv", header = TRUE)


## chemostat environment----
### CO2----
CO2_plot <- ggplot(monitor_var, group = chemostat_ID, aes(x = time, y = carbon_dio_perc)) +
  geom_line(aes(color = dilution_rate, linetype = chemostat_ID), linewidth = 1.0) +
  geom_vline(xintercept = c(15.71,20.79), color = "black", linetype = "dashed") +
  scale_color_manual(values = c("#2B5B6C", "#E34F33","#FFC87E")) +
  pretty.theme() +
  labs (x = NULL, y = expression(CO[2] ~ "(%)"), linetype = "Chemostat ID", color = "Dilution Rate") +
  theme(axis.text.x = element_blank())
CO2_plot

### O2----
O2_plot <- ggplot(monitor_var, group = chemostat_ID, aes(x = time, y = oxygen_perc)) +
  geom_line(aes(color = dilution_rate, linetype = chemostat_ID), linewidth = 1.0) +
  geom_vline(xintercept = c(15.71,20.79), color = "black", linetype = "dashed") +
  scale_color_manual(values = c("#2B5B6C", "#E34F33","#FFC87E")) +
  scale_y_continuous(labels = number_format(accuracy = 0.1)) +  # 1 decimal place only
  pretty.theme() +
  labs(x = NULL, y = expression(O[2] ~ "(%)"), linetype = "Chemostat ID", color = "Dilution Rate") +
  theme(axis.text.x = element_blank())
O2_plot

### Diss O2----
diss_O2_plot <- ggplot(dissolved_oxy, group = chemostat_ID, aes(x = time, y = diss_o2)) +
  geom_line(aes(color = dilution_rate, linetype = chemostat_ID), linewidth = 1.0) +
  geom_vline(xintercept = c(15.71,20.79), color = "black", linetype = "dashed") +
  scale_color_manual(values = c("#2B5B6C", "#E34F33","#FFC87E")) +
  pretty.theme() +
  labs (x = NULL, y = expression(O[2] ~ "(mg" ~ L^{-1} * ")"), linetype = "Chemostat ID", color = "Dilution Rate") +
  theme(axis.text.x = element_blank())
diss_O2_plot

### pH----
ph_plot <- ggplot(monitor_var, group = chemostat_ID, aes(x = time, y = ph)) +
  geom_line(aes(color = dilution_rate, linetype = chemostat_ID), linewidth = 1.0) +
  geom_vline(xintercept = c(15.71,20.79), color = "black", linetype = "dashed") +
  scale_color_manual(values = c("#2B5B6C", "#E34F33","#FFC87E")) +
  pretty.theme() +
  labs (x = NULL, y = "pH", linetype = "Chemostat ID", color = "Dilution Rate") +
  theme(axis.text.x = element_blank())
ph_plot

## combine----
fig2 <- (CO2_plot / O2_plot | diss_O2_plot / ph_plot) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right",
        legend.justification = "center")
fig2

## nutrients----
### nitrate----
nitrate_plot <- ggplot(nuts_data, aes(x = days_after, y = nitrate_conc)) +
  geom_point(aes(shape = chemostat_ID, color = dilution_rate), size = 2) +
  geom_line(aes(color = dilution_rate, linetype = chemostat_ID)) +
  geom_vline(xintercept = c(15.71,20.79), color = "black", linetype = "dashed") +
  scale_color_manual(values = c("#2B5B6C", "#E34F33","#FFC87E")) +
  pretty.theme() +
  labs (x = NULL, y = expression(NO[3]^"-" ~ "(" * mu * "mol L"^{-1} * ")"), shape = "Chemostat ID", color = "Dilution Rate", linetype = NULL) +
  theme(axis.text.x = element_blank())
nitrate_plot

### ammonium----
ammonium_plot <- ggplot(nuts_data, aes(x = days_after, y = ammonium_con)) +
  geom_point(aes(shape = chemostat_ID, color = dilution_rate), size = 2) +
  geom_line(aes(color = dilution_rate, linetype = chemostat_ID)) +
  geom_vline(xintercept = c(15.71,20.79), color = "black", linetype = "dashed") +
  scale_color_manual(values = c("#2B5B6C", "#E34F33","#FFC87E")) +
  pretty.theme() +
  labs(x = NULL, y = expression(NH[4]^"+" ~ "(" * mu * "mol L"^{-1} * ")"), shape = "Chemostat ID", color = "Dilution Rate", linetype = NULL) +
  theme(axis.text.x = element_blank())
ammonium_plot

### phosphate----
phosphate_plot <- ggplot(nuts_data, aes(x = days_after, y = phos_conc)) +
  geom_point(aes(shape = chemostat_ID, color = dilution_rate), size = 2) +
  geom_line(aes(color = dilution_rate, linetype = chemostat_ID)) +
  geom_vline(xintercept = c(15.71,20.79), color = "black", linetype = "dashed") +
  scale_color_manual(values = c("#2B5B6C", "#E34F33","#FFC87E")) +
  pretty.theme() +
  labs(x = NULL, y = expression(PO[4]^"3-" ~ "( " * mu * "mol L"^{-1} * " )"), shape = "Chemostat ID", color = "Dilution Rate", linetype = NULL) +
  theme(axis.text.x = element_blank())
phosphate_plot


### DOC----
DOC_plot <- ggplot(nuts_data, aes(x = days_after, y = doc_conc_mmol)) +
  geom_point(aes(shape = chemostat_ID, color = dilution_rate), size = 2) +
  geom_line(aes(color = dilution_rate, linetype = chemostat_ID)) +
  geom_vline(xintercept = c(15.71,20.79), color = "black", linetype = "dashed") +
  scale_color_manual(values = c("#2B5B6C", "#E34F33","#FFC87E")) +
  pretty.theme() +
  labs(x = "Days after start", y = expression("DOC (mmol L"^"-1" * ")"), shape = "Chemostat ID", color = "Dilution Rate", linetype = NULL) 
DOC_plot

## Combine plots----
fig3 <- (nitrate_plot / ammonium_plot | phosphate_plot / DOC_plot) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right",
        legend.justification = "center")
fig3

combine_2_3 <- fig2 / fig3
combine_2_3


