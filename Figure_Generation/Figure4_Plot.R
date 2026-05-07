# Figure 4
# Refer to "Final_Phyloseq_Feb2024.R" for Phyloseq object generation

# SET UP ENVIRONMENT----

# Load necessary libraries
library("ggplot2"); library(phyloseq); library(ggrepel); library(patchwork)
library(car); library(rstatix); library(ggpubr); library(lme4); library(emmeans)
library(sjstats); library(lmerTest); library(MuMIn)

# For plot consistency 
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

##

# PLOTTING----
## Panels A&B: PCoA plot for MC1----
### Files----
# Read in phyloseq objects as necessary
# per_MC1 <- readRDS("phyloseq_objects/phy_taxaFilter_relAbund_noPond_MC1.rds")
# per_MC2 <- readRDS("phyloseq_objects/phy_taxaFilter_relAbund_noPond_MC2.rds")
# per_nopond <- readRDS("phyloseq_objects/phy_taxaFilter_relAbund_noPond.rds")

# Bray curtis for MC1 specifically
BC_distance1 <- phyloseq::distance(per_MC1, "bray")
bcOrd1 <- ordinate(per_MC1, "PCoA", BC_distance1)
plot_scree(bcOrd1)
bray_plot1 <- plot_ordination(per_MC1, bcOrd1) +
  geom_point(aes(fill = dilution_rate), size = 3, alpha = 1, pch = 21, color = "black") +
  geom_polygon(aes(fill = dilution_rate), alpha = 0.8) +
  geom_text_repel(aes(label= total_days_after), 
                   color = "black",
                   # vjust = 0.02, nudge_y = 0.01,
                   direction = "both",
                   size = 4,
                   min.segment.length = Inf) +
  scale_fill_manual(values=c("#2B5B6C", "#E34F33", "#FFC87E","black")) +
  scale_color_manual(values=c("#2B5B6C", "#E34F33", "#FFC87E","black")) +
  pretty.theme() +
  theme(legend.title = element_blank(), 
        legend.position = "none") +
  ggtitle("A. PCoA for MC1 - Observed")
bray_plot1

# Bray curtis for MC2 specifically
BC_distance2 <- phyloseq::distance(per_MC2, "bray")
bcOrd2 <- ordinate(per_MC2, "PCoA", BC_distance2)
plot_scree(bcOrd2)
bray_plot2 <- plot_ordination(per_MC2, bcOrd2) +
  geom_point(aes(fill = dilution_rate), size = 3, alpha = 1, pch = 21, color = "black") +
  geom_polygon(aes(fill = dilution_rate), alpha = 0.8) +
  geom_text_repel(aes(label= total_days_after), 
                   color = "black",
                   # vjust = 0.02, nudge_y = 0.01,
                   direction = "both",
                   size = 4,
                   min.segment.length = Inf) +
  scale_fill_manual(values=c("#2B5B6C", "#E34F33", "#FFC87E","black")) +
  scale_color_manual(values=c("#2B5B6C", "#E34F33", "#FFC87E","black")) +
  pretty.theme() +
  theme(legend.title = element_blank()) +
  ggtitle("B. PCoA for MC2 - Observed")
bray_plot2

## Panel C: Alpha diversity for MC1 and MC2 (Observed)----
# This is for alpha diversity - chemostats combined
alpha_diversity <- estimate_richness(phy.f.nopond, measures = "Shannon")
alpha_diversity_rare <- estimate_richness(phy.f.rare.nopond, measures = c("Shannon", "Observed"))
write.csv(alpha_diversity_rare, "intermediate_files/alpha_diversity_rare.csv")

alpha_div_data <- read.csv("input_files/alpha_diversity_rare_data.csv", header = TRUE)

alpha_boxplot_ALL <- ggplot(alpha_div_data, aes(x = dilution_rate, y = Shannon)) +
  geom_boxplot(aes(fill = dilution_rate), color = "black", alpha = 0.9, width = 0.25) +
  # geom_jitter(aes(color = dilution_rate, shape = chemostat_ID), width = 0.2, size = 3, alpha = 0.9) +
  scale_color_manual(values=c("#2B5B6C", "#E34F33", "#FFC87E")) +
  scale_fill_manual(values=c("#2B5B6C", "#E34F33", "#FFC87E")) +
  xlab("Dilution Rate") +
  ylab("Shannon Diversity") +
  pretty.theme() +
  ggtitle("C. Alpha Diversity - Observed") +
  theme(strip.background = element_blank(), 
        strip.text = element_blank(),
        legend.position = "none") +
  ylim(0,7) 
alpha_boxplot_ALL

## Panel D: Dissimilarity for MC1 and MC2 (Observed)----
# Need to add in here how dissimilarity was calculated in sequence
# Load data
dis_box_data_ALL<- read.csv("input_files/BRAY_DISS_updated.csv", header = TRUE)

dis_boxplot_ALL <- ggplot(dis_box_data_ALL, aes(x = DILUTION_RATE, y = SIM), ) +
  geom_boxplot(aes(fill = DILUTION_RATE), color = "black", alpha = 0.9, width = 0.25) +
  # geom_jitter(aes(color = DILUTION_RATE, shape = CHEMOSTAT_ID), width = 0.2, size = 3, alpha = 0.9) +
  scale_color_manual(values=c("#2B5B6C", "#E34F33", "#FFC87E")) +
  scale_fill_manual(values=c("#2B5B6C", "#E34F33", "#FFC87E")) +
  xlab("Dilution Rate") +
  ylab("Similarity Metric") +
  pretty.theme() +
  ggtitle("D. Community Stability - Observed") +
  ylim(0,1) +
  theme(legend.position = "none")
dis_boxplot_ALL

## Panel E: Alpha Diversity (Modeled)
# We ran two iterations of the model to represent variation that may occur
# between chemostats but they are not mean to specifically represent
# chemostats 1 and 2 so will be named "model A, model B" and either uniform or beta

alpha_modeled <- read.csv("input_files/alpha_diversity_modeled.csv", header = T)
head(alpha_modeled)

# This is the correct plot with the two models, A & B
alpha_model_plot <- ggplot(alpha_modeled, aes(x =dilution_rate, y = shannon), ) +
  geom_boxplot(aes(fill = dilution_rate), color = "black", alpha = 0.9, width = 0.5) +
  # geom_jitter(aes(color = dilution_rate, shape = model), width = 0.2, size = 3, alpha = 0.9) +
  scale_color_manual(values=c("#2B5B6C", "#E34F33", "#FFC87E")) +
  scale_fill_manual(values=c("#2B5B6C", "#E34F33", "#FFC87E")) +
  xlab("Dilution Rate") +
  ylab("Shannon Diversity") +
  pretty.theme() +
  ggtitle("E. Alpha Diversity - Modeled") +
  # ylim(0,1) +
  facet_wrap(.~factor(distribution, levels=c("uniform", "beta"))) +
  # this changes the order of the facet wrap
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 14, colour = "black"),
        legend.position = "none") +
  ylim(0,7)
alpha_model_plot

## Panel F: Beta Diversity (Modeled)----
# We ran two iterations of the model to represent variation that may occur
# between chemostats but they are not mean to specifically represent
# chemostats 1 and 2 so will be named "model A, model B" and either uniform or beta

# Updated model beta diversity plots----
model_beta <- read.csv("input_files/model_dissimilarity_data_within_dilutions.csv", header = T)

# Only sequential data
model_beta_seq <- model_beta %>%
  filter(sequential == "YES")

beta_model_plot <- ggplot(model_beta_seq, aes(x =dilution_rate, y = sim), ) +
  geom_boxplot(aes(fill = dilution_rate), color = "black", alpha = 0.9, width = 0.5) +
  # geom_jitter(aes(color = dilution_rate, shape = model), width = 0.2, size = 2, alpha = 0.9) +
  scale_color_manual(values=c("#2B5B6C", "#E34F33", "#FFC87E")) +
  scale_fill_manual(values=c("#2B5B6C", "#E34F33", "#FFC87E")) +
  xlab("Dilution Rate") +
  ylab("Similarity Metric") +
  pretty.theme() +
  ggtitle("F. Community Stability - Modeled") +
  # ylim(0,1) +
  facet_wrap(.~factor(distribution, levels=c("uniform", "beta"))) +  # this changes the order of the facet wrap
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 14, colour = "black"),
        legend.position = "none") +
  ylim(0,1)
beta_model_plot

## Combining plots----
fig4 <- (bray_plot1 | bray_plot2) /
  (alpha_boxplot_ALL | dis_boxplot_ALL) /
  (alpha_model_plot | beta_model_plot)
fig4
