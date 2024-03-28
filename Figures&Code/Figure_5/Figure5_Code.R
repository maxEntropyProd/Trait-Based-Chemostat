# Redoing diversity plots for Joe
# Need to get formatting consistent for publishing
# Original code comes from file: 2020_04_23_SES_Code.R
# Refer to "Final_Phyloseq_Feb2024.R" for Phyloseq object generation

# UPLOAD DATA----
path = "/Users/abulseco/Dropbox/MBL_Postdoc/Experiments/2018_Pilot_SESProject/Microbial_Data/Downloads_from_VAMPS/JJV_MEG_Bv4v5" 
setwd("~/Dropbox/MBL_Postdoc/Experiments/2018_Pilot_SESProject/Microbial_Data/Downloads_from_VAMPS/JJV_MEG_Bv4v5")

# Load necessary libraries
library("ggplot2")
library(phyloseq)
library(ggrepel)
library(patchwork)

# For plot consistency
pretty.theme <- function(){
  theme_bw()+
    theme(axis.text.x=element_text(size=14, color = "black", angle = 0, hjust = 1),
          axis.text.y=element_text(size=14, color = "black"),
          axis.title.x=element_text(size=20, color = "black"),             
          axis.title.y=element_text(size=20, color = "black"),             
          panel.grid.major.x=element_blank(),                                          
          panel.grid.minor.x=element_blank(),
          panel.grid.minor.y=element_blank(),
          panel.grid.major.y=element_blank(),  
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), units = , "cm"),
          plot.title = element_text(size=20),
          legend.text = element_text(size=12, face="italic"),          
          legend.title = element_blank(),                              
          legend.position="none")
}

# Panels A&B: PCoA plot for MC1
#################################################################################
#################################################################################

# Interested in subsetting by chemostat ID
# Reminder: per_nopond is where relative abundance has been calculated 
# And pond samples have been removed
per_MC1 <- subset_samples(per_nopond, chemostat_ID=="MC1")
per_MC2 <- subset_samples(per_nopond, chemostat_ID=="MC2")

# Bray curtis for MC1 specifically
BC_distance1 <- phyloseq::distance(per_MC1, "bray")
bcOrd1 <- ordinate(per_MC1, "PCoA", BC_distance1)
plot_scree(bcOrd1)
bray_plot1 <- plot_ordination(per_MC1, bcOrd1) +
  geom_point(aes(fill = dilution_rate), size = 8, alpha = 1, pch = 21, color = "black") +
  # geom_polygon(aes(fill = dilution_rate), alpha = 0.8) +
  geom_text_repel(aes(label= total_days_after), 
                  color = "black",
                  # vjust = 0.02, nudge_y = 0.01,
                  direction = "both",
                  size = 6,
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
  geom_point(aes(fill = dilution_rate), size = 8, alpha = 1, pch = 21, color = "black") +
  # geom_polygon(aes(fill = dilution_rate), alpha = 0.8) +
  geom_text_repel(aes(label= total_days_after), 
                  color = "black",
                  # vjust = 0.02, nudge_y = 0.01,
                  direction = "both",
                  size = 6,
                  min.segment.length = Inf) +
  scale_fill_manual(values=c("#2B5B6C", "#E34F33", "#FFC87E","black")) +
  scale_color_manual(values=c("#2B5B6C", "#E34F33", "#FFC87E","black")) +
  pretty.theme() +
  theme(legend.title = element_blank()) +
  ggtitle("B. PCoA for MC2 - Observed")
bray_plot2

# Panel C
#################################################################################
#################################################################################

# Then this is for alpha diversity
estimate_richness(phy.f.nopond, measures = "Shannon")
estimate_richness(phy.f.nopond.rare, measures = "Shannon")

alpha_boxp_plot_ALL <- plot_richness(phy.f.nopond.rare, x="dilution_rate", measures=c("Shannon")) +
  geom_boxplot(aes(fill = dilution_rate), color = "black") +
  # geom_jitter(aes(color = dilution_rate, shape = chemostat_ID), width = 0.4, size = 4) +
  # scale_color_manual(values=c("#2B5B6C", "#E34F33", "#FFC87E")) +
  scale_fill_manual(values=c("#2B5B6C", "#E34F33", "#FFC87E")) +
  xlab("Dilution Rate") +
  ylab("Shannon Diversity") +
  pretty.theme() +
  ggtitle("C. Alpha Diversity - Observed") +
  theme(strip.background = element_blank(), 
        strip.text = element_blank()) +
  ylim(0,7) 
alpha_boxp_plot_ALL

# Panel D
#################################################################################
#################################################################################

# Load data
dis_box_data_ALL<- read.csv("BRAY_DISS_updated.csv", header = TRUE)
dis_boxplot_ALL <- ggplot(dis_box_data_ALL, aes(x = DILUTION_RATE, y = SIM), ) +
  geom_boxplot(aes(fill = DILUTION_RATE), color = "black", alpha = 0.9) +
  # geom_jitter(aes(color = DILUTION_RATE, shape = CHEMOSTAT_ID), width = 0.4, size = 4) +
  # scale_color_manual(values=c("#2B5B6C", "#E34F33", "#FFC87E")) +
  scale_fill_manual(values=c("#2B5B6C", "#E34F33", "#FFC87E")) +
  xlab("Dilution Rate") +
  ylab("Similarity Metric") +
  pretty.theme() +
  ggtitle("D. Community Stability - Observed") +
  ylim(0,1) 
dis_boxplot_ALL

# Panel E - Modeled Data (Beta)
#################################################################################
#################################################################################
# alpha diversity based on modeled data
beta_alpha <- read.csv("betaALPHA_forboxplot.csv", header = TRUE)

beta_alpha_boxplot <- ggplot(beta_alpha, aes(x = DILUTION_RATE, y = SHANNON), ) +
  geom_boxplot(aes(fill = DILUTION_RATE)) +
  # geom_jitter(aes(fill = DILUTION_RATE), width = 0.4, size = 4, pch = 21) +
  scale_fill_manual(values=c("#2B5B6C", "#E34F33", "#FFC87E","black")) +
  xlab("Dilution Rate") +
  ylab("Shannon Diversity") +
  pretty.theme() +
  ylim(0,7) + 
  labs(title = "E. Alpha Diversity - Beta Model")
beta_alpha_boxplot

# Panel F - Modeled Data (Beta)
#################################################################################
#################################################################################
beta_beta <- read.csv("betaBray_forboxplot_sequential-NEW.csv", header = TRUE)

beta_beta_boxplot <- ggplot(beta_beta, aes(x = DILUTION_RATE, y = SIM), ) +
  geom_boxplot(aes(fill = DILUTION_RATE)) +
  # geom_jitter(aes(fill = DILUTION_RATE), width = 0.4, size = 4, pch = 21) +
  scale_fill_manual(values=c("#2B5B6C", "#E34F33", "#FFC87E","black")) +
  xlab("Dilution Rate") +
  ylab("Similarity Metric") +
  pretty.theme() +
  ggtitle("F. Community Stability - Beta Model") +
  ylim(0,1)
beta_beta_boxplot

# Supplementary Model Data (Uniform)
uniform_alpha <- read.csv("uniformALPHA_forboxplot.csv", header = TRUE)

uniform_alpha_boxplot <- ggplot(uniform_alpha, aes(x = DILUTION_RATE, y = SHANNON), ) +
  geom_boxplot(aes(fill = DILUTION_RATE)) +
  # geom_jitter(aes(fill = DILUTION_RATE), width = 0.4, size = 4, pch = 21) +
  scale_fill_manual(values=c("#2B5B6C", "#E34F33", "#FFC87E","black")) +
  xlab("Dilution Rate") +
  ylab("Shannon Diversity") +
  pretty.theme() +
  ylim(0,7) +
  labs(title = "A. Alpha Diversity - Uniform Model")
uniform_alpha_boxplot

# Plot (after manipulating the data)
uniform_beta <- read.csv("uniformBray_forboxplot_sequential-NEW.csv", header = TRUE)

uniform_beta_boxplot <- ggplot(uniform_beta, aes(x = DILUTION_RATE, y = SIM), ) +
  geom_boxplot(aes(fill = DILUTION_RATE)) +
  # geom_jitter(aes(fill = DILUTION_RATE), width = 0.4, size = 4, pch = 21) +
  scale_fill_manual(values=c("#2B5B6C", "#E34F33", "#FFC87E","black")) +
  xlab("Dilution Rate") +
  ylab("Similarity Metric") +
  pretty.theme() +
  ggtitle("B. Community Stability - Uniform Model")
uniform_beta_boxplot

# Exporting Plots

# Exporting plots as PDFs for Joe
pdf("Fig.5.pdf", width = 15, height = 15)
(bray_plot1 | bray_plot2) /
  (alpha_boxp_plot_ALL | dis_boxplot_ALL) /
  (beta_alpha_boxplot | beta_beta_boxplot)
dev.off()

# Exporting each individually
#################################################################################
#################################################################################
pdf("Fig.5 - Individual.pdf")
bray_plot1
bray_plot2
alpha_boxp_plot_ALL
dis_boxplot_ALL
beta_alpha_boxplot
beta_beta_boxplot
dev.off()
