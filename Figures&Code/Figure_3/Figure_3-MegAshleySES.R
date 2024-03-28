# Figure 2 in Meg/Ashley SES Manuscript
# Environmental parameters from chemostats

# UPLOAD DATA----
path = "/Users/abulseco/Dropbox/MBL_Postdoc/Experiments/2018_Pilot_SESProject/Microbial_Data/Downloads_from_VAMPS/JJV_MEG_Bv4v5" 
setwd("~/Dropbox/MBL_Postdoc/Experiments/2018_Pilot_SESProject/Microbial_Data/Downloads_from_VAMPS/JJV_MEG_Bv4v5")

# Load necessary libraries
library("ggplot2")
library(patchwork)

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

nuts_data <- read.csv("new_nutrients.csv", header = TRUE)

nitrate_plot <- ggplot(nuts_data, aes(x = days_after, y = nitrate_conc)) +
  annotate("rect", fill = "#2B5B6C", alpha = 0.4, 
           xmin = 0, xmax = 15.71,
           ymin = -Inf, ymax = Inf) +
  annotate("rect", fill = "#E34F33", alpha = 0.4, 
           xmin = 15.71, xmax = 20.79,
           ymin = -Inf, ymax = Inf) +
  annotate("rect", fill = "#FFC87E", alpha = 0.4, 
           xmin = 20.79, xmax = 25,
           ymin = -Inf, ymax = Inf) +
  geom_point(aes(shape = chemostat_ID, color = chemostat_ID), size = 2) +
  geom_line(aes(color = chemostat_ID, linetype = chemostat_ID)) +
  labs(y = expression (""), x = "Days After Start") +
  scale_color_manual(values = c("black", "black")) +
  pretty.theme()
nitrate_plot 

phosphate_plot <- ggplot(nuts_data, aes(x = days_after, y = phos_conc)) +
  annotate("rect", fill = "#2B5B6C", alpha = 0.4, 
           xmin = 0, xmax = 15.71,
           ymin = -Inf, ymax = Inf) +
  annotate("rect", fill = "#E34F33", alpha = 0.4, 
           xmin = 15.71, xmax = 20.79,
           ymin = -Inf, ymax = Inf) +
  annotate("rect", fill = "#FFC87E", alpha = 0.4, 
           xmin = 20.79, xmax = 25,
           ymin = -Inf, ymax = Inf) +
  geom_point(aes(shape = chemostat_ID, color = chemostat_ID), size = 2) +
  geom_line(aes(color = chemostat_ID, linetype = chemostat_ID)) +
  labs(y = expression (""), x = "Days After Start") +
  scale_color_manual(values = c("black", "black")) +
  pretty.theme()
phosphate_plot

ammonium_plot <- ggplot(nuts_data, aes(x = days_after, y = ammonium_con)) +
  annotate("rect", fill = "#2B5B6C", alpha = 0.4, 
           xmin = 0, xmax = 15.71,
           ymin = -Inf, ymax = Inf) +
  annotate("rect", fill = "#E34F33", alpha = 0.4, 
           xmin = 15.71, xmax = 20.79,
           ymin = -Inf, ymax = Inf) +
  annotate("rect", fill = "#FFC87E", alpha = 0.4, 
           xmin = 20.79, xmax = 25,
           ymin = -Inf, ymax = Inf) +
  geom_point(aes(shape = chemostat_ID, color = chemostat_ID), size = 2) +
  geom_line(aes(color = chemostat_ID, linetype = chemostat_ID)) +
  labs(y = expression (""), x = "Days After Start") +
  scale_color_manual(values = c("black", "black")) +
  pretty.theme() 
ammonium_plot

DOC_plot <- ggplot(nuts_data, aes(x = days_after, y = doc_conc_mmol)) +
  annotate("rect", fill = "#2B5B6C", alpha = 0.4, 
           xmin = 0, xmax = 15.71,
           ymin = -Inf, ymax = Inf) +
  annotate("rect", fill = "#E34F33", alpha = 0.4, 
           xmin = 15.71, xmax = 20.79,
           ymin = -Inf, ymax = Inf) +
  annotate("rect", fill = "#FFC87E", alpha = 0.4, 
           xmin = 20.79, xmax = 25,
           ymin = -Inf, ymax = Inf) +
  geom_point(aes(shape = chemostat_ID, color = chemostat_ID), size = 2) +
  geom_line(aes(color = chemostat_ID, linetype = chemostat_ID)) +
  labs(y = expression (""), x = "Days After Start") +
  scale_color_manual(values = c("black", "black")) +
  pretty.theme() 
DOC_plot

