# Figure 2 in Meg/Ashley SES Manuscript
# Environmental parameters from chemostats (pH, CO2, O2, Dissolved O2)

# UPLOAD DATA----
path = "/Users/abulseco/Dropbox/MBL_Postdoc/Experiments/2018_Pilot_SESProject/Microbial_Data/Downloads_from_VAMPS/JJV_MEG_Bv4v5" 
setwd("~/Dropbox/MBL_Postdoc/Experiments/2018_Pilot_SESProject/Microbial_Data/Downloads_from_VAMPS/JJV_MEG_Bv4v5")

# Load necessary libraries
library("ggplot2")
library(patchwork)

# This is a slightly different theme function from others 
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

theme.exp1 <- function(){
  theme_bw()+
    theme(axis.text.x=element_text(size=18, vjust=1, hjust=1, color = "black"),
          axis.text.y=element_text(size=18, color = "black"),
          axis.title.x=element_text(size=14, face="plain"),             
          axis.title.y=element_text(size=14, face="plain"),             
          panel.grid.major.x=element_blank(),                                          
          panel.grid.minor.x=element_blank(),
          panel.grid.minor.y=element_blank(),
          panel.grid.major.y=element_blank(),  
          plot.title = element_text(size=20, vjust=1, hjust=0.5),
          legend.text = element_text(size=12, face="italic"),          
          legend.title = element_blank(), 
          axis.title = element_blank())
}

# Read in the data (dissolved O2 had a different output)
monitor_var <- read.csv("monitor_param_LONG.csv", header = T)
str(monitor_var)
dissolved_oxy <- read.csv("dissolved_oxygen.csv", header = T)

# pH Plot
#############################################################################
#############################################################################
ph_plot <- ggplot(monitor_var, group = chemostat_ID, aes(x = time, y = ph)) +
  annotate("rect", fill = "#2B5B6C", alpha = 0.4, 
           xmin = 0, xmax = 15.71,
           ymin = -Inf, ymax = Inf) +
  annotate("rect", fill = "#E34F33", alpha = 0.4, 
           xmin = 15.71, xmax = 20.79,
           ymin = -Inf, ymax = Inf) +
  annotate("rect", fill = "#FFC87E", alpha = 0.4, 
           xmin = 20.79, xmax = 25,
           ymin = -Inf, ymax = Inf) +
  # geom_vline(xintercept = 15.71, color = "black", linetype = "dashed") +
  # geom_vline(xintercept = 20.79, color = "black", linetype = "dashed") +
  # scale_fill_manual(values = c("#3399CC", "#003399", linetype = "dashed")) +
  scale_color_manual(values = c("black", "gray45", linetype = "dashed")) +
  # facet_wrap(chemostat_ID~.) +
  geom_line(size = 0.75, aes(linetype = chemostat_ID)) +
  pretty.theme() +
  theme(strip.background = element_blank(), # Overriding some of the original theme
        strip.text.x = element_blank(), 
        axis.title.x=element_blank(), 
        axis.title.y=element_blank(), 
        legend.position = "none") 
ph_plot

# CO2 Plot
#############################################################################
#############################################################################
CO2_plot <- ggplot(monitor_var, group = chemostat_ID, aes(x = time, y = carbon_dio_perc)) +
  annotate("rect", fill = "#2B5B6C", alpha = 0.4, 
           xmin = 0, xmax = 15.71,
           ymin = -Inf, ymax = Inf) +
  annotate("rect", fill = "#E34F33", alpha = 0.4, 
           xmin = 15.71, xmax = 20.79,
           ymin = -Inf, ymax = Inf) +
  annotate("rect", fill = "#FFC87E", alpha = 0.4, 
           xmin = 20.79, xmax = 25,
           ymin = -Inf, ymax = Inf) +
  # geom_vline(xintercept = 15.71, color = "black", linetype = "dashed") +
  # geom_vline(xintercept = 20.79, color = "black", linetype = "dashed") +
  # scale_fill_manual(values = c("#3399CC", "#003399", linetype = "dashed")) +
  scale_color_manual(values = c("black", "gray45", linetype = "dashed")) +
  # facet_wrap(chemostat_ID~.) +
  geom_line(size = 0.75, aes(linetype = chemostat_ID)) +
  pretty.theme() +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(), 
        axis.title.x=element_blank(), 
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position = "none")  
CO2_plot

# O2 Plot
#############################################################################
#############################################################################
O2_plot <- ggplot(monitor_var, group = chemostat_ID, aes(x = time, y = oxygen_perc)) +
  annotate("rect", fill = "#2B5B6C", alpha = 0.4, 
           xmin = 0, xmax = 15.71,
           ymin = -Inf, ymax = Inf) +
  annotate("rect", fill = "#E34F33", alpha = 0.4, 
           xmin = 15.71, xmax = 20.79,
           ymin = -Inf, ymax = Inf) +
  annotate("rect", fill = "#FFC87E", alpha = 0.4, 
           xmin = 20.79, xmax = 25,
           ymin = -Inf, ymax = Inf) +
  # geom_vline(xintercept = 15.71, color = "black", linetype = "dashed") +
  # geom_vline(xintercept = 20.79, color = "black", linetype = "dashed") +
  # scale_fill_manual(values = c("#3399CC", "#003399", linetype = "dashed")) +
  scale_color_manual(values = c("black", "gray45", linetype = "dashed")) +
  # facet_wrap(chemostat_ID~.) +
  geom_line(size = 0.75, aes(linetype = chemostat_ID)) +
  pretty.theme() +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(), 
        axis.title.x=element_blank(), 
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position = "none") 
O2_plot

# Dissolved O2 Plot
#############################################################################
#############################################################################
diss_O2_plot <- ggplot(dissolved_oxy, group = chemostat_ID, aes(x = time, y = diss_o2)) +
  annotate("rect", fill = "#2B5B6C", alpha = 0.4, 
           xmin = 0, xmax = 15.71,
           ymin = -Inf, ymax = Inf) +
  annotate("rect", fill = "#E34F33", alpha = 0.4, 
           xmin = 15.71, xmax = 20.79,
           ymin = -Inf, ymax = Inf) +
  annotate("rect", fill = "#FFC87E", alpha = 0.4, 
           xmin = 20.79, xmax = 25,
           ymin = -Inf, ymax = Inf) +
  # geom_vline(xintercept = 15.71, color = "black", linetype = "dashed") +
  # geom_vline(xintercept = 20.79, color = "black", linetype = "dashed") +
  # scale_fill_manual(values = c("#3399CC", "#003399", linetype = "dashed")) +
  scale_color_manual(values = c("black", "gray45", linetype = "dashed")) +
  # facet_wrap(chemostat_ID~.) +
  geom_line(size = 0.75, aes(linetype = chemostat_ID)) +
  pretty.theme() +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(), 
        axis.title.x=element_blank(), 
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position = "none") 
diss_O2_plot

# Then using patchwork, connect all the plots
# After exporting to PDF, use InkScape to add in timeline along the top
# axes labels along the left of each plot, and legend for MC1 and MC2.

pdf("Fig2.pdf")
CO2_plot + O2_plot + diss_O2_plot + ph_plot + plot_layout(ncol = 1)
dev.off()

# CREATING A TIMELINE PLOT
# FROM THIS AWESOME BLOG POST: https://benalexkeen.com/creating-a-timeline-graphic-using-r-and-ggplot2/
#############################################################################
#############################################################################
# (But it is not being used in the figure)
df <- read.csv("sampling_timeline_final.csv")

status_levels <- c("A", "B", "C")
status_colors <- c("#2B5B6C", "#E34F33", "#FFC87E")

df$dilution_rate <- factor(df$dilution_rate, levels=status_levels, ordered=TRUE)

# PLOT
timeline_plot<-ggplot(df,aes(x=time,y=0, col=dilution_rate))
# timeline_plot<-timeline_plot+labs(col="dilution_rate")
timeline_plot<-timeline_plot+scale_color_manual(values=status_colors, labels=status_levels, drop = FALSE)
timeline_plot<-timeline_plot+theme_classic()

# Plot horizontal black line for timeline
timeline_plot<-timeline_plot+geom_hline(yintercept=0, 
                                        color = "black", size=0.3)

# Plot scatter points at zero and date
timeline_plot<-timeline_plot+geom_point(aes(y=0), size=2)

# Don't show axes, appropriately position legend
timeline_plot<-timeline_plot+theme(axis.line.y=element_blank(),
                                   axis.text.y=element_blank(),
                                   axis.title.x=element_blank(),
                                   axis.title.y=element_blank(),
                                   axis.ticks.y=element_blank(),
                                   axis.text.x =element_blank(),
                                   axis.ticks.x =element_blank(),
                                   axis.line.x =element_blank(),
                                   legend.title = element_blank(),                              
                                   legend.position="none"
) + scale_x_continuous(limits=c(0, 25)) +
  scale_y_continuous(limits=c(-1,1))

timeline_plot
