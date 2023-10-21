# © This script is written by Anabel von Jackowski @avonjackowski
#####################################
#setup
#####################################
load("./data_import_v4v5.RData")
source("./Source_Scripts/rarefaction_curve.R")
source("./Source_Scripts/multiplot.R")
source("./Source_Scripts/miseqR.R")

library(phyloseq)
library(ggplot2)
library(igraph)
library(vegan)
library(grid)
library(DESeq2)
library(plyr)
library(dplyr)
#####################################
#Plot Alpha diversity indeces - to estimate sampling biases
#####################################
richness <- estimate_richness(all_data, split = TRUE, measures = NULL)

#Set parameters as factors and define order of plotting
all_data@sam_data$Station <- as.factor(all_data@sam_data$Station)
all_data@sam_data$Depth <- ordered(all_data@sam_data$Depth, levels = c("Surface","200 m", "400 m", "600 m" ,"Bathypelagic","Bottom Water")) #water column
all_data@sam_data$Location <- factor(all_data@sam_data$Location, levels=c("Northern Seamount", "Central Seamount","Karasik Seamount", "Reference Site"))

all_data@sam_data$Location <- factor(all_data@sam_data$Location, levels=c("Polaris Vent", "Northern Seamount", "Central Seamount","Seamount Saddle","Karasik Seamount", "Southern Slope","S-Reference"))
all_data@sam_data$Layer <- factor(all_data@sam_data$Layer, levels=c("Surface (0-1cm)","Subsurface (1-5cm)","Deep Layer (14-16cm)"))

#Depth for Water
#Layer for Sediment
plot_richness(all_data, x="Location", color = "Depth",
              measures = c("Observed", "Chao1", "Shannon","InvSimpson"),
              shape = "Depth") +
  theme_bw()+
  geom_point(size=2)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
        strip.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 12)) +
  xlab("Location at Langseth Ridge")+
  scale_shape_manual(values = c(15,17,18,19,24,25))+ #Water Column
  # scale_shape_manual(values = c(15,17,18,24,25)) +
  scale_color_manual(values = c("black","black","black","black","black","black","black","black","black","black"))
  # scale_shape_manual(values = c(19,18,15))+ #Sediment
  # scale_color_manual(values = c("#0B775E","#46ACC8","#46ACC8","#46ACC8","#46ACC8","#ADADAD","#ADADAD","black"))
  # scale_color_manual(values = c("#0B775E","#46ACC8","#ADADAD","black"))

ggsave("./Figures/v4v5_Richness.pdf", width = 8.3, height = 4)
ggsave("./Figures/v4v5_Richness.pdf", width = 7.3, height = 4)
write.csv(richness, "./v3v5_Richness.csv")

#for Statistics on alpha diversity indices see Statistics.R